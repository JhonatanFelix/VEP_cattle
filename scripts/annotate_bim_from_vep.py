#!/usr/bin/env python3

import argparse
import csv
import gzip
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


DNA_COMPLEMENT = str.maketrans({"A": "T", "T": "A", "C": "G", "G": "C"})
IMPACT_RANK = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1}
FEATURE_RANK = {"Transcript": 3, "RegulatoryFeature": 2, "MotifFeature": 1}
PICK_VALUES = {"1", "YES", "Y", "yes", "y"}


@dataclass
class VcfRecord:
    chrom_raw: str
    chrom_norm: str
    pos: int
    variant_id: str
    ref: str
    alts: list[str]
    filt: str
    info: str


class VcfStream:
    def __init__(self, vcf_path: Path):
        self.vcf_path = vcf_path
        self.proc = subprocess.Popen(
            ["gzip", "-dc", str(vcf_path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
            errors="replace",
            bufsize=1024 * 1024,
        )
        if self.proc.stdout is None:
            raise RuntimeError("Could not open VCF stream.")
        self.stdout = self.proc.stdout
        self.csq_fields: list[str] | None = None
        self.current_line: str | None = None
        self.current_key: tuple[int, str, int] | None = None
        self.previous_key: tuple[int, str, int] | None = None
        self.lines_read = 0
        self._prime()

    def _prime(self) -> None:
        for line in self.stdout:
            if line.startswith("##INFO=<ID=CSQ"):
                self.csq_fields = parse_csq_header(line)
                continue
            if line.startswith("#"):
                continue
            self.current_line = line
            self.current_key = parse_line_key(line)
            self.previous_key = self.current_key
            self.lines_read += 1
            return
        self.current_line = None
        self.current_key = None

    def advance(self) -> None:
        for line in self.stdout:
            self.current_line = line
            self.current_key = parse_line_key(line)
            if self.previous_key is not None and self.current_key < self.previous_key:
                raise RuntimeError(
                    "VCF records are not sorted by chromosome and position, so a streaming merge "
                    "cannot be done safely."
                )
            self.previous_key = self.current_key
            self.lines_read += 1
            return
        self.current_line = None
        self.current_key = None

    def collect_current_position_records(self) -> tuple[tuple[int, str, int], list[VcfRecord]]:
        if self.current_key is None or self.current_line is None:
            raise RuntimeError("No current VCF record to collect.")
        key = self.current_key
        records: list[VcfRecord] = []
        while self.current_key == key and self.current_line is not None:
            records.append(parse_vcf_record(self.current_line))
            self.advance()
        return key, records

    def close(self) -> None:
        if self.stdout:
            self.stdout.close()
        stderr_output = ""
        if self.proc.stderr is not None:
            stderr_output = self.proc.stderr.read()
            self.proc.stderr.close()
        return_code = self.proc.wait()
        if return_code not in (0, -13, 141):
            raise RuntimeError(
                f"gzip exited with code {return_code} while reading {self.vcf_path}.\n{stderr_output}"
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Annotate a PLINK BIM file with VEP annotations from a gzipped VCF while "
            "keeping the original BIM SNP order."
        )
    )
    parser.add_argument("--bim", required=True, type=Path, help="Input BIM file.")
    parser.add_argument("--vcf", required=True, type=Path, help="Input VEP VCF .gz file.")
    parser.add_argument("--output", required=True, type=Path, help="Output TSV or TSV.GZ file.")
    parser.add_argument(
        "--summary",
        required=False,
        type=Path,
        help="Optional summary text file with merge statistics.",
    )
    return parser.parse_args()


def normalize_chrom(chrom: str) -> str:
    chrom = chrom.strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    chrom = chrom.upper()
    if chrom == "M":
        chrom = "MT"
    return chrom


def chrom_rank(chrom: str) -> int:
    if chrom.isdigit():
        return int(chrom)
    if chrom == "X":
        return 100
    if chrom == "Y":
        return 101
    if chrom == "MT":
        return 102
    return 1000


def make_key(chrom: str, pos: int) -> tuple[int, str, int]:
    chrom_norm = normalize_chrom(chrom)
    return chrom_rank(chrom_norm), chrom_norm, pos


def parse_line_key(line: str) -> tuple[int, str, int]:
    first_tab = line.find("\t")
    second_tab = line.find("\t", first_tab + 1)
    chrom = line[:first_tab]
    pos = int(line[first_tab + 1 : second_tab])
    return make_key(chrom, pos)


def parse_vcf_record(line: str) -> VcfRecord:
    fields = line.rstrip("\n").split("\t", 8)
    chrom_raw, pos, variant_id, ref, alt, _qual, filt, info = fields[:8]
    return VcfRecord(
        chrom_raw=chrom_raw,
        chrom_norm=normalize_chrom(chrom_raw),
        pos=int(pos),
        variant_id=variant_id,
        ref=ref.upper(),
        alts=[entry.upper() for entry in alt.split(",")],
        filt=filt,
        info=info,
    )


def parse_csq_header(line: str) -> list[str]:
    marker = "Format: "
    if marker not in line:
        raise RuntimeError("Found CSQ header, but could not parse the declared field order.")
    csq_layout = line.split(marker, 1)[1].rsplit('"', 1)[0]
    return csq_layout.split("|")


def extract_info_value(info: str, key: str) -> str:
    prefix = f"{key}="
    for field in info.split(";"):
        if field.startswith(prefix):
            return field[len(prefix) :]
    return ""


def complement_base(base: str) -> str:
    if len(base) != 1:
        return ""
    return base.translate(DNA_COMPLEMENT)


def allele_match_status(bim_a1: str, bim_a2: str, ref: str, alt: str) -> str | None:
    if len(ref) != 1 or len(alt) != 1:
        return None
    bim_set = {bim_a1.upper(), bim_a2.upper()}
    vcf_set = {ref.upper(), alt.upper()}
    if bim_set == vcf_set:
        return "exact_allele_match"
    ref_comp = complement_base(ref.upper())
    alt_comp = complement_base(alt.upper())
    if ref_comp and alt_comp and bim_set == {ref_comp, alt_comp}:
        return "strand_flip_match"
    return None


def select_candidate(records: list[VcfRecord], bim_a1: str, bim_a2: str) -> dict | None:
    exact_candidates = []
    strand_candidates = []

    for record in records:
        for alt_index, alt in enumerate(record.alts, start=1):
            match_status = allele_match_status(bim_a1, bim_a2, record.ref, alt)
            candidate = {
                "record": record,
                "alt": alt,
                "alt_index": alt_index,
                "position_record_count": len(records),
            }
            if match_status == "exact_allele_match":
                candidate["match_status"] = match_status
                exact_candidates.append(candidate)
            elif match_status == "strand_flip_match":
                candidate["match_status"] = match_status
                strand_candidates.append(candidate)

    if exact_candidates:
        return exact_candidates[0]
    if strand_candidates:
        return strand_candidates[0]

    if len(records) == 1 and len(records[0].alts) == 1:
        return {
            "record": records[0],
            "alt": records[0].alts[0],
            "alt_index": 1,
            "match_status": "position_only_ref_alt_mismatch",
            "position_record_count": 1,
        }
    if len(records) == 1:
        return {
            "record": records[0],
            "alt": records[0].alts[0],
            "alt_index": 1,
            "match_status": "position_only_multiallelic_unresolved",
            "position_record_count": 1,
        }
    return None


def parse_csq_entries(
    info: str, csq_fields: list[str], matched_alt: str, matched_alt_index: int
) -> tuple[list[dict[str, str]], list[str]]:
    raw_csq = extract_info_value(info, "CSQ")
    if not raw_csq:
        return [], []

    allele_num_key = "ALLELE_NUM" if "ALLELE_NUM" in csq_fields else None
    entries: list[dict[str, str]] = []
    raw_entries: list[str] = []

    for raw_entry in raw_csq.split(","):
        values = raw_entry.split("|")
        if len(values) < len(csq_fields):
            values.extend([""] * (len(csq_fields) - len(values)))
        elif len(values) > len(csq_fields):
            values = values[: len(csq_fields)]
        entry = dict(zip(csq_fields, values))
        allele_matches = entry.get("Allele", "").upper() == matched_alt.upper()
        allele_num_matches = allele_num_key and entry.get(allele_num_key, "") == str(matched_alt_index)
        if allele_matches or allele_num_matches:
            entries.append(entry)
            raw_entries.append(raw_entry)

    if entries:
        return entries, raw_entries

    fallback_entries: list[dict[str, str]] = []
    fallback_raw_entries: list[str] = []
    for raw_entry in raw_csq.split(","):
        values = raw_entry.split("|")
        if len(values) < len(csq_fields):
            values.extend([""] * (len(csq_fields) - len(values)))
        elif len(values) > len(csq_fields):
            values = values[: len(csq_fields)]
        fallback_entries.append(dict(zip(csq_fields, values)))
        fallback_raw_entries.append(raw_entry)
    return fallback_entries, fallback_raw_entries


def best_csq_entry(entries: list[dict[str, str]]) -> dict[str, str]:
    def score(entry: dict[str, str]) -> tuple[int, int, int, int, int, int, int, int, int, int]:
        pick = 1 if entry.get("PICK", "") in PICK_VALUES or entry.get("PICK_allele", "") in PICK_VALUES else 0
        canonical = 1 if entry.get("CANONICAL", "") == "YES" else 0
        impact = IMPACT_RANK.get(entry.get("IMPACT", ""), 0)
        feature = FEATURE_RANK.get(entry.get("Feature_type", ""), 0)
        protein_coding = 1 if entry.get("BIOTYPE", "") == "protein_coding" else 0
        non_intergenic = 0 if "intergenic_variant" in entry.get("Consequence", "") else 1
        has_hgvsp = 1 if entry.get("HGVSp", "") else 0
        has_hgvsc = 1 if entry.get("HGVSc", "") else 0
        has_predictions = 1 if entry.get("SIFT", "") or entry.get("PolyPhen", "") else 0
        has_gene = 1 if entry.get("SYMBOL", "") or entry.get("Gene", "") else 0
        return (
            pick,
            canonical,
            impact,
            feature,
            protein_coding,
            non_intergenic,
            has_hgvsp,
            has_hgvsc,
            has_predictions,
            has_gene,
        )

    return max(entries, key=score)


def open_output(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "wt", newline="")
    return open(path, "w", newline="")


def annotate_bim(args: argparse.Namespace) -> dict[str, int]:
    vcf_stream = VcfStream(args.vcf)
    try:
        if not vcf_stream.csq_fields:
            raise RuntimeError("No VEP CSQ header was found in the VCF.")

        stats = {
            "total_bim_rows": 0,
            "exact_allele_match": 0,
            "strand_flip_match": 0,
            "position_only_ref_alt_mismatch": 0,
            "position_only_multiallelic_unresolved": 0,
            "position_multiple_records_no_allele_match": 0,
            "no_vcf_position": 0,
            "with_csq_annotation": 0,
        }

        output_fields = [
            "bim_chr",
            "bim_snp_id",
            "bim_cm",
            "bim_bp",
            "bim_allele1",
            "bim_allele2",
            "match_status",
            "position_record_count",
            "vcf_chr",
            "vcf_pos",
            "vcf_id",
            "vcf_ref",
            "vcf_alt",
            "vcf_filter",
            "vep_annotation_count",
            "vep_matched_csq",
        ] + [f"vep_best_{field}" for field in vcf_stream.csq_fields]

        cached_key: tuple[int, str, int] | None = None
        cached_records: list[VcfRecord] = []
        previous_bim_key: tuple[int, str, int] | None = None

        with args.bim.open("r", encoding="utf-8", errors="replace") as bim_handle, open_output(
            args.output
        ) as output_handle:
            writer = csv.DictWriter(output_handle, fieldnames=output_fields, delimiter="\t")
            writer.writeheader()

            for line_number, line in enumerate(bim_handle, start=1):
                stripped = line.rstrip("\n")
                if not stripped:
                    continue
                fields = stripped.split()
                if len(fields) < 6:
                    raise RuntimeError(
                        f"BIM line {line_number} does not have at least 6 columns: {stripped}"
                    )

                bim_chr, bim_snp_id, bim_cm, bim_bp, bim_a1, bim_a2 = fields[:6]
                bim_pos = int(bim_bp)
                bim_key = make_key(bim_chr, bim_pos)
                if previous_bim_key is not None and bim_key < previous_bim_key:
                    raise RuntimeError(
                        "The BIM file is not sorted by chromosome and position, so the streaming "
                        "merge cannot preserve order safely."
                    )
                previous_bim_key = bim_key
                stats["total_bim_rows"] += 1

                if cached_key != bim_key:
                    while vcf_stream.current_key is not None and vcf_stream.current_key < bim_key:
                        vcf_stream.advance()

                    if vcf_stream.current_key == bim_key:
                        cached_key, cached_records = vcf_stream.collect_current_position_records()
                    else:
                        cached_key = bim_key
                        cached_records = []

                row = {
                    "bim_chr": bim_chr,
                    "bim_snp_id": bim_snp_id,
                    "bim_cm": bim_cm,
                    "bim_bp": bim_bp,
                    "bim_allele1": bim_a1,
                    "bim_allele2": bim_a2,
                    "match_status": "",
                    "position_record_count": 0,
                    "vcf_chr": "",
                    "vcf_pos": "",
                    "vcf_id": "",
                    "vcf_ref": "",
                    "vcf_alt": "",
                    "vcf_filter": "",
                    "vep_annotation_count": 0,
                    "vep_matched_csq": "",
                }
                for field in vcf_stream.csq_fields:
                    row[f"vep_best_{field}"] = ""

                if not cached_records:
                    row["match_status"] = "no_vcf_position"
                    stats["no_vcf_position"] += 1
                    writer.writerow(row)
                    continue

                candidate = select_candidate(cached_records, bim_a1, bim_a2)
                if candidate is None:
                    row["match_status"] = "position_multiple_records_no_allele_match"
                    row["position_record_count"] = len(cached_records)
                    stats["position_multiple_records_no_allele_match"] += 1
                    writer.writerow(row)
                    continue

                record = candidate["record"]
                row["match_status"] = candidate["match_status"]
                row["position_record_count"] = candidate["position_record_count"]
                row["vcf_chr"] = record.chrom_raw
                row["vcf_pos"] = record.pos
                row["vcf_id"] = record.variant_id
                row["vcf_ref"] = record.ref
                row["vcf_alt"] = candidate["alt"]
                row["vcf_filter"] = record.filt
                stats[candidate["match_status"]] += 1

                csq_entries, raw_entries = parse_csq_entries(
                    record.info, vcf_stream.csq_fields, candidate["alt"], candidate["alt_index"]
                )
                if csq_entries:
                    best_entry = best_csq_entry(csq_entries)
                    row["vep_annotation_count"] = len(csq_entries)
                    row["vep_matched_csq"] = ",".join(raw_entries)
                    for field in vcf_stream.csq_fields:
                        row[f"vep_best_{field}"] = best_entry.get(field, "")
                    stats["with_csq_annotation"] += 1

                writer.writerow(row)

        return stats
    finally:
        vcf_stream.close()


def write_summary(summary_path: Path, stats: dict[str, int]) -> None:
    with summary_path.open("w", encoding="utf-8") as handle:
        for key, value in stats.items():
            handle.write(f"{key}\t{value}\n")


def main() -> int:
    args = parse_args()
    stats = annotate_bim(args)
    if args.summary:
        write_summary(args.summary, stats)
    return 0


if __name__ == "__main__":
    sys.exit(main())
