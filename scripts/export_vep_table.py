#!/usr/bin/env python3

import argparse
import csv
import gzip
import sys
from dataclasses import dataclass
from pathlib import Path


IMPACT_RANK = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1}
FEATURE_RANK = {"Transcript": 3, "RegulatoryFeature": 2, "MotifFeature": 1}
PICK_VALUES = {"1", "YES", "Y", "yes", "y"}
COMPRESSED_SUFFIXES = {".gz", ".bgz", ".bgzf"}


@dataclass
class VcfRecord:
    chrom: str
    pos: int
    variant_id: str
    ref: str
    alts: list[str]
    qual: str
    filt: str
    info: str


class VcfReader:
    def __init__(self, vcf_path: Path):
        self.vcf_path = vcf_path
        self.handle = open_text(vcf_path)
        self.csq_fields: list[str] | None = None
        self.current_line: str | None = None
        self._prime()

    def _prime(self) -> None:
        for line in self.handle:
            if line.startswith("##INFO=<ID=CSQ"):
                self.csq_fields = parse_csq_header(line)
                continue
            if line.startswith("#"):
                continue
            self.current_line = line
            return
        self.current_line = None

    def __iter__(self):
        while self.current_line is not None:
            line = self.current_line
            self.current_line = next(self.handle, None)
            yield parse_vcf_record(line)

    def close(self) -> None:
        self.handle.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Flatten a VEP-annotated VCF into a TSV so Bos taurus analyses are not limited to a "
            "small hard-coded subset of annotations."
        )
    )
    parser.add_argument("--vcf", required=True, type=Path, help="Input VEP VCF or VCF.GZ file.")
    parser.add_argument(
        "--output",
        type=Path,
        help="Output TSV/TSV.GZ file. Not required when --show-fields is used.",
    )
    parser.add_argument(
        "--mode",
        choices=("all", "best"),
        default="all",
        help=(
            "Export all matched CSQ entries per allele, or only the best-ranked CSQ entry per "
            "allele. Default: all."
        ),
    )
    parser.add_argument(
        "--fields",
        default="all",
        help=(
            'Comma-separated CSQ fields to export, or "all" to export every CSQ field declared in '
            'the VCF header. Default: all.'
        ),
    )
    parser.add_argument(
        "--include-raw-csq",
        action="store_true",
        help="Include the raw matched CSQ record in an extra output column.",
    )
    parser.add_argument(
        "--drop-unannotated",
        action="store_true",
        help="Skip alleles that have no matched CSQ annotation.",
    )
    parser.add_argument(
        "--show-fields",
        action="store_true",
        help="Print the CSQ fields present in the VCF header and exit.",
    )
    args = parser.parse_args()

    if not args.show_fields and args.output is None:
        parser.error("--output is required unless --show-fields is used.")

    return args


def open_text(path: Path):
    if path.suffix.lower() in COMPRESSED_SUFFIXES:
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return open(path, "r", encoding="utf-8", errors="replace", newline="")


def open_output(path: Path):
    if path.suffix.lower() in COMPRESSED_SUFFIXES:
        return gzip.open(path, "wt", encoding="utf-8", newline="")
    return open(path, "w", encoding="utf-8", newline="")


def parse_csq_header(line: str) -> list[str]:
    marker = "Format: "
    if marker not in line:
        raise RuntimeError("Found a CSQ header, but could not parse its field layout.")
    return line.split(marker, 1)[1].rsplit('"', 1)[0].split("|")


def parse_vcf_record(line: str) -> VcfRecord:
    fields = line.rstrip("\n").split("\t", 8)
    if len(fields) < 8:
        raise RuntimeError(f"VCF line does not have at least 8 columns: {line.rstrip()}")
    chrom, pos, variant_id, ref, alt, qual, filt, info = fields[:8]
    return VcfRecord(
        chrom=chrom,
        pos=int(pos),
        variant_id=variant_id,
        ref=ref,
        alts=alt.split(","),
        qual=qual,
        filt=filt,
        info=info,
    )


def extract_info_value(info: str, key: str) -> str:
    prefix = f"{key}="
    for field in info.split(";"):
        if field.startswith(prefix):
            return field[len(prefix) :]
    return ""


def normalize_requested_fields(requested: str, csq_fields: list[str]) -> list[str]:
    if requested.strip().lower() == "all":
        return list(csq_fields)

    selected = [field.strip() for field in requested.split(",") if field.strip()]
    if not selected:
        raise RuntimeError("No CSQ fields were selected.")

    missing = [field for field in selected if field not in csq_fields]
    if missing:
        available = ", ".join(csq_fields)
        raise RuntimeError(
            "Requested CSQ fields are not present in the VCF header: "
            + ", ".join(missing)
            + f"\nAvailable fields: {available}"
        )

    return selected


def parse_matching_csq_entries(
    info: str, csq_fields: list[str], matched_alt: str, matched_alt_index: int
) -> list[dict[str, object]]:
    raw_csq = extract_info_value(info, "CSQ")
    if not raw_csq:
        return []

    allele_num_key = "ALLELE_NUM" if "ALLELE_NUM" in csq_fields else None
    matched_entries: list[dict[str, object]] = []

    for entry_index, raw_entry in enumerate(raw_csq.split(","), start=1):
        values = raw_entry.split("|")
        if len(values) < len(csq_fields):
            values.extend([""] * (len(csq_fields) - len(values)))
        elif len(values) > len(csq_fields):
            values = values[: len(csq_fields)]

        entry = dict(zip(csq_fields, values))
        allele_matches = entry.get("Allele", "").upper() == matched_alt.upper()
        allele_num_matches = allele_num_key and entry.get(allele_num_key, "") == str(matched_alt_index)

        if allele_matches or allele_num_matches:
            matched_entries.append(
                {
                    "entry_index": entry_index,
                    "raw_entry": raw_entry,
                    "fields": entry,
                }
            )

    return matched_entries


def best_csq_entry(entries: list[dict[str, object]]) -> dict[str, object]:
    def score(payload: dict[str, object]) -> tuple[int, int, int, int, int, int, int, int, int, int]:
        entry = payload["fields"]
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


def emit_row(
    writer: csv.DictWriter,
    record: VcfRecord,
    alt: str,
    alt_index: int,
    matched_entries: list[dict[str, object]],
    selected_fields: list[str],
    include_raw_csq: bool,
    payload: dict[str, object] | None,
) -> None:
    row = {
        "vcf_chrom": record.chrom,
        "vcf_pos": record.pos,
        "vcf_id": record.variant_id,
        "vcf_ref": record.ref,
        "vcf_alt": alt,
        "vcf_alt_index": alt_index,
        "vcf_qual": record.qual,
        "vcf_filter": record.filt,
        "vep_annotation_count": len(matched_entries),
        "vep_entry_index": "" if payload is None else payload["entry_index"],
    }

    if include_raw_csq:
        row["vep_raw_csq"] = "" if payload is None else payload["raw_entry"]

    entry_fields = {} if payload is None else payload["fields"]
    for field in selected_fields:
        row[field] = entry_fields.get(field, "")

    writer.writerow(row)


def export_table(args: argparse.Namespace) -> int:
    reader = VcfReader(args.vcf)
    try:
        if not reader.csq_fields:
            raise RuntimeError("No VEP CSQ header was found in the VCF.")

        if args.show_fields:
            for index, field in enumerate(reader.csq_fields, start=1):
                print(f"{index}\t{field}")
            return 0

        selected_fields = normalize_requested_fields(args.fields, reader.csq_fields)
        output_fields = [
            "vcf_chrom",
            "vcf_pos",
            "vcf_id",
            "vcf_ref",
            "vcf_alt",
            "vcf_alt_index",
            "vcf_qual",
            "vcf_filter",
            "vep_annotation_count",
            "vep_entry_index",
        ]
        if args.include_raw_csq:
            output_fields.append("vep_raw_csq")
        output_fields.extend(selected_fields)

        assert args.output is not None
        args.output.parent.mkdir(parents=True, exist_ok=True)

        with open_output(args.output) as output_handle:
            writer = csv.DictWriter(output_handle, fieldnames=output_fields, delimiter="\t")
            writer.writeheader()

            for record in reader:
                for alt_index, alt in enumerate(record.alts, start=1):
                    matched_entries = parse_matching_csq_entries(
                        record.info, reader.csq_fields, alt, alt_index
                    )

                    if not matched_entries:
                        if not args.drop_unannotated:
                            emit_row(
                                writer,
                                record,
                                alt,
                                alt_index,
                                matched_entries,
                                selected_fields,
                                args.include_raw_csq,
                                None,
                            )
                        continue

                    if args.mode == "best":
                        emit_row(
                            writer,
                            record,
                            alt,
                            alt_index,
                            matched_entries,
                            selected_fields,
                            args.include_raw_csq,
                            best_csq_entry(matched_entries),
                        )
                        continue

                    for payload in matched_entries:
                        emit_row(
                            writer,
                            record,
                            alt,
                            alt_index,
                            matched_entries,
                            selected_fields,
                            args.include_raw_csq,
                            payload,
                        )

        return 0
    finally:
        reader.close()


def main() -> int:
    args = parse_args()
    return export_table(args)


if __name__ == "__main__":
    sys.exit(main())
