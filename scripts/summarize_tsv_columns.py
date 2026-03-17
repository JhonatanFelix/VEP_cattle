#!/usr/bin/env python3

import argparse
import csv
import gzip
import json
import math
import re
from collections import Counter
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create dataset- and column-level summary statistics for a TSV/TSV.GZ file, "
            "including missingness and full value counts for every column."
        )
    )
    parser.add_argument(
        "input_file",
        type=Path,
        help="Input TSV file. Gzipped .tsv.gz files are supported.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Directory where the reports will be written. Defaults to <input>__summary.",
    )
    parser.add_argument(
        "--missing-token",
        action="append",
        default=[],
        help=(
            "Extra token to treat as missing. Repeat this option to add multiple tokens. "
            "Empty strings are always treated as missing."
        ),
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=10,
        help="Number of most frequent values to expose directly in the column summary table.",
    )
    return parser.parse_args()


def default_output_dir(input_file: Path) -> Path:
    name = input_file.name
    for suffix in (".tsv.gz", ".txt.gz", ".csv.gz", ".gz", ".tsv", ".txt", ".csv"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return input_file.parent / f"{name}__summary"


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="")
    return open(path, "r", newline="")


def open_text_write(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "wt", newline="")
    return open(path, "w", newline="")


def sanitize_column_name(column_name: str, index: int) -> str:
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", column_name).strip("._")
    if not safe:
        safe = f"column_{index + 1}"
    return f"{index + 1:02d}_{safe}"


def is_int(value: str) -> bool:
    try:
        int(value)
        return True
    except ValueError:
        return False


def is_float(value: str) -> bool:
    try:
        float(value)
        return True
    except ValueError:
        return False


def truncate_preview(value: str, limit: int = 120) -> str:
    if len(value) <= limit:
        return value
    return value[: limit - 3] + "..."


@dataclass
class ColumnProfile:
    index: int
    name: str
    missing_count: int = 0
    non_missing_count: int = 0
    value_counts: Counter = field(default_factory=Counter)
    example_value: str = ""
    min_length: int | None = None
    max_length: int | None = None
    int_candidate: bool = True
    float_candidate: bool = True
    numeric_count: int = 0
    text_count: int = 0
    numeric_min: float | None = None
    numeric_max: float | None = None
    numeric_mean: float = 0.0
    numeric_m2: float = 0.0

    def update(self, raw_value: str, missing_tokens: set[str]) -> None:
        value = raw_value
        if value == "" or value in missing_tokens:
            self.missing_count += 1
            self.value_counts[value] += 1
            return

        self.non_missing_count += 1
        self.value_counts[value] += 1

        if not self.example_value:
            self.example_value = value

        length = len(value)
        if self.min_length is None or length < self.min_length:
            self.min_length = length
        if self.max_length is None or length > self.max_length:
            self.max_length = length

        value_is_int = is_int(value)
        value_is_float = value_is_int or is_float(value)

        if not value_is_int:
            self.int_candidate = False
        if not value_is_float:
            self.float_candidate = False

        if value_is_float:
            number = float(value)
            self.numeric_count += 1
            if self.numeric_min is None or number < self.numeric_min:
                self.numeric_min = number
            if self.numeric_max is None or number > self.numeric_max:
                self.numeric_max = number

            delta = number - self.numeric_mean
            self.numeric_mean += delta / self.numeric_count
            delta2 = number - self.numeric_mean
            self.numeric_m2 += delta * delta2
        else:
            self.text_count += 1

    def detected_type(self) -> str:
        if self.non_missing_count == 0:
            return "empty"
        if self.numeric_count == self.non_missing_count:
            if self.int_candidate:
                return "integer"
            if self.float_candidate:
                return "float"
        if self.numeric_count > 0 and self.text_count > 0:
            return "mixed"
        return "text"

    def numeric_std(self) -> float | None:
        if self.numeric_count < 2:
            return 0.0 if self.numeric_count == 1 else None
        return math.sqrt(self.numeric_m2 / (self.numeric_count - 1))

    def top_values(self, top_n: int) -> list[tuple[str, int]]:
        return self.value_counts.most_common(top_n)


def profile_tsv(
    input_file: Path, missing_tokens: set[str]
) -> tuple[list[ColumnProfile], dict[str, object], list[str]]:
    with open_text(input_file) as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise RuntimeError(f"{input_file} is empty.") from exc

        profiles = [ColumnProfile(index=i, name=name) for i, name in enumerate(header)]
        row_count = 0
        total_missing_cells = 0
        rows_with_any_missing = 0
        min_fields = len(header)
        max_fields = len(header)

        for row_number, row in enumerate(reader, start=2):
            row_count += 1
            min_fields = min(min_fields, len(row))
            max_fields = max(max_fields, len(row))

            if len(row) < len(header):
                row = row + [""] * (len(header) - len(row))
            elif len(row) > len(header):
                raise RuntimeError(
                    f"Row {row_number} has {len(row)} fields, but the header has {len(header)}."
                )

            row_missing = False
            for profile, value in zip(profiles, row):
                was_missing = value == "" or value in missing_tokens
                profile.update(value, missing_tokens)
                if was_missing:
                    total_missing_cells += 1
                    row_missing = True
            if row_missing:
                rows_with_any_missing += 1

    dataset_info = {
        "input_file": str(input_file.resolve()),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "rows": row_count,
        "columns": len(header),
        "shape": [row_count, len(header)],
        "column_names": header,
        "delimiter": "\\t",
        "compression": "gzip" if input_file.suffix == ".gz" else "none",
        "missing_tokens": [""] + sorted(missing_tokens),
        "total_missing_cells": total_missing_cells,
        "rows_with_any_missing": rows_with_any_missing,
        "rows_without_missing": row_count - rows_with_any_missing,
        "min_fields_per_row": min_fields,
        "max_fields_per_row": max_fields,
    }
    return profiles, dataset_info, header


def write_dataset_info(output_dir: Path, dataset_info: dict[str, object]) -> None:
    json_path = output_dir / "dataset_info.json"
    text_path = output_dir / "dataset_info.txt"

    json_path.write_text(json.dumps(dataset_info, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")

    lines = [
        f"input_file\t{dataset_info['input_file']}",
        f"generated_at_utc\t{dataset_info['generated_at_utc']}",
        f"rows\t{dataset_info['rows']}",
        f"columns\t{dataset_info['columns']}",
        f"shape\t{dataset_info['shape'][0]} x {dataset_info['shape'][1]}",
        f"compression\t{dataset_info['compression']}",
        f"delimiter\t{dataset_info['delimiter']}",
        f"total_missing_cells\t{dataset_info['total_missing_cells']}",
        f"rows_with_any_missing\t{dataset_info['rows_with_any_missing']}",
        f"rows_without_missing\t{dataset_info['rows_without_missing']}",
        f"min_fields_per_row\t{dataset_info['min_fields_per_row']}",
        f"max_fields_per_row\t{dataset_info['max_fields_per_row']}",
        "column_names\t" + ", ".join(dataset_info["column_names"]),
        "missing_tokens\t" + ", ".join(repr(token) for token in dataset_info["missing_tokens"]),
    ]
    text_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_column_summary(
    output_dir: Path,
    profiles: list[ColumnProfile],
    row_count: int,
    top_n: int,
    missing_tokens: set[str],
) -> None:
    summary_path = output_dir / "column_summary.tsv"
    fieldnames = [
        "column_index",
        "column_name",
        "detected_type",
        "non_missing_count",
        "missing_count",
        "missing_pct",
        "unique_value_count_including_missing",
        "unique_non_missing_count",
        "example_value",
        "min_length",
        "max_length",
        "numeric_min",
        "numeric_max",
        "numeric_mean",
        "numeric_std",
    ]
    for rank in range(1, top_n + 1):
        fieldnames.extend(
            [
                f"top_{rank}_value",
                f"top_{rank}_count",
                f"top_{rank}_pct",
            ]
        )

    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for profile in profiles:
            column_type = profile.detected_type()
            all_missing_tokens = {""} | missing_tokens
            non_missing_unique = len(
                [value for value in profile.value_counts if value not in all_missing_tokens]
            )
            row = {
                "column_index": profile.index + 1,
                "column_name": profile.name,
                "detected_type": column_type,
                "non_missing_count": profile.non_missing_count,
                "missing_count": profile.missing_count,
                "missing_pct": round((profile.missing_count / row_count) * 100, 6) if row_count else 0.0,
                "unique_value_count_including_missing": len(profile.value_counts),
                "unique_non_missing_count": non_missing_unique,
                "example_value": truncate_preview(profile.example_value),
                "min_length": "" if profile.min_length is None else profile.min_length,
                "max_length": "" if profile.max_length is None else profile.max_length,
                "numeric_min": ""
                if column_type not in {"integer", "float"} or profile.numeric_min is None
                else profile.numeric_min,
                "numeric_max": ""
                if column_type not in {"integer", "float"} or profile.numeric_max is None
                else profile.numeric_max,
                "numeric_mean": ""
                if column_type not in {"integer", "float"} or profile.numeric_count == 0
                else round(profile.numeric_mean, 6),
                "numeric_std": ""
                if column_type not in {"integer", "float"} or profile.numeric_std() is None
                else round(profile.numeric_std(), 6),
            }

            for rank, (value, count) in enumerate(profile.top_values(top_n), start=1):
                row[f"top_{rank}_value"] = truncate_preview(value)
                row[f"top_{rank}_count"] = count
                row[f"top_{rank}_pct"] = round((count / row_count) * 100, 6) if row_count else 0.0

            writer.writerow(row)


def write_value_counts(
    output_dir: Path, profiles: list[ColumnProfile], row_count: int, missing_tokens: set[str]
) -> None:
    value_counts_dir = output_dir / "value_counts"
    value_counts_dir.mkdir(exist_ok=True)
    all_missing_tokens = {""} | missing_tokens

    for profile in profiles:
        filename = sanitize_column_name(profile.name, profile.index) + ".tsv.gz"
        path = value_counts_dir / filename
        with open_text_write(path) as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["value", "count", "percentage", "is_missing"])
            for value, count in sorted(
                profile.value_counts.items(),
                key=lambda item: (-item[1], item[0]),
            ):
                writer.writerow(
                    [
                        value,
                        count,
                        round((count / row_count) * 100, 6) if row_count else 0.0,
                        "yes" if value in all_missing_tokens else "no",
                    ]
                )


def main() -> int:
    args = parse_args()
    input_file = args.input_file
    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")

    output_dir = args.output_dir or default_output_dir(input_file)
    output_dir.mkdir(parents=True, exist_ok=True)

    missing_tokens = {token for token in args.missing_token if token != ""}
    profiles, dataset_info, _header = profile_tsv(input_file, missing_tokens)

    write_dataset_info(output_dir, dataset_info)
    write_column_summary(output_dir, profiles, dataset_info["rows"], args.top_n, missing_tokens)
    write_value_counts(output_dir, profiles, dataset_info["rows"], missing_tokens)

    print(f"Summary written to: {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
