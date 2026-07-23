from __future__ import annotations

import csv
import math
from collections.abc import Hashable, Mapping
from collections.abc import Sequence
from pathlib import Path


RAW_TRANSLOCATION_REPORT_TSV = "raw_translocation_read_counts.tsv"

_SIDE_COLUMNS = (
    ("chrom_a", "point_a_no_spanning_rawread"),
    ("chrom_b", "point_b_no_spanning_rawread"),
)
_TRUE_VALUES = {"1", "true"}
_FALSE_VALUES = {"0", "false"}


def _parse_report_bool(value: str | None, column: str, line_number: int) -> bool:
    normalized = "" if value is None else value.strip().lower()
    if normalized in _TRUE_VALUES:
        return True
    if normalized in _FALSE_VALUES:
        return False
    raise ValueError(
        f"Invalid {column} value at line {line_number}: {value!r}"
    )


def load_normal_prior_excluded_chromosomes(
    report_path: str | Path,
) -> set[str]:
    """Return both chromosomes from every pair with a no-span raw-read flag."""
    report_path = Path(report_path)
    if not report_path.is_file():
        return set()

    with report_path.open("r", newline="") as report_file:
        reader = csv.DictReader(report_file, delimiter="\t")
        required_columns = {
            column
            for side_columns in _SIDE_COLUMNS
            for column in side_columns
        }
        missing_columns = required_columns - set(reader.fieldnames or ())
        if missing_columns:
            missing_text = ", ".join(sorted(missing_columns))
            raise ValueError(
                f"Raw translocation report is missing required columns: {missing_text}"
            )

        excluded_chromosomes = set()
        for line_number, row in enumerate(reader, start=2):
            no_span_flags = [
                _parse_report_bool(
                    row.get(no_span_column), no_span_column, line_number
                )
                for _, no_span_column in _SIDE_COLUMNS
            ]
            if not any(no_span_flags):
                continue

            for chromosome_column, _ in _SIDE_COLUMNS:
                chromosome = (row.get(chromosome_column) or "").strip()
                if chromosome and chromosome != "*":
                    excluded_chromosomes.add(chromosome)

    return excluded_chromosomes


def select_normal_prior_chromosome_columns(
    chromosome_paths: Mapping[str, Hashable],
    path_columns: Mapping[Hashable, int],
    excluded_chromosomes: set[str],
) -> list[tuple[str, int]]:
    """Map eligible chromosome paths to matrix columns, preserving path order."""
    return [
        (chromosome, path_columns[path])
        for chromosome, path in chromosome_paths.items()
        if path in path_columns and chromosome not in excluded_chromosomes
    ]


def calculate_length_normalized_prior_scales(
    base_column_norm_squares: Sequence[float],
    prior_column_norm_squares: Sequence[float],
    depth_row_count: int,
    strength: float,
) -> tuple[list[float], float, float]:
    """Return per-column scales without redistributing excluded prior mass.

    When every base column is eligible, the sum of squared scales is
    ``strength * depth_row_count``, matching the total strength of the former
    scalar normalization.  Each squared scale is proportional to that
    column's clean-depth squared norm, so its strength relative to the data
    curvature is constant across chromosome lengths.

    The base columns, rather than only the currently eligible prior columns,
    define the normalizer.  Excluding one prior therefore removes its
    contribution without strengthening unrelated surviving priors.
    """
    depth_row_count = int(depth_row_count)
    strength = float(strength)
    base_norms = [float(value) for value in base_column_norm_squares]
    prior_norms = [float(value) for value in prior_column_norm_squares]

    if depth_row_count <= 0:
        raise ValueError("depth_row_count must be positive")
    if not math.isfinite(strength) or strength < 0:
        raise ValueError("strength must be finite and non-negative")
    if not base_norms:
        if prior_norms:
            raise ValueError("prior norms require at least one base column norm")
        return [], 0.0, 0.0
    if any(not math.isfinite(value) or value <= 0 for value in base_norms):
        raise ValueError("base column norm squares must be finite and positive")
    if any(not math.isfinite(value) or value <= 0 for value in prior_norms):
        raise ValueError("prior column norm squares must be finite and positive")

    base_norm_square_sum = math.fsum(base_norms)
    length_normalization_factor = (
        strength * depth_row_count / base_norm_square_sum
    )
    scales = [
        math.sqrt(length_normalization_factor * norm_square)
        for norm_square in prior_norms
    ]

    # This is the RMS of the per-column scales over the complete base set.
    # It is also the old scalar scale with a fixed (pre-exclusion) row count.
    reference_scale = math.sqrt(
        strength * depth_row_count / len(base_norms)
    )
    return scales, reference_scale, length_normalization_factor
