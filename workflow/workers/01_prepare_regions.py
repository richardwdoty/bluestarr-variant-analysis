#!/usr/bin/env python3

"""
Prepare genomic regions for BlueSTARR prediction generation.

Reads a region file, filters intervals by size and TSS distance, pads each
interval for BlueSTARR windowing, and writes:

1. region_keys_file (with metadata)
2. regions_all_file (prediction_region only, no header)
"""

from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd
import yaml


REQUIRED_COLUMNS = ["chrom", "start", "end", "tss_distance"]


def load_config(path: Path) -> dict:
    with open(path, "r") as f:
        return yaml.safe_load(f)


def resolve_column(df: pd.DataFrame, spec, has_header: bool):
    """
    Resolve a column from either index or name.
    """
    if spec is None:
        return None

    if has_header:
        if spec not in df.columns:
            raise ValueError(f"Column '{spec}' not found in input file.")
        return df[spec]

    return df.iloc[:, spec]


def format_region(chrom, start, end):
    return chrom.astype(str) + ":" + start.astype(str) + "-" + end.astype(str)


def prepare_regions(config: dict):

    run_root = Path(config["paths"]["run_root"])
    raw_dir = Path(config["paths"]["raw_data"])

    workflow_cfg = config["workflow"]["prediction_generation"]
    gen_cfg = config["prediction_generation"]

    input_path = run_root / raw_dir / gen_cfg["input_file"]

    intermediate_dir = run_root / workflow_cfg["intermediate_dir"]
    output_keys = intermediate_dir / gen_cfg["region_keys_file"]
    output_regions = intermediate_dir / gen_cfg["regions_all_file"]

    has_header = gen_cfg["has_header"]
    column_spec = gen_cfg["columns"]

    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    df = pd.read_csv(
        input_path,
        sep="\t",
        header=0 if has_header else None,
    )

    resolved = {}

    for col in REQUIRED_COLUMNS:
        if col not in column_spec:
            raise ValueError(f"Missing required column spec: '{col}'")
        resolved[col] = resolve_column(df, column_spec[col], has_header)

    gene_col = None
    if column_spec.get("gene") is not None:
        gene_col = resolve_column(df, column_spec["gene"], has_header)

    chrom = resolved["chrom"]
    start = pd.to_numeric(resolved["start"])
    end = pd.to_numeric(resolved["end"])
    tss_dist = pd.to_numeric(resolved["tss_distance"]).abs()

    region_size = gen_cfg["region_size"]
    min_size = gen_cfg["min_size"]
    min_distance = gen_cfg["min_distance"]

    region_length = end - start
    midpoint = (start + end) // 2
    midpoint_distance = tss_dist + (region_length // 2)

    keep_mask = (
        (region_length >= min_size)
        & (midpoint_distance >= min_distance)
        & ~chrom.astype(str).str.endswith(("X", "Y"))
    )

    chrom = chrom[keep_mask]
    start = start[keep_mask]
    end = end[keep_mask]
    midpoint = midpoint[keep_mask]
    region_length = region_length[keep_mask]
    midpoint_distance = midpoint_distance[keep_mask]

    if gene_col is not None:
        gene_col = gene_col[keep_mask]

    large_region = region_length >= region_size

    region_start = start - 150
    region_end = end + 149

    region_start.loc[large_region] = (
        midpoint.loc[large_region] - region_size // 2 - 150
    )

    region_end.loc[large_region] = (
        midpoint.loc[large_region] + region_size // 2 + 149
    )

    out_df = pd.DataFrame(
        {
            "source_region": format_region(chrom, start, end),
            "prediction_region": format_region(chrom, region_start, region_end),
            "tss_distance": midpoint_distance,
        }
    )

    if gene_col is not None:
        out_df.insert(2, "gene", gene_col.values)

    output_keys.parent.mkdir(parents=True, exist_ok=True)
    output_regions.parent.mkdir(parents=True, exist_ok=True)

    out_df.to_csv(output_keys, sep="\t", index=False)
    out_df["prediction_region"].to_csv(
        output_regions,
        index=False,
        header=False,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        default="config/config.yaml",
        help="Path to configuration file",
    )

    args = parser.parse_args()

    config = load_config(Path(args.config))
    prepare_regions(config)


if __name__ == "__main__":
    main()