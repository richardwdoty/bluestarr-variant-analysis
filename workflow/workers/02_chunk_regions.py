#!/usr/bin/env python3

"""
Extract FASTA sequence for prepared prediction regions and split records into
chunk files for BlueSTARR prediction jobs.
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import yaml


def load_config(path: Path) -> dict:
    with open(path, "r") as f:
        return yaml.safe_load(f)


def run_twobit_to_fa(
    regions_file: Path,
    fasta_file: Path,
    twobit_to_fa: Path,
    twobit_file: Path,
) -> None:
    fasta_file.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        str(twobit_to_fa),
        str(twobit_file),
        "-noMask",
        f"-seqList={regions_file}",
        str(fasta_file),
    ]

    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def read_fasta(fasta_file: Path, allow_ambiguous_bases: bool) -> list[tuple[str, str]]:
    records = []
    skipped = 0

    header = None
    seq_lines = []

    with open(fasta_file, "r") as f:
        for line in f:
            line = line.strip()

            if line.startswith(">"):
                if header is not None:
                    sequence = "".join(seq_lines).upper()
                    if allow_ambiguous_bases or set(sequence).issubset({"A", "C", "G", "T"}):
                        records.append((header, sequence))
                    else:
                        skipped += 1

                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)

    if header is not None:
        sequence = "".join(seq_lines).upper()
        if allow_ambiguous_bases or set(sequence).issubset({"A", "C", "G", "T"}):
            records.append((header, sequence))
        else:
            skipped += 1

    if skipped:
        print(f"Skipped {skipped} FASTA record(s) containing non-ACGT bases.")

    return records


def chromosome_sort_key(chrom: str):
    label = chrom.removeprefix("chr")
    if label.isdigit():
        return (0, int(label))
    return (1, label)


def write_chunks(
    records: list[tuple[str, str]],
    output_dir: Path,
    chunk_size: int,
    chunk_template: str,
    chromosome_index_file: str,
) -> None:
    if not records:
        raise ValueError("No FASTA records available for chunking.")

    output_dir.mkdir(parents=True, exist_ok=True)

    chunk_index = 1
    line_count = 0
    current_chrom = None
    chr_index: dict[str, list[int]] = {}

    chunk_path = output_dir / chunk_template.format(chunk=chunk_index)
    out_fh = open(chunk_path, "w")

    try:
        for header, sequence in records:
            chrom = header.split(":", 1)[0]

            if current_chrom is None:
                current_chrom = chrom
                chr_index[current_chrom] = [chunk_index, chunk_index]

            if chrom != current_chrom or line_count >= chunk_size:
                out_fh.close()

                chunk_index += 1
                line_count = 0
                current_chrom = chrom

                chr_index.setdefault(current_chrom, [chunk_index, chunk_index])

                chunk_path = output_dir / chunk_template.format(chunk=chunk_index)
                out_fh = open(chunk_path, "w")

            out_fh.write(f"{header}\t{sequence}\n")
            line_count += 1
            chr_index[current_chrom][1] = chunk_index

    finally:
        out_fh.close()

    index_path = output_dir / chromosome_index_file
    with open(index_path, "w") as f:
        for chrom in sorted(chr_index, key=chromosome_sort_key):
            first, last = chr_index[chrom]
            f.write(f"{chrom}:{first}-{last}\n")

    print(f"Wrote {chunk_index} chunk file(s) to {output_dir}")
    print(f"Wrote chromosome index to {index_path}")


def resolve_paths(config: dict) -> dict[str, Path]:
    run_root = Path(config["paths"]["run_root"])
    ref_root = Path(config["paths"]["ref_root"])

    workflow_cfg = config["workflow"]["prediction_generation"]
    gen_cfg = config["prediction_generation"]
    tools_cfg = config["tools"]

    intermediate_dir = run_root / workflow_cfg["intermediate_dir"]
    chunks_dir = run_root / workflow_cfg["chunks_dir"]

    return {
        "regions_file": intermediate_dir / gen_cfg["prediction_regions_file"],
        "fasta_file": intermediate_dir / gen_cfg["prediction_fasta_file"],
        "chunks_dir": chunks_dir,
        "twobit_to_fa": Path(config["paths"]["repo_root"]) / tools_cfg["twoBitToFa"],
        "twobit_file": ref_root / tools_cfg["twobit_file"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract FASTA sequences and split prediction regions into chunk files."
    )
    parser.add_argument(
        "--config",
        default="config/config.yaml",
        help="Path to configuration file.",
    )
    parser.add_argument(
        "--overwrite-fasta",
        action="store_true",
        help="Regenerate FASTA even if it already exists.",
    )
    args = parser.parse_args()

    config = load_config(Path(args.config))
    gen_cfg = config["prediction_generation"]
    paths = resolve_paths(config)

    if not paths["regions_file"].exists():
        raise FileNotFoundError(f"Region list not found: {paths['regions_file']}")

    if args.overwrite_fasta or not paths["fasta_file"].exists():
        run_twobit_to_fa(
            regions_file=paths["regions_file"],
            fasta_file=paths["fasta_file"],
            twobit_to_fa=paths["twobit_to_fa"],
            twobit_file=paths["twobit_file"],
        )
    else:
        print(f"Using existing FASTA file: {paths['fasta_file']}")

    records = read_fasta(
        paths["fasta_file"],
        allow_ambiguous_bases=gen_cfg.get("allow_ambiguous_bases", False),
    )

    write_chunks(
        records=records,
        output_dir=paths["chunks_dir"],
        chunk_size=gen_cfg["chunk_size"],
        chunk_template=gen_cfg["chunk_filename_template"],
        chromosome_index_file=gen_cfg["chromosome_index_file"],
    )


if __name__ == "__main__":
    main()
