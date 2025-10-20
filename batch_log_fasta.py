#!/usr/bin/env python3
"""
Simple FASTA parser/batcher for protein sequences.

Features
- "Sterilizes" FASTA headers: keeps only the ID (text after ">" up to first whitespace), **and uppercases it**.
- Splits sequences into batches by sequence count
- Writes a summary CSV with columns: id,length,has_nonstandard.

CLI
----
python fasta_batcher.py \
  --input path/to/input.fasta[.gz] \
  --max-seqs 1000 \
  --summary path/to/summary.csv

Outputs
- Summary CSV.
- Batch FASTA files in a sibling directory derived from the CSV name: <summary_stem>_batches/
  Files are named batch_0001.fasta, batch_0002.fasta, ...

Notes
- Streams the input (no whole-file loading), suitable for large FASTAs.
- Accepts gzipped input (.gz).
- Preserves input order; batches contain up to N sequences each.
- Nonstandard amino acids are any characters outside the 20 standard set: ACDEFGHIKLMNPQRSTVWY (case-insensitive).
"""
from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path
from typing import Generator, Iterable, Tuple, TextIO

# 20 standard amino acids
STD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def open_maybe_gzip(path: Path) -> TextIO:
    """Open a file that may be gzipped, in text mode (utf-8)."""
    if str(path).endswith(".gz"):
        # text mode directly, avoids wrapping a binary stream
        return gzip.open(path, "rt", encoding="utf-8", newline=None)
    return open(path, "r", encoding="utf-8", newline=None)


def sterilize_header(raw_header: str) -> str:
    """
    Return only the ID portion of a FASTA header (after '>' up to first whitespace),
    uppercased.
    """
    h = raw_header.strip()
    if h.startswith(">"):
        h = h[1:]
    sid = h.split()[0] if h else ""
    return sid.upper()


def fasta_iter(handle: TextIO) -> Generator[Tuple[str, str], None, None]:
    """Get (sterilized_id, sequence) from a FASTA stream.

    - Collapses line breaks; strips *all* whitespace (spaces, tabs, unusual unicode spaces).
    - Upper-cases sequences.
    - Skips empty sequences (headers with no residues).
    """
    header: str | None = None
    seq_chunks: list[str] = []
    for line in handle:
        if not line:
            continue
        if line.startswith(">"):
            # flush previous record
            if header is not None:
                # remove all whitespace from sequence
                seq = "".join("".join(seq_chunks).split()).upper()
                if seq:
                    yield sterilize_header(header), seq
                seq_chunks.clear()
            header = line.rstrip("\r\n")
        else:
            seq_chunks.append(line.rstrip("\r\n"))

    # flush last record
    if header is not None:
        seq = "".join("".join(seq_chunks).split()).upper()
        if seq:
            yield sterilize_header(header), seq


def has_nonstandard_aa(seq: str) -> bool:
    """True if sequence contains any residue outside the 20 standard amino acids."""
    return any(ch not in STD_AA for ch in seq)


def wrap_fasta(seq: str, width: int = 60) -> Iterable[str]:
    for i in range(0, len(seq), width):
        yield seq[i : i + width]


def write_batch(batch_idx: int, records: list[Tuple[str, str]], out_dir: Path) -> Path:
    out_path = out_dir / f"batch_{batch_idx:04d}.fasta"
    with open(out_path, "w", encoding="utf-8", newline="\n") as w:
        for sid, seq in records:
            w.write(f">{sid}\n")
            for line in wrap_fasta(seq, 60):
                w.write(line + "\n")
    return out_path


def process(input_path: Path, max_seqs: int, summary_csv: Path) -> None:
    out_dir = summary_csv.with_suffix("").with_name(summary_csv.stem + "_batches")
    out_dir.mkdir(parents=True, exist_ok=True)

    # Prepare CSV
    summary_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(summary_csv, "w", newline="", encoding="utf-8") as csv_f:
        csv_writer = csv.writer(csv_f)
        csv_writer.writerow(["id", "length", "has_nonstandard"])  # header

        batch_records: list[Tuple[str, str]] = []
        batch_idx = 1
        n_in_batch = 0

        def flush_batch() -> None:
            nonlocal batch_records, batch_idx, n_in_batch
            if batch_records:
                write_batch(batch_idx, batch_records, out_dir)
                batch_idx += 1
                batch_records = []
                n_in_batch = 0

        with open_maybe_gzip(input_path) as handle:
            for sid, seq in fasta_iter(handle):
                # summary row
                csv_writer.writerow([sid, len(seq), has_nonstandard_aa(seq)])

                # pack by sequence count only
                if n_in_batch >= max_seqs and batch_records:
                    flush_batch()

                batch_records.append((sid, seq))
                n_in_batch += 1

        # flush last batch
        flush_batch()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Sterilize FASTA headers (uppercase IDs), batch by sequence count, and summarize."
    )
    p.add_argument("--input", required=True, type=Path, help="Input FASTA (.fa/.fasta[.gz])")
    p.add_argument("--max-seqs", required=True, type=int, help="Max number of sequences per batch (e.g., 1000)")
    p.add_argument("--summary", required=True, type=Path, help="Output summary CSV path")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    if args.max_seqs <= 0:
        raise SystemExit("--max-seqs must be a positive integer")
    process(args.input, args.max_seqs, args.summary)


if __name__ == "__main__":
    main()
