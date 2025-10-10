#!/usr/bin/env python3
"""
Efficient pLDDT extractor for AlphaFold/ESMFold-style PDB files (plain .pdb only).

Scans a directory for .pdb files, parses ATOM records, and computes:
- ID: filename without the .pdb suffix
- protein length: number of unique residues with ATOM records
- average pLDDT: mean of per-residue pLDDT (read from the B-factor column)

Designed for very large corpora (millions of files):
- Streaming, fixed-column byte slicing (no regex/split in the hot path)
- Per-residue rolling aggregation (no dict of all residues)
- Multiprocessing with chunked dispatch and incremental CSV writing

Usage:
    python plddt_pdb_dir_to_csv.py \
        --input /path/to/pdb_dir \
        --output results.csv \
        --workers 16 \
        --recursive

Notes:
- pLDDT is read from the PDB B-factor column (cols 61-66, 0-based slice [60:66]).
- Residues are identified by (chainID, resSeq, iCode) and averaged per residue.
- Alternate locations: no filtering (safe because we average per residue).
- Files without ATOM records will be reported with length=0 and avg_pLDDT empty.
"""

from __future__ import annotations
import os
import argparse
import math
from concurrent.futures import ProcessPoolExecutor
from typing import Iterator, Tuple, Optional
import csv

ATOM_PREFIX = b"ATOM"  # match "ATOM" at start (robust to spacing)

# Fixed-column indices (0-based, end-exclusive) according to PDB format
SLICE_CHAIN = slice(21, 22)       # Chain identifier (col 22)
SLICE_RESSEQ = slice(22, 26)      # Residue sequence number (cols 23-26)
SLICE_ICODE = slice(26, 27)       # Insertion code (col 27)
SLICE_ALTLOC = slice(16, 17)      # Alternate location indicator (col 17)  [kept for completeness]
SLICE_OCC = slice(54, 60)         # Occupancy (cols 55-60)
SLICE_BFACTOR = slice(60, 66)     # B-factor (pLDDT) (cols 61-66)


def iter_pdb_paths(root: str, recursive: bool = False) -> Iterator[str]:
    """Yield absolute paths to `.pdb` files under root."""
    root = os.path.abspath(root)
    def is_pdb_name(name: str) -> bool:
        return name.endswith('.pdb')
    if recursive:
        for dirpath, _dirnames, filenames in os.walk(root):
            for fn in filenames:
                if is_pdb_name(fn):
                    yield os.path.join(dirpath, fn)
    else:
        with os.scandir(root) as it:
            for entry in it:
                if entry.is_file() and is_pdb_name(entry.name):
                    yield entry.path


def parse_pdb_plddt(path: str) -> Tuple[str, int, Optional[float]]:
    """
    Parse a PDB to compute protein length and average per-residue pLDDT.
    Returns: (id_without_suffix, length, avg_plddt or None if no residues had parsable pLDDT)
    """
    try:
        with open(path, 'rb', buffering=1024 * 1024) as fh:
            current_res = None  # type: Optional[Tuple[bytes, bytes, bytes]]
            atom_sum = 0.0
            atom_count = 0

            residue_sum = 0.0
            residue_count_with_plddt = 0
            length = 0

            for line in fh:
                # Work in bytes for speed
                if not line.startswith(ATOM_PREFIX):
                    continue

                # Identify residue
                resid = (line[SLICE_CHAIN], line[SLICE_RESSEQ], line[SLICE_ICODE])

                # New residue? finalize previous
                if current_res is None:
                    current_res = resid
                elif resid != current_res:
                    if atom_count:
                        residue_sum += (atom_sum / atom_count)
                        residue_count_with_plddt += 1
                    length += 1
                    current_res = resid
                    atom_sum = 0.0
                    atom_count = 0

                # --- B-factor (pLDDT) extraction ---
                bval = None
                # Primary: strict fixed slice
                try:
                    bbytes = line[SLICE_BFACTOR].strip()
                    if bbytes:
                        bval = float(bbytes)
                except Exception:
                    bval = None

                # Fallback 1: slightly wider slice to tolerate drift
                if bval is None and len(line) >= 70:
                    try:
                        bval = float(line[60:70].strip() or b'')
                    except Exception:
                        bval = None

                # Fallback 2: look at occupancy and b-factor fields; pick value that isn't occupancy-like
                if bval is None:
                    try:
                        occ_bytes = line[SLICE_OCC].strip()
                        bf_bytes = line[SLICE_BFACTOR].strip()
                        occ = float(occ_bytes) if occ_bytes else None
                        bf  = float(bf_bytes) if bf_bytes else None
                        if bf is not None and not math.isnan(bf):
                            bval = bf
                        elif occ is not None and occ > 1.5:  # not a typical occupancy
                            bval = occ
                    except Exception:
                        pass

                # Fallback 3: token parse; take the last numeric token > 1.5 (avoid occupancy ~1.00)
                if bval is None:
                    try:
                        s = line.decode('ascii', 'ignore')
                        toks = s.split()
                        for tok in reversed(toks):
                            try:
                                v = float(tok)
                                if v > 1.5:
                                    bval = v
                                    break
                            except ValueError:
                                continue
                    except Exception:
                        pass

                if bval is not None and not math.isnan(bval):
                    atom_sum += bval
                    atom_count += 1

            # Finalize last residue
            if current_res is not None:
                if atom_count:
                    residue_sum += (atom_sum / atom_count)
                    residue_count_with_plddt += 1
                length += 1

        # Build ID and result
        file_id = os.path.basename(path)
        if file_id.endswith('.pdb'):
            file_id = file_id[:-4]

        if length == 0:
            return (file_id, 0, None)
        avg_plddt = (residue_sum / residue_count_with_plddt) if residue_count_with_plddt else None
        return (file_id, length, avg_plddt)

    except Exception:
        file_id = os.path.basename(path)
        if file_id.endswith('.pdb'):
            file_id = file_id[:-4]
        return (file_id, 0, None)


def main():
    ap = argparse.ArgumentParser(description='Parse pLDDT from PDB files into a CSV.')
    ap.add_argument('--input', required=True, help='Directory containing .pdb files')
    ap.add_argument('--output', required=True, help='Output CSV path')
    ap.add_argument('--workers', type=int, default=os.cpu_count() or 1,
                    help='Number of worker processes (default: number of CPUs)')
    ap.add_argument('--recursive', action='store_true', help='Recurse into subdirectories')
    ap.add_argument('--chunksize', type=int, default=64,
                    help='Task chunk size for multiprocessing (default: 64)')
    ap.add_argument('--failures', default=None,
                    help='Optional path to write a list of IDs that could not be parsed')
    args = ap.parse_args()

    paths = list(iter_pdb_paths(args.input, args.recursive))

    # Prepare CSV writer; write incrementally to avoid large memory use
    os.makedirs(os.path.dirname(os.path.abspath(args.output)) or '.', exist_ok=True)
    failures = []

    with open(args.output, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['ID', 'length', 'avg_pLDDT'])

        if args.workers == 1:
            for p in paths:
                file_id, length, avg = parse_pdb_plddt(p)
                writer.writerow([file_id, length, f"{avg:.6f}" if avg is not None else ""])
                if avg is None:
                    failures.append(file_id)
        else:
            with ProcessPoolExecutor(max_workers=args.workers) as ex:
                for file_id, length, avg in ex.map(parse_pdb_plddt, paths, chunksize=args.chunksize):
                    writer.writerow([file_id, length, f"{avg:.6f}" if avg is not None else ""]) 
                    if avg is None:
                        failures.append(file_id)

    if args.failures and failures:
        with open(args.failures, 'w') as fh:
            for fid in failures:
                fh.write(f"{fid}")


if __name__ == '__main__':
    main()
