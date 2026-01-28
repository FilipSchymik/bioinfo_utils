#!/usr/bin/env python3
import argparse
import os
from Bio.PDB import PDBParser, PPBuilder


def extract_chains_with_biopython(pdb_path):
    """
    Return dict {chain_id: sequence_string} for a PDB, using Biopython's PPBuilder.
    If a chain has multiple polypeptides (breaks), they are concatenated.
    """
    parser = PDBParser(QUIET=True)
    structure_id = os.path.basename(pdb_path)
    structure = parser.get_structure(structure_id, pdb_path)

    ppb = PPBuilder()
    chains = {}

    # usually first model is enough; but we’ll just iterate
    for model in structure:
        for chain in model:
            seq_parts = []
            for pp in ppb.build_peptides(chain):
                seq_parts.append(str(pp.get_sequence()))
            if seq_parts:
                chains[chain.id] = "".join(seq_parts)
    return chains


def write_fasta(records, out_path, line_width=60):
    """
    records: iterable of (header, seq)
    """
    with open(out_path, "w") as handle:
        for header, seq in records:
            handle.write(f">{header}\n")
            for i in range(0, len(seq), line_width):
                handle.write(seq[i:i+line_width] + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Extract protein sequences from PDB files in a directory and write to a FASTA file (Biopython-based)."
    )
    parser.add_argument("pdb_dir", help="Directory with .pdb files")
    parser.add_argument("out_fasta", help="Output FASTA file")
    parser.add_argument(
        "--concat-chains",
        action="store_true",
        help="Concatenate all chains from the same PDB into a single FASTA entry",
    )
    parser.add_argument(
        "--use-filename",
        action="store_true",
        help="Use the PDB filename (without extension) as the base ID (default: use filename anyway, but this flag keeps API stable).",
    )
    args = parser.parse_args()

    pdb_dir = os.path.abspath(args.pdb_dir)
    if not os.path.isdir(pdb_dir):
        raise SystemExit(f"ERROR: {pdb_dir} is not a directory")

    records = []

    for fname in sorted(os.listdir(pdb_dir)):
        if not fname.lower().endswith(".pdb"):
            continue
        pdb_path = os.path.join(pdb_dir, fname)
        base_id = os.path.splitext(fname)[0]  # e.g. 1abc
        chains = extract_chains_with_biopython(pdb_path)

        if args.concat_chains:
            # merge all chains, maintain sorted chain order for determinism
            all_seqs = [chains[cid] for cid in sorted(chains) if chains[cid]]
            if not all_seqs:
                continue
            seq = "".join(all_seqs)
            header = base_id
            records.append((header, seq))
        else:
            # one entry per chain
            for chain_id, seq in sorted(chains.items()):
                if not seq:
                    continue
                header = f"{base_id}"
                records.append((header, seq))

    write_fasta(records, args.out_fasta)


if __name__ == "__main__":
    main()
