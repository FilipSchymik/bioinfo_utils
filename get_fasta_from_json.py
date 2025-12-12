#!/usr/bin/env python3
import json
from pathlib import Path
import argparse


def get_query_from_json(json_path: Path, chain_id: str | None = None) -> tuple[str, str]:
    """
    Extract (header, sequence) for the query protein from an AlphaFold3 JSON file.

    - Uses protein.sequence (ungapped, non-MSA)
    - Picks the first chain by default, or a specific chain if chain_id is given.
    - Header: >{name}|{chain_id}
    """
    with json_path.open() as f:
        data = json.load(f)

    sequences = data.get("sequences", [])
    if not sequences:
        raise ValueError(f"No 'sequences' field in {json_path}")

    # choose chain
    protein = None
    if chain_id is None:
        protein = sequences[0].get("protein", {})
    else:
        for seq in sequences:
            p = seq.get("protein", {})
            if p.get("id") == chain_id:
                protein = p
                break
        if protein is None:
            raise ValueError(f"Chain {chain_id} not found in {json_path}")

    seq = protein.get("sequence")
    if not seq:
        raise ValueError(f"'protein.sequence' not found in {json_path}")

    chain_label = protein.get("id", "A")
    base_name = data.get("name") or json_path.stem

    header = f">{base_name}"
    return header, seq


def main():
    parser = argparse.ArgumentParser(
        description="Extract query sequences from AlphaFold3 JSON into a single FASTA."
    )
    parser.add_argument("input", help="JSON file or directory with JSON files")
    parser.add_argument(
        "-o", "--output",
        default="queries.fasta",
        help="Output FASTA file (default: queries.fasta)"
    )
    parser.add_argument(
        "-c", "--chain",
        help="Protein chain ID to extract (default: first chain in each JSON)"
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    out_path = Path(args.output)

    # collect JSON files
    if input_path.is_dir():
        json_files = sorted(input_path.glob("*.json"))
    else:
        json_files = [input_path]

    if not json_files:
        raise SystemExit("No JSON files found.")

    with out_path.open("w") as out_f:
        for jpath in json_files:
            header, seq = get_query_from_json(jpath, chain_id=args.chain)
            out_f.write(f"{header}\n{seq}\n")

    print(f"Wrote {len(json_files)} sequences to {out_path}")


if __name__ == "__main__":
    main()
