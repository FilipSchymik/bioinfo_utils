#!/usr/bin/env python3
import json
import argparse
from pathlib import Path

def extract_msa_from_json(json_path: Path, chain_id: str | None = None) -> str:
    """Return the unpairedMsa string for a given chain (or first chain if not specified)."""
    with json_path.open() as f:
        data = json.load(f)

    sequences = data.get("sequences", [])
    if not sequences:
        raise ValueError(f"No 'sequences' field in {json_path}")

    if chain_id is None:
        # Take the first sequence's MSA
        protein = sequences[0].get("protein", {})
    else:
        # Pick specific chain, e.g. "A"
        protein = None
        for seq in sequences:
            p = seq.get("protein", {})
            if p.get("id") == chain_id:
                protein = p
                break
        if protein is None:
            raise ValueError(f"Chain {chain_id} not found in {json_path}")

    msa = protein.get("unpairedMsa")
    if msa is None:
        raise ValueError(f"'unpairedMsa' not found in {json_path} for chain {chain_id or sequences[0]['protein'].get('id')}")
    return msa

def main():
    parser = argparse.ArgumentParser(description="Extract AlphaFold3 MSA (unpairedMsa) to .a3m")
    parser.add_argument("input", help="Input JSON file or directory with JSONs")
    parser.add_argument("-o", "--outdir", help="Output directory for .a3m files (default: same as input)")
    parser.add_argument("-c", "--chain", help="Protein chain ID to extract (default: first chain)")
    args = parser.parse_args()

    input_path = Path(args.input)
    outdir = Path(args.outdir) if args.outdir else None

    # Collect all JSON files
    if input_path.is_dir():
        json_files = sorted(input_path.glob("*.json"))
    else:
        json_files = [input_path]

    if not json_files:
        raise SystemExit("No JSON files found.")

    for jpath in json_files:
        msa = extract_msa_from_json(jpath, chain_id=args.chain)

        if outdir:
            outdir.mkdir(parents=True, exist_ok=True)
            out_path = outdir / (jpath.stem + ".a3m")
        else:
            out_path = jpath.with_suffix(".a3m")

        # Write exactly what’s in unpairedMsa
        # It already contains '>query', sequences, and newline characters
        with out_path.open("w") as f:
            # Ensure a trailing newline just in case
            if not msa.endswith("\n"):
                msa += "\n"
            f.write(msa)

        print(f"Wrote {out_path}")

if __name__ == "__main__":
    main()
