#!/usr/bin/env python3
"""
Parse the full chimeric rabbit/human IgG1 constructs (data/raw/sequences/antibody_sequences_full_constructs.fasta),
strip the shared 19-residue signal peptide, and write:
  - data/processed/antibody_sequences_mature.fasta   (all mature HC/LC records)
  - data/processed/by_clone/<clone>.fasta            (paired mature HC+LC per clone, for structure prediction input)
  - data/raw/sequences/antibody_metadata.csv          (clone, chain, lengths)
"""
import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RAW_FASTA = ROOT / "data/raw/sequences/antibody_sequences_full_constructs.fasta"
PROCESSED_DIR = ROOT / "data/processed"
BY_CLONE_DIR = PROCESSED_DIR / "by_clone"
METADATA_CSV = ROOT / "data/raw/sequences/antibody_metadata.csv"

SIGNAL_PEPTIDE = "MDWTWRILFLVAAATGAHS"  # 19 aa, shared across all constructs


def parse_fasta(path):
    records = {}
    header = None
    seq_lines = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header:
                records[header] = "".join(seq_lines)
            header = line[1:]
            seq_lines = []
        else:
            seq_lines.append(line)
    if header:
        records[header] = "".join(seq_lines)
    return records


def main():
    records = parse_fasta(RAW_FASTA)

    clones = {}
    for header, seq in records.items():
        clone, chain = re.match(r"(.+)_(HC|LC)$", header).groups()
        assert seq.startswith(SIGNAL_PEPTIDE), f"{header} does not start with expected signal peptide"
        mature = seq[len(SIGNAL_PEPTIDE):]
        clones.setdefault(clone, {})[chain] = {"full": seq, "mature": mature}

    BY_CLONE_DIR.mkdir(parents=True, exist_ok=True)

    with open(PROCESSED_DIR / "antibody_sequences_mature.fasta", "w") as mature_fasta, \
         open(METADATA_CSV, "w", newline="") as meta_file:

        meta_writer = csv.writer(meta_file)
        meta_writer.writerow([
            "clone", "chain", "target_antigen", "isotype",
            "full_construct_length", "mature_length", "has_LALA_PG",
        ])

        for clone in sorted(clones):
            chains = clones[clone]
            for chain in ("HC", "LC"):
                if chain not in chains:
                    continue
                full = chains[chain]["full"]
                mature = chains[chain]["mature"]
                mature_fasta.write(f">{clone}_{chain}\n{mature}\n")
                meta_writer.writerow([
                    clone, chain, "TBD",
                    "chimeric rabbit/human IgG1" if chain == "HC" else "chimeric rabbit/human (kappa/lambda)",
                    len(full), len(mature),
                    "yes" if chain == "HC" else "n/a",
                ])

            with open(BY_CLONE_DIR / f"{clone}.fasta", "w") as clone_fasta:
                for chain in ("HC", "LC"):
                    if chain in chains:
                        clone_fasta.write(f">{clone}_{chain}\n{chains[chain]['mature']}\n")

    print(f"Parsed {len(clones)} clones ({len(records)} chain records).")
    print(f"Wrote {PROCESSED_DIR / 'antibody_sequences_mature.fasta'}")
    print(f"Wrote {len(clones)} paired FASTA files to {BY_CLONE_DIR}")
    print(f"Wrote {METADATA_CSV}")


if __name__ == "__main__":
    main()
