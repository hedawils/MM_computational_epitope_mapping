#!/usr/bin/env python3
"""
CDR-only sequence similarity analysis, BLOSUM62-scored, for the 15 clones.

Steps:
  1. Load VH/VL variable domains (data/processed/variable_domains.fasta, produced
     by sequence_similarity_analysis.py).
  2. Extract CDR1/2/3 per chain using conserved-framework anchor motifs. This is a
     heuristic tuned to the 15 sequences actually observed in this dataset (all
     motifs below were verified present in every clone) -- it is NOT a validated
     Kabat/Chothia/IMGT numbering tool and should not be assumed correct for other
     antibody sequences. See CDR_ANCHOR_NOTES below.
  3. Concatenate CDR1+CDR2+CDR3 per chain and run pairwise global alignment scored
     with BLOSUM62 (affine gap penalties) via Bio.Align.PairwiseAligner.
  4. Percent identity matrices + hierarchical clustering dendrograms, analogous to
     sequence_similarity_analysis.py but restricted to CDR regions.

CDR_ANCHOR_NOTES
-----------------
VH (2 conserved Cys mark FR1/CDR-H1 and FR3/CDR-H3 boundaries exactly):
  FR1 | CDR-H1 | FR2 | CDR-H2 | FR3 | CDR-H3 | FR4
   ...C   ...   W[VF]RQAPG  ...  LE[WY]IG/LQWIG  ...  R[FSA]T[ILV][ST]  ...C  ...  WG.GT ...

VL (kappa + lambda; 2 conserved Cys mark FR1/CDR-L1 and FR3/CDR-L3 boundaries exactly):
  FR1 | CDR-L1 | FR2 | CDR-L2 | FR3 | CDR-L3 | FR4
   ...C   ...   W[FYI]QQ ... (LLIY|LLIH|LVIY|LLIF|LMQL)  ...  G.P[SD]R  ...C  ...  FG.GT ...
"""
import csv
import re
from pathlib import Path
from itertools import combinations

import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio.Align import PairwiseAligner, substitution_matrices

ROOT = Path(__file__).resolve().parent.parent
VAR_DOMAINS_FASTA = ROOT / "data/processed/variable_domains.fasta"
METADATA_CSV = ROOT / "data/raw/sequences/antibody_metadata.csv"
PROCESSED_DIR = ROOT / "data/processed"
RESULTS_DIR = ROOT / "results/cdr_similarity"

VH_FR2_START = re.compile(r"W[VF]RQAPG")
VH_FR2_END = re.compile(r"LE[WY]IG|LQWIG")
VH_FR3_START = re.compile(r"R[FSA]T[ILV][ST]")
VH_FR4_START = re.compile(r"WG[A-Z]GT")

VL_FR2_START = re.compile(r"W[FYI]QQ")
VL_FR2_END = re.compile(r"LLIY|LLIH|LVIY|LLIF|LMQL")
VL_FR3_START = re.compile(r"G[A-Z]P[SD]R")
VL_FR4_START = re.compile(r"FG[A-Z]GT")


def parse_fasta(path):
    records = {}
    header, seq_lines = None, []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header:
                records[header] = "".join(seq_lines)
            header, seq_lines = line[1:], []
        else:
            seq_lines.append(line)
    if header:
        records[header] = "".join(seq_lines)
    return records


def load_metadata():
    meta = {}
    with open(METADATA_CSV) as f:
        for row in csv.DictReader(f):
            meta[(row["clone"], row["chain"])] = row
    return meta


def extract_vh_cdrs(seq, clone):
    c1 = seq.index("C")
    m2 = VH_FR2_START.search(seq, c1 + 1)
    m3 = VH_FR2_END.search(seq, m2.end())
    m4 = VH_FR3_START.search(seq, m3.end())
    c2 = seq.index("C", m4.end())
    m6 = VH_FR4_START.search(seq, c2 + 1)

    cdr1 = seq[c1 + 1: m2.start()]
    cdr2 = seq[m3.end(): m4.start()]
    cdr3 = seq[c2 + 1: m6.start()]
    for name, cdr in (("CDR-H1", cdr1), ("CDR-H2", cdr2), ("CDR-H3", cdr3)):
        assert 2 <= len(cdr) <= 30, f"{clone} {name} implausible length {len(cdr)}: {cdr!r}"
    return cdr1, cdr2, cdr3


def extract_vl_cdrs(seq, clone):
    c1 = seq.index("C")
    m2 = VL_FR2_START.search(seq, c1 + 1)
    m2tail = VL_FR2_END.search(seq, m2.end())
    m5 = VL_FR3_START.search(seq, m2tail.end())
    c2 = seq.index("C", m5.end())
    m7 = VL_FR4_START.search(seq, c2 + 1)

    cdr1 = seq[c1 + 1: m2.start()]
    cdr2 = seq[m2tail.end(): m5.start()]
    cdr3 = seq[c2 + 1: m7.start()]
    for name, cdr in (("CDR-L1", cdr1), ("CDR-L2", cdr2), ("CDR-L3", cdr3)):
        assert 2 <= len(cdr) <= 30, f"{clone} {name} implausible length {len(cdr)}: {cdr!r}"
    return cdr1, cdr2, cdr3


def build_aligner():
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.mode = "global"
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    return aligner


def percent_identity(aligner, a, b):
    aln = aligner.align(a, b)[0]
    t = aln.aligned  # tuple of (target_blocks, query_blocks)
    matches = 0
    total = 0
    s1, s2 = str(aln[0]), str(aln[1])
    for x, y in zip(s1, s2):
        if x == "-" or y == "-":
            total += 1
            continue
        total += 1
        if x == y:
            matches += 1
    return 100.0 * matches / total


def build_identity_matrix(aligner, seq_dict, clones):
    n = len(clones)
    mat = np.full((n, n), 100.0)
    for (i, ci), (j, cj) in combinations(enumerate(clones), 2):
        pid = percent_identity(aligner, seq_dict[ci], seq_dict[cj])
        mat[i, j] = mat[j, i] = pid
    return mat


def write_matrix_csv(mat, clones, path):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([""] + clones)
        for i, clone in enumerate(clones):
            w.writerow([clone] + [f"{v:.1f}" for v in mat[i]])


def plot_dendrogram(mat, clones, meta, chain_label, path):
    n = len(clones)
    dist = 100.0 - mat
    np.fill_diagonal(dist, 0.0)
    condensed = dist[np.triu_indices(n, k=1)]
    Z = linkage(condensed, method="average")

    targets = [meta[(c, "HC")]["target_antigen"] for c in clones]
    labels = [f"{c} ({t})" for c, t in zip(clones, targets)]

    fig, ax = plt.subplots(figsize=(8, 6))
    dendrogram(Z, labels=labels, orientation="right", ax=ax)
    ax.set_xlabel("Distance (100 - BLOSUM62 %identity)")
    ax.set_title(f"{chain_label} CDR-only similarity clustering (BLOSUM62)")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def summarize_by_target(mat, clones, meta, chain, out_lines):
    targets = {c: meta[(c, "HC")]["target_antigen"] for c in clones}
    same, diff = [], []
    for (i, ci), (j, cj) in combinations(enumerate(clones), 2):
        pid = mat[i, j]
        if targets[ci] in ("TBD", "unknown", "NA") or targets[cj] in ("TBD", "unknown", "NA"):
            continue
        (same if targets[ci] == targets[cj] else diff).append(pid)

    out_lines.append(f"\n[{chain}] Mean %identity, same-target pairs: "
                      f"{np.mean(same):.1f} (n={len(same)})" if same else f"\n[{chain}] no same-target pairs")
    out_lines.append(f"[{chain}] Mean %identity, different-target pairs: "
                      f"{np.mean(diff):.1f} (n={len(diff)})" if diff else f"[{chain}] no cross-target pairs")


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    records = parse_fasta(VAR_DOMAINS_FASTA)
    meta = load_metadata()

    clones = sorted({h.rsplit("_", 1)[0] for h in records})

    vh_cdrs, vl_cdrs = {}, {}
    with open(PROCESSED_DIR / "cdr_regions.fasta", "w") as f:
        for clone in clones:
            h1, h2, h3 = extract_vh_cdrs(records[f"{clone}_VH"], clone)
            l1, l2, l3 = extract_vl_cdrs(records[f"{clone}_VL"], clone)
            vh_cdrs[clone] = h1 + h2 + h3
            vl_cdrs[clone] = l1 + l2 + l3
            f.write(f">{clone}_VH_CDR1\n{h1}\n>{clone}_VH_CDR2\n{h2}\n>{clone}_VH_CDR3\n{h3}\n")
            f.write(f">{clone}_VL_CDR1\n{l1}\n>{clone}_VL_CDR2\n{l2}\n>{clone}_VL_CDR3\n{l3}\n")
            f.write(f">{clone}_VH_CDRs_concat\n{vh_cdrs[clone]}\n")
            f.write(f">{clone}_VL_CDRs_concat\n{vl_cdrs[clone]}\n")

    print(f"Extracted CDRs for {len(clones)} clones -> {PROCESSED_DIR / 'cdr_regions.fasta'}")

    aligner = build_aligner()

    print("Aligning VH CDRs (BLOSUM62, all pairs)...")
    vh_mat = build_identity_matrix(aligner, vh_cdrs, clones)
    print("Aligning VL CDRs (BLOSUM62, all pairs)...")
    vl_mat = build_identity_matrix(aligner, vl_cdrs, clones)

    write_matrix_csv(vh_mat, clones, RESULTS_DIR / "cdr_identity_matrix_VH.csv")
    write_matrix_csv(vl_mat, clones, RESULTS_DIR / "cdr_identity_matrix_VL.csv")

    plot_dendrogram(vh_mat, clones, meta, "VH", RESULTS_DIR / "cdr_dendrogram_VH.png")
    plot_dendrogram(vl_mat, clones, meta, "VL", RESULTS_DIR / "cdr_dendrogram_VL.png")

    summary = ["CDR-only BLOSUM62 similarity summary", "=" * 45]
    summarize_by_target(vh_mat, clones, meta, "VH", summary)
    summarize_by_target(vl_mat, clones, meta, "VL", summary)
    summary_text = "\n".join(summary) + "\n"
    (RESULTS_DIR / "cdr_similarity_summary.txt").write_text(summary_text)
    print(summary_text)

    print(f"Wrote outputs to {RESULTS_DIR}")


if __name__ == "__main__":
    main()
