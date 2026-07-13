#!/usr/bin/env python3
"""
Sequence alignment / similarity analysis across the 15 rabbit/human chimeric clones.

Steps:
  1. Load mature HC/LC sequences (data/processed/antibody_sequences_mature.fasta)
     and clone metadata (data/raw/sequences/antibody_metadata.csv).
  2. Extract the rabbit variable domains (VH, VL) by stripping the shared human
     constant region, found as the longest common suffix within each chain/isotype
     group (HC; LC-kappa; LC-lambda) rather than assumed from prior knowledge.
  3. Pairwise global alignment (Needleman-Wunsch, linear gap penalty) of VH-vs-VH
     and VL-vs-VL for all clone pairs -> percent identity matrices.
  4. Hierarchical clustering + dendrograms of VH and VL identity, labeled by clone
     and target antigen.
  5. Summary stats: mean pairwise identity for clone pairs sharing a target vs.
     pairs with different targets.

No third-party alignment libraries required (pure Python Needleman-Wunsch).
"""
import csv
from pathlib import Path
from itertools import combinations

import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
MATURE_FASTA = ROOT / "data/processed/antibody_sequences_mature.fasta"
METADATA_CSV = ROOT / "data/raw/sequences/antibody_metadata.csv"
PROCESSED_DIR = ROOT / "data/processed"
RESULTS_DIR = ROOT / "results/sequence_similarity"

MATCH, MISMATCH, GAP = 1, -1, -2


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


def longest_common_suffix(seqs):
    if not seqs:
        return ""
    shortest = min(seqs, key=len)
    for i in range(len(shortest), 0, -1):
        suffix = shortest[-i:]
        if all(s.endswith(suffix) for s in seqs):
            return suffix
    return ""


def extract_variable_domains(records, meta):
    """Return dict: (clone, chain) -> variable-domain sequence."""
    hc_seqs = [seq for h, seq in records.items() if h.endswith("_HC")]
    hc_suffix = longest_common_suffix(hc_seqs)

    kappa_clones = {c for (c, ch), r in meta.items() if ch == "LC" and "kappa" in r["isotype"]}
    lambda_clones = {c for (c, ch), r in meta.items() if ch == "LC" and "lambda" in r["isotype"]}

    kappa_seqs = [records[f"{c}_LC"] for c in kappa_clones]
    lambda_seqs = [records[f"{c}_LC"] for c in lambda_clones]
    kappa_suffix = longest_common_suffix(kappa_seqs)
    lambda_suffix = longest_common_suffix(lambda_seqs)

    variable = {}
    for header, seq in records.items():
        clone, chain = header.rsplit("_", 1)
        if chain == "HC":
            variable[(clone, "VH")] = seq[: len(seq) - len(hc_suffix)]
        else:
            suffix = kappa_suffix if clone in kappa_clones else lambda_suffix
            variable[(clone, "VL")] = seq[: len(seq) - len(suffix)]

    return variable, hc_suffix, kappa_suffix, lambda_suffix


def needleman_wunsch_identity(a, b):
    """Global alignment, linear gap penalty. Returns percent identity (0-100)."""
    n, m = len(a), len(b)
    score = np.zeros((n + 1, m + 1), dtype=int)
    score[:, 0] = np.arange(n + 1) * GAP
    score[0, :] = np.arange(m + 1) * GAP

    for i in range(1, n + 1):
        ai = a[i - 1]
        row_prev = score[i - 1]
        row_cur = score[i]
        for j in range(1, m + 1):
            diag = row_prev[j - 1] + (MATCH if ai == b[j - 1] else MISMATCH)
            up = row_prev[j] + GAP
            left = row_cur[j - 1] + GAP
            row_cur[j] = max(diag, up, left)

    # traceback to count identical aligned positions
    i, j = n, m
    matches, aligned_len = 0, 0
    while i > 0 and j > 0:
        cur = score[i, j]
        diag = score[i - 1, j - 1] + (MATCH if a[i - 1] == b[j - 1] else MISMATCH)
        if cur == diag:
            if a[i - 1] == b[j - 1]:
                matches += 1
            aligned_len += 1
            i, j = i - 1, j - 1
        elif cur == score[i - 1, j] + GAP:
            aligned_len += 1
            i -= 1
        else:
            aligned_len += 1
            j -= 1
    aligned_len += i + j  # remaining leading gaps
    return 100.0 * matches / aligned_len


def build_identity_matrix(seq_dict, clones):
    n = len(clones)
    mat = np.full((n, n), 100.0)
    for (i, ci), (j, cj) in combinations(enumerate(clones), 2):
        pid = needleman_wunsch_identity(seq_dict[ci], seq_dict[cj])
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
    ax.set_xlabel("Distance (100 - %identity)")
    ax.set_title(f"{chain_label} variable-domain similarity clustering")
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
    records = parse_fasta(MATURE_FASTA)
    meta = load_metadata()

    variable, hc_suffix, kappa_suffix, lambda_suffix = extract_variable_domains(records, meta)

    clones = sorted({c for c, _ in variable})

    with open(PROCESSED_DIR / "variable_domains.fasta", "w") as f:
        for clone in clones:
            f.write(f">{clone}_VH\n{variable[(clone, 'VH')]}\n")
            f.write(f">{clone}_VL\n{variable[(clone, 'VL')]}\n")

    print(f"HC constant-region suffix ({len(hc_suffix)} aa): {hc_suffix[:40]}...")
    print(f"LC kappa constant-region suffix ({len(kappa_suffix)} aa): {kappa_suffix[:40]}...")
    print(f"LC lambda constant-region suffix ({len(lambda_suffix)} aa): {lambda_suffix[:40]}...")

    vh_seqs = {c: variable[(c, "VH")] for c in clones}
    vl_seqs = {c: variable[(c, "VL")] for c in clones}

    print("Aligning VH domains (all pairs)...")
    vh_mat = build_identity_matrix(vh_seqs, clones)
    print("Aligning VL domains (all pairs)...")
    vl_mat = build_identity_matrix(vl_seqs, clones)

    write_matrix_csv(vh_mat, clones, RESULTS_DIR / "identity_matrix_VH.csv")
    write_matrix_csv(vl_mat, clones, RESULTS_DIR / "identity_matrix_VL.csv")

    plot_dendrogram(vh_mat, clones, meta, "VH", RESULTS_DIR / "dendrogram_VH.png")
    plot_dendrogram(vl_mat, clones, meta, "VL", RESULTS_DIR / "dendrogram_VL.png")

    summary = ["Sequence similarity summary", "=" * 40]
    summarize_by_target(vh_mat, clones, meta, "VH", summary)
    summarize_by_target(vl_mat, clones, meta, "VL", summary)
    summary_text = "\n".join(summary) + "\n"
    (RESULTS_DIR / "similarity_summary.txt").write_text(summary_text)
    print(summary_text)

    print(f"Wrote outputs to {RESULTS_DIR}")


if __name__ == "__main__":
    main()
