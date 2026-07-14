#!/usr/bin/env python3
"""
CDR3-only dendrogram (BLOSUM62 distance, leaf labels colored by target antigen)
with a right-hand panel showing the actual aligned CDR3 sequences.

Per chain (VH, VL): CDR3 sequences are aligned across all 15 clones with MAFFT
(L-INS-i). Clustering (BLOSUM62 identity matrix + linkage) is also computed on
CDR3 alone, since it is the most diagnostic loop for antigen specificity.

Requires the `mafft` binary on PATH (installed via `brew install mafft`).
Reuses CDR extraction + BLOSUM62 identity matrix / linkage from
cdr_similarity_analysis.py.
"""
import csv
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cdr_similarity_analysis as base  # noqa: E402

ROOT = base.ROOT
RESULTS_DIR = base.RESULTS_DIR
ALIGN_DIR = ROOT / "data/processed/cdr_alignments"

TARGET_COLORS = {
    "PTPRG": "tab:blue",
    "GARS": "tab:orange",
    "CADM1": "tab:green",
    "ICAM1": "tab:red",
    "NA": "0.5",
    "unknown": "0.5",
}

RESIDUE_GROUPS = {
    "hydrophobic (A,V,L,I,M,C)": ("AVLIMC", "#7EA6E0"),
    "aromatic (F,W,Y,H)": ("FWYH", "#B07EE0"),
    "polar (S,T,N,Q)": ("STNQ", "#7ED6A0"),
    "positive (K,R)": ("KR", "#E07E7E"),
    "negative (D,E)": ("DE", "#E0A87E"),
    "special (G,P)": ("GP", "#E0D97E"),
}
GAP_COLOR = "#EEEEEE"


def residue_color(aa):
    if aa in ("-", ""):
        return GAP_COLOR
    for members, color in RESIDUE_GROUPS.values():
        if aa in members:
            return color
    return "#DDDDDD"


def run_mafft(seqs):
    """seqs: {clone: ungapped_sequence} -> {clone: aligned_sequence}"""
    with tempfile.TemporaryDirectory() as tmpdir:
        in_path = Path(tmpdir) / "in.fasta"
        with open(in_path, "w") as f:
            for clone, seq in seqs.items():
                f.write(f">{clone}\n{seq}\n")
        result = subprocess.run(
            ["mafft", "--amino", "--localpair", "--maxiterate", "1000", "--quiet", str(in_path)],
            capture_output=True, text=True, check=True,
        )
    aligned = _parse_fasta_str(result.stdout)
    assert set(aligned) == set(seqs), "MAFFT output headers don't match input clones"
    return aligned


def _parse_fasta_str(text):
    records = {}
    header, seq_lines = None, []
    for line in text.splitlines():
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


def build_chain_alignment(clones, cdr3, chain_label):
    ALIGN_DIR.mkdir(parents=True, exist_ok=True)
    aln3 = run_mafft(cdr3)

    with open(ALIGN_DIR / f"{chain_label}_CDR3_aligned.fasta", "w") as f:
        for clone in clones:
            f.write(f">{clone}\n{aln3[clone]}\n")

    return aln3


def plot_combined(mat, clones, meta, chain_label, aligned_rows, path, germline=None):
    n = len(clones)
    dist = 100.0 - mat
    np.fill_diagonal(dist, 0.0)
    condensed = dist[np.triu_indices(n, k=1)]
    Z = linkage(condensed, method="average")

    targets = {c: meta[(c, "HC")]["target_antigen"] for c in clones}
    if germline:
        labels = [f"{c} ({targets[c]})  {germline.get(c, '')}" for c in clones]
    else:
        labels = [f"{c} ({targets[c]})" for c in clones]

    n_cols = len(next(iter(aligned_rows.values())))
    fig_w = 8 + n_cols * 0.16
    fig, (ax_d, ax_a) = plt.subplots(
        1, 2, figsize=(fig_w, 6),
        gridspec_kw={"width_ratios": [3, max(2.0, n_cols * 0.12)]},
    )

    R = dendrogram(Z, labels=labels, orientation="right", ax=ax_d,
                    link_color_func=lambda k: "black")
    ax_d.set_xlabel("Distance (100 - BLOSUM62 %identity)")
    ax_d.set_title(f"{chain_label} CDR clustering")

    ordered_clones = [lbl.split(" (")[0] for lbl in R["ivl"]]
    yticks = ax_d.get_yticks()
    for label_obj, clone in zip(ax_d.get_yticklabels(), ordered_clones):
        label_obj.set_color(TARGET_COLORS.get(targets[clone], "black"))

    row_h = yticks[1] - yticks[0] if len(yticks) > 1 else 10
    for tick, clone in zip(yticks, ordered_clones):
        seq = aligned_rows[clone]
        for j, aa in enumerate(seq):
            if aa == " ":
                continue
            ax_a.add_patch(Rectangle((j, tick - row_h * 0.4), 1, row_h * 0.8,
                                      facecolor=residue_color(aa), edgecolor="none"))
            ax_a.text(j + 0.5, tick, aa, ha="center", va="center", fontsize=6, family="monospace")

    ax_a.set_xlim(0, n_cols)
    ax_a.set_ylim(ax_d.get_ylim())
    ax_a.set_yticks([])
    ax_a.set_xticks([])
    ax_a.set_title(f"{chain_label} CDR3 alignment (MAFFT)")

    color_to_labels = {}
    for t, c in TARGET_COLORS.items():
        color_to_labels.setdefault(c, []).append(t)
    target_legend = [Patch(facecolor="none", edgecolor="none", label="Target antigen (label color):")]
    target_legend += [Patch(facecolor=c, label=" / ".join(labels)) for c, labels in color_to_labels.items()]
    residue_legend = [Patch(facecolor=c, label=k) for k, (_, c) in RESIDUE_GROUPS.items()]
    residue_legend.append(Patch(facecolor=GAP_COLOR, label="gap"))
    fig.legend(handles=target_legend, loc="lower left", ncol=3, fontsize=7, frameon=False,
               bbox_to_anchor=(0.02, -0.02))
    fig.legend(handles=residue_legend, loc="lower right", ncol=3, fontsize=7, frameon=False,
               bbox_to_anchor=(0.98, -0.02))

    fig.tight_layout(rect=[0, 0.06, 1, 1])
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def load_vh_germline():
    path = ROOT / "results/germline/vh_germline_assignment.csv"
    if not path.exists():
        return None
    germline = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            germline[row["clone"]] = f"{row['closest_IGHV_allele']} ({row['pident_pct']}% id)"
    return germline


def main():
    records = base.parse_fasta(base.VAR_DOMAINS_FASTA)
    meta = base.load_metadata()
    clones = sorted({h.rsplit("_", 1)[0] for h in records})
    vh_germline = load_vh_germline()

    for chain, extract_fn in (("VH", base.extract_vh_cdrs), ("VL", base.extract_vl_cdrs)):
        cdr3 = {}
        for clone in clones:
            _, _, c3 = extract_fn(records[f"{clone}_{chain}"], clone)
            cdr3[clone] = c3

        print(f"[{chain}] running MAFFT CDR3 alignment...")
        aligned_rows = build_chain_alignment(clones, cdr3, chain)

        print(f"[{chain}] computing BLOSUM62 CDR3 identity matrix + linkage...")
        aligner = base.build_aligner()
        mat = base.build_identity_matrix(aligner, cdr3, clones)

        out_path = RESULTS_DIR / f"cdr3_dendrogram_alignment_{chain}.png"
        germline = vh_germline if chain == "VH" else None
        plot_combined(mat, clones, meta, chain, aligned_rows, out_path, germline=germline)
        print(f"[{chain}] wrote {out_path}")


if __name__ == "__main__":
    main()
