#!/usr/bin/env python3
"""
Assign the closest-matching rabbit IGHV germline gene/allele to each clone's VH
domain, by BLASTP against IMGT's rabbit (Oryctolagus cuniculus) germline V-REGION
reference set.

This is a nearest-germline lookup (BLASTP top hit by %identity over the aligned
region), not a validated IgBLAST/V(D)J assignment -- no D/J calling, no CDR3
junction analysis, no correction for the germline reference covering only
FR1-FR3 (not CDR3/FR4, which derive from D/J segments and somatic diversity).
Treat the result as "closest known germline", not a certified V-gene call.

Requires `blastp`/`makeblastdb` (installed via `brew install blast`) and the
IMGT rabbit IGHV reference fetched into
data/raw/germline/imgt_rabbit_IGHV_raw.txt (raw HTML/text response from
IMGT/GENE-DB's GENElect endpoint, query=7.2+IGHV, species=Oryctolagus cuniculus).
"""
import csv
import re
import subprocess
import sys
import tempfile
from pathlib import Path

from Bio.Seq import Seq

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cdr_similarity_analysis as base  # noqa: E402

ROOT = base.ROOT
GERMLINE_RAW = ROOT / "data/raw/germline/imgt_rabbit_IGHV_raw.txt"
GERMLINE_DIR = ROOT / "data/processed/germline"
RESULTS_DIR = ROOT / "results/germline"

FASTA_RECORD_RE = re.compile(
    r">(\S+)\|([^|]*)\|Oryctolagus cuniculus[^|]*\|([^|]*)\|V-REGION\|[^\n]*\n([^>]+)"
)


def parse_imgt_dump(path):
    """Returns list of (allele_name, functionality, nt_seq) from IMGT/GENE-DB raw output."""
    text = path.read_text()
    records = []
    for m in FASTA_RECORD_RE.finditer(text):
        accession, allele, functionality, seq_block = m.groups()
        nt_seq = re.sub(r"[^ACGTNacgtn]", "", seq_block)
        if nt_seq:
            records.append((allele, functionality.strip() or "?", nt_seq.upper()))
    return records


def translate_to_protein(nt_seq):
    trimmed = nt_seq[: len(nt_seq) - (len(nt_seq) % 3)]
    return str(Seq(trimmed).translate(to_stop=False)).replace("*", "X")


def build_germline_protein_fasta(records, path):
    with open(path, "w") as f:
        for allele, functionality, nt_seq in records:
            protein = translate_to_protein(nt_seq)
            f.write(f">{allele}|{functionality}\n{protein}\n")


def run_blastp(query_fasta, db_fasta, db_dir):
    subprocess.run(
        ["makeblastdb", "-in", str(db_fasta), "-dbtype", "prot", "-out", str(db_dir / "rabbit_ighv")],
        capture_output=True, text=True, check=True,
    )
    result = subprocess.run(
        ["blastp", "-query", str(query_fasta), "-db", str(db_dir / "rabbit_ighv"),
         "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
         "-max_target_seqs", "3"],
        capture_output=True, text=True, check=True,
    )
    return result.stdout


def best_hits_per_query(blast_tsv):
    best = {}
    for line in blast_tsv.strip().splitlines():
        if not line:
            continue
        fields = line.split("\t")
        qseqid, sseqid, pident = fields[0], fields[1], float(fields[2])
        bitscore = float(fields[-1])
        if qseqid not in best or bitscore > best[qseqid][2]:
            best[qseqid] = (sseqid, pident, bitscore)
    return best


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    GERMLINE_DIR.mkdir(parents=True, exist_ok=True)

    germline_records = parse_imgt_dump(GERMLINE_RAW)
    print(f"Parsed {len(germline_records)} rabbit IGHV germline alleles from IMGT dump.")
    assert germline_records, "No germline records parsed -- check IMGT raw dump format"

    germline_protein_fasta = GERMLINE_DIR / "rabbit_IGHV_germline_protein.fasta"
    build_germline_protein_fasta(germline_records, germline_protein_fasta)

    var_records = base.parse_fasta(base.VAR_DOMAINS_FASTA)
    clones = sorted({h.rsplit("_", 1)[0] for h in var_records if h.endswith("_VH")})
    query_fasta = RESULTS_DIR / "query_VH.fasta"
    with open(query_fasta, "w") as f:
        for clone in clones:
            f.write(f">{clone}\n{var_records[f'{clone}_VH']}\n")

    blast_tsv = run_blastp(query_fasta, germline_protein_fasta, GERMLINE_DIR)
    (RESULTS_DIR / "blastp_vs_rabbit_IGHV_raw.tsv").write_text(blast_tsv)

    best = best_hits_per_query(blast_tsv)
    functionality = {allele: func for allele, func, _ in germline_records}

    with open(RESULTS_DIR / "vh_germline_assignment.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["clone", "closest_IGHV_allele", "functionality", "pident_pct", "bitscore"])
        for clone in clones:
            if clone not in best:
                w.writerow([clone, "no_hit", "", "", ""])
                continue
            sseqid, pident, bitscore = best[clone]
            allele = sseqid.split("|")[0]
            w.writerow([clone, allele, functionality.get(allele, "?"), f"{pident:.1f}", f"{bitscore:.0f}"])

    print(f"Wrote {RESULTS_DIR / 'vh_germline_assignment.csv'}")
    for clone in clones:
        if clone in best:
            sseqid, pident, _ = best[clone]
            print(f"  {clone}: {sseqid.split('|')[0]} ({pident:.1f}% id)")


if __name__ == "__main__":
    main()
