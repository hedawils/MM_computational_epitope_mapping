# Antibody sequences

Source: Cyr MG, Wilson HD, Spierling AL, Chang J, Peng H, Steinberger P, Rader C.
"Concerted Antibody and Antigen Discovery by Differential Whole-cell Phage Display
Selections and Multi-omic Target Deconvolution." J Mol Biol. 2023;435(10):168085.
doi: [10.1016/j.jmb.2023.168085](https://doi.org/10.1016/j.jmb.2023.168085), PMID: 37019174.

Sequences provided directly by the author (H. Wilson), not scraped from the published text.

## Contents

- `antibody_sequences_full_constructs.fasta` — full expressed chimeric rabbit/human
  IgG1 constructs as pasted by the author: signal peptide + rabbit VH/VL + human
  constant domains. Heavy chains carry the L234A/L235A/P329G (LALA-PG) Fc-silencing
  mutations. 15 clones x 2 chains (HC/LC) = 30 records.
- `antibody_metadata.csv` — per-record clone/chain/isotype/length summary generated
  by [../../scripts/process_antibody_sequences.py](../../scripts/process_antibody_sequences.py).

## Clones

HW-1, HW-8, HW-16, HW-17, HW-25, HW-28, HW-42, HW-45, HW-56, HW-70, HW-81, HW-97,
HW-101, HW-113, TH9-022 (15 clones total — the user's original count of 12 did not
match the pasted data; flagged here rather than silently dropped).

**Target antigens** (`target_antigen` in `antibody_metadata.csv`) are from an
author-provided clone/target table (PTPRG, GARS, CADM1, ICAM1). Cross-checked
against light-chain isotype (kappa vs lambda) shown in that table, which matched
the isotype implied by each LC sequence's C-terminal residues (`...FNRGEC` = kappa,
`...TVAPTECS` = lambda) for all 15 clones. **HW-8 and TH9-022 targets are still
`TBD`** — cropped out of the source screenshot; update once confirmed.

## Signal peptide

All 30 records share the same 19-residue N-terminal signal peptide
(`MDWTWRILFLVAAATGAHS`), confirmed programmatically for every record. This is
cleaved during secretion and is not part of the folded protein, so it is stripped
in the derived "mature" sequences used for structure prediction (see
[../../data/processed](../../data/processed)).
