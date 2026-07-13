# MM Computational Epitope Mapping

Structural prediction and epitope mapping of antibody sequences.

## Project structure

- `data/raw/` — input antibody sequences (FASTA, VH/VL, CDR annotations)
- `data/processed/` — cleaned/formatted sequence data
- `structures/predicted/` — predicted antibody/complex structures (e.g. AlphaFold, ABodyBuilder, IgFold output)
- `structures/reference/` — reference/experimental structures (PDB)
- `scripts/` — analysis and pipeline scripts
- `notebooks/` — exploratory analysis notebooks
- `results/` — epitope mapping outputs, figures, summary tables
- `docs/` — notes and documentation

## Scope

Pipeline for predicting antibody structures from sequence and mapping antigen epitopes.
