# MM Computational Epitope Mapping

Computational structural prediction and epitope mapping of the rabbit/human
chimeric antibody panel from:

> Cyr MG, Wilson HD, Spierling AL, Chang J, Peng H, Steinberger P, Rader C.
> **Concerted Antibody and Antigen Discovery by Differential Whole-cell Phage
> Display Selections and Multi-omic Target Deconvolution.**
> *J Mol Biol.* 2023;435(10):168085.
> doi: [10.1016/j.jmb.2023.168085](https://doi.org/10.1016/j.jmb.2023.168085),
> PMID: [37019174](https://pubmed.ncbi.nlm.nih.gov/37019174/), PMCID: PMC10148915

The paper used whole-cell phage display against multiple myeloma cells plus
multi-omic target deconvolution to generate a panel of >50 unique mAbs and
identify three validated cell-surface antigen targets — PTPRG, ICAM1, and
CADM1 — with epitope binning done experimentally (paper Figure 3b). This
repository picks up where that paper's sequence-level characterization left
off, working computationally from the antibody sequences toward their
structures and binding footprints.

## Objectives

1. **Where do these antibodies bind?** Predict antibody structures and
   antibody–antigen complexes to localize each clone's binding footprint
   (epitope) on its target antigen.
2. **How diverse are those epitopes?** Characterize epitope diversity within
   and across target antigens — do same-target clones converge on one epitope
   or several, and how does that relate to the sequence-level diversity
   already observed (see [Status](#status) below)?
3. **Are these clones therapeutically developable?** Profile the panel for
   standard developability liabilities (aggregation propensity, hydrophobicity
   and charge patches, sequence liability motifs — deamidation, oxidation,
   isomerization, glycosylation — CDR length outliers, etc.).
4. **How could they be improved?** Use the structural, epitope, and
   developability findings to propose concrete iterative engineering designs
   (affinity maturation, liability removal, framework optimization) for
   individual clones.

## Status

What's done so far, and what each piece feeds into:

| Objective | Done | Not yet done |
|---|---|---|
| Sequence foundation | 15 clones' full chimeric HC/LC sequences ([data/raw/sequences](data/raw/sequences)), mature/variable-domain/CDR extraction ([data/processed](data/processed)) | — |
| 2. Epitope diversity (sequence proxy) | VH/VL and CDR-only similarity clustering, BLOSUM62-scored, same-target vs. different-target comparison ([results/sequence_similarity](results/sequence_similarity), [results/cdr_similarity](results/cdr_similarity)); nearest rabbit IGHV germline gene per clone ([results/germline](results/germline)) | Structural/epitope-based diversity (see objective 1) |
| 1. Where do they bind | — | Antigen structures (PTPRG, ICAM1, CADM1 extracellular domains — note: the target/clone table also includes **GARS**, a fourth target not mentioned in the paper's abstract, worth resolving before proceeding); antibody structure prediction (e.g. IgFold/ABodyBuilder for Fv, AlphaFold-Multimer/AF3 for antibody-antigen complexes); binding footprint/epitope residue identification from predicted complexes; comparison against the paper's experimental epitope binning (Fig 3b) if that data can be recovered |
| 3. Developability | — | Liability motif scanning, aggregation/hydrophobicity/charge patch prediction (e.g. TAP-style metrics), CDR length/composition outliers |
| 4. Enhancement designs | — | Depends on 1–3 |

Target antigen assignment is incomplete for 2 of 15 clones (HW-8: unknown,
TH9-022: NA — see [data/raw/sequences/README.md](data/raw/sequences/README.md))
and should be resolved before epitope-diversity conclusions are drawn per target.

## Project structure

- `data/raw/` — input antibody sequences, target/germline reference downloads
- `data/processed/` — derived sequence data (mature chains, variable domains, CDRs, alignments)
- `structures/predicted/` — predicted antibody/complex structures (e.g. AlphaFold, ABodyBuilder, IgFold output)
- `structures/reference/` — reference/experimental structures (PDB), incl. target antigens once sourced
- `scripts/` — analysis and pipeline scripts
- `notebooks/` — exploratory analysis notebooks
- `results/` — sequence similarity, CDR/germline analysis, and (pending) epitope mapping and developability outputs
- `docs/` — notes and documentation
