# GPR161 deep mutational scan

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23682577.svg)](https://doi.org/10.5281/zenodo.23682577)

Code and data for the GPR161 saturation mutagenesis study. A library built for this work
with the DIMPLE method carries **11,043 designed variants** across **528 positions**
(missense, synonymous and single-residue deletions), read out in **three parallel FACS
assays** — cell-surface expression, whole-cell abundance, and basal cAMP signalling — in
HEK293T cells.

Everything here rebuilds from the files in `data/`: the variant tables, the position
table, the QC tables and every data panel of the figures.

> **Scope.** This repository currently carries what is needed to reproduce
> GPR161-related work. The pooled human GPCRome basal-cAMP survey that accompanies this
> work has its own repository:
> <https://github.com/matthewkarlhoward/2026_GPCRome_basal_cAMP>.

---

## Repository layout

```
data/
  library/      oligo pool (11,043 designed variants) and the FACS/sequencing sample layout
  counts/       per-sample variant counts, 78 samples, gzipped CSV
  scores/       Lilace output, one table per assay (cre / surface / abundance)
  processed/    merged variant table + position table, built by pipeline/
  annotation/   GPCRdb snake-plot layout, generic numbers, and the 8SMV structural join
  structure/    the 8SMV coordinates (fetched) and the ChimeraX attribute files
pipeline/       counts -> scores -> analysis tables
  dumpling/     the read-processing run: Snakemake config and construct reference
scripts/        figure and analysis code
  qc/           library and sequencing QC (the QC supplement)
  structure/    the ChimeraX renders behind Figure 3
figures/        figure panels written by the scripts (PDF + PNG)
tables/         QC tables written by the scripts (CSV)
docs/           METHODS.md, DATA_SOURCES.md, FIGURES.md, COLUMNS.md, RELEASING.md
```

## Rebuilding everything

`bash run_all.sh` does all of this. Step by step, with manuscript figure numbers:

```bash
pip install -r requirements.txt

python3 pipeline/01_merge_scores.py          # data/processed/merged_variant_scores.tsv
python3 pipeline/02_position_summary.py      # data/processed/position_summary.tsv
python3 pipeline/03_structure_annotation.py  # data/annotation/structure_8smv.tsv

python3 scripts/fig1_landscape.py          # Fig 2b,c  density + substitution heatmaps
python3 scripts/fig2_phenotype_classes.py  # Fig 3b-e  classes, snake plot, Sankey
python3 scripts/fig3_structure.py          # Fig 4b    mean surface effect by segment
python3 scripts/fig4_human_variants.py     # Fig 5a,b  gnomAD and ClinVar alleles
python3 scripts/qc/dms_qc_panel.py         # Supp Fig 3, panels a-d in one row

python3 scripts/qc/dms_qc_coverage.py      # write tables/qc_*.csv; add --figures
python3 scripts/qc/dms_qc_scores.py        # for larger standalone QC panels
```

Two scripts are deliberately **not** in that list, and not in `run_all.sh`:

- `pipeline/04_human_alleles.py` rebuilds `data/annotation/human_alleles.tsv`, the input
  to Figure 5. It needs a clinical-variant crosswalk that is not redistributed here (it
  carries case-level clinical fields), so it cannot run from a fresh clone. The table it
  produces is committed instead, and the script records the derivation.
Manuscript Figure 4 panels a and c are ChimeraX renders made by hand; see below.

Panel e of the QC supplement is the FACS gating, exported from the cytometry software,
and panel a of Figure 1 is the assay schematic — neither is drawn by code. Everything
else in the figures comes from the scripts above.

### Figure 4 panels a and c

Those two panels are ChimeraX renders made by hand. The structure is
`data/structure/8smv.cif`, and the colouring is the per-position surface and cAMP values
in `data/processed/position_summary.tsv` (`surf`, `cre`, `resid`), applied as residue
attributes in a ChimeraX session kept with the manuscript figures rather than here.
Neither that session nor the render images is deposited, so those two panels cannot be
rebuilt from this repository — the structure file and the per-position table are what
they were built from.

Panel b (`figures/segment_effects.pdf`) is drawn by `scripts/fig3_structure.py` and
needs no ChimeraX. Nothing in `run_all.sh` does.

`pipeline/01` and `pipeline/02` regenerate the two tables in `data/processed/` byte for
byte; they are committed so the figure scripts run without rerunning the pipeline.

The two upstream steps — FASTQ to per-variant counts, and the model fit that turns
counts into scores — were run with third-party software rather than code in this
repository. Their outputs are this study's data and are committed here; the tools, their
inputs and what still needs recording are in
[`pipeline/README.md`](pipeline/README.md).

## What is measured

| assay | readout | directory key |
|---|---|---|
| **Cell-surface expression** | anti-FLAG stain of the N-terminally tagged receptor | `surface` |
| **Whole-cell abundance** | split-mNeonGreen complementation | `abundance` |
| **Basal cAMP signalling** | CRE-eGFP transcriptional reporter | `cre` |

Each assay: three biological replicates, four FACS bins per replicate, plus three
unsorted library-input samples. Variant effects were estimated from the bin counts with
[Lilace](https://doi.org/10.1186/s13059-026-03934-1), a Bayesian hierarchical model that
returns a posterior mean effect, its standard error, and a local false sign rate (lfsr)
per variant.

## Numbering

**All positions are UniProt [Q8N6U8](https://www.uniprot.org/uniprotkb/Q8N6U8) numbering**
(canonical human GPR161, 529 aa), which is also the numbering of RefSeq NP_001362812.1
and of chain R in PDB 8SMV. Positions 2-529 are assayed.

> **Trap.** Part of the GPR161 literature numbers against the 549-aa isoform, which runs
> **+20** relative to the canonical sequence, and the mouse work adds its own offset.
> Any literature variant must be mapped before it is compared with these scores.
> See [`docs/METHODS.md`](docs/METHODS.md#6-numbering-and-mapping).

## References

**Library design** — DIMPLE, the method the variants in `data/library/` were designed and
cloned with.

> Macdonald, C. B., Nedrud, D., Grimes, P. R., Trinidad, D., Fraser, J. S. &
> Coyote-Maestas, W. DIMPLE: deep insertion, deletion, and missense mutation libraries for
> exploring protein variation in evolution, disease, and biology. *Genome Biology* **24**,
> 36 (2023). doi:[10.1186/s13059-023-02880-6](https://doi.org/10.1186/s13059-023-02880-6)

**Variant-effect scoring** — Lilace, the model fitted to this study's bin counts to
produce `data/scores/`.

> Freudenberg, J., Rao, J., Howard, M. K., Macdonald, C., Greenwald, N. F.,
> Coyote-Maestas, W. & Pimentel, H. Accurate variant effect estimation in FACS-based deep
> mutational scanning data with Lilace. *Genome Biology* **27**, 48 (2026).
> doi:[10.1186/s13059-026-03934-1](https://doi.org/10.1186/s13059-026-03934-1)

**Topology diagram and generic numbering** — third-party data redistributed in
`data/annotation/`, used for the snake plot in Figure 2 and the `bw` column of the
position table. Cite all three if the generic numbers or the snake plot appear in a
figure; the first two are what GPCRdb asks for.

> Taracena Herrera, L. P. *et al.* GPCRdb in 2025: adding odorant receptors, data mapper,
> structure similarity search and models of physiological ligand complexes. *Nucleic Acids
> Research* **53**, D425-D435 (2025). doi:[10.1093/nar/gkae1065](https://doi.org/10.1093/nar/gkae1065)

> Isberg, V. *et al.* Generic GPCR residue numbers — aligning topology maps while minding
> the gaps. *Trends in Pharmacological Sciences* **36**, 22-31 (2015).
> doi:[10.1016/j.tips.2014.11.001](https://doi.org/10.1016/j.tips.2014.11.001)

> Ballesteros, J. A. & Weinstein, H. Integrated methods for the construction of
> three-dimensional models and computational probing of structure-function relations in G
> protein-coupled receptors. *Methods in Neurosciences* **25**, 366-428 (1995).

**Sequence numbering** — every position in this repository is defined against these

> UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2025. *Nucleic Acids
> Research* **53**, D609-D617 (2025). doi:[10.1093/nar/gkae1010](https://doi.org/10.1093/nar/gkae1010)
> — entry [Q8N6U8](https://www.uniprot.org/uniprotkb/Q8N6U8), 529 aa, sequence version 1

> O'Leary, N. A. *et al.* Reference sequence (RefSeq) database at NCBI. *Nucleic Acids
> Research* **44**, D733-D745 (2016). doi:[10.1093/nar/gkv1189](https://doi.org/10.1093/nar/gkv1189)
> · Morales, J. *et al.* A joint NCBI and EMBL-EBI transcript set for clinical genomics and
> research. *Nature* **604**, 310-315 (2022).
> doi:[10.1038/s41586-022-04558-8](https://doi.org/10.1038/s41586-022-04558-8)
> — MANE Select NM_001375883.1 / NP_001362812.1

**Structure** used to interpret positions, and the source of the published point
mutants this scan is benchmarked against

> Hoppe, N. *et al.* GPR161 structure uncovers the redundant role of sterol-regulated
> ciliary cAMP signaling in the Hedgehog pathway. *Nature Structural & Molecular Biology*
> **31**, 667-677 (2024). doi:[10.1038/s41594-024-01223-8](https://doi.org/10.1038/s41594-024-01223-8)
> — PDB [8SMV](https://www.rcsb.org/structure/8SMV), chain R numbered as UniProt

**Software** — `numpy`, `pandas`, `scipy`, `matplotlib`, `pillow`; versions pinned in
[`requirements.txt`](requirements.txt).

**Read processing** — Dumpling, the pipeline that turned this study's FASTQ into
`data/counts/`. Its configuration and construct reference are committed at
`pipeline/dumpling/`.

> Rao, J., Xin, R., Macdonald, C., Howard, M. K., Estevam, G. O., Yee, S. W., Wang, M.,
> Fraser, J. S., Coyote-Maestas, W. & Pimentel, H. Rosace: a robust deep mutational
> scanning analysis framework employing position and mean-variance shrinkage. *Genome
> Biology* **25**, 138 (2024).
> doi:[10.1186/s13059-024-03279-7](https://doi.org/10.1186/s13059-024-03279-7)
> — the paper that released [Dumpling](https://github.com/odcambc/dumpling), the
> Snakemake pipeline used here

## Citing this repository

Releases are archived on Zenodo. Cite the **concept DOI**, which always resolves to the
latest version:

> Howard, M. K. GPR161 deep mutational scan. Zenodo.
> doi:[10.5281/zenodo.23682577](https://doi.org/10.5281/zenodo.23682577)

Each release also gets its own version DOI (v0.1.0 is
[10.5281/zenodo.23682578](https://doi.org/10.5281/zenodo.23682578)). Citation metadata is
in [`CITATION.cff`](CITATION.cff), which GitHub renders as **Cite this repository**, and
in [`.zenodo.json`](.zenodo.json); the release procedure is in
[`docs/RELEASING.md`](docs/RELEASING.md).

## Correspondence

Matthew K. Howard — <matthewkarlhoward@gmail.com>

## Licence

Copyright © 2026 Matthew K. Howard.

- **Code** (`pipeline/`, `scripts/`) — MIT, see [`LICENSE`](LICENSE).
- **Everything else** — data, figures, tables, documentation and this README —
  [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/), see
  [`LICENSE-DATA`](LICENSE-DATA). Please cite the paper.
- **Third-party content, not covered by that copyright**: `data/annotation/` is
  redistributed from GPCRdb under CC BY 4.0, and two files carry it onward — the `bw`
  column of `data/processed/position_summary.tsv` and the snake-plot panel of
  `figures/phenotype_classes`. Provenance and attribution:
  [`data/annotation/SOURCE.md`](data/annotation/SOURCE.md).
