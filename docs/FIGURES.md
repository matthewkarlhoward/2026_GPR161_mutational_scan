# Figure panels → code

Which script draws which panel, what it reads, and what was done by hand. Panels are
assembled and lettered in Illustrator; the scripts write each data panel as a
transparent-background PDF at its final millimetre size, so nothing is rescaled after
export.

## Manuscript figure → script

Headings below for Figures 4 and 5 use the manuscript's numbers. The first two
sections still carry this repository's older internal numbering, which is one less
than the manuscript's. Use this table to go from the paper to the code.

| manuscript | this repo draws | script | output |
|---|---|---|---|
| **Figure 2** — GPR161 DMS overview | panels b, c | `scripts/fig1_landscape.py` | `figures/landscape_density_heat.pdf` |
| **Figure 3** — variant effects across phenotypes | panels b–e | `scripts/fig2_phenotype_classes.py` | `figures/phenotype_classes.pdf` |
| **Figure 4** — structural interpretation | panel b only | `scripts/fig3_structure.py` | `figures/segment_effects.pdf` |
| **Figure 5** — human variant interpretation | panels a, b | `scripts/fig4_human_variants.py` | `figures/human_variants.pdf` |
| **Supplemental Figure 3** — library and sequencing QC | panels a–d | `scripts/qc/*` | `figures/F_QC_panel.pdf` |

Figure 1 and Supplemental Figures 1–2 are the GPCRome screen and live in the companion
repository, `2026_GPCRome_basal_cAMP`. The cartoon panels (Figure 2a, Figure 3a,
Supplemental Figure 3e) are drawn by hand in Illustrator and have no code here.
Figure 4 panels a and c are ChimeraX renders made by hand from
`data/structure/8smv.cif`, coloured by the per-position values in
`data/processed/position_summary.tsv`; the renders themselves are not deposited.

---

## Figure 1 — the platform and the landscape

| panel | content | built by |
|---|---|---|
| **a** | assay schematic (library → stable lines → FACS → NGS → variant effect map) | **drawn by hand** (Illustrator), no code |
| **b** | variant score distributions, one row per assay: synonymous, missense, single-residue deletions | `scripts/fig1_landscape.py` → `figures/landscape_density_heat.pdf` (left column) |
| **c** | phenotype maps: 528 positions × 20 substitutions + a deletion row, one heatmap per assay | same script (right column) |

Reads `data/processed/merged_variant_scores.tsv`. Output is exactly **112 × 43.9 mm**.

Details that matter for the caption:
- Both columns are in **raw DMS effect units**, not z. The raw scale is already centred
  on the synonymous mean (+1e-4 / +1e-4 / -2e-4), so white = 0 = mean synonymous.
- Heatmap rows are ordered charged → polar → special → hydrophobic
  (`DEKRHNQSTGPAVILMFYWC`); the deletion row is set off by a blank strip. The
  single-letter row labels and the panel letters are added in Illustrator.
- Colour scale ±1.2, **blue = gain, red = loss**, saturating on 2.6 % of values, all of
  them deep loss.
- Missing variants appear as white gaps (94 designed variants were never observed).

## Figure 2 — three phenotypes, five mechanism classes

All four panels come from one script, `scripts/fig2_phenotype_classes.py` →
`figures/phenotype_classes.pdf` (**165 × 92 mm**). It reads
`data/processed/merged_variant_scores.tsv` and `data/processed/position_summary.tsv`, and
the snake-plot layout from `data/annotation/snakeplot_raw.svg`.

| panel | content | note |
|---|---|---|
| **a** | abundance vs surface, all 9,927 missense variants, coloured by class | ρ = +0.63 |
| **b** | cAMP vs surface, same variants | ρ = +0.66 |
| **c** | GPCRdb snake plot coloured by **position** class | layout from GPCRdb, colours from this study |
| **d** | all 528 positions flowing through the three gates | ribbon width = positions |

- Surface is on the x axis in **both** scatters because it is the shared reference the
  classes are cut against. Panel a plots abundance (upstream of surface) on y: it is a
  threshold map, not a regression, so no direction is implied.
- The dashed lines are **two synonymous SDs of each axis's own assay**, so they sit at a
  different raw value on every axis: cAMP -0.372, abundance -0.338, surface -0.442.
- Panels a,b use the **variant-level** class rule; c,d use the **position-level** rule.
  They share names and colours but are not the same rule — see `docs/METHODS.md` §5.
- The script prints the full class and path tables, plus ready-made caption text, when it
  runs.

The chain order of the Sankey (abundance → surface → cAMP) is a claim that was tested
separately rather than asserted by the figure; the result is summarised in
`docs/METHODS.md` §5 and its code is `scripts/analysis/dag_model_test.py`.

---

## Figure 4 (manuscript) — the scan on the structure

**Only panel b is drawn here.** `scripts/fig3_structure.py` →
`figures/segment_effects.pdf` (**50 × 40 mm**). It reads
`data/annotation/structure_8smv.tsv`, built by `pipeline/03_structure_annotation.py`.
268 of the 528 assayed positions have coordinates in 8SMV — 28-220 and 266-340 — so
everything here is about half the protein, and the N-terminus, ICL3 and the C-tail are
absent rather than unaffected.

| panel | content | built by |
|---|---|---|
| **a** | the receptor coloured by mean missense **surface** effect | **ChimeraX, by hand** |
| **b** | the same quantity averaged by GPCRdb segment, N to C | `scripts/fig3_structure.py` |
| **c** | the cytoplasmic face coloured by the **coupling residual**, uncoupled positions as sticks | **ChimeraX, by hand** |

Panels a and c were rendered and labelled by hand from `data/structure/8smv.cif`,
coloured by the per-position surface and cAMP values in
`data/processed/position_summary.tsv` (`surf`, `cre`, `resid`), applied as residue
attributes in a ChimeraX session kept with the manuscript figures. Neither the session
nor the render images is deposited. An earlier draft assembled all three panels in matplotlib from
a `figures/renders/` folder — that code and that folder were removed on 2026-10-28.

Details that matter for the caption:

- **Panel b is the signed mean**, in raw DMS effect units: the plain mean of the missense
  surface effects measured at each position, averaged by segment. Losses are negative and
  a gain at one position cancels a loss at another, so ICL1, ICL2 and H8 come out
  slightly positive. It is **surface only** — cAMP and abundance are computed and printed
  by the script for the text, but not plotted.
- Bars are ordered **N to C along the sequence**, not by effect size, and carry no error
  bars. C-term is the one segment dropped: a single modelled position is not a bar.
  ECL1 has four positions and is kept.
- Panel a saturates at 1.2, matching the Figure 2 heatmaps; 9 of the 268 modelled
  positions (3.4 %) are past it. Panel c's ±4 scale never saturates: the residual runs
  −3.81 to +1.30 over the positions it colours.
- **Panel c masks 94 of the 268 positions in grey.** They are the positions whose own
  surface gate is affected (mean surface z < −2). Too little receptor reaches the
  membrane there for coupling to be tested at all, so the residual is not interpretable
  and the figure declines to colour it. All ten uncoupled positions are outside the mask
  by construction — that class requires the surface gate to be spared.
- Panel c shows **Gαs only**. Gβ, Gγ and the Nb35 nanobody are deleted. The interface
  number this repository can still state is **18 positions within 4.5 Å** (`dist_gs`);
  the buried-area footprint was removed with the dSASA column on 2026-10-28, see
  `DATA_SOURCES.md`.

---

## Figure 5 (manuscript) — human variant interpretation

One script, `scripts/fig4_human_variants.py` → `figures/human_variants.pdf`. Reads
`data/processed/merged_variant_scores.tsv` and `data/annotation/human_alleles.tsv`.
Panels are drawn 120 : 40 mm wide; the saved page is narrower than 165 mm because this
figure is written with `bbox_inches="tight"`, so the ratio is what the script sets, not
the absolute width.

| panel | content | note |
|---|---|---|
| **a** | three columns (abundance, surface, cAMP) sharing the z axis; the full scan, the gnomAD set and the ClinVar set as density-spread swarms | 10,425 / 594 / 111 |
| **b** | % loss of function in each set, stacked by mechanism class | 21.2 / 12.8 / 9.0 % |

Two things a caption has to state:

- **The two panels use different universes.** Panel a is missense **plus single-residue
  deletions** (10,425); panel b is missense only (9,927), because the mechanism classes
  are defined on missense. The per-set n differ accordingly: gnomAD is 594 in a and 587
  in b, ClinVar is 111 in both (its one in-frame deletion, V228del, has no deletion
  measurement — the library scored 498 of 528 positions and 228 is a gap).
- **The swarms are not a true beeswarm.** A non-overlapping pack is not achievable for
  10,425 points in a ~24 mm column at a readable marker size, so each point's horizontal
  offset is scaled by the kernel density at its own y. Each swarm is scaled to its *own*
  peak density, so width shows shape, not sample size — 111 ClinVar points spread as
  wide as 10,425. The n in the legend is the only correction for that.

Solid line at 0 is no effect; the dashed line is the loss-of-function cut at two
synonymous SDs (z = −2). No panel letters are drawn; they are added in Illustrator.

## Supplementary figure — library and sequencing QC

`scripts/qc/dms_qc_panel.py` draws panels a-d as a single **165 × 30 mm** row
(`figures/F_QC_panel.pdf`). It reuses the loaders in `dms_qc_coverage.py` and
`dms_qc_scores.py`, which are therefore kept: run on their own they write the QC tables
in `tables/`, and with `--figures` they also draw larger standalone versions of the same
panels (not committed — the row above is the figure).

| panel | content | built by | table |
|---|---|---|---|
| **a** | total read count per residue position in the library input, log scale, dashed lines at the 11 oligo sub-library boundaries | `dms_qc_coverage.py` | `tables/qc_position_coverage.csv` |
| **b** | designed vs observed variants by class | `dms_qc_coverage.py` | `tables/qc_designed_observed.csv`, `qc_dropout_variants.csv` |
| **c** | Lilace standard error, one violin per assay | `dms_qc_scores.py` | `tables/qc_score_se_summary.csv` |
| **d** | mean sequencing coverage per sample, one row per assay | `dms_qc_scores.py` | `tables/qc_sample_coverage.csv` |
| **e** | FACS gating strategy, three assays | **exported from the cytometry software**, no code |

Panel d counts samples after summing the MiniSeq and NextSeq runs of the same biological
sample: 12 sorted samples per assay (3 replicates × 4 bins) and 3 library-input samples.

---

## House style

One place: `scripts/style.py` (main figures) and `scripts/qc/gpr161_style.py` (QC
figures). Helvetica 6 pt black, 0.5 pt black rules, 165 × 180 mm maximum, vector PDF with
live text, transparent background. `style.save_figure()` refuses to expand the canvas past
the cap and reports clipped text, overlapping text and ink within 2 mm of the edge, so a
panel that silently grew in an earlier draft fails loudly instead.
