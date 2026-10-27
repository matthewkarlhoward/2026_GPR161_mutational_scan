# Figure panels → code

Which script draws which panel, what it reads, and what was done by hand. Panels are
assembled and lettered in Illustrator; the scripts write each data panel as a
transparent-background PDF at its final millimetre size, so nothing is rescaled after
export.

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
`figures/phenotype_classes.pdf` (**165 × 104 mm**). It reads
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
`docs/METHODS.md` §5 and its code arrives with the figure that reports it.

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
