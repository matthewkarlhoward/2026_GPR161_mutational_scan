# Methods — as implemented in this repository

What the code does and what the committed data are, written to be checkable against the
files. Anything not establishable from those files is marked **[CONFIRM]** rather than
guessed.

---

## 1. Library

| | |
|---|---|
| Target | human GPR161, UniProt Q8N6U8, 529 aa |
| Designed variants | **11,043** (`data/library/GPR161_oligos.csv`) |
| Missense | 10,032 |
| Synonymous | 503 |
| Single-residue deletions | 508 |
| Positions covered | 2-529 (528 positions, median 19 substitutions per position) |
| Oligo pool | 11 sub-libraries, 230 nt each. Two naming schemes: `GPR161_DMS-<n>_Ser2Cys` for the 10,535 substitutions and synonymous, `GPR161_delete-<n>_3-<pos>` for the 508 deletions |

Each oligo name carries its own variant, which is how the sub-library boundaries in the
QC coverage panel are derived (`scripts/qc/dms_qc_coverage.py::tile_bounds`). That
derivation reads the `GPR161_DMS-` names only; the deletion oligos use the other scheme
and do not contribute boundaries, which does not change them because the two schemes
tile the same positions.

The library was built for this study following the **DIMPLE** method (Macdonald et al.
2023 — the method paper, which did not make this library); the oligo pool is
assembled by Golden Gate into the receptor backbone and integrated at a single genomic
locus so that each cell carries one variant.

The construct is in `pipeline/dumpling/GPR161_LP_ref.fasta`, the reference reads were
mapped to. Read end to end, it is a single open reading frame:

| element | sequence | what it is for |
|---|---|---|
| signal peptide | `MKTIIALSYIFCLVFA` | surface targeting |
| FLAG tag | `DYKDDDDA` | the epitope the **surface** assay stains |
| linker | `SID` | |
| **GPR161 1-529** | nucleotides 128-1714 | the variable region; the 1,587 nt are exactly 529 codons, so assay numbering is GPR161 numbering |
| linker | `GSGSGS` | |
| mNeonGreen2 11th strand | `TELNFKEWQKAFTDMM` | completes the split reporter the **abundance** assay reads |
| P2A | `ATNFSLLKQAGDVEENPGP` | self-cleaving peptide |
| PuroR | puromycin N-acetyltransferase, to the end of the reference | selection |

So two of the three readouts are built into this construct, and the third (cAMP) comes
from the separate CRE-eGFP reporter. The library was integrated into a landing-pad
HEK293T line so that each cell carries a single variant; the line, the integration and
the complementing mNeonGreen2(1-10) fragment are manuscript-methods detail and are not
described here.

## 2. Assays

Three FACS-based assays on the same pooled library in HEK293T:

| key | assay | readout |
|---|---|---|
| `surface` | cell-surface expression | anti-FLAG antibody stain against the extracellular tag |
| `abundance` | whole-cell abundance | split-mNeonGreen complementation |
| `cre` | basal cAMP signalling | CRE transcriptional reporter driving eGFP |

GPR161 is constitutively Gs-coupled, so the cAMP assay reads **basal** signalling with no
agonist added; loss of signalling is a drop in reporter fluorescence.

Sample structure, from `data/library/GPR161_experiment_layout.csv` (78 sequencing
samples):

- 3 biological replicates × 4 FACS bins (the `time` column, 0-3 = bins A-D) per assay
  = 12 sorted samples per assay
- 3 unsorted library-input samples (`untreated`)
- every sample sequenced on two runs, a MiniSeq run (`*_mini`) and a NextSeq run; the two
  runs are **summed per biological sample** before coverage is reported
  (`scripts/qc/dms_qc_scores.py::load_sample_coverage`)

Median sequencing coverage is 337x (surface), 342x (abundance) and 390x (cAMP) reads
per designed variant; individual sorted samples span 22x to 528x
(`tables/qc_sample_coverage.csv`).

**[CONFIRM]** cytometer, sort gates and cell numbers; the gating plots in the QC
supplement were exported from the FACS software and are not reproduced by this code.

## 3. Reads to counts

`data/counts/` holds one table per sequencing sample, each containing the **complete
designed variant set** (11,043 rows) with a `count` column, so variants never observed in
a sample are present as explicit zeros.

| column | meaning |
|---|---|
| `count` | reads assigned to this variant in this sample |
| `pos` | residue position (UniProt numbering) |
| `mutation_type` | `M` missense, `S` synonymous, `D` single-codon deletion |
| `name` | variant key, e.g. `S2C` |
| `codon` | the designed codon |
| `mutation` | substituted residue; `D_1` for a deletion (the score tables write the same thing as `del0`) |
| `length` | 1 |
| `hgvs` | `p.(S2C)` |

These were produced from paired-end Illumina reads (NovaSeq run
`230321_A00111_0991_AHMWGLDSX5`) with **Dumpling** (<https://github.com/odcambc/dumpling>),
the Snakemake pipeline released with the Rosace paper (Rao et al. 2024): BBDuk for
adapter and contaminant removal, BBMerge to merge read pairs, BBMap to the construct
reference, GATK for variant calling, then tabulation against the designed variant list.

The run's configuration is committed at `pipeline/dumpling/GPR161.yaml` and the construct
reference at `pipeline/dumpling/GPR161_LP_ref.fasta`. The settings that affect how the
counts should be read:

| setting | value | consequence |
|---|---|---|
| `orf` | `128-1714` | 1,587 nt = exactly GPR161's 529 codons, so **counts are in GPR161 numbering** with no offset |
| `min_q` | 30 | minimum base quality in GATK |
| `kmers` / `mem` / `sam` | 15 / 16 GB / 1.3 | BBMap settings |
| `regenerate_variants` | true | the designed variant list was rebuilt from the oligo file, which is why every counts file carries all 11,043 rows |
| `remove_zeros` | true | a zero-removed copy was written alongside; the committed counts are the full tables, zeros included |
| adapters / contaminants | BBDuk defaults + PhiX + sequencing artifacts | |

`data/library/GPR161_experiment_layout.csv` **is** the pipeline's experiment file
(byte-identical to `config/GPR161.csv` in the run directory), so `condition`,
`replicate`, `time` and `tile` mean what `GPR161.yaml` documents them to mean.

**[CONFIRM]** the Dumpling commit or release. Variants rejected by the pipeline's filters
and the zero-removed tables were kept outside this repository (≈440 MB). Dumpling also
offers Enrich2 and Rosace scoring; neither produced the scores used here, which come from
Lilace (§4).

Library QC, from `scripts/qc/dms_qc_coverage.py` over the three unsorted input samples:

| | designed | observed (≥1 read) | % | observed (≥10 reads) | % |
|---|---|---|---|---|---|
| Total | 11,043 | 10,949 | 99.1 | 10,932 | 99.0 |
| Synonymous | 503 | 502 | 99.8 | 502 | 99.8 |
| Missense | 10,032 | 9,941 | 99.1 | 9,926 | 98.9 |
| Deletion | 508 | 506 | 99.6 | 504 | 99.2 |

Per-position pooled coverage: median 21,442.5 reads, 9.5× range across positions.

## 4. Counts to scores — Lilace

The scores in `data/scores/` are this study's measurements of GPR161. They are posterior
estimates produced by fitting **Lilace** (Freudenberg et al. 2026, *Genome Biology*
27:48), a Bayesian hierarchical model, to this study's per-bin counts — the paper is
cited for the estimator, not as a source of the data. Each assay was fitted separately;
the committed tables are that fit's rescaled output:

```
data/scores/cre_variant_scores.tsv
data/scores/surface_variant_scores.tsv
data/scores/abundance_variant_scores.tsv
```

**[CONFIRM]** the Lilace version and the exact call. The run directories were named
`baseline_GPR161_<assay>_rescaled`, i.e. the baseline model with the rescaled output
option, but the invocation script was not preserved with the outputs.

Column meanings, verified numerically against the committed tables:

| column | definition | check |
|---|---|---|
| `effect` | posterior mean variant effect on Lilace's latent scale, centred on the synonymous mean | synonymous means are -2.2e-4 (cAMP), +1.3e-4 (surface), +8.8e-5 (abundance) |
| `effect_se` | posterior SD of the effect | `lfsr ≈ Φ(-|effect|/effect_se)`, max abs deviation 0.032 |
| `lfsr` | local false sign rate, min{P(β≥0), P(β≤0)} | estimated on a 0.00025 grid, which is why 5 of 32,895 rows sit marginally above the 0.5 ceiling the definition implies (max 0.50075) |
| `pos_effect`, `pos_sd` | hierarchical position-level mean and SD over non-synonymous variants; NA for synonymous rows | agree with the empirical position mean at r ≥ 0.99995 (0.99998 cAMP, 0.99996 abundance, 0.99996 surface) |
| `discovery01/05/10` | `sign(effect) · 1[lfsr < 0.01/0.05/0.10]` | verified exactly on all rows |

Sign convention: **negative = less** of the readout than synonymous (less basal cAMP,
less surface display, less protein).

Magnitudes are **not** comparable across assays on the raw scale — each assay has its own
rescaling constant and its own null width. Convert to z first.

## 5. Scores to analysis tables

### `pipeline/01_merge_scores.py` -> `data/processed/merged_variant_scores.tsv`

Inner join of the three assay tables on (position, wildtype, mutation, type):
**10,927 variants** measured in all three assays — 9,927 missense, 502 synonymous, 498
deletions.

Adds, per assay, `<assay>_z = effect / SD(synonymous effects of that assay)`:

```
cAMP      0.185761
abundance 0.169196
surface   0.221162
```

> These are quoted to six figures because they are **exact constants of this dataset**
> and the pipeline must reproduce them, not because the SD of a null is known that
> precisely. Each rests on 502 synonymous variants, so as an estimate it carries about
> 3 % relative standard error (0.186 ± 0.006). The same applies to the regression
> coefficients and the residual SD below. Quote them to 3 figures in prose.

Adds `class`, the **variant-level** mechanism class, cut at z < -2. cAMP and surface are
cut first, and abundance is consulted only to split a variant missing from the surface
into one that was never made and one that was made but not delivered:

| surface | cAMP | abundance | class | n |
|---|---|---|---|---|
| lost | lost | lost | `destabilized (low abundance -> low surface)` | 889 |
| lost | lost | intact | `trafficking-defective (made, not delivered)` | 654 |
| lost | intact | — | `surface-loss, signaling-intact` | 452 |
| intact | lost | — | `signaling-selective (expressed, uncoupled)` | 337 |
| intact | intact | — | `WT-like` | 8,595 |

(Counts over all 10,927 variants; restricted to the 9,927 missense they are 799 / 585 /
430 / 295 / 7,818, which are the numbers plotted in Figure 2a,b.)

### `pipeline/02_position_summary.py` -> `data/processed/position_summary.tsv`

One row per assayed position (528 rows):

- `n` — missense variants measured there (max 19)
- `cre`, `abu`, `surf` — position mean of the missense z scores
- `resid` — the **coupling residual**: position mean of `cre_z - (a + b·surface_z)`,
  where `b = 0.631158`, `a = -0.074695` is the OLS fit of cAMP on surface over all 9,927
  missense variants
- `resid_z` — `resid / (0.8423 / sqrt(n))`, where 0.8423 is the SD of the same residual
  over the 502 synonymous variants
- `bw` — GPCRdb generic number, present at 243 of 528 positions. Mostly
  Ballesteros-Weinstein numbers for helix residues, but not exactly: 40 of the 243 sit
  outside TM1-TM7 (H8, C-term, ICL2, ECL2), 30 positions inside the helices have none,
  and 27 TM2 entries carry the literal string `TM2` rather than a number. Treat the
  column as "GPCRdb's label for this position, where it has one" 

> cAMP against surface is **not** linear (quadratic term p = 8e-37 over the 9,927
> missense variants; p = 5e-4 on the 528 position means), so `resid` is a
> linear summary of a curved relationship. It is reported for completeness; the figures
> in this repository use `cre`, `abu` and `surf`. The robustness checks — LOWESS
> conditioning, surface-matched controls, split-half noise handling, permutation — live
> with the analysis that reports them and are not part of this repository yet.

### Position-level classes are a different rule from variant-level classes

The **position** rule (`scripts/style.py::position_gates`, `position_class`) treats the
three assays as a chain and names the class for the **most upstream affected gate**: a
gate is affected when the position mean z < -2; abundance affected -> `destabilized`;
otherwise surface affected -> `not delivered` if cAMP is also lost, else `surface-loss
only`; otherwise cAMP affected -> `uncoupled`. Applied at variant level the two rules
disagree on 2.2 % of calls (243 of 10,927); applied at position level, on 0.4 % (2 of
528). Figure 2a,b are variant-level; Figure 2c,d are
position-level, and the Sankey is drawn in the chain order abundance -> surface -> cAMP.

The gate rule is a real analysis choice and is kept in **one place**,
`style.POSITION_GATE`. The alternative ("at least 3 variants with z < -2 and lfsr <
0.05") was rejected because the three assays detect at very different rates (5.6 %
abundance vs 15.5 % surface vs 14.1 % cAMP), which made the abundance-dependent split a
function of the threshold rather than of the biology.

### Does the chain hold?

The Sankey in Figure 2d is drawn in the order abundance → surface → cAMP, which is a
claim, not a convention. It was tested separately: conditional on surface, abundance
carries no further information about cAMP (partial rho = +0.03, p = 0.49; adding
abundance to a surface-only model of position-mean cAMP buys -0.02 out-of-fold R²), and
the dependence of cAMP on surface is graded rather than a threshold. That test has to
handle three traps — the nonlinearity above, shared measurement noise across assays
(synonymous cross-assay Pearson r = 0.55 cAMP-surface, 0.48 surface-abundance, 0.40
cAMP-abundance, handled by split-half estimation), and regression
dilution from conditioning on a noisy mediator, which biases the test *towards* finding a
direct path, so a null is the strong result. **The code for it is not in this repository
yet**; it arrives with the figure that reports it.

## 6. Numbering and mapping

- Positions are UniProt **Q8N6U8** (529 aa). The assay wild-type sequence was rebuilt
  from the `position`/`wildtype` columns of the score tables and is identical across the
  three assays; it aligns to Q8N6U8 with no indels at 528/528 assayed positions, and
  matches RefSeq NP_001362812.1.
- **PDB 8SMV** chain R author numbering is identical to UniProt numbering; 268 residues
  are modelled (28-220, 266-340).
- **Literature trap.** Part of the GPR161 literature uses the 549-aa isoform, which is
  **+20** relative to canonical numbering; the mouse literature adds its own offset
  (e.g. mouse L465P). Every literature variant must be mapped before comparison.
- GPCRdb's sequence for `gp161_human` matches the canonical sequence, so its generic
  numbers map directly onto these positions
  (`data/annotation/gpcrdb_generic_numbers.tsv`).

## 7. Figure conventions

Set in `scripts/style.py` and `scripts/qc/gpr161_style.py`:

- Helvetica, **6 pt**, black text everywhere — every text element, with no exceptions;
  panel letters are added in Illustrator, not by the code
- **0.5 pt** black lines for axes, ticks, plot lines and marker edges
- composite figures capped at **165 × 180 mm**; `save_figure()` enforces the cap, refuses
  to silently expand the canvas, and reports clipped text, text-on-text collisions and
  ink within 2 mm of the page edge
- PDFs are saved with `pdf.fonttype = 42` (text stays editable text) and a transparent
  background, because panels are assembled in Illustrator
- the figures were rendered on macOS, where Helvetica is present. The font stack falls
  back to Helvetica Neue, Arial, then DejaVu Sans, so a machine without Helvetica
  rebuilds the same figure with different glyph metrics — the PDFs stay valid and the
  numbers identical, but the PNGs will not be byte-identical to the committed ones
- mechanism-class colours are fixed across every figure: uncoupled `#8F2924`, not
  delivered `#4877B5`, destabilized `#65159E`, surface-loss only `#316B3E`, WT-like
  `#CDCECE`. These live in `scripts/qc/gpr161_style.py` and are repeated explicitly at
  the top of `scripts/fig2_phenotype_classes.py`, which is the figure that draws them —
  two copies, so check both if you change one

The heatmaps in Figure 1c use a separate diverging scale built from `BLUE_GOF #3d6ca8`
(gain) and `RED_LOF #9e2b26` (loss), symmetric at ±1.2 raw effect units, which clips
2.6 % of values — all of them deep loss. Note this is the **reverse** of the blue-is-loss
convention used in some earlier drafts.
