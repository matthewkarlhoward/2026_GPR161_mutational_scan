# Pipeline

```
FASTQ ──[1]──▶ data/counts/        ──[2]──▶ data/scores/      ──[3]──▶ data/processed/
      dumpling  78 samples ×        lilace  3 score tables     here    merged_variant_scores.tsv
                11,043 variants                                        position_summary.tsv
```

Every stage operates on this study's own sequencing data. Steps **1** and **2** were run
with third-party software that is not re-run by this repository: their outputs — this
study's counts and scores — are committed, and what was run is recorded below so the
manuscript methods can state it exactly. Step **3** is here and regenerates its outputs
byte for byte.

---

## 1. Reads → counts (third-party software)

**Tool:** [**Dumpling**](https://github.com/odcambc/dumpling), the Snakemake pipeline
released with [Rosace](https://doi.org/10.1186/s13059-024-03279-7) (Rao et al. 2024).
It trims adapters and contaminants with BBDuk, merges read pairs with BBMerge, maps to
the construct reference with BBMap, calls variants with GATK, and tabulates them against
the designed variant list.

The run is documented by the two files in `dumpling/`, copied from the run directory:

| file | what it fixes |
|---|---|
| `GPR161.yaml` | the whole configuration: ORF coordinates, k-mer size, base-quality cut, adapter and contaminant sets, and which of the pipeline's optional outputs were produced |
| `GPR161_LP_ref.fasta` | the construct reference reads were mapped to — the landing-pad plasmid, 2,422 nt |

Settings that matter for reading the counts:

- `orf: "128-1714"` — 1,587 nt, exactly the 529 codons of GPR161, so **variant numbering
  in the counts is GPR161 numbering** with no offset to apply. Upstream of the ORF the
  construct carries an N-terminal signal peptide (`MKTIIALSYIFCLVFA`), a FLAG tag
  (`DYKDDDDA`, what the surface assay stains) and a short `SID` linker; downstream, a
  `GSGSGS` linker, the mNeonGreen2 11th strand (`TELNFKEWQKAFTDMM`, which the abundance
  assay reads), a P2A self-cleaving peptide and a puromycin-resistance marker. See
  `docs/METHODS.md` §1 for the whole architecture.
- `min_q: 30` — minimum base quality for GATK
- `kmers: 15`, `mem: 16`, `sam: "1.3"` — BBMap settings
- `remove_zeros: true`, `regenerate_variants: true`, `noprocess: false` — the designed
  variant list was regenerated from the oligo file, and variants are filtered against it
- `adapters`, `contaminants` — BBDuk's bundled adapter set, plus PhiX and the common
  sequencing-artifact set

**In:** demultiplexed paired-end FASTQ from an Illumina NovaSeq run (run folder
`230321_A00111_0991_AHMWGLDSX5`), plus the designed variant list built from
`data/library/GPR161_oligos.csv`.
**Out:** `data/counts/<sample>.csv.gz`, 78 files.

`data/library/GPR161_experiment_layout.csv` in this repository **is** the pipeline's
experiment file (`config/GPR161.csv` in the run directory, byte-identical), so the
`condition` / `replicate` / `time` / `tile` columns are read by Dumpling exactly as
documented in `GPR161.yaml`.

> **To record before submission:** the Dumpling commit or release used. The run also
> produced rejected-read and zero-count tables (≈440 MB) that are not committed here.
> Dumpling can emit Enrich2 and Rosace scores as well; those are not this study's scores,
> which come from Lilace (step 2).

Raw FASTQ files will be deposited in SRA; the accession goes here once assigned.

## 2. Counts → scores (third-party software)

**Tool:** [Lilace](https://doi.org/10.1186/s13059-026-03934-1) — a Bayesian hierarchical
model over per-bin counts that returns, per variant, a posterior mean effect, its
standard error, and a local false sign rate. It was applied to this study's counts; the
scores it produced are this study's data.

Fitted **once per assay**, over 3 replicates × 4 FACS bins. The output directories were
named `baseline_GPR161_<assay>_rescaled`, i.e. the baseline model with the rescaled
output, and the committed tables are those `variant_scores.tsv` files renamed (the
source files carried 2025-01-29 modification times in the working directory this
repository was assembled from; copying does not preserve that, so it is recorded here
rather than in the file metadata):

| Lilace run directory | committed as |
|---|---|
| `baseline_GPR161_cre_rescaled/variant_scores.tsv` | `data/scores/cre_variant_scores.tsv` |
| `baseline_GPR161_surface_rescaled/variant_scores.tsv` | `data/scores/surface_variant_scores.tsv` |
| `baseline_GPR161_abundance_rescaled/variant_scores.tsv` | `data/scores/abundance_variant_scores.tsv` |

> **To record before submission:** the Lilace version and the exact call. The invocation
> script was not preserved with the outputs.

Lilace's own definitions were re-derived from the committed tables and all check out —
`lfsr ≈ Φ(-|effect|/effect_se)` to within 0.032, `discoveryXX` exactly
`sign(effect)·1[lfsr < XX/100]`, position-level means agreeing with the empirical means
at r ≥ 0.99995. See `docs/METHODS.md` §4.

## 3. Scores → analysis tables (here)

```bash
python3 pipeline/01_merge_scores.py       # -> data/processed/merged_variant_scores.tsv
python3 pipeline/02_position_summary.py   # -> data/processed/position_summary.tsv
```

`01` inner-joins the three assay tables (10,927 variants measured in all three), converts
each effect to its own assay's synonymous-SD units, and assigns the five variant-level
mechanism classes. `02` collapses to positions, fits the cAMP~surface line used for the
coupling residual, and attaches GPCRdb generic numbers.

Both scripts print the constants they derive — the three synonymous SDs, the regression
coefficients, the class counts — so a rerun that changes anything is visible immediately
rather than silently. A rerun of both reproduces the two committed tables **byte for
byte**.

The analysis choices these scripts encode — where the cut goes, which assay is consulted
when, what a position-level call means — are documented in `docs/METHODS.md` §5, and the
position-level gate rule lives in exactly one place, `scripts/style.py::POSITION_GATE`.
