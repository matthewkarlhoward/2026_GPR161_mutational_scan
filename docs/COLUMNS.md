# Column dictionary

Every committed table, column by column. Positions are UniProt Q8N6U8 numbering
throughout (see `METHODS.md` §6).

---

## `data/library/GPR161_oligos.csv`
No header, `name,sequence`, one row per designed oligo (11,043). The file begins with a
UTF-8 BOM, so read it with `encoding="utf-8-sig"` (as the code here does) or the first
name comes back with a stray prefix; `GPR161_experiment_layout.csv` has one too.
Substitution and synonymous oligos (10,535) are named
`GPR161_DMS-<sublibrary>_<Wt><pos><Mut>` in three-letter code, e.g.
`GPR161_DMS-1_Ser2Cys`. The 508 deletion oligos use a different scheme,
`GPR161_delete-<sublibrary>_3-<pos>`, e.g. `GPR161_delete-1_3-2` — no residue letters.

Every sequence is 230 nt: flanking Golden Gate / amplification arms in upper case around
a lower-case insert, with the mutagenized codon itself as a 3 nt upper-case block inside
that insert (the 508 deletions and the 20 substitutions at position 529 have no internal
upper-case block, since there is no codon to mark).

## `data/library/GPR161_experiment_layout.csv`
One row per sequencing sample (78).

| column | meaning |
|---|---|
| `sample`, `file` | sample name; `file` matches `data/counts/<file>.csv.gz` |
| `condition` | `abundance` / `surface` / `cre` / `untreated`, with the `_mini` suffix marking the MiniSeq run of the same biological sample |
| `replicate` | biological replicate, 1-3 |
| `time` | **FACS bin index, 0-3** (bins A-D) — the column name is inherited from the pipeline's growth-assay layout and does not mean time here |
| `tile` | constant 1; the library was sequenced as one amplicon set |

## `data/counts/*.csv.gz`
One file per sequencing sample; every file carries the complete designed set (11,043
rows), so unobserved variants are explicit zeros.

| column | meaning |
|---|---|
| `count` | reads assigned to this variant in this sample |
| `pos` | residue position |
| `mutation_type` | `M` missense, `S` synonymous, `D` single-codon deletion |
| `name` | variant key, e.g. `S2C` |
| `codon` | designed codon |
| `mutation` | substituted residue; `D_1` for a deletion (the score tables write `del0` for the same thing) |
| `length` | 1 (single-residue variants only) |
| `hgvs` | `p.(S2C)` |

## `data/scores/{cre,surface,abundance}_variant_scores.tsv`
Lilace output, one row per variant **scored in that assay**: 10,973 (cAMP), 10,968
(surface), 10,954 (abundance). A variant is dropped from an assay when it has too few
reads there, which is why the merged table (variants present in all three) is smaller.
See `METHODS.md` §4 for the verification of each definition.

| column | meaning |
|---|---|
| `protein`, `exp` | `GPR161`; assay key |
| `hgvs`, `wildtype`, `position`, `mutation`, `type` | variant identity; `type` ∈ {missense, synonymous, deletion} |
| `effect` | posterior mean effect, centred on the synonymous mean. **Negative = less** readout |
| `effect_se` | posterior SD of `effect` |
| `lfsr` | local false sign rate |
| `pos_effect`, `pos_effect_se` | hierarchical position-level mean over non-synonymous variants, and its posterior SD; NA on synonymous rows |
| `pos_sd`, `pos_sd_se` | position-level SD and its posterior SD |
| `discovery01`, `discovery05`, `discovery10` | `sign(effect)` if `lfsr` < 0.01 / 0.05 / 0.10, else 0 |

## `data/processed/merged_variant_scores.tsv`
Built by `pipeline/01_merge_scores.py`. **10,927 variants** measured in all three assays
(9,927 missense, 502 synonymous, 498 deletions).

| column | meaning |
|---|---|
| `position`, `wildtype`, `mutation`, `type` | variant identity |
| `cre_eff`, `cre_se`, `cre_lfsr`, `cre_pos_eff`, `cre_pos_sd` | the Lilace columns above, for the cAMP assay |
| `abundance_*`, `surface_*` | the same five columns for the other two assays |
| `cre_z`, `abundance_z`, `surface_z` | `effect / SD(synonymous)` of that assay — the cross-assay comparable unit |
| `class` | variant-level mechanism class, five values, cut at z < -2 (`METHODS.md` §5) |

## `data/processed/position_summary.tsv`
Built by `pipeline/02_position_summary.py`. One row per assayed position (528).

| column | meaning |
|---|---|
| `position`, `wt` | residue number and wild-type residue |
| `n` | missense variants measured there (16-19) |
| `cre`, `abu`, `surf` | position mean of missense z, per assay |
| `resid` | coupling residual: position mean of `cre_z - (-0.074695 + 0.631158·surface_z)` |
| `bw` | GPCRdb generic number where it has one (243 of 528 positions); mostly but not only helix Ballesteros-Weinstein numbers, and 27 TM2 entries are the literal string `TM2` — see `METHODS.md` §5 |

## `data/annotation/gpcrdb_generic_numbers.tsv`
| column | meaning |
|---|---|
| `position`, `wt` | residue number and residue |
| `gpcrdb_segment` | N-term, TM1-TM7, ICL1-3, ECL1-3, H8, C-term |
| `gpcrdb_generic` | GPCRdb's structure-based generic number |
| `gpcrdb_bw`, `our_bw` | Ballesteros-Weinstein number as served, and the reconciled value used in `position_summary.tsv` |
| `agrees` | True only where both numbers are present and equal (184 of 529 rows). The 345 False rows are three different things: 43 where both values exist and differ, 286 where `our_bw` is blank, and 16 where `gpcrdb_bw` is blank — so the False count is not a disagreement count |

## `tables/` (written by the QC scripts)
| file | contents |
|---|---|
| `qc_position_coverage.csv` | pooled input reads per position |
| `qc_designed_observed.csv` | designed vs observed variants by class, at ≥1 and ≥10 reads |
| `qc_dropout_variants.csv` | the 111 designed variants under 10 pooled reads (94 at zero) |
| `qc_score_se_summary.csv` | n, median score, median SE and lfsr<0.05 fraction per assay × variant class |
| `qc_sample_coverage.csv` | reads and mean coverage per biological sample |
