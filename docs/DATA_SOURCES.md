# Data sources, external annotations and citations

Everything in this repository that did not come out of the GPR161 experiment itself,
with where it came from and what to cite. Written so the entries can be lifted straight
into a manuscript methods section and reference list.

**Confidence flags.** ✅ = verified against the publisher record while assembling this
repository. ⚠️ = identity certain, one field worth checking against your reference
manager. ❓ = needs you (in-house method, or several papers could serve).

---

## 1. Software and methods used to produce this dataset

The library, counts and scores in `data/` were generated for the accompanying paper;
the tools below are cited as methods, not as sources of that data. The landing-pad cell
line and how the library was integrated are manuscript-methods detail and deliberately
out of scope here.

### Library design and cloning — DIMPLE ✅
The method the 11,043-variant oligo pool in `data/library/GPR161_oligos.csv` was designed
and assembled with.

> Macdonald, C. B., Nedrud, D., Grimes, P. R., Trinidad, D., Fraser, J. S. &
> Coyote-Maestas, W. DIMPLE: deep insertion, deletion, and missense mutation libraries
> for exploring protein variation in evolution, disease, and biology.
> *Genome Biology* **24**, 36 (2023). doi:10.1186/s13059-023-02880-6. PMID 36829241.

### Variant-effect scoring — Lilace ✅
The model fitted to this study's FACS bin counts to produce `data/scores/`.

> Freudenberg, J., Rao, J., Howard, M. K., Macdonald, C., Greenwald, N. F.,
> Coyote-Maestas, W. & Pimentel, H. Accurate variant effect estimation in FACS-based deep
> mutational scanning data with Lilace. *Genome Biology* **27**, 48 (2026).
> doi:10.1186/s13059-026-03934-1. PMID 41593735. Preprint: bioRxiv
> 10.1101/2025.06.24.661380.

### Read processing — Dumpling ✅
The pipeline that turned this study's FASTQ into the per-variant counts in `data/counts/`:
a Snakemake workflow (BBDuk → BBMerge → BBMap → GATK → tabulation against the designed
variants), released with the Rosace paper and public at
<https://github.com/odcambc/dumpling>. This run's configuration and construct reference
are committed at `pipeline/dumpling/`. Cite Rosace for the pipeline:

> Rao, J., Xin, R., Macdonald, C., Howard, M. K., Estevam, G. O., Yee, S. W., Wang, M.,
> Fraser, J. S., Coyote-Maestas, W. & Pimentel, H. Rosace: a robust deep mutational
> scanning analysis framework employing position and mean-variance shrinkage.
> *Genome Biology* **25**, 138 (2024). doi:10.1186/s13059-024-03279-7

Dumpling is cited here for read processing only: the scores in `data/scores/` come from
Lilace, not from the Enrich2 or Rosace scoring the pipeline can also run.

### Structure rendering — UCSF ChimeraX ⚠️

Manuscript Figure 4 panels a and c were rendered and labelled by hand in ChimeraX
1.13.dev (daily build, 2026-08-11), from `data/structure/8smv.cif` coloured by the
per-position values in `data/processed/position_summary.tsv`. That session is kept with
the manuscript figures, not here, and the render images are not deposited. **Nothing in
this repository needs ChimeraX** — everything `run_all.sh` builds, Figure 4 panel b
included, is pure Python.

> Meng, E. C. *et al.* UCSF ChimeraX: Tools for structure building and analysis.
> *Protein Science* **32**, e4792 (2023). doi:10.1002/pro.4792

### Scientific Python stack ⚠️
`numpy`, `pandas`, `scipy` and `matplotlib`, imported directly; `Pillow` is pinned
because matplotlib writes PNGs through it.
Versions are pinned in `requirements.txt`; cite as your journal requires.

---

## 2. External annotation shipped in `data/annotation/`

### `snakeplot_raw.svg` — GPCRdb snake-plot layout ✅
The residue layout of the GPR161 topology diagram in manuscript Figure 3c. Retrieved
**2026-09-08** by fetching <https://gpcrdb.org/protein/gp161_human/> and extracting the
`<svg id="snakeplot">` element verbatim from the returned HTML (the few lines that did
this are not part of this repository; the saved SVG is). Only the *layout* is GPCRdb's —
`scripts/style.py::snake_axes` reads the circle coordinates and ignores the colours. The raw SVG is committed so the figure does not depend on a live fetch.

### `gpcrdb_generic_numbers.tsv` — Ballesteros-Weinstein generic numbers ✅
Used for: the `bw` column of `data/processed/position_summary.tsv` (243 of 528
positions), and for naming helix positions (3.50, 6.48, 7.53 …) in the text. Retrieved
**2026-09-08** from the GPCRdb residue service for `gp161_human`. GPCRdb's sequence for
GPR161 matches canonical Q8N6U8, so its numbering maps onto these positions directly.
Columns: `gpcrdb_segment`, `gpcrdb_generic`, `gpcrdb_bw`, `our_bw` (the one used),
`agrees`.

Cite both the database and the numbering scheme:

> Taracena Herrera, L. P. *et al.* GPCRdb in 2025: adding odorant receptors, data mapper,
> structure similarity search and models of physiological ligand complexes.
> *Nucleic Acids Research* **53**, D425-D435 (2025). doi:10.1093/nar/gkae1065
> *(the current release paper; many reference managers still hold the earlier Kooistra
> et al. 2021, NAR 49, D335-D343 entry — either is defensible, but cite one consistently.)*

> Ballesteros, J. A. & Weinstein, H. Integrated methods for the construction of
> three-dimensional models and computational probing of structure-function relations in
> G protein-coupled receptors. *Methods in Neurosciences* **25**, 366-428 (1995).

GPCRdb's own citing page asks that its generic-numbering work be cited as well, so
include this alongside Ballesteros-Weinstein if the generic numbers appear in a figure:

> Isberg, V. *et al.* Generic GPCR residue numbers — aligning topology maps while minding
> the gaps. *Trends in Pharmacological Sciences* **36**, 22-31 (2015).
> doi:10.1016/j.tips.2014.11.001

---

### `human_alleles.tsv` — gnomAD and ClinVar alleles observed in GPR161 ⚠️

One row per distinct protein-level allele seen in human sequence data, used only by
`scripts/fig4_human_variants.py` (manuscript Figure 5).

| column | meaning |
|---|---|
| `protein_label` | the allele, `<wt><pos><mut>` or `<wt><pos>del`, 529-aa canonical numbering |
| `source` | pipe-joined origins: `gnomAD_v4`, `ClinVar`, `literature_PMID…` |
| `consequence_class` | `missense`, `synonymous`, `inframe_indel` |
| `gnomad_pass_any` | the allele passes gnomAD's own quality filter in at least one record |
| `clinvar_classification` | ClinVar's germline call, blank where there is none |

948 alleles: 654 missense, 284 synonymous, 10 in-frame single-residue deletions.
Applying the figure's filters gives 593 gnomAD-PASS missense and 112 ClinVar missense,
of which 587 and 111 have DMS scores — the n the figure reports.

**Provenance.** This file is a **derived annotation file**, committed in full because it
cannot be rebuilt from anything in this repository. It was extracted by
`pipeline/04_human_alleles.py` from a larger clinical-variant crosswalk, which joins
gnomAD v4 and ClinVar records to the scan. That crosswalk is deliberately **not**
redistributed here: it carries case-level clinical fields (case identifiers, disease
category, zygosity) that are not needed to redraw the figure. The script is committed so
that the derivation is on the record, but it cannot run from a clone.

Only the five columns above were carried over, and within `source` the clinical case
identifiers were dropped: they point at individual patient records. The
`literature_PMID…` tokens were kept, because they are the provenance of alleles the
paper names and discusses. Where one allele appeared under several records, the sources
are unioned, `gnomad_pass_any` is the OR, and the first non-null ClinVar classification
is kept — a plain de-duplication would drop the ClinVar source for the ten alleles that
appear twice.

**Dates.** gnomAD v4 and the ClinVar snapshot used for the crosswalk; see the manuscript
methods for the access dates. ❓ Confirm both before release.

## 3. Reference sequences and structures used for numbering

Every position in this repository is defined against these. The sequences are not
committed as files; the structure is — `pipeline/03_structure_annotation.py` downloads
8SMV from RCSB into `data/structure/8smv.cif` on its first run, and ChimeraX reads the
same file.

| resource | accession | used for | cite |
|---|---|---|---|
| UniProt ✅ | **Q8N6U8** (human GPR161, 529 aa, sequence v1; entry v185 when retrieved — the entry version moves, the sequence version does not) | the position numbering of everything | UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2025. *NAR* **53**, D609-D617 (2025). doi:10.1093/nar/gkae1010 |
| RefSeq / MANE ✅ | NP_001362812.1, NM_001375883.1 (the MANE Select pair; RefSeq labels it "transcript variant 12 / isoform 2" although it is byte-identical to UniProt canonical) | codon-level checks and cross-referencing clinical variants | O'Leary, N. A. et al. *NAR* **44**, D733-D745 (2016), doi:10.1093/nar/gkv1189 · Morales, J. et al. *Nature* **604**, 310-315 (2022), doi:10.1038/s41586-022-04558-8 |
| PDB ✅ | **8SMV** (GPR161-Gs cryo-EM; chain R author numbering = UniProt; residues 28-220, 266-340 modelled) | structural interpretation | cite the entry *and* the paper below |

> Hoppe, N. *et al.* GPR161 structure uncovers the redundant role of sterol-regulated
> ciliary cAMP signaling in the Hedgehog pathway. *Nature Structural & Molecular Biology*
> **31**, 667-677 (2024). *(the source of 8SMV)*

### `structure_8smv.tsv` — the DMS joined to 8SMV ⚠️

Built by `pipeline/03_structure_annotation.py`. Every column is computed here from the
files above: `segment` and `generic` are copied from GPCRdb, `modeled` is whether chain
R has coordinates for that residue, and `dist_gs` is the minimum heavy-atom distance to
the G protein heterotrimer. 18 positions come within 4.5 A.

---

## 4. Resources cited by the paper but not needed to rebuild this repository

None of these is read by any script here. They are listed so the reference list can be
assembled in one pass.

| resource | used for | citation |
|---|---|---|
| gnomAD v4.1 ⚠️ | population variants (Figure 5) | Chen, S. et al. *Nature* **625**, 92-100 (2024) — note this paper describes **v4**; v4.1 is a re-release of the same callset with no paper of its own, so cite the paper and state the version separately |
| ClinVar ⚠️ | clinical variant records | Landrum, M. J. et al. *NAR* **48**, D835-D844 (2020) |
| MaveDB ✅ | where these scores are deposited — experiment set urn:mavedb:00001323, one score set per assay (**-a-1 abundance, -b-1 surface, -c-1 cAMP**), published 2026-10-30 under CC BY 4.0 | Esposito, D. et al. *Genome Biology* **20**, 223 (2019) |
| Prior GPCR DMS ✅ | cited as earlier GPCR mutational scans; the Tedman surface-expression values are the x axis of Supplemental Figure 1f against basal cAMP | β2AR: Jones, E. M. et al. *eLife* **9**, e54895 (2020) · MC4R: Howard, C. J. et al. *eLife* **13** (2025) · GPCRome survey: Tedman, A. et al. *Nat. Commun.* (2026), doi:10.1038/s41467-026-76564-7 |

GPR161 biology background used in the introduction (all verified against the proposal's
reference list): Mukhopadhyay et al. *Cell* **152**, 210-223 (2013); Hwang et al. *eLife*
**10**, e67121 (2021); Matteson et al. *PNAS* **105**, 2088-2093 (2008); Begemann et al.
*JCO* **38**, 43-50 (2020); Kim et al. *Hum. Mol. Genet.* **28**, 200-208 (2019).
