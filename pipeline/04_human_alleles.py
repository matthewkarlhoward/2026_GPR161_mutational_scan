#!/usr/bin/env python3
"""Derive data/annotation/human_alleles.tsv from the clinical-variant crosswalk.

`data/annotation/human_alleles.tsv` is an input to manuscript Figure 5
(`scripts/fig4_human_variants.py`). It ships here as a DERIVED ANNOTATION FILE,
committed in full, and this script records the derivation so that the table is
not an undocumented artefact.

Required input, and NOT part of this release:
  clinical_variant_dms_crosswalk.tsv   the join of gnomAD v4 and ClinVar records
      to this scan. It carries case-level clinical fields (case identifiers,
      disease category, zygosity) that are not needed to redraw the figure and
      are deliberately not redistributed. The script therefore does not run from
      a clone; it has no default location, and the crosswalk has to be named:

          python 04_human_alleles.py --crosswalk /path/to/crosswalk.tsv

The rule: take every `record_type == "human_allele"` row, keep five columns, and
collapse to one row per `protein_label`. Where an allele appears under several
records the sources are UNIONED, `gnomad_pass_any` is the OR, and the first
non-null ClinVar classification is kept. A plain de-duplication would silently
drop the ClinVar source for the ten alleles that appear twice, which is the trap
the figure's own source comments warn about.

One kind of token is dropped from `source`: a clinical case identifier
(`ClinGen_CCID:…`), which points at an individual patient record and does not
belong in a public deposit. The database origins `gnomAD_v4` and `ClinVar` are
kept, and so is a `literature_PMID…` token — that is the provenance of alleles
the paper names and discusses, so the deposit should record that it measured
them.

Expected output: 948 alleles — 654 missense, 284 synonymous, 10 in-frame
single-residue deletions.
"""
import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "annotation" / "human_alleles.tsv"
COLS = ["protein_label", "source", "consequence_class", "gnomad_pass_any",
        "clinvar_classification"]
DROP_SOURCE_PREFIXES = ("ClinGen_CCID:",)      # case-level patient identifiers


def union_sources(s):
    parts = set()
    for v in s.dropna():
        parts.update(str(v).split("|"))
    return "|".join(sorted(p for p in parts
                           if not p.startswith(DROP_SOURCE_PREFIXES)))


def first_non_null(s):
    v = s.dropna()
    return v.iloc[0] if len(v) else ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--crosswalk", type=Path, required=True,
                    help="the clinical-variant crosswalk this table is derived "
                         "from; not part of this release, so there is no "
                         "default (see the module docstring and "
                         "docs/DATA_SOURCES.md)")
    a = ap.parse_args()
    if not a.crosswalk.exists():
        raise SystemExit(
            f"crosswalk not found: {a.crosswalk}\n"
            "It is not part of this repository (see the module docstring and "
            "docs/DATA_SOURCES.md), so this script cannot run from a clone. "
            "The table it writes is committed at "
            "data/annotation/human_alleles.tsv.")

    cw = pd.read_csv(a.crosswalk, sep="\t", low_memory=False)
    h = cw[cw.record_type == "human_allele"].dropna(subset=["protein_label"]).copy()
    h["gnomad_pass_any"] = h["gnomad_pass_any"].fillna(False).astype(bool)

    g = (h.groupby("protein_label", as_index=False)
          .agg(source=("source", union_sources),
               consequence_class=("consequence_class", first_non_null),
               gnomad_pass_any=("gnomad_pass_any", "any"),
               clinvar_classification=("clinvar_classification", first_non_null)))
    g = g[COLS].sort_values("protein_label")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    g.to_csv(OUT, sep="\t", index=False)

    print(f"wrote {OUT.relative_to(ROOT)}  {len(g)} alleles "
          f"(from {len(h)} crosswalk rows)")
    print(g.consequence_class.value_counts().to_string())
    mis = g[g.consequence_class.isin(["missense", "inframe_indel"])]
    print(f"\nfigure sets: gnomAD PASS "
          f"{int((mis.source.str.contains('gnomAD') & mis.gnomad_pass_any).sum())}, "
          f"ClinVar {int(mis.source.str.contains('ClinVar').sum())}")


if __name__ == "__main__":
    main()
