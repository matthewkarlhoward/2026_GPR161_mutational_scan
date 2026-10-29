#!/usr/bin/env python3
"""Derive data/annotation/human_alleles.tsv from the clinical-variant crosswalk.

`human_alleles.tsv` is an input to manuscript Figure 5
(`scripts/fig4_human_variants.py`). It is committed to this repository, because
the crosswalk it comes from is NOT redistributed here: that file carries
case-level clinical fields (case identifiers, disease category, zygosity,
ClinGen points) that are not needed to redraw the figure. This script records
the derivation so the committed table is not an undocumented artefact.

It therefore does not run from a fresh clone. Point it at the crosswalk:

    python 04_human_alleles.py --crosswalk /path/to/clinical_variant_dms_crosswalk.tsv

The rule: take every `record_type == "human_allele"` row, keep five columns, and
collapse to one row per `protein_label`. Where an allele appears under several
records the sources are UNIONED, `gnomad_pass_any` is the OR, and the first
non-null ClinVar classification is kept. A plain de-duplication would silently
drop the ClinVar source for the ten alleles that appear twice, which is the trap
the figure's own source comments warn about.

Expected output: 948 alleles — 654 missense, 284 synonymous, 10 in-frame
single-residue deletions.
"""
import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "annotation" / "human_alleles.tsv"
DEFAULT_CROSSWALK = (ROOT.parent / "gpr161_analysis" / "tables"
                     / "clinical_variant_dms_crosswalk.tsv")
COLS = ["protein_label", "source", "consequence_class", "gnomad_pass_any",
        "clinvar_classification"]


def union_sources(s):
    parts = set()
    for v in s.dropna():
        parts.update(str(v).split("|"))
    return "|".join(sorted(parts))


def first_non_null(s):
    v = s.dropna()
    return v.iloc[0] if len(v) else ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--crosswalk", type=Path, default=DEFAULT_CROSSWALK)
    a = ap.parse_args()
    if not a.crosswalk.exists():
        raise SystemExit(
            f"crosswalk not found: {a.crosswalk}\n"
            "It is not part of this repository (see the module docstring and "
            "docs/DATA_SOURCES.md). Pass --crosswalk to point at it.")

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
