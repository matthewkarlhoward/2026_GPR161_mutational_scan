#!/usr/bin/env python3
"""Collapse the variant table to one row per assayed position.

    data/processed/merged_variant_scores.tsv
    data/annotation/gpcrdb_generic_numbers.tsv
        -> data/processed/position_summary.tsv

Columns
  position, wt   residue number (UniProt Q8N6U8 numbering) and wild-type residue
  n              missense variants measured at this position (max 19)
  cre, abu, surf position mean of the missense z scores, one per assay
  resid          COUPLING RESIDUAL: the position mean of
                     cre_z - (a + b * surface_z)
                 where a, b are the ordinary least-squares fit of cAMP on
                 surface over all 9,927 missense variants (b = 0.631158,
                 a = -0.074695). Negative = less cAMP than this much surface
                 expression predicts.
  bw             Ballesteros-Weinstein generic number from GPCRdb (NA outside
                 the seven helices).

Caveat kept with the table, not buried: cAMP against surface is NOT linear
(quadratic term p = 8e-37 over the 9,927 missense variants, p = 5e-4 on the 528
position means), so `resid` is a linear-model summary of a curved
relationship. The figures built by this repository's scripts use `cre`, `abu`
and `surf`; `resid` is the residue attribute the hand-built ChimeraX panel
(manuscript Figure 4c) is coloured by.

Run:  python3 pipeline/02_position_summary.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
VARIANTS = ROOT / "data" / "processed" / "merged_variant_scores.tsv"
GENERIC = ROOT / "data" / "annotation" / "gpcrdb_generic_numbers.tsv"
OUT = ROOT / "data" / "processed" / "position_summary.tsv"


def main():
    d = pd.read_csv(VARIANTS, sep="\t")
    mis = d[d.type == "missense"].copy()

    # the cAMP ~ surface line, fitted at variant level over missense
    b, a = np.polyfit(mis.surface_z, mis.cre_z, 1)
    mis["rv"] = mis.cre_z - (a + b * mis.surface_z)
    print(f"cAMP ~ surface (missense, n = {len(mis):,}): "
          f"slope {b:.6f}, intercept {a:.6f}")

    g = mis.groupby("position")
    out = pd.DataFrame({
        "wt": g.wildtype.first(),
        "n": g.size(),
        "cre": g.cre_z.mean(),
        "abu": g.abundance_z.mean(),
        "surf": g.surface_z.mean(),
        "resid": g.rv.mean(),
    }).reset_index()

    gen = pd.read_csv(GENERIC, sep="\t")
    out["bw"] = out.position.map(dict(zip(gen.position, gen.our_bw)))

    out.to_csv(OUT, sep="\t", index=False, float_format="%.4f")
    print(f"\nwrote {OUT.relative_to(ROOT)}  {len(out)} positions x {out.shape[1]} columns")
    print(f"positions with a generic number: {out.bw.notna().sum()}")
    for col, label in (("cre", "cAMP"), ("abu", "abundance"), ("surf", "surface")):
        print(f"  {label:>9} position mean z < -2: {(out[col] < -2).sum():>3} positions")


if __name__ == "__main__":
    main()
