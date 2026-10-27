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
  resid_z        resid / (SD_syn / sqrt(n)), where SD_syn is the SD of the same
                 residual over the 502 synonymous variants (0.8423) -- i.e. the
                 position mean residual against the synonymous null, corrected
                 for how many variants the position carries.
  bw             Ballesteros-Weinstein generic number from GPCRdb (NA outside
                 the seven helices).

Caveat kept with the table, not buried: cAMP against surface is NOT linear
(quadratic term p = 8e-37 over the 9,927 missense variants, p = 5e-4 on the 528
position means), so `resid` is a linear-model summary of a curved
relationship. The robustness checks that survive LOWESS, surface-matched
controls and permutation live with the analysis that reports them, which is not
part of this repository yet; the figures here use `cre`, `abu` and `surf`, not
`resid`.

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
    syn = d[d.type == "synonymous"].copy()

    # the cAMP ~ surface line, fitted at variant level over missense
    b, a = np.polyfit(mis.surface_z, mis.cre_z, 1)
    for frame in (mis, syn):
        frame["rv"] = frame.cre_z - (a + b * frame.surface_z)
    sd_syn = syn.rv.std(ddof=1)
    print(f"cAMP ~ surface (missense, n = {len(mis):,}): "
          f"slope {b:.6f}, intercept {a:.6f}")
    print(f"synonymous residual SD: {sd_syn:.4f}  (n = {len(syn)})")

    g = mis.groupby("position")
    out = pd.DataFrame({
        "wt": g.wildtype.first(),
        "n": g.size(),
        "cre": g.cre_z.mean(),
        "abu": g.abundance_z.mean(),
        "surf": g.surface_z.mean(),
        "resid": g.rv.mean(),
    }).reset_index()
    out["resid_z"] = out.resid / (sd_syn / np.sqrt(out.n))

    gen = pd.read_csv(GENERIC, sep="\t")
    out["bw"] = out.position.map(dict(zip(gen.position, gen.our_bw)))

    out.to_csv(OUT, sep="\t", index=False, float_format="%.4f")
    print(f"\nwrote {OUT.relative_to(ROOT)}  {len(out)} positions x {out.shape[1]} columns")
    print(f"positions with a generic number: {out.bw.notna().sum()}")
    for col, label in (("cre", "cAMP"), ("abu", "abundance"), ("surf", "surface")):
        print(f"  {label:>9} position mean z < -2: {(out[col] < -2).sum():>3} positions")


if __name__ == "__main__":
    main()
