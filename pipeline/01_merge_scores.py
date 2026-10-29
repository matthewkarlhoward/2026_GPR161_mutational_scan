#!/usr/bin/env python3
"""Merge the three per-assay Lilace tables into one variant-level table.

    data/scores/{cre,abundance,surface}_variant_scores.tsv   (Lilace output)
        -> data/processed/merged_variant_scores.tsv

One row per variant measured in ALL THREE assays (inner join on
position + wildtype + mutation + type), with three derived quantities:

  <assay>_z   effect / SD(synonymous effects of that assay).
              Lilace centres each assay on its own synonymous mean (the
              observed means are -2.2e-4 / 8.8e-5 / 1.3e-4), so z is the
              effect in units of the assay's own null spread and is directly
              comparable across assays. The three SDs are printed on each run
              and are quoted in the figure captions:

                  cAMP      0.185761
                  abundance 0.169196
                  surface   0.221162

  class       the five mechanism classes, cut at z < -2. That is a 4.6 %
              two-sided false-positive rate under a normal null; the synonymous
              variants themselves land outside |z| = 2 at 5.8 % (cAMP), 5.8 %
              (abundance) and 4.4 % (surface), so the real rate is close to,
              but not exactly, the nominal one. cAMP and
              surface are cut first; abundance is consulted only to split a
              variant that is missing from the surface into one that was never
              made (destabilized) and one that was made but not delivered
              to the surface (trafficking defect):

                  surface & cAMP lost   -> destabilized        (abundance also lost)
                                        -> trafficking defect  (abundance intact)
                  surface lost only     -> surface-loss, signaling-intact
                  cAMP lost only        -> signaling-selective (uncoupled)
                  neither               -> WT-like

              The position-level version of this rule lives in
              scripts/style.py (position_gates / position_class) and differs:
              there the most upstream affected gate names the class.

Run:  python3 pipeline/01_merge_scores.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SCORES = ROOT / "data" / "scores"
OUT = ROOT / "data" / "processed" / "merged_variant_scores.tsv"

ASSAYS = ["cre", "abundance", "surface"]
KEEP = {"effect": "eff", "effect_se": "se", "lfsr": "lfsr",
        "pos_effect": "pos_eff", "pos_sd": "pos_sd"}
KEY = ["position", "wildtype", "mutation", "type"]

CLASSES = {
    "destabilized": "destabilized (low abundance -> low surface)",
    "not_delivered": "trafficking-defective (made, not delivered)",
    "surface_only": "surface-loss, signaling-intact",
    "uncoupled": "signaling-selective (expressed, uncoupled)",
    "wt": "WT-like",
}
CUT = -2.0


def load(assay):
    d = pd.read_csv(SCORES / f"{assay}_variant_scores.tsv", sep="\t")
    d = d[KEY + list(KEEP)].rename(columns={k: f"{assay}_{v}" for k, v in KEEP.items()})
    return d


def main():
    out = load(ASSAYS[0])
    for a in ASSAYS[1:]:
        out = out.merge(load(a), on=KEY, how="inner")

    syn = out[out.type == "synonymous"]
    for a in ASSAYS:
        sd = syn[f"{a}_eff"].std(ddof=1)
        out[f"{a}_z"] = out[f"{a}_eff"] / sd
        print(f"{a:>9}  synonymous SD {sd:.6f}  "
              f"(mean {syn[f'{a}_eff'].mean():+.2e}, n = {len(syn)})")

    s, c, a = (out.surface_z < CUT), (out.cre_z < CUT), (out.abundance_z < CUT)
    out["class"] = np.where(
        s & c, np.where(a, CLASSES["destabilized"], CLASSES["not_delivered"]),
        np.where(s, CLASSES["surface_only"],
                 np.where(c, CLASSES["uncoupled"], CLASSES["wt"])))

    out = out.sort_values(["position", "type", "mutation"], kind="stable")
    out.to_csv(OUT, sep="\t", index=False)
    print(f"\nwrote {OUT.relative_to(ROOT)}  {len(out):,} variants x {out.shape[1]} columns")
    print(out.type.value_counts().to_string())
    print(out["class"].value_counts().to_string())


if __name__ == "__main__":
    main()
