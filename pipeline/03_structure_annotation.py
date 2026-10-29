#!/usr/bin/env python3
"""Join the DMS positions to the active-state GPR161-Gs structure (PDB 8SMV).

    data/annotation/gpcrdb_generic_numbers.tsv
    data/structure/8smv.cif                      (fetched from RCSB if absent)
        -> data/annotation/structure_8smv.tsv

One row per assayed position (528) with the structural facts the figures need:

  segment, generic   GPCRdb's segment assignment and Ballesteros-Weinstein
                     generic number, copied through unchanged. ICL2 positions
                     take 34.xx numbers, NOT TM3 extensions -- P132 is 34.50
                     and M133 is 34.51, and an earlier draft of this analysis
                     had both wrong as 3.57/3.58.
  modeled            the residue has coordinates in 8SMV chain R. 268 of the
                     528 do: 28-220 and 266-340. The N-terminus (1-27), ICL3
                     (221-265) and the C-tail (341-529) are absent, so every
                     structural statement in the paper is about half the
                     protein and the figures say so.
  dist_gs            minimum heavy-atom distance from this residue to the
                     G protein heterotrimer (chains A, B, G), computed here.
"""
from pathlib import Path
import sys
import urllib.request

import gemmi
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
CIF = ROOT / "data/structure/8smv.cif"
OUT = ROOT / "data/annotation/structure_8smv.tsv"

RECEPTOR, PARTNER = "R", ("A", "B", "G")        # 8SMV: R = GPR161, A = Gas
# dSASA at the Gs interface used to be merged in here from a topos run in the
# analysis workspace. The manuscript does not use it, so the dependency -- the
# only path in this repository that pointed outside it -- was dropped 2026-10-28.


def fetch_cif():
    if CIF.exists():
        return
    CIF.parent.mkdir(parents=True, exist_ok=True)
    url = "https://files.rcsb.org/download/8SMV.cif"
    print(f"fetching {url}")
    with urllib.request.urlopen(url, timeout=120) as r:
        CIF.write_bytes(r.read())


def heavy_atoms(chain):
    """{residue number: (n_atoms, 3) array} for the amino acids of one chain."""
    out = {}
    for res in chain:
        info = gemmi.find_tabulated_residue(res.name)
        if info is None or not info.is_amino_acid():
            continue
        out[res.seqid.num] = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in res])
    return out


def main():
    fetch_cif()
    st = gemmi.read_structure(str(CIF))
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    model = st[0]

    rec = heavy_atoms(model[RECEPTOR])
    partner = np.vstack([a for ch in PARTNER for a in heavy_atoms(model[ch]).values()])
    tree = cKDTree(partner)

    gp = pd.read_csv(ROOT / "data/annotation/gpcrdb_generic_numbers.tsv", sep="\t")
    posn = pd.read_csv(ROOT / "data/processed/position_summary.tsv", sep="\t")

    t = posn[["position", "wt"]].merge(
        gp[["position", "gpcrdb_segment", "gpcrdb_generic", "gpcrdb_bw"]],
        on="position", how="left")
    t = t.rename(columns={"gpcrdb_segment": "segment",
                          "gpcrdb_generic": "generic",
                          "gpcrdb_bw": "bw"})
    t["modeled"] = t.position.isin(rec)
    t["dist_gs"] = [round(float(tree.query(rec[p])[0].min()), 3) if p in rec else np.nan
                    for p in t.position]

    # Integrity checks, not decoration: each of these was wrong at some point.
    assert len(t) == 528, len(t)
    assert int(t.modeled.sum()) == 268, int(t.modeled.sum())
    assert sorted(set(range(28, 221)) | set(range(266, 341))) == sorted(rec), "modeled range"
    assert t.loc[t.position == 132, "generic"].iloc[0] == "34.50x50"
    assert t.loc[t.position == 133, "generic"].iloc[0] == "34.51x51"

    OUT.parent.mkdir(parents=True, exist_ok=True)
    t.to_csv(OUT, sep="\t", index=False)
    print(f"wrote {OUT.relative_to(ROOT)}  {len(t)} rows")
    print(f"  modeled      {int(t.modeled.sum())}  (28-220, 266-340)")
    print(f"  within 4.5 A {int((t.dist_gs <= 4.5).sum())} positions")
    print(t[t.dist_gs <= 4.5][["position", "wt", "segment", "generic",
                               "dist_gs"]].to_string(index=False))


if __name__ == "__main__":
    sys.exit(main())
