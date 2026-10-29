#!/usr/bin/env python3
"""Which causal model fits: serial chain, or a direct abundance -> cAMP path?

M1 (serial)  abundance -> surface -> cAMP.  Testable prediction: conditional on
             surface, abundance carries NO further information about cAMP.
M3 (direct)  abundance -> cAMP as well, i.e. receptor signals from a pool the
             surface stain does not see.  Prediction: abundance predicts cAMP
             even after surface is conditioned out.

Three traps this script is built around:
 1. cAMP ~ surface is NONLINEAR (saturating). Conditioning linearly biases the
    residual. -> condition with LOWESS, not a line.
 2. The assays share measurement noise (synonymous cross-assay rho ~ 0.55), so a
    naive partial correlation is inflated. -> split-half: estimate the
    conditioning variables from one half of a position's variants and the
    outcome from the other half, making their errors independent.
 3. Conditioning on a NOISY mediator never fully blocks the path (regression
    dilution), which biases the test TOWARD M3. -> so a null here is strong
    evidence for M1, while a positive is ambiguous. Quantified by the
    synonymous-noise floor and a simulation under a true serial chain.
"""
import pathlib

import numpy as np, pandas as pd
from scipy import stats
from statsmodels.nonparametric.smoothers_lowess import lowess

RNG = np.random.default_rng(20261022)
ROOT = pathlib.Path(__file__).resolve().parents[2]
d = pd.read_csv(ROOT / "data/processed/merged_variant_scores.tsv", sep="\t")
mis = d[d.type == "missense"].dropna(subset=["cre_z", "surface_z", "abundance_z"])
syn = d[d.type == "synonymous"].dropna(subset=["cre_z", "surface_z", "abundance_z"])
print(f"missense {len(mis)}, synonymous {len(syn)}")


def resid_lowess(y, x, frac=0.3):
    """y with its LOWESS dependence on x removed."""
    o = np.argsort(x)
    fit = lowess(y[o], x[o], frac=frac, return_sorted=False)
    r = np.empty_like(y, dtype=float)
    r[o] = y[o] - fit
    return r


def partial(y, x, z, frac=0.3):
    """Spearman of y and x, each with its LOWESS dependence on z removed."""
    return stats.spearmanr(resid_lowess(y, z, frac), resid_lowess(x, z, frac))


# ---------------------------------------------------------------- 1. shared noise
print("\n=== 1. the noise floor: synonymous variants (no true effect) ===")
for a, b in [("surface_z", "cre_z"), ("abundance_z", "cre_z"), ("abundance_z", "surface_z")]:
    r = stats.spearmanr(syn[a], syn[b])
    print(f"  syn rho({a[:-2]}, {b[:-2]:8s}) = {r.statistic:+.3f}  p={r.pvalue:.2g}")
print("  -> any apparent association up to this size can be pure shared measurement noise")

# ---------------------------------------------------------------- 2. position level
print("\n=== 2. position level, naive (shared noise NOT removed) ===")
P = mis.groupby("position").agg(cre=("cre_z", "mean"), sur=("surface_z", "mean"),
                                abu=("abundance_z", "mean"), n=("cre_z", "size"))
P = P[P.n >= 10]
print(f"  {len(P)} positions, median {P.n.median():.0f} variants each")
r_sa = stats.spearmanr(P.sur, P.cre); r_aa = stats.spearmanr(P.abu, P.cre)
print(f"  rho(surface, cAMP)   = {r_sa.statistic:+.3f}")
print(f"  rho(abundance, cAMP) = {r_aa.statistic:+.3f}")
pa = partial(P.cre.values, P.abu.values, P.sur.values)
ps = partial(P.cre.values, P.sur.values, P.abu.values)
print(f"  PARTIAL rho(abundance, cAMP | surface) = {pa.statistic:+.3f}  p={pa.pvalue:.2g}")
print(f"  PARTIAL rho(surface, cAMP | abundance) = {ps.statistic:+.3f}  p={ps.pvalue:.2g}")

# ---------------------------------------------------------------- 3. split-half
print("\n=== 3. split-half: predictors and outcome from DIFFERENT variants ===")
print("    (kills the shared-noise inflation; attenuation still biases toward M3)")


def split_half_once(rng):
    rows = []
    for pos, g in mis.groupby("position"):
        if len(g) < 10:
            continue
        idx = rng.permutation(len(g))
        A, B = g.iloc[idx[: len(g) // 2]], g.iloc[idx[len(g) // 2:]]
        rows.append((pos, A.surface_z.mean(), A.abundance_z.mean(), B.cre_z.mean(),
                     B.surface_z.mean()))
    return pd.DataFrame(rows, columns=["pos", "surA", "abuA", "creB", "surB"])


out = {"sur_cre": [], "abu_cre": [], "abu_given_sur": [], "sur_given_abu": []}
for _ in range(40):
    S = split_half_once(RNG)
    out["sur_cre"].append(stats.spearmanr(S.surA, S.creB).statistic)
    out["abu_cre"].append(stats.spearmanr(S.abuA, S.creB).statistic)
    out["abu_given_sur"].append(partial(S.creB.values, S.abuA.values, S.surA.values).statistic)
    out["sur_given_abu"].append(partial(S.creB.values, S.surA.values, S.abuA.values).statistic)
for k, v in out.items():
    v = np.array(v)
    print(f"  {k:16s} rho = {v.mean():+.3f}  (SD over 40 splits {v.std():.3f})")

# ---------------------------------------------------------------- 4. how biased is it?
print("\n=== 4. simulation: a TRUE serial chain, measured with this data's noise ===")
print("    how large a spurious 'abundance | surface' does attenuation alone create?")
n = len(P)
# measurement SE of a position mean, in z units, from the real data
se_pos = {a: (mis.groupby("position")[f"{a}_z"].std() / np.sqrt(mis.groupby("position").size())).median()
          for a in ["cre", "surface", "abundance"]}
print("    position-mean SE (z):", {k: round(float(v), 3) for k, v in se_pos.items()})
sp = []
for _ in range(200):
    true_abu = RNG.normal(0, P.abu.std(), n)
    true_sur = 0.8 * true_abu + RNG.normal(0, 0.6 * P.sur.std(), n)   # abundance -> surface
    true_cre = 0.8 * true_sur + RNG.normal(0, 0.6 * P.cre.std(), n)   # surface -> cAMP ONLY
    o_abu = true_abu + RNG.normal(0, se_pos["abundance"], n)
    o_sur = true_sur + RNG.normal(0, se_pos["surface"], n)
    o_cre = true_cre + RNG.normal(0, se_pos["cre"], n)
    sp.append(partial(o_cre, o_abu, o_sur).statistic)
sp = np.array(sp)
print(f"    spurious partial rho(abundance, cAMP | surface) under a TRUE serial chain:")
print(f"      mean {sp.mean():+.3f}, 95 % range [{np.percentile(sp,2.5):+.3f}, {np.percentile(sp,97.5):+.3f}]")
obs = np.mean(out["abu_given_sur"])
print(f"    observed split-half value: {obs:+.3f}")
print(f"    -> {'INSIDE' if np.percentile(sp,2.5) <= obs <= np.percentile(sp,97.5) else 'OUTSIDE'}"
      " the serial-chain expectation")

# ---------------------------------------------------------------- 5. where it matters
print("\n=== 5. the subset that actually decides it: low surface, is cAMP alive? ===")
low = mis[mis.surface_z < -2]
for lo, hi, lab in [(-3, -2, "mild"), (-5, -3, "moderate"), (-99, -5, "severe")]:
    s = low[(low.surface_z > lo) & (low.surface_z <= hi)]
    if len(s) < 20:
        continue
    alive = s.cre_z > -2
    print(f"  surface ({hi},{lo if lo>-99 else '-inf'}] n={len(s):5d}  cAMP alive {alive.mean():.3f}"
          f"   their abundance median {s.abundance_z.median():+.2f}"
          f" (alive {s[alive].abundance_z.median():+.2f} vs dead {s[~alive].abundance_z.median():+.2f})")
print("  if signalling came from an intracellular pool, the ALIVE ones should hold"
      " abundance while the dead ones lose it")

# ---------------------------------------------------------------- 6. model comparison
print("\n=== 6. out-of-fold prediction of position-mean cAMP (10-fold, gradient boosting) ===")
from sklearn.model_selection import KFold
from sklearn.ensemble import GradientBoostingRegressor

POS = mis.groupby("position").agg(cre=("cre_z", "mean"), sur=("surface_z", "mean"),
                                  abu=("abundance_z", "mean"), n=("cre_z", "size"))
POS = POS[POS.n >= 10]
yv = POS.cre.values
r2 = {}
for lab, cols in [("M2  surface only", ["sur"]), ("M3  abundance only", ["abu"]),
                  ("M1  surface + abundance", ["sur", "abu"])]:
    Xv, oof = POS[cols].values, np.zeros(len(yv))
    for tr, te in KFold(10, shuffle=True, random_state=0).split(Xv):
        g = GradientBoostingRegressor(random_state=0, max_depth=2, n_estimators=200)
        oof[te] = g.fit(Xv[tr], yv[tr]).predict(Xv[te])
    r2[lab] = 1 - ((yv - oof) ** 2).sum() / ((yv - yv.mean()) ** 2).sum()
    print(f"  {lab:26s} R^2 = {r2[lab]:.4f}")
print(f"  adding abundance to surface buys {r2['M1  surface + abundance'] - r2['M2  surface only']:+.4f} R^2")
print(f"  adding surface to abundance buys {r2['M1  surface + abundance'] - r2['M3  abundance only']:+.4f} R^2")

# ---------------------------------------------------------------- 7. the bin-width trap
print("\n=== 7. why QUANTILE surface bins mislead (keep this; it nearly produced a result) ===")
edges = np.percentile(POS.sur, np.linspace(0, 100, 9))
wid, r_ab, r_su = [], [], []
for lo, hi in zip(edges[:-1], edges[1:]):
    s = POS[(POS.sur > lo) & (POS.sur <= hi)]
    if len(s) < 25:
        continue
    wid.append(hi - lo)
    r_ab.append(stats.spearmanr(s.abu, s.cre).statistic)
    r_su.append(stats.spearmanr(s.sur, s.cre).statistic)
    print(f"  ({lo:+.2f},{hi:+.2f}] width {hi-lo:4.2f}  rho(abu,cAMP)={r_ab[-1]:+.3f}"
          f"  rho(sur,cAMP)={r_su[-1]:+.3f}")
print(f"  rank corr( bin width, rho(abu,cAMP) )        = {stats.spearmanr(wid, r_ab).statistic:+.3f}")
print(f"  rank corr( rho(sur,cAMP), rho(abu,cAMP) )    = {stats.spearmanr(r_su, r_ab).statistic:+.3f}")
print("  -> the apparent abundance signal IS leftover within-bin surface variation")
tot = []
for lo in np.arange(-3.0, 2.0, 0.5):
    s = POS[(POS.sur > lo) & (POS.sur <= lo + 0.5)]
    if len(s) >= 30:
        tot.append((len(s), stats.spearmanr(s.abu, s.cre).statistic))
wv = np.array([t[0] for t in tot]); rv = np.array([t[1] for t in tot])
print(f"  FIXED-width 0.5 z bins, weighted mean rho(abu,cAMP) = {(wv*rv).sum()/wv.sum():+.3f} (n={wv.sum()})")
