#!/usr/bin/env python3
"""Figure 5 -- interpretation of human GPR161 variants.

  a  where human alleles fall in EACH assay: three columns (abundance, surface,
     cAMP) sharing the z axis, with the full scan, the gnomAD-observed set and
     the ClinVar set drawn as density-spread swarms
  b  how much of each set is loss of function, and by which mechanism

Panel a's universe is missense PLUS single-residue deletions; panel b is
missense only, because the mechanism classes are defined on missense. The two
panels therefore carry slightly different n, which the figure states.

A true non-overlapping beeswarm cannot be packed for 10,425 points in a ~24 mm
column at a readable marker size, so each point's horizontal offset is scaled by
the kernel density at its own y. That shows where the mass is, which is what a
beeswarm is read for, without the packing.

Inputs:  data/processed/merged_variant_scores.tsv   (pipeline/01_merge_scores.py)
         data/annotation/human_alleles.tsv          (see docs/DATA_SOURCES.md)
Output:  figures/human_variants.{pdf,png}

Usage: python fig4_human_variants.py
"""
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from style import *                                    # noqa: F403
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy import stats

ROOT = pathlib.Path(__file__).resolve().parents[1]
SCORES = ROOT / "data" / "processed" / "merged_variant_scores.tsv"
ALLELES = ROOT / "data" / "annotation" / "human_alleles.tsv"

# the five mechanism classes, named as in Figures 2 and 3
RED, BLUE, PURPLE, GREEN = "#8F2924", "#4877B5", "#65159E", "#316B3E"
GREY1, GREY3 = "#808180", "#CDCECE"
MCOL = {"WT-like": GREY3, "uncoupled": RED, "trafficking defect": BLUE,
        "destabilized": PURPLE, "surface-loss only": GREEN}
MECH = ["WT-like", "uncoupled", "trafficking defect", "destabilized",
        "surface-loss only"]
SHORT = {"WT-like": "WT-like",
         "signaling-selective (expressed, uncoupled)": "uncoupled",
         "trafficking-defective (made, not delivered)": "trafficking defect",
         "destabilized (low abundance -> low surface)": "destabilized",
         "surface-loss, signaling-intact": "surface-loss only"}
SET_COL = {"all": GREY3, "gnomAD": BLUE, "ClinVar": RED}
ASSAYS = [("abundance", "Abundance"), ("surface", "Surface"), ("cre", "cAMP")]
LOF_Z = -2.0                      # two synonymous SDs, the cutoff used throughout

# ------------------------------------------------------------------- data ---
d = pd.read_csv(SCORES, sep="\t", low_memory=False)
mis = d[d.type == "missense"].copy()
mis["mech"] = mis["class"].map(SHORT)
mis["allele"] = mis.wildtype + mis.position.astype(str) + mis.mutation
dele = d[d.type == "deletion"].copy()
dele["allele"] = dele.wildtype + dele.position.astype(str) + "del"
UNI = pd.concat([mis, dele], ignore_index=True)

h = pd.read_csv(ALLELES, sep="\t").fillna({"source": "", "clinvar_classification": ""})
h = h[h.consequence_class.isin(["missense", "inframe_indel"])]
GNOMAD = set(h[h.source.str.contains("gnomAD") & h.gnomad_pass_any].protein_label)
CLINVAR = set(h[h.source.str.contains("ClinVar")].protein_label)

# panel a: missense + single-residue deletions. panel b: missense only.
SETS_A = [("all", UNI), ("gnomAD", UNI[UNI.allele.isin(GNOMAD)]),
          ("ClinVar", UNI[UNI.allele.isin(CLINVAR)])]
SETS_B = [("all missense", mis), ("gnomAD", mis[mis.allele.isin(GNOMAD)]),
          ("ClinVar missense", mis[mis.allele.isin(CLINVAR)])]

rng = np.random.default_rng(0)


def _sina_x(v, centre, half):
    """Beeswarm-style offsets: horizontal spread follows the local density."""
    v = np.asarray(v, float)
    if len(v) < 8 or np.allclose(v, v[0]):
        return centre + (rng.random(len(v)) - 0.5) * half
    grid = np.linspace(v.min(), v.max(), 256)
    dens = stats.gaussian_kde(v)(grid)
    w = np.interp(v, grid, dens / dens.max())
    return centre + (rng.random(len(v)) * 2 - 1) * half * w


# ----------------------------------------------------------------- layout ---
# a : b drawn 120 : 40 mm wide. The saved page is narrower than 165 mm because
# bbox_inches="tight" trims to content; the ratio is what is set here.
W, H = 165.0, 62.0
fig = plt.figure(figsize=(W * MM, H * MM))
gs = GridSpec(1, 2, figure=fig, width_ratios=[120, 40], wspace=0.34)


def nice(ax):
    ax.tick_params(labelsize=6, pad=1.5, length=2)


# ------------------------------------------------------------------ a -------
host = fig.add_subplot(gs[0, 0])
spec = host.get_subplotspec()
host.remove()
sub = spec.subgridspec(1, 3, wspace=0.18)
axes = [fig.add_subplot(sub[0, i]) for i in range(3)]
for a_, b_ in zip(axes, axes[1:]):
    b_.sharey(a_)

for ax, (pre, name) in zip(axes, ASSAYS):
    for i, (lab, x) in enumerate(SETS_A):
        v = x[f"{pre}_z"].dropna().values
        ax.scatter(_sina_x(v, i, 0.40), v, s=0.7 if lab == "all" else 1.6,
                   color=SET_COL[lab], lw=0,
                   alpha=0.5 if lab == "all" else 0.8, rasterized=True)
    ax.axhline(0, color=INK, lw=0.5)                    # no effect
    ax.axhline(LOF_Z, color=INK, lw=0.5, ls="--")       # the LOF cut
    ax.set_xticks([])
    ax.set_xlim(-0.6, len(SETS_A) - 0.4)
    ax.set_title(name, fontsize=6, pad=2)
    nice(ax)
    # the sets are named in the legend and the solid line at 0 reads as the
    # baseline, so the x axis carries nothing
    ax.spines["bottom"].set_visible(False)
    ax.tick_params(bottom=False)
axes[0].set_ylabel("Effect (z)")
axes[0].set_ylim(-11, 6)
# the three columns share one z scale: draw it once, on the left
for ax in axes[1:]:
    ax.tick_params(labelleft=False, left=False)
    ax.spines["left"].set_visible(False)
axes[-1].legend(handles=[Line2D([], [], ls="", marker="o", ms=2.2,
                                color=SET_COL[lab], label=f"{lab} ({len(x):,})")
                         for lab, x in SETS_A],
                loc="lower right", bbox_to_anchor=(1.06, -0.02), fontsize=5,
                handletextpad=0.3, labelspacing=0.25, borderpad=0.1,
                frameon=False)

# ------------------------------------------------------------------ b -------
# Was a stack to 100 %, where three ~85 % grey columns read as "these all look
# the same" and hid a 2.5-fold spread. Dropping WT-like and letting the bars top
# out near 21 % shows the depletion AND its composition in one object.
axb = fig.add_subplot(gs[0, 1])
bottom = np.zeros(len(SETS_B))
for k in [c for c in MECH if c != "WT-like"]:
    frac = np.array([100 * (x.mech == k).mean() for _, x in SETS_B])
    axb.bar(range(len(SETS_B)), frac, bottom=bottom, color=MCOL[k], width=0.68,
            edgecolor=INK, lw=0.4, label=k)
    bottom += frac
for i, tot in enumerate(bottom):
    axb.text(i, tot + 0.4, f"{tot:.1f}", ha="center", va="bottom", fontsize=5)
axb.set_xticks(range(len(SETS_B)))
axb.set_xticklabels([q[0] for q in SETS_B], rotation=18, ha="right")
axb.set_ylabel("% Loss-of-function")
axb.set_ylim(0, float(max(bottom)) * 1.3)
axb.legend(loc="upper right", ncol=1, fontsize=4.6, handlelength=0.8,
           handletextpad=0.3, labelspacing=0.18, borderpad=0.25,
           borderaxespad=0.3, frameon=False)
nice(axb)

for ext in ("png", "pdf"):
    fig.savefig(FIG / f"human_variants.{ext}", dpi=600, bbox_inches="tight",
                pad_inches=0.02, facecolor="white")
plt.close(fig)
print(f"wrote human_variants  (panel a {len(UNI):,} variants, "
      f"panel b {len(mis):,} missense)")

print("\npanel a sets (missense + single-residue deletions):")
for lab, x in SETS_A:
    print(f"  {lab:8s} n={len(x):6,}  (deletions {int((x.type == 'deletion').sum())})")
print("\npanel b, % loss of function by set:")
for (lab, x), tot in zip(SETS_B, bottom):
    n_lof = int((x.mech != "WT-like").sum())
    print(f"  {lab:18s} n={len(x):5,}  {n_lof:4d} LOF ({tot:.1f} %)")
