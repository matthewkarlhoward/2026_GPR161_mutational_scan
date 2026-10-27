#!/usr/bin/env python3
"""Variant effects across the three phenotypes, and the five mechanism classes.

Panel letters below are the ones on the assembled figure; they are added in
Illustrator, so the axis variables here (axA..axD) are layout names, not letters.

  a  abundance vs surface, every missense variant   (axA)
  b  cAMP vs surface, every missense variant        (axB)
  c  GPCRdb snake, position class, with the colour key for the whole figure (axD)
  d  the chain as a flow: all 528 positions through the three gates         (axC)

Units. Panels a and b are in **DMS effect score**, the quantity the assay
reports. The class cuts are at two synonymous SDs, which is a different raw value
on each axis because each assay has its own synonymous spread:

    cAMP      SD 0.1858   cut at -0.372
    abundance SD 0.1692   cut at -0.338
    surface   SD 0.2212   cut at -0.442

so the dashed lines are drawn per axis rather than at a common number. Panels c
and d are position-level and apply the same cuts to the position mean (the gate
rule is style.POSITION_GATE; docs/METHODS.md section 5 states the alternative
that was rejected and why).

Class rules. Variant (a, b): the `class` column of merged_variant_scores.tsv --
cAMP and surface are cut first and abundance is consulted only to split
destabilized from not-delivered. Position (c, d): a gate is affected when the
position mean is past its cut, and the position takes the most upstream affected
gate (abundance > surface > cAMP). The two rules share the five class names and
colours, which is what the single key under the snake refers to; applied to the
same unit they would disagree on 2.2 % of variants.

LEGEND (the text that used to sit above the Sankey; paste into the caption):
    Variant effects across the three phenotypes. (Top left, top middle) All 9,927
    missense variants in DMS effect units. Dashed lines mark two synonymous SDs,
    which is a different effect value in each assay: cAMP -0.372, abundance
    -0.338, surface -0.442. rho is Spearman. (Bottom left) All 528 assayed
    positions routed through the three assays; a gate counts as affected when the
    position mean falls past that assay's cut. 105 positions are affected at one
    gate or more; none loses abundance alone, and none that loses abundance
    spares surface. Ribbon width is the number of positions, and the eight paths
    through the three gates are the five mechanism classes. (Right) GPCRdb snake
    plot coloured by position class.

Whether the assays really act in series is not argued here. It was tested
separately (that code is not in this repository): abundance acts on cAMP only through
surface (partial rho = +0.03, p = 0.49; adding abundance to a surface-only model
of position-mean cAMP buys -0.02 out-of-fold R2), and the dependence of cAMP on
surface is graded rather than a threshold (cAMP survives surface z -2 to -3 in
68 % of variants, -3 to -5 in 26 %, below -5 in 3 %).
"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from style import *
from scipy.stats import spearmanr

# exact house palette (user swatch; gpr161_analysis/scripts/gpr161_style.py)
RED, BLUE, PURPLE, GREEN = "#8F2924", "#4877B5", "#65159E", "#316B3E"
GREY1, GREY2, GREY3, GREY4 = "#808180", "#979798", "#CDCECE", "#E2E2E2"
COL = {"WT-like": GREY3, "uncoupled": RED, "not delivered": BLUE,
       "destabilized": PURPLE, "surface-loss only": GREEN}
KEY_ORDER = ["uncoupled", "not delivered", "destabilized", "surface-loss only", "WT-like"]
DISP = {"uncoupled": "Uncoupled", "not delivered": "Not delivered",
        "destabilized": "Destabilized", "surface-loss only": "Surface-loss only",
        "WT-like": "WT-like"}
# WT-like's fill is too light to read as 6 pt text, so its LABEL uses a darker grey
# of the same family; the ribbon and swatch keep the fill colour.
TXT = dict(COL); TXT["WT-like"] = GREY1
DRAW = ["WT-like", "surface-loss only", "destabilized", "not delivered", "uncoupled"]
SHORT = {"WT-like": "WT-like", "signaling-selective (expressed, uncoupled)": "uncoupled",
         "trafficking-defective (made, not delivered)": "not delivered",
         "destabilized (low abundance -> low surface)": "destabilized",
         "surface-loss, signaling-intact": "surface-loss only"}

d = load_dms()
m = d[d.type == "missense"].copy()
m["mech"] = m["class"].map(SHORT)
posn = load_positions()

GATES = position_gates(posn=posn, variants=m)
posn["mech"] = [position_class(q in GATES["A"], q in GATES["S"], q in GATES["C"])
                for q in posn.position]
pcls = dict(zip(posn.position, posn.mech))

# the z = -2 cut expressed in each assay's own effect units (z is an exact
# per-assay rescaling, so eff/z recovers the synonymous SD to machine precision)
SD = {a: float((m[f"{a}_eff"] / m[f"{a}_z"]).replace([np.inf, -np.inf], np.nan).median())
      for a in ("cre", "abundance", "surface")}
CUT = {a: -2 * v for a, v in SD.items()}

# ---------------------------------------------------------------- layout (mm)
W, H = 165.0, 104.0
fig = plt.figure(figsize=(W * MM, H * MM))


def add_ax(x, top, w, h):
    return fig.add_axes([x / W, 1 - (top + h) / H, w / W, h / H])


LIM = {"cre": (-1.40, 1.10), "abundance": (-1.60, 0.90), "surface": (-1.90, 1.20)}
TICK = {"cre": [-1, -0.5, 0, 0.5, 1], "abundance": [-1.5, -1, -0.5, 0, 0.5],
        "surface": [-1.5, -1, -0.5, 0, 0.5, 1]}
NAME = {"cre": "cAMP", "abundance": "Abundance", "surface": "Surface"}


def scatter(ax, xa, ya):
    for k in DRAW:
        s = m[m.mech == k]
        ax.scatter(s[f"{xa}_eff"], s[f"{ya}_eff"], s=1.5 if k == "WT-like" else 2.1,
                   color=COL[k], alpha=0.45 if k == "WT-like" else 0.85, lw=0, rasterized=True)
    ax.axhline(CUT[ya], color="black", lw=0.5, ls=(0, (3, 2)))
    ax.axvline(CUT[xa], color="black", lw=0.5, ls=(0, (3, 2)))
    ax.set_xlabel(f"{NAME[xa]} effect", labelpad=1.5)
    ax.set_ylabel(f"{NAME[ya]} effect", labelpad=1.5)
    ax.set_xlim(*LIM[xa]); ax.set_ylim(*LIM[ya])
    ax.set_xticks(TICK[xa]); ax.set_yticks(TICK[ya])
    ax.tick_params(pad=1.5, length=2)
    # A key rather than a sentence: a dash of the real line style, then the value.
    # Bottom-right corner, checked empty in both panels (0 points in the left, 1 in
    # the right -- the uncoupled cloud starts above y = 0.12 in axes fraction). The
    # group starts at x = 0.56, well clear of the vertical cut at x = 0.47.
    ax.plot([0.56, 0.66], [0.055, 0.055], transform=ax.transAxes, color="black",
            lw=0.5, ls=(0, (3, 2)), clip_on=False, solid_capstyle="butt")
    ax.text(0.69, 0.055, "2 × syn. SD", transform=ax.transAxes,
            ha="left", va="center")
    r = spearmanr(m[f"{xa}_eff"], m[f"{ya}_eff"])
    ax.text(0.04, 0.97, f"n = {len(m):,}", transform=ax.transAxes, va="top", ha="left")
    ax.text(0.04, 0.89, f"ρ = {r.statistic:+.2f}", transform=ax.transAxes,
            va="top", ha="left")


# Surface is on x in BOTH panels: it is the shared reference the classes are cut
# against, and the pair order (abundance then cAMP) matches the Sankey's left-to-
# right order. Note this plots abundance, which is UPSTREAM of surface, on y --
# the panel is a threshold map, not a regression, so the direction is not implied.
SW, TOPS = 38.0, 8.0
axA = add_ax(12, TOPS, SW, SW)
scatter(axA, "surface", "abundance")
axB = add_ax(62, TOPS, SW, SW)
scatter(axB, "surface", "cre")

# ------------------------------------------------- d (lettered)  chain flow
TOPC = 56.0
ALLPOS = set(posn.position)
PATHS = gate_paths(GATES, ALLPOS)
assert sum(q["n"] for q in PATHS) == len(ALLPOS)

axC = add_ax(4, TOPC, 94, 44)
sankey_axes(axC, PATHS, ["Abundance", "Surface", "cAMP"], COL, min_label=10,
            row_labels=("Spared", "Affected"),
            label_fmt=lambda q: (f"{DISP[q['cls']]} ({q['n']})", TXT[q["cls"]]))
n_any = len(GATES["A"] | GATES["S"] | GATES["C"])
n_A_only = len(GATES["A"] - GATES["S"] - GATES["C"])
n_A_noS = len(GATES["A"] - GATES["S"])


# ------------------------------------------------- c (lettered)  snake + key
SX, SWID = 99.0, 61.0
axD = add_ax(SX, 4.0, SWID, 84.0)
snake_axes(axD, categorical=pcls, catmap=COL, r=10.5)

# key: two columns, 3 + 2, directly under the snake
axK = add_ax(SX + 3, 90.0, SWID - 2, 10.5)
axK.set_xlim(0, 2); axK.set_ylim(3, 0); axK.axis("off")
for i, k in enumerate(KEY_ORDER):
    cx, cy = divmod(i, 3)
    axK.scatter([cx + 0.04], [cy + 0.5], s=14, color=COL[k], lw=0.4, edgecolors="black")
    axK.text(cx + 0.14, cy + 0.5, DISP[k], va="center")

save_figure(fig, "phenotype_classes")

print("effect-unit cuts (two synonymous SDs):", {a: round(v, 4) for a, v in CUT.items()})
print("\nvariant classes (a, b):")
print(m.mech.value_counts().reindex(KEY_ORDER).to_string())
print(f"\nposition classes (c, d) [{POSITION_GATE} gate rule]:")
print(posn.mech.value_counts().reindex(KEY_ORDER).to_string())
print("\npaths through the three gates (abundance, surface, cAMP affected):")
for q in sorted(PATHS, key=lambda r: r["flags"]):
    print("  ", q["flags"], f"{q['n']:>4}", q["cls"])
print(f"  gates: A {len(GATES['A'])}  S {len(GATES['S'])}  C {len(GATES['C'])}")
print("\nlegend text for the caption:")
print(f"  {len(ALLPOS)} positions; {gate_label()}. {n_any} affected at one gate or more; "
      f"abundance alone {n_A_only}"
      + ("; none that loses abundance spares surface." if n_A_noS == 0 else "."))
