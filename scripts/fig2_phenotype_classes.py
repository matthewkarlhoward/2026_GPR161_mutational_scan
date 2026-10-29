#!/usr/bin/env python3
"""Variant effects across the three phenotypes, and the five mechanism classes.

Panel letters below are manuscript Figure 3's, the ones on the assembled figure;
they are added in Illustrator, so the axis variables here (axA..axD) are layout
names, not letters. Panel a, the phenotype cascade schematic, is drawn by hand
and is not produced here.

  b  abundance vs surface, every missense variant   (axA)
  c  GPCRdb snake, position class, with the colour key for the whole figure (axD)
  d  cAMP vs surface, every missense variant        (axB)
  e  the chain as a flow: all 528 positions through the three gates         (axC)

Units. Panels b and d are in **DMS effect score**, the quantity the assay
reports. The class cuts are at two synonymous SDs, which is a different raw value
on each axis because each assay has its own synonymous spread:

    cAMP      SD 0.1858   cut at -0.372
    abundance SD 0.1692   cut at -0.338
    surface   SD 0.2212   cut at -0.442

so the dashed lines are drawn per axis rather than at a common number. Panels c
and e are position-level and apply the same cuts to the position mean (the gate
rule is style.POSITION_GATE; docs/METHODS.md section 5 states the alternative
that was rejected and why).

Class rules. Variant (b, d): the `class` column of merged_variant_scores.tsv --
cAMP and surface are cut first and abundance is consulted only to split
destabilized from trafficking defect. Position (c, e): a gate is affected when the
position mean is past its cut, and the position takes the most upstream affected
gate (abundance > surface > cAMP). The two rules share the five class names and
colours, which is what the single key under the snake refers to; applied to the
same unit they would disagree on 2.2 % of variants.

LEGEND. The caption lives with the manuscript, not in this repository; keep the two
in step. The one trap it has to spell out: panel e labels PATHS, not classes, and the
destabilized class is split across two of them. The ribbon reading "Destabilized (30)"
is the (1,1,1) path; the other 2 destabilized positions take (1,1,0) -- they lose
abundance and surface but keep cAMP -- and fall under min_label, so they are unlabelled.
A reader who takes 30 as the class size is off by 2 against the position-class table
printed below. The other thing worth stating is that b and d are coloured by the
VARIANT rule while c and e use the POSITION rule; they share names and colours and
would disagree on 2.2 % of variants applied to the same unit.
"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from style import *
from scipy.stats import spearmanr

# exact house palette (user swatch; scripts/qc/gpr161_style.py)
RED, BLUE, PURPLE, GREEN = "#8F2924", "#4877B5", "#65159E", "#316B3E"
GREY1, GREY2, GREY3, GREY4 = "#808180", "#979798", "#CDCECE", "#E2E2E2"
COL = {"WT-like": GREY3, "uncoupled": RED, "trafficking defect": BLUE,
       "destabilized": PURPLE, "surface-loss only": GREEN}
KEY_ORDER = ["uncoupled", "trafficking defect", "destabilized", "surface-loss only", "WT-like"]
DISP = {"uncoupled": "Uncoupled", "trafficking defect": "Trafficking defect",
        "destabilized": "Destabilized", "surface-loss only": "Surface-loss only",
        "WT-like": "WT-like"}
# WT-like's fill is too light to read as 6 pt text, so its LABEL uses a darker grey
# of the same family; the ribbon and swatch keep the fill colour.
TXT = dict(COL); TXT["WT-like"] = GREY1
DRAW = ["WT-like", "surface-loss only", "destabilized", "trafficking defect", "uncoupled"]
SHORT = {"WT-like": "WT-like", "signaling-selective (expressed, uncoupled)": "uncoupled",
         "trafficking-defective (made, not delivered)": "trafficking defect",
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
# Three columns, two rows; the snake spans both rows on the right.
#
#     cartoon   |  surface vs abundance |  snake
#     surface vs cAMP |  chain flow     |  (snake continues)  + colour key
#
# No panel letters are drawn. The flow's NODE BAND is made to span exactly the
# same y range as the cAMP scatter's plot box, so the two bottom panels line up
# top and bottom; its stage labels sit in the row gap above and the axes keep a
# strip below for the class counts.
W, H = 165.0, 92.0
fig = plt.figure(figsize=(W * MM, H * MM))


def add_ax(x, top, w, h):
    return fig.add_axes([x / W, 1 - (top + h) / H, w / W, h / H])


SNAKE_ASPECT = 0.824                 # w/h of the GPCRdb layout at r = 10.5
R1, R2 = 4.0, 54.0                   # row tops
CART = 44.0                          # cartoon slot, placed in Illustrator
SC = 30.0                            # scatter axes, square
X_SC1, X_SC2 = 66.0, 14.0            # scatter axes left edges (col 2, col 1)
FLOW_W = 38.0
LBL, NOTE = 4.0, 5.0                 # strips for the stage labels and the counts
HS = LBL + SC + NOTE                 # flow axes box height
X3 = X_SC1 + FLOW_W + 4              # snake column
SN_W = W - 3 - X3
SN_H = SN_W / SNAKE_ASPECT

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


# --------------------------------------------------------- cartoon slot (a) --
# Reserved and left blank; the schematic is placed in Illustrator.
add_ax(3, R1, CART, CART).set_axis_off()

# ----------------------------------------------------------- the scatters ---
# Surface is on x in BOTH panels: it is the shared reference the classes are cut
# against, and the pair order (abundance then cAMP) matches the flow's left-to-
# right order. Note this plots abundance, which is UPSTREAM of surface, on y --
# the panel is a threshold map, not a regression, so the direction is not implied.
scatter(add_ax(X_SC1, R1, SC, SC), "surface", "abundance")
scatter(add_ax(X_SC2, R2, SC, SC), "surface", "cre")

# ------------------------------------------------------------- chain flow ---
ALLPOS = set(posn.position)
PATHS = gate_paths(GATES, ALLPOS)
assert sum(q["n"] for q in PATHS) == len(ALLPOS)

# Strictly proportional: every ribbon AND every node box is its count times one
# scale, which is what lets a reader check 433:95 against the drawing. An earlier
# draft compressed the WT-like band to make the five classes thicker; that also
# rescaled the node boxes (Surface came out 1.01:1 while labelled 433:95), so it
# was dropped. The cost is that the uncoupled ribbon is only ~0.5 mm.
# top/bot put the node band exactly on the cAMP scatter's plot box.
sankey_axes(add_ax(X_SC1, R2 - LBL, FLOW_W, HS), PATHS,
            ["Abundance", "Surface", "cAMP"], COL, min_label=10,
            row_labels=("WT-like", "LOF"),
            top=(NOTE + SC) / HS, bot=NOTE / HS,
            x0=0.13, x1=0.82, node_w=0.062, label_gap=2.3 / HS,
            label_fmt=lambda q: (str(q["n"]), TXT[q["cls"]]))
n_any = len(GATES["A"] | GATES["S"] | GATES["C"])
n_A_only = len(GATES["A"] - GATES["S"] - GATES["C"])
n_A_noS = len(GATES["A"] - GATES["S"])

# ---------------------------------------------------------- snake + key -----
# The flow labels ribbons with counts only, so this key is the single place that
# maps colour to class name, for every panel.
snake_axes(add_ax(X3, R1, SN_W, SN_H), categorical=pcls, catmap=COL, r=10.5)

axK = add_ax(X3 - 2, R1 + SN_H + 4, SN_W + 4, 10.5)
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
