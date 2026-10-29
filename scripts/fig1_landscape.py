#!/usr/bin/env python3
"""GPR161 landscape: score density + substitution heatmap, one row per phenotype.

Density on the left, heatmap and its colour key on the right. Rows are
abundance / surface / cAMP and share a row height, so both columns sit on the
same baseline for every assay.

Raw DMS effect scores. The raw scale is already centred on the synonymous mean
(+0.0001 abundance / +0.0001 surface / -0.0002 cAMP), so white = 0 is the mean
synonymous variant to within 2e-4.
Diverging scale is BLUE = gain, RED = loss, built from the two palette colours
in style.py (BLUE_GOF / RED_LOF). Note the sense: some earlier drafts of this
project used matplotlib's "RdBu_r", in which blue meant loss.

Final file is exactly 112 x 43.9 mm (the 102 x 40 mm layout scaled up at the
same 2.55:1 aspect). Text stays 6 pt and lines 0.5 pt per the style guide; the
panels, gaps and colour key scale, and the margins that only hold text keep
their mm size.
"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from style import *
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
from scipy.stats import gaussian_kde

W_MM = 112.0
H_MM = W_MM * 40.0 / 102.0              # same aspect as the original 102 x 40 mm
# margins that hold 6 pt text: fixed in mm, independent of the canvas size
L_MM, R_MM, B_MM, T_MM = 2.8, 7.3, 5.4, 1.52
# the 102 mm layout between the text margins, scaled to fill the new width:
# density | gap | heatmap | gap | colour key
_D, _G, _HT, _CG, _CB = 27.44, 6.65, 53.51, 3.30, 1.0
SX = (W_MM - L_MM - R_MM) / (_D + _G + _HT + _CG + _CB)
D_MM, GAP_MM, H_HEAT_MM, CB_GAP_MM, CB_W_MM = (v * SX for v in (_D, _G, _HT, _CG, _CB))
CB_X_MM = L_MM + D_MM + GAP_MM + H_HEAT_MM + CB_GAP_MM
AA = list("DEKRHNQSTGPAVILMFYWC")      # charged -> polar -> special -> hydrophobic
DEL_GAP = 0.5                           # blank strip (in row units) between C and del
NROW = len(AA) + DEL_GAP + 1
POS = range(2, 530)
VMAX = 1.2                              # symmetric; clips 2.6% of values, all deep loss
# legend order is NOT draw order: synonymous is the null and reads first, but it
# must still be drawn last so the missense fill does not cover it
LEG_ORDER = ["Synonymous", "Missense", "Deletion"]
ASSAYS = [("Abundance", "abundance"), ("Surface", "surface"), ("cAMP", "cre")]

# red -> white -> blue, from the palette colours rather than matplotlib's RdBu
CMAP = LinearSegmentedColormap.from_list("lof_gof", [RED_LOF, "#ffffff", BLUE_GOF])
CLS = [("Missense", "missense", BLUE_GOF), ("Deletion", "deletion", RED_LOF),
       ("Synonymous", "synonymous", INK)]

d = load_dms()
mis = d[d.type == "missense"]
dele = d[d.type == "deletion"]
COLS = [f"{a}_eff" for _, a in ASSAYS]
DLIM = (min(d[c].min() for c in COLS) * 1.05, max(d[c].max() for c in COLS) * 1.05)
xs = np.linspace(*DLIM, 512)

fig = plt.figure(figsize=(W_MM * MM, H_MM * MM))
gs = GridSpec(3, 2, figure=fig, width_ratios=[D_MM, H_HEAT_MM], hspace=0.13,
              wspace=GAP_MM / ((D_MM + H_HEAT_MM) / 2),
              left=L_MM / W_MM, right=(L_MM + D_MM + GAP_MM + H_HEAT_MM) / W_MM,
              bottom=B_MM / H_MM, top=1 - T_MM / H_MM)

for i, (lab, arm) in enumerate(ASSAYS):
    col = f"{arm}_eff"

    # -------- right: substitution heatmap
    axH = fig.add_subplot(gs[i, 1])
    M = (mis.pivot_table(index="mutation", columns="position", values=col)
            .reindex(index=AA, columns=POS))
    im = axH.imshow(M.values, aspect="auto", cmap=CMAP, vmin=-VMAX, vmax=VMAX,
                    interpolation="nearest", extent=[2, 529, len(AA), 0],
                    rasterized=True)
    # single-residue deletions as one extra row under C, set off by a thin gap
    Dl = dele.set_index("position")[col].reindex(POS).values[None, :]
    axH.imshow(Dl, aspect="auto", cmap=CMAP, vmin=-VMAX, vmax=VMAX,
               interpolation="nearest", extent=[2, 529, NROW, NROW - 1],
               rasterized=True)
    axH.set_ylim(NROW, 0)
    axH.set_yticks([])
    axH.set_xticks([100, 200, 300, 400, 500])

    # -------- left: density of the same scores, all three variant classes
    axD = fig.add_subplot(gs[i, 0])
    for name, typ, colour in CLS:
        y = gaussian_kde(d.loc[d.type == typ, col].dropna().values)(xs)
        if typ == "missense":
            axD.fill_between(xs, y, color=colour, alpha=0.55, lw=0)
        axD.plot(xs, y, color=colour, lw=0.5)
    axD.set_xlim(*DLIM)
    top = max(l.get_ydata().max() for l in axD.lines)
    axD.set_ylim(0, top * (1.55 if i == 0 else 1.06))
    axD.set_yticks([])
    axD.set_ylabel(lab, fontsize=6, labelpad=1.5)
    axD.set_xticks([-1.5, -1.0, -0.5, 0.0, 0.5, 1.0])

    for ax in (axH, axD):
        ax.tick_params(labelsize=6, pad=1.5, length=1.5)
        if i < 2:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel("Residue" if ax is axH else "DMS effect score",
                          fontsize=6, labelpad=1.5)

    if i == 0:
        axD.legend(handles=[Line2D([], [], color=dict((n, c) for n, _, c in CLS)[k],
                                lw=0.5, label=k) for k in LEG_ORDER],
                   loc="upper left", bbox_to_anchor=(0.035, 1.03), fontsize=6,
                   handlelength=1.0, handletextpad=0.30, labelspacing=0.10,
                   borderpad=0.0, borderaxespad=0.0)

# colour key at the far right, spanning all three rows. The density x-axis is
# the same quantity over its full range, so this bar states where the colour
# scale saturates (2.6% of values, all deep loss)
heats = [a for a in fig.axes if a.images]
hb = heats[0].get_position(); hl = heats[-1].get_position()
cax = fig.add_axes([CB_X_MM / W_MM, hl.y0, CB_W_MM / W_MM, hb.y1 - hl.y0])
cb = fig.colorbar(heats[0].images[0], cax=cax, ticks=[-VMAX, 0, VMAX])
cb.ax.set_yticklabels([f"\u2212{VMAX:g}", "0", f"+{VMAX:g}"], fontsize=6)
cb.ax.tick_params(length=1.5, pad=1.5)
cb.set_label("DMS effect score", fontsize=6, labelpad=-4)
cb.outline.set_linewidth(0.5)

for ext in ("png", "pdf"):
    fig.savefig(FIG / f"landscape_density_heat.{ext}", dpi=600)
plt.close(fig)
print(f"wrote landscape_density_heat  ({W_MM:.1f} x {H_MM:.1f} mm, exact)")
