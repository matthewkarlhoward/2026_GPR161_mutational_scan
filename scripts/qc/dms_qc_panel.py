"""DMS QC: all four QC plots in a single 165 x 30 mm row.

    library coverage per residue position (untreated input)
    designed vs observed variant table
    Lilace score and standard error, one violin per assay
    mean sequencing coverage per sample, one row per assay

Axes are placed in millimetres so the single-row geometry is exact.  Labels are
abbreviated relative to the standalone figures because four panels across
165 mm leaves little room.  Reuses the loaders in dms_qc_coverage.py and
dms_qc_scores.py.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
FIGS = ROOT / "figures"

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gpr161_style as S  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker as mticker  # noqa: E402

import dms_qc_coverage as C  # noqa: E402
import dms_qc_scores as Q  # noqa: E402

W_MM, H_MM = 165.0, 30.0
CLASS_ROWS = [("Total", "Total"), ("S", "Synonymous"), ("M", "Missense"), ("D", "Deletion")]


def mm_axes(fig, x, y, w, h):
    """Axes at (x, y) mm from the TOP-LEFT corner, size w x h mm."""
    return fig.add_axes([x / W_MM, 1 - (y + h) / H_MM, w / W_MM, h / H_MM])


def transposed_table(ax, tab):
    """Variant classes as rows -- narrower than classes-as-columns."""
    ax.axis("off")
    cols = ["", "Designed", "Observed", "%"]
    rows = [[label,
             f"{int(tab.loc['designed', key]):,}",
             f"{int(tab.loc['observed', key]):,}",
             f"{tab.loc['pct', key]:.1f}"]
            for key, label in CLASS_ROWS]
    t = ax.table(cellText=rows, colLabels=cols, cellLoc="left",
                 bbox=[0, 0.08, 1, 0.84])
    t.auto_set_font_size(False)
    t.set_fontsize(S.FONT_PT)
    t.auto_set_column_width(col=list(range(len(cols))))
    for (r, c), cell in t.get_celld().items():
        cell.set_linewidth(S.LW)
        cell.set_edgecolor(S.BLACK)
        cell.PAD = 0.03
        if r == 0 and c == 0:
            cell.set_linewidth(0)
    return ax


def main():
    # ------------------------------------------------------------ data -----
    long, _ = C.load_untreated()
    bounds = C.tile_bounds()
    pooled = long.groupby(["name", "pos", "mutation_type"], as_index=False)["count"].sum()
    cov = pooled.groupby("pos")["count"].sum().sort_index()

    tab = pd.DataFrame(index=["designed", "observed", "pct"],
                       columns=["Total", "S", "M", "D"], dtype=float)
    for c in ["Total", "S", "M", "D"]:
        sub = pooled if c == "Total" else pooled[pooled.mutation_type == c]
        tab.loc["designed", c] = len(sub)
        tab.loc["observed", c] = int((sub["count"] > 0).sum())
    tab.loc["pct"] = 100 * tab.loc["observed"] / tab.loc["designed"]

    scores = Q.load_scores()
    samp = Q.load_sample_coverage()

    # ------------------------------------------------------------ figure ---
    S.apply()
    fig = S.figure(width_mm=W_MM, height_mm=H_MM)

    # --- coverage track ---------------------------------------------------
    ax_a = mm_axes(fig, 11, 5.5, 34, 15)
    C.coverage_axis(ax_a, cov.index.values, cov.values, bounds, ylabel="Total count")
    ax_a.set_xlabel("Residue position", labelpad=1.5)
    ax_a.set_xticks([0, 250, 500])
    sup = {0: "\u2070", 1: "\u00b9", 2: "\u00b2", 3: "\u00b3", 4: "\u2074",
           5: "\u2075", 6: "\u2076"}
    ax_a.yaxis.set_major_formatter(mticker.FuncFormatter(
        lambda v, _: "10" + sup[int(round(np.log10(v)))]))

    # --- designed / observed table ---------------------------------------
    ax_b = mm_axes(fig, 50, 4.0, 36, 22)
    transposed_table(ax_b, tab)

    # --- SE violins ---------------------------------------------------------
    keys = list(Q.ASSAYS)
    labels = ["cAMP", "Surf", "Abund"]
    pos = list(range(1, len(keys) + 1))
    se = [scores.loc[scores.assay == k, "effect_se"].to_numpy() for k in keys]

    ax_c = mm_axes(fig, 96, 4.0, 20, 18)
    Q.strip(ax_c, se, pos)
    Q.violin(ax_c, se, pos)
    ax_c.scatter(pos, [np.median(e) for e in se], s=3, color=S.BLACK, zorder=5)
    ax_c.set_ylabel("SE", labelpad=1.5)
    ax_c.set_ylim(0.2, 0.55)
    ax_c.set_yticks([0.2, 0.3, 0.4, 0.5])
    ax_c.set_xticks(pos)
    ax_c.set_xticklabels(labels, rotation=30, ha="right")
    ax_c.set_xlim(pos[0] - 0.8, pos[-1] + 0.8)

    # --- per-sample sequencing coverage ----------------------------------
    order = [("cre", "cAMP"), ("surface", "Surf"),
             ("abundance", "Abund"), ("untreated", "Input")]
    ax_d = mm_axes(fig, 131, 4.0, 31, 18)
    yticks, ylabels = [], []
    for i, (key, label) in enumerate(order):
        y = len(order) - i
        d = samp.loc[samp.assay == key, "coverage"].to_numpy()
        colour = Q.ASSAYS[key][2] if key in Q.ASSAYS else S.GREY1
        ax_d.scatter(d, y + Q.RNG.uniform(-0.18, 0.18, len(d)), s=2,
                     color=colour, linewidths=0, zorder=3)
        ax_d.plot([np.median(d)] * 2, [y - 0.33, y + 0.33], color=S.BLACK,
                  lw=S.LW, zorder=4)
        yticks.append(y)
        ylabels.append(f"{label} ({len(d)})")
    ax_d.set_xscale("log")
    ax_d.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax_d.xaxis.set_minor_formatter(mticker.NullFormatter())
    ax_d.set_yticks(yticks)
    ax_d.set_yticklabels(ylabels)
    ax_d.set_ylim(0.4, len(order) + 0.6)
    ax_d.set_xlim(10, 1000)
    ax_d.set_xlabel("Mean coverage (×)", labelpad=1.5)
    ax_d.spines["left"].set_visible(False)
    ax_d.tick_params(axis="y", length=0)

    S.save(fig, FIGS / "F_QC_panel")
    plt.close(fig)
    print(f"  {FIGS / 'F_QC_panel.pdf'}")
    print(f"  {FIGS / 'F_QC_panel.png'}")


if __name__ == "__main__":
    main()
