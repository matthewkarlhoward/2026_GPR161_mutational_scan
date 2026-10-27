"""DMS library QC: per-position coverage and designed-vs-observed variant table.

Uses the untreated (library input) samples in
    data/counts/*untreat*.csv.gz

Each counts file carries the complete designed variant set (11,043 rows,
identical across samples) with a `count` column, so unobserved variants are
present as zeros.  Variant classes: M missense, S synonymous, D single-codon
deletion.

Outputs (tables/):
    qc_position_coverage.csv  per-position totals
    qc_designed_observed.csv  the variant-class table
    qc_dropout_variants.csv   designed variants under the read threshold

With --figures it also writes figures/F_QC_coverage.{pdf,png}, a larger
standalone version of panels a-b. The committed QC figure is the single 165 x
30 mm row drawn by dms_qc_panel.py, which imports the loaders below.
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
COUNTS = ROOT / "data" / "counts"
OLIGOS = ROOT / "data" / "library" / "GPR161_oligos.csv"
FIGS = ROOT / "figures"
TABS = ROOT / "tables"

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gpr161_style as S  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402

CLASS_LABEL = {"M": "Missense", "S": "Synonymous", "D": "Deletion"}
MIN_READS = 10  # second "observed" threshold


# ---------------------------------------------------------------- data ------
def load_untreated():
    """Return (long dataframe of counts, per-sample read totals)."""
    files = sorted(COUNTS.glob("*untreat*.csv.gz"))
    if not files:
        raise SystemExit(f"no untreated counts files under {COUNTS}")
    frames = []
    for f in files:
        d = pd.read_csv(f)
        m = re.match(r"(R\d)_gpr161_untreat_\d+(_mini)?_", f.name)
        d["replicate"] = m.group(1)
        d["run"] = "miniseq" if m.group(2) else "nextseq"
        d["sample"] = f.name.removesuffix(".csv.gz")
        frames.append(d)
    long = pd.concat(frames, ignore_index=True)
    totals = long.groupby("sample")["count"].sum().sort_index()
    return long, totals


def tile_bounds():
    """First position of each oligo sub-library (DMS-1 ... DMS-11)."""
    o = pd.read_csv(OLIGOS, header=None, names=["name", "seq"], encoding="utf-8-sig")
    ext = o["name"].str.extract(r"GPR161_DMS-(\d+)_[A-Za-z]{3}(\d+)[A-Za-z]{3}$")
    ext.columns = ["tile", "pos"]
    ext = ext.dropna().astype(int)
    g = ext.groupby("tile")["pos"].agg(["min", "max"]).sort_index()
    return g


# ------------------------------------------------------------- plotting -----
def coverage_axis(ax, pos, total, bounds, ylabel="Total count", title=None):
    ax.bar(pos, total, width=1.0, color=S.GREY2, edgecolor="none", zorder=2)
    ax.set_yscale("log")
    lo = total[total > 0].min() if (total > 0).any() else 1
    bottom = max(1.0, 10 ** np.floor(np.log10(lo)))
    ax.set_ylim(bottom, 10 ** np.ceil(np.log10(total.max())))
    ax.axhline(np.median(total), color=S.BLACK, lw=S.LW, ls=(0, (1, 2)), zorder=4)
    ax.set_xlim(pos.min() - 2, pos.max() + 2)
    for b in bounds["min"].iloc[1:]:
        ax.axvline(b - 0.5, color=S.BLACK, lw=S.LW, ls=(0, (2, 2)), zorder=3)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title, loc="left")
    ax.tick_params(which="minor", left=True)
    return ax


def table_axis(ax, tab, bbox=(0, 0.30, 1, 0.42)):
    ax.axis("off")
    cols = ["", "Total"] + [CLASS_LABEL[c] for c in ("S", "M", "D")]
    rows = [
        ["Designed"] + [f"{int(tab.loc['designed', c]):,}" for c in ["Total", "S", "M", "D"]],
        ["Observed"] + [
            f"{int(tab.loc['observed', c]):,} ({tab.loc['pct', c]:.1f}%)"
            for c in ["Total", "S", "M", "D"]
        ],
    ]
    t = ax.table(cellText=rows, colLabels=cols, cellLoc="left",
                 bbox=list(bbox))
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


# ------------------------------------------------------------------ main ----
def main(figures=False):
    long, totals = load_untreated()
    bounds = tile_bounds()

    print("reads per untreated sample:")
    for k, v in totals.items():
        print(f"  {k:<40s} {v:>10,}")
    print(f"  {'TOTAL':<40s} {totals.sum():>10,}\n")

    # ---- pooled counts per designed variant --------------------------------
    pooled = (long.groupby(["name", "pos", "mutation_type"], as_index=False)["count"]
                  .sum())
    assert len(pooled) == 11043, len(pooled)

    # ---- per-position coverage ---------------------------------------------
    cov = pooled.groupby("pos")["count"].sum().sort_index()
    cov_df = cov.rename("total_count").reset_index()
    cov_df["n_variants"] = pooled.groupby("pos").size().values
    cov_df["n_observed"] = pooled.assign(obs=pooled["count"] > 0).groupby("pos")["obs"].sum().values
    cov_df["median_variant_count"] = pooled.groupby("pos")["count"].median().values
    cov_df.to_csv(TABS / "qc_position_coverage.csv", index=False)

    print(f"positions: {cov.index.min()}-{cov.index.max()} (n={len(cov)})")
    print(f"per-position total reads: median {cov.median():,.0f}, "
          f"min {cov.min():,.0f} (pos {cov.idxmin()}), max {cov.max():,.0f} (pos {cov.idxmax()})")
    print(f"fold-range across positions: {cov.max()/cov.min():.1f}x\n")

    # ---- designed vs observed ----------------------------------------------
    tab = pd.DataFrame(index=["designed", "observed", "pct", "observed10", "pct10"],
                       columns=["Total", "S", "M", "D"], dtype=float)
    for c in ["Total", "S", "M", "D"]:
        sub = pooled if c == "Total" else pooled[pooled.mutation_type == c]
        tab.loc["designed", c] = len(sub)
        tab.loc["observed", c] = int((sub["count"] > 0).sum())
        tab.loc["observed10", c] = int((sub["count"] >= MIN_READS).sum())
    tab.loc["pct"] = 100 * tab.loc["observed"] / tab.loc["designed"]
    tab.loc["pct10"] = 100 * tab.loc["observed10"] / tab.loc["designed"]
    tab.loc[["designed", "observed", "observed10"]] = \
        tab.loc[["designed", "observed", "observed10"]].astype(int)
    tab.rename(columns=CLASS_LABEL).to_csv(TABS / "qc_designed_observed.csv")
    missing = (pooled[pooled["count"] < MIN_READS]
               .assign(mutation_type=lambda d: d.mutation_type.map(CLASS_LABEL))
               .sort_values(["pos", "name"])[["name", "pos", "mutation_type", "count"]])
    missing.to_csv(TABS / "qc_dropout_variants.csv", index=False)
    print(f"{len(missing)} designed variants with <{MIN_READS} pooled reads "
          f"({(missing['count'] == 0).sum()} with zero)\n")
    print(tab.rename(columns=CLASS_LABEL).to_string(float_format=lambda v: f"{v:,.1f}"), "\n")

    if not figures:
        print("wrote:")
        for q in ["tables/qc_position_coverage.csv", "tables/qc_designed_observed.csv",
                  "tables/qc_dropout_variants.csv"]:
            print(f"  {ROOT / q}")
        return

    # ---- figure: coverage + table, side by side ----------------------------
    # Only on request: the committed QC figure is the single row drawn by
    # dms_qc_panel.py, and this is a larger standalone version of panels a-b.
    S.apply()
    fig = S.figure(width_mm=165, height_mm=30)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.05],
                          left=0.062, right=0.995, top=0.86, bottom=0.28, wspace=0.10)
    ax = fig.add_subplot(gs[0])
    coverage_axis(ax, cov.index.values, cov.values, bounds,
                  ylabel="Total count", title="Library coverage (untreated input)")
    ax.set_xlabel("Residue position")
    axt = fig.add_subplot(gs[1])
    table_axis(axt, tab)
    S.save(fig, FIGS / "F_QC_coverage")
    plt.close(fig)

    print("wrote:")
    for p in ["figures/F_QC_coverage.pdf", "figures/F_QC_coverage.png",
              "tables/qc_position_coverage.csv", "tables/qc_designed_observed.csv",
              "tables/qc_dropout_variants.csv"]:
        print(f"  {ROOT / p}")


if __name__ == "__main__":
    main(figures="--figures" in sys.argv)
