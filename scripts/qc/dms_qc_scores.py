"""DMS QC part 2: Lilace score/SE distributions and per-sample sequencing coverage.

Outputs (tables/):
    qc_score_se_summary.csv   n, median score, median SE and lfsr<0.05 fraction
                              per assay x variant class
    qc_sample_coverage.csv    reads and mean coverage per biological sample

With --figures it also writes larger standalone versions of two panels:

F_QC_score_se.{pdf,png}
    Two-row violin panel (score on top, standard error below), one violin per
    assay over all scored variants.  Grey violin + jittered points, black dot
    at the median.

F_QC_seq_coverage.{pdf,png}
    Mean sequencing coverage per sample (reads per designed variant), one row
    per assay, log x axis, vertical line at the row median.

The committed QC figure is the single row drawn by dms_qc_panel.py, which
imports the loaders below.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
COUNTS = ROOT / "data" / "counts"
LAYOUT = ROOT / "data" / "library" / "GPR161_experiment_layout.csv"
FIGS = ROOT / "figures"
TABS = ROOT / "tables"

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gpr161_style as S  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker as mticker  # noqa: E402

N_DESIGNED = 11043
RNG = np.random.default_rng(0)
N_POINTS = 700          # jittered points drawn per violin

# assay key -> (display label, score table, colour)
ASSAYS = {
    "cre":       ("cAMP",      ROOT / "data/scores/cre_variant_scores.tsv",      S.BLACK),
    "surface":   ("Surface",   ROOT / "data/scores/surface_variant_scores.tsv",   S.BLUE),
    "abundance": ("Abundance", ROOT / "data/scores/abundance_variant_scores.tsv", S.GREEN),
}
CLASSES = [("synonymous", "Syn."), ("missense", "Missense"), ("deletion", "Deletion")]


# ------------------------------------------------------------------ data ----
def load_scores():
    frames = []
    for key, (label, path, _) in ASSAYS.items():
        d = pd.read_csv(path, sep="\t")
        d["assay"] = key
        d["assay_label"] = label
        frames.append(d[["assay", "assay_label", "hgvs", "position", "type",
                         "effect", "effect_se", "lfsr"]])
    return pd.concat(frames, ignore_index=True)


def load_sample_coverage():
    lay = pd.read_csv(LAYOUT)
    lay["assay"] = lay["condition"].str.replace("_mini", "", regex=False)
    lay["run"] = np.where(lay["condition"].str.endswith("_mini"), "miniseq", "nextseq")
    # biological sample = assay + replicate + sort bin, summed over the two runs
    reads = {}
    for f in lay["file"]:
        p = COUNTS / f"{f}.csv.gz"
        reads[f] = pd.read_csv(p, usecols=["count"])["count"].sum()
    lay["reads"] = lay["file"].map(reads)
    samp = (lay.groupby(["assay", "replicate", "time"], as_index=False)["reads"].sum())
    samp["coverage"] = samp["reads"] / N_DESIGNED
    return samp


# -------------------------------------------------------------- plotting ----
def violin(ax, data, pos, width=0.72):
    v = ax.violinplot(data, positions=pos, widths=width, showextrema=False,
                      showmedians=False)
    for b in v["bodies"]:
        b.set_facecolor(S.GREY3)
        b.set_edgecolor(S.BLACK)
        b.set_linewidth(S.LW)
        b.set_alpha(1.0)
    return v


def strip(ax, data, pos, width=0.30):
    for d, x in zip(data, pos):
        y = d if len(d) <= N_POINTS else RNG.choice(d, N_POINTS, replace=False)
        jitter = RNG.uniform(-width, width, len(y))
        ax.scatter(x + jitter, y, s=0.6, color=S.GREY1, alpha=0.35,
                   linewidths=0, zorder=1, rasterized=True)


def score_se_figure(scores):
    keys = list(ASSAYS)
    labels = [ASSAYS[k][0] for k in keys]
    pos = list(range(1, len(keys) + 1))
    eff = [scores.loc[scores.assay == k, "effect"].to_numpy() for k in keys]
    se = [scores.loc[scores.assay == k, "effect_se"].to_numpy() for k in keys]

    S.apply()
    fig = S.figure(width_mm=70, height_mm=55)
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 0.72],
                          left=0.17, right=0.98, top=0.97, bottom=0.14, hspace=0.18)

    ax0 = fig.add_subplot(gs[0])
    strip(ax0, eff, pos)
    violin(ax0, eff, pos)
    ax0.scatter(pos, [np.median(e) for e in eff], s=4, color=S.BLACK, zorder=5)
    ax0.axhline(0, color=S.GREY2, lw=S.LW, ls=(0, (2, 2)), zorder=0)
    ax0.set_ylabel("Score")
    ax0.set_ylim(-2.5, 1.5)

    ax1 = fig.add_subplot(gs[1], sharex=ax0)
    strip(ax1, se, pos)
    violin(ax1, se, pos)
    ax1.scatter(pos, [np.median(e) for e in se], s=4, color=S.BLACK, zorder=5)
    ax1.set_ylabel("SE")
    ax1.set_ylim(0.2, 0.55)
    ax1.set_xticks(pos)
    ax1.set_xticklabels(labels)
    ax0.tick_params(labelbottom=False)

    for a in (ax0, ax1):
        a.set_xlim(pos[0] - 0.8, pos[-1] + 0.8)

    S.save(fig, FIGS / "F_QC_score_se")
    plt.close(fig)


def seq_coverage_figure(samp):
    order = [("cre", "cAMP"), ("surface", "Surface"),
             ("abundance", "Abundance"), ("untreated", "Library input")]
    S.apply()
    fig = S.figure(width_mm=78, height_mm=38)
    ax = fig.add_subplot(111)
    fig.subplots_adjust(left=0.30, right=0.97, top=0.95, bottom=0.21)

    yticks, ylabels = [], []
    for i, (key, label) in enumerate(order):
        y = len(order) - i
        d = samp[samp.assay == key]["coverage"].to_numpy()
        colour = ASSAYS[key][2] if key in ASSAYS else S.GREY1
        ax.scatter(d, y + RNG.uniform(-0.18, 0.18, len(d)), s=3,
                   color=colour, linewidths=0, zorder=3)
        ax.plot([np.median(d)] * 2, [y - 0.33, y + 0.33], color=S.BLACK,
                lw=S.LW, zorder=4)
        yticks.append(y)
        ylabels.append(f"{label} ({len(d)})")

    ax.set_xscale("log")
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_formatter(mticker.NullFormatter())
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels)
    ax.set_ylim(0.4, len(order) + 0.6)
    ax.set_xlim(10, 2000)
    ax.set_xlabel("Mean sequencing coverage (×)")
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    S.save(fig, FIGS / "F_QC_seq_coverage")
    plt.close(fig)


# ------------------------------------------------------------------ main ----
def main(figures=False):
    scores = load_scores()
    print("scored variants per assay / class:")
    print(scores.pivot_table(index="type", columns="assay_label",
                             values="effect", aggfunc="size").to_string(), "\n")
    summary = (scores.groupby(["assay_label", "type"])
               .agg(n=("effect", "size"),
                    median_score=("effect", "median"),
                    median_se=("effect_se", "median"),
                    frac_lfsr05=("lfsr", lambda v: (v < 0.05).mean()))
               .round(3))
    summary.to_csv(TABS / "qc_score_se_summary.csv")
    print(summary.to_string(), "\n")

    samp = load_sample_coverage()
    samp.to_csv(TABS / "qc_sample_coverage.csv", index=False)
    cov = samp.groupby("assay")["coverage"].agg(["size", "median", "min", "max"])
    print("per-sample sequencing coverage (reads per designed variant):")
    print(cov.round(1).to_string(), "\n")

    print("wrote:")
    for q in ["tables/qc_score_se_summary.csv", "tables/qc_sample_coverage.csv"]:
        print(f"  {ROOT / q}")
    if not figures:
        return

    # Only on request: larger standalone versions of panels c-d, whose committed
    # form is the single row drawn by dms_qc_panel.py.
    score_se_figure(scores)
    seq_coverage_figure(samp)
    for q in ["figures/F_QC_score_se.pdf", "figures/F_QC_score_se.png",
              "figures/F_QC_seq_coverage.pdf", "figures/F_QC_seq_coverage.png"]:
        print(f"  {ROOT / q}")


if __name__ == "__main__":
    main(figures="--figures" in sys.argv)
