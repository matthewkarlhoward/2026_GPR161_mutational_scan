#!/usr/bin/env python3
"""Figure 4 panel b -- mean missense surface effect by receptor segment.

    data/processed/{merged_variant_scores,position_summary}.tsv
    data/annotation/structure_8smv.tsv
        -> figures/segment_effects.pdf/.png          50 x 40 mm

The MEAN MISSENSE surface effect per GPCRdb segment, N to C along the sequence,
in raw DMS units. Signed, so losses are negative and a gain at one position
cancels a loss at another; ICL1, ICL2 and H8 come out slightly positive.

Figure 4 panels a and c are ChimeraX renders made by hand from
data/structure/8smv.cif, coloured by the per-position surface and cAMP values in
data/processed/position_summary.tsv. They are not drawn here and the renders are
not deposited. Panel letters are added in Illustrator.

Run:  python3 scripts/fig3_structure.py
"""
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from style import *                                                  # noqa: F403


# house swatch, as in fig2_phenotype_classes.py
BLUE, GREEN = "#4877B5", "#316B3E"
GREY1 = "#808180"

ASSAYS = [("abundance", "Abundance", GREEN),
          ("surface", "Surface", BLUE),
          ("cre", "cAMP", GREY1)]
DROP_SEGMENTS = {"C-term"}      # one modelled position; a bar would be a point


# --------------------------------------------------------------------- data
def segment_table():
    """Mean missense effect per segment, modelled positions only.

    ✔ 2026-10-28: the plain mean of the missense effects measured at a position,
    in raw DMS effect units, so losses are negative and a gain at one position
    cancels a loss at another. The structure panels were switched to the same
    quantity on the same day, so panels a and c are once again the same number
    at two grains.
    """
    d = load_dms()
    m = d[d.type == "missense"]
    sa = pd.read_csv(ROOT / "data/annotation/structure_8smv.tsv", sep="\t")
    per = pd.DataFrame({
        a: m.groupby("position")[f"{a}_eff"].mean()
        for a, _, _ in ASSAYS}).reset_index()
    t = sa[sa.modeled].merge(per, on="position")
    g = t.groupby("segment")
    out = g.agg(n=("position", "size"),
                **{f"{a}_m": (a, "mean") for a, _, _ in ASSAYS},
                **{f"{a}_se": (a, "sem") for a, _, _ in ASSAYS})
    out = out.drop(index=[s for s in DROP_SEGMENTS if s in out.index])
    # ✔ 2026-10-28: ordered along the sequence (N to C) rather than by effect
    # size, so the panel reads as a walk through the receptor and can be laid
    # beside the structure panels. Ranking by surface effect is still in the
    # printed table that main() emits.
    out["first_position"] = t.groupby("segment").position.min()
    return out.sort_values("first_position"), t


def panel_segments(fig, W, H, mm_axes, seg):
    """Mean missense SURFACE effect per segment, N to C.

    ✔ 2026-10-28: surface only, black bars, no error bars, 50 x 40 mm. Losses are
    negative, so the bars hang from zero. The three-assay version is in the
    repository history; cAMP and abundance by segment are still printed by main()
    and belong in a supplement if they are wanted.
    """
    ax = mm_axes(fig, 15.0, 3.5, W - 17.5, H - 15)
    x = np.arange(len(seg))
    ax.bar(x, seg["surface_m"], 0.72, color=INK, edgecolor=INK, linewidth=0.4,
           zorder=2)
    # ✔ 2026-10-28: no receptor-wide reference line. The zero line stays, since
    # the quantity is signed and the bars hang from it. main() prints the
    # receptor-wide means for the text.
    ax.axhline(0, color=INK, lw=0.5, zorder=1)
    ax.set_xticks(x)
    ax.set_xticklabels(seg.index, rotation=90, fontsize=6, color=INK)
    ax.tick_params(axis="x", length=0, pad=1.5, labelsize=6, colors=INK)
    ax.tick_params(axis="y", labelsize=6, colors=INK)
    ax.set_xlim(-0.65, len(seg) - 0.35)
    # ICL1, ICL2 and H8 have positive means -- the signed quantity puts them
    # above zero, and clipping the axis at 0 would render them as no bar at all.
    lo = float(seg["surface_m"].min()) * 1.12
    hi = max(float(seg["surface_m"].max()) * 1.35, 0.0)
    ax.set_ylim(lo, hi)
    # Explicit ticks: the automatic locator keeps a label just past the top of
    # the view, which the style check reports as ink off the canvas.
    ax.set_yticks(np.arange(np.ceil(lo / 0.2) * 0.2, hi + 1e-9, 0.2))
    ax.set_ylabel("Mean missense\nsurface effect", fontsize=6, color=INK)


# --------------------------------------------------------------------- main
def mm_axes_for(W, H):
    def f(fig, x, top, w, h):
        return fig.add_axes([x / W, 1 - (top + h) / H, w / W, h / H])
    return f


def build(name, W, H, draw):
    fig = plt.figure(figsize=(W * MM, H * MM))
    draw(fig, W, H, mm_axes_for(W, H))
    save_figure(fig, name)
    print(f"  {name:24s} {W:.0f} x {H:.0f} mm")


def main():
    seg, per = segment_table()

    print("segment means (modelled positions, mean missense effect)")
    print(seg.round(3).to_string())
    print(f"all modelled: surface {per['surface'].mean():.3f}  "
          f"cAMP {per['cre'].mean():.3f}  abundance {per['abundance'].mean():.3f}")

    # Figure 4 panels a and c are ChimeraX renders made by hand from
    # data/structure/8smv.cif coloured by the per-position values in
    # data/processed/position_summary.tsv. The repository therefore
    # contributes only panel b.
    BW, BH = 50.0, 40.0
    build("segment_effects", BW, BH,
          lambda fig, W, H, mm: panel_segments(fig, W, H, mm, seg))


if __name__ == "__main__":
    sys.exit(main())
