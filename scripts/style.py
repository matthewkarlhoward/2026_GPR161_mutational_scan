"""Shared style + loaders for the GPR161 paper figures."""
from pathlib import Path
import json, re, io
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

SURFACE = "#ffffff"
INK = "#000000"                         # house style: all text black
LOSS = "#c0472a"                        # default outline colour in snake_axes()

# ---- two-colour diverging scale used by the landscape heatmaps -------------
# Sampled from the project swatch card. These two hexes are the poles of the
# heatmap colour map in Figure 1c; changing them changes that figure.
BLUE_GOF = "#3d6ca8"        # gain of function / positive
RED_LOF  = "#9e2b26"        # loss of function / negative

# The mechanism-class colours are NOT here: fig2_phenotype_classes.py sets them
# explicitly from the house swatch (scripts/qc/gpr161_style.py), so the figure
# that draws the classes is the one place they are defined.

# ---- house style: Helvetica, 6 pt, black text, 0.5 pt black rules ----------
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Helvetica Neue", "Arial", "DejaVu Sans"],
    "font.size": 6,
    "pdf.fonttype": 42, "ps.fonttype": 42,          # keep text as text in vector output
    # transparent on save: these land in Illustrator over other artwork, so an
    # opaque white rectangle behind every panel is a nuisance. On-screen the
    # figure still renders light; only the written file is clear.
    "figure.facecolor": SURFACE, "axes.facecolor": "none",
    "savefig.facecolor": "none", "savefig.transparent": True,
    "text.color": "#000000", "axes.labelcolor": "#000000",
    "axes.edgecolor": "#000000", "axes.titlecolor": "#000000",
    "xtick.color": "#000000", "ytick.color": "#000000",
    "xtick.labelcolor": "#000000", "ytick.labelcolor": "#000000",
    "axes.linewidth": 0.5, "lines.linewidth": 0.5, "patch.linewidth": 0.5,
    "grid.linewidth": 0.5, "hatch.linewidth": 0.5,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.minor.width": 0.5, "ytick.minor.width": 0.5,
    "xtick.major.size": 2.0, "ytick.major.size": 2.0,
    "xtick.minor.size": 1.2, "ytick.minor.size": 1.2,
    "axes.labelsize": 6, "axes.titlesize": 6,
    "xtick.labelsize": 6, "ytick.labelsize": 6,
    "legend.frameon": False, "legend.fontsize": 6,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 150, "savefig.dpi": 300,
})

# ---- millimetre sizing: every figure is laid out in mm, 165 x 180 mm max ---
# The cap is enforced where it matters, in save_figure() at the bottom.
MM = 1 / 25.4

def load_dms():
    """The variant table (pipeline/01_merge_scores.py)."""
    return pd.read_csv(ROOT / "data/processed/merged_variant_scores.tsv", sep="\t")

def load_positions():
    """The position table (pipeline/02_position_summary.py)."""
    return pd.read_csv(ROOT / "data/processed/position_summary.tsv", sep="\t")


# ---- GPCRdb snake plot as a matplotlib panel -------------------------------
_SNAKE_SVG = ROOT / "data/annotation/snakeplot_raw.svg"
_SNAKE_CACHE = {}

def snake_layout():
    """Parse GPCRdb's published snake-plot coordinates for GPR161.

    Returns {position: (x, y, aa, segment)} in matplotlib orientation
    (SVG y runs downward, so it is negated here).
    """
    if _SNAKE_CACHE:
        return _SNAKE_CACHE
    import re
    svg = io.open(_SNAKE_SVG, encoding="utf-8", errors="replace").read()
    m = re.search(r"translate\(0,\s*([\d.]+)\)", svg)
    dy = float(m.group(1)) if m else 0.0
    pat = re.compile(
        r"<circle[^>]*?class='([^']*)'[^>]*?cx='([\d.-]+)'[^>]*?cy='([\d.-]+)'"
        r"[^>]*?id='(\d+)'[^>]*?title='([^']*)'", re.S)
    for m2 in pat.finditer(svg):
        cls, cx, cy, pid, title = m2.groups()
        seg = cls.replace("rcircle", "").replace("long", "").replace("short", "").strip()
        _SNAKE_CACHE[int(pid)] = (float(cx), -(float(cy) + dy), title.split()[0][0], seg)
    return _SNAKE_CACHE

def snake_axes(ax, values=None, cmap="RdBu_r", vmin=-4, vmax=4, nan="#e4e2dc",
               categorical=None, catmap=None, outline=(), outline_color=None,
               r=11.0, backbone="#000000", lw=0.5, letters=False, letter_fs=3.4):
    """Draw the GPCRdb snake plot into `ax`, coloured by `values` (or classes)."""
    import matplotlib as mpl
    from matplotlib.collections import EllipseCollection
    L = snake_layout()
    order = sorted(L)
    xy = np.array([[L[p][0], L[p][1]] for p in order])
    ax.plot(xy[:, 0], xy[:, 1], color=backbone, lw=lw, zorder=1, solid_capstyle="round")
    norm = mpl.colors.Normalize(vmin, vmax)
    # accept either a registered colormap name or a Colormap object, so callers
    # can pass a palette built from BLUE_GOF / RED_LOF
    cm = cmap if isinstance(cmap, mpl.colors.Colormap) else mpl.colormaps[cmap]
    face = []
    for p in order:
        if categorical is not None:
            face.append(catmap.get(categorical.get(p), "#f4f3ee"))
        else:
            v = np.nan if values is None else values.get(p, np.nan)
            face.append(nan if (v is None or (isinstance(v, float) and np.isnan(v)))
                        else cm(norm(v)))
    ec = EllipseCollection(np.full(len(order), 2 * r), np.full(len(order), 2 * r),
                           np.zeros(len(order)), units="xy", offsets=xy,
                           offset_transform=ax.transData, facecolors=face,
                           edgecolors="#000000", linewidths=lw, zorder=3)
    ax.add_collection(ec)
    if outline:
        oc = outline_color or LOSS
        sel = [p for p in order if p in set(outline)]
        if sel:
            o = np.array([[L[p][0], L[p][1]] for p in sel])
            ax.add_collection(EllipseCollection(
                np.full(len(sel), 2 * r * 1.05), np.full(len(sel), 2 * r * 1.05),
                np.zeros(len(sel)), units="xy", offsets=o,
                offset_transform=ax.transData, facecolors="none",
                edgecolors=oc, linewidths=lw * 3, zorder=4))
    if letters:
        for p in order:
            x, y, aa, _ = L[p]
            ax.text(x, y, aa, ha="center", va="center", fontsize=letter_fs, zorder=5)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_xlim(xy[:, 0].min() - 3 * r, xy[:, 0].max() + 3 * r)
    ax.set_ylim(xy[:, 1].min() - 3 * r, xy[:, 1].max() + 3 * r)
    return mpl.cm.ScalarMappable(norm, cm)


# ---- multi-stage binary-gate Sankey as a matplotlib panel --------------------
def sankey_axes(ax, paths, stages, colmap, label_classes=True, min_label=10,
                row_labels=("spared", "affected"), node_w=0.028, gap=0.05,
                top=0.88, bot=0.07, x0=0.10, x1=None, alpha=0.85, counts=True,
                label_gap=0.052, stage_labels=True, label_fmt=None):
    """Draw a Sankey of items flowing through binary gates.

    `label_fmt(path)` may return (text, colour) for the right-hand class labels;
    the default is black "<cls>  <n>".

    `paths` is a list of {'flags': (0/1 per stage), 'n': int, 'cls': key of
    `colmap`}; every combination of flags is one ribbon, so the paths ARE the
    classes. Ribbons are stacked in one fixed order (sorted by flags) at every
    node, which leaves only the genuine crossings: paths that change side.
    Returns {id(path): y} for the terminal centre of each ribbon.
    """
    from matplotlib.path import Path as MPath
    from matplotlib.patches import PathPatch, Rectangle
    ns = len(stages)
    tot = sum(p["n"] for p in paths)
    x1 = (0.72 if label_classes else 0.95) if x1 is None else x1
    X = np.linspace(x0, x1, ns)
    scale = (top - bot - gap) / tot
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    node = {}
    for i in range(ns):
        n_aff = sum(p["n"] for p in paths if p["flags"][i])
        y = top
        node[(i, 0)] = (y - (tot - n_aff) * scale, y)          # spared, on top
        y = y - (tot - n_aff) * scale - gap
        node[(i, 1)] = (y - n_aff * scale, y)                  # affected, below
        for flag, cnt in ((0, tot - n_aff), (1, n_aff)):
            a, b = node[(i, flag)]
            ax.add_patch(Rectangle((X[i] - node_w / 2, a), node_w, b - a,
                                   fc="white", ec="#000000", lw=0.5, zorder=3))
            if counts and (b - a) >= 0.085:
                ax.text(X[i], (a + b) / 2, str(cnt), ha="center", va="center",
                        rotation=90, zorder=4)
            elif counts:        # node too short for a rotated label inside it
                ax.text(X[i], a - 0.018 if flag else b + 0.018, str(cnt),
                        ha="center", va="top" if flag else "bottom", zorder=4)
        if stage_labels:
            ax.text(X[i], top + 0.035, stages[i], ha="center", va="bottom")
    if row_labels:
        for flag, nm in enumerate(row_labels):
            a, b = node[(0, flag)]
            ax.text(X[0] - node_w / 2 - 0.013, (a + b) / 2, nm, ha="right", va="center")

    off = {k: 0.0 for k in node}
    ends = {}
    for p in sorted(paths, key=lambda q: q["flags"]):
        if p["n"] == 0:
            continue
        h = p["n"] * scale
        z = 2 if p["n"] > 0.1 * tot else 2.5      # thin ribbons ride on top
        for i in range(ns - 1):
            src, dst = (i, p["flags"][i]), (i + 1, p["flags"][i + 1])
            ys, yd = node[src][1] - off[src], node[dst][1] - off[dst]
            xa, xb = X[i] + node_w / 2, X[i + 1] - node_w / 2
            cx = (xa + xb) / 2
            verts = [(xa, ys), (cx, ys), (cx, yd), (xb, yd), (xb, yd - h),
                     (cx, yd - h), (cx, ys - h), (xa, ys - h), (xa, ys)]
            codes = [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4,
                     MPath.LINETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4,
                     MPath.CLOSEPOLY]
            ax.add_patch(PathPatch(MPath(verts, codes), fc=colmap[p["cls"]],
                                   ec="none", alpha=alpha, zorder=z))
            off[src] += h
            if i == ns - 2:
                ends[id(p)] = yd - h / 2
                off[dst] += h
    if label_classes:
        prev = 1e9
        for p in sorted([q for q in paths if q["n"] >= min_label],
                        key=lambda q: -ends[id(q)]):
            y = min(ends[id(p)], prev - label_gap)
            txt, tcol = (label_fmt(p) if label_fmt
                         else (f"{p['cls']}  {p['n']}", "#000000"))
            ax.text(X[-1] + node_w / 2 + 0.012, y, txt, va="center", color=tcol)
            prev = y
    return ends


# ---- the position-level gate rule, in ONE place -----------------------------
# Which rule decides "this assay is affected at this position" is a real degree
# of freedom; docs/METHODS.md section 5 states the alternative that was rejected
# and why. Changing POSITION_GATE here changes every figure that calls
# position_gates().
#
#   "mean"  position mean z < CUT, no variant count          <- project default
#   "count" >= MINVAR variants with z < CUT and lfsr < 0.05
#
# Chosen 2026-10-22. The count rule's flat threshold is not a flat bar: the three
# assays detect at 5.6 % (abundance) vs 15.5 % (surface) and 14.1 % (cAMP), so at
# ~19 variants per position >= 3 is 2.8x the expected count for abundance and
# 1.0x for surface, which made the abundance-dependent split (destabilized vs not
# delivered) a function of the threshold. The mean rule also puts Y208 and Y324
# in "uncoupled", agreeing with the Cosmos signalling-specific calls.
POSITION_GATE = "mean"
CUT, MINVAR = -2.0, 3
_GCOL = {"A": ("abu", "abundance_z", "abundance_lfsr"),
         "S": ("surf", "surface_z", "surface_lfsr"),
         "C": ("cre", "cre_z", "cre_lfsr")}


def position_gates(posn=None, variants=None, rule=None, minvar=MINVAR, cut=CUT):
    """{"A","S","C"} -> set of positions whose gate is affected."""
    rule = rule or POSITION_GATE
    if rule == "mean":
        if posn is None:
            raise ValueError("rule='mean' needs the position table")
        return {k: set(posn.loc[posn[col] < cut, "position"]) for k, (col, _z, _l) in _GCOL.items()}
    if rule == "count":
        if variants is None:
            raise ValueError("rule='count' needs the variant table")
        out = {}
        for k, (_col, z, lf) in _GCOL.items():
            n = variants[(variants[z] < cut) & (variants[lf] < 0.05)].groupby("position").size()
            out[k] = set(n[n >= minvar].index)
        return out
    raise ValueError(f"unknown rule {rule!r}")


def position_class(a, s, c):
    """The class is named for the most upstream AFFECTED gate.

    At the surface step the downstream cAMP gate is consulted as well, to split
    "not delivered" from "surface-loss only". Nothing else looks downstream.
    """
    if a:
        return "destabilized"
    if s:
        return "not delivered" if c else "surface-loss only"
    return "uncoupled" if c else "WT-like"


def gate_paths(gates, universe):
    """The eight paths through the three gates, shaped for sankey_axes()."""
    out = []
    for a in (0, 1):
        for s in (0, 1):
            for c in (0, 1):
                sel = set(universe)
                for flag, k in ((a, "A"), (s, "S"), (c, "C")):
                    sel = sel & gates[k] if flag else sel - gates[k]
                out.append(dict(flags=(a, s, c), n=len(sel), cls=position_class(a, s, c)))
    return out


def gate_label(rule=None, minvar=MINVAR, cut=CUT):
    """One-line description of the rule in force, for captions."""
    rule = rule or POSITION_GATE
    return ("a gate counts as affected when the position mean z < −2"
            if rule == "mean" else
            f"a gate counts as affected at ≥{minvar} variants with z < −2 and lfsr < 0.05")


# ---- save with the size cap ENFORCED, not assumed --------------------------
def save_figure(fig, name, max_w_mm=165.0, max_h_mm=180.0, dpi=600, quiet=False):
    """Write PNG + PDF at exactly the requested figsize, and verify the result.

    Also reports text that overlaps other text, and text that falls off the page.

    Deliberately does NOT pass bbox_inches="tight". That expands the canvas past
    the requested size whenever a label overflows, which is how a figure ends up
    wider than the cap without anyone noticing. The trade-off is that overflowing
    text is CLIPPED instead of growing the page, so this also reports anything
    landing outside the canvas rather than letting it disappear quietly.
    """
    w_mm, h_mm = (v * 25.4 for v in fig.get_size_inches())
    if w_mm > max_w_mm + 1e-6 or h_mm > max_h_mm + 1e-6:
        raise ValueError(f"{name}: {w_mm:.1f} x {h_mm:.1f} mm exceeds "
                         f"{max_w_mm:.0f} x {max_h_mm:.0f} mm")
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    PW, PH = fig.bbox.width, fig.bbox.height
    items = list(fig.texts)
    for ax in fig.axes:
        items += list(ax.texts) + [ax.title, ax.xaxis.label, ax.yaxis.label]
        # tick labels too: a long categorical y-tick is the classic way text
        # walks off the left edge without any of the above noticing. Skip an axis
        # that is hidden -- a twinx's x-axis still reports labels it never draws,
        # which otherwise shows up as a phantom overlap with the real ones.
        # ax.axis("off") sets ax.axison, NOT the individual axes' visibility, so
        # both have to be checked or a hidden axes' default ticks (0.0 ... 1.0)
        # get reported as clipped text that is never drawn.
        if getattr(ax, "axison", True):
            # A locator can create a tick label for a value outside the view
            # limits; it is never drawn but still reports a position, which then
            # shows up as a phantom collision. Keep only labels that sit within
            # their own axes' box (generously padded for the label's own size).
            abb = ax.get_window_extent()
            cand = []
            if ax.xaxis.get_visible():
                cand += list(ax.get_xticklabels())
            if ax.yaxis.get_visible():
                cand += list(ax.get_yticklabels())
            for t in cand:
                try:
                    tb = t.get_window_extent(rend)
                except Exception:
                    continue
                pad = max(tb.width, tb.height) + 2
                if (tb.x1 > abb.x0 - pad and tb.x0 < abb.x1 + pad
                        and tb.y1 > abb.y0 - pad and tb.y0 < abb.y1 + pad):
                    items.append(t)
        if ax.get_legend() is not None:
            items.append(ax.get_legend())
    off = []
    for it in items:
        if it is None or not it.get_visible() or not str(getattr(it, "get_text", lambda: "x")()):
            continue
        try:
            bb = it.get_window_extent(rend)
        except Exception:
            continue
        if bb.x0 < -1 or bb.x1 > PW + 1 or bb.y0 < -1 or bb.y1 > PH + 1:
            lbl = getattr(it, "get_text", lambda: "<legend>")()
            off.append(str(lbl)[:46].replace("\n", " "))
    # text-on-text collisions: cheap to check, and the thing most often missed
    boxes = []
    for it in items:
        if it is None or not it.get_visible():
            continue
        lbl = str(getattr(it, "get_text", lambda: "")()).strip()
        if not lbl:
            continue
        try:
            boxes.append((lbl[:38].replace("\n", " "), it.get_window_extent(rend)))
        except Exception:
            pass
    clash = []
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            a, b = boxes[i][1], boxes[j][1]
            ox = min(a.x1, b.x1) - max(a.x0, b.x0)
            oy = min(a.y1, b.y1) - max(a.y0, b.y0)
            # a pixel or two of adjacency between neighbouring tick labels is
            # not a collision; flag only overlaps big enough to read as one
            if ox > 3 and oy > 3:
                clash.append((boxes[i][0], boxes[j][0]))
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"{name}.{ext}", dpi=dpi)
    plt.close(fig)
    # Ink-based margin check on the SAVED file. The text scan above cannot see
    # drawn collections (a snake plot, ribbons, points), so anything non-text that
    # runs to the edge only shows up here.
    margins = None
    try:
        from PIL import Image
        import numpy as _np
        _a = _np.array(Image.open(FIG / f"{name}.png").convert("L"))
        _ink = _a < 250
        _c = _np.where(_ink.any(axis=0))[0]; _r = _np.where(_ink.any(axis=1))[0]
        if len(_c) and len(_r):
            _px2mm = 25.4 / dpi
            margins = dict(left=_c[0] * _px2mm, right=(_a.shape[1] - 1 - _c[-1]) * _px2mm,
                           top=_r[0] * _px2mm, bottom=(_a.shape[0] - 1 - _r[-1]) * _px2mm)
    except Exception:
        pass
    if not quiet:
        print(f"wrote {name}  ({w_mm:.0f} x {h_mm:.0f} mm)")
        if off:
            print(f"  WARNING: {len(off)} item(s) fall outside the canvas and are clipped:")
            for o in dict.fromkeys(off):
                print(f"    {o!r}")
        if clash:
            print(f"  WARNING: {len(clash)} overlapping text pair(s):")
            for a, b in clash[:8]:
                print(f"    {a!r} <-> {b!r}")
        if margins:
            tight = {k: v for k, v in margins.items() if v < 2.0}
            if tight:
                print("  WARNING: ink within 2 mm of the page edge: "
                      + ", ".join(f"{k} {v:.1f} mm" for k, v in tight.items()))
    return off, clash
