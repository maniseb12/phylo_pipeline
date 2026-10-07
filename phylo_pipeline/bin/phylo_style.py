"""
phylo_style.py -- shared colors, layout maths and drawing helpers.

Every plotting script imports this so all figures look consistent.
To recolor for your own data, edit GROUPS below (species -> group).
Any species not listed falls into "Other" (grey).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ---- Edit me: map each sequence ID (first word of the FASTA header) to a group
GROUPS = {
    "Human": "Primates",
    "Chimpanzee": "Primates",
    "Rhesus_Monkey": "Primates",
    "Horse": "Other mammals",
    "Cow": "Other mammals",
    "Chicken": "Birds",
    "Rattlesnake": "Reptiles",
    "Tuna": "Fish",
    "Yeast_iso1": "Fungi",
}

GROUP_COLORS = {
    "Primates":      "#e63946",
    "Other mammals": "#f4a261",
    "Birds":         "#2a9d8f",
    "Reptiles":      "#7fb069",
    "Fish":          "#3a86ff",
    "Fungi":         "#9d4edd",
    "Other":         "#6c757d",
}
GROUP_ORDER = list(GROUP_COLORS)

BG = "#fbfaf7"      # figure background (warm off-white)
INK = "#2b2d42"     # text color
MIXED = "#8d99ae"   # branches spanning several groups
FONT = "DejaVu Sans"


def group_of(name):
    return GROUPS.get(name, "Other")


def clean(name):
    return (name or "").replace("_", " ")


def clade_color(clade):
    groups = {group_of(t.name) for t in clade.get_terminals()}
    return GROUP_COLORS[groups.pop()] if len(groups) == 1 else MIXED


def present_groups(tree):
    found = {group_of(t.name) for t in tree.get_terminals()}
    return [g for g in GROUP_ORDER if g in found]


def prepare(tree, mode="phylogram"):
    """Return (leaves, x, y). y = row per leaf; x = depth (phylogram) or
    negative step-height so tips line up on the right (cladogram)."""
    tree.ladderize()
    root = tree.root
    leaves = tree.get_terminals()
    y = {leaf: float(i) for i, leaf in enumerate(leaves)}

    def ypos(c):
        if c.is_terminal():
            return y[c]
        ys = [ypos(k) for k in c.clades]
        y[c] = (min(ys) + max(ys)) / 2
        return y[c]
    ypos(root)

    x = {}
    if mode == "phylogram":
        def xpos(c, base):
            x[c] = base
            for k in c.clades:
                xpos(k, base + (k.branch_length or 0.0))
        xpos(root, 0.0)
    else:
        h = {}
        def height(c):
            h[c] = 0 if c.is_terminal() else 1 + max(height(k) for k in c.clades)
            return h[c]
        height(root)
        x = {c: -v for c, v in h.items()}
    return leaves, x, y


def group_legend(fig, tree, loc="lower left", anchor=(0.04, 0.03), extra=None, ncol=1):
    handles = [Line2D([0], [0], marker="o", linestyle="", markersize=9,
                      markerfacecolor=GROUP_COLORS[g], markeredgecolor="white",
                      label=g) for g in present_groups(tree)]
    if extra:
        handles += extra
    fig.legend(handles=handles, loc=loc, bbox_to_anchor=anchor, frameon=False,
               fontsize=10, labelcolor=INK, handletextpad=0.4, labelspacing=0.7,
               ncol=ncol)


def support_color(v):
    if v >= 90:
        return "#2d6a4f"
    if v >= 70:
        return "#e9a23b"
    return "#9aa0a6"


def nice_scale(span):
    for s in (0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5):
        if s >= span * 0.12:
            return s
    return 0.5


def draw_rectangular(tree, out_png, title, subtitle, mode="phylogram",
                     show_support=False, scale_label=None):
    leaves, x, y = prepare(tree, mode)
    n = len(leaves)
    xs = list(x.values())
    xmin, xmax = min(xs), max(xs)
    span = (xmax - xmin) or 1.0

    fig = plt.figure(figsize=(10.5, 6.8), facecolor=BG)
    ax = fig.add_axes([0.05, 0.13, 0.83, 0.70], facecolor=BG)
    ax.set_axis_off()

    # branches
    for clade in tree.find_clades(order="level"):
        col = clade_color(clade)
        if clade.clades:
            ys = [y[k] for k in clade.clades]
            ax.plot([x[clade], x[clade]], [min(ys), max(ys)], color=col,
                    lw=2.4, solid_capstyle="round", zorder=2)
        for k in clade.clades:
            ax.plot([x[clade], x[k]], [y[k], y[k]], color=clade_color(k),
                    lw=2.4, solid_capstyle="round", zorder=2)
    # root stub
    stub = span * 0.035
    ax.plot([x[tree.root] - stub, x[tree.root]], [y[tree.root]] * 2,
            color=MIXED, lw=2.4, solid_capstyle="round", zorder=1)

    # tips
    pad = span * 0.018
    for leaf in leaves:
        c = GROUP_COLORS[group_of(leaf.name)]
        ax.scatter([x[leaf]], [y[leaf]], s=85, color=c, edgecolor="white",
                   linewidth=1.5, zorder=4)
        ax.text(x[leaf] + pad * 2.2, y[leaf], clean(leaf.name), va="center",
                ha="left", fontsize=12, color=INK, fontweight="medium")

    # bootstrap support dots (sit on the branch leading to each node)
    if show_support:
        parent = {k: c for c in tree.find_clades() for k in c.clades}
        for c in tree.get_nonterminals():
            if c is tree.root or c.confidence is None:
                continue
            v = float(c.confidence)
            xd = (x[parent[c]] + x[c]) / 2
            ax.scatter([xd], [y[c]], s=430, color=support_color(v),
                       edgecolor="white", linewidth=1.6, zorder=5)
            ax.text(xd, y[c], f"{v:.0f}", color="white", fontsize=8,
                    fontweight="bold", ha="center", va="center", zorder=6)

    # scale bar (phylogram only)
    if mode == "phylogram":
        s = nice_scale(span)
        yb = n - 0.1
        ax.plot([xmin, xmin + s], [yb, yb], color=INK, lw=2.2,
                solid_capstyle="butt")
        ax.text(xmin + s / 2, yb + 0.55, scale_label or f"{s:g}", ha="center",
                va="top", fontsize=9.5, color=INK)

    ax.set_xlim(xmin - span * 0.05, xmax + span * 0.40)
    ax.set_ylim(n + 0.2, -0.9)

    fig.text(0.05, 0.93, title, fontsize=19, fontweight="bold", color=INK,
             ha="left", va="center")
    fig.text(0.05, 0.885, subtitle, fontsize=10.5, color="#6b7280",
             ha="left", va="center")

    extra = None
    if show_support:
        extra = [Line2D([0], [0], marker="o", linestyle="", markersize=9,
                        markerfacecolor=support_color(v), markeredgecolor="white",
                        label=l) for v, l in ((95, "support ≥ 90%"),
                                              (80, "70–89%"), (50, "< 70%"))]
    group_legend(fig, tree, loc="lower right", anchor=(0.97, 0.03), extra=extra)
    fig.savefig(out_png, dpi=200, facecolor=BG)
    plt.close(fig)


def draw_circular(tree, out_png, title, subtitle):
    leaves, x, y = prepare(tree, "phylogram")
    n = len(leaves)
    rmax = max(x[l] for l in leaves) or 1.0
    r = {c: x[c] / rmax for c in x}
    slot = 2 * np.pi / (n + 1)
    th = {c: y[c] * slot for c in y}

    fig = plt.figure(figsize=(10, 10), facecolor=BG)
    ax = fig.add_axes([0.2, 0.15, 0.6, 0.6], projection="polar", facecolor=BG)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_axis_off()

    for clade in tree.find_clades(order="level"):
        if clade.clades:
            ts = [th[k] for k in clade.clades]
            arc = np.linspace(min(ts), max(ts), 80)
            ax.plot(arc, [r[clade]] * len(arc), color=clade_color(clade),
                    lw=2.4, solid_capstyle="round")
        for k in clade.clades:
            ax.plot([th[k], th[k]], [r[clade], r[k]], color=clade_color(k),
                    lw=2.4, solid_capstyle="round")

    for leaf in leaves:
        c = GROUP_COLORS[group_of(leaf.name)]
        t = th[leaf]
        ax.plot([t, t], [r[leaf], 1.04], color="#cfd3d8", lw=1.1,
                ls=(0, (2, 3)), zorder=0)
        ax.scatter([t], [r[leaf]], s=80, color=c, edgecolor="white",
                   linewidth=1.5, zorder=4)
        ax.bar(t, 0.055, bottom=1.06, width=slot * 0.86, color=c, alpha=0.95,
               linewidth=0)
        deg = np.degrees(t)
        left = 180 < deg < 360
        rot = 90 - deg + (180 if left else 0)
        ax.text(t, 1.16, clean(leaf.name), rotation=rot, rotation_mode="anchor",
                ha="right" if left else "left", va="center", fontsize=12,
                color=INK, fontweight="medium")

    ax.set_ylim(0, 1.12)
    fig.text(0.5, 0.95, title, fontsize=20, fontweight="bold", color=INK, ha="center")
    fig.text(0.5, 0.917, subtitle, fontsize=10.5, color="#6b7280", ha="center")
    group_legend(fig, tree, loc="lower center", anchor=(0.5, 0.02),
                 ncol=len(present_groups(tree)))
    fig.savefig(out_png, dpi=200, facecolor=BG)
    plt.close(fig)
