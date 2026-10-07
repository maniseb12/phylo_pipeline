#!/usr/bin/env python3
"""
make_gallery.py -- Tile all pipeline figures into one overview image.

Looks for the five PNGs in the current directory (Nextflow stages them
here) and writes gallery.png.
Usage:  make_gallery.py [output.png]
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import phylo_style as ps

out = sys.argv[1] if len(sys.argv) > 1 else "gallery.png"
panels = [
    ("tree_nj.png",     "1 · Neighbor-Joining phylogram"),
    ("tree_upgma.png",  "2 · UPGMA phylogram"),
    ("cladogram.png",   "3 · Cladogram (topology only)"),
    ("circular.png",    "4 · Circular tree"),
    ("bootstrap.png",   "5 · Bootstrap support"),
]

fig = plt.figure(figsize=(21, 12.5), facecolor=ps.BG)
fig.text(0.5, 0.965, "Phylogenetic tree pipeline  ·  Cytochrome c across 9 species",
         ha="center", fontsize=26, fontweight="bold", color=ps.INK)
fig.text(0.5, 0.932, "MAFFT alignment  →  Neighbor-Joining / UPGMA  →  "
         "5 views of the same evolutionary story",
         ha="center", fontsize=14, color="#6b7280")

gs = fig.add_gridspec(2, 3, left=0.02, right=0.98, bottom=0.03, top=0.90,
                      wspace=0.04, hspace=0.12)
for i, (fname, label) in enumerate(panels):
    ax = fig.add_subplot(gs[i // 3, i % 3])
    ax.imshow(mpimg.imread(fname))
    ax.axis("off")

# 6th tile: key / reading guide
ax = fig.add_subplot(gs[1, 2])
ax.axis("off")
ax.text(0.04, 0.97, "How to read these figures", fontsize=17,
        fontweight="bold", color=ps.INK, transform=ax.transAxes, va="top")
notes = [
    ("Branch color", "clade / group of the species below"),
    ("Grey branch", "mixes several groups"),
    ("Phylogram", "length = amount of change"),
    ("Cladogram", "shape only, lengths ignored"),
    ("Support dots", "% of bootstrap replicates agreeing"),
]
y = 0.86
for k, v in notes:
    ax.text(0.04, y, k, fontsize=13, fontweight="bold", color=ps.INK,
            transform=ax.transAxes, va="top")
    ax.text(0.04, y - 0.05, v, fontsize=12, color="#6b7280",
            transform=ax.transAxes, va="top")
    y -= 0.105
y -= 0.045
for g, c in ps.GROUP_COLORS.items():
    if g == "Other":
        continue
    ax.scatter([0.06], [y], s=120, color=c, transform=ax.transAxes)
    ax.text(0.11, y, g, fontsize=12, color=ps.INK, transform=ax.transAxes,
            va="center")
    y -= 0.05

fig.savefig(out, dpi=110, facecolor=ps.BG)
print(f"Wrote {out}")
