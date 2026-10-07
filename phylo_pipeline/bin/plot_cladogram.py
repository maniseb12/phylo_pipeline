#!/usr/bin/env python3
"""
plot_cladogram.py -- Topology-only view: every split is one step, tips
aligned on the right, so you read *relationships* without distance clutter.

Usage:  plot_cladogram.py <tree.nwk> <output.png>
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from Bio import Phylo
import phylo_style as ps

tree = Phylo.read(sys.argv[1], "newick")
ps.draw_rectangular(
    tree, sys.argv[2],
    title="Cladogram",
    subtitle="Cytochrome c  ·  topology only  ·  branch lengths ignored, "
             "tips aligned",
    mode="cladogram")
print(f"Wrote {sys.argv[2]}")
