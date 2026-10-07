#!/usr/bin/env python3
"""
plot_circular.py -- Circular (radial) phylogram with a colored group ring.

Usage:  plot_circular.py <tree.nwk> <output.png>
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from Bio import Phylo
import phylo_style as ps

tree = Phylo.read(sys.argv[1], "newick")
ps.draw_circular(
    tree, sys.argv[2],
    title="Circular tree",
    subtitle="Cytochrome c  ·  radius = evolutionary distance from the root")
print(f"Wrote {sys.argv[2]}")
