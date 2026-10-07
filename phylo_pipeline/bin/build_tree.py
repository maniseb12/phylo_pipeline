#!/usr/bin/env python3
"""
build_tree.py -- Build a phylogenetic tree (NJ or UPGMA) from an alignment
and draw a styled rectangular phylogram.

Usage:  build_tree.py <aligned.fasta> <nj|upgma> <output_prefix>
Output: <prefix>.nwk  (Newick tree)   <prefix>.png  (figure)
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import phylo_style as ps

aligned_path, method, prefix = sys.argv[1], sys.argv[2], sys.argv[3]

alignment = AlignIO.read(aligned_path, "fasta")
dm = DistanceCalculator("identity").get_distance(alignment)
builder = DistanceTreeConstructor()
if method == "nj":
    tree, label = builder.nj(dm), "Neighbor-Joining"
elif method == "upgma":
    tree, label = builder.upgma(dm), "UPGMA"
else:
    sys.exit(f"Unknown method: {method}")
tree.root_at_midpoint()

Phylo.write(tree, f"{prefix}.nwk", "newick")
ps.draw_rectangular(
    tree, f"{prefix}.png",
    title=f"{label} phylogram",
    subtitle=f"Cytochrome c  ·  {len(alignment)} species  ·  "
             f"{alignment.get_alignment_length()} aligned columns  ·  "
             f"branch length = fraction of differing sites",
    mode="phylogram", scale_label=None)
print(f"Wrote {prefix}.nwk and {prefix}.png")
