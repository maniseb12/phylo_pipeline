#!/usr/bin/env python3
"""
bootstrap_tree.py -- NJ tree annotated with bootstrap support.

Alignment columns are resampled N times and an NJ tree is built from each
replicate. For every split in the reference NJ tree we report the % of
replicates that contain the same split (colored dots on the branches).

Usage:  bootstrap_tree.py <aligned.fasta> <n_replicates> <output_prefix>
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo.Consensus import bootstrap_trees, get_support
import phylo_style as ps

aligned_path, n_rep, prefix = sys.argv[1], int(sys.argv[2]), sys.argv[3]
alignment = AlignIO.read(aligned_path, "fasta")
calc = DistanceCalculator("identity")

reference = DistanceTreeConstructor().nj(calc.get_distance(alignment))
reference.root_at_midpoint()

replicates = list(bootstrap_trees(alignment, n_rep,
                                  DistanceTreeConstructor(calc, "nj")))
tree = get_support(reference, replicates)

Phylo.write(tree, f"{prefix}.nwk", "newick")
ps.draw_rectangular(
    tree, f"{prefix}.png",
    title="Bootstrap-supported tree",
    subtitle=f"Cytochrome c  ·  NJ tree  ·  {n_rep} bootstrap replicates  ·  "
             f"dots show % of replicates supporting each split",
    mode="phylogram", show_support=True)
print(f"Wrote {prefix}.nwk and {prefix}.png")
