#!/usr/bin/env nextflow
/*
 * Phylogenetic Tree Pipeline
 * ---------------------------
 * Input:  unaligned protein/gene sequences (FASTA)
 * Output: 5 forms of phylogenetic graph:
 *           1. Neighbor-Joining phylogram (rectangular, distance-scaled)
 *           2. UPGMA phylogram          (rectangular, distance-scaled)
 *           3. Cladogram                (topology only, no distance scaling)
 *           4. Circular / radial tree
 *           5. Bootstrap-supported tree (branch support dots)
 *         plus one combined overview image (gallery.png)
 *
 * Usage:
 *   nextflow run main.nf --input data/cytochrome_c.fasta --outdir results
 *
 * Default parameter values (input, outdir, bootstrap_reps) live in
 * nextflow.config. Compatible with the strict syntax parser (Nextflow 25.10+).
 */

/*
 * Step 1: Multiple sequence alignment with MAFFT
 */
process ALIGN {
    publishDir params.outdir, mode: 'copy'
    tag "mafft"

    input:
        path fasta

    output:
        path "aligned.fasta"

    script:
    """
    mafft --auto ${fasta} > aligned.fasta
    """
}

/*
 * Step 2a/2b: Distance-based trees (NJ and UPGMA), each with a phylogram
 */
process NJ_TREE {
    publishDir params.outdir, mode: 'copy'
    tag "nj"

    input:
        path aligned

    output:
        path "tree_nj.nwk", emit: newick
        path "tree_nj.png", emit: png

    script:
    """
    build_tree.py ${aligned} nj tree_nj
    """
}

process UPGMA_TREE {
    publishDir params.outdir, mode: 'copy'
    tag "upgma"

    input:
        path aligned

    output:
        path "tree_upgma.nwk"
        path "tree_upgma.png", emit: png

    script:
    """
    build_tree.py ${aligned} upgma tree_upgma
    """
}

/*
 * Step 3: Alternate renderings of the NJ tree topology
 */
process CLADOGRAM {
    publishDir params.outdir, mode: 'copy'
    tag "cladogram"

    input:
        path nj_newick

    output:
        path "cladogram.png", emit: png

    script:
    """
    plot_cladogram.py ${nj_newick} cladogram.png
    """
}

process CIRCULAR_TREE {
    publishDir params.outdir, mode: 'copy'
    tag "circular"

    input:
        path nj_newick

    output:
        path "circular.png", emit: png

    script:
    """
    plot_circular.py ${nj_newick} circular.png
    """
}

/*
 * Step 4: Bootstrap confidence tree (independent of the NJ tree above --
 * it rebuilds trees from resampled alignment columns)
 */
process BOOTSTRAP_TREE {
    publishDir params.outdir, mode: 'copy'
    tag "bootstrap"

    input:
        path aligned

    output:
        path "bootstrap.nwk"
        path "bootstrap.png", emit: png

    script:
    """
    bootstrap_tree.py ${aligned} ${params.bootstrap_reps} bootstrap
    """
}

/*
 * Step 5: Combine all figures into one overview image
 */
process GALLERY {
    publishDir params.outdir, mode: 'copy'
    tag "gallery"

    input:
        path pngs

    output:
        path "gallery.png"

    script:
    """
    make_gallery.py gallery.png
    """
}

workflow {
    log.info """
             PHYLO PIPELINE
             ==============
             input          : ${params.input}
             outdir         : ${params.outdir}
             bootstrap_reps : ${params.bootstrap_reps}
             """.stripIndent()

    fasta_ch   = channel.fromPath(params.input, checkIfExists: true)

    aligned_ch = ALIGN(fasta_ch)

    NJ_TREE(aligned_ch)
    UPGMA_TREE(aligned_ch)
    CLADOGRAM(NJ_TREE.out.newick)
    CIRCULAR_TREE(NJ_TREE.out.newick)
    BOOTSTRAP_TREE(aligned_ch)

    // wait for all five figures, then tile them into one overview
    all_pngs = NJ_TREE.out.png
        .mix(UPGMA_TREE.out.png, CLADOGRAM.out.png,
             CIRCULAR_TREE.out.png, BOOTSTRAP_TREE.out.png)
        .collect()
    GALLERY(all_pngs)
}
