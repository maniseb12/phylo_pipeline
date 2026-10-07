#  Phylogenetic Tree Pipeline

A beginner-friendly [Nextflow](https://www.nextflow.io/) pipeline that turns a
set of related protein sequences into **five different views** of their
evolutionary tree. Sample data: **Cytochrome c** from 9 species (human, chimp,
rhesus monkey, horse, cow, chicken, rattlesnake, tuna, yeast).

![Pipeline overview](results/gallery.png)

## What it produces

| # | Figure | File | What it shows |
|---|--------|------|---------------|
| 1 | Neighbor-Joining phylogram | `tree_nj.png` | Branch length = amount of sequence change |
| 2 | UPGMA phylogram | `tree_upgma.png` | Assumes a constant mutation rate; compare with NJ |
| 3 | Cladogram | `cladogram.png` | Topology only, tips aligned; easiest way to read relationships |
| 4 | Circular tree | `circular.png` | Radial layout with a colored group ring |
| 5 | Bootstrap-supported tree | `bootstrap.png` | Dots show % of resampled trees agreeing with each split |
| ★ | Gallery | `gallery.png` | All five in one overview image |

Each tree is also saved as a Newick file (`.nwk`) that you can open in
[iTOL](https://itol.embl.de/) or FigTree. The alignment is saved as
`aligned.fasta`.

## Pipeline

```
                           ┌─► NJ_TREE ──┬─► CLADOGRAM ─────┐
                           │             └─► CIRCULAR_TREE ─┤
cytochrome_c.fasta ─► ALIGN├─► UPGMA_TREE ──────────────────┼─► GALLERY
                  (MAFFT)  └─► BOOTSTRAP_TREE ──────────────┘
```

## Project layout

```
phylo_pipeline/
├── main.nf                → pipeline definition (steps + how they connect)
├── nextflow.config        → default settings: input path, output dir, bootstrap replicates
├── data/
│   └── cytochrome_c.fasta → input: unaligned protein sequences
├── bin/                   → worker scripts (auto-added to PATH by Nextflow)
│   ├── phylo_style.py       shared colors, layout maths, drawing helpers
│   ├── build_tree.py        NJ / UPGMA tree + phylogram
│   ├── plot_cladogram.py    topology-only view
│   ├── plot_circular.py     circular layout
│   ├── bootstrap_tree.py    bootstrap support on the NJ tree
│   └── make_gallery.py      combines all figures into gallery.png
└── results/               →  output (regenerated on each run)
```

## Requirements

| Tool | Notes |
|------|-------|
| [Nextflow](https://www.nextflow.io/docs/latest/install.html) | Needs Java 11 or newer. Works with the current strict-syntax parser (Nextflow 25.10+). |
| [MAFFT](https://mafft.cbrc.jp/alignment/software/) | Must be on your `PATH` |
| Python 3 | With `biopython`, `matplotlib`, `numpy` |

## Setup on Linux

### 1. Get the files and make the scripts executable

```bash
unzip phylo_pipeline-1.zip
cd phylo_pipeline
chmod +x bin/*.py
```

The `chmod` step matters. Unzipping often drops the executable bit, and
Nextflow then fails with `exit status 126 ... Permission denied`.

### 2. Install dependencies (pick one option)

**Option A: conda (recommended, installs everything together)**

```bash
conda create -n phylo -c conda-forge -c bioconda nextflow mafft biopython matplotlib numpy
conda activate phylo
```

**Option B: system packages + virtual environment (Ubuntu/Debian)**

```bash
sudo apt update
sudo apt install -y default-jre mafft python3 python3-pip python3-venv

python3 -m venv .venv
source .venv/bin/activate
pip install biopython matplotlib numpy

curl -s https://get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin/      # or move it to ~/bin
```

Use only **one** environment at a time. If your prompt shows both `(conda-env)`
and `(.venv)`, run `deactivate` to leave the venv. Nextflow runs the scripts
with whichever `python3` comes first on your `PATH`.

### 3. Check the setup

```bash
java -version
nextflow -version
mafft --version
python3 -c "import Bio, matplotlib, numpy; print('python libs ok')"
```

## Run it

```bash
nextflow run main.nf --input data/cytochrome_c.fasta --outdir results --bootstrap_reps 100
```

| Option | Meaning | Default (in `nextflow.config`) |
|--------|---------|--------------------------------|
| `--input` | Unaligned protein FASTA | `data/cytochrome_c.fasta` |
| `--outdir` | Where results are written | `results` |
| `--bootstrap_reps` | Number of bootstrap replicates | `100` |

Add `-resume` to skip steps whose inputs haven't changed.

### Running without Nextflow

The scripts also work on their own:

```bash
export PATH="$PWD/bin:$PATH"
mkdir -p results
mafft --auto data/cytochrome_c.fasta > results/aligned.fasta
build_tree.py results/aligned.fasta nj results/tree_nj
build_tree.py results/aligned.fasta upgma results/tree_upgma
plot_cladogram.py results/tree_nj.nwk results/cladogram.png
plot_circular.py results/tree_nj.nwk results/circular.png
bootstrap_tree.py results/aligned.fasta 100 results/bootstrap
cd results && make_gallery.py gallery.png
```

## Use your own sequences

1. Put a protein FASTA in `data/` and pass it with `--input`.
2. Open `bin/phylo_style.py` and edit the `GROUPS` dictionary to map each
   sequence ID (the first word of the FASTA header) to a group. Colors come
   from `GROUP_COLORS`. Unlisted names are drawn grey as "Other".

## Troubleshooting

| Problem | Fix |
|---------|-----|
| `Statements cannot be mixed with script declarations` | You are using an old copy of `main.nf`. Use the current one, where defaults live in `nextflow.config` and the banner is inside `workflow { }`. |
| `exit status 126` / `Permission denied` on a `.py` file | Run `chmod +x bin/*.py` |
| `mafft: command not found` | Install MAFFT (`sudo apt install mafft` or via conda) and check with `which mafft` |
| `ModuleNotFoundError: No module named 'Bio'` | The environment with biopython isn't active, or two environments are stacked. Deactivate extras and re-check with the Python one-liner above. |
| Matplotlib display errors on a headless server | `export MPLBACKEND=Agg` before running |
| Pipeline stopped part-way | Fix the error and rerun the same command with `-resume`. Use `nextflow log` and the `work/` directory to inspect failed tasks. |

To tidy up between runs: `nextflow clean -f` removes old work directories.

## Reading the results

- Human, chimp and rhesus monkey cluster tightly (identical or nearly identical
  cytochrome c); yeast branches off first with the longest branch. This matches
  known evolutionary history and is a good sanity check.
- Bootstrap support is modest for the mammal/bird splits. That is expected:
  with 9 short, highly conserved sequences there is little signal to resample.
  Longer alignments and more sequences give stronger support.

## Ideas to extend

- Swap the `identity` distance for `blosum62` in `DistanceCalculator`
- Try a maximum-likelihood tree with IQ-TREE or FastTree instead of NJ
- Add a Docker/Conda profile in `nextflow.config` for reproducibility
