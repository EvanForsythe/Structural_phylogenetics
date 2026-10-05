# Structural Phylogenetics Pipeline

Notes and scripts for running the structural phylogenetics workflow on the Novus HPC.

The basic workflow is:

1. QC the input protein FASTA
2. Predict protein structures with ColabFold (using AlphaFold2 models)
3. Collect the top-ranked structures and summarize prediction confidence
4. Convert structures to 3Di code and generate a structural alignment with FoldMason
5. Infer:
   - a structural phylogeny from the 3Di alignment
   - a standard amino-acid phylogeny from a MAFFT alignment
6. Plot and compare the structural and standard trees

---

## Environment setup

Two separate Conda environments are used.

### 1. ColabFold environment

ColabFold is kept in a separate environment because of its GPU/JAX/CUDA dependencies.

The working environment was installed on scratch to avoid storage limits in the home directory:

```bash
mkdir -p /scratch/forsythe/conda_envs

conda create -p /scratch/forsythe/conda_envs/colabfold python=3.11

conda activate /scratch/forsythe/conda_envs/colabfold
```

pip cache was also moved to scratch:

```bash
mkdir -p /scratch/forsythe/pip_cache
export PIP_CACHE_DIR=/scratch/forsythe/pip_cache
```

Install ColabFold, JAX, and OpenMM:

```bash
python -m pip install --resume-retries 20 --timeout 120 "colabfold[alphafold,openmm]" "jax[cuda12]" "openmm[cuda12]"
```

On Novus, some of the large NVIDIA dependencies had to be downloaded and installed separately when the main pip installation had trouble completing. The packages that required this during the original installation were:

```text
nvidia-cudnn-cu12==9.27.0.42
nvidia-nccl-cu12==2.32.3
nvidia-nvshmem-cu12==3.8.0
```

After installing those dependencies, rerunning the ColabFold pip command completed the installation.

**Important:** The working ColabFold environment uses the pip-installed CUDA libraries. Do not load the cluster `cuda/12.6` module when running the pipeline.

### Test GPU access

Request an A40 interactively:

```bash
srun -p forsythe.q -A forsythe --gres=gpu:a40:1 --cpus-per-task=4 --mem=16G --pty bash
```

Activate ColabFold:

```bash
conda activate /scratch/forsythe/conda_envs/colabfold
```

Check that JAX sees the GPU:

```bash
python -c "import jax; print('JAX version:', jax.__version__); print('Devices:', jax.devices())"
```

The successful setup returned a CUDA device similar to:

```text
Devices: [CudaDevice(id=0)]
```

### ColabFold model cache

The downloaded ColabFold model files are stored on scratch:

```bash
mkdir -p /scratch/forsythe/Structural_phylogenetics/colabfold_cache

cp -a ~/.cache/colabfold/. \
    /scratch/forsythe/Structural_phylogenetics/colabfold_cache/
```

`02_colabfold.sbatch` points ColabFold to this directory with the `--data` option.

---

### 2. Structural phylogenetics environment

The second environment contains the downstream structural and phylogenetic software.

The environment was initially created by cloning a Foldseek environment:

```bash
conda create --name structural_phylo --clone foldseek
conda activate structural_phylo
```

FoldMason, IQ-TREE, MAFFT, and Biopython were then installed:

```bash
conda install -c conda-forge -c bioconda foldmason iqtree biopython mafft -y
```

The `Q.3Di.AF` substitution matrix, obtained from [here](https://edmond.mpg.de/dataset.xhtml?persistentId=doi:10.17617/3.1MJJBH), was used for 3Di phylogenetic inference.


---

# Running the pipeline

Each gene family gets its own output directory. For example:

```text
RNAase_T2_RNS1_family/
```

with an input FASTA such as:

```text
RNAase_T2_RNS1_family/HOG0015056.fa
```

The scripts themselves are run from the repo folder:

```text
/scratch/forsythe/Structural_phylogenetics/
```


---

## Step 1 — Sequence QC

Activate the structural phylogenetics environment:

```bash
conda activate structural_phylo
```

Run:

```bash
python 01_qc_fasta.py --input_fasta RNAase_T2_RNS1_family/HOG0015056.fa --output_dir RNAase_T2_RNS1_family/
```

This checks sequence lengths, duplicate IDs, and invalid characters.

Main output:

```text
RNAase_T2_RNS1_family/qc/sequence_qc.tsv
```

---

## Step 2 — AlphaFold2 structure prediction

Submit the ColabFold GPU job:

```bash
sbatch 02_colabfold.sbatch /scratch/forsythe/Structural_phylogenetics/RNAase_T2_RNS1_family/HOG0015056.fa /scratch/forsythe/Structural_phylogenetics/RNAase_T2_RNS1_family
```

The SLURM script activates the ColabFold environment itself.

Main output directory:

```text
RNAase_T2_RNS1_family/colabfold_results/
```

Each protein should eventually have a `.done.txt` file.

Progress can be checked with:

```bash
watch -n 30 'find RNAase_T2_RNS1_family/colabfold_results -name "*.done.txt" | wc -l'
```

---

## Step 3 — Collect structures and structure QC

Activate:

```bash
conda activate structural_phylo
```

Run:

```bash
python 03_collect_structures.py --input_fasta /scratch/forsythe/Structural_phylogenetics/RNAase_T2_RNS1_family/HOG0015056.fa --output_dir /scratch/forsythe/Structural_phylogenetics/RNAase_T2_RNS1_family
```

This selects the rank-1 AlphaFold structure for each protein and summarizes prediction confidence.

Important outputs:

```text
RNAase_T2_RNS1_family/structures/
RNAase_T2_RNS1_family/qc/structure_qc.tsv
```

---

## Step 4 — FoldMason structural alignment

Submit:

```bash
sbatch 04_foldmason.sbatch /scratch/forsythe/Structural_phylogenetics/RNAase_T2_RNS1_family/HOG0015056.fa /scratch/forsythe/Structural_phylogenetics/RNAase_T2_RNS1_family
```

FoldMason aligns the predicted structures and generates both amino-acid and 3Di alignments.

Important outputs include:

```text
RNAase_T2_RNS1_family/foldmason/HOG0015056_foldmason_3di.fa
RNAase_T2_RNS1_family/foldmason/HOG0015056_foldmason_aa.fa
```

---

## Step 5 — Structural and standard phylogenies

Submit:

```bash
sbatch 05_phylogenies.sbatch /scratch/forsythe/Structural_phylogenetics/RNAase_T2_RNS1_family/HOG0015056.fa /scratch/forsythe/Structural_phylogenetics/RNAase_T2_RNS1_family
```

This performs two analyses.

### Structural phylogeny

Uses:

```text
FoldMason 3Di alignment
        ↓
Q.3Di.AF substitution model
        ↓
IQ-TREE
        ↓
100 ultrafast bootstrap replicates
```

### Standard phylogeny

Uses:

```text
Original protein FASTA
        ↓
MAFFT --auto
        ↓
IQ-TREE ModelFinder
        ↓
best-fitting amino-acid model
        ↓
100 ultrafast bootstrap replicates
```

The resulting trees are:

```text
RNAase_T2_RNS1_family/phylogenies/structural/HOG0015056_3di.treefile
```

and:

```text
RNAase_T2_RNS1_family/phylogenies/standard/HOG0015056_aa.treefile
```

These trees are inferred from the same original set of proteins and can therefore be compared directly.

---

## Step 6 — Plot and compare trees

The convergent species are specified in a TSV file such as:

```text
RNAase_T2_RNS1_family/convergent_trait.tsv
```

Example:

```text
Convergent_species
CepFol
AldVes
DioMus
DroSpa
NepGra
BybGig
GenAur
PinMor
UtrGib
SarPur
HelTat
```

Run:

```bash
python 06_compare_phylogenies.py \
    --standard_tree RNAase_T2_RNS1_family/phylogenies/standard/HOG0015056_aa.treefile \
    --structural_tree RNAase_T2_RNS1_family/phylogenies/structural/HOG0015056_3di.treefile \
    --trait_tsv RNAase_T2_RNS1_family/convergent_trait.tsv \
    --output_dir RNAase_T2_RNS1_family/tree_comparison
```

This:

- plots the standard tree
- plots the 3Di structural tree
- highlights convergent taxa in red
- displays bootstrap support
- calculates all pairwise topology-only distances between convergent tips
- compares the distance distributions between the two trees

The topological distance is the number of tree edges separating two tips and therefore does not depend on estimated branch lengths.

Main outputs:

```text
tree_comparison/standard_tree_unrooted.pdf
tree_comparison/structural_3di_tree_unrooted.pdf

tree_comparison/standard_convergent_pairwise_distances.tsv
tree_comparison/structural_convergent_pairwise_distances.tsv

tree_comparison/convergent_tip_distance_KDE.pdf
```

---

## Environment reminder

For structure prediction:

```bash
conda activate /scratch/forsythe/conda_envs/colabfold
```

For everything else:

```bash
conda activate structural_phylo
```

