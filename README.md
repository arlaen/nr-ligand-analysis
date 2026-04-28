# nr-ligand-analysis

Cleaned and reproducible repository for the analysis of human nuclear receptor ligands.

## What is in the repo

- `notebooks/nr_ligand_analysis.ipynb`: the canonical end-to-end notebook
- `data/input/`: the curated ligand dataset used by the notebook
- `data/processed/`: lightweight summary tables generated from the tracked input data
- `figures/`: polished overview figures generated from the tracked input data
- `drafts/`: archived notebook variants kept only for reference during cleanup
- `scripts/`: helper scripts to rebuild the notebook and lightweight repo assets

## Project structure

```text
.
|-- data/
|   |-- input/
|   |-- processed/
|-- drafts/
|-- figures/
|-- notebooks/
|-- scripts/
|-- environment.yml
`-- README.md
```

## Recommended environment

Use Conda or Mamba so RDKit installs reliably:

```bash
conda env create -f environment.yml
conda activate nr-ligand-analysis
jupyter lab
```

## How to rebuild the repo artifacts

Generate the cleaned notebook:

```bash
python scripts/build_notebook.py
```

Generate the lightweight CSV summaries and overview figures that do not require RDKit:

```bash
python scripts/generate_overview_assets.py
```

Run the full descriptor and chemical-space workflow inside:

```bash
notebooks/nr_ligand_analysis.ipynb
```

## Current cleanup decisions

- One canonical notebook replaces the earlier overlapping notebook variants.
- Paths are repo-relative instead of environment-specific.
- Output naming is standardized under `data/processed/` and `figures/`.
- Draft notebooks are preserved in `drafts/` instead of mixed into the main deliverable.
- Overcrowded duplicate receptor t-SNE plots are replaced in the final notebook with:
  - one t-SNE colored by nuclear receptor family
  - one t-SNE colored by top receptors, with the rest grouped as `Other`

## Input schema expected by the notebook

The final notebook validates these required columns before analysis:

- `Ligand_ID`
- `Ligand_Name`
- `SMILES`
- `UniProtID`
- `QueryGene`
- `NR_Family`
- `Chemical_Family`
