# nr-ligand-analysis

Cleaned and reproducible repository for the analysis of human nuclear receptor ligands.

## What is in the repo

- `notebooks/nr_ligand_analysis.ipynb`: the canonical end-to-end notebook
- `data/input/`: the curated ligand dataset used by the notebook
- `data/processed/`: lightweight summary tables generated from the tracked input data
- `figures/`: polished overview figures generated from the tracked input data
- `results_summary.md`: quick project findings and current data-quality notes
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

This also refreshes:

- `data/processed/data_integrity_summary.csv`
- `data/processed/missingness_summary.csv`
- `data/processed/duplicate_candidates.csv`
- `data/processed/receptor_family_map.csv`
- `figures/receptor_family_map.png`
- `results_summary.md`

Run the full descriptor and chemical-space workflow inside:

```bash
notebooks/nr_ligand_analysis.ipynb
```

## Start here

If you are opening the repo for the first time, use this order:

1. Read `results_summary.md` for the current high-level findings and QC notes.
2. Open `figures/` for the quick visual overview.
3. Open `notebooks/nr_ligand_analysis.ipynb` for the full reproducible workflow.
4. Check `data/processed/duplicate_candidates.csv` if you want to review likely duplicate ligand entries.

## Current cleanup decisions

- One canonical notebook replaces the earlier overlapping notebook variants.
- Paths are repo-relative instead of environment-specific.
- Output naming is standardized under `data/processed/` and `figures/`.
- Draft notebooks are preserved in `drafts/` instead of mixed into the main deliverable.
- Overcrowded duplicate receptor t-SNE plots are replaced in the final notebook with:
  - one t-SNE colored by nuclear receptor family
  - one t-SNE colored by top receptors, with the rest grouped as `Other`

## Current findings snapshot

- The tracked input dataset has 1,377 ligand rows across 34 receptors, 3 nuclear receptor families, and 9 chemical families.
- `ESR1` is currently the largest receptor group in the tracked dataset.
- `Non-steroid` is currently the largest nuclear receptor family.
- The current data has 1 missing `Ligand_ID` row and 22 duplicate `Ligand_ID` / `InChIKey` candidates flagged for review.

## Known limitations

- The full descriptor, scaffold, and t-SNE workflow still requires a local RDKit-enabled environment to execute end-to-end.
- Duplicate candidates are exported for review, but they are not automatically removed because some may reflect biologically meaningful repeated structures across PDB entries.
- The overview figures in the repo are generated from the tracked input data; the chemistry-heavy figures are produced inside the canonical notebook.

## Input schema expected by the notebook

The final notebook validates these required columns before analysis:

- `Ligand_ID`
- `Ligand_Name`
- `SMILES`
- `UniProtID`
- `QueryGene`
- `NR_Family`
- `Chemical_Family`
