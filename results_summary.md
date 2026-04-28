# Results Snapshot

This summary is generated from the tracked input dataset in `data/input/ligands_with_final_chemical_classification.csv`.

## Dataset scale

- Total ligand rows: 1,377
- Unique receptors: 34
- Nuclear receptor families: 3
- Chemical families: 9

## Current high-level findings

- The largest receptor group is `ESR1` with 406 ligand rows.
- The largest nuclear receptor family is `Non-steroid` with 624 ligand rows.
- The most common chemical family is `Fatty acids` with 543 ligand rows.

## Data-quality notes

- Missing `Ligand_ID` rows: 1
- Duplicate `Ligand_ID` entries: 22
- Duplicate `InChIKey` entries: 22
- Duplicate rows are exported to `data/processed/duplicate_candidates.csv` for manual review.

## Top receptors

- `ESR1`: 406 ligand rows
- `PPARG`: 219 ligand rows
- `RORC`: 139 ligand rows
- `NR1H4`: 80 ligand rows
- `RXRA`: 71 ligand rows

## Top chemical families

- `Fatty acids`: 543 ligand rows
- `Steroids`: 507 ligand rows
- `Orphan ligands`: 176 ligand rows
- `Bile acids`: 48 ligand rows
- `Synthetic / drug-like molecules`: 40 ligand rows

## Generated files

- `data/processed/data_integrity_summary.csv`
- `data/processed/missingness_summary.csv`
- `data/processed/duplicate_candidates.csv`
- `data/processed/receptor_summary.csv`
- `data/processed/family_summary.csv`
- `data/processed/chemical_family_summary.csv`
- `data/processed/receptor_family_map.csv`
- `figures/composition_summary.png`
- `figures/dataset_overview.png`
- `figures/receptor_family_map.png`
