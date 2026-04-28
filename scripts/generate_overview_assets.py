from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
INPUT_PATH = ROOT / "data" / "input" / "ligands_with_final_chemical_classification.csv"
PROCESSED_DIR = ROOT / "data" / "processed"
FIGURES_DIR = ROOT / "figures"
RESULTS_SUMMARY_PATH = ROOT / "results_summary.md"

PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

sns.set_theme(style="whitegrid", context="talk")
plt.rcParams["figure.dpi"] = 130
plt.rcParams["savefig.dpi"] = 300

df = pd.read_csv(INPUT_PATH)

required_columns = [
    "Ligand_ID",
    "Ligand_Name",
    "SMILES",
    "UniProtID",
    "QueryGene",
    "NR_Family",
    "Chemical_Family",
]

missingness_summary = pd.DataFrame(
    {
        "Column": df.columns,
        "Missing_Count": [int(df[column].isna().sum()) for column in df.columns],
        "Missing_Percent": [round(float(df[column].isna().mean() * 100), 3) for column in df.columns],
        "Unique_Non_Null": [int(df[column].nunique(dropna=True)) for column in df.columns],
        "Dtype": [str(df[column].dtype) for column in df.columns],
    }
).sort_values(["Missing_Count", "Column"], ascending=[False, True])

dedup_key = df["InChIKey"].fillna(df["Ligand_ID"]) if "InChIKey" in df.columns else df["Ligand_ID"]
duplicate_mask = dedup_key.duplicated(keep=False)
duplicate_candidates = df.loc[duplicate_mask, required_columns + ["PDB_ID", "InChIKey"]].copy()
duplicate_candidates.insert(0, "DedupKey", dedup_key.loc[duplicate_mask].astype(str))
duplicate_candidates = duplicate_candidates.sort_values(["DedupKey", "PDB_ID", "Ligand_Name"])

receptor_summary = (
    df.groupby("QueryGene", dropna=False)
    .agg(
        Ligand_Count=("QueryGene", "size"),
        Unique_Ligands=("Ligand_Name", "nunique"),
        Unique_Families=("NR_Family", "nunique"),
    )
    .sort_values("Ligand_Count", ascending=False)
    .reset_index()
    .rename(columns={"QueryGene": "Receptor"})
)

family_summary = (
    df.groupby("NR_Family", dropna=False)
    .agg(
        Ligand_Count=("NR_Family", "size"),
        Unique_Receptors=("QueryGene", "nunique"),
        Unique_Chemical_Families=("Chemical_Family", "nunique"),
        Mean_Input_Molecular_Weight=("Molecular_Weight", "mean"),
    )
    .sort_values("Ligand_Count", ascending=False)
    .reset_index()
)

chemical_family_summary = (
    df.groupby("Chemical_Family", dropna=False)
    .agg(
        Ligand_Count=("Chemical_Family", "size"),
        Unique_Receptors=("QueryGene", "nunique"),
        Unique_NR_Families=("NR_Family", "nunique"),
    )
    .sort_values("Ligand_Count", ascending=False)
    .reset_index()
)

receptor_family_map = pd.crosstab(df["QueryGene"], df["NR_Family"]).sort_index()

integrity_summary = pd.DataFrame(
    [
        {
            "Total_Rows": int(len(df)),
            "Unique_Receptors": int(df["QueryGene"].nunique()),
            "Unique_NR_Families": int(df["NR_Family"].nunique()),
            "Unique_Chemical_Families": int(df["Chemical_Family"].nunique()),
            "Missing_Ligand_ID": int(df["Ligand_ID"].isna().sum()),
            "Duplicate_Ligand_ID": int(df["Ligand_ID"].dropna().duplicated().sum()),
            "Duplicate_InChIKey": int(df["InChIKey"].dropna().duplicated().sum()) if "InChIKey" in df.columns else 0,
            "Duplicate_Candidate_Rows": int(duplicate_mask.sum()),
        }
    ]
)

missingness_summary.to_csv(PROCESSED_DIR / "missingness_summary.csv", index=False)
duplicate_candidates.to_csv(PROCESSED_DIR / "duplicate_candidates.csv", index=False)
receptor_summary.to_csv(PROCESSED_DIR / "receptor_summary.csv", index=False)
family_summary.to_csv(PROCESSED_DIR / "family_summary.csv", index=False)
chemical_family_summary.to_csv(PROCESSED_DIR / "chemical_family_summary.csv", index=False)
receptor_family_map.to_csv(PROCESSED_DIR / "receptor_family_map.csv")
integrity_summary.to_csv(PROCESSED_DIR / "data_integrity_summary.csv", index=False)

fig, axes = plt.subplots(1, 2, figsize=(18, 6))

sns.barplot(
    data=receptor_summary.head(10),
    x="Ligand_Count",
    y="Receptor",
    color="#4C78A8",
    ax=axes[0],
)
axes[0].set_title("Top 10 Receptors by Ligand Count")
axes[0].set_xlabel("Ligand count")
axes[0].set_ylabel("Receptor")

sns.barplot(
    data=family_summary,
    x="Ligand_Count",
    y="NR_Family",
    color="#72B7B2",
    ax=axes[1],
)
axes[1].set_title("Ligand Count by Nuclear Receptor Family")
axes[1].set_xlabel("Ligand count")
axes[1].set_ylabel("NR family")

plt.tight_layout()
plt.savefig(FIGURES_DIR / "composition_summary.png", bbox_inches="tight")
plt.close()

fig, axes = plt.subplots(1, 2, figsize=(18, 6))

sns.histplot(df["Molecular_Weight"].dropna(), bins=30, ax=axes[0], color="#4C78A8")
axes[0].set_title("Input Molecular Weight Distribution")
axes[0].set_xlabel("Molecular weight (Da)")
axes[0].set_ylabel("Ligand count")

top_chemical_families = chemical_family_summary.head(12).copy()
sns.barplot(
    data=top_chemical_families,
    x="Ligand_Count",
    y="Chemical_Family",
    color="#F58518",
    ax=axes[1],
)
axes[1].set_title("Top Chemical Families")
axes[1].set_xlabel("Ligand count")
axes[1].set_ylabel("Chemical family")

plt.tight_layout()
plt.savefig(FIGURES_DIR / "dataset_overview.png", bbox_inches="tight")
plt.close()

fig, ax = plt.subplots(figsize=(8, 10))
sns.heatmap(receptor_family_map, cmap="Blues", linewidths=0.5, cbar_kws={"label": "Ligand count"}, ax=ax)
ax.set_title("Receptor-to-Family Coverage")
ax.set_xlabel("NR family")
ax.set_ylabel("Receptor")
plt.tight_layout()
plt.savefig(FIGURES_DIR / "receptor_family_map.png", bbox_inches="tight")
plt.close()

top_receptors = receptor_summary.head(5)[["Receptor", "Ligand_Count"]]
top_chemical = chemical_family_summary.head(5)[["Chemical_Family", "Ligand_Count"]]

results_summary_lines = [
    "# Results Snapshot",
    "",
    "This summary is generated from the tracked input dataset in `data/input/ligands_with_final_chemical_classification.csv`.",
    "",
    "## Dataset scale",
    "",
    f"- Total ligand rows: {len(df):,}",
    f"- Unique receptors: {df['QueryGene'].nunique()}",
    f"- Nuclear receptor families: {df['NR_Family'].nunique()}",
    f"- Chemical families: {df['Chemical_Family'].nunique()}",
    "",
    "## Current high-level findings",
    "",
    f"- The largest receptor group is `{top_receptors.iloc[0]['Receptor']}` with {int(top_receptors.iloc[0]['Ligand_Count'])} ligand rows.",
    f"- The largest nuclear receptor family is `{family_summary.iloc[0]['NR_Family']}` with {int(family_summary.iloc[0]['Ligand_Count'])} ligand rows.",
    f"- The most common chemical family is `{top_chemical.iloc[0]['Chemical_Family']}` with {int(top_chemical.iloc[0]['Ligand_Count'])} ligand rows.",
    "",
    "## Data-quality notes",
    "",
    f"- Missing `Ligand_ID` rows: {int(df['Ligand_ID'].isna().sum())}",
    f"- Duplicate `Ligand_ID` entries: {int(df['Ligand_ID'].dropna().duplicated().sum())}",
    f"- Duplicate `InChIKey` entries: {int(df['InChIKey'].dropna().duplicated().sum()) if 'InChIKey' in df.columns else 0}",
    "- Duplicate rows are exported to `data/processed/duplicate_candidates.csv` for manual review.",
    "",
    "## Top receptors",
    "",
]

for _, row in top_receptors.iterrows():
    results_summary_lines.append(f"- `{row['Receptor']}`: {int(row['Ligand_Count'])} ligand rows")

results_summary_lines.extend(
    [
        "",
        "## Top chemical families",
        "",
    ]
)

for _, row in top_chemical.iterrows():
    results_summary_lines.append(f"- `{row['Chemical_Family']}`: {int(row['Ligand_Count'])} ligand rows")

results_summary_lines.extend(
    [
        "",
        "## Generated files",
        "",
        "- `data/processed/data_integrity_summary.csv`",
        "- `data/processed/missingness_summary.csv`",
        "- `data/processed/duplicate_candidates.csv`",
        "- `data/processed/receptor_summary.csv`",
        "- `data/processed/family_summary.csv`",
        "- `data/processed/chemical_family_summary.csv`",
        "- `data/processed/receptor_family_map.csv`",
        "- `figures/composition_summary.png`",
        "- `figures/dataset_overview.png`",
        "- `figures/receptor_family_map.png`",
    ]
)

RESULTS_SUMMARY_PATH.write_text("\n".join(results_summary_lines) + "\n", encoding="utf-8")

print("Generated overview CSV summaries, QC outputs, figures, and results_summary.md.")
