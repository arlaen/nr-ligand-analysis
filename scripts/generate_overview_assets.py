from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
INPUT_PATH = ROOT / "data" / "input" / "ligands_with_final_chemical_classification.csv"
PROCESSED_DIR = ROOT / "data" / "processed"
FIGURES_DIR = ROOT / "figures"

PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

sns.set_theme(style="whitegrid", context="talk")
plt.rcParams["figure.dpi"] = 130
plt.rcParams["savefig.dpi"] = 300

df = pd.read_csv(INPUT_PATH)

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

receptor_summary.to_csv(PROCESSED_DIR / "receptor_summary.csv", index=False)
family_summary.to_csv(PROCESSED_DIR / "family_summary.csv", index=False)
chemical_family_summary.to_csv(PROCESSED_DIR / "chemical_family_summary.csv", index=False)

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

print("Generated overview CSV summaries and figures.")
