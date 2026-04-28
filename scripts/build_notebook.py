from pathlib import Path

import nbformat as nbf


ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK_PATH = ROOT / "notebooks" / "nr_ligand_analysis.ipynb"


def md(text: str):
    return nbf.v4.new_markdown_cell(text.strip() + "\n")


def code(text: str):
    return nbf.v4.new_code_cell(text.strip() + "\n")


cells = [
    md(
        """
        # Nuclear Receptor Ligand Analysis

        This notebook is the canonical, cleaned workflow for the project. It is designed to run from the GitHub repo after cloning, without manual path edits.

        Workflow:
        1. Load the curated ligand dataset from `data/input/`
        2. Validate the required schema
        3. Calculate molecular descriptors from SMILES
        4. Summarize receptor and family composition
        5. Generate polished overview and chemical-space figures
        6. Export processed tables to `data/processed/`
        """
    ),
    md(
        """
        ## Environment

        Recommended setup:

        ```bash
        conda env create -f environment.yml
        conda activate nr-ligand-analysis
        jupyter lab
        ```

        The notebook expects RDKit, pandas, matplotlib, seaborn, and scikit-learn to be installed.
        """
    ),
    code(
        """
        from pathlib import Path
        import warnings
        import re

        import numpy as np
        import pandas as pd
        import matplotlib.pyplot as plt
        import seaborn as sns
        from IPython.display import display

        from sklearn.manifold import TSNE
        from sklearn.preprocessing import StandardScaler

        from rdkit import Chem
        from rdkit.Chem import Crippen, Descriptors, Lipinski, QED, rdMolDescriptors
        from rdkit.Chem.Scaffolds import MurckoScaffold

        warnings.filterwarnings("ignore", category=FutureWarning)
        sns.set_theme(style="whitegrid", context="talk")
        plt.rcParams["figure.dpi"] = 130
        plt.rcParams["savefig.dpi"] = 300


        def find_repo_root(start: Path | None = None) -> Path:
            start = (start or Path.cwd()).resolve()
            for candidate in [start, *start.parents]:
                if (candidate / "data" / "input" / "ligands_with_final_chemical_classification.csv").exists():
                    return candidate
            raise FileNotFoundError("Could not locate the repo root from the current working directory.")


        ROOT = find_repo_root()
        DATA_DIR = ROOT / "data"
        INPUT_DIR = DATA_DIR / "input"
        PROCESSED_DIR = DATA_DIR / "processed"
        FIGURES_DIR = ROOT / "figures"

        PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
        FIGURES_DIR.mkdir(parents=True, exist_ok=True)

        INPUT_PATH = INPUT_DIR / "ligands_with_final_chemical_classification.csv"
        print(f"Repo root: {ROOT}")
        print(f"Input data: {INPUT_PATH}")
        """
    ),
    md(
        """
        ## Load and validate the dataset

        The canonical input is the curated ligand table with receptor and family annotations. We validate the expected columns up front so the notebook fails early and clearly if the file changes.
        """
    ),
    code(
        """
        df = pd.read_csv(INPUT_PATH)

        REQUIRED_COLUMNS = {
            "Ligand_ID",
            "Ligand_Name",
            "SMILES",
            "UniProtID",
            "QueryGene",
            "NR_Family",
            "Chemical_Family",
        }

        missing_columns = sorted(REQUIRED_COLUMNS - set(df.columns))
        if missing_columns:
            raise ValueError(f"Missing required columns: {missing_columns}")

        print(f"Rows: {len(df):,}")
        print(f"Columns: {len(df.columns)}")
        display(df.head())
        """
    ),
    md(
        """
        ## Dataset quality control

        Before doing descriptor-heavy analysis, check whether the tracked input file has missing identifiers, duplicate ligand records, or inconsistent receptor-to-family assignments. These QC outputs are useful both for debugging and for explaining the limitations of the dataset in the final report.
        """
    ),
    code(
        """
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
        duplicate_candidates = df.loc[
            duplicate_mask,
            ["PDB_ID", "Ligand_ID", "Ligand_Name", "InChIKey", "QueryGene", "NR_Family", "Chemical_Family"],
        ].copy()
        duplicate_candidates.insert(0, "DedupKey", dedup_key.loc[duplicate_mask].astype(str))
        duplicate_candidates = duplicate_candidates.sort_values(["DedupKey", "PDB_ID", "Ligand_Name"])

        receptor_family_map = pd.crosstab(df["QueryGene"], df["NR_Family"]).sort_index()

        data_integrity_summary = pd.DataFrame(
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
        receptor_family_map.to_csv(PROCESSED_DIR / "receptor_family_map.csv")
        data_integrity_summary.to_csv(PROCESSED_DIR / "data_integrity_summary.csv", index=False)

        display(data_integrity_summary)
        display(missingness_summary.head(10))
        display(duplicate_candidates.head(10))
        """
    ),
    md(
        """
        ## Helper functions

        These helpers standardize family naming, calculate molecular descriptors, and keep export logic consistent across the workflow.
        """
    ),
    code(
        """
        def slugify(value: str) -> str:
            value = re.sub(r"[^a-z0-9]+", "_", str(value).strip().lower())
            return value.strip("_")


        def safe_mol_from_smiles(smiles: str):
            if pd.isna(smiles) or not str(smiles).strip():
                return None
            return Chem.MolFromSmiles(str(smiles))


        def calculate_descriptors(smiles: str) -> dict:
            mol = safe_mol_from_smiles(smiles)
            if mol is None:
                return {
                    "Canonical_SMILES": None,
                    "MolWt": np.nan,
                    "LogP": np.nan,
                    "TPSA": np.nan,
                    "HBD": np.nan,
                    "HBA": np.nan,
                    "RotatableBonds": np.nan,
                    "RingCount": np.nan,
                    "AromaticRings": np.nan,
                    "StereoCenters": np.nan,
                    "FractionCSP3": np.nan,
                    "MolMR": np.nan,
                    "QED": np.nan,
                    "Lipinski_Compliant": False,
                    "Veber_Compliant": False,
                    "Scaffold": None,
                }

            canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
            mol_wt = Descriptors.MolWt(mol)
            logp = Crippen.MolLogP(mol)
            tpsa = Descriptors.TPSA(mol)
            hbd = Lipinski.NumHDonors(mol)
            hba = Lipinski.NumHAcceptors(mol)
            rotatable_bonds = Lipinski.NumRotatableBonds(mol)
            ring_count = Descriptors.RingCount(mol)
            aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
            stereo_centers = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
            fraction_csp3 = Lipinski.FractionCSP3(mol)
            mol_mr = Crippen.MolMR(mol)
            scaffold = Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol))

            lipinski_compliant = bool(mol_wt <= 500 and logp <= 5 and hbd <= 5 and hba <= 10)
            veber_compliant = bool(rotatable_bonds <= 10 and tpsa <= 140)

            return {
                "Canonical_SMILES": canonical_smiles,
                "MolWt": mol_wt,
                "LogP": logp,
                "TPSA": tpsa,
                "HBD": hbd,
                "HBA": hba,
                "RotatableBonds": rotatable_bonds,
                "RingCount": ring_count,
                "AromaticRings": aromatic_rings,
                "StereoCenters": stereo_centers,
                "FractionCSP3": fraction_csp3,
                "MolMR": mol_mr,
                "QED": QED.qed(mol),
                "Lipinski_Compliant": lipinski_compliant,
                "Veber_Compliant": veber_compliant,
                "Scaffold": scaffold if scaffold else None,
            }
        """
    ),
    md(
        """
        ## Calculate descriptors

        This step creates the master analysis table used by the downstream summaries and figures.
        """
    ),
    code(
        """
        descriptor_df = pd.DataFrame(df["SMILES"].map(calculate_descriptors).tolist())
        analysis_df = pd.concat([df.reset_index(drop=True), descriptor_df], axis=1)

        analysis_df["Receptor"] = analysis_df["QueryGene"].astype(str)
        analysis_df["Family_Slug"] = analysis_df["NR_Family"].map(slugify)
        analysis_df["Chemical_Family_Slug"] = analysis_df["Chemical_Family"].map(slugify)

        analysis_path = PROCESSED_DIR / "ligands_with_descriptors.csv"
        analysis_df.to_csv(analysis_path, index=False)

        print(f"Processed table saved to: {analysis_path}")
        display(analysis_df.head())
        """
    ),
    md(
        """
        ## Receptor and family summaries

        The project uses one naming scheme for receptor, family, and processed outputs.
        """
    ),
    code(
        """
        receptor_summary = (
            analysis_df.groupby("Receptor", dropna=False)
            .agg(
                Ligand_Count=("Receptor", "size"),
                Unique_Ligands=("Ligand_Name", "nunique"),
                Mean_MolWt=("MolWt", "mean"),
                Mean_LogP=("LogP", "mean"),
                Lipinski_Pass_Rate=("Lipinski_Compliant", "mean"),
                Veber_Pass_Rate=("Veber_Compliant", "mean"),
            )
            .sort_values("Ligand_Count", ascending=False)
            .reset_index()
        )

        family_summary = (
            analysis_df.groupby("NR_Family", dropna=False)
            .agg(
                Ligand_Count=("NR_Family", "size"),
                Unique_Receptors=("Receptor", "nunique"),
                Unique_Chemical_Families=("Chemical_Family", "nunique"),
                Mean_MolWt=("MolWt", "mean"),
                Mean_LogP=("LogP", "mean"),
                Mean_TPSA=("TPSA", "mean"),
                Mean_QED=("QED", "mean"),
                Lipinski_Pass_Rate=("Lipinski_Compliant", "mean"),
                Veber_Pass_Rate=("Veber_Compliant", "mean"),
            )
            .sort_values("Ligand_Count", ascending=False)
            .reset_index()
        )

        receptor_summary.to_csv(PROCESSED_DIR / "receptor_summary.csv", index=False)
        family_summary.to_csv(PROCESSED_DIR / "family_summary.csv", index=False)

        display(receptor_summary.head(10))
        display(family_summary)
        """
    ),
    md(
        """
        ## Scaffold summary

        The original scaffold logic counted distinct frequencies rather than distinct scaffolds. Here we calculate scaffold diversity directly from unique scaffold strings.
        """
    ),
    code(
        """
        scaffold_summary = (
            analysis_df.dropna(subset=["Scaffold"])
            .groupby("NR_Family")
            .agg(
                Ligand_Count=("NR_Family", "size"),
                Unique_Scaffolds=("Scaffold", "nunique"),
                Most_Common_Scaffold=("Scaffold", lambda s: s.value_counts().index[0]),
                Most_Common_Scaffold_Count=("Scaffold", lambda s: int(s.value_counts().iloc[0])),
            )
            .reset_index()
        )

        scaffold_summary["Scaffold_Diversity"] = (
            scaffold_summary["Unique_Scaffolds"] / scaffold_summary["Ligand_Count"]
        )
        scaffold_summary = scaffold_summary.sort_values("Scaffold_Diversity", ascending=False)
        scaffold_summary.to_csv(PROCESSED_DIR / "scaffold_summary.csv", index=False)
        display(scaffold_summary)
        """
    ),
    md(
        """
        ## Final figures

        The figure set is intentionally compact:
        - one overview panel for descriptor distributions
        - one receptor count chart
        - one family count chart
        - one t-SNE colored by **family**
        - one t-SNE colored by the **top receptors**, with all remaining receptors grouped as `Other`

        This avoids the duplicated, overcrowded legends from earlier drafts and makes each figure answer a different question.
        """
    ),
    code(
        """
        top_receptors = receptor_summary["Receptor"].head(12).tolist()
        analysis_df["Receptor_Group"] = analysis_df["Receptor"].where(
            analysis_df["Receptor"].isin(top_receptors), "Other"
        )

        receptor_palette = {
            receptor: color
            for receptor, color in zip(
                sorted(analysis_df["Receptor_Group"].dropna().unique()),
                sns.color_palette("tab20", n_colors=analysis_df["Receptor_Group"].nunique()),
            )
        }

        family_palette = {
            family: color
            for family, color in zip(
                sorted(analysis_df["NR_Family"].dropna().unique()),
                sns.color_palette("Set2", n_colors=analysis_df["NR_Family"].nunique()),
            )
        }

        fig, axes = plt.subplots(2, 2, figsize=(18, 12))

        sns.histplot(analysis_df["MolWt"].dropna(), bins=30, ax=axes[0, 0], color="#4C78A8")
        axes[0, 0].set_title("Molecular Weight Distribution")
        axes[0, 0].set_xlabel("Molecular Weight (Da)")

        sns.histplot(analysis_df["LogP"].dropna(), bins=30, ax=axes[0, 1], color="#F58518")
        axes[0, 1].set_title("Lipophilicity Distribution")
        axes[0, 1].set_xlabel("LogP")

        sns.histplot(analysis_df["TPSA"].dropna(), bins=30, ax=axes[1, 0], color="#54A24B")
        axes[1, 0].set_title("Polar Surface Area Distribution")
        axes[1, 0].set_xlabel("TPSA")

        lipinski_counts = (
            analysis_df["Lipinski_Compliant"]
            .map({True: "Compliant", False: "Non-compliant"})
            .value_counts()
        )
        axes[1, 1].pie(
            lipinski_counts.values,
            labels=lipinski_counts.index,
            autopct="%1.1f%%",
            colors=["#4C78A8", "#E45756"],
            startangle=90,
        )
        axes[1, 1].set_title("Lipinski Compliance")

        plt.tight_layout()
        plt.savefig(FIGURES_DIR / "property_overview.png", bbox_inches="tight")
        plt.show()

        fig, axes = plt.subplots(1, 2, figsize=(18, 6))
        sns.barplot(
            data=receptor_summary.head(10),
            x="Ligand_Count",
            y="Receptor",
            ax=axes[0],
            color="#4C78A8",
        )
        axes[0].set_title("Top 10 Receptors by Ligand Count")
        axes[0].set_xlabel("Ligand count")
        axes[0].set_ylabel("Receptor")

        sns.barplot(
            data=family_summary,
            x="Ligand_Count",
            y="NR_Family",
            ax=axes[1],
            color="#72B7B2",
        )
        axes[1].set_title("Ligand Count by Nuclear Receptor Family")
        axes[1].set_xlabel("Ligand count")
        axes[1].set_ylabel("NR family")

        plt.tight_layout()
        plt.savefig(FIGURES_DIR / "composition_summary.png", bbox_inches="tight")
        plt.show()
        """
    ),
    code(
        """
        tsne_features = analysis_df[
            ["MolWt", "LogP", "TPSA", "HBD", "HBA", "RotatableBonds", "QED", "RingCount", "AromaticRings"]
        ].dropna()

        tsne_index = tsne_features.index
        scaled_features = StandardScaler().fit_transform(tsne_features)
        tsne = TSNE(n_components=2, random_state=42, init="pca", learning_rate="auto", perplexity=30)
        tsne_embedding = tsne.fit_transform(scaled_features)

        tsne_df = analysis_df.loc[tsne_index, ["NR_Family", "Receptor_Group"]].copy()
        tsne_df["TSNE_1"] = tsne_embedding[:, 0]
        tsne_df["TSNE_2"] = tsne_embedding[:, 1]

        fig, axes = plt.subplots(1, 2, figsize=(20, 8))

        sns.scatterplot(
            data=tsne_df,
            x="TSNE_1",
            y="TSNE_2",
            hue="NR_Family",
            palette=family_palette,
            alpha=0.8,
            s=45,
            ax=axes[0],
        )
        axes[0].set_title("Chemical Space of Ligands Colored by NR Family")
        axes[0].legend(title="NR Family", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=True)

        sns.scatterplot(
            data=tsne_df,
            x="TSNE_1",
            y="TSNE_2",
            hue="Receptor_Group",
            palette=receptor_palette,
            alpha=0.8,
            s=45,
            ax=axes[1],
        )
        axes[1].set_title("Chemical Space of Ligands Colored by Top Receptors")
        axes[1].legend(title="Receptor", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=True)

        plt.tight_layout()
        plt.savefig(FIGURES_DIR / "chemical_space_tsne.png", bbox_inches="tight")
        plt.show()
        """
    ),
    md(
        """
        ## Final export summary

        The repo-level outputs are written to:
        - `data/processed/ligands_with_descriptors.csv`
        - `data/processed/receptor_summary.csv`
        - `data/processed/family_summary.csv`
        - `data/processed/scaffold_summary.csv`
        - `figures/property_overview.png`
        - `figures/composition_summary.png`
        - `figures/chemical_space_tsne.png`
        """
    ),
]


notebook = nbf.v4.new_notebook()
notebook["cells"] = cells
notebook["metadata"] = {
    "kernelspec": {
        "display_name": "Python 3",
        "language": "python",
        "name": "python3",
    },
    "language_info": {
        "name": "python",
        "version": "3.11",
    },
}

NOTEBOOK_PATH.parent.mkdir(parents=True, exist_ok=True)
nbf.write(notebook, NOTEBOOK_PATH)
print(f"Wrote notebook to {NOTEBOOK_PATH}")
