import argparse
import os
from pathlib import Path
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw, rdMolDescriptors


def calculate_exact_mass(smiles: str) -> float | None:
    """Calculates the exact molecular weight from a SMILES string."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        return rdMolDescriptors.CalcExactMolWt(mol) if mol else None
    except Exception:
        return None


def process_csv(input_csv: str, label_column: str) -> pd.DataFrame:
    """Loads the CSV file and computes exact mass."""
    df = pd.read_csv(input_csv, encoding="utf-8-sig", sep=",", decimal=",", quotechar='"')

    # Clean SMILES column
    df["SMILES"] = df["SMILES"].astype(str).str.strip()
    df = df[df["SMILES"].notna() & df["SMILES"].str.len() > 0]

    # Ensure label is float
    df[label_column] = df[label_column].astype(float)

    # Calculate exact mass and filter invalid
    df["ExactMass"] = df["SMILES"].apply(calculate_exact_mass)
    df = df[df["ExactMass"].notna()]
    df["ExactMassRounded"] = df["ExactMass"].round(4)

    return df


def generate_html(df: pd.DataFrame, label_column: str, output_html: str, img_dir: Path) -> None:
    """Generates an HTML file with RDKit molecule images grouped by label and exact mass."""
    img_dir.mkdir(exist_ok=True)
    sections = []

    for label_val in sorted(df[label_column].unique()):
        df_label = df[df[label_column] == label_val]
        added_label_header = False

        for mass_val in sorted(df_label["ExactMassRounded"].unique()):
            df_mass = df_label[df_label["ExactMassRounded"] == mass_val]

            # Skip singletons
            if len(df_mass) <= 1:
                continue

            if not added_label_header:
                sections.append(f"<h2>{label_column} = {label_val:.2f}</h2>")
                added_label_header = True

            smiles_list = list(df_mass["SMILES"].unique())
            mols = [Chem.MolFromSmiles(s) for s in smiles_list if Chem.MolFromSmiles(s)]

            img_tag = "<p><i>No valid SMILES</i></p>"
            if mols:
                img = Draw.MolsToGridImage(
                    mols,
                    molsPerRow=5,
                    subImgSize=(300, 300),
                    legends=smiles_list
                )
                img_path = img_dir / f"{label_column}_{label_val:.2f}_mass_{mass_val:.4f}.png"
                img.save(str(img_path))
                img_tag = f'<img src="{img_path.as_posix()}" alt="Mass {mass_val:.4f}">'

            sections.append(f"""
                <h3>Mass = {mass_val:.4f} (n = {len(df_mass)})</h3>
                {img_tag}
                <hr>
            """)

    # Save HTML
    with open(output_html, "w", encoding="utf-8") as f:
        f.write(f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>Grouped Viewer by {label_column}</title>
</head>
<body>
  <h1>Grouped Molecules by {label_column} and Exact Mass</h1>
  {''.join(sections)}
</body>
</html>
""")

    print(f"\n✅ Done! Open the HTML file: {output_html}")
    print(f"📁 Images saved in: {img_dir}/")


def main():
    parser = argparse.ArgumentParser(
        description="Process a CSV of molecules, calculate exact mass, group by a label, and generate an HTML visualizer."
    )
    parser.add_argument(
        "csv_path",
        type=str,
        help="Path to input CSV file containing SMILES and a label column."
    )
    parser.add_argument(
        "label_column",
        type=str,
        help="Column name in CSV that contains the grouping parameter (e.g., 'logP')."
    )
    parser.add_argument(
        "--output_html",
        type=str,
        default="grouped_viewer_with_exact_mass.html",
        help="Output HTML file name."
    )
    parser.add_argument(
        "--img_dir",
        type=str,
        default="img",
        help="Directory where molecule images will be saved."
    )

    args = parser.parse_args()

    df = process_csv(args.csv_path, args.label_column)
    generate_html(df, args.label_column, args.output_html, Path(args.img_dir))


if __name__ == "__main__":
    main()
