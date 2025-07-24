# curator.py

**Curator** is a Python script designed for preprocessing CSV files containing SMILES strings (molecular structures) in cheminformatics workflows. It helps clean and standardize your data before applying machine learning or statistical models, ensuring that you’re not feeding in duplicates, stereoisomers, outliers, or corrupted entries.

This script is especially useful when working with large datasets from chemical databases like ChEMBL, PubChem, or ZINC.

---

## ߔ What Does It Do?

### ✅ 1. Deduplication
Rows are compared by their **InChIKey** (a unique identifier for molecules). You can choose how strict the deduplication should be:
- **Default**: Considers full molecular identity.
- `--nostereo`: Ignores stereochemistry (treats stereoisomers as the same compound).
- `--tauto`: Handles tautomers (different forms of the same molecule) using one of the following engines:
  - `off`: No canonicalization.
  - `weak`: Use Open Babel.
  - `strong`: Use RDKit.
  - `molvs`: Use the MolVS library.

### ߧ 2. Fuzzy Duplicate Removal (Optional)
Use ECFP4 fingerprints to catch similar-but-not-identical compounds.
```bash
--tanimoto-thres 0.95
```
Any compound with a Tanimoto similarity ≥ threshold to another will be discarded.

### ߓ 3. Outlier Detection (Optional)
If your CSV file contains a numeric property (e.g., logP, bioactivity, etc.) in the **third column**, this script can automatically:
- Parse the numbers (supports US `1,234.56` and EU `1.234,56` formats).
- Compute **mean ± k×standard deviation**.
- Filter out statistical outliers.
- Optionally display histograms before/after filtering.

Enable this using:
```bash
--outliers 1.5
```

### ⚙️ 4. Parallel Execution
Speed up the process using:
```bash
--njobs 8
```

---

## ߒ Why Use It?

Chemical datasets are **messy**:
- You’ll often get duplicated compounds under different names or notations.
- Tautomers and stereoisomers can lead to data leakage in ML models.
- Some numeric fields come in European formats or have junk text like `"~1.5 ± 0.2"`.

**Curator** handles all that. It’s fast, parallelized, and outputs:
- A cleaned file: `yourfile_singled.csv`
- A detailed reject log: `yourfile_rejects.csv`, tagging each issue:
  - `Duplicate`, `ParseError`, `IndexError`, `FPError`, `Tanimoto`, `Outlier`, `ValueError`

---

## ߓ Installation & Requirements

### Required:
- Python 3.7+
- [Open Babel](http://openbabel.org/)
- [RDKit](https://www.rdkit.org/)

### Optional (recommended):
- `molvs` (tautomer standardization)
- `tqdm` (progress bars)
- `matplotlib` (histograms)

Install missing packages with:
```bash
pip install molvs tqdm matplotlib
```

### Option 2: Conda

If you prefer Conda, use the provided environment.yml to create an environment named data_curator:

```bash
conda env create -f environment.yml -n data_curator
conda activate data_curator
```

---

## ߧ Example Usage

```bash
python curator.py molecules.csv 2 --nostereo --tauto weak --tanimoto-thres 0.9 --outliers 2.0 --njobs 8
```

---

## ߓ Output Files

- ✅ `molecules_singled.csv` → final curated file
- ⚠️  `molecules_rejects.csv` → detailed log of removed/invalid entries

---

## ߧ Clean Code

- Follows [PEP 8](https://peps.python.org/pep-0008/)
- Readable, modular structure
- Compatible with both Unix and Windows systems

---

## ߧ Ideal For:

- Researchers using cheminformatics or QSAR modeling
- Data preprocessing in drug discovery pipelines
- ML practitioners cleaning molecular datasets

---

Happy curation! ߧߧ
