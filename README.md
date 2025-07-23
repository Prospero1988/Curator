# curator.py

**Deduplicate a CSV containing SMILES strings by chemical identity** ߚ

## Description
`curator.py` helps you clean up your chemical datasets by removing duplicate molecules based on their structural identity. It leverages robust cheminformatics tools (Open Babel, RDKit, MolVS) to ensure accurate deduplication, handling stereochemistry, tautomers, and even fuzzy similarity filters.

---

## Features
- **InChIKey identity**: Prefers Open Babel, falls back to RDKit.
- **Stereo merge**: Use `--nostereo` to merge enantiomers/E,Z isomers.
- **Tautomer handling**: `--tauto [off|weak|strong|molvs]`.
- **Fuzzy filtering**: `--tanimoto-thres THRESH` to reject similar molecules.
- **Parallel execution**: `--njobs N` to speed up processing.
- **Reject logging**: Captures `Duplicate`, `ParseError`, `IndexError`, `FPError`, `Tanimoto`.

---

## Installation

### Option 1: pip

1. Clone the repo or download `curator.py`.
2. Ensure you have Python 3.7+ installed.
3. Install dependencies:
   ```bash
   pip install rdkit-pypi tqdm molvs
   ```
4. (Optional) Install Open Babel for preferred InChIKey/tautomer support:
   ```bash
   sudo apt-get install openbabel
   ```

### Option 2: Conda

If you prefer Conda, use the provided `environment.yml` to create an environment named `data_curator`:

```bash
conda env create -f environment.yml -n data_curator
conda activate data_curator
```

---

## Usage

```bash
./curator.py input.csv SMILES_COLUMN [options]
```

### Positional arguments
- `input.csv` &nbsp;&nbsp;Path to your input CSV file.
- `SMILES_COLUMN` &nbsp;1-based index of the SMILES column.

### Options
| Flag                        | Description                                                                                       |
|-----------------------------|---------------------------------------------------------------------------------------------------|
| `--nostereo`                | Merge enantiomers / E,Z stereoisomers.                                                           |
| `--tauto {off,weak,strong,molvs}` | Tautomer handling:<br>• `off`: none<br>• `weak`: Open Babel<br>• `strong`: RDKit<br>• `molvs`: MolVS |
| `--njobs N`                 | Number of worker processes (default: 1).                                                         |
| `--tanimoto-thres THRESH`   | Reject rows with Tanimoto similarity ≥ THRESH (requires RDKit).                                   |
| `--quiet`                   | Silence RDKit log spam.                                                                          |

---

## Examples

- **Basic deduplication**:
  ```bash
  ./curator.py molecules.csv 2
  ```
- **Merge stereoisomers & tautomers with MolVS**:
  ```bash
  ./curator.py molecules.csv 2 --nostereo --tauto molvs
  ```
- **Fuzzy filter similar molecules**:
  ```bash
  ./curator.py molecules.csv 1 --tanimoto-thres 0.95 --njobs 4
  ```
- **Full combo**:  
  ```bash
  python curator.py "molecules.csv 2 --nostereo --tauto "weak" --njobs 8 --tanimoto-thres 0.9
  ```  
---

## Output

- `input_singled.csv` – deduplicated dataset.
- `input_rejects.csv` – rows rejected with reason codes.

---

## License
MIT © 2025

---

*Enjoy cleaner datasets! And yes, your boss will be impressed.* ߘ
