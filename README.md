# SMILES Curator (`curator.py`)

Curates a CSV file containing SMILES strings by removing structural duplicates before QSPR/QSAR modeling.

## Key Features

- **Chemical identity deduplication** based on InChIKey.
- **Stereochemistry merging**: `--nostereo` flag disables R/S, E/Z distinction.
- **Tautomer canonicalization** (optional):
  - `off` – raw SMILES (no normalization)
  - `weak` – Open Babel 3.x (`--tautomer canonical`)
  - `strong` – RDKit TautomerEnumerator
  - `molvs` – [MolVS](https://github.com/mcs07/MolVS) canonicalizer (requires RDKit)
- **Parallel processing**: multithreading via `--njobs`
- **Graceful fallback**: uses RDKit if Open Babel is unavailable or fails.
- **Verbose / silent mode**: enable progress bar with `tqdm`, or suppress logs with `--quiet`

## Usage

```bash
python curator.py data.csv 2 --nostereo --tauto molvs --njobs 4
```

- `data.csv` – input file with SMILES
- `2` – column number (1-based) with SMILES strings
- `--nostereo` – merge stereoisomers
- `--tauto molvs` – normalize tautomers using MolVS
- `--njobs 4` – use 4 parallel workers

## Output Files

- `data_singled.csv` – curated, deduplicated data
- `data_rejects.csv` – rows with invalid or unparsable SMILES (if any)

## Dependencies

- [RDKit](https://www.rdkit.org/)
- [Open Babel](http://openbabel.org/)
- [MolVS](https://github.com/mcs07/MolVS) *(optional, for `--tauto molvs`)*
- `tqdm` *(optional, for progress bar)*

## Example

Given this input:

```csv
ID,SMILES
1,C1=CC=CC=C1
2,c1ccccc1
3,C1=CC=CN=C1
4,C1=CC=CN=C1
```

Running:

```bash
python curator.py input.csv 2 --nostereo --tauto strong
```

Produces:

- `input_singled.csv` with unique structures.
- `input_rejects.csv` (if any errors or duplicates)

## License

MIT License (your choice). See `LICENSE` file if applicable.

## Author

Arkadiusz Leniak © 2025 – arkadiusz.leniak@gmail.com 
