# maccs_atom_ids

A Python script that computes MACCS keys for molecules in a CSV file and, instead of returning a binary fingerprint, reports **which atom indices** matched each of the 166 bits.

## Output format

The output is a rectangular CSV with 168 columns:

| Column | Content |
|---|---|
| `ID` | Molecule identifier |
| `SMILES` | Canonical SMILES (RDKit-normalised) |
| `MACCS_1` … `MACCS_166` | `0` if bit is off; semicolon-separated atom indices (e.g. `0;1;2;3;4;5`) if bit is on |

Multi-atom matches (e.g. aromatic rings) are packed into a single cell as `"0;1;2;3;4;5"`, keeping the file rectangular and directly importable into R or pandas.

## Requirements

```
rdkit
pandas
```

Install with:

```bash
pip install rdkit pandas
# or via conda:
conda install -c conda-forge rdkit pandas
```

## Usage

```bash
# Minimal — auto-detects ID column, writes <input>_maccs.csv
python maccs_atom_ids.py molecules.csv

# Explicit column names (all lookups are case-insensitive)
python maccs_atom_ids.py molecules.csv --smiles-col Structure --id-col CompoundID

# Custom output path
python maccs_atom_ids.py molecules.csv -o results/maccs_output.csv

# Write to stdout (e.g. for piping)
python maccs_atom_ids.py molecules.csv -s | head -n 3
```

### Input CSV format

The only requirement is a column containing SMILES strings. Any other columns are ignored.

```
Name,SMILES,MW
aspirin,CC(=O)Oc1ccccc1C(=O)O,180.2
benzene,c1ccccc1,78.1
```

Column name resolution order for **ID**: explicit `--id-col` → first non-SMILES column → auto-generated `mol_1, mol_2, …`

## Options

| Flag | Default | Description |
|---|---|---|
| `--smiles-col` | `SMILES` | Name of the SMILES column (case-insensitive) |
| `--id-col` | first other column | Name of the ID column (case-insensitive) |
| `-o / --output` | `<input>_maccs.csv` | Output file path |
| `-s / --stdout` | off | Print CSV to stdout instead of a file |

## Design notes

**Atom-index canonicalisation** — input SMILES atom order is arbitrary and would produce different indices for the same molecule written differently. The script renumbers all atoms into RDKit canonical order with `RenumberAtoms` before any substructure matching, so the reported atom IDs depend only on the molecular graph, not the SMILES string.

**Minimum match count** — RDKit's MACCS SMARTS patterns include a minimum number of required matches per bit (e.g. bit 160 requires at least two aromatic atoms). This threshold is respected: a bit is only considered ON when the pattern matches at least the required number of times.

**Reconciliation pass** — a handful of MACCS bits (e.g. bit 125, isotope) are computed at the molecule level rather than via a substructure match. When such a bit disagrees between our SMARTS walk and RDKit's `GenMACCSKeys`, the official fingerprint wins and the cell is marked `1` (bit on, but no atom mapping available).

**Warnings** — unparseable SMILES rows are skipped with a `[WARN]` message to stderr; the rest of the file is processed normally.

## Reading the output in R

```r
library(tidyverse)

df <- read_csv("molecules_maccs.csv")

# Bits with atom IDs are character columns; 0 is numeric — unify to character
df <- df |> mutate(across(starts_with("MACCS_"), as.character))

# Example: which molecules have bit 166 (any ring) set?
df |> filter(MACCS_166 != "0") |> select(ID, SMILES, MACCS_166)
```

## Reading the output in Python

```python
import pandas as pd

df = pd.read_csv("molecules_maccs.csv", dtype=str)

# Parse atom IDs for a specific bit
def parse_bit(cell):
    if cell == "0":
        return []
    return [int(x) for x in cell.split(";")]

df["MACCS_162_atoms"] = df["MACCS_162"].apply(parse_bit)
```
