"""
maccs_atom_ids.py
Reads a CSV with a SMILES column, outputs MACCS key atom-ID mappings.

Usage:
    python maccs_atom_ids.py input.csv              # writes input_maccs.csv
    python maccs_atom_ids.py input.csv -o out.csv   # custom output path
    python maccs_atom_ids.py input.csv -s           # to stdout
    python maccs_atom_ids.py input.csv --id-col Name --smiles-col Structure
"""

import argparse
import csv
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import MACCSkeys


# ── Core logic (unchanged) ────────────────────────────────────────────────────

def maccs_atom_ids(mol: Chem.Mol) -> dict[int, str | int]:
    result = {}
    for bit in range(1, 167):
        entry = MACCSkeys.smartsPatts.get(bit)
        if entry is None or entry[0] == "?":
            result[bit] = 0
            continue
        smarts, min_count = entry
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            result[bit] = 0
            continue
        matches = mol.GetSubstructMatches(patt)
        if len(matches) < min_count:
            result[bit] = 0
            continue
        atom_ids = sorted({idx for match in matches for idx in match})
        result[bit] = ";".join(map(str, atom_ids)) if atom_ids else 0
    return result


def process_mol(mol_id: str, smiles: str) -> list | None:
    """Returns a fully populated row or None on parse failure."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"[WARN] Cannot parse SMILES for '{mol_id}': {smiles}", file=sys.stderr)
        return None

    ranks = list(Chem.CanonicalRankAtoms(mol))
    new_order = sorted(range(mol.GetNumAtoms()), key=lambda i: ranks[i])
    mol = Chem.RenumberAtoms(mol, new_order)
    canonical_smiles = Chem.MolToSmiles(mol)

    official_bv = MACCSkeys.GenMACCSKeys(mol)
    bit_data = maccs_atom_ids(mol)

    for bit in range(1, 167):
        official_on = bool(official_bv[bit])
        our_on = bit_data[bit] != 0
        if official_on != our_on:
            bit_data[bit] = "1" if official_on else 0

    return [mol_id, canonical_smiles] + [bit_data.get(i, 0) for i in range(1, 167)]


# ── CSV reader ────────────────────────────────────────────────────────────────

def find_column(columns: list[str], target: str) -> str:
    """Case-insensitive column search; raises ValueError if not found."""
    matches = [c for c in columns if c.strip().lower() == target.lower()]
    if not matches:
        raise ValueError(
            f"Column '{target}' not found. Available columns: {list(columns)}"
        )
    return matches[0]


def process_csv(
    input_path: Path,
    smiles_col: str = "SMILES",
    id_col: str | None = None,
    out_stream=None,
) -> None:
    df = pd.read_csv(input_path)

    # Resolve SMILES column
    smiles_col = find_column(df.columns.tolist(), smiles_col)

    # Resolve ID column: explicit > first non-SMILES column > auto-index
    if id_col:
        id_col = find_column(df.columns.tolist(), id_col)
        ids = df[id_col].astype(str)
    else:
        other_cols = [c for c in df.columns if c != smiles_col]
        if other_cols:
            id_col = other_cols[0]
            print(f"[INFO] Using '{id_col}' as ID column.", file=sys.stderr)
            ids = df[id_col].astype(str)
        else:
            ids = pd.Series([f"mol_{i+1}" for i in range(len(df))])

    header = ["ID", "SMILES"] + [f"MACCS_{i}" for i in range(1, 167)]
    writer = csv.writer(out_stream, quoting=csv.QUOTE_MINIMAL)
    writer.writerow(header)

    for mol_id, smiles in zip(ids, df[smiles_col]):
        row = process_mol(str(mol_id), str(smiles))
        if row:
            writer.writerow(row)


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="MACCS key atom-ID mapper from CSV.")
    parser.add_argument("input", type=Path, help="Input CSV file.")
    parser.add_argument("--smiles-col", default="SMILES", help="SMILES column name (case-insensitive). Default: SMILES")
    parser.add_argument("--id-col", default=None,    help="ID column name (case-insensitive). Default: first non-SMILES column.")
    parser.add_argument("-o", "--output", type=Path,  default=None, help="Output CSV path. Default: <input>_maccs.csv")
    parser.add_argument("-s", "--stdout", action="store_true", help="Write to stdout instead of a file.")
    args = parser.parse_args()

    if not args.input.exists():
        print(f"[ERROR] File not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    if args.stdout:
        process_csv(args.input, args.smiles_col, args.id_col, out_stream=sys.stdout)
    else:
        out_path = args.output or args.input.with_name(args.input.stem + "_maccs.csv")
        print(f"[INFO] Writing to {out_path}", file=sys.stderr)
        with open(out_path, "w", newline="") as f:
            process_csv(args.input, args.smiles_col, args.id_col, out_stream=f)


if __name__ == "__main__":
    main()
