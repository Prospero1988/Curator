#!/usr/bin/env python3
"""
curator.py
==========

Deduplicate a CSV containing SMILES strings by chemical identity.

Features
--------
* InChIKey identity (Open Babel preferred, RDKit fallback).
* Optional stereo merge  : ``--nostereo``.
* Optional tautomer mode : ``--tauto off|weak|strong|molvs``.
* Optional fuzzy filter  : ``--tanimoto-thres 0.95``.
* Parallel execution     : ``--njobs``.
* Reject log             : Duplicate, ParseError, IndexError, FPError, Tanimoto.
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import shutil
import subprocess
import sys
from collections import OrderedDict
from concurrent.futures import ProcessPoolExecutor
from typing import Iterator, List, Optional, Tuple

# --------------------------------------------------------------------- #
# Optional dependencies
# --------------------------------------------------------------------- #
try:
    from tqdm import tqdm

    HAVE_TQDM = True
except ImportError:  # pragma: no cover
    HAVE_TQDM = False

try:
    from rdkit import Chem, RDLogger, DataStructs
    from rdkit.Chem import AllChem
    from rdkit.Chem.MolStandardize import rdMolStandardize
    from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator

    HAVE_RDKIT = True
except ImportError:  # pragma: no cover
    HAVE_RDKIT = False
    Chem = None  # type: ignore

if HAVE_RDKIT:
    MORGAN_GEN = GetMorganGenerator(
        radius=2,
        fpSize=2048,
        includeChirality=False,
    )
else:
    MORGAN_GEN = None  # type: ignore

try:
    from molvs import Standardizer

    HAVE_MOLVS = True
    MOLVS_STD = Standardizer()
except ImportError:  # pragma: no cover
    HAVE_MOLVS = False
    MOLVS_STD = None  # type: ignore

OBABEL_BIN = shutil.which("obabel")

# --------------------------------------------------------------------- #
# Utility helpers
# --------------------------------------------------------------------- #
def silence_rdkit() -> None:
    """Mute all RDKit logs (Python + C++ layers)."""
    if HAVE_RDKIT:
        from rdkit import rdBase

        RDLogger.DisableLog("rdApp.*")
        rdBase.DisableLog("rdApp.*")


def canonical_tautomer_rdkit(smiles: str) -> Optional[str]:
    """Return RDKit canonical tautomer (silent)."""
    if not HAVE_RDKIT:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return None
        mol = rdMolStandardize.TautomerEnumerator().Canonicalize(mol)
        Chem.SanitizeMol(
            mol,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
            ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE,
        )
        Chem.RemoveStereochemistry(mol)
        return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
    except Exception:
        return None


def canonical_tautomer_molvs(smiles: str) -> Optional[str]:
    """Return MolVS canonical tautomer."""
    if not (HAVE_RDKIT and HAVE_MOLVS):
        return None
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return None
        mol = MOLVS_STD.canonicalize_tautomer(mol)  # type: ignore
        Chem.RemoveStereochemistry(mol)
        return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
    except Exception:
        return None


def canonical_tautomer_obabel(smiles: str) -> Optional[str]:
    """Return Open Babel canonical tautomer."""
    if OBABEL_BIN is None:
        return None
    cmd = [OBABEL_BIN, f"-:{smiles}", "--tautomer", "canonical", "-ocansmi"]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=5)
        return proc.stdout.strip() if proc.returncode == 0 else None
    except Exception:
        return None


def inchikey_obabel(smiles: str, no_stereo: bool, keep_fixedh: bool) -> Optional[str]:
    """Convert SMILES to InChIKey via Open Babel."""
    if OBABEL_BIN is None:
        return None
    cmd = [OBABEL_BIN, f"-:{smiles}", "-oinchikey"]
    if no_stereo:
        cmd += ["-xT", "/nostereo"]
    if keep_fixedh:
        cmd += ["-xT", "/fixedH"]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=5)
        return proc.stdout.strip() if proc.returncode == 0 else None
    except Exception:
        return None


def inchikey_rdkit(smiles: str, no_stereo: bool, keep_fixedh: bool) -> Optional[str]:
    """Convert SMILES to InChIKey via RDKit."""
    if not HAVE_RDKIT:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        if no_stereo:
            Chem.RemoveStereochemistry(mol)
        opts = "-FixedH" if keep_fixedh else ""
        return Chem.InchiToInchiKey(Chem.MolToInchi(mol, options=opts))
    except Exception:
        return None


def fingerprint_ecfp4(smiles: str):
    """Return an ECFP-4 fingerprint (2048 bits)."""
    if not HAVE_RDKIT or MORGAN_GEN is None:
        return None
    mol = Chem.MolFromSmiles(smiles)
    return MORGAN_GEN.GetFingerprint(mol) if mol else None


def tanimoto(fp1, fp2) -> float:
    """Return Tanimoto similarity between two fingerprints."""
    return DataStructs.TanimotoSimilarity(fp1, fp2)  # type: ignore


# --------------------------------------------------------------------- #
# Worker (runs in a subprocess)
# --------------------------------------------------------------------- #
def _process_row(
    payload: Tuple[int, List[str], int, bool, str],
) -> Tuple[str, List[str]]:
    """Return (identity_key | reason, row)."""
    silence_rdkit()
    _, row, col_idx, nostereo, tauto_mode = payload

    try:
        smiles = row[col_idx].strip()
    except IndexError:
        return "IndexError", row

    # tautomer canonicalisation
    if tauto_mode == "strong":
        smiles = canonical_tautomer_rdkit(smiles) or smiles
    elif tauto_mode == "weak":
        smiles = canonical_tautomer_obabel(smiles) or smiles
    elif tauto_mode == "molvs":
        smiles = canonical_tautomer_molvs(smiles) or smiles

    # InChIKey
    keep_fixedh = tauto_mode == "off"
    key = inchikey_obabel(smiles, nostereo, keep_fixedh) or inchikey_rdkit(
        smiles, nostereo, keep_fixedh
    )
    return (key if key else "ParseError", row)


# --------------------------------------------------------------------- #
# Deduplication core
# --------------------------------------------------------------------- #
def deduplicate(
    csv_path: str,
    smiles_col: int,
    nostereo: bool,
    tauto_mode: str,
    n_jobs: int,
    tanimoto_thr: Optional[float],
) -> None:
    col_idx = smiles_col - 1
    total_rows = _count_lines(csv_path) - 1

    out_path = _out_name(csv_path, "_singled")
    rej_path = _out_name(csv_path, "_rejects")

    with open(csv_path, newline="", encoding="utf-8") as fh:
        rdr = csv.reader(fh)
        header = next(rdr)
        payloads = [
            (idx, row, col_idx, nostereo, tauto_mode)
            for idx, row in enumerate(rdr, start=2)
        ]

    executor_cls = ProcessPoolExecutor if n_jobs > 1 else None
    kept: List[List[str]] = [header]
    rejects: List[List[str]] = []
    seen: OrderedDict[str, None] = OrderedDict()

    def _pbar(it: Iterator):
        return tqdm(it, total=total_rows, unit=" rows", ncols=80) if HAVE_TQDM else it

    # pass 1 – strict InChIKey
    if executor_cls:
        with executor_cls(max_workers=n_jobs) as ex:
            for key, row in _pbar(ex.map(_process_row, payloads, chunksize=500)):
                if key in ("ParseError", "IndexError"):
                    rejects.append([key] + row)
                elif key in seen:
                    rejects.append(["Duplicate"] + row)
                else:
                    seen[key] = None
                    kept.append(row)
    else:
        for key, row in _pbar(map(_process_row, payloads)):
            if key in ("ParseError", "IndexError"):
                rejects.append([key] + row)
            elif key in seen:
                rejects.append(["Duplicate"] + row)
            else:
                seen[key] = None
                kept.append(row)

    # pass 2 – fuzzy duplicates (Tanimoto)
    if tanimoto_thr is not None and HAVE_RDKIT:
        fp_cache, final_kept = [], [header]
        for row in kept[1:]:
            smiles = row[col_idx].strip()
            fp = fingerprint_ecfp4(smiles)
            if fp is None:
                rejects.append(["FPError"] + row)
                continue
            if any(tanimoto(fp, old) >= tanimoto_thr for old in fp_cache):
                rejects.append(["Tanimoto"] + row)
            else:
                fp_cache.append(fp)
                final_kept.append(row)
        kept = final_kept

    # write outputs
    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        csv.writer(fh).writerows(kept)

    if rejects:
        with open(rej_path, "w", newline="", encoding="utf-8") as fh:
            csv.writer(fh).writerows([["REASON", *header], *rejects])

    print(f"\n✓ {out_path} — {len(kept) - 1} rows kept")
    if rejects:
        print(f"⚠ {len(rejects)} rejected rows → {rej_path}")


# --------------------------------------------------------------------- #
# Tiny helpers
# --------------------------------------------------------------------- #
def _count_lines(path: str) -> int:
    with open(path, "rb") as handle:
        return sum(1 for _ in handle)


def _out_name(path: str, suffix: str) -> str:
    stem, ext = os.path.splitext(path)
    return f"{stem}{suffix}{ext}"


# --------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------- #
def main() -> None:
    parser = argparse.ArgumentParser(
        description="Deduplicate SMILES by chemical identity (parallel)."
    )
    parser.add_argument("csv_file", help="path to input CSV")
    parser.add_argument("smiles_col", type=int, help="SMILES column (1-based)")
    parser.add_argument(
        "--nostereo",
        action="store_true",
        help="merge enantiomers / E,Z stereoisomers",
    )
    parser.add_argument(
        "--tauto",
        choices=["off", "weak", "strong", "molvs"],
        default="weak",
        help=(
            "tautomer handling: "
            "off = none, weak = Open Babel, strong = RDKit, molvs = MolVS"
        ),
    )
    parser.add_argument(
        "--njobs",
        type=int,
        default=1,
        help="number of worker processes (default 1)",
    )
    parser.add_argument(
        "--tanimoto-thres",
        type=float,
        metavar="THRESH",
        help="enable fuzzy filter; reject rows with Tanimoto "
        "similarity ≥ THRESH (requires RDKit)",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="silence RDKit log spam",
    )
    args = parser.parse_args()

    if args.quiet:
        silence_rdkit()

    logging.basicConfig(
        level=logging.ERROR if args.quiet else logging.INFO,
        format="%(levelname)s: %(message)s",
    )

    if args.tauto == "molvs" and not HAVE_MOLVS:
        parser.error("MolVS not installed: pip install molvs")

    deduplicate(
        args.csv_file,
        args.smiles_col,
        nostereo=args.nostereo,
        tauto_mode=args.tauto,
        n_jobs=max(1, args.njobs),
        tanimoto_thr=args.tanimoto_thres,
    )


if __name__ == "__main__":
    main()
