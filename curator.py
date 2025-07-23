#!/usr/bin/env python3
"""
curator.py
----------------------

Curates a CSV file containing SMILES strings by removing structural duplicates.

• Chemical identity = InChIKey (Open Babel / RDKit).
• Stereo can be merged (--nostereo).
• Tautomer handling:
    off       - keep raw SMILES
    weak      - Open Babel 3.x  (--tautomer canonical)
    strong    - RDKit TautomerEnumerator
    molvs     - MolVS Standardizer (open-source)
• Parallel execution via --njobs.
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
# Optional deps
# --------------------------------------------------------------------- #
try:
    from tqdm import tqdm
    _HAVE_TQDM = True
except ImportError:
    _HAVE_TQDM = False

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem.MolStandardize import rdMolStandardize
    _HAVE_RDKIT = True
except ImportError:
    _HAVE_RDKIT = False
    Chem = None  # type: ignore

try:
    from molvs import Standardizer
    _HAVE_MOLVS = True
    _MOLVS_STD = Standardizer()
except ImportError:
    _HAVE_MOLVS = False
    _MOLVS_STD = None  # type: ignore

_OBABEL = shutil.which("obabel")

# --------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------- #
def _silence_rdkit() -> None:
    """Mute *all* RDKit logs (Python and C++ layers)."""
    if _HAVE_RDKIT:
        from rdkit import rdBase
        RDLogger.DisableLog("rdApp.*")
        rdBase.DisableLog("rdApp.*")


def _canonical_tautomer_rdkit(smiles: str) -> Optional[str]:
    """RDKit canonical tautomer (silent, no kekulization panic)."""
    if not _HAVE_RDKIT:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return None
        enum = rdMolStandardize.TautomerEnumerator()
        mol = enum.Canonicalize(mol)
        Chem.SanitizeMol(
            mol,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
            ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE,
        )
        Chem.RemoveStereochemistry(mol)
        return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
    except Exception:
        return None


def _canonical_tautomer_molvs(smiles: str) -> Optional[str]:
    """MolVS canonical tautomer (requires RDKit + molvs)."""
    if not (_HAVE_RDKIT and _HAVE_MOLVS):
        return None
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return None
        mol = _MOLVS_STD.canonicalize_tautomer(mol)  # type: ignore
        Chem.RemoveStereochemistry(mol)
        return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
    except Exception:
        return None


def _canonical_tautomer_obabel(smiles: str) -> Optional[str]:
    if _OBABEL is None:
        return None
    cmd = [_OBABEL, f"-:{smiles}", "--tautomer", "canonical", "-ocansmi"]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=5)
        return proc.stdout.strip() if proc.returncode == 0 else None
    except Exception:
        return None


def _inchikey_obabel(smiles: str, no_stereo: bool, keep_fixedh: bool) -> Optional[str]:
    if _OBABEL is None:
        return None
    cmd = [_OBABEL, f"-:{smiles}", "-oinchikey"]
    if no_stereo:
        cmd += ["-xT", "/nostereo"]
    if keep_fixedh:
        cmd += ["-xT", "/fixedH"]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=5)
        return proc.stdout.strip() if proc.returncode == 0 else None
    except Exception:
        return None


def _inchikey_rdkit(smiles: str, no_stereo: bool, keep_fixedh: bool) -> Optional[str]:
    if not _HAVE_RDKIT:
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


# --------------------------------------------------------------------- #
# Worker
# --------------------------------------------------------------------- #
def _process_row(payload: Tuple[int, List[str], int, bool, str]) -> Tuple[str, List[str]]:
    """Executed in subprocess; returns (identity_key | reason, row)."""
    _silence_rdkit()  # mute RDKit in each worker

    lineno, row, col_idx, nostereo, tauto_mode = payload
    try:
        smiles = row[col_idx].strip()
    except IndexError:
        return "IndexError", row

    # Tautomer canonicalization
    if tauto_mode == "strong":
        can_smi = _canonical_tautomer_rdkit(smiles)
    elif tauto_mode == "weak":
        can_smi = _canonical_tautomer_obabel(smiles)
    elif tauto_mode == "molvs":
        can_smi = _canonical_tautomer_molvs(smiles)
    else:  # off
        can_smi = None
    if can_smi:
        smiles = can_smi

    keep_fixedh = tauto_mode == "off"
    cid = _inchikey_obabel(smiles, nostereo, keep_fixedh)
    if cid is None:
        cid = _inchikey_rdkit(smiles, nostereo, keep_fixedh)
    return (cid if cid else "ParseError", row)


# --------------------------------------------------------------------- #
# Deduplicator
# --------------------------------------------------------------------- #
def deduplicate(csv_path: str, smiles_col: int, nostereo: bool,
                tauto_mode: str, n_jobs: int) -> None:
    col_idx = smiles_col - 1
    total = _count_lines(csv_path) - 1

    out_path = _new_name(csv_path, "_singled")
    rej_path = _new_name(csv_path, "_rejects")

    with open(csv_path, newline="", encoding="utf-8") as fh:
        rdr = csv.reader(fh)
        header = next(rdr)
        payloads = [(ln, row, col_idx, nostereo, tauto_mode)
                    for ln, row in enumerate(rdr, start=2)]

    executor_cls = ProcessPoolExecutor if n_jobs > 1 else None
    kept, rejects = [header], []
    seen: OrderedDict[str, int] = OrderedDict()

    def _progress(it: Iterator):
        if _HAVE_TQDM:
            yield from tqdm(it, total=total, unit=" rows", ncols=80)
        else:
            yield from it

    iterator: Iterator[Tuple[str, List[str]]]
    if executor_cls:
        with executor_cls(max_workers=n_jobs) as ex:
            iterator = _progress(ex.map(_process_row, payloads, chunksize=500))
            for cid, row in iterator:
                if cid in ("ParseError", "IndexError"):
                    rejects.append([cid] + row)
                elif cid not in seen:
                    seen[cid] = 1
                    kept.append(row)
    else:
        for cid, row in _progress(map(_process_row, payloads)):
            if cid in ("ParseError", "IndexError"):
                rejects.append([cid] + row)
            elif cid not in seen:
                seen[cid] = 1
                kept.append(row)

    # save
    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        csv.writer(fh).writerows(kept)
    if rejects:
        with open(rej_path, "w", newline="", encoding="utf-8") as fh:
            csv.writer(fh).writerows([["REASON", *header], *rejects])

    print(f"\n✓ {out_path} — {len(kept) - 1} rows kept")
    if rejects:
        print(f"⚠ {len(rejects)} rejected rows → {rej_path}")


# --------------------------------------------------------------------- #
# Utils
# --------------------------------------------------------------------- #
def _count_lines(path: str) -> int:
    with open(path, "rb") as fh:
        return sum(1 for _ in fh)


def _new_name(path: str, suf: str) -> str:
    stem, ext = os.path.splitext(path)
    return f"{stem}{suf}{ext}"


# --------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------- #
def main() -> None:
    parser = argparse.ArgumentParser(
        description="Deduplicate SMILES by chemical identity (parallel)."
    )
    parser.add_argument("csv_file")
    parser.add_argument("smiles_col", type=int,
                        help="SMILES column (1-based)")
    parser.add_argument("--nostereo", action="store_true",
                        help="merge enantiomers / E,Z")
    parser.add_argument("--tauto", choices=["off", "weak", "strong", "molvs"],
                        default="weak",
                        help=("tautomer handling: "
                              "off=none, weak=OpenBabel, "
                              "strong=RDKit, molvs=MolVS"))
    parser.add_argument("--njobs", type=int, default=1,
                        help="parallel worker processes (default 1)")
    parser.add_argument("--quiet", action="store_true",
                        help="silence RDKit logs")

    args = parser.parse_args()
    if args.quiet:
        _silence_rdkit()

    logging.basicConfig(level=logging.ERROR if args.quiet else logging.INFO,
                        format="%(levelname)s: %(message)s")

    if args.tauto == "molvs" and not _HAVE_MOLVS:
        parser.error("MolVS not installed: pip install molvs")

    deduplicate(args.csv_file, args.smiles_col,
                nostereo=args.nostereo,
                tauto_mode=args.tauto,
                n_jobs=max(1, args.njobs))


if __name__ == "__main__":
    main()
