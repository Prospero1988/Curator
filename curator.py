#!/usr/bin/env python3
"""
curator.py
==========

Curate a CSV containing SMILES strings by chemical identity and, optionally,
trim statistical outliers from a numeric column (default: column 3, index 2).

Features
--------
* Strict identity de-duplication by InChIKey (Open Babel preferred, RDKit fallback)
* Optional stereo merge      : ``--nostereo``
* Optional tautomer merge    : ``--tauto off|weak|strong|molvs``
* Optional fuzzy filter      : ``--tanimoto-thres 0.95`` (ECFP-4)
* Optional outlier trimming  : ``--outliers 1.5``  →  keep values μ ± K·o
* Locale-aware number parsing (handles “1 234,56”, “1,23”, etc.)
* Parallel execution with progress bar (tqdm) and auto-tuned chunk size
* Reject log: Duplicate, ParseError, IndexError, FPError, Tanimoto,
              ValueError, Outlier
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import re
import shutil
import subprocess
from collections import OrderedDict
from concurrent.futures import ProcessPoolExecutor
from statistics import mean, stdev
from typing import Iterator, List, Optional, Tuple

# --------------------------------------------------------------------------- #
# Optional dependencies
# --------------------------------------------------------------------------- #
try:
    from tqdm import tqdm

    HAVE_TQDM = True
except ImportError:       # pragma: no cover
    HAVE_TQDM = False

try:
    from rdkit import Chem, RDLogger, DataStructs
    from rdkit.Chem import AllChem            # noqa: F401 (side effects)
    from rdkit.Chem.MolStandardize import rdMolStandardize
    from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator

    HAVE_RDKIT = True
except ImportError:       # pragma: no cover
    HAVE_RDKIT = False
    Chem = None  # type: ignore

if HAVE_RDKIT:
    MORGAN_GEN = GetMorganGenerator(radius=2, fpSize=2048, includeChirality=False)
else:
    MORGAN_GEN = None  # type: ignore

try:
    from molvs import Standardizer

    HAVE_MOLVS = True
    MOLVS_STD = Standardizer()
except ImportError:       # pragma: no cover
    HAVE_MOLVS = False
    MOLVS_STD = None  # type: ignore

try:
    import matplotlib.pyplot as plt

    HAVE_MPL = True
except ImportError:       # pragma: no cover
    HAVE_MPL = False

OBABEL_BIN = shutil.which("obabel")

# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def silence_rdkit() -> None:
    """Mute RDKit log spam (Python and C++ layers)."""
    if HAVE_RDKIT:
        from rdkit import rdBase

        RDLogger.DisableLog("rdApp.*")
        rdBase.DisableLog("rdApp.*")


# ---------- locale-aware float coercion ---------- #
_DECIMAL_RE = re.compile(r"[0-9]+(?:[.,][0-9]+)*(?:[eE][-+]?[0-9]+)?")


def coerce_float(token: str) -> Optional[float]:
    """Convert *token* to float, accepting mixed EU/US formats."""
    s = token.replace("\xa0", " ").strip()
    if not s:
        return None
    # fast path
    try:
        return float(s)
    except ValueError:
        pass
    # single comma, no dot   → decimal comma
    if s.count(",") == 1 and s.count(".") == 0:
        try:
            return float(s.replace(",", "."))
        except ValueError:
            pass
    # grouping dots + comma  → remove dots, comma→dot
    if s.count(",") == 1 and s.count(".") >= 1:
        try:
            return float(s.replace(".", "").replace(",", "."))
        except ValueError:
            pass
    # fallback regex
    m = _DECIMAL_RE.search(s)
    if m:
        try:
            return float(m.group(0).replace(".", "").replace(",", "."))
        except ValueError:
            pass
    return None


# ---------- canonical tautomer helpers ---------- #
def _tautomer_rdkit(smiles: str) -> Optional[str]:
    if not HAVE_RDKIT:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return None
        mol = rdMolStandardize.TautomerEnumerator().Canonicalize(mol)
        Chem.SanitizeMol(
            mol,
            sanitizeOps=(Chem.SanitizeFlags.SANITIZE_ALL
                         ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE),
        )
        Chem.RemoveStereochemistry(mol)
        return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
    except Exception:
        return None


def _tautomer_molvs(smiles: str) -> Optional[str]:
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


def _tautomer_obabel(smiles: str) -> Optional[str]:
    if OBABEL_BIN is None:
        return None
    cmd = [OBABEL_BIN, f"-:{smiles}", "--tautomer", "canonical", "-ocansmi"]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=5)
        return proc.stdout.strip() if proc.returncode == 0 else None
    except Exception:
        return None


# ---------- InChIKey helpers ---------- #
def _inchikey_obabel(smiles: str, no_stereo: bool, keep_fixedh: bool) -> Optional[str]:
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


def _inchikey_rdkit(smiles: str, no_stereo: bool, keep_fixedh: bool) -> Optional[str]:
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


# ---------- fingerprints ---------- #
def _ecfp4(smiles: str):
    if not HAVE_RDKIT or MORGAN_GEN is None:
        return None
    mol = Chem.MolFromSmiles(smiles)
    return MORGAN_GEN.GetFingerprint(mol) if mol else None


def _tanimoto(fp1, fp2) -> float:
    return DataStructs.TanimotoSimilarity(fp1, fp2)  # type: ignore


# ---------- misc tiny helpers ---------- #
def _count_lines(path: str) -> int:
    with open(path, "rb") as handle:
        return sum(1 for _ in handle)


def _out_name(path: str, suffix: str) -> str:
    stem, ext = os.path.splitext(path)
    return f"{stem}{suffix}{ext}"


# --------------------------------------------------------------------------- #
# Worker function (executes in subprocesses)
# --------------------------------------------------------------------------- #
def _process_row(
    payload: Tuple[int, List[str], int, bool, str]
) -> Tuple[str, List[str]]:
    silence_rdkit()
    _, row, col_idx, nostereo, tauto_mode = payload

    try:
        smiles = row[col_idx].strip()
    except IndexError:
        return "IndexError", row

    if tauto_mode == "strong":
        smiles = _tautomer_rdkit(smiles) or smiles
    elif tauto_mode == "weak":
        smiles = _tautomer_obabel(smiles) or smiles
    elif tauto_mode == "molvs":
        smiles = _tautomer_molvs(smiles) or smiles

    keep_fixedh = tauto_mode == "off"
    key = (_inchikey_obabel(smiles, nostereo, keep_fixedh)
           or _inchikey_rdkit(smiles, nostereo, keep_fixedh))
    return (key if key else "ParseError", row)


# --------------------------------------------------------------------------- #
# Deduplication core
# --------------------------------------------------------------------------- #
def deduplicate(
    csv_path: str,
    smiles_col: int,
    nostereo: bool,
    tauto_mode: str,
    n_jobs: int,
    tanimoto_thr: Optional[float],
    outlier_factor: Optional[float],
) -> None:
    col_idx = smiles_col - 1  # SMILES column (0-based)
    num_idx = 2               # numeric column for outliers (0-based)

    total_rows = _count_lines(csv_path) - 1
    out_path = _out_name(csv_path, "_singled")
    rej_path = _out_name(csv_path, "_rejects")

    with open(csv_path, newline="", encoding="utf-8") as fh:
        rdr = csv.reader(fh)
        header = next(rdr)
        payloads = [(idx, row, col_idx, nostereo, tauto_mode)
                    for idx, row in enumerate(rdr, start=2)]

    executor_cls = ProcessPoolExecutor if n_jobs > 1 else None
    kept: List[List[str]] = [header]
    rejects: List[List[str]] = []
    seen: OrderedDict[str, None] = OrderedDict()

    def _pbar(it: Iterator):
        return tqdm(it, total=total_rows, unit=" rows", ncols=80) if HAVE_TQDM else it

    # ---------- pass 1 – strict InChIKey ---------- #
    if executor_cls:
        with executor_cls(max_workers=n_jobs) as ex:
            chunk = max(1, len(payloads) // (n_jobs * 5))
            for key, row in _pbar(ex.map(_process_row, payloads,
                                         chunksize=chunk)):
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

    # ---------- pass 2 – fuzzy duplicates ---------- #
    if tanimoto_thr is not None and HAVE_RDKIT:
        fp_cache, final_kept = [], [header]
        for row in kept[1:]:
            smiles = row[col_idx].strip()
            fp = _ecfp4(smiles)
            if fp is None:
                rejects.append(["FPError"] + row)
                continue
            if any(_tanimoto(fp, old) >= tanimoto_thr for old in fp_cache):
                rejects.append(["Tanimoto"] + row)
            else:
                fp_cache.append(fp)
                final_kept.append(row)
        kept = final_kept

    # ---------- pass 3 – outlier trimming ---------- #
    if outlier_factor is not None:
        numeric_vals: List[float] = []
        parse_errors: List[List[str]] = []

        for row in kept[1:]:
            val = coerce_float(row[num_idx]) if len(row) > num_idx else None
            if val is None:
                parse_errors.append(row)
            else:
                numeric_vals.append(val)

        if parse_errors:
            for r in parse_errors:
                rejects.append(["ValueError"] + r)
            logging.warning("%d rows skipped due to non-numeric column 3",
                            len(parse_errors))

        if numeric_vals:
            mu = mean(numeric_vals)
            sigma = stdev(numeric_vals) if len(numeric_vals) > 1 else 0.0
            lo, hi = mu - outlier_factor * sigma, mu + outlier_factor * sigma

            before_vals = numeric_vals.copy()
            kept_after: List[List[str]] = [header]
            after_vals: List[float] = []

            for row in kept[1:]:
                val = coerce_float(row[num_idx]) if len(row) > num_idx else None
                if val is None:
                    continue
                if val < lo or val > hi:
                    rejects.append(["Outlier"] + row)
                else:
                    kept_after.append(row)
                    after_vals.append(val)

            kept = kept_after

            if HAVE_MPL and after_vals:
                _plot_histograms(before_vals, after_vals,
                                 mu, sigma, outlier_factor)
            elif not HAVE_MPL:
                logging.info("matplotlib not installed — skipping histograms")

    # ---------- write outputs ---------- #
    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        csv.writer(fh).writerows(kept)

    if rejects:
        with open(rej_path, "w", newline="", encoding="utf-8") as fh:
            csv.writer(fh).writerows([["REASON", *header], *rejects])

    print(f"\n✓ {out_path} — {len(kept) - 1} rows kept")
    if rejects:
        print(f"⚠ {len(rejects)} rejected rows → {rej_path}")


# --------------------------------------------------------------------------- #
# Plotting helper
# --------------------------------------------------------------------------- #
def _plot_histograms(before: List[float], after: List[float],
                     mu: float, sigma: float, k: float) -> None:
    if not HAVE_MPL:
        return

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    axes[0].hist(before, bins=30)
    axes[0].set_title("Before outlier removal")

    axes[1].hist(after, bins=30)
    axes[1].set_title("After outlier removal")

    fig.suptitle(f"Outlier threshold: μ ± {k}·σ  (μ={mu:.2f}, σ={sigma:.2f})")
    fig.tight_layout()
    plt.show()


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def main() -> None:
    parser = argparse.ArgumentParser(
        prog="curator.py",
        description="Deduplicate SMILES and optionally remove numeric outliers.",
    )
    parser.add_argument("csv_file", help="Path to the input CSV file")
    parser.add_argument("smiles_col", type=int,
                        help="1-based column index of SMILES strings")
    parser.add_argument("--nostereo", action="store_true",
                        help="Merge enantiomer / E,Z stereoisomers")
    parser.add_argument(
        "--tauto",
        choices=["off", "weak", "strong", "molvs"],
        default="weak",
        help="Tautomer handling engine "
             "(off = none, weak = Open Babel, strong = RDKit, molvs = MolVS)",
    )
    parser.add_argument("--njobs", type=int, default=1,
                        help="Number of worker processes (default 1)")
    parser.add_argument("--tanimoto-thres", type=float, metavar="THRESH",
                        help="Fuzzy filter: reject rows with Tanimoto ≥ THRESH")
    parser.add_argument("--outliers", type=float, metavar="K",
                        help="Trim numeric outliers using μ ± K·o")
    parser.add_argument("--quiet", action="store_true",
                        help="Silence RDKit log spam")

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
        outlier_factor=args.outliers,
    )


if __name__ == "__main__":
    main()
