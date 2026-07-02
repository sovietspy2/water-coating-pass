#!/usr/bin/env python3
"""Integrity test for pipeline/pdb-preprocessor.py.

Goal: verify the preprocessor only *adds* atoms and never *moves* an existing
atom (protein heavy atoms or waters). It:

  1. copies an input PDB to a scratch working copy,
  2. runs pdb-preprocessor.py (in --target or --reference mode) on the copy,
  3. matches every pre-existing atom in the input to its counterpart in the
     output and compares coordinates,
  4. prints a report to STDOUT and exits non-zero if any real coordinate change
     is detected.

Coordinate policy (PDB stores coordinates to 3 decimal places):
  * max per-axis |delta| == 0.000       -> unchanged
  * 0.000 <  max per-axis |delta| <= 0.001 -> WARNING (3rd-decimal / rounding-level
                                              change, tolerated)
  * max per-axis |delta| >  0.001       -> ERROR (moved by more than the 3rd
                                              decimal place)
"""
import argparse
import os
import shutil
import subprocess
import sys
import tempfile

# --- configuration ---------------------------------------------------------

# Water residue names we treat as waters when reporting movement.
WATER_RESNAMES = {"HOH", "WAT", "H2O", "TIP3", "TIP", "SOL", "DOD"}

# A per-axis change of one unit in the 3rd decimal place (0.001 Angstrom) is
# tolerated as rounding noise and reported as a warning; anything larger is an
# error.
WARN_TOLERANCE = 0.001

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DEFAULT_PREPROCESSOR = os.path.join(REPO_ROOT, "pipeline", "pdb-preprocessor.py")
DEFAULT_PYTHON = os.path.join(REPO_ROOT, ".venv", "bin", "python")


# --- PDB parsing ------------------------------------------------------------

class Atom:
    __slots__ = ("record", "name", "alt_loc", "res_name",
                 "chain_id", "res_seq", "i_code", "x", "y", "z", "raw")

    def __init__(self, line):
        self.record = line[0:6].strip()
        self.name = line[12:16].strip()
        self.alt_loc = line[16].strip()
        self.res_name = line[17:20].strip()
        self.chain_id = line[21].strip()
        self.res_seq = line[22:26].strip()
        self.i_code = line[26].strip()
        self.x = float(line[30:38])
        self.y = float(line[38:46])
        self.z = float(line[46:54])
        self.raw = line.rstrip("\n")

    @property
    def is_water(self):
        return self.res_name in WATER_RESNAMES

    @property
    def match_key(self):
        # Identify a pre-existing atom independent of serial number (which the
        # preprocessor renumbers) and of resName (which replaceNonstandardResidues
        # may rewrite). altLoc is ignored so collapsed alternate conformers still
        # match; duplicates are resolved by nearest coordinate at compare time.
        return (self.chain_id, self.res_seq, self.i_code, self.name)

    def coords(self):
        return (self.x, self.y, self.z)

    def label(self):
        alt = f":{self.alt_loc}" if self.alt_loc else ""
        return (f"{self.res_name} {self.chain_id}{self.res_seq}{self.i_code} "
                f"{self.name}{alt}")


def parse_pdb(path, first_model_only=True):
    """Parse ATOM/HETATM records. By default only the first model is read,
    because the preprocessor reduces multi-model (e.g. NMR) inputs to a single
    model (keeping the first); comparing against the whole stack is meaningless.
    """
    atoms = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(("ATOM  ", "HETATM")):
                atoms.append(Atom(line))
            elif first_model_only and line.startswith("ENDMDL"):
                break
    return atoms


def count_models(path):
    """Number of MODEL records; 0 means a single implicit model."""
    n = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("MODEL "):
                n += 1
    return n


# --- comparison -------------------------------------------------------------

def axis_deltas(a, b):
    return (abs(a.x - b.x), abs(a.y - b.y), abs(a.z - b.z))


def run_preprocessor(python, script, mode_flag, pdb_path):
    cmd = [python, script, mode_flag, pdb_path]
    print(f"[INFO] Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.stdout:
        for ln in result.stdout.splitlines():
            print(f"    | {ln}")
    if result.returncode != 0:
        sys.stderr.write(result.stderr)
        raise SystemExit(f"[FATAL] preprocessor exited with {result.returncode}")
    return result.stdout + result.stderr


def build_output_index(out_atoms):
    index = {}
    for atom in out_atoms:
        index.setdefault(atom.match_key, []).append(atom)
    return index


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input_pdb", help="input PDB to feed to the preprocessor")
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--target", action="store_true",
                      help="run preprocessor in --target mode (drops all waters)")
    mode.add_argument("--reference", action="store_true",
                      help="run preprocessor in --reference mode (keeps waters, default)")
    ap.add_argument("--preprocessor", default=DEFAULT_PREPROCESSOR,
                    help="path to pdb-preprocessor.py")
    ap.add_argument("--python", default=DEFAULT_PYTHON if os.path.exists(DEFAULT_PYTHON)
                    else sys.executable,
                    help="python interpreter used to run the preprocessor")
    ap.add_argument("--keep-output", action="store_true",
                    help="keep copies next to the input for manual comparison: "
                         "<base>.ORIGINAL.pdb (untouched input) and "
                         "<base>.MODIFIED.pdb (preprocessor output)")
    args = ap.parse_args()

    mode_flag = "--target" if args.target else "--reference"

    input_pdb = os.path.abspath(args.input_pdb)
    if not os.path.isfile(input_pdb):
        raise SystemExit(f"[FATAL] input PDB not found: {input_pdb}")
    if not os.path.isfile(args.preprocessor):
        raise SystemExit(f"[FATAL] preprocessor not found: {args.preprocessor}")

    # Parse the original before it can be touched (first model only).
    in_atoms = parse_pdb(input_pdb)
    in_models = count_models(input_pdb)

    # Run the (in-place) preprocessor on a scratch copy.
    kept_paths = None
    tmp_dir = tempfile.mkdtemp(prefix="pdb-preprocess-test.")
    work_copy = os.path.join(tmp_dir, os.path.basename(input_pdb))
    try:
        shutil.copyfile(input_pdb, work_copy)
        run_preprocessor(args.python, args.preprocessor, mode_flag, work_copy)
        out_atoms = parse_pdb(work_copy)
        out_models = count_models(work_copy)

        if args.keep_output:
            base = os.path.splitext(os.path.basename(input_pdb))[0]
            out_dir = os.path.dirname(input_pdb)
            original_kept = os.path.join(out_dir, f"{base}.ORIGINAL.pdb")
            modified_kept = os.path.join(out_dir, f"{base}.MODIFIED.pdb")
            shutil.copyfile(input_pdb, original_kept)
            shutil.copyfile(work_copy, modified_kept)
            kept_paths = (original_kept, modified_kept)
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)

    out_index = build_output_index(out_atoms)
    # Track which output atoms got consumed by a match so leftovers are "added".
    consumed = {key: [False] * len(v) for key, v in out_index.items()}

    matched = 0
    unchanged = 0
    warnings = []          # (input atom, output atom, max_delta)
    errors = []            # (input atom, output atom, max_delta)
    removed = []           # input atoms whose identity is absent from the output
    collapsed = []         # input atoms dropped by altLoc collapse (identity kept)
    water_moved = 0        # matched water atoms flagged as errors

    for atom in in_atoms:
        candidates = out_index.get(atom.match_key, [])
        flags = consumed.get(atom.match_key)
        # Choose the nearest not-yet-consumed candidate (handles altLoc collapse).
        best_i, best_delta = None, None
        for i, cand in enumerate(candidates):
            if flags[i]:
                continue
            d = max(axis_deltas(atom, cand))
            if best_delta is None or d < best_delta:
                best_delta, best_i = d, i
        if best_i is None:
            # No free counterpart. If a same-identity atom survived in the
            # output (candidates existed but were all consumed by other
            # conformers), this is an alternate-location conformer that was
            # collapsed away, not a real removal. Otherwise the atom's
            # identity is genuinely absent from the output.
            if candidates:
                collapsed.append(atom)
            else:
                removed.append(atom)
            continue

        flags[best_i] = True
        out_atom = candidates[best_i]
        matched += 1
        max_delta = round(max(axis_deltas(atom, out_atom)), 3)

        if max_delta == 0.0:
            unchanged += 1
        elif max_delta <= WARN_TOLERANCE:
            warnings.append((atom, out_atom, max_delta))
        else:
            errors.append((atom, out_atom, max_delta))
            if atom.is_water:
                water_moved += 1

    added = [cand for key, cands in out_index.items()
             for i, cand in enumerate(cands) if not consumed[key][i]]

    waters_in = sum(1 for a in in_atoms if a.is_water)
    waters_removed = sum(1 for a in removed if a.is_water)
    hetero_removed = sum(1 for a in removed if not a.is_water)
    water_collapsed = sum(1 for a in collapsed if a.is_water)
    other_collapsed = sum(1 for a in collapsed if not a.is_water)

    # --- report -------------------------------------------------------------
    print()
    print("=" * 68)
    print("PDB Preprocessor Integrity Report")
    print("=" * 68)
    print(f"Input      : {input_pdb}")
    print(f"Preprocessor: {args.preprocessor}  (mode {mode_flag})")
    print()
    in_eff = max(1, in_models)
    out_eff = max(1, out_models)
    model_reduction_ok = out_eff == 1
    print("Model reduction (preprocessor must collapse to a single model)")
    print(f"  input models               : {in_eff}"
          f"{' (multi-model / NMR)' if in_eff > 1 else ''}")
    print(f"  output models              : {out_eff}")
    print(f"  reduced to a single model  : {'YES' if model_reduction_ok else 'NO -> FAIL'}")
    print()
    print("Atom counts (first model only)")
    print(f"  input atoms                : {len(in_atoms)}")
    print(f"  output atoms               : {len(out_atoms)}")
    print(f"  matched (pre-existing)     : {matched}")
    print(f"  ADDED (new in output)      : {len(added)}")
    print(f"  removed (identity gone)    : {len(removed)}")
    print(f"      waters actually removed: {waters_removed} (of {waters_in} input water atoms)")
    print(f"      other heterogens       : {hetero_removed}")
    print(f"  alt conformers collapsed   : {len(collapsed)} (identity kept, extra altLoc dropped)")
    print(f"      water alt conformers   : {water_collapsed}")
    print(f"      protein/other alt conf.: {other_collapsed}")
    print()
    print("Coordinate integrity of matched atoms")
    print(f"  unchanged (delta = 0.000 A)          : {unchanged}")
    print(f"  WARNINGS  (0 < delta <= {WARN_TOLERANCE:.3f} A)   : {len(warnings)}")
    print(f"  ERRORS    (delta > {WARN_TOLERANCE:.3f} A)        : {len(errors)}")
    print(f"      of which are waters that moved   : {water_moved}")
    print()

    if warnings:
        print("Warnings (rounding-level coordinate changes, tolerated):")
        for a, b, d in warnings:
            print(f"  [WARN] {a.label():<28} moved {d:.3f} A  "
                  f"({a.coords()} -> {b.coords()})")
        print()

    if errors:
        print("Errors (existing atoms moved by more than the 3rd decimal place):")
        for a, b, d in errors:
            tag = "WATER " if a.is_water else "ATOM  "
            print(f"  [ERROR] {tag}{a.label():<28} moved {d:.3f} A  "
                  f"({a.coords()} -> {b.coords()})")
        print()

    if len(added) == 0 and matched == 0:
        print("[NOTE] No atoms matched or added; preprocessor may not have run "
              "(pdbfixer missing?).")
        print()

    if kept_paths:
        original_kept, modified_kept = kept_paths
        print("Kept files for manual comparison:")
        print(f"  ORIGINAL (untouched input) : {original_kept}")
        print(f"  MODIFIED (preprocessed)    : {modified_kept}")
        print()

    result_ok = len(errors) == 0 and model_reduction_ok
    print("=" * 68)
    print(f"RESULT: {'PASS' if result_ok else 'FAIL'}  "
          f"({len(added)} added, {len(errors)} error(s), {len(warnings)} warning(s), "
          f"models {in_eff}->{out_eff})")
    print("=" * 68)

    sys.exit(0 if result_ok else 1)


if __name__ == "__main__":
    main()
