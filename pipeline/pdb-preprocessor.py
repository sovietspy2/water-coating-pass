#!/usr/bin/env python3
import argparse
import os
import sys
import tempfile

try:
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
except ImportError:
    print("[WARNING] PDBFixer not installed; trying without it!", file=sys.stderr)
    sys.exit(0)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Fix a PDB file (missing residues/atoms, nonstandard residues) in place."
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--target",
        action="store_true",
        help="Target PDB mode: remove existing waters and other heterogens.",
    )
    mode.add_argument(
        "--reference",
        action="store_true",
        help="Reference PDB mode: keep waters, remove other heterogens.",
    )
    parser.add_argument("input_pdb", help="Path to the PDB file to fix in place.")
    return parser.parse_args()


def main():
    args = parse_args()
    keep_water = args.reference

    input_pdb = os.path.abspath(args.input_pdb)
    if not os.path.isfile(input_pdb):
        print(f"Error: file not found: {input_pdb}", file=sys.stderr)
        sys.exit(1)

    input_dir = os.path.dirname(input_pdb)
    input_name = os.path.basename(input_pdb)

    fd, tmp_path = tempfile.mkstemp(
        prefix=input_name + ".",
        suffix=".tmp",
        dir=input_dir,
        text=True,
    )
    os.close(fd)

    print(f"[INFO] Mode      : {'reference' if keep_water else 'target'}")
    print(f"[INFO] Input PDB : {input_pdb}")
    print(f"[INFO] TMP file  : {tmp_path}")

    try:
        fixer = PDBFixer(filename=input_pdb)

        # Find missing residues, but do not build missing residues at chain termini
        fixer.findMissingResidues()
        chains = list(fixer.topology.chains())
        for key in list(fixer.missingResidues.keys()):
            chain = chains[key[0]]
            if key[1] == 0 or key[1] == len(list(chain.residues())):
                del fixer.missingResidues[key]

        # Replace nonstandard residues if present
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()

        # Keep waters only in reference mode; always drop other heterogens
        fixer.removeHeterogens(keepWater=keep_water)

        # Add missing heavy atoms and residues
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()

        with open(tmp_path, "w") as handle:
            PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)
            handle.flush()
            os.fsync(handle.fileno())

        os.replace(tmp_path, input_pdb)

        if keep_water:
            print("[DONE] Waters kept; other heterogens removed.")
        else:
            print("[DONE] Waters and other heterogens removed.")
        print("[DONE] Missing heavy atoms and internal residues repaired where possible.")
        print(f"[DONE] Replaced original PDB with fixed version: {input_pdb}")

    except Exception:
        print(f"[ERROR] Failed while fixing PDB: {input_pdb}", file=sys.stderr)
        print(f"[ERROR] TMP file kept at: {tmp_path}", file=sys.stderr)
        raise


if __name__ == "__main__":
    main()
