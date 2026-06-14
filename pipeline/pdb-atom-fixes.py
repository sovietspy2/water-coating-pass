#!/usr/bin/env python3
import os
import sys
import tempfile
from pdbfixer import PDBFixer
from openmm.app import PDBFile

try:
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
except ImportError:
    print("[WARNING] PDBFixer not installed; trying without it!", file=sys.stderr)
    sys.exit(0)

def main():
    if len(sys.argv) != 2:
        print("Usage: pdb-atom-fixes.py input.pdb", file=sys.stderr)
        sys.exit(1)

    input_pdb = sys.argv[1]

    if not os.path.isfile(input_pdb):
        print(f"Error: file not found: {input_pdb}", file=sys.stderr)
        sys.exit(1)

    input_pdb = os.path.abspath(input_pdb)
    input_dir = os.path.dirname(input_pdb)
    input_name = os.path.basename(input_pdb)

    fd, tmp_path = tempfile.mkstemp(
        prefix=input_name + ".",
        suffix=".tmp",
        dir=input_dir,
        text=True
    )
    os.close(fd)

    print(f"[INFO] Input PDB : {input_pdb}")
    print(f"[INFO] TMP file  : {tmp_path}")

    try:
        fixer = PDBFixer(filename=input_pdb)

        fixer.findMissingResidues()
        chains = list(fixer.topology.chains())
        for key in list(fixer.missingResidues.keys()):
            chain = chains[key[0]]
            if key[1] == 0 or key[1] == len(list(chain.residues())):
                del fixer.missingResidues[key]

        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()

        fixer.removeHeterogens(keepWater=False)

        fixer.findMissingAtoms()
        fixer.addMissingAtoms()

        with open(tmp_path, "w") as handle:
            PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)
            handle.flush()
            os.fsync(handle.fileno())

        os.replace(tmp_path, input_pdb)

        print("[DONE] Waters and other heterogens removed.")
        print("[DONE] Missing heavy atoms repaired where possible.")
        print(f"[DONE] Replaced original PDB with fixed version: {input_pdb}")

    except Exception:
        print(f"[ERROR] Failed while fixing PDB: {input_pdb}", file=sys.stderr)
        print(f"[ERROR] TMP file kept at: {tmp_path}", file=sys.stderr)
        raise

if __name__ == "__main__":
    main()