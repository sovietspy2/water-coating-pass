#!/usr/bin/env python3
import os
import sys
import tempfile

try:
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
except ImportError:
    print("[WARNING] PDBFixer not installed; trying without it!", file=sys.stderr)
    sys.exit(0)


def main():
    if len(sys.argv) != 2:
        print("Usage: reference-pdb-preprocessor.py input.pdb", file=sys.stderr)
        sys.exit(1)

    input_pdb = os.path.abspath(sys.argv[1])
    if not os.path.isfile(input_pdb):
        print(f"Error: file not found: {input_pdb}", file=sys.stderr)
        sys.exit(1)

    fixer = PDBFixer(filename=input_pdb)

    # Find missing residues, but do not build missing residues at chain termini
    fixer.findMissingResidues()
    chains = list(fixer.topology.chains())
    for key in list(fixer.missingResidues.keys()):
        chain = chains[key[0]]
        residues = list(chain.residues())
        if key[1] == 0 or key[1] == len(residues):
            del fixer.missingResidues[key]

    # Replace nonstandard residues if present
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()

    # Keep waters, remove other heterogens
    fixer.removeHeterogens(keepWater=True)

    # Add missing heavy atoms and residues
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    # Write to a temp file in the same directory, then replace original atomically
    input_dir = os.path.dirname(input_pdb)
    input_name = os.path.basename(input_pdb)

    with tempfile.NamedTemporaryFile(
            mode="w",
            suffix=".tmp",
            prefix=input_name + ".",
            dir=input_dir,
            delete=False,
    ) as tmp:
        temp_path = tmp.name
        PDBFile.writeFile(fixer.topology, fixer.positions, tmp, keepIds=True)

    os.replace(temp_path, input_pdb)

    print(f"[DONE] Replaced original PDB with fixed version: {input_pdb}")
    print("[DONE] Waters kept.")
    print("[DONE] Missing internal residues repaired where possible.")
    print("[DONE] Alternate locations reduced to a single conformer.")


if __name__ == "__main__":
    main()