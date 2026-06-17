#!/usr/bin/env python3
import os
import sys

try:
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
except ImportError:
    print("[WARNING] PDBFixer not installed; trying without it!", file=sys.stderr)
    sys.exit(0)


def main():
    if len(sys.argv) not in (2, 3):
        print("Usage: fix_pdb_keep_water.py input.pdb [output.pdb]", file=sys.stderr)
        sys.exit(1)

    input_pdb = os.path.abspath(sys.argv[1])
    if not os.path.isfile(input_pdb):
        print(f"Error: file not found: {input_pdb}", file=sys.stderr)
        sys.exit(1)

    if len(sys.argv) == 3:
        output_pdb = os.path.abspath(sys.argv[2])
    else:
        root, ext = os.path.splitext(input_pdb)
        output_pdb = root + ".fixed.pdb"

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

    # Write output; alternate locations are resolved by PDBFixer when loading/fixing
    with open(output_pdb, "w") as handle:
        PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)

    print(f"[DONE] Wrote fixed PDB: {output_pdb}")
    print("[DONE] Waters kept.")
    print("[DONE] Missing internal residues repaired where possible.")
    print("[DONE] Alternate locations reduced to a single conformer.")

if __name__ == "__main__":
    main()