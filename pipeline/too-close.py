#!/usr/bin/env python3
import math
import argparse
import sys


def infer_element(atom_name, element_field=""):
    ef = element_field.strip().upper()
    if ef:
        return ef
    name = ''.join(c for c in atom_name.upper() if c.isalpha())
    if not name:
        return ''
    two = {"CL", "BR", "NA", "MG", "ZN", "FE", "CA", "MN", "CU", "CO", "NI", "CD", "HG"}
    return name[:2] if len(name) >= 2 and name[:2] in two else name[0]


def read_pdb(path):
    atoms = []
    with open(path, 'r', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec not in {"ATOM", "HETATM"}:
                continue
            try:
                serial = int(line[6:11])
                name = line[12:16].strip()
                altloc = line[16].strip()
                resname = line[17:20].strip()
                chain = line[21].strip() or '-'
                resseq = line[22:26].strip()
                icode = line[26].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                occ = float(line[54:60]) if len(line) >= 60 and line[54:60].strip() else None
                element = infer_element(name, line[76:78] if len(line) >= 78 else '')
            except ValueError:
                continue
            atoms.append({
                'serial': serial,
                'name': name,
                'altloc': altloc,
                'resname': resname,
                'chain': chain,
                'resseq': resseq,
                'icode': icode,
                'x': x,
                'y': y,
                'z': z,
                'occ': occ,
                'element': element,
            })
    return atoms


def atom_label(a):
    alt = f"/{a['altloc']}" if a['altloc'] else ''
    ins = a['icode'] if a['icode'] else ''
    return f"{a['serial']}:{a['name']}{alt} {a['resname']} {a['chain']}{a['resseq']}{ins}"


def same_residue(a, b):
    return (a['chain'], a['resseq'], a['icode'], a['resname']) == (b['chain'], b['resseq'], b['icode'], b['resname'])


def dist(a, b):
    dx = a['x'] - b['x']
    dy = a['y'] - b['y']
    dz = a['z'] - b['z']
    return math.sqrt(dx*dx + dy*dy + dz*dz)


def is_hydrogen(a):
    return a['element'] == 'H'


def is_oxygen(a):
    return a['element'] == 'O'


def is_heavy(a):
    return a['element'] not in {'', 'H'}


def acidic_oxygen(a):
    if a['element'] != 'O':
        return False
    rn = a['resname'].upper()
    an = a['name'].upper()
    return (
            (rn == 'ASP' and an in {'OD1', 'OD2'}) or
            (rn == 'GLU' and an in {'OE1', 'OE2'}) or
            an == 'OXT' or
            (rn in {'PO4', 'DOP', 'ADP', 'ATP'} and an.startswith('O'))
    )


def positive_nitrogen(a):
    if a['element'] != 'N':
        return False
    rn = a['resname'].upper()
    an = a['name'].upper()
    return (
            (rn == 'LYS' and an == 'NZ') or
            (rn == 'ARG' and an in {'NE', 'NH1', 'NH2'}) or
            (rn in {'HIS', 'HID', 'HIE', 'HIP'} and an in {'ND1', 'NE2'})
    )


def pair_allowed(a, b, args):
    if args.ignore_hydrogen and (is_hydrogen(a) or is_hydrogen(b)):
        return False
    if args.skip_same_residue and same_residue(a, b):
        return False
    if args.intermolecular_only and a['chain'] == b['chain']:
        return False
    if args.ignore_altloc and a['altloc'] and b['altloc'] and a['altloc'] != b['altloc']:
        return False
    return True


def main():
    ap = argparse.ArgumentParser(description='Check a PDB for suspicious nonbonded contacts relevant to MD.')
    ap.add_argument('pdb', help='Input PDB file')
    ap.add_argument('--oo-cutoff', type=float, default=2.4, help='Flag O-O contacts shorter than this (A); default 2.4')
    ap.add_argument('--heavy-cutoff', type=float, default=1.6, help='Flag heavy-atom clashes shorter than this (A); default 1.6')
    ap.add_argument('--ionic-cutoff', type=float, default=2.2, help='Flag acidic O to positive N contacts shorter than this (A); default 2.2')
    ap.add_argument('--dup-cutoff', type=float, default=0.05, help='Flag near-duplicate coordinates shorter than this (A); default 0.05')
    ap.add_argument('--skip-same-residue', action='store_true', default=True, help='Ignore atom pairs within the same residue (default on)')
    ap.add_argument('--include-same-residue', dest='skip_same_residue', action='store_false', help='Also check pairs within the same residue')
    ap.add_argument('--intermolecular-only', action='store_true', help='Only compare atoms from different chains')
    ap.add_argument('--ignore-hydrogen', action='store_true', default=True, help='Ignore hydrogens (default on)')
    ap.add_argument('--include-hydrogen', dest='ignore_hydrogen', action='store_false', help='Include hydrogens')
    ap.add_argument('--ignore-altloc', action='store_true', default=True, help='Ignore pairs from different altlocs (default on)')
    ap.add_argument('--include-altloc', dest='ignore_altloc', action='store_false', help='Include different altloc pairs too')
    args = ap.parse_args()

    atoms = read_pdb(args.pdb)
    if not atoms:
        print('No readable ATOM/HETATM records found.', file=sys.stderr)
        sys.exit(1)

    oo_hits = []
    heavy_hits = []
    ionic_hits = []
    dup_hits = []

    n = len(atoms)
    for i in range(n):
        ai = atoms[i]
        for j in range(i + 1, n):
            aj = atoms[j]
            if not pair_allowed(ai, aj, args):
                continue
            d = dist(ai, aj)

            if d < args.dup_cutoff:
                dup_hits.append((d, ai, aj))
            if is_heavy(ai) and is_heavy(aj) and d < args.heavy_cutoff:
                heavy_hits.append((d, ai, aj))
            if is_oxygen(ai) and is_oxygen(aj) and d < args.oo_cutoff:
                oo_hits.append((d, ai, aj))
            if ((acidic_oxygen(ai) and positive_nitrogen(aj)) or (acidic_oxygen(aj) and positive_nitrogen(ai))) and d < args.ionic_cutoff:
                ionic_hits.append((d, ai, aj))

    for hits in (dup_hits, heavy_hits, oo_hits, ionic_hits):
        hits.sort(key=lambda x: x[0])

    print(f'Atoms read: {len(atoms)}')
    print(f'skip_same_residue = {args.skip_same_residue}')
    print(f'intermolecular_only = {args.intermolecular_only}')
    print(f'ignore_hydrogen = {args.ignore_hydrogen}')
    print(f'ignore_altloc = {args.ignore_altloc}')
    print()

    print(f'Near-duplicate coordinates (< {args.dup_cutoff:.2f} A):')
    if dup_hits:
        for d, a, b in dup_hits:
            print(f'{atom_label(a)}  <->  {atom_label(b)}    {d:.3f} A')
    else:
        print('None found.')
    print()

    print(f'Severe heavy-atom clashes (< {args.heavy_cutoff:.2f} A):')
    if heavy_hits:
        for d, a, b in heavy_hits:
            print(f'{atom_label(a)}  <->  {atom_label(b)}    {d:.3f} A')
    else:
        print('None found.')
    print()

    print(f'Suspicious nonbonded O-O contacts (< {args.oo_cutoff:.2f} A):')
    if oo_hits:
        for d, a, b in oo_hits:
            print(f'{atom_label(a)}  <->  {atom_label(b)}    {d:.3f} A')
    else:
        print('None found.')
    print()

    print(f'Very short acidic-O ... positive-N contacts (< {args.ionic_cutoff:.2f} A):')
    if ionic_hits:
        for d, a, b in ionic_hits:
            print(f'{atom_label(a)}  <->  {atom_label(b)}    {d:.3f} A')
    else:
        print('None found.')


if __name__ == '__main__':
    main()