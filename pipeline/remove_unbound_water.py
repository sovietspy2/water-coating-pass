#!/usr/bin/env python3
import re
import math
import os
import sys
import argparse

# ---------------------------
# User settings
# ---------------------------
REP_NAME = "removed_waters.tsv"

WATER_RESNAMES = {"SOL", "HOH", "WAT"}
DCUT = 30.0  # Å: delete water residues whose oxygen is farther than this from ANY protein heavy atom

# ---------------------------
# Input / output paths
# ---------------------------
class Parser(argparse.ArgumentParser):
    def error(self, message):
        self.print_usage(sys.stderr)
        self.exit(1, f"{self.prog}: error: {message}\n")


parser = Parser(description="Remove distant water residues from a PDB and renumber atom serials.")
parser.add_argument("INP", help="input PDB path")
parser.add_argument("OUT", help="output PDB path")
args = parser.parse_args()

INP = args.INP
INP = os.path.abspath(INP)
OUT = os.path.abspath(args.OUT)

if not os.path.isfile(INP):
    raise SystemExit(f"Input file not found: {INP}")

WORKDIR = os.path.dirname(OUT)
REP = os.path.join(WORKDIR, REP_NAME)

# ---------------------------
# Helpers (robust to stuck-together coords like 1132.340-665.460-4394.330)
# ---------------------------
float_pat = re.compile(r'[-+]?(?:\d+\.\d+|\d+\.|\.\d+|\d+)')

def rec(line):   return line[0:6]                  # includes trailing spaces
def atom(line):  return line[12:16].strip() if len(line) >= 16 else ""
def resn(line):  return line[17:20].strip() if len(line) >= 20 else ""
def elem(line):  return line[76:78].strip() if len(line) >= 78 else ""

def resseq(line):
    # Standard PDB: residue sequence number is columns 23–26 (1-based) => [22:26] (0-based)
    if len(line) >= 26:
        s = line[22:26].strip()
        if s and s.lstrip("-").isdigit():
            return int(s)
    # Fallback for malformed spacing
    m = re.search(r'\b(?:SOL|HOH|WAT)\b\s*([0-9]+)', line)
    return int(m.group(1)) if m else None

def xyz(line):
    # Try to pick floats from around the coordinate area first (classic x/y/z live around col 31+)
    s = line[30:80] if len(line) >= 54 else line
    nums = float_pat.findall(s)
    if len(nums) < 3:
        nums = float_pat.findall(line)
    if len(nums) >= 3:
        return float(nums[0]), float(nums[1]), float(nums[2])
    return None

def set_serial_atom(line, serial):
    # ATOM/HETATM serial is columns 7–11 (1-based) => [6:11] (0-based)
    if len(line) < 11:
        line = line.rstrip("\n").ljust(11) + "\n"
    return f"{line[:6]}{serial:5d}{line[11:]}"

def min_dist(point, prot_coords):
    x, y, z = point
    best = float("inf")
    for X, Y, Z in prot_coords:
        dx, dy, dz = x - X, y - Y, z - Z
        d2 = dx*dx + dy*dy + dz*dz
        if d2 < best:
            best = d2
    return math.sqrt(best) if best < float("inf") else float("inf")

# ---------------------------
# Pass 0: read file once
# ---------------------------
with open(INP, "r", errors="ignore") as f:
    raw_lines = f.readlines()

# ---------------------------
# Pass 1: collect protein heavy-atom coordinates (non-water, non-H)
# ---------------------------
prot = []
for line in raw_lines:
    rtype = rec(line).strip()
    if rtype not in {"ATOM", "HETATM"}:
        continue

    if resn(line) in WATER_RESNAMES:
        continue

    p = xyz(line)
    if p is None:
        continue

    is_H = (elem(line) == "H") or atom(line).startswith("H")
    if not is_H:
        prot.append(p)

if not prot:
    raise SystemExit("No protein heavy atoms found (check WATER_RESNAMES / input format).")

# ---------------------------
# Pass 2: identify water residues to delete (use oxygen OW/O)
# ---------------------------
water_O = {}  # (resname, resseq) -> (x,y,z)
for line in raw_lines:
    rtype = rec(line).strip()
    if rtype not in {"ATOM", "HETATM"}:
        continue

    r = resn(line)
    if r not in WATER_RESNAMES:
        continue

    rseq = resseq(line)
    if rseq is None:
        continue

    if atom(line) in {"OW", "O"}:
        p = xyz(line)
        if p is not None:
            water_O[(r, rseq)] = p

to_delete = set()
with open(REP, "w") as rep:
    rep.write("resname\tresseq\tx\ty\tz\tmin_dist_to_protein_A\n")
    for key in sorted(water_O.keys()):
        p = water_O[key]
        d = min_dist(p, prot)
        if d > DCUT:
            to_delete.add(key)
            rep.write(f"{key[0]}\t{key[1]}\t{p[0]:.3f}\t{p[1]:.3f}\t{p[2]:.3f}\t{d:.3f}\n")

# ---------------------------
# Pass 3: write filtered PDB, renumber ATOM/HETATM serials only
# - TER lines are written as plain "TER\n" (no index) per your requirement
# ---------------------------
serial = 1
with open(OUT, "w") as out:
    for line in raw_lines:
        rtype = rec(line).strip()

        # Drop whole water residue if flagged
        if rtype in {"ATOM", "HETATM"}:
            r = resn(line)
            if r in WATER_RESNAMES:
                rseq = resseq(line)
                if rseq is not None and (r, rseq) in to_delete:
                    continue

        if rtype in {"ATOM", "HETATM"}:
            out.write(set_serial_atom(line, serial))
            serial += 1
        elif rtype == "TER":
            out.write("TER\n")
        else:
            out.write(line)

print(f"Wrote: {OUT}")
print(f"Removed waters report: {REP}")
print(f"Deleted water residues: {len(to_delete)} (criterion: min_dist(O, protein_heavy) > {DCUT} Å)")
