#!/usr/bin/env sh
set -eu

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
  echo "Usage: $0 input.pdb [output.pdb]" >&2
  exit 1
fi

infile=$1
outfile=${2:-"${infile%.pdb}_fixed.pdb"}

awk '
function trim(s) {
  sub(/^[[:space:]]+/, "", s)
  sub(/[[:space:]]+$/, "", s)
  return s
}

function atom_field(name,    n, l) {
  n = trim(name)
  l = length(n)

  if (l == 1) return " " n "  "
  if (l == 2) return " " n " "
  if (l == 3) return " " n
  return substr(n, 1, 4)
}

{
  rec = substr($0, 1, 6)

  if (rec == "ATOM  " || rec == "HETATM") {
    line = sprintf("%-80s", $0)

    serial  = trim(substr(line,  7, 5))
    name    = trim(substr(line, 13, 4))
    altLoc  = substr(line, 17, 1)
    resName = trim(substr(line, 18, 3))
    chainID = substr(line, 22, 1)
    resSeq  = trim(substr(line, 23, 4))
    iCode   = substr(line, 27, 1)
    x       = trim(substr(line, 31, 8))
    y       = trim(substr(line, 39, 8))
    z       = trim(substr(line, 47, 8))
    occ     = trim(substr(line, 55, 6))
    temp    = trim(substr(line, 61, 6))
    element = trim(substr(line, 77, 2))
    charge  = trim(substr(line, 79, 2))

    if (serial == "") serial = 0
    if (resSeq == "") resSeq = 0
    if (x == "") x = 0
    if (y == "") y = 0
    if (z == "") z = 0
    if (occ == "") occ = 0
    if (temp == "") temp = 0

    aname = atom_field(name)

    printf "%-6s%5d %4s%1s%-3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f",
      rec, serial, aname, altLoc, resName, chainID, resSeq, iCode,
      x+0, y+0, z+0, occ+0, temp+0

    if (element != "" || charge != "") {
      printf "          %2s%2s", element, charge
    }

    printf "\n"
  } else {
    print $0
  }
}
' "$infile" > "$outfile"