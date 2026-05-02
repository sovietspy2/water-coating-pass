#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<USAGE
Usage: $(basename "$0") INPUT_PDB [LOG_FILE]

Rewrites a PDB file to:
  1. Remove rows that start with CONECT
  2. Rename OW to O in atom-name field
  3. Replace residue names 001 or 002 or SOL with WAT
  4. Rename HETATM to ATOM
  5. Renumber ATOM serials so they increment by 1

Defaults:
  OUTPUT_PDB = overwrite INPUT_PDB
  BACKUP_PDB = <input_basename>.original.pdb
  LOG_FILE   = <input_basename>.log
USAGE
}

if [[ ${1:-} == "-h" || ${1:-} == "--help" ]]; then
  usage
  exit 0
fi

if [[ $# -lt 1 || $# -gt 2 ]]; then
  usage >&2
  exit 1
fi

input=$1
if [[ ! -f "$input" ]]; then
  echo "Error: input file not found: $input" >&2
  exit 1
fi

input_dir=$(dirname -- "$input")
input_name=$(basename -- "$input")

if [[ $input_name == *.pdb ]]; then
  input_stem=${input_name%.pdb}
  backup="$input_dir/${input_stem}.original.pdb"
else
  input_stem=$input_name
  backup="$input.original.pdb"
fi

output="$input"

if [[ $# -ge 2 ]]; then
  logfile=$2
else
  logfile="$input_dir/${input_stem}.log"
fi

tmp_output=$(mktemp)
cleanup() {
  rm -f -- "$tmp_output"
}
trap cleanup EXIT

awk -v log_file="$logfile" '
function trim(s) {
  gsub(/^ +| +$/, "", s)
  return s
}

function pad_right(s, width,    out) {
  out = s
  while (length(out) < width) out = out " "
  return substr(out, 1, width)
}

function pad_left(s, width,    out) {
  out = s
  while (length(out) < width) out = " " out
  if (length(out) > width) out = substr(out, length(out) - width + 1)
  return out
}

function atom_field(name,    n) {
  name = trim(name)
  if (name == "") name = "O"
  n = length(name)
  if (n >= 4) return substr(name, 1, 4)
  if (n == 1) return " " name "  "
  if (n == 2) return " " name " "
  if (n == 3) return " " name
  return pad_right(name, 4)
}

function safe_num(s, fallback) {
  s = trim(s)
  if (s ~ /^-?[0-9]+$/) return s + 0
  return fallback
}

function log_change(kind, input_line, out_serial, msg, before, after) {
  print "[" kind "] input_line=" input_line ", output_serial=" out_serial ": " msg >> log_file
  if (before != "") print "  before: " before >> log_file
  if (after  != "") print "  after : " after >> log_file
}

BEGIN {
  removed_conect = 0
  changed_rows = 0
  hetatm_to_atom = 0
  ow_to_o = 0
  res_to_wat = 0
  serial_fixed = 0
  output_atom_rows = 0
  total_input_rows = 0
  next_serial = ""
  first_atom_seen = 0
  print "PDB rewrite log" > log_file
  print "===============" >> log_file
}

{
  total_input_rows++

  if ($0 ~ /^CONECT/) {
    removed_conect++
    log_change("REMOVED", NR, "-", "row starts with CONECT", $0, "")
    next
  }

  record = substr($0, 1, 6)
  is_atom = (record == "ATOM  " || record == "HETATM")

  if (!is_atom) {
    print $0
    next
  }

  old_line = $0
  old_record = record
  old_serial_raw = substr($0, 7, 5)
  old_serial = safe_num(old_serial_raw, 0)
  old_atom = substr($0, 13, 4)
  old_altloc = substr($0, 17, 1)
  old_res = substr($0, 18, 3)
  old_chain = substr($0, 22, 1)
  old_resseq = substr($0, 23, 4)
  old_icode = substr($0, 27, 1)
  old_x = substr($0, 31, 8)
  old_y = substr($0, 39, 8)
  old_z = substr($0, 47, 8)
  old_occ = substr($0, 55, 6)
  old_temp = substr($0, 61, 6)
  old_tail = (length($0) >= 67 ? substr($0, 67) : "")

  new_record = "ATOM  "
  new_atom_trim = trim(old_atom)
  new_res_trim = trim(old_res)

  row_changes = ""
  row_changed = 0

  if (old_record == "HETATM") {
    hetatm_to_atom++
    row_changed = 1
    row_changes = row_changes "HETATM->ATOM; "
  }

  if (new_atom_trim == "OW") {
    new_atom_trim = "O"
    ow_to_o++
    row_changed = 1
    row_changes = row_changes "OW->O; "
  }

  if (new_res_trim == "001" || new_res_trim == "002" || new_res_trim == "SOL") {
    new_res_trim = "WAT"
    res_to_wat++
    row_changed = 1
    row_changes = row_changes "resname->WAT; "
  }

  if (!first_atom_seen) {
    if (trim(old_serial_raw) ~ /^[0-9]+$/) next_serial = old_serial
    else next_serial = 1
    first_atom_seen = 1
  }

  new_serial = next_serial
  next_serial++

  if (old_serial != new_serial) {
    serial_fixed++
    row_changed = 1
    row_changes = row_changes "serial " old_serial "->" new_serial "; "
  }

  new_line = sprintf("%-6s%5d %-4s%1s%3s %1s%4s%1s   %8s%8s%8s%6s%6s%s",
                     new_record,
                     new_serial,
                     atom_field(new_atom_trim),
                     (old_altloc == "" ? " " : old_altloc),
                     pad_right(new_res_trim, 3),
                     (old_chain == "" ? " " : old_chain),
                     pad_left(trim(old_resseq), 4),
                     (old_icode == "" ? " " : old_icode),
                     pad_left(trim(old_x), 8),
                     pad_left(trim(old_y), 8),
                     pad_left(trim(old_z), 8),
                     pad_left(trim(old_occ), 6),
                     pad_left(trim(old_temp), 6),
                     old_tail)

  print new_line
  output_atom_rows++

  if (row_changed || new_line != old_line) {
    changed_rows++
    if (row_changes == "") row_changes = "format normalization; "
    sub(/; $/, "", row_changes)
    log_change("CHANGED", NR, new_serial, row_changes, old_line, new_line)
  }
}

END {
  total_changed = removed_conect + changed_rows
  print "" >> log_file
  print "Summary" >> log_file
  print "-------" >> log_file
  print "Total input rows        : " total_input_rows >> log_file
  print "Output atom rows        : " output_atom_rows >> log_file
  print "Removed CONECT rows     : " removed_conect >> log_file
  print "Changed output rows     : " changed_rows >> log_file
  print "HETATM -> ATOM changes  : " hetatm_to_atom >> log_file
  print "OW -> O changes         : " ow_to_o >> log_file
  print "001/002 -> WAT changes  : " res_to_wat >> log_file
  print "Serial fixes            : " serial_fixed >> log_file
  print "Total rows changed      : " total_changed >> log_file
}
' "$input" > "$tmp_output"

if [[ -e "$backup" ]]; then
  rm -f -- "$backup"
fi
mv -- "$input" "$backup"
mv -- "$tmp_output" "$output"

echo "Wrote cleaned PDB to: $output"
echo "Original PDB backed up to: $backup"
echo "Wrote log to: $logfile"