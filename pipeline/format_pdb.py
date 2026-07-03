#!/usr/bin/env python3
"""Canonicalize a PDB file for the WDROP pipeline.

Rewrites a PDB in place to:
  1. Remove rows that start with CONECT.
  2. Rename the OW atom name to O.
  3. Replace residue names 001 / 002 / SOL with WAT.
  4. Re-emit every ATOM/HETATM record in canonical fixed-width PDB columns.

The original input is backed up to ``<stem>.original.pdb`` and a change log is
written to ``<stem>.log`` (overridable via a second CLI argument).

Why this replaces the old ``format-pdb.sh``
-------------------------------------------
PDB coordinates are a FIXED-WIDTH ``%8.3f`` format, so adjacent fields legally
touch with no separating whitespace whenever a value fills its 8-column field
(e.g. GROMACS pseudo-PBC coordinates near ~5000 A print ``   0.131`` + ``5005.730``
as ``0.1315005.730``). A greedy ``[0-9]+\\.[0-9]+`` regex bleeds x into y on such
lines. Conversely, Tinker's ``xyzpdb`` WIDENS the field to ``%9.3f`` for any frame
containing an atom >= 100 A from the origin, shifting every fixed column.

The coordinate parser here handles all three cases:
  * standard, space-separated lines,
  * touching (but standard-width) lines from GROMACS, and
  * widened / shifted lines and stuck-together negatives from Tinker,
by matching decimals anchored to EXACTLY three fractional digits
(``[-+]?[0-9]+\\.[0-9]{3}``), which corresponds to the PDB ``.3f`` convention and
therefore never bleeds one coordinate into the next.
"""

import argparse
import os
import re
import sys
import tempfile

WATER_RESNAMES_TO_WAT = ("001", "002", "SOL")

# PDB coordinates are written as %8.3f -> always exactly 3 fractional digits.
# Anchoring the decimals to {3} is what disambiguates touching fixed-width fields
# (0.131|5005.730) and stuck-together negatives (1132.340|-665.460|-4394.330).
_COORD_RE = re.compile(r"[-+]?[0-9]+\.[0-9]{3}")


def pad_left(s, width):
    """Right-justify to width; if longer than width keep the rightmost width."""
    if len(s) < width:
        return " " * (width - len(s)) + s
    if len(s) > width:
        return s[len(s) - width:]
    return s


def pad_right(s, width):
    """Left-justify to width; truncate to width."""
    if len(s) < width:
        return s + " " * (width - len(s))
    return s[:width]


def atom_field(name):
    """Format an atom name into the canonical 4-character columns 13-16."""
    name = name.strip(" ")
    if name == "":
        name = "O"
    n = len(name)
    if n >= 4:
        return name[:4]
    if n == 1:
        return " " + name + "  "
    if n == 2:
        return " " + name + " "
    if n == 3:
        return " " + name
    return pad_right(name, 4)


def parse_coords(region):
    """Return the first three x/y/z decimals from a coordinate region, or None.

    ``region`` is the substring of an ATOM/HETATM line from column 31 onward.
    Returns a list ``[x, y, z]`` of the raw numeric strings, or None if fewer
    than three three-decimal numbers are present.
    """
    nums = _COORD_RE.findall(region)
    if len(nums) >= 3:
        return nums[:3]
    return None


def format_atom_line(line):
    """Reformat one ATOM/HETATM record into canonical columns.

    Returns ``(new_line, changes)`` where ``changes`` is a list of applied
    transformations (e.g. ``["OW->O"]``).
    """
    record = line[0:6]
    serial_raw = line[6:11]
    atom = line[12:16]
    altloc = line[16:17]
    res = line[17:20]
    chain = line[21:22]
    resseq = line[22:26]
    icode = line[26:27]

    coords = parse_coords(line[30:])
    if coords is not None:
        x, y, z = coords
        # A standard-width line has z in the canonical columns 47-54. If the
        # fixed-column z matches the parsed z, occupancy/tempFactor/tail are
        # trustworthy too. On a widened (shifted) Tinker frame they are not --
        # and such frames carry no occ/temp/tail -- so blank them out rather
        # than copying shifted bytes.
        if line[46:54].strip(" ") == z:
            occ = line[54:60]
            temp = line[60:66]
            tail = line[66:] if len(line) >= 67 else ""
        else:
            occ = temp = tail = ""
    else:
        # Fewer than three parseable coordinates: fall back to raw fixed slices.
        x = line[30:38]
        y = line[38:46]
        z = line[46:54]
        occ = line[54:60]
        temp = line[60:66]
        tail = line[66:] if len(line) >= 67 else ""

    atom_trim = atom.strip(" ")
    res_trim = res.strip(" ")

    changes = []
    if atom_trim == "OW":
        atom_trim = "O"
        changes.append("OW->O")
    if res_trim in WATER_RESNAMES_TO_WAT:
        res_trim = "WAT"
        changes.append("resname->WAT")

    new_line = (
        pad_right(record, 6)
        + pad_left(serial_raw.strip(" "), 5)
        + " "
        + atom_field(atom_trim)
        + (altloc if altloc != "" else " ")
        + pad_right(res_trim, 3)
        + " "
        + (chain if chain != "" else " ")
        + pad_left(resseq.strip(" "), 4)
        + (icode if icode != "" else " ")
        + "   "
        + pad_left(x.strip(" "), 8)
        + pad_left(y.strip(" "), 8)
        + pad_left(z.strip(" "), 8)
        + pad_left(occ.strip(" "), 6)
        + pad_left(temp.strip(" "), 6)
        + tail
    )
    return new_line, changes


def _new_stats():
    return {
        "total_input_rows": 0,
        "output_atom_rows": 0,
        "removed_conect": 0,
        "changed_rows": 0,
        "ow_to_o": 0,
        "res_to_wat": 0,
    }


def _summary_lines(stats):
    total_changed = stats["removed_conect"] + stats["changed_rows"]
    return [
        "",
        "Summary",
        "-------",
        "Total input rows        : %d" % stats["total_input_rows"],
        "Output atom rows        : %d" % stats["output_atom_rows"],
        "Removed CONECT rows     : %d" % stats["removed_conect"],
        "Changed output rows     : %d" % stats["changed_rows"],
        "OW -> O changes         : %d" % stats["ow_to_o"],
        "001/002 -> WAT changes  : %d" % stats["res_to_wat"],
        "Total rows changed      : %d" % total_changed,
    ]


def format_record(line, lineno, stats):
    """Process one input line.

    Returns ``(out_line_or_None, log_entries)`` and mutates ``stats``.
    ``out_line`` is None when the line is dropped (CONECT). ``log_entries`` is a
    (possibly empty) list of log lines to emit for this row.
    """
    stats["total_input_rows"] += 1

    if line.startswith("CONECT"):
        stats["removed_conect"] += 1
        return None, [
            "[REMOVED] input_line=%d, atom_serial=-: row starts with CONECT" % lineno,
            "  before: " + line,
        ]

    record = line[0:6]
    if record not in ("ATOM  ", "HETATM"):
        return line, []

    new_line, changes = format_atom_line(line)
    stats["output_atom_rows"] += 1

    if "OW->O" in changes:
        stats["ow_to_o"] += 1
    if "resname->WAT" in changes:
        stats["res_to_wat"] += 1

    entries = []
    if changes or new_line != line:
        stats["changed_rows"] += 1
        # Only emit a verbose before/after record for SUBSTANTIVE edits
        # (OW->O, resname->WAT). Pure column normalization touches nearly every
        # atom line, so logging it would bloat the log to hundreds of MB on a
        # multi-model trajectory (e.g. Tinker's 1000-frame system_mdl.pdb).
        if changes:
            serial = line[6:11].strip(" ")
            entries.append(
                "[CHANGED] input_line=%d, atom_serial=%s: %s"
                % (lineno, serial, "; ".join(changes))
            )
            entries.append("  before: " + line)
            entries.append("  after : " + new_line)

    return new_line, entries


def _stream(lines, out_fh, log_fh):
    """Stream-process ``lines``, writing output and log incrementally.

    Neither the whole input nor the whole output/log is held in memory, so this
    stays O(1) in memory on large multi-model trajectory PDBs.
    """
    log_fh.write("PDB rewrite log\n===============\n")
    stats = _new_stats()
    for i, raw in enumerate(lines, start=1):
        line = raw.rstrip("\r\n")
        out_line, entries = format_record(line, i, stats)
        if out_line is not None:
            out_fh.write(out_line + "\n")
        for entry in entries:
            log_fh.write(entry + "\n")
    for summary in _summary_lines(stats):
        log_fh.write(summary + "\n")
    return stats


def process(lines):
    """In-memory convenience wrapper (used by tests).

    Returns ``(out_lines, log_lines, stats)``. Production code uses ``_stream``
    against file handles instead, to avoid buffering the whole file.
    """
    out_lines = []
    log_lines = ["PDB rewrite log", "==============="]
    stats = _new_stats()
    for i, raw in enumerate(lines, start=1):
        out_line, entries = format_record(raw.rstrip("\r\n"), i, stats)
        if out_line is not None:
            out_lines.append(out_line)
        log_lines.extend(entries)
    log_lines.extend(_summary_lines(stats))
    return out_lines, log_lines, stats


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Rewrite a PDB into canonical fixed-width columns "
        "(remove CONECT, OW->O, 001/002/SOL->WAT).",
    )
    parser.add_argument("input", metavar="INPUT_PDB", help="PDB file to rewrite in place")
    parser.add_argument(
        "logfile",
        nargs="?",
        default=None,
        help="log file (default: <input_stem>.log next to the input)",
    )
    args = parser.parse_args(argv)

    input_path = args.input
    if not os.path.isfile(input_path):
        print("Error: input file not found: %s" % input_path, file=sys.stderr)
        return 1

    input_dir = os.path.dirname(input_path) or "."
    input_name = os.path.basename(input_path)
    if input_name.endswith(".pdb"):
        stem = input_name[: -len(".pdb")]
        backup = os.path.join(input_dir, stem + ".original.pdb")
    else:
        stem = input_name
        backup = input_path + ".original.pdb"

    logfile = args.logfile if args.logfile else os.path.join(input_dir, stem + ".log")

    fd, tmp_path = tempfile.mkstemp(dir=input_dir)
    try:
        # Stream input -> temp output and log line-by-line (the temp output is a
        # separate file, so reading the input while writing is safe). The input
        # is only renamed to the backup after both are fully written and closed.
        with os.fdopen(fd, "w") as out_fh, open(logfile, "w") as log_fh, open(
            input_path, "r", errors="ignore"
        ) as in_fh:
            _stream(in_fh, out_fh, log_fh)

        if os.path.exists(backup):
            os.remove(backup)
        os.replace(input_path, backup)
        os.replace(tmp_path, input_path)
    except Exception:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        raise

    print("Wrote cleaned PDB to: %s" % input_path)
    print("Original PDB backed up to: %s" % backup)
    print("Wrote log to: %s" % logfile)
    return 0


if __name__ == "__main__":
    sys.exit(main())
