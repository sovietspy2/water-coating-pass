#!/usr/bin/env python3
"""Quantify how much a set of supposedly-identical structure files differ.

Given N output files produced by N identical runs of the same engine command, this
reports whether they are bit-identical and, when they are not, *how far apart* they
are: pairwise RMSD and maximum per-atom deviation, split into protein and water
atoms, plus the atoms that vary most across replicates.

Everything is reported in angstrom. GROMACS formats store nanometres and are scaled
on read, so Tinker and GROMACS numbers are directly comparable.

Supported formats
-----------------
``tinker-xyz``   Tinker cartesian coordinates (``%12.6f`` A). Also reads ``.arc``
                 archives -- pass ``--trajectory`` to treat the extra frames as a
                 time series instead of erroring.
``g96``          GROMOS96 (``%15.9f`` nm). This is what ``gmx trjconv -o x.g96``
                 writes and it is the only text format GROMACS emits with enough
                 digits to see early divergence.
``gro``          GROMACS ``.gro`` (``%8.3f`` nm = 0.01 A resolution). Too coarse for
                 a determinism verdict; supported for eyeballing only.

Precision matters here: TINKER-DETERMINISM.md section 2 measured divergence of
1e-6 A at MD step 1010, which a ``.gro`` file cannot represent at all.
"""

import argparse
import csv
import hashlib
import itertools
import json
import math
import os
import sys

import numpy as np

# Residue names that mean "water" across the formats the pipeline produces.
# format_pdb.py maps 001/002/SOL to WAT; GROMACS writes SOL; PDB standard is HOH.
WATER_RESNAMES = {"WAT", "SOL", "HOH", "TIP3", "T3P", "H2O"}

# Tinker atom names assigned by pdbxyz for amber99 water. pipeline/tinker.sh:113
# uses exactly these two to find where the protein ends.
WATER_TINKER_NAMES = {"OW", "HW"}

NM_TO_ANGSTROM = 10.0


class ParseError(RuntimeError):
    pass


class Frame:
    """One set of coordinates plus the per-atom labels needed to group them."""

    def __init__(self, coords, names, resnames):
        self.coords = np.asarray(coords, dtype=np.float64)
        self.names = names
        self.resnames = resnames

    def __len__(self):
        return len(self.coords)


# ---------------------------------------------------------------------------
# Readers
# ---------------------------------------------------------------------------


def _is_int(token):
    try:
        int(token)
        return True
    except ValueError:
        return False


def read_tinker(path):
    """Read a Tinker .xyz / .xyz_2 / .arc file into a list of Frames.

    Layout per frame:
        <natoms>  <title>
        [box line]                      -- present iff its first token is not an int
        <idx> <name> <x> <y> <z> <type> <connections...>

    The optional box line is detected the same way pipeline/tinker.sh:23-34 does it.
    """
    with open(path) as handle:
        lines = handle.read().splitlines()

    frames = []
    cursor = 0
    total = len(lines)

    while cursor < total:
        header = lines[cursor].split()
        if not header:
            cursor += 1
            continue
        if not _is_int(header[0]):
            raise ParseError(f"{path}: expected an atom count at line {cursor + 1}, got {lines[cursor]!r}")

        natoms = int(header[0])
        cursor += 1

        # Optional box line: a real atom record always starts with its integer index.
        if cursor < total:
            first_token = lines[cursor].split()[:1]
            if first_token and not _is_int(first_token[0]):
                cursor += 1

        if cursor + natoms > total:
            raise ParseError(f"{path}: truncated frame at line {cursor + 1} (wanted {natoms} atoms)")

        coords = np.empty((natoms, 3), dtype=np.float64)
        names = []
        for i in range(natoms):
            fields = lines[cursor + i].split()
            if len(fields) < 5:
                raise ParseError(f"{path}: short atom record at line {cursor + i + 1}")
            names.append(fields[1])
            coords[i] = (float(fields[2]), float(fields[3]), float(fields[4]))
        cursor += natoms

        # Tinker has no residue names; the atom name carries the water identity.
        resnames = ["WAT" if n in WATER_TINKER_NAMES else "PRO" for n in names]
        frames.append(Frame(coords, names, resnames))

    if not frames:
        raise ParseError(f"{path}: no frames found")
    return frames


def read_g96(path):
    """Read a GROMOS96 file (nm, %15.9f) into a list of Frames.

    Only the POSITION / POSITIONRED blocks are used. POSITION records are
    ``resnr resname atomname atomnr x y z``; POSITIONRED is bare ``x y z``.
    """
    frames = []
    coords = None
    names = None
    resnames = None
    block = None

    with open(path) as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            stripped = line.strip()

            if block is None:
                if stripped in ("POSITION", "POSITIONRED"):
                    block = stripped
                    coords, names, resnames = [], [], []
                continue

            if stripped == "END":
                frames.append(Frame(np.asarray(coords) * NM_TO_ANGSTROM, names, resnames))
                block = None
                continue

            fields = stripped.split()
            if block == "POSITIONRED":
                if len(fields) < 3:
                    raise ParseError(f"{path}: short POSITIONRED record: {stripped!r}")
                coords.append([float(fields[-3]), float(fields[-2]), float(fields[-1])])
                names.append("?")
                resnames.append("?")
            else:
                if len(fields) < 7:
                    raise ParseError(f"{path}: short POSITION record: {stripped!r}")
                resnames.append(fields[1])
                names.append(fields[2])
                coords.append([float(fields[-3]), float(fields[-2]), float(fields[-1])])

    if not frames:
        raise ParseError(f"{path}: no POSITION block found")
    return frames


def read_gro(path):
    """Read a (possibly multi-frame) GROMACS .gro file into a list of Frames.

    Fixed layout: %5d resnr, %5s resname, %5s atomname, %5d atomnr, then three
    coordinates. Field width is 8 by default but grows with -ndec, so the
    coordinates are taken by splitting the tail rather than by column.
    """
    with open(path) as handle:
        lines = handle.read().splitlines()

    frames = []
    cursor = 0
    total = len(lines)

    while cursor < total:
        if not lines[cursor].strip():
            cursor += 1
            continue
        cursor += 1  # title
        if cursor >= total:
            break
        natoms = int(lines[cursor].split()[0])
        cursor += 1

        if cursor + natoms > total:
            raise ParseError(f"{path}: truncated frame at line {cursor + 1}")

        coords = np.empty((natoms, 3), dtype=np.float64)
        names = []
        resnames = []
        for i in range(natoms):
            line = lines[cursor + i]
            resnames.append(line[5:10].strip())
            names.append(line[10:15].strip())
            tail = line[20:].split()
            if len(tail) < 3:
                raise ParseError(f"{path}: short atom record at line {cursor + i + 1}")
            coords[i] = (float(tail[0]), float(tail[1]), float(tail[2]))
        cursor += natoms + 1  # atoms + box line

        frames.append(Frame(coords * NM_TO_ANGSTROM, names, resnames))

    if not frames:
        raise ParseError(f"{path}: no frames found")
    return frames


READERS = {"tinker-xyz": read_tinker, "g96": read_g96, "gro": read_gro}


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------


def water_mask(frame):
    """Boolean mask selecting water atoms."""
    by_res = np.array([r.upper() in WATER_RESNAMES for r in frame.resnames])
    by_name = np.array([n.upper() in WATER_TINKER_NAMES for n in frame.names])
    return by_res | by_name


def deviations(a, b):
    """Per-atom euclidean distance between two coordinate sets, in angstrom."""
    return np.linalg.norm(a - b, axis=1)


def rmsd(dev):
    return float(math.sqrt(float(np.mean(dev ** 2)))) if dev.size else 0.0


def group_stats(dev, mask):
    """(rmsd, max deviation) over the masked subset."""
    subset = dev[mask]
    if subset.size == 0:
        return None, None
    return rmsd(subset), float(np.max(subset))


def coord_fingerprint(frame):
    """Hash of the coordinate array -- identifies structurally distinct results
    even when the surrounding file text (timestamps, titles) differs."""
    return hashlib.sha256(np.ascontiguousarray(frame.coords).tobytes()).hexdigest()


def file_sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sample_frame_indices(count, wanted):
    """Log-spaced frame indices, so early (tiny) divergence and late (saturated)
    divergence both appear in the table without printing 1000 rows."""
    if count <= wanted:
        return list(range(count))
    raw = np.unique(np.round(np.geomspace(1, count, num=wanted)).astype(int)) - 1
    return [int(i) for i in raw]


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------


def read_labels(path):
    """Read an external per-atom label file: one ``<resname> <atomname>`` line per atom.

    Needed because ``gmx trjconv -o x.g96`` emits POSITIONRED blocks -- bare
    coordinates with no atom or residue names -- so the protein/water split has to
    come from the topology instead of from the compared file.
    """
    resnames, names = [], []
    with open(path) as handle:
        for line in handle:
            fields = line.split()
            if not fields:
                continue
            resnames.append(fields[0])
            names.append(fields[1] if len(fields) > 1 else "?")
    return resnames, names


def analyse(paths, fmt, trajectory, frame_samples, labels_path=None):
    reader = READERS[fmt]

    replicates = []
    for path in paths:
        frames = reader(path)
        if not trajectory and len(frames) > 1:
            # A single-structure comparison only ever wants the final state.
            frames = frames[-1:]
        replicates.append({"path": path, "frames": frames, "sha256": file_sha256(path)})

    natoms = len(replicates[0]["frames"][-1])
    for rep in replicates[1:]:
        if len(rep["frames"][-1]) != natoms:
            raise ParseError(
                "replicates have different atom counts "
                f"({natoms} in {replicates[0]['path']}, "
                f"{len(rep['frames'][-1])} in {rep['path']}); "
                "the fixture was not shared correctly"
            )

    reference = replicates[0]["frames"][-1]
    if labels_path:
        resnames, names = read_labels(labels_path)
        if len(resnames) != natoms:
            raise ParseError(
                f"{labels_path} describes {len(resnames)} atoms but the structures have {natoms}"
            )
        reference.resnames = resnames
        reference.names = names
    is_water = water_mask(reference)
    is_protein = ~is_water

    result = {
        "format": fmt,
        "n_replicates": len(replicates),
        "n_atoms": int(natoms),
        "n_water_atoms": int(np.count_nonzero(is_water)),
        "n_protein_atoms": int(np.count_nonzero(is_protein)),
        "files": [{"path": r["path"], "sha256": r["sha256"]} for r in replicates],
    }

    # --- bit-level ---------------------------------------------------------
    result["n_distinct_files"] = len({r["sha256"] for r in replicates})
    result["n_distinct_structures"] = len({coord_fingerprint(r["frames"][-1]) for r in replicates})
    result["bit_identical"] = result["n_distinct_files"] == 1

    # --- pairwise on the final frame ---------------------------------------
    pairs = []
    for (i, a), (j, b) in itertools.combinations(list(enumerate(replicates)), 2):
        dev = deviations(a["frames"][-1].coords, b["frames"][-1].coords)
        prot_rmsd, prot_max = group_stats(dev, is_protein)
        wat_rmsd, wat_max = group_stats(dev, is_water)
        pairs.append({
            "rep_a": i + 1,
            "rep_b": j + 1,
            "identical": a["sha256"] == b["sha256"],
            "rmsd_all": rmsd(dev),
            "max_dev_all": float(np.max(dev)),
            "rmsd_protein": prot_rmsd,
            "max_dev_protein": prot_max,
            "rmsd_water": wat_rmsd,
            "max_dev_water": wat_max,
        })
    result["pairs"] = pairs

    def spread(key):
        values = [p[key] for p in pairs if p[key] is not None]
        if not values:
            return None
        return {
            "mean": float(np.mean(values)),
            "sd": float(np.std(values, ddof=1)) if len(values) > 1 else 0.0,
            "max": float(np.max(values)),
        }

    result["summary"] = {key: spread(key) for key in (
        "rmsd_all", "max_dev_all", "rmsd_protein", "max_dev_protein", "rmsd_water", "max_dev_water"
    )}

    # --- per-atom variability ----------------------------------------------
    stack = np.stack([r["frames"][-1].coords for r in replicates])          # (N, atoms, 3)
    centroid = stack.mean(axis=0)
    per_atom_sd = np.sqrt(((np.linalg.norm(stack - centroid, axis=2)) ** 2).mean(axis=0))
    order = np.argsort(per_atom_sd)[::-1][:10]
    result["most_variable_atoms"] = [
        {
            "index": int(k) + 1,
            "name": reference.names[k],
            "resname": reference.resnames[k],
            "is_water": bool(is_water[k]),
            "sd_angstrom": float(per_atom_sd[k]),
        }
        for k in order
    ]

    # --- divergence over the trajectory ------------------------------------
    if trajectory:
        frame_counts = {len(r["frames"]) for r in replicates}
        if len(frame_counts) > 1:
            result["divergence_warning"] = (
                f"replicates have different frame counts: {sorted(frame_counts)}; "
                "comparing only the first " + str(min(frame_counts))
            )
        usable = min(frame_counts)
        curve = []
        for idx in sample_frame_indices(usable, frame_samples):
            worst_rmsd = 0.0
            worst_max = 0.0
            for a, b in itertools.combinations(replicates, 2):
                dev = deviations(a["frames"][idx].coords, b["frames"][idx].coords)
                worst_rmsd = max(worst_rmsd, rmsd(dev))
                worst_max = max(worst_max, float(np.max(dev)))
            curve.append({"frame": idx + 1, "max_pairwise_rmsd": worst_rmsd, "max_pairwise_dev": worst_max})
        result["divergence_curve"] = curve

    return result


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------


def fmt_num(value, digits=6):
    if value is None:
        return "-"
    if value == 0:
        return "0"
    if abs(value) < 1e-4:
        return f"{value:.2e}"
    return f"{value:.{digits}f}".rstrip("0").rstrip(".")


def verdict(result):
    """One of: BIT-IDENTICAL / NUMERICALLY IDENTICAL / DIVERGENT."""
    if result["bit_identical"]:
        return "BIT-IDENTICAL"
    worst = result["summary"]["max_dev_all"]
    if worst is None or worst["max"] == 0.0:
        return "NUMERICALLY IDENTICAL"
    return "DIVERGENT"


def headline(result, label):
    summary = result["summary"]
    parts = [
        f"{label}: {verdict(result)}",
        f"{result['n_distinct_structures']}/{result['n_replicates']} distinct structures",
    ]
    if not result["bit_identical"]:
        parts.append(f"max pairwise RMSD {fmt_num(summary['rmsd_all']['max'])} A")
        parts.append(f"max per-atom dev {fmt_num(summary['max_dev_all']['max'])} A")
        if summary["rmsd_protein"] and summary["rmsd_water"]:
            parts.append(
                f"(protein {fmt_num(summary['rmsd_protein']['max'])} A / "
                f"water {fmt_num(summary['rmsd_water']['max'])} A)"
            )
    return " | ".join(parts)


def render_markdown(result, label):
    summary = result["summary"]
    out = [
        f"### {label}",
        "",
        f"**Verdict: {verdict(result)}** — "
        f"{result['n_distinct_structures']} distinct structure(s) among "
        f"{result['n_replicates']} identical runs "
        f"({result['n_distinct_files']} distinct file hash(es)).",
        "",
        f"System: {result['n_atoms']} atoms "
        f"({result['n_protein_atoms']} protein, {result['n_water_atoms']} water). "
        f"Format: `{result['format']}`. All distances in angstrom.",
        "",
        "| metric | mean over pairs | SD | max |",
        "|---|---|---|---|",
    ]
    labels = {
        "rmsd_all": "RMSD, all atoms",
        "max_dev_all": "max per-atom deviation",
        "rmsd_protein": "RMSD, protein atoms",
        "max_dev_protein": "max deviation, protein",
        "rmsd_water": "RMSD, water atoms",
        "max_dev_water": "max deviation, water",
    }
    for key, text in labels.items():
        stat = summary.get(key)
        if stat is None:
            out.append(f"| {text} | - | - | - |")
        else:
            out.append(f"| {text} | {fmt_num(stat['mean'])} | {fmt_num(stat['sd'])} | {fmt_num(stat['max'])} |")

    out += ["", "Pairwise detail:", "", "| pair | identical | RMSD all | max dev | RMSD protein | RMSD water |",
            "|---|---|---|---|---|---|"]
    for pair in result["pairs"]:
        out.append(
            f"| {pair['rep_a']} vs {pair['rep_b']} "
            f"| {'yes' if pair['identical'] else 'no'} "
            f"| {fmt_num(pair['rmsd_all'])} "
            f"| {fmt_num(pair['max_dev_all'])} "
            f"| {fmt_num(pair['rmsd_protein'])} "
            f"| {fmt_num(pair['rmsd_water'])} |"
        )

    if not result["bit_identical"]:
        out += ["", "Most variable atoms (SD of position across replicates):", "",
                "| atom # | name | residue | water? | SD |", "|---|---|---|---|---|"]
        for atom in result["most_variable_atoms"]:
            out.append(
                f"| {atom['index']} | {atom['name']} | {atom['resname']} "
                f"| {'yes' if atom['is_water'] else 'no'} | {fmt_num(atom['sd_angstrom'])} |"
            )

    if "divergence_warning" in result:
        out += ["", f"> WARNING: {result['divergence_warning']}"]

    if "divergence_curve" in result:
        out += ["", "Divergence over the trajectory (worst pair at each sampled frame):", "",
                "| frame | max pairwise RMSD | max pairwise per-atom dev |", "|---|---|---|"]
        for point in result["divergence_curve"]:
            out.append(
                f"| {point['frame']} | {fmt_num(point['max_pairwise_rmsd'])} "
                f"| {fmt_num(point['max_pairwise_dev'])} |"
            )

    out += ["", "Output file hashes:", ""]
    for entry in result["files"]:
        # Replicates all produce the same basename, so keep the parent dir for context.
        head, tail = os.path.split(entry["path"])
        shown = os.path.join(os.path.basename(head), tail) if head else tail
        out.append(f"- `{shown}` — `{entry['sha256'][:16]}`")
    out.append("")

    return "\n".join(out)


def write_csv(result, label, path):
    write_header = not os.path.exists(path) or os.path.getsize(path) == 0
    with open(path, "a", newline="") as handle:
        writer = csv.writer(handle)
        if write_header:
            writer.writerow([
                "case", "rep_a", "rep_b", "identical", "rmsd_all", "max_dev_all",
                "rmsd_protein", "max_dev_protein", "rmsd_water", "max_dev_water",
            ])
        for pair in result["pairs"]:
            writer.writerow([
                label, pair["rep_a"], pair["rep_b"], int(pair["identical"]),
                pair["rmsd_all"], pair["max_dev_all"],
                pair["rmsd_protein"], pair["max_dev_protein"],
                pair["rmsd_water"], pair["max_dev_water"],
            ])


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("files", nargs="+", help="structure files from identical runs (>= 2)")
    parser.add_argument("--format", required=True, choices=sorted(READERS),
                        help="input format")
    parser.add_argument("--label", required=True, help="case name used in the report and CSV")
    parser.add_argument("--labels",
                        help="per-atom '<resname> <atomname>' file overriding the labels in the "
                             "compared files; required for GROMACS, whose g96 output carries none")
    parser.add_argument("--trajectory", action="store_true",
                        help="treat files as multi-frame and emit a divergence-vs-frame table")
    parser.add_argument("--frame-samples", type=int, default=20,
                        help="how many frames to sample for the divergence table (default: 20)")
    parser.add_argument("--out-md", help="append a markdown section to this file")
    parser.add_argument("--out-csv", help="append pairwise rows to this CSV")
    parser.add_argument("--out-json", help="write the full result as JSON")
    return parser.parse_args()


def main():
    args = parse_args()

    if len(args.files) < 2:
        print("ERROR: need at least two files to compare", file=sys.stderr)
        return 2

    missing = [p for p in args.files if not os.path.isfile(p)]
    if missing:
        print("ERROR: missing output file(s): " + ", ".join(missing), file=sys.stderr)
        return 2

    try:
        result = analyse(args.files, args.format, args.trajectory, args.frame_samples, args.labels)
    except ParseError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    result["verdict"] = verdict(result)
    result["headline"] = headline(result, args.label)

    if args.out_json:
        with open(args.out_json, "w") as handle:
            json.dump(result, handle, indent=2)

    if args.out_md:
        with open(args.out_md, "a") as handle:
            handle.write(render_markdown(result, args.label) + "\n")

    if args.out_csv:
        write_csv(result, args.label, args.out_csv)

    print(result["headline"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
