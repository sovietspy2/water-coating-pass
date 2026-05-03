#!/usr/bin/env python3

import os
import sys
import shutil

HYDROGEN_NAMES = {"H", "HW", "HW1", "HW2", "H1", "H2"}
OXYGEN_NAMES = {"OW", "O", "OH2"}
COORD_TOL = 1e-6


def is_hydrogen(name):
    return name.strip().upper() in HYDROGEN_NAMES


def is_oxygen(name):
    return name.strip().upper() in OXYGEN_NAMES


def coord_key(x, y, z, tol=COORD_TOL):
    return (
        int(round(x / tol)),
        int(round(y / tol)),
        int(round(z / tol)),
    )


def is_float_token(tok):
    try:
        float(tok)
        return True
    except ValueError:
        return False


def parse_header(line):
    parts = line.rstrip("\n").split(None, 1)
    if not parts:
        raise ValueError("Empty first line")
    natoms = int(parts[0])
    title = parts[1] if len(parts) > 1 else ""
    return natoms, title


def looks_like_box_line(line):
    parts = line.strip().split()
    if len(parts) != 6:
        return False
    return all(is_float_token(p) for p in parts)


def parse_atom_line(line, lineno):
    parts = line.strip().split()
    if len(parts) < 6:
        raise ValueError(
            "Line {}: expected at least 6 fields, got {}: {!r}".format(
                lineno, len(parts), line
            )
        )

    try:
        old_id = int(parts[0])
        name = parts[1]
        x = float(parts[2])
        y = float(parts[3])
        z = float(parts[4])
        atom_type = int(parts[5])
        neighbors = [int(p) for p in parts[6:]]
    except ValueError as exc:
        raise ValueError("Line {}: parse error: {}".format(lineno, exc))

    return {
        "old_id": old_id,
        "name": name,
        "x": x,
        "y": y,
        "z": z,
        "atom_type": atom_type,
        "neighbors": neighbors,
    }


def format_atom_line(new_id, atom, new_neighbors):
    line = "{:6d}  {:<4s}{:14.6f}{:12.6f}{:12.6f}{:6d}".format(
        new_id,
        atom["name"],
        atom["x"],
        atom["y"],
        atom["z"],
        atom["atom_type"],
    )
    for nbr in new_neighbors:
        line += "{:6d}".format(nbr)
    return line.rstrip()


def choose_backup_path(path):
    candidate = path + ".bak"
    if not os.path.exists(candidate):
        return candidate

    i = 2
    while True:
        candidate = path + ".bak_{}".format(i)
        if not os.path.exists(candidate):
            return candidate
        i += 1


def load_tinker_xyz(path):
    with open(path, "r", encoding="utf-8") as fh:
        lines = fh.read().splitlines()

    if not lines:
        raise ValueError("Input file is empty")

    natoms, title = parse_header(lines[0])

    box_line = None
    atom_start = 1

    if len(lines) > 1 and looks_like_box_line(lines[1]):
        box_line = lines[1].rstrip("\n")
        atom_start = 2

    atom_lines = [line for line in lines[atom_start:] if line.strip()]
    atoms = []
    for idx, line in enumerate(atom_lines, start=atom_start + 1):
        atoms.append(parse_atom_line(line, idx))

    if natoms != len(atoms):
        raise ValueError(
            "Header says {} atoms, but parsed {} atom lines".format(
                natoms, len(atoms)
            )
        )

    return title, box_line, atoms


def build_id_map(atoms):
    return {atom["old_id"]: atom for atom in atoms}


def get_water_cluster_from_h(atom_id, id_map):
    atom = id_map.get(atom_id)
    if atom is None or not is_hydrogen(atom["name"]):
        return None

    oxygen_candidates = []
    for nbr in atom["neighbors"]:
        nbr_atom = id_map.get(nbr)
        if nbr_atom is not None and is_oxygen(nbr_atom["name"]):
            oxygen_candidates.append(nbr)

    if not oxygen_candidates:
        return None

    oxygen_id = oxygen_candidates[0]
    oxygen_atom = id_map[oxygen_id]

    h_neighbors = []
    heavy_neighbors = []

    for nbr in oxygen_atom["neighbors"]:
        nbr_atom = id_map.get(nbr)
        if nbr_atom is None:
            continue
        if is_hydrogen(nbr_atom["name"]):
            h_neighbors.append(nbr)
        else:
            heavy_neighbors.append(nbr)

    if heavy_neighbors:
        return None

    if not h_neighbors:
        return None

    cluster = set([oxygen_id] + h_neighbors)
    return cluster


def find_water_residues_to_remove(atoms):
    id_map = build_id_map(atoms)
    remove_ids = set()
    removed_water_oxygen_ids = set()
    removed_messages = []
    unresolved_messages = []

    first_seen_from_bottom = {}

    for idx in range(len(atoms) - 1, -1, -1):
        atom = atoms[idx]
        atom_id = atom["old_id"]
        key = coord_key(atom["x"], atom["y"], atom["z"])

        if key not in first_seen_from_bottom:
            first_seen_from_bottom[key] = atom_id
            continue

        lower_id = first_seen_from_bottom[key]
        lower_atom = id_map[lower_id]

        lower_cluster = get_water_cluster_from_h(lower_id, id_map)
        current_cluster = get_water_cluster_from_h(atom_id, id_map)

        if lower_cluster:
            oxygen_id = None
            for cid in lower_cluster:
                if is_oxygen(id_map[cid]["name"]):
                    oxygen_id = cid
                    break

            if oxygen_id not in removed_water_oxygen_ids:
                remove_ids.update(lower_cluster)
                removed_water_oxygen_ids.add(oxygen_id)
                removed_messages.append(
                    "Removed lower water centered at atom {} because duplicate coordinate "
                    "was detected at atom {} ({}) matching atom {} ({}), coords=({:.6f}, {:.6f}, {:.6f}); "
                    "removed atom IDs={}".format(
                        oxygen_id,
                        lower_id, lower_atom["name"],
                        atom_id, atom["name"],
                        atom["x"], atom["y"], atom["z"],
                        sorted(lower_cluster)
                    )
                )
            first_seen_from_bottom[key] = atom_id
            continue

        if current_cluster:
            oxygen_id = None
            for cid in current_cluster:
                if is_oxygen(id_map[cid]["name"]):
                    oxygen_id = cid
                    break

            if oxygen_id not in removed_water_oxygen_ids:
                remove_ids.update(current_cluster)
                removed_water_oxygen_ids.add(oxygen_id)
                removed_messages.append(
                    "Removed upper water centered at atom {} because duplicate coordinate "
                    "was detected at atom {} ({}) matching lower atom {} ({}), coords=({:.6f}, {:.6f}, {:.6f}); "
                    "removed atom IDs={}".format(
                        oxygen_id,
                        atom_id, atom["name"],
                        lower_id, lower_atom["name"],
                        atom["x"], atom["y"], atom["z"],
                        sorted(current_cluster)
                    )
                )
            continue

        unresolved_messages.append(
            "Duplicate retained: atom {} ({}) and atom {} ({}) share ({:.6f}, {:.6f}, {:.6f})".format(
                atom_id, atom["name"],
                lower_id, lower_atom["name"],
                atom["x"], atom["y"], atom["z"]
            )
        )

    return remove_ids, removed_messages, unresolved_messages


def rebuild_atoms(atoms, remove_ids):
    kept_atoms = [atom for atom in atoms if atom["old_id"] not in remove_ids]

    old_to_new = {}
    for i, atom in enumerate(kept_atoms, start=1):
        old_to_new[atom["old_id"]] = i

    new_lines = []
    for atom in kept_atoms:
        new_neighbors = []
        for nbr in atom["neighbors"]:
            if nbr in old_to_new and nbr != atom["old_id"]:
                new_neighbors.append(old_to_new[nbr])
        new_lines.append(format_atom_line(old_to_new[atom["old_id"]], atom, new_neighbors))

    return kept_atoms, new_lines


def write_fixed_file(path, title, box_line, atom_lines):
    header = str(len(atom_lines))
    if title:
        header += "  " + title

    tmp_path = path + ".tmp_fix"
    with open(tmp_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(header + "\n")
        if box_line is not None:
            fh.write(box_line + "\n")
        for line in atom_lines:
            fh.write(line + "\n")

    os.replace(tmp_path, path)


def main():
    if len(sys.argv) != 2:
        print("Usage: {} /full/path/to/file.xyz".format(os.path.basename(sys.argv[0])), file=sys.stderr)
        return 1

    xyz_path = os.path.abspath(os.path.expanduser(sys.argv[1]))

    if not os.path.exists(xyz_path):
        print("ERROR: file not found: {}".format(xyz_path), file=sys.stderr)
        return 1
    if not os.path.isfile(xyz_path):
        print("ERROR: not a regular file: {}".format(xyz_path), file=sys.stderr)
        return 1

    try:
        title, box_line, atoms = load_tinker_xyz(xyz_path)
    except Exception as exc:
        print("ERROR while reading {}: {}".format(xyz_path, exc), file=sys.stderr)
        return 1

    backup_path = choose_backup_path(xyz_path)
    shutil.copy2(xyz_path, backup_path)

    remove_ids, removed_messages, unresolved_messages = find_water_residues_to_remove(atoms)
    kept_atoms, new_atom_lines = rebuild_atoms(atoms, remove_ids)
    write_fixed_file(xyz_path, title, box_line, new_atom_lines)

    print("Backup created: {}".format(backup_path))
    print("Updated original: {}".format(xyz_path))
    print("Original atom count: {}".format(len(atoms)))
    print("New atom count: {}".format(len(kept_atoms)))
    print("Removed atoms: {}".format(len(remove_ids)))

    if removed_messages:
        print("\nRemoved waters:")
        for msg in removed_messages:
            print("  - " + msg)
    else:
        print("\nNo water residues were removed.")

    if unresolved_messages:
        print("\nUnresolved duplicates:")
        for msg in unresolved_messages:
            print("  - " + msg)

    return 0


if __name__ == "__main__":
    sys.exit(main())
