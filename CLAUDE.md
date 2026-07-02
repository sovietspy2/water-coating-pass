# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

WDROP is a molecular dynamics preprocessing pipeline. The C program (`wdrop`) places water molecule oxygens (HOH) around a protein surface using a PASS-like sphere coating algorithm, then the pipeline scripts drive molecular mechanics (MM) and molecular dynamics (MD) minimization using GROMACS or Tinker as backends.

## Build

```bash
cd src && make          # produces src/wdrop
make install            # copies wdrop to ~/bin
```

CI runs `make` from `src/` on every push/PR (`.github/workflows/c-cpp.yml`).

## Running wdrop directly

```bash
wdrop --file <input.pdb> --sigma <float> --weed-dist <float> --layers <int>
```

Output is written as `<input_base>_<N>WAT.pdb` in the same directory as the input. A log entry is also written alongside it.

## Running the full pipeline

```bash
./pipeline/wdrop.sh <INPUT_PDB> <MODE> [REFERENCE_PDB] [--iterations N] [--layers L]
# MODE: gromacs | tinker
# --iterations N: number of deposit+minimize cycles (default 1)
# --layers L: total water layers across the run (default 5); each cycle deposits L/N
```

- Default (`--iterations 1`, `--layers 5`): all 5 layers in one wdrop run + one MM+MD step → output directory `<base>_i1_l5/`
- `--iterations 5`: 5 iterative cycles (1 layer each), feeding each minimized output into the next; only the final cycle runs MD + MobyWat → output directory `<base>_i5_l5/{1..5}/`
- `--layers` must be an exact multiple of `--iterations` (per-cycle layers = `L / N`).

The input PDB **must be in its own working directory** outside the project folder.

## Testing

```bash
# Full integration test (downloads a real PDB from RCSB, runs wdrop.sh):
./testing/test.sh [OPTIONS]
# -m gromacs|tinker|all   -t 1|5|all   -n <count>   -l (loop forever)

# Determinism test (verifies the C program output is bit-identical across N runs):
./testing/pass_deterministic_test.sh [-n <runs>]
```

## Installation of dependencies

```bash
cd installer && ./install.sh
```

Bootstraps: Python 3, pdbfixer (in `.venv`), GROMACS, Tinker, wdrop, MobyWat. Idempotent.

## Code architecture

### C program (`src/`)

- **`macros.h`** — all compile-time constants and helper print macros. PDB field widths, default tolerances, program name strings.
- **`global.h`** — all PDB-related structs: `ap` (atom), `tp` (topology), `pd` (pairwise distance), `hs`/`hd` (hydrogen-bond network edge/node), `re`/`pl`/`pf` (pool/occupancy lists), `wr`/`mr` (water matching). Struct fields reflect the PDB column layout.
- **`input.c`/`input.h`** — `read_in_pdb()` parses PDB files into `ap[]` arrays; `parse_program_args()` handles CLI flags (`--file`, `--sigma`, `--weed-dist`, `--layers`) into `WdropProgramOptions`.
- **`wdrop.c`/`wdrop.h`** — core algorithm. `three_point_sphere_geometry()` computes candidate water positions from every triple of nearby atom spheres; `pass_like_coating()` iterates over atom triples, filters clashes, and appends accepted HOH atoms. Uses a linked-list spatial grid for neighbor lookup. `sigma_from_atom()` maps element symbols to van-der-Waals radii.
- **`print.c`/`print.h`** — PDB file output and per-run log entry.
- **`main.c`** — entry point: parse args → read PDB → `pass_like_coating()` → write output PDB + log.

Global state lives in `main.c` as `g_`-prefixed variables (e.g. `g_pdb_ref`, `g_ref_name`). The structs in `global.h` carry extra fields inherited from MobyWat that are not used by `wdrop` itself.

### Pipeline scripts (`pipeline/`)

- **`pipeline_common.sh`** — sourced by all pipeline scripts. Provides `log()`, `run_step()`, `setup_logging()`, `normalize_input_pdb()`, `make_output_dir()`, `validate_script_dir_not_input_dir()`, `run_mm_step()`.
- **`wdrop.sh`** — orchestrates the `--iterations` deposit+minimize loop: model reduction → PDB fix (Python) → iterative wdrop + MM/MD → file collection.
- **`gromacs.sh`** / **`tinker.sh`** — backend-specific MM+MD steps; write `next_step.pdb` on success.
- **`pdb-preprocessor.py`** — fixes missing residues/atoms and nonstandard residues in X-ray PDBs (uses pdbfixer). Mode is selected by a required flag: `--target` strips all waters and heterogens (used by `wdrop.sh`); `--reference` keeps waters but removes other heterogens (used by `gromacs.sh`/`tinker.sh` for MobyWat validation).
- **`remove_unbound_water.py`** — reference PDB utility for MobyWat validation mode.

The pipeline auto-creates `.venv` in the repo root if absent, activates it, and uses it for all Python calls. All shell scripts use `set -euo pipefail`.

### Shell coding conventions

- Variable names: `SCREAMING_SNAKE_CASE` throughout.
- Local variables declared with `local`.
- All file paths validated before use with `require_file` or explicit `[[ -f ... ]]`.
