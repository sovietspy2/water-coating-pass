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
./pipeline/wdrop.sh <INPUT_PDB> <MODE> [REFERENCE_PDB] [--layers L] [--refinement default|per_layer] [--intermediate-md-ns NS] [--final-md-ns NS]
# MODE: gromacs | tinker
# --layers L: number of water layers (default 5); one layer is deposited+minimized per iteration, so L is also the iteration count
# --refinement MODE: 'default' runs MD only on the final iteration; 'per_layer' runs MD after each iteration (default 'default')
# --intermediate-md-ns NS: MD length in ns for intermediate iterations in per_layer mode (default 0.1); 0 disables intermediate MD
# --final-md-ns NS: MD length in ns for the final iteration (default 0.5); must be > 0
```

- Iterations are no longer a separate flag: `--layers L` deposits exactly one layer per iteration, so the run is always `L` deposit+minimize iterations, feeding each minimized output into the next → output directory `<base>_l<L>_<refinement>/{1..L}/` (single directory when `L == 1`).
- `--refinement default` (`--layers 5`): 5 iterations of 1 layer each; iterations 1–4 are minimize-only, iteration 5 adds MD (`--final-md-ns`) + MobyWat.
- `--refinement per_layer`: every iteration runs MD; intermediate iterations use `--intermediate-md-ns` (skipped when it is `0`), the final iteration uses `--final-md-ns`.
- **MobyWat runs strictly on the final iteration only**, regardless of refinement mode.
- Backends receive two independent controls: MD is gated on the per-cycle duration being `> 0`, and MobyWat is gated on a separate `run_mobywat` flag (`1` only on the final iteration). This lets per_layer intermediate iterations run MD without MobyWat.

The input PDB **must be in its own working directory** outside the project folder.

## Testing

Tests live in per-purpose subdirectories under `testing/`:

```bash
# Load/integration test (downloads a real PDB from RCSB, runs wdrop.sh repeatedly):
./testing/load-test/test.sh [OPTIONS]
# -m gromacs|tinker|all   -t 1|5|all   -n <count>   -p <path>   -l (loop forever)

# Determinism test (verifies the C program output is bit-identical across N runs):
./testing/algo-deterministic-test/wdrop_deterministic_test.sh [-n <runs>]

# format_pdb.py unit tests / pdb-preprocessor.py tests (pytest):
python testing/format-pdb-test/test_format_pdb.py
python testing/pdb-preprocess-test/test_pdb_preprocessor.py

# Dockerized end-to-end run:
./testing/e2e-test/start.sh [OPTIONS]

# Quick local smoke wrappers (each makes a tmp dir, wgets 1R6J, runs one combo):
./testing/fast/{tinker,gromacs,tinker_gpu}.sh
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
- **`wdrop.sh`** — orchestrates the per-layer deposit+minimize loop (`--layers` iterations, one layer each): model reduction → PDB fix (Python) → iterative wdrop + MM/MD → file collection. `--refinement` selects whether MD runs every iteration or only the last; `--intermediate-md-ns`/`--final-md-ns` set the two MD lengths.
- **`gromacs.sh`** / **`tinker.sh`** — backend-specific MM+MD steps; write `next_step.pdb` on success. `tinker.sh` supports Tinker9-GPU: export `TINKER_GPU=1` to run the MD step with `dynamic9` instead of `dynamic` (file-format utilities `pdbxyz`/`arcedit`/`xyzpdb` are unchanged); GPU mode also injects `INTEGRATOR VERLET` + `REMOVE-INERTIA 0` into `md.key`.
- **`pdb-preprocessor.py`** — fixes missing residues/atoms and nonstandard residues in X-ray PDBs (uses pdbfixer). Mode is selected by a required flag: `--target` strips all waters and heterogens (used by `wdrop.sh`); `--reference` keeps waters but removes other heterogens (used by `gromacs.sh`/`tinker.sh` for MobyWat validation).
- **`format_pdb.py`** — rewrites a PDB into canonical fixed-width columns (removes CONECT, renames OW→O, maps residue names 001/002/SOL→WAT), backing the original up to `<stem>.original.pdb`. Called by both backends before/after MD. Its coordinate parser anchors decimals to exactly 3 fractional digits (`%8.3f`), so it correctly reads touching fixed-width fields (GROMACS pseudo-PBC coords ~5000 Å) and widened/stuck fields (Tinker `xyzpdb` `%9.3f` frames) alike. Per-row change detail in the rewrite log is capped at `_DETAIL_LOG_CAP` (1000 rows) so multi-thousand-frame trajectories don't produce huge logs; summary counts remain exact. Tests in `testing/format-pdb-test/`.
- **`research.sh`** — batch benchmark + report generator. Takes `<INPUT_DIR> <OUTPUT_DIR>`, and for every `.pdb` in the input dir runs all four pipeline combos (`tinker i1`, `tinker i5`, `gromacs i1`, `gromacs i5`), each in its own working dir, distilling runtime, MobyWat SR values, and C-alpha RMSD into a shared `research_report.md`. Downloads nothing (files are provided/reviewed manually); honors `TINKER_GPU=1`.
- **`add-mobywat-analysis-params.sh`** / **`remove-ter-id.sh`** — small helpers for the reference PDB: the first prepends MobyWat `REMARK mobywat_*` tuning params; the second normalizes `TER` records and renumbers atoms.

The pipeline auto-creates `.venv` in the repo root if absent, activates it, and uses it for all Python calls. All shell scripts use `set -euo pipefail`.

`python/` and `research_test/` at the repo root are experimental/scratch areas (e.g. `python/pdbxyz.py`, a WIP pure-Python `pdbxyz` replacement with hardcoded paths) and are **not** wired into the pipeline — the pipeline uses the real Tinker `pdbxyz` binary.

### Shell coding conventions

- Variable names: `SCREAMING_SNAKE_CASE` throughout.
- Local variables declared with `local`.
- All file paths validated before use with `require_file` or explicit `[[ -f ... ]]`.
