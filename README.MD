# WDROP - Water Droplet-based Molecular Dynamics Pipeline via MobyWat

## Software components

The pipeline ties together a C binary, Bash orchestration, a Python virtual
environment, and external MM/MD engines. Every box below must be available:

```text
   wdrop.sh   (Bash orchestrator; sources pipeline_common.sh)
      |
      |   loops L iterations (one water layer each):  deposit water  -->  minimize
      |   per_layer refinement also runs a short MD each iteration;
      |   the final iteration always runs   MD  -->  MobyWat
      |
      +-- [1] preprocess ----> pdb-preprocessor.py ----> Python venv (.venv) + pdbfixer
      |
      +-- [2] place water ---> wdrop (C binary)   <---- built with: C compiler + make
      |
      +-- [3] minimize / MD --> (MODE selects one backend)
      |          gromacs ---> gromacs.sh ---> GROMACS (gmx)
      |          tinker  ---> tinker.sh  ---> Tinker CPU (minimize, dynamic,
      |                                                   pdbxyz, xyzpdb, arcedit)
      |                                   \
      |                                    `--> Tinker9-GPU (dynamic9)  [optional, TINKER_GPU=1]
      |
      +-- [4] validate ------> mobywat   (runs on the final iteration)
```

| Component | Role | Provided by |
| --- | --- | --- |
| `wdrop.sh` + `pipeline_common.sh` | Pipeline orchestration (Bash) | this repo (`pipeline/`) |
| `wdrop` | Places HOH water layers on the protein surface | this repo (`src/`, C binary) |
| C compiler + `make` | Builds the `wdrop` binary | system toolchain |
| `pdb-preprocessor.py` | Fixes/reduces the input PDB | this repo, needs Python venv |
| Python venv (`.venv`) + `pdbfixer` | PDB fixing dependency | auto-created by the pipeline |
| GROMACS (`gmx`) | MM + MD backend (`MODE=gromacs`) | external install |
| Tinker (`minimize`, `dynamic`, `pdbxyz`, `xyzpdb`, `arcedit`) | MM + MD backend (`MODE=tinker`) | external install |
| Tinker9-GPU (`dynamic9`) | GPU-accelerated MD (optional, `TINKER_GPU=1`) | external install |
| `mobywat` | Water prediction/validation on the final iteration | external install |

## Compile
- `cd src` and then `make`

## Install
- `make install`

## Pipeline scripts

The project includes user-facing pipeline script in `pipeline/`:

# wdrop.sh

`wdrop.sh` is a Bash pipeline wrapper for running `wdrop` and molecular mechanics minimization on a PDB file. It deposits water layers around the protein one layer at a time and relaxes the structure over a deposit-and-minimize iteration per layer, controlled by the `--layers` and `--refinement` flags.

## Usage

```bash
./wdrop.sh <INPUT_PDB> <MODE> [REFERENCE_PDB] [--layers L] [--refinement default|per_layer] [--intermediate-md-ns NS] [--final-md-ns NS]
```

Examples:

```bash
# Default: 5 layers, one deposited+minimized per iteration; MD + MobyWat only on the last
./wdrop.sh /tmp/working_dir/1PSV.pdb gromacs
```

```bash
# Refine after every layer: MD runs each iteration (0.1 ns), the final iteration runs 0.5 ns + MobyWat
./wdrop.sh /tmp/working_dir/1PSV.pdb tinker --refinement per_layer --intermediate-md-ns 0.1 --final-md-ns 0.5
```

```bash
# Turn intermediate MD off but keep per-layer minimization; only the last iteration runs MD (0.5 ns) + MobyWat
./wdrop.sh /tmp/working_dir/1PSV.pdb tinker --refinement per_layer --intermediate-md-ns 0 --final-md-ns 0.5
```

```bash
# Validation mode with a reference PDB (waters already present)
./wdrop.sh /tmp/working_dir/1PSV.pdb gromacs /tmp/working_dir/1PSV_cryst.pdb --layers 5
```

Note: Place the input PDB in its own working directory located outside the project folder. The working directory should contain only that PDB file.
## Arguments

- `INPUT_PDB` — Absolute or relative path to an input `.pdb` file.
- `MODE` — Minimization backend, either `gromacs` or `tinker`.
- `REFERENCE_PDB` — Absolute or relative path to a `.pdb` file used as the reference structure (it already contains waters), enabling MobyWat validation mode. Optional.
- `--layers L` — Number of `wdrop` water layers to deposit (default `5`). One layer is deposited and minimized per iteration, so `L` is also the iteration count.
- `--refinement MODE` — `default` runs MD only on the final iteration; `per_layer` runs MD after every iteration (default `default`).
- `--intermediate-md-ns NS` — MD length in nanoseconds for the intermediate iterations in `per_layer` mode (default `0.1`). May be `0` to disable intermediate MD (minimize only). Ignored in `default` mode.
- `--final-md-ns NS` — MD length in nanoseconds for the final iteration (default `0.5`). Must be a positive number (`> 0`).

## Layers, refinement, and MD

Each run deposits `--layers` water layers, one per deposit-and-minimize iteration. Each iteration's minimized structure is stored in a numbered subdirectory and fed into the next iteration.

- **`--refinement default`** — intermediate iterations run `wdrop` + energy minimization only (**no MD**); the final iteration runs `wdrop`, minimization, MD (`--final-md-ns`), and MobyWat.
- **`--refinement per_layer`** — every iteration additionally runs MD. Intermediate iterations use `--intermediate-md-ns` (skipped when it is `0`), and the final iteration uses `--final-md-ns`.
- **MobyWat always runs strictly on the final iteration**, in both refinement modes.

The two MD lengths are controlled independently, so intermediate MD can be tuned or switched off (`--intermediate-md-ns 0`) without affecting the final production MD.

Examples:

- Default (`--layers 5`, `--refinement default`): 5 iterations of 1 layer each; iterations 1–4 are minimize-only, iteration 5 adds MD + MobyWat.
- `--layers 5 --refinement per_layer`: iterations 1–4 run minimize + a short MD, iteration 5 runs the longer MD + MobyWat.

## Output

- Log file: `application.LOG` in the input file directory.
- Output directory: `<input_base>_l<layers>_<refinement>` (e.g. `1PSV_l5_default`, `1PSV_l5_per_layer`).
- Each iteration is stored in a numbered subdirectory: `<input_base>_l5_default/1` ... `<input_base>_l5_default/5` (a single directory is used when `--layers 1`).

## Workflow

1. Validate arguments and normalize the input PDB path.
2. Fix the input PDB (add missing atoms/residues), strip existing waters, and reduce it to a single model.
3. For each iteration: run `wdrop` to deposit one water layer, then minimize with the selected backend. In `per_layer` mode each iteration also runs an intermediate MD; the final iteration always runs MD and MobyWat.
4. Collect and organize each iteration's output files into its run directory.
5. Write progress and runtime information to `application.LOG`.

## Requirements

This script expects the following helper files in the same script directory:

- `pipeline_common.sh`
- `pdb-preprocessor.py` (PDB fixing, water stripping, single-model reduction)
- `gromacs.sh` / `tinker.sh` (MM/MD backends)
- `format_pdb.py` (used by the backends to canonicalize PDB columns)

It also expects the `wdrop` executable to be available in `PATH`.

## Notes

- The script uses strict Bash mode: `set -euo pipefail`.
- `--layers`, `--refinement`, `--intermediate-md-ns`, and `--final-md-ns` may appear in any order, before or after the positional arguments.
- Existing `application.LOG` is reused and appended to, not truncated.

## Exit behavior

The script exits immediately on error and writes a failure message to the log with the line number and exit code. On completion, it logs the total runtime in seconds.

## requirements:
- python with pdbfixer installed (`pip install pdbfixer`; the pipeline also auto-creates a `.venv` with it)
- linux or WSL
- C compiler
- tinker installed and on path
- gromacs installed and on path

## Understanding the output

### Gromacs mode
- Gromacs mode will produce .gro and .xtc files
- load both files to Pymol and run `zoom mobywat_input` command to find it
- use `show spheres` to see the water molecules

### Tinker mode
- Tinker mode will produce PDB file (multiple models)
- load the PDB file to Pymol
- use `show spheres` to see the water molecules


## Acknowledgments

This project reuses code from **MobyWat** http://www.mobywat.com.






