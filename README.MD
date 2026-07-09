# WDROP - Water Droplet-based Molecular Dynamics Pipeline via MobyWat

## Software components

The pipeline ties together a C binary, Bash orchestration, a Python virtual
environment, and external MM/MD engines. Every box below must be available:

```text
   wdrop.sh   (Bash orchestrator; sources pipeline_common.sh)
      |
      |   loops L iterations (one water layer each):  deposit water  -->  minimize
      |   with --intermediate-md-ns > 0 each iteration also runs a short MD;
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

`wdrop.sh` is a Bash pipeline wrapper for running `wdrop` and molecular mechanics minimization on a PDB file. It deposits water layers around the protein one layer at a time and relaxes the structure over a deposit-and-minimize iteration per layer, controlled by the `--layers`, `--intermediate-md-ns`, and `--final-md-ns` flags.

## Usage

```bash
./wdrop.sh <INPUT_PDB> <MODE> [REFERENCE_PDB] [--layers L] [--intermediate-md-ns NS] [--final-md-ns NS]
```

Examples:

```bash
# Default: 5 layers, one deposited+minimized per iteration; MD + MobyWat only on the last
./wdrop.sh /tmp/working_dir/1PSV.pdb gromacs
```

```bash
# Refine after every layer: a short MD (0.1 ns) runs each intermediate iteration, the final iteration runs 0.5 ns + MobyWat
./wdrop.sh /tmp/working_dir/1PSV.pdb tinker --intermediate-md-ns 0.1 --final-md-ns 0.5
```

```bash
# No intermediate MD (the default): only the last iteration runs MD (0.5 ns) + MobyWat
./wdrop.sh /tmp/working_dir/1PSV.pdb tinker --intermediate-md-ns 0 --final-md-ns 0.5
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
- `--intermediate-md-ns NS` — MD length in nanoseconds run after each intermediate iteration (default `0`). `0` means no intermediate MD (minimize only); any value `> 0` runs an MD of that length on every intermediate iteration. This single knob decides whether intermediate iterations run MD — there is no separate mode flag.
- `--final-md-ns NS` — MD length in nanoseconds for the final iteration (default `0.1`). Must be a positive number (`> 0`).

## Layers and MD

Each run deposits `--layers` water layers, one per deposit-and-minimize iteration. Each iteration's minimized structure is stored in a numbered subdirectory and fed into the next iteration.

- **Intermediate iterations** run `wdrop` + energy minimization, and additionally an MD of `--intermediate-md-ns` when that value is `> 0` (no MD when it is `0`, the default).
- **The final iteration** always runs `wdrop`, minimization, MD (`--final-md-ns`), and MobyWat.
- **MobyWat always runs strictly on the final iteration.**

The two MD lengths are controlled independently, so intermediate MD can be tuned or switched off (`--intermediate-md-ns 0`) without affecting the final production MD. The output directory name encodes both MD lengths (see Output below).

Examples:

- Default (`--layers 5`, `--intermediate-md-ns 0`): 5 iterations of 1 layer each; iterations 1–4 are minimize-only, iteration 5 adds MD + MobyWat.
- `--layers 5 --intermediate-md-ns 0.1`: iterations 1–4 run minimize + a short (0.1 ns) MD, iteration 5 runs the longer MD + MobyWat.

## Output

- Log file: `application.LOG` in the input file directory.
- Output directory: `<input_base>_l<layers>_int<intermediate_md_ns>_fin<final_md_ns>` — e.g. `1PSV_l5_int0_fin0.1`, `1PSV_l5_int0.1_fin0.5`.
- Each iteration is stored in a numbered subdirectory: `<input_base>_l5_int0_fin0.1/1` ... `<input_base>_l5_int0_fin0.1/5` (a single directory is used when `--layers 1`).

## Workflow

1. Validate arguments and normalize the input PDB path.
2. Fix the input PDB (add missing atoms/residues), strip existing waters, and reduce it to a single model.
3. For each iteration: run `wdrop` to deposit one water layer, then minimize with the selected backend. When `--intermediate-md-ns > 0`, each intermediate iteration also runs an MD; the final iteration always runs MD and MobyWat.
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
- `--layers`, `--intermediate-md-ns`, and `--final-md-ns` may appear in any order, before or after the positional arguments.
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






