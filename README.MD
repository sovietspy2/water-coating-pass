# WDROP - Water Droplet-based Molecular Dynamics Pipeline via MobyWat

## Software components

The pipeline ties together a C binary, Bash orchestration, a Python virtual
environment, and external MM/MD engines. Every box below must be available:

```text
   wdrop.sh   (Bash orchestrator; sources pipeline_common.sh)
      |
      |   loops N cycles:   deposit water  -->  minimize
      |   the final cycle additionally runs   MD  -->  MobyWat
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
      +-- [4] validate ------> mobywat   (runs on the final cycle)
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
| `mobywat` | Water prediction/validation on the final cycle | external install |

## Compile
- `cd src` and then `make`

## Install
- `make install`

## Pipeline scripts

The project includes user-facing pipeline script in `pipeline/`:

# wdrop.sh

`wdrop.sh` is a Bash pipeline wrapper for running `wdrop` and molecular mechanics minimization on a PDB file. It deposits water layers around the protein and relaxes the structure over one or more deposit-and-minimize cycles, controlled by the `--iterations` and `--layers` flags.

## Usage

```bash
./wdrop.sh <INPUT_PDB> <MODE> [REFERENCE_PDB] [--iterations N] [--layers L]
```

Examples:

```bash
# Default: 5 water layers deposited in a single cycle
./wdrop.sh /tmp/working_dir/1PSV.pdb gromacs
```

```bash
# 5 iterative cycles (1 layer per cycle)
./wdrop.sh /tmp/working_dir/1PSV.pdb gromacs --iterations 5
```

```bash
# Validation mode with a reference PDB (waters already present)
./wdrop.sh /tmp/working_dir/1PSV.pdb gromacs /tmp/working_dir/1PSV_cryst.pdb --iterations 5
```

Note: Place the input PDB in its own working directory located outside the project folder. The working directory should contain only that PDB file.
## Arguments

- `INPUT_PDB` — Absolute or relative path to an input `.pdb` file.
- `MODE` — Minimization backend, either `gromacs` or `tinker`.
- `REFERENCE_PDB` — Absolute or relative path to a `.pdb` file used as the reference structure (it already contains waters), enabling MobyWat validation mode. Optional.
- `--iterations N` — Number of deposit+minimize cycles (default `1`).
- `--layers L` — Total number of `wdrop` water layers deposited across the whole run (default `5`). Each cycle deposits `L / N` layers, so `L` must be an exact multiple of `N`.

## Iterations and layers

Each run deposits `--layers` water layers in total, split evenly across `--iterations` deposit-and-minimize cycles (`--layers / --iterations` layers per cycle).

- **Intermediate cycles** run `wdrop` followed by energy minimization only — **no MD**. Each cycle's minimized structure is stored in a numbered subdirectory and fed into the next cycle.
- **The final cycle** runs `wdrop`, minimization, an MD simulation, and MobyWat.

Examples:

- Default (`--iterations 1`, `--layers 5`): one `wdrop` run of 5 layers, then minimize + MD + MobyWat.
- `--iterations 5` (layers default 5): 5 cycles of 1 layer each; cycles 1–4 are minimize-only, cycle 5 adds MD + MobyWat.

## Output

- Log file: `application.LOG` in the input file directory.
- Output directory: `<input_base>_i<iterations>_l<layers>` (e.g. `1PSV_i1_l5`, `1PSV_i5_l5`).
- For multi-cycle runs (`--iterations > 1`), each cycle is stored in a numbered subdirectory: `<input_base>_i5_l5/1` ... `<input_base>_i5_l5/5`.

## Workflow

1. Validate arguments and normalize the input PDB path.
2. Fix the input PDB (add missing atoms/residues), strip existing waters, and reduce it to a single model.
3. For each cycle: run `wdrop` to deposit water layers, then minimize with the selected backend. The final cycle also runs MD and MobyWat.
4. Collect and organize each cycle's output files into its run directory.
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
- `--iterations` and `--layers` may appear in any order, before or after the positional arguments.
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






