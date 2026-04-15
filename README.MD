# PASS like coating

## Compile
- `cd src` and then `make`

## Install
- `make install`

## Pipeline scripts

The project includes user-facing pipeline script in `pipeline/`:

# pass.sh

`pass.sh` is a Bash pipeline wrapper for running `pass` and molecular mechanics minimization on a PDB file. It supports two execution modes: a single short run or a five-step iterative long run.

## Usage

```bash
./pass.sh <INPUT_PDB> <MODE> <RUN_TYPE>
```

Example:

```bash
./pass.sh /abs/path/FFF.pdb gromacs LONG
```

Note: The input PDB must be placed in a separate work directory outside of the project. 

## Arguments

- `INPUT_PDB` — Absolute or relative path to an input `.pdb` file.
- `MODE` — Minimization backend, either `gromacs` or `tinker`.
- `RUN_TYPE` — Pipeline mode, either `SHORT` or `LONG`.

## Run types

### SHORT

`SHORT` performs one `pass` run with 5 water layers, then one MM and MD step. The generated files are collected into a single output directory named after the input file with the `5x1` suffix.

### LONG

`LONG` performs 5 iterative cycles. Each cycle runs `pass` with 1 layer of water, then MM then MD, stores outputs in a numbered subdirectory, and uses the minimized structure as the input for the next cycle.

## Output

- Log file: `application.LOG` in the input file directory.
- `SHORT` output directory: `<input_base>5x1`
- `LONG` output directory: `<input_base>5x5`
- `LONG` iteration subdirectories: `<input_base>5x5/1` ... `<input_base>5x5/5`

## Workflow

1. Validate arguments and normalize the input PDB path.
2. Remove extra PDB models, keeping only the first one.
3. Run `pass`.
4. Run minimization using the selected backend.
5. Save and organize output files.
6. Write progress and runtime information to `application.LOG`.

## Requirements

This script expects the following helper files in the same script directory:

- `pipeline_common.sh`
- `model-reducer.sh`
- `gromacs.sh`
- `tinker.sh`

It also expects the `pass` executable to be available in `PATH`.

## Notes

- The script uses strict Bash mode: `set -euo pipefail`.
- `RUN_TYPE` is case-insensitive because it is internally normalized to uppercase.
- Existing `application.LOG` is reused and appended to, not truncated.

## Exit behavior

The script exits immediately on error and writes a failure message to the log with the line number and exit code. On completion, it logs the total runtime in seconds.

## requirements:
- python
- linux
- C compiler
- tinker installed and on path
- gromacs installed and on path

## Acknowledgments

This project reuses code from **MobyWat** (http://www.mobywat.com).]()






