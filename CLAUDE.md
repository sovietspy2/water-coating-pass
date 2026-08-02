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
./pipeline/wdrop.sh <INPUT_PDB> <MODE> [REFERENCE_PDB] [--layers L] [--intermediate-md-ns NS] [--final-md-ns NS]
# MODE: gromacs | tinker
# --layers L: number of water layers (default 5); one layer is deposited+minimized per iteration, so L is also the iteration count
# --intermediate-md-ns NS: MD length in ns run after each intermediate iteration (default 0); 0 = no intermediate MD (minimize only)
# --final-md-ns NS: MD length in ns for the final iteration (default 0.1); must be > 0
```

- Iterations are no longer a separate flag: `--layers L` deposits exactly one layer per iteration, so the run is always `L` deposit+minimize iterations, feeding each minimized output into the next → output directory `<base>_l<L>_int<INT>_fin<FIN>/{1..L}/` (single directory when `L == 1`), where `<INT>`/`<FIN>` are the `--intermediate-md-ns`/`--final-md-ns` values.
- **`--intermediate-md-ns` is the single knob for intermediate MD**: intermediate iterations run MD iff `--intermediate-md-ns > 0`. `0` (the default) = intermediate iterations are minimize-only.
- The final iteration always runs MD (`--final-md-ns`, must be `> 0`) + MobyWat, regardless of the intermediate setting.
- **MobyWat runs strictly on the final iteration only.**
- Backends receive two independent controls: MD is gated on the per-cycle duration being `> 0`, and MobyWat is gated on a separate `run_mobywat` flag (`1` only on the final iteration). This lets intermediate iterations run MD without MobyWat.

The input PDB **must be in its own working directory** outside the project folder.

## Testing

Tests live in per-purpose subdirectories under `testing/`:

```bash
# Load/integration test (downloads a real PDB from RCSB, runs wdrop.sh repeatedly):
./testing/load-test/test.sh [OPTIONS]
# -m gromacs|tinker|all   -t 1|5|all   -n <count>   -p <path>   -l (loop forever)

# Determinism test (verifies the C program output is bit-identical across N runs):
./testing/algo-deterministic-test/wdrop_deterministic_test.sh [-n <runs>]

# Engine determinism suites — how much do repeated identical runs differ?
./testing/minimize-deterministic-test/run.sh [-e gromacs|tinker|all] [-n N] [--arms pipeline|hardened|all]
./testing/md-deterministic-test/run.sh [-e gromacs|tinker|tinker-gpu|all] [-n N] [--md-ns NS] [--arms ...]

# format_pdb.py unit tests / pdb-preprocessor.py tests (pytest):
python testing/format-pdb-test/test_format_pdb.py
python testing/pdb-preprocess-test/test_pdb_preprocessor.py

# Dockerized end-to-end run:
./testing/e2e-test/start.sh [OPTIONS]

# Quick local smoke wrappers (each makes a tmp dir, wgets 1R6J, runs one combo):
./testing/fast/{tinker,gromacs,tinker_gpu}.sh
```

### Determinism suites (`testing/{minimize,md}-deterministic-test/`)

These two suites quantify the stochastic component documented in `TINKER-DETERMINISM.md`.
Unlike everything else under `testing/`, they drive the **engine binaries directly**
(`minimize`, `dynamic`, `dynamic9`, `gmx mdrun`) instead of `pipeline/*.sh` — the exact
commands, key files and mdp blocks are copied out of `pipeline/tinker.sh` and
`pipeline/gromacs.sh` into `testing/determinism-common/determinism_lib.sh`, each with a
source-line comment. **Changing a key or mdp setting in the pipeline means changing it there
too**, or the suites stop measuring the shipped configuration.

- Both suites share `testing/determinism-common/`: `determinism_lib.sh` (fixture build, key/mdp
  writers, work dir, report headers), `compare_structures.py` (RMSD / max deviation split by
  protein vs water, most-variable atoms, divergence-vs-frame), `parse_scalars.py` (engine
  energies and iteration counts).
- The fixture is built once per run — `pdb-preprocessor.py --target` → `wdrop --layers 1` →
  `format_pdb.py` — and copied byte-identically into each replicate directory. The input PDB
  (`testing/test_pdbs/1PSV_cryst.pdb`) is never modified; its sha256 is verified on exit.
- **Precision is load-bearing.** GROMACS comparisons go through `gmx trjconv -o traj.g96`
  (`%15.9f` nm); `.gro`'s `%8.3f` nm = 0.01 Å cannot represent the 1e-6 Å divergence these
  suites look for. Tinker MD additionally hashes `.dyn` (16 digits) because `.arc` has only 6.
- **`--md-ns` must stay above ~0.002 ns.** Divergence is invisible below ~500 MD steps, so a
  shorter run reports a false pass. Default is 0.02 ns = 20 000 steps at the pipeline's 1 fs.
- Each case runs one or two *arms*: `pipeline` (as shipped) and `hardened`
  (`OPENMP-THREADS 1` for Tinker, `gmx mdrun -reprod -ntmpi 1 -ntomp 1` for GROMACS).
- `tinker-gpu` reports SKIPPED and exits 0 when `dynamic9` is absent, so `-e all` still completes.
- Reports land in `testing/<suite>/results/<UTC-timestamp>_<case>/{report.md,results.csv}`, plus a
  `summary_<timestamp>.md` per driver run. Results are only comparable **within one machine**;
  each report records host/CPU/threads/GPU.

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

- **`pipeline_common.sh`** — sourced by all pipeline scripts. Provides `log()`, `run_step()`, `setup_logging()`, `normalize_input_pdb()`, `make_output_dir()`, `validate_script_dir_not_input_dir()`, `run_mm_step()`, plus the MobyWat target-selection helpers `protein_chain_ids()`, `mobywat_target_spec()` (see "MobyWat target selection" below).
- **`wdrop.sh`** — orchestrates the per-layer deposit+minimize loop (`--layers` iterations, one layer each): model reduction → PDB fix (Python) → iterative wdrop + MM/MD → file collection. `--intermediate-md-ns` sets the intermediate MD length (`0` = none, the default) and `--final-md-ns` the final one; the output-dir name encodes both (`_l<L>_int<INT>_fin<FIN>`).
- **`gromacs.sh`** / **`tinker.sh`** — backend-specific MM+MD steps; write `next_step.pdb` on success. `tinker.sh` supports Tinker9-GPU: export `TINKER_GPU=1` to run the MD step with `dynamic9` instead of `dynamic` (file-format utilities `pdbxyz`/`arcedit`/`xyzpdb` are unchanged); GPU mode also injects `INTEGRATOR VERLET` + `REMOVE-INERTIA 0` into `md.key`. Both backends also honor `MOBYWAT_DEBUG=1`, which runs an extra MobyWat pass in Analysis mode (`-m Analysis`) after the Prediction one; Analysis mode needs a reference structure, so it is skipped in prediction-only runs. It writes the per-frame success rates to `O_system_succ.txt` at Silent verbosity, leaving the Prediction outputs (`O_system_mt*.lst`, `O_system_rmst.txt`) untouched.
- **`pdb-preprocessor.py`** — fixes missing residues/atoms and nonstandard residues in X-ray PDBs (uses pdbfixer). Mode is selected by a required flag: `--target` strips all waters and heterogens (used by `wdrop.sh`); `--reference` keeps waters but removes other heterogens (used by `gromacs.sh`/`tinker.sh` for MobyWat validation).
- **`format_pdb.py`** — rewrites a PDB into canonical fixed-width columns (removes CONECT, renames OW→O, maps residue names 001/002/SOL→WAT), backing the original up to `<stem>.original.pdb`. Called by both backends before/after MD. Its coordinate parser anchors decimals to exactly 3 fractional digits (`%8.3f`), so it correctly reads touching fixed-width fields (GROMACS pseudo-PBC coords ~5000 Å) and widened/stuck fields (Tinker `xyzpdb` `%9.3f` frames) alike. Per-row change detail in the rewrite log is capped at `_DETAIL_LOG_CAP` (1000 rows) so multi-thousand-frame trajectories don't produce huge logs; summary counts remain exact. Tests in `testing/format-pdb-test/`.
- **`research.sh`** — batch benchmark + CSV history generator. Takes `<INPUT_DIR> <OUTPUT_DIR>`, and for every `.pdb` in the input dir runs each combo in its `COMBOS` array, each in its own working dir. Each `COMBOS` entry is a full `wdrop.sh` invocation — the engine (`tinker`|`gromacs`) followed by any `wdrop.sh` flags (`--layers`, `--intermediate-md-ns`, `--final-md-ns`); omitted flags fall back to `wdrop.sh` defaults. The combo string is slugged into a unique per-combo working dir, and the analysis dir is derived from the resolved `--layers`/`--intermediate-md-ns`/`--final-md-ns` (mirroring wdrop.sh's `_l<L>_int<INT>_fin<FIN>` tag). Each run appends **one row per PDB × combo** to an **append-only, committed CSV history** — `testing/research_history.csv` by default, overridable via the `RESEARCH_CSV` env var. The file is never truncated (it's the long-term performance record for later DB import); the header is written only when the file is missing/empty. Columns: `id` (auto-increment across runs), `pdb_name`, `engine`, `layers`, `intermediate_md_ns`, `final_md_ns`, `total_runtime_seconds`, `rmsd_avg`/`rmsd_max`, `sr_mer`/`sr_ida`/`sr_ide`/`sr_pos`, `sr_frame_min`/`sr_frame_max`/`sr_frame_avg`/`sr_frame_std`, `failure` (`1` when any of the six result metrics is empty = no usable result, else `0`), `run_timestamp` (ISO 8601 UTC), `commit_hash`, and `cpu`/`gpu` (from the `CPU`/`GPU` env vars, `UNKNOWN` when unset). The `sr_frame_*` columns summarize the **per-frame** success rates in `O_system_succ.txt` (min, max, mean, population standard deviation) — that file is only written by MobyWat's Analysis mode, so `research.sh` exports `MOBYWAT_DEBUG=1` by default (set `MOBYWAT_DEBUG=0` to opt out, leaving the four columns empty). They are deliberately not part of the `failure` test, which keeps its original six-metric meaning across the whole history. `total_runtime_seconds` **excludes** the `MOBYWAT_DEBUG` Analysis pass: each backend logs `MobyWat Analysis runtime: N seconds` and `append_csv_row()` subtracts it from `wdrop.sh`'s `Total runtime was` figure, so the column reflects pipeline work rather than the diagnostic extra (rows recorded before this change still include it). Missing metrics are left empty. Downloads nothing (files are provided/reviewed manually); honors `TINKER_GPU=1`.
- **`add-mobywat-analysis-params.sh`** / **`remove-ter-id.sh`** — small helpers for the reference PDB: the first prepends MobyWat `REMARK mobywat_*` tuning params, taking the target spec as an optional `$2` (`add-mobywat-analysis-params.sh <pdb> [AB]`, default `[A]`) which it writes into `REMARK mobywat_reference_target`; the second normalizes `TER` records and renumbers atoms.

### MobyWat target selection (chain IDs)

MobyWat picks the target molecule by chain ID (`-t x-y/[xy...]`), and that selection appears in **three** places that must all name the same protein:

1. the `-t` switch — the target inside the trajectory,
2. `REMARK mobywat_reference_target` in `system_ref.pdb` — the target inside the reference,
3. the chain IDs physically present in the trajectory topology (`system_tpy.pdb` for GROMACS, `system_mdl.pdb` for Tinker).

All three used to be pinned to `[A]`, and `gromacs.sh` additionally ran `gmx editconf -label A`, which stamped chain `A` on *every* atom. On a multi-chain protein that made the trajectory's `[A]` the whole multimer while the reference's `[A]` was only the first chain — Prediction mode tolerated it (`O_system_rmst.txt` is frame-vs-first-frame RMSD), but Analysis mode aborted with `Numbers of CA atoms are not equal ...` **while still exiting 0**, so the pipeline reported success and `O_system_succ.txt`/the `sr_frame_*` CSV columns were silently empty.

Both are documented: MobyWat manual §4.3.2 defines `–t [xy…]` precisely for when *"the target is stored in several chains"*, and §3.1.5.2 says the reference REMARK exists because *"ranges in the reference file may be different from those of the trajectory"* — so each spec is derived from the file it addresses, and the two need not be spelled the same. `mobywat_target_spec()` returns e.g. `[AB]` from a PDB's CA-bearing residues (waters and heterogens drop out; only the first `MODEL` is read). `-t` is quoted at the call sites, since `[AB]` would otherwise glob against files named `A` or `B`.

Single-chain PDBs still resolve to `[A]`, so their results are unchanged; multi-chain PDBs now validate against every chain's surface waters (16GS: 245 reference waters instead of 133), which makes their `sr_*` values **not comparable to rows recorded before this change**.

Two things to know:

- **`gromacs.sh: ensure_chain_labels()`** — `pdb2gmx` only assigns chain letters when the system has 2+ chains, so a **monomer comes out of `md.tpr` with a blank chain ID** that no `[xy...]` target can select. That is what the old unconditional `editconf -label A` was for: the manual's worked GROMACS recipe (§9.1.2.2, the source of this block) says a chain ID *"can be added"* that way, and GROMACS documents `-label` as *"Add chain label for **all residues**"* — a remedy for a missing ID, not something to apply always. It now runs only when the topology has no chain IDs.
- **Nothing verifies that a MobyWat pass produced its outputs.** MobyWat exits 0 even when it aborts, so a run reports success either way and `O_system.log` is where the reason lives — note the Analysis pass overwrites that log, so only the last pass's copy survives. Without a reference MobyWat writes no `O_system_mt*.lst` at all, so empty `sr_*` columns are expected for a prediction-only run. The Analysis pass itself is a `MOBYWAT_DEBUG` extra and is non-fatal, so a finished MD + Prediction run is never discarded over it (`if …; then` bodies keep `errexit`, hence the invocation sits in the `if` condition). The pipeline assumes each run starts in a clean directory and is never resumed.

The pipeline auto-creates `.venv` in the repo root if absent, activates it, and uses it for all Python calls. All shell scripts use `set -euo pipefail`.

`python/` and `research_test/` at the repo root are experimental/scratch areas (e.g. `python/pdbxyz.py`, a WIP pure-Python `pdbxyz` replacement with hardcoded paths) and are **not** wired into the pipeline — the pipeline uses the real Tinker `pdbxyz` binary.

### Shell coding conventions

- Variable names: `SCREAMING_SNAKE_CASE` throughout.
- Local variables declared with `local`.
- All file paths validated before use with `require_file` or explicit `[[ -f ... ]]`.
