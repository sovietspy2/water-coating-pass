# minimize-deterministic-test

How much do two identical energy-minimisation runs differ?

Runs the minimisation step N times on byte-identical inputs and reports whether the
outputs are bit-identical — and when they are not, by how much, split into protein
and water atoms.

The engine binaries are driven directly rather than through `pipeline/*.sh`; see
[`../determinism-common/README.md`](../determinism-common/README.md) for why, and for
the fixture, precision and arm conventions shared with the MD suite.

## Usage

```bash
./run.sh                                  # both engines, all arms, 5 replicates
./run.sh -e tinker --arms all -n 3        # one engine
./run.sh -e gromacs --keep                # keep the temp work dir for inspection

./TINKER/run.sh  -n 5 --arms all          # a single case, directly
./GROMACS/run.sh -n 5 --arms all
```

Options: `-e gromacs|tinker|all`, `-n <replicates>` (default 5),
`-i <input.pdb>` (default `testing/test_pdbs/1PSV_cryst.pdb`),
`--arms pipeline|hardened|all` (default `all`), `--keep`.

## What each case runs

**TINKER** — `minimize sys.xyz -k minimize.key`, gradient criterion `0.01`
(`pipeline/tinker.sh:140`). One stage.

**GROMACS** — two stages, measured separately so you can see which one is stochastic
(`pipeline/gromacs.sh:203-213`):

- `em` — steepest descent, `gmx mdrun -v -s em -o em.trr -c after_em.gro -g em.log`
- `cg` — conjugate gradient, same command against `cg.tpr`

`grompp` runs once per stage, and the `cg` replicates all start from one canonical
`after_em.gro`, so `cg` variance cannot inherit `em` variance.

On a small system `cg` often converges in ≤1 step, because `em` already drove Fmax
below the `cg` `emtol` of 750. The report says so explicitly when it happens — a
stage that did no work is trivially deterministic and its verdict means nothing.

## Output

`results/<UTC-timestamp>_<engine>/`:

- `report.md` — verdict per arm, pairwise RMSD / max-deviation tables, protein-vs-water
  split, most variable atoms, and the engine's own scalars (final energy, iteration
  count) with their spread
- `results.csv` — one row per replicate pair
- `*_structures.json`, `*_scalars.json` — the same data unrounded

`results/summary_<timestamp>.md` collects the verdicts when `run.sh` runs more than
one case.

## Interpreting it

`iterations` and `function_value` matter as much as RMSD. L-BFGS amplifies a low-bit
force difference into a different search path, so two runs can land on structures that
are close in RMSD but took a visibly different number of iterations to get there — and
occasionally on genuinely different minima.

Because the protein is position-restrained (`RESTRAIN-POSITION ... 2.0` for Tinker,
`-DPOSRES` for GROMACS) and the waters are not, water atoms carry most of the
variance. The protein-vs-water split in the report makes that visible; a large water
RMSD next to a near-zero protein RMSD is the expected shape, not a bug.

Results are only comparable **within one machine** — the Tinker user guide states
reproducibility "will not hold across different machine types". Every report records
the host, CPU, thread count and GPU for that reason.
