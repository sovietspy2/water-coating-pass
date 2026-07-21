# md-deterministic-test

How much do two identical MD runs differ, and when does the difference appear?

Runs the molecular-dynamics step N times from a byte-identical starting state and
reports how far the trajectories diverge — including a divergence-vs-frame table, so
you can see the point at which two runs stop agreeing rather than only the end state.

The engine binaries are driven directly rather than through `pipeline/*.sh`; see
[`../determinism-common/README.md`](../determinism-common/README.md) for why, and for
the fixture, precision and arm conventions shared with the minimize suite.

## The three cases

| case | binary | runs on |
|---|---|---|
| `gromacs` | `gmx mdrun` | CPU |
| `tinker` | `dynamic` | CPU |
| `tinker-gpu` | `dynamic9` | GPU |

`tinker-gpu` reports **SKIPPED** and exits 0 when `dynamic9` is not in `PATH`, so
`run.sh -e all` still completes on a machine without a Tinker9-GPU build.

## Usage

```bash
./run.sh                                       # all three cases, 5 replicates, 0.02 ns
./run.sh -e tinker -n 3 --md-ns 0.005 --arms all
./run.sh -e tinker-gpu                         # skips cleanly without a GPU build

./TINKER/run.sh                                # a single case, directly
./TINKER/run.sh --gpu
./GROMACS/run.sh
```

Options: `-e gromacs|tinker|tinker-gpu|all`, `-n <replicates>` (default 5),
`--md-ns <ns>` (default `0.02`), `-i <input.pdb>`,
`--arms pipeline|hardened|all` (default `pipeline` — MD is expensive), `--keep`.

## MD length — do not shorten this casually

`--md-ns 0.02` is 20 000 steps at the pipeline's 1 fs timestep, and that length is
chosen, not arbitrary. `TINKER-DETERMINISM.md` §2 measured, between two identical runs:

| step | max per-atom deviation |
|---|---|
| ≤ 500 | 0.000000 Å |
| 1010 | 0.000001 Å |
| 1510 | 0.001490 Å |
| 1860 | Å-scale |

**A run shorter than ~2000 steps reports a false pass.** Both engines use a 1 fs
timestep here (`pipeline/tinker.sh:50`, `pipeline/gromacs.sh:51`), so the same
`--md-ns` gives both the same simulated time and the same step count.

## What is held fixed

Everything except the engine:

- **Tinker** — one canonical thread-pinned `minimize` run produces `md_start.xyz`,
  shared by every replicate. `md.key` carries the pipeline's fixed
  `RANDOMSEED 28480426`.
- **GROMACS** — `md.tpr` is built once. Since `gen_vel = yes` is resolved at `grompp`
  time with `gen_seed 28480426`, the starting velocities are literally the same bytes
  in every replicate. The prep also runs `em` and `cg` once each.

So minimisation variance cannot leak into these numbers; what is reported is the MD
step alone.

## Output

`results/<UTC-timestamp>_<case>/`:

- `report.md` — verdict, pairwise RMSD / max-deviation split by protein vs water,
  most variable atoms, **divergence over the trajectory** (20 log-spaced frames), and
  the engine's final energies with their spread
- for Tinker, a separate line for the `.dyn` restart file: velocities and
  accelerations at 16 significant digits, which detects divergence still invisible in
  the `.arc` trajectory's 6 decimals
- `results.csv` — one row per replicate pair
- `*_structures.json`, `*_scalars.json` — the same data unrounded

`results/summary_<timestamp>.md` collects the verdicts across cases.

## Interpreting it

**Check the CPU-utilisation line before believing a BIT-IDENTICAL verdict.** Each
report states how many cores the engine actually averaged. An engine that never left
one thread had no OpenMP reduction order to vary, so it is deterministic for a trivial
reason and the result does not generalise. Measured on this codebase: Tinker `dynamic`
averages ~1.0 cores at both 755 and 1644 atoms, while Tinker `minimize` engages more
and *is* non-deterministic without `OPENMP-THREADS 1`. If you get a bit-identical MD
verdict at ~1.0 cores, re-run with a larger input (`-i`) before concluding anything.

Read the divergence table before the headline RMSD. Two runs that are identical for
the first 500 frames and then separate are telling you about chaotic amplification of
a low-bit force difference; two runs that differ from frame 1 are telling you about a
seed or an input that is not actually fixed.

For the GPU case, bitwise reproducibility is not available at all. The Tinker9 manual
states that evaluating forces twice on GPU is not expected to give identical answers,
and that fixed-point accumulation only "elongates the process of the inevitable
divergence". The mitigation is already active in our build
(`installer/tinker/tinker-gpu-install.sh` uses `-DPREC=mixed`, which defaults
`DETERMINISTIC_FORCE` to on). The GPU case therefore measures *how fast* it diverges,
not whether it can be stopped.

Results are only comparable **within one machine**. Every report records the host,
CPU, thread count and GPU.
