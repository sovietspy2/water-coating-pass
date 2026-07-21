# determinism-common

Shared machinery for `testing/minimize-deterministic-test/` and
`testing/md-deterministic-test/`. Nothing here is run directly.

## Why these suites drive the engine binaries, not the pipeline

`pipeline/tinker.sh` and `pipeline/gromacs.sh` run minimize and MD back to back in one
process. Measuring them as a unit cannot tell you which step is stochastic, and their
per-step outputs are overwritten in place. So these suites extract the exact commands,
key files and mdp files from those scripts and invoke `minimize` / `dynamic` /
`dynamic9` / `gmx mdrun` directly. Everything the engine reads is prepared once and
copied byte-identically into each replicate directory, so the only variable left is
the engine.

Anything copied out of the pipeline carries a source reference in a comment
(e.g. `pipeline/tinker.sh:129-133`). If you change the pipeline's key or mdp
settings, change them here too or the suites stop measuring the real thing.

## Files

| file | role |
|---|---|
| `determinism_lib.sh` | CLI parsing, temp work dir + cleanup, input-integrity guard, fixture build, Tinker key writers, GROMACS mdp writers, report headers |
| `compare_structures.py` | structural comparison: bit-identity, pairwise RMSD / max deviation split by protein vs water, most-variable atoms, divergence-vs-frame curve |
| `parse_scalars.py` | engine scalars: Tinker final energy and iteration count, GROMACS potential energy and step count |

## The fixture

Every suite builds the same starting structure, once per run:

```
cp <input.pdb> $WORK/fixture/fixture.pdb
pipeline/pdb-preprocessor.py --target fixture.pdb    # single model, repair atoms, strip waters
wdrop --file fixture.pdb --sigma 1.8 --weed-dist 3.5 --layers 1
pipeline/format_pdb.py fixture_1WAT.pdb              # canonical columns, OW->O, SOL->WAT
```

`wdrop` is deterministic (`testing/algo-deterministic-test/`), so building the fixture
once and sharing it is what keeps the measured variance attributable to the engine.

The input PDB is never written to. Its sha256 is recorded before the run and checked
again on exit; a mismatch fails the run.

## Precision, and why it decides the design

`TINKER-DETERMINISM.md` §2 measured divergence of **1e-6 Å at MD step 1010**. A format
that cannot represent that number will report a false pass:

| format | resolution | used for |
|---|---|---|
| Tinker `.xyz` / `.arc` | `%12.6f` Å = 1e-6 Å | Tinker structures and trajectories |
| Tinker `.dyn` | 16 significant digits | the definitive Tinker MD bit-check |
| GROMOS96 `.g96` | `%15.9f` nm = 1e-8 Å | all GROMACS comparisons |
| GROMACS `.gro` | `%8.3f` nm = 1e-2 Å | **not used** — 4 orders of magnitude too coarse |

GROMACS output is therefore always passed through
`gmx trjconv -o traj.g96` before comparison, never read from `after_*.gro`.

One consequence: `trjconv` writes g96 `POSITIONRED` blocks — bare coordinates with no
atom or residue names — so the protein/water split cannot come from the compared file.
`gmx_write_atom_labels` derives it once from `pdb2gmx`'s `conf.gro` (same atom count
and order as the `.tpr` and every trajectory from it) and passes it as
`compare_structures.py --labels`. Tinker needs none of this: its `.xyz` carries the
`OW`/`HW` atom names directly.

## Arms

Each case can run the command under two configurations, so the report explains the
variance rather than just reporting it:

| case | `pipeline` | `hardened` |
|---|---|---|
| Tinker minimize | `minimize.key` **with** `OPENMP-THREADS 1` (as shipped) | **without** it — the control that is expected to diverge |
| Tinker MD (CPU) | `md.key` as shipped, no thread pin | `md.key` + `OPENMP-THREADS 1` |
| Tinker MD (GPU) | `md.key` + `INTEGRATOR VERLET` + `REMOVE-INERTIA 0` | *(none — upstream offers no lever)* |
| GROMACS, all stages | `gmx mdrun` as the pipeline calls it | `gmx mdrun -reprod -ntmpi 1 -ntomp 1` |

Note the direction differs by case. For Tinker minimize the pipeline already ships the
deterministic setting, so `hardened` is the negative control. For Tinker MD it does
not, so `hardened` is the candidate fix.

`-reprod` and `-ntmpi/-ntomp` change the *numerical result*, not only its
reproducibility — a `hardened` arm answers "is this reproducible", not "does this
agree with the pipeline".

## Using the comparison scripts on their own

```bash
testing/determinism-common/compare_structures.py \
  --format tinker-xyz --label "my run" a/sys.xyz_2 b/sys.xyz_2 c/sys.xyz_2

testing/determinism-common/compare_structures.py \
  --format g96 --trajectory --frame-samples 20 --label "my md" \
  --labels prep/atom_labels.txt */traj.g96

testing/determinism-common/parse_scalars.py --kind tinker-minimize --label "my run" */minimize.log
```

Formats: `tinker-xyz` (also reads `.arc`), `g96`, `gro`.
Scalar kinds: `tinker-minimize`, `tinker-dynamic`, `gromacs-minimize`, `gromacs-md`.
