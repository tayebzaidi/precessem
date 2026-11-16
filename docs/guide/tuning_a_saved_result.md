# Tuning a saved result: the mode capacitance with no second solve

## What you will compute

The **mode capacitance** `C_mode` of the four-port F-small birdcage (exact
port rectangles, κ off, degree 1, perfectly conducting rings and legs cut out
of the mesh) around its saline phantom, read from a 4×4 S-matrix that was
solved once and saved to a Touchstone file. A **fresh** `Study` with the same
inputs reads that file with `Study.load_result` and calls `Study.tune(1)` on
it. The fresh Study never solves: `tune` runs on the loaded matrix.

This is the same read as [Tuning a study](tuning_a_study.md), split in two:
the solve, which needs the mesh, and the tune, which needs only the S-matrix.
The same tune from a study file through `femem run` is
[Tuning a saved result from a YAML file](tuning_a_saved_result_from_yaml.md).

## Before you start

You need the complex DolfinX build (`source /usr/local/bin/dolfinx-complex-mode`)
and **2 MPI ranks**; the block below takes ≈ 1 min, nearly all of it the one
solve. It writes the Touchstone file into the directory you run it from.

## Solve once, save, then tune the saved file

```python loaded-tune
from fem_em_solver import FourPortBirdcageGeometry, LumpedPortDrive, Study
from fem_em_solver.materials import GelledSalinePhantomMaterial


def f_small_study():
    return Study(
        FourPortBirdcageGeometry(as_hole=True, strip_fraction=0.5),
        GelledSalinePhantomMaterial(0.5, 78.0, 64.0e6),
        LumpedPortDrive(),
        frequency_hz=64.0e6,
        degree=1,
    )


# Solve once and save the 4x4 S-matrix.
solved = f_small_study()
result = solved.solve()
path = result.write_touchstone("f_small_pec_64mhz.s4p")

# Later, or in another session: a fresh Study reads the file and tunes it.
fresh = f_small_study()
fresh.load_result(path)
tuned = fresh.tune(1)
print(f"k = {tuned.k}, C_mode = {tuned.c_mode_farad:.4e} F")
```

`load_result` reads a file holding one frequency equal to the Study's own, or
several whose span contains it, and as many ports as the birdcage has legs;
any other file is refused before anything is stored. On a multi-frequency
file, `tune` needs `at_hz`: it picks the file's nearest frequency, with no
interpolation, and refuses a value outside the span. It is defined on a four-port Study driven
by `LumpedPortDrive` only. A loaded result carries an S-matrix and nothing
else, so the field reads (`sar`, `b1_plus`, `b1_plus_at`) raise `TypeError` on
it, naming `load_result`.

## What you should see

The block prints the mode-1 capacitance. It matches the capacitance the same
Study gives when it is tuned straight after its own solve, to the band below:

```text loaded-tune-bands
Saved-file C_mode vs the solved-Study record: 1e-03
```

The test, not the block, computes the reading and checks it against this
band.

## A file holding several frequencies

A Touchstone file can hold a sweep. [A frequency sweep from YAML](sweep_from_yaml.md)
writes one from [f_small_sweep.yaml](../../tests/package/wf1/f_small_sweep.yaml):
produce your own that way and name it `f_small_sweep.s4p`. This repository does
not carry the file. From there on nothing is meshed or solved, so each route
below takes about a second. `tune` needs `at_hz` on such a file, and the
capacitance is read with the selected frequency's own ω.

```python swept-tune
from fem_em_solver import FourPortBirdcageGeometry, LumpedPortDrive, Study
from fem_em_solver.materials import GelledSalinePhantomMaterial

study = Study(
    FourPortBirdcageGeometry(as_hole=True, strip_fraction=0.5),
    GelledSalinePhantomMaterial(0.5, 78.0, 64.0e6),
    LumpedPortDrive(),
    frequency_hz=64.0e6,
    degree=1,
)
study.load_result("f_small_sweep.s4p")
rows = [study.tune(1, at_hz=f) for f in (1.0e7, 6.4e7, 1.28e8)]
for tuned in rows:
    c = float(tuned.c_mode_farad)
    print(f"{tuned.frequency_hz:.6e} Hz: C_mode = {c:.6e} F ({c.hex()})")
```

The same read from a study file, next to the Touchstone file, is one
`at_hz` per file. Run it with `mpiexec -n 2 femem run <file> --output-dir <dir>`;
the `tune` block of the summary in that directory holds the selected
`frequency_hz` and `c_mode_farad`.

```yaml swept-tune-file
femem: 1
study:
  frequency_hz: 64000000.0
  degree: 1
geometry:
  kind: four_port_birdcage
  preset: f_small
  as_hole: true
material:
  sigma: 0.5
  epsilon_r: 78.0
excitation:
  kind: lumped_port_drive
load_result: f_small_sweep.s4p
tune:
  k: 1
  at_hz: 64000000.0
```

On U71's sweep file, both routes give these three lines, the block by printing
them and the study file through its summary at each `at_hz`. The hexadecimal
figure is the exact value:

```text swept-tune-lines
1.000000e+07 Hz: C_mode = 2.472676e-09 F (0x1.53d78952a75e9p-29)
6.400000e+07 Hz: C_mode = 6.166035e-11 F (0x1.0f2f6316ff1dep-34)
1.280000e+08 Hz: C_mode = 1.513381e-11 F (0x1.0a3c9bbb69f6ep-36)
```

Each value matches the single-frequency tune of the same column to the band
below:

```text swept-tune-bands
Swept-file C_mode vs the single-frequency tune of its column: 1e-03
```

`tune` selects the file's nearest frequency and never interpolates between two
of them. An `at_hz` outside the file's span is refused: `tune` raises
`ValueError`, and `femem run` refuses the field `tune.at_hz`, exits with status
2 and writes no output directory. `tune` with no `at_hz` raises `ValueError`
naming `at_hz`. `femem run` with no `tune.at_hz` on a swept file refuses that
field with status 2.

## The human-scale birdcage

On the tuned copper F-human birdcage (32 ring-gap ports, Leontovich copper)
the tune reads the ring mode with `n_slots=16`. `load_result` accepts this
Study since `WF-16` step 15 (a 32-port Touchstone file, one per ring-gap
port, tuned with `n_slots=16`); the block below still solves first and tunes
the solved Study:

```python
from fem_em_solver import FHumanBirdcageGeometry, LeontovichWall, LumpedPortDrive, Study
from fem_em_solver.materials import GelledSalinePhantomMaterial
from fem_em_solver.study import COPPER_SIGMA

study = Study(
    FHumanBirdcageGeometry(as_hole=True),
    GelledSalinePhantomMaterial(0.5, 78.0, 64.0e6),
    LumpedPortDrive(),
    conductor=LeontovichWall(COPPER_SIGMA),
    frequency_hz=64.0e6,
    degree=1,
)
tuned = study.tune(1, n_slots=16)
```

This block is **not run** by this page's test: it is a heavy solve. Its
recorded reading is 1.181185799070e-10 F, at 8 ranks, logged at
`docs/testing/logs/20261015T213534Z_WF-16-s12-W1.log.gz` line 7375 and quoted
at `tests/package/test_wf1_yaml_tune_f_human.py` lines 75–76.

## What this does not show

- A resonance, match, Q or `R_in` statement: `C_mode` is a lossless
  circuit-model read of one solve.
- Loading a saved human-scale (32-port) result: [Tuning a saved human-scale result](tuning_a_saved_human_scale_result.md).
- A value between a file's frequencies: `tune` selects the nearest one.
- A saved human-scale file holding several frequencies: it is refused.
- Degree 2, or the copper (Leontovich) rings on F-small.
- Any external-solver comparison.

## How this page is tested

[tests/package/test_wf17_loaded_tune_tutorial.py](../../tests/package/test_wf17_loaded_tune_tutorial.py)
execs the tagged block at 2 ranks in the complex build, checks that the fresh
Study never solves, compares its capacitance to the recorded one at the band
above, and checks that the same block with `tune(2)` misses it. It also checks
that `load_result` on the human-scale Study refuses a missing file and leaves
no result, and that the human-scale block is never executed. For the several-frequency section it execs the tagged
block and runs the tagged study file at each `at_hz` on U71's file (absent, the
test fails; it never re-solves), with every mesh build and solve made to raise.
It checks both routes against the three lines above and the band, checks that
each `at_hz` moved one frequency along misses, and checks the refusals named
above on the page's own block and study file.
