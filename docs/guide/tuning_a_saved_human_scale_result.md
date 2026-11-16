# Tuning a saved human-scale result: the 32-port birdcage with no second solve

## What you will compute

The **mode capacitance** `C_mode` of the tuned copper F-human birdcage
(32 ring-gap ports, Leontovich copper rings and legs, degree 1) around its
saline phantom, read from a 32×32 S-matrix that was solved once and saved to
a Touchstone file. A **fresh** `Study` with the same inputs reads that file
with `Study.load_result` and calls `Study.tune(1, n_slots=16)` on it. The
fresh Study never solves: `tune` runs on the loaded matrix.

This is [Tuning a saved result](tuning_a_saved_result.md) at human scale:
the solve needs the mesh, the tune needs only the S-matrix.

## Before you start

You need the complex DolfinX build (`source /usr/local/bin/dolfinx-complex-mode`)
and **8 MPI ranks**; the block below takes ≈ 4 minutes, nearly all of it the
one solve. This is a **heavy** run. It writes the Touchstone file
f_human_cu_64mhz.s32p into the directory you run it from.

## Solve once, save, then tune the saved file

```python loaded-tune-human
from fem_em_solver import FHumanBirdcageGeometry, LeontovichWall, LumpedPortDrive, Study
from fem_em_solver.materials import GelledSalinePhantomMaterial
from fem_em_solver.study import COPPER_SIGMA


def f_human_study():
    return Study(
        FHumanBirdcageGeometry(as_hole=True),
        GelledSalinePhantomMaterial(0.5, 78.0, 64.0e6),
        LumpedPortDrive(),
        conductor=LeontovichWall(COPPER_SIGMA),
        frequency_hz=64.0e6,
        degree=1,
    )


# Solve once and save the 32x32 S-matrix.
solved = f_human_study()
result = solved.solve()
path = result.write_touchstone("f_human_cu_64mhz.s32p")

# Later, or in another session: a fresh Study reads the file and tunes it.
fresh = f_human_study()
fresh.load_result(path)
tuned = fresh.tune(1, n_slots=16)
print(f"k = {tuned.k}, C_mode = {tuned.c_mode_farad:.4e} F")
```

On this ring `load_result` reads a file holding exactly one frequency, equal
to the Study's own (a multi-frequency file loads on the four-port only), and
one port per ring gap (32 here); any other file is refused before anything is
stored. The ring tune reads the port names from the file.
A loaded result carries an S-matrix and nothing else, so the field reads
(`sar`, `b1_plus`, `temperature`) raise `TypeError` on it.

## What you should see

The block prints the mode-1 capacitance. It matches the capacitance the same
Study gives when it is tuned straight after its own solve, to the band below:

```text loaded-tune-human-bands
Saved-file C_mode vs the solved-Study record: 1e-03
```

The test, not the block, computes the reading and checks it against this
band.

## What this does not show

- A resonance, match, Q or `R_in` statement: `C_mode` is a lossless
  circuit-model read of one solve.
- The same load from a study file through `femem run`.
- A 64-port (F-clinical) file, or any field computed from a loaded file.
- Degree 2.
- Any external-solver comparison.

## How this page is tested

[tests/package/test_wf17_loaded_tune_human_tutorial.py](../../tests/package/test_wf17_loaded_tune_human_tutorial.py)
execs the tagged block at 8 ranks in the complex build, checks that it solves
exactly once and that the fresh Study never solves, compares its capacitance
to the recorded one at the band above, and checks that the same file with
ports 1 and 2 swapped misses it.
