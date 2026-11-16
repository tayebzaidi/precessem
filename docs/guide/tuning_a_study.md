# Tuning a study: the mode capacitance from Python and from a study file

## What you will compute

The **mode capacitance** `C_mode` of the four-port F-small birdcage (exact
port rectangles, κ off, degree 1) around its saline phantom, with its rings
and legs cut out of the mesh as perfectly conducting holes. Every gap carries a
50 Ω lumped port; the study is solved once for its 4×4 S-matrix, and
`Study.tune(k)` reads the capacitance that, placed in each of the four gaps,
puts mode k of the four-gap structure at resonance in the lossless circuit
model at the study frequency. `tune` reuses the solved result; if the Study
has not been solved yet, it solves it once first.

The same read runs two ways on this page: from Python, with `Study.tune(1)`,
and from a study file whose top-level `tune:` section asks for k = 1, through
`fem_em_solver.run_study_file`.

## Before you start

You need the complex DolfinX build (`source /usr/local/bin/dolfinx-complex-mode`)
and **2 MPI ranks**; each block below takes ≈ 1 min on its own. The study file
is [f_small_hole_pec_tune.yaml](../../tests/package/wf1/f_small_hole_pec_tune.yaml);
the second block reads it by its repository path, so run it from the
repository root.

## Run it from Python

```python tune-python
from fem_em_solver import FourPortBirdcageGeometry, LumpedPortDrive, Study
from fem_em_solver.materials import GelledSalinePhantomMaterial

study = Study(
    FourPortBirdcageGeometry(as_hole=True, strip_fraction=0.5),
    GelledSalinePhantomMaterial(0.5, 78.0, 64.0e6),
    LumpedPortDrive(),
    frequency_hz=64.0e6,
    degree=1,
)
tuned = study.tune(1)
print(f"k = {tuned.k}, C_mode = {tuned.c_mode_farad:.4e} F")
```

`tune` is defined on a four-port Study driven by `LumpedPortDrive` only, and
k is 1 or 2. The returned `TuneResult` also carries the modal impedances
`z_k` and the off-diagonal residual of the modal diagonalisation. To tune a
matrix saved earlier without solving again, see
[Tuning a saved result](tuning_a_saved_result.md).

## Run it from a study file

```python tune-yaml
from fem_em_solver import run_study_file

run = run_study_file("tests/package/wf1/f_small_hole_pec_tune.yaml")
print(f"C_mode = {run.summary['tune']['c_mode_farad']:.4e} F")
```

The file is the same Study as the Python block, written as YAML, with
`tune: {k: 1}`. The run solves it once, calls `Study.tune(1)` on that solved
Study, and writes the result into its summary under the key `tune`.
[Tuning from a YAML file](tuning_from_yaml.md) shows the file and walks
through every key of `summary["tune"]`.

## What you should see

Both blocks print the mode-1 capacitance. The Python block's value matches
the gate's own tuned capacitance, and the study file's value matches the
study-file record, each to the first or second band below. The untuned 4×4
that both reads start from ties to the degree-1 perfect-conductor record
(relative distance of its self / adjacent / opposite class means) to the
third:

```text tune-bands
Python C_mode vs the gate's tuned capacitance: 1e-03
Study-file C_mode vs the study-file record: 1e-03
Untuned 4×4 vs the perfect-conductor record: 1e-09
```

The test, not the blocks, computes the readings and checks them against these
bands.

## Why the tune is lossless

A uniform series resistance on the four capacitors shifts only the real part
of every modal impedance, so `C_mode` does not depend on it and `tune` takes
no resistance argument. The gate's identity that places `C_mode` on the gap
sheets and recovers the circuit reduction of the untuned 4×4 is lossless (no
series resistance); the production route with 0.2 Ω on every capacitor is
[capacitor sheets with ESR](capacitor_sheets_with_esr.md).

## What this does not show

- A resonance, match, Q or `R_in` statement: `C_mode` is a lossless
  circuit-model read of one solve.
- The tuned capacitance on the copper (Leontovich) rings.
- Degree 2.
- The human-scale birdcage's ring read.
- Capacitor sheets from a study file together with `tune:` (refused).
- Any external-solver comparison.

## Worked examples

- [mri:24, a Study tuned, then capacitor sheets](../../examples/mri/24_study_tune_then_capacitor_sheets.md)

## How this page is tested

[tests/package/test_wf17_tune_tutorial.py](../../tests/package/test_wf17_tune_tutorial.py)
execs both blocks at 2 ranks in the complex build and compares their
capacitances at the bands above; it also places the Python block's `C_mode`
on three gap sheets and checks the lossless circuit-reduction identity, with
the k = 2 capacitance as the control that misses, the way the gate
[test_wf16_study_tune_capacitors.py](../../tests/package/test_wf16_study_tune_capacitors.py)
does.
