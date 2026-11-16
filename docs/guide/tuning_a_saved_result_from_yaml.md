# Tuning a saved result from a YAML file

## What you will compute

The **mode capacitance** `C_mode` of the four-port F-small birdcage (exact
port rectangles, κ off, degree 1, perfectly conducting rings and legs cut out
of the mesh) around its saline phantom, read by `femem run` from a 4×4
S-matrix that was solved once and saved to a Touchstone file. The study file
names that file under a top-level `load_result:` key beside its `tune:`
section, so `femem run` calls `Study.load_result` and `Study.tune(1)` and
never meshes or solves.

This is [Tuning a saved result](tuning_a_saved_result.md) with the second
half moved into a study file, and [Tuning from a YAML file](tuning_from_yaml.md)
with the solve replaced by a saved file.

## Before you start

You need the complex DolfinX build (`source /usr/local/bin/dolfinx-complex-mode`)
and **2 MPI ranks**. The first step takes ≈ 1 min, nearly all of it the one
solve; the second takes about a second. Run every step from the same
directory.

## Step 1: solve once and save the S-matrix

Run this block as a Python script on 2 ranks. It writes
`f_small_pec_64mhz.s4p` into the directory you run it from.

```python loaded-yaml-solve
from fem_em_solver import FourPortBirdcageGeometry, LumpedPortDrive, Study
from fem_em_solver.materials import GelledSalinePhantomMaterial

solved = Study(
    FourPortBirdcageGeometry(as_hole=True, strip_fraction=0.5),
    GelledSalinePhantomMaterial(0.5, 78.0, 64.0e6),
    LumpedPortDrive(),
    frequency_hz=64.0e6,
    degree=1,
)
result = solved.solve()
path = result.write_touchstone("f_small_pec_64mhz.s4p")
print(f"wrote {path}")
```

## Step 2: the study file

Save this as `f_small_pec_load.yaml` in the same directory. It is
[`f_small_hole_pec_tune.yaml`](../../tests/package/wf1/f_small_hole_pec_tune.yaml)
(the study of step 1, plus `tune: {k: 1}`) with one line added, `load_result:`.

```yaml loaded-yaml-file
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
load_result: f_small_pec_64mhz.s4p
tune:
  k: 1
```

The study sections must still be there: `load_result` reads a file with four
ports, one per leg, holding one frequency equal to `study.frequency_hz` or
several frequencies. A file holding several needs `tune.at_hz`; without it
`run_study_file` stops with a `ValueError` naming `at_hz`, and `femem run`
refuses the field `tune.at_hz` with status 2. A relative `load_result:` path is
resolved against the directory of the YAML file, not the directory you run
from.

`load_result:` is defined on a `four_port_birdcage`, an `f_clinical_birdcage`
or an `f_human_birdcage`, with `excitation.kind: lumped_port_drive` only, and
it needs a `tune:` section. A loaded result has an S-matrix and no fields, so
`conductor`, `thermal`, `study.sweep_hz` and the field outputs (`sar`,
`b1_plus`, `probes`) are refused when the file is parsed, each by name.
`conductor` is admitted on an `f_human_birdcage` alone, where it records what
the saved file was solved with and nothing is solved.

## Step 3: check and run it

```bash loaded-yaml-cli
source /usr/local/bin/dolfinx-complex-mode
mpiexec -n 2 femem check f_small_pec_load.yaml
mpiexec -n 2 femem run f_small_pec_load.yaml --output-dir femem_out
```

`femem check` prints the resolved path of the saved file, so you can see
which file the run will read before it runs. With the directory you work in
in place of `/your/work/dir`, its `load_result` line is:

```text loaded-yaml-check
load_result: '/your/work/dir/f_small_pec_64mhz.s4p'
```

## What you should see

The run's JSON summary, written into `femem_out`, carries the same `tune` block the solve route
writes (its keys are listed on [Tuning from a YAML file](tuning_from_yaml.md)),
plus `load_result`, the resolved path that was read. It carries no `cells` and
no `touchstone`: nothing was meshed or solved. Its `c_mode_farad` matches the
capacitance the solved Study gives when it is tuned straight after its own
solve, to the band below:

```text loaded-yaml-bands
Saved-file C_mode from femem run vs the solved-Study record: 1e-09
```

The test, not the page, computes the reading and checks it against this
band.

## What this does not show

- A resonance, match, Q or `R_in` statement: `C_mode` is a lossless
  circuit-model read of one solve.
- Loading a saved human-scale (32-port) result: [Tuning a saved human-scale result](tuning_a_saved_human_scale_result.md).
- Any field output (SAR, `B1+`, probes) from a loaded file; these are refused
  by design.
- A worked file holding several frequencies: [A file holding several frequencies](tuning_a_saved_result.md#a-file-holding-several-frequencies).
- Degree 2, or the copper (Leontovich) rings on F-small.
- Any external-solver comparison.

## How this page is tested

[tests/package/test_wf17_loaded_tune_yaml_tutorial.py](../../tests/package/test_wf17_loaded_tune_yaml_tutorial.py)
execs step 1 at 2 ranks in the complex build, writes step 2's file verbatim
next to it, and runs step 3's two `femem` commands through the console entry
point with every mesh build and solve made to raise. It compares the run's
capacitance to the recorded one at the band above, checks that the run's
summary equals the one `run_study_file` writes for the same file, and checks
that the same file with `k: 2` misses the record. It also checks that step 2's
file is the tracked study file plus `load_result:`, and that `femem check`
prints the line shown above. It also checks the statements above on a
two-frequency file with and without `tune.at_hz`, and on an `f_human_birdcage`
file with `conductor` beside `load_result:`.
