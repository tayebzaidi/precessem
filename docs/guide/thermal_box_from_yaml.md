# Thermal-only studies from YAML: the gel box and the ball probes

## What you will compute

The F2182-style gelled-saline box (and the `MAT-19` gel ball) heated by a
prescribed source, with no field solve: a uniform source, the
transverse-quadratic source 3Q(y² + z²)/(B² + C²), or either one switched on
and off by a duty cycle. The walls cool by Robin convection into a bath; the
solve is steady or Crank–Nicolson in time; some files also read ball-average
temperature probes. Each file runs through one call,
`fem_em_solver.run_study_file`, which writes the JSON summary next to the
results.

## Before you start

You need the complex DolfinX build (`source /usr/local/bin/dolfinx-complex-mode`)
and **2 MPI ranks**; each block below takes ≈ 2 min. Copy the files a block
lists next to your script. They are the tracked study files of three worked
examples:

- `mri:29`: [thermal_box_steady_m4.yaml](../../tests/package/wf1/thermal_box_steady_m4.yaml),
  [thermal_box_wrong_wall_m4.yaml](../../tests/package/wf1/thermal_box_wrong_wall_m4.yaml),
  [thermal_box_duty_m6.yaml](../../tests/package/wf1/thermal_box_duty_m6.yaml);
- `mri:31`: [thermal_box_quadratic_steady_m4.yaml](../../tests/package/wf1/thermal_box_quadratic_steady_m4.yaml),
  [thermal_box_quadratic_cn_m4.yaml](../../tests/package/wf1/thermal_box_quadratic_cn_m4.yaml),
  [thermal_box_uniform_steady_m4_p2.yaml](../../tests/package/wf1/thermal_box_uniform_steady_m4_p2.yaml),
  [thermal_box_steady_m4_probes.yaml](../../tests/package/wf1/thermal_box_steady_m4_probes.yaml),
  [thermal_box_wrong_wall_m4_probes.yaml](../../tests/package/wf1/thermal_box_wrong_wall_m4_probes.yaml),
  [thermal_box_duty_m4_probes.yaml](../../tests/package/wf1/thermal_box_duty_m4_probes.yaml);
- `mri:32`: [thermal_box_quadratic_duty_m4.yaml](../../tests/package/wf1/thermal_box_quadratic_duty_m4.yaml),
  [thermal_ball_m19_probes.yaml](../../tests/package/wf1/thermal_ball_m19_probes.yaml),
  [thermal_ball_duty_a_probes.yaml](../../tests/package/wf1/thermal_ball_duty_a_probes.yaml).

## The study files

| File | Example | What it sets | Record it carries |
| --- | --- | --- | --- |
| [thermal_box_steady_m4.yaml](../../tests/package/wf1/thermal_box_steady_m4.yaml) | [mri:29](../../examples/mri/29_femem_run_f2182_gel_box_files.md) | uniform source, steady | `MAT-21` step 1 (P26) |
| [thermal_box_wrong_wall_m4.yaml](../../tests/package/wf1/thermal_box_wrong_wall_m4.yaml) | [mri:29](../../examples/mri/29_femem_run_f2182_gel_box_files.md) | uniform source, steady, wrong wall coefficient (a control) | `MAT-21` step 1 (P26) |
| [thermal_box_duty_m6.yaml](../../tests/package/wf1/thermal_box_duty_m6.yaml) | [mri:29](../../examples/mri/29_femem_run_f2182_gel_box_files.md) | uniform source, transient, duty cycle | `MAT-21` step 1 (P26) |
| [thermal_box_quadratic_steady_m4.yaml](../../tests/package/wf1/thermal_box_quadratic_steady_m4.yaml) | [mri:31](../../examples/mri/31_femem_run_gel_box_source_and_probe_files.md) | transverse-quadratic source, steady, P2 | P36 (S1-quad) |
| [thermal_box_quadratic_cn_m4.yaml](../../tests/package/wf1/thermal_box_quadratic_cn_m4.yaml) | [mri:31](../../examples/mri/31_femem_run_gel_box_source_and_probe_files.md) | transverse-quadratic source, Crank–Nicolson, P2 | P36 (T-quad) |
| [thermal_box_uniform_steady_m4_p2.yaml](../../tests/package/wf1/thermal_box_uniform_steady_m4_p2.yaml) | [mri:31](../../examples/mri/31_femem_run_gel_box_source_and_probe_files.md) | uniform source, steady, P2 | `MAT-21` step 3 |
| [thermal_box_steady_m4_probes.yaml](../../tests/package/wf1/thermal_box_steady_m4_probes.yaml) | [mri:31](../../examples/mri/31_femem_run_gel_box_source_and_probe_files.md) | uniform source, steady, probes | U60 probe record |
| [thermal_box_wrong_wall_m4_probes.yaml](../../tests/package/wf1/thermal_box_wrong_wall_m4_probes.yaml) | [mri:31](../../examples/mri/31_femem_run_gel_box_source_and_probe_files.md) | wrong wall, probes | U60 probe record |
| [thermal_box_duty_m4_probes.yaml](../../tests/package/wf1/thermal_box_duty_m4_probes.yaml) | [mri:31](../../examples/mri/31_femem_run_gel_box_source_and_probe_files.md) | uniform source, duty cycle, probes | U60 probe record |
| [thermal_box_quadratic_duty_m4.yaml](../../tests/package/wf1/thermal_box_quadratic_duty_m4.yaml) | [mri:32](../../examples/mri/32_femem_run_ball_probes_and_duty_quadratic_files.md) | transverse-quadratic source, duty cycle, P2 | P38 (T-qduty) |
| [thermal_ball_m19_probes.yaml](../../tests/package/wf1/thermal_ball_m19_probes.yaml) | [mri:32](../../examples/mri/32_femem_run_ball_probes_and_duty_quadratic_files.md) | gel ball, steady, probes | U60 probe record |
| [thermal_ball_duty_a_probes.yaml](../../tests/package/wf1/thermal_ball_duty_a_probes.yaml) | [mri:32](../../examples/mri/32_femem_run_ball_probes_and_duty_quadratic_files.md) | gel ball, duty cycle, probes | U60 probe record |

## Run them

One block per example. Each runs that example's files, one
`run_study_file` call per file, and prints the source class from the JSON
summary.

```python thermal-box-mri29
from fem_em_solver import run_study_file

FILES = ["thermal_box_steady_m4.yaml", "thermal_box_wrong_wall_m4.yaml", "thermal_box_duty_m6.yaml"]
runs = {name: run_study_file(name, output_dir=f"femem_out/{name[:-5]}") for name in FILES}
for name, run in runs.items():
    print(name, run.summary["temperature"]["source"])
```

```python thermal-box-mri31
from fem_em_solver import run_study_file

FILES = [
    "thermal_box_quadratic_steady_m4.yaml",
    "thermal_box_quadratic_cn_m4.yaml",
    "thermal_box_uniform_steady_m4_p2.yaml",
    "thermal_box_steady_m4_probes.yaml",
    "thermal_box_wrong_wall_m4_probes.yaml",
    "thermal_box_duty_m4_probes.yaml",
]
runs = {name: run_study_file(name, output_dir=f"femem_out/{name[:-5]}") for name in FILES}
for name, run in runs.items():
    print(name, run.summary["temperature"]["source"])
```

```python thermal-box-mri32
from fem_em_solver import run_study_file

FILES = [
    "thermal_box_quadratic_duty_m4.yaml",
    "thermal_ball_m19_probes.yaml",
    "thermal_ball_duty_a_probes.yaml",
]
runs = {name: run_study_file(name, output_dir=f"femem_out/{name[:-5]}") for name in FILES}
for name, run in runs.items():
    print(name, run.summary["family"], run.summary["temperature"].get("source"))
```

The ball files' summary carries no `source` key (a probe-file summary), so
this block prints the family and `None` for them.

`femem check` reads a file without solving it. Three malformed files show
what a refusal looks like; each exits with status 2 and names the field it
refuses (`thermal.duty`, `outputs.probes[1].radius`,
`outputs.probes[1].radius`):

```bash thermal-box-check
femem check malformed_thermal_box_duty_without_transient.yaml
femem check malformed_thermal_box_probe_radius_zero.yaml
femem check malformed_thermal_ball_probe_radius_zero.yaml
```

## What you should see

The temperature rise against the closed-form series, as an interior relative
error: the steady transverse-quadratic file, its Crank–Nicolson file at t\*,
and the duty-cycled transverse-quadratic file at t\*. The test computes these
from the runs (the series lives test-side); the blocks above only print the
source class.

```text thermal-box-records
steady transverse-quadratic interior error = 7.219878e-05
Crank-Nicolson transverse-quadratic interior error at t* = 2.694340e-04
duty-cycled transverse-quadratic interior error at t* = 2.543310e-04
```

## Start from a gradient, inside an acrylic wall

Two more files put the gel box inside a meshed acrylic wall
(`geometry.wall`: thickness, conductivity, volumetric heat capacity, and the
number of cells across it). The wall is on the top and bottom faces only
(the z faces): the box's four side faces are insulated, and the Robin cooling
acts on the outer faces of the wall alone. One file also starts the gel from a
linear profile in z (`thermal.initial.gradient_k_per_m`) instead of a uniform
temperature; the other is the same file without `initial`, starting uniform.

```python thermal-box-wall-gradient
from fem_em_solver import run_study_file

FILES = ["thermal_box_wall_gradient_m4.yaml", "thermal_box_wall_m4.yaml"]
runs = {name: run_study_file(name, output_dir=f"femem_out/{name[:-5]}") for name in FILES}
for name, run in runs.items():
    print(name, run.summary["temperature"].get("initial_gradient_k_per_m"))
```

The gradient file's summary carries `initial_gradient_k_per_m`; the uniform
file's does not, so the block prints `None` for it.

A wall needs its thickness. This file leaves it out; `femem check` exits
with status 2 and names `geometry.wall.thickness_m`:

```bash thermal-box-wall-check
femem check malformed_thermal_box_wall_without_thickness.yaml
```

The test reads the gradient run's temperature at one gel point and checks it
equals a direct solve of the same walled box, built and solved without the
YAML layer; the uniform file's reading misses it by the recorded effect of
the initial profile.

## Read the thermal dose at two points

Two more files are the walled-box files above with a `thermal.dose` section
added: [thermal_box_wall_gradient_m4_dose.yaml](../../tests/package/wf1/thermal_box_wall_gradient_m4_dose.yaml)
starts from the gradient and
[thermal_box_wall_m4_dose.yaml](../../tests/package/wf1/thermal_box_wall_m4_dose.yaml)
starts uniform. `thermal.dose.points_m` lists the gel points whose
temperature history is read at every time step, and `window_s` sets the
trailing window of the running mean (six minutes when omitted). The run's
summary then carries a `dose` entry: for each point, `cem43_s` is `MAT-20`'s
∫ R^(43 − T) dt over the history, in seconds (divide by 60 for CEM43
minutes), beside the peak of the trailing-window mean.

```python thermal-box-dose
from fem_em_solver import run_study_file

FILES = ["thermal_box_wall_gradient_m4_dose.yaml", "thermal_box_wall_m4_dose.yaml"]
runs = {name: run_study_file(name, output_dir=f"femem_out/{name[:-5]}") for name in FILES}
for name, run in runs.items():
    print(name, run.summary["dose"]["points_m"])
```

The dose is read on the walled files only, and never beside a duty cycle. A
dose section needs the wall; this file leaves out its thickness, and
`femem check` exits with status 2 and names `geometry.wall.thickness_m`:

```bash thermal-box-dose-check
femem check malformed_thermal_box_dose_wall_without_thickness.yaml
```

The test checks the gradient run's dose at both points against the dose of
a direct solve of the same walled box, and checks each history's dose and
trailing mean against a quadrature of the same samples; the uniform file's
dose misses the gradient file's at the first point. There is no dose limit,
no damage threshold and no compliance statement here: `cem43_s` is a
number read off a temperature history, nothing more.

## What this does not show

No electromagnetic solve, coil or SAR; no perfusion; P2 for the quadratic
files and the uniform P2 file, P1 for the others; m = 4 and m = 6 meshes
only, so no h-statement. The probe readings are real-build records and are
only printed beside a complex-build run here. The dose section reads CEM43 on
the walled files only: no dose limit or threshold, no dose under a duty
cycle or beside probes, no compliance statement, and no external-solver
comparison.

## Worked examples

- [mri:29, the F2182 gel box files](../../examples/mri/29_femem_run_f2182_gel_box_files.md)
- [mri:31, the gel box source and probe files](../../examples/mri/31_femem_run_gel_box_source_and_probe_files.md)
- [mri:32, the ball probes and the duty-cycled quadratic file](../../examples/mri/32_femem_run_ball_probes_and_duty_quadratic_files.md)

## How this page is tested

[tests/package/test_wf17_thermal_box_tutorial.py](../../tests/package/test_wf17_thermal_box_tutorial.py)
extracts the three blocks, checks each runs exactly its example's file set,
execs them at 2 ranks in the complex build, and recomputes the three records
above at their bands, the way the gates
[test_wf1_yaml_thermal_box_quadratic.py](../../tests/package/test_wf1_yaml_thermal_box_quadratic.py)
and
[test_wf1_yaml_thermal_box_quadratic_duty.py](../../tests/package/test_wf1_yaml_thermal_box_quadratic_duty.py)
do.
[tests/package/test_wf17_thermal_box_wall_tutorial.py](../../tests/package/test_wf17_thermal_box_wall_tutorial.py)
execs the walled-box block at 2 ranks in the real build and compares its
gradient run with the direct walled-box solve of
[test_wf1_yaml_thermal_box_wall.py](../../tests/package/test_wf1_yaml_thermal_box_wall.py).
[tests/package/test_wf17_thermal_box_dose_tutorial.py](../../tests/package/test_wf17_thermal_box_dose_tutorial.py)
execs the dose block at 2 ranks in the real build and checks its doses
against the direct walled-box history of
[test_wf1_yaml_thermal_box_dose.py](../../tests/package/test_wf1_yaml_thermal_box_dose.py).
