# Heating the gel sphere from YAML: temperature rise and probes

## What you will compute

A gelled-saline sphere sits in an imposed quasi-static 64 MHz field. The
field solve gives the SAR inside the sphere; scaled by `source_scale`, that
SAR is the heat source of a steady conduction solve on the sphere's own
submesh, with the sphere surface Robin-cooled to the ambient temperature.
The output is the temperature rise and the power balance in the JSON
summary. The second file adds ball-averaged temperature probes. Each file
runs through one call, `fem_em_solver.run_study_file`. For the same study
from Python, see [the same study from Python](temperature.md).

## Before you start

You need the complex DolfinX build (`source /usr/local/bin/dolfinx-complex-mode`;
the field solve is complex) and **2 MPI ranks**; each block below takes
≈ 30 s. Copy both study files next to your script:
[sphere_mat10_thermal.yaml](../../tests/package/wf1/sphere_mat10_thermal.yaml) and
[sphere_mat10_probes.yaml](../../tests/package/wf1/sphere_mat10_probes.yaml).

## The study files

| File | What it sets | Summary keys it adds | Record it carries |
| --- | --- | --- | --- |
| [sphere_mat10_thermal.yaml](../../tests/package/wf1/sphere_mat10_thermal.yaml) | `thermal: {k_w_per_m_k, h_w_per_m2_k, t_inf_c, source_scale}` on the sphere in its box | `temperature` | `WF-1` step 3: `MAT-10`'s SAR-to-heat record through the file |
| [sphere_mat10_probes.yaml](../../tests/package/wf1/sphere_mat10_probes.yaml) | the same `thermal:` section plus `outputs.probes` (ball centres and radii) | `probes` | `WF-1` step 6: the whole-ball probe identity |

Probes need a `thermal:` section, and a probe radius must be positive; each
refusal names the field it refuses.

## Temperature rise

The block runs the thermal file and prints the temperature section of the
summary.

```python sphere-heating
from fem_em_solver import run_study_file

run = run_study_file("sphere_mat10_thermal.yaml", output_dir="femem_out/sphere_mat10_thermal")
print(run.summary["temperature"])
```

## Temperature probes

The block runs the probes file and prints the ball-averaged temperature rise
of each probe.

```python sphere-probes
from fem_em_solver import run_study_file

run = run_study_file("sphere_mat10_probes.yaml", output_dir="femem_out/sphere_mat10_probes")
print(run.summary["probes"]["values_k"])
```

The same heating run from the command line:

```bash sphere-heating-cli
femem run sphere_mat10_thermal.yaml --output-dir femem_out/sphere_mat10_thermal
```

## What you should see

The temperature section of the summary has three keys, in this order:

```text sphere-heating-keys
source_power_w
convected_power_w
theta_max_k
```

The interior temperature rise is compared with the closed-form steady
solution for the uniform-sphere SAR, both for the field-driven source and
for the closed-form source solved directly on the same submesh; the source
power on the submesh is compared with the field's dissipated power; and the
whole-sphere probe is compared with the volume mean of the temperature rise.
The bands are:

```text sphere-heating-bands
Interior error vs the closed form: 1.1e-01
Interior error, direct P1 comparator: 1e-02
Source power vs the EM dissipation: 1e-10
Whole-ball probe vs the volume mean: 1e-10
```

The test, not the blocks, computes the readings from the runs and checks
them against these bands. A source twice as strong (the factor ½ of the
time-averaged SAR dropped) misses the first band, and a probe ball moved
0.1 mm misses its own reading; the test runs both controls itself.

## What this does not show

- A coil-driven source (`MAT-15`, `WF-10`).
- A transient, a duty cycle or a container wall
  ([duty_cycle_and_container_wall.md](duty_cycle_and_container_wall.md)).
- Perfusion.
- Degree 2.
- A CEM43 dose.
- Tissue.
- Any safety or compliance figure.
- Any external-solver comparison.

## Worked examples

- [mri:13, a temperature from a config file](../../examples/mri/13_temperature_from_yaml.md)

## How this page is tested

[tests/package/test_wf17_sphere_heating_tutorial.py](../../tests/package/test_wf17_sphere_heating_tutorial.py)
extracts the two blocks, execs them at 2 ranks in the complex build, and
recomputes the records above at their bands, the way the gates
[test_wf1_yaml_thermal.py](../../tests/package/test_wf1_yaml_thermal.py)
and
[test_wf1_yaml_probes.py](../../tests/package/test_wf1_yaml_probes.py)
do.
