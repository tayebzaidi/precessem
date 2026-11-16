# Four-port EM outputs from YAML: 10 g SAR, |B₁⁺| and point probes

## What you will compute

The F-small four-port birdcage (exact port rectangles, κ off, degree 1,
64 MHz) around its saline phantom, driven one port at a time at 1 V. From the
one field-keeping solve the file's `outputs` share, a study file can ask for
three readings per port drive: the 10 g mass-averaged SAR in a ball at one
centre per port (the `{ring_radius_m: r}` form puts each centre at radius r on
its port's azimuth), `|B₁⁺|` at listed points, and point SAR and `|E|` at
listed points inside the phantom. Each file runs through one call,
`fem_em_solver.run_study_file`, which writes the JSON summary (and the
Touchstone file of the S-matrix) next to the results.

## Before you start

You need the complex DolfinX build (`source /usr/local/bin/dolfinx-complex-mode`)
and **2 MPI ranks**; each block below takes ≈ 2 min. Copy the two study files
next to your script:
[f_small_sar_b1.yaml](../../tests/package/wf1/f_small_sar_b1.yaml) and
[f_small_em_probes.yaml](../../tests/package/wf1/f_small_em_probes.yaml).

## The study files

| File | What it sets | Summary keys it adds | Record it carries |
| --- | --- | --- | --- |
| [f_small_sar_b1.yaml](../../tests/package/wf1/f_small_sar_b1.yaml) | `outputs.sar` (10 g, centres on a ring) and `outputs.b1_plus` (listed points) | `sar`, `b1_plus` | `WF-1` step 5 (item 344): the 10 g C4 cyclic pairs and the `|B₁⁺|` C4 covariance |
| [f_small_em_probes.yaml](../../tests/package/wf1/f_small_em_probes.yaml) | `outputs.probes` (five points in the phantom) | `probes` | `WF-1` step 20 (U70): the summary is `point_sar` on the run's own fields |

Probe points outside the phantom are refused, because the σ and ρ that turn
`|E|` into SAR are the phantom's. On this family `probes` takes
`{points: […]}`; on the thermal files it takes a list of balls instead.

## Run them

Each block runs one file and prints which port drives its summary holds.

```python four-port-sar-b1
from fem_em_solver import run_study_file

run = run_study_file("f_small_sar_b1.yaml", output_dir="femem_out/f_small_sar_b1")
print(sorted(run.summary["sar"]["values"]))
print(sorted(run.summary["b1_plus"]["values"]))
```

```python four-port-probes
from fem_em_solver import run_study_file

run = run_study_file("f_small_em_probes.yaml", output_dir="femem_out/f_small_em_probes")
print(sorted(run.summary["probes"]["values"]))
```

`femem check` reads a file without solving it. Two malformed files show what
a refusal looks like; each exits with status 2 and names the field it refuses
([malformed_sar_mass_zero.yaml](../../tests/package/wf1/malformed_sar_mass_zero.yaml)
names `outputs.sar.mass_g`;
[malformed_sar_on_sphere.yaml](../../tests/package/wf1/malformed_sar_on_sphere.yaml)
names `outputs.sar`, because SAR outputs need a port study):

```bash four-port-check
femem check malformed_sar_mass_zero.yaml
femem check malformed_sar_on_sphere.yaml
```

## What you should see

The C4 symmetry of the birdcage, read from the run. The 10 g SAR at port k's
own centre under port k's drive, compared cyclically between neighbours; and
`|B₁⁺|` under ports 2, 4 and 3 on a point set rotated by +Δ, −Δ and 2Δ (Δ the
port spacing), against port 1 on the unrotated set. The test, not the block,
computes these from the run; the pairing reads the solved sheets' azimuths,
not the port names.

```text four-port-records
10 g C4 cyclic pairs (%): 0.8702 3.5284 2.0337 4.5024
|B1+| C4 covariance (%): 1.8695 1.9494 1.6405
```

## What this does not show

These are per-1 V port-drive records, not absolute or compliance SAR or
`|B₁⁺|`; there is no spatial-peak search and no homogeneity figure. Degree 1
only (the degree-2 record is a wider-rank record and is not run here); the
coil is the solid σ = 800 S/m conductor, not copper; no F-human; no
external-solver comparison.

## Worked examples

- [mri:14, the same outputs through the Python `Study`](../../examples/mri/14_study_sar_and_b1.md)
- [mri:09, a study from YAML](../../examples/mri/09_study_from_yaml.md)

## How this page is tested

[tests/package/test_wf17_four_port_outputs_tutorial.py](../../tests/package/test_wf17_four_port_outputs_tutorial.py)
extracts the two blocks, execs them at 2 ranks in the complex build, and
recomputes the records above from each run at their bands and printed digits,
the way the gates
[test_wf1_yaml_sar_b1.py](../../tests/package/test_wf1_yaml_sar_b1.py)
and
[test_wf1_yaml_em_probes.py](../../tests/package/test_wf1_yaml_em_probes.py)
do.
