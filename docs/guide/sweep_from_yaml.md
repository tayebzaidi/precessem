# A frequency sweep from YAML: one mesh, one Touchstone file

## What you will compute

The F-small four-port birdcage (exact port rectangles, κ off, degree 1)
around its saline phantom, meshed once and solved at every frequency listed
under `study.sweep_hz`, one port at a time at 1 V. The output is one
multi-frequency Touchstone file of the 4×4 S-matrix plus the JSON summary,
which holds the S- and Z-matrices per frequency. The whole sweep runs through
one call, `fem_em_solver.run_study_file`.

## Before you start

You need the complex DolfinX build (`source /usr/local/bin/dolfinx-complex-mode`)
and **2 MPI ranks**; the block below takes ≈ 1 min. Copy the study file
next to your script:
[f_small_sweep.yaml](../../tests/package/wf1/f_small_sweep.yaml).

## The study file

| File | What it sets | Summary keys it adds | Record it carries |
| --- | --- | --- | --- |
| [f_small_sweep.yaml](../../tests/package/wf1/f_small_sweep.yaml) | `study.sweep_hz` (an increasing list of frequencies) on the four-port birdcage | `sweep_hz`, `s_matrices`, `z_matrices`, `touchstone` | `WF-1` step 21 (U71): the S-matrix at each frequency reproduces the `ANS-4` degree-1 records on one mesh, and the Touchstone file reads back the solved S |

A sweep file refuses a conductor, `as_hole`, capacitor sheets, SAR,
`|B₁⁺|` and probe outputs, a `tune:` section, and the Palace backend; each
refusal names the field it refuses.

## Run it

The block runs the file and prints the swept frequencies and the name of the
Touchstone file it wrote.

```python four-port-sweep
from fem_em_solver import run_study_file

run = run_study_file("f_small_sweep.yaml", output_dir="femem_out/f_small_sweep")
print(run.summary["sweep_hz"], run.paths["touchstone"].name)
```

The same run from the command line:

```bash four-port-sweep-cli
femem run f_small_sweep.yaml --output-dir femem_out/f_small_sweep
```

## What you should see

The S-matrix at each frequency against the degree-1 record at that frequency
(relative Frobenius distance), and the Touchstone file read back against the
solved S-matrices and frequencies. The bands are:

```text four-port-sweep-bands
S vs the P1 records, per frequency: 1e-09
Touchstone round trip, S: 1e-12
Touchstone round trip, frequency (Hz): 1e-03
```

The test, not the block, computes the readings from the run and checks them
against these bands.

## What this does not show

- Degree 2 (its record is a wider-rank record and is not run here).
- Copper, hole, capacitor or F-human sweeps (refused).
- SAR, `|B₁⁺|` or probes per frequency (refused).
- A resonance or tuning statement.
- The Palace backend (refused).
- Any external-solver comparison.

## Worked examples

- [mri:34, the sweep file through `femem run`](../../examples/mri/34_femem_run_sweep_file_touchstone.md)
- [mri:09, a study from YAML](../../examples/mri/09_study_from_yaml.md)

## How this page is tested

[tests/package/test_wf17_sweep_tutorial.py](../../tests/package/test_wf17_sweep_tutorial.py)
extracts the block, execs it at 2 ranks in the complex build, and recomputes
the per-frequency ties and the round trip from the run at the bands above,
the way the gate
[test_wf1_yaml_sweep.py](../../tests/package/test_wf1_yaml_sweep.py)
does.
