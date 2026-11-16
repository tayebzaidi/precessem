# Tuning from a YAML file

## What you will compute

The **mode capacitance** of the four-port F-small birdcage, from a study file:
the capacitance `C_mode` that, placed in each of the four gaps, puts mode k of
the four-gap structure at resonance in the lossless circuit model at the
study frequency. A top-level `tune:` section in the file asks for it.
`femem run` (or `run_study_file` from Python) solves the Study once, calls
`Study.tune(k)` on that solved Study — no second solve — and writes the
result into the run's JSON summary under the key `tune`, beside the
Touchstone file.

## Before you start

- The complex DolfinX build: `source /usr/local/bin/dolfinx-complex-mode`.
- 2 MPI ranks.
- The file sets `degree: 1`. An omitted `degree` means the `Study` default,
  degree 2 (the production order); this page reads degree 1.

## The study file

This is
[`f_small_hole_pec_tune.yaml`](../../tests/package/wf1/f_small_hole_pec_tune.yaml),
verbatim:

```yaml f-small-hole-pec-tune
# WF-1 step 11 (§9 item 442): f_small_hole_pec_control.yaml plus a top-level tune section (Study.tune, k = 1).
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
tune:
  k: 1
```

The coil is the F-small four-leg birdcage with PEC legs cut out of the air
(`as_hole: true`) and four exact-rectangle lumped ports, one per gap. The
`tune:` section is defined on a `four_port_birdcage` geometry only; any other
`geometry.kind` is refused when the file is parsed. `k` is 1 or 2; `tune: {}`
(no `k`) means `Study.tune`'s own default, k = 1.

## Run it

Save the file above as `f_small_hole_pec_tune.yaml` and, from the same
directory, run this block as a Python script on 2 ranks:

```python tune-yaml
from fem_em_solver import run_study_file

# save the study file above as f_small_hole_pec_tune.yaml next to this script
run = run_study_file("f_small_hole_pec_tune.yaml", output_dir="femem_out")
print(f"C_mode = {run.summary['tune']['c_mode_farad']:.4e} F")
```

or make the same call through the console script:

```bash tune-yaml-cli
source /usr/local/bin/dolfinx-complex-mode
PYTHONPATH=/workspace/src mpiexec -n 2 femem run f_small_hole_pec_tune.yaml --output-dir femem_out
```

## What you should see

```text tune-yaml-output
C_mode = 6.6268e-11 F
```

## What `summary["tune"]` carries

```text tune-summary-keys
k                      the mode index that was tuned (1 or 2)
frequency_hz           the Study frequency the capacitance is read at, in Hz
c_mode_farad           C_mode = 1 / (omega Im z_k), omega = 2 pi frequency_hz, in F
z_k                    the four modal impedances, each as [re, im], in ohm
off_diagonal_residual  how far the DFT of the port impedance matrix is from diagonal
```

`z_k` is the diagonal of `F Z F^H`, with `Z` the 4×4 port impedance matrix of
the solved Study and `F` the unitary four-point DFT. All five values are plain
floats, so the JSON summary carries them without loss.

## What this does and does not claim

- **A mode capacitance, not a resonance.** `C_mode` is the circuit-model
  capacitance that would resonate mode k at the study frequency (`TH-17`'s
  construction). It is not a resonance found in the field, not an `f₀`, not a
  Q, not a match and not an absolute `S₁₁`. The tuned resonance of this coil
  is a separate study (`ANS-7`).
- **Lossless.** A uniform series resistance on the four capacitors moves
  `Re z_k` only, so `C_mode` does not depend on it and `tune:` takes no
  resistance.
- **Degree 1 on this page.** The production order is degree 2; the degree-2
  `Study.tune` readings are records of their own tests (listed under
  [`Study`](validation.md#study) on the validation page), not of this file.
- **One mesh.** No h-convergence statement.
- **Two routes through YAML.** On a `four_port_birdcage` the `tune:` section
  reads the four-port route of `Study.tune`. On an `f_human_birdcage` with
  `excitation.kind: lumped_port_drive` it reads the 32-port ring route
  (degree 1, the copper hole fixture; its nodes are registered under
  `run_study_file` on the validation page).

## Worked examples

- [`mri:19`](../../examples/mri/19_study_tune_mode_capacitance.md) —
  `Study.tune` from Python on the same coil.
- [`mri:21`](../../examples/mri/21_femem_run_tune_section.md) — this same
  study file through `femem run`, with the summary's `tune` block read back.

To tune a saved S-matrix with no solve, add `load_result:` to this file:
see [Tuning a saved result from a YAML file](tuning_a_saved_result_from_yaml.md).

## How this page is tested

This page is not executed. Its test,
[`test_wf17_tuning_yaml_page.py`](../../tests/package/test_wf17_tuning_yaml_page.py),
checks that the study file above equals the tracked file byte for byte, that
the key list above is exactly the set of keys `run_study_file` writes, that
the printed `C_mode` is the `TH-17` record computed from imported constants,
and that every test node named below is registered in the validation
envelope. The solve itself is tested by these nodes, all listed on the
[validation page](validation.md) under [`Study`](validation.md#study) and
[`run_study_file`](validation.md#run_study_file):

- `tests/package/test_wf16_study_tune.py::test_wf16_tune_a_c_mode_equals_th17_record`
- `tests/package/test_wf1_yaml_tune.py::test_wf1_s11_yaml_equals_same_study_api`
- `tests/package/test_wf1_yaml_tune.py::test_wf1_s11_c_mode_equals_th17_record`

API reference: [`Study`](api/Study.md), [`run_study_file`](api/run_study_file.md).
