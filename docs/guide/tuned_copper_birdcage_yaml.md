# The tuned copper birdcage from a YAML file

## What you will compute

The same coil as [the tuned copper birdcage](tuned_copper_birdcage.md) — the
human-scale 16-leg copper birdcage, tuned at 64 MHz, a 0.2 Ω series resistance
on every one of its 32 capacitor gaps, driven from two feeds 90° apart — but
described in a study file instead of Python. `femem run` (or
`run_study_file` from Python) solves the file once and writes a JSON summary
holding the node S-matrix, the active admittances `Y_in` for each sense of
rotation, the 10 g mass-averaged SAR at 19 ball centres, and |B₁⁺| at two
points.

The first 18 centres are a fixed lattice in the phantom (three radii, four
azimuths at the mid-plane, and two off-plane rows). The 19th centre is the
image of the lattice's hottest point c\* under the plane bisecting the two
feeds: that mirror maps the feed azimuth 11.25° onto 101.25°, so it lands
exactly on lattice row 6 (r = 45 mm, a quarter turn from c\*). It is **not an
independent centre**; it is listed so that the summary carries the
counter-rotating (`cw`) reading at the mirror image of c\* beside c\*'s
co-rotating (`ccw`) one.

## Before you start

- The complex DolfinX build: `source /usr/local/bin/dolfinx-complex-mode`.
- 8 MPI ranks, ≈ 4–5 min.
- The file sets `degree: 1`. An omitted `degree` means the `Study` default,
  degree 2, which this coil cannot afford at this width.

## The study file

This is
[`f_human_tuned_copper.yaml`](../../tests/package/wf1/f_human_tuned_copper.yaml),
verbatim:

```yaml f-human-tuned-copper
# WF-1 step 7a (§9 item 390a): item 382 / 382b's tuned copper F-human Study (hole route), degree 1, 64 MHz.
femem: 1
study:
  frequency_hz: 64000000.0
  degree: 1
geometry:
  kind: f_human_birdcage
  as_hole: true
material:
  sigma: 0.5
  epsilon_r: 78.0
excitation:
  kind: two_port_quadrature
  c_in_f: 1.1824370000000001e-10
  esr_ohm: 0.2
conductor:
  kind: leontovich
  sigma_s_per_m: 58000000.0
outputs:
  sar:
    mass_g: 10.0
    centres:
      - [0.014711779206048456, 0.0029263548302419238, 0.0]
      - [0.02942355841209691, 0.0058527096604838475, 0.0]
      - [0.044135337618145365, 0.00877906449072577, 0.0]
      - [-0.002926354830241923, 0.014711779206048456, 0.0]
      - [-0.005852709660483846, 0.02942355841209691, 0.0]
      - [-0.008779064490725769, 0.044135337618145365, 0.0]
      - [-0.014711779206048456, -0.0029263548302419255, 0.0]
      - [-0.02942355841209691, -0.005852709660483851, 0.0]
      - [-0.044135337618145365, -0.008779064490725776, 0.0]
      - [0.0029263548302419246, -0.014711779206048456, 0.0]
      - [0.005852709660483849, -0.02942355841209691, 0.0]
      - [0.008779064490725774, -0.044135337618145365, 0.0]
      - [0.014711779206048456, 0.0029263548302419238, -0.05]
      - [0.02942355841209691, 0.0058527096604838475, -0.05]
      - [0.044135337618145365, 0.00877906449072577, -0.05]
      - [0.014711779206048456, 0.0029263548302419238, 0.05]
      - [0.02942355841209691, 0.0058527096604838475, 0.05]
      - [0.044135337618145365, 0.00877906449072577, 0.05]
      - [-0.008779064490725769, 0.044135337618145365, 0.0]
  b1_plus:
    points:
      - [0.0, 0.0, 0.0]
      - [0.044135337618145365, 0.00877906449072577, 0.0]
```

The copper is not meshed: `as_hole: true` cuts the coil out of the air and
`conductor: leontovich` makes its surface a surface-impedance wall.
`c_in_f` is the in-model Larmor capacitance on copper, an input; this page
does not tune the coil.

## Run it

Save the file above as `f_human_tuned_copper.yaml` and, from the same
directory, run this block as a Python script on 8 ranks:

```python tuned-copper-yaml
from fem_em_solver import run_study_file

# save the study file above as f_human_tuned_copper.yaml next to this script
run = run_study_file("f_human_tuned_copper.yaml", output_dir="femem_out")
ccw = run.summary["sar"]["values"]["ccw"][:18]  # the 18 lattice centres
sar_star = max(entry["averaged_sar_w_per_kg"] for entry in ccw)
print(f"SAR10g(c*) = {sar_star:.4e} W/kg")
```

or make the same call through the console script:

```bash tuned-copper-yaml-cli
source /usr/local/bin/dolfinx-complex-mode
PYTHONPATH=/workspace/src mpiexec -n 8 femem run f_human_tuned_copper.yaml --output-dir femem_out
```

The test below shows the console-script line makes the same
`run_study_file` call as the Python block, without solving it; the solve
itself is tested through the Python block only.

## What you should see

```text tuned-copper-yaml-output
SAR10g(c*) = 4.8529e-03 W/kg
```

The test also reads the run's node S, its `Y_in` and all 18 lattice values
against their records, and the mean |B₁⁺| (`ccw`) over 41 mesh sample points
read on the block's own `Study` — those points are mesh cell centroids, so
they stay out of the study file.

## What this does not show

No tuning-accuracy, `f₀`, Q or match claim (`c_in_f` is an input record); no
compliance or absolute SAR claim. c\* is the argmax of a fixed lattice of ball
centres, not a global SAR maximum. Degree 1 only. No Touchstone file: the
node S is not a reference-impedance network.

## How this page is tested

The study file above must equal the tracked file byte for byte, and the
Python block is extracted verbatim and executed at 8 ranks by
[`test_wf17_tuned_copper_yaml_tutorial.py`](../../tests/package/test_wf17_tuned_copper_yaml_tutorial.py);
its readings must reproduce those of
[`test_wf1_yaml_f_human_solve.py`](../../tests/package/test_wf1_yaml_f_human_solve.py)
and the records of
[`test_wf16_study_tuned_f_human_b1_sar.py`](../../tests/package/test_wf16_study_tuned_f_human_b1_sar.py).
