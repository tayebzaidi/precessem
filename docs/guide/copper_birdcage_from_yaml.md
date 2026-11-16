# Copper rings from YAML: the F-small four-port with a Leontovich conductor

## What you will compute

The F-small four-port birdcage (exact port rectangles, κ off, degree 1)
around its saline phantom, with its rings and legs cut out of the mesh
(`as_hole: true`). The cut surfaces carry a Leontovich surface impedance for
copper instead of a perfect conductor. The study is solved at 64 MHz, one port
at a time, and the output is the 4×4 S-matrix in a Touchstone file plus the
JSON summary. The whole study runs through one call,
`fem_em_solver.run_study_file`.

## Before you start

You need the complex DolfinX build (`source /usr/local/bin/dolfinx-complex-mode`)
and **2 MPI ranks**; the block below takes ≈ 1 min. Copy the study file
next to your script:
[f_small_copper.yaml](../../tests/package/wf1/f_small_copper.yaml).

## The study file

| File | What it sets | Summary keys it adds | Record it carries |
| --- | --- | --- | --- |
| [f_small_copper.yaml](../../tests/package/wf1/f_small_copper.yaml) | `conductor: {kind: leontovich, sigma_s_per_m}` on the hole geometry (`geometry.as_hole: true`) of the four-port birdcage | `conductor` | `WF-1` step 4: the S-matrix class means reproduce the `ANS-6` degree-1 copper record |

A conductor needs `geometry.as_hole: true`, is refused on the sphere, and is
refused together with `study.sweep_hz` and with `--backend palace`; each
refusal names the field it refuses.

## Run it

The block runs the file and prints the conductor section of the summary and
the name of the Touchstone file it wrote.

```python copper-yaml
from fem_em_solver import run_study_file

run = run_study_file("f_small_copper.yaml", output_dir="femem_out/f_small_copper")
print(run.summary["conductor"], run.paths["touchstone"].name)
```

The same run from the command line:

```bash copper-yaml-cli
femem run f_small_copper.yaml --output-dir femem_out/f_small_copper
```

## What you should see

The S-matrix class means (self, adjacent, opposite) against the degree-1
copper record (relative distance per class), and the reciprocity, passivity
and class-spread gates on the run's own S-matrix. The bands are:

```text copper-yaml-bands
S class means vs the copper record: 1e-09
Reciprocity, relative: 1e-03
Largest singular value above 1: 1e-09
Class spread: 5e-03
```

The test, not the block, computes the readings from the run and checks them
against these bands. The same file without its `conductor:` section (a
perfectly conducting hole wall) misses the copper record; the test runs that
control itself.

## What this does not show

- Degree 2 (its copper records are wider-rank records and are not run here).
- Solid copper, which needs no `conductor:` section.
- A copper frequency sweep (refused).
- SAR, `|B₁⁺|` or probes on copper: the files
  `tests/package/wf1/f_small_copper_sar_b1.yaml` and
  `tests/package/wf1/f_small_copper_em_probes.yaml` exist and are gated
  elsewhere.
- A tuning, resonance or Q statement.
- The power identity on the copper surface (gated in `WF-1` step 4, not re-run
  here).
- Any external-solver comparison.

## Worked examples

- [mri:16, the same copper birdcage from Python](../../examples/mri/16_study_copper_birdcage.md)
- [mri:09, a study from YAML](../../examples/mri/09_study_from_yaml.md)

## How this page is tested

[tests/package/test_wf17_copper_yaml_tutorial.py](../../tests/package/test_wf17_copper_yaml_tutorial.py)
extracts the block, execs it at 2 ranks in the complex build, and recomputes
the class-mean tie and the three gates from the run at the bands above, the
way the gate
[test_wf1_yaml_copper.py](../../tests/package/test_wf1_yaml_copper.py)
does.
