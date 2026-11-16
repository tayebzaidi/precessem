# Capacitor sheets with ESR from YAML on the copper F-small birdcage

## What you will compute

The same study as [capacitor sheets with ESR](capacitor_sheets_with_esr.md),
written as a study file instead of Python: the F-small four-port birdcage
(exact port rectangles, κ off, degree 1) around its saline phantom, with its
rings and legs cut out of the mesh and carrying a Leontovich surface impedance
for copper. One leg port, P2, is replaced in the model by a series R–C sheet,
a capacitor with its equivalent series resistance (ESR). The remaining sheets
stay lumped ports and are driven in turn, so the study's S-matrix is 3×3 over
P1, P3 and P4. The whole study runs through one call,
`fem_em_solver.run_study_file`.

## Before you start

You need the complex DolfinX build (`source /usr/local/bin/dolfinx-complex-mode`)
and **2 MPI ranks**; the block below takes ≈ 1 min. The study file is
[f_small_copper_capacitors_esr.yaml](../../tests/package/wf1/f_small_copper_capacitors_esr.yaml);
the block reads it by its repository path, so run it from the repository root.

## The study file

```yaml capacitor-esr-yaml
# WF-17 step 15 (§9 item 627): f_small_copper.yaml with one 0.2 ohm ESR capacitor sheet on P2; record: WF-17 step 14's (U84) gated per-sheet circuit-reduction identity.
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
  kind: capacitor_sheets
  c_f: 6.629e-11
  esr_ohm: 0.2
  sheets: [P2]
conductor:
  kind: leontovich
  sigma_s_per_m: 58000000.0
```

The `excitation:` section is `capacitor_sheets`: `c_f` in farads, `esr_ohm`
in ohms, and `sheets`, the port ids whose sheets become capacitors. The key is
`sheets`, not `on`: YAML reads a bare `on` as `true`, so the loader refuses it.
Everything else is the copper file of
[copper rings from YAML](copper_birdcage_from_yaml.md).

## Run it

```python capacitor-esr-run
from fem_em_solver import run_study_file

run = run_study_file("tests/package/wf1/f_small_copper_capacitors_esr.yaml")
print(run.summary["kept_ports"], run.summary["capacitors"])
```

A capacitor-sheet study writes no Touchstone file: the JSON summary carries the
reduced S-matrix over `kept_ports` and a `capacitors` block.

## What you should see

The run's 3×3 S-matrix equals the circuit reduction of the unterminated 4×4
(every sheet a 50 Ω lumped port) with Z = r + 1/(jωC) terminating P2, to the
first band below. The unterminated 4×4 itself ties to the degree-1 copper
record (relative distance of its self / adjacent / opposite class means) to
the second. The bands, and the resistance on the sheet, are:

```text capacitor-esr-yaml-bands
Study S vs the circuit reduction with Z = r + 1/(jωC): 1e-03
Unterminated 4×4 vs the copper record: 1e-09
Series resistance r (Ω): 2e-01
```

The test, not the block, computes the readings from the run and checks them
against these bands. The same file with `esr_ohm: 0.4` misses the 0.2 Ω
reduction; the test runs that control itself, which shows the ESR in the file
reaches the solve.

## What this does not show

- More than one capacitor sheet.
- A tuned capacitance: `c_f` is an input.
- A resonance, match, Q or `R_in` statement.
- The power identity with the sheet's dissipation counted.
- Degree 2.
- Capacitor sheets on the human-scale birdcage (refused).
- The other capacitor file,
  `tests/package/wf1/f_small_hole_pec_capacitors.yaml`: it is **lossless**
  (`esr_ohm: 0.0`) and on the perfectly conducting hole wall, not on copper,
  and is not run here.
- Any external-solver comparison.

## Worked examples

- [mri:26, capacitor sheets from `femem run`](../../examples/mri/26_femem_run_capacitor_sheets.md)
- [mri:24, a Study tuned, then capacitor sheets](../../examples/mri/24_study_tune_then_capacitor_sheets.md)

## How this page is tested

[tests/package/test_wf17_capacitor_esr_yaml_tutorial.py](../../tests/package/test_wf17_capacitor_esr_yaml_tutorial.py)
compares the YAML block with the tracked file byte for byte, execs the run
block at 2 ranks in the complex build, and recomputes the circuit-reduction
identity and the copper tie at the bands above, the way the gate
[test_wf16_study_capacitor_sheets.py](../../tests/package/test_wf16_study_capacitor_sheets.py)
does.
