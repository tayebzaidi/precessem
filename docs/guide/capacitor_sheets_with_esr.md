# Capacitor sheets with ESR on the copper F-small birdcage

## What you will compute

The F-small four-port birdcage (exact port rectangles, κ off, degree 1) around
its saline phantom, with its rings and legs cut out of the mesh and carrying a
Leontovich surface impedance for copper. One leg port, P2, is replaced in the
model by a series R–C sheet: a capacitor with its equivalent series resistance
(ESR). The remaining sheets stay lumped ports and are driven in turn, so the
study's S-matrix is 3×3 over P1, P3 and P4. Everything is built from Python
with one `Study` and one `CapacitorSheets` excitation.

## Before you start

You need the complex DolfinX build (`source /usr/local/bin/dolfinx-complex-mode`)
and **2 MPI ranks**; the block below takes ≈ 1 min.

## Run it

The block builds the study, solves it, and prints the kept ports and the
sheet's impedance Z = r + 1/(jωC).

```python capacitor-esr
from fem_em_solver import (CapacitorSheets, FourPortBirdcageGeometry,
                           LeontovichWall, Study)
from fem_em_solver.materials import GelledSalinePhantomMaterial

study = Study(
    FourPortBirdcageGeometry(as_hole=True, strip_fraction=0.5),
    GelledSalinePhantomMaterial(0.5, 78.0, 64.0e6),
    CapacitorSheets(6.629e-11, esr_ohm=0.2, on=("P2",)),
    frequency_hz=64.0e6,
    degree=1,
    conductor=LeontovichWall(5.8e7),
)
result = study.solve()
print([p.port_id for p in result.port_definitions], result.capacitor_impedance_ohm)
```

`CapacitorSheets(c_f, esr_ohm, on)` takes the capacitance in farads, the
series resistance in ohms and the port ids whose sheets become capacitors; the
other sheets keep the lumped-port drive. The capacitance is an input here, not
a tuning result.

## What you should see

The study's 3×3 S-matrix equals the circuit reduction of the unterminated 4×4
(every sheet a 50 Ω lumped port) with Z = r + 1/(jωC) terminating P2, to the
first band below. The unterminated 4×4 itself ties to the degree-1 copper
record (relative distance of its self / adjacent / opposite class means) to
the second. The bands, and the resistance on the sheet, are:

```text capacitor-esr-bands
Study S vs the circuit reduction with Z = r + 1/(jωC): 1e-03
Unterminated 4×4 vs the copper record: 1e-09
Series resistance r (Ω): 2e-01
```

The test, not the block, computes the readings from the solved study and checks
them against these bands. The same study with the resistance doubled on P2
misses the 0.2 Ω reduction; the test runs that control itself, which shows the
ESR reaches the solve.

## Why every capacitor carries 0.2 Ω

Real RF capacitors are lossy. A lossless capacitor overstates the Q of a copper
coil, whose own conductor loss is small enough that the capacitor's ESR
dominates it. The project therefore puts r = 0.2 Ω on every capacitor in new
studies; lossless results stay only as labelled records. A uniform r adds r to
every mode's eigen-impedance, so the mode capacitance `C_mode` is unchanged.

On the human-scale copper birdcage the same ESR is applied at the circuit layer
instead of as an in-model sheet:
[tuned copper birdcage](tuned_copper_birdcage.md). The four-port S-matrix
without capacitors is on [four-port S](four_port_s.md).

## What this does not show

- More than one capacitor sheet (the three-sheet route is printed, not
  asserted, in the gate).
- A tuned capacitance: `c_f` is an input.
- A resonance, match, Q or `R_in` statement.
- The power identity with the sheet's dissipation counted (gated in `WF-16`
  step 10, not re-run here).
- Degree 2.
- Capacitor sheets on the human-scale birdcage (refused).
- The YAML route: the file
  `tests/package/wf1/f_small_hole_pec_capacitors.yaml` is lossless and on the
  perfectly conducting hole wall, not on copper.
- Any external-solver comparison.

## Worked examples

- [mri:26, capacitor sheets from `femem run`](../../examples/mri/26_femem_run_capacitor_sheets.md)
- [mri:24, a Study tuned, then capacitor sheets](../../examples/mri/24_study_tune_then_capacitor_sheets.md)

## How this page is tested

[tests/package/test_wf17_capacitor_esr_tutorial.py](../../tests/package/test_wf17_capacitor_esr_tutorial.py)
extracts the block, execs it at 2 ranks in the complex build, and recomputes
the circuit-reduction identity and the copper tie at the bands above, the way
the gate
[test_wf16_study_capacitor_sheets.py](../../tests/package/test_wf16_study_capacitor_sheets.py)
does.
