# An implant rod in the coil-driven phantom

## What you will compute

The 4×4 **S-matrix** of the F-small birdcage (exact port rectangles, κ off,
degree 1, 64 MHz) with a thin **conducting rod** inside its saline phantom,
through a `Study`, and the 1 g mass-averaged SAR at the two rod tips. The rod is
a hole in the mesh with a perfectly conducting wall. You will do it twice, from
Python and from a study file, and both give the same S-matrix.

## Before you start

You need the complex DolfinX build (`source /usr/local/bin/dolfinx-complex-mode`)
and **2 MPI ranks**. The Python block takes ≈ 2.5 min and the study file ≈ 1 min.

## Run it from Python

```python implant-rod
from fem_em_solver import FourPortBirdcageGeometry, LumpedPortDrive, Study
from fem_em_solver.materials import GelledSalinePhantomMaterial

rod = {"radius": 1.0e-3, "length": 0.030, "centre": (0.015, 0.0, 0.0), "axis": (0.0, 0.0, 1.0)}
study = Study(
    FourPortBirdcageGeometry(implant_rod=rod),
    GelledSalinePhantomMaterial(0.5, 78.0, 64.0e6),
    LumpedPortDrive(),
    frequency_hz=64.0e6,
    degree=1,
)
result = study.solve()
print(f"cells = {result.cells}")

tips = [(0.015, 0.0, 0.016), (0.015, 0.0, -0.016)]
sar = study.sar(mass_g=1.0, centres=tips)
print(sorted(sar), [len(per_centre) for per_centre in sar.values()])
```

`radius` and `length` are in metres, and `length` is the cylinder between the
centres of the two rounded caps. The two `tips` are the caps' apexes, 1 mm past
those centres. `Study.sar` is per 1 V drive, one entry per port.

## Run it from a study file

Save this as `implant_rod.yaml`:

```yaml implant-rod-file
femem: 1
study:
  frequency_hz: 64000000.0
  degree: 1
geometry:
  kind: four_port_birdcage
  preset: f_small
material:
  sigma: 0.5
  epsilon_r: 78.0
excitation:
  kind: lumped_port_drive
implant:
  rod:
    radius_m: 0.001
    length_m: 0.03
    centre_m: [0.015, 0.0, 0.0]
    axis: [0.0, 0.0, 1.0]
```

Then run it:

```text
mpiexec -n 2 femem run implant_rod.yaml
```

The run writes a JSON summary carrying the S-matrix, and `implant_rod.s4p`,
into `femem_out/implant_rod/`. All four `rod` keys are required. `implant:` is
refused on any other geometry and beside `load_result:`, and the Palace backend
refuses a rod.

## What you should see

After the mesher's own output, at 2 ranks, the Python block prints:

```text implant-rod-expected
cells = 168507
```

The cell count can change with the mesher's version. The S-matrix from the
Python block, from the study file, and from the fixture built directly without a
`Study` agree to this relative band:

```text implant-rod-bands
Python vs YAML vs the direct build: 1e-12
```

The test, not the block, computes these readings. The rod changes the S-matrix
by far more than that band, and the test checks that too.

## With capacitor sheets

`CapacitorSheets` closes one port sheet with a series capacitor and its ESR —
here sheet `P2`, 66.29 pF with 0.2 Ω — and the solve returns the 3×3 S-matrix of
the other three ports. The rod is the same **perfectly conducting** rod as
above, on the same mesh:

```python implant-rod-caps
from fem_em_solver import CapacitorSheets, FourPortBirdcageGeometry, Study
from fem_em_solver.materials import GelledSalinePhantomMaterial

rod = {"radius": 1.0e-3, "length": 0.030, "centre": (0.015, 0.0, 0.0), "axis": (0.0, 0.0, 1.0)}
study = Study(
    FourPortBirdcageGeometry(implant_rod=rod),
    GelledSalinePhantomMaterial(0.5, 78.0, 64.0e6),
    CapacitorSheets(6.629e-11, esr_ohm=0.2, on=("P2",)),
    frequency_hz=64.0e6,
    degree=1,
)
result = study.solve()
print(f"cells = {result.cells}; ports = {[p.port_id for p in result.port_definitions]}")
```

At 2 ranks it prints:

```text implant-rod-caps-expected
cells = 168507; ports = ['P1', 'P3', 'P4']
```

Its 3×3 agrees with the same capacitor-terminated fixture built directly,
without a `Study`, to this relative band:

```text implant-rod-caps-bands
Python vs the direct build: 1e-12
```

This route is **degree 1 only**, and only `Study.solve` on the perfectly
conducting rod takes it. At degree 2 the same `Study` is refused before any
mesh:

```python implant-rod-caps-refused
from fem_em_solver import CapacitorSheets, FourPortBirdcageGeometry, Study
from fem_em_solver.materials import GelledSalinePhantomMaterial

rod = {"radius": 1.0e-3, "length": 0.030, "centre": (0.015, 0.0, 0.0), "axis": (0.0, 0.0, 1.0)}
study = Study(
    FourPortBirdcageGeometry(implant_rod=rod),
    GelledSalinePhantomMaterial(0.5, 78.0, 64.0e6),
    CapacitorSheets(6.629e-11, esr_ohm=0.2, on=("P2",)),
    frequency_hz=64.0e6,
    degree=2,
)
study.solve()  # raises ValueError
```

The two refusals read, first at degree 2 and then on the copper hole route
(`as_hole=True` with a `LeontovichWall` conductor):

```text implant-rod-caps-refusals
solve: implant_rod at degree 2 is not gated; it needs degree 1 (WF-15 step 12)
solve: implant_rod with CapacitorSheets is not gated (WF-15 step 12); only Study.solve on the PEC four-port at degree 1 is (WF-15 step 13)
```

## Heating the phantom around the rod

`Study.temperature` with `drive="P1"` heats the gel phantom with port `P1`
driven at 1 V on the rod `Study` above: the SAR of that one drive is the heat
source, the gel conducts it and cools through its outer boundary into air at
20 °C, and the rod's wall is adiabatic. It steps 40 Crank–Nicolson steps to
t\* = 1254 s and returns the volume mean and the 1 g ball means around the two
tips. Temperatures are rises above 20 °C, per watt dissipated in the
phantom: the block reads that power P from the EM solve and scales the heat
source by 1/P.

```python implant-rod-heating
from fem_em_solver import FourPortBirdcageGeometry, LumpedPortDrive, RobinCooling, Study, TransientSpec
from fem_em_solver.materials import GelledSalinePhantomMaterial
from fem_em_solver.post.sar import mean_sar

rod = {"radius": 1.0e-3, "length": 0.030, "centre": (0.015, 0.0, 0.0), "axis": (0.0, 0.0, 1.0)}
study = Study(
    FourPortBirdcageGeometry(implant_rod=rod),
    GelledSalinePhantomMaterial(0.5, 78.0, 64.0e6),
    LumpedPortDrive(),
    frequency_hz=64.0e6,
    degree=1,
)
transient = TransientSpec(dt_s=31.35, n_steps=40, rho_c_j_per_m3_k=4.18e6, theta=0.5,
                          startup_half_steps=2)
tips = [[0.015, 0.0, 0.016], [0.015, 0.0, -0.016]]
study.sar(mass_g=1.0, centres=tips, drives=["P1"])
fields = study.fields_result.fields["P1"]
P = mean_sar(fields.e_complex, sigma=fields.sigma_field, rho=study.density_kg_m3,
             cell_tags=study.fields_result.cell_tags,
             subdomain_ids=study.geometry.phantom_tag)["dissipated_power_w"]
cooling = RobinCooling(k_w_per_m_k=0.6, h_w_per_m2_k=10.0, t_inf_c=20.0, source_scale=1.0 / P)
heat = study.temperature(cooling, transient=transient, balls={"mass_g": 1.0, "centres": tips},
                         drive="P1")
print(f"cells = {heat.cells}")
for ball in heat.balls:
    print(ball["centre"], f"{ball['theta_k'][-1]:.6e} K")
```

The same run from the tracked study file
[tests/package/wf15/rod_heating.yaml](../../tests/package/wf15/rod_heating.yaml),
copied into the working directory. The tracked file omits
`thermal.source_scale`, so as tracked it reads per 1 V on `P1`; add
`source_scale: <1/P>` under `thermal:` to read per watt. The test runs exactly
that copy with the Python block's P:

```python implant-rod-heating-file
from fem_em_solver.study_yaml import run_study_file

run = run_study_file("rod_heating.yaml")
for ball in run.summary["temperature"]["balls"]:
    print(ball["centre"], f"{ball['theta_k'][-1]:.6e} K")
```

The two blocks' tip temperatures at t\* agree, and the run's energy balance
and the volume heating slope of one backward-Euler step hold, to these bands:

```text implant-rod-heating-bands
Python vs the study file, tip temperature at t*: 1e-12
energy balance (E): 1e-09
volume slope (S2): 1e-04
```

The test, not the block, computes these readings, and checks each tip against
the recorded value of the gated rod-heating test to its printed last digit.

## One step from the study file

The heating slopes are read from the study file too. `rod_heating_step.yaml`
is the tracked file with the `source_scale` line added as above, and its
`transient:` replaced by one backward-Euler step: `dt_s: 0.0015675`,
`n_steps: 1`, `theta: 1.0`, `startup_half_steps: 0`. The block prints the
volume mean rise and each tip ball's mean rise after that one step:

```python implant-rod-heating-step
from fem_em_solver.study_yaml import run_study_file

run = run_study_file("rod_heating_step.yaml")
t = run.summary["temperature"]
print(f"volume mean {t['volume_mean_k'][1]:.6e} K")
for ball in t["balls"]:
    print(ball["centre"], f"{ball['theta_k'][1]:.6e} K")
```

The test prints the step's tip and volume heating slopes beside the Python
route's one step (`study.temperature` with the same transient) and checks each
against the volume slope band above. Their agreement is printed, not gated: one
small step resolves the slopes only to the arithmetic of the step itself. The
study-file gate is the tip temperature at t\* from the file block above, which
the test checks equal to the Python route's to the identity band.

The rod file has no duty cycle: a `thermal.duty` key is refused as an unknown
key, and `Study.temperature(..., duty=...)` on the rod study raises
`ValueError`. [The malformed rod file](../../tests/package/wf15/malformed_rod_heating_duty.yaml)
carries one, and the check names that field and fails:

```bash implant-rod-heating-check
femem check malformed_rod_heating_duty.yaml
```

## What this does not tell you

- **No two-port drive:** every reading is one port driven at 1 V, so this is not
  the quadrature drive of a clinical coil.
- **No 1 g rung:** the 1 g average is read on this one mesh, with no refined
  rung beside it.
- **No published comparand:** the tip SAR, and its ratio to the same point with
  no rod, have no external value to agree with. The test prints the ratio as a
  reading with no record.
- **No tissue:** the phantom is uniform saline, and the rod is a perfect
  conductor with no coating.
- **No other route:** a rod `Study` at any degree but 1 is refused before any
  mesh, and so is a rod with `CapacitorSheets` anywhere but `Study.solve` on the
  perfectly conducting rod, because no test gates them. The capacitor section
  is not a tuned, rod-loaded coil, and the rod there is not copper. The
  heating section is degree 1, one port (`P1`) scaled to 1 W in the phantom,
  in gel, and makes no compliance claim.

## How this page is tested

[tests/package/test_wf17_implant_rod_tutorial.py](../../tests/package/test_wf17_implant_rod_tutorial.py)
execs the Python block and runs the study file through `femem run` at 2 ranks in
the complex build, and compares both S-matrices and the direct build at the band
above. It also runs the Python block with the rod removed and checks that it
misses, and checks the refusals named above. It also execs the capacitor block
and compares its 3×3 with the direct build at its band, and runs the degree-2
block and the copper hole route and checks that each raises the text quoted.
[tests/package/test_wf17_implant_rod_heating_tutorial.py](../../tests/package/test_wf17_implant_rod_heating_tutorial.py)
execs both heating blocks at 2 ranks and checks the bands of the heating
section, and checks that the heating block with the rod removed is refused
before any mesh.
[tests/package/test_wf17_implant_rod_heating_yaml_tutorial.py](../../tests/package/test_wf17_implant_rod_heating_yaml_tutorial.py)
runs the one-step study file and prints its slopes beside the Python route's,
re-checks the study file's tip temperature at t\* against the Python route's,
and checks that the rod route refuses a duty cycle in YAML and in Python.
