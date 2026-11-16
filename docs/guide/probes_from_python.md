# Temperature probes on the gel sphere from Python

## What you will compute

The gelled-saline sphere in its air box (the `MAT-12` fixture) is driven by the
quasi-static uniform field and solved once. Its absorbed power heats the sphere
under Robin cooling to the steady temperature rise θ = T − T∞. Four virtual
probes then read the ball-averaged θ: one ball that encloses the whole sphere,
and three small balls at the centre, off-axis, and near the surface. The
whole-sphere ball must read the volume mean of θ, which is what the probe's
quadrature is checked against.

This is the Python route of
[heating the gel sphere from YAML](gel_sphere_heating_from_yaml.md), whose
probes section runs the same study from a file.

## Before you start

You need the complex DolfinX build and **2 MPI ranks**; the run takes
≈ 30 s. The steady temperature this route reuses is described on
[the temperature page](temperature.md).

## Run it

```python probe-python
from fem_em_solver import QuasiStaticSphereField, RobinCooling, SphereInBoxGeometry, Study
from fem_em_solver.materials import GelledSalinePhantomMaterial
from fem_em_solver.thermal import ProbeBall

study = Study(
    SphereInBoxGeometry(
        sphere_radius=0.01,
        box_half_width=0.02,
        resolution_sphere=0.001,
        resolution_far=0.002,
        sphere_tag=1,
    ),
    GelledSalinePhantomMaterial(0.57, 78.0, 64.0e6),
    QuasiStaticSphereField(e0=1.0),
    frequency_hz=64.0e6,
    density_kg_m3=1000.0,
    degree=1,
)
result = study.solve()
probes = [
    ProbeBall((0.0, 0.0, 0.0), 0.02),
    ProbeBall((0.0, 0.0, 0.0), 0.002),
    ProbeBall((0.0035, 0.0, 0.0), 0.002),
    ProbeBall((0.0037527767497325675, 0.0037527767497325675, 0.0037527767497325675), 0.002),
]
readings = study.probe(probes, RobinCooling(0.6, 100.0, 20.0, source_scale=25000000.0))
print(readings.values, readings.volumes)
```

`study.probe` reuses the solved field: it computes the steady temperature once
and averages θ over each ball on the mesh cells this rank owns.

## What you should see

The block prints one θ value (kelvin) and one mesh volume (m³) per probe. The
test compares them with these bands:

```text probe-python-bands
Whole-ball probe vs the volume mean: 1e-10
Displaced probe must miss by more than: 1e-09
```

The test computes the readings, not the block. The first line bounds the
whole-sphere probe's value against ∫θ dx / V and its volume against V. The
second line is the control: the third small probe moved 0.1 mm outward must
read a different θ by more than that relative margin.

## What this does not show

- A probe on a coil-driven field: the field here is imposed, not a birdcage's.
- A transient or duty-cycled temperature (`study.probe` accepts `transient` and
  `duty`; this page does not run them).
- Degree 2, a perfused tissue model, or a container wall.
- The `-n 1` against `-n 2` width check, which the gate runs and this page
  does not.
- Any statement against a commercial solver.

## Worked examples

- [mri:13, a temperature from a config file](../../examples/mri/13_temperature_from_yaml.md)

## How this page is tested

[tests/package/test_wf17_probe_python_tutorial.py](../../tests/package/test_wf17_probe_python_tutorial.py)
runs the block above verbatim at 2 ranks and checks it against the probe
gate
[tests/package/test_wf15_study_probe.py](../../tests/package/test_wf15_study_probe.py).
