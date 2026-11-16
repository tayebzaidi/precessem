# Four-port S-parameters of a small birdcage

## What you will compute

The 4×4 scattering matrix S of the F-small four-leg birdcage at 64 MHz, with a
gelled-saline phantom (σ = 0.5 S/m, εᵣ = 78) inside it. Each of the four legs
carries one lumped port on an exact rectangular sheet; every port is driven in
turn and terminated in 50 Ω otherwise. One object,
[`Study`](../../src/fem_em_solver/study.py), composes the geometry, the
material, the port drive and the time-harmonic sweep, and returns S and Z.

## Before you start

- The `fem-em-solver` Docker service must be up
  (`docker compose -f docker/docker-compose.yml up -d`), and you work inside it
  with the repository at `/workspace`.
- The solve is complex-valued: switch the container shell to the complex
  DolfinX build first with `source /usr/local/bin/dolfinx-complex-mode`.
- Cost: ≈ 1 min wall time on 2 MPI ranks.
- **Set `degree=1`.** `Study` defaults to degree 2, which on this coil needs
  ≈ 8 ranks; the recorded S below is the degree-1 matrix on 2 ranks.

## The study

```python four-port-s
from fem_em_solver import (
    FourPortBirdcageGeometry,
    GelledSalinePhantomMaterial,
    LumpedPortDrive,
    Study,
)

geometry = FourPortBirdcageGeometry.f_small()  # exact port rectangles
phantom = GelledSalinePhantomMaterial(0.5, 78.0, 64e6)  # sigma S/m, eps_r, Hz
drive = LumpedPortDrive()  # 50 ohm on every port

study = Study(geometry, phantom, drive, frequency_hz=64e6, degree=1)
result = study.solve()
s = result.s_matrix  # 4 x 4 complex, port order P1..P4
print(f"|S11| = {abs(s[0, 0]):.4f}")
print(f"|S21| = {abs(s[1, 0]):.4f}")
```

`result.z_matrix` holds Z, and `result.write_touchstone(path)` writes a
Touchstone file of S (pass a path of your choosing).

## Run it

Save the block as a Python script and run it inside the container, in complex
mode, on 2 ranks (`<script>` is the file you saved):

```bash
source /usr/local/bin/dolfinx-complex-mode
PYTHONPATH=/workspace/src mpiexec -n 2 python3 <script>
```

## What you should see

```text four-port-s-output
|S11| = 0.5901
|S21| = 0.4230
```

The full matrix reproduces the recorded degree-1 S to ≈ 1e-10 relative (the
test holds it to 1e-9). By the coil's four-fold symmetry the other diagonal
entries match |S11| and the other adjacent-port entries match |S21|.

## What this does not show

No tuning, no resonance and no absolute S₁₁ accuracy claim: this is one mesh at
degree 1, and S here is checked against its own record and against
reciprocity / passivity / symmetry identities elsewhere, not against an
external solver at this degree.

## How this page is tested

The block above is extracted from this page verbatim and executed through the
logging harness by
[`test_wf17_more_tutorials.py`](../../tests/package/test_wf17_more_tutorials.py);
its S must reproduce the degree-1 record of
[`test_wf16_study_four_port.py`](../../tests/package/test_wf16_study_four_port.py)
to 1e-9 relative.
