# Estimating memory and time before a solve

## What you will compute

For one study file: the mesh it builds and the number of unknowns its
function space carries, plus a predicted peak memory and wall time for the
solve at the current launch width. Nothing is solved; the mesh and the
function space are built and then discarded.

## Before you start

You need the complex DolfinX build (`source /usr/local/bin/dolfinx-complex-mode`)
and **2 MPI ranks**; each block below takes ≈ 1 min (the mesh is built, the
solve is not). Copy the study file next to your script:
[f_small_four_port.yaml](../../tests/package/wf1/f_small_four_port.yaml).

## Run it

The block loads the study file and asks for its estimate on the ranks the
script was launched with.

```python estimate-python
from mpi4py import MPI
from fem_em_solver import load_study

study = load_study("f_small_four_port.yaml")
est = study.estimate(MPI.COMM_WORLD)
print(est.cells, est.unknowns, est.degree, est.ranks)
```

The same estimate from the command line:

```bash estimate-cli
mpiexec -n 2 femem run --estimate f_small_four_port.yaml
```

## What you should see

The estimate carries these fields, in this order:

```text estimate-fields
cells
unknowns
degree
ranks
predicted_peak_gib
predicted_seconds
rows_used
not_priced
```

The command line prints one `name: value` line per field, on rank 0 only.
The test, not the page, checks the values against the recorded mesh and
unknown counts.

## What the estimate is built from

The prices come from measured solves. For each element degree a power law is
fitted to those measurements: peak memory against the unknown count, and wall
time multiplied by the rank count against the unknown count. The predicted
time is then divided by the current launch width. A study at a degree with
fewer than two measured solves is refused by name, and so is a thermal-only
study, which has no electromagnetic solve to price. `--ranks` together with
`--estimate` is refused, because the estimate prices the width the command
was launched with.

## What this does not show

- How accurate an estimate is on a study unlike the measured solves; the
  leave-one-out bracket is the gate in
  [test_wf13_study_estimate.py](../../tests/package/test_wf13_study_estimate.py).
- Degree 2 on this file.
- Human-scale studies (priced; see
  [test_wf13_estimate_f_human.py](../../tests/package/test_wf13_estimate_f_human.py)).
- The Palace backend (refused with `--estimate`).
- Any external-solver comparison.

## Worked examples

- [mri:12, a study estimate](../../examples/mri/12_study_estimate.md)

## How this page is tested

[tests/package/test_wf17_estimate_tutorial.py](../../tests/package/test_wf17_estimate_tutorial.py)
extracts the Python block, execs it at 2 ranks in the complex build with the
solve disabled, runs the command line in the same process, and checks that the
two agree field by field, the way
[test_wf13_cli_estimate.py](../../tests/package/test_wf13_cli_estimate.py)
does.
