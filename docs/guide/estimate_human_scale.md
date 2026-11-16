# Estimating a human-scale study before you solve it

## What you will compute

For the copper F-human birdcage study file: the mesh it builds and the number
of unknowns its function space carries, at element degree 1 and at element
degree 2, plus the predicted peak memory and wall time of each electromagnetic
solve at the current launch width. Nothing is solved; the mesh and the
function space are built and then discarded.

## Before you start

You need the complex DolfinX build (`source /usr/local/bin/dolfinx-complex-mode`)
and **2 MPI ranks**. The human-scale mesh is large: building it takes
≈ 1–2 min at 2 ranks, and each block below builds it once. Copy the study file
next to your script:
[f_human_hole_copper_tune.yaml](../../tests/package/wf1/f_human_hole_copper_tune.yaml).

The F-small version of this page is
[Estimating memory and time before a solve](estimate_before_a_solve.md); it
explains the fields of the estimate.

## Run it

The first block prices the study file as it is, at degree 1.

```python estimate-f-human-d1
from mpi4py import MPI
from fem_em_solver import load_study

study = load_study("f_human_hole_copper_tune.yaml")
est = study.estimate(MPI.COMM_WORLD)
print(est.cells, est.unknowns, est.degree, est.rows_used, est.not_priced)
```

The second block writes a copy of the same file with the one line
`  degree: 1` changed to `  degree: 2`, and prices the copy.

```python estimate-f-human-d2
from pathlib import Path

from mpi4py import MPI
from fem_em_solver import load_study

comm = MPI.COMM_WORLD
copy = Path("f_human_hole_copper_tune_d2.yaml")
if comm.rank == 0:
    lines = Path("f_human_hole_copper_tune.yaml").read_text().splitlines(keepends=True)
    hits = [i for i, line in enumerate(lines) if line.rstrip("\n") == "  degree: 1"]
    assert len(hits) == 1
    lines[hits[0]] = "  degree: 2\n"
    copy.write_text("".join(lines))
comm.Barrier()

study = load_study(copy)
est = study.estimate(comm)
print(est.cells, est.unknowns, est.degree, est.rows_used, est.not_priced)
```

The same estimate from the command line, on either file:

```bash estimate-f-human-cli
mpiexec -n 2 femem run --estimate f_human_hole_copper_tune.yaml
```

`--ranks` together with `--estimate` is refused: the estimate prices the width
the command was launched with.

## What you should see

Each block prints the cell count, the unknown count, the degree, the measured
rows the price was fitted from, and what is not priced. The two blocks build
the same mesh, so the cell counts agree; the degree-2 unknown count is larger.
The test, not the page, checks the values against the recorded human-scale
mesh and unknown counts and against the measured row each degree is priced
from. The predicted seconds and peak memory are the price table's fit, not a
measurement: no 2-rank human-scale solve exists to compare them against.

## The clinical coil at degree 2

The 0.70 m bore clinical birdcage at 1.5 T is priced differently: no mesh is
built, and the numbers are read from one measured run, not fitted. Copy its
study file next to your script:
[f_clinical_15t.yaml](../../tests/package/wf13/f_clinical_15t.yaml). This block
writes a copy at degree 2, in the same way as the block above, and prices it in
seconds at any width.

```python estimate-f-clinical-d2
from pathlib import Path

from mpi4py import MPI
from fem_em_solver import load_study

comm = MPI.COMM_WORLD
copy = Path("f_clinical_15t_d2.yaml")
if comm.rank == 0:
    lines = Path("f_clinical_15t.yaml").read_text().splitlines(keepends=True)
    hits = [i for i, line in enumerate(lines) if line.rstrip("\n") == "  degree: 1"]
    assert len(hits) == 1
    lines[hits[0]] = "  degree: 2\n"
    copy.write_text("".join(lines))
comm.Barrier()

study = load_study(copy)
est = study.estimate(comm)
print(est.cells, est.unknowns, est.degree, est.ranks, est.predicted_seconds,
      est.predicted_peak_gib, est.rows_used, est.not_priced)
```

The block prints the cell count, the unknown count, the degree, the 16 ranks
the row was measured at, its seconds with the mesh build excluded, its peak
memory, the row's name and what is not priced. `femem run --estimate` prints the
same on the copy. Without `--estimate`, `femem run` refuses the copy before any
mesh and names the tier a solve needs: `xl` for a root grid of up to 3 nodes,
`xxl` for a column set.

## What this does not show

- A solve, or tuning the coil: the solve at 8 ranks is
  [The tuned copper birdcage from a YAML file](tuned_copper_birdcage_yaml.md).
- A thermal price, or any price for the 3 T clinical coil.
- How accurate the prediction is at 2 ranks; the leave-one-out checks are in
  [test_wf13_estimate_f_human.py](../../tests/package/test_wf13_estimate_f_human.py)
  and
  [test_wf13_estimate_f_human_d2.py](../../tests/package/test_wf13_estimate_f_human_d2.py).
- The thermal transient of a human-scale study (stated as not priced).
- Any output (SAR, B1+) of the human-scale geometry: those are refused on it.
- Any external-solver comparison.

## How this page is tested

[tests/package/test_wf17_estimate_human_tutorial.py](../../tests/package/test_wf17_estimate_human_tutorial.py)
extracts both Python blocks, execs them at 2 ranks in the complex build with
every solve disabled, runs the command line on both files in the same
process, and checks that the two agree field by field, the way
[test_wf13_cli_estimate.py](../../tests/package/test_wf13_cli_estimate.py)
does. It also execs the clinical block, checks every printed field against the
measured row, checks that `femem run` refuses the copy with both tiers named,
and checks that the same block left at degree 1 misses the row's unknowns.
