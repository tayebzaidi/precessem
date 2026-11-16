# Backends: DolfinX by default, Palace on request

A **backend** is the field solver behind a `Study`. The `Study`, its YAML file
and `femem run` describe the problem once; the backend turns that description
into fields and port quantities.

- **DolfinX (the default).** Every gated record in this repository was produced
  on it. It is the reference, and it is the coupling host: SAR, `|B₁⁺|`,
  probes, tuning and the thermal hand-off all read DolfinX fields.
- **Palace (opt-in, for scale).** Palace (AWS Labs) is an open-source frequency-domain FEM solver. The package writes a Palace
  case and reads Palace's port output back; it never imports or links Palace.
  Today the Palace backend returns the four-port S-matrix of the F-small
  birdcage, nothing else.

## The three stages

The `Study` does not launch Palace. A Palace-backed solve runs in three
stages: the `Study` writes the case, a separate run executes Palace on it, and
the same call reads the result back after checking the case's checksum.

From Python, the first call writes the case into `case_dir` and raises
`PalaceRunRequired`; the exception carries the Palace command and the
directory to run it in. Calling `solve` again with the same `case_dir`, once
Palace has written its output, returns a `PalaceStudyResult`:

```python palace-backend
from mpi4py import MPI

from fem_em_solver.backends.palace import PalaceRunRequired

# `study` is a four-port F-small Study (see api/Study.md).
try:
    result = study.solve(MPI.COMM_WORLD, backend="palace", case_dir="palace_case")
except PalaceRunRequired as pending:
    print(pending.cwd, pending.argv)
```

From the command line, the same three stages are: `femem run` writes the case
and exits with code 3 ("external run pending"); Palace runs in the case
directory; the same `femem run` reads the output and writes its summary into
the output directory:

```bash palace-backend-cli
femem run tests/package/wf1/f_small_four_port.yaml --output-dir out --backend palace
cd out/palace_case    # run the Palace command below here, on the Palace service
femem run tests/package/wf1/f_small_four_port.yaml --output-dir out --backend palace
```

The Palace command the package asks for, by default, is:

```text palace-command
palace -np 8 case.json
```

`Study.solve(..., palace_ranks=N)` changes the rank count in that command.

## What is refused, and why

- **A parallel communicator.** The Palace backend writes and reads on one rank;
  `Study.solve(backend="palace")` on more than one rank raises. Palace's own
  parallelism is the rank count in its command.
- **A missing `case_dir`**, an unknown backend name, and any geometry,
  conductor, excitation or element degree the writer does not support (only
  the gapped four-port birdcage with lumped ports, degree 1 or 2) raise before
  anything is written.
- **An edited case.** The writer never overwrites a case file. If a file in
  `case_dir` differs from what today's writer would produce, the call raises
  `FileExistsError` naming the file, so Palace output is never read against a
  case it was not computed on.
- **`--backend palace` with `--ranks` or `--estimate`** on `femem run`.
- **A study file asking for an output Palace has no route for.** On
  `femem run --backend palace`, the `outputs.sar`, `outputs.b1_plus` and
  `outputs.probes` sections, a `tune:` section, `study.sweep_hz` and a
  `thermal:` section are each refused by name, as are a `conductor:` and
  `capacitor_sheets` excitation. The error names the key; nothing is
  dropped silently.
- **A study the Palace writer cannot take, on the command line.** On
  `femem run --backend palace`, a hole-wall birdcage (`geometry.as_hole`) and
  a launch on more than one rank are refused by name with exit code 2, before
  any mesh is built and before anything is written. An element degree other
  than 1 or 2 is refused when the file is read (`study.degree`). The checks
  are in `tests/package/test_arch1_u79_palace_inputs.py`.
- **Launching Palace from the `Study`.** Palace lives in its own container
  image, and the solver container has no Docker socket. A `Study` that starts
  Palace itself would need either Palace inside the solver image or a socket
  mount; both are deployment decisions, not package code.

## What the Palace backend carries

The YAML route honours one `outputs` key and returns one result; the lists
below are the code's own (`tests/package/test_arch1_u76_palace_outputs.py`
compares them with `fem_em_solver.study_yaml`):

```text palace-outputs
honoured: outputs.directory
returns: s_matrix
refused: outputs.b1_plus outputs.probes outputs.sar study.sweep_hz thermal tune conductor capacitor_sheets
```

`femem run --backend palace` refuses every other output by name, with exit
code 2, before anything is written.

## The binding condition

Every Palace-backed result names a **reduced cross-check at matched element
order on a shared mesh** against the DolfinX backend. A Palace number without
that cross-check is not a result of this package.

## The cross-check of record

The Palace backend is held to the same input giving the same answer. On the
×1 exact-strip F-small birdcage at 64 MHz, the S-matrix read back through
`Study.solve` (degree 2) and through `femem run` (degree 1) must equal the
recorded Palace read of the same case, Frobenius-relative, within the
tolerance below (rendered from `LEG_D0_REPRODUCTION_BAND` and checked by
`tests/package/test_wf17_backends_page.py`). The degree-2 read also ties per
class to our own degree-2 DolfinX record, and must differ from Palace's
degree-1 read, so the read is not blind to the case.

```text palace-tolerance
Frobenius-relative S tolerance: 1e-09
```

The nodes that hold it:

- `tests/package/test_arch1_palace_backend.py::test_u16_backend_refusals` —
  the refusals above.
- `tests/package/test_arch1_palace_backend.py::test_u16_case_equals_ops61` —
  the written case equals the matched-order Palace inputs of record (config
  keys; mesh byte for byte).
- `tests/package/test_arch1_palace_backend.py::test_u16_read_path_on_recorded_output`
  — the reader on recorded Palace output.
- `tests/package/test_arch1_palace_backend.py::test_u16_palace_read_back` — the
  degree-2 read back through `Study.solve`, the tie and the order control.
- `tests/package/test_arch1_palace_backend.py::test_u16_cli_read_path_on_recorded_output`
  — the CLI reader on recorded Palace output.
- `tests/package/test_arch1_palace_backend.py::test_u16_cli_palace_read_back`
  — the degree-1 read back through `femem run`.

They are registered under [`Study`](validation.md#study) and
[`run_study_file`](validation.md#run_study_file) on the
[validation-envelope page](validation.md); see also the API pages for
[`Study`](api/Study.md) and [`run_study_file`](api/run_study_file.md). The
matched-order Palace comparison these nodes reproduce is the public
"Palace at matched order" section of
[`ans:4`'s comparison](../../examples/ansys_benchmarks/ans4_birdcage_four_port_10_64_128MHz/COMPARISON.md).

**Not claimed:** a `Study` that launches Palace, any output other than S, any
fixture other than the F-small birdcage, or accuracy beyond the cross-check
above.
