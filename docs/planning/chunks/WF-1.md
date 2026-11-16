# WF-1 — MRI example CLI/config: `femem run study.yaml`

Item briefs and execution history (operator directive 2026-10-01). The §7 row is the status of record.

### §9 item 256 — 2026-10-05 (chain U5)

**Outcome: ✅ DONE 2026-10-06 00:41Z** (19:30 CDT slot, attempt 1). Five windows green on first run; no in-slot fix. Execution record at the end of this file.

*Drafted by a general-purpose agent for the 2026-10-05 18:00 daily review; the review reviews it.* Licence: §10 chain U5 (PROJECT_PLAN.md:4739–4741), a weekly-enumerated step (§9 rule (4)). It is serial on U4, and U4 is met: `WF-16` ✅ (`ee9be842`), `WF-15` ✅ (`4239ba88`). The `WF-1` §7 row (PROJECT_PLAN.md:2876) has **no Done-when**, and the 10-03 weekly left it for the daily to write ("the daily writes the Done-when at rubric grain"). The Done-when under "Status it can move" below is that text. **The review copies it into the row in its own commit (rule (j)).** No band, tolerance, record or validation target moves. Every band below is imported by name from the file:line cited, and none is set here.

**`WF-1` step 1: a thin YAML layer over `Study`, and a `femem` console script.** `femem run <file>.yaml` parses a YAML study file into the same `Study` object that `WF-15` / `WF-16` build in Python. It solves the study and writes the outputs: a JSON summary, plus a Touchstone file for a port study.
- Two committed files drive the two gated fixtures at degree 1, `-n 2`, to the records those chunks already assert.
- Three malformed files fail with a named field.
- The installed console script refuses the same three.
- **Executor:** implementer.
- **What it adds:**
  - two new `src/` modules (YAML loader, CLI) plus a `__main__.py`;
  - one additive `[project.scripts]` table in `pyproject.toml`;
  - three `__all__` exports;
  - two new test modules and seven YAML files;
  - the registry and snapshot entries.
- **What it does not touch:** `study.py`, any existing gate module, helper or record.
- **Tier:** W1 standard (180 s); W2 heavy by ceiling (400 s, the ceiling `WF-16` W1 had on the same fixture at the same width); W3 / W4 the `WF-14` / `WF-13` windows verbatim; W5 smoke (installed console script, real build).
- **Ranks:** `-n 2` for every solve window, which is the width of both records. Complex build for W1 / W2.
- **Class:** chain step (solve, degree 1 only). Measured median ≈ 9 min (the 10-05 restock table).
- **Cost:** ≈ 50 min predicted:
  - 15 fixed + ≈ 25 code;
  - windows ≈ 1.7 + 2.5 + 0.1 + 0.6 + 0.5 min;
  - ≈ 5 spare.
- **Depends on:** nothing open. Nothing named `cli`, `femem`, `load_study` or `StudyConfigError` exists in `src/`. `pyproject.toml` has no `[project.scripts]` (read 2026-10-05). `pyyaml` is both a declared dependency (`pyproject.toml`) and installed in the image (`docker/Dockerfile:15`).
- **Shared files.** Items 259 / 260 may also edit `src/fem_em_solver/__init__.py`, `tests/package/public_surface.txt` and `tests/package/validation_envelope.toml`. Make every change there **by edit, adding lines only. On a rebase, keep both sides** (both sets of names, both sets of tables). This item does **not** edit `study.py`, so 259 / 260's `study.py` edits cannot conflict with it.

**Settled from the code (read before designing).**
- **The Python-built comparands, imported and not rebuilt by hand:**
  - Sphere: `tests/package/test_wf15_study_lossy_sphere.py`, `_study(degree=1)` (`:111–119`) on `GEOMETRY` (`:98–104`). That is σ 0.57 / εᵣ 78 (`_material_057`, `:107`), `QuasiStaticSphereField(e0=E0)`, `frequency_hz=FREQUENCY_HZ` and `density_kg_m3=RHO_KG_M3`.
  - Birdcage: `Study(FourPortBirdcageGeometry.f_small(), GelledSalinePhantomMaterial(SALINE_SIGMA, SALINE_EPSILON_R, FREQUENCY_64_HZ), LumpedPortDrive(), frequency_hz=FREQUENCY_64_HZ, degree=1)`, which is `test_wf16_study_four_port.py:139–145`.
- **The records and bands:**
  - Sphere: `MAT4_DEG1_RECORD` 8.291558128e-08 and `MAT4_DEG1_PRINT` `"8.291558e-08"` (`test_wf15…:87–88`), `L0_BAND` 1e-9 (`test_mat5_em_coupled_sphere.py:104`), and `MAT4_BAND` (`test_phantom_material_model.py:222`).
  - Birdcage: `ANS4_EXACT_STRIP_S_P1_RECORDS["64"]` / `ANS4_EXACT_STRIP_S_P2_X1["64"]` and `ANS4_EXACT_STRIP_RECORD_RANK_WIDTH` (`test_ans4_exact_strip_records.py`), `LEG_D0_REPRODUCTION_BAND` (`test_port_birdcage_four_port.py:165`), `ANCHOR_TOL` / `FREQ_ATOL_HZ` (`tests/io/test_touchstone_skrf_crosscheck.py:51–52`), and the `PORT-11` bands and helpers exactly as `test_wf16…:73–89` imports them.
- **What `WF-15` / `WF-16` measured on these comparisons today:**
  - Sphere, `-n 2`, degree 1: (I) Study against direct, in-process, **0.000e+00**; (R1) **6.222e-11**; (R2) the string matched (`20261005T170210Z_WF-15-r2-W1-deg1.log`).
  - Birdcage, `-n 2`, degree 1: (R1) **9.557e-11** against `P1`; (T1) **6.873e-14**; (T2) **9.557e-11**; (I) in-process **9.110e-16**; (c1) cross-degree **0.0491**; (c2) P1↔P2 swap **0.2492** (`20261005T202623Z_WF-16-W1-deg1.log`).
  - Degree 2 at `-n 8` (`20261005T203015Z_WF-16-W2-deg2.log`) is **not** re-run here (Scope).
- **"Bit-identical" cannot mean bitwise solver output.** `WF-16`'s in-process Study-against-direct read 9.110e-16, not 0. `-n 2` MUMPS is also not bitwise across processes (known-issues.md:69, `OPS-43`; `SOLVE_FLOOR_RECORD` at known-issues.md:81). The reading is `WF-15`'s rule-(j) reading (PROJECT_PLAN.md:2890):
  - **bitwise up to the solve**: the YAML-built `Study`'s spec equals the Python-built one field by field, with floats compared by `float.hex`, and the DG0 material fields assembled from it are equal;
  - **the solve outputs within the imported bands**: SAR at `L0_BAND`; S and the Touchstone file at `LEG_D0_REPRODUCTION_BAND`;
  - **plus the registered records the U3 / U4 tests already assert.**

  The Done-when states this as "(read by the 2026-10-05 18:00 review, rule (j); the weekly ratifies)".
- **The JSON round trip is exact.** `json.dump` writes a Python float by `repr`, and `json.load` reads it back to the same double. So the summary file's SAR and S equal the in-memory values bitwise, and the test may compare the file's numbers with `==` where it says so.
- **PyYAML's float resolver (YAML 1.1) needs a dot *and* a signed exponent.** Measured on the host's PyYAML 6.0.1, `yaml.safe_load` gives:
  - `64.0e6` → `str` and `1e-9` → `str`;
  - `6.4e+07`, `64000000.0` and `0.5700000000000001` → `float`;
  - `78` → `int` and `true` → `bool`;
  - for `a: 1\na: 2` → `{'a': 2}`: the duplicate key is dropped silently.

  The container's `pyyaml` is unpinned, so the test prints `yaml.__version__`.
- **The literals the good files carry.** The files write each float as its Python `repr`:
  - `repr(0.01 / 10.0)` is `0.001` and `repr(0.01 / 5.0)` is `0.002`, both `==` the decimal literals (host arithmetic, 2026-10-05);
  - `repr(64.0e6)` is `64000000.0`;
  - `repr(np.nextafter(0.57, 1))` is `0.5700000000000001`.

  (K) checks every bit, so a hand-typed literal that is one ULP off fails there, in the parse node, before any solve.
- **The material frequency is the study's.** `GelledSalinePhantomMaterial(sigma, epsilon_r, frequency_hz, mu_r=1.0)` is a frozen dataclass (`materials/phantom.py:21–40`). The YAML `material` section has **no** frequency key: the loader passes `study.frequency_hz`. A mismatch would raise at `core/time_harmonic.py:462–465`, and that is not one of the three gated failures.
- **`femem` will not be on `PATH` in the dev container.** The image's `pip install -e .` (`docker/Dockerfile:36`) predates the entry point, so in-container windows call `fem_em_solver.cli.main([...])` in-process. The console script is exercised only in `WF-13`'s fresh venv. That venv's `pip install --no-deps --no-build-isolation .` (`scripts/package/wf13_setup.sh:84`) builds a wheel from the copied `pyproject.toml` and does install `[project.scripts]`. **Never `pip install` into the live container.**

**Code.**
- (1) **New `src/fem_em_solver/study_yaml.py`** with three public names:
  - **`StudyConfigError(ValueError)`**, with attributes `.field` (the dotted path, e.g. `"geometry.sphere_raduis"`) and `.source` (the file path). Its message is `f"{source}: {field}: {reason}"`, so the field name appears in `str(exc)`.
  - **`load_study(path) -> Study`.** It parses with a `yaml.SafeLoader` subclass whose `construct_mapping` raises `StudyConfigError` on a duplicate key. Never use `yaml.load`. `import yaml` happens **inside** the function, so `import fem_em_solver` does not import `yaml`.
  - **`run_study_file(path, *, output_dir=None, comm=MPI.COMM_WORLD) -> StudyRun`.** It calls `load_study`, then `study.solve(comm)`, then writes the outputs. `StudyRun` (study, result, summary dict, written paths) is a dataclass and is **not exported**, as `PortStudyResult` was not.
  - **Schema v0. Units are SI and given by the key name; a unit string is a type error.**
    - `femem` (int, must be 1; required).
    - `study`: `frequency_hz` (number, required), `degree` (int in {1, 2}, optional, **default 2**, the Study's own default), `density_kg_m3` (number, optional, default the Study's).
    - `geometry`: `kind` (required).
      - `sphere_in_box`: `sphere_radius`, `box_half_width`, `resolution_sphere`, `resolution_far` (numbers, required), `sphere_tag` (int, optional).
      - `four_port_birdcage`: `preset: f_small` (optional), plus any `FourPortBirdcageGeometry` field name as an override. A list is converted to a tuple for `port_box_size`. Any other key, including `width_correction_kappa`, is an unknown key.
    - `material`: **either** `preset` (`low` / `mid` / `high`, which goes through `Study`'s own string route) **or** `sigma` + `epsilon_r` (+ optional `mu_r`), building `GelledSalinePhantomMaterial(sigma, epsilon_r, study.frequency_hz, mu_r)`. Both forms at once is an error on `material`.
    - `excitation`: `kind` (required).
      - `quasi_static_sphere_field`: `e0` (optional).
      - `lumped_port_drive`: `z0_ohm`, `port_impedance_ohm` (optional numbers), `drive_voltage_v` (optional: a number or a `[re, im]` pair).
    - `outputs` (optional): `directory`.
  - **Numbers.** A number field accepts `int` or `float`, converted by `float()`. It **rejects `bool`** (it is an `int` subclass) and `str`. The `str` message names the YAML-1.1 trap: "YAML reads `64.0e6` as a string; write `6.4e+07` or `64000000.0`".
  - **Validation order is fixed, and the first error raises:**
    - at top level: the version, then unknown top-level keys, then missing required sections;
    - then section by section in the order study → geometry → material → excitation → outputs;
    - within a section: unknown keys, then missing required keys, then types and values.
  - **No solver code, no mesh before the parse completes, and no import from `tests/`.** A `StudyConfigError` is raised identically on every rank, because every rank parses the same file. There is no collective before it.
  - **Outputs.** `run_study_file` resolves `output_dir` as the argument, else `outputs.directory`, else `./femem_out/<stem>`. After every reduction, **rank 0** creates it and writes `summary.json`. Then `comm.Barrier()`, and every rank returns. The summary holds:
    - `femem` 1, `source` (absolute path), `sha256` of the file, `family`, `degree`, `frequency_hz`, `ranks`, `cells`, `elapsed_s`;
    - sphere: the `mean_sar` dict with every value `float()`-ed;
    - ports: `s_matrix` and `z_matrix` as nested `[re, im]` pairs, and `touchstone` (the filename).

    For a port study it calls `result.write_touchstone(<output_dir>/<stem>.s4p)` (`study.py:405–421`, rank 0 plus `Barrier`).
- (2) **New `src/fem_em_solver/cli.py`.** It defines `main(argv=None) -> int` with `argparse` subcommands:
  - `run <file> [--output-dir DIR]` calls `run_study_file` and returns 0;
  - `check <file>` calls `load_study` only and prints, on rank 0, one line per spec field (family, degree, frequency, the material's σ / εᵣ, the excitation kind). It returns 0.

  Behaviour:
  - A `StudyConfigError` prints `femem: error: <message>` to stderr on rank 0 and returns **2** on every rank.
  - `main` **returns** an exit code and never calls `sys.exit`. The console-script wrapper does that.
- (3) **New `src/fem_em_solver/__main__.py`:** `import sys; from .cli import main; sys.exit(main())`. It is not exercised by a window, so it is disclosed as unexercised.
- (4) **`pyproject.toml`: add one table, `[project.scripts]` with `femem = "fem_em_solver.cli:main"`.** Nothing else changes.
- (5) **`src/fem_em_solver/__init__.py`:** import `StudyConfigError`, `load_study` and `run_study_file` from `.study_yaml`, and append them to `__all__` (30 → 33). `cli` is **not** imported by `__init__`.
- (6) **`tests/package/public_surface.txt`:** add the three names, keeping the sort order: `StudyConfigError` after `Study`, `load_study` before `load_touchstone`, `run_study_file` after `run_single_port_excitation_case`.
- (7) **`tests/package/validation_envelope.toml`:**
  - add `[names.StudyConfigError]` (the malformed node and the installed node), `[names.load_study]` (the good-file and malformed nodes) and `[names.run_study_file]` (the two CLI nodes);
  - append `test_wf1_cli_sphere_reproduces_wf15` to `[names.Study]`;
  - node IDs only, no number (`WF-14` (N)).
- (8) **Seven YAML files under `tests/package/wf1/`:**
  - **`sphere_mat4.yaml`:** `femem: 1`; `study: {frequency_hz: 64000000.0, degree: 1, density_kg_m3: 1000.0}`; `geometry: {kind: sphere_in_box, sphere_radius: 0.01, box_half_width: 0.02, resolution_sphere: 0.001, resolution_far: 0.002, sphere_tag: 1}`; `material: {sigma: 0.57, epsilon_r: 78.0}`; `excitation: {kind: quasi_static_sphere_field, e0: 1.0}`.
  - **`f_small_four_port.yaml`:** `femem: 1`; `study: {frequency_hz: 64000000.0, degree: 1}`; `geometry: {kind: four_port_birdcage, preset: f_small}`; `material: {sigma: 0.5, epsilon_r: 78.0}`; `excitation: {kind: lumped_port_drive}`.
  - **`control_sigma_ulp.yaml`:** `sphere_mat4.yaml` with `sigma: 0.5700000000000001` and no other change.
  - **`malformed_unknown_key.yaml`:** `sphere_mat4.yaml` with `sphere_radius` misspelled `sphere_raduis`. Expected `.field` = `geometry.sphere_raduis`.
  - **`malformed_wrong_type.yaml`:** `sphere_mat4.yaml` with `frequency_hz: 64.0e6`, which YAML reads as a `str`. Expected `.field` = `study.frequency_hz`.
  - **`malformed_missing_section.yaml`:** `f_small_four_port.yaml` with no `excitation:` section. Expected `.field` = `excitation`.
  - **`malformed_duplicate_key.yaml`:** `sphere_mat4.yaml` with `sigma` given twice (0.57, then 0.75). Expected `.field` = `material.sigma`. This is the loader's guard, **printed and asserted in the parse node, and not one of the row's three**.

  Each file has a one-line comment header citing this item. Each malformed file differs from its parent in exactly the one stated place.
- (9) **New `tests/package/test_wf1_study_yaml.py`.** `@complex_only` goes on the two solve nodes, so `WF-14`'s real-build collection imports the module cleanly.
  - **Imports from `test_wf15…` / `test_wf16…`:** helpers and constants only, namely `_study`, `GEOMETRY`, `_closed`, `_record`, `_rel` and the constants. **Never a `test_*` function or a fixture (`deg1`, `leg`)**: an imported test name is collected again in this module.
  - **Helper `_spec_diff(a, b) -> list[str]`.** It returns the dotted names that differ between two `Study` objects. It walks `geometry`, `material` and `excitation` dataclass fields plus `frequency_hz`, `degree` and `density_kg_m3`. It compares `type()` first, then floats by `float.hex`, complex numbers by the hex of their real and imaginary parts, tuples elementwise, and everything else by `==`.
  - **Nodes:**
    - `test_wf1_good_files_build_the_python_studies` (no solve; (K), (Deg), (c1) parse half, duplicate key);
    - `test_wf1_malformed_files_name_their_field` (no solve; (X));
    - `test_wf1_cli_sphere_reproduces_wf15` (W1);
    - `test_wf1_cli_f_small_reproduces_wf16` (W2).
  - **Width guard on both solve nodes, asserted and never skipped:** `MPI.COMM_WORLD.size == ANS4_EXACT_STRIP_RECORD_RANK_WIDTH` (2).
  - **Output paths.** Each CLI node first has rank 0 delete `/workspace/logs/wf1/<stem>/` if it exists, then `Barrier`; `/logs/` is gitignored, `.gitignore:82`. It then calls `cli.main(["run", <file>, "--output-dir", "/workspace/logs/wf1/<stem>"])` in-process and reads `summary.json` on every rank after `main` returns. A stale summary therefore cannot pass for a fresh one.
- (10) **New `tests/package/test_wf1_installed_cli.py`.** It uses stdlib only (`subprocess`, `sys`, `pathlib`), with **no import from `tests/` and none of `fem_em_solver`**. It is skipped unless `WF13_INSTALLED_CHECK=1`, and it locates `femem` at `Path(sys.prefix) / "bin" / "femem"`. Node: `test_wf1_installed_femem_refuses_the_malformed_files` (W5, (X′)).

**Anchors (asserted).**
- **(K) The YAML builds the Python `Study`, bitwise and field by field (no solve).**
  - `_spec_diff(load_study("sphere_mat4.yaml"), _study(degree=1)) == []`.
  - `_spec_diff(load_study("f_small_four_port.yaml"), <the WF-16 Python Study at degree=1>) == []`.
  - The preset route: the sphere mapping with `material: {preset: mid}`, built through the private mapping entry point that `load_study` uses (no extra file), gives `_spec_diff(..., _study(material="mid", degree=1)) == []`.
  - Print each list, and the `float.hex` of `geometry.resolution_sphere` / `resolution_far` / `study.frequency_hz` on both sides.
- **(Deg) Degree 2 when omitted, structurally.** Both good mappings with `study.degree` deleted give `.degree == 2`. Exact.
- **(X) Three malformed files fail with a named field (no solve).** For each malformed file:
  - `load_study` raises **`StudyConfigError`**, not a bare `KeyError`, `TypeError` or YAML error;
  - `exc.field ==` the expected dotted path above, and that path `in str(exc)`;
  - `cli.main(["run", <file>, "--output-dir", "/workspace/logs/wf1/malformed"])` returns **2** on every rank (allreduced `MPI.MIN` == `MPI.MAX` == 2), and the output directory **does not exist** afterwards, which shows the parse fails before any mesh or write.
  - The three `.field` values are pairwise distinct.
  - Duplicate key (printed and asserted, beside the three): `malformed_duplicate_key.yaml` raises `StudyConfigError` with `.field == "material.sigma"`.
- **(U) Upstream identity through the CLI's own parse, sphere (W1).**
  - DG0 σ / εᵣ / μᵣ come from `build_material_fields` / `build_mu_r_field`, built on the Python Study's solved `mesh` / `cell_tags` once with `load_study(file).material` as `phantom_material` and once with `_study(degree=1).material`, both with `phantom_tag=SPHERE_TAG` and the same background.
  - Each reads max |Δ| == 0.0 with 0 differing on owned dofs, using `_dg0_owned_diff` (`test_phantom_material_model.py`, as `WF-15` imports it).
- **(S) The CLI's outputs are the run's.** In both CLI nodes:
  - `rc == 0` on every rank;
  - `summary.json` exists, with `femem == 1`, `degree == 1`, `ranks == 2`, and `family` equal to the file's kind;
  - `cells ==` the Python-built Study's global cell count in the same process, exactly; the version-tagged 74 020 / 117 206 are printed and never asserted at equality.
- **(R) The records, through the CLI.**
  - **Sphere:** (R1) |SAR_yaml − `MAT4_DEG1_RECORD`| / record ≤ `L0_BAND`, where SAR_yaml is the summary's `mean_sar_w_per_kg`; (R2) `f"{SAR_yaml:.6e}" == MAT4_DEG1_PRINT`. Backed by `WF-15` W1's 6.222e-11.
  - **Birdcage:** (R1) ‖S_yaml − `P1["64"]`‖_F / ‖P1‖_F ≤ `LEG_D0_REPRODUCTION_BAND`, where S_yaml comes from the summary. Backed by 9.557e-11 (`WF-16` W1 and items 212 / 225).
- **(T) The Touchstone file (birdcage).**
  - (T1): `load_touchstone(<out>/f_small_four_port.s4p)` and `skrf.Network` against the summary's S give max |ΔS| ≤ `ANCHOR_TOL`; frequency within `FREQ_ATOL_HZ` of 64 MHz; `nports` 4; ports P1–P4; z0 50.
  - (T2): the loaded matrix against `P1["64"]` ≤ `LEG_D0_REPRODUCTION_BAND`.
  - These are `WF-16`'s (T1) / (T2) on the CLI's file.
- **(I) The YAML run equals the Python `Study`, in the same process.**
  - Sphere: |SAR_yaml − SAR_py| / SAR_py ≤ `L0_BAND`, where SAR_py comes from `_study(degree=1).solve(); .mean_sar()`.
  - Birdcage: ‖S_yaml − S_py‖ / ‖S_py‖ ≤ `LEG_D0_REPRODUCTION_BAND`.
  - Both are printed at `.3e`. Backed by `WF-15` W1 (I) 0.000e+00 and `WF-16` W1 (I) 9.110e-16, the same comparison in-process on the same fixtures. Not asserted bitwise (traps).
- **(G) The gates the U3 / U4 tests assert, on the CLI's outputs.**
  - Sphere: SAR_yaml against `_closed(SIGMA_HIGH, EPSILON_R_SPHERE)` ≤ `MAT4_BAND`.
  - Birdcage: `PORT-11` (i) reciprocity ≤ `RECIPROCITY_BAND`; (ii) σ_max(S) and every column power ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE`; (iii′) every Z-class spread ≤ `ADJACENT_SPREAD_BAND` and separation ≥ `POOLED_SEPARATION_FLOOR`. These are evaluated on the summary's S and Z with `_circulant_classes` / `_class_spread`, imported as `test_wf16…:73–80` imports them.
- **(X′) The installed `femem` (W5, real build, no `mpiexec`).**
  - `femem` exists under the venv's `bin/`, and its shebang names the venv's interpreter.
  - `femem check sphere_mat4.yaml` returns 0, and its stdout names `sphere_in_box` and degree 1.
  - `femem check control_sigma_ulp.yaml` returns 0.
  - For each of the three malformed files, `femem check` returns **2**, and its stderr contains the expected dotted field.

**Negative controls (rule (e)).**
- **(c1) One value changed, through the parser. Asserted by construction.**
  - The check: `load_study("control_sigma_ulp.yaml")` does **not** raise, and `_spec_diff(it, _study(degree=1)) == ["material.sigma"]` exactly. So the parser reads the field, (K)'s comparator resolves one ULP, and nothing else moved.
  - Backing: `float("0.5700000000000001") == np.nextafter(0.57, 1)` is an IEEE round-trip identity (host arithmetic 2026-10-05).
  - Downstream, in W1 on the solved mesh: (U)'s DG0 σ compare with this material reads exactly 2⁻⁵³ = 1.1102230246251565e-16 on > 0 owned dofs, **asserted**. The backing is the same comparison on the same fixture and code line: `MAT-1` (c1) 2⁻⁵³ on 19 496 dofs (`20260930T003229Z_MAT-1-reanchor-W1.log:491`) and `WF-15` W1 (c1), the same reading. It is the smallest nonzero reading, so it shows that (U)'s `== 0.0` is at full resolution.
- **(c2) Sphere, the record comparator against the wrong degree. Asserted.**
  - The check: |SAR_yaml / `WF15_DEG2_RECORD` − 1| > `L0_BAND`, with `WF15_DEG2_RECORD` = 8.0116577245e-08. The record is restated with its citation (`20261005T170537Z_WF-15-r2-W2a-deg2-study.log`, PROJECT_PLAN.md:2890), the way `WF-15` restates `MAT4_DEG1_RECORD`.
  - Expected reading: ≈ 3.494 % (= 8.291558128e-08 / 8.0116577245e-08 − 1, host arithmetic).
  - Backing: the same two numbers on the same fixture. W2a measured SAR_d2 / `MAT4_DEG1_RECORD` − 1 = −3.3757 %, which is ≈ 3.5e7 × the band.
- **(c3) Birdcage, cross-degree record swap. Asserted.**
  - The check: ‖S_yaml − `P2_X1["64"]`‖ / ‖P2_X1‖ > `LEG_D0_REPRODUCTION_BAND`.
  - Backing: `WF-16` W1 (c1) measured **0.0491**, the same comparator on the same fixture at the same degree and width (`20261005T202623Z_WF-16-W1-deg1.log`).
- **(c4) Touchstone P1↔P2 transposition. Asserted.**
  - The check: (T2)'s comparator on the loaded matrix with rows and columns 1 ↔ 2 swapped must exceed `LEG_D0_REPRODUCTION_BAND`.
  - Backing: `WF-16` W1 (c2) measured **0.2492** on the same comparison. The transpose control stays blind on a reciprocal S and is not reused (`WF-16`'s trap).
- **(c5) The malformed comparator is not a constant. Asserted by construction.**
  - The check: the two good files and the control file return 0 from `cli.main(["check", …])` and from the installed `femem check`, while the three malformed files return 2. A refuse-everything parser fails (c5), and an accept-everything one fails (X).
  - The three expected fields are distinct strings, so a parser that always names the same field fails (X).

**Printed (not asserted).**
- `yaml.__version__`.
- Per-leg seconds: the CLI run (parse, mesh, solve or sweep, write) and the Python-built run.
- Summed peak `ru_maxrss` via `report_peak_rss`, never cgroup `memory.peak`.
- The `sha256` of each YAML file.
- The summary's `cells` beside the version-tagged 74 020 / 117 206.

**Windows.** Foreground; the harness call is the whole Bash command, bare and repo-relative (`bash_guard` rule 3); Bash `timeout: 660000`; `-s`; node IDs, never `-k`; `tests/environment` first; complex windows `source /usr/local/bin/dolfinx-complex-mode` + `FEM_EM_REQUIRE_COMPLEX=1`; `PYTHONPATH=/workspace/src:/workspace`; no pipe or redirect. Every window is well under 10 min, so no durable capture is needed.
- **W1, parse nodes and the sphere through the CLI.** Standard, `-n 2`, complex.
  - Predicted 60–110 s: two degree-1 sphere solves including mesh, at ≈ 14.2 s each (`WF-15` W1); DG0 builds; environment and imports. `WF-15` W1, with three solves, (F)'s mesh and four operator assemblies, took 115 s.
  - `scripts/testing/run_and_log.sh WF-1-W1-parse-sphere "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/package/test_wf1_study_yaml.py::test_wf1_good_files_build_the_python_studies tests/package/test_wf1_study_yaml.py::test_wf1_malformed_files_name_their_field tests/package/test_wf1_study_yaml.py::test_wf1_cli_sphere_reproduces_wf15 -s -v --tb=short'"`
- **W2, the F-small four-port through the CLI.** Heavy by ceiling (400 s, `WF-16` W1's ceiling on the same fixture at `-n 2`), complex.
  - Predicted 110–170 s: two meshes at ≈ 26 s each and two 64 MHz sweeps at ≈ 9.7 s each. `WF-16` W1 had the same shape (Study plus direct call) and took 138 s.
  - `scripts/testing/run_and_log.sh WF-1-W2-f-small "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 400 mpiexec -n 2 python3 -m pytest tests/environment tests/package/test_wf1_study_yaml.py::test_wf1_cli_f_small_reproduces_wf16 -s -v --tb=short'"`
- **W3, `WF-14`'s module** (the envelope clause; rule (c) for the registry and snapshot). Real build, no `mpiexec`, ≈ 5 s.
  - `scripts/testing/run_and_log.sh WF-1-W3-wf14 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 python3 -m pytest tests/package/test_wf14_validation_envelope.py -v -s --tb=short'"`
- **W4, `WF-13` re-run (rule (c)).** `__init__.py` and `pyproject.toml` change what the installed package is. Run the command verbatim from line 12 of `20261005T171823Z_WF-15-r2-W4-wf13.log`, under tag `WF-1-W4-wf13`. It took 37 s.
- **W5, the installed console script.** Smoke by prediction (≈ 30 s: a 9 s setup plus five `femem check` subprocesses at ≈ 3 s each); real build. It re-runs `wf13_setup.sh real`, which is unchanged, and copies only this item's stdlib-only module and YAML files into the run directory.
  - `scripts/testing/run_and_log.sh WF-1-W5-installed-femem "docker compose exec -T fem-em-solver bash -lc 'timeout -k 10 50 bash /workspace/scripts/package/wf13_setup.sh real && mkdir -p /tmp/wf13/run/tests/package/wf1 && cp /workspace/tests/package/test_wf1_installed_cli.py /tmp/wf13/run/tests/package/ && cp /workspace/tests/package/wf1/*.yaml /tmp/wf13/run/tests/package/wf1/ && cd /tmp/wf13/run && env -u PYTHONPATH WF13_INSTALLED_CHECK=1 timeout -k 10 90 /tmp/wf13/venv/bin/python -m pytest -p no:cacheprovider tests/package/test_wf1_installed_cli.py -v -s --tb=short'"`
- **No orphan sweep is needed.** There is no degree-2 compile. The degree-1 sphere and four-port forms are warm from today's `WF-15` / `WF-16` windows. After any kill, run `pgrep -c python3` in the container first.

**Traps already paid for.**
- **No solve output is compared bitwise.** `MAT-1` item 81 went red at 1e-12 on an identical re-run (known-issues.md:81). `MAT-5` item 237's (R0c) read 2.747e-12 in one process, and `WF-16`'s in-process (I) read 9.110e-16. Bitwise belongs upstream of the solve only: (K), (U), (c1), and the JSON round trip.
- **PyYAML (YAML 1.1) typing** (measured above):
  - `64.0e6` and `1e-9` are strings, and a strict loader must say why (the hint in the message);
  - `true` is a `bool`, and `isinstance(True, int)` is `True`, so a number check must exclude `bool` first;
  - duplicate keys keep the last value silently, unless the loader refuses them;
  - YAML lists are `list`, so `port_box_size` must become a `tuple`, or (K) reads a type difference.
- **`yaml.safe_load` / `SafeLoader` only.** `import yaml` stays inside the loader function, so `WF-14`'s real-build collection and `WF-13`'s venv do not depend on it at import time.
- **`main` returns, never exits.** An in-process `sys.exit` ends the pytest run. `argparse`'s own errors raise `SystemExit`, so pass only valid argv in-process.
- **Rank-safety.**
  - The summary and Touchstone are written on rank 0 after every reduction, then `Barrier`. There is no `comm.` call under `if comm.rank == 0:` (carried ruling, 10-04).
  - `cell_tags.values` is rank-local; DG0 compares use the owned slice (`_dg0_owned_diff`); the CLI's return code is allreduced before it is asserted.
- **A stale `summary.json` must not pass for a fresh one.** Delete the output directory at node start (rank 0, `Barrier`), as `WF-15`'s hand-off did.
- **A record holds at its own width** (`OPS-41`). Both CLI nodes assert `-n 2`, and the width guard never skips.
- **`TimeHarmonicSolver` and `run_n_port_sparameter_sweep` default to degree 1, and the Study to 2.** The good files set `degree: 1` explicitly, and (Deg) checks the omitted case. A loader that injects its own default instead of deferring to `Study`'s fails (Deg) if that default ever drifts. Say so in a comment.
- **`femem` is not on `PATH` in the dev container** (the image's editable install predates the entry point). Never `pip install` into the live container; W5's fresh venv is the console script's only test.
- **The `src/` package never imports `tests/`, and the installed module imports neither `tests/` nor `fem_em_solver`.** W4 / W5 catch it.
- **Importing a `test_*` name or a fixture from `test_wf15…` / `test_wf16…`** re-collects it in this module (a duplicated solve, or a fixture leak). Import helpers and constants only.
- **`WF-14` (N): node IDs only.** The snapshot is sorted. Adding to `__all__` without both edits makes W3 red by design.
- **`complex_only` on the two solve nodes,** or `WF-14`'s real-build collection errors out (rc 2).
- **Exact rectangles, κ off** (binding ruling 2026-10-02). The YAML offers no κ key: `width_correction_kappa` is an unknown key, which (X)'s mechanism covers.
- **No pipe, redirect or `; date` in the harness command** (the footer records the pipe's status, and `bash_guard` denies the rest). Status 124 is a finding, not a retry.
- **`-k a or b` splits inside the quoted container command,** so use node IDs.
- **Rule (i):** no edit to `study_yaml.py`, `cli.py`, the test modules, the YAML files, the registry or the snapshot after their last green window.
- **Write ≈, never `~`, in tracked Markdown.** Quote no AED number; this item touches no AED file.

**Scope.**
- **Covered:**
  - YAML schema v0 (`femem: 1`) for the two `Study` families that exist, at degree 1, `-n 2`, 64 MHz;
  - the lossy sphere with an imposed field and mean SAR;
  - the F-small four-port on exact rectangles, κ off, with S, Z and Touchstone;
  - the default degree 2 when omitted, structurally;
  - three named-field refusals, plus the duplicate-key guard;
  - the installed console script's parse and refusals.
- **Not claimed:**
  - **any degree-2 solve through the YAML**: `WF-16`'s degree-2 record holds at `-n 8` and `WF-15`'s degree-2 solve costs 282 s at `-n 2`. That leg is left to the review (open question 2);
  - `femem run` through the installed script with a solve;
  - `-n 8`; 10 / 128 MHz; multi-frequency sweeps;
  - a stable or semver'd schema; unit strings; YAML for any fixture beyond these two;
  - `WF-4`'s presets; docs or a tutorial (`WF-17`);
  - the `__main__.py` route (unexercised, disclosed).

**Status it can move.** `WF-1` 🧪 → ✅ on this Done-when, which the review writes into the row (rule (j)):

> **Done-when (§4):** `femem run <file>.yaml` (console script `fem_em_solver.cli:main`) drives `WF-15`'s lossy sphere and `WF-16`'s F-small four-port, each at degree 1 and `-n 2`, to bit-identical outputs *(read by the 2026-10-05 18:00 review, rule (j); the weekly ratifies: the YAML-built `Study` equals the Python-built one field by field, floats by `float.hex`, degree 2 when omitted; the solve outputs — SAR, S and the Touchstone file — equal `MAT-4`'s / item 192's records and the Python `Study`'s within `L0_BAND` / `LEG_D0_REPRODUCTION_BAND`)*; three malformed files (unknown key, wrong type, missing section) raise `StudyConfigError` naming the field, in-process and through the installed `femem`, while a one-value variant parses and changes only that field; entry points registered in `WF-14`'s envelope; elapsed recorded ⇒ ✅.

The row's tier column reads "smoke", which predates the step. The review should amend it to "heavy by ceiling (W2 400 s)". The audit maps the Done-when onto the anchors:
- YAML-built equals Python-built: (K) and (U);
- degree 2 when omitted: (Deg);
- the solve outputs within the bands: (R), (T), (I) and (G);
- three malformed files: (X), and (X′) for the installed `femem`;
- the one-value variant: (c1);
- the envelope: W3;
- elapsed: every window's footer.

**Negative result.** In every case: report the measurement, park on `attempt/WF-1-*`, mark item 256 BLOCKED in the same commit (rule (d)) and journal it in attempts.md.
- **(K), (Deg) or (U) red:** the YAML builds a different `Study`. This is a loader bug in the new code: fix it, and re-run the window in-slot only if the fix is confined to the new files. A YAML literal that is one ULP off is fixed in the file, which is in-letter.
- **(X) red:** a malformed file is accepted, raises the wrong type, or names the wrong field. The fix is in-letter loader work. Re-run W1.
- **(S) or (T1) red:** a writer or layout defect in the CLI's outputs. Stop if the defect is in `export_touchstone` itself (that is `PORT-7`'s module's business).
- **(R) or (T2) red with (I) and (K) green:** the environment has drifted since today's `WF-15` / `WF-16` windows. This is a finding for the review / `record-reconciler`. No band moves.
- **(I) red with (K) green:** the in-process solve floor exceeds the band through a second `Study` instance. That is a finding for known-issues.md:80–81. Stop.
- **(G) red:** a real finding, because the same S / SAR is green through `WF-15` / `WF-16`. Stop.
- **(c1)–(c5) at or under their floor:** the comparator is blind. Stop.
- **W3 red:** a registry or snapshot omission. Fix it (that is in-letter) and re-run W3.
- **W4 red:** the export or the `[project.scripts]` table broke the installed package. Revert (4) / (5); BLOCKED.
- **W5 red with W4 green:** the console script is not generated or not runnable in the venv. Report the `pip install` lines from the log. `WF-1` holds 🟡, with the in-process half green.
- **A window killed at its ceiling:** BLOCKED with the reading. Never lengthen the timeout and never change `-n`, because the widths are the records'. After a kill, run `pgrep -c python3` first.

No `PREDICTED RED:` line. A red anywhere stops the slot.

#### Execution record — item 256 (2026-10-06, 19:30 CDT slot)

**Landed as briefed.** `src/fem_em_solver/study_yaml.py` (`StudyConfigError`, `load_study`, `run_study_file`, the private `_study_from_mapping`, the non-exported `StudyRun`), `src/fem_em_solver/cli.py` (`main`, `run` / `check`), `src/fem_em_solver/__main__.py`, `[project.scripts] femem = "fem_em_solver.cli:main"`, three `__all__` names (30 → 33), the snapshot and registry entries, seven YAML files under `tests/package/wf1/`, and two test modules. `study.py` is untouched.

**Two implementation choices to disclose.**
- **Duplicate-key guard.** The loader is a `yaml.SafeLoader` subclass whose `construct_mapping` refuses a repeated key. To name the field by its dotted path, it composes the node tree first, then constructs it, and maps the offending mapping node back to its path. It reads `material.sigma` on the guard file.
- **Order inside `geometry` / `excitation`.** `kind` is checked first (missing, then unknown value), because the allowed key set depends on it. Unknown keys, missing keys and types follow in the briefed order. A `geometry.kind` / `excitation.kind` mismatch (a sphere driven by ports, or the reverse) is refused as `excitation.kind`. None of the three gated files reaches either branch.

**Windows (all Status 0, first run, `-n 2` for the solves).**

| Window | Log | Elapsed | Result |
|---|---|---|---|
| W1 parse + sphere (standard, 180 s) | `20261006T003504Z_WF-1-W1-parse-sphere.log` | 100 s | 15 passed |
| W2 F-small (heavy by ceiling, 400 s) | `20261006T003650Z_WF-1-W2-f-small.log` | 142 s | 13 passed |
| W3 `WF-14` envelope | `20261006T003919Z_WF-1-W3-wf14.log` | 5 s | 6 passed |
| W4 `WF-13` re-run (verbatim) | `20261006T003931Z_WF-1-W4-wf13.log` | 37 s | 3 + 3 passed |
| W5 installed `femem` | `20261006T004014Z_WF-1-W5-installed-femem.log` | 15 s | 1 passed |

**Measured.**
- Environment: `yaml.__version__` 6.0.3 in the container.
- (K): spec diff `[]` for the sphere, F-small and "mid" preset routes. `float.hex` is equal on both sides: resolution_sphere `0x1.0624dd2f1a9fcp-10`, resolution_far `0x1.0624dd2f1a9fcp-9`, frequency `0x1.e848000000000p+25`, σ `0x1.23d70a3d70a3dp-1`.
- (Deg): `[2, 2]`.
- (X): the three files name `geometry.sphere_raduis`, `study.frequency_hz` and `excitation`. `run` and `check` both return (2, 2) (min, max over ranks), and no output directory is created. The duplicate-key file names `material.sigma`.
- (c5): `check` returns 0 on the two good files and on the control file.
- Sphere (W1):
  - CLI run 15.0 s, Python Study 14.9 s, 74 020 cells on both;
  - (U) DG0 σ / εᵣ / μᵣ (0.0, 0, 74 020);
  - (c1) 1.1102230246251565e-16 = 2⁻⁵³ on 19 496 owned dofs;
  - (R1) 6.159e-11, and (R2) `8.291558e-08` == `8.291558e-08`;
  - (I) 8.244e-13;
  - (G) 3.536 % against the 10 % `MAT4_BAND`;
  - (c2) 3.4937 % (3.49e+07 × `L0_BAND`).
- F-small (W2):
  - CLI run 36.8 s, Python Study 37.2 s (mesh 28.3 s, sweep 8.5 s), 117 206 cells on both;
  - (R1) 9.557e-11, and (c3) 0.0491 against `P2_X1`;
  - (T1) 6.873e-14 for both readers; frequency 64000000.0; ports P1–P4;
  - (T2) 9.557e-11, and (c4) P1↔P2 swap 0.2492;
  - (I) 9.110e-16;
  - (G) reciprocity 4.509e-16, σ_max 0.999760093664, max column power 0.820520164533, spreads 0.0747 / 0.0691 / 0.0897 %, separation 451.71×.
- (X′) (W5): the shebang is `#!/tmp/wf13/venv/bin/python`. `femem check` returns 0 on the good and control files, with `family: sphere_in_box` / `degree: 1` on stdout. It returns 2 on each malformed file, with the dotted field on stderr. Each call takes ≈ 0.8 s.

**Readings against the backing.** (R1) moved from `WF-15` W1's 6.222e-11 to 6.159e-11, and (I) read 8.244e-13 where `WF-15` (I) read 0.000e+00. Both are inside the `-n 2` cross-process MUMPS floor (`SOLVE_FLOOR_RECORD`, known-issues `MAT-1` entry). The YAML run and the Python run are two `Study` instances, each solved once, so this is the expected class of difference. The F-small readings repeat `WF-16` W1's digits exactly.

**Not claimed (unchanged from Scope):** a degree-2 solve through YAML; `femem run` with a solve through the installed script; `-n 8`; 10 / 128 MHz; the `__main__.py` route (unexercised); a stable schema. *(2026-10-07: the `__main__.py` route is retired by reference — `WF-13` step 3, §9 item 294, launches every `femem run --ranks N` child through `python -m fem_em_solver` and gates its SAR; `20261007T050705Z_WF-13-s3-W1-launch.log`.)*

### §9 item 293 — 2026-10-06 (chain U18, provisional)

**Outcome: ✅ DONE 2026-10-07 05:04Z** (00:00 CDT slot, attempt 1). Three windows green on first run, no in-slot fix. (P): the venv imports `fem_em_solver` from `/tmp/wf13/venv/lib/python3.12/site-packages`, shebang `#!/tmp/wf13/venv/bin/python`, summaries' `source` under `/tmp/wf13/run/u18/` with sha256 equal to the repo files', `ranks` 2, `degree` 1. Sphere (W1, `20261007T050136Z_WF-1-s2-W1-installed-sphere.log`, 41 s; summary elapsed 19.9 s): SAR 8.291558128489977e-08, 5.909e-11 from `MAT4_DEG1_RECORD`, `.6e` string equal, closed form 3.536 %, 74 020 cells, (c2) 3.4937 %. F-small (W2, `20261007T050226Z_WF-1-s2-W2-installed-f-small.log`, 65 s; summary elapsed 47.1 s): S 9.557e-11 from `P1`, (c3) 0.0491; `.s4p` 6.855e-14 from the summary for `load_touchstone` and scikit-rf, 64000000.0 Hz, P1–P4, file-vs-`P1` 9.556e-11, (c4) 0.2492; reciprocity 4.377e-16, σ_max 0.999760093664, spreads 0.0747 / 0.0691 / 0.0897 %, separation 451.71×; 117 206 cells. W3 `WF-14` (`20261007T050336Z_WF-1-s2-W3-wf14.log`, 5 s): 6 passed, 87 node refs, 0 violations.

*Drafted by a general-purpose agent for the 2026-10-06 daily review; the review reviews it.* Licence: §9 rule (6), a provisional step of §10 chain U toward (9a) / (9c) and the Phase 9 exit gate's clause "importing the *installed* package … reproduces a gated result within its band" (PROJECT_PLAN.md ≈ 4287). **No band, tolerance, Done-when or validation target is set or moved; every band and record is imported by name.** Overlaps none of U12–U17 or `EX-92` / `EX-94` / `EX-95`.

**`WF-1` step 2: `femem run` with a solve, through the installed console script, outside the repository.** `WF-1` step 1 drove the solves in-process from `/workspace`; its row lists "`femem run` with a solve through the installed script" as **not claimed**, and `WF-1` W5 exercised the installed `femem` on `check` only. This step runs the installed `/tmp/wf13/venv/bin/femem run` under `mpiexec -n 2`, cwd outside `/workspace`, no `PYTHONPATH`, on the two committed `WF-1` YAML files, then asserts the written `summary.json` and `.s4p` against the records `WF-1` step 1 already asserts.
- **Executor:** implementer. **Adds:** one test module `tests/package/test_wf1_installed_solve.py` (two nodes, below); `tests/package/validation_envelope.toml` gains the two nodes under `main` / `run_study_file` (adding lines only). **No `src/` change, no new YAML file** (the two committed files are copied, sha256-checked). `public_surface.txt` unchanged.
- **Tier / ranks:** W1, W2 standard by prediction, heavy by ceiling (inner timeouts sum to 275 s / 355 s); `-n 2` for both solves (the width of both records); complex build. **Cost:** class median ≈ 9 measured min; **≈ 35 predicted**: 15 fixed + ≈ 12 code + W1 ≈ 1.5 + W2 ≈ 2 + W3 ≈ 0.1 + ≈ 4 spare. Priced from: `WF-1` W5 (setup + installed `femem`) 15 s (`20261006T004014Z_WF-1-W5-installed-femem.log`), `WF-1` W1 CLI sphere run 15.0 s at `-n 2` (`20261006T003504Z_WF-1-W1-parse-sphere.log:281`), `WF-1` W2 CLI F-small run 36.8 s (mesh 28.3 s), `WF-13` W1 installed complex `mpiexec -n 2` route 35 s whole window.
- **Depends on:** nothing open — `WF-13` ✅ (`scripts/package/wf13_setup.sh` on main), `WF-1` ✅ (console script `femem`, the two YAML files).

**Anchors (asserted).** Records and bands imported, none restated: `MAT4_DEG1_RECORD`, `MAT4_DEG1_PRINT`, `WF15_DEG2_RECORD` (`tests/package/test_wf1_study_yaml.py:99`), `L0_BAND` (`tests/validation/test_mat5_em_coupled_sphere.py:104`), `MAT4_BAND`; `ANS4_EXACT_STRIP_S_P1_RECORDS["64"]` / `ANS4_EXACT_STRIP_S_P2_X1["64"]` through `test_wf16_study_four_port._record` / `_rel`, `LEG_D0_REPRODUCTION_BAND` (`tests/validation/test_port_birdcage_four_port.py:165`), `ANCHOR_TOL` / `FREQ_ATOL_HZ` (`tests/io/test_touchstone_skrf_crosscheck.py:51–52`), and the `PORT-11` bands exactly as `test_wf1_study_yaml.py:74–82` imports them.
- **(P) Provenance:** `env -u PYTHONPATH /tmp/wf13/venv/bin/python -c 'import fem_em_solver; print(fem_em_solver.__file__)'` run by `subprocess` with cwd `/tmp` prints a path under `/tmp/wf13/venv` and not under `/workspace/src` (`WF-13`'s (c2) reading); `/tmp/wf13/venv/bin/femem`'s first line is `#!/tmp/wf13/venv/bin/python`; each `summary.json` has `source` under `/tmp/wf13/run/u18/`, `sha256` equal to the repo file's sha256, `ranks == 2`, `degree == 1`.
- **(R-S) Sphere** (`test_wf1_installed_sphere_reproduces_mat4`): `summary["mean_sar"]["mean_sar_w_per_kg"]` (the key `WF-1`'s test reads, `test_wf1_study_yaml.py:324`) within `L0_BAND` of `MAT4_DEG1_RECORD`, its `.6e` string `== MAT4_DEG1_PRINT`, and within `MAT4_BAND` of the closed form; cells printed (74 020 in `WF-1` W1).
- **(R-F) F-small** (`test_wf1_installed_f_small_reproduces_p1`): `summary["s_matrix"]` (pairs → complex) within `LEG_D0_REPRODUCTION_BAND` of `P1` at 64 MHz; `cells == 117 206` printed against `ANS4_EXACT_STRIP_CELL_RECORD`; the installed-written `f_small_four_port.s4p` read back by `load_touchstone` and `scikit-rf` equals the summary S within `ANCHOR_TOL`, frequency within `FREQ_ATOL_HZ`, ports P1–P4; `PORT-11`'s reciprocity / passivity / class-spread / separation gates on the summary S, imported helpers, green.

**Negative controls (asserted, same fixture, prior measurement).**
- **(c2) degree confusion, sphere:** `|SAR / WF15_DEG2_RECORD − 1| > L0_BAND` — `WF-1` W1 read 3.4937 % (`20261006T003504Z_WF-1-W1-parse-sphere.log:288`); ceiling n/a (separation 3.49e+07 × band).
- **(c3) degree confusion, F-small:** `‖S − P2_X1‖/‖P2_X1‖ > LEG_D0_REPRODUCTION_BAND` — 0.0491 (`20261006T003650Z_WF-1-W2-f-small.log:3712`). **(c4) port swap:** P1↔P2 swapped S misses `P1` — 0.2492 (`:3714`).

**Windows** (foreground, Bash `timeout: 660000`; whole bare command; `-s`; node IDs, never `-k`; no pipes; each window re-runs `wf13_setup.sh real`, which wipes `/tmp/wf13`, so each reads only its own outputs).
- **W1, installed sphere.**
  `scripts/testing/run_and_log.sh WF-1-s2-W1-installed-sphere "docker compose exec -T fem-em-solver bash -lc 'timeout -k 10 50 bash /workspace/scripts/package/wf13_setup.sh real && source /usr/local/bin/dolfinx-complex-mode && timeout -k 10 15 bash /workspace/scripts/package/wf13_setup.sh complex && mkdir -p /tmp/wf13/run/u18 && cp /workspace/tests/package/wf1/sphere_mat4.yaml /tmp/wf13/run/u18/ && cd /tmp/wf13/run/u18 && env -u PYTHONPATH FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 120 mpiexec -n 2 /tmp/wf13/venv/bin/femem run sphere_mat4.yaml --output-dir /tmp/wf13/run/u18/out/sphere && cd /workspace && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 WF1_INSTALLED_OUT=/tmp/wf13/run/u18/out timeout -k 30 90 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_installed_solve.py::test_wf1_installed_sphere_reproduces_mat4 -s -v --tb=short'"`
- **W2, installed F-small.**
  `scripts/testing/run_and_log.sh WF-1-s2-W2-installed-f-small "docker compose exec -T fem-em-solver bash -lc 'timeout -k 10 50 bash /workspace/scripts/package/wf13_setup.sh real && source /usr/local/bin/dolfinx-complex-mode && timeout -k 10 15 bash /workspace/scripts/package/wf13_setup.sh complex && mkdir -p /tmp/wf13/run/u18 && cp /workspace/tests/package/wf1/f_small_four_port.yaml /tmp/wf13/run/u18/ && cd /tmp/wf13/run/u18 && env -u PYTHONPATH FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 200 mpiexec -n 2 /tmp/wf13/venv/bin/femem run f_small_four_port.yaml --output-dir /tmp/wf13/run/u18/out/f_small && cd /workspace && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 WF1_INSTALLED_OUT=/tmp/wf13/run/u18/out timeout -k 30 90 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_installed_solve.py::test_wf1_installed_f_small_reproduces_p1 -s -v --tb=short'"`
- **W3, `WF-14` registry + snapshot** (the envelope gained nodes):
  `scripts/testing/run_and_log.sh WF-1-s2-W3-wf14 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 python3 -m pytest tests/package/test_wf14_validation_envelope.py -v -s --tb=short'"`

**Traps already paid for.**
- The test **fails, never skips**, when `WF1_INSTALLED_OUT` or a summary / `.s4p` is missing — a skipped node reads green in the footer.
- `env -u PYTHONPATH` must sit on the `mpiexec` line itself; the complex-mode `source` re-sets the environment, which is why `wf13_setup.sh complex` regenerates the stack `.pth` after it (`WF-13` item 220).
- YAML float resolver (`WF-1` step 1): the copied files are byte-identical (sha256 asserted); never retype a literal.
- `-n 2` MUMPS is not bitwise across processes (known-issues `MAT-1` entry, `SOLVE_FLOOR_RECORD`): compare at the imported bands, never `==`, except the JSON round trip (`WF-1` step 1: exact by `repr`).
- Complex build: `source /usr/local/bin/dolfinx-complex-mode` + `FEM_EM_REQUIRE_COMPLEX=1`; `tests/environment/test_complex_mode.py` only (the full `tests/environment` migration scan cost `GEO-37` step 2 a 30 s ceiling, item 267). A killed window may leave an orphan FFCx `.c`: sweep per the rubric's one-liner only when `pgrep -c python3` reads 0.
- No AED number anywhere. Write ≈, never `~`.

**Scope.** Covered: the installed console script solves both `WF-1` fixtures at degree 1, `-n 2`, outside the repository, to their records. **Not claimed:** a STEP geometry (U12 / U13), degree 2 through YAML, other widths, a published image or PyPI, the exit gate itself (the operator's run, from the docs, on a STEP).

**What a negative result means.** Report and stop; park on `attempt/WF-1-step2-*`; mark item 293 BLOCKED; journal. (P) red ⇒ the venv imports the repo tree (packaging defect, a known-issues entry). (R) red with (P) green ⇒ the installed package diverges from the in-repo one (missing package data, a stale build) — a §7 annotation with the reading; never widen a band. A control not failing ⇒ the instrument is blind; stop.

No `PREDICTED RED:` line.

**Status it can move.** No glyph: `WF-1` stays ✅. On green the row's "Not claimed: … `femem run` with a solve through the installed script" is retired by a **step-2 record** sentence (readings, elapsed, logs); §10 U18 marked DONE; `WF-14`'s envelope gains two nodes; the Phase 9 exit-gate rehearsal on built-in geometry is recorded (the STEP half stays U13's).

### §9 item 295 — 2026-10-06 (chain U20, provisional)

*Drafted by a general-purpose agent for the 2026-10-06 daily review; the review reviews it.* Licence: §9 rule (6), a provisional step of §10 chain U toward (9a) "`femem run study.yaml` drives the same API" (PROJECT_PLAN.md ≈ 4287; results include "temperature via `thermal/`"). **No band, tolerance, Done-when or validation target is set or moved.** Not `EX-94` (a Python example of `Study.temperature`); this is the YAML entry point, and touches no example.

**`WF-1` step 3: a `thermal:` section in the YAML schema, driving `Study.temperature(RobinCooling)`.** `WF-15` step 2 (item 260) made the steady Robin-cooled temperature a `Study` result; schema v0 cannot ask for it. This step adds one optional top-level section and reproduces `MAT-10`'s recorded (A) / (P) digits and `WF-15` step 2's (IT) identity through `femem run`.
- **Executor:** implementer. **Adds (additive):** `src/fem_em_solver/study_yaml.py` — `thermal` in `_TOP_KEYS` (optional); `_parse_thermal(sec, geometry_kind, source)` → `RobinCooling` with keys `k_w_per_m_k`, `h_w_per_m2_k`, `t_inf_c` (required), `source_scale`, `degree` (optional; omitted ⇒ `RobinCooling`'s own defaults, never a loader default — `WF-1`'s (Deg) rule); validation order fixed after `excitation`, before `outputs`; a `four_port_birdcage` geometry with `thermal` raises `StudyConfigError` at field `thermal` (mirrors `Study.temperature`'s `TypeError`); a non-positive value raises at its dotted field **before** `RobinCooling` is constructed (no bare `ValueError` escapes). `run_study_file`: after `mean_sar`, when a cooling is present, `study.temperature(cooling)`; `summary["temperature"] = {"source_power_w", "convected_power_w", "theta_max_k"}` (θ max over owned dofs of Re(T − t_inf), `MPI.MAX`-reduced, before rank 0 writes); `StudyRun` gains `thermal=None` (non-exported). Keep `load_study`'s return a `Study`; carry the cooling inside the private parse tuple and update every caller of `_parse_mapping` in the same diff. Schema stays `femem: 1` (additive key; `WF-1` claims no stable schema). Two YAML files: `tests/package/wf1/sphere_mat10_thermal.yaml` (= `sphere_mat4.yaml` + `thermal: {k_w_per_m_k: 0.6, h_w_per_m2_k: 100.0, t_inf_c: 20.0, source_scale: 25000000.0}`) and `…_control_2g.yaml` (`source_scale: 50000000.0`); two malformed files (`thermal` on the F-small file; `thermal.k_w_per_m_k: 0.0`). New `tests/package/test_wf1_yaml_thermal.py`; envelope nodes under `run_study_file` / `load_study`.
- **Tier / ranks:** W1 standard, `-n 2`, complex. **Cost:** class median ≈ 9 measured min; **≈ 50 predicted**: 15 fixed + ≈ 22 code + W1 ≈ 2 + W2 ≈ 1.7 + W3 ≈ 2.4 + W4 ≈ 0.1 + W5 ≈ 0.3 + ≈ 6 spare. Priced from `WF-15` step 2 W1 on this fixture: degree-1 EM solve + mesh 13.7 s, Study thermal 0.4 s, 2.39 GiB summed, window 96 s (`20261006T035826Z_WF-15-s2-W1-deg1.log:244`, `:251`, `:317`); two EM solves here (YAML route + Python comparand).
- **Depends on:** nothing open — `WF-15` step 2 (`Study.temperature`, `RobinCooling` exported) and `WF-1` step 1 on main. Shares `validation_envelope.toml` with items 293 / 294: edit, adding lines only; on a rebase keep both sides.

**Anchors (asserted).** Imported, none restated: `K`, `H`, `T_INF`, `G`, `A_BAND`, `P_BAND`, `B_BAND`, `C_BAND` and the helpers `_interior_indicator`, `_interior_rel`, `_theta_cf` (`tests/validation/test_mat10_sar_to_heat_sphere.py:52–62`), `MAT10_A_PRINT` / `MAT10_P_PRINT` / `MAT10_SUB_CELLS` and `_theta_owned_rel` (`tests/package/test_wf15_study_temperature.py:81–83`, `:100`), `L0_BAND`, `_study`.
- **(K) bitwise up to the solve:** the YAML-built `RobinCooling` equals `RobinCooling(K, H, T_INF, source_scale=G)` field by field by `float.hex`; the `Study` spec equals `_study(degree=1)`'s (`WF-1`'s `_spec_diff` reads `[]`); omitted `degree` / `source_scale` give `RobinCooling`'s defaults (1, 1.0).
- **(X) malformed, no solve:** the two files raise `StudyConfigError` naming `thermal` and `thermal.k_w_per_m_k`, in-process, rc 2 on both ranks through `cli.main`, nothing written.
- **(IT) YAML vs Python:** the `run_study_file` θ against `Study.temperature(RobinCooling(K, H, T_INF, source_scale=G))` on an independently solved `_study(degree=1)`, owned max rel ≤ `L0_BAND` (`WF-15` (IT) read 5.603e-14 on one Study, `:247`; two Study instances read 8.244e-13 on SAR in `WF-1` W1 `:286`).
- **(R) records on the YAML θ:** (A) and (P) at `.4e` `== MAT10_A_PRINT` / `MAT10_P_PRINT`; sub cells `== MAT10_SUB_CELLS`; (B) ≤ `B_BAND`, (C) ≤ `C_BAND`; `summary["temperature"]` equals the in-memory values exactly (JSON round trip).

**Negative control (asserted, same fixture, prior measurement).** **(c1) missing ½:** the control file's parsed `RobinCooling` (2G), applied to the already-solved YAML Study (no second EM solve; (K) makes the EM inputs identical), gives an (A) comparator `> A_BAND` — `WF-15` step 2 read 1.0649e+00 on this fixture (`20261006T035826Z_WF-15-s2-W1-deg1.log:250`); ceiling n/a (9.7 × band).

**Windows** (foreground, Bash `timeout: 660000`; bare whole command; `-s`; node IDs, never `-k`; no pipes).
- **W1, thermal through YAML.**
  `scripts/testing/run_and_log.sh WF-1-s3-W1-thermal "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_thermal.py -s -v --tb=short'"`
- **W2 / W3, `WF-1` step-1 re-gate (rule (c): `study_yaml.py` is its module).** Item 256's W1 and W2 commands verbatim (this file, `WF-1-W1-parse-sphere` and `WF-1-W2-f-small` lines of item 256's Windows list) under tags `WF-1-s3-W2-wf1-w1` / `WF-1-s3-W3-wf1-w2`; 100 s / 142 s measured. W3 may be skipped only if the review rules the port branch untouched (disclose).
- **W4, `WF-14`:** `scripts/testing/run_and_log.sh WF-1-s3-W4-wf14 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 python3 -m pytest tests/package/test_wf14_validation_envelope.py -v -s --tb=short'"`
- **W5, installed `femem check`:** item 256's W5 command verbatim under tag `WF-1-s3-W5-installed-femem` (15 s measured).

**Traps already paid for.**
- **PyYAML 1.1 floats need a dot and a signed exponent** (`WF-1` step 1): `25000000.0` / `50000000.0` / `0.6` / `100.0` / `20.0`, never `2.5e7`; and the strict loader refuses duplicate keys.
- `ufl.max_value` / ordering comparisons do not compile on complex operands (`OPS-22`): θ max is taken in NumPy on owned dofs, reduced, never in UFL.
- Rank-safety: `theta_max_k` reduced with `MPI.MAX` before rank 0 writes; owned dofs only (`size_local`).
- `Study.temperature` reuses `self.result` (one EM solve); (c1) must not re-solve.
- `src/` never imports `tests/`; `WF-14` (N): node IDs only. A killed compile leaves an orphan `.c` — sweep per the rubric only when `pgrep -c python3` reads 0. Write ≈, never `~`.

**Scope.** Covered: the steady Robin-cooled temperature of the sphere through YAML at degree 1, `-n 2`. **Not claimed:** degree 2 through YAML, transient / Pennes / duty cycle, thermal on a port Study, a stable schema.

**What a negative result means.** Report and stop; park on `attempt/WF-1-step3-*`; item 295 BLOCKED; journal. (K) red ⇒ a parser bug (fix in-slot if confined to new code). (IT) / (R) red with (K) green ⇒ the YAML route assembles a different thermal problem — a §7 annotation, never a wider band. W2 / W3 red ⇒ revert `study_yaml.py`, BLOCKED. (c1) not failing ⇒ blind instrument; stop.

No `PREDICTED RED:` line.

**Status it can move.** No glyph (`WF-1` ✅): the row gains a **step-3 record** (temperature through YAML; readings, elapsed); §10 U20 DONE; `WF-14`'s envelope gains nodes; the (9a) clause "`femem run study.yaml` drives the same API" covers the thermal result.

**Outcome: ✅ DONE 2026-10-07 05:21Z** (00:00 CDT slot, attempt 1, third item). Five windows green on first run, no in-slot fix. `study_yaml.py` additive: `thermal` in `_TOP_KEYS`, `_parse_thermal` (validated after `excitation`, before `outputs`; non-positive `k` / `h` / `source_scale` and `degree < 1` raise at the dotted field before `RobinCooling` is built; `four_port_birdcage` raises at `thermal`), `_parse_mapping` returns `(Study, kind, outputs, cooling)` (callers `_study_from_mapping`, `_load`, `load_study`, `run_study_file` updated), `StudyRun.thermal`. W1 (`20261007T051444Z_WF-1-s3-W1-thermal.log`, 39 s, 6 passed): (K) cooling diff `[]`, spec diff `[]`, control diff `['source_scale']`, omitted keys give 1.0 / 1; (X) fields `thermal` / `thermal.k_w_per_m_k`, run / check rc (2, 2), nothing written; `run_study_file` 17.2 s; `summary["temperature"]` source 0.008651525642063975 W, convected 0.008651525642062842 W, θ max 0.12682242431550605 K, JSON equal; (R) sub cells 19 496, (A) 3.2463961444e-02 [3.2464e-02], (P) 3.1923343273e-03 [3.1923e-03], (B) 6.817e-15, (C) 1.309e-13; (c1) 1.0649e+00, no re-solve; (IT) 6.079e-12 (Python Study + thermal 17.8 s). W2 (`20261007T051531Z_WF-1-s3-W2-wf1-w1.log`, 115 s, 15 passed), W3 (`20261007T051729Z_WF-1-s3-W3-wf1-w2.log`, 164 s, 13 passed), W4 (`20261007T052017Z_WF-1-s3-W4-wf14.log`, 5 s, 6 passed, 34 registry keys, 91 node refs, 0 violations), W5 (`20261007T052023Z_WF-1-s3-W5-installed-femem.log`, 17 s, 1 passed). Envelope: one node under `[names.load_study]`, one under `[names.run_study_file]` (lines added only; no `[names.main]` exists — `main` is not in `__all__`).

### §9 item 343 — 2026-10-08 (chain step U32, daily provisional)

*Drafted by a general-purpose agent for the 2026-10-08 18:00 daily review; the review reviews it.* Licence: §9 rule (6), a provisional step of §10 chain U ("each step reproduces a gated result or asserts a gated identity through the API, within its unmoved band"), tagged "(daily, 2026-10-08 — provisional)"; the 10-10 weekly ratifies or reverts. **No band, tolerance, Done-when, record or validation target is set or moved; every band and record is imported by name.** Owner: `WF-1` (the YAML loader `src/fem_em_solver/study_yaml.py` and its tests are `WF-1`'s; steps 1–3 = items 256 / 293 / 295).

**`WF-1` step 4: a top-level `conductor:` section and the `geometry.as_hole` flag, so `femem run` drives U26's copper Study.** `WF-16` step 7 (item 324) made `Study(…, conductor=LeontovichWall(…))` on `FourPortBirdcageGeometry(as_hole=True, strip_fraction=GATED_WIDTH_FRACTION)` a gated API path; schema v0 cannot reach it (`_parse_geometry` refuses every `bool` field as "not settable from YAML", `study_yaml.py:216–217`, and there is no `conductor` key). This step reproduces item 324's record through `run_study_file`.
- **Executor:** implementer. **Adds (additive, `study_yaml.py` only in `src/`):**
  - `conductor` in `_TOP_KEYS` (optional), validated **after `excitation`, before `thermal`**: `kind` (required; one value, `leontovich`), `sigma_s_per_m` (optional number; omitted ⇒ `LeontovichWall`'s own default `COPPER_SIGMA`, never a loader default — the (Deg) rule), non-positive ⇒ `StudyConfigError` at `conductor.sigma_s_per_m` before `LeontovichWall` is built. `conductor` on `sphere_in_box` ⇒ error at `conductor`; `conductor` without `geometry.as_hole: true` ⇒ error at `conductor` **at parse, before any mesh** (mirrors `_conductor_setup`'s `ValueError`, `study.py:1562–1564`).
  - `geometry.as_hole`: the **only** `bool` field made settable; it must be a YAML `true` / `false` (an `int` ⇒ error at `geometry.as_hole`). Every other `bool` field keeps the existing refusal.
  - `_parse_mapping` passes `conductor=` to `Study` (omitted ⇒ `None`, the unchanged default path). If the private parse tuple grows, update every caller in the same diff (item 295's rule).
  - `run_study_file`: `summary["conductor"] = {"kind": "leontovich", "sigma_s_per_m": float}` when present (absent otherwise; existing summaries byte-unchanged in key set).
  - **YAML files** (one-line header citing this item; floats as Python `repr`): `tests/package/wf1/f_small_copper.yaml` = `f_small_four_port.yaml` + `geometry.as_hole: true` + `conductor: {kind: leontovich, sigma_s_per_m: 58000000.0}`; `f_small_hole_pec_control.yaml` = the same with **no** `conductor` section (one-place difference); three malformed files: `malformed_conductor_on_sphere.yaml` (`sphere_mat4.yaml` + the conductor section ⇒ `conductor`), `malformed_conductor_without_hole.yaml` (`f_small_copper.yaml` minus `as_hole` ⇒ `conductor`), `malformed_conductor_sigma_zero.yaml` (`sigma_s_per_m: 0.0` ⇒ `conductor.sigma_s_per_m`).
  - New `tests/package/test_wf1_yaml_copper.py` (`@complex_only` on the solve nodes); envelope nodes appended under `[names.load_study]`, `[names.run_study_file]` and `[names.LeontovichWall]` (adding lines only). No `cli.py`, `study.py`, `__init__.py` or `public_surface.txt` change.
- **Tier / ranks:** W1 standard by prediction, heavy by ceiling (container 300 s, item 324's ceiling on this fixture), `-n 2` (= `ANS6_EXACT_STRIP_RECORD_RANK_WIDTH`, asserted), degree 1, complex build. Re-gates W2–W6 below.
- **Cost:** class median ≈ 9 measured min; **≈ 55 predicted**: 15 fixed + ≈ 20 code + W1 ≈ 2 + W2 ≈ 2 + W3 ≈ 3 + W4 ≈ 1 + W5 / W6 ≈ 0.5 + ≈ 10 spare. Priced from item 324 W1 (`20261008T190119Z_WF-16-s7-W1.log`, 69 s: two Studies — Cu and PEC — plus the identity field solve, 13 passed); here three Study solves (YAML Cu, YAML PEC control, Python Cu comparand) plus the identity ⇒ predicted 90–120 s.
- **Depends on:** nothing open — `LeontovichWall`, `as_hole`, `GATED_WIDTH_FRACTION` / `COPPER_SIGMA` (`study.py:1519–1530`) and the loader are on `main`. Shares `validation_envelope.toml` with any concurrent item: edit, adding lines only; on a rebase keep both sides.

**Anchors (asserted).** Imported, none restated: `ANS6_EXACT_STRIP_CLASS_MEANS_P1`, `ANS6_EXACT_STRIP_HOLE_CELL_RECORD`, `ANS6_EXACT_STRIP_RECORD_RANK_WIDTH`, `CLASSES`, `_class_means` (`tests/validation/test_ans6_exact_strip.py`), `LEG_D0_REPRODUCTION_BAND`, `PASSIVITY_SIGMA_TOLERANCE`, `DISCRETE_IDENTITY_RTOL`, `RECIPROCITY_BAND`, `ADJACENT_SPREAD_BAND`, and **the helpers `_identity`, `_rel_means`, `_geometry`, `NEGATIVE_CONTROL_FLOOR` from `tests/package/test_wf16_study_copper.py`** (import; never copy; never import its `copper` fixture — it re-solves).
- **(K) bitwise up to the solve:** `load_study("f_small_copper.yaml")`'s spec equals `Study(_geometry(), GelledSalinePhantomMaterial(SALINE_SIGMA, SALINE_EPSILON_R, FREQUENCY_64_HZ), LumpedPortDrive(), frequency_hz=FREQUENCY_64_HZ, degree=1, conductor=LeontovichWall(sigma_s_per_m=COPPER_SIGMA))`'s field by field by `float.hex` (`WF-1`'s `_spec_diff` extended to `conductor`; diff `[]`), including `as_hole is True` and `strip_fraction == GATED_WIDTH_FRACTION`; the control file's diff is exactly `['conductor']`; a file omitting `sigma_s_per_m` parses to `LeontovichWall()`'s default.
- **(R) record through YAML:** `run_study_file("f_small_copper.yaml")` class means within `LEG_D0_REPRODUCTION_BAND` of `ANS6_EXACT_STRIP_CLASS_MEANS_P1["64"]["cu"]` per class (item 324 read 4.093e-11 / 8.673e-11 / 1.376e-10, same fixture, same width); `cells == ANS6_EXACT_STRIP_HOLE_CELL_RECORD` (80 798); `summary["s_matrix"]` equals `run.result.s_matrix` exactly (JSON round trip, `WF-1` step 1).
- **(b) copper power identity on the YAML-solved result:** `_identity(run.result, run.study.omega, comm)` gives `p_surf > 0` and residual ≤ `DISCRETE_IDENTITY_RTOL` (item 324: 8.507e-13).
- **(g) `TH-14` gates** on the YAML Cu and PEC sweeps (`th14._sweep_record(run.result.sweep)`), imported bands (item 324: green both columns).
- **(I) YAML vs Python Study:** class means of an independently solved Python-built Cu Study within `LEG_D0_REPRODUCTION_BAND` of the YAML run's. *Rule (j) reading (`WF-1` Done-when, ratified 10-05): "bitwise" stops at the solve — `-n 2` MUMPS is not bitwise across solves (`WF-1` step 1 F-small in-process (I) read 9.110e-16, not 0; known-issues `MAT-1` entry, `SOLVE_FLOOR_RECORD`), so (I) is asserted at the imported band, never `==`.*
- **(X) malformed, no solve:** the three files raise `StudyConfigError` naming `conductor` / `conductor` / `conductor.sigma_s_per_m`, in-process; `cli.main(["run", f])` and `["check", f]` return 2 on both ranks (min / max reduced), nothing written.

**Negative control (asserted, same fixture, prior measurement).** **(c1) conductor dropped:** `run_study_file("f_small_hole_pec_control.yaml")` reproduces the record's `"pec"` column within `LEG_D0_REPRODUCTION_BAND` and misses the Cu record by `> NEGATIVE_CONTROL_FLOOR` per class — item 324 read 6.420e-4 / 2.811e-4 / 4.861e-4 on this fixture at `-n 2`; ceiling not binding (≥ 28× the floor). This is the blind-parser control: a loader that silently drops `conductor:` passes (K)-free tests and fails here.

**Windows** (foreground, Bash `timeout: 660000`; each a bare whole harness call; `-s`; node IDs, never `-k`; no pipes; no edit after the last window).
- **W1, copper through YAML:** `scripts/testing/run_and_log.sh WF-1-s4-W1-copper "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 300 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_copper.py -s -v --tb=short'"`
- **W2 / W3, step-1 re-gate (rule (c): `study_yaml.py` is `WF-1`'s module):** item 256's `WF-1-W1-parse-sphere` and `WF-1-W2-f-small` commands verbatim (this file, item 256 Windows list) under tags `WF-1-s4-W2-wf1-w1` / `WF-1-s4-W3-wf1-w2`; 115 / 164 s measured at step 3.
- **W4, step-3 re-gate:** item 295's W1 command verbatim under tag `WF-1-s4-W4-thermal` (39 s measured).
- **W5, API pages (regenerate, then gate):** `scripts/testing/run_and_log.sh WF-1-s4-W5-wf17-gen "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 60 python3 scripts/docs/generate_api_reference.py --out docs/guide/api'"`, then `scripts/testing/run_and_log.sh WF-1-s4-W6-wf14-wf17 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 python3 -m pytest tests/package/test_wf14_validation_envelope.py tests/package/test_wf17_api_reference.py -v -s --tb=short'"` (the new envelope nodes move `load_study.md` / `run_study_file.md` / `LeontovichWall.md`; commit the regenerated pages).
- **W7, installed `femem check`:** item 256's `WF-1-W5-installed-femem` command verbatim under tag `WF-1-s4-W7-installed-femem` (15–17 s measured; it copies `wf1/*.yaml`, so the new files ride along).

**Traps already paid for.**
- PyYAML 1.1 floats need a dot (and a signed exponent if any): `58000000.0`, never `5.8e7` (`WF-1` step 1); YAML `true` is `bool`, `1` is `int` — the `as_hole` check must reject the `int`.
- The hole route has **no conductor cells**: `material_map` is `{}` when `as_hole` (`study.py:483–486`, item 324's first W1 raised on the missing tag 1). The loader passes nothing that re-adds it.
- `FEM_EM_ANS6_*` is refused by `ans:7`'s `main()`: set no such variable; do not run `ans:6` / `ans:7`.
- `cell_tags.values` / `assemble_scalar` are rank-local — `_identity` already reduces; do not re-implement it.
- Complex build: `source /usr/local/bin/dolfinx-complex-mode` + `FEM_EM_REQUIRE_COMPLEX=1`, `tests/environment/test_complex_mode.py` first (W1); `-k 30` on every container `timeout`. A killed window may leave an orphan FFCx `.c`: run the rubric's cache-sweep one-liner under tag `WF-1-s4-cache-sweep` only when `pgrep -c python3` reads 0.
- The cell count is mesher-dependent (0.11 image): on an image change it is re-recorded under (1\*), never re-banded. No AED number anywhere; write ≈, never `~`.

**Scope.** Covered: the 64 MHz Leontovich-copper exact-hole four-port at degree 1, `-n 2`, through YAML and `femem run` in-process. **Not claimed:** degree 2 (the `_P2` record at `-n 8`), 10 / 128 MHz, the installed console script solving the copper file, F-human, tuning, any AED statement; `ANS-6` stays 🟡; no stable schema.

**What a negative result means.** Report and stop; park on `attempt/WF-1-step4-*`; item 343 BLOCKED; journal in attempts.md. (K) red ⇒ a parser bug (fix in-slot only inside the new lines). (R) / (I) red with (K) green ⇒ the YAML route assembles a different problem — a §7 annotation with the per-class gap, never a wider band, never a re-recorded constant. (c1) not failing ⇒ the conductor key is not reaching the Study; stop. W2–W4 red ⇒ revert `study_yaml.py`, BLOCKED.

No `PREDICTED RED:` line.

**Status it can move.** No glyph (`WF-1` ✅): the row gains a **step-4 record** (copper through YAML; readings, elapsed, logs); §10 U32 DONE; `WF-14`'s envelope gains nodes; Phase 9 (9a) "`femem run study.yaml` drives the same API" covers the copper conductor.

**Outcome: ✅ DONE 2026-10-09 00:48Z** (2026-10-08 19:30 CDT slot, attempt 1). Seven windows green on first run, no in-slot fix. `study_yaml.py` additive: `conductor` in `_TOP_KEYS`, `_parse_conductor` (validated after `excitation`, before `thermal`; `sphere_in_box` and missing `as_hole: true` raise at `conductor` before any mesh, then unknown keys, `kind`, `sigma_s_per_m > 0`), `geometry.as_hole` the one settable `bool` (an `int` raises at `geometry.as_hole`), `Study(..., conductor=...)`, `summary["conductor"]` when present; the parse tuple did not grow. `_spec_diff` (`test_wf1_study_yaml.py`) now also compares `conductor`. W1 (`20261009T004108Z_WF-1-s4-W1-copper.log`, 94 s, 11 passed): (K) copper diff `[]`, PEC-control diff `['conductor']`, sigma hex equal to `COPPER_SIGMA`, `as_hole True`, `strip_fraction 0.5`; omitted `sigma_s_per_m` → `LeontovichWall()`'s default; `as_hole: 1` → `geometry.as_hole`; (X) `conductor` / `conductor` / `conductor.sigma_s_per_m`, run / check rc (2, 2), nothing written; (R) YAML Cu 4.093e-11 / 8.673e-11 / 1.376e-10 (band 1e-09), cells 80 798, JSON `s_matrix` equal exactly, `run_study_file` 29.1 s; (b) `p_surf` 1.517e-06 W > 0, residual 8.507e-13; (g) Cu reciprocity 6.478e-16, max σ 0.999851545, spreads ≤ 6.575e-4; PEC 8.681e-16, 0.999904317; (I) Python Cu vs YAML Cu 0.000e+00 on all three classes (same process; asserted at the band, rule (j)); (c1) PEC control vs `"pec"` record 5.052e-11 / 3.878e-11 / 5.875e-11, vs Cu record 6.422e-4 / 2.812e-4 / 4.863e-4 (predicted 6.420e-4 / 2.811e-4 / 4.861e-4; floor 1e-05). Re-gates: W2 (`20261009T004249Z_WF-1-s4-W2-wf1-w1.log`, 99 s, 15 passed), W3 (`20261009T004433Z_WF-1-s4-W3-wf1-w2.log`, 142 s, 13 passed), W4 (`20261009T004659Z_WF-1-s4-W4-thermal.log`, 33 s, 6 passed), W5 pages regenerated (`20261009T004737Z_WF-1-s4-W5-wf17-gen.log`, 2 s; `LeontovichWall.md` / `load_study.md` / `run_study_file.md` moved), W6 (`20261009T004743Z_WF-1-s4-W6-wf14-wf17.log`, 8 s, 12 passed, 40 registry keys, 131 node refs, 0 violations), W7 installed `femem check` (`20261009T004757Z_WF-1-s4-W7-installed-femem.log`, 15 s, 1 passed). Envelope: one node under `[names.load_study]`, three under `[names.run_study_file]`, four under `[names.LeontovichWall]` (lines added only).

### §9 item 344 — 2026-10-08 (chain step U33, daily provisional)

*Drafted by a general-purpose agent for the 2026-10-08 18:00 daily review; the review reviews it.* Licence and tag as item 343. **No band, tolerance, Done-when, record or validation target is set or moved.** **Depends on item 343** (same module and test family; serial).

**`WF-1` step 5: `outputs.sar` and `outputs.b1_plus` in the YAML schema, driving `Study.sar` and `Study.b1_plus_at`.** `WF-16` steps 4 / 5 (items 321 / 322) made the 10 g mass-averaged SAR and `|B₁⁺|` port-Study results; schema v0's `outputs` section holds only `directory`. This step reproduces U23's 10 g C4 identity and U24's `|B₁⁺|` C4 covariance identity on a YAML-built, YAML-run Study, and gates the written summary against the API.
- **Executor:** implementer. **Adds (additive, `study_yaml.py` only in `src/`):**
  - `_parse_outputs` gains two optional keys (port geometry only; on `sphere_in_box` ⇒ `StudyConfigError` at `outputs.sar` / `outputs.b1_plus`):
    - `sar: {mass_g: <number > 0>, centres: …}` — `centres` is **either** a list of `[x, y, z]` triples **or** the mapping `{ring_radius_m: <number > 0>}`, meaning one centre per port at `(r·cos az, r·sin az, 0.0)` with `az = radians(sheet["azimuth_deg"])` read off the solved Study's own sheets (`fields_result.sheets`) — the exact expression `tests/package/test_wf16_study_sar.py:155–160` uses, so it reproduces those centres bitwise. `mass_g ≤ 0` ⇒ error at `outputs.sar.mass_g` before any mesh.
    - `b1_plus: {points: [[x, y, z], …]}` — a non-empty list of triples (a bad entry ⇒ error at `outputs.b1_plus.points[i]`).
  - `run_study_file`, port branch, after the S block: if `sar` ⇒ `study.sar(mass_g=…, centres=<all centres>, comm=comm)` (every drive at every centre; 4 × 4 = 16 operator calls) and `summary["sar"] = {"mass_g", "radius_m", "centres": [[x,y,z],…], "values": {port: [mass_averaged_sar dict, float()-ed, per centre]}}`; if `b1_plus` ⇒ `study.b1_plus_at(points, comm=comm)` and `summary["b1_plus"] = {"points", "values": {port: [float]}, "valid": {port: [bool]}}`. Both are collective; rank 0 writes after them (the existing ordering). `Study.solve` is **not** rerouted: the S block keeps `study.solve`, and `sar` / `b1_plus` use the cached `_solve_port_fields` (two solves; disclosed, not optimised here).
  - `StudyRun` keeps its fields; the summary is the carrier.
  - **YAML files:** `tests/package/wf1/f_small_sar_b1.yaml` = `f_small_four_port.yaml` + `outputs: {sar: {mass_g: 10.0, centres: {ring_radius_m: 0.015}}, b1_plus: {points: [[0.0, 0.0, 0.0]]}}`; `malformed_sar_on_sphere.yaml` (`sphere_mat4.yaml` + the `sar` block ⇒ `outputs.sar`); `malformed_sar_mass_zero.yaml` (`mass_g: 0.0` ⇒ `outputs.sar.mass_g`).
  - New `tests/package/test_wf1_yaml_sar_b1.py`; envelope nodes under `[names.load_study]` / `[names.run_study_file]` (adding lines only).
- **Tier / ranks:** W1 standard by prediction, heavy by ceiling (container 400 s, item 321's W1 ceiling on this fixture), `-n 2`, degree 1 (the YAML file's `degree: 1`), complex build.
- **Cost:** class median ≈ 9 measured min; **≈ 60 predicted**: 15 fixed + ≈ 22 code + W1 ≈ 3 + re-gates ≈ 8 + ≈ 12 spare. Priced from item 321 W1 (`20261008T184215Z_WF-16-s4-W1-deg1.log`, 114 s; field-keeping solve 36.41 s, 9 operator calls 21.00 s ⇒ ≈ 2.3 s per call) and item 322 W1 (`20261008T185141Z_WF-16-s5-W1-deg1.log`, 43 s; 4 CG1 projections 2.69 s, evaluations 0.89 s); here `study.solve` ≈ 37 s + field solve ≈ 36 s + 16 calls ≈ 37 s + B1 ≈ 4 s + the test's API comparands ≈ 12 s ⇒ predicted 120–160 s.

**Anchors (asserted).** Imported, none restated: from `tests/package/test_wf16_study_sar.py` / `test_wf16_study_four_port.py` the helpers `_port_study`, `_record`, `_rel`; `CENTRE_RADIUS_M`, `EXACT_IDENTITY_RTOL` (`test_birdcage_sar_mass_averaged.py:149–170`); `C4_COVARIANCE_BAND`, `_sample_points`, `_rotate_z` (`test_birdcage_b1_plus_map.py`); `TEN_GRAM_KG`; `LEG_D0_REPRODUCTION_BAND`. Never import a module-scope solving fixture.
- **(K) bitwise up to the solve:** the YAML Study's spec diff against `_port_study(degree=1)` is `[]`; parsed `mass_g == 10.0` and `1000 · TEN_GRAM_KG` by `float.hex`; parsed `ring_radius_m == CENTRE_RADIUS_M` by `float.hex`; after the run, `summary["sar"]["centres"]` equals the test's own sheet-derived centres (`test_wf16_study_sar.py:155–160` expression on `run.study.fields_result.sheets`) by `float.hex`.
- **(R) same fixture:** `summary["s_matrix"]` within `LEG_D0_REPRODUCTION_BAND` of `_record(1)` (`P1`; item 321 (R) 9.557e-11).
- **(ii) 10 g C4 identity from the summary:** the four cyclic pairs of the diagonal `SAR(c_{k+1}; P_{k+1})` vs `SAR(c_k; P_k)`, read from `summary["sar"]["values"]`, each ≤ `C4_COVARIANCE_BAND` — item 321 read 0.8702 / 3.5284 / 2.0337 / 4.5024 % at degree 1, `-n 2`, same fixture (worst 0.90× band; the `-n 2` cross-process floor is ≈ 1e-10 relative, `SOLVE_FLOOR_RECORD`, so the reading cannot move by a visible amount).
- **(iii) YAML = API:** for each drive k, `run.study.sar(mass_g=10.0, centres=[c_k], drives=[P_k])` on the run's own (cached) fields equals the summary's entry with relative difference ≤ `EXACT_IDENTITY_RTOL` (item 321 (iii) read 0.0 for API vs operator on one field set); the JSON summary equals the in-memory dicts exactly.
- **(B1) `|B₁⁺|` C4 covariance on the YAML-solved Study:** `_sample_points({"mesh", "cell_tags"} of run.study.fields_result)`, rotated by Δ read off the sheets; `run.study.b1_plus_at(...)` for P2@+Δ / P4@−Δ / P3@2Δ against P1, each ≤ `C4_COVARIANCE_BAND` (item 322 read 1.8695 / 1.9494 / 1.6405 %); (B3) 51 of 51 valid.
- **(B2′) YAML = API on the written points:** `summary["b1_plus"]["values"]` equals `run.study.b1_plus_at(points)` on the same Study, max |Δ| / max |B₁⁺| ≤ `EXACT_IDENTITY_RTOL` (item 322 (B2) read 0.0 at degree 1); `valid` all true.
- **(X) malformed, no solve:** the two files raise `StudyConfigError` naming `outputs.sar` / `outputs.sar.mass_g`, in-process; `cli.main` `run` / `check` rc 2 on both ranks; nothing written.

**Negative controls (asserted, same fixture, prior measurement).**
- **(c1) mis-paired far-side ball, from the summary:** `SAR(c_{k+2}; P_k)` vs `SAR(c_k; P_k)` > `C4_COVARIANCE_BAND` at every k — item 321 read 88.12 / 87.94 / 88.66 / 88.93 % at degree 1; ceiling 100 % (non-negative ratio), not binding. It proves the summary's centre ↔ drive indexing is not blind.
- **(c2) mis-rotation:** P2 read at the −Δ points against P1 > `C4_COVARIANCE_BAND` — item 322 read 35.0040 % at degree 1 (its 80–100 % size prediction missed and was reported; the sign was asserted and green).

**Windows** (foreground, Bash `timeout: 660000`; bare whole commands; `-s`; no pipes; no edit after the last window).
- **W1:** `scripts/testing/run_and_log.sh WF-1-s5-W1-sar-b1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 400 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_sar_b1.py -s -v --tb=short'"`
- **W2 / W3 / W4, re-gates:** item 256's two solve commands and item 295's W1 verbatim under tags `WF-1-s5-W2-wf1-w1` / `WF-1-s5-W3-wf1-w2` / `WF-1-s5-W4-thermal`.
- **W5, step-4 re-gate:** item 343's W1 command verbatim under tag `WF-1-s5-W5-copper`.
- **W6 / W7, API pages:** item 343's W5 / W6 commands under tags `WF-1-s5-W6-wf17-gen` / `WF-1-s5-W7-wf14-wf17`.
- **W8, installed `femem check`:** item 256's W5 command verbatim under tag `WF-1-s5-W8-installed-femem`.

**Traps already paid for.**
- Keep `SAR_QUADRATURE_DEGREE` (16) pinned through `Study.sar`'s default — the loader passes no quadrature degree (an unpinned `SpatialCoordinate` integral can compile > 9 min, `POST-5`); the ball is a UFL conditional — compare squares, never `sqrt` (`WF-6` 4c).
- `post.faraday.b1_plus` is a DG0 operator and is **not** applied to the CG1 projection; `|B₁⁺|` is formed after point evaluation (item 322 disclosure). Point evaluation only through `evaluate_vector_field_parallel`; `points` identical on every rank (the parsed list is).
- The 10 g column at a_10g / h = 0.891 on the 0.015 m phantom (item 321): **no `phantom_resolution` key** is added to make anything pass (ruling 321-1).
- A literal centre list would need the sheet azimuths to the last bit; the `ring_radius_m` form exists so (K) can be bitwise. Never hand-type centre literals.
- Orphan `.c` sweep only at `pgrep -c python3` = 0 (tag `WF-1-s5-cache-sweep`); complex build + `test_complex_mode.py` first; `-k 30`.

**Scope.** Covered: 10 g SAR at the four port-azimuth centres and `|B₁⁺|` at listed points through YAML, degree 1, `-n 2`, 64 MHz, item 192's fixture. **Not claimed:** degree 2, a spatial-peak search, C95.3 / IEC compliance, absolute SAR or `|B₁⁺|` (the printed values stay per-1 V records), homogeneity, F-human, the copper fixture's SAR, a stable schema.

**What a negative result means.** Report and stop; park on `attempt/WF-1-step5-*`; item 344 BLOCKED; journal. (ii) red with (iii), (R) and (c1) green would contradict item 321 on the same fixture — a §7 annotation plus a known-issues entry, never a wider band, never a re-mesh. (K) / (iii) / (B2′) red ⇒ plumbing in the new lines (fix in-slot only there). Controls not failing ⇒ blind instrument; stop.

No `PREDICTED RED:` line.

**Status it can move.** No glyph: the `WF-1` row gains a **step-5 record**; §10 U33 DONE; envelope nodes; (9a) "`femem run` drives the same API" covers SAR and B₁⁺.

**Outcome: ✅ DONE 2026-10-09 01:02Z** (2026-10-08 19:30 CDT slot, take-next after item 343, attempt 1). Eight windows green on first run, no in-slot fix. `study_yaml.py` additive: `_parse_outputs(sec, geometry_kind, source)` accepts `sar` / `b1_plus` (sphere ⇒ `outputs.<key>`; `sar`: unknown keys, missing `mass_g` / `centres`, `mass_g > 0`, `centres` a triple list or `{ring_radius_m > 0}`; `b1_plus.points` a non-empty triple list, a bad entry named by index); `run_study_file` port branch, after the S block, `_sar_summary` (every drive at every centre, centres in `fields_result.sheets` order, `radius_m` = `averaging_ball_radius`) and `_b1_plus_summary`. W1 (`20261009T005219Z_WF-1-s5-W1-sar-b1.log`, 119 s, 13 passed; `run_study_file` 106.9 s incl. both solves, 117 206 cells): (K) spec diff `[]`, `mass_g` `0x1.4p+3` = `1000·TEN_GRAM_KG`, `ring_radius_m` hex = `CENTRE_RADIUS_M`, summary centres = the test's sheet-derived centres by `float.hex`; (R) 9.557e-11; (ii) 0.8702 / 3.5284 / 2.0337 / 4.5024 % (item 321's readings to 4 digits); (iii) 0.000e+00, JSON = memory; (B1) 1.8695 / 1.9494 / 1.6405 %, (B3) 51 / 51 on every image; (B2′) 0.000e+00, centre `|B₁⁺|` printed per 1 V; (X) `outputs.sar` / `outputs.sar.mass_g`, run / check rc (2, 2), nothing written; (c1) 88.1162 / 87.9366 / 88.6600 / 88.9315 %; (c2) 35.0040 %. Re-gates: W2 (`20261009T005434Z_WF-1-s5-W2-wf1-w1.log`, 101 s, 15 passed), W3 (`20261009T005618Z_WF-1-s5-W3-wf1-w2.log`, 142 s, 13 passed), W4 (`20261009T005843Z_WF-1-s5-W4-thermal.log`, 34 s, 6 passed), W5 (`20261009T005921Z_WF-1-s5-W5-copper.log`, 93 s, 11 passed), W6 pages regenerated (`20261009T010056Z_WF-1-s5-W6-wf17-gen.log`, 1 s; `load_study.md` / `run_study_file.md` moved), W7 (`20261009T010101Z_WF-1-s5-W7-wf14-wf17.log`, 8 s, 12 passed, 40 registry keys, 140 node refs, 0 violations), W8 installed `femem check` (`20261009T010109Z_WF-1-s5-W8-installed-femem.log`, 14 s, 1 passed). Envelope: one node under `[names.load_study]`, eight under `[names.run_study_file]` (lines added only; the `{ring_radius_m}` convenience disclosed in the comment).

### §9 item 345 — 2026-10-08 (chain step U34, daily provisional)

*Drafted by a general-purpose agent for the 2026-10-08 18:00 daily review; the review reviews it.* Licence and tag as item 343. **No band, tolerance, Done-when, record or validation target is set or moved.** **Depends on item 344** (same module; serial).

**`WF-1` step 6: `outputs.probes` in the YAML schema, driving `Study.probe(probes, cooling)`.** `WF-15` step 4 (item 326, U30) made `MAT-19`'s ball-average probe a `Study` result on `MAT-12`'s sphere under the steady Robin temperature; schema v0 cannot ask for it. This step reproduces U30's (A) identity through `run_study_file`.
- **Executor:** implementer. **Adds (additive, `study_yaml.py` only in `src/`):**
  - `_parse_outputs` gains `probes: [{centre: [x, y, z], radius: <number > 0>}, …]` (non-empty) ⇒ a tuple of `ProbeBall`; it needs a `thermal:` section (absent ⇒ `StudyConfigError` at `outputs.probes`) and a `sphere_in_box` geometry (the thermal section already refuses a port geometry); radius ≤ 0 ⇒ error at `outputs.probes[i].radius` **before** `ProbeBall` is built (no bare `ValueError` escapes). `levels` is not exposed (the API default 3 applies).
  - `run_study_file`, sphere branch, after the temperature block: `pr = study.probe(probes, cooling)`; `summary["probes"] = {"centres", "radii", "values_k", "volumes_m3"}` (floats); `StudyRun` gains `probes=None` (non-exported; holds the `ProbeStudyResult`).
  - **YAML files:** `tests/package/wf1/sphere_mat10_probes.yaml` = `sphere_mat10_thermal.yaml` + `outputs: {probes: [whole ball (centre [0.0, 0.0, 0.0], radius 0.02 = 2R), then the three `CENTRES` with radius `0.002` = `A_PROBE`]}`, every float the Python `repr` of the imported value (`CENTRES[2]` = `6.5e-3 * np.ones(3) / np.sqrt(3.0)` — IEEE basic operations, so the host `repr` is the container's); `malformed_probes_without_thermal.yaml` (`sphere_mat4.yaml` + `outputs.probes` ⇒ `outputs.probes`); `malformed_probe_radius_zero.yaml` (⇒ `outputs.probes[1].radius`).
  - New `tests/package/test_wf1_yaml_probes.py`; envelope nodes under `[names.load_study]` / `[names.run_study_file]`.
- **Tier / ranks:** W1 standard (container 180 s), `-n 2`, EM degree 1, P1 thermal, complex build.
- **Cost:** class median ≈ 9 measured min; **≈ 45 predicted**: 15 fixed + ≈ 15 code + W1 ≈ 1 + re-gates ≈ 9 + ≈ 5 spare. Priced from item 326 W1 (`20261008T201316Z_WF-15-s4-W1-n2.log`, 22 s; EM + mesh 13.9 s, temperature + probe 1.8 s, 2.84 GiB summed) and item 295 W1 (39 s, two EM solves); here the YAML run plus one independent Python comparand ⇒ predicted 40–60 s.

**Anchors (asserted).** Imported, none restated: `B_BAND`, `K`, `H`, `T_INF`, `G`, `R`, `_allsum` (`tests/validation/test_mat10_sar_to_heat_sphere.py`), `A_PROBE`, `CENTRES`, `N_SUB`, `WIDTH_BAND` (`test_mat19_virtual_temperature_probe.py`), `WHOLE`, `PROBES`, `C2_OUT` (`tests/package/test_wf15_study_probe.py`), `L0_BAND`, `_study`.
- **(K) bitwise up to the solve:** the parsed probes equal `[WHOLE] + PROBES` centre by centre and radius by radius by `float.hex`; the cooling equals `RobinCooling(K, H, T_INF, source_scale=G)` (diff `[]`); the Study spec diff against `_study(degree=1)` is `[]`.
- **(A0)** `run.probes.thermal.submesh` global cells `== N_SUB` (19 496).
- **(A) P8's whole-ball identity through YAML:** `summary["probes"]["values_k"][0]` against ∫θ dx / V_sub on the run's own θ (`run.probes.theta`, `_allsum`), and the whole-ball volume against V_sub, each ≤ `B_BAND` (item 326 W1 read 6.040e-16 / 1.221e-15 at `-n 2`).
- **(U′) YAML = API on the run's θ:** each summary value equals `ball_average` called directly on `run.probes.theta` exactly (item 326 (U): equal); the JSON summary equals the in-memory values exactly.
- **(IT) YAML vs an independent Python Study:** `_study(degree=1).probe([WHOLE] + PROBES, RobinCooling(K, H, T_INF, source_scale=G)).values` within `L0_BAND` relative of the summary's, per probe (two Study instances: item 295 (IT) read 6.079e-12 on θ for this fixture; `-n 2` MUMPS not bitwise across solves, so never `==`).
- **(X) malformed, no solve:** the two files name `outputs.probes` / `outputs.probes[1].radius`; `cli.main` rc 2 on both ranks; nothing written.

**Negative control (asserted, same fixture, prior measurement).** **(c1) probe displaced 0.1 mm outward:** `ball_average` at `C2_OUT` on the run's θ misses the summary's c₂ value by `> WIDTH_BAND` relative — item 326 read 7.5333e-03 (`-n 2`, same fixture); ceiling not binding (7.5e6× the band). It shows a centre one parse-error away from the written one is visible.

**Windows** (foreground, Bash `timeout: 660000`; bare whole commands; `-s`; no pipes; no edit after the last window).
- **W1:** `scripts/testing/run_and_log.sh WF-1-s6-W1-probes "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_probes.py -s -v --tb=short'"`
- **W2 / W3 / W4, re-gates:** item 256's two solve commands and item 295's W1 verbatim under tags `WF-1-s6-W2-wf1-w1` / `WF-1-s6-W3-wf1-w2` / `WF-1-s6-W4-thermal`.
- **W5 / W6, steps 4–5 re-gates:** item 343's W1 and item 344's W1 verbatim under tags `WF-1-s6-W5-copper` / `WF-1-s6-W6-sar-b1`.
- **W7 / W8, API pages:** item 343's W5 / W6 commands under tags `WF-1-s6-W7-wf17-gen` / `WF-1-s6-W8-wf14-wf17`.
- **W9, installed `femem check`:** item 256's W5 command verbatim under tag `WF-1-s6-W9-installed-femem`.

**Traps already paid for.**
- `Study.probe` reuses `self.result` (one EM solve, item 326 (Once)); `run_study_file` must not call `study.solve` twice. Calling `study.temperature` and then `study.probe` re-solves the cheap thermal problem only; disclose which θ the summary's `temperature` block and the probe values come from (they may differ at the solver floor, never at `L0_BAND`).
- `ufl.max_value` / complex ordering do not compile (`OPS-22`): nothing new in UFL here; θ max stays NumPy.
- (c2) of item 326 (ghost-included sum) was **inapplicable** at `-n 2` (G_w = 0) and is not re-posed here.
- Write every literal as the imported value's `repr`; (K) catches a one-ULP error before the solve. Orphan `.c` sweep only at `pgrep -c python3` = 0 (tag `WF-1-s6-cache-sweep`). Write ≈, never `~`.

**Scope.** Covered: the steady Robin-cooled ball-average probes on `MAT-12`'s sphere through YAML at `-n 2`, EM degree 1. **Not claimed:** width invariance (U30's `-n 1` leg and (W)), transient / duty-cycled probes, the closed-form probe means (printed only in U30), any probe on a coil-driven phantom, a stable schema.

**What a negative result means.** Report and stop; park on `attempt/WF-1-step6-*`; item 345 BLOCKED; journal. (K) red ⇒ parser bug. (A) / (IT) red with (K) green ⇒ the YAML route builds a different probe problem — a §7 annotation, never a wider band. (c1) not failing ⇒ blind instrument; stop. A re-gate red ⇒ revert `study_yaml.py`, BLOCKED.

No `PREDICTED RED:` line.

**Status it can move.** No glyph: the `WF-1` row gains a **step-6 record**; §10 U34 DONE; envelope nodes; (9a) "`femem run` drives the same API" covers the virtual probe.

**Outcome: ✅ DONE 2026-10-09 02:12Z** (2026-10-08 21:00 CDT slot, attempt 1). Nine windows green on first run, no in-slot fix. `study_yaml.py` additive: `_parse_probes` (non-sphere ⇒ `outputs.probes`; no `thermal:` ⇒ `outputs.probes`; non-empty list of `{centre, radius}` mappings, unknown / missing keys named at `outputs.probes[i].<key>`, `radius > 0` checked at `outputs.probes[i].radius` before `ProbeBall` is built), `_parse_outputs(..., has_thermal=)` (the parse tuple did not grow), `run_study_file` sphere branch after the temperature block: `study.probe(probes, cooling)` → `summary["probes"] = {centres, radii, values_k, volumes_m3}`, `StudyRun.probes`. Disclosed: the probe values come from `Study.probe`'s own thermal solve, a second cheap thermal solve beside the temperature block's (one EM solve; the two θ differ by 8.404e-14 owned max rel, printed). YAML literals are the imported values' `repr` (`CENTRES[2]` component `0.0037527767497325675`). W1 (`20261009T020231Z_WF-1-s6-W1-probes.log`, 40 s, 9 passed; `run_study_file` 16.9 s): (K) four probes hex-equal to `[WHOLE] + PROBES`, cooling diff `[]`, spec diff `[]`; (X) `outputs.probes` / `outputs.probes[1].radius`, run / check rc (2, 2), nothing written; (A0) 19 496; (A) 1.812e-15 / 1.221e-15 (band 1e-10; item 326 read 6.040e-16 / 1.221e-15); (U′) summary values = direct `ball_average` exactly, JSON = memory; (IT) per probe 6.777e-12 / 6.826e-12 / 4.405e-12 / 1.401e-11 (band 1e-09; Python Study + probe 17.4 s); (c1) 7.5333e-03 (item 326: 7.5333e-03). Values: whole ball 0.0919062 K, c0 0.1250640 K, c1 0.1179883 K, c2 0.1006763 K (per the fixture's G). Re-gates: W2 (`20261009T020327Z_WF-1-s6-W2-wf1-w1.log`, 102 s, 15 passed), W3 (`20261009T020513Z_WF-1-s6-W3-wf1-w2.log`, 142 s, 13 passed), W4 (`20261009T020735Z_WF-1-s6-W4-thermal.log`, 33 s, 6 passed), W5 (`20261009T020810Z_WF-1-s6-W5-copper.log`, 95 s, 11 passed), W6 (`20261009T020945Z_WF-1-s6-W6-sar-b1.log`, 122 s, 13 passed), W7 pages regenerated (`20261009T021148Z_WF-1-s6-W7-wf17-gen.log`, 1 s; `load_study.md` / `run_study_file.md` moved), W8 (`20261009T021150Z_WF-1-s6-W8-wf14-wf17.log`, 9 s, 12 passed, 40 registry keys, 145 node refs, 0 violations), W9 installed `femem check` (`20261009T021200Z_WF-1-s6-W9-installed-femem.log`, 16 s, 1 passed). Envelope: one node under `[names.load_study]`, four under `[names.run_study_file]` (lines added only).

## §7 row moved verbatim 2026-10-05 (rotate_plan_archive.py chunks)

 MRI example CLI/config *(2026-10-03 weekly: returns as §10 chain U5 — `femem run study.yaml` drives `WF-15` / `WF-16` to bit-identical outputs and three malformed files fail with a named field; serial on `WF-16`; the daily writes the Done-when at rubric grain)* *(deferred to Phase 9 (9a / 9d) by the 2026-09-30 weekly — no validation target, unmoved over a month; returns as package work when Phase 9 opens; §10 epitaph 2026-09-30; not a lane-1 item)* *(history and §9 briefs: `docs/planning/chunks/WF-1.md`)* | ✅ 2026-10-06 (§9 item 256): `femem run` on `tests/package/wf1/sphere_mat4.yaml` / `f_small_four_port.yaml` (degree 1, `-n 2`) builds Studies whose spec equals the Python-built one field by field by `float.hex` (diff `[]` on both, and on the "mid" preset route; degree 2 when omitted), DG0 σ / εᵣ / μᵣ max \|Δ\| 0.0 on 74 020 owned dofs; sphere SAR 6.159e-11 from `MAT4_DEG1_RECORD` (`.6e` string equal), 8.244e-13 from the in-process Python Study, closed form 3.536 %; F-small S 9.557e-11 from `P1`, Touchstone (T1) 6.873e-14 / (T2) 9.557e-11, in-process (I) 9.110e-16, `PORT-11` gates green; the three malformed files raise `StudyConfigError` naming `geometry.sphere_raduis` / `study.frequency_hz` / `excitation` in-process (rc 2 on both ranks, nothing written) and through the installed `femem` in `WF-13`'s venv; controls (c1) exactly 2⁻⁵³ on 19 496 dofs, (c2) 3.4937 %, (c3) 0.0491, (c4) 0.2492. Elapsed: W1 100 s, W2 142 s, W3 5 s, W4 37 s, W5 15 s (`20261006T003504Z`…`004014Z_WF-1-*.log`). Not claimed: a degree-2 solve through YAML, `femem run` with a solve through the installed script, the `__main__` route (unexercised), a stable schema. *(History: WF-1.md.)* **Done-when (§4):** `femem run <file>.yaml` (console script `fem_em_solver.cli:main`) drives `WF-15`'s lossy sphere and `WF-16`'s F-small four-port, each at degree 1 and `-n 2`, to bit-identical outputs *(read by the 2026-10-05 18:00 review, rule (j); the weekly ratifies: the YAML-built `Study` equals the Python-built one field by field, floats by `float.hex`, degree 2 when omitted; the solve outputs — SAR, S and the Touchstone file — equal `MAT-4`'s / item 192's records and the Python `Study`'s within `L0_BAND` / `LEG_D0_REPRODUCTION_BAND`)*; three malformed files (unknown key, wrong type, missing section) raise `StudyConfigError` naming the field, in-process and through the installed `femem`, while a one-value variant parses and changes only that field; entry points registered in `WF-14`'s envelope; elapsed recorded ⇒ ✅. 

## §7 row moved verbatim 2026-10-07 (rotate_plan_archive.py chunks)

 MRI example CLI/config *(2026-10-03 weekly: returns as §10 chain U5 | ✅ 2026-10-06 (§9 item 256, chain U5): `femem run` on the sphere and F-small YAMLs builds Studies equal to the Python-built ones by `float.hex`, reproduces `MAT4_DEG1_RECORD` and `P1` (SAR 6.159e-11, S 9.557e-11), and the three malformed files raise `StudyConfigError` naming the field in-process and through the installed `femem`. Rule (j) reading of "bit-identical" (bitwise up to the solve, imported bands after it) ratified by the 2026-10-05 weekly; not claimed: a degree-2 solve through YAML or a stable schema. **Step-2 record (§9 item 293, chain U18, 2026-10-07):** the installed console script `/tmp/wf13/venv/bin/femem run` under `mpiexec -n 2`, cwd `/tmp/wf13/run/u18`, no `PYTHONPATH` (venv imports `fem_em_solver` from its own site-packages; sha256 of each copied YAML equal to the repo file's) reproduces sphere SAR 5.909e-11 from `MAT4_DEG1_RECORD` (`.6e` string equal, closed form 3.536 %) and F-small S 9.557e-11 from `P1` (installed-written `.s4p` 6.855e-14 from the summary for both readers, `PORT-11` gates green, 117 206 cells); controls (c2) 3.4937 %, (c3) 0.0491, (c4) 0.2492. Elapsed W1 41 s, W2 65 s, W3 (`WF-14`) 5 s (`20261007T050136Z` / `050226Z` / `050336Z_WF-1-s2-*.log`) — "`femem run` with a solve through the installed script" is retired from Not claimed; still not claimed: a STEP geometry, degree 2 through YAML, other widths, a published image or PyPI. The `__main__` route (step-1 "unexercised") is retired by reference: `WF-13` step 3 (§9 item 294) launches every child through `python -m fem_em_solver` and gates its output. **Step-3 record (§9 item 295, chain U20, 2026-10-07 00:00 slot):** an optional `thermal:` section (`k_w_per_m_k`, `h_w_per_m2_k`, `t_inf_c`; optional `source_scale`, `degree` — omitted ⇒ `RobinCooling`'s own defaults) drives `Study.temperature(RobinCooling)` through `run_study_file`; `summary["temperature"]` carries source / convected power and the `MPI.MAX`-reduced owned θ max (0.126822 K). On `tests/package/wf1/sphere_mat10_thermal.yaml` at degree 1, `-n 2`: (K) the parsed `RobinCooling` equals `RobinCooling(K, H, T_INF, source_scale=G)` by `float.hex` (diff `[]`), Study spec diff `[]`; (R) sub cells 19 496, (A) 3.2464e-02 / (P) 3.1923e-03 equal to `MAT-10`'s `.4e` records, (B) 6.817e-15, (C) 1.309e-13, JSON summary equal to the in-memory values exactly; (IT) vs an independently solved Python Study 6.079e-12 ≤ `L0_BAND`; (X) `thermal` on the F-small file and `thermal.k_w_per_m_k: 0.0` raise `StudyConfigError` naming `thermal` / `thermal.k_w_per_m_k`, rc 2 on both ranks, nothing written; control (c1) the 2G control file on the solved Study reads 1.0649e+00 > `A_BAND`, no re-solve. W1 39 s (`20261007T051444Z_WF-1-s3-W1-thermal.log`, 6 passed); rule (c) re-gates step-1 W1 / W2 15 / 13 passed (115 / 164 s), `WF-14` 6 passed (91 node refs; nodes under `load_study` / `run_study_file`), installed `femem check` 1 passed (17 s). Not claimed: degree 2 through YAML, transient / Pennes / duty cycle, thermal on a port Study, a stable schema. **Done-when:** `femem run <file>.yaml` (console script `fem_em_solver.cli:main`) drives `WF-15`'s lossy sphere and `WF-16`'s F-small four-port, each at degree 1 and `-n 2`, to bit-identical outputs *(read by the 2026-10-05 18:00 review, rule (j); the weekly ratifies: the YAML-built `Study` equals the Python-built one field by field, floats by `float.hex`, degree 2 when omitted; the solve outputs — SAR, S and the Touchstone file — equal `MAT-4`'s / item 192's records and the Python `Study`'s within `L0_BAND` / `LEG_D0_REPRODUCTION_BAND`)*; three malformed files (unknown key, wrong type, missing section) raise `StudyConfigError` naming the field, in-process and through the installed `femem` *History: `docs/planning/chunks/WF-1.md`.* 

**Review rulings (2026-10-08 18:00 review) on items 343–345 — these outrank the briefs above.** These are chain-U steps U32–U34 written under §9 rule (6), **provisional until the 10-10 weekly ratifies or reverts them**; no band, tolerance or target is set by them. (1) "Bitwise class means" reads as bitwise up to the solve, then YAML vs Python at the imported band (the ratified rule-(j) reading of `WF-1` step 1) — confirmed. (2) `sar` / `b1_plus` / `probes` under `outputs:`, `conductor:` top-level, as drafted. (3) The `centres: {ring_radius_m: …}` form is accepted as an API convenience, disclosed in the commit and the envelope; literal centres stay accepted too. (4) The two solves in step 5 are accepted and disclosed; no rerouting of the existing S path. (5) Keep every re-gate window as drafted (rule (c)); `study_yaml.py` is shared by both branches. Each item edits `study_yaml.py` / `cli.py`, so it regenerates the API pages and re-runs `test_wf17_api_reference.py` and `test_wf14_validation_envelope.py` in its own slot.

## §7 row moved verbatim 2026-10-10 (rotate_plan_archive.py chunks)

 `femem run`: the MRI study CLI and YAML config — Phase 9 (9a), §10 chain U5 | ✅ 2026-10-06 (§9 item 256, U5): `femem run` builds Studies equal to the Python-built ones by `float.hex` and reproduces `MAT4_DEG1_RECORD` and `P1` at their bands, with malformed files raising `StudyConfigError`. Steps 2 / 3 (items 293 / 295, 2026-10-07, U18 / U20, ratified) are records: the installed console script under `mpiexec -n 2` with no `PYTHONPATH` reproduces both records, and a YAML `thermal:` section reproduces `MAT-10`'s (A) / (P) at `.4e` and `WF-15` step 2's (IT) at `L0_BAND`. **Step-4 record (§9 item 343, U32, provisional, 2026-10-08 19:30 slot):** a top-level `conductor:` (`kind: leontovich`, optional `sigma_s_per_m`) and `geometry.as_hole` (the one settable bool) drive `WF-16` step 7's copper Study through `run_study_file` at degree 1, `-n 2`: Cu class means 4.093e-11 / 8.673e-11 / 1.376e-10 from `ANS6_EXACT_STRIP_CLASS_MEANS_P1["64"]["cu"]`, 80 798 cells, copper identity 8.507e-13, `TH-14` gates green, Python Study (I) 0.0; the PEC-control file misses Cu by 6.422e-4 / 2.812e-4 / 4.863e-4 (W1 94 s, `20261009T004108Z_WF-1-s4-W1-copper.log`); re-gated green after `4689d8bc` on `63bccf11`, every reading equal (§9 item 360, 11 passed, 98 s, `20261010T171139Z_OPS-70-item360-W21-copper.log`). **Step-5 record (§9 item 344, U33, provisional, 2026-10-08 19:30 slot):** `outputs.sar: {mass_g, centres}` (a triple list, or `{ring_radius_m}` — one centre per port at the solved sheets' azimuths, an API convenience) and `outputs.b1_plus: {points}` drive `Study.sar` / `Study.b1_plus_at` through `run_study_file` on item 192's fixture at degree 1, `-n 2` (two solves, disclosed): S 9.557e-11 from `P1`; the 10 g C4 identity from the summary 0.8702 / 3.5284 / 2.0337 / 4.5024 % and `|B₁⁺|` C4 covariance 1.8695 / 1.9494 / 1.6405 % (51 / 51 valid), both ≤ `C4_COVARIANCE_BAND`; summary = API 0.0 for both; controls (c1) 87.94–88.93 % (range corrected 2026-10-09 18:00 audit), (c2) 35.0040 % (W1 119 s, `20261009T005219Z_WF-1-s5-W1-sar-b1.log`). Not claimed: degree 2, a peak search, compliance, absolute SAR / `|B₁⁺|`. **Step-6 record (§9 item 345, U34, provisional, 2026-10-08 21:00 slot):** `outputs.probes: [{centre, radius}, …]` (sphere with `thermal:` only) drives `Study.probe(probes, cooling)` through `run_study_file` on `MAT-12`'s sphere at EM degree 1, P1 thermal, `-n 2` (one EM solve; the probe θ is `Study.probe`'s own thermal solve, 8.404e-14 owned max rel from the temperature block's — disclosed): parsed probes = `[WHOLE] + PROBES` by `float.hex`; submesh 19 496 cells; P8's whole-ball identity 1.812e-15 (value) / 1.221e-15 (volume) ≤ `B_BAND`; summary = `ball_average` on the run's θ exactly; vs an independent Python Study max 1.401e-11 ≤ `L0_BAND`; control (c1) 7.5333e-03 > `WIDTH_BAND` (W1 40 s, `20261009T020231Z_WF-1-s6-W1-probes.log`). Not claimed: width invariance, transient probes, closed-form probe means, a coil-driven phantom. *History: `docs/planning/chunks/WF-1.md`.* 

### §9 item 390a — 2026-10-11 (chain U43 part 1, daily provisional under §9 rule (6); chain step (no solve))

**Outcome 2026-10-12 (implementer, slot 12:30 UTC): DONE, green.** W0 `20261012T123618Z_WF-1-s7a-W0.log` 6 passed, 4 s: (K) `_spec_diff` `[]`, `c_in_f` `0x1.0405491074fcbp-33` = `_c_f()`, `esr_ohm` = `ESR_OHM` by hex, generator kwargs = `_f_human_kwargs(as_hole=True)`, `femem check` rc (0, 0); (K-c) 19 of 19 centres and both points by hex; (X) 3 files + 17 of 17 mutations name their field, run / check rc (2, 2), nothing written; (c1) `['excitation.esr_ohm']`, (c2) `['excitation.c_in_f']`; (Run-0) fake called (1, 1), summary = fake, JSON = summary. Literals: `20261012T123540Z_WF-1-s7a-W00-literals.log` (4 s). Re-gates (rule (c)) all green: RG1 15 passed 110 s (`…123643Z…RG1-wf1-w1`), RG2 13 / 150 s (`…123841Z…RG2-wf1-w2`), RG3 6 / 35 s (`…124119Z…RG3-thermal`), RG4 11 / 96 s (`…124201Z…RG4-copper`), RG5 13 / 128 s (`…124345Z…RG5-sar-b1`), RG6 9 / 37 s (`…124605Z…RG6-probes`), RG7 6 / 47 s (`…124650Z…RG7-wf2-sphere`), RG8 5 / 79 s (`…124747Z…RG8-wf2-four-port`); `pgrep -c python3` read 0 before RG2 / RG4 / RG5 / RG8. Pages (ruling Q3): Wgen 44 files 2 s (`…124913Z`), Wval 157 nodes 2 s (`…124919Z`), Wapi 22 passed 15 s (`…124928Z`), Wdoc 0 dead 1 s (`…124948Z`); W9 installed check 1 passed 16 s (`…124958Z`). Ruling (3) grep: `test_wf1_installed_solve.py`, `test_wf13_launch_helper.py`, `test_wf13_cli_estimate.py`, `examples/mri/09` / `13` import only public `load_study` / `run_study_file` / `cli` and `test_wf1_study_yaml` helpers — no private `_parse_*` name, and no F-human file, so the new branches are unreachable from them. Findings, not fixed (no `cli.py` edit in scope): `femem check` prints `family: four_port_birdcage` for the F-human file (`cli.py:158` treats every non-sphere geometry as four-port; W0 log line 68); M(c\*) has the same `repr` as the lattice row `r=45 dz=+0 az17+90` (the 56.25° mirror maps az 11.25° onto 101.25°), so `centre[18]` repeats `centre[5]` bit for bit (W0 log lines 80 / 93) — for 390b's reading.

*Drafted by a general-purpose agent for the 2026-10-11 18:00 daily review; the review reviews it.* Licence: §9 rule (6), a provisional step of §10 chain U ("each step reproduces a gated result or asserts a gated identity through the API, within its unmoved band"), tagged "U43 (daily, 2026-10-11 — provisional)"; the 10-12 weekly ratifies or reverts. **Target: reproduction only** — `femem run` drives item 382 / 382b's tuned copper F-human `Study`. **No band, tolerance, Done-when, record or validation target is set or moved; every band and record is imported by name.** Owner: `WF-1` (the loader `src/fem_em_solver/study_yaml.py` and its tests; steps 1–6 = items 256 / 293 / 295 / 343 / 344 / 345). **This is step 7, split in two** (below): 390a lands every `src/` change, the tracked YAML file, the no-solve gates and the rule (c) re-gates; 390b (serial) adds only a test module and runs the one heavy `-n 8` window.

**Why split.** One item would price at ≈ 15 fixed + ≈ 40 code + ≈ 11 re-gate windows + ≈ 6 heavy window ≈ 75 min — past the 1-hour box. Item 382 (61 predicted) was declined at minute 29.8 with ≈ 35 min of code unwritten (attempts.md 2026-10-11T20:31Z) and item 383 parked with its test module unwritten (2026-10-11T21:59Z). Splitting at the `src/` boundary makes 390b test-only: no `src/` edit ⇒ no rule (c) re-gates in the heavy slot, and 390a's re-gates are cheap `-n 2` windows. The cost of the split: 390a writes the `run_study_file` F-human branch without a real solve — it is exercised there through a no-mesh fake (Run-0) and for real only in 390b.

**`WF-1` step 7a: the `f_human_birdcage` geometry kind, the `two_port_quadrature` excitation and the F-human branch of `run_study_file`, gated without a mesh.**
- **Executor:** implementer. **Edits (`study_yaml.py` only in `src/`; additive where possible, every touched existing function named for rule (c)):**
  - **`geometry.kind: f_human_birdcage`** → `FHumanBirdcageGeometry`. Add `"f_human_birdcage": "two_port_quadrature"` to `_EXCITATION_FOR` (`study_yaml.py:100–103`; geometry kinds and excitation kinds are both read off it, so one entry adds both). In `_parse_geometry` (`:216–268`) a new `f_human_birdcage` branch, **copied from the four-port loop, not factored out of it** (item 373's rule: a shared helper would edit the four-port path): defaults `FHumanBirdcageGeometry()`, overrides by `dataclasses.fields`; **no `preset` key** (unknown ⇒ `geometry.preset`); `as_hole` the only settable `bool` (`_SETTABLE_BOOLS`; `emit_leg_ports` keeps "not settable from YAML"); and, in this schema version, **`as_hole` is required and must be `true`** — the YAML F-human route is the copper hole route only: missing ⇒ `geometry.as_hole` "required key is missing", `false` ⇒ `geometry.as_hole` "must be true", an `int` ⇒ `geometry.as_hole` (step 4's rule). Order inside the section as everywhere: unknown keys, missing keys, types, values.
  - **`excitation.kind: two_port_quadrature`** → `TwoPortQuadratureDrive`. `_parse_excitation` gains the parsed `geometry` argument (its one caller, `_parse_mapping`, updated in the same diff — item 295's rule). Keys: `c_in_f` (number, required, > 0), `esr_ohm` (number ≥ 0, optional), `feeds` (list of two distinct strings, optional), `z0_ohm` (number > 0, optional). **Omitted optional keys are not filled in** — the dataclass defaults apply (the (Deg) rule): `esr_ohm` 0.0, `feeds` `("P17", "P21")`, `z0_ohm` 50.0. Every value is checked here so `TwoPortQuadratureDrive.__post_init__`'s `ValueError` never escapes bare: `c_in_f ≤ 0` / non-finite ⇒ `excitation.c_in_f`; `esr_ohm < 0` ⇒ `excitation.esr_ohm`; `feeds` not a two-element list or two equal ids ⇒ `excitation.feeds`; a feed not in `{f"P{i}" for i in geometry.ring_ordinals()}` ⇒ `excitation.feeds[i]` (mirrors `_study_solve_tuned_quadrature`'s own refusal, `study.py:2094–2097`, at parse); `z0_ohm ≤ 0` ⇒ `excitation.z0_ohm`. `lumped_port_drive` on `f_human_birdcage` and `two_port_quadrature` on `four_port_birdcage` both raise at `excitation.kind` through the existing mismatch check (no new code). The four-port and sphere branches of `_parse_excitation` keep their statements and error texts byte for byte.
  - **`conductor:` on `f_human_birdcage`:** `_parse_conductor`'s geometry test (`:342`) admits `f_human_birdcage` beside `four_port_birdcage` (the sphere message unchanged); its `as_hole` test is already generic. **`two_port_quadrature` needs a conductor:** in `_parse_mapping`, after the conductor parse, `kind == "f_human_birdcage"` with no `conductor:` section ⇒ `StudyConfigError` at `conductor` ("required section is missing: two_port_quadrature needs a leontovich conductor") — before any mesh, mirroring `study.py:2092–2093`.
  - **`thermal:`** on `f_human_birdcage` ⇒ error at `thermal` (a new `elif` beside the four-port refusal in `_parse_thermal`, `:371–374`; the four-port text unchanged). **`outputs`:** `sar` / `b1_plus` are admitted on `f_human_birdcage` (`_parse_outputs`, `:485–492`, the geometry test widened); the `centres: {ring_radius_m}` form is **refused on F-human** at `outputs.sar.centres` (it reads `_solve_port_fields`' four-port sheets, `:679–688`; the F-human SAR lattice is a literal list); `probes` keeps its existing refusal (`outputs.probes`).
  - **`run_study_file` (`:714–791`) — the new branch is mandatory, not optional:** today an F-human Study would fall into the sphere `else` (`:763`) and call `result.mesh` / `mean_sar` **after** a full `-n 8` solve. Insert `elif isinstance(study.geometry, FHumanBirdcageGeometry):` between the four-port `if` (`:750`) and the `else`; the existing branches' bodies are untouched. Body: `summary["cells"] = int(result.cells)`; `summary["excitation"] = {"kind": "two_port_quadrature", "c_in_f", "esr_ohm", "feeds", "z0_ohm", "gap_impedance_ohm": [re, im], "feed_impedance_ohm": [re, im], "source_voltage_v": [re, im]}` (from the drive and `result`); `summary["az_feeds"]`, `summary["az_mirror"]` (floats); `summary["s_node"] = _pairs(result.s_node)`; `summary["y_in"] = {sense: [[re, im], …] for sense in result.y_in}`; then the **existing** `_sar_summary` / `_b1_plus_summary` calls if requested (both are drive-key agnostic: `study.sar` / `study.b1_plus_at` return `{sense: …}` through item 382b's hooks and reuse `study.result` — `_tuned_result`, `study.py:2233–2237` — so the file costs one solve). **No Touchstone** is written (`TunedQuadratureResult` has no `write_touchstone`, `study.py:2037–2064`; a node `s_node` is not a 50 Ω network — item 133's trap), so `paths` has no `"touchstone"` key and `summary` no `"touchstone"` / `"z_matrix"`; disclosed in the docstring. `summary["conductor"]` is already generic.
  - Module docstring: the schema list gains the three items (one line each, citing this item). **Name no artifact-shaped token (`.json`, `.s2p`, `.yaml` file names) in any new docstring** — `docrefs` reads them on the rendered API pages (item 274 lost a window to `summary.json` in a docstring).
- **Tracked YAML (one-line header citing this item; floats as Python `repr`):** `tests/package/wf1/f_human_tuned_copper.yaml` — `femem: 1`; `study: {frequency_hz: 64000000.0, degree: 1}`; `geometry: {kind: f_human_birdcage, as_hole: true}`; `material: {sigma: 0.5, epsilon_r: 78.0}`; `excitation: {kind: two_port_quadrature, c_in_f: <repr(_c_f())>, esr_ohm: 0.2}` (`_c_f()` is item 382's `C_IN_F_HUMAN_PF * 1e-12`; its repr read on the node JSON is `1.1824370000000001e-10` — W00 prints it; a literal `118.2437e-12` may differ by one ulp); `conductor: {kind: leontovich, sigma_s_per_m: 58000000.0}`; `outputs: {sar: {mass_g: 10.0, centres: <19 triples>}, b1_plus: {points: [[0.0, 0.0, 0.0], <c* triple>]}}`. **The 19 centres**, in this order: the 18 points of `_lattice(AZ17_RECORD_DEG, Z_C_RECORD_M)` (`tests/validation/test_wf10_f_human_tuned_sar10g.py:45, :58`, its own row order), then **M(c\*)** = `_mirror_xy(np.array([c*]), AZ_MIRROR_RECORD_DEG)[0]` with c\* the lattice row labelled `H9P_CSTAR_LABEL`. The examples tree has no `presets/` directory yet (U42 creates it); the file lives with `WF-1`'s fixtures and U42 may copy it.
- **New `tests/package/test_wf1_yaml_f_human.py`** (no mesh anywhere in it), module-level: `FH_FILE`; two record constants transcribed with `log:line` comments (never bands) — `Z_C_RECORD_M = 0.0` (382b's z_c, printed `0.000000e+00`, which `.6e` prints only for ±0.0: `20261011T214947Z_WF-16-s9b-W1.log:6383`) and `AZ_MIRROR_RECORD_DEG = 56.25` (382's `az_mirror 56.250000000000 vs JSON 56.250000000000 (|d| 0.000e+00)`, `20261011T213654Z_WF-16-s9a-W1.log:6369`); `_expected_centres()` (the 19 triples above) and `_literal_lines()` (the YAML text of `c_in_f`, the 19 centres and the two points, every float `repr`). Two nodes: `test_wf1_f_human_yaml_parse_and_refusals` and `test_wf1_f_human_yaml_run_branch_no_mesh`. Three malformed files (below). Envelope: the first node under `[names.load_study]`, the second under `[names.run_study_file]` (lines added only). No `cli.py`, `study.py`, `__init__.py` or `public_surface.txt` change; `study_yaml.__all__` unchanged (no new public name ⇒ no W13).
- **Tier / ranks:** W0 standard (container 180 s), `-n 2` (the width at which a rank-local rc or a non-collective raise shows), complex build; re-gates as measured. **Class:** chain step (no solve; the re-gates solve). **Cost:** ≈ 12 measured min (class median; `WF-1` steps 4–6 landed in ≈ 14–18 min each); **≈ 60 predicted**: 15 fixed + ≈ 28 code + W00 / W0 ≈ 1.5 + RG1–RG8 ≈ 9.5 (measured: 102 + 142 + 33 + 95 + 122 + 40 + 46 + 73 s) + Wgen / Wapi / Wdoc / W9 ≈ 0.6 + ≈ 5 spare. **If W0 is not green by minute 38, park** (rule (d)) — the re-gates need ≈ 11 min plus the commit.
- **Depends on:** nothing open — items 382 / 382b are on `main` (`e114da39`, `806145f3`). Shares `validation_envelope.toml` with any concurrent item: add lines only; on a rebase keep both sides.

**Anchors (asserted; all imported, none restated).** From `tests/package/test_wf1_study_yaml.py`: `_spec_diff`, `_rc_all`, `_p0`, `LOGS`, `WF1`; from `tests/package/test_wf16_study_tuned_f_human.py`: `_study`, `_c_f` (module-level helpers — 382b imports the same two); from `tests/package/test_wf16_study_tuned_f_human_b1_sar.py`: `H9P_CSTAR_LABEL`; `_lattice`, `AZ17_RECORD_DEG` (`test_wf10_f_human_tuned_sar10g.py`); `_mirror_xy` (`tests/validation/test_birdcage_b1_quadrature.py:203`); `_f_human_kwargs`, `ESR_OHM` as item 382 imports them. Never import a `test_*` name or a module-scope solving fixture.
- **(K) bitwise up to the solve:** `_spec_diff(load_study(FH_FILE), _study())` is `[]` — geometry, material, excitation, frequency, degree, density and conductor field by field, floats by `float.hex` (the `_cmp` walk handles both new dataclasses, `test_wf1_study_yaml.py:120–154`); printed: `c_in_f.hex()` against `_c_f().hex()`, `esr_ohm.hex()` against `ESR_OHM`; also `load_study(FH_FILE).geometry.generator_kwargs() == _f_human_kwargs(as_hole=True)` (item 382 (K0)'s comparand) and `cli.main(["check", FH_FILE])` rc `(0, 0)` on both ranks.
- **(K-c) the literal lattice:** the parsed `outputs.sar.centres` equal `_expected_centres()` triple by triple by `float.hex` (19 of 19); `mass_g == 10.0` by hex; the parsed `b1_plus.points` equal `[(0.0, 0.0, 0.0), c*]` by hex. This pins the transcription; whether the transcribed points sit where 382b's run-derived centres sat is 390b's printed (Δc).
- **(X) planted refusals, no mesh, each `type(e) is StudyConfigError` with the named `.field`** (never a bare `ValueError` / `TypeError` from a dataclass or the solve): three **tracked files**, run in-process and through `cli.main(["run", f])` / `["check", f]` with rc `(2, 2)` on both ranks and nothing written — `malformed_f_human_without_hole.yaml` (`as_hole` removed ⇒ `geometry.as_hole`), `malformed_quadrature_on_four_port.yaml` (`f_small_copper.yaml` with this item's excitation block ⇒ `excitation.kind`), `malformed_quadrature_without_conductor.yaml` (⇒ `conductor`); and **in-process one-key mutations** of `FH_FILE`'s mapping through `_study_from_mapping` (no files): unknown geometry key `ring_raduis` ⇒ `geometry.ring_raduis`; `emit_leg_ports: true` ⇒ `geometry.emit_leg_ports`; `as_hole: false` and `as_hole: 1` ⇒ `geometry.as_hole`; `kind: lumped_port_drive` ⇒ `excitation.kind`; unknown excitation key `c_in_pf` ⇒ `excitation.c_in_pf`; `c_in_f` removed ⇒ `excitation.c_in_f`; `c_in_f: 0.0` and `c_in_f: "1.18e-10"` (a string) ⇒ `excitation.c_in_f`; `esr_ohm: -0.1` ⇒ `excitation.esr_ohm`; `feeds: [P17]` and `feeds: [P17, P17]` ⇒ `excitation.feeds`; `feeds: [P17, P1]` ⇒ `excitation.feeds[1]`; `z0_ohm: 0.0` ⇒ `excitation.z0_ohm`; a `thermal:` section ⇒ `thermal`; `outputs.sar.centres: {ring_radius_m: 0.045}` ⇒ `outputs.sar.centres`; `outputs.probes` ⇒ `outputs.probes`. Every existing malformed file keeps its field (RG1 re-asserts the three Done-when files).
- **(Run-0) the F-human branch of `run_study_file`, no mesh:** inside `pytest.MonkeyPatch.context()`, `fem_em_solver.study._study_solve_tuned_quadrature` is replaced by a counting fake that sets and returns `study.result = SimpleNamespace(cells=369387, s_node=<a fixed 2×2 complex array>, y_in={"ccw": <2 complex>, "cw": <2 complex>}, feeds=("P17", "P21"), az_feeds=(11.25, 101.25), az_mirror=56.25, gap_impedance_ohm=…, feed_impedance_ohm=…, source_voltage_v=…)`; `Study._solve_port_fields` is patched to raise if called. `run_study_file` on an `outputs`-free copy of `FH_FILE` (rank 0 writes it under `LOGS / "f_human_run0"`, then `Barrier`): fake called exactly once per rank; `summary["family"] == "f_human_birdcage"`; `s_node` / `y_in` / `az_*` / `cells` / the `excitation` block equal the fake's values exactly; no `"touchstone"` in `summary` or `paths`; `summary["conductor"]` present; on rank 0 the written JSON equals the in-memory summary exactly. (The fake is the monkeypatch target item 383's journal names: `Study.sar` on the tuned drive reaches the solve through `_tuned_result`, bypassing `Study.solve`.)

**Negative controls (asserted; deterministic parse, so no factor to back — rule (e) has nothing to label).**
- **(c1) dropped key:** `FH_FILE`'s mapping minus `excitation.esr_ohm` parses, its `_spec_diff` against `_study()` is exactly `['excitation.esr_ohm']`, and its `esr_ohm` is `TwoPortQuadratureDrive`'s own default by hex — a loader that silently ignored `esr_ohm` passes (X) and fails (K) or (c1).
- **(c2) one ULP:** `c_in_f` set to `repr(np.nextafter(_c_f(), 1.0))` parses and diffs exactly `['excitation.c_in_f']` — a loader that converts through pF or rounds would pass a looser check and fail here.

**Windows** (foreground, Bash `timeout: 660000`; each a bare whole harness call, repo-relative, nothing chained or piped; `-s`; node IDs, never `-k`; no edit to `study_yaml.py`, the test module, the YAML files or the registry after their last green window — rule (i)). Before RG2, RG4, RG5 and RG8: `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reads 0.
- **W00, the literals (no solve, ≈ 5 s; run once the test module exists, before the YAML is written):** `scripts/testing/run_and_log.sh WF-1-s7a-W00-literals "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace timeout -k 30 60 python3 -c \"from tests.package.test_wf1_yaml_f_human import _literal_lines; print(_literal_lines())\"'"` — paste its output into the YAML verbatim. If the guard refuses the quoting, write the file from (K-c)'s failure print in W0 (the node prints each expected `repr` before asserting) and re-run W0 — one extra ≈ 40 s window.
- **W0 (standard, `-n 2`):** `scripts/testing/run_and_log.sh WF-1-s7a-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_f_human.py -s -v --tb=short'"` (predicted 20–45 s).
- **Standing rule (c), named in full.** `study_yaml.py` is reached by `test_wf1_study_yaml.py`, `test_wf1_yaml_thermal.py`, `test_wf1_yaml_copper.py`, `test_wf1_yaml_sar_b1.py`, `test_wf1_yaml_probes.py`, `test_wf2_bundle_reproducible.py` (`run_study_file`'s bundle, landed after items 343–345, so new to this list), `test_wf1_installed_cli.py` / `test_wf1_installed_solve.py`, `test_wf13_launch_helper.py` / `test_wf13_cli_estimate.py`, and the examples `mri:9` / `mri:13`. The edits touch `_parse_geometry`, `_parse_excitation`, `_parse_conductor`, `_parse_thermal`, `_parse_outputs`, `_parse_mapping` and `run_study_file`, every one on every route — the 10-08 ruling (5) applies ("keep every re-gate window"). Proposed in-slot set, each command **verbatim** from its last green log under a new label:
  - **RG1** `WF-1-s7a-RG1-wf1-w1` — `20261009T020327Z_WF-1-s6-W2-wf1-w1.log:12`: `scripts/testing/run_and_log.sh WF-1-s7a-RG1-wf1-w1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/package/test_wf1_study_yaml.py::test_wf1_good_files_build_the_python_studies tests/package/test_wf1_study_yaml.py::test_wf1_malformed_files_name_their_field tests/package/test_wf1_study_yaml.py::test_wf1_cli_sphere_reproduces_wf15 -s -v --tb=short'"` (102 s, 15 passed).
  - **RG2** `WF-1-s7a-RG2-wf1-w2` — `20261009T020513Z_WF-1-s6-W3-wf1-w2.log:12`: `scripts/testing/run_and_log.sh WF-1-s7a-RG2-wf1-w2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 400 mpiexec -n 2 python3 -m pytest tests/environment tests/package/test_wf1_study_yaml.py::test_wf1_cli_f_small_reproduces_wf16 -s -v --tb=short'"` (142 s, 13 passed).
  - **RG3** `WF-1-s7a-RG3-thermal` — `20261009T020735Z_WF-1-s6-W4-thermal.log:12`: `scripts/testing/run_and_log.sh WF-1-s7a-RG3-thermal "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_thermal.py -s -v --tb=short'"` (33 s, 6 passed).
  - **RG4** `WF-1-s7a-RG4-copper` — `20261009T020810Z_WF-1-s6-W5-copper.log:12`: `scripts/testing/run_and_log.sh WF-1-s7a-RG4-copper "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 300 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_copper.py -s -v --tb=short'"` (95 s, 11 passed).
  - **RG5** `WF-1-s7a-RG5-sar-b1` — `20261009T020945Z_WF-1-s6-W6-sar-b1.log:12`: `scripts/testing/run_and_log.sh WF-1-s7a-RG5-sar-b1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 400 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_sar_b1.py -s -v --tb=short'"` (122 s, 13 passed).
  - **RG6** `WF-1-s7a-RG6-probes` — `20261009T020231Z_WF-1-s6-W1-probes.log:12`: `scripts/testing/run_and_log.sh WF-1-s7a-RG6-probes "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_probes.py -s -v --tb=short'"` (40 s, 9 passed).
  - **RG7** `WF-1-s7a-RG7-wf2-sphere` — `20261010T024009Z_OPS-70-W24a-wf2-sphere.log:12`: `scripts/testing/run_and_log.sh WF-1-s7a-RG7-wf2-sphere "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf2_bundle_reproducible.py::test_wf2_comparator_flags_planted_changes tests/package/test_wf2_bundle_reproducible.py::test_wf2_sphere_bundle_reproduces -s -v --tb=short'"` (46 s, 6 passed).
  - **RG8** `WF-1-s7a-RG8-wf2-four-port` — `20261010T024102Z_OPS-70-W24b-wf2-four-port.log:12`: `scripts/testing/run_and_log.sh WF-1-s7a-RG8-wf2-four-port "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 400 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf2_bundle_reproducible.py::test_wf2_four_port_bundle_reproduces -s -v --tb=short'"` (73 s, 5 passed).
  - **Not in the slot (drafter's proposal; the review rules):** `test_wf1_installed_solve.py`, `test_wf13_launch_helper.py`, `test_wf13_cli_estimate.py` and `mri:9` / `mri:13` reach only the parse and the sphere / four-port run paths that RG1–RG8 cover node for node; adding them is ≈ 4 min more and would push 390a past the box — if required, the review queues them as a serial lane-1 re-gate (item 360's precedent).
- **Wgen / Wapi / Wdoc (the new envelope nodes move `load_study.md` / `run_study_file.md`; commit the regenerated pages):** `scripts/testing/run_and_log.sh WF-1-s7a-Wgen "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 python3 scripts/docs/generate_api_reference.py --out docs/guide/api'"`, then `scripts/testing/run_and_log.sh WF-1-s7a-Wapi "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 python3 -m pytest tests/package/test_wf17_api_reference.py tests/package/test_wf14_validation_envelope.py -v -s --tb=short'"` (12 passed or more), then `scripts/testing/run_and_log.sh WF-1-s7a-Wdoc "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && timeout -k 30 60 python3 scripts/testing/check_example_doc_references.py --docs-root docs/guide --check-links'"` (`0 dead`) — commands as `20261011T214232Z` / `214234Z` / `214244Z_WF-16-s9a-W*.log:12`.
- **W9, installed `femem check`:** `scripts/testing/run_and_log.sh WF-1-s7a-W9-installed-femem "docker compose exec -T fem-em-solver bash -lc 'timeout -k 10 50 bash /workspace/scripts/package/wf13_setup.sh real && mkdir -p /tmp/wf13/run/tests/package/wf1 && cp /workspace/tests/package/test_wf1_installed_cli.py /tmp/wf13/run/tests/package/ && cp /workspace/tests/package/wf1/*.yaml /tmp/wf13/run/tests/package/wf1/ && cd /tmp/wf13/run && env -u PYTHONPATH WF13_INSTALLED_CHECK=1 timeout -k 10 90 /tmp/wf13/venv/bin/python -m pytest -p no:cacheprovider tests/package/test_wf1_installed_cli.py -v -s --tb=short'"` (16 s; `20261009T021200Z_WF-1-s6-W9-installed-femem.log:12`; the new files ride along).

**Traps already paid for.**
- PyYAML 1.1 floats need a dot, and a signed exponent if any (`WF-1` step 1): `1.1824370000000001e-10` and `58000000.0` are fine, but a `repr` with no dot (`5e-05`, `1e-17`) is read as a **string** — write such a value as `%.17e` (`5.00000000000000000e-05`); (K-c) checks the bits. YAML `true` is `bool`, `1` is `int`.
- **Never fill a default in the loader** (the (Deg) rule): `esr_ohm`, `feeds`, `z0_ohm` and `study.degree` omitted mean the dataclass / `Study` default. An F-human file that omits `degree` builds a **degree-2** F-human Study — the tracked file carries `degree: 1`, and (K) catches its absence.
- **Never call `_study_solve_tuned_quadrature`, `build_tuned_quadrature` or any mesh in 390a.** The fake is the patch target for the tuned solve (item 383's 2026-10-11 journal hypothesis); patching `Study.solve` alone would not intercept `Study.sar`.
- `TwoPortQuadratureDrive.feeds` is coerced to a tuple in `__post_init__`; the parser passes a tuple anyway so `_cmp`'s `type(x) is type(y)` holds.
- A new field on an existing frozen dataclass broke a field-set pin once (`WF-16` step 7, first W0 red): this item adds **no** field to any `study.py` dataclass, and `StudyRun` keeps its fields.
- `cell_tags.values` / rc are rank-local — `_rc_all` reduces min / max. On status 124 or a kill: `pgrep -c python3` first; at 0 only, the rubric's orphan-stem sweep as one bare call labelled `WF-1-s7a-cache-sweep`. Never raise a timeout or change `-n`. Write ≈, never a tilde. No AED number anywhere.

**Scope.** Covered: the parse of `f_human_birdcage` (hole route only) + `two_port_quadrature` + `conductor:` into item 382's `Study` bitwise, the literal SAR lattice and points, the refusals, and the F-human branch of `run_study_file` on a fake result. **Not claimed:** any solve (390b); the solid σ = 800 S/m F-human or its 32-port `lumped_port_drive` through YAML (refused); a Touchstone of the node S; `femem run --estimate` / `--ranks` or the installed console script on the F-human file; `capacitor_sheets` (below); a stable schema.

**What a negative result means.** Report and stop; park on `attempt/WF-1-step7a-<UTC>`, mark item 390a BLOCKED with its unblock condition in the same commit (rule (d)), journal in `docs/testing/attempts.md`. (K) / (K-c) red ⇒ a parser or transcription defect (fix in-slot only inside the new lines or the YAML file, then re-run W0). (X) red ⇒ a refusal missing or mis-named (same). (Run-0) red ⇒ the new branch (same). Any RG red ⇒ the edit leaked into an existing route: revert `study_yaml.py`, BLOCKED, known-issues entry if the cause is not the diff. No `PREDICTED RED:` line.

**Status it can move.** `WF-1` stays ✅ (rule (j): its Done-when is step 1's). On green: the `WF-1` row's history gains a step-7a clause (the F-human schema, parse-only), `WF-14`'s envelope gains two nodes, and 390b becomes runnable; §10 U43 is **not** done until 390b lands.

**Executor spawn prompt (verbatim).** "Execute §9 item 390a — `WF-1` step 7a (§10 U43 part 1, provisional) — per `docs/planning/chunks/WF-1.md` §9 item 390a. Read the brief and any review rulings beneath it first; the rulings hold. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Pass timeout: 660000 on every Bash call that runs a window; run pytest with -s. Run W00, W0, RG1–RG8, Wgen, Wapi, Wdoc, W9 in order, foreground. Before RG2, RG4, RG5 and RG8, run `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` and require 0. If W0 is not green by minute 38, park. No mesh and no solve in the new test module; never call `build_tuned_quadrature` or `_study_solve_tuned_quadrature`; never import a `test_*` name; never add a field to a `study.py` dataclass. No edit to `study_yaml.py`, the test module, the YAML files or the registry after the window that closes on them. Never delete a log. You are the executor — do not spawn agents, never background a harness window. If the docker socket is denied, stop and report the exact command string. On a red, park on `attempt/WF-1-step7a-*`, mark item 390a BLOCKED in the same commit, and report."

### §9 item 390b — 2026-10-11 (chain U43 part 2, daily provisional under §9 rule (6); chain step (solve))

**Outcome 2026-10-12 (implementer, slot 12:30 UTC): DONE, green; no `src/` edit.** W0 `20261012T125618Z_WF-1-s7b-W0.log` 5 passed / 5 skipped (knob unset), 4 s. `pgrep -c python3` read 0, then W1 `20261012T125631Z_WF-1-s7b-W1.log` 10 passed on each of 8 ranks, 243 s (heavy; `run_study_file` 238.6 s, build 113.8 s, node pair 26.6 s): (A0) cells 369 387 (rel 0), solve count 1 at every node, `run.result is run.study.result`; (Δc) 1.735e-18 m over 19 (printed); (S_rep) 4.287849e-12 (prior 4.265317e-12); (Y) 2.217e-12 ccw / 3.216e-12 cw; (S10) ccw lattice max 2.857e-10, c\* = `r=45 dz=+0 az17+0`, SAR10g(c\*) 4.852895959e-03 W/kg (rel 4.310e-11, ball 10.010236 g), cw(M c\*) rel 4.714e-11; (iii-a) 0.6413 %; (B) 41 of 41 valid, rel 1.101e-07; (Sense) 21.5155; (M_B) 1.9647e-03; (B2′) 0.000e+00; (J) rank 0 differing keys `[]` (ranks 1–7 `['elapsed_s']` only), no Touchstone; controls (N1) 9.5476e-01 = 19.1× and (N2) 4.266875e-02 = 4.27e+04× the band, both asserted and green. M(c\*) is bitwise lattice row 5 (`r=45 dz=+0 az17+90`; the 56.25° mirror maps 11.25° onto 101.25°), so index 18's cw read is that row's cw value (cw[18] − cw[5] = 0, printed) and still the cw value at M(c\*). Printed only: the z = z_c twelve-point mirror cw(Mx) vs ccw(x), rel L2 1.4734 %. RGa (03:00 ruling Q6, 390a's module re-run after the registry edit) `20261012T130057Z_WF-1-s7b-RGa.log` 6 passed, 3 s. Pages (ruling Q3): Wgen 44 files 2 s (`…130101Z…Wgen`), Wval 166 nodes 2 s (`…130103Z…Wval`), Wapi 22 passed 15 s (`…130105Z…Wapi`), Wdoc 0 dead 1 s (`…130121Z…Wdoc`).

*Drafted by a general-purpose agent for the 2026-10-11 18:00 daily review; the review reviews it.* Licence and tag as item 390a. **No band, tolerance, Done-when, record or validation target is set or moved; every band and record is imported by name.** **Depends on item 390a landing on `main` green. If it did not, skip this item.** No `src/` edit is in scope (below).

**`WF-1` step 7b: `run_study_file("f_human_tuned_copper.yaml")` at `-n 8` reproduces items 382 / 382b's records from its own summary.**
- **Executor:** implementer. **Adds:** `tests/package/test_wf1_yaml_f_human_solve.py` only (plus envelope lines and regenerated pages). Opt-in knob `FEM_EM_WF1_S7B=1` (unset ⇒ the solve nodes skip, stated; any other value raises — item 382b's `_knob` pattern); width asserted `MPI.COMM_WORLD.size == RECORD_RANK_WIDTH` (8) when set; `@complex_only` on the solve nodes. **One module-scoped fixture** runs, inside `pytest.MonkeyPatch.context()`, a **counting call-through wrapper** around `fem_em_solver.study._study_solve_tuned_quadrature`, then `run = run_study_file(FH_FILE, output_dir=LOGS / "f_human_tuned", comm=MPI.COMM_WORLD)`; the count is read after every test-side read and asserted `== 1` (one solve for the file's S, SAR and B — `_tuned_result` reuses `study.result`). Nodes: `test_wf1_fh_solve_knob_and_imports` (no mesh), `test_wf1_fh_yaml_cells_centres_once`, `test_wf1_fh_yaml_node_reproduces_h9p`, `test_wf1_fh_yaml_sar_reproduces_h9p`, `test_wf1_fh_yaml_b1_reproduces_h9p`, `test_wf1_fh_yaml_summary_json`. Envelope: the five solve nodes under `[names.run_study_file]`; the node / SAR / B nodes also under `[names.TwoPortQuadratureDrive]`; the cells node under `[names.FHumanBirdcageGeometry]` (lines added only).
- **Inputs (gitignored, checked present on the host at drafting):** `logs/wf10/h9p-node-118.2437.json` (5 474 bytes, mtime 2026-10-11 20:09:20 UTC — item 381's rewrite, against which item 382 read (S_rep) 4.265317e-12) and `logs/wf10/h17-node-118.2000.json`. **The node JSON absent ⇒ stop, BLOCKED; never re-run 358 or call `build_tuned_quadrature`** (it rewrites the JSON). Read on rank 0 and broadcast, exactly as item 382's `tuned` fixture does (`tests/package/test_wf16_study_tuned_f_human.py:186–210`; copy those lines with a provenance comment — they are inside a fixture, not an importable helper).
- **Tier / ranks:** heavy, `-n 8` (the records' width), complex, container `timeout -k 30 590` (item 382b's ceiling). **Class:** chain step (solve). **Cost:** ≈ 12 measured min (class median); **≈ 42 predicted**: 15 fixed + ≈ 15 code + W0 ≈ 0.5 + W1 ≈ 4.5–5.5 + Wgen / Wapi / Wdoc ≈ 0.3 + ≈ 5 spare. **W1 predicted 230–330 s** (≤ 435 s at the 1.32× box reading): solve + two CG1 projections 146.5 s (build 114.1 s, node pair 27.5 s; `20261011T214947Z_WF-16-s9b-W1.log:6385`), 38 balls (19 centres × 2 senses) at 69.2 s / 36 ≈ 1.92 s each ≈ 73 s, interpreter start + `test_complex_mode.py` + the test-side reads on the cached CG1 ≈ 25 s. **If W1 has not started by minute 40, park.**

**Anchors (asserted; all imported).** From `tests/package/test_wf16_study_tuned_f_human_b1_sar.py` (registered there for downstream import): `H9P_MEAN_B1_PLUS_CCW_RECORD_T`, `H9P_SAR10G_TABLE_RECORD`, `H9P_SAR10G_CSTAR_RECORD_W_PER_KG`, `H9P_CSTAR_LABEL`, `H9P_N_SAMPLE_POINTS`, `RECORD_RANK_WIDTH`; from `tests/package/test_wf16_study_tuned_f_human.py`: `_cx`, `H17_NODE_JSON`, `NODE_JSON` (as it imports it); from 390a's module: `FH_FILE`, `_expected_centres`; `REUSE_REPRODUCTION_RTOL` (`tests/validation/test_port_birdcage_ring_matrix.py:171`, = 1e-6), `_residual` (`tests/validation/test_port_circuit_layer_field.py:142`), `C4_COVARIANCE_BAND` / `_relative_l2` / `_sample_points` (`tests/validation/test_birdcage_b1_plus_map.py:120, :600, :246`), `_mirror_xy` / `_read_senses` (`tests/validation/test_birdcage_b1_quadrature.py:203, :214`), `_lattice` (`test_wf10_f_human_tuned_sar10g.py:58`), `_phantom_geometry` (`tests/validation/test_wf7_f_human_sar10g_hotspot.py:95`), `PHANTOM_CELL_TAG`, `HOLE_CELL_RECORD`, `CELL_COUNT_BAND`, `EXACT_IDENTITY_RTOL` (`tests/validation/test_birdcage_sar_mass_averaged.py:163`), `LOGS` (`test_wf1_study_yaml.py:86`).
- **(A0) cells and one solve:** `summary["cells"] == HOLE_CELL_RECORD` (369 387) within `CELL_COUNT_BAND` and `== run.result.cells`; the solve count `== 1`; `run.result is run.study.result`. **(Δc), printed only:** max \|Δ\| between the YAML's 19 centres and the run-derived ones — `_lattice(float(run.result.az_feeds[0]), z_c)` with z_c from `_phantom_geometry(run.result.mesh, run.result.cell_tags, PHANTOM_CELL_TAG, comm)`, and `_mirror_xy(·, run.result.az_mirror)` for M(c\*) — **predicted ≤ 1e-15 m** (382's feed azimuths agree with 11.25° / 101.25° to 12 printed decimals; never asserted).
- **(S_rep) the node S from the summary:** `_residual(np.array(summary["s_node"] as complex), s_json) ≤ REUSE_REPRODUCTION_RTOL`; prior on the same Study route and file 4.265317e-12 (`20261011T213654Z_WF-16-s9a-W1.log`, the (S_rep) lines `:6381–6386`); predicted ≤ 1e-11 (390a's (K) makes this item 382's Study). **(Y)** per sense, ‖y_in − y_json‖ / ‖y_json‖ from the summary ≤ `REUSE_REPRODUCTION_RTOL` (priors 2.358e-12 ccw / 3.088e-12 cw, `:6387–6388`).
- **(S10) the 10 g SAR from the summary:** `summary["sar"]["values"]["ccw"]` at the 18 lattice centres (indices 0–17, labels from `_lattice`'s order) against the 18 ccw records, each ≤ `REUSE_REPRODUCTION_RTOL` (382b max 2.857e-10, `20261011T214947Z_WF-16-s9b-W1.log:6383`); the argmax over those 18 is `H9P_CSTAR_LABEL`; SAR10g(c\*) against `H9P_SAR10G_CSTAR_RECORD_W_PER_KG` ≤ the band (382b 4.306e-11); `summary["sar"]["values"]["cw"]` at index 18 (M(c\*)) against the c\* row's cw(M) record (4.821772251e-03) ≤ the band (382b 4.70e-11, `:6367`). **(iii-a)** the mirror at c\*, \|cw(M c\*) − ccw(c\*)\| / ccw(c\*) from the summary, ≤ `C4_COVARIANCE_BAND` (record 0.6413 %, `:6384`). Printed, never asserted: the z = z_c row's twelve-point mirror (on the lattice, the feed-bisector reflection maps az₁₇+d to az₁₇+90−d, so the twelve z = z_c points are mirror images of one another to the reflection's round-off — the dz = ±50 mm rows' images are not in the file, so (iii-b) is not claimed), the ccw value at M(c\*) and the cw values at the 18 lattice points; `summary["sar"]["radius_m"]` and each ball's `mass_kg` (382b printed 10.010236 g at c\*).
- **(B) mean \|B₁⁺\|_ccw on the YAML-solved Study:** `points = _sample_points({"mesh": run.result.mesh, "cell_tags": run.result.cell_tags})`; `run.study.b1_plus_at(points, drives=("ccw",))` (the cached CG1 — no projection, no solve); 41 of 41 valid (`H9P_N_SAMPLE_POINTS`); the mean against `H9P_MEAN_B1_PLUS_CCW_RECORD_T` ≤ `REUSE_REPRODUCTION_RTOL` (382b 1.101e-07, of which the `.6e` print contributes ≤ 1.7e-7, `:6361`). **(Sense)** mean\|B₁⁺\|_ccw / mean\|B₁⁻\|_ccw > 1 through `_read_senses(run.study.b1_plus()["ccw"], points)` (record 21.5155, `:6362`). **(M_B)** `_relative_l2(|B₁⁻|_cw(Mx), |B₁⁺|_ccw(x), mask)` ≤ `C4_COVARIANCE_BAND`, M = `_mirror_xy(·, run.result.az_mirror)` (record 1.9647e-03, `:6363`). **(B2′) summary = API:** `summary["b1_plus"]["values"]` at the file's two points equals `run.study.b1_plus_at(those points)` per sense, max \|Δ\| / max \|B₁⁺\| ≤ `EXACT_IDENTITY_RTOL` (item 344 read 0.0 on its fixture), `valid` all true. *Why (B) is read on the run's Study, not from the summary:* 358's 41 sample points are tag-3 cell centroids of the built mesh (`test_birdcage_b1_plus_map.py:246–280`) — transcribing them into a tracked file would make the file a mesh record re-recordable under every image change; the drafter flags this for the review (a `points: phantom_sample` YAML convenience would be new API).
- **(J) the written bundle:** on rank 0 the JSON file equals `run.summary` exactly (other ranks print their differing keys; WF-2 measured `['elapsed_s']` only); `"touchstone"` absent from `summary` and `run.paths`; `summary["family"] == "f_human_birdcage"`; `summary["excitation"]["c_in_f"].hex() == _c_f().hex()` and `summary["conductor"]["sigma_s_per_m"] == 5.8e7`.

**Negative controls (rule (e); both backed by a prior of the same comparison on the same fixture — 390a's (K) makes the YAML Study item 382's).**
- **(N1) mis-paired sense — asserted, backed:** `_relative_l2(|B₁⁺|_cw(Mx), |B₁⁺|_ccw(x), mask) > C4_COVARIANCE_BAND`; 382b read 9.5476e-01 = 19.1× the band (`20261011T214947Z_WF-16-s9b-W1.log:6364`). The factor is printed, not asserted; the ceiling (a relative L2 of two non-negative maps) does not bind a > 1× bar.
- **(N2) the wrong node — asserted, backed:** the summary's `s_node` against `H17_NODE_JSON`'s (118.2000 pF) misses by > `REUSE_REPRODUCTION_RTOL`; 382 read 4.266875e-02 on this comparison (`20261011T213654Z_WF-16-s9a-W1.log:6389`, there labelled "predicted, printed only" because no prior existed; that reading is now the prior). If `H17_NODE_JSON` is absent, print "not run" — (N1) carries the rule (e) control.

**Windows** (foreground, Bash `timeout: 660000`; bare whole harness calls, repo-relative; `-s`; never `-k`; under 10 min each, so no durable capture). Before W1: `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reads 0; never beside another heavy window (item 383's `-n 8` window included).
- **W0 (no mesh, knob unset; ≈ 20 s):** `scripts/testing/run_and_log.sh WF-1-s7b-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_f_human_solve.py -s -v --tb=short'"` (the solve nodes skip, stated).
- **W1 (heavy, `-n 8`; predicted 230–330 s):** `scripts/testing/run_and_log.sh WF-1-s7b-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_WF1_S7B=1 timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_f_human_solve.py -s -v --tb=short'"`.
- **Wgen / Wapi / Wdoc** — 390a's three commands verbatim under labels `WF-1-s7b-Wgen`, `WF-1-s7b-Wapi`, `WF-1-s7b-Wdoc` (the new nodes move `run_study_file.md`, `TwoPortQuadratureDrive.md`, `FHumanBirdcageGeometry.md`; commit the pages).
- **No rule (c) re-gate:** no `src/` file changes. If W1 shows a defect in the F-human branch of `run_study_file`, a fix confined to that branch's body (lines no pre-existing route executes) is allowed in-slot **only if** 390a's W0 and RG1–RG8 can all run again before minute 50 after W1's re-run; otherwise park.

**Traps already paid for.**
- **One solve.** `run_study_file` calls `study.solve` first, so `study.sar` / `study.b1_plus` reuse `study.result` (`_tuned_result`); a test-side `Study` built afresh, or a cleared `result`, would re-solve (≈ 150 s) through `_tuned_result` without passing `Study.solve` — the counting wrapper is on `_study_solve_tuned_quadrature` for that reason (item 383's journal).
- Every `mass_averaged_sar` call is collective on every rank (item 227's deadlock); `_sar_summary` already calls it on every rank — never filter centres per rank. Point evaluation only through `evaluate_vector_field_parallel` (inside `b1_plus_at`); `|B₁±|` from the evaluated vector, never a projected scalar; `points` identical on every rank.
- `SAR_QUADRATURE_DEGREE` (16) stays `Study.sar`'s default — the loader passes none (an unpinned `SpatialCoordinate` integral can compile > 9 min, `POST-5`); the ball is a UFL conditional — compare squares, never `sqrt` (`WF-6` 4c).
- Never import a `test_*` name or call `build_tuned_quadrature`; κ off everywhere; the node `s_node` is not a 50 Ω S (item 133) — no Touchstone. A record row off its last printed digit with every assertion green is reported, not a stop (382b's cw(M) row at r = 45 az₁₇+90, rel 1.24e-10).
- Status 124 or a kill: `pgrep -c python3` first; at 0 only, the orphan-stem sweep labelled `WF-1-s7b-cache-sweep`. Never raise the timeout or change `-n`. Rule (i). Write ≈, never a tilde. No AED number.

**Scope.** Covered: item 382 / 382b's tuned copper F-human (degree 1, 64 MHz, κ off, 369 387 cells, r = 0.2 Ω on all 32 gaps, the two-port 90° feed at C_in = 118.2437 pF) through `femem run`'s Python entry point at `-n 8`: node S, active `Y_in`, the 18-point ccw 10 g lattice with c\*, cw at M(c\*), the c\* mirror, and mean \|B₁⁺\|_ccw with its sense and mirror gates. **Not claimed:** the console script or `--ranks` on this file; (iii-b) over the full lattice; the 18 cw(M) rows; homogeneity, an absolute `|B₁⁺|` or SAR, a global SAR maximum (c\* is a lattice argmax), compliance; `f₀`, Q, a match or any tuning figure (C_in is an input record); degree 2; any AED comparison.

**What a negative result means.** Report and stop; park on `attempt/WF-1-step7b-<UTC>`, mark item 390b BLOCKED (rule (d)), journal. (A0) count > 1 ⇒ the run solves twice — a finding against the run branch (known-issues), not a re-band. (S_rep) / (Y) red with (A0) green ⇒ the YAML Study is not item 382's despite (K) — print the per-entry residuals; known-issues; never re-record. (S10) / (B) red with the mirror gates green ⇒ print (Δc) first (a transcription off by more than round-off is a 390a defect), then the per-row readings beside the records. A mirror gate red ⇒ a finding (known-issues; never re-band). A control inside the band ⇒ a blind comparator: stop. No `PREDICTED RED:` line.

**Status it can move.** `WF-1` stays ✅ (rule (j)). On green of both parts: the `WF-1` §7 row gains one **step-7 record** clause (the tuned copper F-human through YAML at `-n 8`, readings and logs); §10 U43 is marked done (provisional until the 10-12 weekly); `WF-14`'s envelope gains nodes; the Phase 9 (9a) clause "`femem run` drives the same API" covers the tuned copper F-human.

**Executor spawn prompt (verbatim).** "Execute §9 item 390b — `WF-1` step 7b (§10 U43 part 2, provisional) — per `docs/planning/chunks/WF-1.md` §9 item 390b. Skip it unless item 390a is DONE on `main`. Read the brief and any review rulings beneath it first; the rulings hold. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Pass timeout: 660000 on every Bash call that runs a window; run pytest with -s. Run W0, W1, Wgen, Wapi, Wdoc in order, foreground. Before W1, run `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` and require 0. If W1 has not started by minute 40, park. If `logs/wf10/h9p-node-118.2437.json` is absent, stop. Never call `build_tuned_quadrature`, never import a `test_*` name, never edit `src/` except as the brief's W1 clause allows. No edit to the test module or the registry after the window that closes on them. Never delete a log. You are the executor — do not spawn agents, never background a harness window. If the docker socket is denied, stop and report the exact command string. On a red, park on `attempt/WF-1-step7b-*`, mark item 390b BLOCKED in the same commit, and report."

**`capacitor_sheets` in YAML — separable, not in item 390 (for the review).** `TwoPortQuadratureDrive` carries its own ESR and never uses `CapacitorSheets`, so U43's F-human target does not need it. Adding it is a different change: `CapacitorStudyResult` (item 373) has no Touchstone (`write_touchstone` raises) and no `z_matrix`, so it would need a sub-branch **inside** `run_study_file`'s four-port branch — a path every F-small YAML route executes — and its comparand is not a registered constant (item 373's (I) compares against an in-process `reduce_terminated_ports(S0, …)`, so the YAML test would solve a second Study). If the review wants it, a fourth provisional step could read: `excitation: {kind: capacitor_sheets, c_f, esr_ohm, on: [P2]}` on `f_small_copper.yaml`'s geometry and conductor; the summary's reduced 3×3 against `reduce_terminated_ports(S0, z0, {2: _rc(ESR_OHM)})` at `REDUCTION_BAND` (item 373 (I) read 3.036163e-05 on P2), S0 an in-test `LumpedPortDrive` Study; planted ESR (`PLANT_OHM`) beyond the band (backed, 1.095344e-02); standard `-n 2`, two Studies ≈ 60–90 s; the same rule (c) re-gate set as 390a.

**Review ruling 2026-10-11 18:00 (daily) on items 390a / 390b.** (1) The split is accepted; 390a carries every `src/` change and the re-gates, 390b is test-only and runs after 390a lands. (2) The \|B₁⁺\| record is read on the run's Study, not written into the tracked YAML or summary (358's 41 points are mesh cell centroids; a tracked copy would be a mesh record). (3) The re-gate list stands as drafted (RG1–RG8, `WF-2` included); the installed-solve, `WF-13` and example consumers are not re-gated by this item because they do not import the edited parser branches — the executor confirms that by grep before W0 and discloses it. (4) `capacitor_sheets` in YAML is not part of U43. (5) 390b never runs beside another heavy window (383 and 387 precede it in §9). Provisional under §9 rule (6), for the 10-12 weekly.

### §9 item 394 — 2026-10-12 (chain U44, daily provisional under §9 rule (6); lane 1 / chain step (solve))

**Outcome 2026-10-12 (implementer, slot 14:00 UTC): DONE, green; no `src/` edit (`git diff --stat main -- src/` empty).** W0 `20261012T140246Z_WF-1-s8-W0.log` 5 passed, 5 s: (K) both diffs `[]`, outputs parse as written (`mass_g` hex `0x1.4000000000000p+3`); (X=) YAML `StudyConfigError` at `conductor` naming `as_hole`, cli run / check rc (2, 2), nothing written, 0 build calls; control `as_hole: true` diff `[]`, outputs deleted still refused at `conductor`; Python `sar` / `b1_plus_at` / `solve()` refuse with `as_hole` (stub build called 3 times, each returning no mesh — printed). `pgrep -c python3` read 0, then W1 `20261012T140300Z_WF-1-s8-W1.log` 11 passed, 151 s (`run_study_file` 67.9 s = two builds; Study A 39.5 s, Study B 32.9 s; API re-read 6.7 s reusing `fields_result`): (K1) 80 798 cells, one conductor term; (S=Y) 6.442e-14; (R) 4.093e-11 / 8.673e-11 / 1.376e-10; (I-SAR) max 4.902e-13 over four drives and every nonzero key, (iii) 0.000e+00; (I-B) max 5.397e-14, all valid, (B2′) 0.000e+00; (P) residual 8.501e-13, P_surf 1.516873149e-06 W; controls (N1) SAR 4.149e-04 = 414.9×, `|B₁⁺|` 3.576e-04 = 357.6× the band, (N2) 3.1559e-04 = 315.6×, all asserted, green. Printed only: SAR10g(P1, origin) 2.408073397e-06 W/kg and `|B₁⁺|`(P1, origin) 1.833178311e-08 T, equal to item 388's at every printed digit. Pages (ruling Q3): Wgen 44 files 1 s (`…140547Z…Wgen`), Wval 174 nodes 2 s (`…140552Z…Wval`), Wapi 22 passed 14 s (`…140602Z…Wapi`), Wdoc 0 dead 1 s (`…140622Z…Wdoc`). RGa (03:00 ruling Q6, registry shared with 390a / 390b) `20261012T140638Z_WF-1-s8-RGa.log` 6 s: 390a's module 6 passed at `-n 2`, 390b's module 5 passed / 5 skipped with the knob unset (no mesh); 390b's heavy `-n 8` W1 was **not** re-run — no solve code moved, registry lines only.

*Drafted by a general-purpose agent for the 2026-10-12 03:00 daily review; the review reviews it.* Licence: §9 rule (6), a provisional step of §10 chain U ("each step reproduces a gated result or asserts a gated identity through the API, within its unmoved band"), tagged "U44 (daily, 2026-10-12 — provisional)"; the next weekly ratifies or reverts. Origin: item 388's own scope line ("the YAML `conductor` + `outputs.sar` / `outputs.b1_plus` combination (reachable after the fix, untested — a `WF-1` follow-up for the review)", `WF-16.md`, item 388, "Scope"). **Target: reproduction only. No band, tolerance, Done-when, record or validation target is set or moved; every band is imported by name.** Owner: `WF-1` (steps 1–6 = items 256 / 293 / 295 / 343 / 344 / 345; step 7 = items 390a / 390b, queued ahead). The 2026-10-08 18:00 rulings on items 343–345 (this file) hold, in particular (4): the two solves of an `outputs` file are accepted and disclosed.

**`WF-1` step 8: `run_study_file` on the F-small copper hole Study with `outputs.sar` / `outputs.b1_plus` returns the 10 g SAR and `|B₁⁺|` of the Python copper Study's field-keeping route, and a solid-route `leontovich` file refuses as the Python `Study` refuses.**
- **Executor:** implementer. **Test-only: no `src/` edit.** Adds one test module, two tracked YAML files, registry lines (added only) and the regenerated pages. No public name, no docstring of a public object, no `__init__` / `public_surface.txt` change.
- **Tier:** heavy by ceiling (container `timeout -k 30 400`, the ceiling of item 388's W1 and of the `WF-1` four-port windows), `-n 2` throughout, complex. **Class:** chain step (solve); class median **13 measured** min. **Cost: ≈ 38 predicted min** — 15 fixed + ≈ 12 code (module ≈ 200 lines, two YAML files, registry lines) + W0 ≈ 0.7 + W1 ≈ 3 + Wgen / Wval / Wapi / Wdoc ≈ 0.7 + ≈ 4 row / journal + ≈ 3 spare. **If W1 has not started by minute 35, park** (rule (d)).
- **Depends on:** nothing open. Item 388's fix is on `main` (`39fd8a01`). Items 390a / 390b (ahead in §9) edit `study_yaml.py`'s F-human branches and add registry lines; this item edits no `src/` file and shares only `tests/package/validation_envelope.toml` and the generated pages with them — on a rebase keep both sides of the registry and re-run Wgen / Wval / Wapi / Wdoc after it. If 390a lands **after** this item, 390a's rule (c) list gains this module's W1 (it reaches `_parse_conductor`, `_parse_outputs` and `run_study_file`; ≈ 3 min) — for the review.
- **Family:** `WF-1` step 8, attempt 1 of 4.

**What the code says (read by this drafter on `1b8c3532`, not run).**
- **No schema addition is needed.** `_parse_conductor` (`src/fem_em_solver/study_yaml.py:333–357`) admits `kind: leontovich` on a `four_port_birdcage` with `geometry.as_hole: true`; `_parse_outputs` (`:479–495`) admits `sar` / `b1_plus` on any `four_port_birdcage` and never looks at the conductor; `_parse_mapping` (`:498–548`) passes the parsed `conductor` into the `Study`. `run_study_file`'s four-port branch (`:750–762`) calls `study.solve` (the conductor-applying route since step 7) and then `_sar_summary` (`:667–700`) / `_b1_plus_summary` (`:703–711`), which call `Study.sar` / `Study.b1_plus_at` and so `_study_solve_port_fields` (`src/fem_em_solver/study.py:1592–1670`) — since `39fd8a01` that route runs `_conductor_setup` (`:1618`), the hole-route empty `material_map` (`:1625`), `**cond_kwargs` (`:1633`) and the conditional `extra_bilinear_terms` (`:1645`). Before the fix this file would have raised item 388's (F1) `material_map … [1]` inside `_sar_summary` (signature: `20261012T033148Z_WF-16-s11-W0r.log:942`). No test has run it since.
- **The YAML never reaches (F2).** Sections are validated study → geometry → material → excitation → **conductor** → thermal → **outputs** (`:521–542`), so a solid-route file with a `leontovich` conductor raises at `conductor` whether or not it has `outputs`. `StudyConfigError` subclasses `ValueError` (`:106–113`); its text is `"<source>: conductor: a leontovich conductor needs `geometry.as_hole: true` (the hole route)"` (`:346–349`), so it carries the token `as_hole`. The Python `Study` refuses the same combination later, at its first solve-bearing call: `_conductor_setup` raises `ValueError` matching `"as_hole"` before it touches `built` (item 388 (K): `sar`, `b1_plus`, `solve()` on a `{"mesh": None}` stub, `20261012T033238Z_WF-16-s11-W1.log`, 12 passed). **Parity is therefore read as: both refuse with a `ValueError` whose text contains `as_hole`, and neither builds a mesh** — not the same exception type, not the same call (for the review, question 4 of the return).
- **Two solves (disclosed, ruling 2026-10-08 (4)):** `run_study_file` builds and solves through `Study.solve`, then the field-keeping route builds its own mesh and sweeps with `keep_fields=True`. The Python comparand below uses the field route only (no `solve()`), as item 388's records did.

**Files.**
- **`tests/package/wf1/f_small_copper_sar_b1.yaml`** (one-line header citing this item): `f_small_copper.yaml` verbatim plus `outputs: {sar: {mass_g: 10.0, centres: [[0.0, 0.0, 0.0]]}, b1_plus: {points: [[0.0, 0.0, 0.0]]}}` — the origin is item 388's `CENTRE`, where its W1 printed the copper records.
- **`tests/package/wf1/malformed_conductor_solid_with_outputs.yaml`** (header): the same file with `geometry.as_hole: false` written explicitly (the existing `malformed_conductor_without_hole.yaml` omits the key and has no `outputs`).
- **New `tests/package/test_wf1_yaml_copper_sar_b1.py`.** Its docstring's first paragraph backticks **only** `WF-1` and `WF-16` (both §7 rows; `generate_validation_page.py` refuses a derived ID with no row, (F3) — never backtick `U44` or an item number). `@complex_only` on every solve node. Width asserted (never skipped) in the fixture: `MPI.COMM_WORLD.size == ANS6_EXACT_STRIP_RECORD_RANK_WIDTH` (2). **Imports, module-level names only, never a `test_*`:** `LOGS`, `WF1`, `_fresh_dir`, `_mapping`, `_rc_all`, `_spec_diff` (`tests/package/test_wf1_study_yaml.py:85–170`); `COPPER_FILE`, `_py_copper` (`tests/package/test_wf1_yaml_copper.py:84, :98`); `CENTRE`, `_phantom`, `_kept_identity` (`tests/package/test_wf16_study_conductor_fields.py:66, :78, :85`); `_geometry`, `_rel_means` (`tests/package/test_wf16_study_copper.py:80, :76`); `th14` (`tests/validation/test_th14_birdcage_copper.py`); `ANS6_EXACT_STRIP_CLASS_MEANS_P1`, `ANS6_EXACT_STRIP_HOLE_CELL_RECORD`, `ANS6_EXACT_STRIP_RECORD_RANK_WIDTH`, `CLASSES`, `_class_means` (`tests/validation/test_ans6_exact_strip.py:118, :112, :144, :103`, `_class_means` as item 388 imports it); `REUSE_REPRODUCTION_RTOL` (`tests/validation/test_port_birdcage_ring_matrix.py:171`, = `STEP2_COLUMN_POWER_RTOL` = 1e-6); `EXACT_IDENTITY_RTOL` (`tests/validation/test_birdcage_sar_mass_averaged.py:163`, 1e-10); `LEG_D0_REPRODUCTION_BAND` (`tests/validation/test_port_birdcage_four_port.py:165`, 1e-9); `DISCRETE_IDENTITY_RTOL` (`tests/validation/test_birdcage_power_identity.py:227`, 1e-6); `BLIND_FACTOR` (`tests/validation/test_wf12_f_human_copper_matrix.py:61`, 10.0); `FREQUENCY_64_HZ` (`tests/validation/test_lossy_sphere_fullwave.py:94`); `fem_em_solver.study.COPPER_SIGMA`.
- **Module fixture `runs`, in this order:** (1) `_fresh_dir(LOGS / stem)`; `run = run_study_file(CU_SAR_FILE, output_dir=…, comm=comm)`; read `summary.json` back. (2) The same-Study API re-read: `fr = run.study.fields_result`; `api_sar = run.study.sar(mass_g=10.0, centres=[CENTRE], comm=comm)`; `api_b1 = run.study.b1_plus_at([CENTRE], comm=comm)`; assert `run.study.fields_result is fr` (no third build). (3) Python Cu Study A = `_py_copper()` (item 388's Study A, argument for argument); `py_sar = A.sar(mass_g=10.0, centres=[CENTRE], comm=comm)`, `py_b1 = A.b1_plus_at([CENTRE], comm=comm)` — field route only, every drive. (4) Python PEC-hole Study B = `_py_copper("pec")`; `pec_sar = B.sar(mass_g=10.0, centres=[CENTRE], drives=["P1"], comm=comm)`, `pec_b1 = B.b1_plus_at([CENTRE], drives=["P1"], comm=comm)`. (5) `power = _kept_identity(fr, run.study.omega, comm)`. Print every build's `mesh_time` / `sweep_time`, `run_study_file`'s wall time and `summary["elapsed_s"]`. Every `sar` / `b1_plus_at` call is made on both ranks; every power is reduced inside `_kept_identity`.
- **Nodes** (assertions inside the node bodies, never only in the fixture — the validation page traces bands in node bodies only): `test_wf1_s8_parse_and_refusal_parity` (no mesh, uses `monkeypatch`, does not request `runs`); `test_wf1_s8_yaml_field_route_applies_conductor`; `test_wf1_s8_yaml_sar_reproduces_python`; `test_wf1_s8_yaml_b1_reproduces_python`; `test_wf1_s8_kept_field_power_identity`; `test_wf1_s8_negative_pec_column`; `test_wf1_s8_negative_surface_blind`.
- **Registry (lines added only, each block with a `# WF-1 step 8 (§9 item 394)` comment):** the parity node under `[names.load_study]`; the field-route, SAR, `|B₁⁺|` and PEC-column nodes under `[names.run_study_file]`; the field-route, power-identity and surface-blind nodes under `[names.LeontovichWall]`.

**Anchors (asserted; bands imported).**
- **(K) bitwise up to the solve, no mesh:** `_spec_diff(load_study(CU_SAR_FILE), _py_copper()) == []` and `_spec_diff(load_study(CU_SAR_FILE), load_study(COPPER_FILE)) == []` (the `outputs` section adds nothing to the Study); through `study_yaml._load`, `outputs["sar"] == {"mass_g": 10.0, "centres": [(0.0, 0.0, 0.0)]}` with `mass_g` by `float.hex`, `outputs["b1_plus"]["points"] == [(0.0, 0.0, 0.0)]`.
- **(X=) refusal parity, no mesh.** YAML side: `load_study(MALFORMED)` raises `StudyConfigError` with `.field == "conductor"`, `isinstance(e, ValueError)`, `"as_hole" in str(e)`; `cli.main(["run", f, "--output-dir", out])` and `["check", f]` give rc `(2, 2)` through `_rc_all`, and `out` does not exist on any rank (`MPI.LOR`). Python side, inside the node's `monkeypatch.setattr(FourPortBirdcageGeometry, "build", stub)` where `stub` returns `{"mesh": None}` and counts its calls: `solid = Study(dataclasses.replace(_geometry(), as_hole=False), _phantom(), LumpedPortDrive(), frequency_hz=FREQUENCY_64_HZ, degree=1, conductor=LeontovichWall(study_mod.COPPER_SIGMA))`; `solid.sar(mass_g=10.0, centres=[CENTRE])`, `solid.b1_plus_at([CENTRE])` and `solid.solve()` each raise `ValueError` matching `"as_hole"`; `getattr(solid, "fields_result", None) is None` after each. The YAML refusal never calls `build` (stub count 0 after `load_study`).
- **(K1) the YAML Study's field route carries the conductor:** `fr.cells == summary["cells"] == ANS6_EXACT_STRIP_HOLE_CELL_RECORD` (80 798); `fr.conductor_terms` has length 1 and `fr.conductor_facet_tags is not None`; `summary["conductor"] == {"kind": "leontovich", "sigma_s_per_m": COPPER_SIGMA}`.
- **(S=Y) the YAML Study's field route = its `solve()`:** ‖`fr.s_matrix` − `run.result.s_matrix`‖_F / ‖`run.result.s_matrix`‖_F ≤ `LEG_D0_REPRODUCTION_BAND`. Prior, the same comparison on the Python Study A, same fixture: 7.820e-16 (`20261012T033238Z_WF-16-s11-W1.log:2699`). **(R)** `_class_means(fr.s_matrix)` against `ANS6_EXACT_STRIP_CLASS_MEANS_P1["64"]["cu"]`, per class ≤ `LEG_D0_REPRODUCTION_BAND` (prior on the Python field route 4.093e-11 / 8.673e-11 / 1.376e-10, `…W1.log:2703`).
- **(I-SAR) YAML summary SAR = the Python copper Study's field route.** For every drive `P1`…`P4`, every key of the operator dict with a nonzero Python value: |summary − Python| / |Python| ≤ `REUSE_REPRODUCTION_RTOL` (the band 382b / 390b use for a SAR reproduction across Studies). Predicted ≤ 1e-12 (independent builds and factorisations; (S=) on the same fixture reads 7.8e-16), printed at `.3e`. **(iii)** the summary against `api_sar` on the **same** Study ≤ `EXACT_IDENTITY_RTOL` (item 344 read 0.0 on its fixture), and the JSON file's `sar` / `b1_plus` blocks equal `run.summary`'s exactly. **Printed, not asserted:** summary SAR10g(P1, origin) at `.9e` beside item 388's printed Cu record 2.408073397e-06 W/kg (`…W1.log:2694`); predicted equal at every printed digit (the `.9e` print alone contributes ≤ 2.1e-10 relative); a digit difference with every assertion green is reported, not a stop.
- **(I-B) YAML summary `|B₁⁺|` = the Python copper Study's field route.** Per drive, |summary − Python| / |Python| ≤ `REUSE_REPRODUCTION_RTOL`; `summary["b1_plus"]["valid"]` all true and `py_b1` valid. **(B2′)** the summary against `api_b1` on the same Study, max |Δ| / max |B₁⁺| ≤ `EXACT_IDENTITY_RTOL` (item 344 read 0.0). **Printed only:** P1 at the origin beside item 388's 1.833178311e-08 T per 1 V (`…W1.log:2694`), predicted equal at `.9e`.
- **(P) the exact power identity on the YAML Study's kept field:** `power["p_surf"] > 0` and `power["residual"] ≤ DISCRETE_IDENTITY_RTOL`. Prior, the same identity on the Python Study's kept field, same fixture, same drive: residual 8.505e-13, P_src 4.806472101e-03 W, P_surf 1.516873149e-06 W, P_phantom 1.870050049e-06 W (`…W1.log:2692`), printed beside at `.9e`.

**Negative controls (rule (e)).**
- **(N1) the PEC column — asserted, backed.** The YAML Cu summary's SAR10g(P1, origin) misses the PEC Study B's by > `REUSE_REPRODUCTION_RTOL`, and the same for `|B₁⁺|`(P1, origin). Backing: item 388's W1 printed both columns of this comparison on this fixture at this point and drive — SAR 2.408073397e-06 vs 2.409072917e-06 W/kg, `|B₁⁺|` 1.833178311e-08 vs 1.833834143e-08 T (`20261012T033238Z_WF-16-s11-W1.log:2694–2695`), a relative gap of 4.151e-04 (SAR) and 3.578e-04 (`|B₁⁺|`), so ≈ 415× and ≈ 358× the band; the factors are printed. The ceiling is that physical Cu-vs-PEC gap itself; a 1× bar is reachable with more than two decades of margin. This is the control that separates "the YAML file's conductor reached the field route" from a PEC-walled field. (Item 388's backing pair was Python-vs-Python; here the Cu side is the YAML Study, equal to the Python one by (I-SAR). The review may relabel it *predicted*, as it might have 388's (N1); the assertion would then be dropped and the factors printed.)
- **(N2) surface-blind identity — asserted, backed.** The kept-field residual with P_surf omitted is ≥ `BLIND_FACTOR` × `DISCRETE_IDENTITY_RTOL` (= 1e-5). Backing: 3.1559e-04 = 315.6× on the same kept-field route, fixture and drive (`…W1.log:2718`); the ceiling is P_surf / P_src ≈ 3.16e-04, so ≈ 316× is the most this fixture can show and the 10× bar holds with ≈ 32× margin; the factor is printed.
- **(X=)'s control is deterministic, asserted:** the malformed file's mapping (`_mapping`) with `geometry.as_hole` set to `True` parses through `_study_from_mapping` to a Study with `_spec_diff` against `load_study(CU_SAR_FILE)` of `[]` (the refusal is that one bit's alone), and with `outputs` deleted it still raises at `conductor` (the refusal is the conductor's, not the outputs').

**Windows.** Foreground. **The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped.** Bash `timeout: 660000` on every window; pytest with `-s`; `tests/environment/test_complex_mode.py` first on the path; node IDs, never `-k`. Every window is under 10 min, so no durable capture. Before W1: `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reads 0. Windows run serially, never beside another heavy window.
- **W0 — the no-mesh node (standard; predicted 20–40 s):** `scripts/testing/run_and_log.sh WF-1-s8-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_copper_sar_b1.py::test_wf1_s8_parse_and_refusal_parity -v -s --tb=short'"`.
- **W1 — the module (heavy by ceiling; predicted 140–190 s, ≤ 250 s at the 1.32× box reading):** `scripts/testing/run_and_log.sh WF-1-s8-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 400 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_copper_sar_b1.py -v -s --tb=short'"`. Price, from item 388's W1 (`…W1.log:2691`, 99 s for three builds): four builds of the 80 798-cell hole mesh at ≈ 30 s each (mesh 23.2 s + sweep 6.3–6.7 s) — `run_study_file`'s `solve()` and field route, Study A's and Study B's field routes ≈ 120 s; 4 + 4 + 1 balls and three CG1 projections ≈ 12 s (item 388: two balls and two projections with evaluation 4.8 s, `:2696`); identity 0.4 s; start-up, `test_complex_mode.py` and the no-mesh node ≈ 15 s ⇒ ≈ 150 s.
- **Wgen / Wval / Wapi / Wdoc (the registry lines move the API pages and `docs/guide/validation.md`; commit both):** `scripts/testing/run_and_log.sh WF-1-s8-Wgen "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 python3 scripts/docs/generate_api_reference.py --out docs/guide/api'"` (2 s, `20261011T214232Z_WF-16-s9a-Wgen.log:12`); then `scripts/testing/run_and_log.sh WF-1-s8-Wval "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace timeout -k 30 60 python3 scripts/docs/generate_validation_page.py --out docs/guide/validation.md'"` (`20261012T034338Z_WF-14-s2-W0-gen.log:12`); then `scripts/testing/run_and_log.sh WF-1-s8-Wapi "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace timeout -k 30 170 python3 -m pytest tests/package/test_wf14_validation_envelope.py tests/package/test_wf17_api_reference.py tests/package/test_wf14_validation_page.py -v -s --tb=short'"` (22 passed, 14 s, `20261012T034345Z_WF-14-s2-W1-page.log:12, :220`); then `scripts/testing/run_and_log.sh WF-1-s8-Wdoc "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && timeout -k 30 60 python3 scripts/testing/check_example_doc_references.py --docs-root docs/guide --check-links'"` (`0 dead`; `20261011T214244Z_WF-16-s9a-Wdoc.log:12`).
- **No rule (c) re-gate:** no `src/` file changes. The executor confirms `git diff --stat main -- src/` is empty before the commit and says so. If the module shows a defect in `src/`, it is a finding — known-issues, park — not an in-slot fix.

**Traps already paid for.**
- **Two `REUSE_REPRODUCTION_RTOL`s exist:** `tests/validation/test_port_birdcage_ring_matrix.py:171` (1e-6, the one 382b / 390b import) and `tests/validation/test_port19_factor_reuse.py:82` (1e-12). Import the first by its full module path.
- **The stub is function-scoped.** The parity node uses the `monkeypatch` fixture (restored after the node) and never requests `runs`; a module-scoped stub would hand `{"mesh": None}` to the fixture's builds. `_conductor_setup` runs after `build` and before the problem, which is why the stub needs only `"mesh"` (item 388 trap).
- **One field build per Study.** `run_study_file` already populated `run.study.fields_result` through `_sar_summary`; the API re-read must reuse it (asserted by identity). Never clear `fields_result` or `_b1_cg1`, never build a fresh Study for (iii) / (B2′) — that is a third ≈ 30 s build and a different comparison.
- PyYAML 1.1 floats need a dot (`0.0`, `10.0`, `58000000.0`; `WF-1` step 1); `false` is a `bool`, `0` an `int` (refused at `geometry.as_hole`, step 4).
- Every `sar` / `b1_plus_at` call is collective on both ranks (item 227's deadlock); `points` and `centres` identical on every rank; point evaluation only inside `b1_plus_at` (`evaluate_vector_field_parallel`). `cell_tags.values` and `assemble_scalar` are rank-local — `_kept_identity` reduces each power. rc is rank-local — `_rc_all`.
- Never import a `test_*` name or a module-scope solving fixture; helpers and constants only. No solve output compared bitwise (`MAT-1` item 81; `MAT-5` item 237; rule (j)'s "bitwise up to the solve").
- The validation page (item 375, after items 390a / 390b were drafted) reads chunk IDs from the module docstring's first paragraph and bands from node bodies: backtick only §7 row IDs there, and keep each asserted band in the node that asserts it.
- On status 124 or a kill: `pgrep -c python3` first; at 0 only, the rubric's orphan-stem sweep as one bare call labelled `WF-1-s8-cache-sweep`. Never raise a timeout or change `-n`.
- Rule (i): no edit to the test module, the YAML files or the registry after the window that closes on them. Write ≈, never a tilde. No AED number anywhere.

**Scope.** One fixture (the F-small exact-rectangle copper hole mesh, 80 798 cells, κ off, `strip_fraction = GATED_WIDTH_FRACTION`), 64 MHz, degree 1, `-n 2`, `LeontovichWall(COPPER_SIGMA)`, one point (the origin), four drives. **Not claimed:** any SAR or `|B₁⁺|` value on copper as a physical figure (per 1 V records, compared with themselves through two front ends); the `{ring_radius_m}` centre form on the hole route; the PEC-hole YAML file with `outputs` (its Study is reached only as the Python control B); `femem run` with a solve through `cli.main` or the installed console script on this file (the CLI is exercised on the refusal only); degree 2; other frequencies or widths; F-human (items 390a / 390b); a C4 statement at the origin; any AED statement.

**What a negative result means.** Report and stop; park on `attempt/WF-1-step8-<UTC>`, mark item 394 BLOCKED with its unblock condition in the same commit (rule (d)), journal in `docs/testing/attempts.md`. (K) red — the parse moved: a `study_yaml` finding, known-issues, never a test edit to match. (X=) red on the YAML side — the refusal order or text moved; on the Python side — item 388's fix regressed: known-issues, stop. (K1) / (S=Y) red with (R) green — the YAML Study's field route does not carry the conductor although the parse does: a defect against item 388's fix, known-issues, never a widened band. (I-SAR) / (I-B) red with (K), (S=Y) and (iii) / (B2′) green — two equal Study specs give different fields across builds: print per-drive residuals, known-issues, never re-band. (iii) / (B2′) red — the summary is not the API's value on its own Study: a `run_study_file` finding. (P) red with (N2) green — kept-field accounting; (N2) or (N1) inside its bar — the YAML field carries no wall term (the conductor did not reach the field route): stop, known-issues. A Wapi red after Wgen / Wval — regenerate once after a rebase, else park. No `PREDICTED RED:` line: a red anywhere stops the slot.

**Status it can move.**
- **§7 `WF-1` row:** stays ✅ (rule (j): its Done-when is step 1's, which this item does not touch). On green the row's history gains one step-8 record clause, exact text from the review's return, of the form "Step-8 record (§9 item 394, U44, provisional): `outputs.sar` / `outputs.b1_plus` with a `leontovich` conductor on the F-small hole Study return the Python copper Study's field-route 10 g SAR and `|B₁⁺|` at the origin within `REUSE_REPRODUCTION_RTOL` (max <x> / <y> over four drives), the YAML field route equals its `solve()` to <z> (`LEG_D0_REPRODUCTION_BAND`), the kept field's copper identity closes at <w> (`DISCRETE_IDENTITY_RTOL`); a solid-route `leontovich` file refuses at `conductor` as `Study.sar` / `b1_plus_at` / `solve()` refuse (`as_hole`), with no mesh".
- **§10 chain U:** U44 marked DONE (provisional until the next weekly ratifies).
- **`WF-16.md`, item 388's Scope:** the "YAML `conductor` + `outputs.sar` / `outputs.b1_plus` combination … untested" clause is answered by this record (a pointer line under item 388's outcome; the clause itself is not edited). The `WF-16` row does not move.
- **`WF-14`'s registry** gains the seven node IDs; no known-issues entry is filed on green. No §2 sentence moves.

**Executor spawn prompt (verbatim).** "Execute §9 item 394 — `WF-1` step 8 (§10 U44, provisional: copper SAR / |B₁⁺| through the YAML front end, and solid-route refusal parity) — per `docs/planning/chunks/WF-1.md` §9 item 394. Read the brief and any review rulings beneath it first; the rulings hold. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Pass timeout: 660000 on every Bash call that runs a window; run pytest with -s. Write the two YAML files and the test module, then run W0, W1, Wgen, Wval, Wapi, Wdoc in order, foreground. Before W1 run `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` and require 0. If W1 has not started by minute 35, park. No `src/` edit: if the module exposes a `src/` defect, file it in known-issues and park. Import `REUSE_REPRODUCTION_RTOL` from `tests/validation/test_port_birdcage_ring_matrix.py`; never import a `test_*` name. No edit to the test module, the YAML files or the registry after the window that closes on them. Never delete a log. You are the executor — do not spawn agents, never background a harness window. If the docker socket is denied, stop and report the exact command string. On a red, park on `attempt/WF-1-step8-*`, mark item 394 BLOCKED in the same commit, and report."

**Review ruling 2026-10-12 03:00 (daily) on the item 394 draft, and on items 390a / 390b.**
(Q1) Agreed: item 388's printed SAR / \|B₁⁺\| values stay printed beside the YAML values, never asserted (that would register a record under rule (6)). (Q2) `REUSE_REPRODUCTION_RTOL` accepted, as 382b / 390b. (Q4) Refusal parity as drafted — same `ValueError` family, the `as_hole` token in both messages, no mesh built on either side — accepted; it is a parity of outcome, not of call site, and the brief says so. (Q5) (N1) stays **asserted**: it is the same cu-vs-pec comparison on the same fixture, measured at `20261012T033238Z_WF-16-s11-W1.log:2694–2695`; the route that produces the cu column differs, the comparison does not. (Q6) Order in §9 is 390a, 390b, 394; whichever of them lands second re-runs the other's closing module in its own re-gate window if it touched the same registry file.
**(Q3) — binding on items 390a, 390b and 394: any item that adds lines to `WF-14`'s registry regenerates the API pages and `docs/guide/validation.md` and runs the page gates in the same slot** — the four windows of item 394's "Wgen / Wval / Wapi / Wdoc" bullet, with the item's own label prefix (`WF-1-s7a-…`, `WF-1-s7b-…`), committed with the regenerated pages. Item 375 (`cabfed7e`) made `docs/guide/validation.md` byte-checked against the registry; 390a / 390b's briefs predate it, and without these windows their commit would leave `test_wf14_validation_page.py`'s (G) red on `main`. This adds ≈ 2 minutes to each and moves no band.

### §9 item 401 — 2026-10-12 (`WF-1` step 9: `femem check` prints the parsed geometry kind; lane 1, ops)

*Written by the 2026-10-12 10:30 daily review, from item 390a's finding (1).*

**What.** `cli.py:158` labels every study that is not a sphere `four_port_birdcage`. So `femem check tests/package/wf1/f_human_tuned_copper.yaml` prints `family: four_port_birdcage` (`20261012T123618Z_WF-1-s7a-W0.log:68`). `run_study_file` records the same file as `"family": "f_human_birdcage"` (`study_yaml.py:864`, the parsed `kind`). Tracked as known-issues "🟡 OPEN 2026-10-12 (found by `WF-1` step 7a …)".

**Edit (`src/fem_em_solver/cli.py` only, the `check` branch).** Print the parsed `geometry.kind`, the same `kind` that `study_yaml._load(path)` returns and `run_study_file` writes. Take it from `_load`, or from a small public wrapper in `study_yaml` that returns `(study, kind)` if the slot prefers not to import a private name; disclose which. Nothing else changes: not the other five printed lines, the return codes or the `StudyConfigError` path. `study_yaml.py` gains at most that wrapper. The parse and solve paths are not touched.

**Anchors (new `tests/package/test_wf1_check_family.py`, `-n 2`, no mesh, no solve).**
- (F) Take every tracked YAML under `tests/package/wf1/` that is not malformed (glob `*.yaml`, exclude `malformed_*`; print the list, ≥ 10 files). For each, `cli.main(["check", f])` returns 0 on every rank, and its printed `family:` line (captured on rank 0 with `capsys`) equals `study_yaml._load(f)[1]`. Print each pair.
- (F-h) Three literals, the letter of the known-issues entry:
  - `f_human_tuned_copper.yaml` prints `family: f_human_birdcage`;
  - `sphere_mat4.yaml` prints `sphere_in_box`;
  - `f_small_four_port.yaml` prints `four_port_birdcage`.
- (U) For `sphere_mat4.yaml` and `f_small_four_port.yaml`, the other five lines of `check`'s output are unchanged. Compare them against the output of the **pinned pre-change** `cli.py` (`b01dc790`, never `HEAD:`). The host writes that copy under the gitignored `logs/` with a plain `git show b01dc790:src/fem_em_solver/cli.py > logs/cli_b01dc790.py`, and the test loads it by path. Loaded that way, the copy's relative imports need its package context. If that does not work cleanly, compare against the five lines printed in `20261012T123618Z_WF-1-s7a-W0.log` and the RG logs of item 390a instead, and say which.

**Negative control (asserted; backed by the measurement above, the same comparison on the same file).** The pre-change `check` on `f_human_tuned_copper.yaml` printed `family: four_port_birdcage` (`…WF-1-s7a-W0.log:68`), which differs from `_load(...)[1]`. Assert that this logged line differs from the string (F) now prints. Without this control, (F) could pass on a file set that never reaches the bug.

**Re-gates (standing rule (c): the pre-existing gates on the edited module).** These modules call `cli.main(["check", …])`:
- `tests/package/test_wf1_study_yaml.py`
- `tests/package/test_wf1_yaml_sar_b1.py`
- `tests/package/test_wf1_yaml_probes.py`
- `tests/package/test_wf1_yaml_f_human.py`

Grep `cli.main` under `tests/package/` and add any others, plus `test_wf13_launch_helper.py`, which asserts `s["family"]`. Run them in one `-n 2` window at standard tier, skipping and naming any that is solve-gated behind an opt-in variable. If item 390a's W9 window exercises the installed console script, re-run it too.

**Registry.** This item adds no `WF-14` registry line, so (Q3)'s page windows are not owed. If the slot does add one, it runs Wgen / Wval / Wapi / Wdoc as item 394 did.

**Tier, cost.** Smoke for W1. The re-gate window is standard, ≈ 110–150 s at `-n 2` (item 390a's RG1–RG8). ≈ 3 measured (ops / docs median) / 20 predicted min. Depends on nothing.

**Windows (bare, whole, foreground, Bash `timeout: 660000`):**

`scripts/testing/run_and_log.sh WF-1-s9-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src timeout -k 30 60 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_check_family.py -s -v --tb=short'"`

After W1, run one re-gate window `WF-1-s9-RG` (standard, `timeout -k 30 180`) over the modules listed above. Copy the exact paths and env from item 390a's RG logs.

**Traps already paid for.**
- `capsys` sees only this rank's output: assert on rank 0, and reduce the rc over ranks (item 390a's `_rc_all`).
- The pinned control commit is `b01dc790`, never `HEAD:`.
- `-k a or b` splits inside the quoted container command.
- Pytest prints need `-s`.
- Rule (i): no edit after the closing window.

**Scope.** A printed label. Not claimed: anything about the parse, the Study or the solve, none of which is touched. The `--estimate` path is out of scope. If it also prints a family, report that and do not fix it.

**Status it can move.** The known-issues entry is retired in the same commit. `WF-1` stays ✅, and its row gains nothing; the outcome goes in this file.

**What a negative result means.**
- A red re-gate means the edit changed more than the label: revert, mark item 401 BLOCKED, and leave the entry open.
- If (F) finds a second file whose kind disagrees for some other reason, report it, file it in known-issues and stop.

There is no `PREDICTED RED` line.

**Executor.** The slot runs this item directly (ops), or spawns `implementer` with the §9 "Every item" bare-command sentence verbatim in the prompt.

**Outcome 2026-10-12 (12:00 CDT slot, run directly, fourth item of the slot): DONE.**

**The edit.** `cli.py`'s `check` branch now takes `study, kind = _load(args.file)[:2]` and prints `family: {kind}`. It uses the private `study_yaml._load`, imported inside the package; no wrapper was added. The unused `g` binding is dropped. The other five lines, the return codes and the `StudyConfigError` path are unchanged.

**W1, the window of record** (`20261012T172853Z_WF-1-s9-W1.log`, 6 passed, `-n 2`, 2 s). New `tests/package/test_wf1_check_family.py`:
- (F): 11 files; the printed family equals `_load(f)[1]` on each, and the rc is (0, 0).
- (F-h): `f_human_birdcage` / `sphere_in_box` / `four_port_birdcage`.
- (U): the other five lines on `sphere_mat4.yaml` and `f_small_four_port.yaml` equal the pinned `b01dc790` `cli.py`. It is loaded by path as `fem_em_solver._cli_b01dc790`, so its relative imports resolve, and the logged-line fallback was not needed.
- (N): the logged pre-change `family: four_port_birdcage` differs from the new `family: f_human_birdcage`.

**The first W1** (`20261012T172825Z_WF-1-s9-W1.log`, 6 passed) is superseded. `capsys` swallowed the record prints into the next `readouterr()`. The assertions were unaffected, and the prints now go through `capsys.disabled()`.

**Re-gate** (`20261012T172909Z_WF-1-s9-RG.log`, 13 passed, 3 s). It covers every parse / refusal node that calls `cli.main(["check", …])` in `test_wf1_study_yaml`, `_yaml_probes`, `_yaml_sar_b1`, `_yaml_f_human` (two nodes), `_yaml_copper`, `_yaml_thermal` and `_yaml_copper_sar_b1`. The solve-gated nodes of those modules and `test_wf13_launch_helper.py` were not run. They assert `run_study_file`'s summary family, which this edit does not touch. The installed console script (item 390a's W9 command) passed 1 of 1 (`20261012T172921Z_WF-1-s9-W9-installed-femem.log`).

**Status.** No registry line was added. `--estimate` prints no family. The known-issues entry is retired.

### §9 item 431 — 2026-10-13 (§10 chain U46; chain step, solve) — drafted as 431a (no solve) + 431b (solve)

**Outcome 431a 2026-10-14 (implementer, slot 00:30 UTC, third item): DONE, green.** W00 `20261014T005051Z_WF-1-s10a-W00-literals.log` (2 s): literals printed from the series JSON, YAML written from them. W0 `20261014T005102Z_WF-1-s10a-W0.log` 6 passed, 4 s (`-n 2`): (K) both `_spec_diff` `[]`; `source_scale` `0x1.1352fb069c3cep+5` = `G_per_w`, `dt_s` `0x1.f59999999999ap+4` = T\*/N, ρc, θ, c\*_θ and `window_s` = `W_B2` by hex; `femem check` rc (0, 0), `family: f_human_birdcage`; (X) 9 + 390a's k/h/t_inf mutation name their field, cli run / check on a transient-less file rc (2, 2), nothing written; (c1) `['source_scale']` = default 1.0, (c2) one ULP `['source_scale']`; (Run-0) tuned solve (1, 1), temperature hook (1, 1), summary `temperature` = fake, JSON = summary; without `thermal:` hook 0 calls, no key. Re-gates (rule (c)) all green: RG1 15 / 116 s, RG2 13 / 160 s, RG3 6 / 35 s, RG4 11 / 95 s, RG5 13 / 120 s, RG6 9 / 37 s, RG7 6 / 47 s, RG8 5 / 73 s, RG9 6 / 4 s, RG10 6 / 3 s, RG11 14 / 102 s, RG12 11 passed + 4 skipped / 200 s; `pgrep -c python3` 0 before RG2 / RG4 / RG5 / RG8 / RG12. **RG13 as written went red** (`20261014T010844Z…RG13-wf15-probe.log`: "(W) study.py or the test module changed since the -n 2 npz was written") — a sha pin on `study.py` in the `-n 2` npz, not a value: the `-n 1` node needs its `-n 2` writer re-run after any `study.py` edit (OPS-70 ran W1 before W2). Re-run as RG13a (the writer nodes of OPS-70's W1, `…010906Z…RG13a-wf15-probe-n2.log`, 6 passed, 23 s) then RG13b (RG13 verbatim, `…010931Z…RG13b-wf15-probe.log`, 5 passed, 26 s) — **for the review: future briefs that edit `study.py` must list the writer before RG13**. W9 installed check 1 passed 15 s (`…011009Z`); Wgen 44 files (`…011025Z`), Wval 182 nodes (`…011027Z`), Wapi 22 passed 15 s (`…011030Z`), Wdoc 0 dead (`…011045Z`). Deviations: `_study_temperature` also gained keyword-only `sense=None` and `comm=` (the hook takes `sense` / `comm`; `run_study_file` passes them), refused with `balls` / `window_s` on a non-tuned Study; omitted `window_s` means the API's 360 s (`_SIX_MINUTE_WINDOW_S`). Item 431b is runnable.

*Drafted by a general-purpose agent for the 2026-10-13 03:00 daily review; the review reviews it.* Licence: §10 chain U, "U46 (9a). Thermal through the F-human YAML: `tests/package/wf1/f_human_tuned_copper.yaml` gains a `thermal:` section; `femem run` reproduces P17's ball-mean θ(t\*) and `MAT-20` step 3's six-minute figure at `REUSE_REPRODUCTION_RTOL` and writes the series. Heavy (`-n 8`). Serial on nothing." (PROJECT_PLAN.md:4782–4786), weekly-enumerated (§9 rule (4)). **Target: reproduction. No band, tolerance, Done-when, record or validation target is set or moved; every band and record is imported by name.** The rulings of items 361 (`MAT-15.md`), 380 (`MAT-20.md`) and 390a / 390b (this file), and the capacitor-ESR and sheets-and-κ rulings, bind.

**What exists, and what the step must add (read on `59809ba7`).**
- **The comparands are not module constants.** P17's ball-mean θ̄(c\*_θ, t\*) and the series it came from live in the gitignored `logs/mat15_s5/theta-series-ccw.json`, written by item 380 (4 422 bytes, present on the host 2026-10-13). Its keys:
  - `times_s` (41 samples, Δt 31.35 s, t\* 1254.0 s) and `G_per_w` = 34.41551785625744 1/W (ruling 361-1: 1/P_acc,ccw);
  - `c_star_m` = [0.044135337618145365, 0.008779064490725769, 0.0], `c_star_label` `r=45 dz=+0 az17+0`, `a10_m` 0.013365046175719762, `ball_volume_m3` 1.0010236026341093e-05;
  - `theta_ball10g_k` (last 0.049261108496515285), `theta_mean_k` (last 0.035530244441182766);
  - `cells` 369 387, `n_sub` 3 306, `ranks` 8.

  The logged prints match: θ̄(t\*) 4.926110850e-02 K and mean 3.553024444e-02 K (`20261011T183420Z_MAT-20-step3-W1.log`); the six-minute maximum 4.277005e-02 K at t\* = 1 254 s (`20261011T183717Z_MAT-20-step3-W2.log`). That log prints 7 digits, so the 1e-6 comparison reads the JSON.
- **The API refuses today.**
  - `Study.temperature` calls `_refuse_f_human` (`src/fem_em_solver/study.py:526`), which raises on an F-human Study.
  - It allows `transient=` only on a thermal-only ball (`:1110–1130`).
  - The YAML refuses `thermal:` on `f_human_birdcage` at field `thermal` (`src/fem_em_solver/study_yaml.py:477–479`).
  - 390a's module asserts that refusal: mutation `{k_w_per_m_k: 0.6, h_w_per_m2_k: 10.0, t_inf_c: 20.0}` ⇒ field `thermal` (`tests/package/test_wf1_yaml_f_human.py:215–216`, `:233`).
- **The building blocks already in `src/`.**
  - `_tag_submesh`, `_exterior_robin_tags`, `_dg0_cell_average` and `_dg0_to_submesh` in `study.py:1040–1105` (MAT-10's route, the same route P17's `_cell_average_dg0_tag` / `_to_sub` take).
  - `TransientSpec` (`:1239`).
  - `solve_transient_conduction` (`src/fem_em_solver/thermal/transient.py:94`) and `running_average` (`src/fem_em_solver/thermal/dose.py:108`).
  - `averaging_ball_radius` (`src/fem_em_solver/post/sar.py:214`) and `SAR_QUADRATURE_DEGREE` = 16 (`study.py:1589`, equal to the test-side `QUADRATURE_DEGREE`, `tests/validation/test_mass_averaged_sar_standard_masses.py:121`).
  - **Missing:** the tuned-F-human thermal branch, P17's 10 g ball-mean forms (indicator `ufl.lt(ufl.real(dot(off, off)), a10²)` at `QUADRATURE_DEGREE`, `tests/validation/test_mat15_f_human_tuned_heating.py:426–455`), the per-step series, and the YAML / `run_study_file` wiring.
  - **Trap:** `thermal.probe.ball_average` (`thermal/probe.py:168`) is a lattice-refined estimator, **not** P17's. Using it would not reproduce θ̄ at 1e-6. The 5 mm probe read 4.987727e-02 K against the ball's 4.926e-02 K at a different radius (item 380).
- **Window price.** Solve + thermal fit one 590 s window. 390b's W1 was 243 s (build 113.8 s, node pair 26.6 s, 38 SAR balls ≈ 73 s), and P17's thermal half was 4.3 s plus 40 balls in 0.8 s (`20261011T003520Z_MAT-15-step5-W1.log`). A thermal-only file (no SAR / B outputs) prices ≈ 170 s, ≤ 230 s on the 1.32× box.
- **Why the split anyway: code size, the 390a / 390b precedent.** The `src/` work (≈ 30 min) plus the rule (c) re-gates (≈ 15 min of windows) plus 15 fixed is past the hour on its own. Item 382 was declined at minute 29.8 for the same reason. **431a** lands every `src/` change, the tracked YAML, the no-mesh gates and the re-gates. **431b** (serial) adds one test module and runs the `-n 8` window. U46 is done only when both land.

**File choice (deviation from the weekly's letter, for the review).** The brief writes a **new** `tests/package/wf1/f_human_tuned_copper_thermal.yaml` = `f_human_tuned_copper.yaml` + a `thermal:` section, and leaves out the file's `outputs:` (−73 s). It does not edit the existing file. Item 295's precedent is a new `sphere_mat10_thermal.yaml` beside `sphere_mat4.yaml`.
- **Why.** Editing the existing file changes its sha256 and summary, under 390a's (K-c) and 390b's closed record. It would force 390b's heavy window to run on a different file, and two solving modules in one `-n 8` window (≈ 2 × 150 s) price past 590 s on the slow box.
- **If the review holds the letter,** 431b's W1 runs 390b's module and the new one in two serial windows instead, and nothing else changes.

---

#### 431a — `WF-1` step 10a: the tuned-F-human thermal API, the schema and the run branch, gated without a mesh (chain step, no solve)

- **Executor:** implementer.
- **Tier:** W0 standard (`-n 2`, complex); re-gates as measured.
- **Class:** chain step (no solve), median **7 measured** min. **Cost: ≈ 60 predicted** — 15 fixed + ≈ 30 code + W00 / W0 ≈ 1 + re-gates ≈ 12.5 + pages / W9 ≈ 1 + ≈ 1 spare. **If W0 is not green by minute 32, park** (rule (d)).
- **Family:** `WF-1` step 10, attempt 1 of 4.
- **Depends on:** nothing.

**`src/` edits** (every touched existing function is named for rule (c)).
- **`study.py`, the hook.**
  - `TwoPortQuadratureDrive` gains `_temperature_port_study(study, cooling, *, transient, balls, window_s, sense="ccw", comm)` beside its `_sar_port_study` / `_b1_plus_port_study` hooks (the existing pattern). These are new lines.
  - `_study_temperature` gains **two inserted lines at its top**: `hook = getattr(self.excitation, "_temperature_port_study", None)` and `if hook is not None: return hook(self, cooling, transient=…, balls=…, window_s=…)`.
  - Its signature gains keyword-only `balls=None, window_s=None`. When they are not `None` on a non-tuned Study, it raises `ValueError` before any mesh.
  - Every existing statement is unchanged. `_refuse_f_human` stays for the solid route.
- **`study.py`, the branch body — P17's route, lifted with provenance comments** (`test_mat15_f_human_tuned_heating.py` at `59809ba7`, line spans named), never imported:
  1. `res = _tuned_result(study, comm)`, so the field is reused and there is no second solve.
  2. The source: q = `cooling.source_scale` × the DG0 average of ½σ\|E_sense\|² over `PHANTOM_CELL_TAG`, through `_dg0_cell_average` then `_dg0_to_submesh`. The scale is applied **after** the average (ruling 361-4: one compiled form).
  3. Robin on every exterior submesh facet (`_exterior_robin_tags`).
  4. `solve_transient_conduction(sub, ft, k, rho_c, q, h, t_inf, t_init=t_inf, robin_tag, dt=transient.dt_s, n_steps=transient.n_steps, theta=transient.theta, startup_half_steps=transient.startup_half_steps, degree=cooling.degree, snapshot_steps=tuple(range(n_steps + 1)))`.
  5. Per snapshot: the volume mean ⟨θ⟩ₙ and, per ball centre, θ̄ₙ through P17's indicator forms, with radius `averaging_ball_radius(mass_kg=mass_g·1e-3, rho=study.density_kg_m3)` and `SAR_QUADRATURE_DEGREE`.
  6. Per series, `running_average(times, series, window=window_s)` → (peak, t_max).
- **`study.py`, the result.** A new non-exported `TunedTemperatureResult` carrying `times_s`, `volume_mean_k`, `balls` (centre, radius, volume, `theta_k`), `six_minute` (peak, t_max per series), `source_power_w` (∫q), `max_balance_residual`, `cells`, `n_sub`, `sense`, `cooling` and `transient`.
- **Collectives.** Every ball assembly and reduction is collective on every rank. Never guard a collective with "owns phantom cells" (item 227; `MAT-19` D4's empty ranks).
- **`study_yaml.py`.**
  - `_parse_thermal`'s F-human `elif` (`:477–479`) no longer refuses outright. It calls a new `_parse_thermal_f_human(sec, source)` that returns a private `_FHumanThermal(cooling, transient, balls, window_s, sense)`. The keys:
    - `k_w_per_m_k`, `h_w_per_m2_k`, `t_inf_c`, `source_scale`;
    - **`transient: {dt_s, n_steps, rho_c_j_per_m3_k, theta?, startup_half_steps?}`** — required;
    - **`balls: {mass_g, centres: [[x, y, z], …]}`** — required;
    - `window_s` — optional; omitted means the API default, never a loader default (the (Deg) rule);
    - `sense` — optional, `ccw` / `cw`.
  - **A section without `transient` raises at field `thermal`** with "f_human_birdcage needs `thermal.transient` (a transient Study)". This mirrors 390a's "`two_port_quadrature` needs a conductor ⇒ `conductor`" rule, and it keeps 390a's mutation assertion green **unedited**. A malformed key inside `transient` or `balls` raises at its dotted field. Order inside the section as everywhere: unknown keys, missing, types, values.
  - `_parse_mapping` carries the result in the cooling slot of its return tuple for F-human only. Every caller (`_study_from_mapping`, `_load`, `load_study`, `run_study_file`) is updated in the same diff (item 295's rule).
  - **`run_study_file`'s F-human branch** gains, after the existing S / SAR / B lines and before the write: `if fh_thermal is not None: thermal = study.temperature(...)`. It then writes `summary["temperature"] = {"sense", "times_s", "volume_mean_k", "balls": [{"centre", "radius_m", "volume_m3", "theta_k"}], "six_minute": {...}, "source_power_w", "max_balance_residual", "source_scale"}`. The series *is* the summary, and the bundle writes it. `StudyRun.thermal` is set.
  - The four-port and sphere branches are untouched.
  - **Name no artifact-shaped token in any docstring** (item 274's trap).
- **Not touched:** `cli.py`, `__init__.py`, `public_surface.txt`. No new public name, since `TunedTemperatureResult` is returned, not exported.

**Tracked YAML: `tests/package/wf1/f_human_tuned_copper_thermal.yaml`.** It is `f_human_tuned_copper.yaml`'s `femem` / `study` / `geometry` / `material` / `excitation` / `conductor` sections verbatim (the existing file is not edited), with no `outputs:`, plus:
- `thermal: {k_w_per_m_k: 0.6, h_w_per_m2_k: 10.0, t_inf_c: 20.0, source_scale: <repr(G_per_w)>, transient: {dt_s: <repr(T_STAR / N_MAIN)>, n_steps: 40, rho_c_j_per_m3_k: <repr(RHO_C)>, theta: 0.5, startup_half_steps: 2}, balls: {mass_g: 10.0, centres: [<repr triple of c_star_m>]}, window_s: 360.0}`.
- **Literals come from W00's print, never typed.** The `.json` is read on the host-mounted `/workspace/logs/mat15_s5/`.
- The constants are step 1's `K_GEL`, `H_GEL`, `T_INF`, `T_STAR`, `N_MAIN` (`tests/validation/test_mat15_coil_driven_heating_phantom.py:82–89`) and `MAT-12`'s `RHO_C` (`tests/validation/test_mat12_transient_sar_heating_sphere.py:68`).
- It carries a one-line header citing this item and its source (the JSON's `source` string). Floats are Python `repr` with a dot (PyYAML 1.1).

**New `tests/package/test_wf1_yaml_f_human_thermal.py` (no mesh anywhere).**
- `FHT_FILE` and `SERIES_JSON` (= P17 module's `SERIES_JSON`, `:134`, imported).
- `_literal_lines()` reads the JSON on rank 0, broadcasts it, and prints the YAML literals.
- **The docstring's first paragraph backticks only `WF-1`, `MAT-15` and `MAT-20`.**
- **Imports (module-level names only; never a `test_*`):** from `test_wf1_study_yaml`, `_spec_diff`, `_rc_all`, `_p0`, `LOGS`, `WF1`; from `test_wf16_study_tuned_f_human`, `_study`, `_c_f`; from 390a's module, `FH_FILE`; plus the step-1 / `MAT-12` / `MAT-20` constants above, `W_B2` (`tests/validation/test_mat20_running_average.py:36`).
- **Nodes:** `test_wf1_fht_parse_and_refusals` and `test_wf1_fht_run_branch_no_mesh`.
- **Registry (lines added only):** the first node under `[names.load_study]`, the second under `[names.run_study_file]`.

**Anchors (asserted).**
- **(K) bitwise up to the solve.**
  - `_spec_diff(load_study(FHT_FILE), _study()) == []` and `_spec_diff(load_study(FHT_FILE), load_study(FH_FILE)) == []` (the thermal section adds nothing to the Study).
  - The parsed `RobinCooling` equals `RobinCooling(K_GEL, H_GEL, T_INF, source_scale=json["G_per_w"])` field by field by `float.hex`.
  - `TransientSpec` equals `TransientSpec(dt_s=T_STAR / N_MAIN, n_steps=N_MAIN, rho_c_j_per_m3_k=RHO_C, theta=0.5, startup_half_steps=2)` by hex.
  - The ball centre equals `json["c_star_m"]` by hex, `mass_g` is 10.0, and `window_s` equals `W_B2` by hex.
  - `cli.main(["check", FHT_FILE])` returns rc (0, 0) and prints `family: f_human_birdcage` (item 401's comparator).
  - **The JSON absent ⇒ FAIL, never skip**: a skipped node reads green.
- **(X) planted refusals, no mesh**, each `type(e) is StudyConfigError` with the named `.field`, through in-process one-key mutations of `FHT_FILE`'s mapping (no files):
  - `transient` removed ⇒ `thermal`;
  - `transient.dt_s: 0.0` ⇒ `thermal.transient.dt_s`;
  - `transient.n_steps: 0` ⇒ `thermal.transient.n_steps`;
  - `transient.theta: 1.5` ⇒ `thermal.transient.theta`;
  - `balls` removed ⇒ `thermal`;
  - `balls.mass_g: 0.0` ⇒ `thermal.balls.mass_g`;
  - `balls.centres: []` ⇒ `thermal.balls.centres`;
  - `sense: up` ⇒ `thermal.sense`;
  - an unknown `thermal.transient.dt` ⇒ that field.
  - **390a's own mutation (k / h / t_inf only) still raises at `thermal`**; RG9 re-runs that module unedited.
  - Through `cli.main` `run` / `check`: rc (2, 2) and nothing written.
- **(Run-0) the branch on a fake, no mesh.**
  - Inside `pytest.MonkeyPatch.context()`, `fem_em_solver.study._study_solve_tuned_quadrature` is replaced by 390a's counting fake (the same `SimpleNamespace`), and `TwoPortQuadratureDrive._temperature_port_study` by a counting fake returning a fixed `TunedTemperatureResult`-shaped namespace (three times, two balls).
  - `run_study_file(FHT_FILE)` calls each fake once per rank.
  - `summary["temperature"]` equals the fake's values exactly, and the written JSON equals the in-memory summary on rank 0.
  - With `thermal:` removed, the temperature fake is not called and the summary has no `"temperature"` key, so 390b's file path is unchanged.

**Negative controls (asserted; deterministic parse — rule (e) has no factor to label).**
- **(c1) dropped key:** `FHT_FILE`'s mapping minus `thermal.source_scale` diffs exactly in `source_scale`, which equals `RobinCooling`'s default 1.0 by hex. A loader that ignored `source_scale` passes (X) and fails (K) or (c1).
- **(c2) one ULP:** `source_scale` set to `repr(np.nextafter(G, 2G))` diffs exactly in `source_scale`.

**Windows.** Foreground. **The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped.** Bash `timeout: 660000`; `-s`; node IDs, never `-k`. Run `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` (it must read 0) before RG2, RG4, RG5, RG8 and RG12. Rule (i): no edit to `study.py`, `study_yaml.py`, the module, the YAML or the registry after the window that closes on them.
- **W00, the literals** (≈ 5 s; once the module exists, before the YAML is written; if the guard refuses the quoting, write the YAML from W0's (K) failure prints, which print each expected `repr` before asserting, and re-run W0):
  `scripts/testing/run_and_log.sh WF-1-s10a-W00-literals "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace timeout -k 30 60 python3 -c \"from tests.package.test_wf1_yaml_f_human_thermal import _literal_lines; print(_literal_lines())\"'"`
- **W0** (standard, `-n 2`, predicted 20–45 s):
  `scripts/testing/run_and_log.sh WF-1-s10a-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_f_human_thermal.py -s -v --tb=short'"`
- **Standing rule (c), named in full.** The edits reach `_parse_mapping` (every YAML route), `_parse_thermal` (every `thermal:` file), `run_study_file` and `_study_temperature` (every temperature route). Every command is verbatim from its last green log under a new label:
  - **RG1–RG8:** 390a's eight commands (this file, item 390a, "Standing rule (c)") verbatim, labelled `WF-1-s10a-RG1-wf1-w1` … `WF-1-s10a-RG8-wf2-four-port`. Measured 110 / 150 / 35 / 96 / 128 / 37 / 47 / 79 s (390a's outcome).
  - **RG9** (390a's module, unedited; `20261012T130057Z_WF-1-s7b-RGa.log:12`): `scripts/testing/run_and_log.sh WF-1-s10a-RG9-f-human "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_f_human.py -s -v --tb=short'"` (3 s).
  - **RG10** (item 401's glob over `wf1/*.yaml`, which now includes the new file; `20261012T172853Z_WF-1-s9-W1.log`): `scripts/testing/run_and_log.sh WF-1-s10a-RG10-check-family "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src timeout -k 30 60 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_check_family.py -s -v --tb=short'"` (2 s; (F) now reads 12 files).
  - **RG11** (`WF-15` steady temperature; `20261010T023137Z_OPS-70-W17-temperature-deg1.log:12`, 88 s): `scripts/testing/run_and_log.sh WF-1-s10a-RG11-wf15-temperature "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/package/test_wf15_study_temperature.py::test_wf15_temperature_spec_and_refusals tests/package/test_wf15_study_temperature.py::test_wf15_temperature_degree1_reproduces_mat10_record -v -s --tb=short'"`
  - **RG12** (`WF-15` transient, real build; `20261010T022152Z_OPS-70-W13-transient-real.log:12`, 198 s): `scripts/testing/run_and_log.sh WF-1-s10a-RG12-wf15-transient "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 400 mpiexec -n 2 python3 -m pytest tests/environment tests/package/test_wf15_study_transient.py -v -s --tb=short'"`
  - **RG13** (`WF-15` probe; `20261010T020736Z_OPS-70-W2-probe-n1.log:12`, 27 s): `scripts/testing/run_and_log.sh WF-1-s10a-RG13-wf15-probe "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 1 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf15_study_probe.py::test_wf15_probe_width_n1 -v -s --tb=short'"`
  - **Not in this slot:** 390b's heavy `-n 8` window, which is the only route through the edited F-human branch with a real solve. It runs in 431b as RGb.
- **W9, installed `femem check`** (it copies `wf1/*.yaml`, so the new file rides along): 390a's W9 command verbatim, labelled `WF-1-s10a-W9-installed-femem` (16 s).
- **Wgen / Wval / Wapi / Wdoc** (the new nodes move `load_study.md`, `run_study_file.md` and `docs/guide/validation.md`; commit them; the §9 "Every item" bullet): item 394's four commands verbatim, labelled `WF-1-s10a-Wgen` / `-Wval` / `-Wapi` / `-Wdoc`.

**Traps already paid for.**
- **`thermal.probe.ball_average` is not P17's estimator.** Lift the indicator forms; never substitute the lattice probe.
- **Never fill a default in the loader** (the (Deg) rule). An omitted `theta` / `startup_half_steps` / `window_s` / `sense` means the API's default.
- **PyYAML 1.1 floats need a dot** (`31.35`, `10.0`, `4180000.0`); (K) checks the bits.
- **The JSON is gitignored.** Read it on rank 0 and broadcast, and FAIL (never skip) when it is absent in a closing window.
- **Never call `_study_solve_tuned_quadrature`, `build_tuned_quadrature` or any mesh in 431a.** The fakes are the patch targets: `Study.sar` and `Study.temperature` reach the solve through `_tuned_result`, bypassing `Study.solve` (item 383's journal).
- **No field is added to an existing `study.py` dataclass** (`WF-16` step 7's field-set pin went red once). `TunedTemperatureResult` is new; `StudyRun.thermal` already exists.
- **A Python float inside a UFL form is a kernel literal.** Scale after the DG0 average (ruling 361-4).
- `cell_tags.values` and `assemble_scalar` are rank-local; `_rc_all` reduces rc.
- **On status 124 or a kill:** run `pgrep -c python3` first; at 0 only, the orphan-stem sweep labelled `WF-1-s10a-cache-sweep`. Never raise a timeout or change `-n`.
- Write ≈, never a tilde. No AED number.

**Scope.** Covered: the API and schema for a transient Robin-cooled temperature on the tuned copper F-human's own phantom, the per-step series, 10 g ball means and the six-minute maximum, all gated without a mesh. **Not claimed:**
- any solve (431b);
- a steady F-human temperature;
- perfusion, σ(T), the duty cycle or CEM43 through YAML;
- the solid σ = 800 F-human (refused);
- the ASTM box;
- degree 2;
- a stable schema.

**What a negative result means.** Report and stop. Park on `attempt/WF-1-step10a-<UTC>`, mark item 431a BLOCKED in the same commit (rule (d)), and journal it.
- (K) / (X) / (Run-0) red: fix in-slot inside the new lines only, then re-run W0.
- Any RG red: the edit leaked into an existing route. Revert `study.py` / `study_yaml.py`, mark BLOCKED, and file a known-issues entry if the diff is not the cause.

No `PREDICTED RED:` line.

**Status it can move.** `WF-1` stays ✅ (rule (j): its Done-when is step 1's). On green the row gains a step-10a clause (the tuned-F-human thermal schema, parse-only), `WF-14`'s envelope gains two nodes, and 431b becomes runnable. §10 U46 is **not** done until 431b lands.

**Spawn line:** "Use the implementer agent to execute `WF-1` §9 item 431a per `docs/planning/chunks/WF-1.md` (read the latest Review ruling paragraph first). The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run W00, W0, RG1–RG13, W9, Wgen, Wval, Wapi, Wdoc in order, foreground, `timeout: 660000`, pytest with `-s`, container `timeout -k 30`; before RG2, RG4, RG5, RG8 and RG12 run `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` and require 0. If W0 is not green by minute 32, park. No mesh and no solve in the new module; never call `build_tuned_quadrature` or `_study_solve_tuned_quadrature`; never import a `test_*` name; never add a field to an existing `study.py` dataclass; never use `thermal.probe.ball_average` for the 10 g mean. You are the executor — do not spawn agents, never background a harness window. On a red, park on `attempt/WF-1-step10a-*`, mark item 431a BLOCKED in the same commit, and report."

---

#### 431b — `WF-1` step 10b: `run_study_file("f_human_tuned_copper_thermal.yaml")` at `-n 8` reproduces P17's θ̄(t\*) and `MAT-20` step 3's six-minute maximum, and writes the series (chain step, solve)

**Outcome 431b 2026-10-14 (implementer, slot 02:00 UTC): DONE, green; no `src/` edit.** Inputs present (series JSON 4 422 bytes, node JSON 5 481 bytes, npz 20 812 bytes). W0 `20261014T020242Z_WF-1-s10b-W0.log` 5 passed / 5 skipped (knob unset), 3 s: source_scale `0x1.1352fb069c3cep+5` = G, centre = c\*_θ. `pgrep -c python3` 0, then W1 `20261014T020254Z_WF-1-s10b-W1.log` 10 passed on each of 8 ranks, 147 s (heavy; `run_study_file` 140.5 s, build 107.7 s, node pair 28.6 s): (A0) cells 369 387 (rel 0), solve count 1 at every node; (S_rep) 4.663415e-12; (V) 7.749e-14, n_sub 3 306; (B) 3.935e-16; (E) 1.608e-12; (T\*) θ̄(c\*_θ, 1 254 s) 4.926110849651e-02 K, rel 8.005e-14 (ball volume equal to the JSON's at 16 digits, rel 0); (Ser) times exact, ball series 2.588e-13, volume mean 3.998e-13; (6min) 4.277005459539e-02 K at 1 254 s, rel 1.423e-13, t_max exact, API vs operator 0.000e+00; (J) rank 0 differing keys `[]`. Controls, all asserted and green: (N1) M(c\*_θ) 1.321059e-01 = 1.321e+05× band; (N2) window 328.65 s 1.348831e-02 = 1.349e+04×; (N3) 2G miss 1.000000 (predicted 1.000). Printed: η 0.282134; volume-mean six-minute peak 3.070252443e-02 K. `pgrep` 0, then RGb `20261014T020529Z_WF-1-s10b-RGb-s7b.log` 10 passed, 222 s: every 390b reading equal at its printed digits ((S_rep) 4.434473e-12, (S10) 2.857e-10, (B) 1.101e-07, (iii-a) 0.6413 %, (M_B) 1.9647e-03, (N1) 9.5476e-01, (N2) 4.266875e-02). Pages: Wgen 1 s (`…020920Z`), Wval 190 nodes 2 s (`…020922Z`), Wapi 22 passed 15 s (`…020929Z`), Wdoc 0 dead 1 s (`…020946Z`). Deviations: (V)'s V_sub is assembled in-test on `study._tag_submesh` of the run's mesh (the result carries no volume); `HOLE_CELL_RECORD` is imported from `test_wf12_f_human_copper_matrix.py` (as 390b), `CELL_COUNT_BAND` from the prerequisite module. U46 is done.

**Depends on item 431a landing on `main` green; if it did not, skip.** It is test-only, with no `src/` edit (as 390b).
- **Executor:** implementer.
- **Tier:** heavy, `-n 8` (the series' width, `json["ranks"]`, and 390b's `RECORD_RANK_WIDTH`), complex, container `timeout -k 30 590`.
- **Class:** chain step (solve), median **12 measured** min. **Cost: ≈ 45 predicted** — 15 fixed + ≈ 12 code + W0 ≈ 0.3 + W1 ≈ 3–4 + RGb ≈ 4–5.5 + pages ≈ 0.6 + ≈ 6 spare.
- **Price.**
  - W1 predicted 150–230 s: build ≈ 114 s and node pair ≈ 27 s (`20261012T125631Z_WF-1-s7b-W1.log`), the thermal half ≈ 5 s (`20261011T003520Z_MAT-15-step5-W1.log`: thermal 4.3 s, 40 balls 0.8 s), and start-up ≈ 20 s; there are no SAR / B outputs in the file.
  - RGb (390b's W1) measured 243 s.
  - Two heavy windows run serially, never side by side.
  - **If W1 has not started by minute 35, park.**
- **Inputs (gitignored; checked present on the host 2026-10-13):** `logs/mat15_s5/theta-series-ccw.json`, `logs/wf10/h9p-node-118.2437.json` (5 474 bytes) and `logs/wf12/h12c-cu-s32-a2.npz`. **Absent ⇒ stop, BLOCKED.** Never run the step-5 module's series knob, `build_tuned_quadrature` or item 358 in-slot.

**New `tests/package/test_wf1_yaml_f_human_thermal_solve.py`.**
- **Knob and width.** Opt-in knob `FEM_EM_WF1_S10B=1`: unset ⇒ the solve nodes skip, stated; any other value raises (390b's `_knob` pattern). The width is asserted `== 8` when the knob is set. `@complex_only`.
- **One module fixture.** It runs, inside `pytest.MonkeyPatch.context()`, a counting call-through wrapper on `fem_em_solver.study._study_solve_tuned_quadrature`, then `run = run_study_file(FHT_FILE, output_dir=LOGS / "f_human_tuned_thermal", comm=MPI.COMM_WORLD)`. The JSONs are read on rank 0 and broadcast **before** the run (390b / 361's order).
- **Docstring.** The first paragraph backticks only `WF-1`, `MAT-15` and `MAT-20`.
- **Imports (names only; never a `test_*`):**

  | Name | Source |
  |---|---|
  | `REUSE_REPRODUCTION_RTOL` | `tests/validation/test_port_birdcage_ring_matrix.py:171` (1e-6; **not** `test_port19_factor_reuse.py:82`'s 1e-12) |
  | `_residual` | `tests/validation/test_port_circuit_layer_field.py:142` |
  | `NODE_JSON`, `_cx` | as 390b imports them |
  | `HOLE_CELL_RECORD`, `CELL_COUNT_BAND` | `tests/mesh/test_birdcage_port_sheet_prerequisite.py:57` |
  | `V_PHANTOM_H9P_RECORD_M3`, `SERIES_JSON` | `tests/validation/test_mat15_f_human_tuned_heating.py:139`, `:134` |
  | `V_BAND`, `D_CONTROL_FLOOR` | `tests/validation/test_mat15_coil_driven_heating_phantom.py:98`, `:102` |
  | `B_BAND` | `tests/validation/test_mat10_sar_to_heat_sphere.py:61` |
  | `E_BAND` | `MAT-12`, `:76` |
  | `W_B2`, `BAR` as `AVG_BAR` | `tests/validation/test_mat20_running_average.py:36`, `:34` |
  | `PHANTOM_CELL_TAG` | `tests/validation/test_port_birdcage_lumped_column.py:101` |
  | `PHANTOM_RHO_KG_PER_M3` | `tests/validation/test_birdcage_b1_plus_map.py:139` |
  | `mean_sar` | `src/fem_em_solver/post/sar.py:57` |
  | `running_average` | `src/fem_em_solver/thermal/dose.py:108` |
  | `FHT_FILE` | 431a's module |

- **Nodes:**
  - `test_wf1_fht_solve_knob_and_inputs` (no mesh);
  - `test_wf1_fht_cells_once_and_source`;
  - `test_wf1_fht_ball_reproduces_p17`;
  - `test_wf1_fht_series_and_six_minute`;
  - `test_wf1_fht_controls`;
  - `test_wf1_fht_summary_json`.
- **Registry (lines added only):** the solve nodes under `[names.run_study_file]` and `[names.TwoPortQuadratureDrive]`.

**Anchors (asserted; all imported).**
- **(A0) cells and one solve.** `summary["cells"] == HOLE_CELL_RECORD` (369 387) within `CELL_COUNT_BAND`. The solve count is 1 after every read (the thermal half reuses `study.result`). `run.result is run.study.result`.
- **(S_rep) source wiring.** `_residual(summary["s_node"], node_json s_node) ≤ REUSE_REPRODUCTION_RTOL`. Prior on the same route and fixture: 4.287849e-12 (`20261012T125631Z_WF-1-s7b-W1.log`).
- **(V) fixture tie.** \|V_sub / `V_PHANTOM_H9P_RECORD_M3` − 1\| ≤ `V_BAND`, with V_sub from the result's submesh. Prior 7.749e-14 (`20261011T003520Z_MAT-15-step5-W1.log`). `n_sub == json["n_sub"]` (3 306).
- **(B) source transfer.** \|`source_power_w` − G·P_phantom,ccw\| / (G·P_phantom,ccw) ≤ `B_BAND`. P_phantom,ccw is `mean_sar(e_ccw, sigma, rho=PHANTOM_RHO_KG_PER_M3, cell_tags, subdomain_ids=PHANTOM_CELL_TAG)["dissipated_power_w"]` on the run's field, and G = `json["G_per_w"]`. P17's prior on the builder route: 0 (`…MAT-15-step5-W1.log`, (B) ccw).
- **(E)** `max_balance_residual ≤ E_BAND`. P17 read ≤ 3.3e-10 over seven runs.
- **(T\*) P17's ball mean.** \|θ̄(c\*_θ, t\*) / `json["theta_ball10g_k"][-1]` − 1\| ≤ `REUSE_REPRODUCTION_RTOL`. Predicted ≤ 1e-10: the field reproduces at ≈ 4e-12 by (S_rep), the thermal arithmetic is P17's lifted, and the ball volume should equal `json["ball_volume_m3"]` (printed at `.12e`). The logged 4.926110850e-02 K is printed beside it.
- **(Ser) the written series.** The summary's `times_s` equal `json["times_s"]` to 1e-12 s, and the ball series and the volume-mean series each equal the JSON's elementwise within `REUSE_REPRODUCTION_RTOL` of the series maximum (sample 0 is θ = 0 exactly). The bundle's JSON equals `run.summary` exactly on rank 0.
- **(6min) `MAT-20` step 3's six-minute figure.**
  - \|peak / peak_json − 1\| ≤ `REUSE_REPRODUCTION_RTOL`, where peak_json is `running_average(json["times_s"], json["theta_ball10g_k"], window=W_B2)`'s maximum, and t_max == t_max_json exactly. The logged 4.277005e-02 K at 1 254 s is printed beside it.
  - The summary's peak equals an in-test `running_average` on the summary's own series ≤ `AVG_BAR` (an API-vs-operator identity).

**Negative controls (rule (e)).**
- **(N1) the wrong ball — asserted, backed.** θ̄_ccw at the feed-bisector image M(c\*_θ) misses θ̄(c\*_θ) by > `REUSE_REPRODUCTION_RTOL`. M(c\*_θ) is `_mirror_xy(c*, run.result.az_mirror)` (`tests/validation/test_birdcage_b1_quadrature.py:203`), computed through `run.study.temperature(...)` with that one centre; the result is reused, so there is no re-solve. Backing: P17's "mis-paired" reading, the same comparison on the same fixture, was 13.21 % (`20261011T003520Z_MAT-15-step5-W1.log`), ≈ 1.3e5× the band. The factor is printed. The ceiling is not binding.
- **(N2) the window one sample short — asserted, backed.** `running_average(..., window=W_B2 − Δt)` on the summary's ball series misses peak_json by > `REUSE_REPRODUCTION_RTOL`. `MAT-20` step 3's (Sh) read 4.335e-02 K against 4.277e-02 K on the same series (`20261011T183717Z_MAT-20-step3-W2.log`): 1.36 %, ≈ 1.4e4× the band.
- **(N3) the doubled source — the miss asserted against a floor, the factor predicted.** `run.study.temperature` with `source_scale` = 2G, no re-solve. θ̄(t\*) is predicted, by linearity, to miss the JSON by exactly 1.000; this comparison has not been measured on this fixture, so 1.000 is printed. Asserted only > `D_CONTROL_FLOOR` (1e-3).

**Printed only:** η = P_phantom / P_acc (P17 0.2821); seconds per phase; the summed RSS; the ball volume against the JSON's; the volume-mean six-minute peak against MAT-20's 3.070e-02 K.

**Windows.** Foreground. **The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped.** Bash `timeout: 660000`; `-s`; never `-k`. Under 10 min each, so no durable capture. Before W1 and before RGb, `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0.
- **W0** (no mesh, knob unset, ≈ 20 s): `scripts/testing/run_and_log.sh WF-1-s10b-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_f_human_thermal_solve.py -s -v --tb=short'"`
- **W1** (heavy, `-n 8`, predicted 150–230 s): `scripts/testing/run_and_log.sh WF-1-s10b-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_WF1_S10B=1 timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_f_human_thermal_solve.py -s -v --tb=short'"`
- **RGb** (standing rule (c) for 431a's edit to the F-human branch of `run_study_file`; this is the only real-solve route through it): 390b's W1 verbatim (this file, item 390b, "Windows"; `20261012T125631Z_WF-1-s7b-W1.log:12`, 243 s), labelled `WF-1-s10b-RGb-s7b`: `scripts/testing/run_and_log.sh WF-1-s10b-RGb-s7b "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_WF1_S7B=1 timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_f_human_solve.py -s -v --tb=short'"`. Every 390b reading must equal its prior at the printed digits, with residual figures at rounding (item 380's standard).
- **Wgen / Wval / Wapi / Wdoc:** item 394's four commands verbatim, labelled `WF-1-s10b-Wgen` / `-Wval` / `-Wapi` / `-Wdoc`. Commit the pages.

**Traps already paid for.**
- **One solve.** Never build a fresh `Study` for a control. The controls call `run.study.temperature` on the kept result.
- **Every ball and every reduction is collective on all 8 ranks.** The phantom submesh may be owned by only some ranks (`MAT-19` D4), and a guarded collective deadlocks (item 227).
- **Never apply G twice:** the YAML's `source_scale` *is* G.
- `cem43` is not in scope (no dose through YAML). Passing θ for T is the trap item 380 named.
- **Status 124 or a kill:** run `pgrep -c python3` first; at 0 only, the orphan-stem sweep labelled `WF-1-s10b-cache-sweep`. A killed window is not a result; never raise the timeout or change `-n`.
- Rule (i). Write ≈, never a tilde. No AED number.

**Scope.** Covered: H9′'s tuned copper F-human (degree-1 EM, 369 387 cells, κ off, r = 0.2 Ω on 32 gaps, the two-port 90° drive at 118.2437 pF, 64 MHz) heating its 0.015 m cylindrical gel phantom through `femem run`'s Python entry point at `-n 8`: the ccw ball mean at c\*_θ and the volume mean per step, per accepted W (ruling 361-1, through the transcribed G), and the 360 s running maximum. **Not claimed:**
- the console script or `--ranks` on this file;
- the cw sense, or the lattice;
- CEM43;
- an absolute temperature or a compliance figure;
- the ASTM box;
- degree 2 (P22);
- perfusion;
- an `ANS-9` reading.

**What a negative result means.** Report and stop. Park on `attempt/WF-1-step10b-<UTC>`, mark item 431b BLOCKED in the same commit (rule (d)), and journal it.
- (A0) count > 1: the thermal half re-solves — a finding against the branch.
- (S_rep) / (V) red: not P17's fixture. Read nothing thermal; file a known-issues entry.
- (B) / (E) red: the lifted source or thermal path. Known-issues.
- **(T\*) / (Ser) / (6min) red with (S_rep), (V), (B) and (E) green:** the lifted ball forms differ from P17's. Print the ball volume against the JSON's first; never re-band and never re-record the JSON.
- RGb red: 431a's branch edit broke 390b's path. Revert is the review's call; known-issues.
- A control inside its bar: a blind comparator. Stop.

No `PREDICTED RED:` line.

**Status it can move.** `WF-1` stays ✅ (rule (j)). On green of both parts:
- the row gains a **step-10 record** clause: "`thermal:` on the tuned copper F-human YAML: `femem run` at `-n 8` reproduces P17's 10 g ball-mean θ(t\*) (<x>) and `MAT-20` step 3's six-minute maximum (<y>) at `REUSE_REPRODUCTION_RTOL` and writes the 41-step series; one solve";
- §10 U46 is done;
- `WF-14`'s envelope gains nodes;
- (9a)'s "temperature via `thermal/`" covers the tuned F-human.

`MAT-15` and `MAT-20` stay ✅ and gain nothing.

**Spawn line:** "Use the implementer agent to execute `WF-1` §9 item 431b per `docs/planning/chunks/WF-1.md` — skip it unless item 431a is DONE on `main` (read the latest Review ruling paragraph first). The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run W0, W1, RGb, Wgen, Wval, Wapi, Wdoc in order, foreground, `timeout: 660000`, pytest with `-s`, container `timeout -k 30`; before W1 and before RGb run `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` and require 0; two heavy windows never run at once. If `logs/mat15_s5/theta-series-ccw.json`, `logs/wf10/h9p-node-118.2437.json` or `logs/wf12/h12c-cu-s32-a2.npz` is absent, stop (BLOCKED). No `src/` edit; never call `build_tuned_quadrature`; never import a `test_*` name. You are the executor — do not spawn agents, never background a harness window. On a red, park on `attempt/WF-1-step10b-*`, mark item 431b BLOCKED in the same commit, and report."

### Review ruling 2026-10-13 03:00 (daily) on the item 431a / 431b drafts

Queued as drafted, split in two (the 390a / 390b precedent): 431a the schema and branch with no mesh, 431b the one `-n 8` solve, serial on 431a. **A new file `tests/package/wf1/f_human_tuned_copper_thermal.yaml`** instead of a `thermal:` section in `f_human_tuned_copper.yaml` — accepted provisionally (editing the closed file would reopen items 390a / 390b's records); the weekly ratifies. 390a's thermal-mutation gate (`test_wf1_yaml_f_human.py:215–233`) stays unedited and green. Trap: `thermal.probe.ball_average` is not P17's ball-mean estimator — import P17's.

### §9 item 442 — 2026-10-14 (§10 chain U49, daily — provisional)

**Outcome 2026-10-14 (implementer, slot 11:00 UTC): DONE, green.** `src/` edit in `study_yaml.py` only (`_TOP_KEYS` gains `tune`; `_parse_tune`; `_tune_summary`; one line each in `_parse_mapping` and the four-port branch of `run_study_file`; docstring bullet); `_load`'s 4-tuple unchanged; no κ; the YAML says `degree: 1`. W0 `20261014T110309Z_WF-1-s11-W0.log` 7 passed, 3 s: (K) both `_spec_diff` `[]`, parsed outputs `{'tune': {'k': 1}}` (int), control has no `tune`, `tune: {}` → `{"tune": {}}`; (X) the tracked file names `tune.k`, cli run / check rc (2, 2), nothing written; the seven mapping variants name `tune.k` ×3 / `tune.mode` / `tune` ×3; 0 build calls; unknown-top-level message lists `tune`; `test_wf1_check_family.py` green with the new good file. `pgrep -c python3` 0, then W1 `20261014T110326Z_WF-1-s11-W1.log` 10 passed, 62 s (YAML build 29.7 s, mesh 22.7 s, sweep 6.5 s; Python build 28.6 s): (R0) cells 80 798; class means 5.052e-11 / 3.878e-11 / 5.875e-11 (prior equal at every printed digit); (I=) C `0x1.2373183ff8634p-34` both, every `z_k` part, residual, JSON read-back by hex; k = 1, f = 64 MHz by hex; no second solve; (I-Py) rel 0.000e+00 (hex equal, printed only); (a) `C_mode` 6.626796421e-11 F, gap −2.976552e-04 (predicted −2.977e-4, inside [−3.0e-4, −2.95e-4]; item 429's −2.976552e-4); (c1) k = 2 9.873468900e-11 F, +48.95 %, 489.5× the band. RG1 `20261014T110438Z_WF-1-s11-RG1.log` (the nine `[names.load_study]` nodes the registry lists — the brief's "ten" counted one more than are registered) 13 passed, 3 s; `pgrep` 0, RG2 `20261014T110453Z_WF-1-s11-RG2.log` 5 passed, 72 s. Pages: Wgen 2 s (`…110628Z`), Wval 203 nodes 3 s (`…110636Z`), Wapi 22 passed 16 s (`…110645Z`), Wdoc 0 dead 1 s (`…110707Z`). Not re-run (untouched branches): the sphere / thermal / F-human solve nodes. U49 is done (provisional).

*Drafted by a general-purpose agent for the 2026-10-14 03:00 daily review; the review reviews it.* Licence: §9 rule (6), a provisional step of §10 chain U ("each step reproduces a gated result or asserts a gated identity through the API, within its unmoved band"), tagged "U49 (daily, 2026-10-14 — provisional)"; the next weekly ratifies or reverts. Origin: §9 item 429 (`WF-16` step 6(a), `ac947c9c`) added `Study.tune(k, *, comm)` and `TuneResult` (`src/fem_em_solver/study.py:998–1060`). `femem run` cannot reach it yet. **Target: reproduction only. No band, tolerance, Done-when, record or validation target is set or moved; every band is imported by name.** Owner: `WF-1` (steps 1–10 are items 256 / 293 / 295 / 343 / 344 / 345 / 390a–b / 394 / 401 / 431a–b). Rule (j)'s reading of "bit-identical" binds: bitwise up to the solve, imported bands after it (ratified by the 2026-10-05 weekly).

**`WF-1` step 11: a `tune:` section makes `femem run` call `Study.tune` on the F-small PEC hole Study and write the result. The written `C_mode` equals `TH-17`'s record within `REDUCTION_BAND`, equals the same Study's `tune()` bitwise, and equals an independent Python Study's `tune()` within `LEG_D0_REPRODUCTION_BAND`.**
- **Executor:** implementer.
- **Tier:** heavy by ceiling (container `timeout -k 30 300`, item 429's W1 ceiling), standard by prediction. `-n 2` (`ANS6_EXACT_STRIP_RECORD_RANK_WIDTH`), complex.
- **Class:** chain step (solve), ≈ 18 measured min. **Cost: ≈ 42 predicted** = 15 fixed + ≈ 14 code (≈ 30 `src` lines, module ≈ 180 lines, two YAML files, registry lines) + W0 ≈ 0.3 + W1 ≈ 1.7 + RG1 ≈ 0.3 + RG2 ≈ 2.5 + Wgen / Wval / Wapi / Wdoc ≈ 0.7 + ≈ 7 row / journal / spare. **If W1 has not started by minute 35, park** (rule (d)).
- **Family:** `WF-1` step 11, attempt 1 of 4.
- **Depends on:** nothing. Item 429 is on `main` (`ac947c9c`). Items 431a / 431b also edit `study_yaml.py`, in the F-human branches. On a rebase keep both sides, and re-run Wgen / Wval / Wapi / Wdoc after it.

**What the code says (read by this drafter on `b6d366b5`, not run).**
- **The top-level key list.** `_TOP_KEYS` (`src/fem_em_solver/study_yaml.py:103–105`) lists the eight top-level keys. `_parse_mapping` (`:716–775`) refuses any other key at that key's own name (`:731–735`). It validates study → geometry → material → excitation → conductor → thermal → outputs (`:739–769`) and returns `(study, kind, outputs, cooling)`.
- **The 4-tuple must not change.** At least 12 call sites under `tests/package/` unpack `_load(...)` / `_parse_mapping(...)` as four names (`test_wf1_yaml_probes.py:75`, `test_wf1_yaml_copper_sar_b1.py:109`, `test_wf1_yaml_f_human.py:155`, `test_wf1_yaml_thermal.py:111`, …), as do `run_study_file` (`:985`) and `cli.py`'s `check` branch (`_load(...)[:2]`).
- **So the parsed `tune` spec travels in the `outputs` dict, under the key `"tune"`.** The YAML key is top-level `tune:`, as commissioned. This is disclosed, and it is question 1 below.
- **`run_study_file`'s four-port branch** (`:1011–1023`) calls `study.solve(comm)` (`:993`), then writes the Touchstone and the `sar` / `b1_plus` blocks. `Study.tune` reuses `self.result` (`study.py:1040–1043`); `_study_solve_ports` sets it (`study.py:655`). So a `tune` block read after `solve` costs no second solve.
- **The PEC hole file already exists.** `tests/package/wf1/f_small_hole_pec_control.yaml` (item 343's control) parses to `_py_copper("pec")`: `test_wf1_yaml_copper.py:113–118` asserts that its diff against `_py_copper()` is exactly `["conductor"]`.

**Files.**
- **`src/fem_em_solver/study_yaml.py`**, insertions plus one edited line:
  - The edited line: `_TOP_KEYS` gains `"tune"`, last.
  - A new `_parse_tune(sec, geometry_kind, source) -> dict`. Its order follows `_parse_thermal`'s (`:465–496`): geometry first, then unknown keys, then types and values.
    - A `sphere_in_box` or `f_human_birdcage` geometry raises `StudyConfigError(source, "tune", "tune is defined on a four_port_birdcage only, not …")`. This mirrors `Study.tune`'s `TypeError` (`study.py:1030–1032`).
    - `_unknown(sec, ("k",), "tune", source)`.
    - `k` is optional. Omitted means `Study.tune`'s own default: the key is not filled in, as for `study.degree`. When present it is `_integer(sec["k"], "tune.k", source)`, which refuses a bool and the float `1.0`. A value outside {1, 2} raises at `tune.k` with "must be 1 or 2".
    - It returns `{}` or `{"k": k}`.
  - In `_parse_mapping`, after the `outputs` block: `if "tune" in mapping: outputs = {**outputs, "tune": _parse_tune(_section(mapping, "tune", source), kind, source)}`. A non-mapping `tune: 1` raises at `tune` through `_section`.
  - A new `_tune_summary(study, spec, comm) -> dict`: `t = study.tune(**spec, comm=comm)`, returning `{"k": t.k, "frequency_hz": t.frequency_hz, "c_mode_farad": t.c_mode_farad, "z_k": [[re, im] × 4], "off_diagonal_residual": t.off_diagonal_residual}`, all plain floats.
  - In `run_study_file`'s four-port branch, after `b1_plus`: `if "tune" in outputs: summary["tune"] = _tune_summary(study, outputs["tune"], comm)`.
  - The module docstring's schema list gains one `tune` bullet citing `WF-1` step 11.
  - Nothing else changes. `StudyRun`, `load_study`'s signature and `cli.py` are untouched.
- **`tests/package/wf1/f_small_hole_pec_tune.yaml`**: `f_small_hole_pec_control.yaml` verbatim, with its header line replaced by a header citing this step, plus `tune:` / `  k: 1`.
- **`tests/package/wf1/malformed_tune_k_zero.yaml`** (header line): the same file with `k: 0`.
- **New `tests/package/test_wf1_yaml_tune.py`.**
  - The docstring's first paragraph backticks **only** `WF-1`, `WF-16` and `TH-17` (§7 rows), never `U49` or an item number (the (F3) rule).
  - `@complex_only` on every solve node. The width is asserted in the fixture, never skipped: `MPI.COMM_WORLD.size == ANS6_EXACT_STRIP_RECORD_RANK_WIDTH`.
  - **Imports: module-level names only, never a `test_*` name.**

    | Name | Source |
    |---|---|
    | `LOGS`, `WF1`, `_fresh_dir`, `_mapping`, `_rc_all`, `_spec_diff` | `tests/package/test_wf1_study_yaml.py:85–170` |
    | `PEC_CONTROL_FILE`, `_py_copper` | `tests/package/test_wf1_yaml_copper.py:85, :98` |
    | `TH17_STEP3_C_MODE_RECORD_F` | `tests/package/test_wf16_study_tune.py:71`, imported, never re-transcribed |
    | `ANS6_EXACT_STRIP_CLASS_MEANS_P1`, `ANS6_EXACT_STRIP_HOLE_CELL_RECORD`, `ANS6_EXACT_STRIP_RECORD_RANK_WIDTH`, `CLASSES`, `_class_means` | `tests/validation/test_ans6_exact_strip.py`, as item 429 imports them |
    | `_rel_means` | `tests/package/test_wf16_study_copper.py:76` |
    | `REDUCTION_BAND` | `tests/validation/test_port_lumped_rlc_termination.py:88` |
    | `LEG_D0_REPRODUCTION_BAND` | `tests/validation/test_port_birdcage_four_port.py:165` |
    | `EXACT_IDENTITY_RTOL` | `tests/validation/test_birdcage_sar_mass_averaged.py:163` |
    | `FREQUENCY_64_HZ` | `tests/validation/test_lossy_sphere_fullwave.py:94` |

  - **The module fixture `runs`, in this order:**
    1. `_fresh_dir(LOGS / "f_small_hole_pec_tune", comm)`, then `run = run_study_file(TUNE_FILE, output_dir=…, comm=comm)`, then `summary.json` read back.
    2. Record `res0 = run.study.result`, then `api1 = run.study.tune(k=1, comm=comm)`. Assert `run.study.result is res0` (no second solve).
    3. `py = _py_copper("pec")`, `py.solve(comm)`, `pyt = py.tune(k=1, comm=comm)`. This is item 429's route, argument for argument.
    4. Print each build's `mesh_time` / `sweep_time`, the wall time of `run_study_file` and `summary["elapsed_s"]`.
  - **The nodes.** Assertions go in the node bodies; the validation page traces bands in node bodies only.
    - `test_wf1_s11_tune_parse_and_refusals`: no mesh, `monkeypatch`, does not request `runs`.
    - `test_wf1_s11_r0_fixture_tie`
    - `test_wf1_s11_yaml_equals_same_study_api`
    - `test_wf1_s11_yaml_equals_python_route`
    - `test_wf1_s11_c_mode_equals_th17_record`
    - `test_wf1_s11_c1_mode_index_control`
- **The registry.** Lines are added only, each block under a `# WF-1 step 11 (§9 item 442)` comment. The parse node goes under `[names.load_study]` and `[names.StudyConfigError]`. The five solve nodes go under `[names.run_study_file]`.

**Anchors (asserted; bands imported).**
- **(K) The parse is bitwise and needs no mesh.**
  - `_spec_diff(load_study(TUNE_FILE), _py_copper("pec")) == []`.
  - `_spec_diff(load_study(TUNE_FILE), load_study(PEC_CONTROL_FILE)) == []`: `tune` adds nothing to the Study.
  - `study_yaml._load(TUNE_FILE)[2] == {"tune": {"k": 1}}`, with `type(k) is int`.
  - `_load(PEC_CONTROL_FILE)[2]` has no `"tune"` key.
  - The mapping with `tune: {}` parses to `{"tune": {}}`.
- **(X) Every refusal names its field and builds no mesh.** `FourPortBirdcageGeometry.build` is stubbed and counted, and the count stays 0.
  - The tracked file `malformed_tune_k_zero.yaml` raises `StudyConfigError` with `.field == "tune.k"`. `cli.main(["run", f, "--output-dir", out])` and `["check", f]` both give rc `(2, 2)` through `_rc_all`, and `out` does not exist on any rank (`MPI.LOR`).
  - The mapping variants raise in-process, through `_study_from_mapping`:

    | Variant | Field named |
    |---|---|
    | `k: 3` | `tune.k` |
    | `k: 1.0` | `tune.k` |
    | `k: true` | `tune.k` |
    | `mode: 1` | `tune.mode` |
    | `tune: 1` | `tune` |
    | `tune: {k: 1}` added to `sphere_mat4.yaml`'s mapping | `tune` |
    | `tune: {k: 1}` added to `f_human_tuned_copper.yaml`'s mapping | `tune` |

  - The unknown-top-level-key message now lists `tune`. Print it; no existing test asserts that text (the drafter grepped `tests/package/` for "allowed").
- **(R0) The fixture tie.**
  - `summary["cells"] == ANS6_EXACT_STRIP_HOLE_CELL_RECORD` (80 798).
  - `_class_means(run.result.s_matrix)` against `ANS6_EXACT_STRIP_CLASS_MEANS_P1["64"]["pec"]`, per class, ≤ `LEG_D0_REPRODUCTION_BAND`. Prior on the Python route: 5.052e-11 / 3.878e-11 / 5.875e-11 (`20261014T004318Z_WF-16-s6-W1.log`).
- **(I=) The summary is the same Study's `tune()`, bitwise (up to the solve: the same S, the same arithmetic).**
  - `summary["tune"]["c_mode_farad"].hex() == api1.c_mode_farad.hex()`.
  - Every `z_k` entry is equal by `float.hex` on both parts.
  - `off_diagonal_residual` is equal by hex.
  - The JSON read-back's `tune` block equals `run.summary["tune"]` by hex (`json.dump` writes `repr`, which round-trips).
  - `summary["tune"]["k"] == 1` and `frequency_hz.hex() == FREQUENCY_64_HZ.hex()`.
- **(I-Py) The YAML route equals item 429's Python route across two builds:** |summary C / `pyt.c_mode_farad` − 1| ≤ `LEG_D0_REPRODUCTION_BAND`.
  - Predicted ≤ 1e-14: item 256's (I) on S across two builds read 9.110e-16.
  - **Hex equality is printed, never asserted.** Rule (j) does not compare a solve output bitwise across builds (`MAT-1` item 81; `MAT-5` item 237).
- **(a) The record:** |summary C / `TH17_STEP3_C_MODE_RECORD_F` − 1| ≤ `REDUCTION_BAND`.
  - Predicted −2.977e-4, inside [−3.0e-4, −2.95e-4]. Item 429 measured −2.976552e-4 on this fixture at this degree and width, printed beside.
  - Predicted green, so there is no `PREDICTED RED` line.

**Negative controls (rule (e)).**
- **(c1) The mode index through the YAML route — asserted, backed.**
  - The mapping variant `tune: {k: 2}` parses (via `_parse_mapping`) to `{"tune": {"k": 2}}`. `_tune_summary(run.study, {"k": 2}, comm)` then runs on the already solved YAML Study: the same function `run_study_file` calls, with no re-solve, and `run.study.result is res0` is asserted after it.
  - Its C against the record must exceed `REDUCTION_BAND`.
  - Backing: item 429's (c1) measured this comparison on this fixture, degree and width: 9.873468900e-11 F, +48.95 %, 489.5× the band (`20261014T004318Z_WF-16-s6-W1.log`).
  - No ceiling below 100 % binds. The factor is printed.
  - This is the control for a summary that ignores the parsed `k`.
- **(c0) The silently dropped section is caught by construction.**
  - (K)'s `{"tune": {"k": 1}}` and (I=)'s `"tune" in summary` together catch a parse that keeps the key but never reaches `run_study_file`.
  - The file without `tune:` has no `"tune"` in its parsed outputs.
- **Why there is no wrong-fixture control.** Item 429's (c2) already separates fixtures on the same `tune` lines: the solid σ = 800 F-small reads −6.98 %. The code new here is the parse and the summary, and (c1), (I=) and (K) cover their blind modes. A third build (≈ 56 s) would buy nothing new. For the review: question 2.

**Windows.** Foreground. **The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped.** Bash `timeout: 660000`; pytest with `-s`; `tests/environment/test_complex_mode.py` first; node IDs, never `-k`. Before W1 and before RG2, `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0. Two heavy windows never run at once.
- **W0** (no mesh, `-n 2`, ≈ 10 s). It covers the new parse node and `test_wf1_check_family.py`: that module's (F) globs every non-`malformed_` YAML in `wf1/` (`:32`), so the new good file joins it, and `check` must print `family: four_port_birdcage` with rc (0, 0).
  `scripts/testing/run_and_log.sh WF-1-s11-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_tune.py::test_wf1_s11_tune_parse_and_refusals tests/package/test_wf1_check_family.py -v -s --tb=short'"`
- **W1** (the module; heavy by ceiling, `-n 2`; predicted 60–120 s):
  `scripts/testing/run_and_log.sh WF-1-s11-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 300 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_tune.py -v -s --tb=short'"`
  - **Price.** Item 429's W1 took 66 s for two builds, the PEC hole build alone 29.1 s (mesh 22.0 s, sweep 6.6 s, `20261014T004318Z_WF-16-s6-W1.log:944`). Here there are two PEC hole builds (the YAML one and the Python one), ≈ 58 s, plus start-up and `test_complex_mode.py` ≈ 20 s, ≈ 80 s in all.
- **RG1 and RG2: standing rule (c)** (an existing line of `study_yaml.py` moves, `_TOP_KEYS`, and `_parse_mapping` / `run_study_file` gain lines).
  - **RG1** (no mesh, `-n 2`, standard, `timeout -k 30 180`): every `[names.load_study]` node in `tests/package/validation_envelope.toml:238–255`, the existing ten, by node ID, under the label `WF-1-s11-RG1`. Item 401's re-gate ran 13 such nodes in 3 s (`20261012T172909Z_WF-1-s9-RG.log`).
  - **RG2** (the four-port branch of `run_study_file`; heavy by ceiling, `-n 2`, `timeout -k 30 400`): `tests/package/test_wf1_study_yaml.py::test_wf1_cli_f_small_reproduces_wf16`, under the label `WF-1-s11-RG2`. Its W2 took 142 s in item 256 (`20261006T003504Z`…); it is predicted 70–150 s now.
  - The sphere / thermal / F-human branches are not touched. Disclose that their solve nodes were not re-run.
- **Wgen / Wval / Wapi / Wdoc** (ruling Q3 of 2026-10-12: the registry lines move the API pages and `docs/guide/validation.md`; commit both). Item 394's four commands verbatim (this file, item 394, "Windows"), labelled `WF-1-s11-Wgen` / `-Wval` / `-Wapi` / `-Wdoc`; `0 dead`.

**Traps already paid for.**
- **Never compare a solve output bitwise across two builds** (rule (j)). Only (I=) is bitwise, because it compares one S with itself.
- **`_load`'s 4-tuple is load-bearing** for 12+ callers. Never widen it.
- **`Study` defaults to degree 2.** The YAML file says `degree: 1`, and (K) checks it against `_py_copper("pec")` (`degree=1`). Without it, (a) reads ≈ −4.0 % (item 374, measured at degree 2 by `ANS-7` item 217).
- **κ is never passed** (sheets-and-κ ruling). The PEC route is `conductor` absent, with `as_hole: true`.
- **PyYAML 1.1** (`WF-1` step 1): `k: 1` is an int, `k: 1.0` a float (refused), `true` a bool (refused). Floats in the file need a dot.
- **Rank-local values.** rc is reduced (`_rc_all`); `out`'s existence is reduced (`MPI.LOR`). `summary["cells"]` comes from the Study's reduced path. `tune` is NumPy on the gathered S, the same on every rank.
- **The validation page reads the docstring's first paragraph for chunk IDs, and bands from node bodies** (item 375). Keep each band in the node that asserts it.
- **Stubs are function-scoped.** The parse node uses `monkeypatch` and never requests `runs` (item 394's trap).
- **On status 124 or a kill:** run `pgrep -c python3` first. At 0 only, run the rubric's orphan-stem sweep as one bare call labelled `WF-1-s11-cache-sweep`. Never raise a timeout or change `-n`.
- **Rule (i):** no edit to `study_yaml.py`, the module, the YAML files or the registry after the window that closes on them.
- Write ≈, never a tilde. No AED number anywhere.

**Scope.** Covered: `femem run` → `Study.tune` on one fixture: the F-small exact-rectangle PEC hole (80 798 cells, κ off, `strip_fraction = GATED_WIDTH_FRACTION`), degree 1, 64 MHz, `-n 2`, k = 1, lossless, through `run_study_file` in-process and `cli.main` on the refusal. **Not claimed:**
- `tune` through the installed console script;
- `tune` at degree 2 (item 445 is the API at degree 2; YAML at degree 2 is not run);
- `tune` on copper, on the solid σ = 800 four-port, or on F-human (refused);
- k = 2 as a gated mode (a control only);
- a tuned `f₀`, Q, match or any absolute `S₁₁` / `Z_in` / C figure;
- U25(b), the tuned in-model 4×4 (with the weekly);
- a stable schema;
- any AED reading.

`ANS-7` stays 🟡.

**What a negative result means.** Report and stop. Park on `attempt/WF-1-step11-<UTC>`, mark item 442 BLOCKED with its unblock condition in the same commit (rule (d)), and journal it in attempts.md.
- (K) or (X) red: a parse defect. Fix it in-slot inside `_parse_tune` / `_parse_mapping` only, then re-run W0.
- (R0) red: the PEC route moved since item 429. File a known-issues entry and stop; never re-record.
- (I=) red: the summary is not the Study's `tune()`, a `_tune_summary` / JSON defect. Fix it in-slot, then re-run W1.
- (I-Py) red with (R0) and (I=) green: two equal specs give different S across builds beyond 1e-9. Print both S, file a known-issues entry, never re-band.
- (a) red with (R0), (I=) and (I-Py) green: impossible without (R0) moving. Report it.
- (c1) inside the band: the summary ignores `k`. Stop.
- RG1 / RG2 red: the edit changed more than the new section. Revert and stop.

No `PREDICTED RED:` line: a red anywhere stops the slot.

**Status it can move.**
- **§7 `WF-1` row:** stays ✅ (rule (j): its Done-when is step 1's). On green, the row gains a step-11 record clause, exact text from the review's return, of the form: "**Step-11 record (§9 item 442, U49, provisional):** a top-level `tune: {k}` section drives `Study.tune` through `run_study_file` on the F-small PEC hole Study at degree 1, `-n 2`: `C_mode` <x> F, <y> from `TH-17`'s record (`REDUCTION_BAND`), equal to the same Study's `tune()` by `float.hex` and to an independent Python Study's within `LEG_D0_REPRODUCTION_BAND` (<z>); k = 2 misses by <p>; bad `k` / unknown key / wrong geometry raise `StudyConfigError` naming `tune.k` / `tune.<key>` / `tune`, rc 2, nothing written."
- **§10 chain U:** U49 is marked DONE (provisional until the next weekly ratifies).
- **`WF-14`:** the envelope gains six nodes; `docs/guide/api/load_study.md` / `run_study_file.md` / `StudyConfigError.md` and `docs/guide/validation.md` are regenerated.
- No §2 sentence moves; no known-issues entry is filed on green.

**For the review.**
1. Top-level `tune:` (as commissioned, carried inside the returned `outputs` dict) or `outputs.tune` (no `_TOP_KEYS` edit; it matches `outputs.sar` / `b1_plus`, the other post-solve reads)? Both keep the 4-tuple. The brief drafts the first.
2. Is dropping a wrong-fixture control acceptable, given that item 429's (c2) measured it on the same `tune` lines?
3. The commission's "equals item 429's Python route by `float.hex`" is redrafted as (I=) bitwise on one Study, plus (I-Py) at `LEG_D0_REPRODUCTION_BAND` across builds, per rule (j). Accept, or rule otherwise.

**Spawn line:** "Use the implementer agent to execute `WF-1` §9 item 442 per `docs/planning/chunks/WF-1.md` (read the latest Review ruling paragraph first). The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run W0, W1, RG1, RG2, Wgen, Wval, Wapi, Wdoc in order, foreground, `timeout: 660000`, pytest with `-s`, container `timeout -k 30`; before W1 and before RG2 run `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` and require 0. Edits in `src/` only in `study_yaml.py`; never change `_load`'s 4-tuple; never pass κ; the YAML says `degree: 1`; never import a `test_*` name; never compare across builds by `float.hex`. You are the executor — do not spawn agents, never background a harness window. On a red, park on `attempt/WF-1-step11-*`, mark item 442 BLOCKED in the same commit, and report."

### Review ruling 2026-10-14 03:00 (item 442)

The brief's two questions are ruled:
- **Schema.** `tune:` is a **top-level** YAML section, as U46's `thermal:` is. It is carried inside the parsed `outputs` dict so `_load`'s four-tuple (more than 12 callers) does not grow.
- **Wrong-fixture control.** Dropping it is accepted. The k = 2 control (item 429's measured +48.95 %) and the refusals are enough.
- **Bitwise scope.** The bitwise check is against the same Study's `tune()` only. Across two builds the bar is `LEG_D0_REPRODUCTION_BAND` (rule (j) of the brief).
- **Re-gates.** RG1 / RG2 are owed in the slot (standing rule (c)), as are the registry pages.

## §7 row moved verbatim 2026-10-14 (rotate_plan_archive.py chunks)

 `femem run`: the MRI study CLI and YAML config | ✅ 2026-10-06 (§9 item 256, U5): `femem run` builds Studies equal to the Python-built ones by `float.hex` and reproduces `MAT4_DEG1_RECORD` and `P1` at their bands. Steps 2–6 (the installed console script, `thermal:`, `conductor:`, `outputs.sar` / `outputs.b1_plus`, `outputs.probes`) are records at their imported bands — U32–U34 ratified by the 2026-10-10 weekly. Step 7a (§9 item 390a, U43 part 1, provisional, 2026-10-12): the `f_human_birdcage` (hole route) + `two_port_quadrature` schema parses bit for bit into item 382's tuned copper F-human Study and 17 planted mutations plus three files name their field, with `run_study_file`'s F-human branch exercised on a fake (no mesh; the real solve is item 390b). Step-7 record (§9 item 390b, U43 part 2, provisional, 2026-10-12): `run_study_file` on that file at `-n 8` solves once and its summary reproduces items 382 / 382b's node S (4.29e-12), `Y_in`, the 18-row ccw 10 g lattice (max 2.86e-10), cw at M(c\*) and mean \|B₁⁺\|_ccw (1.10e-07) inside `REUSE_REPRODUCTION_RTOL`, with the mirror gates and both controls green (`20261012T125631Z_WF-1-s7b-W1.log`, 243 s, heavy). Step-8 record (§9 item 394, U44, provisional, 2026-10-12): `outputs.sar` / `outputs.b1_plus` with a `leontovich` conductor on the F-small hole Study return the Python copper Study's field-route 10 g SAR and `|B₁⁺|` at the origin within `REUSE_REPRODUCTION_RTOL` (max 4.902e-13 / 5.397e-14 over four drives), the YAML field route equals its `solve()` to 6.442e-14 (`LEG_D0_REPRODUCTION_BAND`), the kept field's copper identity closes at 8.501e-13 (`DISCRETE_IDENTITY_RTOL`); a solid-route `leontovich` file refuses at `conductor` as `Study.sar` / `b1_plus_at` / `solve()` refuse (`as_hole`), with no mesh (`20261012T140300Z_WF-1-s8-W1.log`, 151 s, heavy by ceiling, `-n 2`). Step 10a (§9 item 431a, U46 part 1, provisional, 2026-10-14): the tuned-F-human thermal schema, parse-only — `tests/package/wf1/f_human_tuned_copper_thermal.yaml` parses bit for bit into P17's `RobinCooling` (G from the series JSON), `TransientSpec` (Δt 31.35 s, N 40, CN + 2 half-steps), c\*_θ and `W_B2`; 10 mutations name their field; `run_study_file`'s F-human branch calls `Study.temperature` (new `_temperature_port_study` hook) once on counting fakes, no mesh (`20261014T005102Z_WF-1-s10a-W0.log`, 4 s, standard). Step-10 record (§9 item 431b, U46 part 2, provisional, 2026-10-14): `thermal:` on the tuned copper F-human YAML: `femem run` at `-n 8` reproduces P17's 10 g ball-mean θ(t\*) (4.926110850e-02 K per accepted W, rel 8.0e-14) and `MAT-20` step 3's six-minute maximum (4.277005460e-02 K at 1 254 s, rel 1.4e-13) at `REUSE_REPRODUCTION_RTOL` and writes the 41-step series; one solve; controls (mirror ball 13.21 %, short window 1.35 %, doubled source 1.000) green (`20261014T020254Z_WF-1-s10b-W1.log`, 147 s, heavy). U46 is done. **Step-11 record (§9 item 442, U49, provisional):** a top-level `tune: {k}` section drives `Study.tune` through `run_study_file` on the F-small PEC hole Study at degree 1, `-n 2`: `C_mode` 6.626796421e-11 F, −2.977e-4 from `TH-17`'s record (`REDUCTION_BAND`), equal to the same Study's `tune()` by `float.hex` and to an independent Python Study's within `LEG_D0_REPRODUCTION_BAND` (rel 0.000e+00); k = 2 misses by +48.95 %; bad `k` / unknown key / wrong geometry raise `StudyConfigError` naming `tune.k` / `tune.<key>` / `tune`, rc 2, nothing written (`20261014T110326Z_WF-1-s11-W1.log`, 62 s, heavy by ceiling). *History: `docs/planning/chunks/WF-1.md`.* 

### §9 item 484 — 2026-10-15 (§10 chain U56, daily provisional)

*Drafted by a general-purpose agent for the 2026-10-15 18:00 daily review; the review reviews it.* Licence: §9 rule (6), a provisional next step on live chain U (Phase 9, 9a), after U53 (§9 item 472, `WF-16` step 12, `24088335`) and U54 (§9 item 473, `WF-17` step 6). **Target: reproduction plus an identity. No band, tolerance, record, Done-when or validation target is set or moved; every band is imported by name.** Anchors: `WF-12` step 3's H12c `C_mode` 1.181185799070e-10 F at `REDUCTION_BAND` (the comparand item 472 already reproduced through the API), and summary = API at `EXACT_IDENTITY_RTOL`. The sheets-and-κ ruling (2026-10-02) binds. The capacitor-ESR ruling is not engaged (the read is lossless, `study.py:1027–1029`). Item 442 (this file) and item 472 (`WF-16.md`) are the templates.

**Owner.** `WF-1` owns `femem run`'s `tune:` section: item 442 is `WF-1` step 11 and this file holds its brief (`:1108`); `WF-17` owns the guide pages (step 6, item 473, is a page) and `EX-115` is the example of item 442's route. This is `WF-1` step 12.

**Finding at drafting (code reading): the YAML route cannot reach the ring read today. A `src/` change is needed, `src/fem_em_solver/study_yaml.py` only; `study.py` is untouched.** Three walls, in parse order:
- **Excitation.** `_EXCITATION_FOR["f_human_birdcage"] = "two_port_quadrature"` (`study_yaml.py:119–124`); `_parse_excitation` refuses any other kind at `excitation.kind` (`:361–367`). Item 472's Study needs `LumpedPortDrive()`.
- **Tune.** `_parse_tune` refuses every non-four-port geometry at `tune` (`:724–727`).
- **Run.** `run_study_file`'s F-human branch (`:1061–1094`) reads `drive.c_in_f`, `result.s_node`, `result.y_in`, which are the `TwoPortQuadratureDrive` result's attributes. A `LumpedPortDrive` F-human Study returns a `PortStudyResult` (`study.py:574`, through `_study_solve_ports`, `:609`, (S1) of item 472 at `:619–622`). That branch would raise `AttributeError` after the solve.
- **Nothing else is missing.** `Study.tune` already dispatches the F-human `LumpedPortDrive` Study to `_study_tune_ring` (`study.py:1033–1034`). `_tune_summary` (`study_yaml.py:996–1005`) is geometry-blind: `z_k` is serialised by length, so the ring read's 9 entries serialise as the four-port's 4 do. `PortStudyResult.write_touchstone` takes any port count (`ports/touchstone.py:15–18`, `.s{n}p`).

**Pinned refusals the change must keep (found at drafting).**
- `tests/package/test_wf1_yaml_f_human.py:223` (item 390a): `f_human_tuned_copper.yaml` with only `excitation.kind` set to `lumped_port_drive` (so `c_in_f` / `esr_ohm` remain) must still name `excitation.kind`. Admitting the kind and falling through to the four-port `_unknown` check would name `excitation.c_in_f` and turn that gate red. **The drafter's design keeps it green with no test edit** (E1b below). The alternative is to re-point that row's expected field to `excitation.c_in_f`, which edits a pinned assertion; the drafter does not propose it, and the review rules if E1b is declined.
- `tests/package/test_wf1_yaml_tune.py:132–136` (item 442): `tune` on `f_human_tuned_copper.yaml` (quadrature) names `tune`. It stays so (E2 admits `f_human_birdcage` only with a `LumpedPortDrive`).
- `malformed_quadrature_on_four_port.yaml` → `excitation.kind` (`test_wf1_yaml_f_human.py:70`): unchanged (E1 widens only the F-human row).

**The change — `study_yaml.py` only, ≈ 35 lines.** Existing-line edits are E1, E2, E3 and two docstrings, so standing rule (c) applies: disclose them, and re-run the pre-existing gates (RG1, RG2, and W0's re-gate modules).
- **(E1)** `_parse_excitation`'s mismatch test (`:361`) admits one more pair: `kind == "lumped_port_drive" and geometry_kind == "f_human_birdcage"`. Every other pair refuses as now, message unchanged.
- **(E1b)** Additive, in the same function, before the `LumpedPortDrive` `_unknown` (`:376`): on an `f_human_birdcage`, if the section carries any of `c_in_f`, `esr_ohm`, `feeds` (the quadrature keys), raise `StudyConfigError(source, "excitation.kind", "lumped_port_drive takes no c_in_f / esr_ohm / feeds (two_port_quadrature's keys); …")`. A user who switched the kind but kept quadrature keys is told at the kind; item 390a's `:223` row stays green unedited.
- **(E2)** `_parse_tune(sec, geometry_kind, source, excitation=None)`: admitted on a `four_port_birdcage` (as now) or on an `f_human_birdcage` whose `excitation` is a `LumpedPortDrive`; anything else refuses at `tune` (message names both admitted forms). `k` handling unchanged.
- **(E3)** `_parse_mapping`'s call (`:792`) passes `excitation`.
- **(E4)** Additive refusals in `_parse_mapping`, each placed at its section's place in the fixed validation order: on an `f_human_birdcage` + `LumpedPortDrive`, a `thermal:` section raises at `thermal` (after the `cooling` line), and `outputs.sar` / `outputs.b1_plus` raise at `outputs.sar` / `outputs.b1_plus` (after the `outputs` line). Reason: those routes are the `two_port_quadrature` per-sense routes (`Study.sar` / `b1_plus_at` / `temperature` refuse an F-human `LumpedPortDrive` Study, `study.py:526–530`); without E4 they would fail after the solve.
- **(E5)** Additive branch in `run_study_file`, **between** the four-port branch and the quadrature F-human branch: `elif isinstance(study.geometry, FHumanBirdcageGeometry) and isinstance(study.excitation, LumpedPortDrive):` writing `cells`, the `.s32p` Touchstone, `s_matrix`, `z_matrix` and `touchstone` exactly as the four-port branch's `:1047–1054`, then `if "tune" in outputs: summary["tune"] = _tune_summary(study, outputs["tune"], comm)`. Neither neighbouring branch gains or loses a line.
- **Docstrings:** the module's `tune` bullet (`:66–70`) and the `excitation` lines gain the F-human `lumped_port_drive` form; `run_study_file`'s docstring (`:1014–1016`) says the quadrature route writes no Touchstone and the `lumped_port_drive` route writes `<stem>.s32p`. Wgen regenerates the API page.
- **Never:** change `_load`'s 4-tuple; add a `Study.__init__` keyword; touch `study.py`; pass κ; admit an F-human `lumped_port_drive` file without `as_hole: true` and a `leontovich` conductor (both refusals exist: `:279–283` and `:773–778`, unchanged).

**`WF-1` step 12: `femem run` (`run_study_file`) with a top-level `tune: {k: 1}` on the copper F-human hole `lumped_port_drive` YAML at degree 1, `-n 8`, solves once; `summary["tune"]` equals the same Study's `Study.tune(k=1)` and reproduces `WF-12` step 3's H12c `C_mode` within `REDUCTION_BAND`; k = 2 misses.**
- **Executor:** implementer.
- **Tier:** heavy, `-n 8` (the comparand's width, item 472's W1), complex. W0 / RG1 / RG2 `-n 2`.
- **Class:** chain step (solve), median ≈ 10.5 measured.
- **Cost: ≈ 45 predicted min** = 15 fixed + ≈ 12 code (≈ 35 `src` lines, module ≈ 200 lines, one YAML, registry lines, one guide bullet) + W0 ≈ 0.5 + RG1 ≈ 0.3 + RG2 ≈ 1.5 + W1 ≈ 4–5 + Wgen / Wval / Wapi / Wdoc ≈ 1 + ≈ 7 row / journal / spare. **If the `src/` change and W0 / RG1 / RG2 are not green by minute 30, or W1 has not started by minute 35, park** (rule (d)).
- **Depends on:** nothing open. Reads `logs/wf12/h12c-cu-s32-a2.npz` through `NPZ_PATH` (`tests/validation/test_wf10_f_human_copper_esr_root.py:50`); absent ⇒ stop, BLOCKED, never re-solve. Never beside another heavy window. On a rebase over another item touching `study_yaml.py`, `validation_envelope.toml` or the generated pages, keep both sides and re-run W0, RG1, Wgen / Wval / Wapi / Wdoc.
- **Family:** `WF-1` step 12, attempt 1 of 4.

**Tracked input: new `tests/package/wf1/f_human_hole_copper_tune.yaml`.** First line a `# WF-1 step 12 (§9 item 484): …` comment, then `femem: 1`; `study: {frequency_hz: 64000000.0, degree: 1}`; `geometry: {kind: f_human_birdcage, as_hole: true}`; `material: {sigma: 0.5, epsilon_r: 78.0}`; `excitation: {kind: lumped_port_drive}`; `conductor: {kind: leontovich, sigma_s_per_m: 58000000.0}`; `tune: {k: 1}`. No `outputs:` section, so the summary goes to the `output_dir` argument. (K) proves it equals item 472's Study field for field.

**Test module: new `tests/package/test_wf1_yaml_tune_f_human.py`.** First docstring paragraph backticks only `WF-1`, `WF-16` and `WF-12` (the (F3) rule). Opt-in solve: `FEM_EM_WF1_S12=1`, item 472's `_knob()` pattern (unset ⇒ the solve nodes skip, any other value raises); width 8 asserted when set; `@complex_only` on the solve nodes.
- **Imports (names and constants only, never a `test_*` name):** `_study`, `WF12_STEP3_C_MODE_CU_F`, `WF12_STEP3_K2_SEPARATION`, `RECORD_RANK_WIDTH` from `tests/package/test_wf16_study_tune_f_human_copper.py` (`:104`, `:71`, `:74`, `:66`); `_spec_diff`, `_fresh_dir`, `_mapping`, `LOGS`, `WF1`, `_rc_all` from `tests/package/test_wf1_study_yaml.py` (as item 442, `test_wf1_yaml_tune.py:39–46`); `NPZ_PATH` (`test_wf10_f_human_copper_esr_root.py:50`); `HOLE_CELL_RECORD` (`test_wf12_f_human_copper_matrix.py`); `_residual` (`test_port_circuit_layer_field.py:142`); `REUSE_REPRODUCTION_RTOL`, `REDUCTION_BAND`, `EXACT_IDENTITY_RTOL` from item 472's sources (`WF-16.md` item 472's import table); `load_study`, `run_study_file`, `StudyConfigError`, `cli` from `fem_em_solver`; `_parse_mapping`, `_study_from_mapping`, `_tune_summary` from `fem_em_solver.study_yaml`.
- **Printed beside, transcribed with log lines:** item 472's measured `C_mode` 1.181185799070e-10 F and gap −3.544942e-13 (`20261015T213534Z_WF-16-s12-W1.log:7375`), k = 2 gap −0.3327259 (`:7380`).
- **Nodes:** `test_wf1_s12_parse_and_refusals` (ungated, no mesh); `test_wf1_s12_r0_fixture_tie`; `test_wf1_s12_summary_equals_same_study_tune`; `test_wf1_s12_a_c_mode_equals_wf12_h12c`; `test_wf1_s12_c1_mode_index_control`.
- **Module fixture (solve):** `_fresh_dir(LOGS / "f_human_hole_copper_tune", comm)`; `run = run_study_file(FILE, output_dir=d, comm=comm)`; read `summary.json` on every rank (after the `Barrier`); `res0 = run.study.result`; `api1 = run.study.tune(k=1, comm=comm)`; assert `run.study.result is res0` (one solve); print `res0.mesh_time`, `res0.sweep_time`, `run_study_file` wall time and `report_peak_rss(comm, label="[WF-1 s12] summed peak ru_maxrss")`; npz read on rank 0, broadcast.
- **Registry (lines added only):** the parse node under `[names.StudyConfigError]` and `[names.load_study]` beside item 442's (`tests/package/validation_envelope.toml:262`, `:283`); the four solve nodes under `[names.run_study_file]` after item 442's block (`:345–349`), each with `# WF-1 step 12 (§9 item 484)`.
- **Guide:** `docs/guide/tuning_from_yaml.md:106–109`'s "F-small only through YAML" bullet is now false. Replace it with one bullet: on a `four_port_birdcage` the section reads the four-port route; on an `f_human_birdcage` with `excitation.kind: lumped_port_drive` it reads the 32-port ring route (degree 1, the copper hole fixture), citing the new registered nodes only. No number is added; the embedded YAML block is untouched (item 473's (Y=) stays byte equal).

**Anchors (asserted; bands imported).**
- **(K) Parse, no mesh.** Every geometry's `build` monkeypatched to raise (as `test_wf1_yaml_tune.py:80–85`). `_spec_diff(load_study(FILE), _study()) == []` (item 472's Python Study: `FHumanBirdcageGeometry(as_hole=True)`, saline, `LumpedPortDrive()`, `LeontovichWall(COPPER_SIGMA)`, 64 MHz, degree 1); parsed outputs `== {"tune": {"k": 1}}`, `type(k) is int`; `tune: {}` → `{"tune": {}}`.
- **(X) Refusals, no mesh, each naming its field.** On mapping variants of `FILE`: `tune: {k: 3}` → `tune.k`; `excitation.c_in_f` added → `excitation.kind` (E1b); `thermal:` added (item 390a's `_add_thermal` values) → `thermal`; `outputs.sar` (one literal centre) → `outputs.sar`; `outputs.b1_plus` → `outputs.b1_plus`; `as_hole: false` → `geometry.as_hole`; `conductor` removed → `conductor`. On `f_human_tuned_copper.yaml` + `tune: {k: 1}` → `tune` (unchanged). `cli.main(["check", FILE])` rc (0, 0) via `_rc_all`, rank 0 prints `family: f_human_birdcage`. Zero `build` calls.
- **(R0) Fixture tie.** `summary["cells"] == HOLE_CELL_RECORD` (369 387); `_residual(S_run[p][:, p], s32_npz) ≤ REUSE_REPRODUCTION_RTOL` with `p` the name map to the npz `order` (item 472's construction, `test_wf16_study_tune_f_human_copper.py:242–259`); `len(summary["s_matrix"]) == 32`; `summary["touchstone"] == FILE.stem + ".s32p"` and the file exists (rank 0, broadcast). Predicted ≤ 1e-10 (item 472 measured 1.632e-12 on the same build route, `…213534Z_WF-16-s12-W1.log:7346`).
- **(I=) Summary = API, one solve.** `|summary["tune"]["c_mode_farad"] / api1.c_mode_farad − 1| ≤ EXACT_IDENTITY_RTOL`; the same for `off_diagonal_residual`; `len(z_k) == 9` and every real / imaginary part within `EXACT_IDENTITY_RTOL` of `api1.z_k`'s; JSON read-back equal to the in-memory summary; `k == 1`, `frequency_hz == 64e6`. Printed beside: `float.hex` equality (predicted equal: same `self.result`, same arithmetic, item 442's (I=) reading).
- **(a) The comparand.** `|summary C / WF12_STEP3_C_MODE_CU_F − 1| ≤ REDUCTION_BAND`. Predicted ≤ 1e-9 (item 472: −3.545e-13 through the API on the same fixture and width). Printed beside: the gap to item 472's measured C.

**Negative control (rule (e)).**
- **(c1) The mode index — *asserted*.** `t2 = _tune_summary(run.study, {"k": 2}, comm)`; assert `run.study.result is res0` and `|t2["c_mode_farad"] / WF12_STEP3_C_MODE_CU_F − 1| > REDUCTION_BAND`. Backing: the same comparison, fixture and width measured −0.3327259 (332.7× the band) through the API (`…213534Z_WF-16-s12-W1.log:7380`) and 0.332726 in step 3 (`20261006T170222Z_WF-12-step3-a2-W1.log:6397`). Predicted −33.27 %. Ceiling: the assertion is one-sided (`> band`), so no ceiling binds.
- **Printed only:** `WF12_STEP3_K2_SEPARATION` beside (c1)'s gap.

**Windows.** Foreground, in the order W0, RG1, RG2, `pgrep`, W1, Wgen, Wval, Wapi, Wdoc. **The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped.** Bash `timeout: 660000`; pytest with `-s`; node IDs, never `-k`; `tests/environment/test_complex_mode.py` first. Run W0 only after every edit (rule (i)).
- **W0** (no mesh, `-n 2`; the new parse node plus the rule-(c) parse re-gates: item 442's parse node, `test_wf1_check_family.py` (the new good file joins its glob), item 390a's and 431a's no-mesh modules (E1 / E1b / E4 sit in their path; their Run-0 fake nodes cover the untouched quadrature branch), and item 473's page module (the guide bullet)):
  `scripts/testing/run_and_log.sh WF-1-s12-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_tune_f_human.py::test_wf1_s12_parse_and_refusals tests/package/test_wf1_yaml_tune.py::test_wf1_s11_tune_parse_and_refusals tests/package/test_wf1_check_family.py tests/package/test_wf1_yaml_f_human.py tests/package/test_wf1_yaml_f_human_thermal.py tests/package/test_wf17_tuning_yaml_page.py -v -s --tb=short'"`
- **RG1** (rule (c): every `[names.load_study]` node): item 442's RG1 command verbatim (this file, item 442 "Windows", `20261014T110438Z_WF-1-s11-RG1.log:12`), label `WF-1-s12-RG1`; 13 passed, 3 s there.
- **RG2** (rule (c): `_parse_tune` / `_tune_summary` on the four-port route): item 442's W1 verbatim (`20261014T110326Z_WF-1-s11-W1.log:12`), label `WF-1-s12-RG2`; 10 passed, 62 s there.
- **Not re-run (disclosed; the review rules):** item 390b's `-n 8` quadrature solve (`20261012T125631Z_WF-1-s7b-W1.log`, 243 s): no line of the quadrature branch moves and W0's Run-0 fakes exercise it.
- **Before W1:** `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reads 0.
- **W1** (heavy, `-n 8`; predicted 200–260 s, ≤ 345 s at the box's 1.32×):
  `scripts/testing/run_and_log.sh WF-1-s12-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_WF1_S12=1 timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_tune_f_human.py -s -v --tb=short'"`
  - **Price:** item 472's W1 on this Study, mesh and width: 203.4 s, mesh 106.5 s, sweep 93.4 s, 9.765 GiB summed (`20261015T213534Z_WF-16-s12-W1.log:7343–7344`). Added here: one 32-port Touchstone write and the JSON (< 2 s); `tune` is NumPy (< 1 s). One solve only: no independent Python Study (item 442's (I-Py) is dropped — a second build is ≈ 200 s more and item 472 already measured the API on this fixture).
  - **JIT:** the forms compiled in item 472's window; allow ≈ 30 s if they hash differently.
- **Wgen / Wval / Wapi / Wdoc:** item 394's four commands verbatim (this file, item 394 "Windows"), labels `WF-1-s12-Wgen` / `-Wval` / `-Wapi` / `-Wdoc`. Wval's node count is the count before the item plus 5 (229 + 5 = 234 if nothing else landed; print it). Wdoc reads `0 dead`.
- **`git diff --stat main -- src/`** shows `src/fem_em_solver/study_yaml.py` only. Say so in the commit (rule (c), E1–E5 named).

**Traps already paid for.**
- **`degree: 1` in the YAML.** The Study default is 2; a degree-2 copper F-human solve at `-n 8` was killed at 590 s (item 472's trap list). (K) catches a missing `degree`.
- **YAML floats:** write `64000000.0` and `58000000.0`, never `64.0e6` (`_YAML_FLOAT_HINT`, `study_yaml.py:115`).
- **`femem` prints on rank 0 only:** read `check`'s stdout on rank 0 and reduce rc through `_rc_all`; never assert a print on every rank.
- **`-k a or b` splits into stray argv** inside the quoted container command: node IDs only.
- **Registry nodes added ⇒ API pages and `docs/guide/validation.md` regenerated and committed**, Wgen / Wval / Wapi / Wdoc in the slot (ruling Q3 of 2026-10-12; item 472 did the same).
- **Rank-local values:** cells and S come through the Study's reduced paths; the Touchstone exists on rank 0 only — broadcast the existence bool.
- **The npz is a fixture input**, not the comparand; (a)'s comparand is step 3's `:6366`.
- **Never `TwoPortQuadratureDrive` for the read** (its `tune` refuses; item 442's (X) row pins it at parse).
- **Never `float.hex` across builds** (item 442 ruling, `:1290`); here the identity is the same Study, asserted at `EXACT_IDENTITY_RTOL`, hex printed.
- **Rule (i):** no edit to `study_yaml.py`, the module, the YAML or the registry after the window that closes on them; a `src` edit after W0 means W0, RG1, RG2 again.
- **Status 124 / a kill:** `pgrep -c python3` first; at 0 only, the rubric orphan-stem sweep, label `WF-1-s12-cache-sweep`. Never raise a timeout or change `-n`.
- Write ≈, never a tilde. No AED number.

**Scope.** Covered: k = 1 `C_mode` from a study file's `tune:` section on the copper (Leontovich, σ = 5.8e7 S/m) F-human hole, κ off, degree 1, 64 MHz, `-n 8`, one solve, equal to the same Study's `tune()` and to `WF-12` step 3's H12c reading; `run_study_file` is `femem run`'s call (`cli.py`, `run` → `run_study_file`). **Not claimed:** the installed console script or `--ranks` on this file (not run); degree 2; an h-statement; a tuned `f₀`, Q, `R_in`, match, `C_in` or resonance; ESR; the `.s32p` or `s_matrix` as validated quantities beyond (R0)'s tie; SAR / `|B₁⁺|` / thermal on the `lumped_port_drive` F-human file (refused, E4); the PEC or solid σ = 800 F-human from YAML (still refused); `--estimate` (item 480's rc 2 stands); AED.

**What a negative result means.** Report and stop. Park on `attempt/WF-1-step12-<UTC>`, mark item 484 BLOCKED in the same commit (rule (d)), journal it.
- **(K) / (X) red:** a parse defect. Fix in-slot inside `study_yaml.py` only and re-run W0 (rule-(i) clock restarts).
- **W0's re-gate modules, RG1 or RG2 red:** the edit changed more than the new route. Revert and stop; the review rules. Never edit a pinned row to fit.
- **(R0) red:** the YAML-built Study is not item 472's on the same build route; print the worst column and the cell count, file nothing new until the review rules.
- **(I=) red:** `_tune_summary` / the JSON loses the ring read. A defect in this item's `src`; no band move.
- **(a) red with (R0) and (I=) green:** custody of the comparand; print item 472's C beside and report.
- **(c1) inside the band:** a blind comparator. Stop.
- **W1 status 124:** not a result; park with the measured time.

No `PREDICTED RED:` line. A red anywhere stops the slot.

**Status it can move.**
- **`WF-1`** row stays ✅ (Done-when is step 1's). On green it gains a **step-12 record (§9 item 484, U56, provisional)**: "a top-level `tune: {k: 1}` on the copper F-human hole `lumped_port_drive` file drives the ring read through `run_study_file` at degree 1, `-n 8`, one solve: `C_mode` <x> F, <y> from `WF-12` step 3's H12c reading (`REDUCTION_BAND`), equal to the same Study's `tune()` at <r> (`EXACT_IDENTITY_RTOL`); 32×32 tied to the H12c npz at <s>; k = 2 misses by <p>; `src`: `study_yaml.py` only (rule (c), re-gates green)".
- **§10 U56** marked DONE (provisional until the weekly ratifies).
- **`WF-14`:** five nodes; API pages and `validation.md` regenerated.
- **`WF-17`:** step 6's page loses a statement this item falsifies (one bullet, re-gated in W0).
- **For the weekly:** an `EX-115` F-human companion (the 32-port `tune:` from a file) becomes openable once U56 is ratified; not opened here.

**Spawn line:** "Use the implementer agent to execute `WF-1` §9 item 484 (U56) per `docs/planning/chunks/WF-1.md` (read the brief and any review ruling beneath it first; the rulings hold). The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run W0, RG1, RG2, W1, Wgen, Wval, Wapi, Wdoc in order, foreground, `timeout: 660000`, pytest with `-s`, container `timeout -k 30`; before W1 run `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` and require 0. `src` edits: `src/fem_em_solver/study_yaml.py` only, exactly E1–E5 plus the named docstrings, disclosed in the commit; never edit `study.py`, `_load`'s 4-tuple, or any existing test assertion. The YAML says `degree: 1`; never κ; never import a `test_*` name. You are the executor — do not spawn agents, never background a harness window. On a red, park on `attempt/WF-1-step12-*`, mark item 484 BLOCKED in the same commit, and report."

### Review ruling 2026-10-15 18:00 (item 484, U56 — provisional)

**Queued, with three rulings.**
- **The `src/` change is licensed under ruling (c):** `study_yaml.py` only, E1–E5 as briefed (≈ 35 lines), disclosed in the commit. Any other `src/` file in the diff parks the item.
- **E1b, not the re-point.** Item 390a's pinned row (`tests/package/test_wf1_yaml_f_human.py:223`) stays unedited and green. An existing assertion is never re-pointed to admit a new feature.
- **Item 390b's `-n 8` F-human YAML solve is re-run in this slot**, by its node ID (never a `-k a or b` expression). The change touches the F-human run branch, and ruling (c) requires a touched path's pre-existing gate to be re-run green in the same slot. The brief's "not re-run, disclosed" is overruled. Two heavy windows never run at once: W1 first, then 390b. If both are not green by minute 40, park on `attempt/WF-1-s12-*`, BLOCKED in the same commit.

The anchors set no target. They are the H12c record at `REDUCTION_BAND` (U53's comparison, item 472 `20261015T213534Z_WF-16-s12-W1.log:7375`), the API = CLI identity at `EXACT_IDENTITY_RTOL`, and the k = 2 miss asserted on this fixture (`:7380`). Chain U's provisional count is now 1 of 3. The weekly ratifies or reverts.

**Outcome 2026-10-16 (22:30 slot): green, DONE (provisional).** `src`: `study_yaml.py` only, E1–E5 plus the module / `_parse_tune` / `run_study_file` docstrings. W0 19 passed (`20261016T034124Z_WF-1-s12-W0.log`; the first W0, `20261016T034101Z_WF-1-s12-W0.log`, was red on item 473's pinned `len(nodes) == 3` in `test_wf17_tuning_yaml_page.py` because the guide bullet cited a new node; the bullet now points to the `run_study_file` section of the validation page instead of citing node IDs — no assertion edited); RG1 13 passed (`…034140Z_WF-1-s12-RG1.log`); RG2 10 passed in 82.8 s (`…034148Z_WF-1-s12-RG2.log`). W1 `-n 8`: 9 passed in 283.0 s; cells 369 387, mesh 150.5 s, sweep 127.0 s, `run_study_file` 280.4 s, 9.564 GiB summed (host ≈ 1.4× item 472's 203 s); (K) diff `[]`; every (X) row names its field (the `thermal` row is refused by item 431a's F-human `thermal.transient` check before E4 is reached — same field); (R0) 1.616e-12, worst column 2.621e-12, `.s32p` written; (I=) 0.000e+00, `float.hex` equal (`0x1.03bed95bc98ddp-33`); (a) `C_mode` 1.181185799069577e-10 F, gap −3.580469e-13 (item 472 −3.544942e-13); (c1) k = 2 gap −0.3327259 (332.7× band) (`20261016T034323Z_WF-1-s12-W1.log`). Item 390b's `-n 8` solve re-run by module path (its original command): 10 passed in 329.3 s (`20261016T034818Z_WF-1-s12-390b.log`). Wgen / Wval (235 nodes: the brief's 229 + 5 predicted 234; one node landed between drafting and this slot) / Wapi 22 passed / Wdoc 0 dead (`20261016T0353…`–`035426Z_WF-1-s12-W*.log`). Printing note: the module's `ITEM472_C_MODE_F` carries item 472's printed 13-digit C, which equals the H12c record to the print, so the "vs item 472's C" line repeats the (a) gap.

### §9 item 496 — 2026-10-16 (§10 chain U57, daily — provisional; chain step (solve))

*Drafted by a general-purpose agent for the 2026-10-16 10:30 daily review; the review reviews it.* Licence: §9 rule (6) (PROJECT_PLAN.md ≈ :3850), a provisional next step of live §10 chain U toward subgoal 9a (a study-level Python API + YAML/CLI), after U52 (§9 item 463, `WF-16` step 6(b)) and U56 (§9 item 484, this file). Chain U's provisional count becomes 2 of 3. **Target: reproduction plus an identity. No band, tolerance, record, Done-when or validation target is set or moved; every band is imported by name.** Anchors: item 463's Route-A circuit = field identity (`PORT-22` (i)'s construction) at `REDUCTION_BAND`; the summary = the same Study's result by `float.hex`; the YAML Study = the Python `CapacitorSheets` Study on one mesh at `EXACT_IDENTITY_RTOL`. The sheets-and-κ ruling (2026-10-02) and the capacitor-ESR ruling (2026-10-06: lossless records are never re-read with r) bind: the solved file is r = 0. Item 484 (above) and item 463 (`WF-16.md`) are the templates. Provenance of the shape: item 390a's drafter sketched exactly this step as separable (this file, item 390a, "`capacitor_sheets` in YAML — separable"), and the 2026-10-11 18:00 ruling (4) kept it out of U43 only.

**Owner.** `WF-1` owns `femem run`'s schema (steps 1–12 are items 256 … 484). This is `WF-1` step 13.

**Finding at drafting (code reading on `2007d29d`, nothing run): YAML cannot express capacitor sheets today. A `src/` change is needed, `src/fem_em_solver/study_yaml.py` only; `study.py` is untouched.**
- **Excitation.** `_parse_excitation` admits only the three kinds in `_EXCITATION_FOR`'s values (`study_yaml.py:124–129`, `:365`); `capacitor_sheets` refuses at `excitation.kind`. `CapacitorSheets` is not even imported (`:98–108`). The only ESR key in the schema is `two_port_quadrature`'s `esr_ohm` (F-human non-feed gaps, `:428–432`), and E1b refuses it on any other route (`:383–391`).
- **Run.** Even if parsed, `run_study_file`'s four-port branch (`:1089–1103`) would fail after the solve: `CapacitorStudyResult.write_touchstone` raises `NotImplementedError` (`study.py:743–744`) and its `z_matrix` is `None` (`study.py:833`), which `_pairs` cannot serialise.
- **The API needs nothing.** `CapacitorSheets(c_f, esr_ohm=0.0, on=…, drive=LumpedPortDrive())` (`study.py:690–731`) dispatches through `_study_solve_ports`' hook (`:614–617`) to `_study_solve_capacitor_sheets` (`:747–847`), which refuses an F-human geometry, `on` outside `P1…P{leg_count}` and an `on` that leaves no port **before** any mesh (`:755–766`). `Study.tune` / `sar` / `b1_plus` / `mean_sar` / `solve_sweep` refuse a `CapacitorSheets` Study (`TypeError`, pinned at `tests/package/test_wf16_study_capacitor_sheets.py:126–134`).
- **YAML 1.1 trap found at drafting.** PyYAML's `SafeLoader` reads a bare `on` key as the boolean `True` (`on` / `off` / `yes` / `no` are YAML 1.1 booleans). A YAML key spelled like the API field (`on: [P2, P3, P4]`) would reach the parser as `{True: [...]}` and `_unknown` would name `excitation.True`. The schema key is therefore **`sheets`** (it maps to `CapacitorSheets.on`), and a `True` key in the section is refused with a message naming `excitation.on` and the fix. Disclosed; question 1 for the review.

**The change — `study_yaml.py` only, ≈ 50 lines.** Existing-line edits are C2, C3, C7's one keyword and two docstrings, so standing rule (c) applies: disclose them and re-run the pre-existing gates (W0's re-gate modules, RG1, RG2, RG3).
- **(C1)** The `from .study import (…)` block gains `CapacitorSheets` (one inserted line). `__all__` is unchanged (the `WF-13` installed-surface snapshot, `tests/package/public_surface.txt`, must not move).
- **(C2)** `_parse_excitation`'s `_kind` call (`:365`) admits `"capacitor_sheets"` beside `_EXCITATION_FOR`'s values. The unknown-kind message now lists four kinds; no test pins that text (drafter grepped `tests/` for "unknown kind" / "allowed: quasi": none).
- **(C3)** The mismatch test (`:366–368`) admits one more pair: `kind == "capacitor_sheets" and geometry_kind == "four_port_birdcage"`. On a `sphere_in_box` or `f_human_birdcage` it refuses at `excitation.kind` with the message unchanged (mirrors `study.py:755–759`, raised at parse).
- **(C4)** Additive, after the `two_port_quadrature` return (`:381–382`): `if kind == "capacitor_sheets": return _parse_capacitor_sheets(sec, geometry, source)`.
- **(C5)** New `_parse_capacitor_sheets(sec, geometry, source) -> CapacitorSheets`, in `_parse_two_port_quadrature`'s order (every value checked here so the dataclass's own `ValueError` never escapes bare):
  - a `True` key (the YAML 1.1 `on`) raises at `excitation.on`: "YAML 1.1 reads a bare `on` key as true; write `sheets: [...]`";
  - `_unknown(sec, ("kind", "c_f", "esr_ohm", "sheets"), …)`; `_missing(sec, ("c_f", "sheets"), …)`;
  - `c_f`: `_number`, finite and > 0, else at `excitation.c_f`;
  - `esr_ohm`: optional, `_number`, finite and ≥ 0, else at `excitation.esr_ohm`; **omitted is not filled in** (the dataclass default applies, as `study.degree`);
  - `sheets`: a list (a bare string refuses), each entry `_string`, non-empty, distinct, each in `{f"P{i + 1}" for i in range(geometry.leg_count)}` (at `excitation.sheets[i]`), and leaving at least one port to drive (at `excitation.sheets`);
  - no `drive` / `z0_ohm` / `port_impedance_ohm` keys: the kept ports take `LumpedPortDrive()`'s defaults (disclosed; unknown keys refuse).
- **(C6)** Additive refusals in `_parse_mapping`, each at its section's place in the fixed order: on a `CapacitorSheets` excitation, `outputs.sar` / `outputs.b1_plus` raise at that field (after the existing E4 loop, `:823–828`), and a `tune:` section raises at `tune` (before the existing `if "tune" in mapping:`, `:829`) — the API refuses all three with `TypeError` after a mesh would be built.
- **(C7)** `run_study_file`: a new **first** branch `if isinstance(study.excitation, CapacitorSheets):` writing `cells`, `s_matrix` (`_pairs` of the reduced matrix, kept-port order), `kept_ports` (`[p.port_id for p in result.port_definitions]`) and `capacitors = {"sheets": list(result.capacitor_port_ids), "c_f": float, "esr_ohm": float, "impedance_ohm": [re, im]}` — no Touchstone, no `z_matrix` (the reduced network has none, `study.py:743–744`). The existing `if isinstance(study.geometry, FourPortBirdcageGeometry):` (`:1089`) becomes `elif` — **the only edited line in `run_study_file`'s body**; no line of the four-port, F-human or sphere branches moves.
- **Docstrings:** the module's `excitation` bullet gains `capacitor_sheets` (`c_f`, `esr_ohm`, `sheets`; four-port only; the `on` trap); `run_study_file`'s docstring says a `capacitor_sheets` study writes no Touchstone. Wgen regenerates the API pages.
- **Never:** touch `study.py`; widen `_load`'s 4-tuple; add a name to `__all__`; pass κ; admit `capacitor_sheets` on an F-human geometry.

**`WF-1` step 13: a `capacitor_sheets` excitation makes `femem run` (`run_study_file`) solve item 463's Route-A Study — the degree-1 exact-rectangle PEC hole F-small with lossless capacitor sheets at the k = 1 `C_mode` on P2–P4 — from a tracked file; its summary equals the same Study's result by `float.hex`, equals the Python `CapacitorSheets` Study on the same mesh within `EXACT_IDENTITY_RTOL`, and reproduces item 463's circuit = field identity against the tune file's own untuned 4×4 within `REDUCTION_BAND`; the k = 2 reduction misses.**
- **Executor:** implementer.
- **Tier:** heavy by ceiling (container `timeout -k 30 300`, item 463's / 442's W1 ceiling), standard by prediction. `-n 2` (`ANS6_EXACT_STRIP_RECORD_RANK_WIDTH`, the comparand's width), complex.
- **Class:** chain step (solve), median ≈ 17.5 measured min.
- **Cost: ≈ 42 predicted min** = 15 fixed + ≈ 14 code (≈ 50 `src` lines, module ≈ 230 lines, two YAML files, registry lines) + W0 ≈ 0.3 + RG1 ≈ 0.1 + RG2 ≈ 1.5 + RG3 ≈ 1.7 + W1 ≈ 1.3 + Wgen / Wval / Wapi / Wdoc ≈ 0.5 + ≈ 8 row / journal / spare. **If the `src/` change and W0 / RG1 are not green by minute 25, or W1 has not started by minute 35, park** (rule (d)).
- **Family:** `WF-1` step 13, attempt 1 of 4.
- **Depends on:** nothing open. Item 463 is on `main` (its module supplies `ITEM429_C_MODE_F`); item 442's tune file is tracked. On a rebase over another item touching `study_yaml.py`, `validation_envelope.toml` or the generated pages, keep both sides and re-run W0, RG1, Wgen / Wval / Wapi / Wdoc.

**Tracked inputs (two new files).**
- **`tests/package/wf1/f_small_hole_pec_capacitors.yaml`:** `f_small_hole_pec_control.yaml` verbatim with its header replaced by `# WF-1 step 13 (§9 item 496): f_small_hole_pec_control.yaml with lossless capacitor sheets on P2-P4 at the k = 1 C_mode (20261015T201110Z_WF-16-s6b-W1.log:940).` and the excitation section
  ```yaml
  excitation:
    kind: capacitor_sheets
    c_f: 6.626796421e-11
    esr_ohm: 0.0
    sheets: [P2, P3, P4]
  ```
  `c_f` is the 10-digit print of item 463's in-run `C_mode` (`…s6b-W1.log:940`; item 429 / 442 printed the same digits), i.e. `ITEM429_C_MODE_F` in `tests/package/test_wf16_study_tune_capacitors.py:68` — a disclosed transcription whose custody (a′) asserts. It has a dot and a signed exponent, so PyYAML 1.1 reads it as a float. No `outputs:` / `tune:` section.
- **`tests/package/wf1/malformed_capacitor_sheets_all_ports.yaml`** (header line): the same file with `sheets: [P1, P2, P3, P4]` → `excitation.sheets` ("leaves no port to drive").

**Test module: new `tests/package/test_wf1_yaml_capacitors.py`.**
- **Docstring.** First paragraph backticks only `WF-1`, `WF-16`, `PORT-22`, `TH-17` (§7 row IDs; the (F3) rule), never `U57` or an item number. It says the build cache is a test economy licensed by (K)'s spec equality, and that the transcribed `c_f` is a record whose custody (a′) checks.
- `@complex_only` on every solve node. **Width asserted** in the fixture, never skipped: `MPI.COMM_WORLD.size == ANS6_EXACT_STRIP_RECORD_RANK_WIDTH` (2).
- **Imports (module-level names and constants only; never a `test_*` name):**

  | Name | Source |
  |---|---|
  | `LOGS`, `WF1`, `_fresh_dir`, `_mapping`, `_rc_all`, `_spec_diff` | `tests/package/test_wf1_study_yaml.py:85–170` |
  | `PEC_CONTROL_FILE`, `_py_copper` | `tests/package/test_wf1_yaml_copper.py:85`, `:98` |
  | `TUNE_FILE` | `tests/package/test_wf1_yaml_tune.py:60` |
  | `_geometry` | `tests/package/test_wf16_study_copper.py:80` |
  | `ITEM429_C_MODE_F` | `tests/package/test_wf16_study_tune_capacitors.py:68` |
  | `_study` (item 373's copper Study builder), `C_F`, `ESR_OHM` | `tests/package/test_wf16_study_capacitor_sheets.py:99`; `C_F` / `ESR_OHM` as that module imports them (`tests/validation/test_ans7_exact_strip.py`, `tests/validation/test_port24_series_rc_sheet.py`) |
  | `TH17_STEP3_C_MODE_RECORD_F` | `tests/package/test_wf16_study_tune.py:71` |
  | `ANS6_EXACT_STRIP_CLASS_MEANS_P1`, `ANS6_EXACT_STRIP_HOLE_CELL_RECORD`, `ANS6_EXACT_STRIP_RECORD_RANK_WIDTH`, `CLASSES`, `_class_means` | `tests/validation/test_ans6_exact_strip.py:103–144` |
  | `REDUCTION_BAND` | `tests/validation/test_port_lumped_rlc_termination.py:88` |
  | `LEG_D0_REPRODUCTION_BAND` | `tests/validation/test_port_birdcage_four_port.py:165` |
  | `EXACT_IDENTITY_RTOL` | `tests/validation/test_birdcage_sar_mass_averaged.py:163` |
  | `REFERENCE_IMPEDANCE_OHM` | `tests/validation/test_port_package_sparameters.py:140` |
  | `FREQUENCY_64_HZ`, `SALINE_SIGMA`, `SALINE_EPSILON_R` | `tests/validation/test_lossy_sphere_fullwave.py:92–94` |
  | `reduce_terminated_ports` | `fem_em_solver.ports.circuit:81` |
  | `CapacitorSheets`, `FourPortBirdcageGeometry`, `Study`, `StudyConfigError`, `cli`, `load_study`, `run_study_file` | `fem_em_solver` |
  | `_load`, `_study_from_mapping`, `_tune_summary` | `fem_em_solver.study_yaml` |

- **Printed beside, never asserted (transcribed with log lines):** item 463's (I-A) 1.391811e-05 and S₁₁ −2.787323415e-01 + 9.590409206e-01j (`20261015T201110Z_WF-16-s6b-W1.log:954`); (c2) 3.5268e-01 (`:971`); (a′) −6.139e-11 (`:950`).
- **Module fixture `runs` (solve), in this order, inside `pytest.MonkeyPatch.context()` wrapping `FourPortBirdcageGeometry.build` with item 463's build cache copied with a provenance comment (`tests/package/test_wf16_study_tune_capacitors.py:99–109`):**
  1. `_fresh_dir(LOGS / "f_small_hole_pec_tune_s13", comm)`; `run_t = run_study_file(TUNE_FILE, output_dir=…, comm=comm)` — item 442's route; the one real build; gives S0 (4×4) and `run_t.summary["tune"]`.
  2. `t2 = _tune_summary(run_t.study, {"k": 2}, comm)`; assert `run_t.study.result` unchanged (no re-solve).
  3. `_fresh_dir(LOGS / "f_small_hole_pec_capacitors", comm)`; `run_c = run_study_file(CAP_FILE, output_dir=…, comm=comm)`; `summary.json` read back on every rank after the `Barrier`.
  4. `py = Study(_geometry(), GelledSalinePhantomMaterial(SALINE_SIGMA, SALINE_EPSILON_R, FREQUENCY_64_HZ), CapacitorSheets(ITEM429_C_MODE_F, 0.0, on=("P2", "P3", "P4")), frequency_hz=FREQUENCY_64_HZ, degree=1)`; `py.solve(comm)`.
  5. Assert the cache read `real == 1`, `hits == 2`; print them, each `mesh_time` / `sweep_time` and the `run_study_file` wall times.
- **Nodes:** `test_wf1_s13_parse_and_refusals` (no mesh, `monkeypatch`, does not request `runs`); `test_wf1_s13_r0_fixture_tie`; `test_wf1_s13_summary_equals_same_study`; `test_wf1_s13_yaml_equals_python_study`; `test_wf1_s13_route_a_identity`; `test_wf1_s13_k2_control`.
- **Registry (lines added only, each block under `# WF-1 step 13 (§9 item 496)`):** the parse node under `[names.StudyConfigError]` and `[names.load_study]` after item 484's (`tests/package/validation_envelope.toml:267–268`, `:290–291`); the five solve nodes under `[names.run_study_file]` after item 484's block (`:358–362`); `test_wf1_s13_route_a_identity` and `test_wf1_s13_yaml_equals_python_study` also under `[names.CapacitorSheets]` after item 463's block (`:416–421`).

**Anchors (asserted; bands imported).**
- **(K) Parse, no mesh** (every geometry's `build` monkeypatched to raise, as `test_wf1_yaml_tune.py:80–85`):
  - `_spec_diff(load_study(CAP_FILE), py_spec) == []`, `py_spec` the fixture's step-4 Study built without solving; `type(load_study(CAP_FILE).excitation) is CapacitorSheets`, `.on == ("P2", "P3", "P4")`, `.c_f.hex() == ITEM429_C_MODE_F.hex()`.
  - `_spec_diff(load_study(CAP_FILE), load_study(PEC_CONTROL_FILE)) == ["excitation"]` — the file differs from the PEC control in the excitation only (catches a dropped section).
  - **(K-esr)** the mapping with `esr_ohm: 0.2` differs from the file by exactly `["excitation.esr_ohm"]`, and its `excitation.impedance_ohm(FREQUENCY_64_HZ)` equals `CapacitorSheets(ITEM429_C_MODE_F, 0.2, on=…).impedance_ohm(…)` by hex on both parts (ESR reaches the Study; nothing solved).
  - **(K-cu)** a mapping of `f_small_copper.yaml` (`tests/package/wf1/`) with the excitation replaced by `{kind: capacitor_sheets, c_f: C_F, esr_ohm: ESR_OHM, sheets: [P2]}` (the constants inserted in-test as Python floats, never transcribed) gives `_spec_diff(…, _study(CapacitorSheets(C_F, ESR_OHM, on=("P2",)))) == []` — item 373's copper ESR Study expressed from YAML, parse only.
  - The mapping with `esr_ohm` omitted gives `esr_ohm == 0.0` from the dataclass default (not filled by the parser).
  - `cli.main(["check", CAP_FILE])` rc (0, 0) via `_rc_all`; rank 0 prints `family: four_port_birdcage` and `excitation: CapacitorSheets`.
- **(X) Refusals, no mesh, each naming its field** (0 `build` calls):

  | Variant of `CAP_FILE`'s mapping | Field named |
  |---|---|
  | the tracked `malformed_capacitor_sheets_all_ports.yaml` (also `cli.main` `run` and `check`: rc (2, 2) via `_rc_all`, `out` absent on every rank by `MPI.LOR`) | `excitation.sheets` |
  | `sheets: [P5]` | `excitation.sheets[0]` |
  | `sheets: [P2, P2]` | `excitation.sheets` |
  | `sheets: []` | `excitation.sheets` |
  | `sheets: P2` (a string) | `excitation.sheets` |
  | `sheets` removed | `excitation.sheets` |
  | the key `True` in place of `sheets` (what `on:` loads as) — also through a YAML text written to a `_fresh_dir` scratch file | `excitation.on` |
  | `c_f: 0.0` / `c_f: "66e-12"` (string) / `c_f` removed | `excitation.c_f` |
  | `esr_ohm: -0.1` | `excitation.esr_ohm` |
  | `z0_ohm: 50.0` added | `excitation.z0_ohm` |
  | `geometry` of `sphere_mat4.yaml` with this excitation | `excitation.kind` |
  | `f_human_hole_copper_tune.yaml`'s mapping with this excitation (`tune` removed) | `excitation.kind` |
  | `tune: {k: 1}` added | `tune` |
  | `outputs: {sar: {mass_g: 10.0, centres: [[0.0, 0.0, 0.0]]}}` added | `outputs.sar` |
  | `outputs: {b1_plus: {points: [[0.0, 0.0, 0.0]]}}` added | `outputs.b1_plus` |

  Pinned refusals that must stay green unedited (W0): item 442's `tune` rows (`test_wf1_yaml_tune.py:132–136`), item 484's E1b / E4 rows, item 390a's `excitation.kind` rows (`test_wf1_yaml_f_human.py:70`, `:223`).
- **(R0) Fixture tie** (the S0 route): `run_t.summary["cells"] == run_c.summary["cells"] == ANS6_EXACT_STRIP_HOLE_CELL_RECORD` (80 798); `_class_means(run_t.result.s_matrix)` against `ANS6_EXACT_STRIP_CLASS_MEANS_P1["64"]["pec"]` per class ≤ `LEG_D0_REPRODUCTION_BAND` (prior 5.052e-11 / 3.878e-11 / 5.875e-11 on this route, `20261014T110326Z_WF-1-s11-W1.log`, and `…s6b-W1.log:949`); (a) |`run_t.summary["tune"]["c_mode_farad"]` / `TH17_STEP3_C_MODE_RECORD_F` − 1| ≤ `REDUCTION_BAND` (−2.976552e-4 backed, item 442); **(a′) custody of the transcribed C:** |`ITEM429_C_MODE_F` / `run_t.summary["tune"]["c_mode_farad"]` − 1| ≤ `REDUCTION_BAND` — backed: the same comparison on the same fixture and width read −6.139e-11 (`…s6b-W1.log:950`); also printed against `EXACT_IDENTITY_RTOL` (not asserted).
- **(I=) Summary = the same Study's result, bitwise.** Every `run_c.summary["s_matrix"]` part equals `run_c.study.result.s_matrix`'s by `float.hex`; the JSON read-back equals the in-memory summary by hex; `kept_ports == ["P1"]`; `capacitors.sheets == ["P2", "P3", "P4"]`; `capacitors.c_f` / `esr_ohm` equal the excitation's by hex; `capacitors.impedance_ohm` equals `run_c.study.excitation.impedance_ohm(FREQUENCY_64_HZ)` and `run_c.result.capacitor_impedance_ohm` by hex; `"touchstone"` and `"z_matrix"` absent; `len(s_matrix) == 1`.
- **(I-Py) The YAML Study = the Python `CapacitorSheets` Study on the same cached mesh:** |S_yaml[0,0] − S_py[0,0]| / |S_py[0,0]| ≤ `EXACT_IDENTITY_RTOL`. **Predicted 0** (one mesh object, equal specs by (K), the same arithmetic); hex equality printed, never asserted. Not backed by a prior measurement of this comparison.
- **(I-A) Item 463's identity through the file.** `zc = run_c.study.excitation.impedance_ohm(FREQUENCY_64_HZ)`; `pA = reduce_terminated_ports(run_t.result.s_matrix, float(REFERENCE_IMPEDANCE_OHM), {1: zc, 2: zc, 3: zc})` (0-based); |S_yaml[0,0] − pA[0,0]| / |pA[0,0]| ≤ `REDUCTION_BAND`. **Predicted 1.39e-05** (item 463 read 1.391811e-05 on the API route, `…s6b-W1.log:954`; the C differs by ≈ 6e-11 relative, so the digits should repeat to ≈ 4 places — printed beside, not asserted).

**Negative controls (rule (e)).**
- **(c1) The k = 2 capacitance on (I-A)'s comparator — *asserted, backed*.** `zc2 = CapacitorSheets(t2["c_mode_farad"], 0.0, on=("P2", "P3", "P4")).impedance_ohm(FREQUENCY_64_HZ)` (the excitation's own impedance, never a hand-written 1/(jωC)); `pA2` the Route-A reduction of S0 at `zc2`; |S_yaml[0,0] − pA2[0,0]| / |pA2[0,0]| > `REDUCTION_BAND`. Backing: the same comparison on the same fixture and width measured 3.5268e-01 (352.7× the band) through the API (`…s6b-W1.log:971`); (I-Py) ties the YAML S to the API's. One-sided, so no ceiling binds; the factor is printed. This is the control that shows (I-A)'s comparator reads C.
- **(c0) Dropped-key controls — *asserted, by construction, no solve*:** (K)'s `["excitation"]` diff against the PEC control (a parse that drops the section), (K-esr)'s `["excitation.esr_ohm"]` diff (a parse that drops ESR), and (I=)'s `kept_ports == ["P1"]` (a summary that ignores `sheets`).
- **Why no solved ESR control:** the solved file is r = 0 (the ESR ruling; item 463's Study); item 373's planted-ESR control (1.095344e-02, copper) is the API's and is not re-run. A solved ESR route through YAML is a separate step (scope).

**Windows.** Foreground, in the order W0, RG1, `pgrep`, RG2, `pgrep`, RG3, `pgrep`, W1, Wgen, Wval, Wapi, Wdoc. **The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped.** Bash `timeout: 660000`; pytest with `-s`; node IDs, never `-k`; `tests/environment/test_complex_mode.py` first. Run W0 only after every edit (rule (i)). Before RG2, RG3 and W1: `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reads 0. Two heavy windows never run at once.
- **W0** (no mesh, `-n 2`; the new parse node plus the rule-(c) parse / dispatch re-gates: items 442's and 484's parse nodes, `test_wf1_check_family.py` (the new good file joins its glob, `:32`), item 390a's and 431a's no-mesh modules including their `run_study_file` branch fakes (C7's dispatch head sits in their path)):
  `scripts/testing/run_and_log.sh WF-1-s13-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_capacitors.py::test_wf1_s13_parse_and_refusals tests/package/test_wf1_yaml_tune.py::test_wf1_s11_tune_parse_and_refusals tests/package/test_wf1_yaml_tune_f_human.py::test_wf1_s12_parse_and_refusals tests/package/test_wf1_check_family.py tests/package/test_wf1_yaml_f_human.py tests/package/test_wf1_yaml_f_human_thermal.py -v -s --tb=short'"`
  Price: item 484's W0 (same set less the new node) 19 passed, 7 s (`20261016T034124Z_WF-1-s12-W0.log:297–304`).
- **RG1** (rule (c): the `[names.load_study]` parse nodes not in W0) — item 484's RG1 verbatim (`20261016T034140Z_WF-1-s12-RG1.log:12`), relabelled; 13 passed, 4 s there:
  `scripts/testing/run_and_log.sh WF-1-s13-RG1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_study_yaml.py::test_wf1_good_files_build_the_python_studies tests/package/test_wf1_study_yaml.py::test_wf1_malformed_files_name_their_field tests/package/test_wf1_yaml_thermal.py::test_wf1_thermal_yaml_parse_and_refusals tests/package/test_wf1_yaml_copper.py::test_wf1_copper_yaml_parse_and_refusals tests/package/test_wf1_yaml_sar_b1.py::test_wf1_sar_b1_yaml_parse_and_refusals tests/package/test_wf1_yaml_probes.py::test_wf1_probes_yaml_parse_and_refusals tests/package/test_wf1_yaml_f_human.py::test_wf1_f_human_yaml_parse_and_refusals tests/package/test_wf1_yaml_copper_sar_b1.py::test_wf1_s8_parse_and_refusal_parity tests/package/test_wf1_yaml_f_human_thermal.py::test_wf1_fht_parse_and_refusals -v -s --tb=short'"`
- **RG2** (rule (c): the four-port branch, now `elif`, and `_parse_tune` / `_tune_summary`) — item 442's W1 verbatim (`20261014T110326Z_WF-1-s11-W1.log:12`), relabelled; 10 passed, 82.8 s at its last re-run (`20261016T034148Z_WF-1-s12-RG2.log:1853`):
  `scripts/testing/run_and_log.sh WF-1-s13-RG2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 300 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_tune.py -v -s --tb=short'"`
- **RG3** (rule (c): the sphere `else` branch, reached through C7's new head) — item 256's sphere solve node; its window with the full `tests/environment` and two parse nodes took 99 s (`20261009T004249Z_WF-1-s4-W2-wf1-w1.log:12`, `:431`), predicted 60–100 s here:
  `scripts/testing/run_and_log.sh WF-1-s13-RG3 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_study_yaml.py::test_wf1_cli_sphere_reproduces_wf15 -v -s --tb=short'"`
- **Not re-run (disclosed; question 2 for the review):** the three F-human `-n 8` solves (390b, 431b, 484 W1). No F-human branch line moves; C7 adds one `isinstance(study.excitation, CapacitorSheets)` test (False for every F-human Study) ahead of them, and W0's branch fakes exercise that head for the quadrature branches. Item 484's ruling required 390b because that item edited the F-human branch itself. If the review rules otherwise, item 484's W1 command (this file, item 484 "Windows", `FEM_EM_WF1_S12=1`, `-n 8`, `timeout -k 30 590`; 283 s measured, `20261016T034323Z_WF-1-s12-W1.log`) runs after W1, alone, under the label `WF-1-s13-RG4`, and the cost rises by ≈ 6 min.
- **W1** (heavy by ceiling, `-n 2`; predicted 55–110 s):
  `scripts/testing/run_and_log.sh WF-1-s13-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 300 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_capacitors.py -s -v --tb=short'"`
  - **Price.** One real PEC hole build through the YAML route: 29.7 s (mesh 22.7 s, sweep 6.5 s, `20261014T110326Z_WF-1-s11-W1.log:1830`); two reduced 1×1 `CapacitorSheets` Studies on the cached mesh at 6.1–6.6 s each (`…s6b-W1.log:941`, `:946`); start-up and `test_complex_mode.py` ≈ 20 s ⇒ ≈ 63 s. No new form: the PEC hole at degree 1 with lumped sheets compiled in items 429 / 442 / 463.
- **Wgen / Wval / Wapi / Wdoc** (the registry gains lines and two docstrings move: §9 "Every item" `WF-14` bullet; commit `docs/guide/api/` and `docs/guide/validation.md`) — item 484's four commands verbatim (`20261016T035359Z_WF-1-s12-Wgen.log:12`, `…035401Z…-Wval.log:12`, `…035404Z…-Wapi.log:12`, `…035426Z…-Wdoc.log:12`), labels `WF-1-s13-Wgen` / `-Wval` / `-Wapi` / `-Wdoc`:
  - `scripts/testing/run_and_log.sh WF-1-s13-Wgen "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 python3 scripts/docs/generate_api_reference.py --out docs/guide/api'"`
  - `scripts/testing/run_and_log.sh WF-1-s13-Wval "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace timeout -k 30 60 python3 scripts/docs/generate_validation_page.py --out docs/guide/validation.md'"` — predicted node count = the last Wval's 238 (`20261016T041116Z_WF-13-s7-Wval.log:41`) + 6 distinct new nodes = 244 if nothing else lands; print it.
  - `scripts/testing/run_and_log.sh WF-1-s13-Wapi "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace timeout -k 30 170 python3 -m pytest tests/package/test_wf14_validation_envelope.py tests/package/test_wf17_api_reference.py tests/package/test_wf14_validation_page.py -v -s --tb=short'"`
  - `scripts/testing/run_and_log.sh WF-1-s13-Wdoc "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && timeout -k 30 60 python3 scripts/testing/check_example_doc_references.py --docs-root docs/guide --check-links'"` (`0 dead`)
- **`git diff --stat main -- src/`** shows `src/fem_em_solver/study_yaml.py` only; say so in the commit (rule (c), C1–C7 named).

**Traps already paid for.**
- **YAML 1.1 booleans:** a bare `on:` key loads as `True` (found at drafting); the schema key is `sheets`. `c_f: 66e-12` (no dot) loads as a string — `_number` refuses it with `_YAML_FLOAT_HINT`; the tracked file writes `6.626796421e-11`.
- **`degree: 1` in the YAML and on the Python Study.** The `Study` default is 2; the degree-2 hole reads ≈ −4.03 % on (a) (item 445's (c-deg)). (K) catches a missing degree.
- **No `conductor` on the solved file.** The PEC route is `conductor` absent with `as_hole: true`, never a `LeontovichWall` at σ = ∞; κ is never passed.
- **`sheets` ids are `"P1"…"P4"`; `reduce_terminated_ports` indices are 0-based** (item 373 / 463 trap). A 1-based dict terminates the wrong sheets and (I-A) reads large.
- **`zc` from the excitation's own `impedance_ohm`**, never a hand-written 1/(jωC) (item 463 trap).
- **The cache must key on an equal geometry:** the YAML's `f_small` + `as_hole: true` equals `_geometry()` (item 442's (K), `_spec_diff == []`). Assert `real == 1, hits == 2`; a third real build means the specs diverged — stop and print the diff, never widen the cache key.
- **Never `float.hex` across builds** (item 442 ruling); here every solve is on one mesh object, (I=) compares one array with its own serialisation, (I-Py) is at `EXACT_IDENTITY_RTOL` with hex printed.
- **`_load`'s 4-tuple is load-bearing** (12+ callers); **`__all__` must not grow** (`WF-13`'s installed-surface snapshot).
- **Rank-local values:** cells and S come through the Study's reduced paths (`run_lumped_sheet_port_case` responses are global); rc and `out`'s existence are reduced (`_rc_all`, `MPI.LOR`); `femem` prints on rank 0 only.
- **Stubs are function-scoped:** the parse node uses `monkeypatch` and never requests `runs` (item 394's trap). The build-cache patch lives inside the fixture's `MonkeyPatch.context()` only.
- **The validation page reads the docstring's first paragraph for chunk IDs and bands from node bodies** (item 375): keep each band in the node that asserts it. **A guide page that pins a node count** (item 473's `len(nodes) == 3`, which turned item 484's first W0 red) — no guide page is edited here; do not cite the new nodes from `docs/guide/tuning_from_yaml.md`.
- **`-k a or b` splits into stray argv** inside the quoted container command: node IDs only.
- **Status 124 or a kill:** `pgrep -c python3` first; at 0 only, the rubric's orphan-stem sweep as one bare call labelled `WF-1-s13-cache-sweep`. Never raise a timeout or change `-n`.
- **Rule (i):** no edit to `study_yaml.py`, the module, the YAML files or the registry after the window that closes on them; a `src` edit after W0 means W0, RG1, RG2, RG3 again.
- Write ≈, never a tilde. No AED number anywhere.

**Scope.** Covered: `femem run` (`run_study_file`, `cli.py`'s `run` call) on one tracked file — the F-small exact-rectangle PEC hole (80 798 cells, κ off, `strip_fraction = GATED_WIDTH_FRACTION`), degree 1, 64 MHz, `-n 2`, lossless capacitor sheets on P2–P4 at the transcribed k = 1 `C_mode` — reproducing item 463's Route-A identity at `REDUCTION_BAND`, the summary equal to its own Study and to the Python Study; ESR and the copper conductor expressible at parse (item 373's Study, `_spec_diff == []`), not solved. **Not claimed:**
- a solved ESR or copper capacitor route through YAML (parse only; item 373 holds them through the API);
- the per-sheet (I-4) identity through YAML (item 463 holds it through the API);
- a "tune, then place" file — the file's `c_f` is a transcribed record, not derived in-run (one file cannot hold two Studies; a `c_f: {tune: …}` form would need a second solve inside `run_study_file` and is not drafted);
- a reduced-network Touchstone (`study.py:743–744`);
- `femem run --estimate` on a capacitor file (it prices the four-port's four right-hand sides, `study.py:1675–1678`; not checked) and the installed console script / `--ranks` on this file;
- degree 2, other frequencies, the F-human (refused at parse, as the API refuses);
- a tuned `f₀`, Q, `R_in`, match or absolute `S₁₁` / `Z_in` (`ANS-7` stays 🟡), a field-side mode selection (item 463: Route A does not separate k = 1 from k = 2);
- a stable schema; any AED reading.

**What a negative result means.** Report and stop. Park on `attempt/WF-1-step13-<UTC>`, mark item 496 BLOCKED with its unblock condition in the same commit (rule (d)), journal it in attempts.md.
- **(K) / (X) red:** a parse defect. Fix in-slot inside `study_yaml.py` only and re-run W0, RG1 (rule-(i) clock restarts).
- **W0's re-gate modules, RG1, RG2 or RG3 red:** the edit changed more than the new route. Revert and stop; never edit a pinned row to fit.
- **(R0) or (a) red:** the PEC route moved since items 442 / 463. Known-issues entry; never re-record.
- **(a′) red:** the transcribed `c_f` is not the in-run C — custody of the tracked file; print both by hex and stop.
- **(I=) red:** `run_study_file`'s capacitor branch or the JSON loses the reduced S; a defect in this item's `src`, fixed in-slot only if before W1's clock (rule (i)), else park.
- **(I-Py) red with (I=) and (K) green:** two equal specs on one mesh give different S beyond 1e-10 — a finding about the solve path's determinism; print both by hex, known-issues entry, never re-band.
- **(I-A) red with (R0), (I=) and (I-Py) green:** impossible unless the C or the reduction differs from item 463's; print the residual beside item 463's 1.391811e-05 and report.
- **(c1) inside the band:** a blind comparator (most likely `zc2` was not built from `t2`). Stop.
- **W1 status 124:** not a result; park with the measured time.

No `PREDICTED RED:` line: every assertion is predicted green; a red anywhere stops the slot.

**Status it can move.**
- **§7 `WF-1` row** stays ✅ (rule (j): its Done-when is step 1's). On green it gains: "**Step-13 record (§9 item 496, U57, provisional):** a `capacitor_sheets` excitation (`c_f`, `esr_ohm`, `sheets`) drives item 463's Route-A Study through `run_study_file` on the F-small PEC hole at degree 1, `-n 2`: the reduced S₁₁ reproduces the circuit reduction of the tune file's own 4×4 at <x> (`REDUCTION_BAND`), equals the same Study's result by `float.hex` and the Python `CapacitorSheets` Study on one mesh at <y> (`EXACT_IDENTITY_RTOL`); the k = 2 reduction misses by <p>; ESR and copper expressible at parse (item 373's Study, `_spec_diff == []`); bad `sheets` / `c_f` / `esr_ohm` / the YAML 1.1 `on` key / wrong geometry / `tune` / `sar` / `b1_plus` raise `StudyConfigError` naming their field, rc 2, nothing written; `src`: `study_yaml.py` only (rule (c), re-gates green)."
- **§10 U57** marked DONE (provisional until the weekly ratifies).
- **`WF-14`:** six nodes (two also under `[names.CapacitorSheets]`); API pages and `validation.md` regenerated.
- `WF-16`, `PORT-22`, `PORT-24`, `ANS-7` do not move; no §2 sentence moves; no known-issues entry on green.

**For the review.**
1. The schema key `sheets` (mapping to `CapacitorSheets.on`), forced by the YAML 1.1 `on` boolean — accept, or rule another name (a quoted `"on":` would work but fails silently-confusingly when unquoted; the brief refuses the `True` key explicitly either way).
2. The F-human `-n 8` solves are not re-run (C7 adds a head test ahead of untouched branches; W0's fakes cover the quadrature dispatch). Accept, or add RG4 (item 484's W1, ≈ 6 min).
3. The transcribed `c_f` (custody asserted by (a′)) instead of a `c_f: {tune: {k: 1}}` form — the latter is the user-facing "tune then place" flow but needs two solves in one `run_study_file` call; it is left for the weekly to enumerate if wanted.

**Spawn line:** "Use the implementer agent to execute `WF-1` §9 item 496 (U57) per `docs/planning/chunks/WF-1.md` (read the brief and any review ruling beneath it first; the rulings hold). The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run W0, RG1, RG2, RG3, W1, Wgen, Wval, Wapi, Wdoc in order, foreground, `timeout: 660000`, pytest with `-s`, container `timeout -k 30`; before RG2, RG3 and W1 run `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` and require 0; two heavy windows never run at once. `src` edits: `src/fem_em_solver/study_yaml.py` only, exactly C1–C7 plus the named docstrings, disclosed in the commit; never edit `study.py`, `_load`'s 4-tuple, `__all__`, or any existing test assertion. The YAML key is `sheets`, never a bare `on`; the YAML says `degree: 1`; no conductor on the solved file; never κ; `reduce_terminated_ports` indices 0-based; never import a `test_*` name; never compare across builds by `float.hex`. You are the executor — do not spawn agents, never background a harness window. On a red, park on `attempt/WF-1-step13-*`, mark item 496 BLOCKED in the same commit, and report."

### Review ruling 2026-10-16 10:30 (item 496)

**Queued, with question 2 ruled the stricter way.** Chain U's second provisional step (rule (6)); the 10-17 weekly ratifies or reverts.
1. **`sheets` is ratified** as the YAML key for `CapacitorSheets.on`, with the explicit refusal of a stray boolean `True` key naming `excitation.on`.
2. **RG4 is required** (item 484's W1, `-n 8`, ≈ 6 min), run after W1 and before Wgen, with `pgrep -c python3` = 0 first. C7's `if` → `elif` sits on the dispatch every F-human file passes through, and item 484's ruling re-ran 390b on that ground; standing rule (c) re-runs the pre-existing gate whose path a `src` change touches. If RG4 cannot start before minute 45, park the item on `attempt/WF-1-step13-*` and mark it BLOCKED (rule (d)); do not land the `src` change without it.
3. **The transcribed `c_f` is accepted**, with (a′) asserting custody against the in-run tune at `REDUCTION_BAND`. The step does not claim "tune then place" in one file; whether a `c_f: {tune: …}` form is wanted is the weekly's to enumerate.
- The YAML-vs-Python comparison at `EXACT_IDENTITY_RTOL` is predicted, not backed: a red there is reported as a finding, never widened.
- **Measured class:** chain step (solve), median ≈ 17.5, plus RG4's ≈ 6.

**Outcome 2026-10-16 (12:00 slot, CDT): green, DONE (provisional).** `src`: `study_yaml.py` only, C1–C7 plus the module and `run_study_file` docstrings; `study.py`, `_load`'s 4-tuple and `__all__` untouched. W0 13 passed in 3.2 s (`20261016T170354Z_WF-1-s13-W0.log`; the brief's set, which is item 484's W0 without `test_wf17_tuning_yaml_page.py`, plus the new node; the new good file joined `test_wf1_check_family.py`'s glob, 15 files, `family: four_port_birdcage`); RG1 13 passed (`…170420Z_WF-1-s13-RG1.log`); RG2 10 passed in 62.2 s (`…170432Z_WF-1-s13-RG2.log`); RG3 5 passed in 30.6 s (`…170543Z_WF-1-s13-RG3.log`). W1 `-n 2`: 10 passed in 41.1 s; build cache 1 real build, 2 hits; cells 80 798, mesh 21.8 s, sweep 6.3 s (tune file) / 6.5 s (capacitor file), Python Study 4.5 s; (K) diff `[]`, vs PEC control `['excitation']`, (K-esr) `['excitation.esr_ohm']` with impedance hex-equal, (K-cu) `[]`; every (X) row names its field; (R0) class means 5.052e-11 / 3.878e-11 / 5.875e-11 (the priors to the digit); (a) −2.976552e-04; (a′) `ITEM429_C_MODE_F / C_in-run − 1` = +6.139e-11 (item 463 printed the reciprocal sense, −6.139e-11); (I=) hex equal, `kept_ports == ['P1']`, no Touchstone / `z_matrix`; (I-Py) 2.223e-16 — predicted 0, **not** hex-equal (one ULP on one part; printed, not asserted: the two solves share the mesh but each builds its own solver); (I-A) 1.391811e-05, S₁₁ −2.787323414e-01 + 9.590409206e-01j (item 463 1.391811e-05); (c1) k = 2 C 9.873468900e-11 F, 3.5268e-01 (352.7× band) (`20261016T170623Z_WF-1-s13-W1.log`). RG4 (item 484's W1, `-n 8`, `pgrep` 0 first): 9 passed in 202.4 s (`20261016T170723Z_WF-1-s13-RG4.log`; item 484 read 283.0 s). Wgen / Wval (247 nodes, against the brief's predicted 244: three more nodes landed between drafting and this slot) / Wapi 22 passed / Wdoc 968 links, 0 dead (`20261016T171058Z`–`171129Z_WF-1-s13-W*.log`).

### §9 item 500 — 2026-10-16 (§10 chain U58, daily — provisional; chain step (solve))

*Drafted by a general-purpose agent for the 2026-10-16 18:00 daily review; the review reviews it.* Licence: §9 rule (6) (PROJECT_PLAN.md ≈ :3853), a provisional next step of live §10 chain U toward subgoal 9a (a study-level Python API + YAML/CLI), after U56 (§9 item 484) and U57 (§9 item 496), both in this file. **Chain U's provisional count becomes 3 of 3, the last the daily may write before the 10-17 weekly.** Chain U's validation target, quoted: "each step reproduces a gated result or asserts a gated identity through the API, within its unmoved band". **Target: reproduction plus an identity. No band, tolerance, record, Done-when or validation target is set or moved; every band is imported by name.** Anchors: `WF-15` step 3's records (§9 item 330, `docs/planning/chunks/WF-15.md` ≈ :411–471) — `MAT-17` step 2 cases A / B and `MAT-18` 2(b) (T1) at `.4e`, under `MAT-17`'s `INTERIOR_BAND` / `BALANCE_BAND` — read through `femem run`'s call, `run_study_file`, from three tracked files; and the YAML Study = the Python Study (`_spec_diff == []` at parse, θ at `L0_BAND` after the solve). Items 484 and 496 (above) and item 330 (`WF-15.md`) are the templates. No sheet, κ or capacitor-ESR ruling is engaged: there is no port, no capacitor and no field.

**Owner.** `WF-1` owns `femem run`'s schema (steps 1–13 are items 256 … 496). This is `WF-1` step 14. `WF-15` owns the API it drives (step 3, item 330) and does not move.

**Finding at drafting (code reading on `da7e369c`, nothing run): YAML cannot express a thermal-only Study today, and `study_yaml.py` alone cannot make `femem check` accept one. The smallest feasible scope is `study_yaml.py` (≈ 130 lines, all inserted) plus ≈ 8 inserted lines in `cli.py`'s `check`; `study.py` is untouched.** Five walls, in the order a file meets them:
- **Sections.** `_REQUIRED_SECTIONS = ("study", "geometry", "material", "excitation")` (`study_yaml.py:129`, checked at `:847–849`) and `_parse_study` requires `study.frequency_hz` (`:235`). A thermal-only Study has no material (`Study.__init__` raises `TypeError` on anything but `None`, `study.py:117–118`) and no frequency (`self.frequency_hz = None`, `:127`).
- **Geometry / excitation kinds.** `_parse_geometry`'s kinds are `_EXCITATION_FOR`'s keys (`study_yaml.py:134–139`, `:254`): no thermal-only kind; `UniformHeatSource` is not an excitation kind (`:375`) and none of `ThermalBallGeometry`, `UniformHeatSource`, `TransientSpec`, `DutyCycle`, `ContainerWall` is imported (`:107–118`).
- **Thermal.** On a `sphere_in_box`, `thermal:` admits only `k_w_per_m_k`, `h_w_per_m2_k`, `t_inf_c`, `source_scale`, `degree` (`:564`). **Verified: a transient section exists only on `f_human_birdcage`** (`_parse_thermal_f_human`, `:598–683`, keys `dt_s`, `n_steps`, `rho_c_j_per_m3_k`, `theta`, `startup_half_steps`); there is no `duty` key anywhere and no wall key.
- **Run.** `run_study_file` calls `study.solve(comm)` unconditionally (`:1150`) — `Study.solve` raises `TypeError` on a thermal-only Study (`study.py:171–172`) — and then `float(study.frequency_hz)` (`study_yaml.py:1157`).
- **Check (outside `study_yaml.py`).** `cli.main(["check", f])` prints `m.sigma` / `m.epsilon_r` of `study.material` on rank 0 (`cli.py:159–164`): `AttributeError` on `None`. `tests/package/test_wf1_check_family.py:32` globs every non-`malformed_` file in `tests/package/wf1/`, so a good thermal file there turns that re-gate red. **This is the reason for the `cli.py` lines** (question 1 for the review; the alternative is below).
- **Nothing else is missing.** The API does the work: `Study(ThermalBallGeometry(...), None, UniformHeatSource(q))` and `Study.temperature(RobinCooling, transient=TransientSpec, duty=DutyCycle | None)` dispatch to `_study_transient_temperature` (`study.py:1246–1247`, `:1434–1492`), which builds the mesh itself, refuses `source_scale != 1` (`:1443–1444`) and returns a `TransientStudyResult` (`:1401–1418`) carrying `mesh`, `k` / `rho_c` / `q`, `source_multipliers`, `source_power` and `max_balance_residual`.

**The change.** Every line is inserted; no existing line of `study_yaml.py` or `cli.py` is edited (docstring additions are inserted lines). Standing rule (c) still applies because two insertions sit on the path every file takes (T2, T5): disclose them and re-run the gates of every path they touch (Windows).
- **(T1)** `study_yaml.py`'s `from .study import (…)` block gains `ContainerWall`, `DutyCycle`, `ThermalBallGeometry`, `TransientSpec`, `UniformHeatSource` (five inserted lines). `__all__` unchanged; `public_surface.txt` (`WF-13`) must not move.
- **(T2)** In `_parse_mapping`, **after** the unknown-top-key loop (`:842–846`) and **before** the required-sections loop (`:847`): `geo = mapping.get("geometry")`; `if isinstance(geo, dict) and geo.get("kind") == "thermal_ball": return _parse_thermal_ball_mapping(mapping, source)`. The version and top-key checks stay first for every file; any other file falls through to the unchanged path (a missing `geometry` still raises at `geometry`, `:849`).
- **(T3)** New `_parse_thermal_ball_mapping(mapping, source)`, returning `_load`'s 4-tuple `(Study, "thermal_ball", outputs, _BallThermal)`; validation order fixed, first error raises:
  - required sections `geometry`, `excitation`, `thermal` (missing ⇒ at that name); then refused sections `study`, `material`, `conductor`, `tune` (present ⇒ at that name, reason "a thermal_ball study has no frequency / EM degree / density" (`study`), "… no EM material: k is `thermal.k_w_per_m_k`, ρc is `thermal.transient.rho_c_j_per_m3_k`" (`material`), "… no conductor" / "… no tune route");
  - `geometry` (the `ThermalBallGeometry` field names, as the four-port route uses its dataclass's): `_unknown` (`kind`, `radius`, `resolution`, `wall`), `_missing` (`radius`, `resolution`), each `_number`, finite and > 0, else at `geometry.<key>`; `wall` optional, a mapping (else at `geometry.wall`) with exactly `thickness_m`, `k_w_per_m_k`, `rho_c_j_per_m3_k` (`_unknown` / `_missing` at `geometry.wall.<key>`), each finite and > 0;
  - `excitation`: `_kind(sec, ("uniform_heat_source",), "excitation", source)`; `_unknown` (`kind`, `q_w_per_m3`); `_missing` (`q_w_per_m3`); finite and > 0;
  - `thermal`: `_unknown` (`k_w_per_m_k`, `h_w_per_m2_k`, `t_inf_c`, `degree`, `source_scale`, `transient`, `duty`); `source_scale` present ⇒ at `thermal.source_scale` ("applies to an EM source; scale `excitation.q_w_per_m3` instead" — the API's own refusal, `study.py:1443–1444`, raised at parse); `transient` absent ⇒ at `thermal` ("a thermal_ball Study needs `thermal.transient`" — the API's `TypeError`, `study.py:1439–1440`); the base keys through `_parse_thermal(base, "sphere_in_box", source)` (the F-human route's own reuse, `:621–622`); `transient` **copied** from `_parse_thermal_f_human`'s block (`:624–658`, same five keys, same checks, same field names; not factored out — a shared helper would edit that path); `duty` optional, a mapping with exactly `period_s` (finite > 0) and `duty` (in (0, 1]) at `thermal.duty.<key>`;
  - `outputs` optional, `directory` only (`_unknown` names `outputs.sar` / `outputs.b1_plus` / `outputs.probes`);
  - builds `Study(ThermalBallGeometry(radius, resolution, wall=ContainerWall(...) | None), None, UniformHeatSource(q))` with **no** keyword (the Python Study's defaults: degree 2, density 1000, frequency `None`). Omitted optional keys are **not** filled in (the dataclass defaults apply, as everywhere in this module).
- **(T4)** New frozen `_BallThermal(cooling: RobinCooling, transient: TransientSpec, duty: Optional[DutyCycle])` beside `_FHumanThermal` (`:583–595`).
- **(T5)** In `run_study_file`, **after** the `digest` line (`:1148`) and **before** `result = study.solve(comm)` (`:1150`): `if isinstance(study.geometry, ThermalBallGeometry): return _run_thermal_ball(study, kind, cooling, path, out, digest, started, comm)`.
- **(T6)** New `_run_thermal_ball(...)`: `thermal = study.temperature(cooling.cooling, transient=cooling.transient, duty=cooling.duty, comm=comm)`; summary `femem`, `source`, `sha256`, `family` (`"thermal_ball"`), `ranks`, `cells` (`thermal.mesh.topology.index_map(3).size_global`), and `temperature = {dt_s, n_steps, theta, startup_half_steps, rho_c_j_per_m3_k, thermal_degree, q_w_per_m3, duty: {period_s, duty} | null, wall: {thickness_m, k_w_per_m_k, rho_c_j_per_m3_k} | null, source_multipliers: [...], source_power_w, max_balance_residual, theta_max_k}` (`theta_max_k` = owned max of Re(T − t_inf), `MPI.MAX`-reduced — `_temperature_summary`'s four lines, `:1023–1030`, copied; that function needs a `convected_power` a transient result does not have), `elapsed_s`; rank 0 writes `summary.json`, `Barrier`; returns `StudyRun(study=study, result=None, summary=…, paths=…, thermal=thermal)`. No `degree` / `frequency_hz` key (disclosed: the file has none). No Touchstone.
- **(T7, `cli.py`)** In `main`'s `check` branch, **between** `study, kind = _load(args.file)[:2]` (`:157`) and `if comm.rank == 0:` (`:158`): `if kind == "thermal_ball":` — rank 0 prints `family: thermal_ball` and `excitation: UniformHeatSource`, flushes — `return 0`. ≈ 8 inserted lines; no printed line of any other family moves (`test_wf1_check_family.py::test_wf1_check_other_lines_unchanged` pins the sphere and four-port lines against the pre-change `cli.py`).
- **Docstrings (inserted lines):** a `thermal_ball` bullet in the module docstring (sections, keys, the refusals, "the solve legs are gated in the real build"); one sentence in `run_study_file`'s ("a `thermal_ball` study runs `Study.temperature` with no field solve and writes no Touchstone"). Wgen regenerates the API pages.
- **Never:** touch `study.py`; edit an existing line of `study_yaml.py` or `cli.py`; widen `_load`'s 4-tuple; add a name to `__all__`; admit `transient` / `duty` on a `sphere_in_box` (the API refuses, `study.py:1249–1250`); admit `t_init_c` (not in the F-human transient section either; not claimed); list `thermal_ball` in the unchanged `geometry.kind` unknown-kind message (that would edit `:254`; disclosed, question 3).
- **Alternative if the review declines T7:** `study_yaml.py` only, with the three good files under `tests/package/wf1/thermal_ball/` (outside check-family's non-recursive glob) and `femem check` on them not claimed — a user-facing `AttributeError` left in place, filed in known-issues. The drafter does not recommend it.

**`WF-1` step 14: `femem run` (`run_study_file`) on three tracked thermal-ball files reproduces `WF-15` step 3's records — `MAT-17` step 2 cases A and B (duty cycle) and `MAT-18` 2(b) (container wall) — in the real build at `-n 2`; each run's θ equals the same Python Study's at `L0_BAND`; dropping the duty cycle or the wall's ρc through the file misses.**
- **Executor:** implementer.
- **Tier:** standard for the new legs (W1a / W1b, **real build** — the records' build, item 330 ruling (4); container `timeout -k 30 180`), `-n 2`; heavy for RG4 only (`-n 8`, complex). Every re-gate runs in the build its gate was recorded in.
- **Class:** chain step (solve), median ≈ 17.5 measured min.
- **Cost: ≈ 56 predicted min** = 15 fixed + ≈ 22 code (≈ 130 `study_yaml.py` lines, ≈ 8 `cli.py` lines, module ≈ 320 lines, three good YAML files and one malformed, registry lines) + W0 ≈ 0.3 + RG1 ≈ 0.2 + W9 ≈ 0.3 + RG2 ≈ 1.2 + RG3 ≈ 0.6 + RG5 ≈ 0.8 + W1a ≈ 1.8 + W1b ≈ 2 + RG4 ≈ 3.5–5 + Wgen / Wval / Wapi / Wdoc ≈ 0.6 + ≈ 8 row / journal / spare. **Park (rule (d)) if the `src` change and W0 / RG1 / W9 are not green by minute 33, if W1a has not started by minute 40, or if RG4 cannot start by minute 50** — never land the `src` change without RG4 (the 10:30 ruling's terms). Item 496 measured ≈ 15–20 min against its 42 predicted (12:00 slot, `20261016T170354Z` … `171129Z_WF-1-s13-*.log`).
- **Family:** `WF-1` step 14, attempt 1 of 4.
- **Depends on:** nothing open. Item 330's module is on `main` and supplies every helper and record constant. On a rebase over another item touching `study_yaml.py`, `cli.py`, `validation_envelope.toml` or the generated pages, keep both sides and re-run W0, RG1, W9, Wgen / Wval / Wapi / Wdoc.

**Tracked inputs (four new files under `tests/package/wf1/`).** Each starts with a `# WF-1 step 14 (§9 item 500): …` line, then `femem: 1`. Floats are written as Python `repr` digits so (K) can assert them by `float.hex`; **never `4.18e6`** (PyYAML 1.1 loads it as a string; `_YAML_FLOAT_HINT`, `study_yaml.py:130–132`).
- **`thermal_ball_duty_a.yaml`** (`MAT-17` case A): `geometry: {kind: thermal_ball, radius: 0.1, resolution: 0.005}`; `excitation: {kind: uniform_heat_source, q_w_per_m3: 2000.0}`; `thermal: {k_w_per_m_k: 0.6, h_w_per_m2_k: 10.0, t_inf_c: 20.0, transient: {dt_s: 150.0, n_steps: 94, rho_c_j_per_m3_k: 4180000.0}, duty: {period_s: 1200.0, duty: 0.5}}`. (`n_steps` = `_n(T_END_A, 150.0)` = 94, `T_END_A` = 14 100 s, `tests/validation/test_mat17_duty_cycle_fem.py:59`, `:87–90`; the log's `N=94`, `20261008T213834Z_WF-15-s3-W1a-duty.log:91`.)
- **`thermal_ball_duty_b.yaml`** (case B): as A with `transient: {dt_s: 435.4166666666668, n_steps: 20, …}` and `duty: {period_s: 13933.333333333338, duty: 0.25}` — `DT_B = D_B · P_B / 8` (`test_mat17_duty_cycle_fem.py:61`), `P_B = 0.2 · T_R` (`test_mat17_duty_cycle_closed_form.py:43–44`), `_n(T_END_B, DT_B)` = 20 (`…W1a-duty.log:96`). Repr digits computed at drafting; (K) asserts them by hex against the imported constants.
- **`thermal_ball_wall.yaml`** (`MAT-18` 2(b)): `geometry: {kind: thermal_ball, radius: 0.1, resolution: 0.005, wall: {thickness_m: 0.009999999999999995, k_w_per_m_k: 0.19, rho_c_j_per_m3_k: 1700000.0}}`; `excitation: {kind: uniform_heat_source, q_w_per_m3: 2000.0}`; `thermal: {k_w_per_m_k: 0.6, h_w_per_m2_k: 10.0, t_inf_c: 20.0, transient: {dt_s: 348.3333333333334, n_steps: 40, rho_c_j_per_m3_k: 4180000.0}}` (no `duty`). `thickness_m` is `B_ - A_` = 0.11 − 0.1 in binary, **not** `0.01` (`test_wf15_study_transient.py:284`; both give `outer_radius` 0.11, but `0.01` makes (K) read `["geometry.wall.thickness_m"]`); `dt_s` = `T_STAR / N_CN` (`test_mat18_two_layer_sphere_transient_fem.py:61–62`).
- **`malformed_thermal_ball_rho_c_string.yaml`:** `thermal_ball_duty_a.yaml` with `rho_c_j_per_m3_k: 4.18e6` → `thermal.transient.rho_c_j_per_m3_k` (the string-float trap made a fixture).

**Test module: new `tests/package/test_wf1_yaml_thermal_ball.py`.**
- **Docstring.** First paragraph backticks only `WF-1`, `WF-15`, `MAT-17`, `MAT-18` (the (F3) rule; never `U58` or an item number). It says the solve legs are gated in the real build, the build of the records.
- **Opt-in solve:** `FEM_EM_WF1_S14=1`, item 484's `_knob()` pattern (unset ⇒ the two solve nodes skip, any other value raises); when set, assert `not IS_COMPLEX` (`tests/complex_mode.py:23`) and `MPI.COMM_WORLD.size == 2` in the node, never skip on a mismatch.
- **Imports (module-level names and constants only; never a `test_*` name):**

  | Name | Source |
  |---|---|
  | `LOGS`, `WF1`, `_cmp`, `_fresh_dir`, `_mapping`, `_rc_all`, `_spec_diff` | `tests/package/test_wf1_study_yaml.py:85–170` |
  | `RES`, `MAT17_A_PRINT`, `MAT17_B_PRINT`, `MAT18_T1_PRINT`, `MAT18_T1_CELLS`, `MAT18_C1_UNIFORM_RHOC`, `_duty_study`, `_wall_study`, `_exact_duty`, `_on`, `_same_mesh`, `_theta_owned_rel` | `tests/package/test_wf15_study_transient.py:84–179`, `:283–288` |
  | `R`, `K`, `Q`, `H`, `T_INF`, `RHO_C`, `_interior_error` | `tests/validation/test_mat9_transient_conduction_sphere.py:41–46`, `:109` |
  | `CASES` | `tests/validation/test_mat17_duty_cycle_closed_form.py:44` |
  | `DT_B`, `T_END_A`, `T_END_B`, `ON_TIME_A`, `ON_TIME_B`, `INTERIOR_BAND`, `BALANCE_BAND`, `_census`, `_n` | `tests/validation/test_mat17_duty_cycle_fem.py:59–69`, `:87`, `:105` |
  | `A_`, `B_`, `K1`, `K2`, `Q1`, `RC1`, `RC2` | `tests/validation/test_mat18_two_layer_sphere_closed_form.py:37–39` |
  | `CELL_GUARD`, `N_CN`, `T_STAR`, `_errors`, `_fine` | `tests/validation/test_mat18_two_layer_sphere_transient_fem.py:61–63`, `:114`, `:160` |
  | `L0_BAND` | `tests/validation/test_mat5_em_coupled_sphere.py:104` |
  | `_dg0_owned_diff` | `tests/materials/test_phantom_material_model.py:246` |
  | `IS_COMPLEX` | `tests/complex_mode.py:23` |
  | `ContainerWall`, `DutyCycle`, `RobinCooling`, `Study`, `StudyConfigError`, `ThermalBallGeometry`, `TransientSpec`, `UniformHeatSource`, `cli`, `load_study`, `run_study_file` | `fem_em_solver` |
  | `_load`, `_study_from_mapping` | `fem_em_solver.study_yaml` |

- **Printed beside, never asserted (transcribed with log lines):** `WF-15` step 3's (G) readings 4.3869840018e-04 (`20261008T213834Z_WF-15-s3-W1a-duty.log:91`), 1.0682780964e-03 (`:96`), 8.2066325346e-04 (`20261008T214046Z_WF-15-s3-W1b-wall.log:91`); its (c1) 9.7087e-01 (`…W1a-duty.log:98`); its cells 151 548 (`…W1a-duty.log:89`).
- **Nodes:** `test_wf1_s14_parse_and_refusals` (no mesh, `monkeypatch`, both builds); `test_wf1_s14_duty_reproduces_mat17` (W1a); `test_wf1_s14_wall_reproduces_mat18_2b` (W1b).
- **Registry (lines added only, each block under `# WF-1 step 14 (§9 item 500)`):** the parse node under `[names.StudyConfigError]` and `[names.load_study]` after item 496's (`tests/package/validation_envelope.toml:269–270`, `:294–295`); both solve nodes under `[names.run_study_file]` after item 496's block (`:367–372`); the parse and duty nodes under `[names.ThermalBallGeometry]`, `[names.TransientSpec]` and `[names.DutyCycle]` (`:223–247`); the duty and wall nodes under `[names.UniformHeatSource]` (`:230–234`); the wall node under `[names.ContainerWall]` (`:249–252`).

**Anchors (asserted; bands imported).**
- **(K) Parse, no mesh** (`ThermalBallGeometry.build` and every other geometry's `build` monkeypatched to raise, as `test_wf1_yaml_tune.py:80–85`; 0 `build` calls asserted):
  - `_spec_diff(load_study(DUTY_A), Study(ThermalBallGeometry(R, RES), None, UniformHeatSource(Q))) == []`, the same for `DUTY_B` (one Study spec: the transient and the schedule are `temperature()`'s, not the Study's), and `_spec_diff(load_study(WALL), Study(ThermalBallGeometry(A_, RES, wall=ContainerWall(B_ - A_, K2, RC2)), None, UniformHeatSource(Q1))) == []`; `frequency_hz is None` and `material is None` on all three.
  - The fourth `_load` element, through `_cmp` (floats by hex): case A `cooling == RobinCooling(K, H, T_INF)`, `transient == TransientSpec(150.0, _n(T_END_A, 150.0), RHO_C)`, `duty == DutyCycle(*CASES["A"])`; case B `TransientSpec(DT_B, _n(T_END_B, DT_B), RHO_C)`, `DutyCycle(*CASES["B"])`; wall `RobinCooling(K1, H, T_INF)`, `TransientSpec(T_STAR / N_CN, N_CN, RC1)`, `duty is None`.
  - `cli.main(["check", f])` rc (0, 0) via `_rc_all` on the three files; rank 0 prints exactly one `family: thermal_ball` line and `excitation: UniformHeatSource`.
- **(X) Refusals, no mesh, each naming its field** (`StudyConfigError`, `.field` equal; 0 `build` calls):

  | Variant of `DUTY_A`'s mapping unless named | Field named |
  |---|---|
  | the tracked `malformed_thermal_ball_rho_c_string.yaml` (also `cli.main` `run` and `check`: rc (2, 2) via `_rc_all`, `out` absent on every rank by `MPI.LOR`) | `thermal.transient.rho_c_j_per_m3_k` |
  | `study: {frequency_hz: 64000000.0}` added | `study` |
  | `material: {sigma: 0.5, epsilon_r: 78.0}` added | `material` |
  | `conductor: {kind: leontovich}` added / `tune: {k: 1}` added | `conductor` / `tune` |
  | `thermal` removed / `excitation` removed | `thermal` / `excitation` |
  | `geometry.radius: 0.0` / `geometry.resolution` removed / `geometry.radius_m` added | `geometry.radius` / `geometry.resolution` / `geometry.radius_m` |
  | `geometry.wall: [0.01]` / `geometry.wall` with `thickness_m: -0.01` / with `k2` added | `geometry.wall` / `geometry.wall.thickness_m` / `geometry.wall.k2` |
  | `excitation.kind: quasi_static_sphere_field` / `excitation.q_w_per_m3: 0.0` | `excitation.kind` / `excitation.q_w_per_m3` |
  | `thermal.source_scale: 2.0` / `thermal.balls` added | `thermal.source_scale` / `thermal.balls` |
  | `thermal.transient` removed | `thermal` |
  | `thermal.transient.n_steps: 0` / `theta: 1.5` / `startup_half_steps: 1` / `t_init_c: 25.0` added | `thermal.transient.n_steps` / `.theta` / `.startup_half_steps` / `.t_init_c` |
  | `thermal.duty.duty: 0.0` / `1.5` / `thermal.duty.period_s: -1.0` / `thermal.duty` without `duty` | `thermal.duty.duty` / `.duty` / `thermal.duty.period_s` / `thermal.duty.duty` |
  | `outputs: {probes: […]}` / `outputs: {sar: {…}}` added | `outputs.probes` / `outputs.sar` |
  | `sphere_mat10_thermal.yaml`'s mapping with `thermal.duty` added / `thermal.transient` added / `excitation.kind: uniform_heat_source` (the unchanged sphere route) | `thermal.duty` / `thermal.transient` / `excitation.kind` |

  Pinned refusals that must stay green unedited (W0 / RG1): `malformed_thermal_on_port.yaml` → `thermal`, item 431a's F-human `thermal.transient` rows, item 496's `excitation.on` / `excitation.sheets` rows, item 484's E1b / E4 rows, item 390a's `excitation.kind` rows.
- **(UT) Bitwise upstream of the solve, the file's Study against the Python one** (each solve node; the Python Study is item 330's own helper, `_duty_study(dt, n, DutyCycle(P, D))` / `_wall_study(RC2)`, run in the same process): `summary["cells"]` equals both runs' global cell count, `_same_mesh(run.thermal.mesh, py.mesh)` geometry equal; duty: `source_multipliers` `array_equal` between the two and equal to `summary["temperature"]["source_multipliers"]`; wall: `summary["cells"] == MAT18_T1_CELLS` (215 542), `_same_mesh` also against `_fine()`'s mesh (the comparand's), and DG0 k / ρc / q `_dg0_owned_diff` against `_fine()`'s `k_f` / `rc_f` / `q_f` max |Δ| `== 0.0`. **(S0)** `_census(run.thermal, dt, 2, ON_TIME_*, …)` (its own 1e-12 bar, `test_mat17_duty_cycle_fem.py:116`).
- **(I-Py) The file's θ = the Python Study's:** `_theta_owned_rel(_on(run.thermal.mesh, py.temperature), run.thermal.temperature, comm) ≤ L0_BAND` for A, B and the wall. **Predicted 0 to ≈ 2e-14** (same specs by (K), same code path; item 330's API-vs-direct read 8.146e-15 / 2.080e-14 / 1.433e-14 on a different code path, `…W1a-duty.log:90`, `:95`, `…W1b-wall.log:90`). Not backed by a prior measurement of this comparison; hex equality of the owned arrays printed, never asserted (item 496's (I-Py) read one ULP, not hex-equal).
- **(G) The gates through the file:** case A / B interior relative L2 against `_exact_duty(run.thermal.mesh, T_END_*, P, D)` (`_interior_error`) ≤ `INTERIOR_BAND`; wall (T1) interior (`_errors(msh_f, ct, _on(msh_f, run.thermal.temperature), theta)`, the middle value) ≤ `INTERIOR_BAND` with `_fine()`'s guard `nglob <= CELL_GUARD` asserted **before** any solve; `summary["temperature"]["max_balance_residual"] ≤ BALANCE_BAND` on every run.
- **(R) Records at `.4e`:** `f"{e:.4e}"` equals `MAT17_A_PRINT` (`4.3870e-04`), `MAT17_B_PRINT` (`1.0683e-03`), `MAT18_T1_PRINT` (`8.2066e-04`); each printed at `.10e` beside item 330's reading. **Predicted equal to ≈ 1e-13 relative** (item 330 read the API equal to the direct run to all ten digits).
- **(I=) Summary = its own thermal result:** `summary["temperature"]["theta_max_k"]` equals the reduced owned max recomputed from `run.thermal.temperature` by `float.hex`; `max_balance_residual` and `source_power_w` equal `run.thermal`'s by hex; the JSON read back on every rank (after the `Barrier`) equals the in-memory summary by hex; `summary["family"] == "thermal_ball"`; `"frequency_hz"`, `"degree"`, `"touchstone"` absent; `run.result is None` and `run.study.result is None` (no field solve).

**Negative controls (rule (e)).**
- **(c1) The schedule dropped through a file — *asserted, backed*.** Rank 0 writes `DUTY_A`'s mapping without `thermal.duty` (`yaml.safe_dump`, repr floats) to `_fresh_dir(LOGS / "thermal_ball_c1")/c1.yaml`, `Barrier`; `run_study_file` on it; its case-A interior error against case A's duty comparand `> INTERIOR_BAND`, and `summary["temperature"]["duty"] is None`, `source_multipliers == []`. Backing: the same comparison on the same fixture measured 9.7087e-01 through the API (`20261008T213834Z_WF-15-s3-W1a-duty.log:98`; (I-Py) ties the file to the API). Predicted 9.7087e-01; one-sided, so no ceiling binds (the continuous θ is ≈ 1/D times the duty θ, so the miss saturates near 1/D − 1 = 1.0).
- **(c2) The wall's ρc dropped through a file — *asserted, backed*.** The same with `WALL`'s `geometry.wall.rho_c_j_per_m3_k: 4180000.0` (RC1): (T1) interior `> INTERIOR_BAND`. Backing: 5.3140e-02 through the API on the same fixture (`20261008T214046Z_WF-15-s3-W1b-wall.log:93`; `MAT-18`'s own prior, `20261003T184850Z_MAT-18-step2b-W1b.log:94`, is `MAT18_C1_UNIFORM_RHOC`). Predicted 5.3140e-02, printed beside the constant; one-sided.
- **(c0) Dropped-key controls — *asserted, by construction, no solve*:** (K)'s fourth element of the c1 mapping reads `duty is None` against case A's `DutyCycle` (a parser that drops `duty` is caught at parse); `_spec_diff` of `WALL` against its mapping without `geometry.wall` reads `["geometry.wall"]`; `_spec_diff(DUTY_A, DUTY_B) == []` while their `transient` / `duty` differ by `_cmp` in exactly `dt_s`, `n_steps`, `period_s`, `duty` (the schedule rides in the thermal element, not the Study).
- **Not re-run:** item 330's (c3) planted-spacing comparator (deterministic, `…W1b-wall.log:89`); `_dg0_owned_diff` is the same function.

**Windows.** Foreground, in the order W0, RG1, W9, `pgrep`, RG2, `pgrep`, RG3, `pgrep`, RG5, `pgrep`, W1a, W1b, `pgrep`, RG4, Wgen, Wval, Wapi, Wdoc. **The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped.** Bash `timeout: 660000`; pytest with `-s`; node IDs, never `-k`; `tests/environment/test_complex_mode.py` first (in the real build it passes `test_scalar_type_matches_active_dolfinx_build` and skips the rest, `…W1a-duty.log:57–60`). Run W0 only after every edit (rule (i)). `pgrep` = `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reads 0. Two heavy windows never run at once.
- **W0** (complex, no mesh, `-n 2`; the new parse node plus the rule-(c) parse / dispatch re-gates — T2 sits on every file's parse path and T5 on every run's: items 442's, 484's and 496's parse nodes, `test_wf1_check_family.py` (the three good files join its glob, 15 → 18 files; T7 is gated here), items 390a's and 431a's no-mesh modules including their `run_study_file` branch fakes (T5's head sits in their path)) — item 496's W0 (`20261016T170354Z_WF-1-s13-W0.log:12`) plus one node:
  `scripts/testing/run_and_log.sh WF-1-s14-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_thermal_ball.py::test_wf1_s14_parse_and_refusals tests/package/test_wf1_yaml_capacitors.py::test_wf1_s13_parse_and_refusals tests/package/test_wf1_yaml_tune.py::test_wf1_s11_tune_parse_and_refusals tests/package/test_wf1_yaml_tune_f_human.py::test_wf1_s12_parse_and_refusals tests/package/test_wf1_check_family.py tests/package/test_wf1_yaml_f_human.py tests/package/test_wf1_yaml_f_human_thermal.py -v -s --tb=short'"`
  Price: 13 passed in 3.2 s there (`…170354Z_WF-1-s13-W0.log:270`); predicted 14 passed, < 10 s.
- **RG1** (rule (c): every `[names.load_study]` parse node not in W0, including `test_wf1_yaml_thermal.py`'s sphere `thermal:` refusals) — item 496's RG1 verbatim (this file, item 496 "Windows"), label `WF-1-s14-RG1`; 13 passed there (`20261016T170420Z_WF-1-s13-RG1.log`).
- **W9** (rule (c): T7 edits the installed `femem`) — item 401's installed-CLI window verbatim (`20261012T172921Z_WF-1-s9-W9-installed-femem.log:12`), label `WF-1-s14-W9-installed-femem`; 1 passed in 5.0 s there (`:116`):
  `scripts/testing/run_and_log.sh WF-1-s14-W9-installed-femem "docker compose exec -T fem-em-solver bash -lc 'timeout -k 10 50 bash /workspace/scripts/package/wf13_setup.sh real && mkdir -p /tmp/wf13/run/tests/package/wf1 && cp /workspace/tests/package/test_wf1_installed_cli.py /tmp/wf13/run/tests/package/ && cp /workspace/tests/package/wf1/*.yaml /tmp/wf13/run/tests/package/wf1/ && cd /tmp/wf13/run && env -u PYTHONPATH WF13_INSTALLED_CHECK=1 timeout -k 10 90 /tmp/wf13/venv/bin/python -m pytest -p no:cacheprovider tests/package/test_wf1_installed_cli.py -v -s --tb=short'"`
- **RG2** (rule (c): the four-port branch behind T5's head) — item 496's RG2 verbatim (`20261016T170432Z_WF-1-s13-RG2.log:12`), label `WF-1-s14-RG2`; 10 passed in 62.2 s (`:1853`).
- **RG3** (rule (c): the sphere `else` branch behind T5's head) — item 496's RG3 verbatim (`20261016T170543Z_WF-1-s13-RG3.log:12`), label `WF-1-s14-RG3`; 5 passed in 30.6 s (`:214`).
- **RG5** (rule (c): the capacitor branch, now the first branch behind T5's head) — item 496's W1 verbatim (`20261016T170623Z_WF-1-s13-W1.log:12`), label `WF-1-s14-RG5`; 10 passed in 41.1 s (`:992`).
- **W1a** (real build, `-n 2`; predicted 70–110 s: item 330's duty leg was 54.8 s for the API + direct + (c1) runs, `…W1a-duty.log:99`; here the file + Python + (c1) runs, the same three solves; start-up and `test_complex_mode.py` ≈ 15 s — the full `tests/environment`, ≈ 65 s of item 330's 123 s, is not on the path):
  `scripts/testing/run_and_log.sh WF-1-s14-W1a-duty "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace FEM_EM_WF1_S14=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_thermal_ball.py::test_wf1_s14_parse_and_refusals tests/package/test_wf1_yaml_thermal_ball.py::test_wf1_s14_duty_reproduces_mat17 -s -v --tb=short'"`
- **W1b** (real build, `-n 2`; predicted 85–125 s: `_fine()` build 8.7 s and comparand 17.1 s (`…W1b-wall.log:85`, `:87`), three wall solves at ≈ 14.3 s each (`:91`), start-up ≈ 15 s; 2.04 GiB summed there, `:94`):
  `scripts/testing/run_and_log.sh WF-1-s14-W1b-wall "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace FEM_EM_WF1_S14=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_thermal_ball.py::test_wf1_s14_wall_reproduces_mat18_2b -s -v --tb=short'"`
  - **JIT:** the real-build transient forms compiled in item 431a's RG12 (`20261014T010524Z_WF-1-s10a-RG12-wf15-transient.log`); allow ≈ 30 s if they hash differently.
- **RG4 — required, and why.** T5's head runs ahead of `study.solve` for every file, so it sits on the dispatch every F-human file passes through — the ground on which the 10-15 18:00 ruling re-ran item 390b (item 484's E5 was also a pure insertion) and the 10-16 10:30 ruling required RG4 for item 496's C7. The drafter applies the same rule rather than argue the head is inert; the review may waive it, citing W0's F-human branch fakes, but no prior ruling has. Item 484's W1 verbatim (`20261016T170723Z_WF-1-s13-RG4.log:12`), label `WF-1-s14-RG4`, after W1b, alone, `pgrep` 0 first; 9 passed in 202.4 s there (`:7389`), 283.0 s at item 484's own W1.
- **Not re-run (disclosed):** items 390b / 431b's `-n 8` solves (the quadrature and F-human thermal branches; W0's branch fakes cover the head on their path and RG4 the real F-human dispatch, as at U57); item 330's own module (`study.py` is untouched).
- **Wgen / Wval / Wapi / Wdoc** (registry lines added and two docstrings grow: §9 "Every item" `WF-14` bullet; commit `docs/guide/api/` and `docs/guide/validation.md`) — item 496's four commands verbatim (this file, item 496 "Windows"), labels `WF-1-s14-Wgen` / `-Wval` / `-Wapi` / `-Wdoc`. Wval's node count = the last Wval's 247 (item 496 outcome, `20261016T1711…Z_WF-1-s13-Wval.log`) + 3 = 250 if nothing else lands; print it. Wdoc reads `0 dead`. No guide page is edited (no page claims thermal-only studies are Python-only; drafter grepped `docs/guide/`).
- **`git diff --stat main -- src/`** shows `src/fem_em_solver/study_yaml.py` and `src/fem_em_solver/cli.py` only, insertions only; say so in the commit (rule (c), T1–T7 named).

**Traps already paid for.**
- **The build.** The records are real-build `.4e` strings (item 330 ruling (4)); W1a / W1b never `source … dolfinx-complex-mode` and never set `FEM_EM_REQUIRE_COMPLEX`; the module asserts `not IS_COMPLEX` when the knob is set. The re-gates stay in their own (complex) build.
- **YAML floats:** `4.18e6` / `1.7e6` load as strings (write `4180000.0`, `1700000.0`); `B_ - A_` is `0.009999999999999995`; case B's `dt_s` / `period_s` are 17-digit reprs. (K) by hex is the custody check — a hand-rounded value reads as a spec diff, never "close enough".
- **The comparand is evaluated at dof coordinates in chunks of 20 000** (`_exact_duty` from item 330's module, `MAT-18`'s `_fine()`), never `f.interpolate(callable)` — each leg paid one 180 s kill for it (`20261003T170330Z_MAT-17-step2-W1.log`, `20261003T184348Z_MAT-18-step2b-W1.log`).
- **UFL cannot mix meshes:** compare θ across runs only through `_on` (which asserts the P1 dof coordinates equal first) or in NumPy; evaluate case A / B on the file run's own mesh, the wall on `_fine()`'s.
- **`_study_transient_temperature` builds on `MPI.COMM_WORLD` whatever `comm` says** (`study.py:1446`): pass `MPI.COMM_WORLD`; a sub-communicator is not claimed.
- **No `degree` in a thermal-ball file:** `study:` is refused; the thermal degree is `thermal.degree` (`RobinCooling`'s, default 1), never `Study.degree`.
- **`femem` prints on rank 0 only:** read `check`'s stdout on rank 0; reduce rc through `_rc_all`; `out`'s existence by `MPI.LOR`.
- **Stubs are function-scoped:** the parse node uses `monkeypatch` and builds nothing; the solve nodes never patch `build`.
- **A guide page that pins a node count** (item 473's `len(nodes) == 3`, which turned item 484's first W0 red): cite no new node from any guide page.
- **`-k a or b` splits into stray argv** inside the quoted container command: node IDs only.
- **Never `float.hex` across builds** (item 442 ruling): (I=) compares one result with its own serialisation; (I-Py) is at `L0_BAND` with hex printed.
- **`_load`'s 4-tuple is load-bearing** (12+ callers, `cli.py:157` among them); **`__all__` must not grow.**
- **Status 124 or a kill:** `pgrep -c python3` first; at 0 only, the rubric's orphan-stem sweep as one bare call labelled `WF-1-s14-cache-sweep` — **in the build that was killed** (the real and complex caches are separate). Never raise a timeout or change `-n`.
- **Rule (i):** no edit to `study_yaml.py`, `cli.py`, the module, the YAML files or the registry after the window that closes on them; a `src` edit after W0 means W0, RG1, W9, RG2, RG3, RG5 again.
- Write ≈, never a tilde. No AED number anywhere.

**Scope.** Covered: `femem run` (`run_study_file`, `cli.py`'s `run` call) and `femem check` on three tracked thermal-ball files — the uniform-source duty-cycled ball (`MAT-17` step 2 cases A / B) and the two-region walled ball (`MAT-18` 2(b)), CN + Rannacher, P1 thermal, real build, `-n 2` — reproducing `WF-15` step 3's records at `.4e` inside `INTERIOR_BAND` / `BALANCE_BAND`, each equal to the same Python Study at `L0_BAND`. **Not claimed:**
- an EM-driven transient, duty or wall through YAML (`thermal.transient` / `thermal.duty` on a `sphere_in_box` stay refused, as the API refuses them; `MAT-17` step 3 / `MAT-18` step 3 are not `WF-15`'s either);
- `femem run --estimate` on a thermal-ball file — by code reading (not run) `Study.estimate` builds the ball mesh and then raises `AttributeError` at `self.geometry.sphere_tag` (`study.py:1680–1687`), which `cli._estimate` does not catch (`cli.py:107–112`); a `study.py` refusal is a follow-on for the weekly (lane 1), not this item;
- the installed console script running these files (W9 re-gates the installed `check` on the old files only), `--ranks`;
- `t_init_c`, θ ≠ ½, a thermal degree other than 1, a sub-communicator, the complex build for the solve legs;
- perfusion, `MAT-20`'s running maximum or probes on the ball, a time series in the summary;
- any time-averaging, regulatory or compliance statement, a real container's properties, tissue, a stable schema.

**What a negative result means.** Report and stop. Park on `attempt/WF-1-step14-<UTC>`, mark item 500 BLOCKED with its unblock condition in the same commit (rule (d)), journal it in attempts.md.
- **(K) / (X) red:** a parse defect. Fix in-slot inside `study_yaml.py` (or T7 in `cli.py`) only and re-run W0, RG1, W9 (rule-(i) clock restarts).
- **W0's re-gate modules, RG1, W9, RG2, RG3, RG5 or RG4 red:** the insertion changed more than the new route. Revert and stop; never edit a pinned row to fit.
- **(UT) red:** the file builds a different problem from the Python Study. Print the `_spec_diff` and the first differing DG0 entry; stop.
- **(I-Py) red with (UT) and (K) green:** two equal specs on two builds give θ beyond `L0_BAND` — a solve-floor finding (known-issues table), band unmoved; print both by hex.
- **(R) red with (G), (UT) and (I-Py) green:** an environment drift since 2026-10-14's RG12 — for `record-reconciler` (1*), never re-registered in-slot.
- **(G) red:** `MAT-17` / `MAT-18` / `WF-15` no longer reproduce in this build — a known-issues entry; the review re-reads item 330's record before anything moves.
- **(I=) red:** `_run_thermal_ball` or the JSON loses the result; a defect in this item's `src`, fixed in-slot only before W1a's clock (rule (i)), else park.
- **(c1) or (c2) inside the band:** the knob does not reach the solver through the file. Stop.
- **Collection error in the real build** (an import the real image lacks): a finding; park — never move the solve legs to the complex build.
- **W1a / W1b status 124:** not a result; park with the measured time.

No `PREDICTED RED:` line: every assertion is predicted green; a red anywhere stops the slot.

**Status it can move.**
- **§7 `WF-1` row** stays ✅ (rule (j): its Done-when is step 1's, item 256; nothing here re-reads it). On green it gains: "**Step-14 record (§9 item 500, U58, provisional):** a `thermal_ball` geometry (`radius`, `resolution`, `wall`), a `uniform_heat_source` excitation and `thermal.transient` / `thermal.duty` drive `Study.temperature` through `run_study_file` with no field solve, real build, `-n 2`: `MAT-17` step 2 cases A / B and `MAT-18` 2(b) reproduce at <a> / <b> / <c> (`.4e` records, `INTERIOR_BAND`), each run's θ equal to the Python Study's at ≤ <d> (`L0_BAND`), multipliers and DG0 coefficients bitwise; the duty and wall-ρc drops miss by <p> / <q>; `study` / `material` / `conductor` / `tune` / `source_scale` / the string-float ρc and bad wall / duty / transient keys raise `StudyConfigError` naming their field; `femem check` prints `family: thermal_ball`; `src`: `study_yaml.py` + `cli.py`'s `check`, insertions only (rule (c), re-gates green incl. RG4)."
- **§10 U58** marked DONE (provisional until the weekly ratifies). Chain U's daily provisional budget is then spent until the 10-17 weekly.
- **`WF-14`:** three nodes (also under five `WF-15` names); API pages and `validation.md` regenerated.
- `WF-15`, `MAT-17`, `MAT-18`, `WF-13` do not move; no §2 sentence moves; no known-issues entry on green.
- **For the weekly:** a `WF-17` tutorial page "duty cycle and container wall from YAML" and an `EX-*` companion become openable once U58 is ratified; the `--estimate` `AttributeError` above is a lane-1 candidate. Neither is opened here.

**For the review.**
1. **Scope beyond `study_yaml.py`:** T7's ≈ 8 inserted lines in `cli.py`'s `check` (with W9), forced by `cli.py:159–164` reading `material.sigma` and by check-family's glob — accept, or rule the subdirectory alternative (study_yaml.py only, `femem check` on thermal files not claimed, known-issues entry).
2. **RG4** is drafted as required (≈ 4–5 min) on the 10-15 / 10-16 rulings' ground; the review may waive it, saying why.
3. **The unknown-`geometry.kind` message** stays unchanged (it does not list `thermal_ball`), so no existing line moves; listing it would edit `study_yaml.py:254` and is the review's to allow.
4. **Size.** ≈ 56 predicted is near the slot. If the review wants slack, the smallest re-scope that keeps the target is to drop RG4 (if waived) or to split the step as items 390a / 390b were: 500a (T1–T7, the parse node, W1a — the duty leg) and 500b (the wall solve node, W1b), counted as one chain step.

**Spawn line:** "Use the implementer agent to execute `WF-1` §9 item 500 (U58) per `docs/planning/chunks/WF-1.md` (read the brief and any review ruling beneath it first; the rulings hold). The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Run W0, RG1, W9, RG2, RG3, RG5, W1a, W1b, RG4, Wgen, Wval, Wapi, Wdoc in order, foreground, every harness Bash call with `timeout: 660000`, pytest with `-s`, container `timeout -k 30`; before RG2, RG3, RG5, W1a and RG4 run `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` and require 0; two heavy windows never run at once. W1a and W1b run in the **real** build (no `dolfinx-complex-mode`, no `FEM_EM_REQUIRE_COMPLEX`) with `FEM_EM_WF1_S14=1`; every re-gate runs in its own recorded build. `src` edits: `src/fem_em_solver/study_yaml.py` (T1–T6) and `src/fem_em_solver/cli.py` (T7) only, insertions only, disclosed in the commit; never edit `study.py`, an existing line of either file, `_load`'s 4-tuple, `__all__`, or any existing test assertion. Write floats as repr digits (`4180000.0`, `0.009999999999999995`, never `4.18e6`); evaluate comparands in chunks at dof coordinates, never `interpolate(callable)`; never import a `test_*` name; never compare across builds by `float.hex`. You are the executor — do not spawn agents, never background a harness window. On a red, park on `attempt/WF-1-step14-*`, mark item 500 BLOCKED in the same commit, and report."

#### Review ruling 2026-10-16 18:00 (item 500; overrides the brief where they differ)

(1) **The `cli.py` insertion is accepted** (≈ 8 inserted lines in `check`, no existing line edited) under §9 standing rule (c): disclosed in the commit, re-gated by W9 and `test_wf1_check_family.py` in the same slot. The subdirectory alternative is declined — it leaves a user-facing crash. `study.py` stays untouched. (2) **RG4 is required** (item 484's `-n 8` F-human W1), as for U57: the new early return runs before `study.solve` on every file. (3) **One item, the brief's park rules** (code + W0 / RG1 / W9 green by minute 33; W1a started by minute 40; RG4 started by minute 50). On a park, the split the brief names (500a parser + `check` + duty leg; 500b wall leg) is the next review's re-pose, attempt 2 of 4. (4) **The `--estimate` defect on a thermal-ball file is read from source, not observed;** it is not claimed, gets no known-issues entry until a window shows it, and is named for the 10-17 weekly as a lane-1 follow-on. (5) Chain U's third provisional step (U56, U57, U58): the 10-17 weekly ratifies or reverts all three.

**Outcome 2026-10-16 (19:30 slot, CDT): green, DONE (provisional).** `src`: `study_yaml.py` (T1–T6, 233 inserted lines) + `cli.py` (T7, 6 inserted lines), insertions only; `study.py`, `_load`'s 4-tuple and `__all__` untouched. W0 14 passed in 3.4 s (`20261017T003617Z_WF-1-s14-W0.log`; check-family's glob 18 files; re-run after the test edit below, 14 passed, `…004126Z_WF-1-s14-W0.log`); RG1 13 passed (`…003637Z`); W9 1 passed in 5.1 s (`…003641Z`); RG2 10 passed in 62.5 s (`…003710Z`); RG3 5 passed in 32.0 s (`…003824Z`); RG5 10 passed in 40.4 s (`…003905Z`). (K) every `_spec_diff` `[]`, cooling / transient / duty hex-equal (case B `dt_s` `0x1.b36aaaaaaaaadp+8`); every (X) row names its field; (c0) `['geometry.wall']` and A vs B thermal diff exactly `dt_s`, `n_steps`, `period_s`, `duty`; 0 build calls. **W1a, first run red on the test, not the code** (`20261017T004010Z_WF-1-s14-W1a-duty.log:108`): (I=)'s whole-dict JSON compare included `elapsed_s`, which is each rank's own wall clock (only rank 0's summary is written, as on every route); the compare now excludes `elapsed_s` on ranks ≠ 0 and keeps it on rank 0 (measurement recorded in the module comment). W1a re-run (real build, `-n 2`): 3 passed in 65.5 s (`20261017T004140Z_WF-1-s14-W1a-duty.log`): cells 151 548; (S0) rel 0; (I-Py) A 8.146e-15, B 4.622e-14 (predicted ≤ ≈ 2e-14 — above the prediction, far inside `L0_BAND`; not hex-equal, printed); (G) A 4.3869840018e-04, B 1.0682780964e-03 (item 330's ten digits each); balance 2.823e-14 / 9.975e-15; (c1) 9.7087e-01. W1b: 2 passed in 73.1 s (`20261017T004254Z_WF-1-s14-W1b-wall.log`): cells 215 542 = record, DG0 k / ρc / q max |Δ| 0.0; (I-Py) 5.973e-15; (T1) 8.2066325346e-04 (item 330's digits); balance 1.082e-14; (c2) 5.3140e-02. `pgrep` 0, RG4 9 passed in 210.4 s (`20261017T004425Z_WF-1-s14-RG4.log`). Wgen 44 files / Wval 260 nodes (brief predicted 250; more nodes landed between drafting and this slot) / Wapi 22 passed / Wdoc 1 013 links, 0 dead (`20261017T004807Z`–`004837Z_WF-1-s14-W*.log`).

## §7 row moved verbatim 2026-10-17 (rotate_plan_archive.py chunks)

 `femem run`: the MRI study CLI and YAML config | ✅ 2026-10-06 (item 256, U5): `femem run` builds Studies equal to the Python-built ones by `float.hex` and reproduces `MAT4_DEG1_RECORD` and `P1` at their bands; steps 2–8 and 10–11 are records at their imported bands — the installed console script, `thermal:` (now on the tuned copper F-human at `-n 8`, items 431a / 431b, reproducing P17's θ(t\*) and `MAT-20`'s six-minute maximum at `REUSE_REPRODUCTION_RTOL`), `conductor:`, `outputs.sar` / `outputs.b1_plus` / `outputs.probes`, the F-human two-port schema (items 390a / 390b, 394) and a top-level `tune: {k}` section driving `Study.tune` (item 442). Step-12 record (§9 item 484, U56, provisional): a top-level `tune: {k: 1}` on the copper F-human hole `lumped_port_drive` file drives the ring read through `run_study_file` at degree 1, `-n 8`, one solve: `C_mode` 1.181185799070e-10 F, −3.580e-13 from `WF-12` step 3's H12c reading (`REDUCTION_BAND`), equal to the same Study's `tune()` at 0.000e+00 (`EXACT_IDENTITY_RTOL`); 32×32 tied to the H12c npz at 1.616e-12; k = 2 misses by −33.27 %; `src`: `study_yaml.py` only (rule (c), re-gates green). **Step-13 record (§9 item 496, U57, provisional):** a `capacitor_sheets` excitation (`c_f`, `esr_ohm`, `sheets`) drives item 463's Route-A Study through `run_study_file` on the F-small PEC hole at degree 1, `-n 2`: the reduced S₁₁ reproduces the circuit reduction of the tune file's own 4×4 at 1.391811e-05 (`REDUCTION_BAND`), equals the same Study's result by `float.hex` and the Python `CapacitorSheets` Study on one mesh at 2.223e-16 (`EXACT_IDENTITY_RTOL`); the k = 2 reduction misses by 3.5268e-01 (352.7× band); ESR and copper expressible at parse (item 373's Study, `_spec_diff == []`); bad `sheets` / `c_f` / `esr_ohm` / the YAML 1.1 `on` key / wrong geometry / `tune` / `sar` / `b1_plus` raise `StudyConfigError` naming their field, rc 2, nothing written; `src`: `study_yaml.py` only (rule (c), re-gates green incl. item 484's `-n 8` W1 as RG4). **Step-14 record (§9 item 500, U58, provisional):** a `thermal_ball` geometry (`radius`, `resolution`, `wall`), a `uniform_heat_source` excitation and `thermal.transient` / `thermal.duty` drive `Study.temperature` through `run_study_file` with no field solve, real build, `-n 2`: `MAT-17` step 2 cases A / B and `MAT-18` 2(b) reproduce at 4.3870e-04 / 1.0683e-03 / 8.2066e-04 (`.4e` records, `INTERIOR_BAND`), each run's θ equal to the Python Study's at ≤ 4.622e-14 (`L0_BAND`), multipliers and DG0 coefficients bitwise; the duty and wall-ρc drops miss by 9.7087e-01 / 5.3140e-02; `study` / `material` / `conductor` / `tune` / `source_scale` / the string-float ρc and bad wall / duty / transient keys raise `StudyConfigError` naming their field; `femem check` prints `family: thermal_ball`; `src`: `study_yaml.py` + `cli.py`'s `check`, insertions only (rule (c), re-gates green incl. RG4). Items 431a / 431b / 442 / 484 / 496 / 500 are provisional until the weekly ratifies. *History: `docs/planning/chunks/WF-1.md`.* 

### §9 item 521 — 2026-10-18 (§10 U59, femem run on a thermal-only gel box; chain step, solve)

*Drafted by a read-only agent for the 2026-10-18 03:00 daily review; the review reviews it. Code read on `d94320db`; nothing was run.*

**Licence.** §10 chain U, U59 (`PROJECT_PLAN.md:4848–4854`), enumerated by the 2026-10-17 weekly and flagged for this review's queue (`:5025`): "`femem run` on a thermal-only gel-box file (`geometry.kind: thermal_box`, the F2182 dimensions, Robin walls, `uniform_heat_source`, `thermal.transient` / `thermal.duty`): `MAT-21` step 1's m = 4 steady anchors and P26's duty record reproduced at `INTERIOR_BAND` / `E_BAND`; θ = the Python Study at `L0_BAND`; a wrong-wall control misses. A disclosed `study_yaml.py` admission with every YAML gate re-run. Serial on nothing." Chain U's target: each step reproduces a gated result or asserts a gated identity through the API, within its unmoved band. **No band, tolerance, record, Done-when or target is set or moved; every band is imported by name.** Templates: item 500 (U58, this file), item 330 (`WF-15` step 3, `docs/planning/chunks/WF-15.md:411–473`), item 431 (drafted as 431a + 431b, this file). No sheet, κ or ESR ruling is engaged: no port, no capacitor, no field.

**Finding at drafting: the Python `Study` has no box geometry, so U59 as enumerated is two slots, not one.**
- **No thermal box exists anywhere in `src/`.** The only thermal-only geometry is `ThermalBallGeometry` (`src/fem_em_solver/study.py:1325–1354`: a gmsh ball, optional `ContainerWall`). `grep -rn "thermal_box\|ThermalBox" src tests` returns nothing. Every thermal-only guard is an `isinstance(…, ThermalBallGeometry)`: the constructor (`study.py:114–132`), `solve` (`:171–172`, which also guards `mean_sar` / `sar` / `b1_plus` / `tune` because they call `solve`), `_refuse_thermal_only` (`:533–537`, `estimate`) and the `temperature` dispatch (`:1253–1254`). A box passed today falls through every one of them to the EM path.
- **No steady thermal-only solve exists.** On a thermal-only Study `temperature()` goes only to `_study_transient_temperature`, which refuses `transient=None` (`study.py:1446–1447`). `MAT-21` step 1's m = 4 anchors are **steady**. So the Python Study needs a steady branch for the box before any file can reach those anchors.
- **The YAML side is item 500's pattern again:** the `thermal_ball` hook in `_parse_mapping` (`src/fem_em_solver/study_yaml.py:1014–1016`), `_parse_thermal_ball_mapping` (`:732–850`, which refuses a missing `thermal.transient` at `:785–786`), `_run_thermal_ball` (`:1295–1351`, which reads `study.geometry.wall` at `:1303`), the `run_study_file` dispatch (`:1380–1381`) and `femem check`'s `thermal_ball` branch (`src/fem_em_solver/cli.py:158–163`). None of them admits a box.
- **"θ = the Python Study at `L0_BAND`" presupposes a Python box Study.** It does not exist, so the comparison has no right-hand side until the API half lands.
- **Size.** Item 500 was ≈ 56 predicted with `study.py` untouched, and measured ≈ 30 minutes of windows and code. U59 needs that much again plus a `study.py` geometry, a steady branch, five guards, an export, and the `WF-13` / `WF-15` re-gates the export triggers.

**Proposed split (one chain step, two items, as 390a / 390b and 431a / 431b were).**
- **521a — `WF-15` step 5 (the API half).** `ThermalBoxGeometry` and the box through `Study.temperature` (steady, transient, duty). The Python Study reproduces `MAT-21` step 1's m = 4 steady anchors and P26's duty record, bitwise against the direct call up to the solve. Owner: `WF-15`, which owns `Study.temperature` / thermal-only (steps 2–3, items 260 / 330).
- **521b — `WF-1` step 15 (the file half, U59's letter).** The `thermal_box` YAML admission; `femem run` reproduces the same records; θ_file = θ_Python at `L0_BAND`; the wrong-wall control through a file. Serial on 521a.
- If the review declines the split, 521a's text is the first half of a single item whose park rule is "521a's windows green by minute 40, else park the YAML half"; that is the drafter's second choice.

**Build — a deviation from the letter, for the review (question 1).** U59 says "real build". `MAT-21` step 1 and step 2 (P26) are **complex-build** records: every gate node is `@complex_only` (`tests/validation/test_mat21_box_fem.py:311–396`, `tests/validation/test_mat21_box_duty_fem.py:341–347`), and the logs source `dolfinx-complex-mode` (`20261013T142646Z_MAT-21-step1r-W1-steady.log`, `20261014T003506Z_MAT-21-step2-a2-W3-closed-form.log:12`). U58 ran real because *its* records (`MAT-17` / `MAT-18` through item 330) are real-build (item 330 ruling (4)). Item 442's ruling forbids `float.hex` across builds, and item 500's own rule is "every re-gate runs in the build its gate was recorded in". **The drafter poses the solve legs of 521a and 521b in the complex build**, where the `.4e` records were taken. The alternative is the real build with the `.4e` records printed beside and only the bands asserted.

**Tier — a deviation from the letter (question 2).** U59 says "standard". The m = 4 steady leg is standard: step 1's m = 4 rung built in 0.3 s, solved in 0.6 s and interpolated its P2 comparand in 4.3 s (`20261013T142646Z_MAT-21-step1r-W1-steady.log:99–102`). **P26's duty record is not standard.** Its (T-duty) comparand is the ten-run triple-series sum interpolated at 380 725 P2 points: the fixture took 177.8 s and the window 262 s harness (`20261014T003506Z_MAT-21-step2-a2-W3-closed-form.log:92`, `:229`), over the 180 s ceiling. So the duty leg is **heavy by measurement** (container `timeout -k 30 590`, `-n 2`), in both items.

---

#### 521a — `WF-15` step 5: `ThermalBoxGeometry` through `Study.temperature`

- **Executor:** implementer. **Tier:**
  - W0 standard (no mesh, `-n 2`, complex);
  - W1 standard (`-n 2`, complex; the m = 4 steady leg, predicted 25–45 s);
  - W2 **heavy** (`-n 2`, complex; the m = 6 duty leg, predicted 210–290 s against P26 W3's 262 s);
  - the re-gates in their own builds (below).
- **Class:** chain step (solve). **Cost:** ≈ 18 measured (class median) / ≈ 55 predicted min = 15 fixed + ≈ 20 code (≈ 90 `study.py` lines, module ≈ 330 lines, registry and snapshot lines) + windows ≈ 15.5 (W0 0.2, W1 0.7, W2 4.5, RG1 2.1, RG2 2.2, RG3 2.0, RG4 1.8, RG5 0.4, RG6 0.5, RG7 0.6, Wgen / Wval / Wapi / Wdoc 0.4) + ≈ 5 journal. **Park (rule (d)) if code + W0 + W1 are not green by minute 30, or W2 has not started by minute 35.**
- **Depends on:** nothing open. **Family:** `WF-15` step 5, attempt 1 of 4.

**Code state (cited above).** The box solves the drafter reuses, all on `main`: `_exterior_robin_tags` (`study.py:1176–1184`), a verbatim copy of `MAT-10`'s `_robin_tags` (`tests/validation/test_mat10_sar_to_heat_sphere.py:92–98`, `ROBIN_TAG = 1` at `:56`; `study.py`'s `_ROBIN_TAG = 1` at `:1153`), which is the facet tagging `MAT-21`'s `_rung` uses (`tests/validation/test_mat21_box_fem.py:148–162`); `solve_steady_conduction` (the call at `study.py:1272–1281` and at `test_mat21_box_fem.py:159–160`, same keywords); and `_study_transient_temperature` (`study.py:1441–1499`), which already takes scalar k / ρc / q when `geo.wall is None` (`:1455–1457`), the `fem.Constant` path of `MAT-21`'s `_cn` (`test_mat21_box_duty_fem.py:96–101`).

**The change (`src/fem_em_solver/study.py` and `src/fem_em_solver/__init__.py`; insertions only, disclosed under rule (c)).**
- **(A1) `ThermalBoxGeometry(half_widths_m: tuple, divisions: tuple)`**, a frozen dataclass placed after `ThermalBallGeometry` (`:1354`). `__post_init__` coerces both to tuples of length 3 (floats finite > 0; ints ≥ 1, else `ValueError` naming the field).
  - `build(comm)` returns `(mesh, cell_tags, facet_tags)`: `dmesh.create_box(comm, [np.array([-A, -B, -C]), np.array([A, B, C])], list(divisions), cell_type=dmesh.CellType.tetrahedron)`, exactly `_rung`'s call (`test_mat21_box_fem.py:151–153`); cell tag 1 on every local cell; facets from `_exterior_robin_tags(msh)`.
  - A read-only property `wall` returning `None` (docstring: "no container wall on a box in v0"). This lets `_study_transient_temperature` run **unedited** on its scalar branch. Disclosed: the property exists only for that branch.
- **(A2) Four guard insertions**, each a new `if isinstance(…, ThermalBoxGeometry):` block placed **before** the existing ball check, so that no existing line moves and the pinned ball messages stay byte-identical (the estimate message is pinned exactly as `THERMAL_MESSAGE` in `tests/package/test_wf13_cli_estimate.py`, item 507):
  - the constructor (before `:114`): the thermal-only initialisation copied from `:117–131`, with messages naming `ThermalBoxGeometry`;
  - `solve` (before `:171`): `TypeError("a thermal-only Study has no field solve; use temperature()")`, the ball's text;
  - `_refuse_thermal_only` (before `:535`): `TypeError(f"{method} is not defined on a thermal-only Study (ThermalBoxGeometry has no EM solve; use temperature())")`;
  - the `temperature` dispatch (before `:1253`): `return _study_box_temperature(self, cooling, transient, duty)`.
- **(A3) `_study_box_temperature`.**
  - `transient is None` and `duty is not None` ⇒ `ValueError` ("duty needs a transient").
  - `transient is None` ⇒ steady. Check `cooling` as `:1260–1261` does, refuse `source_scale != 1` as `:1450–1451` does, build, then `solve_steady_conduction(msh, ft, k=…, q=float(q), h=…, t_inf=…, robin_tag=_ROBIN_TAG, degree=cooling.degree)`. Return a `ThermalStudyResult` (`:1139–1150`) with `submesh=msh`, `parent_cells=None`, `facet_tags=ft`, `heat_source=q` and the solver's `source_power` / `convected_power`.
  - Otherwise ⇒ `return _study_transient_temperature(self, cooling, transient, duty)`, unchanged.
- **(A4) Export.** `"ThermalBoxGeometry"` joins `study.py`'s `__all__` (`:35–51`) and `fem_em_solver/__init__.py`'s import and `__all__` (beside `ThermalBallGeometry`, `:52`, `:98`). `tests/package/public_surface.txt` gains one line (43 → 44 names), as items 330 / 373 / 382 did. That snapshot is asserted by `test_wf14_validation_envelope.py:163` (S) and `:170` (every public name has a registry entry) and by the installed rehearsal (`tests/package/test_wf13_exit_rehearsal.py:43`), so RG7 re-runs both.
- **Never:** edit an existing line of `study.py` (the four guards are new blocks); change `ThermalBallGeometry`, `UniformHeatSource`, `TransientSpec`, `DutyCycle` or `RobinCooling`; admit a `ContainerWall` on a box; give a ball Study a steady branch (the ball's `transient` refusal at `:1446–1447` and U58's `thermal.transient`-removed row stay pinned).

**Test module: new `tests/package/test_wf15_study_thermal_box.py`.**
- **Docstring.** The first paragraph backticks only `WF-15` and `MAT-21` (the (F3) rule). It says the solve legs are gated in the complex build, the build of the records.
- **Imports (module-level names only; never a `test_*` name):**

  | Name | Source |
  |---|---|
  | `A_BOX`, `B_BOX`, `C_BOX`, `HALF`, `Q`, `T_C` | `tests/validation/test_mat21_box_closed_form.py` (`HALF` `:35`, `T_C` `:36`) |
  | `RUNGS`, `BASE`, `T_STAR`, `_steady_ex`, `_trans_ex`, `_exact`, `_interior_chi`, `_rung`, `_Clock` | `tests/validation/test_mat21_box_fem.py:70–162` |
  | `M_RUNG`, `_chi_box`, `_cn` | `tests/validation/test_mat21_box_duty_fem.py:83–101` |
  | `DT_BOX`, `HALF_T`, `N_STEPS`, `_duty_cf` | `tests/validation/test_mat21_box_duty_closed_form.py:54–56`, `:71` |
  | `K_GEL`, `H_GEL`, `T_INF`, `N_MAIN` | `tests/validation/test_mat15_coil_driven_heating_phantom.py:82–88` |
  | `RHO_C`, `E_BAND` | `tests/validation/test_mat12_transient_sar_heating_sphere.py:68`, `:76` |
  | `INTERIOR_BAND` | `tests/validation/test_mat13_transient_perfusion_sphere.py:64` (`MAT-21`'s own import) |
  | `S2_GAUSS_BAND`, `_rel`, `_max_imag` | `tests/validation/test_mat15_finite_cylinder_fem.py:69`, `:125`, `:135` |
  | `ROBIN_TAG`, `_allsum` | `tests/validation/test_mat10_sar_to_heat_sphere.py:56` |
  | `L0_BAND` | `tests/validation/test_mat5_em_coupled_sphere.py:104` |
  | `_same_mesh`, `_on`, `_theta_owned_rel` | `tests/package/test_wf15_study_transient.py:117`, `:127`, `:145` |
  | `IS_COMPLEX` | `tests/complex_mode.py:23` |

- **Record constants (transcribed with log lines):**
  - `STEP1_M4_WHOLE = "3.4599e-03"` (`20261013T142646Z_MAT-21-step1r-W1-steady.log:112`; ten digits `3.459853e-03` at `:102`);
  - `STEP1_M4_CELLS = 79872` and `STEP1_M4_P1_DOFS = 15741` (`:102`);
  - `P26_T_DUTY = "9.4014e-04"` (`20261014T003506Z_MAT-21-step2-a2-W3-closed-form.log:92`, `9.401399e-04`);
  - `P26_C3 = 5.0282e-01` (`:93`);
  - `STEP1_M6_CELLS = 269568` (step-1 log `:106`).
- **Nodes:**
  - `test_wf15_box_spec_and_refusals` (no mesh; `ThermalBoxGeometry.build` patched to raise; 0 build calls asserted);
  - `test_wf15_box_steady_reproduces_mat21_step1_m4` (W1);
  - `test_wf15_box_duty_reproduces_mat21_p26` (W2).
- **Opt-in:** `FEM_EM_WF15_S5=1` (item 484's `_knob()` pattern). When it is set, assert `IS_COMPLEX` and `MPI.COMM_WORLD.size == 2`; never skip on a mismatch.
- **Registry (lines added, each block under `# WF-15 step 5 (§9 item 521a)`):** a new `[names.ThermalBoxGeometry]` block (all three nodes), and the two solve nodes under `[names.UniformHeatSource]` and `[names.RobinCooling]`; the duty node under `[names.TransientSpec]` and `[names.DutyCycle]` (`tests/package/validation_envelope.toml:223–259`). **The §9 standing rule applies** (`PROJECT_PLAN.md:4058`): Wgen / Wval / Wapi / Wdoc in the same slot, and the regenerated `docs/guide/api/` and `docs/guide/validation.md` are committed.

**Anchors (asserted; bands imported).**
- **(K0) Spec and refusals, no mesh:**
  - `Study(ThermalBoxGeometry(HALF, (52, 32, 8)), None, UniformHeatSource(Q))` has `material is None` and `frequency_hz is None`;
  - `material=` anything, or an excitation other than `UniformHeatSource` ⇒ `TypeError`;
  - `.solve()`, `.mean_sar()` and `.estimate()` ⇒ `TypeError` (the last with the box message), 0 build calls;
  - `temperature(cooling, duty=DutyCycle(...))` without a transient ⇒ `ValueError`;
  - `RobinCooling(..., source_scale=2.0)` ⇒ `ValueError`;
  - `half_widths_m` of length 2, a zero half-width, `divisions` containing 0 ⇒ `ValueError`;
  - the ball's estimate message is still exactly `THERMAL_MESSAGE` (imported from `tests/package/test_wf13_cli_estimate.py` as a module constant, not a `test_*` name).
- **(UT) Bitwise up to the solve (W1 and W2).** The Study's mesh against `_rung(m, clk)`'s mesh: `_same_mesh` geometry equal, global cell count equal to the record (79 872 at m = 4, 269 568 at m = 6), and the facet-tag `indices` / `values` equal per rank (same partition, same call).
- **(I-dir) θ_Study = θ_direct ≤ `L0_BAND`** through `_theta_owned_rel(_on(direct_mesh, study_T), direct_T, comm)`: steady against `_rung`'s `st.temperature`; duty against `_cn(msh, ft, _chi_box).temperature`. Predicted 0 to ≈ 5e-14 (the same calls with the same arguments; item 500's file-vs-API read 8.1e-15 to 4.6e-14). Hex equality of the owned arrays is printed, never asserted.
- **(S1-m4) Steady, m = 4, P1:**
  - whole-box `_rel(msh, θ, _exact(msh, _steady_ex()))`, printed `.4e`, `== STEP1_M4_WHOLE`;
  - interior `_rel(…, _interior_chi(msh))` ≤ **`INTERIOR_BAND`** (predicted ≈ 2e-3 to 4e-3; m = 4 interior was not printed by step 1, so this is a new reading at the step's band, labelled as such);
  - P1 dofs `== 15741`;
  - Gauss `|∮ h θ ds − Q V_h| / (Q V_h)` ≤ **`S2_GAUSS_BAND`** (step 1's m = 4 reading 2.708e-13, `:102`);
  - θ_h(0) printed beside the series 12.003947 K (`:108`).
- **(S0) Duty schedule:** the Study's `source_multipliers` `array_equal` the direct `_cn(…, _chi_box)` run's. `DutyCycle(period_s=2.0 * HALF_T, duty=0.5)` (drafting arithmetic: `HALF_T` = 282.15, period 564.3) reproduces `_chi_box`'s χ at every evaluated time.
- **(T-duty) P26 through the Study, m = 6, P1, CN N = 40 with two BE half-steps** (`TransientSpec(DT_BOX, N_STEPS, RHO_C)`; defaults θ = ½, `startup_half_steps` 2, as `_cn` uses):
  - interior `_rel` against `_exact(msh, lambda x, y, z: _duty_cf(x, y, z, T_STAR))` ≤ **`INTERIOR_BAND`**, and printed `.4e` `== P26_T_DUTY`;
  - `max_balance_residual` ≤ **`E_BAND`** (P26 2.809e-12, `:92`).

**Negative controls (rule (e)).**
- **(c3) Unscheduled comparand on the duty θ — asserted, backed.** `_rel(msh, θ_Study, _exact(msh, _trans_ex(T_STAR)), χ) > INTERIOR_BAND`. It is the same comparison on the same fixture as P26's (c3), 5.0282e-01 (`…W3-closed-form.log:93`), tied to it by (I-dir). Printed beside `P26_C3`. One-sided; the miss sits near 0.5, far from any saturation.
- **(c-wall) The wrong wall through the API — *predicted*.** `RobinCooling(K_GEL, 2.0 * H_GEL, T_INF)` on the m = 4 Study, interior `_rel` against the comparand at the right h. Asserted only `> INTERIOR_BAND`. Predicted ≈ 0.37 (range 0.32–0.42): step 1's (c1) is the mirror comparison (FEM at h against the series at 2h), 5.9734e-01 at m = 6 (`…step1r-W1-steady.log:110`); renormalising by ‖θ(h)‖ instead of ‖θ(2h)‖ gives ≈ 0.597 / 1.6. No prior measurement exists of this exact comparison, so the factor is printed, never asserted. One extra 0.6 s solve.
- **(c0) By construction:** an `__init__` blind to the box branch would raise at `frequency_hz` (`:133–134`); (K0) asserts the box Study constructs.

**Windows.** Foreground; Bash `timeout: 660000`; pytest `-s`; node IDs, never `-k`; `tests/environment/test_complex_mode.py` first. **The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped.** `pgrep` = `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reads 0 before W2, RG1 and RG3. Run W0 only after every edit (rule (i)). Order: W0, W1, `pgrep`, W2, `pgrep`, RG1, RG2, `pgrep`, RG3, RG4, RG5, RG6, RG7, Wgen, Wval, Wapi, Wdoc.
- **W0** (complex, no mesh; the new spec node plus `WF-13` step 9's estimate module, since the estimate guard gained a block):
  `scripts/testing/run_and_log.sh WF-15-s5-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && export FEM_EM_REQUIRE_COMPLEX=1 && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf15_study_thermal_box.py::test_wf15_box_spec_and_refusals tests/package/test_wf13_cli_estimate.py -v -s --tb=short'"`
  Price: item 507's W1 on the estimate module, 8 passed in 31 s (`20261017T093117Z_WF-13-s9-W1.log`).
- **W1** (standard; predicted 25–45 s):
  `scripts/testing/run_and_log.sh WF-15-s5-W1-steady "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_WF15_S5=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf15_study_thermal_box.py::test_wf15_box_steady_reproduces_mat21_step1_m4 -v -s --tb=short'"`
- **W2** (heavy; predicted 210–290 s: P26's W3 fixture 177.8 s for one CN, the ten-run interpolation and the (c3) interpolation; this node adds one CN, ≈ 10 s; start-up ≈ 15 s with only `test_complex_mode.py`):
  `scripts/testing/run_and_log.sh WF-15-s5-W2-duty "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_WF15_S5=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 590 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf15_study_thermal_box.py::test_wf15_box_duty_reproduces_mat21_p26 -v -s --tb=short'"`
- **RG1 / RG2** (real build, rule (c): the constructor and dispatch insertions sit ahead of the ball's path). Item 330's W1a and W1b verbatim (`20261008T213834Z_WF-15-s3-W1a-duty.log.gz:12`, 124 s; `20261008T214046Z_WF-15-s3-W1b-wall.log.gz:12`, 133 s), labels `WF-15-s5-RG1-ball-duty` / `WF-15-s5-RG2-ball-wall`.
- **RG3 / RG4** (complex, rule (c): the constructor block runs ahead of every EM Study). Item 330's W0a (`20261008T213439Z_WF-15-s3-W0a-wf15.log.gz:12`, 119 s) and W0b (`20261008T213645Z_WF-15-s3-W0b-temp.log.gz:12`, 104 s) verbatim, labels `WF-15-s5-RG3-wf15` / `WF-15-s5-RG4-temp`.
- **RG5** (complex; `Study.probe` calls the edited dispatch). Item 326's W1 verbatim (`20261008T201316Z_WF-15-s4-W1-n2.log.gz:12`, 22 s), label `WF-15-s5-RG5-probe`. It rewrites `logs/wf15/probe-n2.npz` with the new `study.py` hash, which is that module's own guard.
- **RG6** (real then complex installed; the `__init__` surface and the snapshot hash). Item 432's rehearsal W1 verbatim (`20261014T021300Z_WF-13-s5-W1-rehearsal.log:12`, 28 s), label `WF-15-s5-RG6-rehearsal`.
- **RG7** (real + complex installed). Item 330's W4 verbatim (`20261008T214317Z_WF-15-s3-W4-wf13.log.gz:12`, 36 s), label `WF-15-s5-RG7-wf13`.
- **Wgen / Wval / Wapi / Wdoc:** item 500's four commands verbatim (`20261017T004807Z_WF-1-s14-Wgen.log:12`, `…004810Z…Wval.log:12`, `…004818Z…Wapi.log:12`, `…004837Z…Wdoc.log:12`), labels `WF-15-s5-Wgen` / `-Wval` / `-Wapi` / `-Wdoc`. Wapi carries (S) and the every-name check. Wdoc reads `0 dead`.
- **Not re-run (disclosed):** item 500's YAML legs (`run_study_file`'s dispatch is not touched in 521a); the `-n 8` F-human windows (no F-human path is touched; the constructor block is an `isinstance` that is false for them, and RG3 covers the constructor on the EM path).
- **`git diff --stat main -- src/`** shows `study.py` and `__init__.py`, insertions only.

**Traps already paid for.**
- **The comparand is evaluated in chunks at dof coordinates** (`_exact`'s `INTERP_CHUNK`, `test_mat21_box_fem.py:75`, `:94–108`), never `interpolate(callable)` on the series: each `MAT-21` leg paid one 180 s kill for it (`20261003T170330Z_MAT-17-step2-W1.log`; item 500's traps).
- **UFL cannot mix meshes.** Compare across runs only through `_on` (which asserts equal P1 dof coordinates) or in NumPy.
- **`_study_transient_temperature` builds on `MPI.COMM_WORLD`** (`study.py:1453`). The box's steady branch does the same; a sub-communicator is not claimed.
- **Floats:** the period is `2.0 * HALF_T`, computed from the imported constant and never retyped.
- **The thermal degree is `RobinCooling.degree`** (default 1), never `Study.degree` (item 260's trap).
- **Complex build:** take `np.real` after `_allsum`; no UFL ordering comparison.
- **Orphan `.c` after a kill:** run `pgrep` first; at 0, run the rubric sweep as `WF-15-s5-cache-sweep` in the killed build. Never raise a timeout or change `-n`.
- **Rule (i):** an edit to `study.py` after W0 re-runs W0, RG1–RG7.
- Write ≈, never a tilde.

**Scope.** Covered: a structured `create_box` thermal-only box Study with Robin cooling on all six faces under a uniform source; steady P1 at m = 4 and CN + Rannacher P1 at m = 6 with `MAT-21` step 2's square-wave duty; complex build, `-n 2`. **Not claimed:** the YAML (521b); P2 (`MAT-21` steps 3–5); non-uniform sources (U64); a container wall, perfusion or EM coupling on the box; `Study.probe` on the box (U60); `femem run --estimate` on a box (refused by A2); a gmsh or CAD box; the box inside a coil; any compliance statement.

**Negative result.** Report and stop. Park on `attempt/WF-15-step5-<UTC>`, mark item 521a BLOCKED with its unblock condition in the same commit (rule (d)), and journal it in attempts.md.
- **(K0) red:** a guard or dispatch defect. Fix in-slot before W1 only.
- **(UT) red:** the Study's box is not `_rung`'s box. Print the first differing coordinate and stop.
- **(I-dir) red with (UT) green:** the Study does not call the solver as the direct route does. Print both kwargs and stop.
- **(S1-m4) or (T-duty) `.4e` red with its band green:** an environment drift since 2026-10-14. This goes to `record-reconciler` (1*), never re-registered in-slot.
- **A band red:** `MAT-21` no longer reproduces. Open a known-issues entry; the review re-reads step 1 / P26.
- **Any RG red:** an insertion changed more than the box path. Revert `src/` and stop.
- **(c3) inside the band:** the duty schedule does not reach the solver. Stop.
- **Status 124:** not a result; park with the measured time.

No `PREDICTED RED:` line: a red anywhere stops the slot.

**Status it can move.** `WF-15` stays ✅ (its Done-when is step 1's). On green the §7 row gains: "**Step-5 record (§9 item 521a, U59 part 1, provisional):** `ThermalBoxGeometry` (`create_box`, Robin on all six faces) drives `Study.temperature` steady and CN with a duty cycle, complex, `-n 2`: `MAT-21` step 1's m = 4 whole-box L2 `3.4599e-03` and P26's (T-duty) `9.4014e-04` reproduce through the API (`INTERIOR_BAND` / `E_BAND` / `S2_GAUSS_BAND`), θ equal to the direct call at ≤ <x> (`L0_BAND`), multipliers bitwise; (c3) misses by <y>, the wrong wall by <z>." §10 U59 gains "part 1 DONE". `MAT-21` does not move.

**Spawn line:** "Use the implementer agent to execute `WF-15` §9 item 521a (U59 part 1) per `docs/planning/chunks/WF-1.md` "§9 item 521" (read the brief and any review ruling beneath it first; the rulings hold). The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Run W0, W1, W2, RG1–RG7, Wgen, Wval, Wapi, Wdoc in order, in the foreground, every harness Bash call with `timeout: 660000`, pytest with `-s`, container `timeout -k 30`. Before W2, RG1 and RG3, run `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` and require 0. W0–W2 run in the complex build with `FEM_EM_WF15_S5=1`; every re-gate runs in its own recorded build. `src` edits: `study.py` and `__init__.py`, insertions only, disclosed in the commit; never edit an existing line, a pinned message, or an existing test assertion. Evaluate comparands in chunks at dof coordinates, never `interpolate(callable)`; never import a `test_*` name; never compare across builds by `float.hex`. You are the executor: do not spawn agents, and never background a harness window. On a red, park on `attempt/WF-15-step5-*`, mark item 521a BLOCKED in the same commit, and report."

---

**Outcome (part 1 = 521a, 2026-10-18 11:38 UTC): DONE, green, attempt 1 of 4.** (K0) 9 passed; (UT) 79 872 / 269 568 cells, geometry and facet tags equal; (I-dir) 5.9e-16 steady / 9.1e-14 duty ≤ `L0_BAND`; (S1-m4) whole `3.4599e-03` = record, interior 3.2236e-03 (new reading), P1 dofs 15741, Gauss 2.708e-13; (S0) 41 multipliers equal; (T-duty) `9.4014e-04` = record, balance 2.809e-12; (c3) 5.0282e-01 = P26's; (c-wall, predicted 0.32–0.42) 3.8048e-01. RG1–RG7 and Wgen / Wval / Wapi / Wdoc green. Disclosed: the `pgrep` check was run before W2 only, not before RG1 / RG3 (each started after its predecessor's footer returned). Logs `20261018T112129Z_WF-15-s5-W0.log` … `20261018T113759Z_WF-15-s5-Wdoc.log`.

---

#### 521b — `WF-1` step 15: the `thermal_box` file through `femem run`

- **Executor:** implementer. **Tier:** W0 / RG1 / RG3 / RG5 standard; W1a standard (complex; the m = 4 steady file); W1b **heavy** (complex; the m = 6 duty file); RG4 heavy (`-n 8`, complex); RG6 standard (real). `-n 2` everywhere else.
- **Class:** chain step (solve). **Cost:** ≈ 18 measured / ≈ 55 predicted min = 15 fixed + ≈ 18 code (≈ 120 `study_yaml.py` lines, ≈ 6 `cli.py` lines, module ≈ 300 lines, three good files and one malformed) + windows ≈ 16 (W0 0.1, RG1 0.1, W9 0.3, RG2 1.1, RG3 0.6, RG5 0.7, W1a 0.7, W1b 4.8, RG4 3.6, RG6a 1.1, RG6b 1.2, W-gen/val/api/doc 0.4) + ≈ 6 journal. Measured precedent: item 500's windows (`20261017T003617Z` … `004837Z_WF-1-s14-*.log`; RG4 213 s). **Park if code + W0 / RG1 / W9 are not green by minute 30, W1b has not started by minute 38, or RG4 cannot start by minute 48.**
- **Depends on:** item 521a (DONE). **Family:** `WF-1` step 15, attempt 1 of 4.

**The change (insertions only, rule (c), the same seams item 500 used).**
- **(B1)** `study_yaml.py`'s `from .study import (…)` block gains `ThermalBoxGeometry`.
- **(B2)** In `_parse_mapping`, directly **after** the `thermal_ball` hook (`:1014–1016`): `if isinstance(geo, dict) and geo.get("kind") == "thermal_box": return _parse_thermal_box_mapping(mapping, source)`.
- **(B3)** New `_parse_thermal_box_mapping(mapping, source)`, returning `(Study, "thermal_box", outputs, _BallThermal | _BoxSteady)`. It follows `_parse_thermal_ball_mapping`'s order and refusals (`:740–745`, with the reasons re-worded to `thermal_box`):
  - `geometry`: `_unknown` (`kind`, `half_widths_m`, `divisions`); `_missing` both; `half_widths_m` three finite > 0 numbers; `divisions` three integers ≥ 1. `wall` is refused by `_unknown`.
  - `excitation`: item 500's lines (`:770–775`) verbatim.
  - `thermal`: `_unknown` (`k_w_per_m_k`, `h_w_per_m2_k`, `t_inf_c`, `degree`, `source_scale`, `transient`, `duty`); `source_scale` refused as at `:780–784`. **`transient` is optional here** (a steady box); `duty` without `transient` ⇒ at `thermal.duty` ("needs `thermal.transient`"). The base keys go through `_parse_thermal(base, "sphere_in_box", source)` (`:788`); the transient and duty blocks are **copied** from `:792–840` (never factored: a shared helper would edit the ball's path).
  - `outputs`: `directory` only.
  - The Study is `Study(ThermalBoxGeometry(half, divisions), None, UniformHeatSource(q))`.
  - The fourth element is `_BallThermal(cooling, transient, duty)` (`:714–720`) when a transient is present, else `_BallThermal(cooling, None, None)`. That reuses the frozen class unedited, since its `transient` field is typed `Any`.
- **(B4)** In `run_study_file`, **before** the ball dispatch (`:1380`): `if isinstance(study.geometry, ThermalBoxGeometry): return _run_thermal_box(...)`.
- **(B5)** New `_run_thermal_box`: `study.temperature(cooling.cooling, transient=cooling.transient, duty=cooling.duty, comm=comm)`. The summary carries `family "thermal_box"`, `cells`, and `temperature` = {`half_widths_m`, `divisions`, `thermal_degree`, `q_w_per_m3`, `transient` (null or {`dt_s`, `n_steps`, `theta`, `startup_half_steps`, `rho_c_j_per_m3_k`}), `duty` (null or {`period_s`, `duty`}), `source_power_w`, `convected_power_w` (steady only), `source_multipliers` (transient only), `max_balance_residual` (transient only), `theta_max_k`} (the `theta_max_k` lines copied from `:1306–1309`, `:1338`). Rank 0 writes; `Barrier`. No `degree` / `frequency_hz` / Touchstone.
- **(B6, `cli.py`)** In `check`, before `:158`: `if kind == "thermal_box":`, then rank 0 prints `family: thermal_box` and `excitation: {type(study.excitation).__name__}` (the ball branch's own line, `cli.py:161`, so U64's new source class needs no further `cli.py` edit) and returns 0 (≈ 6 lines; `test_wf1_check_family.py::test_wf1_check_other_lines_unchanged` pins every other family's lines).
- **Docstrings:** a `thermal_box` bullet in the module docstring; one sentence in `run_study_file`'s.
- **Never:** edit an existing line of either file; widen `_load`'s 4-tuple; touch `_parse_thermal_ball_mapping` or `_run_thermal_ball`; admit `wall` / `t_init_c` / `source_scale`; add a name to `study_yaml`'s `__all__`.

**Tracked inputs (four new files under `tests/package/wf1/`; each opens `# WF-1 step 15 (§9 item 521b): …`, then `femem: 1`; floats as Python `repr` digits, never `4.18e6`).**
- `thermal_box_steady_m4.yaml`: `geometry: {kind: thermal_box, half_widths_m: [0.325, 0.21, 0.045], divisions: [52, 32, 8]}`; `excitation: {kind: uniform_heat_source, q_w_per_m3: 2000.0}`; `thermal: {k_w_per_m_k: 0.6, h_w_per_m2_k: 10.0, t_inf_c: 20.0}`.
- `thermal_box_duty_m6.yaml`: `divisions: [78, 48, 12]`; `thermal` as above plus `transient: {dt_s: 70.5375, n_steps: 40, rho_c_j_per_m3_k: 4180000.0}` and `duty: {period_s: 564.3, duty: 0.5}`. These are drafting reprs of `DT_BOX` and `2.0 * HALF_T` (computed from `T_C = C_BOX**2 * RHO_C / K`, `test_mat21_box_closed_form.py:36`); (K) asserts them by `float.hex` against the imported constants, so a mis-rounded value fails at parse.
- `thermal_box_wrong_wall_m4.yaml`: the steady file with `h_w_per_m2_k: 20.0` (the (c-wall) control as a tracked file).
- `malformed_thermal_box_duty_without_transient.yaml` → `thermal.duty`.

**Test module: new `tests/package/test_wf1_yaml_thermal_box.py`.** It is item 500's module shape (`tests/package/test_wf1_yaml_thermal_ball.py:113–207`; `_knob`, `_write_yaml`, `_hex_tree_equal`). The knob is `FEM_EM_WF1_S15=1`; when it is set, assert `IS_COMPLEX` and size 2. It imports 521a's record constants by name from `tests/package/test_wf15_study_thermal_box.py` (module constants, not `test_*` names), and item 500's `_spec_diff` / `_cmp` / `_rc_all` / `_fresh_dir` (`tests/package/test_wf1_study_yaml.py:120–169`). Nodes: `test_wf1_s15_parse_and_refusals`, `test_wf1_s15_steady_reproduces_mat21_m4`, `test_wf1_s15_duty_reproduces_p26`. Registry: the parse node under `[names.load_study]` / `[names.StudyConfigError]`; the solve nodes under `[names.run_study_file]` and `[names.ThermalBoxGeometry]`, after item 500's blocks. **The §9 standing rule applies** (Wgen / Wval / Wapi / Wdoc, pages committed).

**Anchors (asserted).**
- **(K) Parse, no mesh** (every `build` patched to raise; 0 calls):
  - `_spec_diff(load_study(file), Study(ThermalBoxGeometry(HALF, divs), None, UniformHeatSource(Q))) == []` on the three good files;
  - the fourth element through `_cmp` by hex: `RobinCooling(K_GEL, H_GEL, T_INF)` (wrong wall: `2.0 * H_GEL`); duty file `TransientSpec(DT_BOX, N_STEPS, RHO_C)` and `DutyCycle(2.0 * HALF_T, 0.5)`; steady `transient is None`;
  - `cli.main(["check", f])` rc (0, 0) via `_rc_all`, with exactly one `family: thermal_box` line on rank 0.
- **(X) Refusals naming their field** (`StudyConfigError.field`): the malformed file (also through `run` and `check`: rc (2, 2), `out` absent by `MPI.LOR`); `study` / `material` / `conductor` / `tune` added; `geometry.half_widths_m` of length 2 / a zero entry / a string; `geometry.divisions` with 0 / a float; `geometry.wall` added; `geometry.radius` added; `excitation.q_w_per_m3: 0.0`; `thermal.source_scale: 2.0`; `thermal.transient.n_steps: 0`; `thermal.duty.duty: 1.5`; `outputs.probes` added. **Pinned refusals that must stay green unedited (W0 / RG1):** item 500's whole (X) table, including `thermal.transient` removed → `thermal` on a **ball** file (the box admits a missing transient; the ball still refuses one).
- **(UT) / (I-Py)** on each solve node: the file's mesh equals 521a's Python Study's (`_same_mesh`; cells 79 872 / 269 568); multipliers `array_equal`; θ_file vs θ_Python `_theta_owned_rel` ≤ **`L0_BAND`** (U59's clause), predicted 0 to ≈ 5e-14 (hex printed).
- **(R) The records through the file:** steady whole-box `.4e` `== "3.4599e-03"` and interior ≤ **`INTERIOR_BAND`**; Gauss ≤ `S2_GAUSS_BAND`; duty (T-duty) `.4e` `== "9.4014e-04"` and ≤ **`INTERIOR_BAND`**; `summary["temperature"]["max_balance_residual"]` ≤ **`E_BAND`**.
- **(I=)** The summary equals its own result by hex (item 500's (I=), with `elapsed_s` compared on rank 0 only, its measured fix at `…004140Z_WF-1-s14-W1a-duty.log`); `family == "thermal_box"`; no `degree` / `frequency_hz` / `touchstone` key; `run.result is None`.

**Negative controls (rule (e)).**
- **(c-wall) The wrong wall through a tracked file — *predicted*:** `thermal_box_wrong_wall_m4.yaml`'s interior miss against the m = 4 comparand at the right h is `> INTERIOR_BAND`. It is backed by 521a's (c-wall) on the same fixture **once 521a has run** (cite its log line in the module before W1a; it is then *asserted, backed*, tied by (I-Py)).
- **(c3) through the duty file — *asserted, backed*** by P26's 5.0282e-01 and 521a's reading.
- **(c0) By construction:** (K)'s fourth element of the wrong-wall file differs from the steady file's in exactly `h_w_per_m2_k`, and `_spec_diff` between them is `[]`. The wall rides in the cooling, not the Study.

**Windows** (the harness call is the whole Bash command; foreground; `timeout: 660000`; `-s`; node IDs; `pgrep` 0 before RG2, RG3, RG5, W1b, RG4 and RG6a). Order: W0, RG1, W9, RG2, RG3, RG5, W1a, W1b, RG4, RG6a, RG6b, Wgen, Wval, Wapi, Wdoc.
- **W0** (complex, no mesh, rule (c): B2 sits on every file's parse path, B4 on every run's). These are **the YAML gate modules**: item 500's W0 plus the new node:
  `scripts/testing/run_and_log.sh WF-1-s15-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_thermal_box.py::test_wf1_s15_parse_and_refusals tests/package/test_wf1_yaml_thermal_ball.py::test_wf1_s14_parse_and_refusals tests/package/test_wf1_yaml_capacitors.py::test_wf1_s13_parse_and_refusals tests/package/test_wf1_yaml_tune.py::test_wf1_s11_tune_parse_and_refusals tests/package/test_wf1_yaml_tune_f_human.py::test_wf1_s12_parse_and_refusals tests/package/test_wf1_check_family.py tests/package/test_wf1_yaml_f_human.py tests/package/test_wf1_yaml_f_human_thermal.py -v -s --tb=short'"`
  Price: item 500's W0, 14 passed in 4 s (`20261017T004126Z_WF-1-s14-W0.log`). Check-family's glob grows from 18 to 21 good files.
- **RG1** (every `[names.load_study]` parse node not in W0): item 500's RG1 verbatim (`20261017T003637Z_WF-1-s14-RG1.log:12`, 4 s), label `WF-1-s15-RG1`.
- **W9** (B6 edits the installed `femem`): item 500's W9 verbatim (`20261017T003641Z_WF-1-s14-W9-installed-femem.log:12`, 16 s), label `WF-1-s15-W9-installed-femem`.
- **RG2 / RG3 / RG5** (the four-port, sphere and capacitor branches behind B4's head): item 500's commands verbatim (`20261017T003710Z_WF-1-s14-RG2.log:12`, 64 s; `20261017T003824Z_WF-1-s14-RG3.log:12`, 33 s; `20261017T003905Z_WF-1-s14-RG5.log:12`, 42 s), labels `WF-1-s15-RG2` / `-RG3` / `-RG5`.
- **W1a** (complex, standard; predicted 25–50 s):
  `scripts/testing/run_and_log.sh WF-1-s15-W1a-steady "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_WF1_S15=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_thermal_box.py::test_wf1_s15_steady_reproduces_mat21_m4 -s -v --tb=short'"`
- **W1b** (complex, heavy; predicted 230–310 s: 521a's W2 plus one file run, ≈ 10 s):
  `scripts/testing/run_and_log.sh WF-1-s15-W1b-duty "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_WF1_S15=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 590 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_thermal_box.py::test_wf1_s15_duty_reproduces_p26 -s -v --tb=short'"`
- **RG4 — required** (B4's head runs ahead of every F-human file's `study.solve`; the 10-16 rulings on items 496 / 500 required it): item 500's RG4 verbatim (`20261017T004425Z_WF-1-s14-RG4.log:12`, 213 s), label `WF-1-s15-RG4`, alone.
- **RG6a / RG6b** (real build, knob `FEM_EM_WF1_S14=1`; B2 and B4 sit directly beside the ball's hook and dispatch): item 507's RG2 / RG3 verbatim (`20261017T093213Z_WF-13-s9-RG2-thermal-duty.log:12`, 108 s; `20261017T093406Z_WF-13-s9-RG3-thermal-wall.log:12`, 111 s), labels `WF-1-s15-RG6a-ball-duty` / `-RG6b-ball-wall`.
- **Wgen / Wval / Wapi / Wdoc:** item 500's four verbatim, labels `WF-1-s15-W*`. The Wval node count is the last Wval's plus 3; print it.
- **`git diff --stat main -- src/`** shows `study_yaml.py` and `cli.py`, insertions only; name B1–B6 in the commit.

**Traps already paid for:** item 500's list (this file, item 500 "Traps") applies whole. In particular: YAML 1.1 loads `4.18e6` as a string; never `float.hex` across builds; `femem` prints on rank 0 only; `-k a or b` splits argv; `_load`'s 4-tuple is load-bearing; `__all__` does not grow. **New here:**
- **The box and ball hooks are adjacent.** B2 must come after the ball's `return` and match only `"thermal_box"`, so a ball file never reaches it. W0's ball parse node is the check.
- **`divisions` must reach `create_box` as Python `int`s.** YAML gives `int`; `_integer` (`study_yaml.py:191`) returns `int`; never `float`.
- **The steady file carries no `transient`.** `_run_thermal_box` must not read `transient.dt_s` on it; the summary writes `null`.

**Scope.** Covered: `femem run` / `femem check` on three tracked `thermal_box` files (uniform source, Robin on six faces, steady at m = 4, duty CN at m = 6, P1, complex, `-n 2`), reproducing 521a's records and equal to the Python Study at `L0_BAND`. **Not claimed:** P2 / `thermal.degree: 2` records (U64 carries P2); the transverse source (U64); probes on the box (U60); the installed console script running the new files (W9 re-gates the installed `check` on the old ones); `--estimate` (refused, rc 2, by 521a's guard and item 480's handler — not asserted here); the real build; a container wall; any compliance statement.

**Negative result.** Item 500's dispositions (this file), with these additions:
- **(I-Py) red with (K) / (UT) green:** two equal specs on one build give θ beyond `L0_BAND`. That is a solve-floor finding, band unmoved.
- **RG6 red:** the box hook leaked into the ball route. Revert and stop.

Park on `attempt/WF-1-step15-<UTC>`, item 521b BLOCKED in the same commit.

No `PREDICTED RED:` line.

**Status it can move.** `WF-1` stays ✅ (rule (j): its Done-when is step 1's). On green the row gains: "**Step-15 record (§9 item 521b, U59 part 2, provisional):** a `thermal_box` geometry (`half_widths_m`, `divisions`), `uniform_heat_source` and optional `thermal.transient` / `thermal.duty` drive `Study.temperature` through `run_study_file`, complex, `-n 2`: `MAT-21` step 1's m = 4 `3.4599e-03` and P26's `9.4014e-04` reproduce (`INTERIOR_BAND` / `E_BAND`), θ equal to the Python Study at ≤ <d> (`L0_BAND`); the wrong-wall file misses by <p>; bad box keys raise `StudyConfigError` naming their field; `femem check` prints `family: thermal_box`; `src`: `study_yaml.py` + `cli.py`, insertions only." §10 U59 DONE (provisional). **For the weekly:** the `mri:` thermal-box example the 10-17 weekly named (`PROJECT_PLAN.md:4997`) becomes openable.

**Spawn line:** item 500's spawn line with `521b` / `step 15` / `attempt/WF-1-step15-*` / `FEM_EM_WF1_S15=1`, the window list W0, RG1, W9, RG2, RG3, RG5, W1a, W1b, RG4, RG6a, RG6b, Wgen, Wval, Wapi, Wdoc, **W1a / W1b in the complex build** (`source /usr/local/bin/dolfinx-complex-mode`, `FEM_EM_REQUIRE_COMPLEX=1`) and RG6a / RG6b in the real build with `FEM_EM_WF1_S14=1`.

**For the review.**
1. **Build.** Complex for the solve legs (the `MAT-21` records' build), against U59's "real build". The alternative is the real build with the records printed, not asserted.
2. **Tier.** The duty leg is heavy by measurement (P26 W3, 262 s), against U59's "standard".
3. **The split** (521a `WF-15` step 5 + 521b `WF-1` step 15, one chain step). Declining it makes one ≈ 100-minute item.
4. **`wall` as a `None` property on `ThermalBoxGeometry`**, so `_study_transient_temperature` runs unedited. The alternative is a one-line edit at `study.py:1455` (`getattr(geo, "wall", None)`), disclosed.
5. **"Wrong-wall" read as h doubled through the cooling.** The other reading, insulating a face pair, is a new boundary option that `RobinCooling` does not have; the drafter does not pose it.
6. **The export** (`public_surface.txt` 43 → 44). It turns the U47 / item 432 "43 names" clause into a dated record; RG6 / Wapi re-gate the new snapshot.

**Outcome (part 2 = 521b, 2026-10-18 18:47–19:02 UTC): DONE, green, attempt 1 of 4.** All 15 windows green in order, `pgrep` 0 before RG2, RG3, RG5, W1b, RG4 and RG6a. W0 15 passed (14 + the new parse node; check-family 21 good files); RG1 13 passed; W9 1 passed; RG2 10 passed (62 s); RG3 5 passed (33 s); RG5 10 passed (42 s). W1a (complex, 10 s): (UT) 79 872 cells, geometry equal; (I-Py) 1.479e-15 (`0x1.aa67490b1fa24p-50`; owned arrays not bitwise, printed only); whole L2 3.459853e-03 → `3.4599e-03` = record, interior 3.223613e-03 (= 521a), Gauss 2.708e-13; (c-wall) through the tracked file 3.8048e-01 = 521a's API reading (asserted > `INTERIOR_BAND`, backed). W1b (complex, heavy, 190 s): (UT) 269 568 cells, geometry equal, 41 multipliers equal; (I-Py) 4.792e-14 (`0x1.af9dd76419c4ap-45`); (T-duty) 9.401399e-04 → `9.4014e-04` = record; balance 2.809e-12; (c3) 5.0282e-01 = P26. (I=) both summaries equal their own results by hex. RG4 (`-n 8`) 9 passed (211 s); RG6a (real) 3 passed 3 skipped (63 s), RG6b 2 passed 3 skipped (75 s) — the precedent's tallies. Wgen / Wval / Wapi (22 passed) / Wdoc green. **Wval node count: 282 against the committed page's 276 (521a's Wval) — +6, not the brief's +3:** the stat counts (name, node) rows, and the three new nodes sit under six names (parse node under `load_study` + `StudyConfigError`; each solve node under `run_study_file` + `ThermalBoxGeometry`). (Item 523's Wval log reads 269; not investigated.) `src`: `study_yaml.py` +189, `cli.py` +6, insertions only (B1–B6); `__all__` unchanged. Logs `20261018T184753Z_WF-1-s15-W0.log` … `20261018T190227Z_WF-1-s15-Wdoc.log`.

### §9 item 528 — 2026-10-18 (§10 U64, the transverse source from YAML; chain step, solve)

*Drafted by a read-only agent for the 2026-10-18 03:00 daily review; the review reviews it. Code read on `d94320db`; nothing was run.*

**Licence.** §10 chain U, U64 (`PROJECT_PLAN.md:4861–4863`): "The `thermal_box` file with P36's source as `source: {kind: transverse_quadratic}` reproduces P36's record at its bands. Serial on U59." Chain U's target applies. **No band, tolerance, record, Done-when or target is set or moved.** P36 is `MAT-21` step 9 (§9 item 497, `docs/planning/chunks/MAT-21.md:1392–1666`, ratified by the 10-17 weekly). Templates: items 521a / 521b (above) and 500.

**Owner.** `WF-1` step 16 (the file), with a `WF-15` step 7 clause (the source class the file names). One item: the API addition is ≈ 40 lines and has no record of its own apart from the file's.

**Code state.** Nothing exists yet: `ThermalBoxGeometry`, `_study_box_temperature`, `_parse_thermal_box_mapping` and `_run_thermal_box` are item 521's (521a / 521b), not on `main` at drafting. The only imposed thermal source is `UniformHeatSource` (`src/fem_em_solver/study.py:1357–1364`), a scalar. `_study_transient_temperature` reads it as `float(self.excitation.q_w_per_m3)` (`:1454`), so a spatial source cannot pass through that function unedited. P36 itself passes a UFL expression straight to both solvers: `q_quad = Q * 3.0 * (xs[1]**2 + xs[2]**2) / (B_BOX**2 + C_BOX**2)` (`tests/validation/test_mat21_box_quadratic_fem.py:141–142`), into `solve_steady_conduction` (`:143–144`) and `solve_transient_conduction` (`:180–183`). The mean of that shape over the box is Q (∫3Q(y² + z²)/(B² + C²) dV / V = Q), so `source_power = Q V`; P36 read `source_power / (Q V_h) − 1 = −1.350e-13` (`20261016T171904Z_MAT-21-step9-W2-quadratic-fem.log:73`).

- **Executor:** implementer. **Tier:** standard throughout, `-n 2`, complex (P36's build; its W2 ran 126 s for one uniform steady, one quadratic steady, one CN and five P3 interpolations, `…quadratic-fem.log:115`). This item's solves are split over two windows so that each stays under 180 s with margin. RG-B is heavy (521a's W2), and the review may waive it (question 3).
- **Class:** chain step (solve). **Cost:** ≈ 18 measured / ≈ 50 predicted min = 15 fixed + ≈ 15 code + windows ≈ 12 (W0 0.1, RG1 0.1, W1a 1.5, W1b 2.0, RG-A 0.8, RG-B 4.5, RG-C 0.7, RG-D 0.5, Wgen / Wval / Wapi / Wdoc 0.4) + ≈ 6 journal. **Park if code + W0 are not green by minute 28, or W1b has not started by minute 38.**
- **Depends on:** item 521 (521a and 521b both DONE). If only 521a has landed, the item is not runnable; it stays OPEN with "waits on 521b". **Family:** `WF-1` step 16, attempt 1 of 4.

**The change (insertions only, rule (c), disclosed).**
- **(Q1, `study.py`)** A new frozen `TransverseQuadraticHeatSource(q_w_per_m3: float)` beside `UniformHeatSource`, with a docstring: "q(x) = 3 q (y² + z²)/(B² + C²) on a `ThermalBoxGeometry` with half-widths (A, B, C); its box mean is q (`MAT-21` step 9's transverse eddy-SAR shape)". `_positive` applies. It is exported (`__all__`, `__init__.py`, `public_surface.txt` 44 → 45).
- **(Q2)** In 521a's box constructor block, admit `(UniformHeatSource, TransverseQuadraticHeatSource)`. This is an **edit of 521a's line**, disclosed under rule (c), and re-gated by RG-A. A ball Study still refuses the new class (the ball's check at `:119` is untouched).
- **(Q3)** At the top of 521a's `_study_box_temperature`: `if isinstance(self.excitation, TransverseQuadraticHeatSource): return _study_box_quadratic_temperature(self, cooling, transient, duty)`. The uniform path below it is unchanged.
- **(Q4)** New `_study_box_quadratic_temperature`. It builds the box; `xs = ufl.SpatialCoordinate(msh)`; `q = q0 * 3.0 * (xs[1] * xs[1] + xs[2] * xs[2]) / (B**2 + C**2)`, with P36's operation order (`:142`) and `B`, `C` read from `geometry.half_widths_m`.
  - Steady: `solve_steady_conduction(msh, ft, k=…, q=q, h=…, t_inf=…, robin_tag=_ROBIN_TAG, degree=cooling.degree)`.
  - Transient: `solve_transient_conduction(...)` with `_study_transient_temperature`'s keyword set **copied** from `:1467–1483`, with `q=q`; `duty` passes through as `source_schedule` (no duty record is claimed).
  - It returns 521a's result types.
- **(Q5, `study_yaml.py`)** In 521b's `_parse_thermal_box_mapping`, `excitation` gains one optional key, `source`. It is a mapping with exactly `kind ∈ ("uniform", "transverse_quadratic")` (`_unknown` / `_missing` / `_kind` at `excitation.source`). It is absent or `uniform` ⇒ `UniformHeatSource(q)` (521b's files unchanged), and `transverse_quadratic` ⇒ `TransverseQuadraticHeatSource(q)`. `excitation.source` on a `thermal_ball` file is refused by item 500's `_unknown` (`study_yaml.py:773`), unchanged, and a row in (X) pins it.
- **(Q6)** 521b's `_run_thermal_box` summary gains `"source": "uniform" | "transverse_quadratic"` under `temperature`.
- `cli.py` needs no edit if 521b's B6 prints the excitation's type name (as its brief says). If it landed a literal, the item adds one disclosed line plus W9.
- **Never:** edit `_study_transient_temperature`, `UniformHeatSource`, the ball route, or any `MAT-21` module.

**Tracked inputs (`tests/package/wf1/`, repr floats).** Each is 521b's steady m = 4 file with `divisions: [52, 32, 8]` and `thermal.degree: 2` (P36 is P2, `tests/validation/test_mat21_box_fem_p2.py:77`):
- `thermal_box_quadratic_steady_m4.yaml`: `excitation: {kind: uniform_heat_source, q_w_per_m3: 2000.0, source: {kind: transverse_quadratic}}`;
- `thermal_box_quadratic_cn_m4.yaml`: the same plus `transient: {dt_s: 70.5375, n_steps: 40, rho_c_j_per_m3_k: 4180000.0}`. That is `T_STAR / N_MAIN`, P36's CN (`:180–183`); `T_STAR` is step 1's 0.2 t_C, and (K) asserts it by hex;
- `thermal_box_uniform_steady_m4_p2.yaml`: the steady file at degree 2 with no `source`, the control and (U).

**Test module: new `tests/package/test_wf1_yaml_thermal_box_quadratic.py`** (knob `FEM_EM_WF1_S16=1`; assert `IS_COMPLEX` and size 2 when it is set).
- Imports: `_quad_ex`, `_quad_trans_ex` (`tests/validation/test_mat21_box_quadratic_fem.py:85–98`); `DEGREE`, `CMP_DEGREE`, `QDEG_P`, `_exact_p`, `_rel_p` (`tests/validation/test_mat21_box_fem_p2.py:77–79`, `:91`, `:108`); `M4`, `STEP3_M4` (`tests/validation/test_mat21_box_separable_fem.py:57–58`); `_steady_ex`, `_trans_ex`, `_interior_chi`, `T_STAR`, `V_BOX` (`tests/validation/test_mat21_box_fem.py`); bands as 521a.
- Record constants transcribed from `20261016T171904Z_MAT-21-step9-W2-quadratic-fem.log`:
  - `P36_S1 = "7.219878e-05"` (`:76`);
  - `P36_T = "2.694340e-04"` (`:83`);
  - `P36_E = 2.285e-13` (`:84`);
  - `P36_S2 = 9.775e-14` (`:73`);
  - `P36_DOFS = 116025` (`:70`);
  - `P36_C_SRCP = 5.1471e-01` (`:78`).
- Nodes: `test_wf1_s16_parse_and_refusals`, `test_wf1_s16_quadratic_steady_reproduces_p36`, `test_wf1_s16_quadratic_cn_reproduces_p36`.
- Registry: a new `[names.TransverseQuadraticHeatSource]` block; the solve nodes under `[names.run_study_file]` and `[names.ThermalBoxGeometry]`. **The §9 standing rule applies** (`PROJECT_PLAN.md:4058`): Wgen / Wval / Wapi / Wdoc, pages committed.

**Anchors (asserted; bands imported by name).**
- **(K) Parse, no mesh:**
  - `_spec_diff(load_study(quad_file), Study(ThermalBoxGeometry(HALF, (52, 32, 8)), None, TransverseQuadraticHeatSource(Q))) == []`;
  - the uniform file's Study has a `UniformHeatSource`;
  - cooling `RobinCooling(K_GEL, H_GEL, T_INF, degree=2)` by hex;
  - `check` prints `excitation: TransverseQuadraticHeatSource`.
- **(X)** `excitation.source: {kind: parabolic}` ⇒ `excitation.source.kind`; `source: [1]` ⇒ `excitation.source`; `source` with an extra key ⇒ `excitation.source.<key>`; on a ball file ⇒ `excitation.source`. Pinned and unedited: item 500's and 521b's (X) tables.
- **(Pdeg)** The P2 dof count `== P36_DOFS` on every solve.
- **(I-Py)** θ_file = θ of the same Python Study ≤ **`L0_BAND`**, steady and CN (predicted 0 to ≈ 5e-14).
- **(S1-quad)** The steady interior `_rel_p(msh, θ, _exact_p(msh, _quad_ex()), _interior_chi(msh))` ≤ **`INTERIOR_BAND`** and `f"{e:.6e}" == P36_S1`. The `.6e` match is P36's own print precision; it is asserted, because the build, fixture and code path are P36's, as (U) shows.
- **(S2-quad)** `|h ∮ θ ds − source_power| / source_power` ≤ **`S2_GAUSS_BAND`**, printed beside `P36_S2`.
- **(T-quad)** The CN interior against `_quad_trans_ex(T_STAR)` ≤ **`INTERIOR_BAND`** and `.6e == P36_T`.
- **(E-quad)** `max_balance_residual` ≤ **`E_BAND`**, printed beside `P36_E`.
- **(U)** The uniform P2 file reproduces step 3's m = 4 digits: whole `2.489307e-05` and interior `1.669359e-05` (`STEP3_M4`), P36's own (U) check (`…quadratic-fem.log:70`).

**Negative controls (rule (e)).**
- **(c-src′) The source dropped through a file — *asserted, backed*.** The uniform P2 file's interior miss against the quadratic series is `> INTERIOR_BAND`, printed beside `P36_C_SRCP` 5.1471e-01. It is the same comparison on the same fixture: P36's (c-src′) is "uniform FEM vs quadratic series" (`:78`). One-sided.
- **(c-src) — *asserted, backed*.** The quadratic file's interior miss against the uniform series is `> INTERIOR_BAND`, beside 4.4856e-01 (`:77`).
- **(c-z) — printed only.** The quadratic file against the no-z series. P36 read 4.7893e-02 (`:79`), only 4.79× the band, so the drafter does not assert it.

**Windows** (the harness call is the whole Bash command; foreground; `timeout: 660000`; `-s`; node IDs; complex; `pgrep` 0 before W1a and RG-B). Order: W0, RG1, W1a, W1b, RG-A, RG-B, RG-C, RG-D, Wgen, Wval, Wapi, Wdoc.
- **W0:** 521b's W0 command with `tests/package/test_wf1_yaml_thermal_box_quadratic.py::test_wf1_s16_parse_and_refusals` added after the `test_complex_mode.py` path, label `WF-1-s16-W0` (these are the YAML gate modules, and check-family's glob grows by three).
- **RG1:** item 500's RG1 verbatim (`20261017T003637Z_WF-1-s14-RG1.log:12`), label `WF-1-s16-RG1`.
- **W1a** (steady; predicted 70–110 s: three P2 steady solves at ≈ 10 s, `…step3-W2-p2-transient.log:88`, and two P3 interpolations at ≈ 10 s, `…quadratic-fem.log`, step 9's clock lines):
  `scripts/testing/run_and_log.sh WF-1-s16-W1a-steady "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_WF1_S16=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_thermal_box_quadratic.py::test_wf1_s16_quadratic_steady_reproduces_p36 -s -v --tb=short'"`
- **W1b** (CN; predicted 90–130 s: two P2 CN N = 40 at 33.7 s each, `20261014T095800Z_MAT-21-step3-W2-p2-transient.log:90`, plus one P3 interpolation):
  `scripts/testing/run_and_log.sh WF-1-s16-W1b-cn "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_WF1_S16=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_thermal_box_quadratic.py::test_wf1_s16_quadratic_cn_reproduces_p36 -s -v --tb=short'"`
  If W1a measures over 120 s, run W1b with `timeout -k 30 590` and say so. Never shrink the case to fit.
- **RG-A** (rule (c): Q2 / Q3 edit 521a's code): 521a's W0 and W1 verbatim, labels `WF-1-s16-RG-A0` / `-RG-A1`.
- **RG-B** (521a's W2, heavy, 590): the duty route passes Q3's `isinstance`. **The review may waive it** citing RG-A1, which runs the same head on the uniform path (question 3).
- **RG-C** (521b's W1a verbatim): the uniform box file is unchanged with `source` absent.
- **RG-D** (the export): 521a's RG6 (item 432's rehearsal) verbatim.
- **Wgen / Wval / Wapi / Wdoc:** item 500's four verbatim, labels `WF-1-s16-W*`.
- **Not re-run (disclosed):** RG2–RG5 / RG4 of 521b. No line on any other family's path moves: Q5 and Q6 sit inside the box's own parser and runner.

**Traps already paid for.** 521a's and item 500's lists, plus:
- **P2 comparands are P3** (`CMP_DEGREE`) with `QDEG_P` (quadrature degree 6) on `_rel_p`. Comparing a P2 θ through P1's `_rel` changes the quantity (P28's own trap; `docs/planning/chunks/MAT-21.md:596–669`).
- **Operation order in Q4 is P36's** (`Q * 3.0 * (y·y + z·z) / (B**2 + C**2)`). A reordered expression changes the compiled form and the last digits, and (S1-quad)'s `.6e` would then read a different number.
- **`SpatialCoordinate` in a volume form** compiles at the automatic degree, as P36's does. No facet integral carries it, so the `POST-5` pinned-degree trap is not engaged; never add a facet term with `xs` without `metadata`.
- **Real operands do not save an ordering comparison**: there is none in Q4.

**Scope.** Covered: P36's transverse quadratic source on the m = 4 box at P2, steady and CN to t\*, from a tracked file through `femem run`, equal to the Python Study at `L0_BAND` and reproducing P36's records at `INTERIOR_BAND` / `S2_GAUSS_BAND` / `E_BAND`. **Not claimed:** a duty-cycled transverse source (P38's); other weights or a user-written source expression; the box's EM SAR; a coil; the real build; an h-rate; probes or dose on the box (U60 / P37 / P39); any compliance statement.

**Negative result.** Report and stop; park on `attempt/WF-1-step16-<UTC>`, item 528 BLOCKED in the same commit.
- **(S1-quad) / (T-quad) `.6e` red with the band and (I-Py) green:** Q4's expression differs from P36's in order or constants. Print both forms; fix in-slot only before W1a's clock (rule (i)), else park.
- **A band red:** P36 no longer reproduces. Open a known-issues entry.
- **(U) red:** the P2 box path drifted. Stop.
- **RG red:** revert.

No `PREDICTED RED:` line.

**Status it can move.** `WF-1` stays ✅; the row gains "**Step-16 record (§9 item 528, U64, provisional):** `excitation.source: {kind: transverse_quadratic}` on a `thermal_box` file (`TransverseQuadraticHeatSource`) reproduces P36's steady `7.219878e-05` and CN `2.694340e-04` at P2, complex, `-n 2` (`INTERIOR_BAND`, Gauss and balance at their bands), θ equal to the Python Study at ≤ <x>; dropping the source misses by <p>." `WF-15` stays ✅ and gains a one-clause step-7 note naming the class. §10 U64 DONE. `MAT-21` does not move.

**Spawn line:** 521b's spawn line with `528` / `step 16` / `FEM_EM_WF1_S16=1` / `attempt/WF-1-step16-*`, the window list above, all solve legs in the complex build, and `src` edits limited to `study.py` (Q1–Q4), `__init__.py` and `study_yaml.py` (Q5–Q6), insertions only except Q2's disclosed edit of 521a's admission line.

**For the review.**
1. **Where `source:` sits.** The drafter reads U64's `source: {kind: transverse_quadratic}` as a sub-key of `excitation`, so the excitation kind stays `uniform_heat_source` and means "an imposed source of mean q". The alternative is a new excitation kind `transverse_quadratic_heat_source`, which is cleaner but not the letter.
2. **The class name** `TransverseQuadraticHeatSource`, and its export (45 names).
3. **RG-B** (521a's heavy duty window, ≈ 4.5 min): required, or waived on RG-A1.
4. **`.6e` record equality** in place of item 500's `.4e`. P36 printed six digits, and the build and code path are P36's.

**Outcome (2026-10-19 00:35–00:45 UTC): DONE, green, attempt 1 of 4.** All 13 windows green in order (`pgrep` 0 before W1a and RG-B; RG-B run, no waiver). W0 16 passed (first W0 run red on a test-only bug, `_study_from_mapping(...)[0]` on a `Study`, fixed before W1a's clock; `20261019T003525Z_WF-1-s16-W0.log`); check-family 24 good files. RG1 13 passed. W1a (78 s): (Pdeg) 116 025 ×3; (I-Py) 1.032e-15; (S1-quad) `7.219878e-05` = P36; (S2-quad) 9.775e-14 (P36 9.775e-14); (U) `2.489307e-05` / `1.669359e-05` = step 3; (c-src′) 5.1471e-01, (c-src) 4.4856e-01, (c-z, printed) 4.7893e-02 — all = P36's. W1b (98 s): (Pdeg) 116 025 ×2; (I-Py) 1.213e-14; (T-quad) `2.694340e-04` = P36; (E-quad) 2.288e-13 (P36 2.285e-13); multipliers equal; (c-src-T, printed) 6.9562e-01. RG-A0 9 passed, RG-A1 5 passed, RG-B 5 passed (186 s), RG-C 5 passed, RG-D 3 + 5 passed — the precedents' tallies. Wgen 46 files; Wval 289 nodes (282 + 7: the parse node under one name, each solve node under three); Wapi 22 passed; Wdoc 0 dead. `src`: `study.py` (Q1–Q4, plus Q2's disclosed admission-line edit), `__init__.py`, `study_yaml.py` (Q5–Q6; disclosed edits of two 521b lines: the `excitation` `_unknown` tuple gains `"source"`, and the Study is built from the parsed `heat_source`), `public_surface.txt` 44 → 45. Disclosed: (I-Py) compares P2 owned arrays through a module-local P2 copy of `_on` (the shared helper is P1-only). Logs `20261019T003542Z_WF-1-s16-W0.log` … `20261019T004455Z_WF-1-s16-Wdoc.log`.


### Review ruling 2026-10-18 03:00 (items 521a / 521b U59, 528 U64; 523 in `WF-15.md`, 525 in `WF-13.md`)

Spot-check: `study.py:1681` (`_refuse_f_human(self, "estimate")`) re-read. No `thermal_box` exists in `src/` or `tests/`.
1. **U59 is split.** 521a (`WF-15` step 5) adds `ThermalBoxGeometry` and the steady thermal-only branch through `Study.temperature`. 521b (`WF-1` step 15) adds the `thermal_box` YAML through `femem run`. 521b depends on 521a.
2. **Complex build, and a heavy duty leg.** `MAT-21`'s records are `@complex_only`. P26's comparand window measured 262 s. §10's "real build, standard" is corrected as execution detail, and no target moves.
3. **`wall` is a `None` property on the box**, so `_study_transient_temperature` runs unedited. No edit at `:1455`.
4. **"Wrong wall" means the Robin h doubled.** The control must miss the record by more than the band, labelled *predicted* (rule (e)).
5. **Accepted, disclosed:** the public surface grows from 43 to 45 names (`ThermalBoxGeometry`, `TransverseQuadraticHeatSource`). The WF-14 registry / API-page rule applies (§9 "Every item").
6. **U60:** the helpers raise `StudyConfigError` (a `ValueError`). MAT-19's band is `P_BAND`, which equals `INTERIOR_BAND` at 1e-2. The gate imports `P_BAND`, and only MAT-19's gate is reproduced, not its digits (different mesh).
7. **U64:**
   - `source:` sits under `excitation`, with the kind unchanged;
   - **no RG-B waiver**: every YAML gate is re-run;
   - record equality at P36's printed `.6e` digits is accepted.
8. **U55b (item 525) is not queued; it is returned to the weekly.**
   - `estimate` refuses the F-human first (`study.py:1681`, pinned by `test_wf16_study_f_human.py:155–156`).
   - Item 431b's log carries no unknowns or GiB.
   - A row priced from the file's own log is "inside ×2" by construction.
   - U61 meets the same refusal.
   - The drafted re-pose U55b′ (heavy, `-n 8`) needs the weekly's rulings on Q1 (lift `WF-16` step 8's estimate refusal) and Q2 (a thermal price).


### §9 item 539 — 2026-10-18 (§10 U65, daily — provisional; chain step, solve)

*Drafted by a read-only agent for the 2026-10-18 18:00 daily review; the review reviews it. Code read on `3c143576`; nothing was run.*

**Licence.** §9 rule (6) (`docs/automation/daily-review.md:196–214`): a provisional next step of chain U, tagged "(daily, 2026-10-18 — provisional)", toward chain U's stated target, "each step reproduces a gated result or asserts a gated identity through the API, within its unmoved band" (item 521's licence, this file). **The target is a reproduction of an existing record at its existing band: U60's box probe readings (`WF-15` step 6, §9 item 523) at `INTERIOR_BAND` / `B_BAND`, and θ_file = θ_Python at `L0_BAND`. No band, tolerance, record, Done-when or target is set or moved; every band is imported by name.** Chain U stands at two provisional steps with this item and item 540 (the cap is three). Templates: item 521b (this file) and item 345 (`WF-1` step 6, the sphere's `outputs.probes`). No sheet, κ or ESR ruling is engaged.

**Code state (`3c143576`).**
- `outputs.probes` is parsed only by `_parse_probes` (`src/fem_em_solver/study_yaml.py:1024–1058`), which refuses any `geometry_kind` other than `sphere_in_box` (`:1034–1037`) and any file without `thermal:` (`:1038–1041`).
- The box parser `_parse_thermal_box_mapping` (`:858–973`) admits `outputs.directory` only: `_unknown(osec, ("directory",), "outputs", source)` at `:968`. So a box file with `probes` is refused today at `outputs.probes`, and **521b pins that refusal**: `tests/package/test_wf1_yaml_thermal_box.py:257–258`. This item makes that row false (question 1).
- `_run_thermal_box` (`:1479–1537`) takes no `outputs` argument and returns `StudyRun(…, thermal=thermal)` with `probes` defaulted to `None` (`StudyRun`, `:1295–1304`, a mutable dataclass). The dispatch is at `run_study_file:1623–1624`.
- The sphere route's probe block (`:1727–1736`) is the model for the summary: `centres`, `radii`, `values_k`, `volumes_m3`.
- `Study.probe` (`_study_probe`, `src/fem_em_solver/study.py:1640–1670`) takes keyword `transient` / `duty` since U60 and passes them to `self.temperature`. On the box it re-solves the temperature, as the sphere route does (disclosed at `study_yaml.py:1728–1729`).
- U60's box records are **real-build** readings (`20261018T123416Z_WF-15-s6-W2-box.log:12`, no complex prefix): (R-box) c0…c4 `2.6992e-03`, `2.6936e-03`, `2.7066e-03`, `2.6468e-03`, `2.6443e-03` (`:64–68`); (A-box) 2.495e-13 / 2.509e-13 (`:63`); (c-wall) `3.5712e-01` (`:72`). So this item's solve legs run in the real build, the build of the records, and the YAML gate re-runs run in their own builds (item 500's rule).

- **Executor:** implementer. **Tier:** W0 / RG1 / W9 / W1 / RG-P / RG2 / RG3 / RG5 / RG7a / RG6a / RG6b standard; RG7b (521b's W1b) and RG4 (`-n 8`) **heavy**. `-n 2` everywhere else.
- **Class:** chain step (solve). **Cost:** ≈ 18.0 measured (class median) / ≈ 50 predicted min = 15 fixed + ≈ 15 code (≈ 45 `study_yaml.py` lines, a module of ≈ 220 lines, three good files and one malformed) + windows ≈ 13.4 + ≈ 6 journal. Window prices, all measured on 521b's windows unless named: W0 4 s (`20261018T184753Z_WF-1-s15-W0.log`), RG1 3 s, W9 15 s, W1 predicted 20–60 s (U60's W2 12 s for two `Study.probe` calls on the m = 4 box, `…WF-15-s6-W2-box.log`; this window adds the file runs, the wrong-wall file and one m = 4 CN pair), RG-P 40 s (`20261009T020231Z_WF-1-s6-W1-probes.log.gz`), RG2 62 s, RG3 33 s, RG5 42 s, RG7a 10 s, RG7b 190 s, RG4 211 s, RG6a 63 s, RG6b 75 s, Wgen / Wval / Wapi / Wdoc 2 / 3 / 19 / 1 s (`20261018T190139Z` … `190227Z_WF-1-s15-W*.log`). **Park (rule (d)) if code + W0 / RG1 / W9 are not green by minute 30, RG7b has not started by minute 40, or RG4 cannot start by minute 48.**
- **Depends on:** nothing open (U59 = items 521a / 521b and U60 = item 523 are on `main`). **Family:** `WF-1` step 17, attempt 1 of 4.
- **Item 528 (U64) touches the same parser and runner.** If 528 has landed when the slot starts, read its `_parse_thermal_box_mapping` / `_run_thermal_box` before editing, add its parse node to W0 and its W1a / W1b to the re-gates (labels `WF-1-s17-RG8a` / `-RG8b`, commands from its outcome), and say so. If it has not, item 528's brief (which pins "521b's (X) tables unedited") meets this item's amended row; the review notes it on 528.

**The change (`src/fem_em_solver/study_yaml.py` only; insertions plus one disclosed edited line, rule (c)).**
- **(R1) `_parse_box_probes(value, source)`**, a new function after `_parse_probes`. It is `_parse_probes`'s list-and-entry block (`:1042–1058`) **copied**, not factored (521b's rule: a shared helper would edit the sphere's path), with the same field names (`outputs.probes`, `outputs.probes[i]`, `….centre`, `….radius`). It has no geometry-kind or `has_thermal` check: a box file always carries `thermal:` (`_BALL_REQUIRED`).
- **(R2) The one edited line:** `:968` becomes `_unknown(osec, ("directory", "probes"), "outputs", source)`. Then one inserted pair after the `directory` lines (`:969–970`): `if "probes" in osec: outputs["probes"] = _parse_box_probes(osec["probes"], source)`.
- **(R3) `_run_thermal_box_probes(study, kind, cooling, outputs, path, out, digest, started, comm)`**, a new function. It calls `_run_thermal_box(...)` unedited, then `probes = study.probe(outputs["probes"], cooling.cooling, transient=cooling.transient, duty=cooling.duty)` (collective). It writes `run.summary["probes"]` with the sphere route's four keys (`:1731–1736`, copied), re-stamps `elapsed_s`, rewrites `summary.json` on rank 0, `Barrier`, sets `run.probes = probes` and returns `run`. Disclosed: on a probes file, `summary.json` is written twice, and the second write is the record.
- **(R4)** In `run_study_file`, inserted **before** `:1623`: `if isinstance(study.geometry, ThermalBoxGeometry) and "probes" in outputs: return _run_thermal_box_probes(...)`. A box file without probes reaches 521b's dispatch at `:1623` unchanged.
- **Docstrings:** one clause in the module docstring's `thermal_box` bullet (`:52–53`) and in `run_study_file`'s (`:1612`).
- **Never:** edit `_parse_probes`, `_parse_outputs`, `_run_thermal_box`, `_run_thermal_ball` or `_parse_thermal_ball_mapping`; admit probes on a `thermal_ball` file (its refusal, `tests/package/test_wf1_yaml_thermal_ball.py:317–318`, stays green unedited); touch `study.py` or `cli.py`; add to `study_yaml`'s `__all__`.

**Tracked inputs (`tests/package/wf1/`; each opens `# WF-1 step 17 (§9 item 539): …`, then `femem: 1`; floats as Python `repr` digits).**
- `thermal_box_steady_m4_probes.yaml`: 521b's `thermal_box_steady_m4.yaml` plus `outputs: {probes: [...]}`. That is `WHOLE_BOX` then U60's five `BOX_PROBES` in order: radius 0.4 at the origin, then radius 0.0135 at (0, 0, 0), (0, 0.105, 0), (0, −0.105, 0), (0.1625, 0, 0), (−0.1625, 0, 0). (K) asserts each against `tests/package/test_wf15_study_probe_thermal.py:85–94` by `float.hex`. 0.105 and 0.1625 are `0.5 * B_BOX` / `0.5 * A_BOX`; halving is exact, so the repr of the product parses to the same double.
- `thermal_box_wrong_wall_m4_probes.yaml`: the same with `h_w_per_m2_k: 20.0` and the c0 probe only (U60's (c-wall) used `BOX_PROBES[:1]`, `:273–274`).
- `thermal_box_duty_m4_probes.yaml`: m = 4 (`divisions: [52, 32, 8]`), 521b's duty `transient` / `duty` blocks (`dt_s: 70.5375`, `n_steps: 40`, `rho_c_j_per_m3_k: 4180000.0`, `period_s: 564.3`, `duty: 0.5`), and `BOX_PROBES`. This is the pass-through file. **No record exists for it, and none is claimed.**
- `malformed_thermal_box_probe_radius_zero.yaml` → `outputs.probes[1].radius`.

**Test module: new `tests/package/test_wf1_yaml_thermal_box_probes.py`.** It takes 521b's module shape (`tests/package/test_wf1_yaml_thermal_box.py`; `_knob`, `_load`, `_say`). Knob `FEM_EM_WF1_S17=1`; when set, assert `not IS_COMPLEX` and size 2, never skip on a mismatch. The first docstring paragraph backticks only `WF-1` and `WF-15`.
- **Imports (module-level names only):**
  - `A_BOXPROBE`, `BOX_CENTRES`, `BOX_PROBES`, `WHOLE_BOX`, `BOX_DIVISIONS`, `_series_at`, `_domain_mean` from `tests/package/test_wf15_study_probe_thermal.py:85–138`;
  - `HALF`, `Q` (`test_mat21_box_closed_form.py`); `K_GEL`, `H_GEL`, `T_INF` (`test_mat15_coil_driven_heating_phantom.py:82–84`); `DT_BOX`, `N_STEPS`, `HALF_T` (`test_mat21_box_duty_closed_form.py`); `RHO_C` (`test_mat12_transient_sar_heating_sphere.py:68`);
  - `INTERIOR_BAND` (`test_mat13_…:64`), `B_BAND` (`test_mat10_sar_to_heat_sphere.py:61`), `L0_BAND` (`test_mat5_em_coupled_sphere.py:104`);
  - `_same_mesh`, `_on`, `_theta_owned_rel` (`tests/package/test_wf15_study_transient.py:117–152`); `_spec_diff`, `_cmp`, `_rc_all`, `_fresh_dir` (`tests/package/test_wf1_study_yaml.py:120–169`); `IS_COMPLEX`.
- **Record constants (transcribed with log lines):** `U60_R_BOX = ("2.6992e-03", "2.6936e-03", "2.7066e-03", "2.6468e-03", "2.6443e-03")` (`20261018T123416Z_WF-15-s6-W2-box.log:64–68`); `U60_C_WALL = 3.5712e-01` (`:72`); `U60_A_BOX = 2.495e-13` (`:63`, printed beside, never compared).
- **Nodes:**
  - `test_wf1_s17_parse_and_refusals` (no mesh; every `build` patched to raise; 0 calls);
  - `test_wf1_s17_box_probes_reproduce_u60` (W1, steady and wrong wall);
  - `test_wf1_s17_transient_probes_equal_python` (W1, the pass-through).
- **Registry:** under `# WF-1 step 17 (§9 item 539)`, the parse node under `[names.load_study]` / `[names.StudyConfigError]`, and the solve nodes under `[names.run_study_file]` / `[names.ThermalBoxGeometry]` (`tests/package/validation_envelope.toml`). **The §9 standing rule applies** (`PROJECT_PLAN.md:4060`): Wgen / Wval / Wapi / Wdoc in this slot, with the regenerated `docs/guide/api/` and `docs/guide/validation.md` committed. Wval's row count is the last committed page's plus 6 (521b's finding: the stat counts (name, node) rows, 3 nodes × 2 names). Print it.

**Anchors (asserted; bands imported by name).**
- **(K) Parse, no mesh:**
  - `_spec_diff(load_study(f), Study(ThermalBoxGeometry(HALF, BOX_DIVISIONS), None, UniformHeatSource(Q))) == []` on the three good files;
  - `outputs["probes"]` centres and radii equal `[WHOLE_BOX, *BOX_PROBES]` (wrong wall: `BOX_PROBES[:1]`) by `float.hex`;
  - the cooling through `_cmp` by hex (`2.0 * H_GEL` on the wrong-wall file); duty file `TransientSpec(DT_BOX, N_STEPS, RHO_C)` and `DutyCycle(2.0 * HALF_T, 0.5)` by hex;
  - `cli.main(["check", f])` rc (0, 0) via `_rc_all`, one `family: thermal_box` line on rank 0.
- **(X) Refusals naming their field:**
  - the malformed file (also through `run` and `check`: rc (2, 2), `out` absent by `MPI.LOR`);
  - `outputs.probes: []` → `outputs.probes`; an entry `[1]` → `outputs.probes[0]`; an extra key → `outputs.probes[0].<key>`; a missing `radius` → `outputs.probes[0].radius`; a two-element `centre` → `outputs.probes[0].centre`.
  - **Pinned, unedited, green:** the ball's `outputs.probes` refusal (`test_wf1_yaml_thermal_ball.py:317–318`, in W0) and the sphere-without-`thermal:` refusal (`malformed_probes_without_thermal.yaml`, `test_wf1_yaml_probes.py:56`, in RG1).
- **(UT) / (I-Py), steady:** the file run's probe mesh equals the Python Study's (`_same_mesh`; 79 872 cells). The Python twin is U60's exact call (`test_wf15_study_probe_thermal.py:241–242`). Then:
  - θ_file vs θ_Python: `_theta_owned_rel(_on(py_mesh, run.probes.theta), py.theta, comm)` ≤ **`L0_BAND`** (predicted 0 to ≈ 5e-14; 521b read 1.479e-15);
  - the probe values file vs Python, printed only;
  - `run.thermal` vs `run.probes.thermal` θ (the re-solve) at `L0_BAND`, as step 6's (`test_wf1_yaml_probes.py:128`).
- **(R-box) U60's records through the file:**
  - for each of the five probes, `abs(P_j / (_series_at(c_j) − Q·A_BOXPROBE²/(10·K_GEL)) − 1)` ≤ **`INTERIOR_BAND`**, and `f"{r:.4e}" == U60_R_BOX[j]`. The digits are asserted: same build, same mesh, same `Study.probe` call, tied by (I-Py);
  - **(A-box)** the whole-box probe |P₀ − ∫θ/V_h| / |∫θ/V_h| ≤ **`B_BAND`** and |V₀/V_h − 1| ≤ `B_BAND` (`_domain_mean`), printed beside `U60_A_BOX`.
- **(U′) The summary is its own result:** `summary["probes"]["values_k"]` equals `[float(v) for v in run.probes.values]` exactly. JSON on disk equals the in-memory summary (step 6's (U′), `test_wf1_yaml_probes.py:158–166`). `family == "thermal_box"`; `run.result is None`.
- **(I-tr) The pass-through, an identity only:** the duty m = 4 probes file against `Study(...).probe(BOX_PROBES, cooling, transient=TransientSpec(DT_BOX, N_STEPS, RHO_C), duty=DutyCycle(2.0 * HALF_T, 0.5))`:
  - θ ≤ **`L0_BAND`** (`_theta_owned_rel` through `_on`);
  - `source_multipliers` `array_equal` between the two;
  - the probe values printed. No comparand is posed (U60's scope: "no exact ball mean is posed" on a transient box, `docs/planning/chunks/WF-15.md:588`).

**Negative controls (rule (e)).**
- **(c-wall) The wrong wall through a tracked file — *asserted, backed*.** The wrong-wall file's c0 reading against the right-h ball mean misses: `> INTERIOR_BAND`, printed `.4e` beside `U60_C_WALL`. It is U60's own comparison (`…W2-box.log:72`) on the same mesh in the same build, so the digits are predicted equal; equality is printed, not asserted. One-sided; the 3.6e-1 miss is far from saturation.
- **(c-ball) The refusal that must not move — *asserted*:** the ball file with `outputs.probes` still raises at `outputs.probes` (W0, unedited). The sphere-without-`thermal:` file still raises (RG1). Both show the admission is box-only.
- **(c0) By construction:** R4's `isinstance … and "probes" in outputs` is false for every file without probes. RG7a / RG7b (521b's own nodes, whose (I=) pins the summary keys, `20261018T185143Z_WF-1-s15-W1a-steady.log:67`) show the no-probe box route is unchanged.

**Windows.** Foreground; Bash `timeout: 660000`; pytest `-s`; node IDs, never `-k`; **the harness call is the whole Bash command: bare, repo-relative, nothing chained or piped.** `pgrep` = `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reads 0 before W1, RG2, RG3, RG5, RG7b, RG4 and RG6a. Run W0 only after every edit (rule (i)). Order: W0, RG1, W9, W1, RG-P, RG2, RG3, RG5, RG7a, RG7b, RG4, RG6a, RG6b, Wgen, Wval, Wapi, Wdoc.
- **W0** (complex, no mesh; the YAML gate modules, 521b's W0 plus the new node; check-family's glob grows from 21 to 24 good files, or more if 528 has landed):
  `scripts/testing/run_and_log.sh WF-1-s17-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_thermal_box_probes.py::test_wf1_s17_parse_and_refusals tests/package/test_wf1_yaml_thermal_box.py::test_wf1_s15_parse_and_refusals tests/package/test_wf1_yaml_thermal_ball.py::test_wf1_s14_parse_and_refusals tests/package/test_wf1_yaml_capacitors.py::test_wf1_s13_parse_and_refusals tests/package/test_wf1_yaml_tune.py::test_wf1_s11_tune_parse_and_refusals tests/package/test_wf1_yaml_tune_f_human.py::test_wf1_s12_parse_and_refusals tests/package/test_wf1_check_family.py tests/package/test_wf1_yaml_f_human.py tests/package/test_wf1_yaml_f_human_thermal.py -v -s --tb=short'"`
  If 528 has landed, insert `tests/package/test_wf1_yaml_thermal_box_quadratic.py::test_wf1_s16_parse_and_refusals` after the s15 node.
- **RG1** (complex; every other `[names.load_study]` parse node, including step 6's sphere probes):
  `scripts/testing/run_and_log.sh WF-1-s17-RG1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_study_yaml.py::test_wf1_good_files_build_the_python_studies tests/package/test_wf1_study_yaml.py::test_wf1_malformed_files_name_their_field tests/package/test_wf1_yaml_thermal.py::test_wf1_thermal_yaml_parse_and_refusals tests/package/test_wf1_yaml_copper.py::test_wf1_copper_yaml_parse_and_refusals tests/package/test_wf1_yaml_sar_b1.py::test_wf1_sar_b1_yaml_parse_and_refusals tests/package/test_wf1_yaml_probes.py::test_wf1_probes_yaml_parse_and_refusals tests/package/test_wf1_yaml_f_human.py::test_wf1_f_human_yaml_parse_and_refusals tests/package/test_wf1_yaml_copper_sar_b1.py::test_wf1_s8_parse_and_refusal_parity tests/package/test_wf1_yaml_f_human_thermal.py::test_wf1_fht_parse_and_refusals -v -s --tb=short'"`
- **W9** (real installed; the installed package carries the edited `study_yaml.py`, and the step copies `wf1/*.yaml`):
  `scripts/testing/run_and_log.sh WF-1-s17-W9-installed-femem "docker compose exec -T fem-em-solver bash -lc 'timeout -k 10 50 bash /workspace/scripts/package/wf13_setup.sh real && mkdir -p /tmp/wf13/run/tests/package/wf1 && cp /workspace/tests/package/test_wf1_installed_cli.py /tmp/wf13/run/tests/package/ && cp /workspace/tests/package/wf1/*.yaml /tmp/wf13/run/tests/package/wf1/ && cd /tmp/wf13/run && env -u PYTHONPATH WF13_INSTALLED_CHECK=1 timeout -k 10 90 /tmp/wf13/venv/bin/python -m pytest -p no:cacheprovider tests/package/test_wf1_installed_cli.py -v -s --tb=short'"`
- **W1** (real, standard; predicted 20–60 s):
  `scripts/testing/run_and_log.sh WF-1-s17-W1-probes "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace FEM_EM_WF1_S17=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_thermal_box_probes.py::test_wf1_s17_box_probes_reproduce_u60 tests/package/test_wf1_yaml_thermal_box_probes.py::test_wf1_s17_transient_probes_equal_python -s -v --tb=short'"`
- **RG-P** (complex; the sphere `outputs.probes` route, unchanged by construction; step 6's W1 verbatim, `20261009T020231Z_WF-1-s6-W1-probes.log.gz:12`, 40 s):
  `scripts/testing/run_and_log.sh WF-1-s17-RGP-sphere-probes "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_probes.py -s -v --tb=short'"`
- **RG2** (complex; tune; `20261018T184902Z_WF-1-s15-RG2.log:12`, 62 s):
  `scripts/testing/run_and_log.sh WF-1-s17-RG2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 300 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_tune.py -v -s --tb=short'"`
- **RG3** (complex; the sphere CLI; `…185011Z_WF-1-s15-RG3.log:12`, 33 s):
  `scripts/testing/run_and_log.sh WF-1-s17-RG3 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_study_yaml.py::test_wf1_cli_sphere_reproduces_wf15 -v -s --tb=short'"`
- **RG5** (complex; capacitors; `…185055Z_WF-1-s15-RG5.log:12`, 42 s):
  `scripts/testing/run_and_log.sh WF-1-s17-RG5 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 300 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_capacitors.py -s -v --tb=short'"`
- **RG7a / RG7b** (complex; 521b's W1a / W1b verbatim, `…185143Z_WF-1-s15-W1a-steady.log`, 10 s; `…185201Z_WF-1-s15-W1b-duty.log`, 190 s, heavy):
  `scripts/testing/run_and_log.sh WF-1-s17-RG7a-box-steady "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_WF1_S15=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_thermal_box.py::test_wf1_s15_steady_reproduces_mat21_m4 -s -v --tb=short'"`
  `scripts/testing/run_and_log.sh WF-1-s17-RG7b-box-duty "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_WF1_S15=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 590 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_thermal_box.py::test_wf1_s15_duty_reproduces_p26 -s -v --tb=short'"`
- **RG4 — required, alone** (complex, heavy, `-n 8`; R4's head runs ahead of every F-human file's dispatch; `…185519Z_WF-1-s15-RG4.log:12`, 211 s):
  `scripts/testing/run_and_log.sh WF-1-s17-RG4 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_WF1_S12=1 timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_tune_f_human.py -s -v --tb=short'"`
- **RG6a / RG6b** (real; the ball route beside R4; `…185901Z_WF-1-s15-RG6a-ball-duty.log:12`, 63 s; `…190014Z_WF-1-s15-RG6b-ball-wall.log:12`, 75 s):
  `scripts/testing/run_and_log.sh WF-1-s17-RG6a-ball-duty "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace FEM_EM_WF1_S14=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_thermal_ball.py::test_wf1_s14_parse_and_refusals tests/package/test_wf1_yaml_thermal_ball.py::test_wf1_s14_duty_reproduces_mat17 -s -v --tb=short'"`
  `scripts/testing/run_and_log.sh WF-1-s17-RG6b-ball-wall "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace FEM_EM_WF1_S14=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_thermal_ball.py::test_wf1_s14_wall_reproduces_mat18_2b -s -v --tb=short'"`
- **Wgen / Wval / Wapi / Wdoc** (521b's four verbatim, `20261018T190139Z` … `190227Z_WF-1-s15-W*.log:12`):
  `scripts/testing/run_and_log.sh WF-1-s17-Wgen "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 python3 scripts/docs/generate_api_reference.py --out docs/guide/api'"`
  `scripts/testing/run_and_log.sh WF-1-s17-Wval "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace timeout -k 30 60 python3 scripts/docs/generate_validation_page.py --out docs/guide/validation.md'"`
  `scripts/testing/run_and_log.sh WF-1-s17-Wapi "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace timeout -k 30 170 python3 -m pytest tests/package/test_wf14_validation_envelope.py tests/package/test_wf17_api_reference.py tests/package/test_wf14_validation_page.py -v -s --tb=short'"`
  `scripts/testing/run_and_log.sh WF-1-s17-Wdoc "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && timeout -k 30 60 python3 scripts/testing/check_example_doc_references.py --docs-root docs/guide --check-links'"`
  Wdoc reads `0 dead`.
- **Not re-run (disclosed):** U60's and 521a's API windows (no `study.py` line moves); `WF-13`'s estimate module (no `cli.py` line moves).
- **`git diff --stat main -- src/`** shows `study_yaml.py` only: insertions, plus R2's one edited line, named in the commit.

**Traps already paid for.** Item 521's lists (521a's and 521b's, this file) and item 500's apply whole. In particular:
- never `float.hex` across builds; YAML 1.1 loads `4.18e6` as a string; `femem` prints on rank 0 only; `-k a or b` splits argv; `_load`'s 4-tuple is load-bearing; `__all__` does not grow;
- **UFL cannot mix meshes**: compare through `_on` or in NumPy;
- the thermal degree is `RobinCooling.degree`, never `Study.degree`;
- orphan `.c` after a kill: `pgrep` first, then at 0 the rubric sweep as `WF-1-s17-cache-sweep` in the killed build; never raise a timeout or change `-n`;
- rule (i): an edit to `study_yaml.py` after W0 re-runs W0 and every RG.

**New here:**
- **The probe is an owned-cell integral** (`MAT-19`'s contract; item 227's deadlock). Every rank enters `Study.probe`, including a rank that owns no cell of a probe ball. R3 calls it unconditionally on every rank.
- **The real build here, the YAML gates in theirs.** W1 has no complex prefix, and `_require` asserts `not IS_COMPLEX`. The digits in `U60_R_BOX` are real-build digits.
- **The probe re-solves the temperature.** `run.thermal` and `run.probes.thermal` are two solves. (I-Py) compares `run.probes.theta`, the one the values come from.
- **The duty file's floats** (`70.5375`, `564.3`) are 521b's reprs, already hex-checked against `DT_BOX` / `2.0 * HALF_T` at `test_wf1_yaml_thermal_box.py`'s (K). Copy them, never retype.
- Write ≈, never a tilde.

**Scope.** Covered: `outputs.probes` on `thermal_box` files through `femem run` / `run_study_file`, steady at m = 4 (P1, real, `-n 2`), reproducing U60's five (R-box) readings and the whole-box identity at their bands, equal to the Python Study at `L0_BAND`; the transient / duty pass-through as an identity only. **Not claimed:**
- probes on a `thermal_ball` file (still refused);
- any comparand on a transient or duty box probe;
- probes at P2 or with U64's transverse source (none on `main` at drafting);
- scalar point evaluation (P37);
- the complex build for the probe legs;
- the installed console script running the new files (W9 re-gates the installed `check`);
- a probe's own volume or self-heating; any compliance statement.

**Negative result.** Report and stop. Park on `attempt/WF-1-step17-<UTC>`, mark item 539 BLOCKED with its unblock condition in the same commit (rule (d)), and journal it in attempts.md.
- **(K) / (X) red:** a parse defect. Fix in-slot before W1 only.
- **(I-Py) red with (UT) green:** two equal specs on one build give θ beyond `L0_BAND`. That is a solve-floor finding; the band is unmoved.
- **(R-box) `.4e` red with its band green:** an environment drift since U60's window (`20261018T123416Z`). Send it to `record-reconciler` (1*); never re-register in-slot.
- **A band red:** U60 no longer reproduces through the file. Open a known-issues entry.
- **(c-ball) red:** the admission leaked into the ball route. Revert and stop.
- **Any RG red:** R2–R4 changed more than the box-probes path. Revert `src/` and stop.
- **Status 124:** not a result; park with the measured time.

No `PREDICTED RED:` line: a red anywhere stops the slot.

**Status it can move.** `WF-1` stays ✅ (rule (j): its Done-when is step 1's). On green the §7 row gains: "**Step-17 record (§9 item 539, U65, provisional):** `outputs.probes` on a `thermal_box` file drives `Study.probe` through `run_study_file`, real build, `-n 2`: U60's five box-probe readings reproduce to their printed digits (`INTERIOR_BAND`), the whole-box probe identity ≤ <a> (`B_BAND`), θ_file = θ_Python ≤ <x> (`L0_BAND`), the duty pass-through equal to the Python Study ≤ <y>; the wrong-wall file misses by <p>; `outputs.probes` stays refused on `thermal_ball` and without `thermal:`; `src`: `study_yaml.py`, insertions plus one disclosed edited line." §10 U65 DONE (provisional). `WF-15` stays ✅; the "probes from YAML on thermal-only files" clause of item 523's scope (`docs/planning/chunks/WF-15.md:590`) narrows to the ball. `MAT-19` / `MAT-21` do not move.

**Spawn line:** "Use the implementer agent to execute `WF-1` §9 item 539 (U65) per `docs/planning/chunks/WF-1.md` "§9 item 539" (read the brief and any review ruling beneath it first; the rulings hold). The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Run W0, RG1, W9, W1, RG-P, RG2, RG3, RG5, RG7a, RG7b, RG4, RG6a, RG6b, Wgen, Wval, Wapi, Wdoc in order, in the foreground, every harness Bash call with `timeout: 660000`, pytest with `-s`, container `timeout -k 30`. Before W1, RG2, RG3, RG5, RG7b, RG4 and RG6a, run `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` and require 0. W1 runs in the real build with `FEM_EM_WF1_S17=1`; every re-gate runs in its own recorded build. `src` edits: `study_yaml.py` only, insertions plus the one disclosed `:968` line, named in the commit. Edit `tests/package/test_wf1_yaml_thermal_box.py:257–258` only as the review ruling on this item says, and never any other existing test assertion. Never import a `test_*` name; never compare across builds by `float.hex`. You are the executor: do not spawn agents, and never background a harness window. On a red, park on `attempt/WF-1-step17-*`, mark item 539 BLOCKED in the same commit, and report."

**For the review.**
1. **521b's pinned row must change.** `test_wf1_yaml_thermal_box.py:257–258` asserts that a box file with one valid probe is refused at `outputs.probes`; this item's letter makes it admitted. The drafter proposes a **disclosed one-line edit of that row's input** to `{"probes": []}`, which keeps the expected field `outputs.probes` (R1's empty-list refusal). The alternative is to delete the row and pin the refusal in the new module only. Either is an edit of an existing test assertion and needs this review's ruling. Without one, the item is not runnable.
2. **Build.** The solve legs run in the real build, the build of U60's records, against 521b's complex box legs. The alternative is complex with U60's digits printed only.
3. **R3 as a wrapper**, so `_run_thermal_box` stays byte-identical, at the cost of a second `summary.json` write. The alternative is an edit of `_run_thermal_box`'s signature, call and return line (three disclosed edits).
4. **(I-tr)** admits probes on transient / duty box files with an identity only. The alternative refuses `outputs.probes` with `thermal.transient` (no comparand), which is a narrower admission.

**Outcome (2026-10-19 03:33–03:51 UTC, 22:30 CDT slot): DONE, green, attempt 1 of 4.** All 19 windows green, run one at a time in the brief's order, with 528's two re-gates inserted after RG6b; `pgrep` 0 before W1, RG2, RG3, RG5, RG7b, RG4, RG6a, RG8a and RG8b. **Item 528 had landed**, so its `_parse_thermal_box_mapping` / `_run_thermal_box` were read before editing. The `:968` line of the brief is `:985` on this `main` (528 inserted lines above it). W0 17 passed, 6 s, with 528's `test_wf1_s16_parse_and_refusals` after the s15 node; check-family 27 good files. RG1 13 passed. W9 1 passed. **W1 (real, `-n 2`, 25 s):** (UT) 79 872 / 79 872 cells, geometry equal; (I-Py) 1.184e-15; temperature block vs the probe's re-solve 1.184e-15; (A-box) 2.495e-13 / 2.509e-13 (U60 2.495e-13); (R-box) c0…c4 `2.6992e-03`, `2.6936e-03`, `2.7066e-03`, `2.6468e-03`, `2.6443e-03`, all equal to U60's digits; (c-wall) 3.5712e-01 = U60; (U′) values = `run.probes.values`, JSON = memory on all three runs. (I-tr): 41 multipliers equal, θ ≤ 6.912e-14; the probe values are printed only. RG-P 9 passed; RG2 10 passed (63 s); RG3 5 passed; RG5 10 passed; RG7a 5 passed; RG7b 5 passed (187 s); RG4 `-n 8` 9 passed (209 s); RG6a 3 passed 3 skipped; RG6b 2 passed 3 skipped; **RG8a** (528's W1a, `WF-1-s17-RG8a`) 5 passed (74 s); **RG8b** (528's W1b) 5 passed (98 s). Wgen 46 files; Wval 301 nodes, the committed page's rows + 6, as predicted; Wapi 22 passed; Wdoc 0 dead. **Disclosed:** (1) per the review ruling, 521b's refusal row (`test_wf1_yaml_thermal_box.py:257–258`) moved its input to `outputs: {probes: []}`, which still raises at `outputs.probes` (W0, line 108 of the log). No other existing assertion was edited. 528's "(X) tables unedited" clause now reads "unedited except for 539's disclosed row". (2) (I-Py) and (I-tr) compare **temperatures** (`thermal.temperature` through `_on`), not `theta`. `_theta_owned_rel` normalises by max\|T − T∞\|, so feeding it θ (already a rise) would divide by ≈ 20 K. (3) The (X) table also pins `outputs.probes` given as a mapping, and `outputs.sar` beside probes, which is refused on a box. (4) `WF-15.md:590`'s "probes from YAML on thermal-only files" clause now narrows to the ball. `src`: `study_yaml.py` +67 / −1: R1 `_parse_box_probes`, R2 the one edited `_unknown` line plus the inserted pair, R3 `_run_thermal_box_probes`, R4 the dispatch head, and the `run_study_file` docstring. The module docstring's `thermal_box` bullet was not edited. Logs `20261019T033322Z_WF-1-s17-W0.log` … `20261019T035110Z_WF-1-s17-Wdoc.log`.


### §9 item 540 — 2026-10-18 (§10 U66, daily — provisional; chain step, solve)

*Drafted by a read-only agent for the 2026-10-18 18:00 daily review; the review reviews it. Code read on `3c143576`; nothing was run.*

**Licence.** §9 rule (6), a provisional chain-U step, tagged "(daily, 2026-10-18 — provisional)", toward chain U's stated target (above). **The target is a reproduction of an existing record at its existing band: P38 (`MAT-21` step 11, §9 item 524), the duty-cycled transverse source, (T-qduty) interior `2.543310e-04` ≤ `INTERIOR_BAND` (`20261018T124640Z_MAT-21-step11-W2-quadratic-duty-fem.log:82`), and θ_file = θ_Python at `L0_BAND`. No band, tolerance, record, Done-when or target is set or moved.** Item 528's brief scoped this out: "Not claimed: a duty-cycled transverse source (P38's)" (this file, item 528 "Scope"). Templates: item 528 and 521b.

**Gate: serial on item 528 (U64) landing. If item 528 is not DONE on `main` when the slot starts, skip to the next item, at no attempt cost.** At drafting, 528 is OPEN (`PROJECT_PLAN.md:4096`).

**Code state.**
- P38's fixture (`tests/validation/test_mat21_box_quadratic_duty_fem.py:82–185`) is m = 4 (`M4`), P2 (`DEGREE`), and CN N = `N_MAIN` with `dt = T_STAR / N_MAIN`, θ = ½ and two startup half-steps (`_cn`, `:82–88`). Its source is `Q * 3.0 * (xs[1] * xs[1] + xs[2] * xs[2]) / (B_BOX**2 + C_BOX**2)` (`:121`), the operation order 528's Q4 copies, and its schedule is `_chi_box`.
- P38's readings: (T-qduty) `2.543310e-04` (`:82`), whole 3.319219e-04; (E-qduty) 2.288e-13 (`:83`); (c-step) 5.0333e-01 (`:84`); (S0) `11111000011110000111100001111000011110000` (`:81`; `SCHEDULE`, module `:65`); P2 dofs 116 025 (`:89`).
- **There is no named record constant for (T-qduty)** in the module: it is only printed. This item transcribes it from the log, as 528 transcribed P36's.
- 528's Q4 passes `duty` through as `source_schedule` on the quadratic transient branch (item 528 "(Q4)"), and 521b's parser already admits `thermal.duty` beside `thermal.transient` (`study_yaml.py:951–963`). **If 528 lands as briefed, this item needs no `src/` edit**: one tracked file and one test module. **First action of the slot:** read the landed `_study_box_quadratic_temperature`. If it does not pass `duty` as `source_schedule` (it refuses `duty`, or drops it), stop before any edit and report it; the review re-poses. Nothing has run, so no attempt is consumed.

- **Executor:** implementer. **Tier:** W0 standard (complex, no mesh); W1 **heavy** (complex, `-n 2`; predicted 120–200 s against P38's 170 s window, which ran four solves and three P3 interpolations; this window runs two P2 CNs at ≈ 33.8 s each, the duty and continuous comparands at 27.2 s / 20.7 s, `:72–79`, plus start-up); page gates standard.
- **Class:** chain step (solve). **Cost:** ≈ 18.0 measured (class median) / ≈ 32 predicted min = 15 fixed + ≈ 8 code (module ≈ 170 lines, one good file) + windows ≈ 3.8 (W0 0.1; W1 ≈ 3.3; Wgen / Wval / Wapi / Wdoc 0.4, measured 2 / 3 / 19 / 1 s on 521b's) + ≈ 5 journal. **Park if the module + W0 are not green by minute 25, or W1 has not started by minute 30.**
- **Depends on:** item 528 (DONE on `main`). **Family:** `WF-1` step 18, attempt 1 of 4.

**Tracked input (`tests/package/wf1/`; repr floats).** `thermal_box_quadratic_duty_m4.yaml`: 528's `thermal_box_quadratic_cn_m4.yaml` verbatim (P2, `source: {kind: transverse_quadratic}`, `transient: {dt_s: 70.5375, n_steps: 40, rho_c_j_per_m3_k: 4180000.0}`) plus `duty: {period_s: 564.3, duty: 0.5}` (521b's duty block). Its header is `# WF-1 step 18 (§9 item 540): …`.

**Test module: new `tests/package/test_wf1_yaml_thermal_box_quadratic_duty.py`.** Knob `FEM_EM_WF1_S18=1`; when set, assert `IS_COMPLEX` and size 2.
- **Imports (module-level names only):**
  - `SCHEDULE` (`test_mat21_box_quadratic_duty_fem.py:65`);
  - `_quad_duty` (`test_mat21_box_quadratic_duty_closed_form.py`, as imported at `…_duty_fem.py:61`);
  - `_quad_trans_ex` (`test_mat21_box_quadratic_fem.py`);
  - `DEGREE`, `CMP_DEGREE`, `QDEG_P`, `_exact_p`, `_rel_p` (`test_mat21_box_fem_p2.py`);
  - `T_STAR`, `_interior_chi` (`test_mat21_box_fem.py`); `M4` (`test_mat21_box_separable_fem.py:57`); `DT_BOX`, `N_STEPS`, `HALF_T`, `DUTY` (`test_mat21_box_duty_closed_form.py`);
  - `HALF`, `Q`, `K_GEL`, `H_GEL`, `T_INF`, `RHO_C`, `E_BAND`, `INTERIOR_BAND`, `L0_BAND`, `_same_mesh`, `_theta_owned_rel`, as item 539 and 521a import them;
  - 528's module constants by name if it defines them (`P36_DOFS`).
- **Record constants (transcribed, log `20261018T124640Z_MAT-21-step11-W2-quadratic-duty-fem.log`):** `P38_T_QDUTY = "2.543310e-04"` (`:82`); `P38_E = 2.288e-13` (`:83`); `P38_C_STEP = 5.0333e-01` (`:84`); `P38_DOFS = 116025` (`:89`).
- **Nodes:**
  - `test_wf1_s18_parse_and_refusals` (no mesh);
  - `test_wf1_s18_quadratic_duty_reproduces_p38` (W1).
- **Registry:** under `# WF-1 step 18 (§9 item 540)`, the parse node under `[names.load_study]` / `[names.StudyConfigError]`, and the solve node under `[names.run_study_file]` / `[names.TransverseQuadraticHeatSource]` / `[names.DutyCycle]`. **The §9 standing rule applies:** Wgen / Wval / Wapi / Wdoc in this slot, pages committed. Print the Wval row count, expected +5 rows.

**Anchors (asserted; bands imported by name).**
- **(K) Parse, no mesh:**
  - `_spec_diff(load_study(f), Study(ThermalBoxGeometry(HALF, (52, 32, 8)), None, TransverseQuadraticHeatSource(Q))) == []`;
  - cooling `RobinCooling(K_GEL, H_GEL, T_INF, degree=2)`, `TransientSpec(DT_BOX, N_STEPS, RHO_C)` and `DutyCycle(2.0 * HALF_T, DUTY)` by hex;
  - `check` prints `excitation: TransverseQuadraticHeatSource`, rc (0, 0).
- **(X):** the same file with `transient` removed raises at `thermal.duty` (521b's refusal, now on the quadratic file). 528's and 521b's (X) tables stay pinned, as amended by item 539's ruling if 539 landed first.
- **(Pdeg)** P2 dofs `== P38_DOFS` on both solves. **(UT)** the file's mesh equals the Python Study's (`_same_mesh`; 79 872 cells).
- **(S0)** `"".join(str(int(m)) for m in summary["temperature"]["source_multipliers"]) == SCHEDULE`, and `array_equal` with the Python Study's multipliers.
- **(I-Py)** θ_file vs θ of `Study(ThermalBoxGeometry(HALF, (52, 32, 8)), None, TransverseQuadraticHeatSource(Q)).temperature(RobinCooling(..., degree=2), transient=..., duty=...)` ≤ **`L0_BAND`** (predicted 0 to ≈ 5e-14). At P2, compare the owned arrays after asserting equal P2 `tabulate_dof_coordinates`, or use 528's landed (I-Py) helper (trap below).
- **(T-qduty)** the interior `_rel_p(msh, θ_file, ex_d, _interior_chi(msh))` with `ex_d = _exact_p(msh, lambda x, y, z: _quad_duty(x, y, z, T_STAR))` ≤ **`INTERIOR_BAND`** and `f"{e:.6e}" == P38_T_QDUTY`. The `.6e` equality is P38's print precision, accepted for P36 by the 03:00 ruling (point 7) on the same grounds: same build, fixture and code path, tied by (I-Py).
- **(E-qduty)** `summary["temperature"]["max_balance_residual"]` ≤ **`E_BAND`**, printed beside `P38_E`.
- **(I=)** the summary equals its own result by hex (item 500's (I=), `elapsed_s` on rank 0 only); `summary["temperature"]["source"] == "transverse_quadratic"` (528's Q6); the `duty` block equals `{period_s, duty}` by hex.

**Negative controls (rule (e)).**
- **(c-step) The unscheduled comparand — *asserted, backed*.** The file's duty θ against `_quad_trans_ex(T_STAR)` (interior) is `> INTERIOR_BAND`, printed `.4e` beside `P38_C_STEP`. It is P38's own (c-step) on the same fixture, tied by (I-Py), 50× the band; one-sided, far from saturation.
- **(c0) By construction:** (S0) shows the schedule reaches the solver. A route that dropped `duty` would read 41 ones and fail (S0) before (T-qduty).

**Windows** (foreground; Bash `timeout: 660000`; `-s`; node IDs; **the harness call is the whole Bash command: bare, repo-relative, nothing chained or piped**; `pgrep` 0 before W1). Order: W0, W1, Wgen, Wval, Wapi, Wdoc.
- **W0** (complex, no mesh; 528's W0, which is 521b's W0 plus the s16 parse node, plus this item's; check-family's glob grows by one):
  `scripts/testing/run_and_log.sh WF-1-s18-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_thermal_box_quadratic_duty.py::test_wf1_s18_parse_and_refusals tests/package/test_wf1_yaml_thermal_box_quadratic.py::test_wf1_s16_parse_and_refusals tests/package/test_wf1_yaml_thermal_box.py::test_wf1_s15_parse_and_refusals tests/package/test_wf1_yaml_thermal_ball.py::test_wf1_s14_parse_and_refusals tests/package/test_wf1_yaml_capacitors.py::test_wf1_s13_parse_and_refusals tests/package/test_wf1_yaml_tune.py::test_wf1_s11_tune_parse_and_refusals tests/package/test_wf1_yaml_tune_f_human.py::test_wf1_s12_parse_and_refusals tests/package/test_wf1_check_family.py tests/package/test_wf1_yaml_f_human.py tests/package/test_wf1_yaml_f_human_thermal.py -v -s --tb=short'"`
  If item 539 has landed, insert `tests/package/test_wf1_yaml_thermal_box_probes.py::test_wf1_s17_parse_and_refusals` after the s16 node. If 528's landed parse node has a different ID, use the landed one and say so.
- **W1** (complex, heavy; predicted 120–200 s):
  `scripts/testing/run_and_log.sh WF-1-s18-W1-qduty "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_WF1_S18=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 590 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_thermal_box_quadratic_duty.py::test_wf1_s18_quadratic_duty_reproduces_p38 -s -v --tb=short'"`
- **Wgen / Wval / Wapi / Wdoc:** item 539's four commands with labels `WF-1-s18-Wgen` / `-Wval` / `-Wapi` / `-Wdoc`.
- **Not re-run (disclosed), valid only if no `src/` line moves:** every RG of 521b / 528. The item adds a file and a module, nothing on any route. `git diff --stat main -- src/` must read empty; if it does not, the item has left its brief, so revert and stop.

**Traps already paid for.** Item 528's list (this file) and 521a's / 521b's apply whole, in particular:
- **P2 comparands are P3** (`CMP_DEGREE`) with `QDEG_P` on `_rel_p`; P1's `_rel` changes the quantity (P28's trap);
- **the comparand is evaluated in chunks at dof coordinates** (`_exact_p`), never `interpolate(callable)`;
- **`_on` is P1-only** (`tests/package/test_wf15_study_transient.py:127–142` builds a `("Lagrange", 1)` space). At P2, assert the two P2 spaces' `tabulate_dof_coordinates` and local sizes equal, then compare owned arrays with `_theta_owned_rel`, or reuse 528's landed helper. Never mix meshes in UFL;
- **the period is `2.0 * HALF_T`**, from the imported constant;
- complex build: `np.real` after `_allsum`; no UFL ordering comparison;
- orphan `.c` after a kill: `pgrep` first, then `WF-1-s18-cache-sweep`; never raise a timeout or change `-n`;
- write ≈, never a tilde.

**Scope.** Covered: P38's duty-cycled transverse source on the m = 4 box at P2 (CN to t\* with the square-wave duty) from a tracked `thermal_box` file through `femem run`, equal to the Python Study at `L0_BAND` and reproducing (T-qduty) at `INTERIOR_BAND` to its printed digits, complex, `-n 2`. **Not claimed:**
- probes on this file (item 539's route at P1 only);
- P38's volume-mean and point records (printed by P38, not gated);
- another duty, period or weight; an h-rate; the real build;
- the box's EM SAR; a coil; any compliance statement.

**Negative result.** Report and stop; park on `attempt/WF-1-step18-<UTC>`, item 540 BLOCKED in the same commit.
- **(S0) red:** the duty does not reach the quadratic solve, which is a 528 defect. Open a known-issues entry; the review re-poses.
- **(T-qduty) `.6e` red with the band and (I-Py) green:** environment drift or an operation-order difference in 528's Q4. Print both forms. It goes to `record-reconciler` (1*), never re-registered in-slot.
- **A band red:** P38 does not reproduce through the file. Open a known-issues entry.
- **Status 124:** not a result; park with the measured time.

No `PREDICTED RED:` line.

**Status it can move.** `WF-1` stays ✅. On green the row gains: "**Step-18 record (§9 item 540, U66, provisional):** a `thermal_box` file with `source: {kind: transverse_quadratic}` and `thermal.duty` reproduces P38's (T-qduty) `2.543310e-04` at P2 (`INTERIOR_BAND`), schedule = step 2's, balance ≤ <e> (`E_BAND`), θ equal to the Python Study at ≤ <x> (`L0_BAND`), complex, `-n 2`; the unscheduled comparand misses by <p>; no `src/` change." §10 U66 DONE (provisional). `MAT-21` and `WF-15` do not move.

**Spawn line:** "Use the implementer agent to execute `WF-1` §9 item 540 (U66) per `docs/planning/chunks/WF-1.md` "§9 item 540" (read the brief and any review ruling beneath it first; the rulings hold). If item 528 is not DONE on `main`, stop and report. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. First read the landed `_study_box_quadratic_temperature`; if it does not pass `duty` as `source_schedule`, stop before any edit and report. Run W0, W1, Wgen, Wval, Wapi, Wdoc in order, in the foreground, every harness Bash call with `timeout: 660000`, pytest with `-s`, container `timeout -k 30`; before W1 run `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` and require 0. W0 / W1 run in the complex build with `FEM_EM_WF1_S18=1` on W1. No `src/` edit; never edit an existing test assertion; never import a `test_*` name. You are the executor: do not spawn agents, and never background a harness window. On a red, park on `attempt/WF-1-step18-*`, mark item 540 BLOCKED in the same commit, and report."

**For the review.**
1. **No `src/` edit, so no YAML re-gates.** This holds only if 528 lands as briefed; the first-action check enforces it. If the review wants 521b's re-gate list run anyway, add ≈ 13 min of windows (item 539's) and raise the predicted cost to ≈ 45.
2. **`.6e` equality on (T-qduty)**, as accepted for P36 (03:00 ruling, point 7).
3. **Placement after 528** in the queue; the gate makes an early slot skip at no cost.

### Review ruling 2026-10-18 18:00 (items 539 U65, 540 U66)

- **Both are provisional chain-U steps under rule (6)** (2 of 3 allowed for chain U before the 10-19 weekly). Each reproduces an existing record at its existing band. No band, tolerance, Done-when or target is set or moved.
- **Item 539, the refusal row.** Admitting `outputs.probes` on `thermal_box` files makes 521b's refusal row (`test_wf1_yaml_thermal_box.py:257–258`) false by design. This is a feature change, not a loosened gate. **The row's input moves to a still-refused form** (an empty `probes: []`, or another malformed probes value the parser refuses with field `outputs.probes`), in the same commit, disclosed in the journal. **W0 asserts that the new input still raises with field `outputs.probes`**; if no such input exists, the row is deleted and the deletion disclosed. No other existing assertion is edited.
- **Item 539 runs in the real build**, where U60's records were taken (`20261018T123416Z_WF-15-s6-W2-box.log:12`).
- **Ordering with item 528.** Item 528 sits above 539 in the queue. If 528 has landed, 539 adjusts its W0 and the (X) tables as its brief says. If it has not, 539 still runs, and 528's "(X) tables unedited" clause reads "unedited except for 539's disclosed rows".
- **Item 540 depends on item 528.** If 528 is not DONE on `main`, skip it at no attempt cost. Its first action checks that 528 landed as briefed (duty passed through as the source schedule); otherwise it stops and journals.
- Spot-checked by the review: `tests/package/test_wf1_yaml_thermal_box.py:82–85`, `:156–157`; `20261018T185143Z_WF-1-s15-W1a-steady.log:65–67`.

**Outcome, item 540 (2026-10-19 04:30 slot, take-next): DONE.** First action held: `_study_box_quadratic_temperature` passes `source_schedule=duty` (`src/fem_em_solver/study.py:1559`). W0 `20261019T094458Z_WF-1-s18-W0.log` 18 passed, 4 s (item 539's s17 parse node inserted): (K) spec diff empty, `dt_s` / period / duty by hex, `check` rc (0, 0); (X) refuses at `thermal.duty` without `transient`. W1 `20261019T094516Z_WF-1-s18-W1-qduty.log` 5 passed, Status 0, 121 s (predicted 120–200): (T-qduty) `2.543310e-04` = `P38_T_QDUTY`, within `INTERIOR_BAND`; (I-Py) 2.221e-14 ≤ `L0_BAND`; (S0) = `SCHEDULE`, multipliers equal to the Python Study's; (E-qduty) 2.288e-13 (= `P38_E`); (Pdeg) 116 025 / 116 025; (UT) 79 872 cells, bitwise-equal geometry; (I=) green; (c-step) 5.0333e-01 (= `P38_C_STEP`), > band. Wgen 46 files; Wval 306 rows (+5); Wapi 22 passed (the planted (c2) dead reference reads as designed); Wdoc dead 0. `git diff --stat main -- src/` empty. Disclosed: the P2 (I-Py) comparison uses 528's `_on_p` / `_owned_rel` in place of the P1-only `_theta_owned_rel`, which the trap list permits; `CELLS_M4 = 79872` is a literal from the brief, not an imported constant. Commit `af3d6764`.

### §9 item 545 — 2026-10-19 (chain U, U67, daily — provisional; chain step)

*Drafted by a read-only agent for the 2026-10-19 03:00 daily review; the review reviews it. Code read on `ddfde4aa`; nothing was run.*

**Licence.** §9 rule (6) (`docs/automation/daily-review.md:196–214`): a provisional next step of chain U, tagged "(daily, 2026-10-19 — provisional)", toward chain U's stated target, "each step reproduces a gated result or asserts a gated identity through the API, within its unmoved band" (item 521's licence, this file). **The target is a reproduction of an existing record at its existing band: U60's ball-Study probe readings (`WF-15` step 6, §9 item 523, `MAT-19`'s virtual-probe gate) at `P_BAND` / `B_BAND`, and θ_file = θ_Python at `L0_BAND`. No band, tolerance, record, Done-when or target is set or moved; every band is imported by name.** **Cap check for the review:** with items 539 (U65, DONE) and 540 (U66, OPEN), this is chain U's **third** provisional step since the 10-17 weekly, which is the cap. If the 10-19 weekly has already run when the review reads this, the count restarts from its ratification. It closes item 539's scope line "probes on a `thermal_ball` file (still refused)" (this file, item 539 "Scope") and the clause item 539 narrowed at `docs/planning/chunks/WF-15.md:590`. Templates: item 539 (the box half, this file) and item 500 (`WF-1` step 14, the ball files). No sheet, κ or ESR ruling is engaged.

**Feasibility (checked at drafting, source only).**
- **U60's ball record exists and is a real-build record** (`20261018T123401Z_WF-15-s6-W1-ball.log:50`, `PETSC_ARCH='linux-gnu-real64-32'`):
  - (T\*-ball) c0 / c1 / c2 rel `9.7501e-04` / `5.0491e-04` / `9.4257e-04` ≤ `P_BAND` (`:79–81`);
  - (A-ball) 3.007e-16 / 1.332e-15 ≤ `B_BAND` (`:78`);
  - (c-disp) `5.4208e-02`, comparand separation 5.4663e-02 (`:83`);
  - 19 489 cells (`:77`); the probe call took 2.8 s (`:76`), and the window 8 s (`:113`).
  - It is asserted by `tests/package/test_wf15_study_probe_thermal.py::test_wf15_s6_ball_probe_reproduces_mat19_gate` (`:186–229`). The fixture is `Study(ThermalBallGeometry(R, 1.0e-3), None, UniformHeatSource(_q_cf())).probe([WHOLE_BALL, *BALL_PROBES, DISP_PROBE], RobinCooling(K, H, T_INF), transient=TransientSpec(T_STAR / N_MAIN, N_MAIN, RHO_C))` (`:189–193`).
- **The ball file route can already express that fixture.** `_parse_thermal_ball_mapping` (`src/fem_em_solver/study_yaml.py:741–859`) admits:
  - `geometry.radius` / `resolution`;
  - `excitation.q_w_per_m3` (any finite value > 0, so `_q_cf()`'s repr);
  - `thermal.k_w_per_m_k` / `h_w_per_m2_k` / `t_inf_c`, with `degree` omitted (so `RobinCooling`'s default of 1, as in U60);
  - `thermal.transient.{dt_s, n_steps, rho_c_j_per_m3_k}`.

  The one refusal is `outputs`: `_unknown(osec, ("directory",), "outputs", source)` at `:854`. `Study.probe` already takes `transient` / `duty` (`src/fem_em_solver/study.py:1734–1756`, U60). **So the change is an admission plus a dispatch, with no physics or `study.py` edit.**

**Code state (`ddfde4aa`, after item 539 `7dd35a6f`).**
- `_parse_box_probes(value, source)` (`study_yaml.py:1080–1106`) is item 539's copy of the list-and-entry block. It has no geometry or `thermal:` check, and its field names are `outputs.probes`, `outputs.probes[i]`, `….centre` and `….radius`.
- `_run_thermal_ball` (`:1468–1524`) takes no `outputs` and returns `StudyRun(…, thermal=thermal)` with `probes=None`. `_run_thermal_box_probes` (`:1592–1618`) is the wrapper to copy.
- `run_study_file`'s dispatch heads (`:1706–1713`), in order: box with probes (`:1706`), box (`:1710`), ball (`:1712`). Docstring `:1691–1695`. The module docstring's ball bullet ends "``outputs`` admits ``directory`` only" (`:67`); the box bullet inherits it ("with the ``thermal_ball`` sections below", `:53–54`), so that clause has been stale for the box since 539.
- **The two rows that pin the refusal this item makes false** (rule-(c) edits, see "For the review" 1):
  - (a) 500's (X) row `({**a, "outputs": {"probes": [{"centre": [0.0, 0.0, 0.0], "radius": 0.01}]}}, "outputs.probes")` (`tests/package/test_wf1_yaml_thermal_ball.py:317–318`);
  - (b) item 539's (c-ball) block: `ball_f = _field(lambda: _study_from_mapping({**b, "outputs": {"probes": [ok]}}))` … `assert ball_f == "outputs.probes"` (`tests/package/test_wf1_yaml_thermal_box_probes.py:217–222`). Its module docstring names it at `:17–18`.
- The ball (X) row `outputs.sar` (`test_wf1_yaml_thermal_ball.py:319–320`) stays true, because R2 admits `probes` only.
- `femem check` on a ball file prints `family: thermal_ball` and `excitation: …` (`src/fem_em_solver/cli.py:171–176`). It reads no outputs, so `cli.py` is not touched. `test_wf1_check_family.py:32` globs every non-malformed `wf1/*.yaml`, so the glob grows by the two good files below. W0 of item 539 read 27; this item expects 29, or 30 if item 540 has landed.

- **Executor:** implementer. **Tier:** Wr smoke (real, no mpiexec); W0 / RG1 / W9 / W1 / RG-B / RG2 / RG3 / RG5 / RG7a / RG6a / RG6b standard; **RG4 heavy (`-n 8`)**; `-n 2` everywhere else.
- **Class:** chain step (solve). **Cost:** ≈ 20.0 measured (chain-solve class median, 18:00 restock) / **≈ 28 predicted min**. The base is item 539's measured 25 min, a like-sized admission (`7dd35a6f` committed at minute 25, `docs/testing/attempts.md:1715`). Then:
  - minus ≈ 6.6 min for the four box re-gates not re-run here (RG7b 187 s, RG8b 98 s, RG8a 74 s, RG-P 38 s; `20261019T033816Z` / `034903Z` / `034729Z` / `033444Z_WF-1-s17-*.log`);
  - plus ≈ 4.5 min for Wr, the second retargeted row and the duty leg's aliased imports;
  - plus ≈ 5 contingency.
- **Window prices (measured on item 539's windows unless named):**

  | Window | Price | Source |
  |---|---|---|
  | Wr | ≈ 5 s | U60's W0a, `20261018T123111Z_WF-15-s6-W0a-observe-omega.log` |
  | W0 | 6 s | `20261019T033322Z_WF-1-s17-W0.log:397` |
  | RG1 | 4 s | |
  | W9 | 15 s | |
  | W1 | predicted 30–90 s (below) | |
  | RG-B (item 539's W1) | 25 s | `20261019T033412Z_WF-1-s17-W1-probes.log:117` |
  | RG2 | 63 s | |
  | RG3 | 33 s | |
  | RG5 | 41 s | |
  | RG7a | 10 s | |
  | RG4 | 209 s | |
  | RG6a | 65 s | |
  | RG6b | 74 s | |
  | Wgen / Wval / Wapi / Wdoc | 1 / 3 / 20 / 2 s | |

  Windows total ≈ 10.6 min.
- **W1's prediction:** U60's probe call takes 2.8 s on 19 489 cells. The file run solves twice (the temperature block, then `Study.probe`'s re-solve), and the Python twin once more. The duty leg runs 500's case-A ball (151 548 cells), whose file run took 11.2 s (`20261019T034504Z_WF-1-s17-RG6a-ball-duty.log:106`); here that is two file solves plus a twin. Start-up is ≈ 5 s.
- **Park (rule (d)) if any of these is missed:**
  - code, Wr, W0 and RG1 green by **minute 25**;
  - W1 started by minute 28;
  - RG4 able to start by minute 45.
- **Depends on:** nothing open (U60 = item 523 and U65 = item 539 are on `main`). **Family:** `WF-1` step 19, attempt 1 of 4. Step 18 is item 540 (U66, OPEN); this item does not wait for it.
- **Item 540 (U66) may land first.** It adds no `src/` line (its brief, "Not re-run"). If it is on `main` when the slot starts, insert its `tests/package/test_wf1_yaml_thermal_box_quadratic_duty.py::test_wf1_s18_parse_and_refusals` into W0 after the s17 node, and say so. Its W1 is a box-only route the new lines never reach, so it is not re-run (disclosed). The order of the two items does not matter.

**The change (`src/fem_em_solver/study_yaml.py` only; insertions plus two disclosed edited lines, rule (c)).**
- **(R1) No new parser.** The ball path **calls** `_parse_box_probes` unedited. Its field names are the ones wanted, and calling it edits nothing on the box path. The alternative, a copied `_parse_ball_probes`, is "For the review" 2.
- **(R2) Edited line 1:** `:854` becomes `_unknown(osec, ("directory", "probes"), "outputs", source)`. Then one inserted pair after the `directory` lines (`:855–856`): `if "probes" in osec:  # `WF-1` step 19 (§9 item 545)` / `outputs["probes"] = _parse_box_probes(osec["probes"], source)`.
- **(R3) `_run_thermal_ball_probes(study, kind, cooling, outputs, path, out, digest, started, comm)`**, a new function after `_run_thermal_ball`. It is `_run_thermal_box_probes` (`:1592–1618`) **copied** with `_run_thermal_box` → `_run_thermal_ball`:
  - it calls `_run_thermal_ball(...)` unedited;
  - then `probes = study.probe(outputs["probes"], cooling.cooling, transient=cooling.transient, duty=cooling.duty)`, collective, every rank;
  - it writes `summary["probes"]` with the same four keys, re-stamps `elapsed_s`, and rewrites `summary.json` on rank 0;
  - then `Barrier`, `run.probes = probes`, and it returns `run`.

  Disclosed, as in item 539: `summary.json` is written twice; the second write is the record; the probe re-solves the temperature.
- **(R4)** In `run_study_file`, inserted **before** `:1712`: `if isinstance(study.geometry, ThermalBallGeometry) and "probes" in outputs:  # §9 item 545` / `return _run_thermal_ball_probes(...)`. Box files return at `:1706` / `:1710` before reaching it. A ball file without probes reaches item 500's dispatch at `:1712` unchanged. An EM file evaluates the same `isinstance` predicate `:1712` already evaluates, false.
- **Edited line 2 (docstring):** `:67`'s "``outputs`` admits ``directory`` only" becomes "``outputs`` admits ``directory`` and ``probes`` (`WF-1` steps 17 / 19)". The re-wrap may move the next line; disclose the hunk. Insert one clause in `run_study_file`'s docstring after `:1695`.
- **Never:**
  - edit `_parse_box_probes`, `_parse_probes`, `_parse_outputs`, `_run_thermal_ball`, `_run_thermal_box`, `_run_thermal_box_probes` or `_parse_thermal_box_mapping`;
  - touch `study.py` or `cli.py`;
  - add to `__all__`;
  - refuse or admit anything else in the ball's `outputs` (the `outputs.sar` row stays green unedited).

**Tracked inputs (`tests/package/wf1/`; each opens `# WF-1 step 19 (§9 item 545): …`, then `femem: 1`; floats as Python `repr` digits taken from Wr, never retyped).**
- `thermal_ball_m19_probes.yaml`: U60's ball fixture.
  - `geometry: {kind: thermal_ball, radius: <repr R>, resolution: 0.001}`;
  - `excitation: {kind: uniform_heat_source, q_w_per_m3: <repr _q_cf()>}`;
  - `thermal: {k_w_per_m_k: <repr K>, h_w_per_m2_k: <repr H>, t_inf_c: <repr T_INF>, transient: {dt_s: <repr T_STAR / N_MAIN>, n_steps: 40, rho_c_j_per_m3_k: <repr RHO_C>}}`, with `degree` omitted;
  - `outputs: {probes: [...]}`: U60's five, in U60's order. These are `WHOLE_BALL` (centre 0, radius `<repr 2.0 * R>`), then `A_PROBE` balls at the three `CENTRES`, then `A_PROBE` at `DISPLACED` (`test_wf15_study_probe_thermal.py:81–83`). CENTRES[2] and DISPLACED are NumPy products, so their components come from Wr's `repr(float(x))`.
- `thermal_ball_duty_a_probes.yaml`: 500's `thermal_ball_duty_a.yaml` verbatim plus `outputs: {probes: [{centre: [0.0, 0.0, 0.0], radius: 0.2}, {centre: [0.0, 0.0, 0.0], radius: 0.02}]}`. The first probe is the whole ball (2R with MAT-9's R = 0.1). This is the duty pass-through file. **No record exists for it, and none is claimed.**
- `malformed_thermal_ball_probe_radius_zero.yaml` → `outputs.probes[1].radius` (the m19 file with probe 1's radius 0.0).

**Test module: new `tests/package/test_wf1_yaml_thermal_ball_probes.py`.** It takes item 539's module shape (`tests/package/test_wf1_yaml_thermal_box_probes.py`; `_knob`, `_p0`, `_say`, `_field`, `_probe_hex`, `_check_summary` with `family == "thermal_ball"`). Knob `FEM_EM_WF1_S19=1`; when set, assert `not IS_COMPLEX` and size 2, never skip on a mismatch. The first docstring paragraph backticks only `WF-1`, `WF-15` and `MAT-19`.
- **Imports (module-level names only; never a `test_*` name):**
  - `WHOLE_BALL`, `BALL_PROBES`, `DISP_PROBE`, `_domain_mean` from `tests/package/test_wf15_study_probe_thermal.py:81–83`, `:128–133`;
  - `A_PROBE`, `CENTRES`, `DISPLACED`, `_ball_mean`, `_q_cf` from `tests/validation/test_mat19_virtual_temperature_probe.py:76–97`;
  - `R`, `K`, `H`, `T_INF`, `P_BAND` (`:59`), `B_BAND` (`:61`) from `tests/validation/test_mat10_sar_to_heat_sphere.py`;
  - `N_MAIN`, `RHO_C`, `T_STAR` from `tests/validation/test_mat12_transient_sar_heating_sphere.py:68–71`;
  - `L0_BAND` (`tests/validation/test_mat5_em_coupled_sphere.py:104`);
  - for the duty leg, **aliased**: `R as R9`, `K as K9`, `H as H9`, `Q as Q9`, `RHO_C as RHO_C9` from `tests/validation/test_mat9_transient_conduction_sphere.py:41–46`; `CASES` (`test_mat17_duty_cycle_closed_form.py:44`); `T_END_A`, `_n` (`test_mat17_duty_cycle_fem.py`); `RES` (`tests/package/test_wf15_study_transient.py:84`);
  - `_same_mesh`, `_on`, `_theta_owned_rel` (`test_wf15_study_transient.py:117–152`); `LOGS`, `WF1`, `_cmp`, `_fresh_dir`, `_mapping`, `_rc_all`, `_spec_diff` (`tests/package/test_wf1_study_yaml.py`); `IS_COMPLEX`.
- **Record constants (transcribed with log lines, `20261018T123401Z_WF-15-s6-W1-ball.log`):**
  - `U60_T_BALL = ("9.7501e-04", "5.0491e-04", "9.4257e-04")` (`:79–81`);
  - `U60_C_DISP = 5.4208e-02` (`:83`);
  - `U60_A_BALL = 3.007e-16` (`:78`, printed beside, never compared);
  - `U60_BALL_CELLS = 19489` (`:77`, printed beside).
- **Nodes:**
  - `test_wf1_s19_parse_and_refusals` (no mesh; every `build` patched to raise; 0 calls);
  - `test_wf1_s19_ball_probes_reproduce_u60` (W1);
  - `test_wf1_s19_duty_probes_equal_python` (W1, the pass-through).
- **Registry:** under `# WF-1 step 19 (§9 item 545)`, the parse node under `[names.load_study]` (`tests/package/validation_envelope.toml:345`) and `[names.StudyConfigError]` (`:327`), and the two solve nodes under `[names.run_study_file]` (`:376`) and `[names.ThermalBallGeometry]` (`:232`). **The §9 standing rule applies:** Wgen / Wval / Wapi / Wdoc run in this slot, and the regenerated `docs/guide/api/` and `docs/guide/validation.md` are committed. Wval's row count is the last committed page's plus 6 (3 nodes × 2 names, item 539's arithmetic); print it.

**Anchors (asserted; bands imported by name).**
- **(K) Parse, no mesh:**
  - `_spec_diff(load_study(m19), Study(ThermalBallGeometry(R, 1.0e-3), None, UniformHeatSource(_q_cf()))) == []`;
  - `_spec_diff(load_study(duty_probes), load_study(WF1 / "thermal_ball_duty_a.yaml")) == []`;
  - the parsed probes equal `[WHOLE_BALL, *BALL_PROBES, DISP_PROBE]` (m19) and the two duty probes by `float.hex` (`_probe_hex`);
  - the cooling, transient and duty through `_cmp` by hex: m19 `RobinCooling(K, H, T_INF)`, `TransientSpec(T_STAR / N_MAIN, N_MAIN, RHO_C)`, duty `None`; duty file `RobinCooling(K9, H9, T_INF)`, `TransientSpec(150.0, _n(T_END_A, 150.0), RHO_C9)`, `DutyCycle(*CASES["A"])` (500's want table, `test_wf1_yaml_thermal_ball.py:238–239`);
  - `kind == "thermal_ball"`, `set(outputs) == {"probes"}`;
  - `cli.main(["check", f])` rc (0, 0) via `_rc_all`, one `family: thermal_ball` line on rank 0;
  - **walled ball, parse only:** `_study_from_mapping({**_mapping(WF1 / "thermal_ball_wall.yaml"), "outputs": {"probes": [ok]}})` parses, and `_load`'s outputs carry one probe (see "Scope").
- **(X) Refusals naming their field**, all on m19's mapping:
  - the malformed file through `load`, `run` and `check`: rc (2, 2), `out` absent by `MPI.LOR`;
  - `probes: []` → `outputs.probes`; `probes` as a mapping → `outputs.probes`; `[1]` → `outputs.probes[0]`; an extra key → `outputs.probes[0].<key>`; a missing `radius` → `outputs.probes[0].radius`; a two-element `centre` → `outputs.probes[0].centre`;
  - `outputs: {probes: [ok], sar: {}}` → `outputs.sar`;
  - `{probes: [ok]}` → `None`.
- **(UT) / (I-Py), U60's fixture:**
  - the file run's probe mesh equals the Python Study's (`_same_mesh`; cells printed beside `U60_BALL_CELLS`). The Python twin is U60's exact call (`test_wf15_study_probe_thermal.py:190–193`);
  - temperatures, not θ (item 539's disclosure (2)): `_theta_owned_rel(_on(m_py, run.probes.thermal.temperature), py.thermal.temperature, comm)` ≤ **`L0_BAND`** (predicted 0 to ≈ 5e-14; item 539 read 1.184e-15);
  - `run.thermal` vs `run.probes.thermal` (the re-solve) at `L0_BAND`;
  - the file's and the Python Study's probe values, printed only.
- **(T\*-ball) U60's records through the file:**
  - for each of the three `CENTRES`, `abs(P_{1+j} / _ball_mean(norm(c_j), A_PROBE, T_STAR, _q_cf()) − 1)` ≤ **`P_BAND`**, and `f"{r:.4e}" == U60_T_BALL[j]`. The digits are asserted: same build, same mesh, same `Study.probe` call, tied by (I-Py). This is item 539's (R-box) rule.
  - **(A-ball):** the whole-ball probe gives \|P₀ − ∫θ/V_h\| / \|∫θ/V_h\| ≤ **`B_BAND`** and \|V₀/V_h − 1\| ≤ `B_BAND` (`_domain_mean`), printed beside `U60_A_BALL`.
- **(U′) The summary is its own result:** `summary["probes"]["values_k"] == [float(v) for v in run.probes.values]` exactly; JSON on disk equals the in-memory summary; `family == "thermal_ball"`; `run.result is None` (item 539's `_check_summary`, copied with the family changed).
- **(I-du) The duty pass-through, an identity only.** The duty file is compared against `Study(ThermalBallGeometry(R9, RES), None, UniformHeatSource(Q9)).probe(DUTY_PROBES, RobinCooling(K9, H9, T_INF), transient=TransientSpec(150.0, _n(T_END_A, 150.0), RHO_C9), duty=DutyCycle(*CASES["A"]))`:
  - temperatures ≤ **`L0_BAND`** (`_theta_owned_rel` through `_on`);
  - `source_multipliers` `array_equal` between the two runs and against the summary;
  - the whole-ball identity at **`B_BAND`** (`MAT-19`'s (P8), a gated identity);
  - the inner-probe values printed only. No comparand is posed.

**Negative controls (rule (e)).**
- **(c-disp) The displaced probe through a tracked file — *asserted, backed*.** The m19 file's fifth reading (`DISPLACED`, 1 mm radially inward of c₂), against c₂'s nominal comparand, misses: `> P_BAND`, printed `.4e` beside `U60_C_DISP`.
  - It is U60's own comparison (`…W1-ball.log:83`) on the same mesh in the same build, so the digits are predicted equal; equality is printed, not asserted.
  - **Ceiling:** the comparand separation itself is 5.4663e-02 (`:83`), so the miss cannot exceed ≈ 5.5 × the band on this fixture. A blind solver (θ ≡ 0) saturates at 1.0. One-sided; far from both.
- **(c-sar) The refusal that must not move — *asserted*.** The ball `outputs.sar` row (`test_wf1_yaml_thermal_ball.py:319–320`, W0, unedited) and the sphere-without-`thermal:` probes refusal (`test_wf1_yaml_probes.py:56`, RG1) still raise. Together they show that R2 admits `probes` only, and only on thermal-only files.
- **(c0) By construction:** R4's `isinstance … and "probes" in outputs` is false for every ball file without probes. RG6a / RG6b are 500's own nodes, whose (I=) pins the summary keys with no `probes` (`20261019T034504Z_WF-1-s17-RG6a-ball-duty.log:108`); they show the no-probe ball route is unchanged.

**Windows.**
- Foreground; Bash `timeout: 660000`; pytest `-s`; node IDs, never `-k`. **The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped.**
- `pgrep` means `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3`, and it must read 0 before W1, RG2, RG3, RG5, RG4 and RG6a.
- Wr runs **before** the YAML files are written. W0 runs only after every edit (rule (i)).
- **Order:** Wr, code + files + module + the two retargeted rows, W0, RG1, W9, W1, RG-B, RG2, RG3, RG5, RG7a, RG4, RG6a, RG6b, Wgen, Wval, Wapi, Wdoc.

The windows:
- **Wr** (real, smoke; it prints the reprs the files carry and asserts nothing). The `-c` program has no single quote (item 523's trap), only escaped double quotes:
  `scripts/testing/run_and_log.sh WF-1-s19-Wr-reprs "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 10 60 python3 -c \"from tests.validation.test_mat19_virtual_temperature_probe import A_PROBE, CENTRES, DISPLACED, _q_cf; from tests.validation.test_mat12_transient_sar_heating_sphere import N_MAIN, RHO_C, T_STAR; from tests.validation.test_mat10_sar_to_heat_sphere import H, K, R, T_INF; print(repr(float(_q_cf())), repr(float(T_STAR / N_MAIN)), N_MAIN, repr(float(RHO_C)), repr(float(R)), repr(float(2.0 * R)), repr(float(K)), repr(float(H)), repr(float(T_INF)), repr(float(A_PROBE))); print([[repr(float(x)) for x in c] for c in CENTRES]); print([repr(float(x)) for x in DISPLACED])\"'"`
  - If the guard denies it, write the files from W0's (K) mismatch print (the node prints `repr` of each expected value on failure), then re-run W0. That costs one 6 s window.
- **W0** (complex, no mesh; item 539's W0 plus the new node; if item 540 has landed, insert its s18 parse node after the s17 node):
  `scripts/testing/run_and_log.sh WF-1-s19-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_thermal_ball_probes.py::test_wf1_s19_parse_and_refusals tests/package/test_wf1_yaml_thermal_box_probes.py::test_wf1_s17_parse_and_refusals tests/package/test_wf1_yaml_thermal_box.py::test_wf1_s15_parse_and_refusals tests/package/test_wf1_yaml_thermal_box_quadratic.py::test_wf1_s16_parse_and_refusals tests/package/test_wf1_yaml_thermal_ball.py::test_wf1_s14_parse_and_refusals tests/package/test_wf1_yaml_capacitors.py::test_wf1_s13_parse_and_refusals tests/package/test_wf1_yaml_tune.py::test_wf1_s11_tune_parse_and_refusals tests/package/test_wf1_yaml_tune_f_human.py::test_wf1_s12_parse_and_refusals tests/package/test_wf1_check_family.py tests/package/test_wf1_yaml_f_human.py tests/package/test_wf1_yaml_f_human_thermal.py -v -s --tb=short'"`
  - W0 must show three things:
    - both retargeted rows raising at `outputs.probes`;
    - (K) green in the complex build, which is predicted because `_q_cf` is pure-Python complex arithmetic (`src/fem_em_solver/post/sar.py:381–409`) and no PETSc scalar is involved;
    - check-family's good-file count printed.
- **RG1** (complex): item 539's RG1 verbatim (`20261019T033338Z_WF-1-s17-RG1.log:12`), label `WF-1-s19-RG1`.
- **W9** (real installed): item 539's W9 verbatim (`20261019T033342Z_WF-1-s17-W9-installed-femem.log:12`), label `WF-1-s19-W9-installed-femem`.
- **W1** (real, standard; predicted 30–90 s; the parse node again in the real build, then the two solve nodes):
  `scripts/testing/run_and_log.sh WF-1-s19-W1-ball-probes "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace FEM_EM_WF1_S19=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_thermal_ball_probes.py::test_wf1_s19_parse_and_refusals tests/package/test_wf1_yaml_thermal_ball_probes.py::test_wf1_s19_ball_probes_reproduce_u60 tests/package/test_wf1_yaml_thermal_ball_probes.py::test_wf1_s19_duty_probes_equal_python -s -v --tb=short'"`
- **RG-B** (real; the box-probes route shares `_parse_box_probes` and sits on the dispatch head above R4; item 539's W1 verbatim, `20261019T033412Z_WF-1-s17-W1-probes.log:12`, 25 s):
  `scripts/testing/run_and_log.sh WF-1-s19-RGB-box-probes "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace FEM_EM_WF1_S17=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_thermal_box_probes.py::test_wf1_s17_box_probes_reproduce_u60 tests/package/test_wf1_yaml_thermal_box_probes.py::test_wf1_s17_transient_probes_equal_python -s -v --tb=short'"`
- **RG2 / RG3 / RG5 / RG7a** (complex; every EM file now evaluates R4's predicate): item 539's RG2, RG3, RG5 and RG7a commands verbatim (`20261019T033527Z` / `033637Z` / `033715Z` / `033801Z_WF-1-s17-*.log:12`), labels `WF-1-s19-RG2` / `-RG3` / `-RG5` / `-RG7a-box-steady`.
- **RG4 — required, alone** (complex, heavy, `-n 8`; R4 runs ahead of every F-human file's dispatch, so the 18:00 ruling's reason for item 539's RG4 holds; `20261019T034129Z_WF-1-s17-RG4.log:12`, 209 s):
  `scripts/testing/run_and_log.sh WF-1-s19-RG4 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_WF1_S12=1 timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_tune_f_human.py -s -v --tb=short'"`
- **RG6a / RG6b** (real; the no-probe ball route, the one R2 edits; `20261019T034504Z` 65 s / `20261019T034609Z` 74 s):
  `scripts/testing/run_and_log.sh WF-1-s19-RG6a-ball-duty "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace FEM_EM_WF1_S14=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_thermal_ball.py::test_wf1_s14_parse_and_refusals tests/package/test_wf1_yaml_thermal_ball.py::test_wf1_s14_duty_reproduces_mat17 -s -v --tb=short'"`
  `scripts/testing/run_and_log.sh WF-1-s19-RG6b-ball-wall "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace FEM_EM_WF1_S14=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_thermal_ball.py::test_wf1_s14_wall_reproduces_mat18_2b -s -v --tb=short'"`
- **Wgen / Wval / Wapi / Wdoc:** item 539's four commands verbatim, labels `WF-1-s19-Wgen` / `-Wval` / `-Wapi` / `-Wdoc`. Wdoc must read `0 dead`.
- **Not re-run (disclosed, by construction; valid only if `git diff --stat main -- src/` shows `study_yaml.py` alone):**
  - RG-P (sphere probes; `_parse_probes` and the sphere route are untouched);
  - RG7b, RG8a, RG8b and item 540's W1 (box files return at `:1706` / `:1710`, above every new line);
  - U60's windows (no `study.py` line moves);
  - `WF-13`'s estimate module (no `cli.py` line moves).

**Traps already paid for.** Item 539's list (this file) and items 521 and 500 apply whole. In particular:
- never `float.hex` across builds; `femem` prints on rank 0 only; `-k a or b` splits argv; `_load`'s 4-tuple is load-bearing; `__all__` does not grow;
- **UFL cannot mix meshes**: compare through `_on` (P1; the ball's thermal degree is 1), or in NumPy;
- **(I-Py) compares temperatures, not θ**: `_theta_owned_rel` normalises by max\|T − T∞\| (item 539's disclosure (2));
- **the probe is an owned-cell integral**: every rank enters `Study.probe`, including a rank that owns no cell of a 2 mm probe ball (item 227's deadlock). R3 calls it unconditionally;
- **the real build here, the YAML gates in theirs**: the digits in `U60_T_BALL` are real-build digits;
- orphan `.c` after a kill: `pgrep` first, then at 0 the rubric's sweep as `WF-1-s19-cache-sweep`; never raise a timeout or change `-n`;
- rule (i): an edit to `study_yaml.py` after W0 re-runs W0 and every RG.

**New here:**
- **Two constant sets share names.** MAT-10 / MAT-12 (`R` 0.01, `H` 100.0, U60's) and MAT-9 (`R` 0.1, `H` 10.0, 500's duty ball) both export `R`, `K`, `H`, `RHO_C`. Import MAT-9's under the `…9` aliases only; a shadowed `R` silently builds the wrong ball.
- **YAML 1.1 reads an exponent float without a dot (`1e-05`) as a string.** If a Wr repr has that form, write it with a dot (`1.0e-05`), which is the same double; `_number` would refuse the string at its field anyway. The m19 values are predicted positional, but check each.
- **`CENTRES[2]` / `DISPLACED` are NumPy float64**: take `repr(float(x))`, never `repr(x)` (NumPy 2 prints `np.float64(…)`).
- **Two existing rows move, not one** (question 1). Edit them only as the review ruling says.
- Write ≈, never a tilde.

**Scope.** Covered: `outputs.probes` on `thermal_ball` files through `femem run` / `run_study_file`, real build, `-n 2`. Two legs:
- U60's `MAT-19` fixture (CN N = 40 to t\*): U60's three (T\*-ball) readings and the whole-ball identity reproduced at their bands and digits, equal to the Python Study at `L0_BAND`;
- the duty pass-through on 500's case-A ball, as an identity only.

**Not claimed:**
- any solve on a walled ball with probes: it is admitted by the parser and parse-checked only (question 3);
- any comparand on the duty ball's probes;
- `MAT-19`'s own-mesh 1.7315e-03 or its 40-time series;
- scalar point evaluation (P37);
- the complex build for the probe legs;
- the installed console script running the new files (W9 re-gates the installed `check`);
- `--estimate` on a ball probes file;
- a probe's own volume or self-heating; any compliance statement.

**Negative result.** Report and stop. Park on `attempt/WF-1-step19-<UTC>`, mark item 545 BLOCKED with its unblock condition in the same commit (rule (d)), and journal it in attempts.md.
- **(K) / (X) red:** a parse defect. Fix it in-slot before W1 only. A (K) hex red in the complex W0 alone, green in W1's real parse node, means `_q_cf` is build-dependent. Record it in known-issues, assert (K)'s q in the real build only, and say so.
- **(I-Py) red with (UT) green:** two equal specs on one build give a temperature beyond `L0_BAND`. That is a solve-floor finding; the band is unmoved.
- **(T\*-ball) `.4e` red with its band green:** environment drift since U60's window (`20261018T123401Z`). Send it to `record-reconciler` (1\*); never re-register in-slot.
- **A band red:** U60 no longer reproduces through the file. Open a known-issues entry.
- **(c-sar) red, or RG6a / RG6b red:** R2 / R4 changed more than the ball-probes path. Revert `src/` and stop.
- **Any other RG red:** revert `src/` and stop.
- **Status 124:** not a result; park with the measured time.

No `PREDICTED RED:` line: a red anywhere stops the slot.

**Status it can move.** `WF-1` stays ✅ (rule (j)). On green its §7 row gains:

"**Step-19 record (§9 item 545, U67, provisional):** `outputs.probes` on a `thermal_ball` file drives `Study.probe` through `run_study_file`, real build, `-n 2`: U60's three ball-probe readings reproduce to their printed digits (`P_BAND`), the whole-ball probe identity ≤ <a> (`B_BAND`), θ_file = θ_Python ≤ <x> (`L0_BAND`), the duty pass-through equal to the Python Study ≤ <y>; the displaced probe misses by <p>; `outputs.probes` is admitted on both thermal-only families and stays refused without `thermal:` on a sphere; `src`: `study_yaml.py`, insertions plus two disclosed edited lines."

Also on green:
- §10 U67 DONE (provisional).
- Item 539's ball scope line resolves, and `docs/planning/chunks/WF-15.md:590`'s "probes from YAML on thermal-only files" clause resolves; add a dated note under it.
- `WF-15`, `MAT-19` and `MAT-21` do not move.

**Commits.**
1. `WF-1 step 19 (§9 item 545, U67): outputs.probes on thermal_ball files`. It carries:
   - `src/fem_em_solver/study_yaml.py`;
   - the three `tests/package/wf1/` files and the new module;
   - the two retargeted rows (`test_wf1_yaml_thermal_ball.py`, `test_wf1_yaml_thermal_box_probes.py`, including the latter's (c-ball) docstring and print text);
   - `tests/package/validation_envelope.toml`, the regenerated `docs/guide/api/` and `docs/guide/validation.md`;
   - every `WF-1-s19-*` harness log and its test-results rows;
   - `PROJECT_PLAN.md` (§7 row, §9 item 545 DONE line, §10 U67), this file's outcome paragraph and the `WF-15.md:590` note.

   The message names both edited `study_yaml.py` lines and both moved test rows.
2. `attempts: WF-1 step 19 §9 item 545 journal` (`docs/testing/attempts.md` only).

**Spawn line:** "Use the implementer agent to execute `WF-1` §9 item 545 (U67) per `docs/planning/chunks/WF-1.md` "§9 item 545" (read the brief and any review ruling beneath it first; the rulings hold). The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Run Wr first and write the YAML floats from its reprs. Make the `study_yaml.py` change: R2's `:854` line and pair, R3, R4 and the `:67` docstring line; nothing else in `src/`. Then run W0, RG1, W9, W1, RG-B, RG2, RG3, RG5, RG7a, RG4, RG6a, RG6b, Wgen, Wval, Wapi, Wdoc in order, in the foreground, every harness Bash call with `timeout: 660000`, pytest with `-s`, container `timeout -k 30`. Before W1, RG2, RG3, RG5, RG4 and RG6a, run `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` and require 0. W1 runs in the real build with `FEM_EM_WF1_S19=1`; every re-gate runs in its own recorded build. If item 540 is on `main`, add its s18 parse node to W0 after the s17 node. Edit `tests/package/test_wf1_yaml_thermal_ball.py:317–318` and `tests/package/test_wf1_yaml_thermal_box_probes.py:217–222` only as the review ruling on this item says, and never any other existing test assertion. Never import a `test_*` name; never compare across builds by `float.hex`; import MAT-9's constants under the `…9` aliases only. You are the executor: do not spawn agents, and never background a harness window. Park if code + Wr + W0 + RG1 are not green by minute 25. On a red, park on `attempt/WF-1-step19-*`, mark item 545 BLOCKED in the same commit, and report."

**For the review.**
1. **Two pinned rows must change.** 500's (X) row (`test_wf1_yaml_thermal_ball.py:317–318`) and item 539's (c-ball) (`test_wf1_yaml_thermal_box_probes.py:217–222`) both assert that a ball file with one valid probe is refused at `outputs.probes`, and this item's letter admits it. The drafter proposes the 18:00 ruling's pattern for 521b's row: **move each input to `outputs: {probes: []}`**, which still raises at `outputs.probes` (`_parse_box_probes`'s empty-list check, `study_yaml.py:1090–1091`), and have W0 show both. The (c-ball) print text and docstring (`:17–18`) follow. Without a ruling, the item is not runnable.
2. **R1 calls `_parse_box_probes` rather than copying it.** No line of the box path moves, but a future edit to that function would reach the ball. The alternative is a copied `_parse_ball_probes` (≈ 27 lines, item 539's copy-not-factor rule), at no extra window cost.
3. **The walled ball is admitted, with a parse check only.** The alternative is a refusal of `outputs.probes` beside `geometry.wall` (no record exists for a probe spanning gel and wall). That is a narrower admission and adds one (X) row.
4. **The re-gate set** drops item 539's RG-P, RG7b, RG8a and RG8b (≈ 6.6 min), by construction. If the review wants item 539's full list, add them with their item 539 commands and raise the prediction to ≈ 35 min.
5. **Chain U's provisional count** is three with this item (cap). Queue it only if the count holds when the review reads this.

### Review ruling 2026-10-19 03:00 (item 545, U67)

- **Rule (6) count:** U65, U66 and U67 are chain U's three provisional steps before the 10-19 weekly (the cap). The weekly ratifies or reverts all three. No band, tolerance, Done-when or target is set or moved.
- **Q1, the two pinned rows: as proposed.** Move both inputs to `outputs: {probes: []}`: 500's row at `tests/package/test_wf1_yaml_thermal_ball.py:317–318`, and item 539's (c-ball) at `tests/package/test_wf1_yaml_thermal_box_probes.py:217–222`, with its print text and docstring. W0 asserts that both still raise with field `outputs.probes`. These are the only existing assertions this item edits, and both edits are disclosed in the journal. This is the 18:00 pattern for 521b's row.
- **Q2: R1 stands.** The ball path calls `_parse_box_probes` unedited. Calling an existing parser is neither a copy nor a factor, and one parser is the smaller surface.
- **Q3: refuse the walled ball** (the narrower admission). No record exists for a probe in a gel-plus-wall ball, and an item never admits an ungated configuration. `outputs.probes` beside `geometry.wall` raises `StudyConfigError` with field `outputs.probes`, and one new (X) row asserts it. Drop the brief's "walled ball, parse only" check, and its "Scope" line becomes "a walled ball with probes is refused".
- **Q4: the reduced re-gate set stands**, on one condition. Before commit, `git diff main -- src/fem_em_solver/study_yaml.py` must show no changed line inside `_parse_box_probes`, `_run_thermal_box_probes` or `_run_thermal_box`. Otherwise RG-P, RG7b, RG8a and RG8b run with their item 539 commands.
- **Ordering:** independent of item 540. If 540 has landed, W0 adds its s18 parse node, as the brief says.

**Outcome 2026-10-19 (06:00 slot) — DONE (provisional; U67).** All windows green, run in order, foreground. Wr (`20261019T110057Z`, 2 s) printed the reprs the files carry (q `2002.088514415672`, dt `3.4833333333333334`, CENTRES[2] `0.0037527767497325675`, DISPLACED `0.0031754264805429417`; every value positional). W0 (`110312Z`, complex, 5 s, 19 passed) with item 540's s18 parse node after s17: (K) both new files hex-equal (5 and 2 probes, cooling / transient / duty diffs `[]`), q's hex green in the complex build; (X) nine rows green including the walled-ball refusal at `outputs.probes` (ruling Q3); both retargeted rows (`test_wf1_yaml_thermal_ball.py:317`, `test_wf1_yaml_thermal_box_probes.py` (c-ball)) raise at `outputs.probes` on `probes: []`; check-family reads **30** good files. RG1 (`110339Z`, 4 s, 13 passed); W9 (`110348Z`, 14 s, 1 passed). **W1** (`110412Z_WF-1-s19-W1-ball-probes.log`, real `linux-gnu-real64-32`, `-n 2`, 60 s, 4 passed): (UT) 19 489 / 19 489 cells (U60 19 489), geometry equal; (I-Py) 2.973e-13 ≤ `L0_BAND`, temperature block vs probe re-solve 2.973e-13; (A-ball) 4.511e-16 / 1.332e-15 ≤ `B_BAND` (U60 3.007e-16); (T\*-ball) 9.7501e-04 / 5.0491e-04 / 9.4257e-04, digits equal to U60's; (c-disp) 5.4208e-02, digits equal to U60's, > `P_BAND`; (U′) green; (I-du) 151 548 cells, 8.146e-15 ≤ `L0_BAND`, multipliers bit-equal (95 entries) to the Python Study and the summary, whole-ball identity 8.476e-15 / 6.994e-15 ≤ `B_BAND`. Re-gates: RG-B 24 s, RG2 63 s, RG3 31 s, RG5 43 s, RG7a 10 s, RG4 (`-n 8`) 216 s, RG6a 64 s, RG6b 73 s, all green; `pgrep` read 0 before each required window. Wgen 1 s, Wval 3 s (**312** rows = 306 + 6), Wapi 20 s (22 passed), Wdoc 2 s (0 dead). **Q4 check:** `git diff main -- src/fem_em_solver/study_yaml.py` shows no changed line inside `_parse_box_probes`, `_run_thermal_box_probes` or `_run_thermal_box`. **Disclosed edits:** the ball `_unknown(osec, …)` line; the `:67` docstring line, re-wrapped onto a second line; one `run_study_file` docstring line (the item-539 sentence's closing period became a semicolon to carry the ball clause — a third edited line, docstring only, made before W0). The walled-ball refusal is three inserted lines (ruling Q3); the "Scope" line reads "a walled ball with probes is refused". Not re-run, by construction (`git diff --stat main -- src/` shows `study_yaml.py` alone): RG-P, RG7b, RG8a, RG8b, item 540's W1, U60's windows, `WF-13`'s estimate module.

## §7 row moved verbatim 2026-10-19 (rotate_plan_archive.py chunks)

 `femem run`: the MRI study CLI and YAML config | ✅ (item 256, U5): `femem run` builds Studies equal to the Python-built ones by `float.hex` and reproduces `MAT4_DEG1_RECORD` and `P1` at their bands; steps 2–14 are records at their imported bands — the installed console script, `thermal:`, `conductor:`, the SAR / B₁⁺ / probe outputs, the F-human two-port schema, `tune:`, `capacitor_sheets` (U57) and the thermal-only `thermal_ball` files (U58) — all ratified by the weekly. *History: `docs/planning/chunks/WF-1.md`.* **Step-15 record (§9 item 521b, U59 part 2, provisional, 2026-10-18):** a `thermal_box` geometry (`half_widths_m`, `divisions`), `uniform_heat_source` and optional `thermal.transient` / `thermal.duty` drive `Study.temperature` through `run_study_file`, complex, `-n 2`: `MAT-21` step 1's m = 4 `3.4599e-03` and P26's `9.4014e-04` reproduce (`INTERIOR_BAND` / `E_BAND`), θ equal to the Python Study at ≤ 4.792e-14 (`L0_BAND`); the wrong-wall file misses by 3.8048e-01; bad box keys raise `StudyConfigError` naming their field; `femem check` prints `family: thermal_box`; `src`: `study_yaml.py` + `cli.py`, insertions only (`20261018T185143Z_WF-1-s15-W1a-steady.log`, 10 s; `20261018T185201Z_WF-1-s15-W1b-duty.log`, 190 s). **Step-16 record (§9 item 528, U64, provisional, 2026-10-19):** `excitation.source: {kind: transverse_quadratic}` on a `thermal_box` file (`TransverseQuadraticHeatSource`) reproduces P36's steady `7.219878e-05` and CN `2.694340e-04` at P2, complex, `-n 2` (`INTERIOR_BAND`, Gauss and balance at their bands), θ equal to the Python Study at ≤ 1.213e-14; dropping the source misses by 5.1471e-01 (`20261019T003608Z_WF-1-s16-W1a-steady.log`, 78 s; `20261019T003738Z_WF-1-s16-W1b-cn.log`, 98 s). **Step-17 record (§9 item 539, U65, provisional, 2026-10-19):** `outputs.probes` on a `thermal_box` file drives `Study.probe` through `run_study_file`, real build, `-n 2`: U60's five box-probe readings reproduce to their printed digits (`INTERIOR_BAND`), the whole-box probe identity ≤ 2.495e-13 / 2.509e-13 (`B_BAND`), θ_file = θ_Python ≤ 1.184e-15 (`L0_BAND`), the duty pass-through equal to the Python Study ≤ 6.912e-14; the wrong-wall file misses by 3.5712e-01; `outputs.probes` stays refused on `thermal_ball` and without `thermal:`; `src`: `study_yaml.py`, insertions plus one disclosed edited line (`20261019T033412Z_WF-1-s17-W1-probes.log`, 25 s). **Step-18 record (§9 item 540, U66, provisional, 2026-10-19):** `excitation.source: {kind: transverse_quadratic}` plus `thermal.duty` on a `thermal_box` file reproduces P38's duty-cycled (T-qduty) `2.543310e-04` at P2 to its printed digits (`INTERIOR_BAND`), complex, `-n 2`: schedule (S0) `11111000011110000111100001111000011110000` reaches the solve, 116 025 P2 dofs on 79 872 cells, CN balance 2.288e-13 (`E_BAND`), θ_file = θ_Python ≤ 2.221e-14 (`L0_BAND`); the unscheduled comparand misses by 5.0333e-01; the file without `transient` raises at `thermal.duty`; no `src/` edit (`20261019T094516Z_WF-1-s18-W1-qduty.log`, 121 s). **Step-19 record (§9 item 545, U67, provisional, 2026-10-19):** `outputs.probes` on a `thermal_ball` file drives `Study.probe` through `run_study_file`, real build, `-n 2`: U60's three ball-probe readings reproduce to their printed digits 9.7501e-04 / 5.0491e-04 / 9.4257e-04 (`P_BAND`) on 19 489 cells, the whole-ball probe identity ≤ 4.511e-16 / 1.332e-15 (`B_BAND`), θ_file = θ_Python ≤ 2.973e-13 (`L0_BAND`), the duty pass-through equal to the Python Study ≤ 8.146e-15 (multipliers bit-equal, whole-ball identity 8.476e-15); the displaced probe misses by 5.4208e-02; `outputs.probes` is admitted on both thermal-only families, refused beside `geometry.wall` (review ruling 2026-10-19 03:00 Q3), and stays refused without `thermal:` on a sphere; `src`: `study_yaml.py`, insertions plus disclosed edited lines (the ball `outputs` `_unknown` line, the module-docstring `:67` line re-wrapped, one `run_study_file` docstring line); two refusal rows retargeted to `probes: []` (`20261019T110412Z_WF-1-s19-W1-ball-probes.log`, 60 s). 

### §9 item 564 — 2026-10-20 (§10 chain U70, 9a; lane 2 feature admission)

*Drafted by a read-only agent for the 2026-10-20 03:00 daily review; the review reviews it. Code read on `20a2ba1f`; nothing was run.*

**Licence.** §10 chain U, U70, as the 2026-10-19 weekly enumerated it (`PROJECT_PLAN.md`, "Live chains and lanes", *Chain U*): "`outputs.probes` on an EM file: point SAR / \|E\| probes on item 192's four-port F-small file through `point_sar` under the owner-rule default (3b-i), width-invariant at `-n 1 / 2 / 4` to 1e-9 (the `POST-8` 3a band), θ_file = θ_Python at `L0_BAND`; a displaced probe misses. A disclosed `study_yaml.py` admission with every YAML gate re-run. Serial on nothing." Lane 2 by the same block. **No band, tolerance, record, Done-when or target is set or moved; every band is imported by name:** `WIDTH_RTOL` = 1e-9 (`tests/validation/test_post8_step3a_averaged_row1.py:88`, the `POST-8` 3a band), `L0_BAND` (`tests/validation/test_mat5_em_coupled_sphere.py:104`), `LEG_D0_REPRODUCTION_BAND` (`tests/validation/test_port_birdcage_four_port.py:165`), `EXACT_IDENTITY_RTOL` (`tests/validation/test_birdcage_sar_mass_averaged.py:163`). This is `WF-1` **step 20** (step 19 is item 545). Templates: item 344 (`outputs.sar` / `b1_plus` on the same fixture) and item 545 (the latest admission, its re-gate set).

**Feasibility (source only).**
- `post.sar.point_sar(e_real, e_imag, points, *, sigma, rho, comm, owner="lowest_original_cell")` (`src/fem_em_solver/post/sar.py:332–378`) is `POST-8` 3b-i's owner-rule default (§9 item 526): σ\|E\|²/(2ρ) from the split real / imaginary fields; a point no rank evaluates raises.
- The port Study's field-keeping solve `Study._solve_port_fields` (`src/fem_em_solver/study.py:1952–2033`) is cached on `study.fields_result` and returns `{port_id: TimeHarmonicFields}`, whose `e_real` / `e_imag` (`core/time_harmonic.py:134–135`) are exactly the DG split `point_sar` and the owner rule are written for. `run_study_file`'s four-port branch already reaches it through `_sar_summary` / `_b1_plus_summary` (`src/fem_em_solver/study_yaml.py:1795–1798`).
- The phantom is the cylinder r < `phantom_radius` = 0.03 m, \|z\| < `phantom_height`/2 = 0.04 m, centred at the origin (`FourPortBirdcageGeometry`, `study.py`; `io/mesh.py:3686–3692`); σ is `study.material.sigma` (`materials/phantom.py`), ρ is `study.density_kg_m3` (the value `Study.sar` uses).
- Today `outputs.probes` on a `four_port_birdcage` raises at `outputs.probes` ("probes are defined only on a sphere_in_box", `_parse_probes`, `study_yaml.py:1050–1085`, dispatched at `:1141–1142`). **No existing test pins that refusal for a four-port file** (grep of `tests/package` at drafting); the F-human row (`tests/package/test_wf1_yaml_f_human.py:235–236`) and the sphere / box / ball rows stay true.

**Executor:** implementer. **Tier:** standard by prediction, heavy by ceiling (item 344's container `timeout -k 30 400` on this fixture), complex build, `-n 2` for W0 / W1; the width legs at `-n 1` and `-n 4`. **Class:** chain step (solve), **≈ 16.7 measured / ≈ 60 predicted min**: 15 fixed + ≈ 12 `src` + ≈ 8 file and module + windows ≈ 24 (below) + ≈ 1 spare. **Park (rule (d)) if code + W0 are not green by minute 25, or W1 has not started by minute 28.** **Depends on:** nothing open. Items 562 (`study.py`) and 563 (docs / tests) share no `src` file with this one; if either lands first, rebase, keep both sides of `validation_envelope.toml` and the generated pages, and re-run Wgen / Wval / Wapi / Wdoc. **Family:** `WF-1` step 20, attempt 1 of 4.

**The change (`src/fem_em_solver/study_yaml.py` only; insertions plus one disclosed edited line, rule (c)).**
- **(P1) `_parse_em_probes(value, source) -> tuple`**, a new function after `_parse_box_probes`: `outputs.probes` on a `four_port_birdcage` is a mapping `{points: [[x, y, z], …]}` (the `b1_plus` shape; EM probes are points, thermal probes are balls). Unknown keys ⇒ `outputs.probes.<key>`; missing / empty `points` ⇒ `outputs.probes.points`; a bad entry ⇒ `outputs.probes.points[i]` (`_triple`). Returns a tuple of float triples.
- **(P2) Edited line:** `_parse_outputs`'s `out["probes"] = _parse_probes(...)` (`:1142`) becomes a two-way dispatch — `_parse_em_probes(sec["probes"], source)` when `geometry_kind == "four_port_birdcage"`, else `_parse_probes(...)` unchanged. `_parse_probes`, `_parse_box_probes` and the ball / box `outputs` lines are not edited; the F-human and sphere refusals keep their fields.
- **(P3) The in-phantom check, before any mesh:** inserted in `_parse_mapping` (`:1172–1265`) where the parsed geometry and outputs are both in hand: each point must satisfy `hypot(x, y) < geometry.phantom_radius` and `abs(z) < geometry.phantom_height / 2` (strict; σ and ρ below are the phantom's), else `StudyConfigError` at `outputs.probes.points[i]` ("outside the phantom").
- **(P4) `_em_probe_summary(study, points, comm) -> dict`**, a new function after `_b1_plus_summary`: `res = study._solve_port_fields(comm)` (cached; collective); for each port id, in `res.fields` order, `sar = point_sar(f.e_real, f.e_imag, pts, sigma=float(study.material.sigma), rho=float(study.density_kg_m3), comm=res.comm)` (owner default, **never** the complex field twice), `e_abs = sqrt(2 ρ sar / σ)`. Returns `{"points", "owner": "lowest_original_cell", "sigma_s_per_m", "rho_kg_m3", "values": {port: {"sar_w_per_kg": [...], "e_abs_v_per_m": [...]}}}` (floats). Function-local imports from `.post.sar`; `src` never imports `tests`.
- **(P5) Inserted in `run_study_file`'s four-port branch** after the `b1_plus` lines (`:1797–1798`): `if "probes" in outputs:  # WF-1 step 20 (§9 item 564)` / `summary["probes"] = _em_probe_summary(study, outputs["probes"], comm)`. `StudyRun.probes` stays `None` on EM files (the summary is the carrier, as for `sar`). The module docstring and `run_study_file`'s docstring each gain one clause (disclosed hunks).
- **Never:** touch `study.py`, `cli.py`, `post/`, `__all__`; admit `probes` on `f_human_birdcage`, `capacitor_sheets` or a copper / hole four-port beyond what `_parse_outputs` already lets through (the copper F-small file reaches the same branch — reachable, not measured, not claimed).

**Tracked input:** `tests/package/wf1/f_small_em_probes.yaml` = `f_small_four_port.yaml` verbatim (item 192's fixture: degree 1, 64 MHz, `preset: f_small`, `lumped_port_drive`) with the first comment line `# WF-1 step 20 (§9 item 564): …` and `outputs: {probes: {points: [[0.02, 0.0, 0.0], [0.0, 0.02, 0.0], [-0.02, 0.0, 0.01], [0.0, -0.02, -0.01], [0.014, 0.014, 0.02]]}}` — five off-axis points inside the phantom at r = 0.02 m (the axis is avoided: E is smallest there, so a 1e-9 relative statement is hardest to read). Decimal literals; the module carries the same literals, so (K) is hex-equal by construction.

**Test module: new `tests/package/test_wf1_yaml_em_probes.py`** (item 545's module shape; first docstring paragraph backticks only `WF-1`, `POST-8`, `WF-16`). Knob `FEM_EM_WF1_S20=1` for the solve nodes; when set, assert `IS_COMPLEX` and the width the node is for, never skip on a mismatch.
- **Imports (names only; never a `test_*` name):** `point_sar` (`fem_em_solver.post.sar`); `_port_study` (`tests/package/test_wf16_study_sar.py:87`); `_record`, `_rel` (`tests/package/test_wf16_study_four_port.py:110`, `:115`); `ANS4_EXACT_STRIP_CELL_RECORD` (`tests/validation/test_ops61_exact_strip_ports.py:64`); the four bands above; `_bitwise_equal_across_ranks` (`tests/post/test_point_owner_rule.py:84`); `LOGS`, `WF1`, `_mapping`, `_fresh_dir`, `_rc_all`, `_spec_diff` (`tests/package/test_wf1_study_yaml.py`); `_field` / `_study_from_mapping` as item 545's module reaches them.
- **Nodes:** `test_wf1_s20_parse_and_refusals` (no mesh; every `build` patched to raise, 0 calls), `test_wf1_s20_em_probes_n2` (W1), `test_wf1_s20_em_probes_width` (Wn1, Wn4).
- **Registry (lines added only):** the parse node under `[names.load_study]` and `[names.StudyConfigError]`, the two solve nodes under `[names.run_study_file]`, with `# WF-1 step 20 (§9 item 564)`. **The §9 standing rule applies:** Wgen / Wval / Wapi / Wdoc run in this slot and the regenerated `docs/guide/api/` and `docs/guide/validation.md` are committed; Wval's row count = the committed page's + 4 (print it).

**Anchors (asserted; bands imported by name).**
- **(K) parse, no mesh:** `_spec_diff(load_study(EM_PROBES), load_study(WF1 / "f_small_four_port.yaml")) == []`; the parsed points equal the module's five triples by `float.hex`; `kind == "four_port_birdcage"`; `cli.main(["check", f])` rc (0, 0).
- **(X) refusals naming their field** (on the file's mapping, no mesh): `probes: {}` ⇒ `outputs.probes.points`; `points: []` ⇒ `outputs.probes.points`; a two-element point ⇒ `outputs.probes.points[0]`; an extra key ⇒ `outputs.probes.<key>`; `[0.0, 0.0, 0.05]` (\|z\| > 0.04) and `[0.031, 0.0, 0.0]` ⇒ `outputs.probes.points[i]`; a ball-form list `[{centre, radius}]` on the four-port ⇒ `outputs.probes`; the F-human row (`test_wf1_yaml_f_human.py:235–236`) and the sphere-without-`thermal:` row (`test_wf1_yaml_probes.py:56`) unedited and green in W0 / RG1.
- **(R) same fixture** (W1): `summary["cells"] == ANS4_EXACT_STRIP_CELL_RECORD`; `_rel(summary S, _record(1)) ≤ LEG_D0_REPRODUCTION_BAND` (item 344 read 9.557e-11).
- **(I-Py) θ_file = θ_Python — the weekly's `L0_BAND` target** (W1): an independent `_port_study(degree=1)` Study, `_solve_port_fields(comm)`, then `point_sar(f.e_real, f.e_imag, pts, sigma=…, rho=…, comm=…)` per port; every summary SAR value against it, relative, ≤ `L0_BAND` (predicted ≤ 1e-12: item 256's in-process (I) on this fixture read 9.110e-16).
- **(U′) the summary is the API on its own fields:** the summary values equal `point_sar` called on `run.study.fields_result` exactly; `_bitwise_equal_across_ranks` on the value array; the JSON on disk equals the in-memory summary; `e_abs² σ / (2ρ)` against `sar_w_per_kg` ≤ `EXACT_IDENTITY_RTOL`.
- **(W) width invariance — the weekly's 1e-9 target:** W1 (rank 0, after an `allgather` check) writes `logs/wf1/s20/em-probes-n2.npz` (points, port order, the 4 × 5 SAR and \|E\| arrays, cells); Wn1 (`-n 1`) and Wn4 (`-n 4`) run the file and assert, per value, `|v_n / v_2 − 1| ≤ WIDTH_RTOL`, cells equal, port order equal. The npz absent ⇒ the width node fails with a message (never re-created there).
- **Printed only:** the legacy `owner=None` values at each width beside the owner values (the `POST-8` motivation; never asserted here); per-point \|E\| and SAR per 1 V drive; the width spreads in units of the band.

**Negative control (rule (e)).**
- **(c-disp) a displaced probe misses — *asserted*:** in W1, point 0 displaced 1 mm radially, `(0.021, 0.0, 0.0)`, through `point_sar` on `run.study.fields_result`, against the summary's point-0 value per port: relative difference > `WIDTH_RTOL` (and so > `L0_BAND`). **First measurement on this fixture — predicted 1e-2 to 3e-1** (\|E\| grows ≈ linearly with r near the axis, so δSAR/SAR ≈ 2δr/r ≈ 0.1). The bar is ≥ 1e7× below the prediction and a relative difference has no ceiling here; a value ≤ 1e-6 means the probe ignores its coordinates — stop.
- **(c-ref) the refusals that must not move — *asserted*:** the F-human and sphere probe rows above, unedited.

**Windows.** Foreground; **the harness call is the whole Bash command: bare, repo-relative, nothing chained or piped**; Bash `timeout: 660000`; pytest `-s`; node IDs, never `-k`. `pgrep` (`docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3`) must read 0 before W1, Wn1, Wn4, RG-SB, RG-F, RG-CuSB and RG4. W0 runs only after every edit (rule (i)). **Order:** code + file + module, W0, RG1, W9, W1, Wn1, Wn4, RG-P, RG-SB, RG-F, RG-Cu, RG-CuSB, RG2, RG3, RG5, RG4, Wgen, Wval, Wapi, Wdoc.

| Window | Price | Source |
|---|---|---|
| W0 | ≈ 6 s | item 545 |
| RG1 / W9 | 4 s / 15 s | item 545 |
| W1 (`-n 2`) | predicted 110–170 s: the file run's mesh + S solve + field solve ≈ 70 s (item 344: 106.9 s with 16 SAR calls), plus the Python twin's mesh + field solve ≈ 45 s | `20261009T005219Z_WF-1-s5-W1-sar-b1.log.gz` |
| Wn1 (`-n 1`) | predicted 90–220 s | none measured at `-n 1` on this fixture |
| Wn4 (`-n 4`) | predicted 50–120 s | |
| RG-P sphere probes | 38 s | `20261019T033444Z_WF-1-s17-RGP-sphere-probes.log` |
| RG-SB sar / b1 (the same branch) | 119 s | item 344 |
| RG-F F-small reproduction | 166 s | `20261018T100919Z_ARCH-1-U16a-W3-regate-wf1.log` |
| RG-Cu / RG-CuSB copper | 94 s / 151 s | `20261009T004108Z_WF-1-s4-W1-copper.log.gz`, `20261012T140300Z_WF-1-s8-W1.log` |
| RG2 tune / RG3 sphere / RG5 capacitors | 63 / 31 / 43 s | item 545 |
| RG4 F-human tune (`-n 8`) | 216 s | item 545 |
| Wgen / Wval / Wapi / Wdoc | 1 / 3 / 20 / 2 s | item 545 |

Windows ≈ 24 min.
```
scripts/testing/run_and_log.sh WF-1-s20-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_em_probes.py::test_wf1_s20_parse_and_refusals tests/package/test_wf1_yaml_thermal_ball_probes.py::test_wf1_s19_parse_and_refusals tests/package/test_wf1_yaml_thermal_box_quadratic_duty.py::test_wf1_s18_parse_and_refusals tests/package/test_wf1_yaml_thermal_box_probes.py::test_wf1_s17_parse_and_refusals tests/package/test_wf1_yaml_thermal_box.py::test_wf1_s15_parse_and_refusals tests/package/test_wf1_yaml_thermal_box_quadratic.py::test_wf1_s16_parse_and_refusals tests/package/test_wf1_yaml_thermal_ball.py::test_wf1_s14_parse_and_refusals tests/package/test_wf1_yaml_capacitors.py::test_wf1_s13_parse_and_refusals tests/package/test_wf1_yaml_tune.py::test_wf1_s11_tune_parse_and_refusals tests/package/test_wf1_yaml_tune_f_human.py::test_wf1_s12_parse_and_refusals tests/package/test_wf1_check_family.py tests/package/test_wf1_yaml_f_human.py tests/package/test_wf1_yaml_f_human_thermal.py -v -s --tb=short'"
scripts/testing/run_and_log.sh WF-1-s20-RG1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_study_yaml.py::test_wf1_good_files_build_the_python_studies tests/package/test_wf1_study_yaml.py::test_wf1_malformed_files_name_their_field tests/package/test_wf1_yaml_thermal.py::test_wf1_thermal_yaml_parse_and_refusals tests/package/test_wf1_yaml_copper.py::test_wf1_copper_yaml_parse_and_refusals tests/package/test_wf1_yaml_sar_b1.py::test_wf1_sar_b1_yaml_parse_and_refusals tests/package/test_wf1_yaml_probes.py::test_wf1_probes_yaml_parse_and_refusals tests/package/test_wf1_yaml_f_human.py::test_wf1_f_human_yaml_parse_and_refusals tests/package/test_wf1_yaml_copper_sar_b1.py::test_wf1_s8_parse_and_refusal_parity tests/package/test_wf1_yaml_f_human_thermal.py::test_wf1_fht_parse_and_refusals -v -s --tb=short'"
scripts/testing/run_and_log.sh WF-1-s20-W9-installed-femem "docker compose exec -T fem-em-solver bash -lc 'timeout -k 10 50 bash /workspace/scripts/package/wf13_setup.sh real && mkdir -p /tmp/wf13/run/tests/package/wf1 && cp /workspace/tests/package/test_wf1_installed_cli.py /tmp/wf13/run/tests/package/ && cp /workspace/tests/package/wf1/*.yaml /tmp/wf13/run/tests/package/wf1/ && cd /tmp/wf13/run && env -u PYTHONPATH WF13_INSTALLED_CHECK=1 timeout -k 10 90 /tmp/wf13/venv/bin/python -m pytest -p no:cacheprovider tests/package/test_wf1_installed_cli.py -v -s --tb=short'"
scripts/testing/run_and_log.sh WF-1-s20-W1-em-probes-n2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_WF1_S20=1 timeout -k 30 400 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_em_probes.py::test_wf1_s20_em_probes_n2 -s -v --tb=short'"
scripts/testing/run_and_log.sh WF-1-s20-Wn1-em-probes "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_WF1_S20=1 timeout -k 30 400 mpiexec -n 1 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_em_probes.py::test_wf1_s20_em_probes_width -s -v --tb=short'"
scripts/testing/run_and_log.sh WF-1-s20-Wn4-em-probes "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_WF1_S20=1 timeout -k 30 400 mpiexec -n 4 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_em_probes.py::test_wf1_s20_em_probes_width -s -v --tb=short'"
scripts/testing/run_and_log.sh WF-1-s20-RGP-sphere-probes "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_probes.py -s -v --tb=short'"
scripts/testing/run_and_log.sh WF-1-s20-RGSB-sar-b1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 400 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_sar_b1.py -s -v --tb=short'"
scripts/testing/run_and_log.sh WF-1-s20-RGF-f-small "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 400 mpiexec -n 2 python3 -m pytest tests/environment tests/package/test_wf1_study_yaml.py::test_wf1_cli_f_small_reproduces_wf16 -s -v --tb=short'"
scripts/testing/run_and_log.sh WF-1-s20-RGCu-copper "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 300 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_copper.py -s -v --tb=short'"
scripts/testing/run_and_log.sh WF-1-s20-RGCuSB-copper-sar-b1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 400 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_copper_sar_b1.py -v -s --tb=short'"
scripts/testing/run_and_log.sh WF-1-s20-RG2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 300 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_tune.py -v -s --tb=short'"
scripts/testing/run_and_log.sh WF-1-s20-RG3 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_study_yaml.py::test_wf1_cli_sphere_reproduces_wf15 -v -s --tb=short'"
scripts/testing/run_and_log.sh WF-1-s20-RG5 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 300 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_capacitors.py -s -v --tb=short'"
scripts/testing/run_and_log.sh WF-1-s20-RG4 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_WF1_S12=1 timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_tune_f_human.py -s -v --tb=short'"
```
Wgen / Wval / Wapi / Wdoc: item 545's four commands verbatim (logs `20261019T111510Z` / `111512Z` / `111516Z` / `111536Z_WF-1-s19-*.log:12`), labels `WF-1-s20-Wgen` / `-Wval` / `-Wapi` / `-Wdoc`; Wdoc reads `0 dead`.
- **Not re-run, by construction (valid only if `git diff --stat main -- src/` shows `study_yaml.py` alone and `git diff main -- src/fem_em_solver/study_yaml.py` shows no changed line inside `_parse_probes`, `_parse_box_probes`, `_run_thermal_box*` or `_run_thermal_ball*` — the 2026-10-19 03:00 Q4 condition):** the thermal-only solve gates RG-B, RG6a / RG6b, RG7a / RG7b, RG8a / RG8b, item 540's W1 and item 545's W1 — those files return at the thermal-only dispatch heads above the four-port branch, and their parse nodes run in W0. If the condition fails, they run with item 545's / 539's commands.

**Traps already paid for.**
- **Never hand `point_sar` the complex field twice** — it takes `np.real` of each argument and silently halves SAR (`tests/validation/test_birdcage_sar_map.py:295–302`); pass `f.e_real`, `f.e_imag`.
- **Point evaluation only through `point_sar` / `evaluate_vector_field_parallel`;** `points` identical on every rank (the parsed list is); every rank enters `_em_probe_summary` (collective); never `f.eval(points, …)`.
- **A bitwise claim names its batch** (rubric, `POST-8` item 537): (U′) compares the summary against `point_sar` on the **same** five-point batch, never a re-batched or single-point call; (W) is a band, not a byte claim.
- **One field solve per Study:** `_solve_port_fields` caches on `study.fields_result`; never clear it, never build a second Study except the declared Python twin.
- **The npz lives under the gitignored `logs/`;** rank 0 writes after an `allgather`, `Barrier` before return; the width node reads it on every rank. A stale npz from an earlier attempt is overwritten only by W1.
- **YAML 1.1 reads `1e-05` without a dot as a string** — the five points are plain decimals.
- `-k a or b` splits argv; `-k 30` on every container `timeout`; never pipe pytest; never import a `test_*` name; `femem` prints on rank 0 only.
- Status 124 or a kill: `pgrep` first; at 0 only, the rubric's orphan-stem sweep labelled `WF-1-s20-cache-sweep`. Never raise a timeout or change `-n`. Rule (i): a `study_yaml.py` edit after W0 re-runs W0 and every RG. Write ≈, never a tilde.

**Scope.** Covered: `outputs.probes: {points}` on the four-port F-small file (item 192's fixture, degree 1, 64 MHz, `-n 1 / 2 / 4`): point SAR and \|E\| per 1 V port drive through `point_sar` under the owner-rule default, equal to an independent Python Study at `L0_BAND` and width-invariant at `WIDTH_RTOL`. **Not claimed:** probes on F-human, copper or capacitor files (copper F-small reachable, unmeasured); degree 2; absolute or compliance SAR (per-1 V values); a peak search; points outside the phantom (refused); the installed console script running the file (W9 re-gates `check`); `POST-8`'s census closure (a new site born on the default; the `POST-8` row does not move).

**What a negative result means.** Report and stop; park on `attempt/WF-1-step20-<UTC>`, mark item 564 BLOCKED in the same commit (rule (d)), journal.
- **(K) / (X) red:** parser defect; fix before W1 only.
- **(I-Py) red with (R) green:** two equal Studies at one width disagree pointwise beyond `L0_BAND` — a solve-floor finding; record it, never re-band.
- **(W) red at `-n 1` or `-n 4` with (I-Py) green:** the owner rule is width-invariant but the solve is not, at these points, to 1e-9 — print the per-value spread and each point's \|E\| / max\|E\|, open a known-issues entry naming the points; **never** move the band or the points to pass.
- **(c-disp) not missing:** blind instrument; stop.
- **Any RG red:** revert `src/` and stop. **Status 124:** not a result; park with the measured time.

No `PREDICTED RED:` line: a red anywhere stops the slot.

**Status it can move.** `WF-1` stays ✅ (rule (j)). On green its state line's list gains "point SAR / \|E\| probes on the four-port file (U70)" and the history a **step-20 record**: "`outputs.probes: {points}` on item 192's F-small file drives `post.sar.point_sar` (owner rule) through `run_study_file`: θ_file = θ_Python ≤ <x> (`L0_BAND`), width spread `-n 1` / `-n 4` vs `-n 2` ≤ <w> (`WIDTH_RTOL`), S tie <s>; the 1 mm displaced probe misses by <p>; `src`: `study_yaml.py`, insertions plus one disclosed edited line." §10 U70 DONE; the Phase 9 example "a `mri:` probes-on-EM-file example after U70" (10-19 weekly, step 4) becomes openable.

**Commits.**
1. `WF-1 step 20 (§9 item 564, U70): point SAR / |E| probes on the four-port YAML file` — `src/fem_em_solver/study_yaml.py`; the new YAML file and module; `validation_envelope.toml`, the regenerated `docs/guide/api/` and `docs/guide/validation.md`; every `WF-1-s20-*` log and its test-results rows; `PROJECT_PLAN.md` (§7 `WF-1` state line, §9 item 564 DONE line, §10 U70); this file's outcome paragraph. The message names the edited `_parse_outputs` line.
2. `attempts: WF-1 step 20 §9 item 564 journal` (`docs/testing/attempts.md` only).

**Spawn line:** "Use the implementer agent to execute `WF-1` step 20, §9 item 564 (U70), per `docs/planning/chunks/WF-1.md` "§9 item 564" (read the brief and any review ruling beneath it first; the rulings hold). The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. `src/` change: `study_yaml.py` only — `_parse_em_probes`, the one edited `_parse_outputs` dispatch line, the in-phantom check in `_parse_mapping`, `_em_probe_summary`, the inserted pair in the four-port branch of `run_study_file`; no `study.py` / `cli.py` / `post/` edit. Write the YAML file and the module whole, then run W0, RG1, W9, W1, Wn1, Wn4, RG-P, RG-SB, RG-F, RG-Cu, RG-CuSB, RG2, RG3, RG5, RG4, Wgen, Wval, Wapi, Wdoc in order, foreground, every harness Bash call with `timeout: 660000`, pytest with `-s`, container `timeout -k 30`; `pgrep -c python3` must read 0 before W1, Wn1, Wn4, RG-SB, RG-F, RG-CuSB and RG4. W1 before Wn1 / Wn4 (they read its npz). Pass `f.e_real` / `f.e_imag` to `point_sar`, never the complex field twice. Never import a `test_*` name; never edit an existing test assertion. You are the executor: do not spawn agents, never background a harness window. Park if code + W0 are not green by minute 25 or W1 has not started by minute 28. On a red, park on `attempt/WF-1-step20-*`, mark item 564 BLOCKED in the same commit, and report."

**For the review.**
1. **"Every YAML gate re-run"** is read as every gate the diff can reach (the parse nodes of all families in W0 / RG1, the installed `check`, every four-port / sphere / F-human solve gate) plus item 545's Q4 condition for the thermal-only solve gates. The literal reading adds ≈ 12 min of thermal-only windows (≈ 72 predicted) and would need a split into 564a (code + W0–Wn4 + four-port re-gates) and 564b (the rest).
2. **The schema shape.** `probes` on a four-port file is `{points: [...]}` (the `b1_plus` form), while sphere / box / ball files keep `[{centre, radius}, …]` — one key, two shapes by family, each refused in the other's shape. A separate key (`point_probes`) would avoid that but is not the weekly's `outputs.probes`.
3. **The in-phantom refusal** is the drafter's: σ and ρ in `point_sar` are the phantom's, so a point in air would print a meaningless SAR. Dropping it would admit air points with \|E\| only.

### Review ruling 2026-10-20 03:00 (item 564)

Binds over the brief. Accepted as provisional execution choices under rule (6), disclosed in the commit, the weekly ratifies: on four-port EM files `outputs.probes` takes the `{points: [...]}` shape (thermal files keep the ball list); points outside the phantom are refused (σ and ρ are the phantom's). "Every YAML gate" is read as every gate the diff can reach plus the 10-19 Q4 git-diff condition for the thermal-only solve gates. The width legs at `-n 1 / 4` are the first width measurement on this fixture: a red there is a known-issues finding, reported, never a band or probe move.

**Outcome 2026-10-20 (07:30 CDT slot, 12:30–12:58 UTC) — DONE (provisional; U70), attempt 1 of 4.** All 19 windows green, run in order, foreground; `pgrep` read 0 before W1, Wn1, Wn4, RG-SB, RG-F, RG-CuSB and RG4. W0 (`20261020T123357Z`, 20 passed, 4 s): (K) spec diff `[]` vs `f_small_four_port.yaml`, the five points hex-equal, `kind == "four_port_birdcage"`, `check` rc (0, 0); (X) nine rows green (`probes: {}` and `points: []` ⇒ `outputs.probes.points`; a two-element point ⇒ `points[0]`; an extra key ⇒ `outputs.probes.mass`; `[0, 0, 0.05]` ⇒ `points[1]`; `[0.031, 0, 0]` ⇒ `points[0]`; the ball list ⇒ `outputs.probes`; capacitor sheets ⇒ `outputs.probes`; the good mapping ⇒ none); 0 build calls; the F-human and sphere probe refusals unedited and green (W0 / RG1). RG1 13 passed; W9 1 passed. **W1** (`20261020T123448Z_WF-1-s20-W1-em-probes-n2.log`, complex, `-n 2`, 106 s; run_study_file 68.2 s, the Python twin 34.3 s): (R) cells 117 206 = `ANS4_EXACT_STRIP_CELL_RECORD`, S tie 9.558e-11 ≤ `LEG_D0_REPRODUCTION_BAND`; (U′) the summary equals `point_sar` on `run.study.fields_result` exactly, bitwise across ranks, JSON = memory, \|E\|² σ/(2ρ) vs SAR 2.213e-16 ≤ `EXACT_IDENTITY_RTOL`; **(I-Py) 0.000e+00** ≤ `L0_BAND` (predicted ≤ 1e-12); **(c-disp) asserted**: the 1 mm displaced point 0 misses by 6.3865e-02 / 2.1232e-01 / 8.3482e-02 / 5.6184e-02 (P1…P4), > `WIDTH_RTOL`, inside the predicted 1e-2 to 3e-1. Per-1 V point SAR on the driven port's near point ≈ 2.6e-05 to 4.9e-05 W/kg (printed only). The legacy `owner=None` values equal the owner-rule values to all printed digits at every width (printed only). **Wn1** (`20261020T123649Z`, 78 s) max spread 2.025e-12 and **Wn4** (`20261020T123818Z`, 66 s) 2.414e-12, i.e. ≤ 2.5e-03 × `WIDTH_RTOL`; cells and port order equal. Re-gates: RG-P 9 passed (35 s), RG-SB 13 passed (116 s), RG-F 13 passed (161 s), RG-Cu 11 passed (89 s), RG-CuSB 11 passed (136 s), RG2 10 passed (56 s), RG3 5 passed (29 s), RG5 10 passed (37 s), RG4 `-n 8` 9 passed (199 s). Wgen 1 s; Wval 3 s (319 nodes; the page's table rows 405 → 409, + 4 as predicted); Wapi 22 passed; Wdoc 0 dead. **Q4 condition held:** `git diff --stat main -- src/` shows `study_yaml.py` alone (+73 / −1), and the one removed line is `_parse_outputs`'s `out["probes"] = _parse_probes(...)` dispatch; nothing inside `_parse_probes`, `_parse_box_probes`, `_run_thermal_box*` or `_run_thermal_ball*` changed, so the thermal-only solve gates were not re-run. **Disclosed:** (1) the edited dispatch line now reads `_parse_em_probes(...) if geometry_kind == "four_port_birdcage" else _parse_probes(...)`; (2) the in-phantom block in `_parse_mapping` also refuses `probes` on a `capacitor_sheets` four-port at `outputs.probes`, so a capacitor file is not admitted beyond what `_parse_outputs` let through before (the brief's "Never" line); (3) `import math` was added for `math.hypot`; (4) the module docstring gains one bullet clause and `run_study_file`'s docstring one sentence. The copper F-small file reaches the same branch: reachable, not measured, not claimed.

### §9 item 572 — 2026-10-20 (§10 U71, chain step (rule (6), provisional); lane 2)

*Drafted by a read-only agent for the 2026-10-20 10:30 daily review; the review reviews it. Code read on `5017d1bc`; nothing was run.*

**Licence.** §9 rule (6), a provisional chain-U step tagged *(daily, 2026-10-20 — provisional)*, toward chain U's stated target (Phase 9 (9a), `PROJECT_PLAN.md` §6: "every public entry point declares its validation envelope; gated results reproduced through the API within their unmoved bands"). **The target is a reproduction of existing records at their existing bands:** item 259's (`WF-16` step 2, §10 U9) multi-frequency sweep on item 192's F-small fixture — `ANS4_EXACT_STRIP_S_P1_RECORDS` at 10 / 64 / 128 MHz (`tests/validation/test_ans4_exact_strip_records.py:141`, keys `FREQUENCIES` `:135`), degree 1, `-n 2` (`ANS4_EXACT_STRIP_RECORD_RANK_WIDTH` `:227`). **No band, tolerance, record, Done-when or target is set or moved; every band is imported by name:** `LEG_D0_REPRODUCTION_BAND` (`tests/validation/test_port_birdcage_four_port.py:165`), `L0_BAND` (`tests/validation/test_mat5_em_coupled_sphere.py:104`), `ANCHOR_TOL` / `FREQ_ATOL_HZ` (`tests/io/test_touchstone_skrf_crosscheck.py:51–52`), `RECIPROCITY_BAND` (`tests/validation/test_port_lumped_sheet_sweep.py:86`), `PASSIVITY_SIGMA_TOLERANCE` (`test_port_birdcage_four_port.py:129`), `ANS4_EXACT_STRIP_CELL_RECORD` (`tests/validation/test_ops61_exact_strip_ports.py:64`). This is `WF-1` **step 21**. Templates: item 564 (the latest four-port admission and its re-gate set) and item 259 (the Python sweep's anchors). The Python `Study.solve_sweep` has been gated since 2026-10-06; `femem run` has never driven it ("multi-frequency sweeps" is in item 256's Not-claimed list, above).

**Feasibility (source only).**
- `Study.solve_sweep(frequencies_hz, comm)` (`src/fem_em_solver/study.py:963–1043`) builds the geometry once, solves each frequency on that mesh with `dataclasses.replace(self.material, frequency_hz=f)`, stores `self.sweep_result` (a `PortStudySweepResult`, `:915–945`: `frequencies_hz`, `results`, `mesh_time`, `sweep_times`, `s_matrices()`, `write_touchstone(path)` — rank 0 writes one multi-frequency file, then `Barrier`), and leaves `self.result` / `self.frequency_hz` unchanged. It refuses a non-four-port geometry and a non-`LumpedPortDrive` excitation (`TypeError`) and an empty / unsorted / duplicate / non-positive / non-finite list (`_validate_sweep_frequencies`, `:947–960`, `ValueError`).
- **It has no conductor route:** unlike `_study_solve_ports` / `_study_solve_port_fields`, it never calls `_conductor_setup`, so a `LeontovichWall` on the Study is silently ignored, and its `material_map` always names `geom.conductor_tag`, which a hole mesh lacks (item 324's first W1 raised on that missing tag). The YAML layer must therefore refuse `conductor` and `geometry.as_hole: true` with a sweep; the Python method's own silence is reported, not fixed here (For the review, 2).
- `_parse_study` (`src/fem_em_solver/study_yaml.py:263–280`) admits `frequency_hz`, `degree`, `density_kg_m3`; its return is passed to `Study(**study_kwargs)` in `_parse_mapping` (`:1222`, `:1297–1301`), so `sweep_hz` must not enter those kwargs. `run_study_file` calls `study.solve(comm)` once for every EM family (`:1827`); `_run_palace` refuses only `sar` / `b1_plus` / `tune` (`:1730`).

**Executor:** implementer. **Tier:** standard (every new window predicted ≤ 100 s, container `timeout -k 30 180`), complex build, `-n 2`; the re-gates keep their landed commands (each measured ≤ 166 s at `-n ≤ 2` on 2026-10-20; RG4 is `-n 8` and runs only if its condition below fails). **Class:** chain step (solve), **≈ 16.7 measured / ≈ 60 predicted min**: 15 fixed + ≈ 12 `src` + ≈ 10 file and module + windows ≈ 20 (below) + ≈ 3 spare. **Park (rule (d)) if code + W0 are not green by minute 25, or W1a has not started by minute 28.** **Depends on:** nothing open. Items 573 (test-only, registry lines) and 574 (docs only) share no `src` file with this one; **whichever of 572 / 573 lands second rebases, keeps both sides of `validation_envelope.toml` and the generated pages, and re-runs Wgen / Wval / Wapi / Wdoc** before its commit. **Family:** `WF-1` step 21, attempt 1 of 4.

**The change (`src/fem_em_solver/study_yaml.py` only; insertions plus one disclosed edited line, rule (c)).**
- **(E1) Edited line:** `_parse_study`'s `_unknown(sec, ("frequency_hz", "degree", "density_kg_m3"), p, source)` gains `"sweep_hz"`. `_parse_study`'s return is unchanged (it never puts `sweep_hz` in the kwargs).
- **(S1) `_parse_sweep(sec, kind, geometry, excitation, conductor, outputs, has_tune, source) -> tuple`**, a new function after `_parse_tune`, called from `_parse_mapping` after the `tune` block and before `Study(...)` (inserted lines) when `"sweep_hz" in mapping["study"]`. Field `study.sweep_hz`: a YAML list (else `study.sweep_hz`); each entry a number via `_number` at `study.sweep_hz[i]`, finite and > 0 at `study.sweep_hz[i]`; non-empty and strictly increasing at `study.sweep_hz` (the message names the pair, as `_validate_sweep_frequencies` does). Then, at `study.sweep_hz`, naming the conflict: `kind != "four_port_birdcage"`; `excitation` not a `LumpedPortDrive` (capacitor sheets); `conductor is not None`; `geometry.as_hole`; any of `sar` / `b1_plus` / `probes` in `outputs`; a `tune:` section. Returns a float tuple, carried in the parsed outputs under `"sweep_hz"` (the 4-tuple kept, as `tune` rides).
- **(S2) `_run_four_port_sweep(study, kind, outputs, path, out, digest, started, comm) -> StudyRun`**, a new function after `_run_palace`: `sw = study.solve_sweep(outputs["sweep_hz"], comm)`; `touchstone = out / f"{path.stem}.s4p"`; `sw.write_touchstone(touchstone)` (collective); summary `{femem, source, sha256, family, degree, frequency_hz, sweep_hz, ranks, cells, s_matrices: [_pairs(S_f) …], z_matrices: [_pairs(Z_f) …], touchstone, mesh_time_s, sweep_times_s, elapsed_s}`; rank 0 writes `summary.json`, `Barrier`; returns `StudyRun(study=study, result=sw, summary=…, paths={summary, touchstone})`.
- **(S3) Inserted in `run_study_file`** after the thermal-only dispatch heads and before `result = study.solve(comm)` (`:1826–1827`): `if "sweep_hz" in outputs:  # WF-1 step 21 (§9 item 572)` / `return _run_four_port_sweep(...)`. **(S4) Inserted in `_run_palace`** after the `for key in ("sar", "b1_plus", "tune")` loop: `if "sweep_hz" in outputs: raise StudyConfigError(path, "--backend", "palace: study.sweep_hz is not supported")`. The module docstring (`study:` bullet) and `run_study_file`'s docstring each gain one clause (disclosed hunks).
- **Never:** touch `study.py`, `cli.py`, `ports/`, `__all__`; admit a sweep on copper, hole, capacitor, F-human, sphere or thermal files; change `_parse_study`'s return, `_parse_outputs`, `_parse_tune` or any existing branch of `run_study_file`.

**Tracked input:** `tests/package/wf1/f_small_sweep.yaml` = `f_small_four_port.yaml` verbatim (item 192's fixture: degree 1, 64 MHz, `preset: f_small`, `lumped_port_drive`) with the first comment line `# WF-1 step 21 (§9 item 572): …` and `study.sweep_hz: [10000000.0, 64000000.0, 128000000.0]` — decimal literals with a dot (YAML 1.1 reads `1e7` as a string).

**Test module: new `tests/package/test_wf1_yaml_sweep.py`** (item 564's module shape; first docstring paragraph backticks only `WF-1`, `WF-16`, `ANS-4`). Knob `FEM_EM_WF1_S21=1` for the solve nodes; when set, assert `IS_COMPLEX` and `comm.size == ANS4_EXACT_STRIP_RECORD_RANK_WIDTH`, never skip on a mismatch.
- **Imports (names only; never a `test_*` name):** `_records`, `KEYS`, `_port_study` (`tests/package/test_wf16_study_sweep.py:83`, `:103`, `:112` — `_port_study` returns `(study, phantom)`), `C4_HAND_TABLE` (`:88`); `FREQUENCIES` (`test_ans4_exact_strip_records.py:135`); `_rel` (`tests/package/test_wf16_study_four_port.py:115`); `load_touchstone` (`fem_em_solver.ports.touchstone`); `LOGS`, `WF1`, `_mapping`, `_fresh_dir`, `_rc_all`, `_spec_diff` (`tests/package/test_wf1_study_yaml.py:85–169`); the bands above; `_field` as item 564's module reaches it.
- **Nodes:** `test_wf1_s21_parse_and_refusals` (no mesh; every `build` patched to raise, 0 calls), `test_wf1_s21_sweep_file_n2` (W1a), `test_wf1_s21_sweep_python_twin` (W1b).
- **Registry (lines added only):** the parse node under `[names.load_study]` and `[names.StudyConfigError]`, the two solve nodes under `[names.run_study_file]`, each block with `# WF-1 step 21 (§9 item 572)`. **The §9 standing rule applies:** Wgen / Wval / Wapi / Wdoc run in this slot and the regenerated `docs/guide/api/` and `docs/guide/validation.md` are committed; Wval's table rows grow by 4 (print the count).

**Anchors (asserted; bands imported by name).**
- **(K) parse, no mesh:** `_spec_diff(load_study(SWEEP), load_study(WF1 / "f_small_four_port.yaml")) == []`; the parsed `sweep_hz` equals `tuple(FREQUENCIES.values())` by `float.hex`; `kind == "four_port_birdcage"`; `cli.main(["check", f])` rc (0, 0).
- **(X) refusals naming their field** (on the file's mapping, no mesh): `sweep_hz: []`, `6.4e7` (a scalar), `[64000000.0, 10000000.0]`, `[10000000.0, 10000000.0]` ⇒ `study.sweep_hz`; `[0.0]`, `[-1.0]`, `["10 MHz"]` ⇒ `study.sweep_hz[0]`; the sweep on `f_small_copper.yaml`'s mapping (conductor), on `f_small_hole_pec_control.yaml`'s (`as_hole`), on `f_small_hole_pec_capacitors.yaml`'s, on `sphere_mat4.yaml`'s, on `f_human_hole_copper_tune.yaml`'s, and beside `outputs.sar` / `outputs.b1_plus` / `outputs.probes` / `tune:` ⇒ `study.sweep_hz`; `_run_study_file_palace(SWEEP)` ⇒ `--backend` with nothing written; the Python `solve_sweep` raises `ValueError` on the four bad lists (parity of outcome, printed). Every existing parse node green in W0 / RG1, unedited.
- **(R′) per frequency** (W1a, from the summary): `_rel(S_f, _records(1)[k]) ≤ LEG_D0_REPRODUCTION_BAND` for `k` in `KEYS` (item 259 read 6.601e-11 / 9.560e-11 / 9.449e-11, `20261006T033609Z_WF-16-s2-W1-deg1.log.gz`); `summary["cells"] == ANS4_EXACT_STRIP_CELL_RECORD`; one mesh object across `run.result.results`.
- **(T′) the one Touchstone:** `load_touchstone(<stem>.s4p)` has three points, frequencies equal `sweep_hz` within `FREQ_ATOL_HZ`, S per point against the in-memory S ≤ `ANCHOR_TOL` (item 259 (T1′) read 5.541e-14 / 6.409e-14 / 4.969e-14).
- **(G′) `PORT-11`'s identities per frequency:** reciprocity ≤ `RECIPROCITY_BAND`, σ_max ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE` (item 259: reciprocity 4.593e-15 / 1.416e-15 / 6.624e-16).
- **(U′)** the JSON on disk equals the in-memory summary; `summary["s_matrices"]` equals `run.result.s_matrices()` exactly (`_pairs`); `run.study.result is None` (the scalar solve did not run). W1a (rank 0, after an `allgather` check) writes `logs/wf1/s21/sweep-n2.npz` (frequencies, the 3 × 4 × 4 S, cells).
- **(I-Py) θ_file = θ_Python** (W1b, a separate process at the same width): `_port_study(1)[0].solve_sweep(tuple(FREQUENCIES.values()), comm)`; per frequency `_rel(S_py, S_file) ≤ L0_BAND` against the npz (predicted ≤ 1e-12: two processes on the same partition; item 259 vs item 256's records moved only in the 4th digit of a 1e-10 tie). The npz absent ⇒ the node fails with a message (never re-created there).
- **Printed only:** mesh and per-frequency sweep seconds for both legs; the summary's `frequency_hz` (64 MHz) beside `sweep_hz`.

**Negative control (rule (e)).**
- **(c-swap) the 10 ↔ 128 MHz exchange misses — *asserted, backed*:** in W1a, the file's 10 MHz S against `_records(1)["128"]` and its 128 MHz S against `_records(1)["10"]`, each > `LEG_D0_REPRODUCTION_BAND`. Backing: `C4_HAND_TABLE[1]` holds the record-to-record distances 1.2367 (`("10", "128")`) and 1.2763 (`("128", "10")`), reproduced to two digits by item 259's (c4) on this mesh; the readings are predicted equal at that precision (printed beside, not asserted). ≈ 1e9× the band; a relative distance has no ceiling here. A value ≤ 1e-6 means the frequency labels do not reach the solve — stop.
- **(c-ref) the refusals that must not move — *asserted*:** item 564's nine (X) rows and the F-human / sphere probe rows, unedited, in W0 / RG1.

**Windows.** Foreground; **the harness call is the whole Bash command: bare, repo-relative, nothing chained or piped**; Bash `timeout: 660000`; pytest `-s`; node IDs, never `-k`. `pgrep` (`docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3`) must read 0 before W1a, W1b, RG-SB, RG-F and RG-CuSB. W0 runs only after every edit (rule (i)). **Order:** code + file + module, W0, RG1, W9, W1a, W1b, RG-EMP, RG-P, RG-SB, RG-F, RG-Cu, RG-CuSB, RG2, RG3, RG5, Wgen, Wval, Wapi, Wdoc.

| Window | Price | Source |
|---|---|---|
| W0 | ≈ 6 s | item 564 W0, 4 s (`20261020T123357Z_WF-1-s20-W0.log`) |
| RG1 / W9 | 2 s / 5 s | item 564 (`20261020T123415Z_WF-1-s20-RG1.log`, `…123425Z…W9-installed-femem.log`) |
| W1a (`-n 2`) | predicted 60–100 s: one F-small mesh 27.82 s + three sweeps 8.94 / 9.11 / 8.17 s + start-up | `20261006T033609Z_WF-16-s2-W1-deg1.log.gz` |
| W1b (`-n 2`) | predicted 60–100 s: the same solve in Python | as W1a |
| RG-EMP (item 564's W1, the four-port branch it shares) | 106 s | `20261020T123448Z_WF-1-s20-W1-em-probes-n2.log` |
| RG-P / RG-SB / RG-F | 35 / 116 / 161 s | `…123933Z…RGP-sphere-probes.log`, `…124017Z…RGSB-sar-b1.log`, `…124222Z…RGF-f-small.log` |
| RG-Cu / RG-CuSB | 89 / 136 s | `…124512Z…RGCu-copper.log`, `…124650Z…RGCuSB-copper-sar-b1.log` |
| RG2 / RG3 / RG5 | 56 / 29 / 37 s | `…124915Z…RG2.log`, `…125019Z…RG3.log`, `…125056Z…RG5.log` |
| Wgen / Wval / Wapi / Wdoc | 1 / 3 / 18 / 1 s | `20261020T125521Z` / `125531Z` / `125542Z` / `125609Z_WF-1-s20-*.log` |

Windows ≈ 17 min.
```
scripts/testing/run_and_log.sh WF-1-s21-W1a-sweep-file "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_WF1_S21=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_sweep.py::test_wf1_s21_sweep_file_n2 -s -v --tb=short'"
scripts/testing/run_and_log.sh WF-1-s21-W1b-sweep-python "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_WF1_S21=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_sweep.py::test_wf1_s21_sweep_python_twin -s -v --tb=short'"
```
W0: item 564's W0 command verbatim with `tests/package/test_wf1_yaml_sweep.py::test_wf1_s21_parse_and_refusals` and `tests/package/test_wf1_yaml_em_probes.py::test_wf1_s20_parse_and_refusals` both in the node list, label `WF-1-s21-W0`. RG1, W9, RG-P, RG-SB, RG-F, RG-Cu, RG-CuSB, RG2, RG3, RG5: item 564's commands verbatim (above), labels `WF-1-s21-*`. RG-EMP: item 564's W1 command verbatim, label `WF-1-s21-RGEMP`. Wgen / Wval / Wapi / Wdoc: item 564's four (logs above, `:12`), labels `WF-1-s21-Wgen` / `-Wval` / `-Wapi` / `-Wdoc`; Wdoc reads `0 dead`.
- **Not re-run, by construction (valid only if `git diff --stat main -- src/` shows `study_yaml.py` alone, and `git diff main -- src/fem_em_solver/study_yaml.py` shows exactly one changed existing line, `_parse_study`'s `_unknown` tuple, with no changed line inside `_parse_probes`, `_parse_box_probes`, `_run_thermal_box*`, `_run_thermal_ball*`, `_parse_thermal_ball_mapping`, `_parse_thermal_box_mapping`, or the `FHumanBirdcageGeometry` branches of `run_study_file`):** the thermal-only solve gates (item 564's Q4 list — they return at their dispatch heads above the inserted line and never read `study:`), and **RG4** (the F-human tune at `-n 8`): every F-human file's parse runs in W0, which proves its `Study` unchanged, and no line of its solve branch moves. If the condition fails, they run with items 545 / 564's commands (RG4 at `-n 8`, 590 s ceiling — say so in the journal).

**Traps already paid for.**
- **`solve_sweep` silently ignores `study.conductor`** (source above): refuse before the solve; never "fix" `study.py` here.
- **YAML 1.1 reads `1e7` without a dot as a string** — the file's literals are `10000000.0` etc.; (X)'s `"10 MHz"` row is the string case.
- **One mesh, one Touchstone:** never call `study.solve` beside the sweep (the scalar route would mesh again and write a second file); (U′) asserts `run.study.result is None`.
- **The record holds at its own width** (`OPS-41`): `-n 2` only; never compare a `-n 1` or `-n 4` sweep against `P1`.
- **The npz lives under the gitignored `logs/`;** rank 0 writes after an `allgather`, `Barrier` before return; W1b reads it on every rank; a stale npz is overwritten only by W1a.
- `-k a or b` splits argv; `-k 30` on every container `timeout`; never pipe pytest; never import a `test_*` name; `femem` prints on rank 0 only; `capsys.readouterr()` discards earlier prints (item 562's record defect) — print the per-frequency lines after any capture.
- Status 124 or a kill: `pgrep` first; at 0 only, the rubric's orphan-stem sweep labelled `WF-1-s21-cache-sweep`. Never raise a timeout or change `-n`. Rule (i): a `study_yaml.py` edit after W0 re-runs W0 and every RG. Write ≈, never a tilde.

**Scope.** Covered: `study.sweep_hz` on the solid F-small four-port file (item 192's fixture, degree 1, `-n 2`): one mesh, three frequencies, one `.s4p`, each S equal to item 259's `P1` record at `LEG_D0_REPRODUCTION_BAND` and to an independent Python `solve_sweep` at `L0_BAND`; the conflicting combinations refused by field. **Not claimed:** degree 2 (`P2_X1` at `-n 8`); copper / hole / capacitor / F-human sweeps (refused); SAR, `|B₁⁺|` or probes per frequency (refused); a sweep through `--backend palace` (refused) or the installed console script (W9 re-gates `check` only); `femem run --estimate` on a sweep file (it prices one frequency through `Study.estimate`, unchanged, unmeasured); a resonance or tuning statement; any AED statement.

**What a negative result means.** Report and stop; park on `attempt/WF-1-step21-<UTC>`, mark item 572 BLOCKED in the same commit (rule (d)), journal.
- **(K) / (X) red:** parser defect; fix before W1a only.
- **(R′) red with (T′) / (G′) green:** the file's sweep is not the record's — compare `run.result.results[i].problem` with item 259's per-frequency problem (σ / εᵣ retag) and stop; never re-band.
- **(I-Py) red with (R′) green:** two equal Studies in two processes disagree beyond `L0_BAND` — a determinism finding for known-issues, never a band move.
- **(c-swap) not missing:** the frequency labels do not reach the solve — stop.
- **Any RG red:** revert `src/` and stop. **Status 124:** not a result; park with the measured time.

No `PREDICTED RED:` line: a red anywhere stops the slot.

**Status it can move.** `WF-1` stays ✅ (rule (j)). On green its state line's list gains "multi-frequency sweep from one file (U71)", item 256's Not-claimed "multi-frequency sweeps" clause is retired by a **step-21 record**: "`study.sweep_hz` drives `Study.solve_sweep` through `run_study_file`: one mesh, one `.s4p`; per-frequency S vs `P1` <a> / <b> / <c> (`LEG_D0_REPRODUCTION_BAND`), Touchstone round trip ≤ <t>, θ_file = θ_Python ≤ <x> (`L0_BAND`); the 10 ↔ 128 MHz swap misses by <p> / <q>; `src`: `study_yaml.py`, insertions plus one disclosed edited line." §10 U71 DONE.

**Commits.**
1. `WF-1 step 21 (§9 item 572, U71): femem run drives Study.solve_sweep from study.sweep_hz` — `src/fem_em_solver/study_yaml.py`; the new YAML file and module; `validation_envelope.toml`, the regenerated `docs/guide/api/` and `docs/guide/validation.md`; every `WF-1-s21-*` log and its test-results rows; `PROJECT_PLAN.md` (§7 `WF-1` state line, §9 item 572 DONE line, §10 U71); this file's outcome paragraph. The message names the edited `_parse_study` line.
2. `attempts: WF-1 step 21 §9 item 572 journal` (`docs/testing/attempts.md` only).

**Spawn line:** "Use the implementer agent to execute `WF-1` step 21, §9 item 572 (U71), per `docs/planning/chunks/WF-1.md` "§9 item 572" (read the brief and any review ruling beneath it first; the rulings hold). The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. `src/` change: `study_yaml.py` only — `"sweep_hz"` added to `_parse_study`'s `_unknown` tuple (the one edited line), `_parse_sweep`, its call in `_parse_mapping`, `_run_four_port_sweep`, the inserted head in `run_study_file`, the inserted refusal in `_run_palace`; no `study.py` / `cli.py` edit. Write the YAML file and the module whole, then run W0, RG1, W9, W1a, W1b, RG-EMP, RG-P, RG-SB, RG-F, RG-Cu, RG-CuSB, RG2, RG3, RG5, Wgen, Wval, Wapi, Wdoc in order, foreground, every harness Bash call with `timeout: 660000`, pytest with `-s`, container `timeout -k 30`; `pgrep -c python3` must read 0 before W1a, W1b, RG-SB, RG-F and RG-CuSB. W1a before W1b (W1b reads its npz). RG4 and the thermal-only solve gates run only if the brief's git-diff condition fails. Never import a `test_*` name; never edit an existing test assertion. You are the executor: do not spawn agents, never background a harness window. Park if code + W0 are not green by minute 25 or W1a has not started by minute 28. On a red, park on `attempt/WF-1-step21-*`, mark item 572 BLOCKED in the same commit, and report."

**For the review.**
1. **Key placement.** `study.sweep_hz` keeps the frequency keys together and needs no new top-level key; a top-level `sweep: {frequencies_hz}` (the `tune:` shape) would also need `"sweep"` added to the thermal-only refusal tables (`_BALL_REFUSED` and the box's), or a thermal file would silently ignore it. The brief takes `study.sweep_hz`; `study.frequency_hz` stays required and is **not** required to be one of the sweep points (the Python method does not require it either).
2. **Two defects found while drafting, outside this item's `src` scope:** (a) `Study.solve_sweep` ignores `study.conductor` silently (a copper Study sweeps as if PEC-free σ = 800), and (b) `_run_palace` refuses `sar` / `b1_plus` / `tune` but not `outputs.probes`, so `femem run --backend palace` on U70's file drops the probes silently. Each is a lane-1 refusal item (`study.py` / `study_yaml.py`, one inserted `raise` each); neither is folded in here unless the review rules so.
3. **RG4 by construction.** The brief does not re-run the `-n 8` F-human gate when the git-diff condition holds (W0 parses every F-human file). The literal "every YAML gate" adds RG4 (199 s, `-n 8`).

**Outcome (2026-10-20 15:00 slot, implementer): green.** **Step-21 record (§9 item 572, U71, provisional):** `study.sweep_hz` drives `Study.solve_sweep` through `run_study_file`: one mesh, one `.s4p`; per-frequency S vs `P1` 6.601e-11 / 9.557e-11 / 9.449e-11 (`LEG_D0_REPRODUCTION_BAND`) on 117 206 cells (`ANS4_EXACT_STRIP_CELL_RECORD`), Touchstone round trip ≤ 6.873e-14 (`ANCHOR_TOL`, \|Δf\| = 0), reciprocity ≤ 4.593e-15, σ_max 0.999994 / 0.999760 / 0.999136, θ_file = θ_Python ≤ 9.010e-15 (`L0_BAND`; 10 / 64 / 128 MHz 9.010e-15 / 9.110e-16 / 0); the 10 ↔ 128 MHz swap misses by 1.2367 / 1.2763 (equal to `C4_HAND_TABLE[1]` at four digits); `run.study.result is None`; `src`: `study_yaml.py`, insertions plus one disclosed edited line (`_parse_study`'s `_unknown` tuple) and one re-punctuated module-docstring line (`density_kg_m3 …` now ends in a comma before the new `sweep_hz` bullet text). Refusals (X): 16 rows by field (`study.sweep_hz` / `study.sweep_hz[0]`), `--backend palace` ⇒ `--backend` with nothing written, Python `solve_sweep` raises `ValueError` on the five numeric bad lists. Windows: W1a 53 s (run_study_file 50.3 s: mesh 25.28 s, sweeps 8.11 / 8.02 / 8.34 s), W1b 53 s (Python 49.7 s). **The first W1a (`20261020T201341Z_WF-1-s21-W1a-sweep-file.log`) was red on (U′) "JSON == memory" only:** the new module compared the summary's rank-local timings (`mesh_time_s`, `sweep_times_s`) against rank 0's JSON on rank 1 too; the template excluded only `elapsed_s`. Fixed in the new test (timing keys excluded off rank 0, full equality still asserted on rank 0, a differing-keys print added — it reads `[]`); no anchor or band moved, and every anchor reading was already green in that run. Every YAML gate re-ran per ruling (1), RG4 at `-n 8` included: W0 4 s (21 passed), RG1 3 s, W9 15 s, RG-EMP 107 s, RG-P 37 s, RG-SB 116 s, RG-F 165 s, RG-Cu 90 s, RG-CuSB 142 s, RG2 61 s, RG3 31 s, RG4 206 s (9 passed), RG5 41 s, Wgen / Wval / Wapi / Wdoc 0-status, page node rows 323, 0 dead. Logs `20261020T201259Z` … `20261020T203455Z_WF-1-s21-*.log`.

### §9 item 573 — 2026-10-20 (§10 U72, chain step (rule (6), provisional); lane 1)

*Drafted by a read-only agent for the 2026-10-20 10:30 daily review; the review reviews it. Code read on `5017d1bc`; nothing was run.*

**Licence.** §9 rule (6), a provisional chain-U step tagged *(daily, 2026-10-20 — provisional)*, toward chain U's stated target (above). It measures a path U70 left "reachable, not measured, not claimed" (item 564's Outcome, last sentence): `outputs.probes` on the copper F-small file. **The records reproduced are existing ones at their existing bands:** `ANS6_EXACT_STRIP_CLASS_MEANS_P1["64"]["cu"]` (`tests/validation/test_ans6_exact_strip.py`) at `LEG_D0_REPRODUCTION_BAND` (item 343 / 394 read 4.093e-11 / 8.673e-11 / 1.376e-10, `20261012T140300Z_WF-1-s8-W1.log`), `ANS6_EXACT_STRIP_HOLE_CELL_RECORD` = 80 798 (`:112`) at `ANS6_EXACT_STRIP_RECORD_RANK_WIDTH` = 2 (`:144`), and U70's own target θ_file = θ_Python at `L0_BAND` (`tests/validation/test_mat5_em_coupled_sphere.py:104`). Also imported: `EXACT_IDENTITY_RTOL` (`tests/validation/test_birdcage_sar_mass_averaged.py:163`), `WIDTH_RTOL` (`tests/validation/test_post8_step3a_averaged_row1.py:88`). **No band, tolerance, record, Done-when or target is set or moved.** This is `WF-1` **step 22**. Templates: item 564 (the probe anchors) and item 394 (the copper file through YAML, test-only).

**Feasibility (source only).** `_parse_mapping`'s in-phantom block (`src/fem_em_solver/study_yaml.py:1269–1283`) admits `probes` on any `four_port_birdcage` that is not `capacitor_sheets`; a `conductor:` with `geometry.as_hole: true` passes `_parse_conductor`. `run_study_file`'s four-port branch (`:1856–1870`) calls `_em_probe_summary`, which reads `study._solve_port_fields(comm)` — the field route that applies the conductor since item 388 (`study.py:1968–2050`; item 394's (K1) asserted one conductor term on the YAML Study). σ and ρ are `study.material.sigma` / `study.density_kg_m3`, and the phantom (r < 0.03 m, \|z\| < 0.04 m) is the preset's on the hole mesh too. **No `src` edit is needed or allowed.**

**Executor:** implementer. **Tier:** standard (W1 predicted 100–140 s, container `timeout -k 30 180`), complex build, `-n 2`. **Class:** lane 1 (test-only chain step, solve), **≈ 8.8 measured / ≈ 35 predicted min**: 15 fixed + ≈ 10 file and module + windows ≈ 4 (below) + ≈ 6 spare. **Park (rule (d)) if W0 is not green by minute 20, or W1 has not started by minute 23.** **Depends on:** nothing open. Item 572 also adds registry lines: **whichever of 572 / 573 lands second rebases, keeps both sides of `validation_envelope.toml` and the generated pages, and re-runs Wgen / Wval / Wapi / Wdoc** before its commit. **Family:** `WF-1` step 22, attempt 1 of 4.

**Tracked input:** `tests/package/wf1/f_small_copper_em_probes.yaml` = `f_small_copper.yaml` verbatim (degree 1, 64 MHz, `preset: f_small`, `as_hole: true`, `leontovich` σ = 58000000.0) with the first comment line `# WF-1 step 22 (§9 item 573): …` and `outputs: {probes: {points: […]}}` carrying **U70's five points as the same decimal literals** as `f_small_em_probes.yaml` (so (K) is hex-equal by construction).

**Test module: new `tests/package/test_wf1_yaml_copper_em_probes.py`** (item 564's module shape; first docstring paragraph backticks only `WF-1`, `ANS-6`, `POST-8`). Knob `FEM_EM_WF1_S22=1` for the solve node; when set, assert `IS_COMPLEX` and `comm.size == ANS6_EXACT_STRIP_RECORD_RANK_WIDTH`, never skip on a mismatch.
- **Imports (names only; never a `test_*` name):** `POINTS`, `EM_PROBES` (`tests/package/test_wf1_yaml_em_probes.py:66`, `:72`); `COPPER_FILE`, `_py_copper` (`tests/package/test_wf1_yaml_copper.py:84`, `:98`); `_rel_means` (`tests/package/test_wf16_study_copper.py:76`); `ANS6_EXACT_STRIP_CLASS_MEANS_P1`, `ANS6_EXACT_STRIP_HOLE_CELL_RECORD`, `ANS6_EXACT_STRIP_RECORD_RANK_WIDTH`, `CLASSES`, `_class_means` (`tests/validation/test_ans6_exact_strip.py`, as `test_wf1_yaml_copper.py:68–75` imports them); `point_sar` (`fem_em_solver.post.sar`); `_bitwise_equal_across_ranks` (`tests/post/test_point_owner_rule.py:84`); `LOGS`, `WF1`, `_fresh_dir`, `_rc_all`, `_spec_diff` (`tests/package/test_wf1_study_yaml.py`); the bands above.
- **Nodes:** `test_wf1_s22_parse` (no mesh; every `build` patched to raise, 0 calls), `test_wf1_s22_copper_em_probes_n2` (W1).
- **Registry (lines added only):** the parse node under `[names.load_study]`; the solve node under `[names.run_study_file]` and `[names.LeontovichWall]`, each block with `# WF-1 step 22 (§9 item 573)`. **The §9 standing rule applies:** Wgen / Wval / Wapi / Wdoc run in this slot and the regenerated pages are committed; Wval's table rows grow by 3 (print the count).

**Anchors (asserted; bands imported by name).**
- **(K) parse, no mesh:** `_spec_diff(load_study(CU_PROBES), load_study(COPPER_FILE)) == []`; the parsed points equal `POINTS` by `float.hex`; the conductor parses as `LeontovichWall` with σ hex-equal to the copper file's; `cli.main(["check", f])` rc (0, 0).
- **(K1) the fixture** (W1): `summary["cells"] == ANS6_EXACT_STRIP_HOLE_CELL_RECORD`; `summary["conductor"] == {"kind": "leontovich", "sigma_s_per_m": 58000000.0}`.
- **(R-Cu)** `_rel_means(_class_means(S_summary), ANS6_EXACT_STRIP_CLASS_MEANS_P1["64"]["cu"])` per class ≤ `LEG_D0_REPRODUCTION_BAND` (predicted 4.093e-11 / 8.673e-11 / 1.376e-10, item 394).
- **(U′) the summary is the API on its own fields:** the summary values equal `point_sar(f.e_real, f.e_imag, POINTS, sigma=…, rho=…, comm=…)` on `run.study.fields_result` exactly, on the same five-point batch; `_bitwise_equal_across_ranks` on the value array; JSON = memory; `e_abs² σ / (2ρ)` against `sar_w_per_kg` ≤ `EXACT_IDENTITY_RTOL` (item 564 read 2.213e-16 on the solid fixture).
- **(I-Py) θ_file = θ_Python — U70's target on copper:** an independent `_py_copper()` Study, `_solve_port_fields(comm)` then `point_sar` per port; every summary SAR value against it, relative, ≤ `L0_BAND` (predicted ≤ 1e-12: item 394's (I-SAR) read 4.902e-13 for the 10 g SAR on this fixture through the same two front ends; item 564 read 0.000e+00 on the solid one).
- **Printed only:** per-point SAR and \|E\| per 1 V drive per port; the build seconds of both Studies.

**Negative control (rule (e)).**
- **(c-disp) a displaced probe misses — *asserted*:** point 0 displaced 1 mm radially, `(0.021, 0.0, 0.0)`, through `point_sar` on `run.study.fields_result`, against the summary's point-0 value per port: relative difference > `WIDTH_RTOL`. **First measurement on the hole mesh — predicted 1e-2 to 3e-1** (item 564 read 5.62e-02 … 2.12e-01 on the solid mesh at the same point; δSAR/SAR ≈ 2δr/r ≈ 0.1). The bar is ≥ 1e7× below the prediction; no ceiling binds. A value ≤ 1e-6 means the probe ignores its coordinates — stop.
- **(c-ref) — *asserted*:** item 564's parse node (`test_wf1_s20_parse_and_refusals`) green, unedited, in W0.

**Windows.** Foreground; **the harness call is the whole Bash command: bare, repo-relative, nothing chained or piped**; Bash `timeout: 660000`; pytest `-s`; node IDs, never `-k`. `pgrep -c python3` (in the container) must read 0 before W1. **Order:** file + module, W0, W1, Wgen, Wval, Wapi, Wdoc.

| Window | Price | Source |
|---|---|---|
| W0 | ≈ 6 s | item 564 W0 (`20261020T123357Z_WF-1-s20-W0.log`) |
| W1 (`-n 2`) | predicted 100–140 s: `run_study_file` with two hole-mesh builds 67.9 s, Study A's field route 39.5 s, five-point evaluations ≈ 1 s, start-up | `20261012T140300Z_WF-1-s8-W1.log` (item 394, 151 s with a third build) |
| Wgen / Wval / Wapi / Wdoc | 1 / 3 / 18 / 1 s | `20261020T125521Z` / `125531Z` / `125542Z` / `125609Z_WF-1-s20-*.log` |

Windows ≈ 3 min.
```
scripts/testing/run_and_log.sh WF-1-s22-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_copper_em_probes.py::test_wf1_s22_parse tests/package/test_wf1_yaml_em_probes.py::test_wf1_s20_parse_and_refusals tests/package/test_wf1_yaml_copper.py::test_wf1_copper_yaml_parse_and_refusals -v -s --tb=short'"
scripts/testing/run_and_log.sh WF-1-s22-W1-copper-em-probes "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_WF1_S22=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_copper_em_probes.py::test_wf1_s22_copper_em_probes_n2 -s -v --tb=short'"
```
Wgen / Wval / Wapi / Wdoc: item 564's four commands verbatim (logs above, `:12`), labels `WF-1-s22-Wgen` / `-Wval` / `-Wapi` / `-Wdoc`; Wdoc reads `0 dead`. **No re-gate:** confirm `git diff --stat main -- src/` is empty and say so (rule (c) does not bind).

**Traps already paid for.**
- **Never hand `point_sar` the complex field twice** (`tests/validation/test_birdcage_sar_map.py:295–302`, it silently halves SAR); pass `f.e_real`, `f.e_imag`.
- **A bitwise claim names its batch** (`POST-8` item 537): (U′) compares on the same five-point batch, never re-batched.
- **The copper Study's field route is the one that carries the conductor** (item 388's fix); never compare against `Study.solve()`'s S block for the fields, and never build a PEC twin here (the third build is what took item 394's W1 to 151 s).
- **One field solve per Study:** `_solve_port_fields` caches on `study.fields_result`; build no Study beyond the declared Python twin.
- **YAML 1.1:** the five points and σ stay plain decimals with a dot.
- `-k a or b` splits argv; `-k 30` on every container `timeout`; never pipe pytest; never import a `test_*` name; `capsys.readouterr()` discards earlier prints (item 562).
- Status 124 or a kill: `pgrep` first; at 0 only, the rubric's orphan-stem sweep labelled `WF-1-s22-cache-sweep`. Never raise a timeout or change `-n`. Write ≈, never a tilde.

**Scope.** Covered: `outputs.probes: {points}` on the copper (Leontovich, σ = 5.8e7 S/m) F-small hole file at degree 1, 64 MHz, `-n 2`: the `cu` class means at `LEG_D0_REPRODUCTION_BAND`, point SAR and \|E\| per 1 V drive equal to an independent Python copper Study at `L0_BAND`. **Not claimed:** width invariance on the hole mesh (U70's `-n 1 / 4` legs were on the solid mesh); a Cu-vs-PEC point gap (no PEC twin); degree 2; F-human or capacitor files; absolute or compliance SAR; any AED statement; `POST-8`'s census (no new site — U70's `_em_probe_summary` is the site).

**What a negative result means.** Report and stop; park on `attempt/WF-1-step22-<UTC>`, mark item 573 BLOCKED in the same commit (rule (d)), journal.
- **(K) red:** the file is not the copper file — fix the file before W1 only.
- **(R-Cu) red:** the YAML Study's S is not the record on today's `main` — re-run item 564's RG-CuSB command once as a diagnostic; if it is also red, a `main` regression: known-issues entry, stop.
- **(I-Py) red with (R-Cu) green:** two equal copper Studies disagree pointwise beyond `L0_BAND` — record it, never re-band.
- **(c-disp) not missing:** blind instrument; stop.

No `PREDICTED RED:` line: a red anywhere stops the slot.

**Status it can move.** `WF-1` stays ✅ (rule (j)). On green, U70's "copper F-small reachable, not measured, not claimed" clause is retired by a **step-22 record**: "`outputs.probes` on the copper F-small file: `cu` class means <a> / <b> / <c> (`LEG_D0_REPRODUCTION_BAND`), 80 798 cells, θ_file = θ_Python ≤ <x> (`L0_BAND`); the 1 mm displaced probe misses by <p>; no `src` edit." §10 U72 DONE.

**Commits.**
1. `WF-1 step 22 (§9 item 573, U72): point probes on the copper F-small YAML file` — the YAML file and module; `validation_envelope.toml`, the regenerated `docs/guide/api/` and `docs/guide/validation.md`; every `WF-1-s22-*` log and its test-results rows; `PROJECT_PLAN.md` (§7 `WF-1` state line, §9 item 573 DONE line, §10 U72); this file's outcome paragraph.
2. `attempts: WF-1 step 22 §9 item 573 journal` (`docs/testing/attempts.md` only).

**Spawn line:** "Use the implementer agent to execute `WF-1` step 22, §9 item 573 (U72), per `docs/planning/chunks/WF-1.md` "§9 item 573" (read the brief and any review ruling beneath it first; the rulings hold). The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. No `src/` edit. Write the YAML file and the module whole, then run W0, W1, Wgen, Wval, Wapi, Wdoc in order, foreground, every harness Bash call with `timeout: 660000`, pytest with `-s`, container `timeout -k 30`; `pgrep -c python3` must read 0 before W1. Pass `f.e_real` / `f.e_imag` to `point_sar`, never the complex field twice; build no PEC twin. Never import a `test_*` name; never edit an existing test assertion. If item 572 has landed first, rebase and re-run Wgen / Wval / Wapi / Wdoc. You are the executor: do not spawn agents, never background a harness window. Park if W0 is not green by minute 20 or W1 has not started by minute 23. On a red, park on `attempt/WF-1-step22-*`, mark item 573 BLOCKED in the same commit, and report."

**For the review.**
1. **No PEC twin.** Item 394's (N1) PEC column is the control that separates "the conductor reached the field route"; here (R-Cu) plus (I-Py) against a copper Python Study carry that, and the third build would put W1 at item 394's 151 s against a 180 s ceiling. Adding it as a second window (≈ 70 s) is the alternative.

**Outcome (2026-10-20, 21:30 slot, take-next; implementer).** Green, every anchor, no `src` edit (`git diff --stat main -- src/` empty; rule (c) does not bind). W0 `20261020T215730Z_WF-1-s22-W0.log` (3 s, 7 passed): (K) spec diff `[]`, points hex-equal to `POINTS` and to `f_small_em_probes.yaml`'s, conductor `LeontovichWall` σ hex `0x1.ba81400000000p+25` = the copper file's, `femem check` rc (0, 0), 0 build calls; (c-ref) item 564's and item 343's parse nodes green, unedited. W1 `20261020T215743Z_WF-1-s22-W1-copper-em-probes.log` (93 s, `-n 2`, complex): `run_study_file` 58.8 s; (K1) 80 798 cells, conductor `{"kind": "leontovich", "sigma_s_per_m": 58000000.0}`; (R-Cu) self 4.093e-11 / adjacent 8.673e-11 / opposite 1.376e-10 ≤ `LEG_D0_REPRODUCTION_BAND` (= item 394's digits); (U′) summary = `point_sar` on the run's fields exactly, bitwise across ranks, JSON = memory, \|E\| identity 2.179e-16 ≤ `EXACT_IDENTITY_RTOL`; (I-Py) the Python copper twin (29.9 s, 80 798 cells) ≤ 1.824e-13 ≤ `L0_BAND`; (c-disp) per port 6.5011e-02 / 2.4157e-01 / 8.0900e-02 / 3.6739e-03, min 3.6739e-03 > `WIDTH_RTOL` — **P4's value sits below the predicted 1e-2 to 3e-1 range** (P4's own point-0 SAR is small, 1.19e-07 W/kg per 1 V, so its gradient need not follow 2δr/r), still ≥ 3e6× above the bar and far above the 1e-6 blind-instrument stop. Wgen / Wval / Wapi / Wdoc `20261020T215936Z` / `215938Z` / `215941Z` (22 passed) / `220002Z_WF-1-s22-*.log` (2 / 3 / 20 / 1 s), Wdoc `dead=0`; validation page +3 rows (416 table lines), item 572's lines kept. Step-22 record: "`outputs.probes` on the copper F-small file: `cu` class means 4.093e-11 / 8.673e-11 / 1.376e-10 (`LEG_D0_REPRODUCTION_BAND`), 80 798 cells, θ_file = θ_Python ≤ 1.824e-13 (`L0_BAND`); the 1 mm displaced probe misses by ≥ 3.674e-03; no `src` edit." U72 DONE.

### Review ruling 2026-10-20 10:30 (items 572, 573; rule (6), provisional)

(1) **572:** `study.sweep_hz` is accepted as the key, provisionally, until the 10-21 weekly rules. **Every YAML gate re-runs, RG4 at `-n 8` included** (199 s). The git-diff skip condition in the brief is struck, because the ratified U65–U67 pattern re-runs every gate.
(2) **572 refuses `study.sweep_hz` together with `study.conductor`**, by field, as one inserted refusal. The drafter found that `Study.solve_sweep` never calls `_conductor_setup` (`src/fem_em_solver/study.py:963–1043`), so a copper sweep would silently drop its copper boundary. The refusal does not fix the Python path; that fix is the weekly's to pose.
(3) **573:** accepted as drafted: no PEC twin and no `-n 1 / 4` width legs, both not claimed. The second-window alternative is not taken.
(4) **A finding for the weekly, not queued:** `_run_palace` refuses `sar` / `b1_plus` / `tune` but not `outputs.probes` (`src/fem_em_solver/study_yaml.py:1730`). So `femem run --backend palace` on U70's file drops the probes silently. Found by reading the code, not observed in a run.

### §9 item 588 — 2026-10-22 (§10 U74, lane-2 feature step)

*Drafted read-only for the 2026-10-22 03:00 daily review; the review rules on it.* Code read on `11fc13b0`; nothing was run.

**Licence.** §10 chain U, U74 (`PROJECT_PLAN.md:4947–4952`): "`n_slots` on `ring_mode_rows` (additive, default 16) and `Study.tune(n_slots=…)`: equal to C10a's `ring_c_modes` at N = 16 / 32 to 0 and to the U2 closed form at `REUSE_REPRODUCTION_RTOL`, the two recorded F-human `C_mode` reproduced bitwise through the default; every `Study.tune` / `tune:` gate re-run. Serial on nothing." The 10-21 weekly's ruling (`PROJECT_PLAN.md:4787–4792`; `docs/planning/chunks/GEO-38.md:1560`): "`ring_mode_rows` for N slots — ruled: an additive `src/` keyword, chain U (U74) … default 16, so every F-human consumer is bitwise … A future F-clinical `Study.tune` is the consumer." It answers the 10-21 03:00 ruling (2) on item 581 (`GEO-38.md:1552`). **No band, tolerance, record, Done-when or target is set or moved; every band is imported by name:** `REUSE_REPRODUCTION_RTOL` = 1e-6 (`tests/validation/test_port_birdcage_ring_matrix.py:171`, never the 1e-12 namesake at `tests/validation/test_port19_factor_reuse.py:82`), `EXACT_IDENTITY_RTOL` (`tests/validation/test_birdcage_sar_mass_averaged.py:163`), `REDUCTION_BAND` (`tests/validation/test_port_lumped_rlc_termination.py:88`). This is **`WF-16` step 13** (`Study.tune` is `WF-16`'s: steps 6 / 12; no `WF-16-s13` label or record exists), filed here as assigned. **Class: lane-2 feature step** — it adds a public keyword to `src/`, the item-564 precedent ("lane 2 feature admission"); the weekly priced it at the chain-solve median (`PROJECT_PLAN.md:5012`).

**What exists (`11fc13b0`).**
- `src/fem_em_solver/ports/ring_modes.py` is 16-slot by construction: `N_SLOTS = 16` (`:38`); `ring_port_map_from_names` hard-codes P17…P48 (`:42–53`); `_canonical_permutation` demands 32 pairs (`:56–65`); `_dft_matrix` (`:68–71`), `k_blocks` (`:74–85`) and `k_spectrum` (`:100–113`, folds `range(1, 8)` and bin 8) read the module constant; `ring_mode_rows` (`:164–230`) loops `range(N_SLOTS // 2 + 1)`, slices rings by `N_SLOTS`, and builds `others` from `range(2 * N_SLOTS)`. Last change `9b90d7b3` (`PORT-24` step 1); unchanged since.
- `Study.tune` (`src/fem_em_solver/study.py:1075–1118`) dispatches an F-human `LumpedPortDrive` Study to `_study_tune_ring` (`:1121–1150`), which calls `ring_port_map_from_names`, `_canonical_permutation`, `k_blocks` and `antisymmetric_eigenpair` directly (not `ring_mode_rows`) and loops `range(9)`. `FHumanBirdcageGeometry.leg_count` defaults to 16 (`:455`); `ring_ordinals()` = `range(leg_count + 1, 3·leg_count + 1)` (`:500–501`), so a 32-leg geometry names its ring gaps P33…P96 — the F-clinical naming C10a verified (`GEO-38.md:1554`, (K3) 64/64).
- C10a's test-local reader `ring_c_modes` (`tests/validation/test_geo38_f_clinical_c_mode.py:138–164`) with `dft_k_domain` (`:115`), `canonical_permutation` (`:125`), `map_from_names_64` (`:182`), `_sym_blocks` / `_block_circulant` / `_shuffled` (`:199`, `:214`, `:229`), `C_MODE_CU_D1_RECORD_F` (`:74`), `C_MODE_CU_D2_COMPARAND_F` (`:83`), `H12C_A2_NPZ` / `H16B_D2_NPZ` (`:67–68`). It read (K1) 1.554e-15 / 1.443e-15 at N = 32 / 16, (K2) and (I16) 0.000e+00 against `ring_mode_rows`, `C_mode` 1.181185799070e-10 F (d1) and 1.145221310828e-10 F (d2) (`20261021T093341Z_GEO-38-C10a-W1.log:65–67`, `:73`, `:83–84`).
- The pinned-pre-change bitwise pattern: `tests/unit/test_port24_esr_circuit.py:67–79` (`git show <sha>:…` inside the container, `exec` into a module with `__package__ = "fem_em_solver.ports"`; `_canon` `:120–121`). The injected-result `Study.tune` pattern (no mesh): `tests/package/test_wf16_study_tune_f_human_copper.py:143–202` (`SimpleNamespace(s_matrix=…, port_definitions=…)`, every `build` stubbed to raise).
- `study_yaml.py` calls `study.tune(**spec, comm=comm)` (`:1589`) with `spec` ⊆ `{"k"}` (`_parse_tune`, `:1170–1193`): YAML never reaches `n_slots`, so `study_yaml.py` is not edited.

**Consumers of `ring_modes` (grep of `src tests examples scripts`, `11fc13b0`)** — each re-gated below unless stated: `study.py` (`_study_tune_ring`); `tests/unit/test_ring_mode_rows.py`, `tests/unit/test_port24_esr_circuit.py` (pinned `b2174873` bitwise on `logs/ans8/s32.npz` — stays the stronger default-bitwise gate), `tests/unit/test_ans8_private_comparison.py` (synthetic, scratch copy; never opens `aed_results/`); `tests/validation/test_ans8_copper_degree2_rows.py` (H25), `tests/validation/test_ans8_degree2_rows.py` (`ring_port_map_from_names` only), `tests/validation/test_geo38_f_clinical_c_mode.py` (C10a); `tests/package/test_wf16_study_tune_f_human_copper.py`; the `ans:8` script `08_f_human_quadrature_mode_content.py` (`rm.N_SLOTS`, `rm.ring_mode_rows`; reached through `tests/unit/test_ans8_v3_rows.py` / `test_ans8_v3_definition.py`). Not re-run: `examples/mri/23_f_human_study_tune_ring.py` (an example, not a gate — For the review 5) and `examples/ansys_benchmarks/ans8_*/aed/ans8_private_rows_v2.py` (private, never opened). `tests/validation/test_wf10_f_human_copper_esr_root.py` and `test_ans7_exact_strip.py` name `ring_modes` in comments only.

**`Study.tune` / `tune:` gates (with their last logs):** `test_wf16_study_tune.py` (`20261015T213242Z_WF-16-s12-RG1-tune.log`, `-n 2`, 69 s), `test_wf16_study_tune_capacitors.py` (`20261015T201110Z_WF-16-s6b-W1.log`, `-n 2`, 69 s), `test_wf16_study_tune_deg2.py` (`20261014T111026Z_WF-16-s6c-W1.log`, `-n 8`, 82 s), `test_wf16_study_tune_deg2_copper.py` (`20261014T185251Z_WF-16-s6d-W1.log`, `-n 8`, 83 s), `test_wf16_study_tune_f_human_copper.py` (`20261015T213534Z_WF-16-s12-W1.log`, `-n 8`, 205 s), `test_wf1_yaml_tune.py` (RG2, `20261020T202817Z_WF-1-s21-RG2.log`, 61 s), `test_wf1_yaml_tune_f_human.py` (RG4, `20261020T202959Z_WF-1-s21-RG4.log`, `-n 8`, 206 s), `test_wf1_yaml_capacitors.py` (RG5, `20261020T203329Z_WF-1-s21-RG5.log`, 41 s), the `tune:` parse nodes in item 572's W0 (`20261020T203852Z_WF-1-s21-W0-r2.log`, 4 s), `test_wf17_tuning_yaml_page.py` (`20261015T214647Z_WF-17-s6-W0.log`, 5 s). `test_wf13_price_t49.py` imports a constant only.

**Code (`src/`: `ring_modes.py` and `study.py` only; keyword insertions, disclosed edited lines, rule (c)).**
- **`ring_modes.py`** — every new keyword defaults to `N_SLOTS`, and with the default each function executes the same arithmetic in the same order (bitwise by construction; (B) proves it):
  - `_check_n_slots(n_slots)` (new): an `int` (not `bool`), even, ≥ 4 (k = 2 is read), else `TypeError` / `ValueError` naming `n_slots`.
  - `ring_port_map_from_names(names, n_slots=N_SLOTS)`: P(N+1+j) → (j, bottom), P(2N+1+j) → (j, top); the message names the range. At N = 16 this is P17…P48 as now; at N = 32 it is C10a's `map_from_names_64`.
  - `_canonical_permutation(port_map, n_slots=N_SLOTS)`, `_dft_matrix(n_slots=N_SLOTS)` (`np.arange(n)`, `/ np.sqrt(n)` — the same expressions), `k_blocks(z_canonical, n_slots=N_SLOTS)`, `k_spectrum(ring_currents, n_slots=N_SLOTS)` (folds `range(1, n // 2)` plus bin `n // 2`; a length ≠ n raises — disclosed, it can only fire on a misuse).
  - `ring_mode_rows(…, esr_ohm=0.0, n_slots=N_SLOTS)`: `n` threaded through the permutation, `k_blocks`, the mode loop `range(n // 2 + 1)`, `canon`, the ring slices and `others`. Docstring and module docstring gain one clause each (`F_kj = e^{−2πi jk/N}/√N`).
  - **Not edited:** `ring_mode_rows_v2`, `terminated_currents`, `antisymmetric_eigenpair`, `_esr_canonical`, `_gap_admittance`, `N_SLOTS`, `RINGS`, `__all__`.
- **`study.py`**:
  - `_study_tune(self, k=1, *, comm=…, n_slots=None)`: the F-human dispatch passes `n_slots`; on every other route an explicit `n_slots` raises `TypeError("n_slots: defined on the f_human_birdcage ring read only")` **before** `_refuse_f_human` and any solve (one inserted `if`). The four-port read is otherwise untouched.
  - `_study_tune_ring(self, k=1, *, comm=…, n_slots=None)`: after the existing `k` checks and **before** `self.solve`, `n = self.geometry.leg_count if n_slots is None else n_slots`; `_check_n_slots(n)`; `n != self.geometry.leg_count` ⇒ `ValueError` naming `n_slots`. Then `ring_port_map_from_names(names, n_slots=n)`, `_canonical_permutation(…, n_slots=n)`, `k_blocks(…, n_slots=n)`, `range(n // 2 + 1)` in place of `range(9)`. With `leg_count` 16 every existing call takes the same path.
  - `TuneResult`'s docstring "k = 0…8" → "k = 0…N/2" (one disclosed line).
  - **Never:** touch `study_yaml.py`, `cli.py`, the four-port read's arithmetic, `__all__`, any existing test assertion; add `n_slots` to the YAML schema.

**Test module: new `tests/package/test_wf16_study_tune_n_slots.py`** (first docstring paragraph backticks only `WF-16`, `GEO-38`, `PORT-24`). Knob `FEM_EM_WF16_S13` (unset ⇒ the records node skips; `1` ⇒ it runs and an absent npz **fails**; anything else ⇒ `ValueError`).
- **Imports (names only; never a `test_*` name):** `import tests.validation.test_geo38_f_clinical_matrix as c3` (whole module, C10a's pattern); from `tests.validation.test_geo38_f_clinical_c_mode`: `ring_c_modes`, `map_from_names_64`, `map_from_geometry`, `_sym_blocks`, `_block_circulant`, `_shuffled`, `_load`, `H12C_A2_NPZ`, `H16B_D2_NPZ`, `C_MODE_CU_D1_RECORD_F`, `C_MODE_CU_D2_COMPARAND_F`; `_study`, `_phantom`, `FEEDS`, `DRIVES` (`tests/package/test_wf16_study_tune_f_human_copper.py:104`, `:100`, `:81–82`); `FREQUENCY_HZ`, `DEGREE` (`tests.validation.test_wf7_f_human_b1_quadrature`); `ESR_OHM` (`tests.validation.test_wf10_f_human_copper_esr_root`); `_bitwise_equal_across_ranks` (`tests/post/test_point_owner_rule.py:84`); the three bands. Copied, not imported: PORT-24's `pre` loader and `_canon`, pinned **`PRE_CHANGE_SHA = "9b90d7b3"`** (assert `"n_slots" not in src` and `"esr_ohm" in src`).
- **Nodes:** `test_wf16_s13_n_slots_coding` (ungated, no mesh; W0) and `test_wf16_s13_default_bitwise_records` (`FEM_EM_WF16_S13=1`; W1). Every `build` stubbed to raise; 0 calls asserted in both.
- **Registry (lines added only):** both nodes under `[names.Study]` and `[names.FHumanBirdcageGeometry]`, with `# WF-16 step 13 (§9 item 588)`. The §9 standing rule applies: Wgen / Wval / Wapi / Wdoc in this slot, the regenerated `docs/guide/api/` and `docs/guide/validation.md` committed; Wval's table rows grow by 2 (print the count).

**Anchors (asserted; bands imported by name).**
- **(N16 / N32) the weekly's "equal to C10a's `ring_c_modes` … to 0"** (W0): on (i) C10a's (K1) synthetics (`_shuffled(_block_circulant(_sym_blocks(N, ω, 1.0731e-10), N), N, seed=581 + N)`, N = 16, 32), (ii) C10a's (K2) C16-broken N = 16 matrix (seed 5810 / 17, as `:286–290`) and its N = 32 analogue (seed 5820 / 33), and (iii) C3a's ideal-layout 64×64 (`c3._synthetic_s(*c3._ideal_layout())`, `tests/validation/test_geo38_f_clinical_matrix.py:129`, `:136`, on C10a's geometry map, as `:299–308`): `ring_mode_rows(s, port_map=…, n_slots=N, freq_hz=FREQUENCY_HZ, c_f=1.0731e-10, feeds=[(0, "bottom"), (N // 4, "bottom")], drives={"x": (1.0, 0.0)}, z0=50.0)` against `ring_c_modes(s, …, n_slots=N, z0=50.0)` — λ_anti for k = 0…N/2, `c_mode_f`, `c_k2_f`, `off_block`: **worst relative difference == 0.0** (the same expressions in the same order; C10a's (K2) read 0.000e+00). `len(modes) == N // 2 + 1`.
- **(U2) the U2 closed form at `REUSE_REPRODUCTION_RTOL`** (W0): on (i), \|C_mode/c₀ − 1\| and \|C₂ · ω · 52 − 1\| ≤ `REUSE_REPRODUCTION_RTOL`; anti(k = 1) ≥ 1 − `REUSE_REPRODUCTION_RTOL`. Printed against `EXACT_IDENTITY_RTOL` (predicted ≤ 2e-15: C10a `:65–66`).
- **(T32) `Study.tune(n_slots=32)`** (W0): `_study(geom=FHumanBirdcageGeometry(leg_count=32, as_hole=True))`, with `result` injected as the shuffled N = 32 synthetic and port ids `f"P{o}"` for `geometry.ring_ordinals()` in the same shuffle. `tune(k=1, n_slots=32)` and `tune(k=2, n_slots=32)` meet (U2) at `REUSE_REPRODUCTION_RTOL`; each equals `ring_c_modes(s, map_from_names_64(ids), n_slots=32, …)`'s `c_mode_f` / `c_k2_f` exactly (`==`); `len(z_k) == 17`; `tune(k=1)` (default ⇒ `leg_count`) `==` `tune(k=1, n_slots=32)`; `_bitwise_equal_across_ranks` on the two C values (`-n 2`).
- **(X) refusals naming `n_slots`, no build** (W0): `ring_mode_rows` with `n_slots` ∈ {2, 3, 15, `True`, 16.0} ⇒ `ValueError` / `TypeError`; the 64-port with `n_slots=16` ⇒ `ValueError` (as C10a `:310–311`); `ring_port_map_from_names(["P17"], n_slots=32)` ⇒ `ValueError`; the 32-leg Study with `n_slots=16` ⇒ `ValueError`, `st.result` untouched; a four-port PEC Study and a sphere Study with `n_slots=16` ⇒ `TypeError` with `st.result is None`. The existing refusals (`test_wf16_s12_k0_refusals_and_synthetic`, in W0) green and unedited.
- **(B-rm) default bitwise against the pinned pre-change source** (W1): on `H12C_A2_NPZ` and `H16B_D2_NPZ`, with `kw = dict(port_map=pre.ring_port_map_from_names(order), freq_hz=FREQUENCY_HZ, feeds=FEEDS, drives=DRIVES)`, for c ∈ {`C_MODE_CU_D1_RECORD_F`, `C_MODE_CU_D2_COMPARAND_F`} and r ∈ {default, `ESR_OHM`}: `_canon(rm.ring_mode_rows(s, c_f=c, **kw[, esr_ohm=r]))` `==` `_canon(pre.ring_mode_rows(…))`, and again with `n_slots=16` explicit on the new side; `rm.ring_port_map_from_names(order) == pre.ring_port_map_from_names(order)`; per k, `np.array_equal(rm.k_blocks(zc)[0][k], pre.k_blocks(zc)[0][k])` and the off-blocks `==`.
- **(B-tune) "the two recorded F-human `C_mode` reproduced bitwise through the default"** (W1): `_study()` (16 legs, copper hole) with each npz's S and `order` injected; `tune(k=1)` and `tune(k=2)`, default `n_slots`, against the pre-change `_study_tune_ring` body (`study.py:1140–1150`, copied into the test with `pre.`'s functions and the Study's own `z0_ohm` / `omega`): `c_mode_farad ==`, `np.array_equal(z_k, lam)`, `off_diagonal_residual ==`. **The recorded prints reproduce as strings:** `f"{C:.12e}"` is `"1.181185799070e-10"` (d1; `20261015T213534Z_WF-16-s12-W1.log:7375`, C10a `:73`) and `"1.145221310828e-10"` (d2; C10a `:83`, H25 `:171`), and each is within `REUSE_REPRODUCTION_RTOL` of its constant (C10a read 3.515e-13 / 4.212e-13).

**Negative controls (rule (e)).**
- **(c-k2) asserted, backed:** on the N = 32 synthetic, `tune(k=2, n_slots=32)` misses c₀ by more than `REDUCTION_BAND`. Backing: C10a measured C₂ = 1/(ω·52) on this synthetic to 1.110e-15 (`20261021T093341Z_GEO-38-C10a-W1.log:65`), so the miss is \|1/(ω·52·c₀) − 1\| — **predicted 0.554**, ≈ 550× the bar; a ratio has no ceiling here.
- **(c-bit) asserted, backed:** the (B-rm) comparator can see a change — the new code at r = `ESR_OHM` against the pre-change code at r = 0 gives `_canon` ≠ (`PORT-24` measured the shift Re λ_k(r) − Re λ_k(0) = r on F-human, `20261011T093040Z_PORT-24-step1-T25p.log:96`).
- **(c-N) asserted (a raise, no factor):** a 32-leg read at N = 16 is refused, never computed — (X).
- **Printed only:** the top ring offset by one slot at N = 32 leaves `C_mode` and drops anti(k = 1) to cos²(π/32) ≈ 0.9904 (C10a read 0.961945 at N = 16).

**Predictions (printed, never asserted).** (N16 / N32) 0.0 on every quantity; (U2) ≤ 2e-15; (T32) ≤ 2e-15 vs closed form; (B-rm) / (B-tune) `True`; the re-gates' readings equal their last logs' printed digits (no-solve gates) or sit inside their own bands (re-solves; e.g. `test_wf16_s12_a` printed gap −3.544942e-13, `:7375`).

**Windows.** Foreground; **the harness call is the whole Bash command: bare, repo-relative, nothing chained or piped**; Bash `timeout: 660000`; pytest `-s`; node IDs, never `-k`. `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0 before RG-T, RG-Tcap, RG2, RG5, RG-Td2, RG-Td2cu, RG-Tfh and RG4. W0 runs only after every edit (rule (i)). No window needs durable capture (every ceiling ≤ 590 s).

| Window | Price | Source |
|---|---|---|
| W0 (`-n 2`) / W1 (`-n 2`) | ≈ 6 s / ≈ 6 s | C10a W0 5 s (`20261021T093330Z_GEO-38-C10a-W0.log`); T25p body 8 s |
| RG-U1 (PORT-24, ring rows, private writer, circuit layer, API reference; `-n 1`) | 10 s | `20261011T093040Z_PORT-24-step1-T25p.log` |
| RG-U2 (ring rows, private writer, v3 rows / definition; `-n 1`) | 84 s | `20261013T094942Z_ANS-8-item411-W1.log` |
| RG-C10a / RG-H25 / RG-1c (`-n 1`) | 96 / 96 / 74 s | `20261021T093341Z_GEO-38-C10a-W1.log`, `20261018T095559Z_ANS-8-H25-W1.log`, `20261011T123443Z_ANS-8-step1c-W1.log` |
| RG-Y0 (`tune:` and every parse node) / RG-page | 4 / 5 s | `20261020T203852Z_WF-1-s21-W0-r2.log`, `20261015T214647Z_WF-17-s6-W0.log` |
| RG-T / RG-Tcap (`-n 2`) | 69 / 69 s | `…213242Z_WF-16-s12-RG1-tune.log`, `20261015T201110Z_WF-16-s6b-W1.log` |
| RG2 / RG5 (`-n 2`) | 61 / 41 s | `20261020T202817Z_WF-1-s21-RG2.log`, `…203329Z_WF-1-s21-RG5.log` |
| RG-Td2 / RG-Td2cu (`-n 8`, ceiling 400) | 82 / 83 s | `20261014T111026Z_WF-16-s6c-W1.log`, `20261014T185251Z_WF-16-s6d-W1.log` |
| RG-Tfh / RG4 (`-n 8`, ceiling 590) | 205 / 206 s | `20261015T213534Z_WF-16-s12-W1.log`, `20261020T202959Z_WF-1-s21-RG4.log` |
| Wgen / Wval / Wapi / Wdoc | 1 / 3 / 20 / 2 s | `20261020T203425Z` / `203429Z` / `203434Z` / `203455Z_WF-1-s21-*.log` |

Windows ≈ 20.5 min.
```
scripts/testing/run_and_log.sh WF-16-s13-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf16_study_tune_n_slots.py tests/package/test_wf16_study_tune_f_human_copper.py::test_wf16_s12_k0_refusals_and_synthetic -s -v --tb=short'"
scripts/testing/run_and_log.sh WF-16-s13-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_WF16_S13=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf16_study_tune_n_slots.py -s -v --tb=short'"
```
Every RG window is the `- Command:` line of the log in the table, **verbatim**, under the label `WF-16-s13-<RG name>`. RG-Y0 is item 572's W0-r2 list with `tests/package/test_wf1_yaml_copper_em_probes.py::test_wf1_s22_parse` appended (item 573's parse node). RG-C10a never sets `FEM_EM_GEO38_C10`.

**Traps already paid for.**
- **Bitwise is call by call on the same inputs** (rubric, item 537): (B-tune) compares the Study against the pre-change body on the **same** S in the **same** port order. A shuffled-vs-unshuffled S differs by 1 ulp through `s_to_z` (`test_wf16_study_tune_f_human_copper.py:192–199`, `20261015T213212Z_WF-16-s12-W0.log:63`).
- **Pin the pre-change commit, never `HEAD:`** (rubric, `OPS-43` (a)): `9b90d7b3`. The source assertions catch a wrong pin.
- **`REUSE_REPRODUCTION_RTOL` has two definitions**; import `test_port_birdcage_ring_matrix.py:171`'s.
- **At N = 32, k = 1 is bins 1 and 31,** and the mode loop runs to 16 (C10a trap); the name route depends on N: P33…P48 are *top* at N = 16 and *bottom* at N = 32.
- **The npz are gitignored and never re-solved;** absent ⇒ BLOCKED. Never write under `logs/wf12/`.
- **A Python `float` `n_slots` (16.0) is refused, not cast** — a silent cast is the defect class the 10-21 weekly just ruled on (U75).
- `-k a or b` splits argv; `-k 30` on every container `timeout`; never pipe; never import a `test_*` name; `capsys.readouterr()` discards earlier prints; status 124 or a kill ⇒ `pgrep` first, then only the rubric's orphan-stem sweep labelled `WF-16-s13-cache-sweep`. Never raise a timeout or change `-n`. Write ≈, never a tilde.

**Scope.** Covered: `n_slots` on `ring_mode_rows` and `Study.tune`, with the default proven bitwise on the two copper F-human records and on `PORT-24`'s σ = 800 record (RG-U1). The N = 32 path is proven on synthetics (closed form, C10a's reader) and on C3a's ideal-layout 64×64. **Not claimed:** any F-clinical `C_mode`, node grid, root or tuning figure (C10 still reads the C10a reader, as the weekly ruled); drive rows, `k_spectrum` folds and the feed 2×2 at N ≠ 16 (computed, not anchored); `ring_mode_rows_v2` at N ≠ 16; `n_slots` from YAML; a 32-leg F-human build; any AED value.

**Negative result.** Report and stop. Park on `attempt/WF-16-step13-<UTC>` and mark item 588 BLOCKED in the same commit (rule (d)).
- **(B-rm) / (B-tune) red:** the default is not bitwise — revert, print the first differing key; never relax to a band.
- **(N16 / N32) ≠ 0 with (U2) green:** the two codes differ in operation order. Print the worst key; the review rules (0 is the weekly's letter).
- **(U2) red:** an indexing or normalisation error at N ≠ 16; stop.
- **Any RG red:** revert `src/` and stop. **Status 124:** not a result; park with the measured time.

No `PREDICTED RED:` line: a red anywhere stops the slot.

**Status it can move.** `WF-16` stays ✅ (rule (j): its Done-when is the four-port S record; no glyph moves). Its state line gains a **step-13 record:** "`n_slots` on `ring_mode_rows` / `Study.tune`: default bitwise against `9b90d7b3` on the d1 / d2 copper records (`C_mode` 1.181185799070e-10 / 1.145221310828e-10 F); equal to C10a's reader at N = 16 / 32 to <x>; U2 closed form ≤ <u>; `src`: `ring_modes.py`, `study.py`, keyword insertions." §10 U74 DONE. The `ring_modes` consumer rows (`ANS-8`, `PORT-24`, `GEO-38`) do not move. The 10-21 weekly's "a `ports:` example after U74 if an F-human record shows N-slot tuning at standard tier" (`PROJECT_PLAN.md:5081–5083`) is **not** opened by U74 alone: no N ≠ 16 F-human record exists.

**For the review.**
1. **`Study.tune`'s default.** The drafter takes `n_slots=None` ⇒ `geometry.leg_count`, refusing a mismatch by name. That makes every 16-leg call bitwise and turns a 32-leg Study with a stray `n_slots=16` into a refusal, not a silent wrong read. The weekly's "(additive, default 16)" is read as `ring_mode_rows`'s default. The alternative is a plain `n_slots=16` on `Study.tune` with no tie to the geometry.
2. **Beyond the letter, disclosed:** `ring_port_map_from_names`, `_canonical_permutation`, `_dft_matrix`, `k_blocks` and `k_spectrum` gain the same keyword. `Study.tune` cannot read a 64-port without the name route. `ring_mode_rows_v2` is **not** given `n_slots`.
3. **(U2)'s band** is the weekly's `REUSE_REPRODUCTION_RTOL` (1e-6), although C10a asserted the same closed form at `EXACT_IDENTITY_RTOL` and read ≤ 1.6e-15. The drafter asserts the letter and prints against the tighter band.
4. **"To 0"** is asserted as exact equality. It holds only because both codes run the same expressions in the same order. A future refactor of either breaks it without a physics change.
5. **Examples not re-run:** `mri:19`, `21`, `23`, `24`, `26` call `Study.tune` (items 564 / 572 re-ran gates, not examples). `mri:23` is the one ring-route example.
6. **Length:** ≈ 63 predicted minutes against a one-hour slot. If the review wants margin, split 588a (code, W0, W1, the no-solve RGs, page gates) from 588b (the six solve re-gates, ≈ 12 min).
7. **Filed in `WF-1.md`** as assigned, though it is `WF-16` step 13. The review may move it to `WF-16.md`.

**Cost.** Lane-2 feature step. **≈ 14.6 measured** (predicted × 0.231) / **≈ 63 predicted** slot-minutes: 15 fixed + ≈ 10 `src` + ≈ 12 module + windows ≈ 21 + commit and journal ≈ 5. The weekly's supply line carries 16.7. **Executor:** implementer. **Tier:** standard (W0 / W1, the `-n 1` / `-n 2` RGs); heavy by ceiling for the four `-n 8` re-gates. **Ranks:** `-n 2` (≤ 8 in the re-gates, as landed). **Daytime:** no. **Depends on:** nothing. **Family:** `WF-16` step 13, attempt 1 of 4. **Park (rule (d))** if code + W0 + W1 are not green by minute 30, or RG-Tfh has not started by minute 48.

**Execution order.** Host Glob `logs/wf12/h12c-cu-s32-a2.npz`, `logs/wf12/h16b-cu-s32-d2.npz`, `logs/ans8/s32.npz` (absent ⇒ BLOCKED) → code → module → W0 → W1 → RG-U1 → RG-U2 → RG-C10a → RG-H25 → RG-1c → RG-Y0 → RG-page → RG-T → RG-Tcap → RG2 → RG5 → RG-Td2 → RG-Td2cu → RG-Tfh → RG4 → Wgen → Wval → Wapi → Wdoc → commit 1 (`src`, module, registry, regenerated pages, every `WF-16-s13-*` log with its test-results rows, `PROJECT_PLAN.md` §7 `WF-16` state line, §9 item 588 DONE, §10 U74; this file's outcome paragraph) → commit 2 (journal, `docs/testing/attempts.md` only).

**Executor spawn prompt (verbatim).** "Execute §9 item 588 — `WF-16` step 13 (§10 U74, lane 2) — per `docs/planning/chunks/WF-1.md` § "§9 item 588". Read the brief and any review ruling beneath it first; the rulings hold. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. `src/` change: `ports/ring_modes.py` and `study.py` only — the `n_slots` keyword as the brief lists it, default bitwise; no `study_yaml.py` / `cli.py` edit, no `ring_mode_rows_v2` edit, no YAML schema change. Write the module whole, then run W0, W1, RG-U1, RG-U2, RG-C10a, RG-H25, RG-1c, RG-Y0, RG-page, RG-T, RG-Tcap, RG2, RG5, RG-Td2, RG-Td2cu, RG-Tfh, RG4, Wgen, Wval, Wapi, Wdoc in order, foreground, every harness Bash call with `timeout: 660000`, pytest with `-s`, container `timeout -k 30`. `pgrep -c python3` must read 0 before every solve re-gate. Every RG command is copied verbatim from the log the brief names. Glob the three npz first; if any is absent, stop and mark item 588 BLOCKED; never re-solve and never write under `logs/wf12/`. Never set `FEM_EM_GEO38_C10`. Never import a `test_*` name; never edit an existing test assertion. Never open `aed_results/`, `COMPARISON_private.md`, `aed/` or `docs/private/`. You are the executor — do not spawn agents, never background a harness window. If the docker socket is denied, stop and report the exact command string. Park if code + W0 + W1 are not green by minute 30 or RG-Tfh has not started by minute 48. On a red, park on `attempt/WF-16-step13-*`, mark item 588 BLOCKED in the same commit, and report."

**Outcome (2026-10-22, 11:00 UTC slot, second item): DONE, first try.** `src`: `ports/ring_modes.py` (`_check_n_slots`; `n_slots` on the name route, `_canonical_permutation`, `_dft_matrix`, `k_blocks`, `k_spectrum`, `ring_mode_rows`; `ring_mode_rows_v2` untouched) and `study.py` (`n_slots` on `Study.tune` / `_study_tune_ring`, default `leg_count`, mismatch refused before any solve; `TuneResult` docstring). New `tests/package/test_wf16_study_tune_n_slots.py`; registry +2 nodes under `[names.Study]` and `[names.FHumanBirdcageGeometry]`. W0 (`20261022T111700Z_WF-16-s13-W0.log`, 6 passed 1 skipped): (N16 / N32) worst rel 0.000e+00 on all five cases; (U2) 1.443e-15 / 6.661e-16 (N = 16), 1.554e-15 / 1.110e-15 (N = 32); (T32) tune = reader `==`, len(z_k) 17, bitwise across 2 ranks; (c-k2) 0.5543 (predicted 0.554); offset print anti 0.990393 = cos²(π/32); (X) all refusals, 0 builds. W1 (`20261022T111711Z_WF-16-s13-W1.log`, 6 passed): (B-rm) 8/8, name route and `k_blocks` equal, (c-bit) differs ×4; (B-tune) `1.181185799070e-10` / `1.145221310828e-10` strings equal, 3.515e-13 / 4.212e-13 vs constants. Re-gates all green at their landed pass counts: RG-U1 43 (9 s), RG-U2 36 (97 s), RG-C10a 14+1s (96 s), RG-H25 14 (97 s), RG-1c 14 (96 s), RG-Y0 22 (4 s), RG-page 11 (5 s), RG-T 10 (65 s), RG-Tcap 9 (62 s), RG2 10 (59 s), RG5 10 (40 s), RG-Td2 9 (80 s), RG-Td2cu 9 (73 s), RG-Tfh 9 (202 s; `s12 a` gap −3.549383e-13, inside `REDUCTION_BAND`), RG4 9 (205 s); Wgen 46 files, Wval 330 nodes (+2 rows per class table), Wapi 22 passed, Wdoc dead=0.

### §9 item 589 — 2026-10-22 (§10 U75, chain step (no solve))

*Drafted read-only for the 2026-10-22 03:00 daily review; the review rules on it.* Code read on `11fc13b0`; nothing was run.

**Licence.** §10 chain U, U75 (`PROJECT_PLAN.md:4953–4957`): "`Study.solve_sweep` refuses a Study with a conductor by a named field (the YAML refusal's twin, same message), gated on both routes; `femem check` lists `sweep_hz` on a sweep file (the key printed equals the file's); every sweep gate re-run (U71's records bitwise). Serial on nothing." The 10-21 weekly posed it from two defects found by reading, not by a run (`PROJECT_PLAN.md:4819–4823`): "`Study.solve_sweep` ignores `study.conductor` (the YAML already refuses) … `femem check` on a sweep file prints `frequency_hz`, not `sweep_hz` (U75)." Sources: this file's review ruling 2026-10-20 10:30 (2) (`:3019`, the YAML refusal "does not fix the Python path; that fix is the weekly's to pose") and (4) (`:3021`, `_run_palace` — U76's, not here); `EX-139.md` ruling 18:00 Q5 (`docs/planning/chunks/EX-139.md:214`; the question at `:210`). **No band, tolerance, record, Done-when or target is set or moved; every band is imported by name** (the re-gates import their own). This is **`WF-1` step 23** (step 22 is item 573), with a `WF-16` refusal. **Class: chain step (no solve)** — a defect closure, not a feature. The new nodes build nothing; the re-gates solve as landed. The weekly priced it at the chain no-solve median (`PROJECT_PLAN.md:5012`).

**What exists (`11fc13b0`).**
- `_study_solve_sweep` (`src/fem_em_solver/study.py:963–1043`) refuses a non-four-port geometry and a non-`LumpedPortDrive` excitation (`TypeError`, `:977–980`), validates the list (`:981`), then builds. It never reads `self.conductor`, and its `material_map` always names `geom.conductor_tag` (`:994–998`). **A copper Study sweeps as the σ = 800 S/m proxy, silently.**
- The YAML refusal (`src/fem_em_solver/study_yaml.py:1228–1231`, field `study.sweep_hz`): reason `"a sweep with a conductor is not supported (Study.solve_sweep has no conductor route)"`. It fires before `as_hole` (`:1232–1233`). `StudyConfigError` (`:175–182`) stores `.field` and `.reason`.
- `femem check` (`src/fem_em_solver/cli.py:162–186`): `study, kind = _load(args.file)[:2]` (`:164`). On EM files rank 0 prints six lines — `family`, `degree`, `frequency_hz`, `material.sigma`, `material.epsilon_r`, `excitation` (`:177–185`). `_load` returns `(study, kind, outputs, cooling)` for every family (`study_yaml.py:873`, `:1006`, `:1357`); `outputs["sweep_hz"]` is the parsed float tuple (`:1344–1347`).
- Gates on `check` output: `tests/package/test_wf1_check_family.py` — (F) one `family:` line per good file, which `f_small_sweep.yaml` keeps; (U) a six-line equality on two files without a sweep — against a pinned `cli.py`. The in-container pinned-source pattern: `tests/unit/test_port24_esr_circuit.py:67–79`. Last `cli.py` change: `acf99b1b`.
- `examples/mri/34_femem_run_sweep_file_touchstone.md:41–43` states "**`femem check` does not list the sweep.** … Whether `check` should list the sweep is a question for the weekly" — **false once this lands** (For the review 3).
- Sweep gates: `tests/package/test_wf16_study_sweep.py` (degree 1 `-n 2`, `20261006T033609Z_WF-16-s2-W1-deg1.log.gz`, 199 s; degree 2 `-n 8`, `20261006T034120Z_WF-16-s2-W2-deg2.log.gz`, 303 s); `tests/package/test_wf1_yaml_sweep.py` (W0 parse, W1a `20261020T201454Z_WF-1-s21-W1a-sweep-file.log` 53 s, W1b `20261020T201555Z_WF-1-s21-W1b-sweep-python.log` 53 s). The `solve_sweep` refusal nodes: `test_wf16_study_sweep.py::test_wf16_sweep_api_and_refusals`, `test_wf16_study_capacitor_sheets.py::test_wf16_s10_api_and_refusals` (`:128`), `test_wf16_study_f_human.py::test_wf16_f_human_w0_kwargs_tags_refusals` (`:165`); last run together at `-n 1` in `20261011T201926Z_WF-16-s10-W0.log` (76 s with the full `tests/environment`).

**Code (`src/`: `study.py` and `cli.py`; insertions plus one disclosed edited line, rule (c)).**
- **`study.py`:**
  - (P1) a module constant beside `_validate_sweep_frequencies`: `_SWEEP_CONDUCTOR_REASON = "a sweep with a conductor is not supported (Study.solve_sweep has no conductor route)"` — the YAML reason, byte for byte.
  - (P2) in `_study_solve_sweep`, inserted after the excitation check (`:980`) and **before** `_validate_sweep_frequencies` and `build`: `if self.conductor is not None:  # WF-1 step 23 (§9 item 589)` / `raise TypeError(f"conductor: {_SWEEP_CONDUCTOR_REASON}")`. This is the named field: `TypeError`, as the two refusals above it. One docstring clause.
- **`cli.py`:** the **edited line** `:164` becomes `study, kind, outputs = _load(args.file)[:3]`. After the `frequency_hz:` print (`:181`), inside the rank-0 block, insert `if "sweep_hz" in outputs:  # WF-1 step 23 (§9 item 589)` / `print(f"sweep_hz: {list(outputs['sweep_hz'])!r}")`. The module docstring gains one clause. The thermal-only branches are untouched; they return before the print.
- **Never:** edit `study_yaml.py` (the YAML message stays a literal; the twin is proven by equality), `_run_palace` (U76), the sweep's solve path, `__all__`, any existing test assertion.

**Test module: new `tests/package/test_wf1_sweep_conductor_check.py`** (first docstring paragraph backticks only `WF-1`, `WF-16`). It has no knob: every node builds nothing.
- **Imports (names only; never a `test_*` name):** `SWEEP`, `FOUR_PORT`, `GOOD_SWEEP`, `_with_sweep`, `_field` (`tests/package/test_wf1_yaml_sweep.py:70`, `:71`, `:111`, `:107`, `:99`); `_port_study`, `FREQUENCIES` (`tests/package/test_wf16_study_sweep.py:112`; `tests/validation/test_ans4_exact_strip_records.py:135`); `WF1`, `LOGS`, `_mapping`, `_fresh_dir`, `_rc_all` (`tests/package/test_wf1_study_yaml.py:85–86`, `:157`, `:163`, `:169`); `_study_from_mapping`, `_load` (`fem_em_solver.study_yaml`); `_SWEEP_CONDUCTOR_REASON`, `COPPER_SIGMA`, `LeontovichWall`, `Study`, `FourPortBirdcageGeometry`. Copied, not imported: PORT-24's in-container `git show` loader, pinned **`PRE_CHANGE_CLI_SHA = "acf99b1b"`** (assert `"sweep_hz" not in src`), loaded as `fem_em_solver._cli_acf99b1b`.
- **Nodes:**
  - `test_wf1_s23_python_sweep_refuses_conductor`: (Py), (c-reach).
  - `test_wf1_s23_yaml_twin_and_cli`: (Y).
  - `test_wf1_s23_check_lists_sweep_hz`: (C), (U), (c-check).
- **Registry (lines added only):** the first node under `[names.Study]` and `[names.LeontovichWall]`; the second under `[names.load_study]`, `[names.StudyConfigError]` and `[names.run_study_file]`; the third under `[names.load_study]`; each block with `# WF-1 step 23 (§9 item 589)`. The §9 standing rule applies: Wgen / Wval / Wapi / Wdoc in this slot; Wval's table rows grow by 3 (print the count).

**Anchors (asserted; all in W0, no mesh; every `build` stubbed to raise).**
- **(Py) the Python route refuses by name.** Two copper Studies: the hole, `load_study(WF1 / "f_small_copper.yaml")`; and the solid, `Study(g, m, e, frequency_hz=…, degree=1, conductor=LeontovichWall(COPPER_SIGMA))` from `_port_study(1)[0]`'s parts. On both, `solve_sweep(tuple(FREQUENCIES.values()))` raises `TypeError` with `str(exc) == "conductor: " + _SWEEP_CONDUCTOR_REASON`. Also: 0 build calls, `not hasattr(st, "sweep_result")`, `st.result is None`. With a bad list (`()`), the conductor refusal still fires first (printed). F-human and capacitor Studies keep their earlier `TypeError` (the three refusal nodes in RG-R, unedited).
- **(Y) the YAML twin, same message.** On `_with_sweep(_mapping(WF1 / "f_small_copper.yaml"), GOOD_SWEEP)`: `_field(lambda: _study_from_mapping(m)) == "study.sweep_hz"`, and the caught `StudyConfigError.reason == _SWEEP_CONDUCTOR_REASON` (string `==`). The same mapping written by rank 0 to `LOGS / "s23" / "copper_sweep.yaml"` (`_fresh_dir`, `Barrier`): `cli.main(["check", f])` and `cli.main(["run", f, "--output-dir", …])` both return rc (2, 2) on every rank (`_rc_all`), the output directory absent on every rank.
- **(C) `femem check` lists the sweep, and the key printed equals the file's.** `cli.main(["check", SWEEP])` rc (0, 0). On rank 0, exactly one line starts with `"sweep_hz: "` and it sits immediately after the `frequency_hz:` line. Its key is `sweep_hz`, which is a key of `yaml.safe_load(SWEEP.read_text())["study"]`. `ast.literal_eval` of its value equals that YAML list element by element by `float.hex`, and equals `tuple(FREQUENCIES.values())`.
- **(U) every other line unchanged against the pinned `cli.py`.** Files: every `tests/package/wf1/*.yaml` not prefixed `malformed_`, plus every `examples/mri/*.yaml`. Rc is (0, 0) on both sides. Where the file has no `study.sweep_hz`, new lines `==` old lines. On `f_small_sweep.yaml`, new lines `==` old lines with the one (C) line inserted. Print the file count and the number of sweep files (predicted 1).

**Negative controls (rule (e)).**
- **(c-reach) asserted, backed — the refusal is conductor-specific:** the same F-small Study without a conductor (`_port_study(1)[0]`), with `build` stubbed to raise a sentinel, reaches `build` (1 call). Backing: U71's W1a and W1b meshed and solved this Study's sweep (`20261020T201454Z_WF-1-s21-W1a-sweep-file.log:1934`, mesh 25.28 s). A blanket refusal would fail here.
- **(c-check) asserted (a line present or absent, no factor):** the pinned `acf99b1b` `check` on `f_small_sweep.yaml` prints no `sweep_hz:` line. This is the reading `EX-139.md:100` / `:210` recorded from `cli.py:177–185`.
- No negative control exists for "the old route silently ignored the conductor" without a copper solve (≈ 60 s at `-n 2`). The defect is read from source (`study.py:963–1043`), as the weekly posed it. It is not measured here (For the review 4).

**"U71's records bitwise" — read as follows (For the review 2).** The re-gates must reproduce U71's printed record digits **as strings**:
- (R′) 6.601e-11 / 9.557e-11 / 9.449e-11, cells 117206 (`…W1a-sweep-file.log:1934`);
- (c-swap) 1.2367 / 1.2763 (`:1935`);
- the parse's `float.hex` equality (W0).

(I-Py) is a two-process comparison and has no bitwise contract: it read 9.010e-15 in W1b (`…W1b-sweep-python.log:1934`) and 6.328e-15 in `EX-139`'s run. It stays at `L0_BAND` and is printed beside.

**Windows.** Foreground; **the harness call is the whole Bash command: bare, repo-relative, nothing chained or piped**; Bash `timeout: 660000`; pytest `-s`; node IDs, never `-k`. `pgrep -c python3` (`docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3`) must read 0 before RG-W1a, RG-W1b, RG-S2d1 and RG-S2d2. W0 runs only after every edit (rule (i)).

| Window | Price | Source |
|---|---|---|
| W0 (`-n 2`: the new module + every parse / `check` node) | ≈ 8 s | item 572 W0-r2, 4 s (`20261020T203852Z_WF-1-s21-W0-r2.log`) |
| RG-R (the three `solve_sweep` refusal nodes, `-n 1`) | ≈ 8 s | `20261011T201926Z_WF-16-s10-W0.log` (76 s with the full `tests/environment`; here `test_complex_mode.py` only) |
| RG1 / W9 (installed `femem check`) | 3 / 15 s | `20261020T201312Z_WF-1-s21-RG1.log`, `20261020T201320Z_WF-1-s21-W9-installed-femem.log` |
| RG-W1a / RG-W1b (`-n 2`) | 53 / 53 s | the U71 logs above |
| RG-S2d1 (`-n 2`, ceiling 400) / RG-S2d2 (`-n 8`, ceiling 590) | 199 / 303 s | the `WF-16-s2` logs above |
| Wgen / Wval / Wapi / Wdoc | 1 / 3 / 20 / 2 s | `20261020T203425Z` / `203429Z` / `203434Z` / `203455Z_WF-1-s21-*.log` |

Windows ≈ 11 min.
```
scripts/testing/run_and_log.sh WF-1-s23-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_sweep_conductor_check.py tests/package/test_wf1_yaml_sweep.py::test_wf1_s21_parse_and_refusals tests/package/test_wf1_yaml_copper_em_probes.py::test_wf1_s22_parse tests/package/test_wf1_yaml_em_probes.py::test_wf1_s20_parse_and_refusals tests/package/test_wf1_yaml_thermal_ball_probes.py::test_wf1_s19_parse_and_refusals tests/package/test_wf1_yaml_thermal_box_quadratic_duty.py::test_wf1_s18_parse_and_refusals tests/package/test_wf1_yaml_thermal_box_probes.py::test_wf1_s17_parse_and_refusals tests/package/test_wf1_yaml_thermal_box.py::test_wf1_s15_parse_and_refusals tests/package/test_wf1_yaml_thermal_box_quadratic.py::test_wf1_s16_parse_and_refusals tests/package/test_wf1_yaml_thermal_ball.py::test_wf1_s14_parse_and_refusals tests/package/test_wf1_yaml_capacitors.py::test_wf1_s13_parse_and_refusals tests/package/test_wf1_yaml_tune.py::test_wf1_s11_tune_parse_and_refusals tests/package/test_wf1_yaml_tune_f_human.py::test_wf1_s12_parse_and_refusals tests/package/test_wf1_check_family.py tests/package/test_wf1_yaml_f_human.py tests/package/test_wf1_yaml_f_human_thermal.py -v -s --tb=short'"
scripts/testing/run_and_log.sh WF-1-s23-RGR "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf16_study_sweep.py::test_wf16_sweep_api_and_refusals tests/package/test_wf16_study_capacitor_sheets.py::test_wf16_s10_api_and_refusals tests/package/test_wf16_study_f_human.py::test_wf16_f_human_w0_kwargs_tags_refusals -s -v --tb=short'"
```
RG1, W9, RG-W1a, RG-W1b, RG-S2d1, RG-S2d2, Wgen, Wval, Wapi and Wdoc are the `- Command:` lines of the logs in the table, **verbatim** (the two `.log.gz` read with `zcat`), under labels `WF-1-s23-<name>`.

**Traps already paid for.**
- **Refuse before the frequency list and before `build`:** a conductor Study with a bad list must name `conductor`, and nothing may mesh (0 build calls).
- **The twin is proven by string equality,** not by editing the YAML literal into the constant. `study_yaml.py` stays untouched, so no YAML gate's input moves.
- **`femem` prints on rank 0 only;** `capsys.readouterr()` discards earlier prints (item 562's record defect) — print the (C) / (U) tables after the last capture, through `capsys.disabled()` (`test_wf1_check_family.py:41–45`).
- **Pin the pre-change commit** (`acf99b1b`), never `HEAD:`; the source assertion catches a wrong pin.
- **YAML 1.1 reads `1e7` without a dot as a string;** the copper sweep file reuses `GOOD_SWEEP`'s literals.
- **The record holds at its own width** (`OPS-41`): the degree-1 sweep at `-n 2`, degree 2 at `-n 8`, as landed.
- `-k a or b` splits argv; `-k 30` on every container `timeout`; never pipe pytest; never import a `test_*` name; status 124 or a kill ⇒ `pgrep` first, then only the rubric's orphan-stem sweep labelled `WF-1-s23-cache-sweep`. Never raise a timeout or change `-n`. Rule (i): a `src` edit after W0 re-runs W0 and every RG. Write ≈, never a tilde.

**Scope.** Covered: `Study.solve_sweep` refuses a Study with a conductor by name (hole and solid copper F-small), with the YAML's reason byte for byte, on the Python and the YAML / CLI routes; `femem check` prints `sweep_hz` on a sweep file and nothing new elsewhere; the sweep gates re-run. **Not claimed:** a copper sweep (the fix is a refusal, not a conductor route); `as_hole` without a conductor on the Python route (it reaches `build` and fails loudly on the missing conductor tag, item 324 — not silent, not in the letter); `_run_palace`'s silent drop (U76); the `mri:34` example re-run; any AED statement.

**Negative result.** Report and stop. Park on `attempt/WF-1-step23-<UTC>` and mark item 589 BLOCKED in the same commit (rule (d)).
- **(Py) / (Y) red:** the message or the order differs; fix before the re-gates only.
- **(U) red:** a non-sweep file's `check` output moved. Revert `cli.py` and stop.
- **(c-reach) red (the conductor-free Study is refused):** the refusal is too broad; stop.
- **Any RG red:** revert `src/` and stop. **Status 124:** not a result; park with the measured time.

No `PREDICTED RED:` line: a red anywhere stops the slot.

**Status it can move.** `WF-1` stays ✅ (rule (j): no glyph). Its state line gains a **step-23 record:** "`Study.solve_sweep` refuses a conductor by name (`conductor: <reason>`, the YAML reason byte for byte, hole and solid copper, 0 builds); `femem check` prints `sweep_hz` on a sweep file, every other file's lines equal `acf99b1b`'s; U71's printed records reproduced; `src`: `study.py` (insertions), `cli.py` (insertions plus one edited line)." Also:
- the 10-21 weekly's silent-drop bullet (`PROJECT_PLAN.md:4819–4823`) loses its `solve_sweep` and `check` halves; `_run_palace` stays U76's;
- `EX-139`'s Q5 is answered;
- §10 U75 DONE.

**For the review.**
1. **The Python refusal's shape:** `TypeError("conductor: " + reason)`, as `solve_sweep`'s two refusals above it. The alternative is `ValueError`, or a `StudyConfigError`-like field object. The YAML literal is not rewired to the constant, which keeps `study_yaml.py` untouched.
2. **"U71's records bitwise"** is read as U71's printed digits, string-equal, plus the parse's `float.hex`. A re-solve's S has no bitwise contract across processes ((I-Py) moved 9.010e-15 → 6.328e-15 between W1b and `EX-139`).
3. **`mri:34`'s guide becomes false** (`examples/mri/34_femem_run_sweep_file_touchstone.md:41–43`). Options:
   - (a) this item rewrites that paragraph to state the new line, citing U75, and runs the docrefs check — `scripts/testing/check_example_doc_references.py` over the examples root, `dead=0` — with no script edit and no example re-run;
   - (b) a separate example-lane item;
   - (c) leave it and journal it.

   The drafter's default is (a); it edits a ✅ example's page after its window, so it is the review's call.
4. **No measured negative for the old silent path.** Measuring it needs a copper sweep solve; the drafter does not propose one.
5. **RG-S2d2** (degree 2, `-n 8`, 303 s) is included under "every sweep gate". Dropping it saves ≈ 5 min, and the diff cannot reach the solve path.

**Cost.** Chain step (no solve). **≈ 8.0 measured** (class median) / **≈ 45 predicted** slot-minutes: 15 fixed + ≈ 5 `src` + ≈ 10 module + windows ≈ 12 + commit and journal ≈ 3. **Executor:** implementer. **Tier:** standard; heavy by ceiling for RG-S2d1 (400 s) and RG-S2d2 (590 s). **Ranks:** `-n 2` (`-n 1` for RG-R, `-n 8` for RG-S2d2, as landed). **Daytime:** no. **Depends on:** nothing. **Family:** `WF-1` step 23, attempt 1 of 4. **Park (rule (d))** if code + W0 are not green by minute 25.

**Order against item 588 — either order.** The two items are independent.
- **`study.py` hunks are disjoint:** 588 edits `_study_tune` / `_study_tune_ring` / `TuneResult` (`:1054–1150`); 589 inserts in `_study_solve_sweep` and beside `_validate_sweep_frequencies` (`:947–1043`).
- **No shared module or re-gate:** neither item's new module imports what the other edits, and neither re-gate set reaches the other's hunk.
- **The one shared file is `validation_envelope.toml`** (both add lines under `[names.Study]`). Whichever lands second rebases, keeps both sides, and re-runs Wgen / Wval / Wapi / Wdoc before its commit (the §9 standing rule).

**Execution order.** Code → module → W0 → RG-R → RG1 → W9 → RG-W1a → RG-W1b → RG-S2d1 → RG-S2d2 → Wgen → Wval → Wapi → Wdoc → commit 1 (`src`, module, registry, regenerated pages, every `WF-1-s23-*` log with its test-results rows, `PROJECT_PLAN.md` §7 `WF-1` state line, §9 item 589 DONE, §10 U75; this file's outcome paragraph; plus the guide paragraph only if the review rules (3)(a)) → commit 2 (journal, `docs/testing/attempts.md` only).

**Executor spawn prompt (verbatim).** "Execute §9 item 589 — `WF-1` step 23 (§10 U75, chain step) — per `docs/planning/chunks/WF-1.md` § "§9 item 589". Read the brief and any review ruling beneath it first; the rulings hold. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. `src/` change: `study.py` (the `_SWEEP_CONDUCTOR_REASON` constant and the inserted refusal in `_study_solve_sweep`, before the frequency check and `build`) and `cli.py` (the one edited `_load(...)[:3]` line and the inserted `sweep_hz:` print) only; no `study_yaml.py` edit, no `_run_palace` edit. Write the module whole, then run W0, RG-R, RG1, W9, RG-W1a, RG-W1b, RG-S2d1, RG-S2d2, Wgen, Wval, Wapi, Wdoc in order, foreground, every harness Bash call with `timeout: 660000`, pytest with `-s`, container `timeout -k 30`. `pgrep -c python3` must read 0 before RG-W1a, RG-W1b, RG-S2d1 and RG-S2d2. Every RG command is copied verbatim from the log the brief names. Never import a `test_*` name; never edit an existing test assertion. You are the executor — do not spawn agents, never background a harness window. If the docker socket is denied, stop and report the exact command string. Park if code + W0 are not green by minute 25. On a red, park on `attempt/WF-1-step23-*`, mark item 589 BLOCKED in the same commit, and report."

**Outcome (2026-10-22, 20:00 UTC slot, first item): DONE, first try.** `src`: `study.py` (`_SWEEP_CONDUCTOR_REASON` beside `_validate_sweep_frequencies`; the refusal after the excitation check in `_study_solve_sweep`, before the list check and `build`; one docstring clause) and `cli.py` (the edited `_load(args.file)[:3]` line; the `sweep_hz:` print after `frequency_hz:`; one module-docstring clause). New `tests/package/test_wf1_sweep_conductor_check.py`; registry +3 nodes, item 588's lines kept. W0 (`20261022T200300Z_WF-1-s23-W0.log`, 25 passed, 6 s): (Py) hole (`as_hole` True) and solid copper both raise `TypeError` `'conductor: a sweep with a conductor is not supported (Study.solve_sweep has no conductor route)'`, string-equal, also with `()`, 0 builds, no `sweep_result`, `result` None; (c-reach) the conductor-free F-small Study reaches `build`, 1 call; (c-old, printed per ruling 03:00 (4)) the pinned `6f0c8453` `study.py` on a solid copper Study **reached build (conductor ignored), 1 call**; (Y) field `study.sweep_hz`, reason `==` the constant, `femem check` / `run` rc (2, 2), output directory absent, 0 builds; (C) one `sweep_hz: [10000000.0, 64000000.0, 128000000.0]` line right after `frequency_hz:`, key in the file's `study`, hex-equal to the file and to `FREQUENCIES`; (U) 38 files (wf1 non-malformed + `examples/mri/*.yaml`), 1 sweep file (predicted 1), every other line equal to `acf99b1b`'s (its `check` branch is identical to `6f0c8453`'s, which changed one `_estimate` docstring; pin kept as the brief named it); (c-check) the pinned `check` prints no `sweep_hz` line. Re-gates (commands verbatim from the named logs): RG-R 7 passed (3 s, `-n 1`), RG1 13 passed (3 s), W9 installed 1 passed (15 s), RG-W1a 5 passed (53 s; (R′) `6.601e-11 / 9.557e-11 / 9.449e-11`, cells 117206, (c-swap) `1.2367` / `1.2763` — string-equal to U71's `:1934–1935`), RG-W1b 5 passed (53 s; (I-Py) 0.000e+00 / 9.110e-16 / 7.784e-16 ≤ `L0_BAND`, printed beside), RG-S2d1 24 passed (219 s, `-n 2`), RG-S2d2 17 passed (331 s, `-n 8`); `pgrep -c python3` read 0 before each of the four solve re-gates. Wgen 46 files, Wval nodes 340 (table lines 334 → 340: the 3 new nodes sit in 6 per-name tables; distinct nodes 239 → 242, the brief's +3), Wapi 22 passed (20 s), Wdoc `dead=0` (1 s). `mri:34` guide paragraph corrected (ruling 03:00 (3)), disclosed on the page; the example script asserts only `check`'s rc and was not re-run; docrefs over the examples root `20261022T201755Z_WF-1-s23-docrefs.log`: `dead=0 guide=0 stale=125`, exit 2 from the stale-artifact report only (as `EX-137`'s docrefs logs). U75 DONE.

### Review ruling 2026-10-22 03:00 (items 588 U74, 589 U75)

Binds over the briefs. 588 is `WF-16` step 13 by code but filed here because the weekly's U74 row list names `WF-1`; its §9 line names `WF-16`.
- **588:** (1) `Study.tune` defaults `n_slots` to `geometry.leg_count` and refuses a mismatch — accepted, provided the default path is **bitwise** on every existing consumer (the brief's pin at `9b90d7b3`). (2) The name route and helpers taking `n_slots` go beyond the letter: accepted, disclosed (rule (c)); `ring_mode_rows_v2` untouched. (3) The U2 comparison stays at `REUSE_REPRODUCTION_RTOL`; C10a's 1.6e-15 is printed beside it. (4) "To 0" is exact equality. (5) Tune examples are not re-run (default bitwise). (6) One item, not split; if the re-gates do not fit, the slot parks at the last green window and marks 588 BLOCKED (rule (d)).
- **589:** (1) The exception type is the one the drafter chose unless `Study` already raises a dedicated class for invalid combinations, in which case that class. (2) "U71's records bitwise" means the printed record digits equal as strings; the Python-twin comparison stays at its band. (3) The `mri:34` guide's now-false sentence is corrected in this item, disclosed, docrefs `dead=0`. (4) The old silent path's negative is printed against the pinned pre-change commit, not asserted. (5) The degree-2 `-n 8` sweep re-gate is kept (it is one of U71's records).

### Review ruling 2026-10-22 10:30 (item 589)

Binds over the brief and the 03:00 ruling where they differ.
- (1) **Slot-first.** Two slots declined 589 at take-next, at minutes 23 and 28 (attempts.md 12:30Z, 14:00Z). Its own park rule (code + W0 green by minute 25), with ≈ 16 min of RG windows behind it, cannot be met by a slot that reaches it late. So 589 is taken only as a slot's first item, and a later slot skips it at no attempt cost.
- (2) **Line numbers have moved.** Items 588 (`a1194c06`) and 585 (`6f0c8453`) have since edited `study.py` and `cli.py`, so the brief's line numbers are stale. Locate each insertion by its named anchor (the excitation check in `_study_solve_sweep`; the `frequency_hz:` print in `cli.main`'s `check` branch, now `:182`) rather than by line. 588 landed first, so 589 is the item that rebases `validation_envelope.toml`, keeps both sides and re-runs Wgen / Wval / Wapi / Wdoc (brief, "Order against item 588").

### §9 item 674 — 2026-10-27 (§10 chain U, U99 — daily, 2026-10-27, provisional; `WF-1` step 24: `load_result:` in the YAML, so that `femem run` tunes a saved four-port with no solve; chain step (solve))

*Written by the 2026-10-27 10:30 daily review under rule (6): chain U's next step toward subgoal 9a (a Study API with a thin YAML / CLI of equal power), provisional until the 10-28 weekly ratifies it. U94 (§9 item 661) gave the Python API `Study.load_result(path)` → `LoadedPortResult`, and `tune(k=1)` on it with 0 solves (`study.py:1281`, `:1258`; `tests/package/test_wf16_study_load_result.py`). The YAML has `tune:` (`WF-1` step 11, `study_yaml.py:1178`, `_tune_summary` `:1595`) but no key that loads a result: the fact sheet found 0 hits for `load` / `result` in `study_yaml.py` and `cli.py`. Every earlier YAML step was the CLI twin of a Python route, gated CLI = API. **Depends on nothing** (U94 is on `main`). Independent of item 671 (U90b, a page and a test, with no `src` edit).*

**Route — additive in `src` (rule (c)).**
- `study_yaml.py` takes a new optional top-level key **`load_result: <path>`** (appended to `_TOP_KEYS` at `:157–158`; the path resolves relative to the YAML file, as the file's other paths do).
- **Admitted** only with `geometry.kind: four_port_birdcage`, `excitation.kind: lumped_port_drive` and a `tune:` section: the loader's own type checks (`study.py:1295`, `:1298`) raised at parse time.
- **Refused at parse, by field, before any mesh:** `outputs.sar`, `outputs.b1_plus`, `outputs.probes`, `thermal:`, `conductor:`, `study.sweep_hz` and `--backend palace`. The refusal texts name `load_result`, as `_refuse_loaded` does (`study.py:1270`).
- `run_study_file` (`:1915`) calls `study.load_result(path)` and then `_tune_summary` on it, and writes the same `tune` summary block the solve route writes. `femem check` prints the key.
- Every existing file parses and runs unchanged. Prove it with the `WF-1` step-11 tune test unedited and green (Wreg).

**Anchor (imported, unmoved).**
- **(I=)** The YAML run's `summary["tune"]` equals the same Study's `load_result(path)` + `tune(1)` **by `float.hex`** on `c_mode_farad`, `off_diagonal_residual` and every `z_k` (`tests/package/test_wf1_yaml_tune.py:198–215`'s pattern, copied, not imported).
- **(CM)** That `c_mode_farad` equals `ITEM429_C_MODE_F` at `.9e`, and within `LEG_D0_REPRODUCTION_BAND` (both imported as `test_wf16_study_load_result.py:36–37` does).
- **(S0)** `Study.solve` and `study._solve_ports` are monkeypatched to raise; the YAML route makes 0 calls.
- **(CLI = API)** `femem run <file>` at `-n 2`: rc 0, and its summary JSON diffs ∅ against the API run field by field (`_diff` / `_parse` from `tests/package/test_wf13_cli_estimate.py`, as `test_wf17_estimate_tutorial.py:256` uses them).
- **(Ref)** Each refused field gives rc 2 / `StudyConfigError` naming the field and `load_result`, with 0 mesh builds (`MeshGenerator` entry points and `gmsh.initialize` patched to raise, U96's (NoB) pattern).

**The matrix.** In-test, solve the step-6(b) fixture once (`_study` from `test_wf16_study_tune_capacitors.py`, imported as U94's test does) and `export_touchstone` it to `/workspace/logs/wf1_s24/s6b_pec_d1_n{size}.s4p` (gitignored). The YAML file is written beside it.
- **Synthetic matrices for the coding / refusal tests:** U94's circulant `_write` route (`test_wf16_study_load_result.py:87–89`), with no mesh.

**Negative control (asserted).** A YAML pointing at the same file with P1 ↔ P2 swapped returns a `c_mode_farad` that misses `ITEM429_C_MODE_F` by > `LEG_D0_REPRODUCTION_BAND`. U94's (Sw) measured 0.197 relative.

**Envelope.** `tests/package/validation_envelope.toml` gains the new test nodes under `[keys.load_result]` (or the file's existing YAML-key convention, as written), with the anchor `WF-16` step 14 and this step. Regenerate the API / validation pages before Wapi, as U95 did.

**Windows** (bare, foreground, Bash `timeout: 660000`; complex; `source /usr/local/bin/dolfinx-complex-mode`, `FEM_EM_REQUIRE_COMPLEX=1`):
- **W0** coding + refusals, `-n 1` (smoke).
- **W1** the solve + (I=) / (CM) / (S0) / (CLI = API) / control, `-n 2`. Predicted 40–90 s from U94's W1 of 34 s plus the CLI re-launch.
- **Wreg** the step-11 tune test, `-n 2`.
- **Wgen → Wapi**.
- **Wdoc** `docrefs` 0 dead.
- **W-env.** After `OPS-72` the full `tests/environment` may lead W1 instead (the 10:30 ruling).

**Executor:** implementer. **Tier:** standard, `-n 2`, complex. **Class:** chain step (solve). **Cost:** ≈ 11.5 measured / 40 predicted slot-minutes. **Family:** U99 (`WF-1` step 24), attempt 1 of 4.

**Traps already paid for.**
- `git show` inside the container needs `-c safe.directory=/workspace` (item 662's W1 red).
- Stale generated API pages make Wapi red: regenerate first (item 654's first Wapi).
- `cell_tags.values` and `assemble_scalar` are rank-local.
- A `float.hex` comparison is per call, on the same matrix: never compare across two solves.
- `-k a or b` inside the quoted container command splits into stray argv: select tests by node ID.

**Scope.** The YAML / CLI twin of `load_result` + `tune`. **Not claimed:** a guide page (U90b is the Python page; a YAML page is a later step), Route A on a loaded matrix (returned to the 10-28 weekly), or F-human files (`load_result` refuses a non-four-port geometry, and this item keeps that refusal).

**Negative result.** (I=) red ⇒ the YAML route is not the API route: report the first differing field and park. Wreg red ⇒ the default path changed: revert and park. No `PREDICTED RED:` line.

**Status it can move.** `WF-1` gains step 24 (provisional record; the weekly ratifies). Chain U's 9a "YAML of equal power" covers the saved-result route. The `--estimate` / `check` lines are unchanged.

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 674 — `WF-1` step 24 (U99) — per `docs/planning/chunks/WF-1.md` "§9 item 674". The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. The `src` change is additive; every existing YAML parses unchanged. Spawned executors run every window foreground with Bash `timeout: 660000`."

**Outcome (2026-10-27, 16:30 slot, second item; implementer): DONE, green on the first W1.** `src` (rule (c), additive): `study_yaml.py` — `load_result` appended to `_TOP_KEYS`; `_refuse_with_load_result` (refuses `conductor`, `thermal`, `study.sweep_hz`, `outputs.sar` / `b1_plus` / `probes` by field, and a missing `tune` as field `load_result`, right after `_parse_study`, before geometry); `_parse_load_result` (four-port + `LumpedPortDrive` only; a relative path resolves against the YAML file's directory); a pre-branch refusal for `thermal_ball` / `thermal_box`; `_run_load_result` (no Touchstone written; summary keys `degree`, `family`, `femem`, `frequency_hz`, `load_result`, `port_ids`, `ranks`, `s_matrix`, `sha256`, `source`, `tune`, `elapsed_s`; no `cells`); `_PALACE_OUTPUT_FIELD` gains `load_result`; module docstring bullet. `cli.py`: `check` prints `load_result: '<abs path>'` when the key is present. New `tests/package/test_wf1_yaml_load_result.py`; registry +2 nodes under `[names.run_study_file]`. W0b (`20261027T214337Z_WF-1-s24-W0b.log`, `-n 1`, 5 passed, 2 s): parse spec diff vs `f_small_hole_pec_tune.yaml` `[]`; synthetic (I=) `0x1.7ae3f5a5c44f3p-35` both; 7 refusals + 3 non-four-port files + Palace, 0 builds. W1 (`20261027T214345Z_WF-1-s24-W1.log`, `-n 2`, 6 PASSED per rank, 34 s; solve 29.8 s, 80 798 cells): (RT) 5.869e-14 printed; (I=) `0x1.2373183ff85f7p-34` both, residual and 4 `z_k` equal, differing `[]`; (CM) 6.626796420593e-11 F, `.9e` equal, rel 6.140e-11 ≤ 1e-9; (S0) 0 calls; (CLI = API) rc (0, 0), ∅; control 0.19679. Wreg (`20261027T214427Z_WF-1-s24-Wreg.log`, 12 passed per rank, 66 s): `test_wf1_yaml_tune.py` unedited plus the two no-mesh `test_wf1_study_yaml.py` parse nodes. Wgen / Wval / Wapi (22 passed) / Wdoc (0 dead). **Deviations:** (CLI = API) compares the two `summary.json` files by key equality (`elapsed_s` excepted) rather than through `_diff` / `_parse`, which parse printed dataclass lines, not a JSON summary; the step-11 Wreg was widened by the two parse nodes; the W0 was re-run once (W0b) only so its parse prints were not swallowed by `capsys` — no assertion changed. U99 DONE (provisional; the 10-28 weekly ratifies).

### §9 item 697 — 2026-10-30 (§10 U105, chain U, 9a — daily, 2026-10-30, provisional; `WF-1` step 25: YAML `load_result:` on an `f_clinical_birdcage`; chain step (solve))

*Written by the 2026-10-30 18:00 daily review. **Depends on item 696 landing; if it did not, skip at no cost.** U102 (the same on the F-human) is held for the weekly because `f_human_birdcage` requires `conductor`, which `load_result` refuses (`study_yaml.py:348`, `:1263–1265`). **That design question does not arise here:** `f_clinical_birdcage` takes no key but `kind` (`:308–311`). This is U99's route on the fixture, not a resolution of U102.*

**Route (additive, `src/fem_em_solver/study_yaml.py`).**
- `_parse_load_result` (`:1281`): `kind` in (`four_port_birdcage`, `f_clinical_birdcage`); the message names both.
- `_parse_mapping`'s early refusal (`:1321`) reads the same; `thermal_ball` / `thermal_box` stay refused.
- `run_study_file` keeps the order at `:2061`: the `load_result` branch runs first, so the "xl tier only" refusal (`:1993`, `:2036`) still fires on a file without `load_result:`. `--estimate` is unchanged; if `:2036`'s path sits ahead of `:2061`, move nothing and report it.
- `study_yaml.py`'s module docstring (`:110`) gains the kind.
- `_refuse_with_load_result`: unchanged.

**New module `tests/package/test_wf1_yaml_load_result_f_clinical.py`,** after `tests/package/test_wf1_yaml_load_result.py`. Knob `FEM_EM_WF1_S25=1` for the record node; the file written to `tmp_path` from `logs/geo38/c3-cu-s64-d1.npz` as item 696 does (**absent ⇒ fail**). The YAML sets `study`, `geometry: {kind: f_clinical_birdcage}`, the sections `_REQUIRED_SECTIONS` demands, `excitation.kind: lumped_port_drive`, `tune: {k: 1}`, `load_result: <path>`.

**Anchors.**
- (I=) the run's summary `C_mode` equals the API's `load_result` + `tune(1)` on the same file by `float.hex`;
- (CM-r) `.4f` = "195.0265" (item 696's constant, imported);
- (S0) 0 solves and 0 `build` calls;
- CLI = API: `femem run` with the same file, the `_diff` route (`tests/package/test_wf13_cli_estimate.py`), ∅;
- (K2) `tune: {k: 2}` misses (`LEG_D0_REPRODUCTION_BAND`);
- refusals: a sweep, a field output or a Palace backend beside `load_result:` is refused before any mesh (`_refuse_with_load_result`), and `f_clinical_birdcage` **without** `load_result:` still returns rc 2, "xl tier only" (`test_wf13_estimate_f_clinical.py:152`, unedited, in Wreg).

**Negative control (asserted).** The swapped file misses by > `LEG_D0_REPRODUCTION_BAND`, as item 696.

**Registration.** Envelope, beside step 24's lines (`validation_envelope.toml:550`). A schema addition, not a version bump: `SCHEMA_VERSION` unchanged, as U99.

**Windows.** W0 `WF-1-s25-W0` (`-n 2`, complex, `timeout -k 30 120`); W1 `WF-1-s25-W1` (`-n 2`, complex, `FEM_EM_WF1_S25=1`, `timeout -k 30 180`); Wreg `WF-1-s25-Wreg` (`-n 2`, `timeout -k 30 300`: `test_wf1_yaml_load_result.py`, `test_wf13_estimate_f_clinical.py`, `test_wf13_cli_estimate.py`, unedited); W-env. **Class:** chain step (solve), ≈ 8.3 measured / 30 predicted minutes (item 674's family, no solve run). Tier standard.

**Traps already paid for.** Item 674's: `load_result` rides in `outputs` (`:1419`), and a relative path resolves against the YAML file. `capsys` swallows parse prints (item 674's W0b). Check the literal `run_study_file` order before asserting refusals. `.s64p` is the file extension: `_parse_load_result` does not check it, `load_touchstone` does.

**Scope.** The YAML twin of item 696. **Not claimed:** the F-human YAML route (U102 stays the weekly's), a page (no tutorial: its input is a gitignored `xl` product, so a clean clone cannot run it), a schema-stability statement.

**Negative result.** (I=) red means the YAML route computes through a different path: report and stop. No `PREDICTED RED:` line.

**Status it can move.** §10 U105 → DONE; the `WF-1` row gains step 25.

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 697 — `WF-1` step 25 — per `docs/planning/chunks/WF-1.md` "§9 item 697". Only after item 696 is DONE. Do not resolve U102. Absent npz ⇒ fail, never re-solve. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped; foreground, Bash `timeout: 660000`. Do not spawn agents."

**Outcome (2026-10-31T00:53Z, 19:30 CDT slot): DONE — `WF-1` step 25, with one departure for the weekly.** `femem run` tunes C3's saved 64×64 from YAML with no mesh and no solve.

**The departure: the brief's route could not reach its target, and two clauses beyond its letter landed under standing rule (c).** The brief reads "the `load_result` branch runs first". It does not: `run_study_file`'s `f_clinical_birdcage` refusal (then `:2036`) sits ahead of the branch at `:2061`, and the brief says that in that case, "move nothing and report it". Separately, `_parse_tune` refused `f_clinical_birdcage` outright, which the brief did not name.
- **`run_study_file`:** nothing moved. The refusal's condition gained `and "load_result" not in outputs`, so it fires on exactly the files the brief says it must.
- **`_parse_tune`:** admits `f_clinical_birdcage` with `lumped_port_drive`. So `tune:` on an F-clinical file **without** `load_result:` now parses (`load_study` succeeds) and is refused at `run_study_file` instead of at parse: rc 2, "xl tier only", asserted (W0-r2 `:60`).
- **`_run_study_file_palace`'s refusal is unchanged,** so `--backend palace` beside `load_result:` is refused by "xl tier only", ahead of Palace's own `load_result` refusal: rc 2 and before any mesh, but not by the name U99's route gives.

I read "move nothing" as "do not relocate the refusal", not "stop": the item's target and the refusal it protects both hold this way. If the weekly meant "stop", reverting is the two clauses plus this module. (Status it can move stays this row's step record, provisional under rule (6).)

Anchors at `-n 2` (`20261031T005004Z_WF-1-s25-W1.log:65–70`; 2 passed, Status 0, 4 s):
- **(I=):** summary C `0x1.acde2c7d42a1cp-33` equals the API's; all 17 z_k and the JSON too.
- **(CM-r):** "195.0265" (item 696's constant, imported); 3.531e-14 rel to C10, printed.
- **(S0):** 0 calls with `solve`, `_solve_ports`, `build` and `gmsh.initialize` patched to raise.
- **CLI = API:** rc (0, 0), summary diff [].
- **(K2):** k = 2 reads 1.512684367669e-10 F, a miss of 2.243699e-01.
- **Control (asserted):** the swapped file misses by 2.711891e-02, the same reading as item 696.

W0 was first run green with its readings swallowed by `capsys` (`20261031T004937Z_WF-1-s25-W0.log`, kept, the brief's named trap); the module's prints then went past `capsys`, and W0-r2 (`20261031T004957Z_WF-1-s25-W0-r2.log`, 1 passed, 1 skipped, 3 s) carries them, including (Ref) on a sweep, `outputs.sar`, `outputs.b1_plus`, no `tune`, Palace, and the two no-`load_result:` files. W1 ran on that same module. Wreg (`20261031T005017Z_WF-1-s25-Wreg.log`): the three named modules unedited, 9 passed per rank, 49 s, which re-runs `WF-13` step 12's refusal and U99's gates green after the `src` edit. Wgen / Wval Status 0 (`…005108Z…`, `…005112Z…`; `run_study_file.md` and `validation.md` regenerated); Wapi 22 passed, 24 s (`…005119Z…`); Wdoc 0 dead (`…005145Z…`); W-env 12 passed per rank (`…005150Z…`).

**Also disclosed.** The files live under gitignored `logs/wf1_s25/`, not `tmp_path` (a per-rank path under `-n 2`). CLI = API uses U99's `_summary_diff` on `summary.json`: the `_diff` route the brief cites compares an `estimate` dataclass, which this run does not produce. `femem check` on this file was not exercised. `--estimate` is unchanged. `SCHEMA_VERSION` is unchanged. U102 is untouched.

## §7 row moved verbatim 2026-10-31 (rotate_plan_archive.py chunks)

 `femem run`: the MRI study CLI and YAML config | ✅ (item 256, U5): `femem run` builds Studies equal to the Python-built ones by `float.hex` and reproduces `MAT4_DEG1_RECORD` and `P1` at their bands; steps 2–19 are records at their imported bands — the installed console script, `thermal:`, `conductor:`, the SAR / B₁⁺ / probe outputs, the F-human two-port schema, `tune:`, `capacitor_sheets`, and the thermal-only `thermal_ball` / `thermal_box` files with their sources, duty cycles and probes (U57 … U67) — all ratified by the weekly; step 20 (2026-10-20, provisional under rule (6), the weekly ratifies): point SAR / \|E\| probes on the four-port file (U70); step 21 (2026-10-20, provisional): multi-frequency sweep from one file (U71) — `study.sweep_hz` drives `Study.solve_sweep` through `run_study_file`, one mesh, one `.s4p`, per-frequency S vs `P1` 6.601e-11 / 9.557e-11 / 9.449e-11, θ_file = θ_Python ≤ 9.010e-15; sweep + `conductor` / hole / capacitors / SAR / B₁⁺ / probes / tune / `--backend palace` refused by field; step 22 (2026-10-20, provisional): U70's point probes on the copper F-small file (U72) — `cu` class means 4.093e-11 / 8.673e-11 / 1.376e-10, 80 798 cells, θ_file = θ_Python ≤ 1.824e-13, the 1 mm displaced probe misses, no `src` edit. **Step-23 record (§9 item 589, §10 U75, 2026-10-22):** `Study.solve_sweep` refuses a conductor by name (`conductor: <reason>`, the YAML reason byte for byte, hole and solid copper, 0 builds); `femem check` prints `sweep_hz` on a sweep file, every other file's lines equal `acf99b1b`'s; U71's printed records reproduced (6.601e-11 / 9.557e-11 / 9.449e-11, 117206 cells; c-swap 1.2367 / 1.2763); `src`: `study.py` (insertions), `cli.py` (insertions plus one edited line) (`20261022T200300Z_WF-1-s23-W0.log`). **Step-24 record (§9 item 674, §10 U99, 2026-10-27, provisional under rule (6), the weekly ratifies):** `load_result: <path>` tunes a saved four-port with no mesh and no solve; the summary's `tune` is `float.hex`-equal to `load_result` + `tune(1)`, `C_mode` `.9e` = `ITEM429_C_MODE_F`, `femem run` = API field by field, field outputs / sweep / conductor / thermal / Palace refused by name; the P1 ↔ P2 swapped file misses by 0.19679 (`20261027T214345Z_WF-1-s24-W1.log`). **Step-25 record (§9 item 697, §10 U105, 2026-10-30, provisional under rule (6), the weekly ratifies):** the same on an `f_clinical_birdcage`: C3's saved 64×64 through `femem run`, (I=) `float.hex`-equal, `.4f` = 195.0265 pF, 0 solves, CLI = API, k = 2 and the swapped file miss by 0.2244 / 0.02712, and a file without `load_result:` still returns rc 2 (`20261031T005004Z_WF-1-s25-W1.log`). *History: `docs/planning/chunks/WF-1.md`.* 

### §9 item 714 — 2026-11-01 (§10 U102, chain U, 9a — `WF-1` step 26: YAML `load_result:` on an `f_human_birdcage`; chain step (solve), standard, complex)

*Enumerated by the 10-31 weekly, which ruled the design question: `conductor` is admitted beside `load_result:` on this kind (it documents what the saved column was solved with; nothing is solved), and every solve-time key is refused as U99 refuses it. Depends on nothing. Next: U108. Family U102, attempt 1 of 4.*

**Route (additive, `src/fem_em_solver/study_yaml.py`).**
- `_refuse_with_load_result` (`:1263`) refuses `conductor` by name. Its `for name in ("conductor", "thermal")` needs the kind: refuse `conductor` unless `geometry.kind == "f_human_birdcage"`, read defensively as `_parse_mapping` already reads `geo` (`:1325`). `thermal`, `study.sweep_hz`, `outputs.sar`, `outputs.b1_plus` and `outputs.probes` stay refused on every kind; the Palace refusal stays.
- `_parse_load_result` (`:1292`) and the early refusal's message (`:1330`): admit `f_human_birdcage`, naming all three kinds. Its `LumpedPortDrive` check stands: the file's excitation must be `lumped_port_drive`, which the F-human admits since `WF-1` step 12 (`:428`). The `two_port_quadrature` default is refused by that check, by field.
- `_run_load_result` (`:1961`): pass `n_slots=16` to `study.tune` on this kind. U45's route (`test_wf1_yaml_tune_f_human.py`) shows how the YAML's `tune:` carries it; follow that, no new key.
- Module docstring (`:110`) gains the kind.

**New module `tests/package/test_wf1_yaml_load_result_f_human.py`,** after `test_wf1_yaml_load_result_f_clinical.py`. Knob `FEM_EM_WF1_S26=1` for the record node. The file is **U103's saved `logs/wf16/s15_f_human_cu_d1_n8.s32p`** (present today; gitignored; **absent ⇒ fail, never re-solve**), copied to `tmp_path`.

**Anchors.**
- (I=) the run's `tune` block equals `Study.load_result` + `tune(1, n_slots=16)` by `float.hex`;
- (CM) `C_mode` against `WF12_STEP3_C_MODE_CU_F` at `REDUCTION_BAND`, as `test_wf16_load_result_f_human.py` reads it (item 685: 3.421e-13 rel);
- (S0) 0 solves and 0 `build` calls, with `solve`, `_solve_ports`, `build` and `gmsh.initialize` patched to raise;
- CLI = API: `femem run` with the same file, the `_diff` route, ∅;
- (Ref) `study.sweep_hz`, `thermal`, `outputs.sar`, `outputs.b1_plus`, `outputs.probes` and `--backend palace` each refused by field, before any mesh, beside a present `conductor`;
- (Ref-c) `conductor` beside `load_result:` on a `four_port_birdcage` is **still** refused: U99's gate, in Wreg, unedited.

**Negative control (asserted).** The swapped file `logs/wf16/s15_swap12_n8.s32p` misses by > `LEG_D0_REPRODUCTION_BAND`.

**Registration.** Envelope line beside step 25's. `SCHEMA_VERSION` unchanged, as U99.

**Windows.** W0 `WF-1-s26-W0` (`-n 2`, complex, `timeout -k 30 120`); W1 `WF-1-s26-W1` (`-n 2`, complex, `FEM_EM_WF1_S26=1`, `timeout -k 30 180`); Wreg `WF-1-s26-Wreg` (`-n 2`, `timeout -k 30 300`: `test_wf1_yaml_load_result.py`, `test_wf1_yaml_load_result_f_clinical.py`, `test_wf1_yaml_tune_f_human.py`, unedited); Wgen, Wval, Wapi, Wdoc; W-env. **Price:** item 697's W1 read 4 s and its Wreg 49 s; a 32-port file is smaller than its 64-port one. **Class:** chain step (solve), ≈ 7.6 measured / 30 predicted minutes. **Tier:** standard.

**Traps already paid for.** Item 697: the literal order of `run_study_file` differs from a brief's reading, so check it before asserting a refusal. `f_human_birdcage` requires `conductor` (`:349`) and `as_hole: true`, and `_REQUIRED_SECTIONS` still applies. Item 674: `load_result` rides in `outputs`, and a relative path resolves against the YAML file. `.s32p` is `load_touchstone`'s check, not `_parse_load_result`'s. `capsys` swallows parse prints: print past it, as item 697's W0-r2 did. `_parse_conductor`'s `sigma` is parsed and unused here: assert it reaches no solver through (S0), and do not add a check that the file matches it.

**Scope.** The YAML twin of item 685. **Not claimed:** that the saved file was solved with the stated conductor, a sweep, the two-port drive, a schema-stability statement.

**Negative result.** (I=) red ⇒ the YAML computes through a different path: report and stop. A departure from this brief's route is reported as item 697 reported its own. No `PREDICTED RED:` line.

**Status it can move.** §10 U102 → DONE; the `WF-1` row gains step 26.

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 714 — `WF-1` U102 — per `docs/planning/chunks/WF-1.md` "§9 item 714". Absent `.s32p` ⇒ fail, never re-solve. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped; foreground, Bash `timeout: 660000`. Do not spawn agents."

**Outcome (2026-11-01T12:32Z, 06:00 CST slot): DONE — `WF-1` step 26, a record, with one landed test line edited against the brief (below).** `src/fem_em_solver/study_yaml.py`: `_refuse_with_load_result`, `_parse_load_result`, the early refusal's message and the module docstring. New module `tests/package/test_wf1_yaml_load_result_f_human.py`. Envelope lines added; `SCHEMA_VERSION` unchanged.

Asserted (W1-r2 `20261101T122847Z_WF-1-s26-W1-r2.log`, `-n 2`, 2 passed per rank, 4 `PASSED` tokens, Status 0, 4 s):
- **(I=):** the summary and its JSON are `float.hex`-equal to `Study.load_result` + `tune(1, n_slots=16)`; differing [] on both; 32 ports, 9 z_k.
- **(CM):** 3.421e-13 against `REDUCTION_BAND` 1e-3, item 685's reading exactly.
- **(S0):** 0 calls.
- **CLI = API:** rc (0, 0), summary diff [].
- **(Ref):** `study.sweep_hz`, `thermal`, `outputs.sar`, `outputs.b1_plus` and `outputs.probes` each refused by their own field, beside a present `conductor`, rc (2, 2), no output directory. `--backend palace` rc (2, 2).
- **(Ref-c):** U99's gate on a `four_port_birdcage` still green in Wreg-r2.
- **Control:** the swapped file misses by 4.254907e-01 against `LEG_D0_REPRODUCTION_BAND` 1e-9.

Other windows: W0 `20261101T122435Z_WF-1-s26-W0.log` (1 passed, 1 skipped). Wreg-r2 `20261101T122855Z_WF-1-s26-Wreg-r2.log` (4 passed, 5 skipped per rank, 33 s). Wgen, Wval (`run_study_file.md`, `validation.md` regenerated), Wapi (22 passed, `…122946Z…`), Wdoc (dead=0, `…123011Z…`), W-env (12 passed, `…123014Z…`), all Status 0.

**The departure, for the 11-02 weekly.** The brief lists `test_wf1_yaml_load_result.py` as unedited in Wreg, but its line 216 asserted that the F-human file beside `load_result:` is refused at parse, which is precisely what the weekly's ruling removes. Both cannot hold. Wreg ran first unedited and read Status 1 at that loop and nowhere else (`20261101T122750Z_WF-1-s26-Wreg.log`: `DID NOT RAISE`, 1 failed, 3 passed, 5 skipped). The F-human then left that list, one line with a comment citing this item. Item 709 set the precedent for a one-line edit disclosed against "unedited". This is a behaviour change the weekly ruled, not a loosened gate, and the other two kinds on that list keep their refusal. Because the old loop raised before the Palace refusal below it, Wreg-r2 is also the first run of those lines since the `src` edit.

**One re-run, under the 10-30 ruling.** W1 (`20261101T122449Z_WF-1-s26-W1.log`, Status 124, 181 s) failed both nodes in 1.6 s and then hung until the timeout, a rank-0 raise leaving the other rank at a `Barrier`. It is a crash in non-gate code before any gate value, not a size kill: `pgrep -c python3` read 0 after it. The two causes were the module's own: `_fresh_dir` does not recreate the directory `shutil.copy` wrote into, and the assertion below. W1-r2 ran at the same `-t`; the first log is kept.

**Also disclosed.**
- **The brief's refusal route for the quadrature default is not the one that fires.** A `two_port_quadrature` file with `tune:` is refused by the `tune` parser, field `tune` (item 697's "literal order" trap), ahead of `_parse_load_result`'s `LumpedPortDrive` check, so that check is not reached by this file. The module asserts the field that fires.
- **`_run_load_result` was not edited.** `tune: {k: 1}` reaches `Study.tune` with `n_slots` unset, which is `leg_count` = 16, and (I=) against `n_slots=16` is equal by `float.hex`. That is U45's route, and passing it explicitly would have added a kind branch for nothing.
- **`--backend palace` is refused by its kind check** (`palace supports four_port_birdcage only`), field `--backend`, so the message does not name `load_result`.
- W0's synthetic node carries no closed-form check: `_synthetic()` returns no capacitance, and the record node holds (CM).
- The files live under gitignored `logs/wf1_s26/`, not `tmp_path` (per-rank under `-n 2`), as item 697.
- W0 ran before the last edits; W1-r2 runs its node on the module as committed.

**Status moved.** §10 U102 DONE; the `WF-1` row gains step 26 by this outcome.

### §9 item 730 — 2026-11-01 (chain U, 9a, U111 — daily, provisional, rule (6))

**Why.** Item 722 gave `Study.tune` an `at_hz`, and YAML `load_result:` on a multi-frequency file now raises the `at_hz` refusal, with no way to supply one: `_parse_tune` admits `k` only. The API and the file route differ, and "YAML is a thin layer" (9a) is the subgoal. **Provisional.**

**Anchor.** On U71's file, `logs/wf16/wf16_sweep_f_small_p1_n2.s4p` (absent ⇒ fail, never re-solve), for each f in the file: (I=) the summary and its JSON are `float.hex`-equal to `Study.load_result` + `tune(k=1, at_hz=f)`; (C=) `C_mode` at `REDUCTION_BAND` against item 722's record column (`2.472675544637e-09` / `6.166035421775e-11` / `1.513380716932e-11` F); (S0) 0 solves; CLI = API, rc (0, 0), summary diff [].

**Negative control.** The 128 MHz summary read against the 64 MHz `C_mode` misses by ≈ 0.75 (the ratio from the file's own columns; compute the ceiling from them, as item 722 did, and assert against `LEG_D0_REPRODUCTION_BAND`). (Ref): `tune.at_hz` outside the span, NaN, a non-number, `at_hz` on a one-frequency file, and `at_hz` with no `load_result:` each refuse by field `tune.at_hz`, rc (2, 2), no output directory.

**Execution.** `study_yaml.py`, additive under rule (c): `_parse_tune` admits `at_hz` (`_unknown(sec, ("k", "at_hz"), p, source)`), `_tune_summary` forwards it, the summary records the selected frequency. New module `tests/package/test_wf1_yaml_load_result_sweep.py`, knob-free. One envelope line; the F-clinical multi-frequency refusal, which item 722 disclosed as unexercised by name, gets one assertion here. `SCHEMA_VERSION` unchanged (additive key).

**Tier, ranks, wall clock.** Standard, `-n 2`, complex, `timeout -k 30 170`. W1 ≈ 5 s (item 722: 3 s). Wreg: `tests/package/test_wf1_yaml_load_result.py tests/package/test_wf1_yaml_load_result_f_human.py tests/package/test_wf1_yaml_tune.py tests/package/test_wf16_load_result_sweep.py`, ≈ 60 s. Wgen, Wval, Wapi, Wdoc, W-env. **Class:** chain step (solve), by precedent (items 714, 722 were classed so). **Cost:** ≈ 5.5 measured / 20 predicted min. **Depends on:** nothing.

**Traps already paid for.** `summary["frequency_hz"]` reads `study.frequency_hz`, which is not the selected frequency on a multi-frequency load: do not assert it equal to f unless you make it so. The record is `C_mode` and the selected frequency, printed separately. A rank-0 raise leaves the other rank at a `Barrier` (item 714): raise on every rank. `_fresh_dir` does not recreate a directory `shutil.copy` writes into (item 714). `tmp_path` differs per rank: files live under gitignored `logs/wf1_s27/`. `--backend palace` refuses by its own kind check, and the message does not name `load_result` (item 714).

**Scope.** The YAML key. **Not claimed:** interpolation between frequencies (`tune` selects the nearest), a multi-frequency F-human path, a resonance.

**Negative result.** (I=) red is a defect in the thin layer: stop, BLOCKED. No `PREDICTED RED:` line.

**Status it can move.** `WF-1` gains step 27; §10 U111 closed; item 737 unblocks. `WF-1` stays as it is.

**Outcome 2026-11-01 (15:00 slot, item 730): DONE — `WF-1` step 27, a record.** Run inline, no executor; every gate green on the first W1. `src` (rule (c), additive): `study_yaml.py`, `_parse_tune` admits `at_hz` and takes `has_load_result`; `_run_load_result` refuses a one-frequency file and an out-of-span `at_hz` by field, on every rank, before any write; `_tune_summary` needed no edit, because it already forwards the spec and writes the selected frequency. New module `tests/package/test_wf1_yaml_load_result_sweep.py`, knob-free, one envelope line. W1-r2 (`20261101T210829Z_WF-1-s27-W1-r2.log`, `-n 2`): 3 passed per rank, 4 s, Status 0.

| f (Hz) | selected | `C_mode` (F) | rel to item 722's record | (I=) differing | CLI rc, diff |
|---|---|---|---|---|---|
| 1.0e7 | 1.0e7 | 2.472675544637e-09 | 1.935e-13 | [] | (0, 0), [] |
| 6.4e7 | 6.4e7 | 6.166035421775e-11 | 2.265e-14 | [] | (0, 0), [] |
| 1.28e8 | 1.28e8 | 1.513380716932e-11 | 1.912e-13 | [] | (0, 0), [] |

(`20261101T210812Z_WF-1-s27-W1.log:53–61`, the same module's gates.) (S0) 0 calls (`:52`). **Control:** 128 MHz against the 64 MHz `C_mode` misses by 0.754562, equal to the ceiling from the file's own columns, against `LEG_D0_REPRODUCTION_BAND` 1e-09 (`:62`). (Ref), all field `tune.at_hz`, rc (2, 2), no output directory: below and above the span, NaN, a non-number, a one-frequency file, no `load_result:` (W1-r2). F-clinical: a multi-frequency `.s64p` raises "exactly one frequency", nothing stored. Wreg (`20261101T210833Z_WF-1-s27-Wreg.log`, `-n 2`): 11 passed, 3 skipped per rank, 93 s, Status 0. Wgen, Wval Status 0; Wapi 22 passed, 25 s, whose one `dead=1` line is that test's own fixture; Wdoc Status 0; W-env 12 passed, 48 s.

**Four disclosures, for the 11-02 weekly.** (1) **Two W1 logs.** W1 passed with five of six (Ref) readings swallowed by `capsys`, item 697's named trap; those prints moved past the capture and W1-r2 re-ran the module, so rule (i) holds on r2. No assertion changed, and both logs are kept. (2) **`test_wf16_load_result_sweep.py` has one comment edited,** the now-false "YAML adds no route", ahead of Wreg; its docstring's "Not claimed: a YAML route" is left as that step's record. YAML with no `at_hz` still raises `tune`'s own `ValueError`, unchanged. (3) **Wreg is one module wider than the brief,** adding `test_wf1_yaml_load_result_f_clinical.py`, because `_parse_tune` serves that route. (4) **`at_hz` alone admits** and takes `Study.tune`'s default k = 1, as `tune: {}` does. `summary["frequency_hz"]` stays the `Study`'s; the selected frequency is `summary["tune"]["frequency_hz"]`. **Not claimed:** interpolation between frequencies, a multi-frequency F-human path, a resonance. **Provisional:** the 11-02 weekly ratifies or reverts.

**Status moved.** `WF-1` gains step 27 as a record; §10 U111 closes for the weekly; item 737 unblocks.

### Review ruling 2026-11-05 03:00 — U125 is step 28; the walled box's public keys

The 11-04 weekly wrote U125 as "`WF-1` step 10". That number is stale: step 10 is used, and the latest step is 27 (item 730). **U125 is step 28.**

Reproducing P63's R1 "on the same mesh" needs two public-surface choices. Both are written here as execution detail, with the target kept, and **named in the dashboard for the 11-07 weekly to ratify or revert** (as `field_strength_t` was).

**1. `thermal.initial` placement and the wall keys.**
- `thermal.initial` sits under `thermal:` as a sibling of `transient`, so `_TOP_KEYS` is unchanged.
- `geometry.wall` uses the ball's `_WALL_KEYS` names (`thickness_m`, `k_w_per_m_k`, `rho_c_j_per_m3_k`), **plus `divisions`**.
- `divisions` is the number of cells across the wall. P63 uses N_w = 2 (`test_mat21_wall_fem.py:74–75`), and it cannot be derived from the other keys.

**2. A walled box is P48b's 1D-reduced form.**
- The wall sits on ±z only.
- Robin cooling applies on z = ±(C + d) only.
- The four side faces are Neumann (`test_mat21_wall_fem.py:12–15`; P63's `_robin_facets(..., (2,))`).

**"The test's own record".** P63's W1b prints only R1 − R0 at P5, not θ_R1 itself. So the new module computes P63's R1 directly in the same window and compares against that, and it checks the drift against W1b's print.

### §9 item 830 — 2026-11-05 (§10 U125, chain U, 9a — `WF-1` step 28, the thermal-box YAML gains an initial profile and a wall; chain step (solve))

**Executes** §10 U125 under the ruling above.
- `thermal.initial: {gradient_k_per_m: g}` sets θ₀ = g·z (P47).
- `geometry.wall` adds P48b's meshed acrylic.
- A YAML file reproduces P63's R1: θ at P5 equals the test's own record at `IDENTITY_BAND` on the same mesh and Δt, (E) holds, and a wall without a thickness is refused.
- Family U125, attempt 1 of 4. Next: item 838.

**Code** (`src`, rule (c), additive; the defaults leave every existing path byte-identical).

*`ThermalBoxGeometry` (`study.py:1737–1771`).*
- The `wall` property (`:1756–1759`) becomes a field, `wall: Optional[ContainerWall] = None`.
- New field `wall_divisions: Optional[int] = None`, required to be ≥ 1 exactly when `wall` is set.
- **Walled `build`:**
  - `create_box` over [−A, A] × [−B, B] × [−(C + d), C + d], with z cells = `divisions[2] + 2·N_w`.
  - `_remap`, copied from `test_mat21_wall_fem.py:78–88` with its provenance, applied to z only.
  - Cell tags: 1 for gel (by the midpoint test), 2 for the wall.
  - Facet tag 1 on |z| = C + d only (`_robin_facets`' marker, `:128–139`, copied).
- **Unwalled `build`:** unchanged, so `test_wf15_study_thermal_box.py:109`'s `geo.wall is None` holds.

*`TransientSpec` (`:1908`).*
- Gains `initial_gradient_k_per_m: Optional[float] = None`.
- `_study_transient_temperature` (`:1981`) passes `t_init_expr=lambda x: t_init + g * x[2]` when the gradient is set.
- A gradient is refused on `ThermalBallGeometry` and in `_tuned_quadrature_temperature` (`:3220`). It is never silently ignored (item 709's lesson).
- The walled box reuses `_dg0_two_region` (`:1995–2004`).

*Refusals.*
- `_study_box_temperature` (`:1774`) refuses a wall on the steady path and on `TransverseQuadraticHeatSource`.
- `Study.probe` refuses a walled box.

*`study_yaml.py`, `_parse_thermal_box_mapping` (`:968–1098`).*
- `geometry` admits `wall`.
- **Order matters.** With no `thermal.transient`, refuse at the field `geometry.wall` *before* parsing the wall's keys. This keeps the pinned row at `test_wf1_yaml_thermal_box.py:250–251` green, unedited.
- Then `_unknown` / `_missing` run on the wall keys. A missing thickness is reported at `geometry.wall.thickness_m`.
- `thermal.initial` needs `transient`.
- A wall beside `outputs.probes`, `thermal.duty` or `excitation.source` is refused at that field.

*`_run_thermal_box` (`:1844–1906`).*
- Adds the summary keys `wall` / `initial_gradient_k_per_m` **only when set**.

*New files.*
- `tests/package/wf1/thermal_box_wall_gradient_m4.yaml` (R1), with these literals:
  - half widths [0.325, 0.21, 0.045]; divisions [52, 32, 8];
  - wall {0.006, 0.19, 1749300.0, divisions 2};
  - q 2000.0, k 0.6, h 10.0, t_inf 20.0;
  - transient: dt_s = `repr(T_STAR / N_MAIN)`, n_steps 40, ρc as `_rho_c_field` reads it, θ 0.5, startup 2;
  - initial: gradient 10.0.
- `thermal_box_wall_m4.yaml` (R0, the control): the same file without `initial`.
- `malformed_thermal_box_wall_without_thickness.yaml`.
- The module `tests/package/test_wf1_yaml_thermal_box_wall.py`, with knob `FEM_EM_WF1_S28`.
- Envelope lines under `[names.ThermalBoxGeometry]` and `[names.TransientSpec]`.

**Windows** (real build; `pgrep -c python3` reads 0 before each solve window).
- **W0** (`-n 1`, `timeout -k 30 120`, no solve): the parse and the refusals. Each YAML literal equals its imported constant by `float.hex`, including dt == `T_STAR / N_MAIN` exactly. The pinned row's field.
- **W1:**
  `scripts/testing/run_and_log.sh WF-1-s28-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace FEM_EM_WF1_S28=1 timeout -k 30 300 mpiexec -n 2 python3 -m pytest tests/package/test_wf1_yaml_thermal_box_wall.py -v -s --tb=short'"`
  Predicted 40–90 s. P63's four runs took 48 s, R1 alone 2.7 s (`20261104T060527Z_MAT-21-P63-W1b.log:60`, `:94`).
- **Wreg-r** (real, `-n 2`, `timeout -k 30 600`, knobs unset):
  - the parse and refusal nodes of `test_wf1_yaml_thermal_box*.py`, `test_wf1_yaml_thermal_ball.py`, `test_wf15_study_thermal_box.py` and `test_wf17_thermal_box_tutorial.py`;
  - plus `test_wf15_study_transient.py::test_wf15_wall_reproduces_mat18_2b` (159 s).
- **Wreg-c** (complex, `FEM_EM_WF15_S5=1`, `timeout -k 30 600`): `test_wf15_study_thermal_box.py::test_wf15_box_duty_reproduces_mat21_p26` (191 s).
- **Then:** Wgen, Wval, Wapi, Wdoc, in that order.

**Anchors** (bands imported, unmoved).
- **(I):** θ at P5 from the YAML R1 file against P63's direct R1 in the same window: |Δ| / |θ| ≤ `IDENTITY_BAND`.
  - The direct R1 imports `_walled_box(RUNG, (False, False, True))`, `_materials`, `_rho_c_field`, `_robin_facets` and `specs["R1"]` by name (`test_mat21_wall_pre_equilibration_fem.py:202–215`).
  - Evaluation uses the same `READ_PTS` batch, call for call.
  - The Python `Study` equals the YAML at the same band.
- **(E):** `max_balance_residual` ≤ `E_BAND`, and `source_power` equals the direct route's.
- **(Cells):** 119 808 on both builds.
- **(Ref):** the malformed file is refused at `geometry.wall.thickness_m`, and `femem check` returns rc [2, 2].

**Control (asserted).**
- **(c-init):** the R0 file misses R1 at P5 by more than `IDENTITY_BAND`.
- **(Drift):** that miss equals W1b's printed drift, +1.828983e-01 K (`W1b:76`), within 1e-7 (the print's last digit).

**Tier.** Heavy by label, standard by the measured W1. `-n 2`, real. **Class:** chain step (solve). **Cost:** ≈ 11.2 measured / 50 predicted min.

**Traps already paid for.**
- The pinned `geometry.wall` row and its field order (above).
- The box puts Robin on all six faces today (`_exterior_robin_tags`, `:1771`). The walled box must not inherit that.
- `src` imports nothing from `tests/`. `_remap` and the facet marker are copied with file:line.
- A YAML dt literal can miss `T_STAR / N_MAIN` by one ULP. Write the `repr`, and assert it.
- `.eval` bytes depend on the batch (item 537).
- Wgen runs before Wapi (item 653).

**Scope.** The 1D-reduced walled box: one rung, one Δt, from the YAML and from the API. **Not claimed:** the six-wall box, a wall under duty, probes beside a wall, a gradient on a ball, or a cited PMMA ρc (W3).

**Negative result.**
- **(I) red with (E) green:** print the per-ingredient comparison (cells, Σk, Σρc, facet count, θ₀ at P5), park on `attempt/WF-1-s28-*`, and mark BLOCKED.
- **Wreg-c red:** the edit reached the unwalled path. Revert it.
- No `PREDICTED RED:` line.

**Status it can move.**
- `WF-1` gains step 28.
- §10 U125 is DONE.
- Item 838 is unblocked.
- The public keys go to the 11-07 weekly.

**Depends on:** nothing.

**Outcome 2026-11-05 (implementer, slot 21:00 UTC): DONE, green.** W0 (`-n 1`) `20261105T210548Z_WF-1-s28-W0.log` 1 passed / 1 skipped, 5 s: 16 literals equal their imported constants by `float.hex` (dt `0x1.1a26666666666p+6` = `T_STAR / N_MAIN`), spec diffs `[]`, `femem check` rc (0, 0), the pinned row refused at `geometry.wall`, 16 of 16 variant rows and 10 of 10 API refusals, build calls 0. W1 `20261105T210601Z_WF-1-s28-W1.log` (16 s) read every brief anchor green but failed on a test-only gate, summary = `summary.json` on rank 1, because `elapsed_s` is each rank's own (item 390b's record). That gate now leaves `elapsed_s` out and prints the differing keys. No brief anchor or band changed. W1-r2 `20261105T210638Z_WF-1-s28-W1-r2.log` 2 passed, 16 s (`pgrep` 0 before each solve window): direct R1 119 808 cells (52, 32, 12), 6 656 Robin facets on z = ±(C + d); (I) θ(P5) +1.224895450933e+00 K, file vs direct rel 0.000e+00 and Python `Study` vs direct 0.000e+00 (`IDENTITY_BAND` 1e-12; the whole `READ_PTS` batch max rel 7.993e-15, printed; W1 r1 read 2.900e-15); (E) direct 2.306e-13, R1 file 2.283e-13, R0 file 2.364e-13, Python 2.312e-13 (`E_BAND` 1e-9); source_power 4.914000000001e+01 W, rel 0 on both files; (Cells) 119 808 on all builds; (c-init) R0 file vs R1 at P5 rel 1.493e-01 (asserted > 1e-12); (Drift) R1 − R0 at P5 +1.828982665e-01 K vs W1b +1.828983e-01, |Δ| 3.345e-08 ≤ 1e-7; summary keys `wall` / `initial_gradient_k_per_m` present only when set. Wreg-r `20261105T210709Z_WF-1-s28-Wreg-r.log` 12 passed, 71 s (the parse / refusal nodes plus `test_wf15_wall_reproduces_mat18_2b`). Wreg-c `20261105T210828Z_WF-1-s28-Wreg-c.log` 1 passed, 204 s (P26 box duty, unwalled path unchanged). Pages: Wgen 47 files (`20261105T211203Z`), Wval 382 nodes (`20261105T211205Z`), Wapi 22 passed, 28 s (`20261105T211213Z`; its `dead=1` line is a planted control), Wdoc 0 dead (`20261105T211241Z`). **Disclosures for the 11-07 weekly.** (1) The public keys `thermal.initial.gradient_k_per_m` and `geometry.wall.{thickness_m, k_w_per_m_k, rho_c_j_per_m3_k, divisions}` are named for ratification. (2) The API refuses a wall with a duty too, and the YAML refuses the same at `thermal.duty`. A gradient is refused with a transverse-quadratic source, at `thermal.initial` and in the API. (3) The `_tuned_quadrature_temperature` gradient refusal is in `src` but has no node that exercises it, because reaching it needs a tuned Study.

### §9 item 873 — 2026-11-08 (§10 U134, chain U — 9a, WF-1 step 29, the YAML twin of U131; chain step (solve), standard, complex)

*Written by the 2026-11-08 03:00 daily review from the 11-07 weekly's enumerated §10 U134 (rule (4); not provisional). The target, anchor and control are the weekly's; the execution detail is the review's. U131 (§9 item 862, `WF-15.md` "§9 item 862") put `MAT-12`'s transient on the EM sphere in the `Study`, Python only: the sphere's `thermal:` section still refuses `transient` as an unknown key (`study_yaml.py:701`). Family U134, attempt 1 of 4. Independent of U135 / U139. Next: item 887 (U136).*

**Executes** §10 U134. The EM-sphere study file gains a `thermal: {transient: {t_end_s, n_steps, …}}` section. `femem run` reproduces U131's (A_t) and (S2) records at the imported bands, on the same mesh and Δt, with CLI = API by `repr`. Duty, gradient and transient probes are refused from YAML with U131's messages, and a four-port file that names a transient is refused. Envelope entry.

**Code** (`src`: `study_yaml.py` only, rule (c), additive; every existing path byte-identical; `study.py` untouched).
- **`_parse_thermal`** (`study_yaml.py:687–717`):
  - A `sphere_in_box` section carrying `transient`, `duty` or `initial` goes to a new `_parse_thermal_sphere_transient(sec, source)`. A section without them takes the landed lines.
  - The ball, box and F-human parsers call `_parse_thermal(base, "sphere_in_box", …)` with those keys stripped (`:896–897`, `:1083–1084`, `:758–759`), so they never reach the new branch.
  - The four-port check (`:695–698`) stays first, so a four-port file that names `thermal.transient` is refused at `thermal` with the landed reason.
- **`_parse_thermal_sphere_transient`**, in this order:
  - **`duty` present:** `StudyConfigError` at `thermal.duty`. The reason is U131's text, copied literally: "duty is not defined on the EM sphere (WF-15 step 15 gates transient only; MAT-17's duty gate is the thermal-only ball's)" (`study.py:1653–1654` at `9e077d56`).
  - **`initial` present:** refused at `thermal.initial` with U131's text, "initial_gradient_k_per_m is defined on a ThermalBoxGeometry only" (`study.py:1659`).
  - **The `transient` block:** copied from `_parse_thermal_box_mapping`'s (`:1089–1128`, gradient lines dropped; the copy cites the span), with the weekly's keys.
    - `t_end_s`, `n_steps` and `rho_c_j_per_m3_k` are required; `theta` and `startup_half_steps` are optional.
    - `dt_s` is not a key here: it is refused as unknown at `thermal.transient`.
    - `t_end_s` must be finite and > 0 (checked at `thermal.transient.t_end_s`). Build `TransientSpec(dt_s=t_end_s / n_steps, …)`.
  - **Cooling:** `_parse_thermal(<sec without transient>, "sphere_in_box", source)`. Return `_BallThermal(cooling, transient, None)` (`:823–829`; no new carrier class).
- **`_parse_mapping`** (`:1538–1544`): after `outputs` parses, an `outputs.probes` beside a sphere `_BallThermal` is refused at `outputs.probes` with U131's text, "transient probes are not defined on the EM sphere (WF-15 step 15 gates temperature only)" (`study.py:2236–2237`).
- **`run_study_file`**'s sphere branch (`:2352–2369`):
  - If `isinstance(cooling, _BallThermal)`: `thermal = study.temperature(cooling.cooling, transient=cooling.transient, comm=comm)`, which reuses `study.result`.
  - Its `summary["temperature"]` copies `_run_thermal_box`'s transient lines (`:1923–1941`): `transient {t_end_s, dt_s, n_steps, theta, startup_half_steps, rho_c_j_per_m3_k}`, `source_power_w`, `source_multipliers`, `max_balance_residual`, `theta_max_k`.
  - Otherwise the landed lines run, unchanged.
- **Module docstring:** the `thermal` bullet (`:88–92`) names the sphere's `transient` (`WF-1` step 29) and the duty / initial / probes refusals.
- **New files** in `tests/package/wf1/`:
  - `sphere_mat12_transient_a.yaml`: `sphere_mat10_thermal.yaml` with `t_inf_c: 0.0` and `transient: {t_end_s: <repr(T_STAR)>, n_steps: 40, rho_c_j_per_m3_k: <repr(RHO_C)>, theta: 0.5, startup_half_steps: 2}`.
  - `sphere_mat12_transient_s.yaml`: the same with `t_end_s: <repr(DT_1)>`, `n_steps: 1`, `theta: 1.0`, `startup_half_steps: 0`.
  - `malformed_sphere_transient_duty.yaml`: the A file plus a `thermal.duty`.
  - `malformed_transient_on_port.yaml`: `malformed_thermal_on_port.yaml` plus a `thermal.transient`.
  - The other refusals are in-memory variants (`_with` / `_without`, `test_wf1_yaml_thermal_box_wall.py:118–133`).
- **New module `tests/package/test_wf1_yaml_sphere_transient.py`**, knob `FEM_EM_WF1_S29`. The solve node skips without the knob; the parse node runs without it.
  - **Imports, by name, never restated:**
    - `_study` from `tests.package.test_wf15_study_lossy_sphere`;
    - `MAT12_A_T_PRINT`, `MAT12_S2_PRINT` and `_owned_rel` from `tests.package.test_wf15_study_sphere_transient`, plus U131's bands, constants and helpers exactly as that module imports them (`:61–92`);
    - `LOGS`, `WF1`, `_spec_diff`, `_mapping`, `_fresh_dir`, `_rc_all` from `tests.package.test_wf1_study_yaml`;
    - `_with`, `_without`, `_field` from `tests.package.test_wf1_yaml_thermal_box_wall`.
  - **Docstring:** backticks `WF-1`, `WF-15`, `MAT-12` only, and quotes U131's records with log:line.
- **Envelope** (`tests/package/validation_envelope.toml`), as U131's lines (`:365–368`), each with a `# WF-1 step 29 (§9 item 873)` comment:
  - the parse node under `[names.load_study]` (`:433`);
  - the solve node under `[names.run_study_file]` (`:477`) and `[names.TransientSpec]` (`:350`).
  - Wval regenerates `docs/guide/validation.md`.

**Gates (asserted, bands imported, none set).**
- **W0, no mesh** (`test_wf1_s29_parse_and_refusals`):
  - **(K)** On both files, `_spec_diff(load_study(f), _study(degree=1)) == []`.
    - The parsed `RobinCooling`'s `repr` equals `repr(RobinCooling(K, H, 0.0, source_scale=G))`.
    - The parsed `TransientSpec` `repr`s equal U131's A, `TransientSpec(T_STAR / N_MAIN, N_MAIN, RHO_C)`, and S, `TransientSpec(DT_1, 1, RHO_C, theta=1.0, startup_half_steps=0)` (`test_wf15_study_sphere_transient.py:159–161`). `repr` carries every float exactly, so this is the spec half of CLI = API.
    - The parsed dt of each file equals U131's by `float.hex`.
    - `femem check` reads rc (0, 0) on both files.
  - **(Ref)** Each refusal raises `StudyConfigError` at its field, with **the reason equal to the API's `str(exc)`**. The API message is captured in the same node from U131's no-solve calls on `_study(degree=1)` (`:122–135`). The refusals:
    - `thermal.duty` (the malformed file);
    - `thermal.initial` (a variant);
    - `outputs.probes` (a variant beside the A file's transient);
    - `thermal` on `malformed_transient_on_port.yaml` (the landed four-port reason);
    - a `dt_s` in place of `t_end_s`, at `thermal.transient`.
    - On the two tracked files, `cli.main(["run", …])` and `cli.main(["check", …])` read rc (2, 2) on both ranks and nothing is written (`test_wf1_yaml_thermal.py:142–162`'s form).
    - Build calls 0; `study.result is None`.
- **W1, the solve** (`test_wf1_s29_sphere_transient_reproduces_u131`): `run_study_file` on A and on S.
  - **(UT)** Sub cells = `MAT10_SUB_CELLS` on both. **(Once)** `run.study.result is run.result` on both.
  - **(A_t)** On A, `_interior_rel(sub, chi, θ − θ_cf, θ_cf)` ≤ `A_BAND`. Here θ_cf = `_exact_theta(sub, T_STAR, q_cf)`, with q_cf built as U131 builds it (`:181–185`). **(R-A_t)** Its `.4e` == `MAT12_A_T_PRINT`.
  - **(S2)** On S, |ρc·θ̄¹/(Δt₁·q̄_h) − 1| ≤ `S2_BAND` (the volume mean, `:194–200`). **(R-S2)** Its `.4e` == `MAT12_S2_PRINT`.
  - **(E)** `summary["temperature"]["max_balance_residual"]` ≤ `E_BAND` on both.
  - **(IT)** In the same window, solve `_study(degree=1)` independently and call `temperature(RobinCooling(K, H, 0.0, source_scale=G), transient=<U131's spec>)`. `_owned_rel` against the file θ is ≤ `L0_BAND` on A and on S (step 3's (IT), `test_wf1_yaml_thermal.py:230–238`).
  - **(CLI = API by `repr`)**
    - `summary.json` on disk equals `run.summary` key by key, `elapsed_s` excluded. The JSON round trip is exact by `repr` (`WF-1.md:50`).
    - `summary["temperature"]`'s `max_balance_residual`, `source_power_w` and `source_multipliers` equal the run's own `TransientStudyResult` fields by `repr`.
    - `cli.main(["run", <A>, "--output-dir", <dir>])` returns rc (0, 0). Its `summary.json` has the API run's key set, with every number within `L0_BAND`. This is a second solve, so it is never compared with `==` (`WF-1.md:359`).

**Control (asserted, backed).** **(c-qs)** On the solved A Study, the steady θ from `run.study.temperature(cooling.cooling)` (no re-solve; `once` still True) stands in for θ(t\*) and misses `A_BAND`. U131 read 1.0108e+00 for this control (`20261107T013321Z_WF-15-s15-W1.log:141`), ≈ 9× the band. Assert `err_qs > A_BAND`, the gate's own band (rule (e)).

**Printed (no band).**
- Both `.10e` errors beside U131's 3.4708293773e-02 (`…s15-W1.log:138`) and ratio 0.999950006 (`:139`).
- Each file's parsed dt in `float.hex`, beside `T_STAR / N_MAIN` and `DT_1`.
- Per-run wall time.

**Windows.** Complex build. Every window is bare and foreground, Bash timeout 660000 ms. `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0 before W1 and before each Wreg.
- **W0:** `scripts/testing/run_and_log.sh WF-1-s29-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 120 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf1_yaml_sphere_transient.py -v -s --tb=short'"`. The solve node skips.
- **W1:** W0's command with label `WF-1-s29-W1`, `FEM_EM_WF1_S29=1` added to the env and `timeout -k 30 300`. Wrap it in the durable capture of `20261104T030134Z_IMP-8-I55-W2-item803.log:12` (redirect to `/workspace/logs/wf1-s29-w1-raw.log`, the `[capture] rc=` line, `cat`), including the trailing `; exit $rc`.
  - **Cost from measurement:** U131's W1 took 30 s for one EM solve and five thermal runs (EM 15.8 s; `…s15-W1.log:142`, `:151–152`). This window runs four EM solves (A, S, the Python twin, the CLI re-run of A). Predicted 70–120 s, inside `-t 300`.
  - A kill parks the item. Never raise the timeout.
- **Wreg-c** (complex, `-n 2`, `timeout -k 30 180`, knobs unset):
  - `tests/package/test_wf1_yaml_thermal.py` whole. The steady sphere YAML must be unchanged; it read 6 passed in 33 s (`20261009T020735Z_WF-1-s6-W4-thermal.log.gz:12`, `:230–231`).
  - `tests/package/test_wf15_study_sphere_transient.py::test_wf15_s15_sphere_transient_refusals`.
- **Wreg-p** (real, `-n 2`, `timeout -k 30 180`): the parsers that call `_parse_thermal(base, "sphere_in_box", …)`, in the s28 Wreg-r form (`20261105T210709Z_WF-1-s28-Wreg-r.log:12`):
  - `tests/package/test_wf1_yaml_thermal_ball.py::test_wf1_s14_parse_and_refusals`
  - `tests/package/test_wf1_yaml_thermal_box.py::test_wf1_s15_parse_and_refusals`
  - `tests/package/test_wf1_yaml_thermal_box_wall.py::test_wf1_s28_parse_and_refusals`
  - `tests/package/test_wf1_yaml_f_human_thermal.py::test_wf1_fht_parse_and_refusals`
- **Then Wgen, Wval, Wapi, Wdoc,** in that order: `WF-1.md:2295–2298`'s four commands, labelled `WF-1-s29-W…` (Wgen before Wapi, item 653).

**Tier.** Standard, `-n 2`, complex. **Class:** chain step (solve). **Cost:** ≈ 6.8 measured / 30 predicted min.

**Traps already paid for.**
- **Degree.** The `Study`'s default EM degree is 2. Both files say `study.degree: 1`, because U131's records are degree 1 (U132's degree-2 records are not this step's).
- **T∞.** `MAT-12` runs T∞ = T(0) = 0. Both files carry `t_inf_c: 0.0`, never `MAT-10`'s 20.0 from `sphere_mat10_thermal.yaml`.
- **dt.** `t_end_s / n_steps` is the same division U131 makes (`T_STAR / N_MAIN`; `DT_1 / 1`), so dt is bitwise U131's. Write `t_end_s` as the `repr` and assert the parsed dt by `float.hex` (U125's one-ULP trap, item 830).
- **Plumbing.** Do not edit `_load`'s 4-tuple; the carrier is `_BallThermal`. `src` imports nothing from `tests/`.
- **Probes.** The API refuses `Study.probe` with `transient` on the sphere (U131). The YAML refusal must fire at parse, before any mesh.
- **U135 lands first.** If U135 (`WF-15` step 17) has already lifted the API's duty refusal on the sphere, the YAML still refuses `thermal.duty` with the text quoted above (no sphere duty key is commissioned). The outcome names the split for the weekly.
- **MUMPS.** `-n 2` MUMPS is not bitwise across solves. (IT) and the CLI re-run compare at `L0_BAND`; `repr` equality applies only to specs and to one run's own JSON round trip.
- **`elapsed_s`.** Rank 0 writes `summary.json`, and `elapsed_s` is each rank's own. Exclude it (item 830's first W1, `20261105T210601Z_WF-1-s28-W1.log`).
- **Estimate table.** No new YAML lands under `examples/` or `docs/guide/`. The estimate table discovers study files there (`scripts/docs/generate_estimate_page.py:56`, `:107–117`) and is pinned at 14 entries (`tests/package/test_wf17_estimate_table_page.py:47`, `:67–72`).
- Count verdicts per rank, by token; use `-s`; reduce before asserting.

**Scope.** The YAML route to U131's two runs: degree-1 EM, the `MAT-4` / `MAT-10` sphere, CN N = 40 to t\* and one BE step. **Not claimed:** degree 2 through YAML, duty / gradient / probes on the sphere, the tip ball, any band, U135's duty lift.

**Negative result.**
- **(IT) red with (E) green:** the YAML route does not build U131's transient. Print the per-ingredient comparison (cells, facets, dt `float.hex`, θ, startup, q_h max\|Δ\|), park on `attempt/WF-1-s29-*`, BLOCKED.
- **(R-A_t) / (R-S2) red with (IT) green:** the sphere's mesh or source moved since U131. Print both, park, BLOCKED: a finding for the weekly, never a re-record.
- **Wreg red:** the edit reached a landed path. Revert it.
- No `PREDICTED RED:` line.

**Status it can move.**
- `WF-1` gains step 29, a record on its ✅ row.
- §10 U134 becomes DONE.
- `MAT-12` joins `run_study_file`'s envelope (`docs/guide/validation.md` regenerated).
- The public key `thermal.transient` (`t_end_s`, `n_steps`, …) on a `sphere_in_box` goes to the 11-09 weekly to ratify.
- Item 887 (U136) is unblocked.

**Depends on:** nothing (U131 is on `main`).

**Outcome (06:00 slot, 2026-11-08): BLOCKED — the brief's negative result "(IT) red with (E) green"; parked on `attempt/WF-1-s29-20261108T1224Z` (`ba30efd2`).** W0 `20261108T122049Z_WF-1-s29-W0.log`: 5 passed, 1 skipped per rank (2.2 s) — (K) the parsed `RobinCooling` / `TransientSpec` `repr`s and dt `float.hex` equal U131's (A `0x1.bdddddddddddep+1`, S `0x1.c891457ed6e76p-8`), `femem check` rc (0, 0); (Ref) duty / initial / probes refused with reasons equal to the API's `str(exc)`, the four-port file at `thermal` with the landed reason, `dt_s` at `thermal.transient.dt_s` (`_unknown`'s field form), rc (2, 2), nothing written, build calls 0. W1 `20261108T122102Z_WF-1-s29-W1.log` (73 s) and W1b with the per-ingredient print `20261108T122244Z_WF-1-s29-W1b.log` (73 s): 5 passed, 1 failed per rank. Green: (UT) 19 496 on both, (Once), (A_t) 3.4708293773e-02 (W1) / 3.4708293770e-02 (W1b) beside U131's 3.4708293773e-02, (R-A_t), (S2) 4.9993987957e-05 (ratio 0.999950006, U131 0.999950006), (R-S2), (E) ≤ 4.8e-15, (c-qs) 1.0108e+00 > `A_BAND`, CLI = API (JSON round trip by `repr`, summary fields = the run's own by `repr`, CLI rc (0, 0), re-run within `L0_BAND`). **Red: (IT) S 4.500e-09 (W1) / 1.122e-09 (W1b) > `L0_BAND` 1e-9; A 1.130e-10 / 3.781e-11 green.** Per ingredient (W1b `:310–311`): cells 19 496 / 19 496, Robin facets 3 086 / 3 086, dt `float.hex`, θ and startup identical on A and S; only q_h differs, owned max rel 3.47e-09 (A) / 9.52e-10 (S) between the two independent `-n 2` MUMPS EM solves, so the one BE step's θ inherits that pointwise noise (A's 40 CN steps diffuse it). The two YAML runs' own (A_t) `.10e` differ in the last digit for the same reason. This reads as the comparand of (IT) rather than the YAML route: the route builds U131's transient ingredient for ingredient. Band not touched. **Unblock:** a review re-poses (IT) — e.g. against `run.study.temperature(cooling, transient=<U131's spec>)` on the YAML's own solved Study (no second EM solve), or with a band backed by this measurement — then W1 and the Wreg / Wgen / Wval / Wapi / Wdoc windows are re-run (none ran). Disclosed on the branch: `summary["temperature"]["transient"]["t_end_s"]` is `dt_s · n_steps` (the `_BallThermal` carrier holds no t_end); the CLI re-run compares `max_balance_residual` as both ≤ `E_BAND`, not relatively. U135 has not lifted the API's duty refusal (W0 captured it live).

### §9 item 878 — 2026-11-08 (§10 U139, chain U — 9a, WF-1 step 30, the dose readout through YAML; chain step (solve), standard, real)

*Written by the 2026-11-08 03:00 daily review from the 11-07 weekly's enumerated §10 U139 (rule (4); not provisional). The target, anchor and control are the weekly's; the execution detail is the review's. Family U139, attempt 1 of 4. Independent of U134 / U135. Next: item 888 (U140).*

**Executes** §10 U139. A `thermal: {dose: …}` readout (`MAT-20`'s `cem43` and six-minute window) on U125's thermal-box files. It asserts (Q-wf), P68's dose gate, at `DOSE_BAR` on P68's mesh and Δt (the walled box with heat capacity), plus (E). CLI = API by `repr`. A wall named without a thickness is refused. Envelope entry.

**How the target is read (execution detail; the target is unchanged).**
- **(Q-wf) is P68's gate:** `cem43` on the FEM point history against `_q_check_g`'s quadrature at `DOSE_BAR` (`tests/validation/test_mat21_wall_pre_equilibration_duty25_dose_fem.py:287–292`).
- **P68's mesh is item 682's walled box:** 119 808 cells, (52, 32, 12) (`20261105T224651Z_MAT-21-P68-W1.log:58`). **P68's Δt** is 70.5375 s (`:53`).
- **Why U125's files:** P68's own runs are D = 0.25, and a wall beside a duty is a ratified refusal (11-07 weekly; `study.py:1905–1906`, `study_yaml.py:1061–1062`). So the gate runs on U125's continuous walled files: R1 is P63's R1, R0 its uniform-start twin. They use that mesh and Δt (`20261105T210638Z_WF-1-s28-W1-r2.log:90–92`). P68's D = 0.25 dose values are not reachable through YAML and are not claimed.
- **The dose form is new.** No landed YAML key carries `cem43`. The only six-minute form in YAML is the F-human `thermal.window_s` (`study_yaml.py:808–813`, U46's), not U57 / U58's: U57 is `capacitor_sheets` (item 496) and U58 the `thermal_ball` duty / wall files (item 500).

**Code** (`src`, rule (c), additive; with `dose=None` every path is byte-identical).
- **`study.py`:**
  - **New `@dataclass(frozen=True) class DoseReadout`**, next to `TransientSpec` (`:2041`):
    - fields `points_m: tuple` (a non-empty tuple of finite (x, y, z)) and `window_s: float = 360.0`;
    - the comment cites `running_average`'s default (`thermal/dose.py:108`) and `_SIX_MINUTE_WINDOW_S` (`:3352`);
    - `__post_init__` validates both fields;
    - exported in `fem_em_solver.__all__` (`src/fem_em_solver/__init__.py:54`, `:61`).
  - **`_study_temperature`** (`:1601`) gains `dose=None`. A `dose` with the tuned two-port hook, or on any geometry except `ThermalBoxGeometry`, raises `ValueError("dose is defined on a transient ThermalBoxGeometry only (WF-1 step 30)")` before any solve. Otherwise the dose passes to `_study_box_temperature`.
  - **`_study_box_temperature`** (`:1898`): a `dose` with `transient=None`, with a `TransverseQuadraticHeatSource`, or beside a `duty` raises `ValueError` naming the case. Otherwise it passes to `_study_transient_temperature`.
  - **`_study_transient_temperature`** (`:2125`) gains `dose=None`. With a dose:
    - `extra["snapshot_steps"] = tuple(range(n_steps + 1))`, through `solve_transient_conduction` (`thermal/transient.py:110`, `:288–294`);
    - a dose point outside the gel box [−A, A] × [−B, B] × [−C, C] raises `ValueError`;
    - after the solve, read each snapshot on **one fixed batch** (the dose points, in order) with `post.evaluation.evaluate_vector_field_parallel` (`post/evaluation.py:224`; P68's `_read` route). A point not found raises `ValueError`. Then clear the snapshots (P68's `:193–200`);
    - per point: `cem43(times, T[:, j])` and `running_average(times, T[:, j], window=window_s)` (`thermal/dose.py:45`, `:108`), with `times = dt * np.arange(n_steps + 1)`.
  - **`TransientStudyResult`** (`:2090`) gains a last, defaulted field `dose: Optional[dict] = None`, in plain floats:
    - `points_m`, `times_s`;
    - `temperature_c` (n_t × n_pts);
    - `cem43_s` (one per point);
    - `six_minute: {window_s, peak_c, t_max_s}` (per point).
  - **Docstrings** of `_study_temperature`, `TransientStudyResult` and `DoseReadout` name `WF-1` step 30 and `MAT-20`. No dose limit, threshold or compliance wording (`thermal/dose.py:18`).
- **`study_yaml.py`:**
  - **`_parse_thermal_box_mapping`** (`:968–1158`):
    - `"dose"` joins the allowed `thermal` keys (`:1052`). The section is `thermal.dose: {points_m: [[x, y, z], …], window_s: w}`;
    - `points_m` is required, parsed with `_triples` (`:1167`). `window_s` is optional, parsed with `_positive_finite` (`:832`); omitted means `DoseReadout`'s own default;
    - refused at `thermal.dose` without `thermal.transient`, beside `thermal.duty`, or with a `transverse_quadratic` source;
    - a point outside the gel is refused at `thermal.dose.points_m[i]`;
    - beside `outputs.probes`, refused at `outputs.probes`.
    - The parsed readout rides in `outputs["dose"]`, as `tune` does; `_load`'s 4-tuple and `_BallThermal` stay unedited.
    - Geometry parses first, so a wall without a thickness is refused at `geometry.wall.thickness_m` before the dose is read.
  - **Other files:** on `thermal_ball` and `sphere_in_box` files, `thermal.dose` stays an unknown key, refused by the landed `_unknown`.
  - **`run_study_file`** (`:2242–2247`): a `ThermalBoxGeometry` with `"dose" in outputs` goes to a new `_run_thermal_box_dose`. It copies `_run_thermal_box`'s lines (`:1904–1976`), passes `dose=outputs["dose"]` to `study.temperature`, and adds `summary["dose"] = thermal.dose`. `_run_thermal_box` stays unedited.
  - **Module docstring:** the `thermal_box` bullet (`:63–69`) names `thermal.dose`.
- **New files** in `tests/package/wf1/`:
  - `thermal_box_wall_gradient_m4_dose.yaml` and `thermal_box_wall_m4_dose.yaml`: U125's R1 / R0 files byte for byte, plus `thermal.dose: {points_m: [P5, ADJ_PT], window_s: 360.0}`. Each coordinate is the `repr` of P68's `HIST_PTS` (`…duty25_dose_fem.py:104`; `ADJ_PT`'s z is `0.9 * C_BOX`, `tests/validation/test_mat21_pre_equilibration_closed_form.py:59`).
  - `malformed_thermal_box_dose_wall_without_thickness.yaml`: the R1 dose file without `geometry.wall.thickness_m`.
- **New module `tests/package/test_wf1_yaml_thermal_box_dose.py`**, knob `FEM_EM_WF1_S30`, real build.
  - **Imports, by name:**
    - from `tests.package.test_wf1_yaml_thermal_box_wall`: `R1_FILE`, `R0_FILE`, `CELLS`, `DT`, `_field`, `_with`, `_without`, `_py_study`, and the direct-R1 names it imports (`_walled_box`, `_materials`, `_rho_c_field`, `_robin_facets`, `ROBIN`, `RUNG`, `A_BOX`, `B_BOX`, `C_BOX`, `D_WALL`, `H_GEL`, `K`, `Q`, `RHO_C`, `T_INF`, `G_GRAD`, `N_MAIN`, `IDENTITY_BAND`, `E_BAND`, `I_P5`);
    - from P68's module: `HIST_PTS`, `HIST_LABELS`, `GRID_BAR`;
    - `DOSE_BAR` and `T_GRID` from `tests.validation.test_mat21_box_dose_closed_form`;
    - `_q_check_g` and `_win_check_g` from `tests.validation.test_mat21_box_dt_duty125_dose_closed_form` (`:131`, `:139`);
    - `AVG_BAR` and `W_B2` from `tests.validation.test_mat20_running_average`;
    - `_dose` from `tests.validation.test_mat21_pre_equilibration_duty_dose_closed_form` (`:56`);
    - `_read` from `tests.validation.test_mat21_wall_pre_equilibration_fem` (`:160`);
    - `LOGS`, `WF1`, `_spec_diff`, `_mapping`, `_fresh_dir`, `_rc_all` from `tests.package.test_wf1_study_yaml`;
    - `cem43` and `running_average` from `fem_em_solver.thermal.dose`.
  - **Helpers for item 888:** the direct-R1 history and the (Q-wf) / (Win-wf) collection are module-level helpers (`_direct_r1_history`, `_dose_gates`), never `test_*` functions, so item 888 can import them.
  - **Docstring:** backticks `WF-1`, `MAT-20`, `MAT-21`, and cites P68's and U125's logs.
- **Envelope:**
  - a new `[names.DoseReadout]` entry with both nodes;
  - the W1 node under `[names.run_study_file]`, `[names.ThermalBoxGeometry]` and `[names.TransientSpec]`;
  - the W0 node under `[names.load_study]`;
  - each with `# WF-1 step 30 (§9 item 878)`.
  - Wval regenerates `docs/guide/validation.md` (`MAT-20` becomes a gating chunk). Wgen writes `docs/guide/api/DoseReadout.md`.

**Gates (asserted, bands imported, none set).**
- **W0, no mesh** (`test_wf1_s30_parse_and_refusals`, `-n 2`):
  - **(K)**
    - With their `dose` lines removed, both dose files equal U125's files' text (custody).
    - `_spec_diff` against U125's parsed Studies is `[]`.
    - The parsed `DoseReadout`'s `repr` equals `repr(DoseReadout(tuple(map(tuple, HIST_PTS)), W_B2))`, and each coordinate matches by `float.hex`.
    - `femem check` reads rc (0, 0).
  - **(Ref)**
    - The malformed file is refused at `geometry.wall.thickness_m`; `femem run` and `femem check` read rc (2, 2) and nothing is written.
    - Variants, each refused at its field: a dose without a transient; a dose beside a duty (unwalled-box variant); with the quadratic source; with a point outside the gel; beside probes; on a `thermal_ball` file.
    - API: `dose=` raises `ValueError` on the steady box, on the ball, on the EM sphere (`_study(degree=1)` from `tests.package.test_wf15_study_lossy_sphere`, no solve), beside a duty, and with the quadratic source.
    - Build calls 0.
- **W1, the solve** (`test_wf1_s30_dose_reproduces_p68_gate`). Four runs:
  - **direct R1:** as `test_wf1_yaml_thermal_box_wall.py:330–338`, with `snapshot_steps=tuple(range(N_MAIN + 1))`; each snapshot read by `_read` on `HIST_PTS` (one batch); history `T_INF + θ`;
  - `run_study_file` on both dose files;
  - **the Python twin:** `_py_study().temperature(RobinCooling(K, H_GEL, T_INF), transient=TransientSpec(DT, N_MAIN, RHO_C, initial_gradient_k_per_m=G_GRAD), dose=DoseReadout(…))`.
  
  The gates:
  - **(Grid)** `summary["dose"]["times_s"]` against `T_GRID`: max \|t − T_GRID\| / DT ≤ `GRID_BAR` (P68's W0 form; it read 3.553e-15 steps, `…P68-W1.log:53`). **(Cells)** 119 808 on every build.
  - **(Q-wf)** On R1, R0 and the Python run, at P5 and `ADJ_PT`: `_q_check_g(T_hist[:, j], times)` rel ≤ `DOSE_BAR`. P68's six histories read ≤ 5.374e-16 (`…P68-W1.log:78–83`).
  - **(Win-wf)** `_win_check_g(T_hist[:, j], times, window=W_B2)`: rel ≤ `AVG_BAR`, and lo ≤ peak ≤ hi (item 708's form).
  - **(I-h) / (I-d)**
    - The R1 file's history against the direct one: max rel ≤ `IDENTITY_BAND` (U125's (I) read 0.000e+00, `…s28-W1-r2.log:94`).
    - Its `cem43_s` at each point against `_dose` of the direct history: rel ≤ `IDENTITY_BAND`.
    - The same two checks on the Python run.
  - **(R)** The R1 history's last row at P5, minus `T_INF`, equals U125's θ(P5) +1.224895450933e+00 K (`…s28-W1-r2.log:91`) within 5e-13: half the print's last digit (§9's standing ruling on I41, `IMP-3.md`).
  - **(E)** `max_balance_residual` ≤ `E_BAND` on the direct run, both files and the Python run.
  - **(CLI = API by `repr`)**
    - Recompute `cem43(times, T[:, j])` and `running_average(…)` in the test from the run's own `thermal.dose` history. They equal `summary["dose"]["cem43_s"][j]` and the `six_minute` peak / t_max by `repr`.
    - `summary.json` on disk equals `run.summary`, `elapsed_s` excluded.
    - Only the dose files' summaries carry `dose`.

**Control (asserted at the gate's band, rule (e)).** **(c-init-d)** The R0 file's `cem43_s` at P5 misses R1's by more than `DOSE_BAR`. This is P68's control form (`…duty25_dose_fem.py:293–299`); its D = 0.25 miss read 4.6009e-01 (`…P68-W1.log:85`). No log backs a size on the continuous files, so the miss is printed and asserted only > `DOSE_BAR`. **Predicted, printed, not asserted:** positive, and at least of order U125's θ(P5) rel miss 1.493e-01 (`…s28-W1-r2.log:95`), since the dose is exponential in T.

**Printed (records, no band).**
- D at P5 and `ADJ_PT` for R1 and R0, in P68's `.6e` on its peak scale (`_dose` is `cem43` of T_INF + θ on `T_GRID`, no normalisation).
- The R-swap miss, as P68 prints it.
- The six-minute peak and t_max.
- The read time per run.

**Windows.** Real build. Every window is bare and foreground, Bash timeout 660000 ms. `pgrep -c python3` in the container must read 0 before W1, Wreg-r and Wreg-c.
- **W0:** `scripts/testing/run_and_log.sh WF-1-s30-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 120 mpiexec -n 2 python3 -m pytest tests/package/test_wf1_yaml_thermal_box_dose.py -v -s --tb=short'"`.
- **W1:** W0's command with label `WF-1-s30-W1`, `FEM_EM_WF1_S30=1` added to the env and `timeout -k 30 300`. Wrap it in the durable capture of `20261104T030134Z_IMP-8-I55-W2-item803.log:12`, including the trailing `; exit $rc`.
  - **Price:** U125's W1-r2 took 16 s for four runs (`…s28-W1-r2.log:106–107`). P68's runs took 2.7–3.5 s to solve plus 4.4–4.7 s for 41 snapshot reads each (`…P68-W1.log:60–64`). Four runs with histories: predicted 40–90 s.
- **Wreg-r1** (real, `-n 2`, `timeout -k 30 600`): `20261105T210709Z_WF-1-s28-Wreg-r.log:12`'s command, verbatim (12 passed, 71 s, `:267–268`).
- **Wreg-r2** (real): U125's module with `FEM_EM_WF1_S28=1`, the `20261105T210638Z_WF-1-s28-W1-r2.log:12` command (16 s). U125's runs must read unchanged.
- **Wreg-c** (complex, `-n 2`, `timeout -k 30 600`): `20261105T210828Z_WF-1-s28-Wreg-c.log:12`'s command (204 s, `:83–84`), plus `tests/package/test_wf15_study_sphere_transient.py::test_wf15_s15_sphere_transient_refusals`.
- **Then Wgen, Wval, Wapi, Wdoc,** in that order: `WF-1.md:2295–2298`, labelled `WF-1-s30-W…`.

**Tier.** Standard (W1), `-n 2`, real (Wreg-c complex). **Class:** chain step (solve). **Cost:** ≈ 6.8 measured / 35 predicted min.

**Traps already paid for.**
- **Real build.** The walled-box records are real-build (U125, P63, P68). The W1 node asserts `not tests.complex_mode.IS_COMPLEX`, as `test_wf17_thermal_box_wall_tutorial.py:171` does.
- **The window default.** `_SIX_MINUTE_WINDOW_S` is defined at `study.py:3352`, after `TransientSpec`. A dataclass default that names it at class creation raises `NameError`. Write `360.0` with the citation.
- **Batch bytes.** `.eval` bytes depend on batch composition (item 537). The direct run and the Study read the same two-point batch, call by call, per snapshot.
- **Snapshots.** Hold them, read them on one batch, then clear them (P68's `:193–200`).
- **`cem43` input.** It refuses non-increasing times and a non-zero imaginary part (`thermal/dose.py:65–80`). Pass `np.real` values.
- **P5's position.** P5 sits on the gel–wall interface z = C, a mesh plane. Read it as P68 did; never perturb it.
- **The duty refusal.** A wall beside a duty is a ratified refusal. Do not lift it to reach P68's D = 0.25 runs.
- **Not a probe.** `Study.probe` refuses a walled box (`study.py:2229–2230`). The dose readout is a point history, not a probe, and does not route through `_study_probe`.
- **Estimate table.** No new YAML lands under `examples/` or `docs/guide/` (the estimate-table discovery trap, item 873).
- **Bookkeeping.** `src` imports nothing from `tests/`; the copied lines cite their span. `elapsed_s` is per rank. Count verdicts per rank, by token; use `-s`. Run Wgen before Wapi.

**Scope.** The 1D-reduced walled box (R1 / R0, continuous source): two points, one rung, one Δt, from the YAML and from the API. **Not claimed:** P68's D = 0.25 dose values (wall + duty is refused); a dose on the ball, the EM sphere, F-human, the unwalled box or under a duty; a dose band, limit or threshold; compliance; a cited PMMA ρc (W3).

**Negative result.**
- **(Q-wf) / (Win-wf) red:** `cem43` / `running_average` on this FEM history disagree with their quadratures, a finding about `MAT-20`'s arithmetic. Print, park on `attempt/WF-1-s30-*`, BLOCKED.
- **(I-h) / (I-d) red with (E) green:** the Study's history is not the direct one. Print the per-ingredient comparison (cells, Σk, Σρc, facet count, θ₀ at P5, the batch, the `times` hex), park, BLOCKED.
- **(R) red with (I-h) green:** U125's record moved, a finding for the weekly, never a re-record.
- **Wreg red:** revert the edit.
- No `PREDICTED RED:` line.

**Status it can move.**
- `WF-1` gains step 30.
- §10 U139 becomes DONE.
- `DoseReadout`, `Study.temperature(…, dose=)` and the YAML `thermal.dose` keys go to the 11-09 weekly to ratify.
- `MAT-20` joins the `run_study_file` / `DoseReadout` envelope.
- Item 888 (U140) is unblocked.

**Depends on:** nothing (U125's files and P68's module are on `main`).

**Outcome (09:00 slot, 2026-11-08): DONE — `WF-1` step 30, a record on the ✅ row.** W0 `20261108T150630Z_WF-1-s30-W0.log`: 1 passed, 1 skipped per rank (1.8 s) — (K) custody (dose files = U125's text plus the three dose lines), `_spec_diff` [], readout `repr` = `DoseReadout(tuple(map(tuple, HIST_PTS)), W_B2)` with coordinates by `float.hex`, `femem check` rc (0, 0); (Ref) the malformed file at `geometry.wall.thickness_m`, run / check rc (2, 2), nothing written; eight YAML variants at their fields; eight API refusals; build calls 0. W1 `20261108T150754Z_WF-1-s30-W1.log` (durable capture, 33 s, `[capture] rc=0`): 2 passed per rank. (Grid) 0.0 steps on all three; (Cells) 119 808 on every build; (Q-wf) rel ≤ 5.715e-16 and (Win-wf) rel ≤ 6.713e-16, inside [lo, hi], on R1 / R0 / Python at P5 and `ADJ_PT`; (I-h) 5.079e-16 (file) / 6.716e-16 (Python), (I-d) ≤ 3.441e-15; (R) θ(P5) +1.224895450933e+00 K, \|Δ\| 4.499e-13 ≤ 5e-13; (E) ≤ 2.34e-13 on five runs; CLI = API by `repr` on all three; `summary.json` = `run.summary`; only the dose files carry `dose` (U125's R1 file run alongside, no `dose`). Control (c-init-d): R0 D 9.046411e-11 s vs R1 1.284092e-10 s at P5, miss 2.9550e-01 > `DOSE_BAR`; predicted "of order 1.493e-01": the miss is larger, signed R0/R1 − 1 = −2.955e-01 (R0 colder, so lower dose). D printed: R1 P5 1.284092e-10 / `ADJ_PT` 1.352102e-10; R0 9.046411e-11 / 9.596378e-11 s; six-minute peak 21.1688 / 21.2277 °C (R1), t_max 2821.5 s. Wreg-r1 `20261108T150856Z_WF-1-s30-Wreg-r1.log` 12 passed (67 s); Wreg-r2 `20261108T151015Z_WF-1-s30-Wreg-r2.log` 2 passed, U125's θ(P5) +1.224895450933e+00 / +1.041997184385e+00 and drift unchanged (14 s); Wreg-c `20261108T151040Z_WF-1-s30-Wreg-c.log` 2 passed (187 s); Wgen / Wval regenerated `api/DoseReadout.md` and `validation.md` (`MAT-20` now gating); Wapi `20261108T151545Z_WF-1-s30-Wapi.log` 22 passed; Wdoc `20261108T151631Z_WF-1-s30-Wdoc.log` 0 dead. **Disclosed:** (1) `DoseReadout.__post_init__` stores coordinates and `window_s` as plain floats — the brief's literal `tuple(map(tuple, HIST_PTS))` holds `np.float64`, whose `repr` differs (first W0 `20261108T150611Z_WF-1-s30-W0.log`, red on that `repr` only; same bits); (2) `tests/package/public_surface.txt` gains `DoseReadout` (the first Wapi `20261108T151427Z_WF-1-s30-Wapi.log` read 4 failed on the snapshot alone); (3) the first W1 `20261108T150700Z_WF-1-s30-W1.log` was green; one control print was relabelled (it said "positive" for a non-zero check) and W1 re-ran, so rule (i) holds on the second; (4) `AVG_BAR` is `test_mat20_running_average.BAR`, imported under that alias as precedent does; (5) W1 also runs U125's R1 file (no dose) for the "only the dose files" check; (6) the Study's read time is not separable from its solve in the summary — printed per run as solve + reads (6.6–6.9 s against the direct 3.1 + 4.2 s). Bands untouched.

**Review ruling (2026-11-08 10:30 daily review) on item 873's red — a re-pose under §9 rule (6), provisional; the 11-09 weekly ratifies or reverts.** The red is (IT) on S alone: 4.500e-09 / 1.122e-09 against `L0_BAND` 1e-9 (`20261108T122102Z_WF-1-s29-W1.log`, `20261108T122244Z_WF-1-s29-W1b.log`). Its comparand is a **second, independent `-n 2` EM solve**. W1b's per-ingredient print (`:310–311`) shows cells, Robin facets, dt, θ and startup identical. q_h differs by an owned-max relative 3.47e-09 (A) / 9.52e-10 (S), while the integrated source power agrees to 2.9e-12 / 7.4e-13. One BE step carries the pointwise difference into θ; A's 40 CN steps diffuse it (A green at 1.1e-10 / 3.8e-11).

**The `log-pathologist` (this review, journal §4) REFUTED the slot's reading as stated:**
- **What was asserted.** (c-qs) and CLI = API **were never asserted** at attempt 1, because they sit after the (IT) assert. Only (UT), (Once), dt, (E), (A_t), (R-A_t), (S2) and (R-S2) were asserted green, on both ranks: (A_t) 3.4708293773e-02, and the (S2) ratio 0.999950006, equal to U131's.
- **The cited floor.** `SOLVE_FLOOR_RECORD` records integrated quantities, so it does not cover a pointwise owned-max spread. 3.47e-09 is ≈ 78× its upper end.
- **What stays open.** Whether the spread is a solver floor or a route difference is not established. The best in-log evidence for "floor" is that the YAML route's own A and S solves of one EM spec also differ in integrated power, by 3.6e-12.

Known-issues (the `-n 2` MUMPS observation, new row 2026-11-08) records this.

**The re-pose separates the readings and moves no band:**
- **(IT-run)**, asserted at `L0_BAND` on A and on S: the YAML route's θ against the **same run's** Study, with no second EM solve. This tests the route by itself.
- **(IT-ind)**, the independent Study: asserted at `L0_BAND` on A, as before. **On S it is printed**, beside both q_h owned-max spreads.
- **(Q-rep)**, the decider the pathologist named, printed: `_owned_rel(q_h(A run), q_h(S run))` between the YAML route's two runs. These are two same-route solves of one EM spec in one process.
  - **If (Q-rep) ≥ 1e-9**, the pointwise spread is a floor of the solve, not of the route. Record it.
  - **If (Q-rep) < 1e-10 while (IT-ind) S stays > 1e-9**, the route and the Python twin build different q_h. That is a finding: park.
  - Between those values, record and stop.
- **The weekly's call:** whether S's independent read takes a band backed by these measurements.

**Since attempt 1:**
- U135 (item 877) made `Study.temperature` accept `duty` on the sphere. The API's "duty is not defined on the EM sphere" message no longer exists; `study.py` now raises "duty needs a transient". The YAML refusal of `thermal.duty` can therefore no longer be compared with an API capture.
- U139 (item 878) edited `study_yaml.py` and the envelope.

### §9 item 891 — 2026-11-08 (§10 U134, chain U — 9a, WF-1 step 29 re-posed; chain step (solve), standard, complex)

*Written by the 2026-11-08 10:30 daily review under §9 rule (6) (a re-pose of a blocked step, keeping its target). Family U134, attempt 2 of 4. Base: `attempt/WF-1-s29-20261108T1224Z` (`ba30efd2`). Item 873's brief (above) binds except where the ruling above and this brief say otherwise. Depends on nothing. Next: item 887 (U136).*

**Code onto `main`:**
1. **By path (new files):**
   - `tests/package/test_wf1_yaml_sphere_transient.py`;
   - the four `tests/package/wf1/` YAML files (`sphere_mat12_transient_a.yaml`, `…_s.yaml`, `malformed_sphere_transient_duty.yaml`, `malformed_transient_on_port.yaml`);
   - attempt 1's three logs (`…122049Z_WF-1-s29-W0.log`, `…122102Z_…-W1.log`, `…122244Z_…-W1b.log`). Append their three test-results rows **by hand**, insertions only.
2. **By three-way apply, never by path checkout** (U139 edited both files since): write `git diff 5054e5fd ba30efd2 -- src/fem_em_solver/study_yaml.py tests/package/validation_envelope.toml` to a scratch file in the repo root, then run `git apply -3` on it and delete the scratch file. A conflict that is not a pure adjacency means stop BLOCKED. After the apply, U139's `thermal.dose` lines and U135's envelope lines must still be present (`git diff HEAD` shows insertions only in both files).

**The edits this attempt makes** (test module and `study_yaml.py` only):
- **(IT) → (IT-run) + (IT-ind)**, as the ruling above:
  - **(IT-run):** `_owned_rel(<file θ>, run.study.temperature(RobinCooling(K, H, 0.0, source_scale=G), transient=<U131's spec>).<θ>)` ≤ `L0_BAND` on A and on S, with `run.study.result is run.result` asserted before and after. This is no second EM solve.
  - **(IT-ind):** the independent `_study(degree=1)` solve, as attempt 1 ran it. On A it is asserted ≤ `L0_BAND`. On S it is printed beside `_owned_rel(q_h_file, q_h_twin)`, with the label "pointwise cross-solve spread (known-issues `-n 2` MUMPS observation, 2026-11-08 row); a band is the weekly's". Keep W1b's per-ingredient print.
  - **(Q-rep):** print `_owned_rel(q_h(A run), q_h(S run))` (`.3e`) and the integrated-power difference between the YAML route's two runs. The decision rule is the ruling's.
- **The `thermal.duty` refusal.** In `_parse_thermal_sphere_transient`, replace the quoted U131 text with a YAML-owned literal: "duty is not available from YAML on the EM sphere (WF-1 step 29 parses transient only; the API's DutyCycle is WF-15 step 17)". (Ref) asserts the reason equal to that module constant by exact string, at `thermal.duty`. The `initial` and `outputs.probes` refusals keep their live API capture: `study.py:1672` and `:2373` still carry U131's texts. Print the split for the weekly: the API accepts duty since U135, and YAML refuses it, because a sphere `thermal.duty` key is not commissioned.
- Change nothing else: every other gate, control and print is attempt 1's.

**Gates.** These are item 873's, with (IT) replaced as above:
- (K), (Ref), (UT), (Once);
- (A_t) at `A_BAND` and (R-A_t) `.4e` == `MAT12_A_T_PRINT`;
- (S2) at `S2_BAND` and (R-S2) == `MAT12_S2_PRINT`;
- (E) and CLI = API by `repr`. **CLI = API and (c-qs) run as gates for the first time here**: at attempt 1 they sat after the red (IT) assert, so their printed values were never asserted.

**Control (asserted, backed).** (c-qs): the steady θ misses `A_BAND`. U131 read 1.0108e+00 (`20261107T013321Z_WF-15-s15-W1.log:141`); attempt 1 printed the same value but did not reach the assert.

**Windows.** These are item 873's W0, W1, Wreg-c, Wreg-p, Wgen, Wval, Wapi and Wdoc, labelled `WF-1-s29r-…`. Each is bare and foreground, with Bash timeout 660000 ms and `pgrep -c python3` reading 0 before W1 and each Wreg. W1 measured 73 s at attempt 1. It now runs the same four EM solves plus two in-process thermal calls, so expect ≈ 80 s inside `-t 300`. Add one regression: W0 runs `tests/package/test_wf1_yaml_thermal_box_wall.py::test_wf1_s28_parse_and_refusals` and U139's module with its knob unset (`tests/package/test_wf1_yaml_thermal_box_dose.py`, landed by `7250188c`) beside the new module.

**Tier.** Standard, `-n 2`, complex. **Class:** chain step (solve). **Cost:** ≈ 8.2 measured / 30 predicted min.

**Traps already paid for.** All of item 873's, and also:
- **The apply.** A path checkout of `study_yaml.py` from the branch would revert U139's `thermal.dose`.
- **The duty text moved under U135.** Never compare the YAML's duty reason with a live API call.
- (IT-run) must be read from the same `run` object; a fresh `load_study` would be a second solve.

**Scope.** This is U134's target, unchanged. **Not claimed:**
- a band for the independent S read;
- a sphere duty key in YAML;
- degree 2 through YAML;
- ratification of this re-pose.

**Negative result.**
- **(IT-run) red:** the YAML route does not build U131's transient. Print the per-ingredient comparison, park on `attempt/WF-1-s29r-*`, and mark BLOCKED.
- **(IT-ind) A red, or (R-A_t) / (R-S2) red:** this is attempt 1's finding moving. Park, BLOCKED.
- **(Q-rep) < 1e-10 while (IT-ind) S > 1e-9:** the route and the twin build different q_h, and a printed-only S would hide that. Park, BLOCKED, and print every ingredient.
- **CLI = API or (c-qs) red** (first time asserted): park, BLOCKED.
- **Wreg red:** revert.
- **No `PREDICTED RED:` line.** The S independent read is printed, so it cannot go red.

**Status it can move.**
- `WF-1` gains step 29 (a record on its ✅ row), and §10 U134 is DONE.
- `MAT-12` joins `run_study_file`'s envelope.
- `thermal.transient` on a `sphere_in_box`, and this re-pose, go to the 11-09 weekly to ratify.
- §9 item 873 → DONE via 891.
- Item 887 is unblocked.

**Depends on:** nothing.

**Outcome (12:00 slot, 2026-11-08): BLOCKED — the brief's negative result "Wreg red"; parked on `attempt/WF-1-s29r-20261108T1821Z` (`0db69d87`).** The landing went as written: new files by path, attempt 1's three logs and their test-results rows by hand, `git apply -3` clean on `study_yaml.py` and three pure-adjacency conflicts in `validation_envelope.toml` (end-of-list appends beside U135's / U139's lines), resolved keeping both sides; `git diff HEAD` insertions only on the envelope, and on `study_yaml.py` only attempt 1's own two modified lines. The duty reason is the YAML-owned `_SPHERE_DUTY_REFUSAL`. W0 `20261108T181815Z_WF-1-s29r-W0.log`: 7 passed, 2 skipped per rank (2.1 s, with the s28 parse node and U139's module; the duty split printed). W1 `20261108T181834Z_WF-1-s29r-W1.log` (durable capture, `[capture] rc=0`, 70.6 s): 6 passed per rank. (UT) 19 496 / 19 496, (Once), (A_t) 3.4708293773e-02 (= U131 to ten digits), (S2) ratio 0.999950006 (4.9993987956e-05), (E) ≤ 4.812e-15, **(IT-run) A 5.603e-16, S 3.292e-16**, once before and after, **(IT-ind) A 7.845e-16** asserted, **S 3.958e-09 printed**, q_h owned max rel A 0.000e+00 / S 3.472e-09, **(Q-rep) 3.472e-09** with source-power rel 3.097e-11, so by the ruling's rule this is a **floor of the solve**, recorded. (c-qs) 1.0108e+00 > `A_BAND` and CLI = API (JSON round trip, fields by `repr`, CLI rc (0, 0), re-run with no mismatches) were asserted green for the first time. Wreg-c `20261108T181957Z_WF-1-s29r-Wreg-c.log`: 3 passed per rank (33 s). **Red: Wreg-p `20261108T182051Z_WF-1-s29r-Wreg-p.log`: 1 failed, 3 passed per rank.** `test_wf1_s14_parse_and_refusals` row 322 (`test_wf1_yaml_thermal_ball.py:322`) expects a sphere file carrying the ball's `dt_s`-keyed transient to be refused at `thermal.transient` (the landed unknown-key refusal). Now that the sphere owns `transient`, the refusal is at `thermal.transient.dt_s` (`_unknown` inside the block). It is still a refusal, but the field moved, so the edit reached a landed assertion. Neither brief foresaw it, and moving a landed expectation is not the slot's call. Wgen / Wval / Wapi / Wdoc were not run. **Unblock:** a review rules on row 322: either move its expected field to `thermal.transient.dt_s`, or have `_parse_thermal_sphere_transient` refuse a `dt_s` key at `thermal.transient`. Then the branch lands and W1 (rule (i) after any edit), Wreg-c, Wreg-p, Wgen, Wval, Wapi and Wdoc run.

**Review ruling (2026-11-08 18:00 daily review) on item 891's red — the landed row 322 stays unedited; the parser refuses `dt_s` at `thermal.transient`.** The only red is Wreg-p, `test_wf1_s14_parse_and_refusals` row 322 (`tests/package/test_wf1_yaml_thermal_ball.py:322`, `20261108T182051Z_WF-1-s29r-Wreg-p.log`): a sphere file carrying the ball's `dt_s`-keyed transient is refused at `thermal.transient.dt_s` where the landed row expects `thermal.transient`. Every gate of the re-pose itself was green (W1 `20261108T181834Z_WF-1-s29r-W1.log`: (IT-run) A 5.603e-16 / S 3.292e-16, (Q-rep) 3.472e-09 ⇒ floor of the solve by the 10:30 rule, (c-qs) and CLI = API asserted for the first time). The reading is not disputed (the slot's journal and the branch log agree; no pathologist needed).
- **Option (b), keeping every landed assertion unedited.** Item 873's own brief already fixed the field: "`dt_s` is not a key here: it is refused as unknown at `thermal.transient`" (this file, item 873 brief, the `_parse_thermal_sphere_transient` paragraph) and "(Ref) … a `dt_s` in place of `t_end_s`, at `thermal.transient`". The branch's module diverged from that text by letting `_unknown` raise inside the block. So `_parse_thermal_sphere_transient` gains, **before** its `_unknown(tsec, …)` call, one new block: a `dt_s` key in the transient mapping raises `StudyConfigError(source, "thermal.transient", <reason>)` with a YAML-owned reason constant naming `t_end_s` (dt = t_end_s / n_steps). The new module's (Ref) rows at `test_wf1_yaml_sphere_transient.py:199–204` (branch) expect `thermal.transient` — the brief's text, restored on unlanded code; row 322 is untouched.
- **Not option (a).** Moving row 322's expected field would edit a landed, pinned refusal to fit new code; the review declines it.
- **Family count:** U134 attempt 3 of 4 (873, 891, now 892). A red at 892 that is not the brief's pre-registered reading freezes the family at the next review.

### §9 item 892 — 2026-11-08 (§10 U134, chain U — 9a, WF-1 step 29 re-posed, second landing; chain step (solve), standard, complex)

*Written by the 2026-11-08 18:00 daily review under §9 rule (6) (a re-pose of a blocked step, keeping its target). Family U134, attempt 3 of 4. Base: `attempt/WF-1-s29r-20261108T1821Z` (`0db69d87`). Items 873 and 891's briefs (above) bind except where the ruling above and this brief say otherwise. Depends on nothing. Next: item 887 (U136).*

**Code onto `main`:**
1. **By path (new files):** `tests/package/test_wf1_yaml_sphere_transient.py`; the four `tests/package/wf1/` YAML files; the seven logs on the branch (attempt 1's three, 891's four). Append their seven test-results rows **by hand**, insertions only (`git diff --stat` on `docs/testing/test-results.md` shows insertions only).
2. **By three-way apply, never by path checkout:** write `git diff 1af95234 0db69d87 -- src/fem_em_solver/study_yaml.py tests/package/validation_envelope.toml` to a scratch file in the repo root, `git apply -3` it, delete the scratch file. Since `1af95234`, item 884 (`b05bdc1d`) edited both files (`_parse_thermal_rod`, `_run_rod_thermal`, one envelope line): a conflict that is not pure adjacency ⇒ stop BLOCKED. After the apply, U139's `thermal.dose`, U135's and U127′'s envelope lines and `_parse_thermal_rod` must still be present (`git diff HEAD` insertions only in both files, apart from 891's own two disclosed lines).

**The one edit this attempt makes** (the ruling above): in `_parse_thermal_sphere_transient`, a new block before `_unknown(tsec, …)` refusing a `dt_s` key at field `thermal.transient` with a module constant `_SPHERE_TRANSIENT_DT_REFUSAL` (text: "the sphere transient takes t_end_s (dt = t_end_s / n_steps), not dt_s"); the new module's (Ref) `dt_s` rows expect field `thermal.transient` and reason equal to that constant by exact string. Nothing else changes.

**Gates, control, prints:** item 891's, unchanged — (K), (Ref), (UT), (Once), (A_t) at `A_BAND` with (R-A_t), (S2) at `S2_BAND` with (R-S2), (E), (IT-run) at `L0_BAND` on A and S, (IT-ind) asserted on A and printed on S, (Q-rep) printed with its decision rule, CLI = API by `repr`; control (c-qs) asserted (U131 read 1.0108e+00, `20261107T013321Z_WF-15-s15-W1.log:141`; 891's W1 reproduced it).

**Windows** (labels `WF-1-s29r2-…`; each bare and foreground, Bash timeout 660000 ms, `pgrep -c python3` reading 0 before W1 and each Wreg): W0 (891's W0 set: the new module's parse nodes, `test_wf1_s28_parse_and_refusals`, U139's module; ≈ 2 s), W1 (891's W1, durable capture, `-t 300`; 70.6 s at 891), **Wreg-p first among the regressions** (891's Wreg-p set, which holds row 322; it is the one that went red), then Wreg-c, Wgen, Wval, Wapi, Wdoc as item 873 lists them. Add `tests/package/test_wf15_study_rod_heating.py`'s parse node(s) only (no solve) to W0, since 884 touched the same parser module.

**Tier.** Standard, `-n 2`, complex. **Class:** chain step (solve). **Cost:** ≈ 8.2 measured / 30 predicted min (891's windows summed ≈ 3 min of harness).

**Traps already paid for.** Items 873 and 891's, and also:
- the base of the apply is `1af95234` (891's own base), not `5054e5fd`;
- the new `dt_s` block goes **before** `_unknown` — after it, `_unknown` raises first and row 322 stays red;
- 884's `_parse_thermal_rod` lives in the same module; a path checkout would delete it.

**Scope.** U134's target, unchanged. Not claimed: a band for the independent S read; a sphere duty key; degree 2 through YAML; ratification.

**Negative result.** As 891's. In addition: **row 322 still red after the new block** ⇒ the block is not reached (order) or the field string differs — print both, park on `attempt/WF-1-s29r2-*`, BLOCKED; that is the family's fourth attempt and the next review freezes it. **Any other Wreg red ⇒ revert, BLOCKED.** No `PREDICTED RED:` line.

**Status it can move.** `WF-1` gains step 29 (a record on its ✅ row); §10 U134 DONE; §9 items 873 and 891 → DONE via 892; item 887 unblocked; `MAT-12` joins `run_study_file`'s envelope; the re-poses and `thermal.transient` on a `sphere_in_box` go to the 11-09 weekly to ratify.

**Depends on:** nothing.

**Outcome (19:30 slot, 2026-11-08): DONE — `WF-1` step 29, a record on the ✅ row; 873 and 891 DONE via 892.** Landing as written: new files and the seven branch logs by path, their seven test-results rows by hand (insertions only); `git apply -3` of `1af95234..0db69d87` clean on both files, `git diff HEAD` insertions only apart from 891's two disclosed `study_yaml.py` lines; `thermal.dose`, `_parse_thermal_rod` and the U135 / U127′ envelope lines intact. The one edit: `_SPHERE_TRANSIENT_DT_REFUSAL` raised at `thermal.transient` before `_unknown`; the new module's (Ref) `dt_s` row expects that field and constant by exact string (docstring line updated to match); row 322 untouched. W0 `20261109T013155Z_WF-1-s29r2-W0.log`: 9 passed, 2 skipped per rank (3.3 s; the s28 parse node, U139's module and both `test_wf15_study_rod_heating.py` no-mesh nodes passed; the two skips are the knob-gated s29 solve and s30 dose nodes), `dt_s` at `'thermal.transient'` with the constant's text. W1 `20261109T013212Z_WF-1-s29r2-W1.log` (durable capture, `[capture] rc=0`, 71.8 s): 6 passed per rank. (UT) 19 496 / 19 496, (Once), dt A `0x1.bdddddddddddep+1` / S `0x1.c891457ed6e76p-8`, (A_t) 3.4708293773e-02 (= U131), (S2) ratio 0.999950006 (4.9993987955e-05), (E) ≤ 4.411e-15, (c-qs) 1.0108e+00 > `A_BAND`, **(IT-run) A 6.724e-16 / S 5.486e-16**, once before and after; **(IT-ind) A 9.282e-12** asserted, **S 4.614e-09 printed** (q_h owned max rel A 5.404e-10 / S 3.792e-09); **(Q-rep) 3.814e-09**, source-power rel 2.634e-11 ⇒ floor of the solve, recorded; CLI = API (JSON round trip, fields by `repr`, CLI rc (0, 0), no mismatches). Wreg-p `20261109T013339Z_WF-1-s29r2-Wreg-p.log` (run first): 4 passed per rank, row 322 green. Wreg-c `20261109T013349Z_WF-1-s29r2-Wreg-c.log`: 3 passed per rank (34 s). Wgen `20261109T013433Z_WF-1-s29r2-Wgen.log` / Wval `20261109T013438Z_WF-1-s29r2-Wval.log` regenerated `api/TransientSpec.md`, `api/load_study.md`, `api/run_study_file.md` and `validation.md` (`MAT-12` now gating). Wapi `20261109T013446Z_WF-1-s29r2-Wapi.log`: 22 passed. Wdoc `20261109T013515Z_WF-1-s29r2-Wdoc.log`: 0 dead. `pgrep -c python3` read 0 before W1 and each Wreg. The (IT-ind) S and (Q-rep) values move between solves (891: 3.958e-09 / 3.472e-09), as the 10:30 ruling's reading predicts; the band for the independent S read stays the weekly's. Bands untouched.

## §7 row moved verbatim 2026-11-09 (rotate_plan_archive.py chunks)

 `femem run`: the MRI study CLI and YAML config | ✅ (item 256, U5): `femem run` builds Studies equal to the Python-built ones; every later YAML step, through `load_result:` on the F-clinical / F-human files (steps 25–26) and `tune.at_hz` (step 27, U111 — the public key ratified by the 11-02 weekly), is a record in the chunk file. U117's `load_result:` on the 3 T 64×64 (1.5 T YAML at `frequency_hz` 128 MHz) equals the API by `float.hex` (`WF-16` step 19, item 781, `20261103T134456Z_WF-16-s19-W1.log`). U118's `femem run --estimate` on the 3 T YAML (new key `field_strength_t`, for the 11-04 weekly to ratify) equals `Study.estimate` by `repr` (`WF-13` step 15, item 782, `20261103T135021Z_WF-13-s15-W1.log`). U119's `femem run --estimate` on the rod-in-clinical-box YAML (`implant:` on `f_clinical_birdcage`, entry 48's rod only — the refusal lift for the 11-04 weekly to ratify) equals `Study.estimate` by `repr` (`WF-13` step 16, item 783, `20261103T135850Z_WF-13-s16-W1.log`). U122's capacitor-sheet section on the implant page (`WF-17` step 28, item 798) reproduces item 802's 3×3 on the PEC rod at `IDENTITY_BAND` (6.377e-14, `20261104T031801Z_WF-17-s28-W1.log`), its degree-2 and copper-hole refusals raising the quoted text (`20261104T031743Z_WF-17-s28-W0.log`). U123 (`WF-17` step 29, item 806) puts `mri:33`'s (I-Py) twin on `solve_sweep(fields_at=)`: (I-Py) 1.232e-12 ≤ `L0_BAND` on the closing window (`20261104T032240Z_EX-138.log`, Status 0, 114 s), and no `_solve_port_fields` is left in `examples/`. U125 (`WF-1` step 28, item 830): `thermal.initial` and `geometry.wall` (P48b's 1D-reduced walled box; public keys for the 11-07 weekly) reproduce P63's R1 at P5 at `IDENTITY_BAND` (rel 0.000e+00 on 119 808 cells) and the R0 file misses it by W1b's drift (|Δ| 3.345e-08 ≤ 1e-7) (`20261105T210638Z_WF-1-s28-W1-r2.log`). U139 (`WF-1` step 30, item 878): `thermal.dose` / `DoseReadout` (public keys for the 11-09 weekly) on U125's walled files reproduces P68's (Q-wf) `cem43` gate at `DOSE_BAR` (≤ 5.715e-16) and (Win-wf) at `AVG_BAR`, the history equal to a direct R1 at `IDENTITY_BAND`; the R0 dose misses R1's by 2.955e-01 (`20261108T150754Z_WF-1-s30-W1.log`). U134 (`WF-1` step 29, item 892 via 873 / 891): `thermal.transient` on a `sphere_in_box` (public key and re-poses for the 11-09 weekly to ratify) reproduces U131's (A_t) 3.4708293773e-02 at `A_BAND` and (S2) 4.9994e-05 at `S2_BAND`, (IT-run) ≤ 6.724e-16 at `L0_BAND`, CLI = API by `repr`; the independent S read is printed (4.614e-09, a solve floor by (Q-rep) 3.814e-09, no band) (`20261109T013212Z_WF-1-s29r2-W1.log`). *History: `docs/planning/chunks/WF-1.md`.* 
