# ARCH-1 — The MRI-safety layer becomes solver-agnostic

 **The MRI-safety layer becomes solver-agnostic** — SAR averaging, B₁⁺ maps, loss partition, the thermal solve, the benchmark cases and their gates take fields and S-parameters from *any* field solver (ours, Palace's output, an HFSS field export) through one interface; the FEniCSx engine stays the in-house backend and the multiphysics home (operator directive 2026-09-20). **EMBARGOED: not before the weekly of 2026-10-17** — four weeks out by the operator's instruction, so that `PORT-21`, `OPS-61` (Palace) and the first human-scale physics land first and the interface is designed against two real backends, not one | ⬜ **Not queueable by a daily review** (outside §9 rule (4): it is architecture, not an enumerated chain step) **and not to be scheduled by the weeklies of 09-26, 10-03 or 10-10** — those may only annotate this row with what they learn. **Why:** the project's replaceable part is the curl-curl solve; its hard-to-replace parts are the verification culture, the MRI-safety layer and the multiphysics reach. Today the safety layer is written against DolfinX function spaces, so it can only judge our own fields; decoupled, it is an evidence toolkit that works with the field solver a user already has — the wider audience arrives through "works with Palace / HFSS output", not through replacing the engine — and it is also what lets the engine be dropped later without losing anything else, should the target narrow. **Shape (for the 10-17 weekly to write at execution grain):** (1) a `FieldSolution` interface — an unstructured mesh with cell tags plus `E` (complex, cell- or edge-based) and the material map, read from XDMF/VTU or from a `dolfinx` function directly; and a `PortSolution` — a Touchstone file plus port metadata; (2) `post/sar.py`, `post/b1.py`, the loss partition and `thermal/` consume only those; (3) two backends gated *identically*: the DolfinX one reproduces every existing `MAT-4` / `WF-6` / `TH-14` record to the digit (the refactor moves no number), the Palace one reproduces `OPS-61`'s S and the phantom power on the `ANS-4` fixture; (4) the benchmark cases' runnable halves dispatch on the backend. **Gate of record:** digit-identical reproduction of every imported band and record on the DolfinX backend before any second backend is wired — a refactor that moves a record is a defect. **Not in scope:** a GUI, a new solver, changing any physics. **Prerequisites the embargo waits on:** `OPS-61` step 1 (a second backend's fields to design against), `PORT-21` (so the port interface is not frozen around a known offset), one F-human physics step (so the interface sees a human-scale field once). 

### §9 item 517 — 2026-10-18 (§10 U16, `ARCH-1` step 1; chain step, solve)

*Drafted by a read-only general-purpose agent for the 2026-10-18 03:00 daily review; the review reviews it. U16 is drafted as two items, **517 (U16a, the Study-side backend, the case writer and the read path, no Palace run)** and **517b (U16b, the Palace runs and the S read-back, no code edit)**: the whole step does not fit one 1-hour slot (≈ 25 min of `src/` + test code, four re-gate / build windows ≈ 10 min, two Palace windows ≈ 5 min, fixed + commit ≈ 21 min). This is a split on price and on a hard architectural fact (below), not a re-target. Four rulings are asked for at the end.*

**Licence.** §10 chain U, U16 (PROJECT_PLAN.md:4864–4870), weekly-enumerated (§9 rule (4)), quoted, not reworded: "`Study.solve(backend="palace")` and `femem run --backend palace` on the F-small four-port exact fixture: the Study writes the Palace config and mesh through `OPS-61`'s writer, runs the `fem-em-palace` service, reads S back; target: equality with `OPS-61` F26's recorded Palace read at ≤ 1e-9 (same input, same answer) and the tie to our degree-2 record at the recorded (T5) ≤ 1.9e-10 per class; the DolfinX path bitwise unchanged." The 2026-10-17 weekly lifted `ARCH-1`'s embargo and named U16 its first step, glyph 🟡 when U16 lands (PROJECT_PLAN.md:4705–4709); the row's precondition (digit-identical DolfinX reproduction) is met by `WF-16` step 1 and `WF-14`'s snapshot (same lines). The Palace directive (§5.4 (i)–(v), PROJECT_PLAN.md:1176–1197), the sheets-and-κ ruling (exact rectangles, κ off), the 2026-10-13 03:00 rulings (1) and (4) (`docs/planning/chunks/OPS-61.md:1208–1213`) and the 2026-10-15 03:00 ruling on item 476 (`docs/planning/chunks/ANS-4.md:651–653`) bind. **No band, tolerance, Done-when or validation target is set or moved.**

**Finding that shapes the design: the solver container cannot run Palace.** The `fem-em-solver` service mounts only `..:/workspace` and the `fem-em-cache` volume (`docker/docker-compose.yml:8–10`); there is no docker socket inside it, and Palace lives only in its own image (`fem-em-palace:0.18.1`, `docker/palace.Dockerfile`), which has no DolfinX and no `fem_em_solver`. Every `OPS-61` Palace window was a **separate host-side harness call**, `docker compose -f docker/docker-compose.yml --profile palace run --rm -T fem-em-palace bash -lc 'cd /workspace/logs/ops61 && timeout -k 30 590 palace -np 8 <config>.json'` (e.g. `20261003T110249Z_OPS-61-item194-W4-palace-o2.log:12`), between two solver-container windows (write, then read). Mounting the socket is a mount change, which the §9 compose constraint forbids (PROJECT_PLAN.md:3946–3965: a blocked finding for the operator); putting Palace into the solver image is a Dockerfile change and also the operator's. **So `Study.solve(backend="palace")` cannot literally "run the `fem-em-palace` service".** The honest design, drafted here for the review to rule on (Ruling R2), is a **three-stage backend**:
1. **Write.** `Study.solve(backend="palace", case_dir=D)` builds the Study's mesh, writes `D/case.json` + `D/mesh.msh` + `D/case.sha256` (a manifest of both files), or verifies byte-equality if they already exist. If `D/postpro/port-S.csv` is absent it raises `PalaceRunRequired`, whose `.argv == ["palace", "-np", "8", "case.json"]` and `.cwd == D`. `femem run --backend palace` prints the command and returns **3** (proposed new exit code, "external run pending"; R4).
2. **Run.** A separate harness window runs that command in `fem-em-palace` (host-side, exactly the `OPS-61` shape).
3. **Read.** The same call again rebuilds and verifies the case byte for byte against the manifest, then reads `D/postpro/port-S.csv` through the lifted `OPS-61` reader and returns a result carrying `s_matrix`. `femem run --backend palace` writes `summary.json` and returns 0.

The chain of custody (manifest checked at read time) is what makes "same input, same answer" a statement about the Study's own input, not about whatever sits in the directory.

**Defects in the step text (for the review).**
- **"`OPS-61` F26's recorded Palace read" names no Palace read.** F26 is `ANS-6`'s copper h-rung, `xl` entry 31 (`docs/planning/chunks/ANS-6.md:509`, `:620`); the Palace read on that rung is F37 (§9 item 454: h = 0.0075 *hole* mesh, `Boundaries.Conductivity` on 401, four single-excitation windows of 142.9–144.9 s, `docs/planning/chunks/OPS-61.md:2076`). That is not the Study's F-small default (`FourPortBirdcageGeometry.f_small()` is the ×1 strip route, σ = 800 S/m volume conductor, `src/fem_em_solver/study.py:225–269`), and writing it would need a conductivity-boundary writer, `as_hole` and a resolution knob — three features beyond U16. **Drafted on the fixture that is the Study's `f_small()` verbatim and that `OPS-61` did read:** the ×1 exact strip, `logs/ops61/ops61_ans4_x1_strip.msh` (117 206 tets, F10a), Palace v0.18.1 order 2 at 64 MHz, `f05x_o2.json` → `postpro_f05x_o2/` (F17, §9 item 194, `20261003T110249Z_OPS-61-item194-W4-palace-o2.log`; read `20261003T110728Z_OPS-61-item194-W5-read.log:49–56`). Ruling R1.
- **"≤ 1.9e-10 per class" is the figure the 2026-10-15 03:00 review already ruled false** (item 476: the bought maximum is 2.542e-10; the public ceiling is **2.6e-10**, `docs/planning/chunks/ANS-4.md:613`, `:653`; `20261015T112551Z_ANS-4-item476-ties.log`). On the drafted rung the recorded tie is **2.464e-10 / 9.668e-11 / 9.883e-11** (self / adjacent / opposite, `20261004T214342Z_OPS-61-item229-W4-read.log:41`), so "≤ 1.9e-10" would fail the self class on Palace's own recorded output. Also, no cross-solver band constant exists and ties have been printed, never asserted (2026-10-13 ruling (1), `OPS-61.md:1210`). Ruling R3.
- **"through `OPS-61`'s writer" cannot be an import.** The writer is `export()` in a 5 410-line probe under `scripts/` that imports `tests/` (`scripts/probes/ops61_step1a_palace_ans4.py:106–136`), and `src/` may not import either (`WF-13`'s rule, as `WF-16` step 1 applied it, `docs/planning/chunks/WF-16.md`, item 244). U16a **lifts the writer and reader verbatim** into `src/` — the Gmsh 2.2 writer `:213–302`, the config dict `:315–361`, the reader `_read_port_s` `:400–425` — and pins the copies with equality anchors against the probe's own files and functions (the `WF-16` step 1 precedent). The probe is not edited.
- **§5.4 still says "nothing in `src/` is built on Palace"** (PROJECT_PLAN.md:1180–1181), and both Palace files say "Never the solver — nothing in src/ uses it" (`docker/docker-compose.yml:107`, `docker/palace.Dockerfile:2`). U8 (operator directive 2026-10-05, Palace as an opt-in backend behind the `Study`) supersedes that sentence in part, but it was never amended. U16 adds a Palace *writer and reader* to `src/`: it never imports or links Palace, and the default path never touches it. The plan sentence is the review's / weekly's to amend (Ruling R4). The compose comment is operator-owned (only `environment:` keys may be edited); flag it, do not edit it.
- **The `ARCH-1` row still reads "EMBARGOED … Not queueable by a daily review"** (this file `:3` and the §7 row). The review amends it when it queues 517.
- **No Touchstone and no `z_matrix` on the Palace path.** Palace's `z_matrix_ohm` is a terminated transimpedance (2026-09-29 ruling, never read), and `PortStudyResult.write_touchstone` needs an `SParameterSweepResult` (`study.py:597–613`). U16's summary carries `s_matrix` only; a Palace Touchstone is a follow-on, not U16.

**`ARCH-1` step 1a / U16a: the `palace` backend's write and read stages on the F-small ×1 exact strip, the written case proven equal to `OPS-61`'s, the read path proven on `OPS-61`'s recorded output, and the DolfinX path re-gated.**
- **Executor:** implementer. Additive `src/` (a new private subpackage plus one dispatch branch each in `Study.solve`, `run_study_file` and `cli.main`), disclosed under ruling (c). One new test module. No probe, gate module, record, registry or generated page is edited.
- **Tier:**
  - W0 smoke (`-n 1`, no mesh);
  - W1 heavy by ceiling (`-n 1`, mesh-only builds, no solve; ≈ 6 builds at ≈ 26–31 s each, F10a W3 31 s / mesh 26.4 s, `20261001T220353Z_OPS-61-F10a-W3-export.log`; predicted 200–330 s);
  - W2 standard → heavy by ceiling (`-n 2`; 138 s measured, `20261005T202623Z_WF-16-W1-deg1.log`);
  - W3 standard → heavy by ceiling (`-n 2`; 164 s measured, `20261007T051729Z_WF-1-s3-W3-wf1-w2.log`);
  - W4 smoke (`-n 1`; 17 s measured, `20261015T214657Z_WF-17-s6-Wapi.log`).
- **Class:** chain step (solve; the DolfinX re-gates solve). **Cost:** ≈ 18 measured / ≈ 55 predicted slot-minutes: 15 fixed + ≈ 25 code + windows ≈ 10 (10 + 270 + 138 + 164 + 17 s) + commit and journal ≈ 6. Tight. **If W1 is not green by minute 40, park** (rule (d)).
- **Depends on:**
  - nothing in §9;
  - on disk (gitignored, all present 2026-10-18, checked): `logs/ops61/ops61_ans4_x1_strip.msh` (Oct 1), `logs/ops61/ops61_ans4_x1.msh` (Sep 29, the ragged control), `logs/ops61/f05x.json`, `logs/ops61/f05x_o2.json`, `logs/ops61/postpro_f05x/port-S.csv`, `logs/ops61/postpro_f05x_o2/port-S.csv`, `logs/ops61/postpro_f05x_o2_f10/port-S.csv`.
  - If any is absent: the test that needs it skips with `BLOCKED: <file> absent`, and the item lands as BLOCKED with that line. **Never re-export or overwrite an `OPS-61` file** (their traps: `OPS-61.md:20`, `:98`, `:244`).
- **Family:** `ARCH-1` step 1 (U16); attempt 1 of 4.

**What the code is today.**
- `Study.solve(self, comm=MPI.COMM_WORLD)` (`src/fem_em_solver/study.py:170–197`) has no backend: a port geometry dispatches to `_solve_ports` (`:173–174` → `_study_solve_ports`, `:616`, bound at `:687`). No `palace` / `backend` string exists anywhere in `src/`.
- `FourPortBirdcageGeometry.build` (`study.py:271–392`) returns the mesh, cell tags, the **narrowed** sheet tags and per-sheet records; the full-sheet tags from `_sheet_areas` (`:317–325`) are reassigned and not returned. The writer needs F_i ∖ N_i (attributes 30 + i) to reproduce `OPS-61`'s file, so it recomputes `_sheet_areas` from `ports/_sheet_geometry.py:93` with the geometry's `sheet_iface` / `port_lower` / `port_upper`. `build` is not edited.
- `run_study_file` (`src/fem_em_solver/study_yaml.py:1354`) calls `study.solve(comm)` at `:1383`; the four-port branch writes the Touchstone and `s_matrix` / `z_matrix` (`:1412–1426`).
- `cli.main` (`src/fem_em_solver/cli.py:120–181`) maps `StudyConfigError` to rc 2 (`:178–181`); the `run` parser is `:35–48`.
- The probe's writer runs serial only (`assert comm.size == 1`, `ops61_step1a_palace_ans4.py:139`).

**Code (additive; disclosed under ruling (c)).**
- **New `src/fem_em_solver/backends/__init__.py` (empty) and `src/fem_em_solver/backends/palace.py`** (found by `[tool.setuptools.packages.find] where = ["src"]`, `pyproject.toml:73–74`). Nothing is added to `fem_em_solver.__all__`, so `WF-14`'s snapshot is unchanged.
  - `PALACE_NARROW_BASE = 20`, `PALACE_REMAINDER_BASE = 30`, `PALACE_EXTERIOR_ATTR = 1` (the probe's `:86–88`).
  - `write_gmsh22(path, msh, cell_tags, tri_blocks)`: the probe's `:213–302` verbatim, in operation order (positive-orientation swap, `.17g` nodes, triangles before tets).
  - `palace_config(*, mesh_name, output, materials, ports, order, freq_ghz)`: the probe's `config()` (`:336–361`) verbatim, with `Order` and `Freq` as arguments. `freq_ghz = study.frequency_hz / 1e9`, asserted to `json.dumps` as `0.064` / `0.01`.
  - `read_port_s(csv_path)`: `_read_port_s` (`:400–425`) verbatim, returning the matrix; plus the frequency column. Raises if that column differs from the config's `Freq`.
  - `write_case(study, built, case_dir)`: refuses `comm.size != 1` (`ValueError`); writes `case.json` (`Output "postpro"`, `Model.Mesh "mesh.msh"`) and `mesh.msh`, then `case.sha256`. If the files exist it compares bytes and raises `FileExistsError` naming the differing file; it never overwrites.
  - `PalaceRunRequired(RuntimeError)` with `.argv`, `.cwd`; `PalaceStudyResult` (dataclass: `s_matrix`, `frequency_hz`, `order`, `cells`, `case_dir`, `sha256`, `palace_total_s` read from `postpro/palace.json` `ElapsedTime.Durations.Total`, `backend = "palace"`).
  - Order is `study.degree` (1 → `Order 1`, 2 → `Order 2`); `R` is `LumpedPortDrive.port_impedance_ohm`, and `z0_ohm != port_impedance_ohm` is refused (Palace's S is referenced to `R`); every sheet must have `drive_direction == (0, 0, 1)` → `"+Z"`.
  - **Refused before any mesh** (`TypeError` / `ValueError`, naming the reason): `SphereInBoxGeometry`, `ThermalBallGeometry`, `FHumanBirdcageGeometry`, `as_hole=True`, a `conductor`, `CapacitorSheets`, `TwoPortQuadratureDrive`, `solve_sweep`. Materials are the geometry's conductor σ (`conductor_sigma_s_per_m`, 800) and the Study's phantom σ / εᵣ — as the probe writes them (`:317–334`).
- **`study.py`: `Study.solve(self, comm=MPI.COMM_WORLD, *, backend="dolfinx", case_dir=None, palace_ranks=8)`.** The first new line is `if backend != "dolfinx": return _solve_backend(self, backend, comm, case_dir=case_dir, palace_ranks=palace_ranks)`; an unknown backend raises `ValueError` naming `("dolfinx", "palace")`. **The `dolfinx` body is not touched, and neither is `Study`'s class docstring or `__init__` signature** — the generated API page renders those (`scripts/docs/generate_api_reference.py:74–116`), so editing them makes `test_wf17_api_reference.py` (G) red.
- **`study_yaml.py`: `run_study_file(path, *, output_dir=None, comm=..., backend="dolfinx")`.** `backend == "palace"` branches before `:1383`: it refuses every non-`four_port_birdcage` family, `outputs.sar` / `b1_plus` / `tune`, a conductor and `capacitor_sheets` (`StudyConfigError`, field `--backend`), calls `study.solve(comm, backend="palace", case_dir=out / "palace_case")`, and writes `summary.json` with the existing keys plus `"backend": "palace"`, `"palace": {"case_dir", "order", "sha256", "palace_total_s"}`, `"s_matrix"` — no `z_matrix`, no Touchstone. **The default path writes no new key.**
- **`cli.py`:** `run.add_argument("--backend", choices=["dolfinx", "palace"], default="dolfinx")`, threaded into `run_study_file`. `PalaceRunRequired` is caught: the command is printed on rank 0 as `femem: palace run required: (cd <cwd> && palace -np 8 case.json)`, and the call returns **3**. `--backend palace` with `--ranks` or `--estimate` returns 2 before any launch.
- **New `tests/package/test_arch1_palace_backend.py`.** Docstring paragraph 1 names `ARCH-1` (the validation page's (I) resolves §7 IDs from it). `@complex_only` on the mesh nodes. Constants: `U16_DIR = ROOT / "logs" / "u16"`, `D2_CASE = U16_DIR / "f_small_d2_64"`, `CLI_OUT = U16_DIR / "f_small_d1_64_cli"`, `READ_ENV = "FEM_EM_U16_PALACE_READ"`. The probe is loaded by `importlib` (the `_load` precedent, `tests/validation/test_wf5_q_two_ways.py:71–78`).
  - **(T0) `test_u16_backend_refusals`** (no mesh): every refusal above, each raising before a mesh (wall time < 2 s, printed), and no file written under `tmp_path`; `backend="hfss"` → `ValueError`; `Study.solve` at `comm.size == 2` with `backend="palace"` → `ValueError` (asserted in W2's width only if cheap; otherwise by monkeypatching `comm`).
  - **(T1) `test_u16_case_equals_ops61`** (`-n 1`): the Python Study (`FourPortBirdcageGeometry.f_small()`, phantom 0.5 / 78, `LumpedPortDrive()`, 64 MHz, **degree 2**) solved with `backend="palace", case_dir=D2_CASE` raises `PalaceRunRequired` (D2_CASE has no `postpro/`); `cli.main(["run", "tests/package/wf1/f_small_four_port.yaml", "--output-dir", str(CLI_OUT), "--backend", "palace"])` (degree 1) returns 3. Then (W-cfg), (W-msh) and controls (c1) / (c2) on both written cases.
  - **(T2) `test_u16_read_path_on_recorded_output`** (`-n 1`): copies `D2_CASE`'s three files into `tmp_path` and **copies** (never moves) `logs/ops61/postpro_f05x_o2/` to `tmp_path/postpro/`; `Study.solve(backend="palace", case_dir=tmp_path)` returns a `PalaceStudyResult`. Asserts (Rd) and (P-rec); then (c3) and (c4).
  - **(T3) `test_u16_cli_read_path_on_recorded_output`** (`-n 1`): the same with the degree-1 CLI case and `postpro_f05x/`; `cli.main([... "--backend", "palace"])` returns 0; `summary["backend"] == "palace"`, `"z_matrix" not in summary`, `"touchstone" not in summary`, and `summary["s_matrix"]` equals `read_port_s` on the copied CSV entry for entry (JSON pairs of floats round-trip exactly).
  - **(T4) `test_u16_palace_read_back`** and **(T5) `test_u16_cli_palace_read_back`**: item 517b's nodes, skipped unless `READ_ENV=1` (they are committed here so that 517b edits nothing; rule (i)).

**Anchors (asserted; bands imported, never restated).**
- **(W-cfg)** `json.loads(D2_CASE/"case.json")` against `logs/ops61/f05x_o2.json`: the set of differing key paths is exactly `{"Problem.Output", "Model.Mesh"}`; every other key and value is equal (`Materials` including the air list `[2, 101…104, 201…204]`, `PEC [1]`, four `LumpedPort`s on `[21]…[24]`, `+Z`, `R 50.0`, `Order 2`, `Freq [0.064]`, `Linear` block). Likewise the CLI's degree-1 case against `logs/ops61/f05x.json`.
- **(W-msh)** each written `mesh.msh` against `logs/ops61/ops61_ans4_x1_strip.msh`, parsed: node count, the node coordinate array in file order (exact `==` on the parsed floats), every element's type / attribute / node list in file order, 117 206 tets (= `ANS4_EXACT_STRIP_CELL_RECORD`, `tests/validation/test_ops61_exact_strip_ports.py:64`), per-attribute triangle counts (attributes 1, 21–24, 31–34; 40 triangles on each N_i, F10a W3). **Printed, predicted equal:** the sha256 of both files (the same code at the same width on the same generator kwargs, `WF-16` (K)).
- **(Rd)** `backends.palace.read_port_s(p)` against the probe's `_read_port_s(p)[0]` on `logs/ops61/postpro_f05x_o2/port-S.csv` and `postpro_f05x/port-S.csv`: `np.array_equal` (the same element-wise arithmetic in the same order; no solve or inverse, so the 2026-10-15 18:00 "bit for bit" note, `WF-16.md:1751–1753`, does not apply).
- **(P-rec)** on the planted read: reciprocity ≤ `RECIPROCITY_BAND` (`tests/validation/test_port_lumped_sheet_sweep.py:86`), σ_max and every column power ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE` (`tests/validation/test_port_birdcage_four_port.py:129`), every S-class spread ≤ `ADJACENT_SPREAD_BAND` (`tests/validation/test_port_birdcage_lumped_column.py:121`), pooled ≥ `POOLED_SEPARATION_FLOOR` × the worst (`test_port_birdcage_four_port.py:173`), classes from `_circulant_classes` / `_class_spread` (`:193`, `:176`). σ_max printed beside item 194's 0.999819887373 (`20261003T110728Z_OPS-61-item194-W5-read.log:50`).
- **(Bit) the DolfinX path:** W2 reproduces `WF-16` step 1's degree-1 set (19 passed; (R) 9.557e-11, σ_max 0.999760093664, spreads 0.0747 / 0.0691 / 0.0897 %, `20261005T202623Z_WF-16-W1-deg1.log`); W3 reproduces `WF-1`'s F-small CLI set (13 passed; (R1) 9.557e-11, `20261007T051729Z_WF-1-s3-W3-wf1-w2.log`), with no new key in its `summary.json`; W4 passes the envelope snapshot, the API-page byte check and the validation page unchanged (22 passed, `20261015T214657Z_WF-17-s6-Wapi.log`).

**Negative controls.**
- **Asserted, by construction:**
  - **(c1)** the CLI's degree-1 `case.json` against `f05x_o2.json` differs in exactly `{"Problem.Output", "Model.Mesh", "Solver.Order"}` — the order mapping reaches the config;
  - **(c2)** the written mesh against the ragged `logs/ops61/ops61_ans4_x1.msh` is unequal: 117 206 vs 116 085 tets (`20261001T215942Z_OPS-61-F10a-W1.log:3581`);
  - **(c3)** a case in `tmp_path` holding `postpro_f05x_o2_f10/port-S.csv` (10 MHz) is refused by the frequency check, naming 0.01 vs 0.064;
  - **(c4)** a one-byte edit to the copied `case.json` makes the next `solve(backend="palace")` raise `FileExistsError` naming `case.json` — the manifest check can fail.
- **Predicted, printed:** W1's summed RSS ≤ 3 GiB.

**Windows.** Each runs in the foreground with Bash `timeout: 660000`, with `-s`, in the complex build, never with `-k`. **The harness call is the whole Bash command.**
- **W0** (smoke, `-n 1`, T0 only):
  `scripts/testing/run_and_log.sh ARCH-1-U16a-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 60 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_arch1_palace_backend.py::test_u16_backend_refusals -v -s --tb=short'"`
- **W1** (`-n 1`, the whole module, knob unset; T4 / T5 skip; predicted 200–330 s):
  `scripts/testing/run_and_log.sh ARCH-1-U16a-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 590 mpiexec -n 1 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_arch1_palace_backend.py -v -s --tb=short > /workspace/logs/arch1-u16a-w1-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/arch1-u16a-w1-raw.log; cat /workspace/logs/arch1-u16a-w1-raw.log; exit \$rc'"`
- **Before W2 and W3:** `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0.
- **W2** (the `WF-16` step 1 re-gate; that log's `:12` verbatim, relabelled; 138 s measured):
  `scripts/testing/run_and_log.sh ARCH-1-U16a-W2-regate-wf16 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_WF16_DEGREE=1 timeout -k 30 400 mpiexec -n 2 python3 -m pytest tests/environment tests/package/test_wf16_study_four_port.py -s -v --tb=short'"`
- **W3** (the `WF-1` F-small CLI re-gate; that log's `:12` verbatim, relabelled; 164 s measured):
  `scripts/testing/run_and_log.sh ARCH-1-U16a-W3-regate-wf1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 400 mpiexec -n 2 python3 -m pytest tests/environment tests/package/test_wf1_study_yaml.py::test_wf1_cli_f_small_reproduces_wf16 -s -v --tb=short'"`
- **W4** (envelope, API pages, validation page — unchanged; 17 s measured):
  `scripts/testing/run_and_log.sh ARCH-1-U16a-W4-wapi "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace timeout -k 30 170 python3 -m pytest tests/package/test_wf14_validation_envelope.py tests/package/test_wf17_api_reference.py tests/package/test_wf14_validation_page.py -v -s --tb=short'"`
- No edit to any `src/` file or the test module after W0 starts (rule (i)); any edit means W0 to W4 again.

**Traps already paid for.**
- **Never write into `logs/ops61/`, and never plant output into `logs/u16/`.** T2 / T3 copy recorded `postpro_*` into `tmp_path` only. `logs/u16/f_small_d2_64/` and `logs/u16/f_small_d1_64_cli/palace_case/` must reach item 517b holding `case.json`, `mesh.msh`, `case.sha256` and **no `postpro/`**; Palace writes that directory itself.
- **`src/` imports neither `tests/` nor `scripts/`.** The lifted code is a copy; (W-cfg), (W-msh) and (Rd) pin it to the originals.
- **Freq is in GHz in the config.** `64e6 / 1e9` must serialise as `0.064`; (W-cfg) catches a drift. Never round it by hand.
- **Writer width.** The Gmsh writer is serial (probe `:139`); a `-n 2` call must raise before meshing, not write a rank-local file.
- **The API page renders `Study`'s class docstring and `__init__` signature.** Touch neither; `backend` lives on `solve`. W4 proves it.
- **The default `summary.json` gains no key** (`"backend"` only on the Palace path).
- **Harness logs are gzipped after rotation** (`docs/testing/logs/*.log.gz`); cited `.log` names are read with `zcat`.
- **A killed window:** `pgrep -c python3` first; at 0, run the rubric-4 cache sweep as `ARCH-1-U16a-cache-sweep`.
- Never `-k`; no pipe; never delete a log; two heavy windows never run at once; durable capture ends in `; exit \$rc`.
- Compare S classes only; never read Palace's `z_matrix_ohm`; never open `docs/private/`, `aed_results/` or `COMPARISON_private.md`.

**Scope.** The backend's dispatch, the written case equal to `OPS-61`'s, the read path on `OPS-61`'s recorded output, the CLI's two stages, the DolfinX path re-gated. **Not claimed:** any Palace run, the U16 equality or tie (517b), fields from Palace (SAR, B₁⁺, thermal — `ARCH-1`'s `FieldSolution` interface), sweeps, copper / Leontovich / hole routes, F-human, tuning, any accuracy or AED figure.

**Negative result.** Report and stop. Park on `attempt/ARCH-1-U16a-*`, mark item 517 BLOCKED in the same commit (rule (d)), and withdraw 517b.
- **(W-cfg) red:** the lifted config differs from the probe's. Fix once inside the new files; otherwise park.
- **(W-msh) red with 117 206 tets:** the Study's mesh is not `OPS-61`'s file (an environment move since 2026-10-01). Open a known-issues entry, BLOCKED; 517b cannot claim "same input". This is a `record-reconciler` question, never a re-export of `OPS-61`'s mesh.
- **(Rd) red:** a copy error. Fix once; otherwise park.
- **W2 / W3 red:** the dispatch moved the default path. Revert the `src/` diff on the branch; BLOCKED.
- **W4 red:** a public-surface or docstring leak. Fix in-letter (remove the leak) and re-run W0–W4.
- Status 124 is not a result; park with the measured time.

There is no `PREDICTED RED:` line, so a red anywhere stops the slot.

**Status it can move.** No glyph: the `ARCH-1` row gains "**Step 1a record** (§9 item 517, U16a): `Study.solve(backend="palace")` / `femem run --backend palace` write a Palace case equal to `OPS-61`'s ×1 exact-strip `f05x_o2` / `f05x` inputs (config keys, mesh content), read S through the lifted reader on `OPS-61`'s recorded output; DolfinX path re-gated (<W2 log>, <W3 log>)". §10 U16 records "a DONE, b open".

**Execution order.** Code → W0 → W1 → `pgrep` → W2 → `pgrep` → W3 → W4 → commit (the `backends/` subpackage, the three dispatch edits, the test module, harness logs, the §7 sentence) → journal.

**Spawn line (if delegated):** "Use the implementer agent to execute `ARCH-1` §9 item 517 (U16a) per `docs/planning/chunks/ARCH-1.md` (read the brief and any review ruling beneath it first). The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run W0–W4 in order, in the foreground, with `timeout: 660000`, pytest with `-s`, container `timeout -k 30`. Before W2 and W3 run `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` and require 0. Never write into `logs/ops61/`; never put a `postpro/` into `logs/u16/`. Do not edit `Study`'s class docstring or `__init__`, the probe, any gate module, the registry or any generated page. You are the executor: do not spawn agents and never background a harness window. If the docker socket is denied, stop and report the exact command string. On a red, park on `attempt/ARCH-1-U16a-*`, mark item 517 BLOCKED in the same commit, and report."

### §9 item 517b — 2026-10-18 (§10 U16, `ARCH-1` step 1b; chain step, solve; serial on item 517)

*Drafted with item 517; the review reviews both. No code edit.*

**`ARCH-1` step 1b / U16b: Palace runs on the two cases the `Study` and `femem` wrote, and the `Study` reads S back — equal to `OPS-61`'s recorded read on the same input, tied to our degree-2 record.**
- **Executor:** implementer. **No edit to any file** except the commit's record lines.
- **Tier:**
  - W1 heavy (Palace `-np 8`, order 2, 64 MHz; 268 s measured on the same config, `20261003T110249Z_OPS-61-item194-W4-palace-o2.log` — `Total` 265.6 s `:351`, 4 × 37 GMRES iterations `:164–320`, peak 6.1 G `:368`; predicted 240–400 s on today's box);
  - W2 smoke → standard (Palace `-np 8`, order 1; 26 s measured, `20261002T021838Z_OPS-61-F10b-W2-palace-f05x.log`, Palace 23.35 s, 1.8 G);
  - W3 standard (`-n 1`, two mesh rebuilds and reads; predicted 70–120 s).
- **Class:** chain step (solve). **Cost:** ≈ 18 measured / ≈ 30 predicted slot-minutes: 15 fixed + windows ≈ 7 (≈ 300 + 26 + 100 s) + pre-checks ≈ 2 + commit and journal ≈ 6. **Start only before minute 30; W1 not started by minute 35 ⇒ do not start it** (the slot rule; the Palace window may not cross the timebox).
- **Depends on:**
  - **item 517 DONE on `main`** — skip otherwise;
  - on disk: `logs/u16/f_small_d2_64/{case.json, mesh.msh, case.sha256}` and `logs/u16/f_small_d1_64_cli/palace_case/{…}` with **no `postpro/`**; and `logs/ops61/postpro_f05x_o2/port-S.csv`, `logs/ops61/postpro_f05x/port-S.csv`.
  - If a `logs/u16` case is absent, re-run item 517's W1 verbatim relabelled `ARCH-1-U16b-W0-rewrite` first (it writes them). If a `postpro/` is already there, stop and report: never overwrite a Palace output. If an `OPS-61` CSV is absent: BLOCKED ("comparand absent"); never re-run `OPS-61`'s config to make one.
- **Family:** `ARCH-1` step 1 (U16); attempt 1 of 4.

**Anchors (asserted, in W3, by T4 / T5 as committed in item 517).**
- **(E2) same input, same answer, degree 2:** ‖S_study − S_ops61‖_F / ‖S_ops61‖_F ≤ `LEG_D0_REPRODUCTION_BAND` (`tests/validation/test_port_birdcage_four_port.py:165`; the weekly's "≤ 1e-9"), S_study from `Study.solve(backend="palace", case_dir=D2_CASE)` and S_ops61 from `read_port_s("logs/ops61/postpro_f05x_o2/port-S.csv")` — one reader on both sides. **Predicted ≤ 1e-11** (the CSV prints 13 significant digits; MUMPS / GMRES can drift by an ULP across processes, the `OPS-43` observation, `docs/testing/known-issues.md:79`). Entry-wise max printed.
- **(E1) the CLI, degree 1:** the same against `postpro_f05x/port-S.csv`, from `summary.json` after `cli.main([... "--backend", "palace"])` returns 0.
- **(P)** on S_study and S_cli: (P-rec)'s set, bands imported as in item 517.
- **(Tie2) the tie to our degree-2 record:** per class \|mean_study − mean_P2\| / \|mean_P2\| against `ANS4_EXACT_STRIP_S_P2_X1["64"]` (`tests/validation/test_ans4_exact_strip_records.py:233`), printed in e-notation beside the recorded 2.464e-10 / 9.668e-11 / 9.883e-11 (`20261004T214342Z_OPS-61-item229-W4-read.log:41`). **Asserted ≤ 2.6e-10 only if the review rules R3 (a)**; otherwise printed, never asserted (2026-10-13 ruling (1)). Never against 1.9e-10.
- **(Tie1)** per class against `ANS4_EXACT_STRIP_S_P1_RECORDS["64"]` (`:141`): printed only (item 159's "0.0000 %", `20261002T022338Z_OPS-61-F10b-W6-ours-p1.log:1975`).
- **Palace log check (executor, journalled).** W1 / W2 each show 117 206 elements, `Rs = 4.375e+01 Ω/sq` on attributes 21–24, ω/2π = 6.400e-02 GHz, four `GMRES solver converged` lines, and `ND (p = 2): 753062` (W1) / `ND (p = 1): 140460` (W2) — item 194's numbers (`OPS-61.md:275`).

**Negative controls.**
- **Asserted (backed on the same fixture and comparison):** **(c-order)** S_study's class means against Palace order 1 (`postpro_f05x/port-S.csv`) differ by > 1e-6 in every class. Measured 6.6674 / 2.5442 / 4.5224 % (`20261003T110728Z_OPS-61-item194-W5-read.log:55`); the weakest class clears the bar by ≈ 2.5e4×. The ceiling is the order move itself. A reader blind to the case (reading a stale or wrong-order file) fails it.
- **Predicted, printed:** **(c-freq)** S_study against the 10 MHz order-2 read (`postpro_f05x_o2_f10/port-S.csv`): predicted O(1e-2 … 1e-1) Frobenius-relative (the 10 → 64 MHz change); not measured as a comparison, so never asserted.

**Windows.** Foreground, bare, `timeout: 660000`; `-np` ≤ 12; Palace v0.18.1, image `fem-em-palace:0.18.1`, profile `palace`.
- **Pre-check** (not harness calls): `docker compose -f docker/docker-compose.yml ps` shows `fem-em-solver` Up; `ls logs/u16/f_small_d2_64/ logs/u16/f_small_d1_64_cli/palace_case/` as above. A missing Palace image is a stop (2026-10-13 ruling (4)).
- **W1** (Palace, the Study's degree-2 case):
  `scripts/testing/run_and_log.sh ARCH-1-U16b-W1-palace-d2 "docker compose -f docker/docker-compose.yml --profile palace run --rm -T fem-em-palace bash -lc 'cd /workspace/logs/u16/f_small_d2_64 && timeout -k 30 590 palace -np 8 case.json'"`
- **W2** (Palace, the CLI's degree-1 case):
  `scripts/testing/run_and_log.sh ARCH-1-U16b-W2-palace-d1 "docker compose -f docker/docker-compose.yml --profile palace run --rm -T fem-em-palace bash -lc 'cd /workspace/logs/u16/f_small_d1_64_cli/palace_case && timeout -k 30 590 palace -np 8 case.json'"`
- After any kill: `docker compose -f docker/docker-compose.yml ps` for a lingering `fem-em-palace` container before W3.
- **W3** (the read-back, `-n 1`, knob on, T4 and T5 only):
  `scripts/testing/run_and_log.sh ARCH-1-U16b-W3-read "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_U16_PALACE_READ=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 400 mpiexec -n 1 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_arch1_palace_backend.py::test_u16_palace_read_back tests/package/test_arch1_palace_backend.py::test_u16_cli_palace_read_back -v -s --tb=short'"`

**Traps already paid for.**
- **"Same input" includes `-np 8`.** Item 194 ran `-np 8`; a different width changes Palace's partition and GMRES path. Never change `-np`, `Order`, `Tol` or `MaxIts` in-slot.
- **A `GMRES … did not converge` line is a finding**, not a retry; that run's S is "not measured".
- **An abort on the port set** (axis alignment) is a finding: never change `Direction`.
- **Palace writes `postpro/` relative to its cwd**; the window must `cd` into the case directory, as above.
- The test module is used as committed by item 517 (rule (i)): no edit, not even a comment.
- Compare S classes only; never read Palace's `z_matrix_ohm`. Record row and journal carry Palace's numbers in full (version, config, mesh cited) and our-vs-Palace gaps; **no AED figure and no Palace figure set against an AED one** in any tracked file (§5.4 (iv)).
- Status 124 on W1 is the price finding: record it; never lengthen the timeout.

**Scope.** One fixture (×1 exact strip, κ off), one frequency (64 MHz), orders 1 and 2, S only. **Not claimed:** any other rung or frequency, copper / hole / F-human, fields from Palace (SAR, B₁⁺, the loss partition, thermal — `ARCH-1` (1)–(2)), the benchmark dispatch (`ARCH-1` (4)), an in-container Palace launch, accuracy (a cross-solver identity on our discretisation, not an accuracy level), GPU / AMR.

**Negative result.** Report and stop; mark 517b BLOCKED in the same commit if anything parks.
- **(E2) or (E1) > `LEG_D0_REPRODUCTION_BAND` with (P) green and the log check green:** Palace does not reproduce its own recorded read on the same input at 1e-9. This is a reading, not a fix: DONE as a record with a known-issues entry ("Palace re-run not reproducible at 1e-9 on `f05x_o2`"); `ARCH-1` stays ⬜ and the review / weekly rules.
- **(E2) red and the log check red** (another element count, Rs, ND or frequency): the case is not `OPS-61`'s input despite item 517's (W-msh). That is a writer defect; BLOCKED.
- **(c-order) red:** the read is blind to the case; BLOCKED.
- **(Tie2) above 2.6e-10 with (E2) green:** impossible unless the record moved; print it and report.
- Palace aborts: record verbatim, known-issues entry, BLOCKED.
- No `PREDICTED RED:` line.

**Status it can move.** **`ARCH-1` ⬜ → 🟡** (the 2026-10-17 weekly: "Glyph 🟡 when U16 lands", PROJECT_PLAN.md:4709). Rule (j): the row has no Done-when (it is an embargo description), so the review writes the 🟡 sentence and lifts the "Not queueable" clause in the same commit: "🟡 step 1 (U16, §9 items 517 / 517b): `Study.solve(backend="palace")` and `femem run --backend palace` write the Palace case, a separate `fem-em-palace` window runs it, and the `Study` reads S back — equal to `OPS-61`'s recorded read on the same input to (E2) / (E1) at `LEG_D0_REPRODUCTION_BAND`, tied to our degree-2 record per class at (Tie2), ×1 exact strip at 64 MHz (Palace v0.18.1, orders 2 / 1, <W1 log>, <W2 log>, <W3 log>)". Not ✅: `ARCH-1`'s gate of record is the whole safety layer on two backends. §10 U16 → DONE; U62 (item 529) becomes runnable.

**Spawn line (if delegated):** "Use the implementer agent to execute `ARCH-1` §9 item 517b (U16b) per `docs/planning/chunks/ARCH-1.md`. Skip it unless item 517 is DONE on `main`. The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped; foreground, `timeout: 660000`. Make no edit to any file but the record lines. Never overwrite a `postpro/`; never change `-np`, `Order`, `Tol`, `MaxIts` or `Direction`. After a kill, check `docker compose -f docker/docker-compose.yml ps` for a lingering `fem-em-palace` container. You are the executor: do not spawn agents and never background a harness window. On a red, mark item 517b BLOCKED in the same commit and report."

**Outcome 517b (2026-10-18 14:07 UTC, 09:00 CDT slot): BLOCKED — the gate nodes are stubs.** W1 (`20261018T140040Z_ARCH-1-U16b-W1-palace-d2.log`, Status 0, 264 s): 117 206 elements, Rs = 4.375e+01 Ω/sq on 21–24, ω/2π = 6.400e-02 GHz, ND (p = 2): 753062, four `GMRES solver converged in 37 iterations`, Palace `Total` 261.6 s, peak 6.0 G. W2 (`20261018T140512Z_ARCH-1-U16b-W2-palace-d1.log`, Status 0, 24 s): the same element count, Rs and frequency, ND (p = 1): 140460, four converged (31 / 31 / 32 / 32 iterations), `Total` 22.96 s, 1.8 G. W3 (`20261018T140543Z_ARCH-1-U16b-W3-read.log`, 4 passed / 2 skipped, 3 s): item 517 committed `test_u16_palace_read_back` / `test_u16_cli_palace_read_back` as unconditional `pytest.skip("item 517b: asserts are written by 517b's brief")`. This brief, though, forbids any edit ("No edit to any file except the commit's record lines"; "The test module is used as committed by item 517"). So no anchor was asserted, and a skip closes nothing. **Record only, not a gate:** a gitignored scratch reader (`logs/arch1-u16b-record.py`, `read_port_s` + the module's `_p_rec`, run through the harness as `20261018T140633Z_ARCH-1-U16b-record-read.log`) reads: (E2) Frobenius-rel 0.000e+00, entry max 0.000e+00 (Palace re-run bitwise equal to `postpro_f05x_o2`); (E1) 2.149e-14 / 3.893e-14; (P) green on both: d2 reciprocity 9.848e-11, σ_max 0.999819887373, spreads 0.0371 / 0.0248 / 0.0500 %; d1 1.692e-10, 0.999760093669, 0.0852 / 0.0783 / 0.0921 %; (Tie2) 2.464e-10 / 9.668e-11 / 9.883e-11, equal to the recorded values and ≤ 2.6e-10; (Tie1) 2.242e-10 / 9.601e-11 / 1.823e-10; (c-order) 6.6674 / 2.5442 / 4.5224 %, reproducing the measured values; (c-freq) 8.515e-01. That reading compares CSVs; it does not go through `Study.solve`'s read path. **Side effect for the review:** `logs/u16/*/postpro/` now exist, as the design intends. Item 517's T1 (`D2_CASE` has no `postpro/` ⇒ `PalaceRunRequired`) and `test_u16_cli_write_stage_and_controls` (`assert not (CLI_CASE / "postpro").exists()`) will now go red if re-run. Those nodes need a ruling: move the write stage to `tmp_path`, or accept them as one-shot. **Unblock:** license a slot to write the T4 / T5 bodies to this brief's anchors (`Study.solve(backend="palace", case_dir=D2_CASE)` and the CLI `summary.json`), then run W3 only. Neither Palace window needs re-running.

**Rulings asked of the review (items 517 / 517b).**
- **R1 — the fixture.** "F26's recorded Palace read" is read as the ×1 exact-strip order-2 64 MHz read (F17, `postpro_f05x_o2`), the Study's `f_small()` verbatim; F37's h = 0.0075 hole Cu read would need three features U16 does not name. Ratify, or name the rung.
- **R2 — the three-stage design** (write → separate `fem-em-palace` harness window → read, manifest-checked) as the meaning of "runs the `fem-em-palace` service". The alternative — the `Study` launching Palace itself — needs a docker-socket mount or Palace in the solver image: an operator decision (compose constraint), not a slot's.
- **R3 — the tie.** (a) assert (Tie2) ≤ 2.6e-10, item 476's ratified public ceiling, as a transcription correction of the weekly's 1.9e-10 (rule (f), the 2026-10-15 03:00 precedent); or (b) print it beside the recorded per-class values, as every `OPS-61` tie has been.
- **R4 — the plan text and the CLI.** §5.4's "nothing in `src/` is built on Palace" (PROJECT_PLAN.md:1180–1181) amended to admit an opt-in writer / reader that never imports or links Palace (U8); the new CLI exit code 3 ("external run pending") accepted or replaced; the split 517 / 517b accepted. The stale compose / Dockerfile comments go to the operator.

### Review ruling 2026-10-18 10:30 (item 517b BLOCKED → item 533)

Read: the 517b outcome above, attempts.md 2026-10-18T14:07Z, `tests/package/test_arch1_palace_backend.py:47–90`, `:190–219`, `:283–290` (re-read), and `src/fem_em_solver/backends/palace.py:126–204` (`write_case` never overwrites; a differing file raises `FileExistsError`) and `:242–246` (the read path takes `postpro/port-S.csv` when present).
- **The block is a brief defect, not a solver finding.** Item 517 committed T4 / T5 as stubs whose bodies "are written by 517b's brief", and 517b's brief forbade any test edit. Neither brief owned the bodies. 517b is attempt 1 of 4 in the family; **item 533 is the re-pose (attempt 2)**, same target, under §9 rule (6)'s re-pose clause.
- **The Palace windows stand.** W1 / W2 ran on the item-517 cases with every log check green. Their `postpro/` stays on disk. No Palace re-run in 533.
- **The record read is a record.** (E2) 0.000e+00 compares CSVs, not `Study.solve`'s read path. It moves nothing.
- **T1 and `test_u16_cli_write_stage_and_controls` become repeatable.** Their write stage moves to `tmp_path`. The on-disk `D2_CASE` / `CLI_CASE` are never written again and keep their Palace output. The write is deterministic: `write_case` writes the same bytes, and T4's `Study.solve(case_dir=D2_CASE)` re-checks them against today's writer (a differing file raises). So the tmp case asserts what item 517 asserted, on every re-run. "One-shot" nodes are declined: a node that turns red because the design worked is not a gate.
- **R3 stands as ruled at 03:00.** (Tie2) is asserted at 2.6e-10. (c-order) is asserted at > 1e-6 per class. (c-freq) is printed.
- **Item 529 now depends on item 533.**

### §9 item 533 — 2026-10-18 (§10 U16, `ARCH-1` step 1b re-pose / U16c; chain step, no solve; serial on nothing)

*Written by the 10:30 review. Re-poses item 517b (attempt 2 of 4, family `ARCH-1` step 1).*

**`ARCH-1` step 1b / U16c: write T4 / T5's bodies to item 517b's anchors, move T1's and the CLI write stage to `tmp_path`, and run the module once, whole, with the read knob on.**
- **Executor:** implementer.
- **Edits allowed:** `tests/package/test_arch1_palace_backend.py` only.
  - **T1 (`test_u16_case_equals_ops61`)** takes `tmp_path` and writes to `tmp_path / "d2"` instead of `D2_CASE`. Every assertion stays, with `D2_CASE` replaced by the tmp case.
  - **`test_u16_cli_write_stage_and_controls`** takes `tmp_path` and passes `--output-dir tmp_path / "cli"`. Every assertion stays, on the tmp case.
  - **T4 (`test_u16_palace_read_back`)** and **T5 (`test_u16_cli_palace_read_back`):** the stub `skip` lines are replaced by bodies. The `skipif` on `FEM_EM_U16_PALACE_READ` stays. Each body skips "BLOCKED: … absent" when its case, its `postpro/port-S.csv` or its `OPS-61` comparand is missing. The module docstring's T4 / T5 sentence is updated.
  - Nothing in `src/`, no other test, no registry line. `logs/arch1-u16b-record.py` is the arithmetic to lift (gitignored; it stays where it is).
- **Tier:** heavy (one window, `-n 1`, container `timeout -k 30 590`). Price: item 517's W1 ran the module in 165 s at `-n 1` (`20261018T101308Z_ARCH-1-U16a-W1.log`). T4 and T5 each rebuild the F-small mesh once inside `write_case`'s check, as T2 / T3 do. **Predicted 200–330 s.**
- **Class:** chain step (no solve). **Cost:** ≈ 6.6 measured / ≈ 30 predicted slot-minutes.
- **Depends on:** nothing. On disk: `logs/u16/f_small_d2_64/{case.json, mesh.msh, case.sha256, postpro/port-S.csv}`, `logs/u16/f_small_d1_64_cli/palace_case/{…, postpro/port-S.csv}`, and `logs/ops61/postpro_f05x_o2/`, `postpro_f05x/`, `postpro_f05x_o2_f10/` (all present at 10:30, re-listed). A missing one is BLOCKED ("comparand absent"); never re-run Palace or `OPS-61` to make one.

**Anchors (asserted).**
- **T4 (E2):** `res = _d2_study().solve(MPI.COMM_WORLD, backend="palace", case_dir=D2_CASE)` returns a `PalaceStudyResult` with `cells == 117206` and `order == 2`. ‖S_res − S_ops61‖_F / ‖S_ops61‖_F ≤ `LEG_D0_REPRODUCTION_BAND` (`tests/validation/test_port_birdcage_four_port.py:165`). S_ops61 = `read_port_s(OPS61 / "postpro_f05x_o2" / "port-S.csv")`. Entry-wise max printed. The record read gives 0.000e+00.
- **T4 (P):** `_p_rec(res.s_matrix)`.
- **T4 (Tie2):** per class \|mean_res − mean_P2\| / \|mean_P2\| ≤ 2.6e-10 (the 03:00 R3 ruling; item 476's ratified ceiling). `mean_P2` is from `ANS4_EXACT_STRIP_S_P2_X1["64"]` (`tests/validation/test_ans4_exact_strip_records.py:233`) through `_circulant_classes` (the record script's `means`). The recorded 2.464e-10 / 9.668e-11 / 9.883e-11 (`20261004T214342Z_OPS-61-item229-W4-read.log:41`) are printed beside. Write the 2.6e-10 as a named module constant whose comment cites the ruling.
- **T5 (E1):** a stale `CLI_OUT / "summary.json"` is removed first (never the case or `postpro/`). Then `cli.main(["run", str(YAML), "--output-dir", str(CLI_OUT), "--backend", "palace"])` returns 0, `summary["backend"] == "palace"`, and S from `summary["s_matrix"]` (as T3 parses it) is within `LEG_D0_REPRODUCTION_BAND` of `read_port_s(OPS61 / "postpro_f05x" / "port-S.csv")`, Frobenius-relative. Then `_p_rec` on it.
- **T1 / CLI write on `tmp_path`:** item 517's assertions, unchanged in value.

**Negative controls.**
- **Asserted (c-order), in T4:** each class mean of S_res differs from Palace order 1 (`postpro_f05x`) by > 1e-6 relative. Backed on the same fixture and comparison: 6.6674 / 2.5442 / 4.5224 % (`20261003T110728Z_OPS-61-item194-W5-read.log:55`, reproduced by `20261018T140633Z_ARCH-1-U16b-record-read.log`). The weakest class clears 1e-6 by ≈ 2.5e4×.
- **Predicted, printed (c-freq), in T4:** Frobenius-relative distance to `postpro_f05x_o2_f10`. The record read gives 8.515e-01. Never asserted.
- **Inherited:** T2's (c3) / (c4) stay as committed.

**Window.** One harness call, foreground, the whole Bash command, `timeout: 660000`:

`scripts/testing/run_and_log.sh ARCH-1-U16c-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_U16_PALACE_READ=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 590 mpiexec -n 1 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_arch1_palace_backend.py -v -s --tb=short'"`

Pre-check: `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reads 0 (no heavy window running). **Expected:** every node passes, none skipped (T0, T1, CLI write, T2, T3, T4, T5, Rd ×3).

**Traps already paid for.**
- **Rule (i):** no edit to the module after W1. An edit means W1 again.
- `-n 1` only: the backend refuses a parallel comm by design (T0 asserts it).
- **The test never writes into `D2_CASE` or `CLI_CASE`** except T5's `summary.json` under `CLI_OUT` (gitignored). `write_case`'s own check stages and removes a `mesh.msh.new` and rewrites an identical `case.sha256` (`palace.py:187–205`); that is the package's read path, not a test write. Never delete or overwrite a `postpro/`.
- No `-k` selector and no pipe in the harness command. `-s` is on.
- Compare S classes only. Never read Palace's `z_matrix_ohm`. No AED figure, and no Palace figure set beside an AED one (§5.4 (iv)). If the private-leak guard flags the log, move it unchanged to `logs/leak-held/`, mark the item BLOCKED and stop.

**Scope.** U16 on one fixture: the ×1 exact strip at 64 MHz, orders 2 / 1, S only. Not claimed: anything in 517b's scope list, a Study that launches Palace, or accuracy.

**Negative result.** Report and stop. Mark item 533 BLOCKED in the same commit (rule (d)).
- (E2) / (E1) red with (P) green: a known-issues entry ("Study read path differs from the CSV read"). The record read says the CSVs agree, so the defect is in the read path. `ARCH-1` stays ⬜.
- `FileExistsError` from T4 / T5: the on-disk case differs from today's writer. Something changed `write_case` or the fixture after item 517. BLOCKED, named; never delete the case to clear it.
- (c-order) red: the read is blind to the case. BLOCKED.
- No `PREDICTED RED:` line.

**Status it can move.** **`ARCH-1` ⬜ → 🟡**, with the sentence item 517b's brief wrote ("🟡 step 1 (U16, §9 items 517 / 517b / 533): …"), citing 517b's W1 / W2 logs and this item's W1 log. The row's "Not queueable" clause is lifted in the same commit. §10 U16 → DONE. Item 529 (U62) becomes runnable.

**Spawn line (if delegated):** "Use the implementer agent to execute `ARCH-1` §9 item 533 (U16c) per `docs/planning/chunks/ARCH-1.md` (heading `§9 item 533` and the 2026-10-18 10:30 ruling above it). The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped; foreground, `timeout: 660000`. Edit only `tests/package/test_arch1_palace_backend.py`. Never write into `logs/u16/*` except T5's `summary.json`, and never touch a `postpro/`. You are the executor: do not spawn agents and never background a harness window. On a red, mark item 533 BLOCKED in the same commit and report."

**Outcome 533 (2026-10-18 17:28 UTC, 12:00 CDT slot, take-next): DONE, attempt 2 of 4.** Edits only in `tests/package/test_arch1_palace_backend.py`.
- T1 and `test_u16_cli_write_stage_and_controls` now write to `tmp_path`; every assertion is unchanged in value.
- T4 / T5 carry bodies. `TIE2_BAND = 2.6e-10` is a named constant citing the R3 ruling.
- W1 (`20261018T172228Z_ARCH-1-U16c-W1.log`, Status 0, 224 s, `-n 1`): **14 passed, none skipped**.
  - T4: `cells` 117 206, order 2. (E2) Frobenius-rel 0.000e+00, entry max 0.000e+00, through `Study.solve(backend="palace", case_dir=D2_CASE)`.
  - (P): reciprocity 9.848e-11, σ_max 0.999819887373. (Tie2) 2.464e-10 / 9.668e-11 / 9.883e-11, equal to the record. (c-order) 6.6674 / 2.5442 / 4.5224 %.
  - (c-freq) 8.515e-01, printed only. The print line's "predicted O(1e-2 .. 1e-1)" text came from 517b's brief; the record read had already measured 0.85, and rule (i) left the line unedited.
  - T5: `cli.main` rc 0, `backend` palace, (E1) 2.149e-14 (entry max 3.893e-14). (P): reciprocity 1.692e-10, σ_max 0.999760093669.
- No `postpro/` was touched. `ARCH-1` → 🟡, and §10 U16 is DONE. Item 529 (U62) is runnable.

### §9 item 529 — 2026-10-18 (§10 U62, Palace backend guide page; lane 1, no solve)

*Drafted by a read-only general-purpose agent for the 2026-10-18 03:00 daily review; the review reviews it.*

**Licence.** §10 U62 (PROJECT_PLAN.md:4871–4872), weekly-enumerated (§9 rule (4)): "A 'backends' guide page and the envelope entry for the Palace backend; `docrefs` 0 dead. Serial on U16." §9's "Every item" rule on `WF-14`'s registry (PROJECT_PLAN.md:4058: regenerate the API pages and `docs/guide/validation.md` and run the page gates in the same slot). Templates: §9 item 433 (`docs/planning/chunks/WF-14.md:196–285`, the registry + regeneration windows) and `WF-17` step 6 (`tests/package/test_wf17_tuning_yaml_page.py`, the page-as-test). **No band, tolerance, Done-when or target is set or moved.**

**`ARCH-1` step 1c / U62: `docs/guide/backends.md`, a page test that keeps it true, and the U16 nodes registered under `Study` / `run_study_file`.**
- **Executor:** implementer.
- **Tier:** smoke / standard, no MPI (`-n 1`): the page test ≈ 5 s (`20261015T214647Z_WF-17-s6-W0.log`), Wapi 17 s (`20261015T214657Z_WF-17-s6-Wapi.log`), Wdoc 2 s (`20261015T214715Z_WF-17-s6-Wdoc.log`), the generators ≈ 2–3 s each (`20261014T021649Z_WF-14-u48-Wgen.log`, `…021653Z…Wval.log`).
- **Class:** lane 1. **Cost:** ≈ 3.5 measured / ≈ 25 predicted slot-minutes: 15 fixed + ≈ 6 page and test + ≈ 1 windows + commit ≈ 3.
- **Depends on:** **§9 item 517b DONE on `main`** (§10: serial on U16); skip otherwise. If 517b landed as a record with (E2) red or R3 undecided, the page states what landed and no more; if 517b is BLOCKED, skip.
- **Family:** `ARCH-1` step 1 (U62); attempt 1 of 4.

**Edits allowed.**
- **New `docs/guide/backends.md`.** Sections:
  - what a backend is — the DolfinX default (reference, coupling host) and the opt-in Palace backend for scale, U8's roles;
  - the three stages, as three tagged blocks: `python` `palace-backend` (one `Study.solve(` call with `backend="palace"` and a `case_dir`), `bash` `palace-backend-cli` (`femem run <file> --backend palace`, then the Palace command, then the same `femem run` again), `text` `palace-command` (the `.argv` line);
  - what is refused, and why the Study cannot launch Palace itself (separate image, no socket);
  - the binding condition: every Palace-backed result names a reduced matched-order cross-check on a shared mesh (U8);
  - the cross-check of record: U16's nodes, cited as `` `tests/…::…` ``, with links to `validation.md`, `api/Study.md`, `api/run_study_file.md` and the public `ans:4` `COMPARISON.md` "Palace at matched order" section.
  - **No number is typed into the page.** The tolerance line is rendered from `LEG_D0_REPRODUCTION_BAND` and the command line from `backends.palace` at test time ((E) below). No AED figure, and no Palace figure beside an AED one.
- **New `tests/package/test_wf17_backends_page.py`**, on `test_wf17_tuning_yaml_page.py`'s pattern (`extract_tagged_block` from `tests/package/test_wf17_first_study_tutorial.py`, `run_checker` from `tests/unit/test_doc_reference_exit_codes.py`). Docstring paragraph 1 names `ARCH-1`.
- **`tests/package/validation_envelope.toml`:** one comment line `# U16 Palace backend (§10, weekly 2026-10-17; §9 items 517 / 517b)` plus the U16 nodes under `[names.Study]` (`:151`): `test_arch1_palace_backend.py::test_u16_backend_refusals`, `::test_u16_case_equals_ops61`, `::test_u16_read_path_on_recorded_output`, `::test_u16_palace_read_back`; and under `[names.run_study_file]` (`:316`): `::test_u16_cli_read_path_on_recorded_output`, `::test_u16_cli_palace_read_back`. Function-level IDs, no `[param]` suffix, no number in the comment (`test_registry_restates_no_number`). This attribution is the drafter's proposal; the review rules on it.
- The regenerated `docs/guide/api/*.md` and `docs/guide/validation.md` (generated, never hand-edited).
- Nothing in `src/`, no gate module, no generator.

**Anchors (asserted).**
- **(X)** one block per tag; the `python` block parses with exactly one `Study.solve` call carrying `backend="palace"`; the `bash` block contains `femem run` and `--backend palace` twice.
- **(E)** the `palace-command` block equals `" ".join(PalaceRunRequired(...).argv)` built from `backends.palace` (rendered, not typed); the page's tolerance line equals `f"{LEG_D0_REPRODUCTION_BAND:.0e}"` (`tests/validation/test_port_birdcage_four_port.py:165`).
- **(L)** every cited node is registered under `("Study", "run_study_file")`; links to `validation.md`, `api/Study.md` and `api/run_study_file.md` resolve, with validation anchors on headings.
- **(D)** `run_checker --docs-root docs/guide --check-links` reads 0 dead (the U62 target).
- **Wapi:** `test_wf14_validation_envelope.py` collects the new nodes (S / R / C / N); `test_wf17_api_reference.py` (G) finds the regenerated pages byte-identical; `test_wf14_validation_page.py` passes (I) with `ARCH-1` resolved from the new module's docstring.
- **Wval** prints the node count. Read the current count first; predict +6. Printed only.

**Negative controls (asserted, by construction).**
- **(c1)** the `palace-command` block with `-np 8` → `-np 4` differs from (E) on exactly one line.
- **(c2)** a planted unregistered node is named by (L).
- **(c3)** a tmp copy of the guide tree with one planted dead link reads dead == 1 (the `test_wf17_tuning_yaml_page.py:197–215` construction).
- **(C-api) / (C-val):** after the registry edit and **before** Wgen / Wval, both `--check` runs exit 1, naming `Study.md` and `run_study_file.md` (`generate_api_reference.py`'s differing-files list) and `DIFFERS: validation.md`. An exit 0 means the lines did not reach the pages: stop.

**Windows.** Foreground, bare, repo-relative, Bash `timeout: 660000`, `-s`; item 433's commands with this item's labels.
1. Write the page and the test; edit the registry.
2. **C-api** (predicted Status 1):
   `scripts/testing/run_and_log.sh ARCH-1-U62-ctl-api "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 python3 scripts/docs/generate_api_reference.py --out docs/guide/api --check'"`
3. **C-val** (predicted Status 1):
   `scripts/testing/run_and_log.sh ARCH-1-U62-ctl-val "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace timeout -k 30 60 python3 scripts/docs/generate_validation_page.py --out docs/guide/validation.md --check'"`
4. **Wgen:**
   `scripts/testing/run_and_log.sh ARCH-1-U62-Wgen "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 python3 scripts/docs/generate_api_reference.py --out docs/guide/api'"`
5. **Wval:**
   `scripts/testing/run_and_log.sh ARCH-1-U62-Wval "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace timeout -k 30 60 python3 scripts/docs/generate_validation_page.py --out docs/guide/validation.md'"`
6. **W0** (the page test):
   `scripts/testing/run_and_log.sh ARCH-1-U62-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace timeout -k 30 120 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf17_backends_page.py -s -v --tb=short'"`
7. **Wapi:**
   `scripts/testing/run_and_log.sh ARCH-1-U62-Wapi "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace timeout -k 30 170 python3 -m pytest tests/package/test_wf14_validation_envelope.py tests/package/test_wf17_api_reference.py tests/package/test_wf14_validation_page.py -v -s --tb=short'"`
8. **Wdoc:**
   `scripts/testing/run_and_log.sh ARCH-1-U62-Wdoc "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && timeout -k 30 60 python3 scripts/testing/check_example_doc_references.py --docs-root docs/guide --check-links'"`
9. Commit the page, the test, the registry, the regenerated pages and the logs together.

**Traps already paid for.**
- **A control window's Status 1 is its pre-registered result**, labelled `-ctl-`; it is not a red.
- **Never hand-edit a generated page**; never edit either generator or `test_wf17_api_reference.py`.
- **An uncollectable node is a stop.** The envelope test collects with MPI env stripped; the U16 module must import cleanly on the real build (`@complex_only` on its mesh nodes, as item 517 wrote it).
- **No number in a registry comment**; no number typed into the page outside a rendered block.
- **No AED figure and no Palace-beside-AED sentence** (§5.4 (iv); the leak hook).
- `≈`, never a tilde, in the page. No `-k`, no pipe, never delete a log, no edit after the last window (rule (i)).
- **Registry conflicts:** if another item touched the registry, a rebase keeps both sides and windows 2–8 re-run.

**Scope.** One guide page, its test, six registry lines, regenerated pages. **Not claimed:** that a node is green today (the registry asserts collectability), any figure beyond the rendered constant, any change to `WF-14`'s Done-when or the public surface, a tutorial that runs Palace (a follow-on example, the daily's to open after U16).

**Negative result.** Report and stop; park on `attempt/ARCH-1-U62-*` and mark item 529 BLOCKED in the same commit (rule (d)).
- (C) collection red on a U16 node: BLOCKED for the review's re-attribution.
- (C-api) / (C-val) exit 0: the edit did not reach the pages; stop.
- (G) red after Wgen: regenerate once after a rebase; otherwise park.
- (D) red: a dead link; fix in-letter once, else park.
- No `PREDICTED RED:` line; the control windows' Status 1 is their asserted result.

**Status it can move.** No glyph (`ARCH-1` 🟡 from item 517b; `WF-14` ✅, Done-when met since 2026-10-04, no rule (j) amendment). The `ARCH-1` row gains "step 1c (§9 item 529, U62): `docs/guide/backends.md` and the U16 envelope entry, 0 dead links (<Wdoc log>)"; §10 U62 → DONE; (9d) gains the backends page.

**Spawn line (if delegated):** "Use the implementer agent to execute `ARCH-1` §9 item 529 (U62) per `docs/planning/chunks/ARCH-1.md`. Skip it unless item 517b is DONE on `main`. The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped; foreground, `timeout: 660000`, pytest with `-s`, container `timeout -k 30`. Run C-api and C-val before Wgen / Wval and require Status 1 on both; never hand-edit a generated page; no number in a registry comment or typed into the page; no AED figure. You are the executor: do not spawn agents. On a red, park on `attempt/ARCH-1-U62-*`, mark item 529 BLOCKED in the same commit, and report."

**Outcome 529 (2026-10-19 00:53 UTC, 00:30 UTC slot): DONE, attempt 1 of 4.** Dependency read as item 533 (10:30 ruling).
- New `docs/guide/backends.md` and `tests/package/test_wf17_backends_page.py`. Six U16 nodes registered: four under `[names.Study]`, two under `[names.run_study_file]`, all six names verified in `test_arch1_palace_backend.py`. The comment line names items 517 / 517b / 533 (533 added to the brief's text).
- C-api (`20261019T005017Z_ARCH-1-U62-ctl-api.log`, Status 1, 2 s): differing `['Study.md', 'run_study_file.md']`. C-val (`20261019T005022Z_ARCH-1-U62-ctl-val.log`, Status 1, 2 s): `DIFFERS: validation.md`.
- Wgen (`20261019T005038Z…`, Status 0, 2 s). Wval (`20261019T005044Z…`, Status 0, 3 s): nodes 289 → 295, as predicted (+6, printed only).
- W0, first run (`20261019T005051Z_ARCH-1-U62-W0.log`, Status 1): (D) and (c3) red. The link pass read 0 dead, but the reference pass found 2 dead: `summary.json` (prose) and `case.json` (inside the rendered `palace-command` block, which must equal `PalaceRunRequired.argv`). Fixed in-letter once:
  - the prose was reworded;
  - **`case.json` got an `ALLOWLIST` entry, with its reason, in `scripts/testing/check_example_doc_references.py`.** That is the checker's own exemption mechanism (precedent `lineplot.csv`). It is outside the brief's edit list and disclosed here. Neither generator was touched.
- W0, re-run (`20261019T005151Z_ARCH-1-U62-W0.log`, Status 0, 5 s): **9 passed**. (X) 4 tags, 1 `solve` call, `femem run` ×2 with `--backend palace`. (E) `palace -np 8 case.json` = argv, and the tolerance `1e-09` = `LEG_D0_REPRODUCTION_BAND` at `.0e`. (c1) differs on 1 line. (L) 6 nodes, 0 unregistered. (c2) planted node named. (D) dead 0. (c3) 0 → 1.
- Wapi (`20261019T005208Z_ARCH-1-U62-Wapi.log`, Status 0, 19 s): 22 passed.
- An extra window, Wchk (`20261019T005236Z_ARCH-1-U62-Wchk.log`, Status 0, 15 s), ran `tests/unit/test_doc_reference_exit_codes.py` (20 passed) because the checker was edited.
- Wdoc (`20261019T005254Z_ARCH-1-U62-Wdoc.log`, Status 0, 2 s): `dead=0`, 1 allowlisted, link pass 0 dead.
- No edit after Wdoc. §10 U62 is DONE. No (9d) page list exists in the plan, so none was extended.

### Review ruling 2026-10-18 03:00 (items 517, 517b, 529 — U16 split, U62)

Accepted as split. Spot-check: `docker/docker-compose.yml:8–10` mounts only `..:/workspace` and the cache volume, so there is no docker socket in the solver container (re-read).
- **R1: accepted.** The fixture is the ×1 exact strip at 64 MHz, order 2 (`postpro_f05x_o2`). "F26's recorded Palace read" names no Palace read of the Study's default fixture. The target is unchanged in kind ("same input, same answer" at ≤ 1e-9 against the recorded Palace read), and the swap is named for the weekly to ratify.
- **R2: the three-stage design is accepted:**
  - `Study.solve(backend="palace", case_dir=…)` writes the case and raises `PalaceRunRequired`;
  - a separate harness window runs the `fem-em-palace` service;
  - the same call reads the case back, after a checksum check.
  
  A Study that launches Palace itself would need a socket mount or Palace inside the solver image. Both are the operator's (the §9 compose constraint) and are named on the dashboard. Exit code 3 for `femem run --backend palace` is accepted.
- **R3:**
  - The same-input equality is **asserted** at `LEG_D0_REPRODUCTION_BAND` 1e-9.
  - The tie to our degree-2 record is **asserted at the ratified 2.6e-10 ceiling** (item 476's ruling). The recorded 2.464e-10 / 9.668e-11 / 9.883e-11 (`20261004T214342Z_OPS-61-item229-W4-read.log:41`) are printed beside it.
  - §10's "≤ 1.9e-10" is a stale figure, named for the weekly. No band is set here.
- **R4 — named for the weekly, not edited by a slot:**
  - §5.4's "nothing in `src/` is built on Palace" becomes false when 517b lands;
  - the `ARCH-1` row's "EMBARGOED" text;
  - the stale compose / Dockerfile comments, which are the operator's files.
- **529:** filing the U16 test nodes under `Study` and `run_study_file` in the registry is accepted. The WF-14 page-gate rule applies (§9 "Every item").
- **`ARCH-1` ⬜ → 🟡 when 517b lands**, per the weekly.

**Outcome 517 (2026-10-18 10:16 UTC): DONE.** `src/fem_em_solver/backends/palace.py` (writer / config / reader lifted; `solve_palace`), `Study.solve(backend=, case_dir=, palace_ranks=)`, `femem run --backend palace` (rc 3 on "run required"; via private `_run_study_file_palace`, so `run_study_file`'s signature and body are unchanged — a first W4 had gone red on `run_study_file.md` when `backend` was a `run_study_file` kwarg, fixed in-letter by removing the leak). W1 (`20261018T101308Z_ARCH-1-U16a-W1.log`, `-n 1`, 12 passed / 2 skipped, 165 s): (W-cfg) differing keys exactly {Problem.Output, Model.Mesh} for both cases, (c1) adds Solver.Order; (W-msh) `mesh.msh` byte-identical to `ops61_ans4_x1_strip.msh` (sha256 b77a7648…6c82, 117 206 tets) for both the Study and CLI cases; (c2) ragged 116 085 tets; (Rd) `np.array_equal` on all three recorded CSVs; (P-rec) on f05x_o2 reciprocity 9.848e-11, σ_max 0.999819887373 (= item 194), spreads 0.0371 / 0.0248 / 0.0500 %, pooled 844.9× worst; (c3) 0.01 vs 0.064 refused; (c4) one-byte edit → FileExistsError naming case.json; T0 refusals 0.000 s each, `-n 2` refused by a fake comm. W2 (`20261018T100621Z_ARCH-1-U16a-W2-regate-wf16.log`, 19 passed, 172 s; (R1) 9.559e-11 vs record 9.557e-11 in band 1e-9, σ_max 0.999760093664, spreads 0.0747 / 0.0691 / 0.0897 %); W3 (`20261018T100919Z_ARCH-1-U16a-W3-regate-wf1.log`, 13 passed, 166 s, (R1) 9.557e-11); W4 (`20261018T101247Z_ARCH-1-U16a-W4-wapi.log`, 22 passed). Disclosed: W2 / W3 ran before the `study_yaml.py` / `cli.py` fix that moved the palace branch out of `run_study_file`; the fix returns `run_study_file` to its pre-edit text and leaves `Study.solve` untouched since W2, so they were not re-run (tight clock). Also run: the reader-only smoke `20261018T100033Z_ARCH-1-U16a-Rd-partial.log`, the first W1 `20261018T100133Z_ARCH-1-U16a-W1-partial.log` and `20261018T100319Z_ARCH-1-U16a-W1.log` (pre-fix, green), the red W4 `20261018T101209Z_ARCH-1-U16a-W4-wapi.log`. `logs/u16/f_small_d2_64/` and `logs/u16/f_small_d1_64_cli/palace_case/` hold case.json / mesh.msh / case.sha256, no postpro/ — ready for 517b.

### Review ruling 2026-10-19 03:00 (item 529 deviation — `case.json` in the checker allowlist)

**The global `ALLOWLIST` entry stands.** Item 529 added `case.json` to `scripts/testing/check_example_doc_references.py` `ALLOWLIST` with its reason, after the first W0 (`20261019T005051Z`) read two dead references, and re-gated the checker (`tests/unit/test_doc_reference_exit_codes.py`, 20 passed, `…005236Z_ARCH-1-U62-Wchk`). The allowlist is the checker's only exemption mechanism, and `lineplot.csv` is the precedent: a filename that names a user-run output, never a tracked file. A page-local exemption would be a new checker feature, and nothing here needs one. **Revisit only if** a tracked file named `case.json` appears in the repo (the entry would then mask a dead reference to it). The disclosure and the extra Wchk window satisfy §9 rule (c).

### §9 item 590 — 2026-10-22 (§10 U76, chain step (no solve); chain U, 9e)

**Outcome (2026-10-22, 17:00 UTC slot, second item): DONE.** `_run_palace` now admits only `_PALACE_ROUTED_OUTPUTS = ("directory",)` and guards `cooling`; field `--backend` on every refusal. W0 real `-n 1` 3 passed / 3 skipped, 3 s (`20261022T171648Z_ARCH-1-U76-W0-real.log`): (X) 15 rows, field as predicted, rc (2, 2), nothing written, 0 builds; (Th) `palace: thermal is not supported`; (A) one solve each, EM keys `{s_matrix}`; (c-drop) widened allow-list reaches the spy once with no `probes` key; (Pg) mismatch `[]`; (c-pg) `['refused: missing outputs.probes']`. One read change: `f_small_sar_b1.yaml` now names `outputs.b1_plus` (sorted keys), not `outputs.sar`; field unchanged. RG-parse complex `-n 2` 22 passed, 4 s (`20261022T171704Z_ARCH-1-U76-RGparse.log`). `pgrep -c python3` 0, then RG-U16 complex `-n 1` 14 passed, 220 s (`20261022T171719Z_ARCH-1-U76-RGU16.log`): (E2) 0.000e+00, (Tie2) 2.464e-10 / 9.668e-11 / 9.883e-11, (c-order) 6.6674 / 2.5442 / 4.5224 %, (E1) 2.149e-14, `palace_total_s` 261.619549041. W-bit 2 passed (`20261022T172109Z_ARCH-1-U76-Wbit.log`): (c-bit) `['[U16c T4 Tie2]']` unequal; (Bit) all four lines string-equal to U16c's. RG-U62 9 passed, (L) 6 nodes, (D) 0 dead (`20261022T172119Z_ARCH-1-U76-RGU62.log`). C-api Status 1 `['StudyConfigError.md']`, C-val Status 1 `DIFFERS: validation.md`; Wgen 46 files; Wval nodes 333 → 334 (+1 as predicted); Wapi 22 passed (21 s); Wdoc `dead=0`. No Palace run, no `postpro/` written. Out of scope per ruling (1), unchanged: hole-mesh / degree-3 / `-n 2` Palace calls still exit by `ValueError` traceback in `solve_palace`, not rc 2.

*Drafted read-only for the 2026-10-22 03:00 daily review; the review rules on it.*

**Licence.** §10 chain U, U76, weekly-enumerated 2026-10-21 (§9 rule (4)), quoted: "**U76 (9e; standard, real, `-n 1`).** `_run_palace` refuses `outputs.probes`, `sar`, `b1_plus`, `tune`, `thermal` and every output it has no route for, by a named field, rc 2 from `femem run --backend palace`; the backends page says which outputs Palace carries; every Palace-route gate re-run (U16's (E2) / (Tie2) bitwise). Serial on nothing." (`PROJECT_PLAN.md:4958–4963`). The weekly posed it as one of "two silent-drop defects (read, not observed)" (`PROJECT_PLAN.md:4819–4823`), from the 2026-10-20 10:30 ruling (4) in `docs/planning/chunks/WF-1.md:3021` ("`_run_palace` refuses `sar` / `b1_plus` / `tune` but not `outputs.probes` … drops the probes silently. Found by reading the code, not observed in a run"). Queue position: after e4 / 3b-iii / U74 / U75 (`PROJECT_PLAN.md:5112–5115`); the weekly's Palace line: "U76 makes the backend's refusals explicit" (`:5095–5096`). The newest ruling in this file (2026-10-19 03:00, the `case.json` allowlist) stands and is not touched. **No band, tolerance, Done-when or validation target is set or moved.** This is `ARCH-1` **step 1d**.

**What the code is today (read 2026-10-22, `11fc13b0`).**
- `_run_palace(study, kind, outputs, path, out, digest, started, comm)` (`src/fem_em_solver/study_yaml.py:1782–1820`) refuses, each with field `--backend`: a family other than `four_port_birdcage` (`:1784–1785`); `outputs.sar` / `b1_plus` / `tune` (`:1786–1788`, a `for key in (…)` deny-list); `study.sweep_hz` (`:1789–1790`, item 572); a conductor (`:1791–1792`); `CapacitorSheets` (`:1793–1794`). Then `study.solve(comm, backend="palace", case_dir=out / "palace_case")` (`:1795`) and a summary whose only EM output is `s_matrix` (`:1796–1813`).
- **The silent drop:** a four-port `outputs` dict can carry `directory`, `sar`, `b1_plus`, `probes` (`_parse_outputs`' `_unknown` tuple, `:1142`), `tune` (`:1339–1343`) and `sweep_hz` (`:1344–1351`). `probes` is on none of `_run_palace`'s lists, so `femem run --backend palace tests/package/wf1/f_small_em_probes.yaml` today writes the case and returns 3, and the read stage would return 0 with no `probes` key. A deny-list also drops any output key added later.
- **`thermal` is refused already, but not by `_run_palace`:** `_run_study_file_palace` (`:1855–1869`) unpacks `_load`'s fourth value as `_cooling` and never passes it on. A `thermal:` section on a `four_port_birdcage` is refused at parse with field `thermal` (`_parse_thermal`, `:603–606`; `tests/package/wf1/malformed_thermal_on_port.yaml`). On `f_human_birdcage`, `thermal_ball` and `thermal_box` files the family check refuses it with `--backend`. So `thermal` has no file-reachable silent drop today. It is unguarded only structurally.
- `cli.main` maps `StudyConfigError` to rc 2 and prints `femem: error: …` on rank 0's stderr (`src/fem_em_solver/cli.py:203–206`); the Palace branch (`:189–200`) returns 3 on `PalaceRunRequired`. **No `cli.py` edit is needed.**
- **No existing test asserts `_run_palace`'s sar / b1_plus / tune refusals.** Only the sweep one is asserted: `tests/package/test_wf1_yaml_sweep.py:187–194` asserts field `--backend` and nothing written. U16's T0 (`tests/package/test_arch1_palace_backend.py:117–156`) covers `Study.solve`'s refusals, not `_run_palace`'s. No test asserts any of these message texts (grep `is not supported` under `tests/`: no hit).

**`ARCH-1` step 1d / U76: `_run_palace` admits only what Palace carries and refuses every other output by name, with field `--backend` and rc 2. The backends page states the carried / refused split, and a test renders it from the code. The U16 gates are re-run and their (E2) / (Tie2) / (c-order) / (E1) lines must be byte-identical to U16c's.**
- **Executor:** implementer.
- **Tier:**
  - W0: standard, real build, `-n 1`, no mesh.
  - RG-parse: smoke, complex, `-n 2`, item 572's landed parse bundle.
  - RG-U16: **heavy by ceiling**, complex, `-n 1`, 224 s measured (`20261018T172228Z_ARCH-1-U16c-W1.log`, `Elapsed (s): 224`, container `timeout -k 30 590`). The weekly's "standard" fits the new node, not this re-run. **No Palace window.**
  - The rest are smoke.
- **Class:** chain step (no solve). **Cost:** ≈ 8.0 measured / ≈ 37 predicted slot-minutes: 15 fixed, ≈ 10 code (≈ 15 `src` lines, one module, one page section), ≈ 6 windows (≈ 4.7 min summed below plus harness overhead), ≈ 6 commit and journal.
- **Depends on:** nothing in §9. On disk (gitignored, all present 2026-10-22, listed): `logs/u16/f_small_d2_64/{case.json, mesh.msh, case.sha256, postpro/port-S.csv}`, `logs/u16/f_small_d1_64_cli/palace_case/{…, postpro/port-S.csv}`, and `logs/ops61/postpro_f05x_o2/`, `postpro_f05x/`, `postpro_f05x_o2_f10/`. Any one absent ⇒ BLOCKED ("comparand absent"). Never re-run Palace or `OPS-61` to make one.
- **Family:** `ARCH-1` step 1d (U76); attempt 1 of 4.
- **Daytime:** no.

**Does the re-run need Palace?** No. Every Palace-route gate reads Palace's *recorded* output; none runs Palace:
- T4 reads `logs/u16/f_small_d2_64/postpro/port-S.csv` through `Study.solve(backend="palace", case_dir=D2_CASE)`.
- T5 reads `logs/u16/f_small_d1_64_cli/palace_case/postpro/port-S.csv` through `femem run --backend palace`.
- T2 / T3 read copies of `logs/ops61/postpro_*` in `tmp_path`.

For the record only, a Palace window would be a host-side harness call into the `palace` profile, e.g. item 544's (`EX-134`) `docker compose -f docker/docker-compose.yml --profile palace run --rm -T fem-em-palace bash -lc 'cd <case> && timeout -k 30 590 palace -np 8 case.json'`:
- `20261019T093301Z_EX-134-palace-d2.log:12`, 262 s, order 2;
- U16b W1 264 s (`20261018T140040Z_ARCH-1-U16b-W1-palace-d2.log:373`);
- U16b W2 24 s at order 1 (`20261018T140512Z_ARCH-1-U16b-W2-palace-d1.log:347`).

**This item runs none, and never writes a `postpro/`.**

**Code (`src/fem_em_solver/study_yaml.py` only; disclosed under rule (c)).**
- **New module-level private constants**, beside `_run_palace`:
  - `_PALACE_ROUTED_OUTPUTS = ("directory",)` (the only `outputs` key the Palace route honours);
  - `_PALACE_OUTPUT_FIELD = {"tune": "tune", "sweep_hz": "study.sweep_hz"}` (display names for keys that are not `outputs.*`).
  - Nothing goes into `__all__`, so `WF-14`'s snapshot is unchanged.
- **`_run_palace`:**
  - **Replaced:** the `for key in ("sar", "b1_plus", "tune")` deny-list and the `sweep_hz` line (`:1786–1790`). In their place, one allow-list loop: `for key in sorted(outputs): if key not in _PALACE_ROUTED_OUTPUTS: raise StudyConfigError(path, "--backend", f"palace: {_PALACE_OUTPUT_FIELD.get(key, 'outputs.' + key)} is not supported")`.
  - **New keyword-only argument `cooling=None`**, with `if cooling is not None: raise StudyConfigError(path, "--backend", "palace: thermal is not supported")` placed before the loop.
  - The family, conductor and `CapacitorSheets` lines stay byte-unchanged.
  - The field stays `--backend` on every refusal: item 572's existing assertion `g_pal == "--backend"` (`test_wf1_yaml_sweep.py:193`) must not move, and every existing refusal already uses that field. The message names the key.
- **`_run_study_file_palace`:** `_cooling` → `cooling`, passed as `cooling=cooling`.
- **Never:**
  - touch `run_study_file` (its body or its docstring; the API page renders the docstring), the module docstring, `_parse_*`, `cli.py`, `study.py`, `backends/palace.py`, or any existing test;
  - edit `test_arch1_palace_backend.py` or `test_wf17_backends_page.py`. Their windows re-run them as committed.

**New `tests/package/test_arch1_u76_palace_outputs.py`.** Docstring paragraph 1 names `ARCH-1` and `WF-1` (the validation page's (I) resolves §7 IDs from it). It imports cleanly in both builds and builds no mesh. Every `*Geometry.build` is patched to raise, and 0 calls are asserted (item 572's stub, `test_wf1_yaml_sweep.py:126–136`). Derived YAML files are written from `_mapping` (`tests/package/test_wf1_study_yaml.py:164`) with `yaml.safe_dump`. Rank 0 writes them into the gitignored scratch `logs/arch1/u76/` after `_fresh_dir` (`:168`), followed by a `Barrier`. Every `--output-dir` lives under that scratch too, because the refusal node also runs at `-n 2` and `tmp_path` is per-process (item 572's `LOGS / "s21" / "palace_refused"` precedent). `_field` is defined locally (the `test_wf1_yaml_sweep.py:99–104` pattern). The module has three nodes:
- **`test_u76_palace_refusals_by_field`** (W0 real, and again in RG-parse complex).
  - **Rows.** Each row is one file run through `_run_study_file_palace(f, output_dir=scratch/…, comm)` and through `cli.main(["run", f, "--output-dir", scratch/…, "--backend", "palace"])`:
    - `f_small_four_port.yaml` with exactly one added key: `outputs.sar`, `outputs.b1_plus` (both from `f_small_sar_b1.yaml`'s mapping), `outputs.probes` (from `f_small_em_probes.yaml`) or `tune: {}`;
    - the tracked `f_small_em_probes.yaml`, `f_small_sar_b1.yaml`, `f_small_sweep.yaml`, `f_small_copper.yaml`, `f_small_hole_pec_capacitors.yaml`, `f_small_hole_pec_tune.yaml`, `sphere_mat4.yaml`, `f_human_tuned_copper_thermal.yaml`, `thermal_box_steady_m4.yaml` and `thermal_ball_wall.yaml`;
    - `malformed_thermal_on_port.yaml`.
  - **(X) by field, rc 2, nothing written.** The field is `--backend` on every row except `malformed_thermal_on_port.yaml`, whose field is `thermal` (the parse refusal, unchanged). The single-key rows' messages name `outputs.sar` / `outputs.b1_plus` / `outputs.probes` / `tune`; the sweep row's names `study.sweep_hz`. `cli.main` returns 2 on every rank (`_rc_all`, `tests/package/test_wf1_study_yaml.py:169`). Rank 0's stderr carries `femem: error:` and the key. The output directory does not exist after the call, on every rank.
  - **(Th) the structural guard:** `_run_palace(…, cooling=object())` on the plain four-port Study raises with field `--backend` and message `palace: thermal is not supported`, with 0 builds.
- **`test_u76_palace_admits_and_carries`** (W0 only, `-n 1`). `fem_em_solver.study.Study.solve` is monkeypatched to a spy. It returns a stand-in result (`s_matrix` = `np.eye(4)`, `cells` 0, `order` 1, `case_dir`, `sha256` `{}`, `palace_total_s` `None`) only when called with `backend="palace"`; otherwise it raises.
  - **(A) admitted:** `f_small_four_port.yaml` with no outputs, and the same with only `outputs.directory`, each reach the spy once. Each writes `summary.json` whose EM-output keys (∩ {`s_matrix`, `z_matrix`, `touchstone`, `sar`, `b1_plus`, `probes`, `tune`, `temperature`, `s_matrices`}) are exactly `{"s_matrix"}`.
  - **(Pg) the page says it:** the new ```` ```text palace-outputs ```` block in `docs/guide/backends.md` has three lines:
    - `honoured: outputs.directory`, which must equal `["outputs." + k for k in _PALACE_ROUTED_OUTPUTS]`;
    - `returns: s_matrix`, which must equal (A)'s key set;
    - `refused: …`, which must equal, as a set, the key names that (X) / (Th) messages carry. Extract them with `palace: (?:a )?(\S+) is not supported` over every refusal raised in this module, so `outputs.b1_plus outputs.probes outputs.sar study.sweep_hz thermal tune conductor capacitor_sheets` is predicted.

    The block is read with `extract_tagged_block` (`tests/package/test_wf17_first_study_tutorial.py:78`).
- **`test_u76_u16_readings_bitwise`** (knob `FEM_EM_U76_BIT=1`; unset ⇒ skip; `"1"` ⇒ run with absent inputs **failing**; anything else ⇒ `ValueError`, the item-576 knob pattern).
  - **(Bit)** It reads two logs from `docs/testing/logs/` (`.log` or `.log.gz`):
    - the comparand `20261018T172228Z_ARCH-1-U16c-W1.log`;
    - the newest by name of `*_ARCH-1-U76-RGU16*.log`.

    In each log, the lines starting `[U16c T4 E2]`, `[U16c T4 Tie2]`, `[U16c T4 c-order]` and `[U16c T5 E1]` each occur exactly once. Each must be **string-equal** across the two logs: comparand `:12548`, `:12550`, `:12551`, `:14336`.
  - This is "U16's (E2) / (Tie2) bitwise" as an assertion on the prints of one function, one batch, one call. The lines print `.3e` / `.4f` of exact equalities, and (E2) is 0.000e+00, which a non-zero S difference cannot print.

**Page (`docs/guide/backends.md`).**
- **New section `## What the Palace backend carries`,** after "What is refused, and why" (`:57–75`). It holds the `palace-outputs` block and one sentence: `femem run --backend palace` refuses every other output by name, with exit code 2, before anything is written.
- **"What is refused" gains one bullet** for the YAML refusals. It is prose: name the keys, no number.
- **`:12–13` and the "Not claimed" line (`:121–123`) stay.**
- **The page must cite no new `tests/…::node`.** The U62 test's (L) asserts exactly six cited nodes (`tests/package/test_wf17_backends_page.py:124`). Name the new module by its path only (no `::`), e.g. "`tests/package/test_arch1_u76_palace_outputs.py`".
- Write `≈`, never a tilde. No AED string, no Palace number.

**Registry (lines added only).** One comment line `# U76 Palace output refusals (§10, weekly 2026-10-21; §9 item 590)`, plus `tests/package/test_arch1_u76_palace_outputs.py::test_u76_palace_refusals_by_field` under `[names.StudyConfigError]` (`tests/package/validation_envelope.toml:338`). It has no number in the comment. **The §9 standing rule binds** (`PROJECT_PLAN.md:4087`, `:3965`): C-api / C-val, then Wgen / Wval / Wapi / Wdoc in this slot, and the regenerated pages are committed. Wval: read the current node count first and predict +1. It is printed only.

**Anchors (asserted).**
- **(X)** rc 2 by field, nothing written, on 15 rows (above).
- **(Th)** the structural `thermal` guard.
- **(A)** the admitted path reaches the solve once and carries `s_matrix` only.
- **(Pg)** the page's carried / refused lines equal the code's.
- **(Bit)** U16's four reading lines are string-equal to U16c's.
- **RG-U16:** the U16 module as committed, 14 passed, none skipped (U16c's count). This includes (E2) ≤ `LEG_D0_REPRODUCTION_BAND` (`tests/validation/test_port_birdcage_four_port.py:165`), (Tie2) ≤ `TIE2_BAND` 2.6e-10 (the module constant, 2026-10-18 03:00 R3), (c-order) > 1e-6 and (E1) ≤ `LEG_D0_REPRODUCTION_BAND`.
- **RG-parse:** item 572's 21-node bundle, 21 passed, with the new module's refusal node added (22 predicted). (A) is not run at `-n 2`: `_run_palace` writes `summary.json` on every rank (`:1815–1818`, no rank-0 guard; the route is serial by design).
- **RG-U62:** the backends page test as committed, 9 passed. (L) still reads 6 nodes, and (D) reads 0 dead.

**Negative controls (rule (e)).**
- **Asserted, by construction, backed by the read-only finding (`WF-1.md:3021`) — (c-drop).** With `_PALACE_ROUTED_OUTPUTS` monkeypatched to `("directory", "probes")`, `_run_study_file_palace(f_small_em_probes.yaml)` reaches the spy once. Its summary has no `probes` key. That is today's silent drop, reproduced on the same file, and it shows the allow-list is the only thing that refuses. Placed in `test_u76_palace_admits_and_carries`.
- **Asserted, by construction — (c-pg).** A tmp copy of the page with `outputs.probes` deleted from the `refused:` line fails (Pg)'s comparator, and the comparator names `outputs.probes`.
- **Asserted, by construction — (c-bit).** (Bit)'s comparator on a copy of the U16c `Tie2` line with one digit changed returns unequal. It runs in the knob node before the real comparison.
- **Predicted, printed:** none needed. (c-order) inside RG-U16 is U16c's own asserted control (6.6674 / 2.5442 / 4.5224 %, `20261018T172228Z_ARCH-1-U16c-W1.log:12551`).

**Predictions (printed, never asserted).**
- (X): all 15 rows rc 2 and 0 builds; W0 ≈ 3–8 s.
- RG-U16: 14 passed, ≈ 220–230 s. (E2) 0.000e+00, (Tie2) 2.464e-10 / 9.668e-11 / 9.883e-11, (E1) 2.149e-14 (entry max 3.893e-14), `palace_total_s` 261.619549041, identical to U16c.
- Wval +1 node.
- The `femem_out` directory is never created: every call takes an explicit `--output-dir` under the scratch.

**Windows.** Foreground; **the harness call is the whole Bash command: bare, repo-relative, nothing chained or piped**; Bash `timeout: 660000`; pytest `-s`; node IDs, never `-k`. Run `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` before RG-U16; it must read 0. **Order:** code + module + page + registry, W0, RG-parse, RG-U16, W-bit, RG-U62, C-api, C-val, Wgen, Wval, Wapi, Wdoc.

| Window | Price | Source |
|---|---|---|
| W0 (real, `-n 1`) | ≈ 5 s | item 572 W0-r2 4 s (`20261020T203852Z_WF-1-s21-W0-r2.log:469`) |
| RG-parse (complex, `-n 2`) | ≈ 5 s | the same log |
| RG-U16 (complex, `-n 1`) | 224 s | `20261018T172228Z_ARCH-1-U16c-W1.log:14366` |
| W-bit | ≈ 3 s | log parsing only |
| RG-U62 | 5 s | `20261019T005151Z_ARCH-1-U62-W0.log:76` |
| C-api / C-val / Wgen / Wval | 2 / 2 / 2 / 3 s | `20261019T005017Z` … `005044Z_ARCH-1-U62-*.log` |
| Wapi / Wdoc | 21 / 2 s | `20261021T003644Z_WF-17-s8-Wapi.log:228`, `20261019T005254Z_ARCH-1-U62-Wdoc.log` |

```
scripts/testing/run_and_log.sh ARCH-1-U76-W0-real "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_arch1_u76_palace_outputs.py::test_u76_palace_refusals_by_field tests/package/test_arch1_u76_palace_outputs.py::test_u76_palace_admits_and_carries -s -v --tb=short'"
scripts/testing/run_and_log.sh ARCH-1-U76-RGparse "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_arch1_u76_palace_outputs.py::test_u76_palace_refusals_by_field tests/package/test_wf1_yaml_sweep.py::test_wf1_s21_parse_and_refusals tests/package/test_wf1_yaml_em_probes.py::test_wf1_s20_parse_and_refusals tests/package/test_wf1_yaml_thermal_ball_probes.py::test_wf1_s19_parse_and_refusals tests/package/test_wf1_yaml_thermal_box_quadratic_duty.py::test_wf1_s18_parse_and_refusals tests/package/test_wf1_yaml_thermal_box_probes.py::test_wf1_s17_parse_and_refusals tests/package/test_wf1_yaml_thermal_box.py::test_wf1_s15_parse_and_refusals tests/package/test_wf1_yaml_thermal_box_quadratic.py::test_wf1_s16_parse_and_refusals tests/package/test_wf1_yaml_thermal_ball.py::test_wf1_s14_parse_and_refusals tests/package/test_wf1_yaml_capacitors.py::test_wf1_s13_parse_and_refusals tests/package/test_wf1_yaml_tune.py::test_wf1_s11_tune_parse_and_refusals tests/package/test_wf1_yaml_tune_f_human.py::test_wf1_s12_parse_and_refusals tests/package/test_wf1_check_family.py tests/package/test_wf1_yaml_f_human.py tests/package/test_wf1_yaml_f_human_thermal.py -v -s --tb=short'"
scripts/testing/run_and_log.sh ARCH-1-U76-RGU16 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_U16_PALACE_READ=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 590 mpiexec -n 1 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_arch1_palace_backend.py -v -s --tb=short'"
scripts/testing/run_and_log.sh ARCH-1-U76-Wbit "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace FEM_EM_U76_BIT=1 timeout -k 30 60 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_arch1_u76_palace_outputs.py::test_u76_u16_readings_bitwise -s -v --tb=short'"
scripts/testing/run_and_log.sh ARCH-1-U76-RGU62 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace timeout -k 30 120 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf17_backends_page.py -s -v --tb=short'"
```
- RGparse is item 572's W0-r2 command (`20261020T203852Z_WF-1-s21-W0-r2.log:12`) with the new refusal node prepended. RGU16 is U16c W1's command (`20261018T172228Z_ARCH-1-U16c-W1.log:12`) verbatim, relabelled. RGU62 is item 529's W0 command verbatim, relabelled.
- **C-api** (predicted Status 1, differing `['StudyConfigError.md']`), **C-val** (predicted Status 1, `DIFFERS: validation.md`), **Wgen**, **Wval**, **Wapi** and **Wdoc** are item 529's six commands verbatim (`ARCH-1.md` "§9 item 529", windows 2–5, 7, 8), labelled `ARCH-1-U76-ctl-api` / `-ctl-val` / `-Wgen` / `-Wval` / `-Wapi` / `-Wdoc`. Wdoc reads `dead=0`.
- No durable capture: every window is well under ≈ 10 min. **No module, page or `src` edit after W0 starts** (rule (i)); an edit restarts at W0.

**Traps already paid for.**
- **The U62 page test counts six cited nodes** (`test_wf17_backends_page.py:124`). A `tests/…::…` citation added to the page turns RG-U62 red. Cite the module path only. Never edit that test.
- **Keep the field `--backend`.** Item 572's assertion reads it, and a field change is an existing-assertion edit by proxy.
- **Never write into `logs/u16/*`** except T5's `summary.json`, which the committed module rewrites. Never touch a `postpro/`, never re-run Palace, and never write into `logs/ops61/`. `write_case`'s own check rewrites an identical `case.sha256`; that is the package's read path (item 533's trap).
- **`-n 1` for RG-U16:** the backend refuses a parallel comm by design. RG-parse is `-n 2` because item 572's bundle ran there.
- **A refusal must precede the mesh:** the build stub's 0-call assertion catches a refusal that moved below `study.solve`.
- **`femem_out` is not gitignored:** pass an explicit `--output-dir` under `logs/arch1/u76/` on every call; `git status --short` must read clean before the commit.
- **`capsys.readouterr()` discards earlier prints** (items 507 / 562): print the table after any capture.
- **Leak guard:** if a log is flagged, move it unchanged to `logs/leak-held/`, mark the item BLOCKED and stop.
- **Harness logs may be gzipped after rotation:** (Bit) reads `.log` or `.log.gz`.
- Never `-k`; `-k 30` on every container `timeout`; no pipe; never delete a log; never import a `test_*` name; no AED string; write ≈, never a tilde.
- **Status 124 or a kill:** run `pgrep` first. At 0 only, run the rubric's orphan-stem sweep, labelled `ARCH-1-U76-cache-sweep`.

**Scope.** Covered: `femem run --backend palace` refuses, by field with rc 2 and before any mesh, every `outputs` key but `directory`, plus `tune`, `study.sweep_hz` and a structural `thermal`, with conductor / capacitor sheets / other families as before; the page states the carried / refused split from the code; U16's gates re-run with their reading lines byte-identical. **Not claimed:**
- any new Palace output (SAR, `|B₁⁺|`, probes or fields from Palace are `ARCH-1` (1)–(2));
- refusing `geometry.as_hole` / degree ∉ {1, 2} / a parallel comm through the YAML route with rc 2 (For the review, 1);
- `Study.solve(backend="palace")`'s own refusals, which are unchanged;
- any Palace run, accuracy figure or AED statement.

**Negative result.** Report and stop. Park on `attempt/ARCH-1-U76-*` and mark item 590 BLOCKED in the same commit (rule (d)).
- **(X) red on a pre-existing row** (sar / b1_plus / tune / sweep / conductor / capacitor / family): the allow-list changed an existing outcome. Revert `src/` and park.
- **(A) red:** the allow-list refuses an admissible file, and U16's T3 / T5 would go red too. Fix once in-letter; otherwise park.
- **RG-U16 red, or (Bit) unequal with RG-U16 green:** the Palace read path moved, though this item does not touch it. Open a known-issues entry naming both logs, BLOCKED. Never re-run Palace to "refresh" a comparand.
- **RG-parse red on a pre-existing node:** a parse regression. Revert and park.
- **RG-U62 red:** a page defect. Fix the page in-letter once; the clock restarts at W0.
- **C-api / C-val exit 0:** the registry line did not reach the pages. Stop.

No `PREDICTED RED:` line, so a red anywhere stops the slot.

**Status it can move.** **No glyph:** `ARCH-1` stays 🟡. Under rule (j), the row's Done-when, "digit-identical reproduction of every imported band and record on the DolfinX backend before any second backend is wired" (`PROJECT_PLAN.md:1388`), is neither met nor touched by this item. The row gains: "Step 1d (§9 item 590, U76): `femem run --backend palace` refuses every output the Palace route has no path for — `outputs.sar` / `b1_plus` / `probes`, `tune`, `study.sweep_hz`, `thermal` — by field with rc 2, before any mesh; `docs/guide/backends.md` states the carried / refused split from the code; U16's (E2) / (Tie2) / (c-order) / (E1) lines re-read byte-identical (<RGU16 log>, <Wbit log>)." §10 U76 → DONE. The weekly's silent-drop finding (`PROJECT_PLAN.md:4819–4823`) is closed for its `_run_palace` half; U75 holds the other half. No known-issues entry moves, because none was opened (the weekly posed U75 / U76 instead).

**For the review.**
1. **Non-output refusals still escape as tracebacks.** `femem run --backend palace tests/package/wf1/f_small_hole_pec_control.yaml` passes `_run_palace` and raises `ValueError("backend='palace': as_hole=True is not supported")` inside `solve_palace` (`src/fem_em_solver/backends/palace.py:217–218`), uncaught by `cli.main`. So it exits by traceback, not rc 2. A degree-3 file and `mpiexec -n 2 femem run --backend palace …` do the same (`:231–233`). These are geometry and width, not outputs, so they are outside U76's letter. The options:
   - (a) leave them, print the hole-control row's outcome (`ValueError`), and assert nothing;
   - (b) add `geometry.as_hole` and degree refusals to `_run_palace` with `--backend` (two lines, the same gate).

   The drafter's default is (a); (b) is the review's to add.
2. **Re-gate breadth.** The brief re-runs the Palace-route gates plus item 572's 21-node parse bundle, because `run_study_file` never calls `_run_palace`. The 2026-10-20 10:30 ruling (1) (`WF-1.md:3018`) re-ran every YAML gate for item 572, which edited the parse path. If the review applies it here, add item 572's RG list (≈ 17 min, RG4 at `-n 8`). Cost is then ≈ 55 predicted minutes.
3. **Registration.** `StudyConfigError` only (one line) is the drafter's choice. `run_study_file` would misdescribe a private route. Unregistered (the `WF-17` tutorial precedent) would drop C-api … Wval.
4. **The step's "real, `-n 1`"** fits W0. RG-U16 is complex and heavy by ceiling (224 s against 180 s), so it needs the `pgrep` pre-check, and no other heavy window may run beside it.

**Execution order.** Host listing of the eight on-disk inputs → code → module → page → registry → W0 → RG-parse → `pgrep` → RG-U16 → W-bit → RG-U62 → C-api → C-val → Wgen → Wval → Wapi → Wdoc → commit (`study_yaml.py`, the module, the page, the registry, the regenerated `docs/guide/api/` and `docs/guide/validation.md`, every `ARCH-1-U76-*` log and its test-results rows, the §7 sentence, §9 item 590 DONE, §10 U76) → journal.

**Executor spawn prompt (verbatim).** "Execute §9 item 590 — `ARCH-1` step 1d (§10 U76) — per `docs/planning/chunks/ARCH-1.md` "§9 item 590". Read the brief and any review ruling beneath it first; the rulings hold. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. List `logs/u16/f_small_d2_64/postpro/port-S.csv`, `logs/u16/f_small_d1_64_cli/palace_case/postpro/port-S.csv` and `logs/ops61/postpro_f05x_o2/`, `postpro_f05x/`, `postpro_f05x_o2_f10/` first; if any is absent, stop and mark item 590 BLOCKED; never run Palace and never write a `postpro/`. `src/` edit: `study_yaml.py` only — the two private constants, the allow-list loop and `cooling=None` guard in `_run_palace`, `cooling` passed by `_run_study_file_palace`; keep every refusal's field `--backend`; do not touch `run_study_file`, the module docstring, `cli.py`, `study.py` or `backends/palace.py`. Create `tests/package/test_arch1_u76_palace_outputs.py`; edit `docs/guide/backends.md` without citing any new `tests/…::…` node; add one registry line under `[names.StudyConfigError]`. Never edit `test_arch1_palace_backend.py`, `test_wf17_backends_page.py` or any existing test. Run W0, RG-parse, RG-U16 (after `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reads 0), W-bit, RG-U62, C-api, C-val, Wgen, Wval, Wapi, Wdoc in order, in the foreground, with `timeout: 660000` on every Bash call that runs a window, pytest with -s, container `timeout -k 30`. Pass an explicit output directory under `logs/arch1/u76/` on every `femem run` call; `git status --short` clean before the commit. You are the executor — do not spawn agents, never background a harness window. If the docker socket is denied, stop and report the exact command string. On a red, park on `attempt/ARCH-1-U76-*`, mark item 590 BLOCKED in the same commit, and report."

### Review ruling 2026-10-22 03:00 (item 590, U76)

Binds over the brief. Queued as drafted. (1) Non-output refusals (hole mesh, degree 3, `-n 2` tracebacks) stay out of scope, named in the outcome. (2) Re-gate breadth: the Palace-route gates plus item 572's parse bundle, not every YAML gate. (3) The refusal node is registered under `StudyConfigError` only. No Palace window is run.

### §9 item 601 — 2026-10-22 (§10 U79, chain step (no solve); chain U, 9e — daily, 2026-10-22 — provisional)

*Written by the 2026-10-22 18:00 daily review under §9 rule (6), from item 590's "For the review" 1 (option (b)) and its outcome's out-of-scope line. The weekly ratifies or reverts.*

**Licence.** §9 rule (6): the next numbered step of live chain U toward subgoal (9e) and the Phase 9 exit gate. The target is U76's own, applied to the two inputs U76 left out: a usable package refuses an unsupported study by a named field, never by a traceback. **No band, tolerance, record, Done-when or validation target is set or moved.** This is `ARCH-1` **step 1e**. Every U76 / U16 assertion stays as committed.

**What the code is today (read 2026-10-22, `8f89f85a`).**
- `_run_palace` (`src/fem_em_solver/study_yaml.py:1790–1806`) refuses, in this order: family (`:1793`), `cooling` (`:1795`), every `outputs` key but `directory` (`:1797–1801`), a conductor (`:1802`), `CapacitorSheets` (`:1804`). It then calls `study.solve(comm, backend="palace", …)` (`:1806`).
- `solve_palace` (`src/fem_em_solver/backends/palace.py:217–233`) raises **`ValueError`** for `geo.as_hole` (`:218`), `comm.size != 1` (`:231`) and a degree outside {1, 2} (`:233`). `cli.main` maps only `StudyConfigError` to rc 2. So on the YAML route, two inputs exit by traceback: `femem run --backend palace tests/package/wf1/f_small_hole_pec_control.yaml` (as_hole, PEC hole wall, no conductor, degree 1), and any file under `mpiexec -n 2 femem run --backend palace …`.
- **Degree is already refused at parse** (`study_yaml.py:275–278`: `study.degree` must be 1 or 2, field `study.degree`). Item 590's "For the review" 1 listed "a degree-3 file" as a third traceback. **That reading is wrong for the YAML route:** the file never reaches `_run_palace`. This item records the parse refusal and adds no degree check.

**`ARCH-1` step 1e / U79: `femem run --backend palace` refuses a hole-wall birdcage (`geometry.as_hole`) and a launch wider than one rank. Each refusal has field `--backend` and rc 2, and comes before any mesh and before anything is written. A degree-3 file is shown refused at parse (field `study.degree`, rc 2), unchanged.**
- **Executor:** implementer.
- **Tier:**
  - W0: standard, real build, `-n 1`, no mesh.
  - W0b: smoke, complex, `-n 2`: the width row plus U76's refusal node and item 572's parse bundle, i.e. item 590's RG-parse command with the new nodes prepended.
  - RG-U16: **heavy by ceiling**, complex, `-n 1`, 220 s measured (`20261022T171719Z_ARCH-1-U76-RGU16.log`).
  - RG-U62, C-api, C-val, Wgen, Wval, Wapi, Wdoc: smoke.
  - **No Palace window.**
- **Class:** chain step (no solve). **Cost:** ≈ 11.2 measured / ≈ 35 predicted slot-minutes. Item 590 landed in 11.2 on the same window list.
- **Depends on:** nothing in §9. It needs the same on-disk inputs as item 590 (`logs/u16/…`, `logs/ops61/postpro_*`). If any is absent ⇒ BLOCKED ("comparand absent"); never re-run Palace.
- **Family:** `ARCH-1` step 1e; attempt 1 of 4. **Daytime:** no.

**Code (`src/fem_em_solver/study_yaml.py` only; disclosed under rule (c)).**
- **New private helper** `_palace_input_refusal(study, path, comm)` beside `_run_palace`. It raises `StudyConfigError(path, "--backend", …)`:
  - `if study.geometry.as_hole:` → `palace: geometry.as_hole is not supported`;
  - `if comm.size != 1:` → `palace: a {comm.size}-rank launch is not supported (the case writer is serial; run at mpiexec -n 1)`.
- **It is called from `_run_palace` right after the `CapacitorSheets` line (`:1804–1805`) and before `study.solve` (`:1806`).** Placed after every existing refusal, it leaves every U76 row's outcome and message unchanged. In particular, U76's `f_small_hole_pec_capacitors.yaml` row is an as_hole file refused first by `capacitor_sheets`. It must stay that way, because U76's (Pg) set depends on it.
- Nothing goes into `__all__`. These stay byte-unchanged: `backends/palace.py`, `cli.py`, `study.py`, `run_study_file`, the module docstring and every existing test. `Study.solve(backend="palace")`'s own `ValueError`s (U16's `test_u16_backend_refusals`, `tests/package/test_arch1_palace_backend.py:117–156`) are not touched.

**New `tests/package/test_arch1_u79_palace_inputs.py`.**
- Docstring paragraph 1 names `ARCH-1` and `WF-1`. The module imports cleanly in both builds and builds no mesh.
- It reuses helpers by name (never a `test_*` name): U76's `_field`, `_all`, `_stub_builds`, `_write` (`tests/package/test_arch1_u76_palace_outputs.py:73`, `:82`, `:86`, `:99`), and `WF1`, `_fresh_dir`, `_mapping`, `_rc_all` (`tests/package/test_wf1_study_yaml.py`).
- Scratch is gitignored `logs/arch1/u79/`: rank 0 writes after `_fresh_dir`, then `Barrier`. Every `--output-dir` lives under it.

Nodes:
- **`test_u79_hole_refused_by_field`** (W0, `-n 1`, real; again in W0b at `-n 2`). It runs `f_small_hole_pec_control.yaml` through `_run_study_file_palace` and through `cli.main(["run", f, "--output-dir", …, "--backend", "palace"])`. Asserted: field `--backend`, message `palace: geometry.as_hole is not supported`, rc 2 on every rank, 0 builds, output directory absent. At `-n 2` the as_hole check fires first (it precedes the width check), so the message is the same at both widths.
- **`test_u79_width_refused_by_field`** (W0b only, `-n 2`). At `-n 1` it asserts `comm.size == 1` and skips with the reason printed.
  - `f_small_four_port.yaml` with no outputs: field `--backend`, message starting `palace: a 2-rank launch is not supported`, rc 2 on both ranks, 0 builds, nothing written.
  - The order row: the same file with `outputs.sar` added is still refused by `outputs.sar`, because the allow-list precedes the width check.
- **`test_u79_degree_refused_at_parse`** (W0, `-n 1`). `f_small_four_port.yaml` with `study.degree: 3`, written from `_mapping`. Asserted: `cli.main([..., "--backend", "palace"])` returns rc 2, and `_load` raises `StudyConfigError` with field `study.degree`. This records the existing parse refusal (`study_yaml.py:277–278`); a spy on `_run_palace` reads 0 calls.
- **`test_u79_admitted_unchanged`** (W0, `-n 1`). `Study.solve` is spied with U76's (A) stand-in result (`test_arch1_u76_palace_outputs.py:227–296`). `f_small_four_port.yaml` reaches the spy once, and its summary's EM keys are `{"s_matrix"}`. So the new helper refuses nothing admissible at `-n 1`.

**Anchors (asserted).**
- (H) the hole file is refused by field, rc 2, at both widths.
- (W) a 2-rank launch is refused by field, rc 2, with the order row.
- (D) degree 3 is refused at parse, with 0 `_run_palace` calls.
- (A′) the admissible file reaches the solve once.
- **Re-gates:**
  - U76's module as committed: its two W0 nodes at `-n 1` (3 passed / 3 skipped, as item 590's W0), and its refusal node inside W0b at `-n 2`.
  - Item 572's parse bundle inside W0b, every node passed.
  - RG-U16: 14 passed, with its four reading lines printed. Item 590's W-bit is **not** re-run: U16's read path is not edited, and RG-U16's own (E2) / (Tie2) asserts stand.
  - RG-U62: 9 passed, (L) 6 nodes, (D) 0 dead.

**Negative controls (rule (e)).**
- **(c-helper), asserted, by construction.** With `study_yaml._palace_input_refusal` monkeypatched to a no-op and `Study.solve` spied, `f_small_hole_pec_control.yaml` reaches the spy once at `-n 1`. That is today's route to `solve_palace`, where `:218`'s `ValueError` fires, and it shows the helper is the only new refusal. Backing: the read of `palace.py:217–218` above and item 590's outcome line ("hole-mesh … still exit by `ValueError` traceback in `solve_palace`").
- **(c-old), printed.** Run `_run_study_file_palace(f_small_hole_pec_control.yaml)` with the helper no-op'd, `Study.solve` **not** spied, and `FourPortBirdcageGeometry.build` stubbed to raise. Print the exception type and message. Predicted: `ValueError("backend='palace': as_hole=True is not supported")`, raised before `build`, because `palace.py:217–218` precedes `:236`. Never asserted.

**Page (`docs/guide/backends.md`).**
- Add one prose bullet to "What is refused, and why" (`:57–77`), after the "A study file asking for an output …" bullet. It says that on `femem run --backend palace`, a hole-wall birdcage (`geometry.as_hole`) and a launch on more than one rank are refused by name with exit code 2, before anything is written; an element degree other than 1 or 2 is refused when the file is read.
- **Cite no new `tests/…::node`:** U62's (L) counts six cited nodes. Name the module by path only.
- **Do not edit the `palace-outputs` block.** U76's (Pg) computes its `refused:` set from U76's own rows, and the new refusals are inputs, not outputs.

**Registry (lines added only).**
- Add one comment line, `# U79 Palace input refusals (§10 chain U, daily 2026-10-22 provisional; §9 item 601)`, and `tests/package/test_arch1_u79_palace_inputs.py::test_u79_hole_refused_by_field` under `[names.StudyConfigError]`, beside item 590's U76 line.
- The §9 standing rule binds. Run C-api / C-val (predicted Status 1), then Wgen, Wval, Wapi and Wdoc in this slot, and commit the regenerated pages. For Wval, read the node count first and predict +1 (printed).
- If another landed item regenerated the pages first: rebase, keep both sides, and re-run Wgen / Wval / Wapi / Wdoc (§9 digest).

**Windows.**
- Foreground. **The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped.** Bash `timeout: 660000`; pytest `-s`; node IDs, never `-k`.
- `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0 before RG-U16.
- **Order:** code + module + page + registry → W0 → W0b → `pgrep` → RG-U16 → RG-U62 → C-api → C-val → Wgen → Wval → Wapi → Wdoc → commit → journal.
```
scripts/testing/run_and_log.sh ARCH-1-U79-W0-real "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_arch1_u79_palace_inputs.py tests/package/test_arch1_u76_palace_outputs.py::test_u76_palace_refusals_by_field tests/package/test_arch1_u76_palace_outputs.py::test_u76_palace_admits_and_carries -s -v --tb=short'"
```
- **W0b** is item 590's RG-parse command (`20261022T171704Z_ARCH-1-U76-RGparse.log:12`) verbatim, with `tests/package/test_arch1_u79_palace_inputs.py::test_u79_hole_refused_by_field tests/package/test_arch1_u79_palace_inputs.py::test_u79_width_refused_by_field` inserted after `tests/environment/test_complex_mode.py`. Label `ARCH-1-U79-W0b-n2`.
- **RG-U16** is item 590's RG-U16 command (`20261022T171719Z_ARCH-1-U76-RGU16.log:12`) verbatim, labelled `ARCH-1-U79-RGU16`.
- **RG-U62, C-api, C-val, Wgen, Wval, Wapi, Wdoc** are item 590's landed commands (`20261022T172119Z_ARCH-1-U76-RGU62.log:12` and the `ARCH-1-U76-ctl-api` … `-Wdoc` logs, `:12`), relabelled `ARCH-1-U79-…`.
- No durable capture: every window is well under ≈ 10 min. **No module, page or `src` edit after W0 starts** (rule (i)).

**Traps already paid for.**
- Item 590's list binds verbatim (this file, "§9 item 590", "Traps already paid for"): six cited nodes on the page; field `--backend`; never write into `logs/u16/*` except T5's `summary.json`; `-n 1` for RG-U16; a refusal precedes the mesh; `femem_out` is not gitignored; `capsys.readouterr()` discards earlier prints; the leak guard; gzipped logs; never `-k`; `-k 30` on every container `timeout`; no pipe; never delete a log; never import a `test_*` name; no AED string; `≈`, never a tilde.
- **New here:** at `-n 2` both ranks must raise before any collective. The helper reads only `comm.size` and the geometry, so no rank waits in a `Barrier`. A hang at W0b means a collective slipped in front of the refusal: kill at the ceiling, run `pgrep`, report.

**Scope.** Covered: the YAML Palace route's two remaining traceback exits become named refusals. **Not claimed:**
- any Palace run or output;
- `Study.solve(backend="palace")`'s own `ValueError`s, which are unchanged by design (the Python API raises; the CLI refuses);
- the `--ranks` / `--estimate` combinations, which are refused already (`cli.py:126–131`).

**Negative result.** Report and stop. Park on `attempt/ARCH-1-U79-*` and mark item 601 BLOCKED in the same commit (rule (d)).
- **(H) or (W) red:** the helper's placement or message is wrong. Fix once in-letter and restart at W0.
- **A U76 node red:** the helper moved an existing outcome. Revert `src/` and park.
- **RG-U16 red:** open a known-issues entry naming the log, BLOCKED. This item does not touch the read path.
- **W0b hang:** as in the trap above.

No `PREDICTED RED:` line, so a red anywhere stops the slot.

**Status it can move.** No glyph: `ARCH-1` stays 🟡, because rule (j) leaves its Done-when untouched. The row gains: "step 1e (§9 item 601, U79): `femem run --backend palace` refuses `geometry.as_hole` and a launch wider than one rank by field `--backend`, rc 2, before any mesh; a degree other than 1 or 2 is refused at parse (`study.degree`)". §10 U79 → DONE. Item 590's out-of-scope line closes.

**Executor spawn prompt (verbatim).** "Execute §9 item 601 — `ARCH-1` step 1e (§10 U79) — per `docs/planning/chunks/ARCH-1.md` "§9 item 601". Read the brief and any review ruling beneath it first; the rulings hold. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. List `logs/u16/f_small_d2_64/postpro/port-S.csv`, `logs/u16/f_small_d1_64_cli/palace_case/postpro/port-S.csv` and `logs/ops61/postpro_f05x_o2/`, `postpro_f05x/`, `postpro_f05x_o2_f10/` first; if any is absent, stop and mark item 601 BLOCKED; never run Palace. `src/` edit: `study_yaml.py` only — the `_palace_input_refusal` helper, called after the `CapacitorSheets` line and before `study.solve`. Create `tests/package/test_arch1_u79_palace_inputs.py`; add one bullet to `docs/guide/backends.md` without citing any `tests/…::…` node and without touching the `palace-outputs` block; add one registry line under `[names.StudyConfigError]`. Never edit an existing test. Run W0, W0b, RG-U16 (after `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reads 0), RG-U62, C-api, C-val, Wgen, Wval, Wapi, Wdoc in order, in the foreground, with `timeout: 660000` on every Bash call that runs a window, pytest with -s, container `timeout -k 30`. You are the executor — do not spawn agents, never background a harness window. If the docker socket is denied, stop and report the exact command string. On a red, park on `attempt/ARCH-1-U79-*`, mark item 601 BLOCKED in the same commit, and report."

**Outcome of item 601 (2026-10-23, 00:30 UTC slot): DONE.** `_palace_input_refusal(study, path, comm)` in `study_yaml.py`, called from `_run_palace` after the `CapacitorSheets` line and before `study.solve`; field `--backend` on both refusals; nothing added to `__all__`. W0 real `-n 1` 6 passed / 4 skipped, 5 s (`20261023T003505Z_ARCH-1-U79-W0-real.log`): (H) `f_small_hole_pec_control.yaml` field `--backend`, msg `palace: geometry.as_hole is not supported`, rc (2, 2), nothing written, 0 builds; (D) degree 3: `_load` field `study.degree` (`must be 1 or 2, got 3`), CLI rc (2, 2), 0 `_run_palace` calls; (A′) `f_small_four_port.yaml` 1 solve, EM keys `{s_matrix}`; (c-helper) helper no-op'd → the hole file reaches the spy once; (c-old), printed: `ValueError("backend='palace': as_hole=True is not supported")`, build calls 0, as predicted; U76's two nodes 3 passed / 3 skipped as item 590's W0 ((X) 15 rows, 0 builds). W0b complex `-n 2` 24 passed, 5 s (`20261023T003520Z_ARCH-1-U79-W0b-n2.log`): (H) same message and rc (2, 2) at `-n 2`; (W) `palace: a 2-rank launch is not supported (the case writer is serial; run at mpiexec -n 1)`, rc (2, 2), nothing written; order row `+outputs.sar` still refused by `outputs.sar`; U76 (X) 15 rows, 0 builds; item 572's parse bundle all passed; no hang. `pgrep -c python3` 0, then RG-U16 complex `-n 1` 14 passed, 226 s (`20261023T003534Z_ARCH-1-U79-RGU16.log`): (E2) 0.000e+00, (Tie2) 2.464e-10 / 9.668e-11 / 9.883e-11, (c-order) 6.6674 / 2.5442 / 4.5224 %, (E1) 2.149e-14. RG-U62 9 passed, (L) 6 nodes, (D) 0 dead (`20261023T003926Z_ARCH-1-U79-RGU62.log`). C-api Status 1 `['StudyConfigError.md']`, C-val Status 1 `DIFFERS: validation.md` (both as predicted); Wgen 46 files; Wval nodes 340 → 341 (+1 as predicted; the 340 is item 589's `20261022T201637Z_WF-1-s23-Wval.log`, which landed after item 590's 334); Wapi 22 passed, 23 s; Wdoc `dead=0` (1369 links). No Palace run, no `postpro/` written. Item 590's out-of-scope line closes: on the YAML route, the hole mesh and an `-n 2` launch are now rc 2 by name, and a degree-3 file was already refused at parse. One deviation: the module's scratch YAML files are written by U76's reused `_write`, so they land under the gitignored `logs/arch1/u76/yaml/` (`u79_*.yaml`); every `--output-dir` is under `logs/arch1/u79/`.

## §7 row moved verbatim 2026-10-24 (rotate_plan_archive.py chunks)

 The MRI-safety layer becomes solver-agnostic | 🟡 step 1 (U16, §9 items 517 / 517b / 533): `Study.solve(backend="palace")` and `femem run --backend palace` write the Palace case, a separate `fem-em-palace` window runs it, and the `Study` reads S back. That S equals `OPS-61`'s recorded read on the same input: (E2) 0.000e+00 at degree 2 and (E1) 2.149e-14 at degree 1, both within `LEG_D0_REPRODUCTION_BAND`. It ties to our degree-2 record per class at (Tie2): 2.464e-10 / 9.668e-11 / 9.883e-11, each ≤ 2.6e-10. Fixture: ×1 exact strip at 64 MHz (Palace v0.18.1, orders 2 / 1; `20261018T140040Z_ARCH-1-U16b-W1-palace-d2.log`, `20261018T140512Z_ARCH-1-U16b-W2-palace-d1.log`, `20261018T172228Z_ARCH-1-U16c-W1.log`). Step 1c (§9 item 529, U62): `docs/guide/backends.md` and the U16 envelope entry, 0 dead links (`20261019T005254Z_ARCH-1-U62-Wdoc.log`). Step 1d (§9 item 590, U76): `femem run --backend palace` refuses every output the Palace route has no path for — `outputs.sar` / `b1_plus` / `probes`, `tune`, `study.sweep_hz`, `thermal` — by field with rc 2, before any mesh; `docs/guide/backends.md` states the carried / refused split from the code; U16's (E2) / (Tie2) / (c-order) / (E1) lines re-read byte-identical (`20261022T171719Z_ARCH-1-U76-RGU16.log`, `20261022T172109Z_ARCH-1-U76-Wbit.log`). Step 1e (§9 item 601, U79): `femem run --backend palace` refuses `geometry.as_hole` and a launch wider than one rank by field `--backend`, rc 2, before any mesh; a degree other than 1 or 2 is refused at parse (`study.degree`) (`20261023T003505Z_ARCH-1-U79-W0-real.log`, `20261023T003520Z_ARCH-1-U79-W0b-n2.log`). Not ✅: `ARCH-1`'s gate of record is the whole safety layer on two backends. *Was (embargo, lifted by the 2026-10-17 weekly):* no daily review could queue it; it waited on `OPS-61` step 1, `PORT-21` and one F-human physics step, so that the solver-agnostic interface is designed against two real backends. **Done-when:** digit-identical reproduction of every imported band and record on the DolfinX backend before any second backend is wired — a refactor that moves a record is a defect. *History: `docs/planning/chunks/ARCH-1.md`.* 

### §9 item 631 — 2026-10-25 (§10 U91, `ARCH-1` step 1f; chain step (solve); chain U, 9e)

*Written by the 2026-10-25 03:00 daily review from the 10-24 weekly's enumerated step U91 (§9 rule (4)). **One citation is corrected, not the target.** The weekly cites "`OPS-61` step 2c's recorded Cu read at degree 2 ((E2-Cu))". Step 2c is the **h = 0.0075** hole rung (226 904 tets, `OPS-61.md:1849`, `:2076`). The copper F-small **×1** hole read at order 2 — the fixture `ARCH-1` runs on — is **step 2b** (§9 item 241, `OPS-61.md:569`). "(E2-Cu)" is defined nowhere else, so it is defined below by analogy with U16's (E2). The 10-26 weekly confirms or corrects this reading.*

**Question.** Can `Study.solve(backend="palace")` run the **copper** F-small (Leontovich wall → Palace's `Boundaries.Conductivity`) and reproduce the `OPS-61` probe's own copper read?

**What must change in `src/` (a feature step; U16–U79 edited `backends/palace.py` the same way).** `solve_palace` refuses both inputs today: `as_hole=True` (`src/fem_em_solver/backends/palace.py:217–218`) and any conductor (`:219–220`). `write_case` writes volume materials only (`:175–183`), and `palace_config` writes only PEC on attribute 1 plus `LumpedPort` (`:99–111`).
- **Narrow acceptance:** `as_hole=True` **together with** `conductor=LeontovichWall(σ)` is accepted.
  - Every other combination keeps its refusal by field: `as_hole` without a conductor, a conductor on a solid geometry, and any non-Leontovich conductor.
  - The writer adds the hole-wall facet attribute and a `Boundaries.Conductivity` block with σ, ported from the probe: `scripts/probes/ops61_step1a_palace_ans4.py` (`export_hole_pec` `:1106`, `_hole_wall_attr` `:785`, the Conductivity block `:1533`). Cite those lines in the code.
  - PEC stays on attribute 1 only, as step 2b's config does (`OPS-61.md:1915`).
- **The YAML refusal stays:** `_palace_input_refusal` refuses `geometry.as_hole` (`study_yaml.py:1798–1799`). This step is Python only; the YAML half is the next U step, the weekly's.
- **Existing refusal tests stay green unedited.** If an existing test asserts exactly the newly accepted combination (`as_hole` + Leontovich) is refused, **stop and report**: that refusal is a gate, and changing it is the weekly's.

**Anchors.**
- **(E2-Cu)** The `Study`-written case on `ops61_ans6_x1_hole_strip.msh` with σ = 5.8e7 reproduces the probe's own Palace port S (`logs/ops61/postpro_h05x_cu_o2/port-S.csv`, on disk) at `LEG_D0_REPRODUCTION_BAND`, entry by entry. This is the same mesh, order and boundary: the analogue of U16's (E2), which read 0.000e+00 on the PEC route (`ARCH-1.md:391`).
- **(Tie-Cu)** The Study's class means against `ANS6_EXACT_STRIP_CLASS_MEANS_P2["64"]["cu"]` at the band step 2b's read asserted (T4) at. Import that band by name from the probe or its test, and say which. Step 2b read 1.612e-10 / 1.586e-10 / 2.033e-10 (`20261005T093623Z_OPS-61-item241-W6-read.log.gz:56`).
- **(RG-U16)** The PEC route's (E2) / (Tie2) re-read: (E2) 0.000e+00, (Tie2) 2.464e-10 / 9.668e-11 / 9.883e-11, 14 passed, 220 s (`ARCH-1.md:391`), unchanged.
- **Control (asserted):** the same case written with the conductor dropped (the PEC hole) misses the copper read, \|S(PEC) − S(Cu)\| > `NEGATIVE_CONTROL_FLOOR`. Backing: `WF-1`'s PEC-vs-copper miss of 6.422e-04 / 2.812e-04 / 4.863e-04 on the degree-1 classes (`20261020T202415Z_WF-1-s21-RGCu-copper.log:2716`). At degree 2 the miss is printed beside it.

**Windows (bare harness calls, foreground).**
1. W0: coding and refusal tests (every refusal field still refused, the new combination accepted at write time, no mesh).
2. Wwrite `ARCH-1-s1f-Wwrite`: the case writer at `-n 1` (it must run at `-n 1`, `palace.py:230–231`).
3. Wpalace `ARCH-1-s1f-Wpalace`: `docker compose -f docker/docker-compose.yml --profile palace run --rm -T fem-em-palace bash -lc 'cd <case_dir> && timeout -k 30 590 palace -np 8 case.json'` (`ARCH-1.md:159`'s shape). Step 2b's ×1 hole Cu run took 175.950 s Palace Total for four excitations (`OPS-61.md:1869`).
4. Wread `ARCH-1-s1f-Wread`: the Study reads S back; (E2-Cu), (Tie-Cu), the control.
5. RG-U16: the PEC regression, ≈ 220 s, heavy.

**Tier:** heavy (Wpalace ≈ 180–260 s, RG-U16 220 s). The solver container has no docker socket, so the Study writes the case, a separate window runs Palace, and the Study reads the result back (dashboard item 2). **Class:** chain step (solve). **Cost:** ≈ 6.0 measured / 60 predicted slot-minutes. **Park (rule (d)) if Wpalace has not started by minute 35.** **Executor:** implementer. **Depends on:** nothing. **Family:** `ARCH-1` step 1f, attempt 1 of 4.

**Traps.**
- **No Palace number set against an AED one** in any tracked file. Palace's own numbers go public at full detail, with the version (0.18.1), config and mesh cited.
- The mesh and configs under `logs/ops61/` are gitignored scratch: read them, never commit them.
- A Palace window is one-shot `run --rm`. The `fem-em-palace` service has its own profile and is not the solver service, so `pgrep` in the solver container does not see it.
- `-k 30` on every container `timeout`; never background a window.

**Scope boundary.** Python route only, copper ×1 hole only, degree 2. Not closed: the YAML route, other rungs, the loss partition (`OPS-61` step 2c excluded Palace's surface-loss postprocessing, `:2033`), field export (which H34 would need, returned to the weekly). `ARCH-1` stays 🟡.

**A negative result means:**
- (E2-Cu) red: the writer's case differs from the probe's. Diff the two JSON configs and the attribute census, report, park.
- (Tie-Cu) red with (E2-Cu) green: impossible by construction unless the record or band was misread. Stop.
- RG-U16 red: the edit moved the PEC route, a defect. Revert.

**Status it can move.** §10 U91 → DONE. `ARCH-1`'s state line gains "the copper hole on the Palace backend (step 1f)". The route's refusal list loses `as_hole` for the Leontovich combination only.

**Outcome of item 631 (2026-10-25, 07:30 CDT slot): DONE.** `backends/palace.py` only: `solve_palace` accepts `as_hole=True` together with `LeontovichWall(σ)` and nothing else (`as_hole` without a conductor or with a non-Leontovich one: `ValueError` "as_hole=True is not supported without a LeontovichWall conductor"; a conductor on a solid geometry: the old "a conductor is not supported"); `write_case(..., hole_wall_sigma=σ)` splits the exterior into box → attribute 1 and the wall (`BIRDCAGE_CONDUCTOR_SURFACE_TAG` = 401, read from `built["mesh_facet_tags"]`) → 401, drops the conductor volume material (phantom + air, the probe's order), and `palace_config(conductivity=…)` appends `Boundaries.Conductivity` `[{Attributes [401], Conductivity 5.8e7, Permeability 1.0}]`; PEC stays on attribute 1 only. Probe lines cited in the code. `_palace_input_refusal` untouched; no existing test edited (U16's `test_u16_backend_refusals` passes unedited: it refuses `as_hole` with no conductor). New `tests/package/test_arch1_s1f_palace_copper.py`. W0 complex `-n 1` 6 passed, 3 s (`20261025T124940Z_ARCH-1-s1f-W0.log`): the three other combinations refused, the new one passes every refusal and reaches `build` (stubbed); the written `Boundaries` / `Solver` equal `logs/ops61/h05x_cu_o2.json`'s. Wwrite `-n 1` 25 s (`20261025T124955Z_ARCH-1-s1f-Wwrite.log`): `case.json` `Domains` / `Boundaries` / `Solver` equal the probe's; 19 525 nodes both; triangle census on attributes 1 / 401 / 21–24 (3320 / 19928 / 40 × 4) and the tet census per cell tag equal `ops61_ans6_x1_hole_strip.msh`'s; the one difference is our writer's unreferenced interior remainder blocks 31–34 (48 / 48 / 50 / 50), as on the PEC route. Wpalace (Palace v0.18.1, `palace -np 8 case.json`, config `logs/arch1/s1f/cu_case/case.json`, mesh its `mesh.msh` ≡ the probe's mesh by census, order 2, 64 MHz, four excitations) exit 0, 157 s, Palace Total 155.195 s, peak 4.7 G (`20261025T125029Z_ARCH-1-s1f-Wpalace.log`). Wread 25 s (`20261025T125321Z_ARCH-1-s1f-Wread.log`): **(E2-Cu)** entry-wise max \|S − S_probe\|/\|S_probe\| 0.000e+00 against `postpro_h05x_cu_o2/port-S.csv` (Frobenius 0.000e+00; band `LEG_D0_REPRODUCTION_BAND` 1e-9), cells 80 798; **(Tie-Cu)** 1.612e-10 / 1.586e-10 / 2.033e-10, step 2b's (T4) digits exactly. **Band disclosure:** step 2b's (T4) was *printed* (predicted ≤ 1e-6, `ops61_step1a_palace_ans4.py:1679–1690`), never asserted, so no (T4) band exists to import; (Tie-Cu) is asserted at `TIE2_BAND` = 2.6e-10 imported from `tests/package/test_arch1_palace_backend.py` (U16's ratified PEC-tie ceiling, tighter than the 1e-6 prediction). **Control (asserted):** the probe's PEC-hole run on the same mesh and order (`postpro_h05x_pec_o2`, the wall PEC instead of Conductivity — the probe's own copy of "the same case with the conductor dropped"; no second Palace run, because the dropped-conductor case is refused by design) misses the Study's copper class means by 6.332e-04 / 2.867e-04 / 4.949e-04 (> `NEGATIVE_CONTROL_FLOOR` 1e-5), beside `WF-1`'s degree-1 6.422e-04 / 2.812e-04 / 4.863e-04. `pgrep -c python3` 0, then RG-U16 complex `-n 1` 14 passed, 217 s (`20261025T125409Z_ARCH-1-s1f-RGU16.log`): (E2) 0.000e+00, (Tie2) 2.464e-10 / 9.668e-11 / 9.883e-11, (c-order) 6.6674 / 2.5442 / 4.5224 %, (E1) 2.149e-14, unchanged. No registry line added (no new refusal type), so no page regeneration. Not claimed: the YAML route (still refuses `geometry.as_hole`), other rungs, Palace's surface-loss postprocessing, field export.
