# MAT-19 — The virtual temperature probe

**The virtual temperature probe: a point time series with a stated probe volume — Phase 8, gel** (opened 2026-10-03 by the weekly; §10 chain P7, `PROJECT_PLAN.md:4587–4592`; lane 3). §7 row `PROJECT_PLAN.md:2512`. **Done-when (§4, the row):** an additive probe post-processor in `thermal/` (ball-averaged temperature against time, owner rule as `POST-8`); on `MAT-12`'s SAR-heated sphere the series at three radii is within 1 % of the closed form at every recorded time, and the series is width-invariant to 1e-9 at `-n 1 / 2 / 4`; a displaced-probe control misses; elapsed recorded ⇒ ✅. Not claimed: a real probe's thermal mass or lead artefact.

### §9 item 243 — 2026-10-05 (chain P7)

**Outcome (2026-10-05 09:30 UTC slot): ✅ GREEN, `MAT-19` ✅ (the row was amended under rule (j) by the 03:00 review).** W0 (`20261005T095530Z_MAT-19-step1-W0.log`, 87 s, 14 passed): (U1) 0.0, (U3) 2.0e-14 / volume 0.0, (U3c) ghost miss 8.3e-2 volume, (U4) bitwise; `MAT-12` re-gated at its record ((P_t) 1.6845e-03). W1 `-n 2` (`…095808Z_…W1.log`, 104 s): (N1) 1.168e-14·θ_ss(0), a → 0 1.1e-16; (A0) 19 496 cells, V_h/V − 1 ≤ 6.3e-4; **(T) max 1.7315e-03** (n = 22, c₀); (G) ≤ 1.12e-5, so (T) is FEM error; q_h series +3.47 … +3.73 % at t\* (printed); imag ratio 0; **(C1) miss max 5.4207e-02** (predicted 5.466e-2). W2 `-n 4` (`…100001Z_…W2.log`, 99 s) and W3 `-n 1` (`…100145Z_…W3.log`, 117 s): **(W) `series_cf` 6.661e-16**, cells equal; `series_h` 3.1e-11 / 7.9e-11 (printed). (C2) as printed is FEM-ball / point − 1 (−0.629 %, −0.815 %, −1.353 % at t\*), the sign-inverted form of the brief's point / exact-ball offsets (+0.732 / +0.872 / +1.275 %), shifted by the FEM error; it is a record only.

`MAT-19` step 1: the ball-averaged probe on `MAT-12`'s sphere (implementer; standard by prediction, **heavy by ceiling** (`timeout -k 30 590`), `-n 2` / `-n 4` / `-n 1`, complex; ≈ 16 measured / ≈ 40 predicted min; depends on nothing). This is the family's first step: attempt 1 of 4.

**The owner-rule trap, and how this item routes around it (read first).**
- `MAT-12`'s thermal solve is **on a submesh**: `_sphere_submesh` (`test_mat10_sar_to_heat_sphere.py:74–89`, `dmesh.create_submesh`), called at `test_mat12_transient_sar_heating_sphere.py:159`.
- `evaluate_vector_field_parallel(owner="lowest_original_cell")` raises `RuntimeError: Mixed topology unsupported.` on any `create_submesh` mesh (`docs/testing/known-issues.md:1372–1379`, OPEN). That entry's Watch clause says `MAT-12`'s sphere is "full-mesh" and unaffected. **That clause is wrong**: the thermal field lives on the tag-1 submesh.
- **Route.** A ball average is an integral, not a point sample, so it needs no point owner. Each rank sums its **owned cells only** (`index_map(3).size_local`, never ghosts). It then `comm.allreduce(SUM)`s the pair (∫θ, ∫1) over the probe ball and divides. Every cell is counted exactly once at every width, by construction of the cell index map. That is the ownership discipline the `POST-8` rule exists to enforce, applied to cells instead of points.
- **Never call** `evaluate_vector_field_parallel(..., owner=...)` or `f.eval` on the submesh, and never read `original_cell_index` there.
- **Wording flag for the review (rule (j)).** The row's parenthetical reads "owner rule as `POST-8`". This item does not call `POST-8`'s point-owner rule. It meets the rule's intent, width invariance to 1e-9 at `-n 1 / 2 / 4`, with the row's own gate (W). Before writing "→ ✅", the review should amend the parenthetical in its own commit to: *"(ball-averaged temperature against time: an integral over owned cells, no point owner — width-invariant as `POST-8`)"*. §10 P7's "through `evaluate_vector_field_parallel`'s owner rule" is the weekly's text and the weekly re-reads it. If the row is not amended, the executor closes on the row as written: it holds 🟡 and flags the clause.

**The fixture (imported, never restated).**
- **The EM side** is `MAT-10`'s and `MAT-12`'s: `_solve_lossy_sphere(SIGMA_HIGH, R / 10.0, R / 5.0, return_fields=True)` (`test_mat12_transient_sar_heating_sphere.py:141`). Here R = 0.01 m, h_in = R/10 = 1 mm, ρ 1000, σ 0.57 S/m, 64 MHz.
- **The thermal side** is the tag-1 submesh with every exterior facet Robin (`_robin_tags`, `test_mat10_sar_to_heat_sphere.py:92–98`). Its parameters are `K` 0.6, `H` 100 (Bi 5/3), `RHO_C` 4.18e6, T∞ = T(0) = 0 and `G` (`:52–55`; MAT-12 `:67–72`).
- **q_cf** = `G·RHO_KG_M3·SAR_cf` = 2002.088514 W/m³, as at MAT-12 `:146–150` (log `20260930T110230Z_MAT-12-W1.log:157`).
- **The submesh** is **19 496 cells** with volume 4.173655e-06 m³ (`…MAT-12-W1.log:159–160`; reproduced in `20261003T171523Z_MAT-17-step2-W3-regate.log:159`). Mean cell volume is 2.14e-10 m³, an ≈ 1.2 mm equilateral edge, consistent with h_in = 1 mm.
- **Recorded times.** The run is `MAT-12`'s own `_cn(sub, ft, q, snapshot_steps=range(N_MAIN + 1))` (`test_mat12_transient_sar_heating_sphere.py:115–125`; the `snapshot_steps` pass-through is additive, `EX-67`). It is CN with a Rannacher start, N = 40 to t\* = 0.2 t_R = 139.3333 s, so Δt = 3.4833 s.
  - The recorded times are **t_n = n·Δt, n = 1 … 40** (snapshot n is full step n; the two BE half-steps produce snapshot 1, `thermal/transient.py:323–330`).
  - n = 0 is θ ≡ 0. It is asserted to read exactly 0 and is excluded from relative checks.
- **The run** reproduces `MAT-12`'s (P_t) interior relative L2 1.6845e-03 at t\* (`…MAT-12-W1.log:161`; `…MAT-17-step2-W3-regate.log:161`). Print it beside.

**The probes (pinned).**
- Probe radius **a = 2 mm = 2 h_in = 0.2 R**, ≈ 156 mean cell volumes per ball. The volume is a fixture parameter, not a real probe size.
- The three centres:
  - c₀ = (0, 0, 0);
  - c₁ = 3.5 mm · (1, 0, 0);
  - c₂ = 6.5 mm · (1, 1, 1)/√3 (off-axis, generic).
- **(M), asserted exactly:** \|cⱼ\| + a ≤ 0.9 R for every probe. The ball must lie inside `_interior_indicator`'s r < 0.9 R (`test_mat10_sar_to_heat_sphere.py:140`), which is where `MAT-12`'s (P_t) 1.68e-3 was measured, and away from the faceted surface (mesh volume −0.36 % from 4πR³/3).

**Code.**
- **(1) New `src/fem_em_solver/thermal/probe.py`** (additive), plus exports in `thermal/__init__.py`.
  - `ProbeBall(centre: tuple[float, float, float], radius: float)`. The centre is plain floats.
  - `ball_average(temperature, probe, *, levels=3) -> (value, volume)`.
  - `probe_series(snapshots: Mapping[int, fem.Function], dt: float, probes, *, levels=3) -> ProbeSeries`, with fields `times`, `values` (n_t × n_probes), `volumes`, `nominal_volumes` = 4πa³/3 and `max_imag_ratio`.
  - **Restrictions, which raise:** scalar Lagrange degree 1 on an affine tetrahedral mesh (`geometry.cmap.degree == 1`). The ball must be given in mesh coordinates.
  - **Algorithm (pure numpy, no UFL form, no JIT, no comparison on a UFL operand).**
    - Per owned cell K (`0 … size_local − 1`): the dof coordinates come from `V.tabulate_dof_coordinates()[V.dofmap.list[K]]`. For P1 these are the vertices, ordered like the dofs by construction, so never mix them with `geometry.dofmap`. Take the centroid and the bounding radius ρ_K = max vertex distance. Skip the cell if \|centroid − c\| − ρ_K > a.
    - Each remaining cell is sampled on a fixed barycentric lattice: the centroids of the 8^L equal-volume sub-tetrahedra of L-fold red (Bey) refinement of the reference tetrahedron. L = 3 gives 512 points of weight \|K\|/512, precomputed once.
    - At each lattice point, θ = λ·θ_K (exact for P1). The indicator is the squared distance compared with a² in numpy float64.
    - **Geometry is precomputed once per probe.** Each snapshot is then one weighted sum. Return the real part and record max \|Im θ\|/max \|Re θ\| over the dofs used (the complex build).
  - **Reduction:** one `comm.allreduce(np.array([num, vol]), op=MPI.SUM)` **on every rank**, including ranks with zero candidate cells. A collective entered by a subset of ranks deadlocks (item 227, `eda85f26`).
  - An internal `_ball_average_on_cells(..., cells)` takes the cell list. The public path passes owned cells only. The test uses the internal path for the planted ghost control.
  - The docstring cites this item and states the restriction and the known-issues route.
- **(2) New `tests/solver/test_thermal_probe.py`**, smoke tier: `dolfinx.mesh.create_box` on [−1, 1]³, 8³ × 6 tets, P1 real-valued fields. In the complex build they are complex-typed with Im 0.
  - **(U1)** θ ≡ 3.7 inside a ball of radius 0.5 at (0.1, 0.2, 0.05) reads 3.7 to ≤ `B_BAND` 1e-10 relative (`test_mat10_sar_to_heat_sphere.py:61`, imported). V_h/(4πa³/3) − 1 is printed, predicted \|·\| ≤ 5e-3.
  - **(U2)** The affine θ = 2 + x − 0.5y + 0.25z reads θ(c). This is **printed, not asserted**: the predicted relative deviation is ≤ 1e-3, from the lattice centroid shift ≈ 3h_s²/(2√π a) with h_s = h/8.
  - **(U3)** A ball of radius 10 ⊇ the box, on θ = interpolant of 2 + x² + y. The probe equals `_allsum(θ dx)/_allsum(1 dx)` (MAT-10's `_allsum`, `:70–71`) to ≤ `B_BAND`, and the volume equals 8 to ≤ `B_BAND`. This is the rank-safety identity. It is visible only at `-n 2`, the width where a missing reduction or a double-counted ghost shows.
  - **(U3c) Planted ghost control, at `comm.size > 1` only:** `_ball_average_on_cells` over owned + ghost cells must miss (U3)'s identity by more than `B_BAND`. That miss is asserted. The predicted separation is O(1e-2), printed. Rule (e): no prior measurement, so the factor is predicted.
  - **(U4)** At `-n 2`, the c₀ reading is bitwise identical on every rank (allgather, `==`).
- **(3) New `tests/validation/test_mat19_virtual_temperature_probe.py`** (`complex_only`), with one module-scoped fixture. It imports and never restates:
  - `E0`, `EPSILON_R_SPHERE`, `OMEGA`, `RHO_KG_M3`, `SIGMA_HIGH`, `_solve_lossy_sphere` (`test_lossy_sphere_sar.py`);
  - `R`, `K`, `H`, `G`, `ROBIN_TAG`, `P_BAND`, `A_BAND`, `B_BAND`, `_allsum`, `_say`, `_sphere_submesh`, `_robin_tags`, `_cell_average_dg0`, `_to_sub` (`test_mat10_sar_to_heat_sphere.py`);
  - `RHO_C`, `T_R`, `T_STAR`, `N_MAIN`, `_cn` (`test_mat12_transient_sar_heating_sphere.py:67–72,115`);
  - `robin_sphere_eigenvalues`, `robin_sphere_transient_coefficients`, `robin_sphere_transient_temperature` (`utils/analytical.py:1192,1220,1232`).

  Two runs on the same submesh, **q_cf** (uniform) and **q_h** (`MAT-12`'s DG0 FEM SAR, `:164–166`), each with `snapshot_steps=range(41)`. The probe series is taken on both.
  - **Width artefact.** At `-n 2`, rank 0 writes `/workspace/logs/mat19/step1-n2.npz` (`/logs/` is gitignored, `.gitignore:82`; `os.makedirs(..., exist_ok=True)`) with `series_cf`, `series_h`, `volumes` and `cells`, then `Barrier`. At `-n 4` and `-n 1` the module loads it. If it is missing, `pytest.fail`; never recompute (`POST-8` step 3a's guard, `docs/planning/chunks/POST-8.md:171–173`).

**The comparand: the closed form averaged over the same ball (pinned).**
- θ(r, t) is `robin_sphere_transient_temperature` with θᵢ = 0, T∞ = 0: θ = θ_ss(r) − Σₙ Bₙ e^{−λₙ² αt/R²} Xₙ(r), with Xₙ = sin(κₙr)/(κₙr), κₙ = λₙ/R, and Bₙ from `robin_sphere_transient_coefficients` (`analytical.py:1220–1229`, 400 terms as the default).
- Its exact mean over the ball B(c, a), d = \|c\|, is:

  θ̄(d, a, t) = qR/(3h) + q(R² − d² − 3a²/5)/(6k) − Σₙ Bₙ e^{−λₙ²αt/R²} Xₙ(d) · 3[sin(κₙa) − κₙa cos(κₙa)]/(κₙa)³.

  - The steady part: θ_ss is quadratic, and the ball mean of \|x\|² is d² + 3a²/5.
  - The modes: Xₙ solves ∇²u + κₙ²u = 0 in all of ℝ³, so the Helmholtz mean-value theorem holds. Its spherical mean at radius s is u(c)·sin(κs)/(κs). Integrating with weight 3s²/a³ gives the factor 3j₁(κa)/(κa).
  - A small-argument guard: use 1 − x²/10 for x < 1e-4.
- **(N1), asserted numpy anchor.** θ̄ above against a direct radial Gauss–Legendre quadrature of `robin_sphere_transient_temperature`. The quadrature uses the shell-cap weight: the area of the sphere of radius r inside B(c, a) is πr(a² − (r − d)²)/d on \|d − a\| ≤ r ≤ d + a, plus 4πr² on r < a − d, normalised by 4πa³/3. Use ≥ 200 nodes per sub-interval. The check runs at all 3 probes × 40 times, together with the a → 0 limit against the point value. The band is ≤ 1e-10·θ_ss(0), the series-route agreement band of `MAT-17` step 1 (P3), imported from `test_mat17_duty_cycle_closed_form.py:206–207`, where θ_ss(0) = 0.1223499 K. Predicted ≤ 1e-13.
- **Predicted comparands** (pure-Python evaluation of the formula above while drafting; it reproduces the logged θ(0, t\*) 0.059564 and θ(R, t\*) 0.035599 of `…MAT-12-W1.log:158`):

  | probe | θ̄(t₁) [K] | θ̄(t₁₀) [K] | θ̄(t\*) [K] | point θ(cⱼ, t\*) [K] |
  | --- | --- | --- | --- | --- |
  | c₀, d = 0 | 1.668407e-03 | 1.666485e-02 | 5.913173e-02 | 5.956431e-02 |
  | c₁, d = 3.5 mm | 1.668407e-03 | 1.658213e-02 | 5.678637e-02 | 5.728180e-02 |
  | c₂, d = 6.5 mm | 1.668186e-03 | 1.595905e-02 | 5.025018e-02 | 5.089074e-02 |

  The probe reads these to the printed digits (n = 1 is q_cf t₁/ρc = 1.6684e-3 at the interior probes).

**Anchors (asserted).**
- **(A0)** The submesh has 19 496 cells at every width, compared exactly (`…MAT-12-W1.log:159`). Probe volume ratios V_h/(4πa³/3) − 1 are printed, predicted \|·\| ≤ 2e-3. Snapshot 0 reads exactly 0.0 at every probe.
- **(T), the row's gate.** For every probe j and every n = 1 … 40, on the **q_cf** series: \|P_h,j(t_n) − θ̄_cf(dⱼ, a, t_n)\|/θ̄_cf(dⱼ, a, t_n) ≤ **1 %**. This is the row's band, carried by the imported constant `P_BAND` = 1.0e-2 (`test_mat10_sar_to_heat_sphere.py:59`, "`MAT-9` (i)").
  - **Predicted:** max over (j, n) ≤ 3e-3, ≥ 3× headroom. `MAT-12`'s (P_t) interior L2 on the same run is 1.6845e-3 at t\* (`:161`), and for n ≤ 3 the interior probes sit in the uniform-growth core, which P1 + BE/CN reproduce to round-off.
  - **Why q_cf and not q_h.** On `MAT-12`'s sphere only the uniform closed-form-SAR leg, (P_t), carries a 1 % band. The FEM-SAR leg (A_t) carries `A_BAND` 11 % and measured 3.4708e-2, almost all of it `MAT-4`'s +3.54 % SAR offset (`:162–163`). So the q_h series against the closed form is **printed** beside `A_BAND` (predicted ≈ +3.4–3.6 % at t\*) and not asserted.
- **(G), printed.** The probe's geometric error alone: the same lattice, with θ_cf evaluated in numpy at the lattice points (no FEM), against θ̄ at t\*, per probe. The predicted \|·\| is ≤ 1e-3; the dominant term is centroid shift ≈ 0.01 mm × the outer probe's sensitivity of 5.5 %/mm. The split of (T) into geometry (G) and FEM (T − G) is printed per probe at t\*.
- **(W), the row's gate, at `-n 4` and `-n 1` against the `-n 2` npz.**
  - `cells` must be equal.
  - Every entry of `series_cf`, n = 1 … 40 × 3 probes, must be within **1e-9 relative** of its `-n 2` value. This is the row's band, as `WIDTH_BAND = 1.0e-9  # MAT-19 row, PROJECT_PLAN.md:2512`, a band this item does not set.
  - Predicted ≤ 1e-12: the thermal solve is real-valued, and the probe is a sum whose order alone moves with the width.
  - **Printed:** the same comparison on `series_h`. It carries the complex EM solve's width round-off. Predicted ≤ 3e-10: the worst cross-width figure this project has measured is `POST-8`'s 2.98e-10 (`docs/planning/chunks/POST-8.md:13`), and `MAT-15`'s P₁ spread across six runs is ≤ 5.3e-12.
  - **Printed:** `max_imag_ratio`, predicted 0 for q_cf.

**Negative controls (rule (e)).**
- **(C1) The displaced probe: the row's control, its miss asserted, its factor predicted.** The FEM probe is moved 1 mm radially inward, c₂′ = (5.5/6.5)·c₂ with the same a, and its q_cf series is compared with the **nominal** θ̄_cf(6.5 mm, a, t_n). **Asserted: max over n of the relative miss > `P_BAND`** (it misses (T)).
  - **Predicted, printed beside** (closed form only; no prior measurement on this fixture):

    | n | 1 | 5 | 6 | 10 | 20 | 40 |
    | --- | --- | --- | --- | --- | --- | --- |
    | miss | 1.28e-4 | 8.15e-3 | 1.089e-2 | 2.116e-2 | 3.887e-2 | **5.466e-2** |

    The miss grows monotonically with n. It first exceeds 1 % at n = 6, so the assertion is on the maximum over n, which is n = 40, never on "every n". An early-time series is uniform in r, so no displacement can show there.
  - **Ceiling:** the steady limit, t → ∞, is 6.843e-2. That is the most a 1 mm radial displacement at 6.5 mm can show on this fixture, 6.8× the band. At t\* the predicted miss is 5.5× the band, less the FEM error of ≤ 3e-3.
  - **Printed:** the displaced probe against its **own** closed form θ̄(5.5 mm), predicted within (T)'s 1 %. It is a fourth radius for free and is not asserted.
- **(C2) The volume-omitted comparand: predicted, printed only.** The q_cf ball series is compared with the **point** closed form θ(dⱼ, t_n). At t\* the predicted offsets are +0.732 % (c₀), +0.872 % (c₁) and +1.275 % (c₂). At c₂ the point comparand exceeds 1 % from n ≈ 20 (1.14 %). This is what stating the probe volume buys at this a/R. It is printed as a record because the separation from the band (1.27×) is too thin to assert.
- **(U3c)** is in W0 (see above): the ghost-included sum misses `B_BAND`. The miss is asserted; the factor is predicted.

**Windows** (complex; foreground; Bash `timeout: 660000`; each command bare, whole and repo-relative; `-s`; no `-k` inside the quoted command; no pipes or redirects).
- **First, two checks.** `pgrep -c python3` in the container reads 0. Then run the rubric's orphan sweep as one bare call:

  `scripts/testing/run_and_log.sh MAT-19-step1-cache-sweep "docker compose exec -T fem-em-solver bash -lc 'cd /root/.cache/fenics && for c in libffcx_*.c; do s=\${c%.c}; ls \$s.*.so >/dev/null 2>&1 || { echo orphan \$s; rm -f \$s.c \$s.c.cached \$s.o; }; done; echo sweep done'"`
- **Pricing basis.** `MAT-12`'s own window took 77 s at `-n 2`, test body 30.5 s, EM 15.4 s (`20260930T110230Z_MAT-12-W1.log:171`, footer). Its 2026-10-03 re-gate read EM 14.4 s / body 21.1 s inside a 105 s window together with `MAT-14` (`20261003T171523Z_MAT-17-step2-W3-regate.log:171,396`).
  - Today's box is slow: `MAT-5` step 3 measured its EM at ≈ 2× its own basis (`20261005T050759Z_MAT-5-step3-W1-alpha0.log`; journal in `docs/planning/chunks/MAT-5.md:96`). Hence the 590 ceiling on a standard prediction. **A window killed at 590 is not a result:** stop BLOCKED as a slow box, never lengthen.
- **W0, the probe and the rule (c) re-gate** (`-n 2`; predicted ≈ 100 s, up to 200 s on a slow box).
  - The `thermal/__init__.py` edit is imported by `MAT-12`'s module, so that module re-runs green in the same window. Its (P_t) / (A_t) / (S1) / (S2) / (E) / (∞) / controls must reproduce `…MAT-12-W1.log:161–170` at their bands.
  - A red W0 stops the item before W1: revert the `src/` edit and park.

  `scripts/testing/run_and_log.sh MAT-19-step1-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 590 mpiexec -n 2 python3 -m pytest tests/environment tests/solver/test_thermal_probe.py tests/validation/test_mat12_transient_sar_heating_sphere.py -s -v --tb=short'"`
- **W1, `-n 2`, which writes the npz** (predicted ≈ 85 s: environment ≈ 45 s, EM ≈ 15 s, two CN runs ≈ 1–4 s each, numpy comparands and (G) ≈ 5 s):

  `scripts/testing/run_and_log.sh MAT-19-step1-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 590 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_mat19_virtual_temperature_probe.py -s -v --tb=short'"`
- **W2, `-n 4`** (predicted ≈ 80 s): the same command with `-n 4`, label `MAT-19-step1-W2`.
- **W3, `-n 1`** (unmeasured; predicted ≈ 110 s, with the EM ≤ 2× its `-n 2` time): the same command with `-n 1`, label `MAT-19-step1-W3`.
- **Order:** W0 → W1 → W2 → W3. W2 and W3 need W1's npz from the **same** module bytes.

**Traps.**
- **The owner rule on a submesh raises** (known-issues `:1372`). Use no point evaluation at all, with or without an owner. The probe is an owned-cell integral.
- **`cell_tags.values` and `assemble_scalar` are rank-local.** The probe's partial sums are rank-local until the one allreduce. Every rank enters that allreduce (item 227 deadlocked on a rank-0-only collective). (U3) / (U3c) at `-n 2` are the check.
- **Ghosts.** `V.dofmap.list` has rows for owned + ghost cells. Iterate `range(size_local)` only. The ghost **dof values**, however, are needed, and the snapshots copy `T_new.x.array` after `scatter_forward` (`transient.py:265–269,302`).
- **Complex build.** θ arrays are complex-typed. Take `np.real` after the weighted sum and print the imaginary ratio. Put no ordering comparison on any UFL operand: the indicator is numpy-side on float64 squared distances, never `ufl.conditional` / `lt` / `sqrt` (`OPS-22`; `WF-6` step 4c). `_allsum` returns complex: take `np.real`.
- **Dof coordinates, not the geometry dofmap,** for the P1 vertices of a cell. The ordering matches by construction only that way.
- **The ball centre is plain floats**, not a numpy scalar of a complex dtype (`POST-8` 3a's trap, `sar.py:257`).
- **Snapshot indices are full steps.** With `startup_half_steps=2`, snapshot 1 is after both BE half-steps, at t = Δt. Times come from n·Δt, never from counting balance entries (41 balance residuals for 40 steps, `…MAT-12-W1.log:161`).
- **No cross-run comparison is bitwise** (`OPS-43`). (W) is relative 1e-9; only (U4) is bitwise, and it is in-run.
- **Orphan `.c`** sweep before W0. After any kill, run `pgrep -c python3` first.
- **Rule (i):** no edit to `probe.py` or either test module after its last window. An edit after W1 re-runs W1, W2 and W3, because the npz must come from the committed bytes. Check file times against the log starts before committing.
- **The durable-capture idiom is not needed:** every window is ≪ 10 min. Never pipe pytest; the footer must be pytest's status.

**Scope.** This step covers one fixture (`MAT-12`'s 1 cm sphere, one mesh at R/10, P1, CN N = 40), three probes of one stated volume (a = 2 mm), the q_cf series gated and the q_h series printed, and recorded times only (no interpolation between steps). **Not covered:**
- a real probe's thermal mass, self-heating or lead artefact (the row);
- a fibre-optic probe's actual sensing volume;
- an h-convergence statement on the probe;
- the phantom or coil fixtures (`MAT-15`, `ANS-9`);
- tissue or perfusion;
- any compliance or ASTM F2182 reading;
- the submesh owner-rule fix itself. The known-issues entry stays OPEN; it resolves only with a `post/evaluation.py` fix.

**Status it can move.**
- **`MAT-19` ⬜ → ✅** per the row: (T) green at three radii and 40 times, (W) green at `-n 1 / 2 / 4`, (C1) misses, elapsed recorded in the row. This is conditional on the review's rule (j) amendment of the "owner rule as `POST-8`" parenthetical (see the flag above). Without it the executor holds **🟡** and flags the clause.
- §10 **P7** can then be read closed by the weekly.
- In the closing commit, the executor corrects the known-issues Watch clause at `:1378`: `MAT-12`'s thermal is on a submesh too, and `MAT-19` routes around the rule by integration. The entry stays OPEN.

**Negative result.** In every case: report the measurement, park on `attempt/MAT-19-step1-*`, mark the item BLOCKED in the same commit (rule (d)) and journal it in attempts.md. Never widen 1 % or 1e-9, and never change a, L or the radii after reading.
- **(N1) red:** the mean-value comparand or the shell-cap quadrature is wrong. This is a numpy bug: stop before any FEM reading is interpreted.
- **(U1) / (U3) / (U4) red:** a probe or reduction bug. Stop at W0.
- **(A0) red:** the mesh has drifted from `MAT-12`'s record. Stop; a `record-reconciler` ruling is needed, never a new count.
- **(T) red with (G) ≤ 3e-3:** the FEM probe misses the closed form at a radius or time where `MAT-12`'s interior L2 says it should not. Print the full 3 × 40 table and the (G) / (T − G) split. The review rules.
- **(T) red with (G) > 3e-3:** the lattice does not resolve the ball. Report it. The review re-poses L or a; the executor does not raise L in-slot.
- **(W) red on `series_cf`:** the probe or the thermal solve depends on the width. Print the entries and the cell counts, and stop.
- **(C1) does not miss:** the FEM series is insensitive to position, which contradicts (T). Stop.
- **A slow-box kill at 590:** BLOCKED. W0 and W1 records stand, and the unblock is a re-run on a box whose `MAT-12` EM reads ≤ 20 s.

No `PREDICTED RED:` line. A red anywhere stops the slot.

### §9 item 253 — 2026-10-05 (chain P8, daily provisional)

*Drafted by a general-purpose agent for the 2026-10-05 10:30 daily review; the review reviews it.* Licence: §9 rule (6), a provisional next step on chain P after P7 (`MAT-19` ✅, item 243). The weekly ratifies or reverts it. No band, tolerance, Done-when or validation target is invented: every asserted anchor is an identity at an imported band.

`MAT-19` step 2: step 1's probe, unchanged, read on `MAT-15` P1's coil-driven gel phantom at `ANS-9`'s SPEC probe centres (implementer; heavy by ceiling (`timeout -k 30 590`), `-n 2` / `-n 4` / `-n 1`, complex; ≈ 12 measured / ≈ 40 predicted min; depends on nothing). Family `MAT-19` step 2: attempt 1 of 4.

**Why this step (read first).**
- `ANS-9`'s SPEC (`examples/ansys_benchmarks/ans9_coil_driven_gel_heating_10MHz/SPEC.md`, row 4) asks for the **point** rise at t = 1254 s at five probes, metres: A (0, 0, 0), B (0.020, 0, 0), C (0, 0.020, 0), D (−0.020, 0, 0), E (0, 0, 0.030). Its readout puts A–D at 10 % and E as a record. `ans:9` produces row 4 by point evaluation through the legacy `owner=None` path. θ is continuous P1, so that is sound (`09_coil_driven_gel_heating_10MHz.py:286–291`). Our side reads, at `20261004T142653Z_ANS-9-step1a-W1.log:1943`: A 0.905798, B 3.564854, C 0.685305, D 0.697444, E 0.958639 K.
- A real gel-phantom measurement is a fibre-optic probe with a sensing volume, placed to a tolerance. The 09-30 weekly's sampling ruling says "public comparison rows quote cell or mass averages" (`POST-8.md:24`). Row 4 is the one SPEC row that quotes neither.
- This step therefore measures two things on our own side, before the operator's AED run: **the probe-volume sensitivity of row 4** (a stated ball against the point) and **its placement sensitivity** (a 0.1 mm displacement at B). It also asserts the identities that make a ball-averaged row 4 trustworthy on this fixture: whole-phantom = volume mean, width invariance with the EM in the loop, and the σ_y mirror image of the one-port drive.
- Which probes: **the SPEC's centres, as named**, plus C′ = (0, −0.020, 0), the σ_y image of C. MAT-19 step 1's three-radii analogue is not used, because the SPEC names its points.
- **What is new against step 1.** Step 1's (W) was already on a `create_submesh` mesh (`MAT-12`'s tag-1 sphere; known-issues `:1378` as corrected). It asserted (W) on the closed-form-source series and only **printed** the EM-sourced `series_h` (3.1e-11 / 7.9e-11). Here the asserted series carries a complex MUMPS solve on 199 920 cells at every width. (A) is step 1's (U3) identity, until now asserted only on a `create_box` mesh, on a coil-driven submesh partition.

**The fixture (imported, never restated).**
- The fixture is `MAT-15` step 1's: `_run_heated(return_fields=True)` (`tests/validation/test_mat15_coil_driven_heating_phantom.py:180–351`). That is `ANS-2`'s 10 MHz loaded F-small birdcage at `phantom_resolution` 0.0025, 199 920 / 58 866 cells, N1curl degree 1, ragged f = 0.5 sheets with κ none, the four single-leg drives, G = 1 W / P₁, P1 thermal on the tag-3 submesh with Robin on every face, and gel k 0.6, ρc 4.18e6, h 10, T∞ = T(0) = 20.
  - This is an **existing** fixture. The 2026-10-02 exact-rectangle ruling binds new fixtures, and ragged-sheet records stay valid as records. This step builds nothing new and its readings are records on the fixture `ANS-9`'s SPEC names.
- **Drive k = 0 only** (`out["by_k"][0] == "P1"`, asserted). Take `sub`, `ft`, `q_h = out["q_h"][0]`, `v_sub` and `T_star[0]`.
- **The series.** Re-run CN with every full step kept: `_transient(sub, ft, q_h, dt=T_STAR / N_MAIN, n_steps=N_MAIN, theta=0.5, startup_half_steps=2, snapshot_steps=range(N_MAIN + 1))`. This is `ans:9`'s `snap_run` (`09_…py:222–224`), with `_transient` imported from the `MAT-15` module. Δt = 31.35 s, and t_n = n·Δt for n = 1 … 40. The SPEC times 313.5 / 627 / 940.5 / 1254 s are n = 10 / 20 / 30 / 40.
  - Also run one `_be_one(sub, ft, q_h)` (Δt₁ = 2.5e-6 t_R) for (S).
- **Imports.**
  - From the `MAT-15` module: `_run_heated`, `_transient`, `_be_one`, `A`, `N_MAIN`, `T_STAR`, `T_INF`, `DT_1`, `RHO_C`, `E_BAND`, `B_BAND`, `V_BAND`, `V_SUB_RECORD_M3`, `F_CONTROL_FLOOR`, `C4_COVARIANCE_BAND`.
  - `PHANTOM_HEIGHT` from `tests.mesh.test_birdcage_port_tags`.
  - `WIDTH_BAND` from `test_mat19_virtual_temperature_probe.py:85`.
  - `ProbeBall`, `probe_series` and `ball_average` from `fem_em_solver.thermal`, plus the internal `_ball_average_on_cells`, `_lattice` and `_owned_cells` from `thermal.probe`.
  - The SPEC centres come from `ans:9`'s `SPEC_THERMAL["probes"]` (`09_…py:91–97`), loaded by `importlib` as `ans:9` loads `ans:2` (`_load_example`, `:102–107`; the stem starts with a digit). Never restate them.
- **No edit** to `src/`, to the `MAT-15` module, to the `ans:9` script, or to its `metrics.json` / `COMPARISON.md`. Each owes a gate's window under rule (a). Do not add a drive selector to `_run_heated`: the three unused drives cost ≈ 60 s per window, which is cheaper than two re-gates.

**The probes (pinned).**
- Radius **a = 5 mm = 2 h_p**, as step 1's 2 h_in. The mean cell volume is V_sub/58 866 = 3.84e-9 m³, a ≈ 3.2 mm regular-tet edge, so the ball holds ≈ 136 mean cell volumes against step 1's 156. a is a fixture parameter, not a fibre-optic sensing volume.
- **Centres:** A, B, C, D and E from the SPEC, plus C′ = (0, −0.020, 0).
- **(M), asserted exactly:** √(x² + y²) + a ≤ `A` and \|z\| + a ≤ `PHANTOM_HEIGHT`/2 for every centre. The closest are B / C / C′ / D, at 25 mm against a 30 mm radius, and E, at 35 mm against a 40 mm half-height. Lateral faceting sags inward by ≤ h²/(8a) ≈ 2.6e-5 m. V_h/(4πa³/3) − 1 is printed, predicted \|·\| ≤ 5e-3.
- **Ω, the whole-phantom ball:** centre (0, 0, 0), radius 0.1 m. Its distance to the farthest phantom vertex is √(0.03² + 0.04²) = 0.05 m.
- **B⁺, the displaced control:** (0.0201, 0, 0), B moved 0.1 mm radially toward the driven leg.

**Code.** One new module, `tests/validation/test_mat19_probe_coil_phantom.py` (`complex_only`), with one module-scoped fixture and no `src/` change.
- θ Functions: θ_n.x.array = snapshot_n.x.array − `T_INF`, then `scatter_forward`. **Probe θ, never T.** A ball mean of T = 20 + θ dilutes every relative comparison by θ/20; (W) would then be 20/θ ≈ 60× too lenient at n = 1.
- `probe_series` over {A, B, C, C′, D, E, Ω, B⁺} on θ_0 … θ_40, at the default `levels=3`.
- **(S) helper, test-only numpy.** For each probe and each owned cell K, the lattice-inside volume v_K(B) is computed with the same `_lattice(3)` mask as `probe._geometry`. The q_h ball mean is then Σ_K q_K·v_K(B)/Σ_K v_K(B), one `allreduce` on every rank, so numerator and denominator share one quadrature. q_K = q_h.x.array[Wq.dofmap.list[K, 0]].
- **Width artefact.** At `-n 2`, rank 0 writes `/workspace/logs/mat19/step2-n2.npz`: the series for the eight probes, the volumes, the cell count, and the sha256 of this module's bytes and of `thermal/probe.py`. Then `Barrier`.
  - At `-n 4` / `-n 1` the module loads it. It runs `pytest.fail` if the file is missing **or either hash differs**. Never recompute.
  - The hash is new against step 1: it makes rule (i)'s "npz from the committed bytes" mechanical.
- Wrap the phases in `time.perf_counter()` and print each with `flush=True`: `_run_heated`, the snapshot re-run, probe geometry, the series sums and the (S) helper.

**Anchors (asserted; every band imported).**
- **(A0) Fixture identity.**
  - `by_k[0] == "P1"`.
  - Submesh cells == 58 866 at every width (`20261001T123309Z_MAT-15-step1-W1.log:1926`).
  - \|V_sub/`V_SUB_RECORD_M3` − 1\| ≤ `V_BAND`.
  - (M) as above.
  - θ_0 reads exactly 0.0 at every probe.
- **(R0) The re-run is the gated run.** ‖θ⁴⁰_rerun − θ\*_gated‖/‖θ\*_gated‖ ≤ `E_BAND`, in `ans:9`'s (R0) form (`09_…py:246–252`). Predicted ≤ 1e-14.
- **(A) The whole-phantom probe is the volume mean.** This is step 1's (U3) on the coil-driven submesh.
  - For n = 0 … 40 and for the BE θ¹: \|P_Ω(θ) − ∫θ dx/V_sub\| ≤ `B_BAND`·\|∫θ dx/V_sub\| (n ≥ 1), and \|V_Ω,h/V_sub − 1\| ≤ `B_BAND`.
  - **Why it is exact, not banded by luck:** the 8^L equal-volume red children's centroids average to the parent centroid at every L, so each vertex weight is \|K\|/4. Σ_K \|K\|·mean(θ_K) is ∫θ for P1, and FFCx's degree-1 rule for `θ·dx` is exact for P1 too. Predicted ≤ 1e-14.
  - Printed beside: the four SPEC-time values against row 2's 0.318096 / 0.614315 / 0.892104 / 1.153314 K (`…142653Z_ANS-9-step1a-W1.log:1936`), predicted equal to the printed digits.
- **(U4)** At `-n 2`, every probe reading at n = 40 is bitwise identical on every rank (allgather, `==`). This is in-run.
- **(W) Width invariance, with the EM in the loop** (`-n 4` and `-n 1` against the `-n 2` npz).
  - Cells are equal.
  - Every entry of the series, 8 probes × n = 1 … 40, is within `WIDTH_BAND` 1e-9 relative of its `-n 2` value.
  - **Predicted ≤ 1e-10.** Backing, measured on this EM fixture (same 199 920-cell mesh, same `_solve_driven`): the P1 phantom power moves 4.9e-12 / 8.9e-13 / 1.1e-11 across `-n 2 / 4 / 1`, and raw point SARs move ≤ 2.98e-10 pairwise, on the smallest-magnitude entries (`docs/planning/chunks/POST-8.md:12–13`). G = 1/P₁ cancels the global scale, and θ at a probe is a smoothed integral of q, so it should sit at the integral class (≈ 1e-11), not the pointwise one.
  - `max_imag_ratio` is printed, predicted 0.
- **(Σ) The σ_y mirror image of the one-port drive at t\*.**
  - \|θ̄_C(t₄₀) − θ̄_C′(t₄₀)\|/θ̄_C(t₄₀) ≤ `C4_COVARIANCE_BAND` (5 %), the band `MAT-15` step 1 imported for 10 g θ ball means and step 3 applied to σ_y pairs.
  - Precondition: P1 sits at 360.0° (`…123309Z…W1.log:1924`, printed by `_run_heated`), and its ragged sheet is σ_y-symmetric to O(h) only.
  - **Predicted ≤ 0.5 %.** Step 1's C4 θ 10 g pairs read 0.020–0.049 % (`:1924–1954`). Step 3's σ_y θ pair at φ₁ ± 90° read 1.36 % on only 537 phantom cells; this fixture has 58 866.
  - **Asserted at n = 40 only.** Early-time θ is local q, and the q_h pair at ±90° can be large: step 3 read 22.7 % at the q minimum. The full C / C′ series and the q_h C / C′ pair are printed.

**Printed, not asserted (with the reason).**
- **(P4) Probe volume against point, row 4.** θ̄_j(t₄₀) / θ_point,j − 1 for A–E, against the logged row-4 values above. (R0) makes the logged θ⁴⁰ the same field to 1e-14.
  - Predicted: B in [−8 %, 0] (B sits on the steep flank of the leg hot spot; the secant from B to the row-5 maximum, 4.457206 K at r = 27.5 mm, is 0.119 K/mm). A, C, D and E within ±3 %.
  - This is the number the weekly needs to judge whether row 4's 10 % is robust to a probe volume.
- **(P5) Placement sensitivity at B.** \|θ̄_B⁺ − θ̄_B\|/θ̄_B at n = 40 is predicted ≈ 3.3e-3 (0.1 mm × 0.119 K/mm / 3.565 K), range [1e-3, 1e-2]. Equivalently ≈ 3 %/mm against the 10 % readout.
- **(S) Local calorimetric slope per probe.** ρc·θ̄_j(θ¹)/(Δt₁·q̄_j) − 1, where q̄_j is the lattice q_h ball mean.
  - **Why print only.** `MAT-15`'s `S2_BAND` 1e-4 is a bound on the **volume-mean** slope. There, testing with v = 1 removes the consistent-mass projection exactly, and the only deficit is the Robin leak (measured 5.754e-6 on P1: row 3, 1.058563285e-3 against 1.058569376e-3 K/s, `…142653Z…W1.log:1938`).
  - The physical premise holds locally: the Fourier number αΔt₁/a² = 1.435e-7 × 0.015675 / 2.5e-5 = 9.0e-5, so conduction is negligible at the probe.
  - The discrete premise fails locally. θ¹/Δt₁ → P_h q/ρc, the L2 projection of the DG0 source onto P1, and its ball mean differs from q_h's by a projection layer at the ball surface on a ball ≈ 3 cells across. That deficit does not vanish as Δt₁ → 0, and no imported band covers it.
  - Magnitude unpredicted; order estimate 1e-3–3e-2. The Ω row reads `MAT-15`'s (S2) 5.754e-6 by (A).
- The C / C′ q_h pair, the series table at n = 10 / 20 / 30 / 40 for every probe, the ball volumes, the phase clocks and the pairwise width spread at `-n 1` vs `-n 4`.

**Negative controls (rule (e)).**
- **(U3c) Ghost-included sum, asserted, factor predicted.** At `comm.size > 1` only: `_ball_average_on_cells` on Ω over owned + ghost cells (`range(size_local + num_ghosts)`) must miss (A)'s volume identity by more than `B_BAND`.
  - **Ceiling first:** at `-n 2` a cell is ghost on at most one other rank, so the ghost-included volume is ≤ 2 V_sub and the miss is ≤ 1. The floor is > 0 whenever the partition has ghosts.
  - Predicted ≈ 3e-2, range [1e-2, 1e-1], at `-n 2`. That is step 1's 8.3e-2 on a 3 072-cell box, scaled by the interface fraction (3 072/58 866)^{1/3} = 0.374. Predicted [2e-2, 2e-1] at `-n 4`.
- **(C1″) Displaced probe against (W)'s band, asserted, factor predicted.** \|θ̄_B⁺(t₄₀) − θ̄_B(t₄₀)\|/θ̄_B(t₄₀) > `WIDTH_BAND`, at every width.
  - Predicted ≈ 3.3e-3, so (W)'s 1e-9 resolves a 0.1 mm probe-geometry change by ≈ 3e6×. This is the control that (W) would see a rank-dependent probe geometry.
  - No prior measurement exists on this fixture.
- **(f″) Off-mirror pair, asserted, factor predicted.** \|θ̄_C(t₄₀) − θ̄_B(t₄₀)\|/θ̄_B(t₄₀) > `F_CONTROL_FLOOR` (= `C4_COVARIANCE_BAND`), in `MAT-15` step 3's (f′) form (θ̄(φ₁ + 90°) against θ̄(φ₁)).
  - Predicted ≈ 0.81 (1 − 0.685305/3.564854 from row 4's points), range [0.6, 0.9].
  - **Ceiling:** < 1 while θ̄_C ≤ θ̄_B, so the 5 % floor sits ≈ 16× under the prediction.
  - C against D is **not** a control: their points differ by only 1.8 % (0.685305 / 0.697444).

**Windows** (complex; foreground; Bash `timeout: 660000`; each command bare, whole and repo-relative; `-s`; no `-k` inside the quoted command; no pipes or redirects; every window ≪ 10 min, so the durable-capture idiom is not needed).
- **First:** `docker compose exec -T fem-em-solver pgrep -c python3` reads 0. Then the any-orphan sweep as one bare call:

  `scripts/testing/run_and_log.sh MAT-19-step2-cache-sweep "docker compose exec -T fem-em-solver bash -lc 'cd /root/.cache/fenics && for c in libffcx_*.c; do s=\${c%.c}; ls \$s.*.so >/dev/null 2>&1 || { echo orphan \$s; rm -f \$s.c \$s.c.cached \$s.o; }; done; echo sweep done'"`
- **Pricing basis (measured).** Per the 03:00 pathologist lesson, divide by the number of solves in the cited interval.
  - `MAT-15` step 1 W1 took **218 s** at `-n 2` (`20261001T123309Z_MAT-15-step1-W1.log:2132`). Its fixture total was 171.4 s, split as mesh + rung 61.3 s, EM 81.3 s for **four** `_solve_driven` (≈ 20 s each), and thermal 28.8 s for nine runs (`:1924,1954`). Re-gates read 208 s and 220 s (`20261002T112729Z_EX-79-mat15-regate.log:2132`, `20261004T141433Z_ANS-9-step1a-mat15-regate.log:2132`).
  - The snapshot re-run and the BE step cost ≈ 3 s: `ans:9`'s three thermal re-runs took 4.2 s (`…142653Z…W1.log:1956`).
  - The probe geometry is dominated by Ω, which materialises (m × 512 × 3) float64 over all owned cells: ≈ 0.36 GB per rank at `-n 2` and ≈ 0.72 GB at `-n 1`, a few seconds.
  - There is no `MAT-15` window at `-n 4` or `-n 1`. `POST-8` step 2 ran the same 199 920-cell EM (four drives) at 195 / 191 / 228 s for `-n 2 / 4 / 1` (`20261002T200331Z_POST-8-step2-W1.log:2137`, `…200654Z…W2.log:2341`, `…201013Z…W3.log:2034`). The width ratios are 0.98 and 1.17.
  - Measured-basis total: 218 + 218 × 0.98 + 218 × 1.17 ≈ 687 s ≈ 11.5 min.
- **W1, `-n 2`, writes the npz** (predicted 225–245 s; ≤ 345 s on item 214's ×1.4 slow box):

  `scripts/testing/run_and_log.sh MAT-19-step2-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 590 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_mat19_probe_coil_phantom.py -s -v --tb=short'"`
- **W2, `-n 4`** (predicted 220–245 s): the same command with `-n 4`, label `MAT-19-step2-W2`.
- **W3, `-n 1`** (predicted 270–320 s; ≤ 450 s slow): the same command with `-n 1`, label `MAT-19-step2-W3`. The thermal side is ≈ 2× its `-n 2` 28.8 s.
- **Order:** W1 → W2 → W3. A red W1 stops the item. Total predicted window time ≈ 13 min; the slot fits in ≈ 40 min including the module.

**Traps.**
- **No point evaluation on the submesh with an owner** (known-issues `:1372`, OPEN). The probe is an owned-cell integral. Row 4's point values are read from the `ans:9` log, never re-evaluated here.
- **θ, not T** (see Code). (W) on T is 60× too lenient.
- **Every rank enters every allreduce,** including the (S) helper on ranks with no candidate cells (item 227 deadlocked on a rank-0-only collective). (U3c) at `-n 2` is the check on the probe path.
- **Ghosts.** Iterate `range(size_local)` for owned cells. Ghost **dof values** are needed: snapshots copy `x.array` after `scatter_forward`.
- **Complex build.** `np.real` before every `float`. Make no ordering comparison on a UFL operand: the probe's indicator is numpy-side, and this module compiles no new indicator form. `_run_heated`'s own ball forms are unchanged.
- **Snapshot indices are full steps.** Snapshot 1 is after both Rannacher half-steps, at t = Δt; times are n·Δt.
- **No cross-run comparison is bitwise** (`OPS-43`; complex MUMPS is not bit-reproducible across processes). (W) is relative 1e-9; only (U4) is bitwise, and it is in-run.
- **Stale npz.** `logs/` persists between sessions. The hash check in the npz is the guard. Never delete-and-recompute at `-n 4` / `-n 1`.
- **Orphan `.c`:** any stem without its `.so`, not only 0-byte (item 160's 44 KB partial). Sweep only when `pgrep -c python3` reads 0.
- **Rule (i):** no edit to the module after W1. An edit re-runs W1–W3, and the hash enforces it.
- **The SPEC and `ans:9` are not this item's to edit.** Changing row 4 to a ball average is the weekly's decision and the operator's export.

**Scope.** Covered: `MAT-15` P1's fixture (ragged sheets, degree-1 EM, P1 thermal, one mesh), drive P1 only, the SPEC's five centres plus C′, one stated volume (a = 5 mm), recorded times only. **Not covered:**
- `ANS-9`'s row 4 as specified (a point), or any adjudication, AED reading or AED number;
- a real probe's sensing volume, thermal mass or lead artefact;
- exact rectangles or degree 2: this is a record on the SPEC's fixture, not a new one;
- an h-statement on the probe readings;
- the copper / tuned fixture (`MAT-15` step 3);
- tissue, perfusion, ASTM F2182;
- the submesh owner-rule fix (known-issues stays OPEN);
- any edit to `ans:9`'s metrics or COMPARISON.

**Status it can move.**
- `MAT-19` is ✅ on step 1 and stays ✅ whatever this step reads, unless (A) is red (see below).
- **On green:** the §7 row gains a **step 2 record** sentence: the probe on the coil-driven phantom at the SPEC centres; (A), (W) with the EM in the loop and (Σ) asserted; and the (P4) / (P5) sensitivities printed. **Rule (j):** the reviewing daily adds the Done-when clause "step 2 (§10 P8): identities green ⇒ a record, no glyph", in the pattern of `MAT-15` step 3.
- §10 **P8** is marked DONE.
- **(8c):** the weekly gets the measured probe-volume and placement sensitivity of row 4 at the 10 % AGREE level. It can then decide, **before the operator's AED run**, whether `ANS-9`'s SPEC should also ask for a ball-averaged row 4 under the `POST-8` ruling ("public comparison rows quote cell or mass averages"). That decision and any SPEC edit are the weekly's. If the AED run has already been made, the weekly judges from the private side whether a re-export is worth asking for. This brief does not read that side.

**Negative result.** In every case: report the measurement, park on `attempt/MAT-19-step2-*`, mark the item BLOCKED in the same commit (rule (d)) and journal it in attempts.md. Never widen a band, and never change a, L or the centres after reading.
- **(A0) or (R0) red:** the fixture moved under `MAT-15`'s record. Stop; a `record-reconciler` ruling is needed.
- **(A) or (U4) red:** the probe's ownership contract fails on a coil-driven submesh partition. Open a known-issues entry. The review re-reads `MAT-19`'s ✅, because the row's Done-when is width invariance as `POST-8`.
- **(W) red with (A0) / (R0) / (A) green:** the complex solve's width round-off reaches 1e-9 in a probe reading, contradicting `POST-8`'s ≤ 2.98e-10 on this EM. Print the 8 × 40 spread, P₁ and G per width, and stop. The review rules whether a probe readout on an EM-driven series carries a width caveat. The band is not moved.
- **(Σ) red:** print the q_h C / C′ pair and the series. If the q_h C / C′ ball pair is also beyond 5 %, the ragged P1 sheet breaks σ_y. Record it in known-issues as a ragged-sheet consequence. Stop.
- **A control does not miss** ((U3c), (C1″) or (f″)): stop. That contradicts (A), (W) or (Σ) respectively.
- **A slow-box kill at 590:** BLOCKED. W1's record stands, and the unblock is a re-run on a box whose `MAT-15` EM reads ≤ 30 s per solve (the 03:00 re-price).

No `PREDICTED RED:` line. A red anywhere stops the slot.

### §9 item 258 — 2026-10-05 (chain P8, daily provisional; re-pose of item 253, attempt 2 of 4)

*Drafted by a general-purpose agent for the 2026-10-05 18:00 daily review; the review reviews it.* Licence: §9 rule (6), the re-pose of a blocked step inside the four-attempt family cap. The weekly ratifies or reverts it. **Unchanged:** the step's target (a whole-phantom probe = ⟨θ⟩ at `B_BAND`; width invariance at `WIDTH_BAND` 1e-9 with the EM in the loop; the σ_y mirror pair within `C4_COVARIANCE_BAND` at t\*), every band, the centres, a, L and the family count. **Changed:** only *where* the planted control (U3c) and the "partition" reading of (A) / (U4) are asserted, with the rule pre-registered below, plus one optional speed-up of a printed helper.

`MAT-19` step 2, re-pose of item 253: step 1's probe, unchanged, on `MAT-15` P1's coil-driven gel phantom at `ANS-9`'s SPEC centres (implementer; heavy by ceiling (`timeout -k 30 590`), `-n 2` / `-n 4` / `-n 1`, plus `-n 8` only by the decision table below, complex; ≈ 9 measured / ≈ 45 predicted min; depends on nothing). Family `MAT-19` step 2: attempt 2 of 4.

**What attempt 1 measured (read first).** Journal: `docs/testing/attempts.md` 2026-10-05T18:34Z. Code and logs are parked on `attempt/MAT-19-step2-20261005T1834Z` (`e4ffe9a4`).
- W1 `-n 2` (`20261005T183423Z_MAT-19-step2-W1.log`, 326 s, Status 1) was red **only** on (U3c). The ghost-included Ω volume missed by 9.770e-15 and the value by 1.577e-14, against "> `B_BAND` 1e-10" (`:1977`).
- Every other W1 reading was green (`:1956–1980`):
  - (A0) 58 866 cells, V_sub rel 3.447e-13;
  - (R0) 2.525e-15;
  - (A) 2.386e-14, BE 7.351e-15, volume 9.770e-15, and row 2 reproduced to the printed digits;
  - (U4) bitwise;
  - (Σ) 5.919e-4;
  - (C1″) 4.398e-3;
  - (f″) 0.8018.
- Phase clocks (`:1980`): `_run_heated` 168.2 s (EM 144.9 s for four solves, `:1955`), snapshot re-run + BE 3.1 s, probe geometry + sums 17.0 s, (S) helper 59.2 s, fixture total 258.7 s.
- **The cause was measured** by the diagnostic `mat19_step2_ghost_probe.py` at `-n 2` (`20261005T184012Z_MAT-19-step2-ghost-diag.log:1871–1873`, 82 s):
  - parent owned 101 934 / 97 986 cells, ghosts 869 / 883;
  - **tag-3 submesh owned 58 866 / 0, ghosts 0 / 0.**
- So at `-n 2` the phantom is not partitioned:
  - (U3c) iterates `range(size_local + num_ghosts)` on the **submesh** index map, so it has nothing to double-count;
  - (A) and (U4) are true there, but they do not exercise a submesh partition.

**The correct "something to double-count" condition (checked against the parked module, lines 193–198).** (U3c) calls `_ball_average_on_cells(θ⁴⁰, Ω, arange(size_local + num_ghosts))` on the tag-3 submesh `sub.topology.index_map(3)`. It then allreduces. Every ghost cell on a rank is owned by another rank, and Ω contains every phantom cell. The planted volume miss is therefore exactly Σ_ranks (volume of that rank's submesh ghost cells) / V_sub, up to round-off: for Ω every lattice point of every cell is inside, so each cell weighs its exact \|K\|. It follows that:
- **(U3c) applies iff G_w := Σ_ranks `num_ghosts`(submesh) > 0.** Ownership by ≥ 2 ranks is neither necessary nor sufficient for this. A rank can hold ghosts of phantom cells while owning none. Whether the parent's ghost phantom cells reach the submesh at all depends on whether `cell_tags.indices` covers parent ghost cells; `_tag_submesh` keeps tagged cells over `size_local + num_ghosts`, `test_mat15_coil_driven_heating_phantom.py:106–121`.
- **The partition reading of (A) / (U4) applies iff O_w := #{ranks with submesh `size_local` > 0} ≥ 2.**

These two conditions are what this item pre-registers. The review may tighten them to the conjunction "O_w ≥ 2 and G_w > 0" for (U3c) if it prefers; see the open question in the return.

**Pre-registered applicability rule (the slot does not decide).**
- **D4 comes first.** Run the parked diagnostic at `-n 4` (window D4 below). It prints the per-rank submesh owned / ghost counts, from which O_4 and G_4 are read.
- **D8 runs only if G_4 = 0 or O_4 < 2.** It is the same diagnostic at `-n 8`.
- **Decision table.** It is fixed now and is applied to the D logs as printed:

  | D4 reads | D8 | (U3c) asserted at | partition reading of (A) / (U4) at | W4 (`-n 8`) window |
  | --- | --- | --- | --- | --- |
  | G_4 > 0 and O_4 ≥ 2 | not run | `-n 4` (W2) | `-n 4` (W2) | not run |
  | G_4 > 0, O_4 < 2 | run | `-n 4` (W2) | `-n 8` if O_8 ≥ 2, else not exercised | run iff O_8 ≥ 2 |
  | G_4 = 0, O_4 ≥ 2 | run | `-n 8` if G_8 > 0, else **inapplicable** | `-n 4` (W2) | run iff G_8 > 0 |
  | G_4 = 0, O_4 < 2 | run | `-n 8` if G_8 > 0, else **inapplicable** | `-n 8` if O_8 ≥ 2, else not exercised | run iff G_8 > 0 or O_8 ≥ 2 |

- **"Inapplicable" / "not exercised" is a record, not a red.** In that case the item closes on the remaining anchors. The chunk-file outcome and the journal state, in these words: "(U3c) inapplicable on this fixture: the tag-3 submesh carries no ghost cells at `-n 2 / 4 / 8` (G_w = 0, `<D logs>`)". For the partition: "(A) / (U4) asserted on a submesh owned by one rank at every measured width; the submesh partition is not exercised on this fixture". No band moves, and no other width is tried in-slot.
- **Predicted (no prior measurement at `-n 4`; rule (e), predicted, printed):**
  - **O_4 ≥ 2 by counting.** At `-n 2` the parent split was 101 934 / 97 986, a 2 % imbalance. At `-n 4` each rank owns ≈ 49 980 × [0.97, 1.03] < 58 866 cells, so the phantom cannot fit on one rank unless the partitioner's imbalance exceeds ≈ 17 %.
  - **G_4 is genuinely unmeasured.** The `-n 2` reading of 0 / 0 is consistent with either answer to the `cell_tags`-covers-ghosts question, because the rank-0/1 cut did not touch the phantom. The predicted G_4 is in [500, 6000] if ghosts propagate and 0 if they do not.
  - The expected row is therefore row 1 or row 3.

**Carried verbatim from item 253's brief** (above; read it there, it is not restated):
- "Why this step";
- "The fixture (imported, never restated)", including the Imports list and the no-edit clause on `src/`, the `MAT-15` module, `ans:9`, its metrics and its COMPARISON;
- "The probes (pinned)", including (M), Ω and B⁺;
- "Code", except the edits E1–E6 below;
- the anchors (A0), (R0), (A) and (Σ) as stated;
- (W) as stated, asserted at `-n 4` and `-n 1` against the `-n 2` npz at `WIDTH_BAND`;
- "Printed, not asserted", meaning (P4), (P5), (S) and the tables;
- the negative controls (C1″) and (f″) as stated. Both misses are asserted. Their factors are now **predicted equal to attempt 1's measured 4.398e-3 and 0.8018** (W1 `:1972`, `:1975`), printed beside;
- every trap in item 253's list;
- the scope list;
- every negative-result disposition in item 253's list. The (U3c) line, "A control does not miss", is amended below.

**Landing the parked code (by path; never a record file by path).**
1. `git checkout attempt/MAT-19-step2-20261005T1834Z -- tests/validation/test_mat19_probe_coil_phantom.py mat19_step2_ghost_probe.py docs/testing/logs/20261005T183255Z_MAT-19-step2-cache-sweep.log docs/testing/logs/20261005T183423Z_MAT-19-step2-W1.log docs/testing/logs/20261005T184012Z_MAT-19-step2-ghost-diag.log`
2. `git mv mat19_step2_ghost_probe.py scripts/probes/mat19_step2_ghost_probe.py`. Keep its bytes identical: no edit, so its `-n 2` log stays a record of the same bytes. The repo root holds no tracked `.py`, and `scripts/probes/` is where the project's probes live. The `-n 2` log's command names the old root path; the journal says it was moved unedited.
3. Append the attempt's three test-results rows **by hand** to `docs/testing/test-results.md` (copy them from `git show e4ffe9a4 -- docs/testing/test-results.md`). `git diff --stat` must show only insertions on that file.
4. The stale `/workspace/logs/mat19/step2-n2.npz` from attempt 1 carries the old module hash. W2 / W3 / W4 refuse it by construction. Never delete it by hand; W1 overwrites it.

**Code: edits to the parked module, all made before W1 (rule (i)).** Every constant, centre, band and import stays as parked.
- **E1, docstring.** Cite item 258 and the applicability rule: (U3c) is asserted iff G_w > 0, and the partition reading is iff O_w ≥ 2.
- **E2, the partition table in the fixture.** Every rank enters `comm.allgather((imap.size_local, imap.num_ghosts))` on `sub.topology.index_map(3)`, at every width including `-n 1`. Compute O_w and G_w from it. Print one line, `MAT19-2 (part) -n w: sub owned [...] ghosts [...] O_w G_w; partition exercised: <O_w ≥ 2>`.
- **E3, (U3c).** At `comm.size > 1` the ghost-included Ω reading is computed and printed **always**, as parked. **Asserted (volume miss > `B_BAND`) iff G_w > 0.** When G_w = 0 it prints `MAT19-2 (U3c) inapplicable at -n w: G_w = 0` and asserts nothing. This makes the module's assertion mechanical at every width. The item-level table above says which width must carry it.
- **E4, (U4).** Assert the in-run bitwise allgather at **every** `comm.size > 1`; it was only at `== 2`. It costs nothing.
  - Predicted True at `-n 4` / `-n 8`: it was measured True at `-n 2` on this fixture (`:1978`). Each reading is a quotient of identical allreduced sums, and float addition is commutative, so the power-of-two pairwise exchanges give the same bits on every rank.
- **E5, (W) widths.** Assert (W) (cells equal; series ≤ `WIDTH_BAND`) at `comm.size ∈ {4, 1}` only, as the row's `-n 1 / 2 / 4` reads. At `-n 8`, (W) against the `-n 2` npz is **printed** with the hash guard intact, and no `-n 8` npz is written. The `-n 4` npz write and the `-n 1` vs `-n 4` print stay as parked.
- **E6, the optional speed-up of the (S) helper (printed-only quantity), applied before W1.** In `_lattice_q_mean`, pre-filter `owned` by bounding radius exactly as `thermal/probe.py`'s `_geometry` does (`probe.py:115–119`): take `X = xdof[V.dofmap.list[owned]]`, `cen = X.mean(1)`, `rho = sqrt(max ‖X − cen‖²)`, and keep `‖cen − c‖ − rho ≤ a`. Then run the parked chunk loop over the kept cells.
  - **Why it cannot change a reading beyond round-off:** every lattice point is a convex combination of the cell's vertices, so ‖P − c‖ ≥ ‖cen − c‖ − ρ_K > a for a dropped cell. Its v_K is exactly 0, and dropping it changes only the pairwise-summation grouping of `np.sum`.
  - For Ω (r = 0.1 m against a farthest vertex at 0.05 m), every cell is kept and the path is byte-identical.
  - **Predicted:** the (S) phase drops from 59.2 s to ≈ 8 s (Ω's ≈ 1/8 share), and every (S) value equals attempt 1's (`:1976`) to the printed 5 digits.
  - If the slot does not apply E6, nothing else changes; the pricing below then adds ≈ 50 s per window.

**Anchors and controls, re-posed parts (everything else carries verbatim, see above).**
- **(U3c), the planted ghost control, asserted where the table says, factor predicted.** At the width the table names: the ghost-included Ω volume misses (A)'s identity by > `B_BAND`. The miss is asserted; the factor is predicted.
  - **Floor:** one ghost cell is ≈ 3.84e-9 / 2.26e-4 ≈ 1.7e-5 of V_sub, so any G_w ≥ 1 clears 1e-10 by ≥ 1e4× unless a ghost is a sliver below 1e-5 of the mean cell volume.
  - **Ceiling:** ≤ w − 1, since each cell is ghost on at most w − 1 other ranks.
  - **Predicted factor:** item 253's range [2e-2, 2e-1] at `-n 4` is kept and printed. The sharp prediction, printed beside, is G_w / 58 866 within ×[0.5, 2], with G_w from the D log at that width. The phantom mesh is uniform at h_p = 2.5 mm, so ghost cells are ≈ mean-sized.
  - The value miss is printed.
- **The partition reading of (A) / (U4).** (A) stays asserted at every width as in item 253. (U4) is asserted at every width > 1 (E4). The reading "(A) and (U4) hold on a **coil-driven submesh partition**" is licensed only from the window at the width the table names, whose own (part) line shows O_w ≥ 2.
- **Cross-check, D against the in-run partition (printed; a mismatch stops the item).** At the width that carries (U3c) or the partition reading, the window's (part) line must match its D log's per-rank owned / ghost counts. The mesh and partitioner path is the same `build_four_port_sweep(PHANTOM_RESOLUTION_1G_RUNG)`.

**Windows** (complex; foreground; Bash `timeout: 660000`; each command bare, whole and repo-relative; `-s`; no `-k` inside the quoted command; no pipes or redirects; every window ≤ 590 s container-side, so the §5.1 durable-capture idiom is not needed).
- **First:** `docker compose exec -T fem-em-solver pgrep -c python3` reads 0. Then the any-orphan sweep as one bare call:

  `scripts/testing/run_and_log.sh MAT-19-step2-a2-cache-sweep "docker compose exec -T fem-em-solver bash -lc 'cd /root/.cache/fenics && for c in libffcx_*.c; do s=\${c%.c}; ls \$s.*.so >/dev/null 2>&1 || { echo orphan \$s; rm -f \$s.c \$s.c.cached \$s.o; }; done; echo sweep done'"`
- **Pricing basis (measured: attempt 1's W1 at `-n 2`, 326 s).** The fixture took 258.7 s, so environment tests, collection and teardown are ≈ 67 s. The width ratios are `POST-8` step 2's 0.98 (`-n 4`) and 1.17 (`-n 1`) on the same EM (item 253's basis). This box's EM read 144.9 s for four solves (`:1955`), above item 253's ≤ 30 s / solve unblock bar, and W1 still finished at 326 s, so the 590 ceiling holds ≈ 1.8× headroom. There is no `-n 8` window on this fixture.
- **D4, ownership at `-n 4`** (predicted 70–120 s; the `-n 2` diagnostic took 82 s, mesh build dominated):

  `scripts/testing/run_and_log.sh MAT-19-step2-a2-D4 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace timeout -k 30 590 mpiexec -n 4 python3 scripts/probes/mat19_step2_ghost_probe.py'"`
- **D8, ownership at `-n 8`.** Run only per the table (predicted 70–150 s): the same command with `-n 8`, label `MAT-19-step2-a2-D8`.
- **Then** make the E1–E6 edits. No edit after W1.
- **W1, `-n 2`, writes the npz.** This is the rule (i) re-run on the landed bytes (predicted 240–330 s; ≈ 275 s with E6, ≈ 326 s without):

  `scripts/testing/run_and_log.sh MAT-19-step2-a2-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 590 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_mat19_probe_coil_phantom.py -s -v --tb=short'"`

  - Expected (part) line: O_2 = 1, G_2 = 0, so (U3c) prints inapplicable at `-n 2` and is not asserted.
  - If the in-run line at `-n 2` instead shows G_2 > 0, the module asserts (U3c) there too. That is harmless, and it is recorded.
- **W2, `-n 4`** (predicted 230–330 s): the same command with `-n 4`, label `MAT-19-step2-a2-W2`. It asserts (W), and (U3c) / (U4) per E3 / E4.
- **W3, `-n 1`** (predicted 260–360 s; the thermal side at `-n 2` was already serial on rank 0, so only the EM ratio applies): the same command with `-n 1`, label `MAT-19-step2-a2-W3`. It asserts (W).
- **W4, `-n 8`.** Run only per the table (unmeasured on this fixture; predicted 230–400 s): the same command with `-n 8`, label `MAT-19-step2-a2-W4`. (W) is printed there (E5); (U3c), (U4) and (A) are asserted per E3 / E4.
- **Order:** pgrep → sweep → land → D4 → [D8] → E1–E6 → W1 → W2 → W3 → [W4].
  - A red W1 stops the item. W2 / W3 / W4 need W1's npz from the same bytes, and the hash enforces it.
  - **Total window time:** ≈ 16.5 min on table row 1 (81 + 82 + 275 + 260 + 290 s), ≈ 23 min with D8 + W4. The slot fits its hour with the edits.

**Traps (added to item 253's list).**
- **Rule (i) and the npz hash.** E1–E6 all land before W1. Any edit afterwards re-runs W1 and everything after it. Check file mtimes against the W1 log start before committing.
- **Every rank enters E2's allgather** and every allreduce of the (S) helper, including ranks that own no phantom cells. At `-n 2` rank 1 owns none (item 227's rank-subset deadlock).
- **The diagnostic is not a gate.** It is moved unedited, never edited, and runs through the harness like any window. Its `GHOSTPROBE` lines are the D reading.
- **Never pick the width after seeing W-window readings.** The table is applied to the D logs only.
- **The bash guard** denies any `docker` / `run_and_log.sh` call that is not bare and whole. The spawn prompt carries the every-item harness sentence verbatim.
- **Landing.** Record files are never checked out by path (the `8d4cf58` trap). Only the module, the diagnostic and the three logs are.

**Scope (added to item 253's list).** Not covered:
- a submesh partition at widths other than those the table runs;
- any statement about DolfinX's ghost propagation through `create_submesh` beyond the counts printed;
- (W) at `-n 8` as a gate (printed only).

**Status it can move** (restated, because the record wording now depends on the table row).
- `MAT-19` is ✅ on step 1 and stays ✅ whatever this step reads, unless (A) or (U4) is red (item 253's disposition: a known-issues entry, and the review re-reads the ✅).
- **On green:** the §7 row gains a **step 2 record** sentence. It reads: the probe on the coil-driven phantom at the SPEC centres; (A) and (U4) asserted, "on a submesh partition at `-n w`" only if the table's partition width ran with O_w ≥ 2, else "on a phantom owned by one rank (partition not exercised)"; (W) with the EM in the loop at `-n 4 / -n 1`; (Σ) asserted; (U3c) "missed by x at `-n w`" or "inapplicable on this fixture (G_w = 0 at `-n 2 / 4 / 8`)"; the (P4) / (P5) sensitivities printed. Elapsed times and logs are cited.
- **Rule (j):** if the row still lacks it, the reviewing daily adds the Done-when clause "step 2 (§10 P8): identities green ⇒ a record, no glyph", as item 253's brief asked; it is not in the row as of this draft.
- §10 **P8** is marked DONE. (8c) is as item 253's brief states.

**Negative result (amended; item 253's list otherwise carries).** In every case: report the measurement, park on `attempt/MAT-19-step2-*`, mark the item BLOCKED in the same commit (rule (d)), and journal it in attempts.md as attempt 2 of 4. Never widen a band, and never change a, L, the centres or the table after reading.
- **(U3c) asserted (G_w > 0) and does not miss:** ghost cells were present and not double-counted. That contradicts the probe's cell loop. Stop, print the (part) line and both volumes, and open a known-issues entry. The review re-reads step 1's (U3c).
- **(U3c) inapplicable per the table:** a record, not a red. Continue on the remaining anchors.
- **D / in-run mismatch** at the width that carries (U3c) or the partition reading: the partition is not reproducible across runs. Stop BLOCKED with both counts. The review rules; the slot does not reassign the width.
- **D4 or D8 fails to run** (import error or kill): stop. W windows do not start without the D reading that sets the table.
- **(U4) red at `-n 4` / `-n 8`:** item 253's "(A) or (U4) red" disposition.
- **A slow-box kill at 590:** BLOCKED. Every window that finished stands as a record. The unblock is a re-run, never a longer timeout.

No `PREDICTED RED:` line. A red anywhere stops the slot.

**Outcome (2026-10-06 01:11Z, 19:30 CDT slot, third item): ✅ GREEN, attempt 2 of 4.**
- **Landed** the parked module, the diagnostic (moved unedited to `scripts/probes/`) and attempt 1's three logs by path from `e4ffe9a4`. Attempt 1's three test-results rows were appended by hand. Sweep `20261006T005419Z_MAT-19-step2-a2-cache-sweep.log`.
- **D4** (`20261006T005546Z_MAT-19-step2-a2-D4.log:1871–1875`, 68 s): sub owned 23 401 / 35 465 / 0 / 0, ghosts 701 / 677 / 0 / 0, so O_4 = 2 and G_4 = 1 378. That is **table row 1**: (U3c) and the partition reading go on `-n 4`, with no D8 and no W4. The row-1/row-3 prediction held, and ghosts do propagate through `_tag_submesh` at `-n 4`.
- **E1–E6 applied before W1.** No edit came after W1, and the npz hash held across W2 / W3.
- **W1 `-n 2`** (`20261006T005727Z_MAT-19-step2-a2-W1.log`, 275 s, 13 passed):
  - (part) [58 866, 0] / [0, 0], O 1, G 0, so (U3c) printed "inapplicable at -n 2: G_w = 0";
  - (A0) 58 866 cells, V_sub rel 3.447e-13; (R0) 2.531e-15;
  - (A) 2.190e-14, BE 2.654e-15, volume 9.770e-15, with row 2 reproduced to the printed digits;
  - (U4) True;
  - phase clocks: `_run_heated` 175.0 s; (S) helper 6.8 s with E6, against 59.2 s before. The (S) values equal attempt 1's (`:1976`) to all printed digits.
- **W2 `-n 4`** (`20261006T010208Z_MAT-19-step2-a2-W2.log`, 248 s, 13 passed):
  - (part) matches D4 exactly;
  - (A) 6.825e-15, volume 1.332e-15; (U4) True;
  - **(U3c) volume miss 2.6047e-2** (asserted > 1e-10), value miss 3.080e-3. That is inside the predicted [2e-2, 2e-1]. The sharp prediction G_4/58 866 = 2.341e-2 is within ×1.11;
  - **(W) 3.714e-11** (band 1e-9, predicted ≤ 1e-10), volumes 4.4e-16.
- **W3 `-n 1`** (`20261006T010623Z_MAT-19-step2-a2-W3.log`, 290 s, 13 passed): (A) 1.780e-14; **(W) 5.735e-11**; the pairwise `-n 1` vs `-n 4` print is 7.287e-11.
- **Width-independent readings (all three windows):**
  - (Σ) 5.919e-4;
  - (C1″) 4.398e-3 and (f″) 0.8018, equal to attempt 1;
  - (P4) A +7.18e-3, B −1.238e-2, C +1.807e-2, D −4.35e-3, E +2.84e-3. B is inside its predicted [−8 %, 0], and A / C / D / E are inside ±3 %;
  - (P5) 4.398e-3, inside [1e-3, 1e-2];
  - the q_h C / C′ ball pair, printed: 9.55e-3.
- **Record wording:** (A) and (U4) hold on a submesh partition at `-n 4`; (U3c) misses by 2.6e-2 at `-n 4`. No band moved.

**Spawn prompt (for the slot).** "Execute §9 item 258, `MAT-19` step 2, re-pose of item 253 (attempt 2 of 4), per `docs/planning/chunks/MAT-19.md` `### §9 item 258` and the parts of `### §9 item 253` it carries verbatim. The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Every window runs foreground with Bash `timeout: 660000`. Land the module, the diagnostic and the three logs by path from `attempt/MAT-19-step2-20261005T1834Z`; append record rows by hand. Apply the decision table to the D logs only. Never widen a band. A red outside the table's 'inapplicable' rows stops the item."

### §9 item 265 — 2026-10-06 (chain P10)

**Outcome (2026-10-06 09:00 CDT slot, attempt 1 of 4): BLOCKED — (Σ) red at W1; parked on `attempt/MAT-19-step3-20261006T1405Z` (`a5b3c631`).** W1 `20261006T140518Z_MAT-19-step3-W1.log` (168 s, `-n 8`, summed ru_maxrss 12.96 GiB; clocks mesh+tie 48.6 s, EM drive 41.5 s, transient 0.8 s, probe 0.1 s). **Red: (Σ) \|C − C′\|/C at t\* = 5.7365e-2 against `C4_COVARIANCE_BAND` 5e-2.** Per the brief's named disposition, printed: q_h ball pair C 340.33 / C′ 264.85 W/m³ (22.18 %, cf. item 236's 22.7 % 10 g q_h pair), ball volumes C 5.0774e-7 / C′ 5.1643e-7 m³ (V_h/nominal − 1: C −3.03e-2, C′ −1.37e-2). Green before it: (A0) 80 798 cells, tag-3 537 (record :1926), (Φ) \|sin φ₁\| = 0, (M) true, snapshot 0 = 0.0; (Src) (W) 1.116e-12 (`NODE_S_RTOL` 1e-9), (P) 7.393e-13 (1e-6); (A) max 4.549e-16, BE 0, V_Ω/V_sub − 1 = 4.441e-16 (`B_BAND` 1e-10), **on a submesh partition at `-n 8`** (owned 100/29/45/57/130/30/79/67, ghosts 40/19/25/25/51/19/28/22, O_w 8, G_w 229); (U4) bitwise; (U3c) asserted (G_w > 0) and missed by 0.5474 (predicted G_w/537 = 0.4264, within ×[0.5, 2]). Computed but not reached as asserts (the (Σ) assert stopped the test): (C1″) 1.9498e-3 (predicted [1e-4, 1e-2]), (f″) 0.5474 (predicted [0.1, 1.0]). Printed: ⟨θ⟩(t\*) 0.687806 K per accepted watt vs item 236's 0.687806 (rel 1.6e-8, 6-digit record); (P4′) ball/point − 1 at A +2.35 %, B +0.02 %, C +0.10 %, D −1.20 %, E −2.04 % (A and E above the predicted ≤ 2 %); ratio (i) (mixed normalisation) A 0.6926, B 0.3195, C 0.7296, D 1.5990, E 0.6799, C′ 0.6874, Ω 0.5964; ratio (ii) (per deposited watt at both) A 1.1989, B 0.5530, C 1.2631, D 2.7680, E 1.1769, C′ 1.1899, Ω 1.0324. W2 / W3 not run, so (W) is not measured. Other logs: `20261006T140242Z_MAT-19-step3-cache-sweep.log` (no orphans), `20261006T140407Z_MAT-19-step3-import.log` (knob-unset collection smoke, killed at its own 60 s ceiling during `tests/environment`; the module collected, `pgrep` 0 afterwards). **For the review:** the brief's (Σ) disposition applies: on ≈ 15 mm cells the 5 mm probe pair inherits the unsymmetric phantom mesh (q_h pair 22 %). Re-pose on a phantom-refined exact hole mesh (a new fixture with its own C_mode), not in-slot. No band moved.

*Drafted by a general-purpose agent for the 2026-10-06 03:00 daily review; the review reviews it.* Licence: §10 chain P10 (grep "P10. `MAT-19` step 3" in PROJECT_PLAN.md), quoted, not changed: "the probe series on the tuned copper fixture at 64 MHz (P5's source, degree 2, exact rectangles) at `ANS-9`'s six centres, 40 steps. Asserted: P8's identities (whole-phantom probe = ⟨θ⟩ at `B_BAND`, width invariance at 1e-9 with the EM in the loop, the σ_y mirror pair at `C4_COVARIANCE_BAND`); printed, ball/point at A–E and the 64 MHz / 10 MHz ratio of the probe rises per accepted watt. The first probe-level thermal figure on a tuned coil; a record. Target: identities. Serial on nothing." The 2026-10-02 sheets-and-κ ruling binds (exact rectangles, κ off). §10 ruling P8 (item 258's ghost / width rule) binds. **No band, tolerance, Done-when or validation target is set or moved.** Every band is imported by name.

`MAT-19` step 3: step 2's probe series (item 258), unchanged, on `MAT-15` step 3's tuned copper source (item 236).
- **Executor:** implementer. One new test module, `tests/validation/test_mat19_probe_tuned_copper.py`. There is no `src/` change. Do not edit `thermal/probe.py`, item 236's module (`test_mat15_tuned_copper_heating.py`), item 258's module (`test_mat19_probe_coil_phantom.py`), step 1's `MAT-15` module, the `ans:7` / `ans:9` scripts, or any `ans7_exact` test. An edit to any of them owes those gates' windows under rule (a).
- **Tier:** heavy, `-n 8` / `-n 4` / `-n 2`, complex build.
- **Cost:** ≈ 17 measured / ≈ 50 predicted min. That is 15 fixed + ≈ 12 code + sweep 0.5 + W1 ≈ 3.3 + W2 ≈ 5 + W3 ≈ 8 + ≈ 6 spare.
- **Depends on nothing on `main`.** It needs T21's two gitignored JSONs, `logs/ans7_exact/tune_p2.json` and `tune_inmodel_p2.json`, which item 236's `_read_t21` reads. If either is absent, stop and take the next item. Never re-create them.

**The fixture (imported, never restated).**
- **The source** is item 236's, exactly. Use the `_exact_geo` ×1 hole mesh (80 798 cells), Leontovich copper, and C_mode,Cu read from `tune_p2.json`. Apply Z_C on P2–P4 and C₁ across P1. Drive one P1 at 64.0 MHz through `_route_a_solve(..., degree=2)`, with G = 1 W / P_acc.
- **Reuse before you copy.** Read `test_mat15_tuned_copper_heating.py:142–310` (`_run`) first.
  - If `_run()` returns the submesh, Robin tags, q_h, G, P_acc and η, call it.
  - If it does not, rebuild the same steps in the same order from the **imported** helpers: `_read_t21`, `_fixture_tie`, `_exact_geo`, `_route_a_solve`, the `TH-14` (b) power block with the Z_C-replaced specs, `_tag_submesh`, `_robin_tags`, `_cell_average_dg0_tag` and `_to_sub`. Disclose the mapping.
  - Never edit item 236's module to make it return more.
- **The series** is item 258's: `_transient(sub, ft, q_h, dt=T_STAR / N_MAIN, n_steps=N_MAIN, theta=0.5, startup_half_steps=2, snapshot_steps=range(N_MAIN + 1))`, `_be_one`, then `probe_series` over θ = T − `T_INF`.
- **The probes** are item 258's `_probe_set()`, imported: `ANS-9`'s SPEC centres A–E plus C′ (σ_y image of C), so six centres, plus Ω (r = 0.1 m) and the control ball B⁺. Keep a = `A_PROBE` 5 mm.
- **Imports:** `B_BAND`, `C4_COVARIANCE_BAND`, `E_BAND`, `F_CONTROL_FLOOR`, `DT_1`, `N_MAIN`, `T_STAR`, `T_INF`, `RHO_C`, `_transient`, `_be_one`, `_allsum` and `_say`, from `test_mat15_coil_driven_heating_phantom` exactly as item 258's module imports them (`test_mat19_probe_coil_phantom.py:56–74`). Also import `WIDTH_BAND` from `test_mat19_virtual_temperature_probe`, and `_probe_set`, `A_PROBE`, `OMEGA_RADIUS` and `_lattice_q_mean` from item 258's module. Restate none of them.
- **(Φ) is reused.** Item 236 measured φ₁ = 0.000000000000° on this mesh (`20261005T040135Z_MAT-15-step3-W1.log`), so C′ = (x, −y, z) is the drive's σ_y image here as on step 2's fixture. Re-assert item 236's (Φ) in-run.
- **Opt-in knob** `FEM_EM_MAT19_STEP3=1`, in item 236's pattern. When it is unset, every test skips. Any other value raises `ValueError`.
- **npz:** `/workspace/logs/mat19/step3-n8.npz`. Guard it with the SHA-256 of the new module, `thermal/probe.py` and item 236's module.

**Anchors (asserted; bands imported).**
- **(A0)** `_assert_tie` (80 798 cells), (Φ) \|sin φ₁\| ≤ item 236's `PHI_BAND`, and the tag-3 cell count == 537. That count is item 236's printed record (`…040135Z_MAT-15-step3-W1.log`; cite the line). Snapshot 0 reads exactly 0.0, and every probe ball lies inside the phantom ((M), item 258's test).
- **(Src)** item 236's (W), \|S11 − S11_ref\|/\|S11_ref\| ≤ `NODE_S_RTOL`, and (P) ≤ `DISCRETE_IDENTITY_RTOL`, re-asserted on this run's field. This ties the source to item 236's gated run.
- **(A)** whole-phantom probe P_Ω = ∫θ/V_sub at n = 1 … 40 and for the BE θ¹, plus V_Ω/V_sub − 1, all ≤ `B_BAND`.
- **(U4)** the rank rows agree bitwise at n = 40, at every width > 1.
- **(W) width invariance with the EM in the loop.** The `-n 4` and `-n 2` series against the `-n 8` npz: cells equal, and max rel ≤ `WIDTH_BAND`. The hash guard stays as item 258's.
- **(Σ)** \|θ̄_C − θ̄_C′\|/θ̄_C at t\* ≤ `C4_COVARIANCE_BAND`.

**§10 ruling P8 (item 258's rule), applied mechanically.** Every rank enters `comm.allgather((size_local, num_ghosts))` on the tag-3 submesh's `index_map(3)`, at every width. Print item 258's (part) line with O_w and G_w.
- (U3c) is asserted iff G_w > 0 (as item 258's E3).
- The reading "(A) / (U4) on a submesh partition" is licensed iff O_w ≥ 2.
- The widths are fixed now, before any reading: 8 / 4 / 2. No width is added or chosen after a reading. "Inapplicable" and "not exercised" are records, not reds, in item 258's words.
- **Predicted, not measured:** the phantom is 537 of 80 798 cells (0.66 %). It may sit on one rank at every width, giving O_w = 1 and G_w = 0. In that case, record the outcome in item 258's wording.

**Negative controls (rule (e)).**
- **(U3c) ghost-included Ω** misses (A)'s volume identity by > `B_BAND`. *Asserted iff G_w > 0.* The factor is *predicted*: ≈ G_w/537 within ×[0.5, 2]. A single ghost cell is ≈ 1/537 ≈ 1.9e-3 of V_sub, ≥ 1e7× the band. Backing on another fixture: item 258's W2 measured 2.605e-2 against a sharp prediction of 2.341e-2 (`20261006T010208Z_MAT-19-step2-a2-W2.log`).
- **(C1″) B⁺ vs B > `WIDTH_BAND`**, *predicted*. This fixture has no prior measurement. Step 2 read 4.398e-3 on 2.5 mm cells. Here the cells are ≈ 15 mm (V_sub/537 ≈ 4.1e-7 m³), so the 0.1 mm shift moves the lattice set less. The prediction is [1e-4, 1e-2], with a ceiling far above the band. Assert only the miss.
- **(f″) C vs B > `F_CONTROL_FLOOR`**, *predicted*. The 64 MHz tuned field is not the 10 MHz one. Item 236's off-mirror (f′) read 35.78 % on this fixture with 10 g balls. The prediction is [0.1, 1.0]. Assert only the miss.
- **The width controls are not re-bought.** Item 258's (U3c) / (C1″) record stands.

**Printed (not asserted).**
- **(P4′) ball/point − 1 at t\* at A–E.** "Point" means the P1 θ at the centre. Compute it from the owned cell that contains the centre (barycentric λ ≥ −1e-12) with `allreduce((Σ value, count))`, and average over the containing cells. No `eval` and no `evaluate_vector_field_parallel` on the submesh: the owner rule raises there (known-issues 2026-10-04). On 15 mm cells a 5 mm ball is smaller than a cell, so the predicted magnitudes are ≤ 2 %. Also print each ball's V_h/(4πa³/3) − 1 (level-3 lattice, sub-tet edge ≈ 1.9 mm, predicted ≤ 3e-2).
- **The 64 / 10 MHz ratio** of the probe rises at t\* at A–E, C′ and Ω, printed twice:
  - **(i)** this run's rise per accepted watt ÷ step 2's rise at 10 MHz. Step 2 was normalised per *deposited* watt, G = 1 W / P₁, from the `MAT19-2   n=40` row of `20261006T005727Z_MAT-19-step2-a2-W1.log` (cite the line). Hard-code those values as printed records with their log line, never as comparands;
  - **(ii)** the same with the 64 MHz rise divided by η (per deposited watt at both frequencies, a shape-only ratio).
  - Label the mixed normalisation of (i) in the print. No 10 MHz per-accepted-watt figure exists.
- Also printed: the series table at n = 1, 10, 20, 30, 40; ⟨θ⟩(t\*) beside item 236's 0.687806 K per accepted watt (`R0`-style, ≤ `E_BAND` expected, printed); (S) via `_lattice_q_mean`; the q_h C / C′ pair; per-phase clocks; summed `ru_maxrss`.

**Windows.** Foreground. The harness call is the **whole** Bash command, bare and repo-relative. Bash `timeout: 660000`; `-s`; `tests/environment` first; no pipe or redirect. Each window is ≤ 590 s, so the durable-capture idiom is not needed.
- **Price.** Item 236's W1 measured 156 s at `-n 8` (fixture 90.2 s, EM drive 41.8–50.3 s, mesh 24.4 s). The probe adds the 40-step transient on 537 cells plus the probe sums, ≈ 5–20 s (item 258's series on 58 866 cells took 17 s). The `-n 4` / `-n 2` degree-2 EM on this mesh is **unmeasured**. It is predicted at ×1.6–2 and ×3–4 of the `-n 8` drive. Memory is predicted at ≈ 10–14 GiB (`ANS-4` item 225: 36.3 GiB on 259 524 cells at degree 2, scaled by cells).
- **First:** `docker compose exec -T fem-em-solver pgrep -c python3` must read 0. Then run the sweep:
  `scripts/testing/run_and_log.sh MAT-19-step3-cache-sweep "docker compose exec -T fem-em-solver bash -lc 'cd /root/.cache/fenics && for c in libffcx_*.c; do s=\${c%.c}; ls \$s.*.so >/dev/null 2>&1 || { echo orphan \$s; rm -f \$s.c \$s.c.cached \$s.o; }; done; echo sweep done'"`
- **W1, `-n 8`, writes the npz** (predicted 170–230 s):
  `scripts/testing/run_and_log.sh MAT-19-step3-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_MAT19_STEP3=1 timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment tests/validation/test_mat19_probe_tuned_copper.py -v -s --tb=short'"`
- **W2, `-n 4`** (predicted 220–330 s): the same command with `-n 4`, label `MAT-19-step3-W2`.
- **Box-speed gate, pre-registered.** If W2's EM drive clock exceeds 120 s, do not start W3. In that case (W) is asserted at `-n 4` only, and the record says "`-n 2` not run: W2 drive > 120 s". This is a time gate, not a reading.
- **W3, `-n 2`** (predicted 330–480 s): the same command with `-n 2`, label `MAT-19-step3-W3`.
- **Order:** pgrep → sweep → W1 → W2 → [gate] → W3. A red W1 stops the item. No edit after W1 (rule (i)); the npz hash enforces it.
- **The fallback is parking, never lengthening.** A kill at 590 is a finding. Never raise the timeout, add ranks or split a solve differently in-slot.

**Traps already paid for.**
- **Every rank enters every collective.** That covers the (part) allgather, each `_allsum` / `allreduce`, the point-value allreduce and the `Barrier` after the npz write, **including ranks that own no phantom cells**. A comm call under a `if comm.rank == 0` branch, or under `if owned.size`, deadlocks (item 227's rank-subset deadlock; at `-n 8` most ranks own no tag-3 cell). Only the `np.savez` itself sits under rank 0.
- **Rank-local values:** `cell_tags.values`, `assemble_scalar`, `size_local` and owned-slice max/min are rank-local. Reduce before asserting. `np.max` on an empty owned slice raises, so use `initial=`.
- **Orphan `.c`** (item 160's 44 KB partial): sweep only when `pgrep -c python3` reads 0, and delete orphan stems only.
- **Complex build:** `ufl.real` inside `ufl.lt`, and no `sqrt` in a comparison (`WF-6` 4c). Apply `np.real` before every `float`.
- **No point evaluation on the submesh** (known-issues 2026-10-04).
- **The (P) sheet terms use the Z_C-replaced specs** (item 236's trap).
- **Inherited from items 195 / 196 / 217 / 236:** κ is `None`; `z_from_s` is never applied; there is no `ans:7` `_stage_*` call; nothing is written under `logs/ans7_exact/`.
- **Never pick a width after a reading.** The stale `step2-n2.npz` is not this item's and is never read.
- **Bash guard:** harness calls are bare and whole; no `-k` selection inside the quoted command; status 124 is a finding, not a retry.
- **Write ≈, never `~`.** No AED file is opened or quoted.

**Scope.** Covered: one fixture (item 236's tuned copper F-small, 80 798 cells, 537 phantom cells, degree 2 EM, P1 thermal), one drive at 64.0 MHz at C_mode,Cu, the a = 5 mm probe at A–E, C′ and Ω, and widths 8 / 4 / 2.

Not claimed:
- any `ANS-9` agreement, adjudication or comparison (`ANS-9` is the 10 MHz ragged-sheet SPEC, and no AED value enters);
- that the coil is tuned to 64.000 MHz (T21 / T21′ are open);
- a real probe's sensing volume, thermal mass or lead artefact;
- probe-level accuracy on a 537-cell phantom (no h-statement);
- tissue, perfusion or compliance.

**Negative result.** In every case: report the measurement, park on `attempt/MAT-19-step3-*`, mark item 265 BLOCKED in the same commit (rule (d)), and journal it in attempts.md as attempt 1 of 4. Never widen a band, and never change a, the centres, the widths or the gate after a reading.
- **(A0) / (Src) red:** the source is not item 236's run. Stop.
- **(A) or (U4) red:** item 253's disposition. Open a known-issues entry; the review re-reads `MAT-19`'s ✅.
- **(W) red** with (A) green: the EM-in-the-loop series differs across widths above 1e-9 on degree 2. That is a finding for known-issues' solve-floor table, with the per-probe maxima printed.
- **(Σ) red:** on 15 mm cells, the probe pair inherits the unsymmetric phantom mesh. Item 236's 10 g q_h pair at φ₁ ± 90° read 22.7 %. Print the q_h C / C′ pair and both ball volumes. The review may re-pose on a phantom-refined exact hole mesh, a new fixture with its own C_mode. Not in-slot.
- **(U3c) asserted (G_w > 0) and not missing:** item 258's disposition.
- **A control not missing:** stop.
- **Kill at 590:** BLOCKED. Every finished window stands as a record.

No `PREDICTED RED:` line. A red anywhere stops the slot.

**Status it can move.** No glyph: `MAT-19` stays ✅. On green, the §7 row gains a **step 3 record** sentence: the probe on the tuned copper fixture at 64 MHz, degree 2, exact rectangles; (A) / (U4), with "on a submesh partition at `-n w`" only if O_w ≥ 2; (W) at the widths run; (Σ); (U3c) "missed by x" or "inapplicable (G_w = 0 at `-n 8 / 4 / 2`)"; ball/point at A–E and both 64 / 10 MHz ratios printed; elapsed times and logs. Rule (j): if the row lacks a step-3 Done-when clause, the reviewing daily adds "step 3 (§10 P10): identities green ⇒ a record, no glyph". §10 P10 is marked DONE.

**Spawn prompt (for the slot).** "Execute §9 item 265, `MAT-19` step 3, per `docs/planning/chunks/MAT-19.md` `### §9 item 265`. The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Every window runs foreground with Bash `timeout: 660000`. Every rank enters every collective. Widths 8 / 4 / 2 are fixed; apply the W2 time gate only. Never widen a band. A red stops the item."

**Review disposition of item 265 (2026-10-06 10:30 daily review): not re-queued; P10 goes to the 10-07 weekly.** This review has read the W1 log (`20261006T140518Z_MAT-19-step3-W1.log`) and lands it on `main` (by path, with its rows by hand) as a record.
- **The red is the brief's named (Σ) case, not a defect.** Every identity before it is green: (A) 4.549e-16, (U4) bitwise on a real submesh partition at `-n 8` (O_w 8, G_w 229), and (U3c) firing and missing by 0.547.
- **What goes red is the mesh's symmetry, not the probe's.** The σ_y pair of 5 mm balls sits in ≈ 15 mm cells on a 537-cell phantom that the mesher did not build mirror-symmetric. The source term itself differs 22.18 % between the two balls (q_h C / C′), against item 236's 22.7 % on its 10 g pair. On step 2's 58 866-cell phantom, the same probe read (Σ) 5.919e-4.
- **The brief's re-pose needs a phantom-refined *tuned copper* fixture with its own `C_mode`.** No such fixture exists. `FEM_EM_ANS7_EXACT_RESOLUTION` is a global resolution knob, not a phantom-only one, and every rung needs its own T-chain tuning window. That is a new fixture, which rule (6) does not let a daily open on a chain the weekly has not posed.
- **The alternative is a target question.** P10 could hold (Σ) on this 537-cell phantom as a printed reading instead of an assertion. That moves a target, so it is the weekly's.
- **For the weekly:** choose (a) a phantom-refined tuned fixture as a new P-step, serial on a T-step that tunes it, or (b) re-pose P10's mirror identity as printed on the record mesh, keeping (A), (U4) and (W) as its asserted identities.
- **Family:** `MAT-19` step 3 stands at attempt 1 of 4. The module stays on `attempt/MAT-19-step3-20261006T1405Z` (`a5b3c631`).

## §7 row moved verbatim 2026-10-05 (rotate_plan_archive.py chunks)

 **The virtual temperature probe: a point time series with a stated probe volume — Phase 8, gel** (opened 2026-10-03 by the weekly; §10 chain P7; lane 3) | ✅ Done 2026-10-05 (§9 item 243, 09:30 UTC slot): `thermal/probe.py` (owned-cell lattice ball average, one allreduce); (T) max 1.7315e-03 over 3 probes × 40 times (band 1 %); (W) `series_cf` 6.661e-16 at `-n 4` and `-n 1` vs `-n 2` (band 1e-9), cells equal; (C1) displaced probe misses at 5.4207e-02; (N1) 1.168e-14·θ_ss(0); (A0) 19 496 cells; elapsed W0 87 s / W1 104 s / W2 99 s / W3 117 s, standard tier (`20261005T095530Z_MAT-19-step1-W0.log`, `…095808Z_…W1.log`, `…100001Z_…W2.log`, `…100145Z_…W3.log`). **Step 2 record (2026-10-06, §9 item 258, attempt 2 of 4, §10 P8; attempt 1, item 253, was blocked because its (U3c) could not fire at `-n 2`):** the a = 5 mm probe on `MAT-15` P1's coil-driven gel phantom at `ANS-9`'s SPEC centres A–E plus C′. (A) whole-phantom probe = ∫θ/V_sub to ≤ 2.19e-14 (band `B_BAND` 1e-10), and (U4) is bitwise. Both hold **on a submesh partition at `-n 4`** (owned 23 401 / 35 465 / 0 / 0, ghosts 701 / 677 / 0 / 0). (W) with the EM in the loop is 3.714e-11 at `-n 4` and 5.735e-11 at `-n 1` (band 1e-9). (Σ) C / C′ at t\* is 5.919e-4 (band 5 %). (U3c) misses by 2.605e-2 at `-n 4`; it is inapplicable at `-n 2`, where G = 0. (C1″) reads 4.398e-3 and (f″) 0.8018. Printed only: (P4) ball/point − 1 at t\* is A +0.72 %, B −1.24 %, C +1.81 %, D −0.43 %, E +0.28 %, and (P5) the 0.1 mm displacement of B moves it 0.44 %. Elapsed: D4 68 s, W1 275 s, W2 248 s, W3 290 s, heavy by ceiling (`20261006T005546Z_MAT-19-step2-a2-D4.log`, `…005727Z_…a2-W1.log`, `…010208Z_…a2-W2.log`, `…010623Z_…a2-W3.log`). The fixture has ragged sheets at degree 1, so this is a record, no glyph. **Done-when (§4):** an additive probe post-processor in `thermal/` (ball-averaged temperature against time: an integral over owned cells with no point owner, width-invariant as `POST-8`. *Amended 2026-10-05 03:00 review, rule (j), item 243, because `MAT-12`'s thermal field lives on a `create_submesh` mesh, where the point-owner rule raises; the weekly ratifies*); on `MAT-12`'s SAR-heated sphere the series at three radii is within 1 % of the closed form at every recorded time, and the series is width-invariant to 1e-9 at `-n 1 / 2 / 4`; a displaced-probe control misses; elapsed recorded ⇒ ✅. Step 2 (§10 P8, provisional): its identities green ⇒ a record in this row, no glyph (added 2026-10-05 18:00 review, rule (j), item 258). Not claimed: a real probe's thermal mass or lead artefact. 

### §9 item 314 — 2026-10-08 (chain P10′)

**Outcome (2026-10-08 07:30 CDT slot, take-next): green — a record, no glyph.** Module `tests/validation/test_mat19_probe_refined_copper.py` (local `_geo_p` / `_z4_cu` copies, ruling (1); the local 4×4 ran the Cu column only — no PEC column, as item 310 on the h rung — because no anchor reads C_mode,PEC). W0: `_build_ladder` and the local default build both 80 798 cells (asserted); refined 93 956 cells, 8 349 tag-3 (predictions met; under the 130 000 stop). W1: Cu 4×4 sweep 65.7 s, recip 2.0e-15; C_mode,Cu 63.586252 pF; C_in 63.587076 pF (−4.74e-5 vs ×1, predicted ≈ 1e-3); \|g\| 3.5e-17 S in 2 iterations. W2 (`-n 8`): drives 49–60 s; (Φ) 0; (AB′) 1.257e-5; in-model zero f₀/64 MHz − 1 = −2.534e-6 = 0.20 × the AB′ resolution (printed); (P) 4.891e-13; η 0.5887; ⟨θ⟩(t\*) 0.680969 K per accepted W; (A) ≤ 3.95e-15; (U4) bitwise on a two-owner partition; (Σ) 1.2148e-3 (predicted [2e-3, 3e-2]: below range); (U3c) misses at factor 1.14; (C1″) 4.06e-3; (f″) 0.653. W3 (`-n 4`): (W) 1.956e-11, drive 87.4 s < 240 s, so W4 (`-n 2`) ran: (W) 1.694e-11. Printed predictions missed: q_h C/C′ ball pair 10.4 % (predicted ≤ 5 %) — the temperature pair agrees to 0.12 % because diffusion smooths it; ball/point at C +2.75 % (predicted ≤ 2 %). Logs: `20261008T124112Z_MAT-19-step3p-cache-sweep.log`, `20261008T124243Z_…-W0.log` (173 s), `…124546Z_…-W1.log` (166 s), `…124839Z_…-W2.log` (273 s), `…125322Z_…-W3.log` (194 s), `…125645Z_…-W4.log` (210 s). Attempt 1's module lands by path (its three logs and rows were already on main). Branch `attempt/MAT-19-step3-20261006T1405Z` deleted.

*Drafted by a general-purpose agent for the 2026-10-08 03:00 daily review; the review reviews it.* Licence: §10 P10′ (weekly 2026-10-07, PROJECT_PLAN.md:4648–4656) and its ruling (`:4499–4504`), quoted, not changed: "`MAT-19` step 3 re-pose — the probe series on a phantom-refined tuned copper fixture. `ANS-7`'s exact hole fixture with the phantom at 0.005 m, degree 2, retuned on C_in (three nodes, the in-model zero within its grid's resolution, printed); then P8's identities at `ANS-9`'s six centres over 40 steps: whole-phantom probe = ⟨θ⟩ at `B_BAND`, width invariance at 1e-9, the σ_y mirror pair at `C4_COVARIANCE_BAND`; printed, ball/point at A–E and the 64 MHz / 10 MHz probe-rise ratio per accepted watt. Target: identities. Heavy (predicted 300–450 s at `-n 8` from I14's 1 g-rung sweep). Serial on nothing." No band, tolerance, Done-when or target moves. **Attempt 2 of 4** on `MAT-19` step 3 (item 265 was attempt 1, parked `attempt/MAT-19-step3-20261006T1405Z`, `a5b3c631`).

- **Executor:** implementer. **Tier:** heavy, `-n 8` (record width) / `-n 4` (width), complex. **Time:** ≈ 30 measured / ≈ 55 predicted min (sweep 1.5 + W0 1.5 + W1 ≈ 6 + W2 ≈ 8 + W3 ≈ 6 + code ≈ 15 + spare). Tight for one hour: if W2 has not started by minute 35, park after W1 (its JSON is the save point) and journal.
- **Depends on nothing.** Reads no T21 JSON: this fixture is new and carries its own C_in.

**One window cannot hold it — split, never a longer timeout.** On the ×1 mesh the 4×4 at 64 MHz (Cu + PEC column, degree 2) cost 215 s (`20261004T140347Z_ANS-7-item217-W1-cmode.log`) and three C drives cost 261 s with a 171.4 s fixture (`…140738Z…-W2-inmodel.log:2052`); item 265's single drive + 40-step probe cost 168 s (`20261006T140518Z_MAT-19-step3-W1.log`, drive 41.5 s). The refined phantom adds ≈ 8 000 tag-3 cells (`GEO-27`'s 0.005 rung: 8 497 on the solid mesh, `20260908T033218Z_GEO-27.log.gz:3573`) to 80 798, so cells ≈ 90 000–115 000 and degree-2 cost ×1.2–1.8 (I14's 201 072-cell 1 g rung: sweep 246 s, 44 GiB at `-n 8`, `IMP-3.md:912–916`). Retune + 40 steps in one window ≈ 570–860 s > 590 s. So:
- **W0 (`-n 8`, build only, ≈ 90 s):** mesh, print cells and tag-3 cells. **Cost gate, pre-registered:** cells > 130 000 ⇒ stop before W1 (W2 projects past 590 s); record and park.
- **W1 (`-n 8`, ≈ 260–390 s):** the 4×4 at 64 MHz (Cu, PEC column) and C_in by `_c_in(z4, 64e6, c_mode)` (`test_ans7_exact_strip.py:1199`, secant to \|g\| ≤ 1e-15 S); writes `logs/mat19/step3p-cin.json` (C_in, C_mode, z4, cells, tag-3 cells; open mode `"x"`).
- **W2 (`-n 8`, ≈ 330–490 s):** reads the JSON; three Route A drives at C_in at 63.9 / 64.0 / 64.1 MHz (T21′'s construction, item 262), in-model zero and its grid resolution printed; at the 64.0 MHz node, the TH-14 (b) power block, G = 1 W / P_acc, P1 thermal 40 steps to t\*, the probe series at A–E, C′ and Ω; writes `logs/mat19/step3p-n8.npz`.
- **W3 (`-n 4`, ≈ 250–400 s):** one 64.0 MHz drive at the stored C_in + thermal + probe; (W) against the npz. **Time gate (item 265's):** if W3's drive clock > 240 s, `-n 2` is not run and the record says so; otherwise W4 at `-n 2` is the same command. A time gate, not a reading.

**Code.** New module `tests/validation/test_mat19_probe_refined_copper.py`, `@complex_only`, opt-in knob `FEM_EM_MAT19_STEP3P` ∈ {unset ⇒ skip, `build`, `cin`, `series`; anything else raises}. Template: the parked module `test_mat19_probe_tuned_copper.py` (`a5b3c631`, `run` at `:153`, `_point_values` at `:126`); take it **and its three logs by path** as the attempt-1 record (new files), append their `test-results.md` rows by hand, never path-checkout `test-results.md`. The mesh: a local `_exact_geo_p` = `_exact_geo`'s body (`test_ans7_exact_strip.py:409–421`) with `phantom_resolution=0.005` passed to `th14._sheets_build(True, as_hole=True, …)` (`tests/mesh/test_birdcage_port_sheets.py:80` accepts it). The 4×4: `th14._build_ladder` (`test_th14_birdcage_copper.py:263`) has **no** phantom knob — **review question 1**; the default plan is an additive `phantom_resolution=None` kwarg forwarded to `_sheets_build`, unset-identical by construction (item 214's `resolution` precedent), with rule (a)'s owed windows named by the review. Do not edit `test_ans7_exact_strip.py`, `test_mat15_tuned_copper_heating.py` (`ans:9` / `ans:7` import them) or the parked module.

**Anchors** (asserted; bands imported by name — P8's, unmoved):
- **(A0)** cells and tag-3 cells printed; tag-3 > 3 000 asserted as the precondition that the phantom *is* refined (item 265 had 537). **(Φ)** \|sin φ₁\| ≤ 1e-9 (item 236).
- **(Src)** (AB′) circuit = field at the 64.0 MHz C_in node ≤ `REDUCTION_BAND` (`test_ans7_exact_strip.py:81`; T21′ read 1.258e-5 on ×1); (P) the accepted-power identity ≤ `DISCRETE_IDENTITY_RTOL` (`test_birdcage_power_identity.py:227`).
- **(A)** whole-phantom probe = ∫θ/V_sub at every step ≤ `B_BAND` (`test_mat10_sar_to_heat_sphere.py:61`); (U4) bitwise.
- **(W)** series at `-n 4` (and `-n 2` if the gate passes) against the `-n 8` npz ≤ `WIDTH_BAND` (`test_mat19_virtual_temperature_probe.py:85`), npz hash checked as item 265.
- **(Σ)** \|C − C′\|/C at t\* ≤ `C4_COVARIANCE_BAND` (`test_birdcage_b1_plus_map.py:120`).

**Negative controls (rule (e)).** (U3c) ghost-included Ω misses (A) by > `B_BAND`, *asserted iff G_w > 0*, factor *predicted* G_w/N_sub within ×[0.5, 2] (item 265: 0.547 vs 0.426; item 258: 2.605e-2 vs 2.341e-2). (C1″) B⁺ vs B > `WIDTH_BAND`, *predicted* [1e-3, 1e-2] (item 258 4.398e-3 on 2.5 mm cells; item 265 1.95e-3 on 15 mm). (f″) C vs B > `F_CONTROL_FLOOR`, *predicted* [0.3, 0.8] (item 265 0.547). Each asserted only as a miss.

**Predictions (printed).** Tag-3 ≈ 6 000–10 000. In-model zero: f₀,A/64 MHz − 1 within the new grid's resolution (max AB′ node residual), order 1e-5 — printed, not asserted (the weekly's word). C_in within ≈ 1e-3 of the ×1 C_in 63.590091 pF. (Σ) predicted 2e-3 … 3e-2 (step 2's 5.919e-4 on 2.5 mm cells, item 265's 5.74e-2 on 15 mm); q_h C/C′ pair ≤ 5 %. η and ⟨θ⟩(t\*) per accepted watt beside item 236's 0.5777 / 0.687806 K (different mesh, reference only). (P4′) ball/point − 1 at A–E ≤ 2 %. The 64 / 10 MHz probe-rise ratio per accepted watt at A–E against `STEP2_N40_PER_DEP_W` (parked module `:108`), printed with its normalisation stated (item 265 printed two: mixed and per deposited watt).

**Commands** (bare, whole, foreground; Bash `timeout: 660000`; `-s`; `tests/environment` first, which runs `test_complex_mode.py`). First `docker compose exec -T fem-em-solver pgrep -c python3` must read 0; then `scripts/testing/run_and_log.sh MAT-19-step3p-cache-sweep "docker compose exec -T fem-em-solver bash -lc 'cd /root/.cache/fenics && for c in libffcx_*.c; do s=\${c%.c}; ls \$s.*.so >/dev/null 2>&1 || { echo orphan \$s; rm -f \$s.c \$s.c.cached \$s.o; }; done; echo sweep done'"`. Then, for STAGE/label/ranks = build/W0/8, cin/W1/8, series/W2/8, series/W3/4 (series/W4/2 if gated in):

`scripts/testing/run_and_log.sh MAT-19-step3p-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_MAT19_STEP3P=cin timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment tests/validation/test_mat19_probe_refined_copper.py -v -s --tb=short'"`

Every window ≤ 590 s, so §5.1's durable-capture idiom is not needed. No edit after W1 (rule (i)); the npz hash enforces it for W3/W4.

**Traps already paid for.** Item 265's list stands: every rank enters every collective, including ranks with no phantom cell (item 227's deadlock); only `np.savez` under rank 0. No point evaluation on the submesh (known-issues 2026-10-04): points by owned-cell barycentrics + allreduce. Never write under `logs/ans7_exact/` (T21′'s JSONs are another fixture's). Orphan `.c` of any size (item 160). The sheet terms of (P) use the Z_C-replaced specs (item 236). κ `None`, exact rectangles (ruling 2026-10-02). `ufl.real` inside `ufl.lt`, no `sqrt` in a comparison. A killed window is a record, never a retry with a longer timeout; check `pgrep -c python3` after any kill.

**Scope.** Covered: one phantom-refined exact hole mesh, degree-2 EM, copper by Leontovich, C_in at 64 MHz on that mesh, one drive, P1 thermal, the a = 5 mm probe at A–E, C′, Ω. **Not claimed:** a tuning claim (the in-model zero is printed), an h-statement (one phantom rung), absolute SAR or temperature, a real probe's thermal mass, `ANS-9`, tissue, compliance.

**Negative result.** Report, park on `attempt/MAT-19-step3p-*`, BLOCKED (rule (d)), attempt 2 of 4. W0 gate tripped: the cell count is the finding; the review re-prices (`xl` or a coarser phantom). (Src) red: the retune is not T21′'s construction on this mesh — stop before thermal. (A)/(U4) red: item 253's disposition. (W) red: known-issues' solve-floor table. **(Σ) red again on ≈ 5 mm cells:** print the q_h pair and ball volumes; the review decides between probe radius and cell size — no band moves. Kill at 590: every finished window stands. No `PREDICTED RED:` line.

**Status it can move.** None: `MAT-19` stays ✅ (Done-when, PROJECT_PLAN.md:2537, is step 1's). On green the §7 row gains a "step 3 record" sentence and §10 P10′ closes. **Rule (j):** the reviewing daily adds "step 3 (§10 P10′): identities green ⇒ a record, no glyph" (step 2's clause is the precedent). Item 265's parked branch is deleted on the landing.

**Review rulings (2026-10-08 03:00 daily review, item 314).** (1) **A local copy of the 4×4 build** in the new module, with the `phantom_resolution` argument; `test_th14_birdcage_copper.py` is not edited, so no rule (a) window is owed. W0 asserts the local build at the default resolution yields `_build_ladder`'s cell count (build only, no solve). (2) The W0 cap of 130 000 cells and the W3 240 s gate are **resource stops**, not bands: accepted as stops that park the item with the measurement, never as pass/fail of physics. (3) Accepted: T21′'s circuit = field at `REDUCTION_BAND` and (P) are asserted as T21′'s construction already asserts them; the in-model zero stays printed. (4) No `PREDICTED RED` line: (Σ) red parks the item as attempt 2 of 4 with the measurement, and the probe-radius-vs-cell-size call goes to the weekly. (5) The minute-35 park after W1 is accepted; W1's JSON is the save point and the next slot resumes from it. (7) No Done-when amendment: no glyph moves, so rule (j) does not apply; the outcome is recorded here and on the row's state line.
