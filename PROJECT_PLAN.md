# PROJECT_PLAN.md — FEM Electromagnetics Solver for MRI Coil Simulation

**Single source of truth for scope, status, and sequencing.** Resolved defects and
the reasoning behind past decisions live in [`docs/project-history.md`](docs/project-history.md);
nothing there is a task.

---

## 1. Mission *(rescoped 2026-08-04)*

A FEniCSX/DolfinX finite element toolkit that reproduces, for the slice of
electromagnetics relevant to **MRI RF safety**, the workflow an engineer runs
in Ansys Electronics Desktop — HFSS plus the circuit solver first, with the
Pennes bioheat equation for thermal simulation as the long-term extension.
The canonical workflow the tool must support end to end:

1. **Construct** a birdcage coil + gelled saline phantom simulation from
   parametric geometry — often with an implant inside the phantom.
2. **Tune** the birdcage at 64 MHz (1.5 T) and 128 MHz (3 T): EM solve plus
   circuit co-simulation to pick capacitor values, verify the mode spectrum,
   and match the ports.
3. **Drive** the tuned coil and extract the safety quantities: B1+ maps,
   SAR (whole-phantom and local, including near-implant hot spots), and port
   S-parameters.
4. Long term: **couple to thermal** via the Pennes bioheat equation.

Commercial solvers are expensive and black-box; open-source alternatives lack
the MRI-specific workflow. The scope is deliberately **the MRI-safety slice of
HFSS, not HFSS**: general-purpose 3-D full-wave parity is out of reach for
this project and is not the goal. Parity claims are made per-workflow ("tunes
a shielded 8-rung birdcage at 128 MHz to within X of AED"), never
per-product.

**Cross-validation against Ansys is part of the method.** `examples/` carries
runnable examples with XDMF outputs for visual review in ParaView, and
periodically a benchmark case specified precisely enough for the human
operator to replicate in Ansys Electronics Desktop; the returned numbers
become gates (§5.4).

---

## 2. Current state — read before planning work

What is validated, to what tolerance, and what must not be trusted.
(The history of how each claim was earned lives in §7's result blocks,
`docs/planning/plan-archive.md`, and `docs/project-history.md`.)

### 2.1 Validated, with the number that licenses it

- **Magnetostatics** (`core/solvers.py`, N1curl + gauge penalty, default
  1.0): Helmholtz **0.04%** centre / 0.83% mean; circular loop 7.07%;
  straight wire **`E_Ω` = 10.6172%** at h = 0.0025, falling
  25.2868 → 10.6172 → **6.6458%** on the h = 0.004/0.0025/0.0018 ladder at
  fitted rate **1.69** (`MAG-18`, 2026-08-22 — the annulus-restricted
  domain L2 that replaced the 10-point radial sample; the retired sampled
  statistic read 12.75% at n_points = 10 and 15.80% at 8 on the *same*
  field, so it never licensed a number); PEC cavity modes **0.0436%**
  (re-gated on the 0.11 image 2026-08-25 by `OPS-24` — `core/cavity.py` had
  been missed by the `OPS-18` migration and the `TH-9` gates were
  non-executing 2026-08-23 → 08-25; the migrated calls reproduce the 0.7.2
  per-mode errors 0.0123 / 0.0153 / 0.0201 / 0.0436% to the printed digit,
  known-issues entry retired).
  The `E_Ω` digits quoted here are the **0.11 image** ones, re-gated on
  `main` 2026-08-23 (§9 item 1, ruling (3\*)): the ladder is monotone at
  fitted rate **1.6854**, the record rung reproduces
  1.0617170177e-01 / 1.0617170374e-01 across two in-slot `-n 2` runs and
  1.0617175341e-01 at `-n 4` (4.86e-07 relative, inside the re-registered
  1e-6), and the natural-BC wall stays strictly worse (32.3155% vs
  10.6172%, ratio 0.3285). The 0.7.2 digits it replaces were
  25.3787 → 10.7288 → 6.6708% at rate 1.6842; no band moved. The
  `OPS-18` ✅ scope caveat is discharged.
- **Time-harmonic complex solve** (`TH-1`, complex build mandatory — real
  mode raises): MMS rate 0.9929; lossy plane wave α **0.019%** / β 0.059%
  (`TH-6`); evanescent waveguide γ **0.006%** (`TH-7`); quasi-static
  sphere 2.443% (`TH-8`); resonance guard calibrated (`TH-1` step 5).
- **The Larmor regime, on an imposed field** (`TH-10`): interior field vs
  the Mie series **3.643% / 1.826%** at 64 / 128 MHz (the 128 MHz figure
  re-recorded **1.769%** on the 0.11 image *with its mesh*,
  55 251 → 55 241 cells — `OPS-18` step 3 attempt 1, green log
  `20260822T123746Z`; 64 MHz reproduces bit-identically), and the
  SAR-relevant ½∫σ|E|² to **3.629%** — quasi-statics is the *wrong answer* there
  (under-predicts absorbed power 2.4× at 1.5 T). Degree-2 N1curl is gated
  on the same fixture (`TH-12` step 1, 2026-08-18): **0.1405%** interior
  relL2 on the coarse 5 866-cell rung — 25.9× the degree-1 fine-rung
  accuracy at 3.01× fewer cells, for 4.3× wall and 2.7× memory.
- **Coil loading, eddy-current regime only** (`MAT-6`): ΔR vs Dodd–Deeds
  **0.280%** against the filament closed form at 10 MHz, σ = 100 S/m, on
  the production projected drive and the slab-refined **418 888-cell**
  fixture promoted by step 11 (2026-09-06; 1.58% was the pre-refinement
  138 619-cell record); ΔX is reported, never gated (not box-converged).
  The closed form is a filament; the fixture is a 5 mm wire, and that
  modelling term is quantified (`MAT-8` ✅ 2026-09-02): **+0.115%** on ΔR,
  **+0.145%** on ΔX at the fixture. It raises the closed form, and the
  refined FEM sits *below* the filament value, so the discrepancy against
  the corrected form is **≈ −0.40%**, and 0.115% is the floor under any
  sub-0.5% ΔR claim on this fixture.
- **SAR machinery on an imposed uniform field** (`MAT-4`): mean SAR vs the
  lossy-sphere closed form to 3.5%; the 1 g/10 g averaging operator exact.
- **Mass-averaged 10 g SAR on the coil-driven field** (`MAT-4` step 4 ✅
  2026-09-06): the same C95.3 operator applied to the four single-drive
  F-small birdcage solves on the finer-phantom rung passes a **C4
  covariance identity** at 0.3303 / 0.0756 / 0.0574 / 0.3132% against the
  unmoved 5% band, with the whole-phantom coverage identity exact to
  1.58e-14 and the far-side ball control at ≈86%. **A symmetry identity on
  one fixture at fixed `h` at 10 MHz — not an absolute SAR number, not a
  C95.3 compliance or limit claim, not Larmor, not converged in `h`; the
  1 g column misses the same band at 6.84% on that rung and is a printed
  record — but on the 2.5 mm phantom rung (`MAT-4` step 5b, 2026-09-08,
  199 920 cells) **both columns are gated at that same unmoved 5% band**:
  the 1 g pairs read 0.0957 / 0.1199 / 0.1305 / 0.1065%, asserted, with
  the mis-paired 1 g control at 87.01–87.06% and the rung's mesh pinned as
  a version-tagged record at the 1% band. Still one fixture, one
  frequency, two rungs — not a convergence rate.**
- **S-parameters, two-torus validation fixture only** (`PORT-1` ✅
  2026-08-15): field-derived through `run_n_port_sparameter_sweep`,
  reciprocity `‖S−Sᵀ‖/‖S‖ = 2.5494e-05` vs the 1e-3 gate — carrying the
  **two named systematics** of 3b-xviii (PEC box, an effective-range
  extrapolation; gap-generator feed model, Jin 3e §10.4.2.1). The retired
  heuristic is reachable only behind a `DeprecationWarning`.
  *Caveat lifted 2026-08-24 (`PORT-9` leg (d3), §9 item 2):* the sweep's
  S no longer comes via the terminated `Z` — the gated routes assemble it
  from power waves (`S_ij = b_i/a_j`), and on the asymmetric two-torus at
  the matched drive `Z_p = z0` the per-pair asymmetry fell from
  2.831857978e-03 (old conversion, same run) to **2.972992845e-15**. The
  `Z` this bullet's records once went through is retained as a documented
  *terminated transimpedance* diagnostic and is never reciprocity-gated.
  The two moved records were re-recorded route-tagged under (1\*) with no
  band moved (σ_max 0.864809457, `‖S−Sᵀ‖/‖S‖` 4.758625e-05); the
  gap-voltage route's undriven ports are not terminated in `z0`, so its
  residual is the route's own and stays gated at 1e-3. The two named
  systematics and the two-torus-only scope are unchanged.
  *Corrected 2026-09-20 (`PORT-20`, from the `ANS-3` adjudication):* **S on
  the current route is `z_to_s(Z)` since `PORT-20`** — on the gap-voltage /
  impressed-current route the undriven ports are *open*, so `b_i/a_j` was
  never a column of S there; that branch now converts its open-circuit `Z`
  (agreeing with `ports/circuit.py`'s independent implementation at
  9.4e-17 relative, and with the closed form `(Z−z0)(Z+z0)⁻¹` exactly on
  synthetic 2-/3-ports). Only this bullet's two S records moved —
  σ_max 0.864809457 → **0.861356894**, `‖S−Sᵀ‖/‖S‖` 4.758625e-05 →
  **3.112129e-05** — with every band unmoved and every `Z`-derived record
  (the mutual ratio above included) digit-identical. The **lumped-sheet**
  routes of the two bullets below are untouched and re-ran byte-identical.
- **Four lumped-sheet ports on the gapped 4-leg birdcage, at 10 MHz
  only** (`PORT-9` ✅ 2026-08-25, step 3 legs (c)/(d0)/(d)/(d1′)). Four
  driven solves on `GEO-18`'s gapped, sheeted fixture (116 085 cells,
  `f = 0.5`, `w = A/h`, `Z_p = z0 = 50 Ω`) assemble a 4×4 that passes
  three pre-stated gates: **(i)** `‖S−Sᵀ‖/‖S‖` ≈ 1e-14 vs 1e-3 on the
  power-wave assembly; **(ii)** `σ_max(S)` = 0.999992805 ≤ 1 + 1e-9 with
  column power sums ≤ 0.793823974; **(iii′)** every C4 circulant class of
  `Z` spreads ≤ **0.5%** (measured 0.0553 / 0.0353 / 0.0214%), the band
  tightened from 5% with this closure. The gate has its **geometric
  negative control**: rotating one leg 22.5° breaks all three classes by
  two orders (6.2219 / 7.1142 / 2.8474%) while reciprocity holds at
  2.259e-14, and the pre-fix terminated-`Z` route read 5.57e-03 on that
  same displaced fixture — a 2.466e+11× separation. **Scope, hard:**
  10 MHz is the port model's frequency, not a Larmor frequency; this
  licenses no resonance, tuning, B1+ or SAR claim. See §2.2.
- **The same three gates pass at both Larmor frequencies on the same
  fixture** (`PORT-11` ✅ 2026-08-26, steps 2 and 3; audited COMPLIANT
  18:00 review). Frequency is demonstrably the only knob turned — the
  in-run 10 MHz rung reproduces leg (d)'s 4×4 to 1.158e-10 vs 1e-6 in
  both modules. **64 MHz:** reciprocity 2.581325834e-14 vs 1e-3, σ_max
  0.999721388 ≤ 1 + 1e-9, C4 class spreads 0.0573 / 0.0599 / 0.0370% vs
  0.5%. **128 MHz:** 7.030990825e-15, σ_max 0.998974779, spreads
  0.1012 / 0.0916 / 0.0654% — on a phantom that has crossed to
  displacement-dominated (loss tangent 0.9002, cells/λ 12.5024 against a
  pre-stated floor of 10, enforced before any gate is read). The 22.5°
  displaced control breaks (iii′) at 12.9 / 27.8% (64) and 16.7 / 34.7%
  (128) while (i) holds at ≈1e-15. **Scope, hard:** this is a
  self-consistency identity set (reciprocity, passivity, C4 symmetry)
  on one gated fixture — it is *not* an absolute-accuracy, resonance,
  tuning, B1+ or SAR claim, and the feed systematics remain the
  two-torus ones (`PORT-1` 3b-xviii, `PORT-10`). The absolute comparison
  is the AED benchmark's job (`ANS-4` — runnable half ✅ 2026-08-30; the
  operator's AED replication at Zero and First Order **landed 2026-09-04**,
  numbers private — **adjudicated 2026-09-06 (weekly review, run
  interactively): AGREE at 10 MHz** (couplings and the terminated-drive
  input impedance agree to the few-percent level a converged AED pair
  brackets, so the port model is not wrong by a constant factor), and
  **INCONCLUSIVE at 64 / 128 MHz** on 09-06 — **re-adjudicated 2026-09-13
  (weekly review): AGREE at 128 MHz** on `ANS-4` step 2d's order-matched
  rung (our degree 2 is HFSS First Order by `ANS-5`; the degree-2 sequence
  converges, successive change 1.19 → 0.62 % on `S₁₁`, and the residual
  against AED's converged column is in the class of AED's own two orders'
  mutual agreement) and **AGREE by mechanism at 64 MHz** (the rung ran at
  128 MHz only; `ANS-4` step 3, `xl`, is the 64 MHz evidence — **its window
  ran 2026-09-16 02:00**: 16 passed, 1658 s, `memory.peak` 280.8 GiB, every
  imported gate green on the 592 744-cell degree-2 rung, and the degree
  1 → 2 move at 64 MHz on the driven column is 6.51 / 2.60 / 4.11 % by C4
  class (`20260916T070008Z_ANS-4-step3.log:2427, 4596`); the private miss
  against AED and the pre-registered decision rule are the 09-19 weekly's,
  so this sentence does not move until then). The 09-06
  "not element order" inference was wrong — AED was order-insensitive
  because it had converged, ours moved 6.09 / 5.38 / 6.70 % because it had
  not. **What this means for every degree-1 S figure at 128 MHz: the gate
  fixture's own entries sit 5–7 % from their order-matched value and every
  identity gate passes on them unchanged** — identities do not see this
  class of error. No absolute figure is quoted, nothing at Larmor beyond
  the verdict word is promoted, numbers private, §10 2026-09-13).
- **A package-level multi-port drive exists, gated on the 4-leg
  identities** (`POST-6` step 1, 2026-09-04; feature ladder A1):
  `ports.superpose_drives` combines N stored single-drive solves under a
  complex weight vector on the drives' own N1curl space, and through it the
  quadrature weights reproduce `WF-6` step 2's C4 / mirror records
  (0.9818% / 0.8087%, package path vs fixture path to 1.2e-15), linearity
  holds at 1e-12 and `w = e_k` returns drive k bit for bit
  (`20260905T004202Z_POST-6.log:1911–1914`). **The drive-level power
  statement is now gated in its exact form** (`PORT-16` step 2,
  2026-09-07): on the superposed quadrature field,
  `P_src,exact(w) = P_vol(w) + P_sheet,exact(w)` closes at **5.877e-15 /
  1.435e-14** against 1e-6 and the attribution
  `P_acc(w) − P_vol(w) = C(w) − sheets_terminal(w)` at **9.881e-14 /
  1.581e-13** against 1e-1
  (`20260907T110826Z_PORT-16.log:1942–1946`), with the `w = e_k`
  generalisation reproducing step 1's single-drive shares bit for bit
  (`:1948–1951`). **What it is not:** the *terminal* form
  `½aᴴ(I − SᴴS)a = ½∫σ|E_w|²` still reads **11.648%** — that figure is a
  **record, not a conservation band** (03:00 review ruling, 2026-09-07):
  the exact identity closes at 1e-14 on every drive and the ≈6.7e-05 W the
  terminal form misses is its own Cauchy–Schwarz deficit (reproduced to
  1e-13), a fixed 1.06% sheet-field non-uniformity at this `h` over a
  drive-dependent denominator, not a leak. `POWER_BALANCE_BAND` keeps its
  1e-2 value and its single-drive asserts; the deliberate red retired with
  step 2's commit.
  Beside it, **the lumped RLC sheet and the circuit layer are gated on the
  same fixture** (`PORT-14` ✅ 2026-09-13: the termination-reduction
  identity holds under 1e-3 at 10 MHz on the registered floor and at 64 MHz
  on the κ-derived width `w/(1 + κ)`, κ = 1.06 % the sheet's named
  systematic — 5.359129e-05 / 1.998404e-06 for C / L,
  `20260914T003252Z_PORT-14-step3-64mhz.log:1937–1959`; `PORT-15` ✅
  2026-09-14: `C_tuned` = 15.570 pF from the stored 4×4 zeroes `Im Z_in`
  to 2.1e-15 and one in-model solve at that C reproduces the circuit's
  tuned `S₁₁` to 8.26e-5 and the reduced 2×2 to 6.12e-5,
  `20260914T020419Z_PORT-15.log:1967–1968, 3732–3742`) — **self-consistency
  between the field solve and the circuit reduction on one F-small fixture
  at 64 MHz; "tuned" is a series resonance (`R_in` 6.77 Ω, `|S₁₁|` 0.761),
  not a match, 128 MHz is a printed reading, and no mode frequency is
  claimed** (`TH-17`).
  **The F-small tuned resonance is located by a driven sweep consistent
  with the circuit layer across 44–84 MHz** (`PORT-22` step 1 ✅
  2026-09-20): on the 116 085-cell gate mesh at degree 1, `Im Z_in(f)`
  from the in-model tuned drive runs −34.62 → +25.60 Ω over the grid with
  **exactly one sign change, in the [60, 64] MHz bracket** (interpolated
  zero 63.9986 MHz, `R_in` 6.7725 Ω, loaded Q 6.93 — printed, not gated),
  and at **every** one of the 11 grid frequencies the circuit reduction of
  that frequency's *untuned* 4×4 reproduces the in-model tuned `S₁₁`
  within the imported `REDUCTION_BAND` = 1e-3 (worst 8.256069e-05 at
  64 MHz, 0.083× band; κ computed in-run, 1.0600e-02 → 1.0618e-02 across
  the band), while the same prediction at `0.5 × C_tuned` misses by
  1.116902 = **1117× the band** (`20260920T110503Z_PORT-22.log:2136–2153`,
  `20260920T110833Z_PORT-22.log:2082–2108`). **Circuit-vs-field identity
  on one fixture, degree 1, series resonance — no absolute `S₁₁` or
  `Z_in` figure (known-issues 2026-09-19) and no mode-spectrum claim
  (`TH-17`'s).**
  **The circuit layer also runs through a vetted engine, `scikit-rf`
  2.1.0** (`PORT-23` step 1, 2026-09-23; `OPS-60` ✅ pins it): an
  `skrf.circuit.Circuit` netlist built in our code terminates the stored
  4×4 and reproduces the raw (C2) reduction to ≤ 3.5e-16 at 64 and 10 MHz,
  its sweep zero is `C_tuned` exactly, and a mis-ordered port map misses by
  201.9× the floor (`20260923T093456Z_PORT-23.log`). The raw layer stays
  as the gate; this is an engine identity, not a new physics claim.
- **A PEC interior body is gated on the sphere's exterior dipole
  coefficient at 10 MHz** (`TH-15` step 1, 2026-09-05): the conductor
  solved as a *hole* with `n × E = 0` on its cavity facets fits
  β = 1.019746 against the perfect-conductor β = 1, **1.97% vs the 4.886%
  band**, with the field's relative-L2 exterior miss converging at
  **+1.19 in h** over the `TH-8` rungs (`20260905T170225Z_TH-15.log:290–296`).
  The natural-cavity control misses by 29.0× the band. **Not** loss, Q,
  copper (`TH-14`), the two-torus (step 2) or a coil (step 3); the
  pointwise field miss on the gate rung is 46.07%, a first-order N1curl
  floor recorded and asserted nowhere. The two-torus and the 4-leg
  birdcage both **mesh** as holes through `MeshGenerator` (steps 2a / 3a,
  2026-09-06: 161 461 and 80 181 cells, every sheet and cavity-area
  identity asserted, `20260906T050829Z_TH-15.log:1463–1477`,
  `20260906T213913Z_TH-15.log:2628–2650`). **Since `TH-15` step 3
  (2026-09-13) the birdcage hole also *solves* with its four lumped-sheet
  ports:** the 4×4 passes `PORT-9`/`PORT-11`'s imported identity gates at
  10 / 64 / 128 MHz (reciprocity ≤ 1.6e-14, σ_max ≤ 0.999994, class spread
  ≤ 0.0734 % vs 0.5 %) with the solid control reproducing `LEG_D_S_MATRIX_10MHZ`
  in the same window (`20260913T202954Z_TH-15.log:951–954, 2768`; 64 / 128
  MHz `…202718Z_TH-15.log:951–954`, `…202828Z_TH-15.log:956–959`) —
  identities on the hole route, **not** an accuracy claim: the hole-vs-solid
  `max\|ΔS\|` is a printed record, and the terminal-power sum's excess over
  the volume loss (7.6937e-05 W, 1 207× the phantom loss) **is the sheets'
  own Cauchy–Schwarz deficit `C − terminal`, equal to every printed digit on
  the hole and on the solid** (step 3c,
  `20260914T005452Z_TH-15-step3c.log:969–970, 2790–2791`; hole/solid
  `C/terminal − 1` 0.9997) — the registered power sentence is
  `Σ½Re(V I*) = P_vol + (C − terminal)` (re-registered 2026-09-14 03:00
  review under rule (h); the 2026-09-13 known-issues entry is retired). The
  two-torus hole carries no gated port current (step 2, killed by the
  2026-10-07 weekly). `TH-15` ✅ (2026-10-19 weekly): this paragraph is its
  conductor-model line — PEC interior bodies gated as above; loss, Q and
  copper are `TH-14`'s, below.
  **Copper as a Leontovich surface (`TH-14`, ✅ 2026-09-14):** the third-kind term
  `jωμ₀/Z_s ∫_Γc (n × E)·(n × W) dS` is gated on the lossy-wall `TH-9`
  cavity's Q against Pozar's closed form (+0.010 % at σ = 1e4, `Q ∝ √σ` to
  0.05 %, the PEC pencil damping-free to 4.8e-19,
  `20260914T004807Z_TH-14.log`) and on the F-small birdcage as a copper
  hole — the three port identities at 10 / 64 / 128 MHz, the surface-loss
  identity to 3e-13, a σ-ladder 5.8e7 → 5.8e11 whose distance from the PEC
  4×4 falls 10× per 100× σ, coil-loss share 0.929 / 0.448 / 0.219 printed
  (`20260914T021500Z_TH-14.log:1037–1364`), and — the closed form at
  copper σ — on the **Dodd–Deeds copper slab as a Leontovich floor** under
  `MAT-6`'s loop at 10 MHz: ΔR(5.8e7) −0.299 % from `coil_impedance_change`
  (band 2 %), the surface-loss identity to 1.31e-8, ΔR(5.8e7)/ΔR(5.8e9) =
  9.9910 (linear in `R_s`, band 1 %), the thin-skin identity (ΔX − ΔX_PEC)/ΔR
  = 1.00097 vs the closed form's 1.00095 (band 2 %), PEC-floor control 0
  (`20260914T094224Z_TH-14-step2.log:283–293`, 340 s; the slab is removed
  by `create_submesh` so the floor is a true domain boundary). One loop,
  one liftoff, one frequency; the absolute S on a copper coil is `ANS-6`'s,
  below (AGREE within 0.33 % of AED First Order, 2026-10-12 weekly).
  **On a resonance, externally (`ANS-7`, adjudicated by the 2026-09-30
  weekly):** on the tuned copper F-small birdcage both codes find the one
  k = ±1 parallel resonance, and the conductor's part of 1/Q, `Q_Cu/Q_PEC`
  and the Cu − PEC shift of `f₀` agree with AED to within 5 % — the second
  external check of `TH-14`, the first on a resonance. Degree 1, one rung;
  the absolute `f₀` / Q / `R_in` are §2.2's (numbers private). **The tuned
  `f₀`, with its h-statement (T13, decided by the 2026-10-10 weekly):** on
  exact rectangles at degree 2, `f₀` sits within 0.11 % of AED First Order at
  ×1 and **within 0.03 % at h = 0.005, converging first order in h**; Q and
  `R_in` converge in h too (T41, 2026-10-17: successive-change ratios 0.717 /
  0.705 on the steps 0.015 / 0.01 / 0.005 m, where first order gives 1.00 and
  second 0.60; p ≈ 1.65) to a level ≈ 3 % below AED First Order — an
  **INCONCLUSIVE** absolute Q, not an agreement, until AED's own pass history
  is read (T48); no Q, `Z_in` or match figure is licensed.
- **The F-human birdcage rung is a gated *mesh*** (`GEO-25` step 2,
  2026-09-05): the 16-leg, 32-ring-port coil at `ring_radius` 0.15 m
  meshes to **504 642 cells** at fixed absolute sizing with the volume
  partition exact, all 32 terminal ratios inside [0.95, 1.0] and
  meshed/CAD conductor mass 0.965414 ≥ 0.95, while the same coil with the
  sizing scaled to the radius reads 0.893028 and fails the gate
  (`20260905T183654Z_GEO-25.log:9728–9735, 20017–20027`). CAD identities
  only — and since **`WF-7` step 0 (2026-09-13) one degree-1 single-drive
  lumped-sheet solve has run on the longitudinal-sheet variant of this rung
  at 64 MHz** (507 266 cells, 607 039 unknowns, 37 s solve, 10.93 GiB summed
  `ru_maxrss` at `-n 8`, `20260913T190102Z_WF-7-step0.log:10419–10427`) —
  a price, not a gate: nothing about that field is asserted, and no Phase 6
  date rests on it.
- **Human scale, symmetry-gated (`WF-7` ✅ 2026-09-25, audited PASS
  2026-09-26; written by the 2026-09-26 weekly).** On the F-human 32-port
  birdcage at 64 MHz, degree 1: the 32×32 passes `PORT-13`'s imported
  identities, the ccw quadrature `|B₁⁺|` is C16-invariant to 0.15 %, and the
  10 g SAR hotspot is located (r ≈ 45 mm) with its 15 rotated images agreeing
  to 0.49 % and its mirror to 0.05 % against the imported 5 % band
  (`20260925T093315Z_WF-7-step3.log:11615–11617`; `EX-62`). **What this is
  not:** the coil is **untuned** (32 ports in 50 Ω), one phantom, no external
  check, no compliance or absolute figure; its S carries a 5.5 % degree 1 → 2
  order move, and the degree-2 terminal accounting residual (≈ 1 % on every
  column) is unattributed until `POST-7`. **The tuned coil (`WF-10` ✅
  2026-10-11, audited PASS by the 03:00 daily; line written by the 2026-10-10
  weekly; degree-2 clause §10 H31, 2026-10-22):** the tuned copper F-human (in-model root 118.2437 pF, 0.2 Ω ESR,
  two-port 90° drive) has `|B₁⁺|` and 10 g SAR readings gated on the
  feed-bisector mirror at degree 1 (0.1965 % and 0.6413 % against the 5 %
  band; `20261010T170356Z_WF-10-step2-W1.log:6335`,
  `20261010T171705Z_WF-10-step2s-W1.log:6373–6374`) — a record. **At the
  degree-2 root** (`xl` entry 41, 2026-10-22; degree 2 at C_in(d2) =
  115.5209 pF, every gate green — (M_B) 1.7093e-03 and the 10 g mirror
  (iii-a) 0.3052 % against the 5 % band): mean `|B₁⁺|` 2.792653e-07 T and
  SAR10g(c\*) 4.547175546e-03 W/kg (5.402616e-02 W/kg per W of P_src) under
  the two-port drive; the root-to-root order moves −3.17 % / −6.30 % (per
  phantom watt +0.60 / +1.13 %) fire the pre-registered rule's second branch,
  so the human-scale tuned figures are quoted **at degree 2 only, at
  C_in(d2)** (`20261022T070008Z_WF-10-step4-H28.log:6847`, `:6852`,
  `:13586–13587`) — read literally by the 10-22 03:00 review and **ruled
  by the 10-24 weekly: Phase 6 exited on this reading** (§10 rulings;
  §6). No absolute-external, compliance or tissue claim.
- **F-clinical, degree-1 records and the degree-2 price** (`GEO-38` chain C,
  flagged by the 10-31 weekly; internal statements, not gates on a tuned
  coil). 1.5 T: the untuned 64×64 is reciprocal and passive, and the in-model
  Larmor root is C_in = 193.5996 pF at 0.25 pF resolution (`xxl` entry 44,
  `xl` entry 51).
  3 T: the untuned 64×64 at 128 MHz is reciprocal, passive and
  copper-identity green, with one class asymmetric beyond the small-entry
  reading, physics or mesh open (`xl` entry 47; known-issues 2026-10-27); the
  in-model degree-1 root is **C_in = 44.0615 pF at the grid's 0.28 pF
  resolution**, one crossing on all four series (`xl` entry 49, read first
  branch by the 11-02 weekly; a record — the sign of root − `C_mode` differs
  from the 1.5 T twin's and is C26's to read), **refined to 44.0525 pF at
  0.05 pF** (`xl` entry 57, 2026-11-08, every gate green, one crossing on
  P33 (1, +1j) at 44.052487 pF by the bracketing pair; read first branch by
  the 03:00 review, the 11-09 weekly reads it — the node of entries 58 / 61). **The rod in the 1.5 T ASTM box
  under the two-port drive is a record** (`xl` entry 48, read by the 11-02
  weekly): every gate that ran green, tip 1 g SAR ratio rod / gel twin
  4.34–4.46, η_rod(T1) 3.642, the rod's share of the volume power 3.7e-5, the
  port quantities perturbed at ≈ 1e-7; box containment not asserted on this
  route. **The box-loaded 1.5 T coil at its degree-1 root is the first
  clinical-scale `|B₁⁺|` / 10 g record** (`xl` entry 52, C6, read first branch
  by the 11-04 weekly, every assertion green): mean `|B₁⁺|`_ccw 2.090e-08 T
  over 400 tag-3 centroids and ≈ 1.498e-08 T over the registered 18-point
  lattice (each the record under its own definition); SAR10g(c\*)
  1.379e-04 W/kg, 2.274e-03 W/kg per accepted watt, with c\* the argmax of the
  18 lattice balls, not a box-wide maximum (C6b reads that). **Degree 2:** the
  1.5 T root is bracketed at 185.763346 pF by the chord at 1.25 pF node
  spacing (`xl` entry 53, read first branch, a record; −4.05 % from the
  degree-1 root), **read at 0.25 pF as C_in(d2) = 185.7629 pF** (`xl` entry
  56, read first branch by the 11-07 weekly; the node of `xxl` entry 62 and
  `ANS-11` SPEC v2). **The degree-2 1.5 T 64×64 is a record** (`xxl` entry 54,
  read second branch by the 11-07 weekly: reciprocal, passive,
  copper-consistent; one class above the spread band at 0.6424 % against
  0.9262 % at degree 1, so the asymmetry is not degree-1 discretisation;
  `C_mode`(d2) 185.4048 pF, C18). **The rod in the box heated by the coil**
  (`xl` entry 55, I32): the EM record stands; the tip calorimetric gate was
  red on its comparand's form and is re-posed on the consistent-mass
  projection (11-07 weekly; I32b, then entry 63). **The box heated by the
  coil is a record at both field strengths** (`xl` entries 59 / 60, P33 /
  P33-3T, read first branch by the 11-11 weekly — every executed assertion
  green): θ(t\*) at c\* per accepted watt 1.611905e-03 K at 1.5 T and
  6.516151e-03 K at 3 T (ccw); the 3 T / 1.5 T lever 4.0425 ccw / 4.5397 cw
  (C39) beside the SAR10g lever 4.1014 and the closed-form pictures' 3.6935;
  `cem43` per accepted watt is the 20 °C baseline, not a dose — degree 1,
  20 °C gel, no container wall, the 3 T gel's constants provisional;
  `ANS-13` (commissioned 11-11) is the external reading. The degree-2 column is priced at 472.42 GiB and ≈ 51 min
  per node pair (`xl` entry 50), so root grids run on `xl` at ≤ 3 nodes and
  column sets on `xxl`. **Not claimed:** an absolute or external tuning
  figure, an absolute |B₁⁺| or SAR figure, a homogeneity or compliance
  figure, any external check — every clinical figure above is an internal
  statement (the roots internal tuning statements, the matrices and maps
  records; `ANS-11`, commissioned 11-04, is the external reading, on the
  degree-2 column of `xxl` entry 62).

- **Heat conduction in the gel regime, on a prescribed source** (`MAT-9`
  ✅ 2026-09-29, audited PASS by the 18:00 daily; the lane-3 opening awaits
  the weekly's ratification). In a Robin-cooled homogeneous ball, the
  steady solve meets the closed form to 0.083 % (P1 rate 2.04,
  `20260929T034247Z_MAT-9-step1.log:88–91`), and the transient solve meets
  the eigen-series to 0.046 % (CN rate 2.001), with the Bi = 5/3 roots
  accurate to 2.2e-16 relative over n ≤ 4000
  (`20260929T112146Z_MAT-9-step2.log:113–129`,
  `20260929T170037Z_MAT-9-step2b.log:105`). **What this is not:** no
  perfusion (Pennes), no SAR-driven source, no temperature-dependent
  property, one geometry.
- **The FEM SAR as a heat source, on the uniform-field sphere** (`MAT-10`
  ✅ 2026-09-30, audited PASS by the 03:00 daily; lane-3 ratification the
  weekly's). `MAT-4`'s SAR, cell-averaged onto the sphere's own submesh,
  heats the Robin ball to 3.25 % of the closed form against a composed 11 %
  band (`MAT-4`'s 10 % + `MAT-9`'s 1 %); the tight assertions are the
  uniform-source premise (0.319 % against 1 %), the power transfer
  (2.8e-15) and the surface balance (1.2e-13)
  (`20260930T021802Z_MAT-10-step1-W1.log:160–164`). **What this is not:**
  the absolute SAR accuracy is still `MAT-4`'s 10 %; an imposed uniform
  field, never a coil; this entry is the steady state without perfusion or
  σ(T) — those are `MAT-11`, `MAT-12` and `MAT-5`, below.
- **Layered impedance boundary** (`TH-18` ✅ 2026-09-30, audited PASS by
  the 03:00 daily; **its gate's σ = 10 S/m amendment awaits the weekly's
  ratification**). A volumetric εᵣ = 4 slab terminated by a conductor layer
  carried as a surface term reproduces the stack's closed-form reflection
  to 1.2e-6 in \|Γ\| and 1.4e-5 rad at 127.74 MHz, normal incidence,
  degree 2, inside a band (3e-6, 5e-5 rad) frozen from a bulk-Leontovich
  run of the same line; bulk-vs-thin control 0.093
  (`20260930T021508Z_TH-18-step2-W2.log:100–116`). **What this is not:** a
  copper gate (copper's \|Γ\| is within 3e-5 of 1 at any thickness here),
  oblique incidence, or any birdcage use.
- **First-order absorbing boundary** (`TH-5` ✅ 2026-10-08, step 2c, §9
  item 320; the weekly 2026-10-07's line, verbatim): the first-order ABC is
  valid for k₀R ≳ π (`TH-5` 1a / 1b); on an electrically small box it is a
  reflector whose box term is opposite in sign to PEC's and 0.23–0.48 of it,
  and a ⅓ PEC + ⅔ natural mixed wall cancels it to 0.25 pp (two-torus,
  10 MHz). (`20261008T183105Z_TH-5-step2c-W1.log:938–940`,
  `20261008T183455Z_TH-5-step2c-W2.log:606`.) **What this is not:** the
  band k₀R ∈ (0.04, π), which holds the 64 / 128 MHz birdcage box; a
  second-order ABC or PML; any open-boundary birdcage, SAR or HFSS
  *Radiation* parity figure.
- **Gel presets and the two boundary options** (`MAT-1` ✅, `TH-3` ✅
  2026-09-30, both audited PASS by the 03:00 daily). The low / mid / high
  σ gelled-saline presets read SAR within 3.52–3.56 % of the lossy-sphere
  closed form (`MAT-4`'s 10 % band) on a phantom path that assembles
  bitwise the same material and operator as the explicit map
  (`20260930T003229Z_MAT-1-reanchor-W1.log:485–498`); zero-data PEC and
  `NATURAL` each converge to their own manufactured field at rate 0.994 /
  0.996 with exactly 18n² constrained edges under PEC and none under
  `NATURAL` (`20260930T022051Z_TH-3-W1.log:280–289`). **What this is not:**
  one frequency and one rung for the presets (64 MHz, degree 1); a box at
  degree 1 with homogeneous data for the boundary options.
- **Perfusion, the rise curve and σ(T), each on the gel ball** (`MAT-11`,
  `MAT-12`, `MAT-5` ✅ 2026-09-30, audited PASS by the 10:30 daily; lane-3
  ratification the weekly's). The steady Pennes sink meets its sinh closed
  form to 0.183 % (P1 rate 2.06, balance ≤ 3.4e-15;
  `20260930T095245Z_MAT-11-W1.log:115–121`). The FEM SAR through the
  θ-scheme meets the eigen-series to 3.47 % against the composed 11 % band,
  with the calorimetric one-step slope ρc·θ̄/Δt = q̄ to 5.0e-5 and the
  per-step energy balance at 7e-15 (`20260930T110230Z_MAT-12-W1.log:161–170`).
  With σ = σ₀[1 + αθ̄] fed back per cell, eight undamped Picard iterates of
  EM solve → SAR → Robin ball land 2.32 % from the lumped-root profile
  against a composed 11.7 % band, the frozen-σ₀ solve missing by 29.2 %
  (`20260930T124701Z_MAT-5-step2-W2.log:175–183`). **What this is not:**
  the 3.5 % and 2.3 % figures mostly re-measure `MAT-4`'s +3.5 % SAR offset
  — those bands would pass a ≲ 5 % scale error, and the tight statements
  are the identities and the slope; the sink is real build, uniform
  source; the σ(T) loop lives in a test, not in `src/` (its per-cell
  transfer identity sees a planted absent / shifted / reversed transfer at
  1.4e-3 / 4.4e-4 / 9.9e-4 against 1e-10, not a missing ghost scatter,
  `20260930T172753Z_MAT-5-Dprime-controls.log:158–164`); an imposed uniform
  field in a 1 cm sphere, never a coil; no tissue, dose or ASTM F2182 claim.
  **The transient sink** (`MAT-13` ✅ 2026-09-30, audited PASS by the 18:00
  daily): the full Pennes equation in the θ-scheme meets the shifted
  eigen-series to 0.164 % (CN rate 2.002, balance ≤ 4.6e-14), its exponent
  checked independently by a time-convolution of `MAT-9`'s series to
  1.8e-15 K (`20260930T172102Z_MAT-13-W1.log:92–111`); real build, uniform
  source, one w — the SAR-driven perfused case is not gated.
- **Coil-driven heating of the gel phantom** (`MAT-15` ✅ 2026-10-02,
  audited PASS by the 10:30 daily; Phase 8 subgoal (8b)). Step 1 uses the
  coil's own SAR on `ANS-2`'s 10 MHz fixture as the source:
  - the θ-scheme's discrete energy identity closes to round-off;
  - the first-step slope matches ⟨SAR⟩ / c;
  - the four drives' maps are C4-gated
    (`20261001T123309Z_MAT-15-step1-W1.log:1924–1954`).

  Step 2(b) solves on the phantom's own finite-cylinder mesh with a uniform
  source and Robin walls, against the double eigen-series:
  - LSQ rate 2.15 and finest interior 7.4e-4;
  - Crank–Nicolson transient 3.5e-4;
  - backward Euler reaches the steady state at 2.2e-14
    (`20261002T093152Z_MAT-15-step2b-W4-steady.log:7174–7180`,
    `…093610Z_MAT-15-step2b-W5-transient.log:1878–1883`).

  On the tuned copper F-human (`xl` entry 42, read by the 2026-10-19
  weekly) every step-6 identity is green at the degree-2 root, and the
  ball-mean θ(t\*) per accepted watt moves **−1.27 %** from degree 1
  (3.507742e-02 against 3.553024e-02 K/W,
  `20261019T070008Z_MAT-15-step6-P22-d2root.log:6883` against
  `20261018T094815Z_MAT-15-step6b-W3-regate.log:6392`), so P17's
  thermal figures are quotable at degree 1 with that move beside them; the
  dose at degree 2 is P40.

  **What this is not:** a P1 gel-only result with no perfusion and no
  tissue. There is no absolute temperature claim on a tuned coil. **External
  check (`ANS-9` step 1, adjudicated by the 2026-10-10 weekly against AEDT
  HFSS First Order + Icepak FEA on the 10 MHz gel fixture): AGREE on the EM →
  thermal hand-off** — ⟨θ⟩ within 0.8 %, initial slope within 0.5 %, surface
  heat within 1.4 % (bands 5 %), probes A–C within 2 % — and **DISAGREE at
  probe D** (ours 11 % below AED against a 10 % band; known-issues
  2026-10-10, diagnosis `ANS-9` step 1b). **Re-read by the 2026-10-12 weekly
  (steps 1b / 1c): the probe-D gap is discretisation on our side — θ_D rises
  13 % with phantom refinement at degree 1 and a further 8 % with the
  element order. Re-adjudicated by the 2026-10-14 weekly on the re-recorded
  column (step 1d: degree 2, exact rectangles, h_p = 0.0025 m — the
  production convention): every row AGREES — probe D within 2.5 %, probes
  A–C within 2.3 % (band 10 %, B > D > C held), ⟨θ⟩ within 0.8 %, slope
  within 0.5 %, surface heat within 1.3 % (bands 5 %); `ANS-9` ✅, the
  DISAGREE retired (known-issues 2026-10-10). Still one gel fixture at
  10 MHz, one-way coupling, no tissue, no absolute temperature claim.**
- **A metal rod in the coil-driven phantom** (`IMP-3` ✅ 2026-10-02,
  audited PASS by the 18:00 daily; Phase 7 subgoal (7b), on `IMP-1` /
  `IMP-2`'s closed forms and `GEO-36`'s generator). A PEC rod inside the
  F-small phantom of the 64 MHz four-port birdcage, with and without:
  reciprocity ≤ 1.5e-15, σ_max < 1, the exact power partition per drive
  ≤ 7e-15, the rod wall's n × E = 0 on 7 437 dofs, on the record mesh and
  the 1 g rung (`20261002T170843Z_IMP-3-step1-W2.log:3798–3839`). The tip
  1 g SAR ratio is printed (P1 4.50, quadrature 3.96) beside `IMP-2`'s
  estimate, never asserted. **What this is not:** identities on one
  untuned fixture at degree 1; the σ_y mirror does not discriminate the
  rod (the phantom takes ≈ 3e-4 of the power), and the ratio's accuracy is
  `ANS-10`'s to adjudicate — no heating, dose or ASTM F2182 claim.
- **The consistency diagnostics that can fail** (`POST-3` ✅ 2026-09-30,
  adjudicated by the 18:00 daily on the 0.11 image). Poynting real-power
  balance 4.13 % at 24³ with the σ-blind solve missing by 95 %, piecewise σ
  and μᵣ through both legs, the total-current divergence at rate 0.94 with
  a CG1 vacuity control (`20260930T200032Z_POST-3-closing-W1.log:88–150`,
  `…200244Z_POST-3-closing-W2.log:89–106`). **What this is not:** a
  structured-box fixture at first-order rates; reciprocity is `PORT-1`'s.
- **Symmetry planes on the lumped-port path** (`TH-16` ✅ 2026-09-30,
  audited PASS by the 10:30 daily; **the gate was re-anchored by the 03:00
  daily — the weekly ratifies**). Per-face PEC / PMC on cut planes with the
  port rescaled reproduces, from three quarter solves, the 4×4 S of the
  quarter's own mirror-built full birdcage mesh to 1.5e-12 (band 1e-9) at
  10 MHz, degree 1, with the quarter sweeps costing 0.36 GiB against the
  full's 1.24 (`20260930T093451Z_TH-16-step3c-W1.log:951,1019–1024`).
  **What this is not:** an identity on one discretisation — against an
  *independently meshed* full birdcage the same quarter differs by 1–4 %
  (the feed neighbourhood's known degree-1 sensitivity, printed `:1023`),
  so no two-mesh agreement, no 64 / 128 MHz or degree-2 statement, and no
  F-human use yet; the control is arithmetic on the class means, not a
  mis-set plane in the solver.
- **Plumbing re-anchored on closed forms** (`PORT-5`, `WF-3`, `GEO-5` ✅
  2026-09-30, audited PASS by the 10:30 daily). The N-port assembler
  returns a typed S to 1.1e-16 including a non-reciprocal circulator — the
  `[recv, drive]` index convention is not transposed
  (`20260930T095825Z_PORT-5-W1.log:53–105`); the quick-look report's eight
  printed numbers sit within 0.72 % of the plane wave's closed form and
  equal its JSON (`20260930T100020Z_WF-3-W1.log:141–150`); the phantom's
  requested mesh size is realised, its inscribed-volume deficit falling at
  rate 1.97 (`20260930T123217Z_GEO-5-W1.log:2464–2468`). **What this is
  not:** the assembler test fakes the field solve (the real-solve half is
  `PORT-20`'s); the report's `status` line and balance row are outside the
  ✅; the mesh policy is gated on the phantom leg only, to h_p = 0.005.
- **The package, first Phase 9 measurements** (`WF-13`, `WF-14`, `GEO-37` ✅
  2026-10-04, audited PASS by the 2026-10-05 03:00 daily).
  - **Install.** `pip install .` into a fresh venv runs off `/workspace` with no `PYTHONPATH`. It reproduces `TH-6`'s α / β (0.017 / 0.060 %, bit-identical to the record) and `MAG-14` (`20261005T022512Z_WF-13-W1.log:258,305`).
  - **Registry.** Every top-level `fem_em_solver.__all__` name (25) maps to collectable gate nodes (`20261005T022839Z_WF-14-W1.log:44–52`).
  - **Facet tags.** A CAD-drawn interior rectangle reaches DolfinX's facet tags complete at `-n 1 / 2 / 4`, area to 1.1e-16 (`20261005T033214Z_GEO-37-W2.log:48–52`).
  - **What this is not:**
    - an offline, dependency-resolving install;
    - a check that a registered node exercises its name (collectability only);
    - the sub-package surfaces;
    - a `dS` assembly over a native tag;
    - a `Study` API (`WF-15` next).
- `post/evaluation.py` point location (all point evaluation goes through
  it), gmsh generation + tag QA, the runner, and the logging harness.

### 2.2 Not validated — do not trust, do not extend

- **One external finding is open (the `f₀` h-rung, 2026-09-30); the 09-19
  self-class disagreement closed 2026-10-05 and the 09-19 current-drive one
  closed 2026-09-20 — read them before quoting any S, `f₀` or Q entry.**
  (1) On the two-torus **current-drive**
  route (`PORT-1`, `EX-20`, `ans:3`) the reported `S` is **not the 50 Ω
  S-matrix**: the undriven port is open, so `S_ij = b_i/a_j` is not a column
  of S. The same run's `Z` is the open-circuit matrix, is what `PORT-1`'s
  gate of record reads, and agrees with the external code (`ANS-3` AGREE on
  the primary row). Every lumped-sheet-route matrix is unaffected.
  **FIXED 2026-09-20 by `PORT-20` (✅, steps 1–3): that route reports
  `S = z_to_s(Z)`, the two S-derived records are re-recorded and the
  known-issues entry is retired.** One residue, known-issues 2026-09-20:
  `ans:3`'s *generated* `metrics.json` / `COMPARISON.md` still show the
  pre-correction S table, because re-running that case in a slot would
  rewrite its gitignored private comparison. (2) The **order-matched**
  `ANS-4` comparison (our degree 2
  against AED First Order) AGREES at 128 MHz and on the coupling classes at
  all three frequencies, and **DISAGREES on the self class `S₁₁` at 64 MHz
  (the pre-registered rule fired) and at 10 MHz**, where the earlier
  degree-1 AGREE is re-read as partly coincidental. Neither rung is yet
  shown `h`-converged (`xl` 3e / 3g). Until that is resolved **no absolute
  `S₁₁`, `Z_in`, match or tuned-`S₁₁` figure is licensed at 10 or 64 MHz at
  either order**; the circuit-layer identities (`PORT-14`/`PORT-15`) are
  unaffected because both of their sides are ours. Diagnosis: `PORT-21`,
  `ANS-6`'s PEC column. Numbers private; known-issues 2026-09-19; §10.
  **The leading explanation is a hypothesis, not a result** (2026-09-30
  weekly): the κ-corrected sheet width at degree 2 on the ×1 cut closes the
  10 MHz self-class gap to 0.8 % — one frequency, one mesh; `PORT-21`
  steps 4–5 and `xl` entries 22 / 23 test it (§10 chain F, rule T6).
  **Re-read by the 2026-10-03 weekly (F9): degree 2 on exact SPEC
  rectangles, κ off, agrees with AED First Order to ≤ 0.73 % on every class
  at 10 / 64 / 128 MHz and with Palace at matched order (0.0000 %), on one
  mesh** — the self-class gap was the ragged port sheet, which κ had been
  correcting. **CLOSED 2026-10-05 (weekly, T6 decided on F18a / F18b): the
  h = 0.0075 exact rung (`xl` entry 27) and the conductor-halved rung
  (items 212 / 225) stay within 1 % of AED First Order on every class
  (0.97 % / 0.77 % worst), so the licence now reads — degree 2 on exact CAD
  rectangles, κ off, is the convention for every quoted lumped-sheet S
  figure, and on F-small that figure is within 1 % of AED First Order on
  every class at 10 / 64 / 128 MHz on three meshes and equals Palace at
  matched order; the "no absolute `S₁₁` / `Z_in`" sentence above is
  superseded for degree 2 on exact rectangles.** **The h-statement, written
  by the 2026-10-12 weekly (F27, `xl` entry 30, `20261011T070008Z_ANS-4-step6c.log`):
  on the third rung (h = 0.005, 593 478 cells) the self-class change from
  h = 0.0075 is 0.014 / 0.057 / 0.066 % at 10 / 64 / 128 MHz (adjacent and
  opposite ≤ 0.11 %), so the degree-2 exact F-small four-port S is
  h-converged to ≤ 0.1 % on three rungs, and the h = 0.005 rung is still
  within 1 % of AED First Order on every class — the 10 MHz self class is the
  worst and does not move with h, so that ≈ 1 % is a level, not a resolution
  effect.** Still not claimed: any degree-1 absolute; any `Z_in` or match
  figure at human scale. `PORT-21` ✅; known-issues 2026-09-19 retired
  2026-10-05. **Copper, the same way (`ANS-6`, F26 / F19, 2026-10-12
  weekly): the exact hole route at degree 2 agrees with AED First Order's
  Finite-Conductivity and PEC columns within 0.33 % on every class at 10 /
  64 / 128 MHz on the ×1 mesh and within 0.21 % at h = 0.0075 at 64 MHz,
  the h-move toward the column; `ANS-6` ✅ — the loss-partition AGREE of
  09-26 and the absolute copper S now both stand externally.**
  Palace at matched order (`OPS-61`, v0.18.1, order 2 on the shared meshes)
  ties ours to ≤ 2.6e-10 per class on every rung bought (×1, 0.0075 and
  0.005 strip; ×1 and 0.0075 hole at 64 MHz) — a cross-solver identity on
  our discretisation, not an accuracy level (the tracked `COMPARISON.md`
  sections; F40, `20261015T112551Z_ANS-4-item476-ties.log`).
  (3) **No tuned `f₀`, Q, `Z_in` or match claim** (`ANS-7`, known-issues
  2026-09-30): the mode and the copper contribution agree (§2.1). At degree
  1 on one rung our `f₀` sat 1.2–1.9 % above AED's columns; **the 2026-10-03
  weekly (T13) re-read it as degree-1 discretisation: at degree 2 on exact
  rectangles `f₀` is within 0.11 % of AED First Order and Q within 1.9 %,
  on the order-matched ×1 mesh — an agreement on one mesh, not a tuning
  claim; the h-rung is owed (T12′)** and known-issues 2026-09-30 retires at
  T13, not before (§10 chain T). **T13 decided 2026-10-10 (weekly): `f₀`
  AGREE with its h-statement (§2.1); `R_in` and Q move away from AED with h
  (Q −1.9 → −3.1 %). T41 decided 2026-10-17 (weekly): Q and `R_in` converge
  in h (ratios 0.717 / 0.705 on the steps 0.015 / 0.01 / 0.005 m, p ≈ 1.65;
  Richardson limits 0.24 % / 0.23 % below the h = 0.005 values) to ≈ 3 %
  below AED First Order, the level Palace reproduces to six digits — `ANS-7`
  ✅ on its Done-when with the absolute Q INCONCLUSIVE pending T48; no Q,
  `Z_in` or match figure is licensed.**
  (4) **The lumped feed sheet's terminal accounting misses the field's
  sheet dissipation by the sheet-field variance — a documented
  sheet-resolution property with a measured width law** (known-issues
  2026-10-11, re-titled 2026-10-22 by §10 H32 on H30's reading, **retired
  as a record 2026-10-26 on H33**: on the F-small exact-rectangle
  four-port ΣX_s / P_acc falls monotonically with the sheet width at both
  orders — 4.274e-04 / 9.958e-05 / 4.135e-05 at w₀/4 / w₀/2 / w₀ at degree
  1, exponent −1.685; 8.826e-04 / 1.683e-04 / 4.838e-05 at degree 2,
  exponent ≈ −2.09 — larger at degree 2 at every width,
  `20261025T095010Z_ANS-4-H33a-W1.log:5662–5665` and the H33b TABLE lines;
  `ports:23`): on the tuned copper F-human each driven sheet
  dissipates X_s more in the field than its lumped Norton term
  ½ Re(1/Z_th)\|V\|², where X_s = CS_s, the Cauchy–Schwarz excess
  ½ Re(Y_s)(∫\|n̂×E\|² dA − \|∫E·ĥ dA\|²/A_s) of a field that varies
  across the sheet (CS_s / X_s = 1.000000000 at degree 1,
  `20261011T170627Z_ANS-8-v3b-V3.log:6204–6208`). ΣX_s / P_acc reads
  8.316822e-03 at degree 1 and 2.149715e-02 at degree 2 at the root
  (`20261022T093054Z_ANS-8-H30-W1.log:89`) — two points, one mesh, not a
  trend. It is a property of the lumped-sheet model, not of the solve or
  the ESR term; every per-watt figure states its normaliser — the lumped
  line cut P_acc,SPEC (weekly ruling 2026-10-12), with the field partition
  printed beside it — and no assertion equates field and terminal feed
  power. The 10-24 condition's first clause (H33's width law at both
  orders) was met 2026-10-25 and the entry retired 2026-10-26; H34 (the
  Palace twin) was withdrawn — no field export exists (§10 epitaph).

- **The birdcage's port S-matrix at 10 / 64 / 128 MHz is a gated
  self-consistency identity set, not an absolute-accuracy result; no
  coil is tuned or resonant, and no Larmor-frequency figure has been
  compared against anything outside this code.** *(Rewritten 2026-08-26
  18:00 review — `PORT-11` ✅ at both Larmor frequencies; the 08-25 head
  "ports at 10 MHz only" is retired, see §2.1. What remains unvalidated
  and is what this bullet now covers: the **absolute accuracy** of any
  S/Z/coupling figure at 64/128 MHz (the gates are reciprocity, passivity
  and C4 symmetry — a wrong-by-a-constant-factor port model passes all
  three; the independent check is the AED benchmark `ANS-4` — its runnable
  half landed 2026-08-30, the operator's AED replication landed 2026-09-04,
  **adjudicated 2026-09-06: AGREE at 10 MHz, INCONCLUSIVE at 64/128 MHz;
  re-adjudicated 2026-09-13 (weekly): AGREE at 128 MHz on the
  order-matched degree-2 rung, AGREE by mechanism at 64 MHz pending
  `ANS-4` step 3** — so the *order-matched* answer is externally checked,
  while the **degree-1 production fixture's own 128 MHz entries are
  measured 5–7 % from it** (our own order move, 2b/2d) and every identity
  gate passes on them: this bullet stands for degree-1 figures exactly as
  written, and the production-order question is `TH-19` step 3), any
  resonance, mode-spectrum or
  tuning claim (Phase 6, no lumped capacitors exist), and B1+/SAR on a
  solved coil field. Two readings on record, gated by nothing: the C4
  spreads grow ≈1.7× per Larmor step (0.055 → 0.057 → 0.101%) and
  `|Im P|/Re P` at the driven port rises 0.34 → 1.76 → 2.66 — stored
  energy, physics not noise, but the trend is what a resolution study at
  128 MHz would have to explain before any tighter band is written.
  Older text, kept for the record: the old head of this bullet on 08-25 was
  "The birdcage has ports at 10 MHz only; nothing is validated at a
  Larmor frequency", itself replacing "No coil or birdcage has ports".
  The 16-leg *mesh* is gated (`GEO-19` ✅ 2026-08-25), the 32-ring-port
  layout is built on both sheet orientations (`GEO-20` step 2 / `EX-35`
  transverse; `GEO-26` ✅ 2026-09-04 longitudinal, 270 728 cells, all 32
  sheets exact to 1e-12), and **`PORT-13` step 1 (2026-09-04) put the first
  field on it: one column at 10 MHz, degree 1, P17 driven with 31 ports
  terminated at 50 Ω — power accounting residual 9.68e-3 inside the imported
  1e-2, the two opposite ports agreeing to 0.35%, one solve 28 s at `-n 8`
  (`20260904T050538Z_PORT-13.log:10750–10763`), and **step 2 (2026-09-04,
  04:30 slot) made it a 4×4 sub-block** — four drives over the one mesh (P17,
  its top-ring mirror P33, the two opposites P25 / P41): reciprocity of the
  sub-block 4.118219e-13 vs the unmoved 1e-3 (machine-level because the
  operator is complex-symmetric; the 1%-column control at 7.045× the band is
  what gives it teeth), column passivity 0.9158–0.9160 under 1, the
  top/bottom mirror identity worst pair 0.0308% vs 5% on two independently
  solved columns, the power residual 9.68e-3 on all four columns to six
  digits (`20260904T093638Z_PORT-13.log:10786–10854`), and **step 3
  (2026-09-04, 17:00 slot) assembled the full 32×32**: the 32×32 exists and
  is reciprocal (‖S − Sᵀ‖_F/‖S‖_F = 5.446798e-13 vs the unmoved 1e-3, with
  the 1%-column control at 2.496× the band), passive (σ_max = 0.999999452 ≤
  1) and C16 × mirror symmetric (all 18 measured classes inside 5%, worst
  0.4426%) to those digits, over 32 drives whose power residuals all sit at
  9.33–9.68e-3 inside the imported 1e-2
  (`20260904T171419Z_PORT-13.log:77–119`). That is a **network-level
  self-consistency reading on one fixture at 10 MHz, degree 1** — no σ_max
  record, no absolute-accuracy claim, no resonance, tuning or mode-spectrum
  claim exists at 16 legs.** The 10 MHz result is a **port-model validation
  on a gated fixture**, not an MRI-regime result; do not quote it as one.)* The history below is kept
  because every band and convention in the port model was set in it. The
  birdcage-port direction is scoped —
  `PORT-9` (lumped/circuit-element port BC, Jin ch. 11) — and as of
  2026-08-16 both named prerequisites have **executed and closed**:
  `PORT-10` (the two systematics compose additively, cross-term
  −0.0604 pp) and `GEO-15` (graded conductor sizing reaches 0.967 of CAD
  mass; `PORT-9` budgets from 98 k cells — measured on 0.7.2; the gate
  behind that figure was non-executing on `main` from the 0.11 merge until
  `GEO-21` ✅ 2026-08-26, which greened it on a **coarse-graded** control —
  0.966977 at 98 666 cells now reproduces the 0.967 / 98 k figure, but the
  live gate measures fine-vs-coarse grading, and the graded-vs-ungraded
  answer to this prerequisite question remains the 0.7.2 close). `PORT-9` **step 1 closed
  2026-08-17**: the mesh prerequisite `GEO-16` landed (longitudinal sheet
  on an opt-in kwarg, area = CAD to roundoff), the parked formulation
  branch was merged with its six exact identities green, and the first
  lumped-port `Z` was solved on the two-torus fixture at 10 MHz —
  `Im Z₁₂ = 0.829782 × ωM₁₂` against the gated gap route's `0.894310`,
  a cross-route deviation of **7.7095%**. **Step 2 executed 2026-08-17:
  both pre-stated bands MISS** (cross-route 7.7095% vs 5%, lumped mutual
  12.6931% vs 10%; the gap route stays inside at 6.0391%), neither
  widened, and the miss is **diagnosed** — it is the transverse average
  over the full-width sheet (7.7783 pp of it; path/projection residual
  only 0.0763 pp), a property of the two feed definitions on this box,
  not of the solver or the mesh. **Step 2b executed 2026-08-17 (12:00
  slot): the band HOLDS at the narrowed width** — ladder 7.7095% →
  3.6730% → **1.8333%** at f = 0.5 against the unmoved 5% band, f = 1.0
  reproducing step 2's record to < 1e-4. Step 2's gate is closed **at
  the narrowed definition**: a lumped port sheet is specified by its
  interior width fraction f and its width is measured as `w = A/h` on
  the filtered facet set — the convention is part of the port model's
  spec. **Step 2c executed 2026-08-18 (22:30 slot): the reciprocity leg
  is closed** — the lumped-sheet route landed in
  `run_n_port_sparameter_sweep` and the two-torus sweep is reciprocal at
  `‖S−Sᵀ‖/‖S‖ = 2.574249e-11` against the unmoved 1e-3 band, cross-route
  1.6079% / 1.5950% inside the 5% band, with **0.23 pp of drive
  dependence** off step 2b's impressed-gap 1.8333% — a lumped reading
  should be quoted with its drive stated. Step 3's gate (i) prerequisite
  is discharged. **Step 3 executed its preflight 2026-08-19/20 and is
  🚫-blocked on the mesh**: the birdcage fixture has no port-sheet facet
  (global facet set exactly `{1}`) and its port boxes have **no
  terminals** — conductor-facing area exactly 0.000000e+00 m² on all
  four, because the coil is uncut and the boxes float in air outside it.
  The prerequisite was `GEO-18` (cut the legs, boxes straddling the cuts;
  commissioned 2026-08-20) and it **closed 2026-08-22**: the gapped
  birdcage now has planar disk terminals (0.988616 of `2·π·r_leg²`) and
  a port sheet per leg gap (area = `dx·g` to 1.000000000000, C4 spread
  8.470e-16) — a mesh with somewhere to put a port. **Step 3 closed
  2026-08-25** across legs (c)/(d0)/(d)/(d1′): see §2.1's birdcage-port
  bullet for what it licenses. B1+ remains §10 subgoal 4; its port
  prerequisite `PORT-11` (the same three gates at 64/128 MHz) **closed
  2026-08-26**, so B1+ on the solved birdcage field is no longer blocked
  on ports — it is `WF-6` (steps 1–2b ✅ by 2026-08-31: `|B₁⁺|` maps and
  five symmetry identities at 10/64/128 MHz on the loaded F-small fixture,
  CG1 covariance ≈ 1.9–2.2% vs 5%; **identities on one unconverged mesh —
  no homogeneity, absolute, tuning or SAR claim**). Its B₁⁺ target is **a
  convergence statement, registered by step 5 (✅ 2026-09-13)**: on the
  unloaded F-small birdcage at 10 MHz the worst-radius C4 four-copy spread
  of the CG1 `|B₁⁺|` falls **5.2506 % → 2.0719 %** (ratio 0.3946) as the
  global resolution goes 0.015 → 0.012 m, the C4 covariance falls
  3.6159 % → 1.6815 % alongside it, and the power residual stays inside the
  unmoved 1 % band on both rungs — spreads reproduced at rtol 1e-3 and both
  falls asserted (`20260913T170259Z_WF-6-step5.log:1987, 3831, 3836`). The
  closed-form comparand (odd-order cube sum `S₁₁`) misses the ×1 rung by
  ≈ 3.2 % at the centre and ≈ 5.3 / 7.9 % at 0.4R / 0.5R — a record; the
  third rung (×0.0095) stalls near 2 % and its ≈ 2 % common-mode power
  residual is a banked, unexplained negative (known-issues). **No
  closed-form, homogeneity, absolute, C95.3 or Larmor B₁⁺ claim exists**
  (§10 epitaph, 2026-09-13). Its step 3
  (2026-08-31) read the coil-driven **point-SAR** map off the primal
  N1curl `E` and measured it **25–41%** off the same five identities
  against the same 5% band — the pointwise-`E` estimator's own floor,
  filed as five deliberate reds on `main` (known-issues). Step 3b
  (2026-08-31) put the CG1-`E` estimator beside it and read **worse**
  (152 / 110 / 170 / 53 / 41%) with the primal column reproducing to every
  digit — and its own diagnostics (`‖E_cg1 − E‖/‖E‖` **1876%** over the
  phantom, CG1 phantom power **+35 199%**) say `project_to_cg1` does not
  fit an N1curl `E` at all, so that is now the open finding. **No SAR
  claim exists**; the `|B₁⁺|` gates project `B` and are untouched.
- **Coil loading at the Larmor frequencies is measured flat in f on one
  XL window — a record, not a gate** *(head re-worded 2026-09-13 weekly
  review; it read "is an extrapolation" until `TH-11` step 5d's 64 MHz
  bracket landed, see the end of this bullet)*. The
  apparent frequency trend is now attributed: `TH-11` step 4's fixed-f
  h-ladders read **flat in f** — the h → 0 brackets [−2.15, −0.91]% at
  10 MHz and [−3.37, −0.38]% at 30 MHz overlap at ≈−1%, so the monotone
  three-point set 1.58 / 5.59 / 10.27% was the resolution term, not
  physics. That is printed evidence adjudicated by the 2026-08-17 review,
  not a gate, and 64 MHz itself still has **no h → 0 bracket** (finest
  rung +2.81% at 2.52 cells/δ). **The degree-1 h-ladder to that bracket
  is closed as a measured negative** (`TH-11` step 5, adjudicated
  2026-08-18 10:30 review): the memory wall is superlinear in cells —
  0.42 M cells comfortable, **0.99 M pegs `memory.peak` at
  `memory.max` = 64.00 GiB**, 2.81 M OOMs at every legal rank count
  (MUMPS factor fill-in) — and an affordable third rung would need a
  refinement ratio ≈ 1.2 whose difference signal sits at the 0.01 pp
  run-to-run floor, so the fit would be noise. **The degree-2 axis is
  now measured too** (`TH-12` step 2, 2026-08-18): on the 138 619-cell
  coil fixture at 10 MHz, degree 2 walks the coarse-rung ΔR deviation
  +1.5834% → **−0.8508%** — h → 0 quality on an unrefined mesh — but at
  **61.94 GiB summed peak RSS (96.8% of `memory.max`)**, i.e. against
  the same wall, and the 64 MHz case needs ≈2.5× the cells at fixed
  cells/δ. **Adjudicated 2026-08-18 18:00 review: no affordable
  (order, h) route to a gated 64 MHz h → 0 bracket exists on this box**
  — there is no rung swap to scope, and the degree-2 reading stands as
  corroborating evidence from the order axis, not a gate. This bullet
  moves only when a review adjudicates a gated 64 MHz bracket (which
  now requires either more memory or an out-of-core/iterative solver
  path, neither scoped).
  **Pointer, 2026-09-09 18:00 daily review — the premise above is now
  measured false, and this review does not move the bullet.** "More
  memory" arrived: `TH-11` step 5d ran the 2 808 204-cell 64 MHz third
  rung to completion on the 512 GiB XL service — `Status 0`, 4838 s at
  `-n 8`, 17 passed / 1 skipped, 263.4 GiB `memory.peak` — a container-lifetime
  maximum spanning the killed first attempt, **not this run's attributable
  peak** (ledger note; `OPS-43` (c) measured `memory.peak` unable to separate
  two runs, `20260910T050613Z_OPS-43c.log:90–95`) (`docs/testing/xl-ledger.md`,
  `20260909T153910Z_TH-11-step5d.log`) — giving the ladder
  **+10.2698% → +2.8063% → +0.3824%** at h 0.005 / 0.0025 / 0.00125 and a
  two-rung bracket **[−2.0415%, −0.4256%]** that overlaps step 4's 10 MHz
  [−2.15, −0.91] and 30 MHz [−3.37, −0.38]. So the sentence "no affordable
  (order, h) route exists **on this box**" was a statement about a 64 GiB
  ceiling. Whether that bracket is *gated* — and therefore whether this
  bullet moves — is **the 2026-09-13 weekly's ruling**, assigned there by
  the operator in `f5071f6` and not taken here: the daily review does not
  own §2's capability claims, and three caveats stand in the way (Dodd–Deeds
  is quasi-static and is the *comparison* rather than the reference at
  64 MHz; the two-rung `p_eff` = 2.876 disagrees with the three-rung
  `p` = 1.623, so the bracket is the honest object and the point estimate is
  not; and `d₀` = −0.7834% sits just outside the 10 MHz bracket's upper end).
  Recorded here so no reader takes the paragraph above as still current.
  **Ruled 2026-09-13 (weekly review): the bracket is a record, not a
  gate**, for exactly those three caveats plus a fourth — a 4838 s XL
  window is not a CI assertion. What it licenses: ΔR against the
  quasi-static Dodd–Deeds kernel is **measured flat in f across
  10–64 MHz to the bracket width (≈ 1.6 pp)**, confirming the 08-23
  epitaph's physics half on the box that could afford the rung; the
  bullet's head moves from "extrapolation" to "measured, not gated" and
  nothing asserts it (§10 2026-09-13).
- **An absolute coil-driven SAR is now externally checked at 10 MHz on
  F-small (`ANS-2` AGREE, 2026-09-18): whole-phantom dissipated power to
  5 % of HFSS (AED lower), pointwise SAR at the four ball centres within
  ±10 % (scatter on our side); normalised to accepted power the codes agree
  to ≈ 1–2 %, so the residual is the lumped-port feed systematic, not the
  SAR operator. The *compliance* SAR number — the IEEE C95.3 *limit*
  claim — is still open** (it needs a body model, a real input power and a
  real coil), and remains open after
  `MAT-4` step 4 (2026-09-06) moved that row to ✅: what step 4 gates is
  the **10 g mass-averaged C4 covariance identity** on the coil-driven
  field (0.3303 / 0.0756 / 0.0574 / 0.3132% at the unmoved 5% band,
  `tests/validation/test_birdcage_sar_mass_averaged.py`), a
  self-consistency statement on one fixture at fixed `h` at 10 MHz. No
  drive normalisation, no `h` convergence, no Larmor rung and no
  compliance comparison exists; the 1 g column misses the same band at
  6.84% on that rung and is a printed record there, while on the 2.5 mm
  phantom rung (`MAT-4` step 5b, 2026-09-08,
  `tests/validation/test_birdcage_sar_1g_rung.py`) the 1 g column **is
  gated** at the same unmoved 5% band, reading 0.0957–0.1305% — a second
  self-consistency identity on the same fixture at the same frequency, and
  still not an absolute or compliance number. The first
  coil-driven SAR *symmetry* readings exist (`WF-6` step 3, 2026-08-31)
  and miss their band by 5–8× read **pointwise** off the primal `E`. Since
  `WF-6` step 3h (2026-09-02) the repo gates exactly one coil-driven SAR
  quantity: **a C4 symmetry identity of quadrant powers on one fixture at
  10 MHz at fixed `h`** (the twelve integral pairs ≤ 1.52% against the
  imported 5% band, `tests/validation/test_birdcage_sar_integral.py`).
  Since step 3i (2026-09-03) the mirror identity through each drive's own
  azimuth is gated on the same four solves — the four single-drive flank
  pairs ≤ 1.76% against the same band, the flank-vs-opposite control at
  38% — a single-drive statement, so independent evidence rather than a
  re-reading of the rotation identity. Those are self-consistency
  identities and nothing else — no absolute SAR, no homogeneity, no
  C95.3, no Larmor, no convergence claim.
- **`⚠️` chunks** (`TH-2`;
  `WF-2` left the list ✅ 2026-10-09 (§9 item 346, provisional — the 10-10 weekly ratifies),
  `PORT-8`, `PORT-4` and `MAT-1` left the list ✅ 2026-09-29,
  `TH-3` ✅ 2026-09-30, `POST-2` ✅ 2026-09-30 (§9 item 92),
  `PORT-5` ✅ 2026-09-30 (§9 item 98), `WF-3` ✅ 2026-09-30 (§9 item 99);
  `POST-2` added 2026-09-29 18:00 review — §7 carried it ⚠️, this list
  had omitted it) carry code whose tests assert too little to protect
  them; revalidate before building on them. `OPS-17` (§7, ✅ 2026-08-21)
  removed or replaced the finiteness-only tests themselves and reconciled
  the complex-suite baseline (216 of 232 validation tests observed in
  completed complex runs; the two absentees are formally deferred with
  reasons); the `⚠️` glyph stays until each chunk's own physics is
  revalidated.

### 2.3 Before you debug a failing test

**Check [`docs/testing/known-issues.md`](docs/testing/known-issues.md)
first.** It lists every currently-failing test with symptom, the commit it
was verified failing at, and the diagnosed cause. Several tests fail on
`main` for reasons unrelated to any change you are making.
---

## 3. Status legend

| Symbol | Meaning |
|---|---|
| ✅ VERIFIED | Test executed, assertion is quantitative, passing |
| 🟡 IN PROGRESS | Actively being implemented |
| 🧪 UNVERIFIED | Code landed; test has never actually executed anywhere |
| ⚠️ PLACEHOLDER-BACKED | Implemented and "green", but the green rests on tests that assert too little (§2.2) |
| ⬜ NOT STARTED | |
| 🚫 BLOCKED | Cannot proceed; blocker named in the chunk |

`⚠️` is not a bug report against the chunk's own code — it means the chunk
cannot be trusted until revalidated against the real solve (`OPS-17`, ✅
2026-08-21, replaced the tests; the glyph itself retires per chunk as each
one's physics is revalidated). Do not "fix" a `⚠️` chunk by
loosening its test.

**A measurement-only step is `🧪`, never `✅`** *(clarified 2026-08-02, 18:00
review, after an audit found two probe steps carrying `✅`)*. Several chunks are
split into a probe step whose stated product is "measurements, assert nothing"
and a gate step that turns those measurements into bounds. The probe step lands
real code and real logs, but §4.3 asks for an executed quantitative *assertion*,
and a script that only prints has none — so the probe step is `🧪` until its
gate lands, however good its numbers are. This is a clarification of §4, not an
exception to it: the fix is to write the gate, never to widen §4. Demoted on
this ground 2026-08-02: `PORT-1` step 1, `MAT-6` step 2a. (`PORT-1` step 1 was
restored to ✅ the same day when its gate, step 2, landed — which is exactly
the intended route back.)

---

## 4. Definition of done

A chunk is `✅ VERIFIED` when **all** of the following hold:

1. **Code and docs are committed.**
2. **The agent executed the verification command itself** and recorded the result
   via `scripts/testing/run_and_log.sh`. A chunk does not park on a human.
3. **The assertion is quantitative** — at least one check compares against:
   - a closed-form analytic solution, with a stated tolerance
   - a measured convergence rate under h- or p-refinement
   - a conservation, reciprocity, or symmetry identity
   - a documented reference value from literature or a prior validated run

   > **Finiteness-only gate:** a chunk may **not** be `✅` if every assertion it
   > adds is "is finite", "is non-zero", "has shape N", or "did not raise". Those
   > are welcome as *additional* assertions and insufficient as the *only* ones.
4. **Its runtime tier is declared and its measured elapsed time recorded.**
5. **Any dependency on a `⚠️` chunk is stated explicitly** in the chunk entry.

Never loosen a failing assertion to make a test pass. A failing analytic
comparison is evidence about the test as much as about the code.

### When human verification *is* required

Reserve it for judgments a test cannot make, and name the judgment: "does this B1+
map look physically plausible in ParaView?", "is this refinement acceptable near
the port faces?", "does this S11 curve resemble published birdcage behavior?"
Never for routine pass/fail.

---

## 5. Execution policy

### 5.1 Compute budget — shared machine

The development server is a **shared 36-core box; this project may use at most 12
cores.** Every verification command declares a tier and must not exceed it:

| Tier | Ceiling | Use for |
|---|---|---|
| `smoke` | 30 s | Imports, pure-Python logic, config validation |
| `standard` | 3 min | Coarse meshes, single small solves — the default |
| `heavy` | 20 min | Convergence studies, sweeps — must be labeled `heavy` |
| `xl` | **4 h, 512 GiB, 16 ranks — six runs per trailing 7 days, i.e. every night Sun–Fri** (4 h since 2026-09-13, was 2 h: step 2d's 7225 s window is the measured reason; a full window now ends 06:00 and overlaps the 04:30 slot, accepted) (operator directive 2026-09-13 — the windows cost no tokens and the box is quietest at 02:00; three from 2026-09-10; was one, raised because the measured constraint turned out to be wall-clock rather than memory) | The convergence rung the heavy tier cannot hold (operator directive 2026-09-05). Runs against the separate `fem-em-solver-xl` compose service (profile `xl`), **commissioned by the weekly planning review** — pre-registered in `docs/testing/xl-pending.md` with a named chunk and step, the exact command, a *measured* price and the readout — **and queued by the daily review as clerk** into `docs/testing/<tier>-queue.env` whenever the queue is empty and the ledger budget allows (operator directive 2026-09-13, daily-review.md step 6b; before that only the weekly could queue, so at most one of the three windows a week ever ran), and recorded in `docs/testing/xl-ledger.md`; the bash guard denies a second XL exec inside 7 days, any XL exec outside `run_and_log.sh`, any `timeout` above 7200 s, and `-n` above 16. First slot reserved for the `ANS-4` 128 MHz refinement rung |
| `xxl` | **8 h, 754 GiB, 16 ranks — one run per trailing 7 days, Saturday 02:00** (operator directive 2026-09-10) | The window a 2 h `xl` slot cannot hold: the human-scale coil, or a many-drive sweep before the factorisation is reused. Runs against the separate `fem-em-solver-xxl` service (profile `xxl`), from `docs/testing/xxl-queue.d/` via `xl-run.sh xxl`, recorded in `docs/testing/xxl-ledger.md`. 754 GiB is this WSL VM's **whole** allocation — half the machine's 1.5 TB — so the kernel OOM killer reaches a runaway before the cgroup does; that is the operator's stated choice and is recorded in `docker-compose.yml` rather than silently softened. Saturday is its own cron day so the two tiers cannot both start at 02:00 |

- Wrap commands in `timeout -k 30 <s>` at the tier ceiling — the `-k` is
  mandatory: a plain TERM does not reliably stop an `mpiexec` job, and an
  overrun can wedge the container (MAT-6 step 10, 2026-08-12; known-issues has
  the recovery recipe). **If a run overruns, kill it and
  redesign the case smaller.** Never re-run with a longer timeout. 20 minutes is a
  hard ceiling for any single compute command regardless of tier.
- **`mpiexec -n 12` is the hard rank ceiling.** Use the smallest count that fits
  the tier — more ranks on a fixed mesh stop paying quickly, and a bigger mesh is
  usually a better use of the budget than more ranks on a small one. Keep `-n 2`
  for anything a rank-local bug could hide in (`MAG-11` was a missing allreduce
  visible only under MPI), and note that CI runs at `-n 2`, so a test that only
  passes at wider ranks is not CI-portable.
- **Memory ceiling: the container is capped at 128 GiB** (`docker-compose.yml`
  `deploy.resources.limits.memory`; raised from 64 GiB by operator directive
  2026-08-24, host has 754 GiB). Read it from
  `/sys/fs/cgroup/memory.max` rather than assuming — `memory.peak` does not
  exist at container level on this box, so track `memory.current` between
  commands. **Memory is a tier dimension like time:** a case that pegs the
  ceiling is redesigned smaller, never granted more. Raising the limit is an
  operator decision, not a chunk's (§9 standing constraint).
  **⚠️ Several recorded negatives were measured against the *old* 64 GiB wall
  and are not evidence about this box any more** — `TH-11` step 5 (2.81 M
  cells OOM, 0.99 M pegged), `TH-12`'s degree-2 wall (61.94 GiB on the coarse
  rung), `OPS-17`'s `coil_loading_degree2` deferral, and the §10 epitaph
  killing the coil-loading-trend target. None of them is automatically
  reopened: re-pricing one is a **review** decision, and any revival needs its
  finest rung priced first (the epitaph's own lesson). Do not cite "does not
  fit the box" from a pre-2026-08-24 measurement without re-measuring.
- **XL windows run at 02:00 from cron, not from a session (operator directive
  2026-09-09) — and that is the whole answer to box contention.** The box is
  shared with work this sandbox cannot see, so scheduling a quiet hour beats
  detection that cannot work. `scripts/automation/xl-run.sh` at 02:00 daily
  runs the first `*.env` file in `docs/testing/<tier>-queue.d/` (a FIFO
  since 2026-09-15 — the one-slot `<tier>-queue.env` it replaces lost every
  window after a review-less night, and is honoured only as a legacy
  fallback), checks the ledger budget before taking it, and deletes the
  file afterwards; an empty directory is the normal state and the entry
  then costs a second. **No Claude session is involved**, which is also what makes
  it possible: a foreground harness call is capped at 660 s and an implementer
  slot is killed at 65 min, so no scheduled *session* can hold a 2 h window,
  while a cron script has neither limit. It takes its own lock rather than the
  automation flock, so a long window cannot starve the 02:15 weekly or the
  03:00 review, and it restarts the service first so `memory.peak` belongs to
  the run. A failed window still consumes its queue entry, for the same reason
  a started run has spent the week: re-queue deliberately or not at all.
- **This box cannot be measured from inside the sandbox (measured 2026-09-09) —
  which is *why* the schedule above exists rather than a detector.** This is WSL2: `/proc/loadavg`,
  `nproc` and `free` describe **this Linux VM only**. Work on the Windows host,
  in another WSL distro, or in another VM is invisible, and it cannot be probed
  — there is no `powershell.exe` on PATH and `/mnt` is unreadable. The
  measurement that settles it: the operator's Task Manager showed **100 % CPU
  with two solves running** while this VM reported **load 0.84** and both
  project containers **0.00 %**. A load-average check would have called the box
  free and started a 16-rank, 2-hour window straight into it. Consequences:
  (a) `scripts/testing/box_check.sh` reports containers, VM load and VM memory
  and prints its own blind spot, and its green verdict means only "nothing *I*
  can see is using the box"; (b) a mandatory per-run confirmation token was
  tried and **removed the same day** — it blocked ordinary work, since any file
  or command merely *containing* the XL command matches the guard's trip-wire,
  so it prevented editing the guard and even writing the queue file; the 02:00
  schedule replaces it and `FEM_EM_XL_BOX_OK` survives only as something the
  launcher sets; (c) the
  docker socket is reachable from a direct agent command but **not from inside
  a script** here, so `box_check.sh`'s container section goes blank exactly when
  run the convenient way and says "unavailable" rather than reporting zero.
- **An XL window waits for the box (operator directive 2026-09-09).** It holds
  up to 16 of 36 cores for up to two hours, so starting one onto a machine
  somebody else is already using is the rudest thing this project can do.
  `bash_guard.py` now refuses an XL command when the 1-minute load average
  leaves fewer free cores than the window's rank count plus 4, and the denial
  says what the load is and what it must fall below. Two properties worth
  knowing: it **fails open** — a guard that cannot read `/proc/loadavg` must
  not become the reason nothing can run — and load average is a decaying mean
  that cannot tell our own load from anyone else's, so for a minute or two
  after one of our runs ends it will postpone the next one. That is the safe
  direction; `FEM_EM_XL_IGNORE_LOAD=1` is the operator's escape hatch when the
  load is known to be ours and decaying, and never a way past somebody else's
  work. Verified across the threshold: at 36 cores and 16 ranks it allows at
  load 12, postpones at 16.1, and the same load allows an 8-rank window.
- **The big-compute slots are a budget, not a loophole — and, since
  2026-09-15, one that is meant to be spent every night.** Six `xl` runs and
  one `xxl` run per trailing 7 days, commissioned by the weekly review on the
  measurements the heavy tier cannot hold, with the daily review holding a
  backlog floor of ≥ 4 `xl` and ≥ 1 `xxl` entries ahead in
  `docs/testing/xl-pending.md` and a licence to fill a shortfall with
  **priced-family variants** and **cost probes** only (daily-review.md step
  6b; operator directive 2026-09-15, after four of five nightly windows
  fired empty in the first wind-down week); the
  ledger row is written by the harness *when the run starts* — since
  2026-09-19 into the gitignored sidecar `logs/<tier>-ledger.pending`,
  folded into the tracked ledger when the window ends (a killed wrapper's
  row is folded in by the next window, elapsed blank), so that a window
  overlapping the 03:00 review or an implementer slot no longer leaves
  `main` dirty for its whole length — so a killed
  run has still spent its slot **if it consumed box time**. The budget counts
  elapsed seconds, not attempts: a row whose elapsed column is 0 or empty
  never started and is not charged (2026-09-10 — two of five rows that day
  were 0-second environment failures, and pricing those like a two-hour solve
  would be absurd). That cannot be gamed the way a retry could, because the
  harness measures the elapsed time itself. Bring the service up for the slot
  (`docker compose -f docker/docker-compose.yml --profile xl up -d
  fem-em-solver-xl`) and stop it afterwards. Every other rule above (kill and
  shrink, `timeout -k 30`, rank-local bugs at `-n 2` first) applies to it
  unchanged.
- **A long window must survive its client (measured 2026-09-09).** `TH-11`
  step 5d's first attempt was killed at ≈40 min *on the wrapper side*: the
  harness process died, so the log got no `## Exit` footer and nothing in it
  was countable — while the container-side `timeout`, `mpiexec` and all eight
  ranks kept running at ≈85% CPU on 260 GiB with **nobody consuming their
  output**. Two rules follow, and they apply to any window over ≈10 minutes,
  not only XL ones:
  * **Redirect container-side to a file first, then echo it back** —
    `... pytest ... > /workspace/logs/<name>-raw.log 2>&1; rc=$?; echo
    "[capture] rc=$rc" >> /workspace/logs/<name>-raw.log; cat
    /workspace/logs/<name>-raw.log; exit $rc`. `/logs/` is gitignored, so it
    never dirties a preflight. The `[capture] rc=` line puts the exit status
    *in the file*, not only on the dead client's stdout. If the client dies,
    the file still holds everything; once the container command has finished,
    `scripts/testing/run_and_log.sh --capture-orphan <ID>
    logs/<name>-raw.log` turns it into a normal log with a real `## Exit`
    (`Status:` from the last rc line, one test-results row, exit = that
    status; **76** `NO_RC_CAPTURE` and `Status: unknown (no rc line)` when the
    file has none). **Mechanised 2026-09-10 (`OPS-43` (a))**; the escaped
    host-side form is in `run_and_log.sh`'s header, and the gate is
    `scripts/testing/test_durable_capture.sh`. The killed wrapper's own log
    has no `## Exit`: cite the recovered one. **In normal mode the footer
    honours a final rc line (`OPS-45`, 2026-09-11):** if the command exits 0
    but its last non-blank output line is `[capture] rc=<n>` with n ≠ 0, the
    Status, row Exit and harness exit are n and `## Exit` carries a
    `Capture note:` — so an idiom that drops `; exit $rc` can no longer
    record green; a non-zero command exit is never overwritten (gate
    `scripts/testing/test_capture_status.sh`). Still write the trailing
    `exit $rc`. A bare pipe or `tee` does **not** do this: when the client
    goes, `tee` takes SIGPIPE and the run dies with it.
  * **After any killed window, check for orphaned ranks before doing anything
    else** — `docker exec <service> ps -eo pid,etime,pcpu,comm`. A dead
    wrapper does not stop the compute. Kill the `timeout` supervisor, verify
    the rank count is zero, and only then re-run. An unnoticed orphan holds
    cores and memory on a shared box indefinitely. **Mechanised 2026-09-10
    (`OPS-43` (b)):** `run_and_log.sh` now lists the target service's
    processes before any `docker compose … exec` command and exits **75**
    (`ORPHAN_REFUSAL`), printing PID / elapsed / args and writing no log or
    row, while any `mpiexec|hydra_pmi_proxy|pytest` survives. A 75 means
    "clean up first", never "retry".
- **`memory.peak` is per *container lifetime*, not per run, and cannot be
  reset on this kernel** (WSL2 6.6; `/sys/fs/cgroup/memory.peak` is
  read-only, measured 2026-09-09). A reading taken after a second run is the
  **maximum over every run since the container started**, so either restart
  the service immediately before the window — which zeroes it — or have the
  module print its own `ru_maxrss`. Label any figure that cannot be
  attributed to one run as what it is; `ANS-4` step 2b's ledger row carries
  exactly that caveat.
- **No chunk is marked `xl` without a *measured* memory price.** An
  extrapolated one is how the 2026-09-09 slot went to a case that measured
  249 s and ≈16.6 GiB — an eighth of the ordinary service's own limit — while
  the genuinely memory-bound case waited (§7 `ANS-4` step 2b). A censored
  reading is not a measurement either: `TH-11` step 5's "0.99 M cells pegged
  at `memory.max` = 64.00 GiB" was clipped at the limit, so it bounded demand
  from below and licensed no extrapolation at all. The honest price for a
  case that has only ever been killed is **unmeasured**, and the first XL
  window on it is a cost probe whose readout says so. A **priced-family
  variant** (same script and mesh as a measured row; frequency, degree,
  drive or port set, rank count or one env knob varied) inherits the
  family's measured price with the scaling stated, and is the one class
  the daily review may commission without a new measurement (2026-09-15).
- Record real elapsed time in `docs/testing/test-results.md`.
- **A tier is a measurement, not an intention.** A chunk whose runtime has never
  been measured is `unmeasured`. Cost-probe first: build the mesh, print the cell
  count, solve a tiny case, extrapolate, then size the real case to fit.

Cost is dominated by mesh size: ≈8×10³ cells solve in seconds, ≈4×10⁵ cells take
minutes, and CI cannot host the `heavy` tier.

### 5.2 Agent autonomy and loop hygiene

- Agents implement **and verify**. Human-gated completion is prohibited (§4.2).
- **No-op guard:** if a work cycle produces only documentation edits and executes
  no verification command, stop and escalate rather than commit an audit note.
- **Do not append duplicate status blocks.** Status lives in §7 tables. (The
  legacy human-gated queue `docs/testing/pending-tests.md` and its
  `AWAITING-HUMAN-TEST` status were removed 2026-08-04 — verification is
  agent-executed per §4; the old queue survives in git history.)

### 5.3 Verification environment — Docker

All verification runs inside the `fem-em-solver` container. Preflight
`docker compose -f docker/docker-compose.yml ps` — STATUS must be "Up"; `exec`
fails with `service ... is not running` otherwise, and that is a *setup* error,
not a test failure. Do **not** use the `cd docker && ...` form: a `cd` inside a
compound command prompts for permission regardless of the allowlist, which fails
in scheduled sessions.

One-time setup: `docker compose -f docker/docker-compose.yml build` (≈4 min), then
`... up -d`. The repo mounts at `/workspace`, so source edits need no rebuild.

| Property | Value |
|---|---|
| Image / base | `fem-em-solver:latest` on `dolfinx/dolfinx:v0.11.0` |
| dolfinx / Python | `0.11.0.post0` / 3.12.3 |
| numpy / gmsh | 2.4.6 / 4.15.2-git-657c8e9 |
| h5py / HDF5 | 3.16.0 / 2.1.1 (built from source against the image's HDF5) |
| petsc4py / mpi4py | 3.25.1 / 4.1.2 |
| Default PETSc scalar | `numpy.float64` — **real mode** |
| Complex build | `/usr/local/dolfinx-complex` |
| Memory cap | 64 GB (raised from 16 GB by `MAT-6` step 7, operator-approved 2026-08-10) |

Upgraded 0.7.2 → 0.11.0 by `OPS-18` (2026-08-23, step 3b), measured on the
built image in `20260823T200740Z_OPS-18-step3b-env-probe.log`. The compose
`PYTHONPATH` literal (`.../lib/python3.12/dist-packages`) is the **only**
version-encoded path in the project — the `dolfinx-{real,complex}-mode`
wrappers and `src/sitecustomize.py` derive the tag themselves. The gate
`tests/environment/test_dolfinx_version.py` asserts the adopted version, and
derives the expected Python tag from `sys.version_info`, so the next upgrade
fails loudly rather than silently importing the wrong tree.

**Complex scalars**, required by `TH-1`:
`source /usr/local/bin/dolfinx-complex-mode` (→ `numpy.complex128`);
`dolfinx-real-mode` switches back. The chunk commands set
`PYTHONPATH=/workspace/src`, which drops the container's dolfinx path;
`src/sitecustomize.py` re-appends the dist-packages directory matching the active
`PETSC_ARCH`, so imports resolve in **either** mode. Set
`FEM_EM_REQUIRE_COMPLEX=1` on any run that is supposed to be complex, so the
real-mode skips in `tests/environment/test_complex_mode.py` become failures, and
put `tests/environment` first in the pytest path list.

Run everything through the harness so results land in
`docs/testing/test-results.md` and `docs/testing/logs/`:

```bash
scripts/testing/run_and_log.sh <CHUNK-ID> "docker compose exec -T fem-em-solver bash -lc '...'"
```

Use `exec -T` — without it `exec` allocates a TTY and can hang under an agent. The
harness exports `COMPOSE_FILE` itself; a bare `docker compose` outside `docker/`
needs `-f docker/docker-compose.yml`.

Logs are retained under `docs/testing/retention-policy.md` (adopted
2026-09-03): the harness collapses gmsh optimisation chatter at write time,
every log is gzipped at 7 days, and logs not cited from this plan, its archive,
known-issues, or an example guide are deleted at 14 days by
`scripts/maintenance/housekeeping.py --apply`. Cite the log basename in the §7
closure annotation — that citation is what keeps it. Read `.log.gz` with `zcat`.

### 5.4 Examples and Ansys cross-validation

- **A second independent reference (operator directive 2026-09-20): Palace.**
  Beside the private AED columns, the open-source Palace solver (`OPS-61`) is
  a reference the project *may publish* — its numbers can live in tracked
  `COMPARISON.md` columns and in gate bands, which AED numbers never can.
  The rule that a benchmark is a judge and never the solver applies to it
  unchanged: nothing in `src/` is built on Palace. Commissioning a Palace
  case is the weekly's, like an AED case; the first is the `ANS-4` fixture,
  to adjudicate the ≈ 3 % self-class offset `ANS-6` isolated to the port
  model.

- **Palace is the public external reference, long term (operator directive
  2026-10-02, interactive session).** Extends the 2026-09-20 bullet above:
  Palace (open source, Apache-2.0) is the external check the project
  **publishes at full detail**, the open counterpart to AED's private
  columns, and the comparison the repository leads with once it is public
  (the operator intends it public before 2027-01-10). Rules:
  (i) **full numbers are publishable** — Palace values, our values and their
  gaps go in tracked `COMPARISON.md` columns, example guides, gate bands and
  external talks; each record cites the Palace version, the config and the
  mesh it ran on.
  (ii) **matched order, shared mesh** — a published our-vs-Palace figure
  states both element orders and is at matched order on one mesh unless it
  is labelled an order or mesh study (`OPS-61` step 1c: 0.0000 % at matched
  order; the 6.7 / 2.5 / 4.5 % degree-2-vs-order-1 gaps are an order
  mismatch, not a disagreement).
  (iii) **reproducible** — a published Palace case keeps its Palace config,
  the shared mesh (or its generator) and our run script tracked together,
  so a reader can rerun both sides.
  (iv) **never a bridge to AED** — no tracked file sets a Palace number
  against an AED number, states a Palace-vs-AED gap, or otherwise lets an
  AED value be recovered; the three-way reading (ours / Palace / AED) is
  done only in the gitignored `COMPARISON_private.md` / `docs/private/`,
  and its public trace is the qualitative verdict alone.
  (v) the benchmark-is-a-judge rule stands — nothing in `src/` is built on
  Palace.
  **For the weekly:** keep Palace as a standing lane — where an `ANS-*`
  case's physics is in Palace's reach (driven ports, lumped RLC,
  surface-impedance conductors, eigenmodes), commission its Palace twin
  beside the AED request, so every externally checked claim has a
  publishable column. Commissioning and pacing stay the weekly's.

- **`examples/` is a maintained product surface, not a scratch area.** Each
  phase keeps **at least five** clean, runnable examples demonstrating its
  capability from distinct angles (different geometries, materials, drives, or
  output quantities — five trivial variations of one case do not count),
  executed via `./run_examples.sh` and producing combined-XDMF output that
  opens in ParaView — this is how the human operator reviews progress
  independently of the test suite. When a chunk changes what an example
  demonstrates, the same commit updates the example. A broken example is a
  defect (known-issues discipline applies). **Examples accrue with gate
  closures, not at phase end**: each time a chunk closes a quantitative gate
  (§4 ✅), the next daily review enqueues a standalone example chunk
  demonstrating that newly gated capability (daily-review.md step 5). The
  audited bar for an in-progress phase is therefore a ramp —
  `examples ≥ min(5, gating chunks closed ✅)` — checked by the weekly review
  (weekly-review.md step 4); the flat five binds once the phase completes.
  *(Tightened 2026-09-09, weekly review: the count is **gating chunks or
  gated steps closed ✅**, whichever granularity the phase's §7 entries use.
  Reason — the 2026-09-09 ramp check measured Phase 5's requirement as
  `min(5, 0) = 0`, because its only live gating chunk `WF-6` is 🟡 while
  carrying five ✅ gated steps (1, 2, 2b, 3h, 3i); a bar of zero cannot bind
  anything, which is the opposite of what this rule is for. Under the
  tightened count Phase 5 owes 5 and has 6 — still no shortfall, and now a
  number that binds. This raises the bar and never lowers it: a phase whose
  chunks are all ✅ counts exactly as before.)*
  *(Clarified 2026-08-09, weekly review: the ramp binds the physics phases
  1+. Phase 0's infrastructure capability — Docker, CI, the harness, the
  runner — is exercised by every example run and does not owe a separate
  example set; its meshing slice is covered by the `mesh:` group.)*
  Example chunks are their own §7 entries sized for one implementer run,
  never riders on physics chunks, and never target ungated capability.
  *(Operator directive 2026-08-10:)* **every runnable example also ships
  with a same-stem step-by-step guide page** — what it demonstrates, how
  to run it, and how to analyze the output step by step so the operator
  can understand what is going on without reading the source (required
  sections and the mechanical gate: §7 `EX-15`). A missing or stale guide
  is a defect like a broken example; the doc-reference checker enforces
  presence and structure, and the weekly review audits guides with the
  ramp. Example chunks scoped after 2026-08-10 include the guide page in
  their done-when; the pre-existing examples are backfilled by `EX-15`.
- **Ansys benchmark cases** live in `examples/ansys_benchmarks/<case>/`, each
  containing: `SPEC.md`, precise enough to replicate in Ansys Electronics
  Desktop with no judgement calls (geometry with dimensions, materials,
  boundary conditions, port definitions, frequencies, mesh guidance, and
  exactly which quantities to export); the runnable script; our results
  (metrics JSON + XDMF); and `COMPARISON.md` with our numbers filled in and
  blank columns for the AED numbers.
- **Cadence is the weekly planning review's call**
  (docs/automation/weekly-review.md) — roughly one case per phase milestone,
  and only on gated capability. A benchmark on ungated physics wastes a
  licence-hour measuring noise.
- **Returned AED numbers are adjudicated by the next weekly review**: recorded
  in the case's `COMPARISON.md`, promoted into §7 gates where they agree, and
  opened as known-issues/chunks where they disagree — a disagreement with the
  commercial solver is a finding to diagnose, never to explain away.

---


**Tier label for host-runner example legs (ruled 2026-09-02 weekly review).** An `EX-*` leg executed through `./run_examples.sh` is tiered by the per-window ceiling its §7 entry states, not by §5.1's 180 s pytest figure — write it as "standard (host-runner, ≤ 500 s window)"; a window past its stated ceiling is an overrun, a window past 180 s is not. Measured: `EX-30` 105 / 447 / 935 s, `EX-36` 141 / 228 / 182 s under the same label.

## 6. Phase map

| Phase | Goal | Gating chunks | State |
|---|---|---|---|
| 0 | Infrastructure, packaging, CI, meshing | `OPS-1`, `OPS-2` | Done |
| 1 | Magnetostatics + analytic validation | `MAG-1`…`MAG-6` | **Complete and trustworthy** |
| 2 | Time-harmonic Maxwell, complex materials, ABC/PML | `TH-1`…`TH-9` | In progress — every analytic gate closed (`TH-1`/`TH-6`/`TH-7`/`TH-8`/`TH-9` ✅); Larmor sphere `TH-10` ✅; coil trend `TH-11` ✅ closed on a measured negative (no 64 MHz h → 0 bracket fits the box, 2026-08-18); degree-2 `TH-12` ✅ (closed 2026-09-02 on the re-affirmed production-order clause: degree 1 coil-fed, degree 2 imposed-field), `TH-13` ✅ 2026-08-31 (the injector is the degree-1-only source projection; coil degree-2 identity reds stay open in known-issues); **`TH-19` steps 1–2 ✅ 2026-09-10, outcome (a) — the matched projection clears the degree-2 coil identity by five orders on both σ-halves; production order unchanged until step 3 tests the birdcage's sheet drive (§10 2026-09-13)**; **`TH-19` step 3 ✅ 2026-09-13 (both degree-2 identities green on the sheet drive at 10 / 128 MHz); `TH-14` ✅ 2026-09-14 (copper via Leontovich: cavity Q +0.010 %, Dodd–Deeds floor −0.299 %, the copper birdcage's identities); 2026-09-19 weekly: degree 2 is the target order and the default is not flipped — decision pre-registered for 09-26 on `xl` 3e / 3g (§10)**; **2026-09-26 weekly: not flipped — the 10 MHz ladder is not converged by the rule's own letter, and the rule's ratio test is withdrawn as mis-posed (§10 epitaph); re-tested on chain F7 / F8**; **2026-09-30 weekly: `TH-3` ✅ and `TH-4` ✅ 2026-09-29 (PEC and natural walls each converge to their own exact field; rates 0.994 / 0.996), `POST-2` ✅, `POST-3` ✅ on the 0.11 image; ladder Tier B — `TH-16` ✅ (symmetry planes: three quarter solves equal the 4×4 of the quarter's mirrored mesh to 1.5e-12), `TH-18` ✅ (layered impedance at σ = 10 S/m, ratified), `TH-5` 🟡 (first-order ABC: closed forms at k₀R = π to 1.8e-3; a reflector for electrically small boxes — step 2 re-posed, PML unopened); `TH-2` ⚠️ leaves on one logged window (§10 rulings)**; **2026-10-07 weekly: `TH-5` 2b met the mixed-wall bar (0.25 pp against 1.95 pp); 2a's (R) retired as a sphere-model band on a cube ladder (r(W) falls 0.48 → 0.23 as the natural wall's does), its (S) stands; ✅ on step 2c (module landed, §2.1 line)**; **2026-10-17 weekly: `TH-19` ✅ ruled — steps 1–3 green with asserted identities, the order question it reopened answered by T6 (the daily writes the row with its banked Done-when)** |
| 3 | Material models, phantoms, SAR | `MAT-1`…`MAT-6` | `MAT-2` ✅; **`MAT-4` externally anchored 2026-09-18 — `ANS-2` AGREE at 10 MHz on F-small (5 % integral / ±10 % pointwise, residual = feed; no compliance claim)**; `MAT-6` ✅ (ΔR to **0.2747%** on the production projected drive and the 418 888-cell slab-refined fixture promoted by step 11, 2026-09-06 — 1.5834% was the pre-refinement record; eddy-current regime; **externally checked 2026-09-02 — `ANS-1` adjudicated AGREE against Maxwell 3D, numbers private; the weekly re-checks that verdict against the moved column**); SAR gated on an **imposed** uniform field only (`MAT-4` steps 1+3: lossy-sphere closed form 3.5%, mass-averaging exact at 1 g/10 g) — **`MAT-4` ✅ 2026-09-06 — step 4 puts the C95.3 mass-averaging operator on the *coil-driven* field: 10 g C4 identity 0.3303 / 0.0756 / 0.0574 / 0.3132% at the unmoved 5% band on F-small at 10 MHz, whole-phantom coverage identity 1.58e-14; the 1 g column is a printed record at 6.84% and the absolute / compliance claim is still open**; **step 5 + 5b close the 1 g column as a gate 2026-09-08 — the four 1 g C4 pairs 0.0957 / 0.1199 / 0.1305 / 0.1065% asserted at the same unmoved 5% band on `GEO-27`'s `h_p` = 0.0025 m rung (199 920 / 58 866 cells, a version-tagged record at the imported `CELL_COUNT_BAND`), both steps audited PASS** *(step 5 by the 2026-09-08 10:30 review; step 5b by the 2026-09-09 03:00 review — this cell was written by the 02:15 weekly 35 min before 5b's audit actually ran, and the audit then returned PASS, so the claim is true as of 03:00 and was premature when written)* — still a symmetry identity at one frequency on one fixture, still no absolute, C95.3-compliance, homogeneity or Larmor-SAR claim; **2026-09-30 weekly: `MAT-1` ✅ 2026-09-29 (gel presets on the lossy sphere within 3.6 %), `MAT-5` ✅ 2026-09-30 (σ(T) fed back to a fixed point, 2.3 % against the lumped root), `MAT-3` 🟡 (Debye through the plane wave; Cole–Cole tissues wait for the Gabriel table, operator) — all on the imposed field** |
| 4 | Coil modeling, lumped elements, ports, S-params | `PORT-1`…`PORT-11` | `PORT-1` ✅ 2026-08-15 (field-derived S through the package, two-torus fixture only, two named systematics); `PORT-10` ✅ 08-16; **`PORT-9` ✅ 2026-08-25 at 10 MHz on the gapped 4-leg birdcage** — leg (d1′)'s geometric negative control passed on the power-wave route (displaced classes 6.2219 / 7.1142 / 2.8474% vs the tightened (iii′) 0.5%, reciprocity 2.259e-14 vs 1e-3, 2.466e+11× from the pre-fix 5.57e-03), no Larmor/resonance/tuning claim; history: steps 1–2c ✅ on the two-torus (lumped-sheet BC, 1.8333% cross-route, reciprocity 2.6e-11), step 3 on the gapped birdcage has two gated legs (c)/(d0) at 10 MHz (C4 spread 0.0152–0.0159% vs 5%, 50 Ω termination separates the classes 598× — re-recorded image-tagged on the 0.11 image 2026-08-24 by leg (d3c): 0.0359%, 253.2002×) and **leg (d) closed 2026-08-23 — the 4×4 passes all three gates** (reciprocity 2.495292352e-05 vs 1e-3, σ_max 0.862659137 ≤ 1, class spreads 0.0199 / 0.0180 / 0.0108% vs 5%, gate (iii) since tightened to 0.5%); leg (d1)'s geometric control ran 2026-08-23 and **found the route loses reciprocity (5.57e-03 vs 1e-3) once the fixture is asymmetric**; leg (d2) (asymmetric two-torus, 13:30 slot) traced it to the assembly — the readout *is* the source's adjoint (1.33e-10), the asymmetry is the terminated-`Z` per-column normalisation — and the 18:00 review ruled the power-wave S fix (leg (d3)) with the class re-record (d3b), (d1′) serial on (d3b); **`PORT-11` ✅ 2026-08-26 — the same three gates at 64 and 128 MHz on the same fixture** (64: 2.581325834e-14 / σ_max 0.999721388 / spreads 0.0573 / 0.0599 / 0.0370%; 128: 7.030990825e-15 / 0.998974779 / 0.1012 / 0.0916 / 0.0654%, cells/λ 12.5024 ≥ 10 enforced; audited COMPLIANT 18:00 review) — self-consistency identities only, no absolute-accuracy/resonance/tuning claim; **externally checked — `ANS-4` adjudicated AGREE at 10 MHz against HFSS (2026-09-06) and AGREE at 64/128 MHz (2026-09-13 weekly: on the order-matched degree-2 rung at 128 MHz, by mechanism at 64 MHz pending `ANS-4` step 3; numbers private; the degree-1 gate fixture's own 128 MHz entries sit 5–7 % from their order-matched value, invisible to the identity gates)**; `PORT-13` ✅ 2026-09-04 — the 32×32 on the 16-leg / 32-ring-port longitudinal rung passes the same three gates (reciprocity 5.4e-13, σ_max 0.9999995, 18 C16 × mirror classes ≤ 0.45%), self-consistency only; **2026-09-19 weekly — two external findings, numbers private: (1) `ANS-3` AGREE on the mutual-coupling row, but `PORT-1`'s current-drive route reports an `S` that is not the 50 Ω S-matrix (its `Z` is right; `PORT-20`, known-issues); (2) the order-matched `ANS-4` rungs: 128 MHz AGREE stands (3b), the 64 MHz rule fired on the self class (the "by mechanism" AGREE is withdrawn for `S₁₁`, kept for the couplings) and the 10 MHz degree-1 AGREE is re-read as partly coincidental (3c) — no absolute `S₁₁` / `Z_in` claim at 10 or 64 MHz at either order (`PORT-21`, known-issues)**; **2026-09-30 weekly: `PORT-4` … `PORT-8` ✅ 2026-09-29 / 30 — the 🧪 / ⚠️ port rows re-gated against the real solve (Y symmetry 8.5e-15, typed S to 1.1e-16, the sweep planner bitwise, Touchstone through `scikit-rf` to 6.3e-14, `PORT-8`'s S′ = DSD identity exact); `ANS-7` adjudicated — the tuned resonance agrees on the mode and the copper contribution, `f₀` 1.2–1.9 % above AED at degree 1 (known-issues 2026-09-30); the κ-corrected sheet width is the leading hypothesis for the `S₁₁` gap, under test (§10 chain F)**; **2026-10-05 weekly: T6 decided — degree 2 on exact CAD rectangles, κ off, is the convention; within 1 % of AED First Order on every class at 10 / 64 / 128 MHz on three meshes (×1, h = 0.0075, conductor ×0.5) and = Palace at matched order; `PORT-21` ✅, `OPS-61` ✅; the h = 0.005 rung (`xl` entry 30) owes the convergence statement** |
| 5 | Full MRI system: loaded birdcage, B1+, SAR maps | `WF-5`…`WF-8` | Blocked on Phases 2–4 for excitation; both meshes (coil+phantom, birdcage) generate and are identity-gated in CI (`GEO-9`, 2026-08-03); the birdcage fixture is loaded (phantom inside) and since `GEO-18` ✅ 2026-08-22 has terminals and port sheets — **`WF-6` steps 1–2b ✅ (2026-08-30/31): `\|B₁⁺\|` maps on the loaded F-small birdcage symmetry-gated at CG1 at 10, 64 and 128 MHz, in ParaView (`EX-38`/`39`/`40`), no homogeneity/absolute/tuning claim; coil-driven SAR is measured, not gated — steps 3–3e: the packaged phantom-restricted `E` estimator is honest (best-approximation and power anchors) and the five SAR identities still miss 5% at 6–9.5%, verdict (c), the ≈1 cm phantom cells; **step 3f (2026-09-02) halves the phantom's `h` and all five land inside the band at 2.5–3.5%, clause (a) — verdict (c) confirmed, no gate registered in-slot, and the `\|B₁⁺\|` identities turn out not to be mesh-converged either (2.19% → 0.62%, three deliberate reds + known-issues)**; **step 3g (same day) reads the C4 identities as cell integrals of the primal `σ\|E\|²` on the coarse mesh and lands all twelve pairs at ≤ 1.52% — the construction, not `h`, was the binding mechanism**; **step 3h ✅ (2026-09-02) registers the repo's first coil-driven SAR gate — a C4 symmetry identity of quadrant powers on one fixture at 10 MHz at fixed `h` (twelve integral pairs ≤ 1.52% against the imported, unmoved 5% band), with the five pointwise asserts retired to records that still exceed the band; no mirror identity, no absolute SAR, no homogeneity, no C95.3, no Larmor, no convergence claim, and `WF-6` stays 🟡**; **2026-09-09 weekly: the SAR half of the F-small exit is DONE and the B₁⁺ half is now the single blocker.** The 09-06 review's "one item away — `MAT-4` step 2" landed 2026-09-06 as **step 4** (the step was renumbered, not skipped) and step 5b extended it to 1 g on 09-08, so the 09-06 watch condition ("if `MAT-4` step 2 is not ✅ by the 09-09 weekly, re-plan the SAR route") is **met, three days early**. What remains is subgoal 4's other half — the B₁⁺ closed-form gate, `WF-6` step 4 — which has taken **six attempts in five days (4, 4a, 4b, 4c, 4d, 4e) and landed two records and four parked negatives, zero gates**; the 09-08 10:30 ruling finished the comparand (odd-order cube sum `S_11`, bar ≈ 0.5–1.1%) and put the residual miss on the FEM, so the item is now **`WF-6` step 4f — the `h`-ladder**, priced by `GEO-29` (§10 Phase-5 assessment 2026-09-09); **2026-09-13 weekly: the 09-09 rule fired — no closed-form gate landed, and what the `h`-ladder measured is a monotone fall of the C4 four-copy spread (5.2506 → 2.0719 % on `main`, step 4g; the third rung stalls and its power residual is a banked negative), so subgoal 4's B₁⁺ target is re-scoped to that convergence statement (epitaph in §10) and Phase 5 exits on F-small when `WF-6` step 5 registers it — one tests-only item, ≈ 2026-09-14…15 if its first window is green**; **`WF-6` step 5 ✅ 2026-09-13 (12:00 slot) — see the §7 row**; **2026-09-19 weekly: Phase 5 has exited on F-small (`WF-6` ✅; `ANS-2` AGREE 2026-09-18, ratified) on the scope §10 states — no homogeneity, compliance, Larmor-SAR or human-scale claim; the human-scale maps are Phase 6's chain H**; **2026-09-26 weekly: `WF-7` ✅ 2026-09-25 — the first human-scale safety quantity (10 g SAR hotspot on the 32-port F-human birdcage at 64 MHz, C16 0.49 %), degree 1, one phantom, *untuned* coil, symmetry-gated only; the tuned coil is `WF-10`** |
| 6 | Birdcage tuning at 64/128 MHz: mode spectrum, lumped capacitors, circuit co-simulation (the HFSS + Circuit split); production target: **32-port high-pass birdcage at 1.5 T** (§10 operator directive 2026-08-17); **fixture scale re-directed 2026-08-25 — two fixtures, F-small (today's 0.07 m gate fixture, records frozen) and F-human (≈ 0.15 m radius / 0.30 m long high-pass, the deliverable fixture); the `N ≤ 25` ceiling was arithmetic on the wrong radius and dissolves at human scale. Full directive in §10 Phase 6 — the 2026-08-30 weekly review must dispose of it, cost probe first** | **Gating chunks, listed explicitly 2026-09-09 weekly so §5.4's ramp is mechanically checkable on this phase (it was the one row with no `X-1…X-n` range): `GEO-19`, `GEO-20`, `GEO-25`, `GEO-26`, `PORT-13`, `PORT-14`, `PORT-15`, `TH-17`** — six closed ✅ at chunk-or-gated-step granularity, so the ramp is 5 and the corpus carries 8. Subgoals owned by the weekly review (§10); mesh prerequisites `GEO-19` (16 legs, cost rung) + `GEO-20` (ring-gap ports) scoped 2026-08-23, + `GEO-26` (longitudinal ring sheets — opened 2026-09-03 after `PORT-13` step 1 measured the `GEO-20` sheets transverse, `h = 0` for the lumped port model; the first ring-port solve waits on it) | **Started 2026-09-04 on the feature ladder** — mesh prerequisites all ✅ (`GEO-19`, `GEO-20`, `GEO-26`, `GEO-25`: the F-human rung is a gated 504 642-cell fixture, cost exponent 0.84 not r³), `PORT-13` ✅ (the 32×32), `PORT-15` step 1 ✅ (ladder-network closed form + S/Z reduction at machine precision), `POST-6` step 1b ✅; `PORT-14` 🟡 (lumped RLC sheet single-mode to 3.4e-3 — ruled a fixture record 2026-09-06, §7); conductor lineage `TH-15` steps 1/2a ✅; **2026-09-09 weekly: `TH-15` step 3a ✅ 2026-09-06 (`birdcage_port_domain(as_hole=True)` — the coil as a PEC hole) and step 2d ✅ 2026-09-07 (the gap-displacement port current on the open-circuit anchor), so `TH-17`'s mesh-side prerequisite is discharged and only `PORT-14` step 2 remains serial ahead of it; `PORT-16` ✅ 2026-09-07 (the exact discrete power identity on the 4-leg birdcage, single and superposed drives) closes the 1%-of-supplied accounting gap `POST-6` step 1b opened. `TH-15` step 2 itself is still open after seven sub-steps — 2f attributed the 2% `Z` asymmetry to the point-sampled `_path_voltage` and 2g measured the calibrated gap average as reading B; the `src/` replacement is specified but unwritten. `PORT-14` step 1e (the 09-06 ruling's `REDUCTION_FLOOR_F_SMALL` re-registration) has not been queued in three days and is the phase's oldest unstarted owed item** *(landed 2026-09-11)*; **2026-09-13 weekly: `PORT-19` ✅ (one factorisation per sweep, 11.2× on the 32×32), `TH-15` step 3b ✅ (`gap_cell_tags` landed, default unflipped; step 2's unitarity gate and step 3 still open), `PORT-14` step 2 frozen on a positive reading (a (1 + κ)-corrected width takes both 64 MHz residuals under 1e-3, unregistered — step 3 is the κ-derived registration); the tuned-birdcage chain is enumerated in §10 as ten numbered steps — internal milestone `TH-17` mode 1 at 64 MHz ≈ 2026-09-25…10-01 at the measured 2–3 days per numbered step, the AED-matched tuned `S₁₁` undated on the operator's AED queue**; **2026-09-19 weekly: seven of the ten chain steps landed 09-13/14 (`PORT-14` ✅, `PORT-15` ✅, `TH-15` step 3, `TH-19` step 3, `TH-14` ✅, `POST-6` ✅); `TH-17` step 1(i) found no eigenvalue at 64 MHz (1b queued; `PORT-22`, the driven-sweep resonance, opened as the independent route); `WF-7` 0b: human-scale order sensitivity 5.50 %, degree 2 priced at 106.1 GiB / 348 s per drive; `ANS-6` commissioned; three chains (T tuned, F feed model, H human scale) enumerated in §10 with 11 slot-executable steps — 2026-09-30 … 11-19 at the interval's two measured rates; the AED-matched tuned `S₁₁` stays undated, serial on chain F**; **2026-09-26 weekly: 10 of the 17 chain steps landed — `PORT-22` ✅ (the tuned resonance by a driven sweep, circuit = field at 11 frequencies), `ANS-6` adjudicated AGREE on the copper loss partition (within 0.4–4.2 %) with the conductor interior cleared, `TH-17` 1c green on a closed-form LC loop and its gate (a) re-posed (the four-capacitor spectrum starts at 129.6 MHz); the 10 / 64 MHz `S₁₁` gap is re-read as unconverged on our side (5.7 % / 2.5 %, the ladder moving toward AED); production order not flipped; `ANS-7` (tuned copper F-small) commissioned — adjudication at the 10-10 weekly if `PORT-22` step 2 lands by 09-30; human scale (`WF-10`, tuned F-human) opened, undated**; **2026-09-30 weekly: zero chunks, twelve gated steps in 3.98 days — `TH-17` steps 2–3 (the four-capacitor spectrum equals the circuit prediction; the k = 1 pair tuned to 64 MHz by the eigen route), `ANS-7`'s runnable half (Routes A = B to 9.3e-5), `WF-10` step 1 (F-human's circuit ↔ field identity; its first zero is not the k = 1 mode, and the mode-tuned fixture went red and was re-posed by the operator); `ANS-7` adjudicated — AGREE on the mode and the copper contribution (5 %), `f₀` 1.2–1.9 % above AED at degree 1, INCONCLUSIVE as a tuning claim (known-issues 2026-09-30); the κ-corrected sheet width is the leading hypothesis for the `S₁₁` gap (15.0 % → 0.8 % at 10 MHz, degree 2, one mesh); production order not flipped; F-human on copper (`WF-12`) and the quarter F-human (`TH-20`) opened; exit defined in §10 — earliest the 10-10 weekly, 10-14 … 10-21 at the measured re-pose rate**; first physics target `TH-17` eigenmodes. **2026-10-03 weekly: `TH-17` ✅ (10-01); on exact port rectangles at degree 2, κ off, the F-small four-port S agrees with AED First Order within 0.73 % on every class at 10 / 64 / 128 MHz and with Palace at matched order, and the tuned `f₀` within 0.11 % — one mesh, h-rungs owed (F18, T12′), so the 09-30 `f₀` finding was degree-1 discretisation; `ANS-8` mode content AGREE within 0.6 %, `C_mode` INCONCLUSIVE at degree 1; H8b banked, tuned human scale moves to copper (H12c) and has no date; F-small side decided at the 10-10 weekly if the rungs hold**. **2026-10-05 weekly: T6 decided (the absolute `S₁₁` convention on F-small — the second exit item); T13 waits on `xl` 10-08; H14 discharged (κ = 4.5e-3 on exact copper sheets is not a sheet defect), so H12c is unblocked and human scale is dated only if it is green by 10-07; T21 / T22 re-posed on C_in with a resolution band (the 1e-6 band was an identity band on a frequency target)**. **2026-10-07 weekly: T21′ / T22′ landed — tuned on C_in, the degree-2 in-model zero sits at 64.000 / 127.74 MHz within the grid's own resolution (records on a ×1 mesh); H12c fired its pre-registered red (the copper F-human's circuit ↔ field identity at its `C_mode` reads 1.58 while the 50 Ω closure is exact) — human scale has no date and its route is the in-model Larmor root (H17, `xl`) then H9′; H16 killed at the ordinary ceiling → `xl` 10-13; capacitor ESR 0.2 Ω adopted for new steps (operator 2026-10-06; `PORT-24`); T13 at the 10-10 weekly on `xl` 10-08**. **2026-10-10 weekly: T13 decided — `ANS-7`'s tuned `f₀` AGREES with AED First Order with an h-statement (0.11 % at ×1, 0.03 % at h = 0.005, converging first order in h); `R_in` / Q move away with h (−2.7 % / −3.1 %, a stated level, T30 / F31 own the loss's h-sequence). Human scale: the in-model root on copper with ESR exists (H17, C_in,F-human = 118.2437 pF, heavy tier) and the tuned two-port drive's `\|B₁⁺\|` and 10 g hotspot are gated on the feed-bisector mirror (H9′; `WF-10` ✅ on audit) — degree 1, records, no tuning claim; degree 2 priced 10-13 (`xl` 32), the degree-2 readings `xl` 10-16 (entry 35); `ANS-8` step 2 commissioned on that node**. **2026-10-12 weekly: F27 decided — the degree-2 exact F-small S is h-converged to ≤ 0.1 % on three rungs and within 1 % of AED First Order on every class; `ANS-6` ✅ (the copper S within 0.33 % of First Order on every class at three frequencies, 0.21 % on the 64 MHz h-rung) — subgoal (e)'s external copper check is complete; `ANS-8` rows 1 / 4 AGREE at degree 2 (`C_mode` / `f₀` within 0.7 % / 0.4 % of First Order on the saved σ = 800 32×32, no solve) — human-scale tuning quantities are quoted at degree 2 only; Q's h-statement waits on `xl` 10-15 (F31) → the 10-17 weekly; `PORT-24` ✅ (ESR identities), `WF-10` ✅ on audit; the F-small side of the exit is a decision at the 10-17 weekly either way**. **2026-10-14 weekly: Palace at matched order ties ours on every exact rung bought (h = 0.0075 at 10 / 64 / 128 MHz, h = 0.005 at 10 MHz, the copper hole rung) and on the ESR-tuned F-small's driven reads at ×1 and h = 0.01 at both Larmor frequencies (T42–T44, Q's ×1 → 0.01 move equal to five / six digits), so Q's h-sequence is the discretisation's — read on F31 (`xl` 10-15) → T41; the degree-2 copper F-human column priced at 514 s / 89 GiB (`xl` 10-13), so the degree-2 32×32, thermal column and quarter are `xl` 10-18 / 10-19 / 10-20 beside entry 35 (10-16); F-clinical (operator 2026-10-14, chain C) does not move the exit**. **2026-10-17 weekly: T41 — Q and `R_in` converge in h (p ≈ 1.65 on named steps) to ≈ 3 % below AED First Order, `ANS-7` ✅ on its Done-when with the absolute Q INCONCLUSIVE (T48, the operator's pass history); the F-small side of the exit is closed. Human scale: at the degree-1 node the degree-2 copper F-human is detuned (entry 35: Im Y_in +1.35e-2 S, per phantom watt ≤ 1.6 %), and its degree-2 root exists at C_in(d2) = 115.5209 pF (H26, `xxl` 10-17, one crossing on all four readings); the tuned figures at that root are `xl` 10-22 (entry 41) — Phase 6 exits at the first weekly that reads it green (10-24 at the earliest) as an internal tuning statement, `ANS-8` step 2 re-posed at the root (H29) for the operator's column**. **2026-10-19 weekly: entry 36 — the degree-2 copper 32×32 is passive, reciprocal and C16-symmetric at degree 2, and the circuit layer's `C_mode`(d2, Cu) sits 0.87 % below the in-model root (H25); `ANS-8` step 2 is posed at the root (SPEC v4, the operator's W2); `TH-15` ✅ ruled (interior PEC bodies gated on three fixtures; the hole default stays a fixture kwarg, epitaph); the exit still waits on entry 41 (10-22) → the 10-24 weekly at the earliest**. **2026-10-21 weekly: entry 37 — the four degree-2 quarter sectors reassemble against the mirror-built column at 1.7e-11 (`TH-20` ✅, the quarter is the licensed degree-2 human-scale route); entry 40 prices F-clinical-1.5T's degree-1 column (317 + 199 s, 40 GiB at `-n 8`; `GEO-38` ✅ as a probe) and C3 is commanded for `xxl` 10-24; H30 narrowed to what entry 41's record carries; the exit still waits on entry 41 (tonight) → the 10-24 weekly at the earliest**. **2026-10-24 weekly: EXITED.** Entry 41 (`xl` 10-22) solved the tuned two-port drive on the copper F-human at its degree-2 root (C_in(d2) = 115.5209 pF) with every gate green; the root-to-root order moves (mean `\|B₁⁺\|` −3.17 %, SAR10g −6.30 %) fire the rule's second branch, so **the human-scale tuned figures are quoted at degree 2 only, at that root** (§2.1; H31 carried them into `WF-10`, `mri:17`, `ans:8`). Licensed: the F-small degree-2 exact-rectangle S / `f₀` / Q as h-converged levels (T6, T13, T41) and, on the copper F-human with 0.2 Ω ESR, the in-model Larmor roots at degree 1 and 2, the passive reciprocal degree-2 32×32, the degree-2 `\|B₁⁺\|` / 10 g hotspot, thermal column and quarter route at the root — **internal tuning statements**. Not claimed: any absolute or external tuning figure at human scale (`ANS-8` step 2, W2, is adjudicated when the operator's column lands and does not hold the exit), compliance, C95.3, homogeneity, tissue, F-clinical. One finding travels with the exit: the lumped-sheet excess ΣX_s / P_acc 2.15e-02 at degree 2 (known-issues 2026-10-11, chain H's H33 / H34), so every per-watt human-scale figure states its normaliser. CLAUDE.md's exit sentence is the operator's to update. **2026-10-26 weekly: the travelling finding is a record** — H33's width law (ΣX_s / P_acc falls monotonically with the sheet width at degrees 1 and 2, exponents −1.69 / ≈ −2.09 on the F-small exact-rectangle four-port), known-issues 2026-10-11 retired, §2.2 (4) re-filed; the normaliser rule stands; H34 withdrawn (no Palace field export). **2026-11-04 weekly: the clinical twin's degree-2 root is bracketed** — one − → + crossing on P33 (1, +1j) at 185.763 pF by the chord (entry 53, −4.05 % from the degree-1 root), every gate green at 480 GiB; the 0.25 pF reading is entry 56 (11-05); an internal statement, no tuning claim. **2026-11-07 weekly: the clinical twin's degree-2 root is read at 0.25 pF — 185.7629 pF (entry 56, one crossing on P33 (1, +1j) by the bracketing pair) — and its degree-2 64×64 is a record (entry 54: reciprocal, passive, copper-consistent; `C_mode`(d2) 185.4048 pF; the asymmetric class survives the order move ×1.44, so it is not degree-1 discretisation); internal statements on F-clinical, outside the exit; the degree-2 column `ANS-11` adjudicates on is `xxl` 11-14 (entry 62, READY).** **2026-11-09 weekly: the 3 T twin's degree-1 root is read at 0.05 pF — 44.0525 pF (entry 57, first branch, all four series inside the grid) — and the first 128 MHz clinical-scale `|B₁⁺|` / 10 g figure is a record at it (entry 58: every gate green; 3 T / 1.5 T per accepted watt — lattice `|B₁⁺|` 1.41, SAR10g(c\*) 4.10, of which the box's power share is 3.69); `ANS-12` (the 3 T box-loaded case) commissioned with its node TBD — the 3 T degree-2 root is `xl` 11-15 (entry 64) then a 0.25 pF refinement; internal statements, outside the exit.** **2026-11-14 weekly: the 1.5 T clinical degree-2 column exists (entry 62, every gate green on both calls at 486.57 GiB) and its per-watt order moves are ≤ 1.26 % with the node change included — the order rule does not fire for the box-loaded figures, and `ANS-11`'s `COMPARISON.md` carries the degree-2 column (C37); internal statements on F-clinical, outside the exit.** |
| 7 | Implants: parametric implant geometry in the phantom, local SAR / near-implant hot spots | `IMP-1`, `IMP-2`, `GEO-36`, `IMP-3`, `WF-9` (listed 2026-09-30 weekly; subgoals (7a)–(7d) in §10), `IMP-7` (added 2026-10-24 weekly, subgoal (7e)) | **Open 2026-09-28, in parallel with Phase 6** (operator directive; §9 rule (5) lane 3). **Started 2026-09-29:** `IMP-1` — a PEC sphere in lossy gel under an imposed quasi-static field, dipole coefficient to 0.30 % and shell SAR within 0.8 % of the closed form at three rungs; bands-only closure ratified 2026-09-30 (the residual floors near 0.5 %, a record). Imposed field, one sphere — no coil, no rod, no heating, no published-data comparison. **2026-10-03 weekly: four of five gating chunks ✅** — `IMP-1`, `IMP-2` (spheroid tip enhancement), `GEO-36`, `IMP-3` (the rod in the coil-driven F-small phantom, identities only; the tip ratio is unadjudicated) — plus the lead-shaft canonicals `IMP-4` ✅ / `IMP-5` 🟡; `WF-9` 🟡 on exact adjoint identities, its ✅ gate written (chain I8). Still no heating, no published-data comparison (7c waits on the operator's reference). Chain I in §10; no date. **2026-10-05 weekly: `WF-9` ✅ as a fixture (I8) and `IMP-5` ✅ — all five gating chunks closed; Palace twins at orders 1 / 2 (I9, I12, I13) agree with ours to ≤ 1.5e-8 on near-tip E, and the tip ratio moves −4.961 % with order in both codes — not p-stable at 1 %; I14 / I15 carry it on the exact-rectangle fixture**. **2026-10-07 weekly: I14 — on the 1 g rung η_rod(T1) moves −6.0 % with order while the 1 g ball ratio moves 0.28 %; I15's ladder runs (r/(4√2) heavy, r/8 `xl` 10-14); the TF at degree 2 (I16) is a record; the insulated lead's TF (I18) opened; (7c) still waits on the operator's reference**. **2026-10-10 weekly: I18 landed (the coated rod and the insulated lead's TF, `imp:7` rows); I15's r/8 rung runs `xl` 10-14; the rod moves to human scale — `GEO-36` step 6 (the rod in the F-human phantom, I21a) and `IMP-3` step 8 (the rod under the tuned copper F-human's two-port drive, I21b), identities and a printed tip ratio, no published comparand until (7c)**. **2026-10-12 weekly: I21a / I21b landed — the rod under the tuned copper F-human's two-port drive, every identity green, tip 1 g ratio 3.80 printed (1.38× the `IMP-2` estimate); `IMP-6` ✅ (the PEC oblate spheroid — a disc edge — rim enhancement 1/L_x to 0.08 / 0.63 % at aspect 2 / 5, opened by the daily under lane 3 and ratified), `imp:8`; the r/8 tip rung runs `xl` 10-14; the disc moves into the phantom generator and the coil-driven phantom (I24 / I25); (7c) still waits on the operator's reference**. **2026-10-14 weekly: the disc is in the coil-driven F-small phantom (I24 / I25, rim 1 g Q 1.592, 0.962 of the closed-form estimate); entry 33's cap-ball ladder reversed at r/8 (η_rod(T1) 3.482 → 3.556 → 3.530) — the pointwise tip ratio is not resolved by the cap ball alone (known-issues 2026-10-14; I17 / I30); the disc rim floor is the comparand's finite-k error with a ≈ 0.9 % k → 0 residual that rises with rim refinement (I27 / I28; I29, then the family freezes); new human-scale implant work moves to F-clinical (chain C)**. **2026-10-17 weekly: I17 landed (at 4 mm the near-T1 cells sit inside the ball; I30 is over the rod cell guard → `xl` 10-23); I29 froze the disc family (body-corrected miss +4.1–4.4 %, E3-not); I31's composed disc-in-F-clinical build died in gmsh at the heavy ceiling — re-posed as a daytime window with a gmsh stop rule, then `xl`; (7c) still waits on the operator's reference**. **2026-10-19 weekly: the composed clinical builds fit a daytime window (the disc 441 s, the rod 419 s of gmsh at `-n 4`, both under the 450 s rule) — `IMP-6` step 4 and `GEO-36` step 8 landed: the aspect-5 disc and the F2182-style rod each sit in F-clinical's ASTM box with identities green, no solve; the rod's first solve is a priced `xl` cost probe (I34, entry 45, after entry 40 prices the fixture); I30 runs `xl` 10-23; (7c) still waits on the operator's reference**. **2026-10-21 weekly: I35 landed (the rod-in-clinical-box column module builds in 457 s at loadavg 0.36 — the 10-20 kill was load, not geometry) and I34, the rod's first priced clinical solve, is commanded for `xl` 10-25; `IMP-7` opened as subgoal (7e) — implant-driven gel heating in closed form (a PEC sphere's dipole SAR heating a Robin-cooled shell, a Carslaw & Jaeger derivation on the `MAT-21` precedent), FEM step after; (7c) still waits on the operator's reference**. **2026-10-26 weekly: `IMP-7` ✅ with P2 records on curved cells (I42 / I43 — the affine P2 rate ≈ 2 was the polyhedral-boundary floor; the curved rate 3.04–3.23 on both boundary types; the rate gate is I44, the transient I45 / I46); the rod's own SAR heats the F-small phantom on a matched mesh pair (I41, tip 1 g θ ratio 1.10, a record); `xl` entry 45 priced the rod-in-ASTM-box clinical column (773 s, 41.76 GiB at `-n 8`) and I36 — the rod under the clinical two-port drive — is commissioned for `xl` 10-30 once its body (I36b) lands; (7c) still waits on the operator's reference**. **2026-10-31 weekly (the 10-28 session was lost to the host outage): implant heating has its transient closed form and FEM — the hollow-sphere series (I45), the CN transient on the curved P2 hole within 1.8e-4 / 7.7e-5 of it (I46′), its duty-cycled form (I47) and dose (I48a / I48b, `imp:18`), the curved P2 rate gated ≥ 2.6 (I44); I36's body and its stub proofs landed (I36b–I36d) and entry 48 runs Monday 11-02 at the ratified 193.5996 pF node; I32a (the clinical implant-heating body) and entry 55 enumerated; (7c) still waits on the operator's reference**. **2026-11-02 weekly: entry 48 read — the rod in the clinical ASTM box under the two-port drive is a record (every gate green; tip 1 g Q 4.34–4.46, η_rod(T1) 3.642, the rod's power share 3.7e-5 at degree 1; box containment not asserted on this route), I32 queued at its price (entry 55, Friday 11-06); the pole's dose gap (I48c–I54) closed as a record — not h, Δt, their pairing or truncation; `IMP-8` ratified (the finite line source and its two-tip twin in closed form, `imp:21`) with its FEM half posed (I55 / I56) and the like-for-like tip-Q reading on F-small (I57); (7c) still waits on the operator's reference**. **2026-11-04 weekly: `IMP-8`'s FEM half landed — the line source, the two-tip twin and each under a duty cycle, (I) ≤ 4.4e-3 against 1e-2 on an infinite medium (I55–I63, `imp:22`); I57 read the like-for-like tip Q (F-small 3.90–3.97 against the clinical box's 4.34–4.46, both kept-cell rods) and a ≈ 40 % rod-model lever on η_rod(T1) (kept cells 2.158 against the hole's 3.48–3.56 — I67 reads it on one definition); entry 55 (the rod's heating in the clinical box) runs Friday 11-06; the rod's own thermal properties enter with I64–I66 (the composite cylinder, Carslaw & Jaeger); (7c) still waits on the operator's reference**. **2026-11-07 weekly: the rod's own thermal properties are in — the composite cylinder in closed form and FEM with the source in the shell and in the core (I64, I65, I70, I71; `imp:24`, `imp:25`), the gel / Ti / hole bracket on I55's box (I66, I69, I73: Ti +55.5 %, hole −8.5 % at the near-tip point); the first clinical implant-heating window ran (entry 55) with every EM gate green and the tip calorimetric gate red on its comparand's form — ruled option (b) (the consistent-mass projection at the unmoved band) and re-run as entry 63 (Friday 11-13); the rod-model lever is ≈ 3 %, not ≈ 40 % (I67 — the rest was the mesh); the F2182-shaped readout (probe 1 mm beyond the tip, 900 s) is enumerated (I75, I76, I78); (7c) still waits on the operator's reference**. **2026-11-09 weekly: the F2182-shaped readout is built on the bracket fixture in closed form and FEM (I75, I78–I84; `imp:26–28`) — a 1 mm placement step moves the line-source θ by 29 %, and 1 mm beyond a 1 mm-radius rod the zero-radius line is not the like-for-like form (I83's uniform-cylinder Φ_a; the gel arm agrees with it at `INTERIOR_BAND`, ratified); the re-posed tip gate is green on the F-small rod four-port (U127′) and its stub (I32b), so entry 63 (Friday 11-13) buys the clinical figure under it; the 3 T twin's body (I85) and window (entry 67) are enumerated; (7c) still waits on the operator's reference**. **2026-11-11 weekly: the 3 T heating body (I85) and entry 67's reader (I88) are stub-proved, so the 3 T implant window waits on windows only (61 Thursday, 63 Friday); in closed form the tip weighting (× 6.56) dominates the rod's radius (× 0.956) at the probe 1 mm beyond the tip (I86); I87 (the rod's own thermal properties as a knob) re-posed with its anchor — the default arm is the adiabatic exclusion, the ρc-weighted tip form rules on balls holding rod cells; entry 63 (Friday 11-13) still buys the first clinical implant-heating figure; (7c) still waits on the operator's reference**. **2026-11-14 weekly: the 3 T rod-in-box record is read (entry 61: tip 1 g Q 4.75–4.88, η_rod(T1) 3.41, rod share 7.4e-05; Q ≈ 1.08–1.10× the 1.5 T record's, C40); the first clinical implant-heating window (entry 63) ran red — every EM gate and θ record green, two one-step thermal readings above band, one of them a gate I32b folded into a shared maximum below the route's floor — so the route is repaired by I32c (a slot step) and the figure is bought by entry 71 (11-17 at the earliest), the 3 T twin by entry 67 after it; the rod's degree-2 reading (entry 62, call 2) is C56's; (7c) still waits on the operator's reference** |
| 8 | Thermal: Pennes bioheat driven by SAR | `MAT-9` … `MAT-14`, `MAT-5`, `MAT-15`, `ANS-9` (listed 2026-09-30 weekly; subgoals (8a)–(8e) in §10; `MAT-17`, `MAT-18`, `MAT-19` added under (8a) 2026-10-03, not gating) | **Open 2026-09-28, in parallel with Phase 6** (operator directive; §9 rule (5) lane 3). **Subgoal (8a) done 2026-09-29 / 30:** conduction, FEM SAR as the source, the Pennes sink steady and transient, the SAR-driven transient and σ(T) fed back, each against a closed form (`MAT-9` … `MAT-14`, `MAT-5` ✅) — **an imposed uniform field in a small sphere: no coil, no phantom-in-coil, no tissue, no external check.** Next: coil-driven heating of the gel phantom (`MAT-15`), then `ANS-9` (HFSS + Icepak) on `ANS-2`'s fixture. **2026-10-03 weekly: (8b) done — `MAT-15` ✅ 2026-10-02** (coil-driven heating of the F-small gel phantom at 10 MHz: energy identity, calorimetric slope, C4, and the phantom's own cylinder against the Robin eigen-series; σ = 800 coil); (8a) extended by `MAT-17` ✅ (duty cycle) and `MAT-18` ✅ (container wall); **`ANS-9` commissioned — SPEC written, our column is chain P4; no external check has landed.** Chain P in §10; no date. **2026-10-05 weekly: our `ANS-9` column landed (P4) and the case is ready for the operator; (8d) done at F-small (`MAT-5` step 3 on the coil-driven phantom, contraction 1.83e-4); the first thermal figure on a tuned coil (`MAT-15` step 3, 64 MHz copper, identities only); `MAT-19` ✅ (the virtual probe) with its series on the coil-driven phantom — probe volume moves a 10 % readout by ≤ 1.9 %**. **2026-10-07 weekly: P11 / P12 landed as records (the container wall raises ⟨θ⟩(t\*) by 4.63 %; the 50 % duty cycle reaches 0.986 of D·P1, Duhamel at 3.4e-14); P10 red on the 537-cell tuned phantom (a 5 mm probe reads the mesh), re-posed on a refined phantom; `MAT-20` opened — thermal dose (CEM43) and the six-minute window against closed forms; `ANS-9` still waits on the operator (W1)**. **2026-10-10 weekly: (8c) adjudicated — `ANS-9` step 1 AGREES on the EM → thermal hand-off (⟨θ⟩ within 0.8 %, slope 0.5 %, surface heat 1.4 % of AEDT HFSS + IcepakFEA) and DISAGREES at the far-side probe D (ours 11 % below, band 10 %; known-issues, step 1b = P21); (8e)'s first step P17 ratified — the tuned copper F-human heats its own cylindrical phantom, identities green; the ASTM-size phantom does not fit F-human's bore (operator, dashboard W5); `MAT-20` steps 1–2 landed, the FEM clause is P19**. **2026-10-12 weekly: `MAT-20` ✅ (dose and the six-minute window on closed forms and on the tuned F-human's own heating series, `mat:18`); `ANS-9`'s probe-D DISAGREE re-read as discretisation on our side (phantom resolution 13 % + element order 8 %; the degree-2 refined candidate within 2.5 % of AED) — the column is re-recorded at degree 2 on exact rectangles with the 0.0025 m phantom (P24) and re-adjudicated 10-14; `MAT-21` opened under lane 3 (the F2182-sized gel box's triple Robin series, step 0 green; step 1's precondition bar ruled) as (8e)'s closed-form thermal half; the first thermal figure at degree 2 is `xl` entry 38 (P22, pending the 10-13 price)**. **2026-10-14 weekly: `ANS-9` ✅ — every row AGREES with AEDT HFSS + IcepakFEA on the degree-2 exact-rectangle column (probe D within 2.5 %, integral rows within 1.3 %): Phase 8's external check is complete at F-small, one-way, gel only; `MAT-21` ✅ with degree-2 and duty-cycle records (P28–P30); (8e) is C6's thermal half on F-clinical (P33); entry 38 (P22, degree 2) runs `xl` 10-19; P27 (two-way) not posed until the operator asks**. **2026-10-17 weekly: P32 / P34–P36 landed — the gel box under the separable, DG0, P2-interpolated and transverse-quadratic (CP eddy-SAR) sources, records against the triple Robin series at degree 2; entry 38 is superseded by entry 42 — the degree-2 thermal column is posed at the degree-2 root (P22b first) so it does not heat a detuned coil; P37–P39 (probe points, duty cycle, dose on the box) are the closed-form supply**. **2026-10-19 weekly: entry 42 — the tuned copper F-human's heating at the degree-2 root moves −1.27 % per accepted watt from the degree-1 root with every identity green, so P17's thermal figures are quotable at degree 1 with the move beside them (first branch); its dose at degree 2 is P40; P37 / P38 landed (the F2182 box's probe points and duty cycle, records against the triple Robin series); P39's anchor is written (the dose's sensitivity bound) and P41 (the duty-cycled dose) enumerated**. **2026-10-21 weekly: P39′–P44 — the F2182 box's dose at all six probe points under the continuous, D = 0.5 and D = 0.25 sources, the (Bd) sensitivity bound asserted on every cell of the matrix with its headroom falling as the duty halves (5.04 → 3.05 → 1.72× at P5; P45 checks Δt); P40 — the dose at the degree-2 root moves −1.3 % of the excess over baseline beside θ's −1.27 %; no new phase claim**. **2026-10-26 weekly: P46 closed as a measured negative (the CN Δt-rate of the F2182 box's dose against the series is floor-bound by the P2 spatial error; CN is second order against itself — `mat:27`'s record); no step landed this interval; P47 (a non-uniform initial temperature, the pre-equilibration drift) and P48 (the acrylic wall as an effective Robin coefficient, then meshed) enumerated as closed-form + FEM pairs; (8e) is C6's thermal half (P33, after `xl` entry 52)**. **2026-10-31 weekly: P47–P50a landed as closed-form + FEM pairs — the drift from θ₀ = g·z (P47a / P47b′, `mat:29`), the acrylic wall as h_eff and meshed (P48a / P48b, `mat:28`), the walled box's transient (P48a″ / P48c, `mat:30`), duty-cycled (P49a′ / P49b, `mat:31`) and its dose (P50a); records at `MAT-21` step 1's bands, no new phase claim; P50b–P54 enumerated; (8e) is still C6's thermal half (P33, after `xl` entry 52, Tuesday 11-03)**. **2026-11-02 weekly: P50b–P60 landed as closed-form + FEM pairs — the walled box's dose (`mat:32`), the pre-equilibration bias in temperature and in dose (`mat:33`, `mat:34`; a skipped wait costs ≈ 2.5× more in dose than in temperature), the dose gap on and off a mesh plane (`mat:35`), the wait under D = 0.25 (`mat:36`); records, no new phase claim; (8e)'s FEM half is P33 = entry 59, its body P61 a slot step tonight, its window priced by the 11-04 weekly from entry 52 — ≈ 11-09 … 11-12 if 52 is green; P62–P65 (the drift and the wait inside the walled box) enumerated**. **2026-11-04 weekly: P61 (the clinical box's heating body) landed and entry 59 is READY on one knob step (P61b: both senses), priced ≈ 950–1 150 s, Tuesday 11-10 at the earliest — the first dated Phase 8 step at clinical scale; P62–P64 landed as records; P65 found the massless-wall drift comparand short by ≈ 3× its band under D = 0.25 while its massless twin is green — the comparand's wall capacity, whose closed form is P67 (the two-layer composite slab) with P68 / P69 re-reading P65 / P63 against it; no new phase claim**. **2026-11-07 weekly: P67–P72 landed — the composite-slab drift reproduces P63's measured bias at P5 and P68 closes P65 at (I) 5.5e-03 against 1e-2 (the massless comparand's gap was the wall's heat capacity), the dose under the capacitive wall (+3.73 % at P5, D = 0.25) and two `mat:` examples; entry 59 (P33, the clinical box heated by the coil, both senses) is QUEUED for Tuesday 11-10 and entry 60 (the 3 T twin) PENDING on 58; the dose vs the wall's thickness (P73 / P74) and the uniform-source comparand at entry 59's probe points (P75 / P76) enumerated; no new phase claim**. **2026-11-09 weekly: P73–P82 landed — the dose under the capacitive wall across the wait in closed form and FEM (composite / massless − 1 falls +3.73 → +1.01 %; FEM / composite − 1 +0.476 → +0.007 % over t_w = 0 … 2/λ₁), the three closed-form SAR pictures per deposited watt at entry 59's probe points and on C6's lattice (U3, eddy, end-weighted — only the end-weighted one reproduces entry 52's x+0mm shape; `mat:40`, `mat:41`); entry 59 (P33) runs Tuesday 11-10 and entry 60 (the 3 T twin) is QUEUED for 11-11 on entry 58's reading; the 1.5 T `ANS` thermal case is posed by the 11-11 weekly if 59 is green; no new phase claim**. **2026-11-11 weekly: (8e)'s FEM half landed at both field strengths — entries 59 / 60 (every executed assertion green, records): θ(t\*) at c\* per accepted watt 1.61e-03 K at 1.5 T and 6.52e-03 K at 3 T, the long-side wall probes 2–3× hotter, the 3 T / 1.5 T lever 4.04 ccw against the SAR10g lever 4.10 and the closed-form pictures' 3.69 (the pictures' shape part was 1 by construction — epitaph; P89 computes it); `cem43` per watt is the 20 °C baseline, so no dose is per-watt (P88); **`ANS-13` commissioned — the first external thermal check at clinical scale** (HFSS + IcepakFEA on `ANS-11` / `ANS-12`'s projects, per deposited watt, `ANS-9`'s bands), its adjudication columns at degree 2 (`xxl` entries 69 / 70 on P87); records on a 20 °C gel with no container wall, no compliance claim, no external check yet**. **2026-11-14 weekly: `ANS-13`'s runnable half landed (C53 / C54, Done-when (1) met) and its degree-2 body is stub-proved (P87); in closed form a dose is not per-watt (P88) and the 100 mm lattice resolves the lever's ≈ 9 % excess (P89′: predicted 4.078 against the measured 4.0425); the 1.5 T degree-2 thermal column is queued (`xxl` entry 69, 11-28), the 3 T one waits on the 3 T degree-2 root (68 → 66 → 70); records, no compliance claim, no external check yet** |
| 9 | **A usable package for MR-safety research** (redefined 2026-09-30, operator directive): a study-level Python API + YAML/CLI, the external-CAD STEP → solve path, packaging and docs; plus `POST-9` adjoint sensitivities, `ARCH-1`, AMR, MPI scaling | `WF-13` … `WF-17`, `WF-1`, `GEO-37`, `GEO-34` (listed 2026-10-03 weekly; chain U in §10) | **Defined 2026-09-30; opened 2026-10-03 in parallel with Phase 6 for (9a) API, (9c) packaging, (9d) docs — Phase 6 steps take precedence for a slot; `POST-9` / AMR / `ARCH-1` wait for Phase 6's exit** (operator directive 2026-10-03, §10). Previously: **`ARCH-1` (solver-agnostic MRI-safety layer) is embargoed until the 2026-10-17 weekly by operator instruction, §7**; **AMR is feature ladder C4** (operator directive 2026-09-04): not a boundary condition, but the largest remaining HFSS workflow gap (adaptive passes are what make an HFSS answer trustworthy without a mesh study); revisit after Phase 6; **`POST-9` adjoint sensitivities added 2026-09-28 (operator directive) — the capability HFSS does not offer, basis of the worst-case implant-heating research direction**. **2026-10-03 weekly: first steps enumerated as chain U (U1 installable … U7 the tutorial as a test), each reproducing a gated result inside its unmoved band; nothing has run, so no date is given — the operator's target is a usable v0 by the end of November**. **2026-10-05 weekly: seven chunks ✅ in 1.99 days (`WF-13`, `WF-14`, `GEO-37`, `WF-15`, `WF-16`, `WF-1`, `WF-17`) — a `Study` + YAML + installed package reproducing two gated records inside their unmoved bands; U9–U17 enumerated (sweep, temperature, native `dS`, STEP port sheet, STEP F-small, resource estimate, API reference, U8's Palace backend after 10-17, `GEO-35`); the v0 date is bounded by the CAD kernel and the docs-hosting choice (operator), not by slot pace; Phase 9 owes five examples (`EX-91` … `EX-95`)**. **9d progress: the generated API reference landed 2026-10-06 (U15, §9 item 286; `docs/guide/api/`, 34 pages + index).** **2026-10-07 weekly: U11–U20 all landed by 10-07 — the installed `femem run --ranks N` with a YAML `thermal:` section, `Study.estimate()`, the API reference; seven examples (`mri:7` … `mri:13`), ramp met; (9b) is blocked on the operator's CadQuery kernel (Gmsh's own STEP writer keeps no names, known-issues 2026-10-07); U21–U30 enumerated over the gated coil identities (10 g SAR, `\|B₁⁺\|`, tuning, copper, F-human through the `Study`; `WF-8` placed under 9d with a target); v0 still bounded by the operator's two items**. **9a results: 10 g SAR through the API landed 2026-10-08 (U23, §9 item 321; `Study.sar`, symmetry identity on item 192's fixture, no compliance claim).** **9a results: CG1 `\|B₁⁺\|` through the API landed 2026-10-08 (U24, §9 item 322; `Study.b1_plus` / `b1_plus_at`, C4 covariance identity on the same fixture, no absolute or homogeneity claim).** **10 g SAR and B₁⁺ through the Study: example present (`mri:14`, `EX-106` ✅ 2026-10-09, §9 item 340).** **2026-10-10 weekly: 15 U steps in 2.98 days — copper, SAR / B₁⁺ and probes from YAML (`WF-1` 4–6), the `Study` on copper and on F-human (`WF-16` 7–8), figures (`WF-8` ✅), two more tutorials; 14 examples; next: `Study.tune` (U25), capacitor sheets with ESR (U37), the tuned copper F-human through the API (U36), the validation-envelope page (U41), presets (U42); `WF-18` (PrecessEM) placed as an interactive atomic commit after the 10-12 / 10-14 weekly (W6); 9b still on the operator's CadQuery**. **2026-10-12 weekly: 8 U steps in 1.98 days — the tuned copper F-human runs from YAML alone (`f_human_birdcage` + `two_port_quadrature`, copper SAR / `\|B₁⁺\|` outputs) and its tutorial is a test (U43–U45); `CapacitorSheets(C, esr_ohm)` and the F-human `Study` (U36 / U37); the validation-envelope page generated from the registry (U41); `GEO-35` ✅ (the leg-box clean-up, all consumers re-recorded); 108 examples; next: `Study.tune` (U25), presets (U42), thermal from the F-human YAML (U46), the installed-package rehearsal (U47); v0 still bounded by the operator's three items (CadQuery, docs hosting, the rename)**. **2026-10-14 weekly: 8 U steps in 1.98 days — `Study.tune` on PEC and copper at degree 1 / 2 (U25(a), U50, U51), a `tune:` YAML section (U49), thermal from the F-human YAML (U46), the installed-package rehearsal (U47), the envelope page's h-statement entries (U48); 111 examples; U52–U55 enumerated; v0 still bounded by the same three operator items**. **2026-10-17 weekly: 9 U steps in 2.98 days — presets ✅ (`WF-4`), `capacitor_sheets` and thermal-only `thermal_ball` files from YAML (U57 / U58), `--estimate` refusals with rc 2, the tuning page; 121 examples (+10); `ARCH-1`'s embargo lifted — U16, the Palace backend through the `Study`, opens 9e; U59 / U60 / U64 (thermal-box files, `Study.probe` on thermal-only Studies) enumerated; v0 still bounded by the operator's three items**. **2026-10-19 weekly: 10 U steps in 1.98 days — thermal-only `thermal_box` files with sources, duty cycles and probes from YAML (U59–U67, ratified), the Palace backend through the `Study` and the CLI (U16 ✅, `ARCH-1` 🟡, `mri:30`), the backends guide page; 127 examples (+6); U68–U70 enumerated (the F-human estimator branch, the thermal-box tutorial, probes on EM files); v0 still bounded by the operator's three items**. **2026-10-21 weekly: 6 U steps in 1.99 days — the F-human estimator branch at degree 1 (U68), the thermal-box tutorial and the four-port outputs page as tests (U69 / U73), point probes on EM files, PEC and copper (U70 / U72), the multi-frequency sweep from YAML to one `.s4p` (U71); 131 examples (+4); U74–U78 enumerated (N-slot `ring_mode_rows` / `Study.tune`, the sweep + conductor and Palace-output refusals, the degree-2 estimator row, the sweep page); v0 still bounded by the operator's three items**. **2026-10-26 weekly: 8 U steps in 1.98 days — the ESR, tune, four-port-outputs, EM-probe and human-scale-estimate pages as tests (U85–U88, U90a), copper on the Palace backend (U91, (E2-Cu) 0), a Cole–Cole medium through the `Study` gated against the lossy-sphere closed form (U92, 0.0421 %) and its page (U93); 148 examples (+6); U94 (a public saved-result loader), U90b, U95 (`ColeColeMaterial` public with its envelope), U96 (F-clinical through `Study.estimate`), U97 / U98 (the rod in the phantom through the `Study` and its page) enumerated; v0 still bounded by the operator's three items**. **2026-10-31 weekly: 12 U steps in 4.99 days — a saved S-matrix is tuned with no field solve through the `Study` and the YAML on the F-small four-port, the copper 32-port F-human and the F-clinical 64×64 (U94, U99, U100′, U104, U105; `mri:35`, `mri:36`), `ColeColeMaterial` public (U95), F-clinical through `Study.estimate` (U96), two more pages as tests; 156 examples (+8); U97 re-posed on the `Study`'s own fixture (U97′), U98 / U102 / U106–U108 enumerated; v0 still bounded by the operator's three items**. **2026-11-02 weekly: 12 U steps in 1.99 days — the implant rod through the `Study` and YAML at one frequency and on a sweep (U97′, U109, U110; `mri:37`, `mri:38`), the saved multi-frequency file tuned at a chosen frequency from the API and YAML (U108, U111; `mri:39`), the clinical degree-2 estimate (U106, U113), four pages as tests (U98, U107, U112′, U114); `OPS-76` fixed a landed refusal the same day; 173 examples (+17); U115–U122 enumerated (a public field route on a sweep, the 3 T twin's saved tune and estimate, the rod fixture's estimate, the installed-package rehearsal re-run, the rod with capacitor sheets as U110's gate); v0 still bounded by the operator's three items**. **2026-11-04 weekly: 9 U steps in 1.99 days — `solve_sweep(…, fields_at=)` (the public field route on a sweep), the 3 T twin's saved tune and estimate, the rod fixture's estimate, the installed-package rehearsal re-run, the rod with ESR capacitor sheets through the `Study` (on the PEC fixture, ratified) and its page, `mri:33` / `mri:38` on the public route; 175 examples (+2); U124–U130 enumerated (the measured degree-2 clinical estimate, the thermal-box initial profile and wall from YAML, the rod's heating through the `Study`, the estimate table as a page); v0 still bounded by the operator's three items**. **2026-11-07 weekly: 8 U steps in 2.99 days — the thermal box's initial profile and wall from YAML (U125 / U126), the estimate table pinned on every shipped file (U130), `MAT-12`'s SAR-driven transient through `Study.temperature` on the EM sphere at degrees 1 and 2 and on a page that is a test (U131–U133, `mri:40` — the one narrow new workflow of the interval), the degree-2 clinical `C_mode` through `load_result` + `tune` (U129); the rod's heating through the `Study` (U127) is re-posed on tonight's tip-gate comparand (U127′); 182 examples (+7); U134–U136 (the YAML twin of the sphere transient, its duty cycle, its page), U139 / U140 (the dose readout from YAML) enumerated; v0 still bounded by the operator's three items**. **2026-11-09 weekly: 13 U steps in 1.98 days — the implant rod's heating through the `Study` and the rod four-port `thermal:` YAML (U127′ / U128, the re-posed tip gate green), the sphere transient's YAML twin, its duty cycle at degrees 1 and 2 and the `thermal.dose` readout from YAML (U134–U136, U139 / U140, U143–U145), the thermal routes re-run green after six merges (U141 / U142 / U146); 187 examples (+5); U147 (the tuned two-port `duty` refusal) and U148–U152 (the guide pages' solve nodes re-run on the current tree, the implant page's YAML heating section) enumerated; v0 still bounded by the operator's three items**. **2026-11-11 weekly: 8 U steps in 1.98 days — the tuned two-port `duty` refusal (U147) and the guide pages' solve nodes re-run on the current tree in six batches with every record at its prior (U148–U151, U154, U155; the record / reading and "constant is the record" rulings now standing; the solves are not run-to-run reproducible at the eleventh digit, 2e-12 at degree 1); 190 examples (+3); U152 re-posed on the identity the route can hold (θ(t\*) at `IDENTITY_BAND`, the slope rows as readings); U156–U158 enumerated (two page sections as tests, the determinism measurement); v0 still bounded by the operator's three items**. **2026-11-14 weekly: 9 U steps in 2.99 days — the CadQuery STEP import with a named port sheet (U12′: 9b's kernel is in), the implant page's one-step section (U152′), the determinism reading (U158: the degree-2 solve differs run to run at 1e-10 with the figure flag unset), the clinical page's rod-in-box and closed-form sections, the estimates page's measured 3 T and 1.5 T sections; 192 examples (+2); two of the three operator items settled on 11-14 (CadQuery installed; docs on GitHub Pages — U164 / U165 pose the build-as-a-test and the Pages config); U13 re-posed on identities the import route can hold (U13a / U13b); v0 waits on the rename (W6) and the manual push** |

Phases 6–8 are the 2026-08-04 scope adjustment (§1). Their phase goals and
subgoals live in §10's long-horizon roadmap and are **owned by the weekly
planning review**; the daily review breaks current-phase subgoals into
implementer-sized items and does not restructure phases.

### Critical path

```
TH-1 (real complex time-harmonic formulation)
   ├─> TH-6/7/8 (analytic validation gates)
   ├─> MAT-2 ──> MAT-6 (materials actually affect fields; Dodd–Deeds)
   └─> PORT-1 (real port excitation) ──> PORT-2…8 ──> WF-5…8
```

**`TH-1` landed 2026-07-31 and `PORT-1` closed 2026-08-15; the constraint
moves down the chain.** The `⚠️` backlog still may not be extended until
revalidated against the real solve, and nothing S-parameter-shaped beyond
the two-torus fixture grows except through the scoped birdcage-port lineage
(`PORT-9` at 10 MHz → `PORT-11` at 64/128 MHz, 2026-08-23).

---

## 7. Chunk backlog

IDs are prefixed by subsystem and **globally unique and stable**. Legacy
`A1`/`B2`/`C3` IDs in commit messages and old logs map here via §8.
**Row budget (operator directive 2026-10-01): every row ≤ 2 KB** — ID, short
title, status glyph, a ≤ 2-sentence current state, the **Done-when** and the
tier; the step history lives in `docs/planning/chunks/<ID>.md` (each row's
*History* pointer). The 2026-10-01 migration moved 72 rows' text verbatim
(`rotate_plan_archive.py chunks`, spec
`docs/planning/s7-row-migration-2026-10-01.spec`, log
`20261002T041822Z_S7-ROW-MIGRATION.log`); the weekly keeps the budget
(weekly-review.md step 6) and hand-condenses rows the migration condensed
mechanically.
These tables are the authoritative *status*. Closed chunks and steps keep
compressed result blocks here (gated numbers, dates, log IDs, live
carry-forwards); their full plans, execution journals and audit narratives are
archived verbatim in `docs/planning/plan-archive.md` — grep there before
re-deriving a closed step's diagnosis. (The older per-chunk log,
`docs/testing/pending-tests.md`, was removed 2026-08-04; see git history.)

### OPS — Infrastructure & testing operations

| ID | Title | Status | Tier |
|---|---|---|---|
| `OPS-35` | **Repository entry-point and supported-environment consolidation** — interactive maintenance pass requested by the operator after the 2026-09-03 repository-health review. Replaced the stale README (including its contradictory phase table and nonexistent `CircularLoop` API example) with a concise alpha-status entry point linked to the live dashboard; added `CONTRIBUTING.md` with the repository's Docker, evidence, MPI, and private-data rules; pointed MkDocs at the live dashboard while retaining the 2026-07-31 page as an explicitly historical snapshot; aligned all CI solver jobs and package metadata with the already-supported DolfinX 0.11 image; updated `setup-python` v4 → v5 and Codecov v3 → v4; and installed the declared documentation toolchain in the development image. No solver source, physics assertion, tolerance, record, or capability claim changed. Verification: `20260903T183816Z_OPS-CLEANUP-SMOKE.log`, **16 passed**, Status 0, 1 s; `run_examples.sh --list` confirmed every README selector; TOML and both YAML files parse; shell syntax and `git diff --check` pass. The preceding combined attempt, `20260903T183801Z_OPS-CLEANUP.log`, records the smoke suite green followed by Status 127 because the running pre-change container lacks `mkdocs`; the Dockerfile fix addresses that missing dependency, but the strict documentation build remains for CI or the next image rebuild. | ✅ *(documentation/operations change; smoke assertions unchanged and green)* | smoke |
| `OPS-36` | Log retention policy, harness chatter filter, weekly housekeeping sweep | ✅ Closed as documentation and operations: the log retention policy, harness chatter filter and weekly housekeeping sweep landed, with gating logs exempt from the volume budget by operator decision. *History: `docs/planning/chunks/OPS-36.md`.* | smoke |
| `OPS-37` | Re-base the `PORT-1` step-4 gate module's two mutual-ratio records onto the 0.11 image and tighten its reproduction control to 1e-6 | ✅ Closed: the `PORT-1` step-4 gate module's two mutual-ratio records are re-based onto the 0.11 image under a 1e-6 reproduction band, with every physics band unmoved. *History: `docs/planning/chunks/OPS-37.md`.* | standard (179 s + 173 s measured) |
| `OPS-38` | `write_xdmf_with_tags` learns facet tags, and the three improvising examples move onto it | ✅ Closed: `write_xdmf_with_tags` writes facet tags under a closed-form round-trip gate, and the three improvising examples now call the helper. *History: `docs/planning/chunks/OPS-38.md`.* | smoke + standard (host-runner; measured 6 / 65 / 103 / 113 / 2 s) |
| `OPS-39` | The `TH-15` step-0 probe's `_census` made rank-safe, gated on the count identity at `-n 2` | ✅ Closed: the `TH-15` step-0 probe's census is rank-safe and gated on the count identity at `-n 2`, and its known-issues entry is retired. *History: `docs/planning/chunks/OPS-39.md`.* | smoke + standard |
| `OPS-40` | `evaluate_vector_field_parallel` refuses a point list that differs across ranks | ✅ Closed: `evaluate_vector_field_parallel` refuses a point list that differs across ranks, gated by executed asserts in a new test module. *History: `docs/planning/chunks/OPS-40.md`.* | smoke + standard (measured 4 / 2 / 57 / 47 s) |
| `OPS-41` | Declare `test_port_package_sparameters.py`'s three digit-reproduction records `-n 2` records, and attribute the 1e-4 width sensitivity to `V` or `I` | ✅ Closed and audited PASS 2026-09-08 10:30: the three digit-reproduction records are declared `-n 2` records, and the width sensitivity is pinned as a strict xfail under its override variable. *History: `docs/planning/chunks/OPS-41.md`.* | standard × 3 (181 + 136 + 142 s) |
| `OPS-42` | The corpus census's staleness window reports the run cadence, not staleness | ✅ Closed 2026-09-09: the corpus census's staleness window no longer reports the run cadence, and the exit-code contract is kept. One unrelated red on `main` was filed to known-issues, not edited. *History: `docs/planning/chunks/OPS-42.md`.* | smoke |
| `OPS-43` | Long-window robustness: durable capture, orphan-rank cleanup, and per-run memory instrumentation | ✅ Closed at the 2026-09-11 03:00 daily review: durable capture, orphan-rank refusal, per-run memory instrumentation and the solver-progress inertness gate all landed. The tier was corrected to standard, not the evidence. *History: `docs/planning/chunks/OPS-43.md`.* | standard *(re-declared 2026-09-11 03:00 review from *smoke*: measured 6–61 s under 120–180 s wrappers)* |
| `OPS-44` | Re-pin `COMMITTED_EXAMPLE_ARTIFACTS` to the five artifacts git actually tracks | ✅ Closed 2026-09-09: `COMMITTED_EXAMPLE_ARTIFACTS` is re-pinned to the five artifacts git tracks, and the exemption still means tracked by git rather than anything under `examples/`. *History: `docs/planning/chunks/OPS-44.md`.* | smoke |
| `OPS-45` | The harness footer must not call a red durable-capture window green | ✅ Closed and audited PASS at the 2026-09-11 18:00 review: the harness footer no longer calls a red durable-capture window green. *History: `docs/planning/chunks/OPS-45.md`.* | smoke |
| `OPS-46` | Move the heavy §7 chunk histories out of `PROJECT_PLAN.md` | ✅ Closed 2026-09-12 21:00 slot: all 73 listed rows are moved byte for byte with anchors (i), (ii) and (iv) green (plan 806 962 B, under 850 000 B); (iii) landed by the operator (`2a0ca4d`) and exercised by the 2026-09-13 10:30 review. **Audited PASS 2026-09-13 10:30** (log lines re-cited in the history file; two caveats: (iv) was a printed `wc -c`, not a scripted assert — an `OPS-47` step 1 rider — and "74 rows" was the census discrepancy, 73 were listed and moved). The narrative blocks are `OPS-47`. *History: `docs/planning/chunks/OPS-46.md`.* | smoke |
| `OPS-47` | Move the open chunks' blockquote narratives out of `PROJECT_PLAN.md` | ✅ Closed 2026-09-13: the four open-chunk narratives (`POST-6`, `PORT-14`, `TH-15`, `WF-6`) were moved byte for byte to docs/planning/chunks/ with all anchors green, shrinking the plan to 5 546 lines. The 4 000-line guide is not met, a miss recorded with the guide not widened. *History: `docs/planning/chunks/OPS-47.md`.* | smoke |
| `OPS-48` | Restore the `EX-14` VTX read-back gate on the 0.11 image | ✅ Closed 2026-09-16: the `EX-14` VTX reader was ported to the adios2 2.12 bindings and both examples now execute the read-back gate at relative difference 0 with the wrong-comparison control at 1. Audited PASS 2026-09-18; the writer-side hole that survey found went to `OPS-50`. **Done-when:** the reader ported to the 2.12 bindings (`adios2.bindings.ADIOS` or `FileReader`, the VTX local-block walk kept), both examples run through the harness at `-n 2` with `relative difference ≤ VTX_ROUNDTRIP_RTOL` (1e-10, unchanged) **asserted and executed** (the log shows the `✓ written .bp reproduces` line, not the `⚠ unavailable` line), a wrong-comparison control (read-back vs 0.5 × in-memory) printed at rel ≈ 1, the guides' read-back sections re-cited to the new logs, censuses `broken=0` / docrefs `exit != 1`, elapsed recorded, the known-issues entry retired in the same commit. *History: `docs/planning/chunks/OPS-48.md`.* | standard (`mag:1` 7 s, `mag:2` 136 s at `-n 2`, measured) |
| `OPS-49` | A time-series XDMF writer beside `write_xdmf_with_tags` | ✅ Closed 2026-09-16: `write_xdmf_time_series` landed beside the byte-identical `write_xdmf_with_tags`, with the one-collection count identity, an exact h5py round trip and the collapse control asserted at -n 2. Audited PASS 2026-09-18; `EX-56` cannot adopt it and no follow-up is queued. **Done-when:** additive `write_xdmf_time_series(filename, mesh, cell_tags, steps, comm)` with `steps = [(t, {name: Function}), …]`, one temporal collection whose `n` children each carry every field and CellTags; unit test at `-n 2`: (a) the XDMF parses to exactly one collection with `n` `<Time>` values equal to the `t`s written and every attribute on every child (in = out), (b) each step's arrays read back through `h5py` equal the gathered function arrays at rel ≤ 1e-12, (c) negative control — the same steps through the existing `write_xdmf_with_tags` at the pre-change commit (pinned `<sha>^`) collapse to ≤ 1 time value; `write_xdmf_with_tags` / `consolidate_xdmf_grids` byte-identical; elapsed recorded; *History: `docs/planning/chunks/OPS-49.md`.* | smoke (2 s at `-n 2`, measured) |
| `OPS-50` | Close the writer-side hole in the `EX-14` VTX read-back gate, and port the last dead-API `adios2` reader | ✅ Closed 2026-09-21: (a) a B writer failure now raises in both examples with the exit-status control asserted, and (b) the post4 step-5 probe runs on adios2.bindings and prints PROBE_RESULT PASS after the licensed 0.11-image re-record and the two-sided `STEP4_SEP` ruling. The B-median nondeterminism known-issues entry stays open. **Done-when:** (a) a `B` writer failure **raises** in both examples (the `A` writer may stay tolerant — no gate reads it); asserted control: a harness window that runs `mag:1` with `dolfinx.io.VTXWriter` monkeypatched to raise exits **non-zero**, and the same wrapper on the example file taken from the pinned pre-change commit exits **0** (the defect); the unpatched `mag:1` and `mag:2` reproduce `OPS-48`'s printed digits (`relative difference = 0.000e+00`, control `1.000e+00`) at `-n 2`; `PARAVIEW_GUIDE.md:166–170` reworded (both failure paths now exit non-zero). *History: `docs/planning/chunks/OPS-50.md`.* | standard (`mag:1` 8 s, `mag:2` 141 s, controls 6 / 7 s, probe 5 s, censuses 1 s — all at `-n 2`, measured 2026-09-18; step 3b probe 4 s + control 4 s at `-n 2`, 2026-09-21) |
| `OPS-51` | Log labels that state the run's real rank width and frequency | ✅ Closed 2026-09-18: the six hardcoded rank-width and frequency log labels now format the real communicator size and frequency, with a unit guard counting 6 to 0 against the pre-change commit and the step-3a flag-off ladder control green at -n 8. No band, record or assertion was touched. **Done-when:** the six literals replaced by the communicator's size and the rung's actual frequency; a `tests/unit` guard (real build, `-n 2`) asserting the count of the two literal patterns in `tests/validation/*.py` is **0**, with the asserted control that the same count on the files at the pinned pre-change commit is **6** (count identity, in = out); the ladder's step-3a flag-off control re-run at `-n 8` (its record width — `20260914T093325Z_ANS-4-step3a-w1.log`'s command verbatim) green with the 2a″ record assert **executed** at `STEP3A_RECORD_RTOL` (1e-6, unchanged) and the rung header reading `-n 8`; *History: `docs/planning/chunks/OPS-51.md`.* | standard (guard 2 s at `-n 2`, collect-only 3 s, ladder 68 s at `-n 8`, lumped_column 36 s at `-n 2`, measured 2026-09-18) |
| `OPS-52` | The F-human cost probe's labels say what ran | ✅ Closed 2026-09-19: the WF-7 cost probe's phase-4 label is now a pure function that says flag-off control only at degree 1, the bracket dicts are re-dated from measured points, and the -n 8 flag-off control window passed with the S_driven record asserted. **Done-when:** the phase-4 label is a pure function that says `flag-off control` only at degree 1 and otherwise prints the relative move beside the degree-1 record, asserted by a `tests/unit` test (including 5.50 % on the 09-19 values) with a source-text count control 1 → 0 against the pinned pre-change sha; the two non-gating bracket dicts re-dated from the measured points with their logs; the `-n 8` flag-off control window green with the probe's own executed `S_driven` assert against an untouched `S_DRIVEN_STEP0_RECORD`; the 2026-09-19 known-issues entry retired in the same commit. *History: `docs/planning/chunks/OPS-52.md`.* | heavy by record (one 154 s control window at `-n 8`, measured 2026-09-19) |
| `OPS-53` | `write_setup_figure` refuses a title it cannot draw | ✅ Closed 2026-09-19: `write_setup_figure` now raises a ValueError on every rank for titles over `MAX_TITLE_CHARS` = 130, with the unit test, call-site scan and pre-change control green and the census unchanged. The clipped-legends half of the known-issues entry stays open for other items. **Done-when:** a `MAX_TITLE_CHARS` constant set by measurement (bracket [126 whole, 150 clipped] ⇒ 130) and a `ValueError` raised before plotting on every rank; a `tests/unit` test asserting 131 raises / 130 does not, an `ast` scan asserting every constant `title=` at the 18 `write_setup_figure(` call sites under `examples/` is inside the limit, and the absent-guard control against the pinned pre-change sha; setup-figure census unchanged (`ok=18 missing=36 broken=0`) and one flagged `mesh:5` re-run digit-identical; the 2026-09-19 known-issues entry retired in the same commit. *History: `docs/planning/chunks/OPS-53.md`.* | standard (unit 1–2 s; census 1 s; one 12 s `mesh:5` window at `-n 2`, measured 2026-09-19) |
| `OPS-54` | **XL / XXL windows run against an immutable checkout** — opened 2026-09-19 (interactive session, from an external code review). The three services bind-mount the live working tree (`docker/docker-compose.yml`), and a window holds its own lock, so a slot that commits to `src/` at 05:00 changes `/workspace/src` underneath ranks launched at 02:00. The commit-side half is fixed (2026-09-19: `xl-run.sh` stages and commits only the window's own paths and refuses a HEAD that has left `main`; the harness keeps the ledger row and the in-flight log under gitignored `/logs/` until the window ends — `test_xl_window_isolation.sh` 23/23, `20260920T004915Z_OPS-XL-ISOLATION.log`). The source-side half is open and is **not yet observed to have bitten**: every window from 09-16 on committed before 03:00, and a window is one Python invocation whose imports resolve at start, so the exposure is late imports and data files. **Sketch:** `xl-run.sh` adds a `git worktree` at the 02:00 HEAD under a gitignored path, the xl / xxl services mount *that* at `/workspace`, results are copied out and the worktree removed; touches compose, `bash_guard.py`'s path assumptions, the durable-capture raw-file paths and `--capture-orphan`. **Done-when (§4):** a scratch-repo test in which a concurrent commit to a tracked source file mid-window leaves the window's view of that file byte-identical to its 02:00 content, plus one real XL window through the new path with its record committed as today. For the weekly to schedule or decline; do not queue until a measured window has actually overlapped a slot, or the weekly decides the risk warrants it. | ⬜ | smoke (scratch-repo test) + one ordinary XL window |
| `OPS-55` | `evaluate_vector_field_parallel` refuses equal-count point lists whose coordinates differ across ranks | ✅ Closed 2026-09-21: `evaluate_vector_field_parallel` now folds a point digest into its allgather and raises on every rank when equal-count point lists differ, with the pre-change control showing the silent defect and the caller subset unchanged. No current caller was passing rank-local points. **Done-when:** in `tests/post/test_evaluation_collective_guard.py`, at `-n 1` and `-n 2`: equal counts with coordinates shifted by `0.01 * rank` inside the unit cube raise on every rank (allgathered flags `== [True] * size`), the message naming each rank's first differing row and deviation, and at `-n 1` the same construction does **not** raise (negative control); a pre-change control — the function at the pinned pre-change sha, `git show` into a scratch module under `/logs/` as `OPS-53` did — returns an all-valid mask on that input **and** values deviating from the closed form `(x, 2y, 3z)` at rank 0's points by the amount the shift predicts (quantitative: the defect was silent, and by how much); *History: `docs/planning/chunks/OPS-55.md`.* | smoke (unit, `-n 1` / `-n 2`) + standard (caller subset at `-n 2`) |
| `OPS-56` | Sweep-planner endpoint inclusion must not depend on the carrier frequency | ✅ Closed 2026-09-21: the sweep planner now decides endpoint inclusion and de-duplication relative to the step rather than the carrier frequency and rejects non-finite inputs, with the pre-change control reproducing all three failures. An auditor DEMOTE at the 2026-09-23 review was overruled and the closure stands. **Done-when:** unit tests asserting exact grids for the three reproduced cases (endpoint present, length as computed by hand), a nondivisible span, the float-noise case `(0.1, 0.3, 0.1)` → 3 points, an on-grid endpoint bitwise equal to `stop_hz`, a coarse/refined near-duplicate merging to one point, and the intended `ValueError` for each non-finite input; a pre-change control showing the three reproduced cases fail on the pinned pre-change sha; the two existing planner tests green unchanged. *History: `docs/planning/chunks/OPS-56.md`.* | smoke |
| `OPS-57` | S ↔ Z conversions: let the solve report singularity; one impedance-to-S implementation | ✅ Closed 2026-09-21: the three determinant singularity tests became LinAlgError-to-ValueError wraps, `s_to_z` and `z_to_s` refuse non-finite input, and `sparameters_from_impedance` delegates to `z_to_s`. The unit round trip and the PORT gate re-runs at -n 2 are green with no tolerance touched. **Done-when:** unit tests — round trip `z_to_s(s_to_z(S))` to ≤ 1e-12 on a 16-port matrix scaled so its determinant underflows (the pre-change control: the pinned pre-change sha raises on it); an exactly singular `I − S` still raises `ValueError` with the same message at all three sites; `sparameters_from_impedance` and `z_to_s` agree to ≤ 1e-14 on a random passive Z. *History: `docs/planning/chunks/OPS-57.md`.* | smoke (unit) + standard (PORT gate re-run at `-n 2`) |
| `OPS-58` | Private-leak check: commit messages, fail-closed on inspection errors, "skipped" is not "clean", and a positive control | ✅ Closed 2026-09-22: the private-leak checker now scans commit messages, fails closed on inspection errors, reports SKIPPED distinct from clean, and has an exit-code test (24/24) plus a narrow shared-input allowlist for the one false positive. Audited PASS 2026-09-23; the operator still needs to re-run `install_git_hooks.sh` once. **Done-when:** that test asserts exit codes exactly — planted figure staged → 1; planted figure in the message → 1; a truncated spelling of it → 1; clean → 0 with `clean:`; no reference data → 0 with `SKIPPED`; broken git → 2 — plus the pre-change control (the pinned pre-change script passes the message plant and the broken-git case); the CI `private-leak` job runs the same test, which becomes the one part of that job that tests something on a runner; `--audit` on this clone still reports clean (the count of figures printed, never a figure). *History: `docs/planning/chunks/OPS-58.md`.* | smoke (host-side, no container) |
| `OPS-59` | Two count pins that the `EX-57` figure task breaks once per figure: make them identities | ✅ Closed 2026-09-21: (b) replaced the literal call-site count with the call-sites equal setup-PNGs equal census-ok identity, and (a) now admits setup figures by rule plus a declared list after the ans6 `metrics.json` was declared on the review's ruling. Both modules are green at -n 2 with the artifact identity reading 28 and the negative controls passing. **Done-when:** (a) the artifact pin admits exactly the files matching `examples/*/figures/*_setup.png` that have a same-stem example script, by rule rather than by list, and nothing else — asserted negative control: a non-setup PNG and an orphan `*_setup.png` injected into the listed set both still fail the pin; count identity extras-before = committed setup PNGs, extras-after = 0. *History: `docs/planning/chunks/OPS-59.md`.* | smoke |
| `OPS-60` | `scikit-rf` in the image, pinned | ✅ Closed 2026-09-23: `scikit-rf` 2.1.0 is pinned in the image and its version test is green on the ordinary service at -n 2, unblocking `PORT-23`. An auditor DEMOTE on the tier word was overruled at the 2026-09-25 review and the tier relabelled standard. *History: `docs/planning/chunks/OPS-60.md`.* | standard (re-labelled 2026-09-25 review; was smoke — the closing window measured 39 s) |
| `OPS-61` | Palace as a second independent reference | ✅ Palace at matched order ties ours to the printed digit on every exact rung bought so far — the ×1 strip and hole meshes at 10 / 64 / 128 MHz, the h = 0.0075 strip rung at all three frequencies (F34, items 418 / 434 / 435), the h = 0.005 strip rung at 10 / 64 MHz (F36 / F39, items 437 / 458 / 458b) and the h = 0.0075 hole rung, PEC and copper (F37, items 454 / 457) — and its driven reads of the ESR-tuned F-small at ×1 and h = 0.01 (T39–T44, items 420 / 421 / 436 / 444 / 448) read f₀ and Q within the zero resolution of ours at both Larmor frequencies. F36, F37 Cu and T42–T44 were ratified by the 2026-10-14 weekly; F37 PEC and F39 landed as §9 items 457 / 458 / 458b (2026-10-15). **Done-when:** Palace `-np 8` solves `f05x` at 64 MHz: 117 206 elements, Rs 43.75 Ω/sq, port i = N_i, `+Z`, 50 Ω. Its 4×4 meets `RECIPROCITY_BAND`, σ_max and every column power ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE`, and every S-class spread ≤ `ADJACENT_SPREAD_BAND` with pooled ≥ `POOLED_SEPARATION_FLOOR` × the worst — each asserted with the band imported. *History: `docs/planning/chunks/OPS-61.md`.* | heavy (step 0 cost probe first; step 1 predicted minutes on the 116 k-cell rung) |
| `OPS-62` | The `ANS-4` ladder's refinement guard counts unknowns, not cells | ✅ Closed 2026-09-22 out of band: the `ANS-4` ladder's refinement guard now compares estimated unknowns instead of cells, strictly stronger, with pure-function controls green at -n 2. Audited PASS 2026-09-23; step 3e was ruled not re-run and 3f ran green under the new guard. *History: `docs/planning/chunks/OPS-62.md`.* | smoke (pure-function controls, `-n 2`, 7 s) |
| `OPS-64` | `xl-run.sh` establishes "no orphans" itself and handles a window's declared outputs | ✅ Closed 2026-09-28: `xl-run.sh` now records container rank counts before and after each window with an orphan verdict, and stages declared `XL_OUTPUTS` or refuses bad ones before the window, with the extended isolation test at 40/40 PASS. The tier was relabelled standard by the 2026-09-28 review. **Done-when:** `xl-run.sh` prints the container-namespace rank count before and after every window into the log it commits, independent of the queued command; and either stages the output paths a queue file declares or refuses *before* the window starts when a declared output is neither ignored nor declared — asserted by extending `test_xl_window_isolation.sh` with one positive and one negative control each, run through the harness, elapsed recorded. *History: `docs/planning/chunks/OPS-64.md`.* | standard *(relabelled from smoke by the 2026-09-28 18:00 review on the auditor's DEMOTE(8), overruled as `OPS-60`'s was: the closing window measured 47 s, `20260928T173030Z_OPS-64.log:126`)* |
| `OPS-65` | Log retention deletes §4 evidence cited only from the chunk history files, review journals and XL ledgers (opened 2026-09-27 by the 03:00 daily review, from … | ✅ Closed 2026-09-27: the housekeeping citing set now includes chunk files, review journals, XL ledgers and xl-pending, 210 deleted logs were restored, and cited-but-missing fell from 284 to 74, all unrecoverable and listed. Audited PASS at the 2026-09-28 review. **Done-when:** a unit test asserts that `cited_basenames()` contains a basename cited only from a chunk file and one cited only from the `xl` ledger; its control, pinned to the pre-change commit (`<sha>^`), is red; the dry-run census reads cited-but-missing = 0 (bar the unrecoverable, listed in `docs/testing/housekeeping.md`); elapsed recorded. *History: `docs/planning/chunks/OPS-65.md`.* | smoke |
| `OPS-66` | *(history and §9 briefs: `docs/planning/chunks/OPS-66.md`)* The known-issues rotation finds undated RETIRED headings | ✅ Steps 1 and 2 closed: `KI_ENTRY_RE` admits the undated-retired headings (audited PASS 2026-10-01) and `chunks` strips a state line's leading row glyph, refusing a different one (2026-10-11). **Done-when:** the entry pattern becomes a module constant `KI_ENTRY_RE` that also admits `## [✅ ]RETIRED\|RESOLVED\|FIXED\|CLOSED <date>`; a unit test shows (a) the four headings, as verbatim literals, are entries and are closed; (b) on the live file and `known-issues-archive.md` every heading the old pattern matched is still matched, new − old is exactly the closed undated-retired headings, and no section header becomes an entry; (c) the pinned pre-change tool's dry-run on a copy of the live file reports none of the four and the new one all four; no rotation is run (the weekly's) ⇒ ✅ (§9 item 139). *History: `docs/planning/chunks/OPS-66.md`.* | smoke |
| `OPS-67` | *(history and §9 briefs: `docs/planning/chunks/OPS-67.md`)* **known-issues: retire three live entries whose resolutions already landed** (opened 2026-10-07 10:30 daily review, §9 rule (5) lane 1) | ✅ 2026-10-07: all three met their own Resolves-with and are re-marked `✅ RETIRED` in place, bodies verbatim (`GEO-26` `2445b9ed` / `20260904T034441Z_GEO-26.log`; `OPS-28` `fd1479ca` / `20260828T033055Z_OPS-28-gate.log`; `TH-13` 3a‴ `f494dc75` / `20260901T0334*Z_TH-13-step3a3-*.log`); W0 `20261007T185729Z_OPS-67-W0.log` lists 1 closed, W1 `20261007T185849Z_OPS-67-W1.log` lists 4, W1 − W0 = exactly the three (set equality); live-vs-archive diff `20261007T185802Z_OPS-67-diff.log` (not byte-identical, nothing deleted); no rotation run. **Done-when:** each entry is either re-marked `✅ RETIRED` with a retired-by line citing a landed log that meets its own "Resolves with", or annotated with what is still owed; the rotation dry-run's closed list grows by exactly the re-marked set (W1 − W0); no rotation run ⇒ ✅. | smoke |
| `OPS-68` | *(history and §9 briefs: `docs/planning/chunks/OPS-68.md`)* **known-issues: retire four more live entries whose resolutions already landed** (opened 2026-10-07 18:00 daily review, §9 rule (5) lane 1) | ✅ 2026-10-07 (19:30 slot): all four entries re-marked `✅ RETIRED 2026-10-07 (OPS-68)`, each condition met by a cited green log and ✅ row; dry-run closed list W0 4 → W1 8, W1 − W0 = exactly the four (`20261008T005516Z_OPS-68-W0.log`, `20261008T005620Z_OPS-68-W1.log`, Status 0); no rotation run. Was: `MAG-18` anchor (ii), the `OPS-19` docrefs exit-code contract, the four `OPS-17` step 2 defects and the `EX-49` control-floor observation each record their resolution in the body but keep live headings (sweep, 2026-10-07 18:00 journal §5). **Done-when:** each entry is either re-marked `✅ RETIRED` with a retired-by line citing a landed log and a ✅ row that meet its own resolution condition, or annotated with what is still owed; the rotation dry-run's closed list grows by exactly the re-marked set (W1 − W0); no rotation run ⇒ ✅. | smoke |
| `OPS-69` | docrefs `dead=2`: the `ans:9` SPEC names gitignored private results files *(brief: `docs/planning/chunks/OPS-69.md`; opened 2026-10-08 18:00 review)* | ✅ The known-issues entry 2026-10-08 (found by `EX-89`) traces both dead references to `examples/ansys_benchmarks/ans9_coil_driven_gel_heating_10MHz/SPEC.md` :36 / :38 (operator commit `4d4d93c5`), which name two gitignored JSONs as filenames; :39 names a private script a fresh clone lacks. **✅ 2026-10-08 (§9 item 336, 21:00 slot):** :36–39 reworded to descriptions (no filename token, no new digit; disclosed in commit and journal only, ruling Q3); census `dead=2 … exit=1` (W0, `20261009T021416Z_OPS-69-W0.log:98`) → `dead=0 guide=0 stale=56 … exit=2` (W1, `20261009T021429Z_OPS-69-W1.log:96`, 2 s); `test_doc_reference_exit_codes.py` 20 passed, tree anchor and `[dead-artifact|dead-script]` controls green (W2, `20261009T021433Z_OPS-69-W2.log`, 13 s); entry RETIRED in place, W3 dry-run lists it (`20261009T021502Z_OPS-69-W3.log`). **Done-when (met):** :36–39 reworded to form no filename token (meaning kept, no number added, change disclosed in the commit and journal); census W1 reads `dead=0`, exit ≠ 1 against W0's `dead=2 … exit=1`; `test_doc_reference_exit_codes.py` green including the planted-dead cases; entry re-marked RETIRED in place ⇒ ✅. | smoke |
| `OPS-70` | Package regression re-gate after the 2026-10-08 `study.py` / `cli.py` edits (`e763fa93`, `24728153`, `df25001f`, `ee87e99b`, `b9ad9e94`, `384b0a88`, … | ✅ 2026-10-09 (§9 item 348): step 1 green whole on `H` = `1edad4ac` — 25 windows, 278 passed / 0 failed / 10 skipped, 1 993 s summed, controls at priors. Step 2 (degree-2 records) is written, not queued. **Done-when:** step 1 — every `tests/package/` module plus `tests/post/test_wf8_figures.py` run through the harness on one recorded HEAD (every log header `Commit:` equal), each degree-1 / non-degree-2 node at its gated width and build, all green with each module's own records, identities and asserted controls at their imported bands (no assertion added, no band touched), elapsed recorded ⇒ ✅, and the `WF-8` / `WF-15` step 3 auditor caveats closed in their chunk files. Any red ⇒ known-issues entry plus a §7 note on the owning row. Not claimed: degree-2 records (step 2, written, not queued). *History: `docs/planning/chunks/OPS-70.md`.* | heavy |
| `OPS-72` | Make `tests/environment`'s DolfinX API-migration sweep cheap enough for a standard window again (opened 2026-10-26 weekly, lane 1; owner of known-issues 2026-10-26) | ✅ 2026-10-27 (§9 item 665): full `tests/environment` 44.90 s at `-n 1`, 44.73 s at `-n 2` (harness footers `20261027T124941Z_OPS-72-W1-n1.log`, `…125038Z_OPS-72-W2-n2.log`), verdicts and CLI `RESULT:` unchanged, negative control green; known-issues 2026-10-26 retired; the environment-first ruling is the next daily's to lift. Was: `tests/environment` costs ≈ 145 s at `-n 1` because `scripts/testing/check_dolfinx_api_migration.py` parses and resolves every `.py` under `src` + `tests` twice (baseline and the negative control's temp copy) and sweeps `examples` + `scripts`, on every rank. **Done-when:** the sweep runs once (rank 0, result broadcast) and skips nothing it checked before — the planted negative control still red, every `tests/environment` verdict unchanged, the checker's `RESULT:` line byte-identical on a clean tree — and the full `tests/environment` at `-n 1` and `-n 2` measures ≤ 60 s (harness footer), with the 145 s run cited as the before; known-issues 2026-10-26 retires on that measurement and the environment-first ruling (03:00 daily, 2026-10-26) is then lifted by the daily. Not a change to any gate; the checker's scope (`src`, `tests`, `examples`, `scripts`) is kept. *History: `docs/planning/chunks/OPS-72.md`.* | standard by measurement (harness 46 s at both widths; re-declared from smoke by the 2026-10-27 10:30 audit), `-n 1` / `-n 2`, real |
| `OPS-71` | Re-record `mri:2`'s guide table to the 0.11 image under the `(1*)` licence (opened 2026-10-19 18:00 daily review, lane 1; §9 item 552) | ✅ 2026-10-20, run of record `20261020T004018Z_OPS-71-mri2.log` (Status 0, 10 s; cells 74 218, DG0 2.81e-15, 10 g kernel 0.0041 %; guide updated, 2026-08-08 figures kept as history). Was: the guide's "On record" table is the 2026-08-08 run. Item 549's re-gate (`20261019T170830Z_POST-8-step3b-ii-e1b-mri2-regate.log`, Status 0) reads three rows differently: cells 74 218 against 74 216, the DG0 round-off 2.81e-15 against 1.32e-15, and the 10 g kernel 0.0041 % against 0.0044 %; every band still holds. **Done-when:** a fresh run of record of the unedited `examples/mri/02_mass_averaged_sar.py` at `-n 2` passes all of its own assertions at unchanged bands; every digit in the guide's record table and printed blocks equals that log, with a guide-line → log-line trace in the journal and the 2026-08-08 figures kept as history; docrefs `exit != 1` on both sides ⇒ ✅. Any red ⇒ a known-issues entry and no re-record. Not claimed: other guides, or any band or gate. *History: `docs/planning/chunks/OPS-71.md`.* | standard |
| `OPS-73` | A retention-safe reader for tests that read a dated harness log (opened 2026-11-01 10:30 daily review, lane 1; owner of known-issues 2026-11-01) | ✅ 2026-11-01 (§9 item 724): `tests/log_records.py` reads a record in either form; the node passes with only the `.gz` on disk, (H) and (A) asserted (`20261101T182231Z_OPS-73-Wreg.log`, `20261101T182228Z_OPS-73-W1.log`). `test_ans8_v3_d2_block` stays red on (R2), a separate entry in known-issues. *History and brief: `docs/planning/chunks/OPS-73.md`.* **Done-when:** `tests/log_records.py::read_log_text` is the only route those modules use; the node passes on `main` with no plain `.log` present; (H) plain = `.gz` text and (A) both-absent refusal asserted; the known-issues entry retired in the same commit. | standard (`-n 2`) |
| `OPS-74` | The four example closures' owed windows, and the first real `./run_examples.sh` through the harness (opened 2026-11-01 18:00 daily review; the 10-31 weekly's closing-window rule) | ✅ 2026-11-01 §9 item 738: four Status 0 windows on `./scripts/run_examples.sh` through the harness (`20261102T013021Z_EX-169-close.log`, `…013146Z_EX-166-close.log`, `…013342Z_EX-168-close.log`, `…013649Z_EX-167-close.log`), the bad selector Status 2 (`…013018Z_OPS-74-Wc.log`), no `-T` needed; the 10-31 rule's form is observed to run. *History and brief: `docs/planning/chunks/OPS-74.md`.* **Done-when:** each of the four has a Status 0 harness log whose Command line is `./scripts/run_examples.sh …`, cited in its §7 row beside the existing log; or the failure is recorded and the rule returned to the weekly with the measurement. | heavy (`-n 2` / `-n 4`) |
| `OPS-76` | `femem run` on a swept file with no `tune.at_hz` ends in an uncaught `ValueError` (opened 2026-11-02 03:00 daily review, lane 1; item 744's disclosure, known-issues 2026-11-02) | ✅ 2026-11-02 §9 item 758 (`20261102T180047Z_OPS-76-item758-W1.log`, `-n 2`, complex, Status 0, 87 s, 26 passed, 3 skipped per rank: rc (2, 2) with `tune.at_hz` in stderr, 0 calls, no output, the API's plain `ValueError` unchanged, the one-frequency control hex-equal to the API; item 754's pin now asserts the return code). *History and brief: `docs/planning/chunks/OPS-76.md`.* **Done-when:** `femem run` on that file exits 2 and names `tune.at_hz` on every rank with no solve and no output; the API's `ValueError` unchanged; a single-frequency file still tunes; the known-issues entry retired. | standard (`-n 2`) |
| `OPS-77` | How close a 1e-12 cross-solve gate sits to the known MUMPS drift (opened 2026-11-02 03:00 daily review, lane 1; `OPS-74`'s disclosure, known-issues' `OPS-43` observation) | ✅ 2026-11-02 as a record, the tail branch: five samples 1.295e-15 … 1.147e-13, the largest 0.115 × `IDENTITY_BAND`, none red (`20261102T122110Z_OPS-77-reader.log:36–44`); no band moved, §9 item 755. **The 11-02 weekly ruled the design stands** (no band, no gate redesigned) and made these five samples the first reading of any future separate-solve (S) red: inside the drift ⇒ the observation firing, filed there; outside ⇒ a defect. *History and brief: `docs/planning/chunks/OPS-77.md`.* **Done-when:** a table of every recorded (S) reading plus two identical new windows, largest read against the imported `IDENTITY_BAND`; the observation gains the row; the pre-registered branch named; no band moved. | heavy (`-n 2`) |
| `OPS-78` | Which identity gates compare two separate solves (opened 2026-11-02 10:30 daily review, lane 1; `OPS-77`'s tail branch, known-issues' `OPS-43` observation) | ✅ 2026-11-02, §9 item 767, attempt 2 of 4: 369 lines classified, 0 unread, counts equal `grep` per name, three controls as required (92 separate solve, 51 shared route, 217 not a solve output, 9 out of scope); 21 separate-solve gates tabled, none above its band, tightest `MAT-21` (H0) at 0.38 × band, 2.6 × headroom (`20261103T013749Z_OPS-78-census.log`, `Status: 0`, 8 s); no band moved. **The 11-02 weekly ruled the design stands** (`IDENTITY_BAND` across two factorisations sits ≈ 9× inside the `OPS-43` exposure line on the worst sample) and made this table, beside `OPS-77`'s samples, the reading of any future separate-solve red; `mri:37` vs `mri:38` stays a recorded difference with no claim. **Refreshed at HEAD 2026-11-04 (§9 item 815):** 402 lines, 0 unread, exit 0 (`20261104T180148Z_OPS-78-census-r2.log`, 10 s). Three new separate-solve gates (`WF-17 s28` (Id) 6.377e-14, `WF-15 s13` (Id) 4.722e-16, `WF-16 s18` (I) 0), none above band. Four pins re-pointed. **`mri:38` (S) now reads 1.076e-13 on one window (`20261103T225531Z_EX-170.log:5646`), so it is 6 of 8 zero, not 6 of 6**: a record for the weekly. *History and brief: `docs/planning/chunks/OPS-78.md`.* | smoke |
| `OPS-79` | The two tracked closing-window texts amended to the 11-02 ruling | ✅ 2026-11-03 §9 item 792: both texts carry the three clauses citing the ruling; the script's diff is comment-only (30 lines, `20261103T223106Z_OPS-79-W2.log`), its `-e mat:2 --dry-run` byte-identical before/after (`…223107Z_OPS-79-W3.log`), `bash -n` Status 0 (`…223105Z_OPS-79-W1.log`), a one-character control differs (`…223108Z_OPS-79-W4.log`); `example-runner.md` untouched (W4). Opened 2026-11-03 03:00 review from the 11-02 weekly's closing-window ruling (§10): `scripts/run_examples.sh:30–43` and `docs/automation/implementer-run.md:133–150` still say emit-then-harness; `example-runner.md` waits for the operator's worktree session (W4). **Done-when:** both texts state the ruling's three clauses — (i) `./scripts/run_examples.sh -e <id> -n <n> -t <s>` through the harness is the closing window, (ii) a setup figure has its own earlier flagged window that closes nothing and the closing window is last, (iii) a guide edited after its closing window admits only digits in the closing log or arithmetic on them, said so — each citing the 2026-11-02 weekly ruling; no run; `.claude/agents/example-runner.md` untouched. | smoke (docs, no compute) |
| `OPS-80` | Which callers of the two-tag power-identity helper (`_exact_shares` / `_loss_power_w`) run on a mesh with lossy cells outside tags 3 and 1 (opened 2026-11-04 18:00 daily review, lane 1; known-issues 2026-11-03's "Resolves with") | ✅ 2026-11-05, §9 item 817 (`20261105T013028Z_OPS-80-grep.log`, Status 0, 22 s): all 46 lines in 16 files classified (`_exact_shares` 35, `_loss_power_w` 11; 3 definitions, 9 imports, 14 prose, 20 calls; 0 unclassified). Both controls hold (entry48 is (K) all-cells with a two-tag control; `test_imp3_rod_birdcage.py` is (U), a `_build_rod` hole); **no other (K) two-tag gate**. Not widened (11-04 weekly; `OPS-81` documents the scope). **Audited PASS 2026-11-05 03:00 review.** *History and brief: `docs/planning/chunks/OPS-80.md`.* **Done-when:** every line classified (definition / import / prose / call; each call (U) / (K) two-tag / (K) all-cells / (?) with file:line), totals equal `grep` per name, 0 unclassified; `test_imp3_tip_q_entry48_def.py` reads (K) all-cells with its two-tag control (K), `test_imp3_rod_birdcage.py` reads (U); the known-issues Scope gains the result; no test, helper, band or `src/` edit. | smoke (no compute) |
| `OPS-81` | **The two-tag power-identity helper states its scope and refuses a lossy tag outside (1, 3)** (`_loss_power_w` / `_exact_shares`; opened 2026-11-05 03:00 daily review, lane 1, on the 2026-11-04 weekly's ruling; retires known-issues 2026-11-03) | ✅ 2026-11-05, audit PASS (18:00 review) (§9 item 828, 12:00 slot; W0 12 passed per rank, W1 / W2 / W3 verdicts unchanged, heavy, 302 + 310 + 134 s; the ack key awaits the 11-07 weekly's ratification). The case carries its material map (`case["problem"]`, `test_port_birdcage_four_port.py:469–487`), so the guard applies; the two (K) wrappers that map 403 lossy (`_arm_rod`, `_arm_with_403`) acknowledge it on the case (review ruling 2026-11-05 03:00, for the weekly to ratify). *History and brief: `docs/planning/chunks/OPS-81.md`.* **Done-when:** both docstrings state the two-tag scope and name the all-cells form; `_loss_power_w` raises `ValueError` naming the tag and the all-cells form on a synthetic case with 403 at σ > 0 unacknowledged, and does not raise with 403 at σ = 0, acknowledged, or without `problem` (asserted); the existing callers' gates unchanged — `test_imp3_tip_q_entry48_def.py` 14 passed / 3 skipped, `test_imp3_rod_heating_matched.py` 16 passed, `test_birdcage_power_identity.py` 0 failed; no band, gate or `src/` edit; the known-issues entry retires in the same commit. | heavy (`-n 2`, ≤ 600 s a window) |
| `OPS-82` | **The `IMP-8` DG0 source helpers build in both builds** (`_rod_source` / `_tips_source`; opened 2026-11-05 18:00 daily review, lane 1, from `EX-182`'s finding; retires known-issues 2026-11-05 "real-build only"). *History and brief: `docs/planning/chunks/OPS-82.md`.* | ✅ 2026-11-06 (§9 item 849, 21:00 slot): both helpers assemble `one * ufl.conj(w) * ufl.dx`; `tests/unit/test_ops82_dg0_source_helpers_both_builds.py` green in both builds at `-n 2` — (V) rel 3.614e-16, (S0) rod / balls ≤ 7.772e-16, complex refusal and real bitwise controls green (`20261106T030402Z_OPS-82-W2-complex-r2.log`, `20261106T030411Z_OPS-82-W1-real-r2.log`, 3 s each; audit PASS 2026-11-06 03:00 review). **Done-when:** both assemble `one * ufl.conj(w) * ufl.dx`; a unit module gates (V) the owned DG0 volume sum = the box volume (1e-12) and (S0) ∫q = P, each ball P/2 (`B_BAND`) in the real and complex builds, with the old form's complex refusal and the real-build vector's bitwise equality as controls; one line per helper module, no `src/` edit; the known-issues entry retired. | smoke (`-n 2`, both builds) |
| `OPS-83` | **The (S2-tip) comparand's form on a rough DG0 source: a stub measurement** (opened 2026-11-06 03:00 daily review, lane 1, from the known-issues entry of that date — `xl` entry 55's and U127's tip-slope reds). *History and brief: `docs/planning/chunks/OPS-83.md`.* | ✅ 2026-11-06 (§9 item 859, 06:00 slot; audit PASS, 10:30 review): all gates green at `-n 2` complex, 13 s standard (`20261106T120232Z_OPS-83-W1.log`); reading **form CONFIRMED** — on the checkerboard R_DG0 1.58e-03 / 2.11e-03 Δt-independent, R_Π 6.8e-06 / 7.9e-06 scaling ≈ 9.3–9.9× with Δt; known-issues Cause line updated, entry stays open for the weekly's re-pose. **Done-when:** `tests/validation/test_ops83_s2tip_form_on_stub.py` green at `-n 2` (complex) on I32a's stub — (E), the global (S2) at `S2_BAND` on a uniform and a checkerboard DG0 source at Δt₁ and Δt₁/10, K10-S's tip 8.0400e-06 reproduced, the projection identity on the uniform source — with R_DG0 / R_Π / R_L and their Δt-ratios printed and the pre-stated reading line; the known-issues entry's Cause line updated with the measured reading ⇒ ✅. Not claimed: a re-pose, a band, a clinical number. | standard (`-n 2`, complex) |
| `OPS-84` | **The (S2-tip) comparand's form on a field-derived source: I41's matched rod pair** (opened 2026-11-06 10:30 daily review, lane 1, from the known-issues entry of 2026-11-06 — its Cause line's "Not measured" clause after `OPS-83`). *History and brief: `docs/planning/chunks/OPS-84.md`.* | ✅ 2026-11-06 (12:00 CST slot, §9 item 860; audit PASS, 18:00 review; `20261106T180958Z_OPS-84-W1.log`, 14 passed, 267 s at `-n 2`, heavy): on I41's matched rod pair every asserted gate is green (both arms' ccw (S2) reproduce `OPS-81` W2 to 2.2e-12 / 2.3e-13) and the tip gap is the comparand's form, R_DG0 1.6e-04 … 9.4e-03, Δt-independent on the seven classified rows (the INCONCLUSIVE row's ratio is 2.2), with R_Π ≤ 1.52e-04 at Δt₁ scaling O(Δt) and ≤ 1.53e-05 at Δt₁/10 (no RESIDUAL row). The known-issues entry stays open for the 11-07 weekly's re-pose; not claimed on U127's or entry 55's own mesh. | heavy (`-n 2`, complex) |
| `OPS-85` | **An exercising node for the `_tuned_quadrature_temperature` gradient refusal** (opened 2026-11-08 03:00 daily review, lane 1, from U125's disclosure (3) as ratified by the 2026-11-07 weekly: `src/fem_em_solver/study.py:3396–3397` refuses `TransientSpec.initial_gradient_k_per_m` on a tuned two-port Study, and no node reaches it). *History and brief: `docs/planning/chunks/OPS-85.md`.* | ✅ 2026-11-08 (§9 item 879, 09:00 slot; audited PASS 10:30: exact message, frame and call counts meet §4 in its operative sense — the first refusal-only closure, a precedent the 11-09 weekly ratifies; the refusal now sits at `study.py:3533–3534`): `tests/package/test_ops85_tuned_gradient_refusal.py` raises the exact message from the `_tuned_quadrature_temperature` frame for g = 10.0 and 0.0 with `_tuned_result` / builder calls 0, and the gradient-free control reaches `_tuned_result` once (`20261108T152052Z_OPS-85-W1.log`, complex, 8 passed per rank, 3 s; `20261108T152105Z_OPS-85-W-real.log`, real, 5 passed / 3 skipped per rank). **Done-when:** one node at `-n 2` in the complex build, with no mesh. The tuned F-human Study's `temperature(…, transient=TransientSpec(…, initial_gradient_k_per_m=g))` raises `ValueError` with the exact message, from the `_tuned_quadrature_temperature` frame (the same message also stands at `:1659` and `:2139`), before `_tuned_result`. The gradient-free spec on the same call reaches `_tuned_result`, and builder calls are 0 ⇒ ✅. The audit rules whether that call-count identity meets §4's quantitative clause; if it does not, the row stays 🟡. Not claimed: any solve, a gradient on a tuned Study, an envelope line. | smoke (`-n 2`, complex) |
| `OPS-75` | The two example scripts that read a dated log by its plain name (opened 2026-11-01 18:00 daily review, lane 1; `OPS-73`'s disclosure 3) | ✅ 2026-11-02 (`20261102T062832Z_OPS-75-W1.log`, `-n 1`, complex, Status 0, 3 s, 10 passed: (H′) dict and text equal on `.gz`-only copies, (A′) refused, the bare `open` control raises, and `H_TABLE_LOGS` read directly reproduces `h_table_s.json` byte for byte with two of its three listed names absent on disk; `…062812Z_OPS-75-Wreg.log`, `…062835Z_OPS-75-Wreg-n2.log`). *History and brief: `docs/planning/chunks/OPS-75.md`.* **Done-when:** both read a record in either form; (H′) equal text and (A′) both-absent refusal asserted; the tests' scratch copies untouched. | standard (`-n 1`) |
| `OPS-63` | The `ans:2` example says ASSERTED only where it asserts, and carries no AED figure in tracked text | ✅ Done 2026-09-25 with audit PASS 2026-09-26: the `ans:2` example now asserts every quantity it labels ASSERTED, importing each band from the gate modules, and the AED figure is replaced by a pointer to the private adjudication file. The knob-route run reproduced the tracked metrics within the identity tolerance, the phrase census reads 0 after the change, and the `OPS-58` checker's failure to flag the figure is recorded as a finding. **Done-when:** (i) the knob route at the tracked rung (`FEM_EM_ANS2_PHANTOM_RESOLUTION=0.0025`, `-n 4`, the route step 4a proved digit-identical to the unset run) runs green with the new asserts, and its `metrics_h0.0025.json` reproduces the tracked `metrics.json` on every primary leaf within the imported `EXACT_IDENTITY_RTOL` (step 4a measured worst 8.19e-11) — no number moved; *History: `docs/planning/chunks/OPS-63.md`.* | heavy (the example at the 0.0025 rung ≈ 200 s at `-n 4`) + smoke |
| `ARCH-1` | The MRI-safety layer becomes solver-agnostic | 🟡 Step 1 (U16) and steps 1c–1e landed: `Study.solve(backend="palace")` / `femem run --backend palace` write the Palace case, a separate window runs it and the `Study` reads S back, tied to our degree-2 record; the route refuses every unrouted output, `as_hole` and a multi-rank launch by field before any mesh; the copper hole on the Palace backend (step 1f, §9 item 631: `as_hole` + `LeontovichWall` → `Boundaries.Conductivity`, (E2-Cu) 0.000e+00 against `OPS-61` step 2b's read). Not ✅: the gate of record is the whole safety layer on two backends. **Done-when:** digit-identical reproduction of every imported band and record on the DolfinX backend before any second backend is wired — a refactor that moves a record is a defect. *History: `docs/planning/chunks/ARCH-1.md`.* | unmeasured — the 10-17 weekly costs it |
| `OPS-1` | Executable verification environment (Docker) | ✅ | smoke |
| `OPS-2` | CI runs the real test suite, not just `tests/unit` | ✅ | standard |
| `OPS-3` | Deterministic test tolerance policy | ✅ | smoke |
| `OPS-4` | Lightweight smoke matrix | ✅ | smoke |
| `OPS-5` | Testing status dashboard | ✅ | smoke |
| `OPS-6` | Expanded run-and-log metadata | ✅ | smoke |
| `OPS-7` | Guided pending-test queue helper *(retired 2026-08-04 — queue tooling removed; verification is agent-executed per §4; 🚫 by the 2026-09-30 weekly — superseded by §9's On-deck queue, nothing left to gate)* | 🚫 | smoke |
| `OPS-8` | v1 milestone acceptance checklist | 🧪 | smoke |
| `OPS-9` | Prune duplicate/stale entries from `pending-tests.md` | ✅ | smoke |
| `OPS-10` | Complex-mode CI job for the frequency-domain gates | ✅ | smoke |
| `OPS-11` | Put `tests/mesh` in CI — the directory no job runs | ✅ | smoke |
| `OPS-12` | Adjudicate the residual-trend classifier (known-issues 2) and return `test_convergence_diagnostics.py` to CI | ✅ 2026-08-08 | standard |
| `OPS-13` | Land the rank-safe `_validate_material_map_tags` fix on `main` with its own gate | ✅ 2026-08-08 | standard |
| `OPS-14` | Diagnose the rank-dependence of `test_single_port_excitation` (known-issues 6) | ✅ | standard |
| `OPS-15` | Retire the checker's standing freshness tax: default `--max-age-s` 1 h → 48 h | ✅ 2026-08-10 | smoke |
| `OPS-16` | Retry-on-529 in the three automation launchers (two review slots lost 2026-08-13; rubric in the §9 item) | ⬛ **WON'T FIX — operator decision 2026-08-22** | smoke |
| `OPS-17` | Delete or replace the finiteness-only test suites (operator directive 2026-08-16) | ✅ Closed at the 2026-08-21 18:00 review: every runnable validation test was observed in a completed complex run, with two files formally deferred and the denominator re-based. `OPS-18` was queued per the commitment. *History: `docs/planning/chunks/OPS-17.md`.* | standard |
| `OPS-18` | DolfinX version upgrade, recurring (0.7.2 → newest qualifying; operator directive 2026-08-16) | ✅ Closed 2026-08-23: `main` boots the 0.11 image with every §2.1 family re-gated green on it, and the real-mode `MAG` scope caveat was discharged the same day. *History: `docs/planning/chunks/OPS-18.md`.* | heavy |
| `OPS-19` | Doc-reference checker: staleness must not own the exit code (2 runs flagged the masked signal 2026-08-16) | ✅ (2026-08-16: exit 0/1/2 split + `--stale-severity {fail,report}` default `report`; on `main` the checker now reads `dead=0 guide=0 stale=24 exit=2` where it read exit 1, guide pass green 21/21; 8 tests, 1.91 s, smoke) | smoke |
| `OPS-20` | Disposition the coil-phantom `ComplexComparisonError`: localize with `--tb=long`, then fix the form or mark `@real_only` (known-issues 2026-08-18; commissioned 2026-08-18 10:30 review) | ✅ *(2026-08-19, 06:00 slot — fixed, no `@real_only`; the commissioned `ComplexComparisonError` was already dead, killed by `OPS-22` through the **imported** drive callable, and a free grep found it. Only the predicted second layer remained: the complex build now **passes the same 30% gate at the same 17.1233%**, both ranks identical, real-mode digits unmoved across control and re-run; collect count unchanged at 49. *Audited COMPLIANT 2026-08-19 10:30 review — all five footers, the 17.1233% in all three counted runs and the line-142 `Im`-assertion verified; the skipped cold-cache clear is journalled in three places and adjudicated sound; cosmetic journal error on record: the exit-124 batch log does carry a footer (124 / 481 s), the uncounted disposition stands*)* | standard |
| `OPS-21` | Make the combined-XDMF test scalar-type-aware and rank-deterministic (known-issues 2026-08-18, two defects in one test; commissioned 2026-08-18 10:30 review) | ✅ | standard |
| `OPS-22` | Make the three magnetostatic loop-drive fixtures complex-safe: replace the `ufl.max_value` / `<=` predicates in their `current_density` callables (known-issues 2026-08-19; commissioned 2026-08-19 03:00 review from the `OPS-17` leg-(b2) attempt-2 diagnosis; unblocks 5 tests in leg (b2)) | ✅ *(2026-08-19, 04:30 slot — all three files fixed, no `@real_only` needed; real-mode digits unmoved to the last printed figure across three runs, and the complex build now runs all three files to a footer: **5 passed, 412.12 s, exit 0**, both ranks identical. *Audited COMPLIANT 2026-08-19 10:30 review — footers, closed-form assertions and the new `Im`-bound idiom verified against all five logs; one caveat on record: `test_helmholtz_v2.py`'s complex coverage rests on a silenced `ComplexWarning` `float()` cast, not an assertion — fold an `Im`-bound in whenever that file is next touched*)* | standard |
| `OPS-23` | Sweep the `OPS-21` rank-0-return defect pattern (4 measured sites in 3 test files) + the `test_helmholtz_v2.py` Im-bound (commissioned 2026-08-20 03:00 review from the 00:00 slot's grep survey) | ✅ | smoke-to-standard | 3 real sites (all in `test_csv_export_stats_parity.py`) + the Im-bound fixed; 2 of the commissioned sites were print-only false positives and 1 exempted site was a real defect; 12 passed both ranks, 5.00 s. *Audited COMPLIANT 2026-08-21 18:00 review — red-baseline byte-identity re-verified against the log's rank blocks; benign omission: the first exit-0 helmholtz-real log is in test-results.md but uncited in the annotation* |
| `OPS-24` | Migrate `core/cavity.py`'s two `assemble_matrix(..., diagonal=)` sites to the 0.11 signature — `TH-9`'s cavity gate + resonance guard have been **non-executing on `main` since the 0.11 merge** (known-issues 2026-08-24; found by `EX-30` leg (th); commissioned 2026-08-24 18:00 review; closed 2026-08-25 — `diagonal=` → `diag=`, all four green, 0.0436% worst-mode reproduced to the printed digit) | ✅ | standard |
| `OPS-25` | Re-join `th:7` to its gate: hoist the series-interior interpolation into the gate module and import it, migrating the repo's only `interpolate(cells=)` site (known-issues 2026-08-24; ruled hoist-not-repair by the 2026-08-24 18:00 review) | ✅ (2026-08-25: hoisted to `series_interior_function` in `test_lossy_sphere_fullwave.py`; `th:7` green in 14 s with both element-order records reproducing — degree 1 8.1541% / 8.3869%, degree 2 0.1405% / 0.0058%, drifts ≤ 1.48e-03 in a 1% band — and the gate's `P_series(meshed)` **bit-identical to all ten printed digits** across the refactor; `13 passed in 25.28s`) | standard |
| `OPS-26` | Systematic dolfinx-0.11 migration completeness sweep | ✅ Closed: `OPS-17`'s completed-run census was re-run on the 0.11 image and `src/` was swept for un-migrated call sites. The numbered findings and the rules they set are in the history file. *History: `docs/planning/chunks/OPS-26.md`.* | heavy (split across ≥ 2 slots) |
| `OPS-27` | Re-record the 0.7.2-era exact records the `OPS-18` re-record did not reach, version-tagged on the `GEO-16` precedent, and sweep for siblings | ✅ Closed: the 0.7.2-era exact records the `OPS-18` re-record did not reach are re-recorded version-tagged, and the prose residue is left as a coupled-constant job for whichever chunk re-prices those fixtures. *History: `docs/planning/chunks/OPS-27.md`.* | heavy (both steps — step 1 measured 256 s on its `mesh_cache` window, over the 180 s standard ceiling; label corrected 2026-08-28 03:00 audit) |
| `OPS-28` | **Give `tests/ports/test_port_orientation_sensitivity.py`'s `_DummyComm` the `allgather` that `OPS-14`'s rank-safety reduction calls, then read the module's real assertions back against known-issues entry 3** — census leg (b) finding 12: a correct reduction outgrew a test double, a class step 1's static sweep cannot see. The reduction stays; the deprecated placeholder route stays runnable (`PORT-1` step 4's negative control). Commissioned 2026-08-27 10:30 review; full rubric in §9 item 3. ***✅ 2026-08-28, 22:30 implementer slot** — one added `staticmethod allgather(value) -> [value]` on `_DummyComm`, nine lines including its comment, `src/` untouched (`git diff -- src/` empty). Bracketed by measurement on the identical command (`tests/ports`, `-n 2`, real, smoke, `-k 30 120`): red baseline `3 failed, 14 passed in 1.50s` / Status 1 / 3 s, gate `2 failed, 15 passed in 0.79s` / Status 1 / 2 s. The sign-flip anchor is **green** — `V(P2) = +5.000000e-02 V` aligned vs `−5.000000e-02 V` flipped, magnitudes equal to `rel=1e-12`, coupling factor `+1.0e-01 → −1.0e-01`. The S-matrix name reaches its assertion for the first time since `OPS-14` and is **red there**, so entry 3 is re-dated, not retired — with a correction it measured: on that 2-port fake the diagonal is **not** zero (`S11 = S22 = 9.047e-01 − 1.289e-02j`); the **off-diagonal** is, because the undriven port is the matched one (`V = 5.000000e-02 = Z₀I` at `Z₀ = 50 Ω` ⇒ `b = 0` exactly). Entry 3's mechanism is confirmed, its old title was imprecise for this name, and its disposition is unchanged (`PORT-0`/`PORT-1`). The leg (b) `allgather` known-issues entry retires whole. Negative control: the other three `tests/ports` modules unchanged — `sparameter_assembly` still 3 passed / 1 failed (entry 3's other name), planner 3 and `port_definition` 8 green in both runs. Logs `20260828T033037Z_OPS-28-red-baseline.log`, `20260828T033055Z_OPS-28-gate.log`.* | ✅ | smoke |
| `OPS-29` | **Rank-safe the `phantom_material` empty-tag check in `build_material_fields`** — the `OPS-13` defect survived 20 lines below its own fix; measured breaking `examples/mri/01_coil_phantom_fields.py` at `-n 12` (interactive session, 2026-08-28) | ✅ 2026-08-28 | smoke |
| `OPS-30` | Migrate the two filed `scripts/probes/` survivors to dolfinx 0.11 | ✅ Closed: both filed `scripts/probes/` survivors are migrated to dolfinx 0.11, and the audit re-tiered the row to standard with no demotion. *History: `docs/planning/chunks/OPS-30.md`.* | standard (re-tiered from smoke 2026-09-02, 37 s measured) |
| `OPS-32` | Private-mode comparison writers for `ANS-3` and `ANS-4` | ✅ Closed: the private-mode comparison writers for `ANS-3` and `ANS-4` landed with no AED value in any tracked table, as ratified at audit. The control half's follow-up is `OPS-33`. *History: `docs/planning/chunks/OPS-32.md`.* | standard (measured 171 s + 172 s + 1 s) |
| `OPS-33` | Re-base `ans:3`'s four `RECORDED_*` constants to the 0.11 image under the (1*) licence, then register the 1e-6 run-to-run control `OPS-32` could not | ✅ Closed: `ans:3`'s four records are re-based to the 0.11 image and the 1e-6 run-to-run control is registered with a negative control that bites. The audit re-tiered the row to heavy with no demotion. *History: `docs/planning/chunks/OPS-33.md`.* | heavy (re-tiered from standard 2026-09-03, 206 s measured; 165 s + 206 s + 1 s) |
| `OPS-34` | Measure the `ports:1` terminated-`Z` records against the 0.11 image; re-base the example's four `RECORDED_*` only if stale | ✅ Closed: the `ports:1` terminated-`Z` records were measured against the 0.11 image under the 1e-6 band. The stale digits in `test_port_package_sparameters.py` and the `gap_voltage.py` docstring are left to a review. *History: `docs/planning/chunks/OPS-34.md`.* | standard (143 / 142 / 143 s) |
| `OPS-31` | Re-record the `ports:3` cross-route narrative to the 0.11 image | ✅ Closed: the `ports:3` cross-route narrative is re-recorded to the 0.11 image, and the audit re-tiered the row to heavy with no demotion. *History: `docs/planning/chunks/OPS-31.md`.* | heavy (re-tiered from standard 2026-09-02, 235 s measured) |

**`OPS-29` — rank-safe the `phantom_material` empty-tag check in
`build_material_fields`** ✅ *(2026-08-28; full 122-line entry archived in
`docs/planning/plan-archive.md`, 2026-09-06 weekly review).* Result: the
`OPS-13` defect survived 20 lines below its own fix; measured breaking
`examples/mri/01_coil_phantom_fields.py` at `-n 12` (red baseline), green
at `-n 12` and `-n 2` after the reduction, example green at `-n 12`. Live
carry-forward: none. Logs: `20260828T165319Z_OPS-29-red-baseline.log`,
`20260828T165646Z_OPS-29-green-n12.log`,
`20260828T165709Z_OPS-29-example-n12.log`.

**`OPS-24` — migrate `core/cavity.py` to 0.11; turn `TH-9`'s gates back on** ✅
*(commissioned 2026-08-24 18:00 review from `EX-30` leg (th)'s finding 1;
closed 2026-08-25, 21:00 CDT slot.)*
> **Closed 2026-08-25.** The break was a pure keyword rename, established by
> introspecting the installed 0.11 `dolfinx.fem.petsc.assemble_matrix` rather
> than assumed: `diagonal=` → **`diag=`**, docstring "Rows/columns that are
> constrained by a Dirichlet boundary condition are zeroed, with the diagonal
> to set to `diag`" — semantics unchanged, so the constrained-DOF eigenvalues
> still land at `bc_diagonal`/1.0 and `solve_pec_cavity_modes`'s
> `spurious_cutoff = 0.5 · bc_diagonal` reasoning holds verbatim. Two lines
> changed in `core/cavity.py` (`:129`, `:131`) plus a migration comment; no
> test, band, tolerance or recorded eigenfrequency touched, and no solver
> path altered. **Red baseline reproduced in-slot first** — `4 failed, 9
> passed in 1.83s`, all 9 `tests/environment` green, Status 1
> (`20260825T020052Z_OPS-24-red-baseline.log`, 4 s) — matching the
> commissioning probe's `4 failed, 9 passed in 2.11s` exactly. **After the
> fix, `13 passed` twice**: 32.11 s / Status 0 / 33 s harness
> (`20260825T020111Z_OPS-24-green.log`) and 29.71 s / Status 0 / 31 s
> harness with `-s` for the printed diagnostics
> (`20260825T020157Z_OPS-24-green-quoted.log`), both `-n 2`, complex,
> `FEM_EM_REQUIRE_COMPLEX=1`, `tests/environment` first. **Every recorded
> figure reproduces the pre-0.11 record to the printed digit** — the
> closed-form eigenfrequency comparison on 720 cells / 5330 dofs prints
> 239.9805 / 291.3904 / 312.3465 / 346.5469 MHz against analytic
> 239.9510 / 291.3459 / 312.2838 / 346.3958 MHz, i.e. **0.0123 / 0.0153 /
> 0.0201 / 0.0436%**, worst-mode **0.0436%** equal to the `TH-9` record and
> to the test module's own header table; `null_modes_in_band = 0`;
> refinement h 0.1667 → 0.1143 takes max err 0.0436% → **0.0102%** at fitted
> rate **3.85** (gate > 2.0); the gradient-mode zero cluster is 8/8 below
> 2.529e-07 with max |λ| = **5.560e-14** against k₁² = 25.2909 (gate 1e-8
> relative); the energy-continuity guard fires **137.554** near-resonant and
> reads **21.951** clear against the 50.0 threshold. That digit-for-digit
> agreement is itself the evidence the rename was semantics-preserving.
> Retires the cavity known-issues entry and the §2.1 non-executing caveat;
> `th:2` / `th:5` are unblocked for `EX-30` leg (th) (§9 item 4), which owns
> re-running them — this chunk did not.
`cavity.py:129`/`:131` still pass `diagonal=` to
`dolfinx.fem.petsc.assemble_matrix`, dropped in 0.11 — `OPS-18` step 2's
migration missed the module because nothing scheduled runs
`test_cavity_resonances.py` / `test_resonance_guard.py`, so the 0.11 merge
landed with a dead subsystem inside a green-looking tree. Gate probe on
`main`: `4 failed, 9 passed in 2.11s`, all 9 `tests/environment` green
(`20260824T213908Z_EX-30-th-cavity-gate-probe.log`). Full rubric (anchors,
red-baseline negative control, traps, stop rules) in the §9 item. Done-when:
all four tests green at `-n 2` complex with the closed-form eigenfrequency
comparison quoted (record 0.0436% worst-mode), no band or recorded
eigenfrequency touched; retires the cavity known-issues entry and removes
the §2.1 non-executing caveat.

**`OPS-26` — systematic dolfinx-0.11 migration completeness sweep** ✅
*(step 1 ✅ 2026-08-25 static sweep — 434 call sites over 29 APIs clean,
`20260825T201054Z_OPS-26.log`; step 2 ✅ 2026-08-27 over seven census legs
(a)–(h), the 19:30 → 16:30 slots; audited COMPLIANT 2026-08-27 18:00 review.)*
**Gated result: 452 of 478 collected tests observed in a footered run on the
0.11 image — 436 green, 16 red in three owned families, 26 deferred with
measured reasons, zero `not reached in slot` and zero `deferred — no
footer`.** The 16 reds are ten stale 0.7.2-era exact records over 8 modules /
5 meshes (`OPS-27`, closed), three gmsh "overlapping facets" reds (`GEO-23`),
three test-double-drift names (`OPS-28`, closed). Chunk-level reconciliation:
`docs/testing/attempts.md` 2026-08-27T22:20Z. **Live carry-forwards:**
`GEO-23`; finding 36's permanently deferred `reactance_box_truncation` and
the sixth mesh suspected behind its fixture (a `MAT-6` fixture-pricing
question, still open). The table row above carries the per-leg evidence.
**Full narrative archived 2026-09-09 → `docs/planning/plan-archive.md`.**

**`OPS-27` — re-record the stale 0.7.2-era exact records, version-tagged;
sweep for siblings** ✅ *(step 1 2026-08-27 19:30 slot, step 2 2026-08-27
21:00 slot, step 3 2026-08-28 16:30 slot; commissioned 2026-08-27 10:30
review from `OPS-26` step 2 findings 19 and 23.)* **Re-recorded — each an
exact equality version-tagged with its 0.7.2 digit and its census log
in-comment, no band introduced or moved, `git diff -- src/` empty:**
`RECORD_128_RELL2` 0.01826 → 0.017686 and `RECORD_128_SEPARATION` 57.31 →
59.16; `NCELLS_BASELINE` 138 619 → 138 490; `NCELLS_THIRD` 2 807 309 →
2 808 204; `NCELLS_FINE` 417 914 → 418 888; `NCELLS_COMBINED` 697 401 →
697 926; `wire_resolution:268` 366 207 → 365 970. Logs: step 1 four windows /
604 s, step 2 three windows / 1 523 s, step 3
`20260828T213049Z_OPS-27-step3-larmor-resolution.log` +
`20260828T213807Z_OPS-27-step3-thirdrung.log` (719 s). **Live carry-forward:**
finding 45's coupled-constant residue belongs to whichever chunk re-prices
those fixtures. **Full narrative archived 2026-09-09 → `docs/planning/plan-archive.md`.**

**`OPS-28` — restore the `_DummyComm` double behind `OPS-14`'s reduction**
✅ *(closed 2026-08-28, 22:30 implementer slot)* *(commissioned 2026-08-27 10:30 review from `OPS-26` step 2 finding
12.)* `excitation.py:262-268` reduces `cell_tags.values` with
`comm.allgather` so the deprecated placeholder route "stays runnable"; the
only module that runs that route stubs its comm with `rank` + `allreduce`
and has been red since — undetected because nothing scheduled runs
`tests/ports`. Done-when: the double gains `allgather(value) -> [value]`,
the induced-voltage sign-flip test is green from `main`, the S-matrix
sign-flip test's outcome is *recorded* and known-issues entry 3 is either
retired for that name or re-dated with the reason it still fails, no
`src/` change. Full rubric is §9 item 3.
> **Closed 2026-08-28.** All four done-when clauses met, and the outcome
> of the S-matrix name is the branch the rubric anticipated: it is red at
> its own assertion, so entry 3 is **re-dated with the reason**, not
> retired. The finding worth keeping is that entry 3's one-line statement
> was wrong for this name — the placeholder's matched port on the 2-port
> orientation fixture is the *undriven* one, so the identically-zero wave
> lands on the **off-diagonal** (`S21 = S12 = 0`) while the diagonal is a
> healthy `9.047e-01 − 1.289e-02j`. Same `b = (V − Z₀I)/(2√Z₀) = 0`
> mechanism entry 3 diagnosed for the 3-port fake's diagonal, different
> matrix entry; the entry's title now says "power wave", not "diagonal".
> **Finding 44 (method, for the review):** this is the second time in two
> days that a test double, not the physics, produced a red — and the
> module is still in no scheduled command, so nothing prevents the next
> drift. The durable fix is coverage, not another one-line double repair;
> a review may want to price `tests/ports` (17 names, **2 s** at `-n 2`)
> into a scheduled command. Quantitative anchor, tier and elapsed times
> are in the §7 table row above.

**`OPS-25` — re-join `th:7` to its gate (hoist the series-interior
interpolation)** ✅ *(2026-08-25; full 50-line entry archived in
`docs/planning/plan-archive.md`, 2026-09-06 weekly review).* Result:
hoisted to `series_interior_function` in `test_lossy_sphere_fullwave.py`;
`th:7` green in 14 s with both element-order records reproducing (degree 1
8.1541% / 8.3869%, degree 2 0.1405% / 0.0058%), the gate's
`P_series(meshed)` bit-identical across the refactor, `13 passed`. Live
carry-forward: none. Logs: `20260825T033114Z_OPS-25-red-baseline.log`,
`20260825T033152Z_OPS-25-th7-green.log`,
`20260825T033221Z_OPS-25-gate-green.log`.

**`OPS-18` — DolfinX version upgrade, recurring** ✅ *(closed 2026-08-23;
glyph reconciled 2026-09-01; full 519-line entry archived in
`docs/planning/plan-archive.md`, 2026-09-06 weekly review).* Result:
`0.11.0.post0` builds and boots in both modes; `418 collected / 0 errors`
against the red baseline's `124 / 75`; the migration was one module
(`io/mesh.py`: `gmshio` → `dolfinx.io.gmsh`, 11 call sites, one shim). Step 3
re-recorded the mesher-dependent records on the 0.11 image under the (1*)
licence (`TH-10` 128 MHz 1.769%, one new-gmsh volume drift 4.251e-04 filed
and disposed); step 3b closed the chunk. Live carry-forward: none — the
recurring-upgrade trigger is the row's own "newest qualifying image" clause.
Logs: `20260822T093934Z_OPS-18-step1-real.log`,
`20260822T093943Z_OPS-18-step1-complex.log`,
`20260822T094005Z_OPS-18-step2-census-real.log` (11 in the archived entry).

**`OPS-19` — doc-reference checker: staleness must not own the exit code** ✅ *(commissioned 2026-08-16, 10:30 review; step 1 closed 2026-08-16, 16:30 implementer slot; archived 2026-08-23.)*
> **Result.** `scripts/testing/check_example_doc_references.py` now exits `EXIT_OK`/`EXIT_HARD`/`EXIT_STALE_ONLY` = 0/1/2, with `--stale-severity {fail,report}` defaulting to `report` (`fail` reproduces the pre-split reading bit for bit) and a machine-readable `RESULT: dead=… guide=… stale=… stale_severity=… exit=…` line; `--max-age-s` (`OPS-15`'s 48 h) unchanged, no example re-run or refreshed.
> **Gated:** `tests/unit/test_doc_reference_exit_codes.py`, 8 tests, 1.91 s, smoke, `-n 1` (`20260816T213312Z_OPS-19-step1-rerun.log`) — on the committed tree `dead=0 guide=0 stale=24 stale_severity=report exit=2`, guide pass green at **21/21 examples, 0 pending**; each fixture asserts the exit code against the literal and against arithmetic over the printed counts. **Negative controls:** a dead artifact reference and a missing `.py` both still exit 1 (`dead=1 stale=0`); boundary on the untouched default: 47 h → `stale=0 exit=0`, 49 h → `stale=1 exit=2`.
> **Bug fixed in passing:** `collect_references` called `doc.relative_to(REPO_ROOT)` unconditionally, so any `--docs-root` outside the repo raised `ValueError` (first run `20260816T213248Z_OPS-19-step1.log`, 7 failed / 1 passed, 2 s); now `display_path()`.
> **Carry-forwards:** the checker has exactly one call site class (ad hoc harness commands; `run_examples.sh` has none) — if a docrefs call is ever added to the runner, `0` and `2` are both pass. The stale `paraview_output/` artifacts remain a standing backlog, not this chunk's — and the figure is **55**, not 24: `EX-29` closed 2026-08-24 and the pre-fix count was a census of the 5 examples that write to the repo-root directory (`EX-30` refreshes the set).
> Full narrative: `docs/planning/plan-archive.md`, entry «§7 OPS-19 full narrative — archived 2026-08-23 (weekly review)».

**`OPS-16` — retry-on-529 in the automation launchers** ⬛ **WON'T FIX** *(commissioned 2026-08-13, 10:30 review; blocked 2026-08-14, 21:00 run; **declined by the human operator 2026-08-22, interactive session**; archived 2026-08-23. Do not re-commission, do not re-raise the permission ask, and do not queue it in §9.)*
> **Result.** Blocked by the permission layer, not the work: every file lives under `scripts/automation/`, and `Edit(scripts/automation/**)` is in the **`ask`** section of `.claude/settings.json`, which a headless `claude -p --permission-mode acceptEdits` run cannot answer — no scheduled session can execute this chunk in any scoping. **Operator decision 2026-08-22:** the rule **stays under `ask`** (a session that can edit its own launcher can change its own model, effort, timeout, `--permission-mode` and `--disallowedTools`), and the avoided cost — a few implementer slots lost to an out-of-credits, API-500 or 529 launch failure (worst observed: the 2026-08-19 18:00 review, `logs/automation/20260819T230001Z_daily-review.log`, one review plus effectively four implementer slots) — is **acceptable**. Settled, not open. The full design is recorded in the 2026-08-14T02:03Z `docs/testing/attempts.md` entry for reference only.
> **Standing instructions, in force until the operator says otherwise:** (1) do not queue `OPS-16` in §9, as a fallback or spare, or toward queue depth; (2) do not re-escalate the unblock on the dashboard's Waiting-on-you, in attempts.md, or in a review summary; (3) journal launch failures as ordinary weather — record the cost, never attribute it to a missing `OPS-16`.
> **Live carry-forwards:** `.gitignore:13` is a bare `lib/` (no leading slash), so any `*/lib/` in the repo is **ignored at any depth** — still true, still worth fixing if such a directory is ever added (use another name or a `!scripts/automation/lib/` negation). The *supervisor* framing (retry logic above the launcher, in an ungated directory called by cron with its own outer `timeout`) is recorded so it is not re-derived, **not** as a queued alternative.
> Full narrative: `docs/planning/plan-archive.md`, entry «§7 OPS-16 full narrative (won't fix) — archived 2026-08-23 (weekly review)».

**`OPS-17` — delete or replace the finiteness-only test suites** ✅ *(commissioned 2026-08-16, operator session; step 1 ✅ 2026-08-17, step 2 ✅ 2026-08-17, step 3 legs (a) 2026-08-17 / (b1) 2026-08-18 / (b2) closed by the 2026-08-21 18:00 review — **chunk ✅ 2026-08-21**; archived 2026-08-23.)*
> **Step 1 — sweep (`scripts/testing/finiteness_sweep.py`, AST):** `20260817T020244Z_OPS-17-step1-sweep.log`, smoke, 2 s, exit 0 — **306 test functions in 89 files: 225 `QUANT`, 22 `pytest.raises`-only, 59 candidates** (11 assert nothing); every candidate read; disposition **10 replace / 4 delete / 45 keep**. Finding: no `⚠️` chunk was propped up by a swept row.
> **Step 2 — dispositions (2026-08-17, 06:00 slot):** 4 deletes + 10 replacements landed, all `-n 2`: `20260817T111036Z_OPS-17-step2-collect.log` (359 collected, 6 s), `20260817T111054Z_OPS-17-step2-mesh-n2.log` (15 s), `20260817T111217Z_OPS-17-step2-solver-n2.log` (41 s), `20260817T112448Z_OPS-17-step2-th-smoke2-n2.log`, `20260817T113031Z_OPS-17-step2-portgap-n2.log` (1 passed, 448 s), `20260817T113806Z_OPS-17-step2-xfail-n2.log` (10 passed, 2 xfailed, 202 s). Anchors: `solver/test_cylinder.py` straight-wire `μ₀I/2πr` **13.2751%** L2 (band 25%); `solver/test_coil_phantom_magnetostatics.py` on-axis two-loop Biot–Savart **17.1233%** L2 (30%); `solver/test_two_torus.py` and `mesh/test_mesh_tag_integrity.py` volume-partition ratio **1.000000000000** (1e-9); `mesh/test_birdcage_port_tags.py` layout diagnostics exact (1e-12); `validation/test_straight_wire.py` fitted h-rate in `[0.7, 1.5]`; `validation/test_port_gap_voltage_impedance.py` 3b-x record pinned (1%). Four defects surfaced as `xfail(strict=True)`, **no band loosened**: region-resolution policy shrinks coil volumes **−21.68% / −22.62%** (→ `GEO-17`); Coulomb multiplier spread **7.836781e+00** on a divergence-free source (→ `MAG-17`); Poynting imbalance **116.7465%** vs 25% with wrong-sign flux, and `poynting_power_balance` raising on scalar `sigma=0.0` (→ `POST-5`). `⚠️`-retirement clause: confirmed, nothing to retire.
> **Step 3 leg (a) — real mode, closed 2026-08-17:** sweep control **56 candidates, reconciled** against the prescribed 45 (two `replace` rows kept finiteness bodies via a new sibling test; 10 newcomers keep-class; `20260817T200056Z_OPS-17-step3-sweep.log`); `tests/ --ignore=tests/validation` **3 failed, 134 passed, 32 skipped, 2 xfailed** / 218 s (`20260817T201248Z_OPS-17-step3-real-nonvalidation-n2.log`); real validation 33 + 1 + 5 = **39 passed, 167 skipped = 206** (`20260817T213419Z_OPS-17-step3b-real-validation-remainder.log`, `20260817T213843Z_OPS-17-step3b-real-mesh-cache.log`, `test_convergence.py` priced at 119.61 s); collect **377 = 171 + 206** exactly (`20260817T214141Z_OPS-17-step3b-collect-real-unpiped.log`). The 3 failures are named: two `_DummyComm` `allgather` regressions from `PORT-1` step 4 and known-issues entry 3's zero-diagonal.
> **Step 3 leg (b1) — complex non-validation, closed 2026-08-18:** `3 failed, 122 passed, 1 xfailed` / 392.76 s (`20260818T123045Z_OPS-17-step3d-complex-nonsolver.log`); `tests/solver` **46 passed, 2 xfailed** / 111.22 s warm (`20260818T141104Z_OPS-17-step3e-complex-solver-warm.log` — defect 3's Poynting xfail and `MAG-17`'s observed in a completed leg); `test_coil_phantom_magnetostatics` FAILED in its own completed log (`20260818T124712Z_OPS-17-step3d-coilphantom-complex.log`, `ComplexComparisonError` → `OPS-20`); collect 49 (`20260818T141312Z_OPS-17-step3e-collect-solver.log`); 126 + 45 = **171**, the same 171 as leg (a).
> **⚠️ Memory-premise caveat (2026-08-24).** The `coil_loading_degree2` deferral below rests on the `TH-12` degree-2 memory wall measured at the **old 64 GiB** ceiling; §5.1 now records **128 GiB**. If that file becomes runnable the **216 denominator moves** — which is a **review** decision and an explicit re-based close, never a silent edit. The deferral stands until then; `port_gap_voltage_padding` is unaffected (its reason is padded-record-only, not memory).
> **Step 3 leg (b2) — complex validation, closed 2026-08-21 18:00 review at 216 of 216:** nine attempts 2026-08-19 → 2026-08-21 across 25 exit-0 complex logs (`20260819T020055Z_OPS-17-step3f-complex-portgap-impedance.log`, 24 passed / 488.37 s, through `20260821T033534Z_OPS-17-step3j-sar-padding-group.log`, 14 passed / 247.68 s); every printed physics figure bit-identical to its `MAT-6`/`TH-11`/`PORT-10` record (e.g. `dR` 1.5763% / 1.0562% / +1.5834% / +5.5912%, `ΔR = +1.3838746e+00 Ω`, 417914 cells), only wall-clock moved (≤ +9.28%); denominator re-based 206 → 225 → 227 → **232**, audited `20260821T050352Z_OPS-17-step3l-collect-audit2.log` (236 = 4 + 232), observed file set 49 of 51 with the complement exactly the two deferred files. **Formally deferred:** `test_coil_loading_degree2.py` (14 — record `5 passed, 13 skipped`, the skips *are* the `TH-12` memory wall) and `test_port_gap_voltage_padding.py` (2, since attempt 3); attempting either is new-chunk work with a memory prescription.
> **Live carry-forwards:** the 3 `tests/ports/` failures remain named expected failures owned by `PORT-1` / known-issues entry 3; the 10 post-sweep candidate newcomers (`TH-11` step 4/5a ×8, `PORT-9` step 2b ×1, `OPS-17`'s deliberate `isnan` half ×1) await a review disposition; the `OPS-18` deferral commitment was executed 2026-08-21. **Standing rules learned here:** never pipe pytest inside a harness command (the footer reports the pipe's exit); complex ≈ 2.6–2.7× real on a warm FFCx cache, a family's first cold command ≈ ×3 (size as a throwaway warm-up); size each file from *all* its own recorded logs (width and elapsed), a padded record is an upper bound; rung/mode env vars select meshes, not test partitions — confirm splits against collect IDs; preflight `find /root/.cache/fenics -name '*.c' -size 0` and delete stubs only; coverage cannot be re-derived from node IDs at `-n 2` (footer arithmetic + file presence is the sound route); a non-collective complex raise hangs `mpiexec` ≈300 s on exit; `memory.peak` is pinned at `memory.max` and read-only — use `memory.current`.
> Full narrative: `docs/planning/plan-archive.md`, entry «§7 OPS-17 full narrative — archived 2026-08-23 (weekly review)».

**`OPS-20` — disposition the coil-phantom `ComplexComparisonError`** ✅ *(closed 2026-08-19, 06:00 implementer slot; disposition (a) — fixed, not marked; no `@real_only`, complex collect stays at 49)*.
- The commissioned `ComplexComparisonError` was dead on arrival: the test imports `azimuthal_current_density` from `tests/validation/test_circular_loop.py`, already repaired by `OPS-22`; no cold-cache window spent. The real defect was the second layer (`ValueError: Unknown format code '%' for object of type 'complex'` at `test_coil_phantom_magnetostatics.py:145`, rank-split because only rank 0 prints); fixed with the `OPS-22` idiom — assert `max|Im B_z| ≤ 1e-12·max|B_z|`, compare on `np.real`. The ≈300 s non-collective exit hang died with the raise.
- Gated numbers (`-n 2`, standard): real control before edit L2 **17.1233%** vs the 30% band (`OPS-17` step-2 record to the digit); complex after fix **1 passed / 5.11 s, L2 17.1233%**, both ranks identical; real re-run **17.1233%** unmoved. Stub sweep clean before and after.
- Logs: `20260819T110051Z` (real control, 7 s), `20260819T110111Z` (complex diagnosis, `--tb=long`, user frame at :145), `20260819T110144Z` (complex fixed), `20260819T110156Z` (real re-run); uncounted whole-`tests/solver` complex batch `20260819T110220Z` timed out at 89% (exit 124, 481 s) — no count claim.
- Carry-forwards: complex `tests/solver` no longer fits a 480 s window (111.22 s warm on 2026-08-18; `POST-5` step 2 is the candidate for the added cold forms); the two examples journaled by `OPS-22` (`02_circular_loop.py:173`, `04_helmholtz_analytic_comparison.py:79`) will carry this second layer too.

Full narrative: `docs/planning/plan-archive.md`, entry «§7 OPS-20 full narrative — archived 2026-08-23 (weekly review)».

**`OPS-21` — make the combined-XDMF test scalar-type-aware and
rank-deterministic** ✅ *(closed 2026-08-19, 16:30 implementer slot; test-side only, no writer change)*.
- The rank split was never a tmp-path race (the fixture has broadcast rank 0's path since `8c6ac03`): the mechanism was the test's own `if comm.rank != 0: return`, so non-zero ranks passed unconditionally and real-mode coverage was rank-0-only. Fix: rank 0 parses light data + every referenced heavy array (`_read_combined`), `comm.bcast`s, every rank runs every assertion. `SCALAR_IS_COMPLEX` selects `EXPECTED_NAMES`; the complementary spelling is `FORBIDDEN_NAMES`, asserted disjoint; imaginary parts asserted identically zero.
- Gated: exact set identity in both builds at `-n 2`. Real `{CellTags, F, G}` with the six split names absent — 1 passed / 3 s (`20260819T213140Z_OPS-21-step1-real.log`). Complex `{real_F, imag_F, real_G, imag_G, real_CellTags, imag_CellTags}` with the three bare names absent — 5 passed / 2 s (`20260819T213153Z_OPS-21-step1-complex.log`). Both ranks' summary lines identical in each run.
- Red baseline (predicate inverted): 1 failed on **both** ranks, byte-identical message, exit 1 / 2 s (`20260819T213221Z_OPS-21-step1-redbaseline.log`); reverted and re-confirmed green, 1 passed (`20260819T213234Z_OPS-21-step1-real-final.log`).
- Carry-forwards: `OPS-17` leg (b1) may count this file in complex; known-issues entry removed. Follow-up, not forced: the `G` field still has no value assertion in either build (presence + zero imaginary part only). The rank-0-return pattern was swept repo-wide by `OPS-23`.

Full narrative: `docs/planning/plan-archive.md`, entry «§7 OPS-21 full narrative — archived 2026-08-23 (weekly review)».

**`OPS-22` — make the three magnetostatic loop-drive fixtures
complex-safe** ✅ *(closed 2026-08-19, 04:30 implementer slot; fixed, not marked — all three files, no `@real_only` anywhere)*.
- Two defects: (i) the commissioned `ufl.max_value` / `<=` predicates in the fixtures' `current_density` callables — regularised inside the `sqrt` (`+ 1e-24`), wire predicates rewritten as `ufl.le(ufl.real(r²), a²)`; (ii) a second layer behind it — `evaluate_vector_field_parallel` returns the complex scalar type for a real magnetostatic solution (`ValueError: Unknown format code '%'`), fixed by asserting `max|Im B_z| ≤ 1e-12·max|B_z|` then comparing on `np.real` (a new complex-mode assertion, no-op in real mode).
- Gated (`-n 2`): real baseline before any edit 5 passed / 223.24 s — loop relL2 **7.0658%**, max **13.8212%**, |B_z|max 2.974560e-05 T; Helmholtz centre **0.728%**, mean **0.644%**, CV **0.1602%** (`20260819T093105Z`); after predicate fix (`20260819T093529Z`) and after real-part fix (`20260819T095414Z`, 5 passed / 199.91 s) every digit identical. Complex build, all three files (`20260819T094710Z`, `FEM_EM_REQUIRE_COMPLEX=1`): **5 passed / 412.12 s / exit 0**, both ranks identical, digits match the real record to the last figure. Interim `20260819T093933Z` (1 failed / 3 passed) is the between-fixes state. Stub sweep clean before and after.
- Costs for sizing: `test_circular_loop` is the sink in complex — **289.41 s** (on-axis) + **102.46 s** (symmetry); magnitude file 18.99 s; `_v2` 0.74 s.
- Carry-forwards: the same `max_value` idiom stands in `examples/magnetostatics/02_circular_loop.py:173` and `04_helmholtz_analytic_comparison.py:79`, unexercised in complex mode (expect layer (ii) there too); `OPS-17` leg (b2) may draw its 5 blocked tests; `test_helmholtz_v2.py`'s `float()` cast caveat was closed by `OPS-23`.

Full narrative: `docs/planning/plan-archive.md`, entry «§7 OPS-22 full narrative — archived 2026-08-23 (weekly review)».

**`OPS-13` — land the rank-safe `_validate_material_map_tags` fix** ✅
*(2026-08-08; full narrative in `docs/planning/plan-archive.md`)*. The one
hunk from the 3b-xiii branch: the tag set is reduced with
`mesh.comm.allgather` before it is tested. New gate
`tests/materials/test_material_map_rank_safety.py` — worst-case fixture
(exactly one tagged cell of 162), exact set identity + volume identity to
1e-12, every digit identical at `-n 2`/`-n 4`; red baseline reproduced the
3b-xiii hang (exit 124) with the hunk stashed. In CI at both widths.
Observed and left alone (worth a review's scoping):
`build_material_fields`'s `phantom_cells.size == 0` check reads the same
rank-local array — unmeasured, not known-broken.

**`OPS-14` — diagnose the rank-dependence of `test_single_port_excitation`
(known-issues 6)** ✅ *(closed as a diagnosis 2026-08-08; full narrative in
`docs/planning/plan-archive.md`)*. Two defects, both wholly inside the
`PORT-0` placeholder: the fixture tags over **rank-local** indices (the
global tag set is rank-count dependent — required tags genuinely absent at
`-n 4/8`), and `excitation.py` handed rank-local `cell_tags.values` to the
validator (non-collective raise — the 3b-xiii hang family). Separated by
counterfactuals; neither fix suffices alone. Pre-registered not-to-fix
branch taken: known-issues 6 re-pointed at `PORT-1` (whose step 4 later
fixed the validator defect); only a docstring hazard warning landed under
`src/`. Survey: no other non-collective tag read remains in `src/`.

**`OPS-11` — `tests/mesh` in CI** ✅ *(created 2026-08-02, closed 2026-08-03,
12:00 run; full narrative archived in `docs/planning/plan-archive.md`)*
> The directory no job ran — which is why known-issues 7 (mesh generators
> failing outright) sat undiscovered. The `validation` job's
> `Mesh generation suite` step now runs the **whole directory** with exactly
> one exclusion left: the known-issues-5 `--deselect` (off-centre sizing
> arithmetic) — **removed 2026-08-06 by `GEO-4` step 1; the step now excludes
> nothing and runs 27 passed 1 skipped in 85.3 s**. **20 passed 1 skipped
> 1 deselected in 42.15 s, exit 0**
> (`20260803T200504Z_GEO-9-step2b-gate.log`; the birdcage `--ignore` was
> removed by `GEO-9` step 2b as required). The "those and only those" control
> was executed, not quoted (`20260803T170132Z_OPS-11-fullsweep.log`); the
> §4.3 assertion is the volume-partition identities (`1e-9`) executing in CI
> across three files. The remaining exclusion is annotated at its
> known-issues entry and must be removed by the commit that fixes it.

**`OPS-10` — complex-mode CI job** ✅ *(2026-07-31; full narrative archived)*
> Before this, CI executed **no** time-harmonic solve at all (`@complex_only`
> skips). The `validation-complex` job sources
> `/usr/local/bin/dolfinx-complex-mode` and runs `tests/environment` first
> (so an environment regression is not blamed on the formulation), then the
> frequency-domain gates, under `FEM_EM_REQUIRE_COMPLEX=1` — which converts
> skips into failures so the job cannot pass by skipping. Verified with the
> CI-fidelity invocation (no `PYTHONPATH` override, `pip install -e`) plus a
> real-mode negative control that fails rather than skips.
>
> **This job has never executed on a GitHub runner.** Local `main` is well
> ahead of `origin/main` — nothing pushed since 2026-07-27 — so every "in CI"
> claim in this file is verified by local reproduction of the CI invocation
> only. The runner-environment caveat settles on the first push, which is a
> human action, not a scheduled-session one.
> CI notes for anyone editing `.github/workflows/ci.yml`: MPI here is
> **MPICH/Hydra**, so `--allow-run-as-root` is not a valid flag and will break the
> job. The `validation` job timeout is 45 min (`MAG-13` is heavy tier). Two tests
> are explicitly `--deselect`ed, both downstream of the time-harmonic proxy and to
> be revisited with `TH-1`: `test_phantom_material_assignment_and_time_harmonic_pipeline_wiring`
> and `test_phantom_field_metrics_and_exports_are_finite`.

**`OPS-12` — adjudicate the residual-trend classifier** ✅ *(2026-08-08;
retires known-issues 2; full narrative in
`docs/planning/plan-archive.md`)*. The code was wrong on all three counts,
not the test: undocumented asymmetric thresholds (labels now partition by
the sign of `f − 0.5`, table in the docstring); the second failure's
recorded symptom was wrong (an under-resourced `ksp_max_it` cap — the cap
moved, not the assertion); and the classifier was **unreachable in
production** (the time-harmonic path never called
`setConvergenceHistory()`, so the membership assertion had passed
vacuously) — now armed and tied to the unit identity. Gate: exact discrete
identity on an 11-row family spanning both threshold sides; in
`validation-complex`.

**Open follow-up in OPS — the `lint` CI job is red on `main`, and that is
adjudicated as expected-red for now** *(surfaced by `OPS-12`'s regression
log, decided 2026-08-08, 03:00 review)*. `flake8`/`black --check`/`isort
--check-only` fail on `src` and `tests` from **pre-existing** debt (W293
throughout `solvers.py`, E501 in `time_harmonic.py`, etc.) — no recent chunk
introduced any of it (verified for `OPS-12`: zero findings on its added
lines). A repo-wide reformat is deliberately **deferred until the `PORT-1`
lineage lands**: `attempt/PORT-1-step3bxiii-…` carries a 2000+-line test
file plus `src/` edits, and reformatting `main` underneath it would turn the
eventual landing into a conflict festival for zero behavioral gain. Whoever
lands that branch should scope the reformat as the next OPS chunk, in its
own commit, immediately after. Until then: red `lint` is known and expected;
do not "fix in passing", and do not read it as a chunk failure.

**`OPS-15` — retire the checker's standing freshness tax** ✅ *(2026-08-10;
full narrative in `docs/planning/plan-archive.md`)*. The doc-reference
checker's default `--max-age-s` went 3600 → **172800 (48 h)** — the 1 h
window was shorter than the 90-min slot grid and taxed every
example-touching slot an 80–200 s refresh. Two-sided anchor met: default
exits 0 on day-old artifacts with zero refresh solves; `--max-age-s 1`
still flags all 14 (the branch is retuned, not disabled — the 158-h
`EX-14` catch stays 3.3× over the new limit).

**`OPS-23` — sweep the `OPS-21` rank-0-return defect pattern + the
`test_helmholtz_v2.py` Im-bound** ✅ *(step 1 closed 2026-08-20, 09:00 implementer slot — and the chunk closes; test-side only, no `src/` change, no gate value moved)*.
- The commissioned census was wrong in both directions: `test_degree2_energy_mechanism.py:237` and `test_lossy_sphere_degree2.py:249` are guards inside a `_print_table` helper, not defects (left untouched, not run); the exempted `test_csv_export_stats_parity.py:252` is a real instance (rank-0-only negative control of `POST-1` step 6). Three real sites, all in `test_csv_export_stats_parity.py` (`:143`, `:192`, `:252`), fixed with the `OPS-21` template (rank 0 parses + `bcast`, every rank asserts, prints stay rank-0). `test_helmholtz_v2.py` gets `max|Im B_z| ≤ 1e-12·max|B_z|` before the `float()` casts, then explicit `np.real`.
- Gated (`-n 2`, smoke, both ranks identical every run): csv complex 11 passed / 5.51 s — 5 184 default rows / 4 896 guarded / **288** drops per tag, worst round-trip disagreement **3.808e-16** vs unmoved 1e-12 (`20260820T140248Z_OPS-23-step1-csv-green.log`); helmholtz real 2 passed / 3 skipped / 0.82 s (`20260820T140344Z_OPS-23-step1-helmholtz-real2.log`) and complex 5 passed / 1.09 s (`20260820T140330Z_OPS-23-step1-helmholtz-complex.log`), mean B_z **4.219228e-09 T**, CV **0.1873%** vs unmoved 1% gate, `max|Im B_z| = 0.000e+00` vs 4.231e-21 bound.
- Red baseline, all four fixed predicates inverted: 8 failed / 4 passed / exit 1 / 5.13 s, the eight `AssertionError` lines byte-identical between ranks (`20260820T140405Z_OPS-23-step1-redbaseline.log`); final 12 passed / exit 0 / 5.00 s (`20260820T140438Z_OPS-23-step1-final.log`).
- Standing nuance: round-off-scale nondeterminism in the iterative solves (round-trip 3.808e-16 vs 3.822e-16; helmholtz `std B_z` moves in the 6th significant digit, 7.902679 / 7.902639 / 7.902744e-12), four–six orders below the gates; gated digits bit-stable.

Full narrative: `docs/planning/plan-archive.md`, entry «§7 OPS-23 full narrative — archived 2026-08-23 (weekly review)».

### MAG — Magnetostatics (Phase 1)

| ID | Title | Status | Tier | Result |
|---|---|---|---|---|
| `MAG-1` | Vector-potential formulation, N1curl, gauge penalty | ✅ | standard | 0.04% centre vs closed form |
| `MAG-2` | Straight-wire analytic validation | ✅ | standard | |
| `MAG-3` | Circular-loop analytic validation | ✅ | standard | |
| `MAG-4` | Helmholtz analytic validation | ✅ | standard | 0.04% centre / 0.83% mean |
| `MAG-5` | h-refinement convergence study | ✅ | standard | |
| `MAG-6` | Coil+phantom B-field symmetry metric strategy | ✅ | landed on DG0 at h = 0.010: mirror-symmetry `max_rel_diff` **0.323844 / 0.302661 / 0.308407** at `-n 1/2/4` vs the untouched 0.350, three-way rank spread **7.00%** (gate ≤ 10%); `-n 1` byte-reproduces the on-record 0.323844 | steps 1–3 ✅ 2026-08-08 — CG1 interpolation owned it (rank-dependent 3.03×, non-convergent under `h`); boundary and gauge exonerated; ≈0.53 was discretisation (p ≈ 1.07); step 3 re-pointed **both** sampled metrics at DG0 and refined one rung, no tolerance touched; known-issues 4 retired. Gates discretisation symmetry, **not** phantom physics (uniform μ) |
| `MAG-7` | Fix point evaluation in validation tests | ✅ | standard | |
| `MAG-8` | Restrict straight-wire current density to the wire | ✅ | standard | |
| `MAG-9` | Re-size validation meshes to fit the tier budget | ✅ | standard | |
| `MAG-10` | Gauge penalty was below the safe window | ✅ | standard | default now 1.0 |
| `MAG-11` | Parallel energy was rank-local (missing allreduce) | ✅ | smoke | |
| `MAG-12` | `evaluate_at_points` used the MAG-7 broken pattern | ✅ | smoke | |
| `MAG-13` | Analytic-Dirichlet outer boundary for wire/loop | ✅ | heavy | wire 12.75%, loop 7.07%, rate 1.10; 167 s + 196 s. **2026-08-25: its rate gate (`test_h_refinement_straight_wire`, band [0.7, 1.5]) is red on `main` on the 0.11 image — fitted 1.9038, the finest rung's error collapsed 9.26% → 4.4605% (known-issues 2026-08-25). Disposition is `MAG-19`; the ✅ here is the 0.7.2 close and stands.** **Disposed 2026-08-25 by `MAG-19` step 2 (ruling (i)): the sampled two-sided band is retired with its basis, the rate duty now belongs to `MAG-18`'s one-sided `E_Ω` ≥ 0.7 gate, and `test_h_refinement_straight_wire` gates monotone decay and is green.** |
| `MAG-14` | Helmholtz magnitude comparison in the test suite | ✅ | smoke | 0.671 % vs closed form on the 0.11 image (`20261007T112319Z_MAG-14-rerecord-W1.log`, 2026-10-07); 0.728 % on the pre-0.11 image (`20260819T093105Z`, 2026-08-19); 1.731 % before `GEO-8`; 11 s, in CI |
| `MAG-15` | Lagrange-multiplier Coulomb gauge (cross-check) | ✅ | smoke | 7 passed, 13 s |
| `MAG-16` | Complex-build-safe magnetostatic energy | ✅ 2026-08-05 | smoke | 10 passed complex `-n 2` in 4.9 s; cross-build pin 2.9e-07, `Im W` exactly 0; retires known-issues 8 |
| `MAG-17` | Coulomb-gauge multiplier does not vanish for a divergence-free source: h-ladder discriminator (`OPS-17` step-2 defect 2, known-issues 2026-08-17; commissioned 2026-08-17 10:30 review) | ✅ *(audited COMPLIANT 2026-08-21 18:00 review — rate 2.4476 vs the pre-registered ≥ 0.7 verified in `…final2.log`; one nuance on record: the cited "ladder" log is exit 1 — a sign-convention fix in the fit sits between it and the record run, same spreads, band unmoved)* | standard |
| `MAG-18` | Sampler-independent straight-wire gate: annulus-restricted domain L2 of `\|B_h\| − \|B_ana\|` with a pre-registered rate band (`OPS-18` step 3 attempt 5 finding, … | ✅ Closed: the annulus-restricted domain L2 gate is sampler-independent, with the `-n 2` / `-n 4` readings inside the re-registered band and the natural-BC control strictly worse. No record moved and no band moved. *History: `docs/planning/chunks/MAG-18.md`.* | heavy |
| `MAG-19` | Dispose of the red straight-wire rate gate on 0.11 (fitted 1.9038 vs [0.7, 1.5]; the finest rung's sampled error collapsed on the image): anomalous rung vs … | ✅ Closed: the red straight-wire rate gate on 0.11 is disposed of, with the retired band's basis stated in-comment and `RATE_MIN` / `RATE_MAX` kept report-only. The residual is commissioned as `MAG-20`. *History: `docs/planning/chunks/MAG-19.md`.* | standard |
| `MAG-20` | Dispose of the residual two-sided sampled rate band in `test_straight_wire_convergence` — a two-rung 8-point sampled fit still gated on [0.7, 1.5], green at 0.7900, on an instrument `MAG-19` measured swinging 34% under its own sampler (commissioned 2026-08-26 03:00 review from `MAG-19` step 2's filed residual; measure-first, own decision rule — ruling (i) is **not** inherited) | ✅ **2026-08-28** *(step 1 measured and disposed in the 00:00 slot: the pre-stated sweep returns **no crossing** — fitted two-rung rate **0.7900 / 0.7246 / 0.9934** at n_points 8 / 10 / 20, all inside [0.7, 1.5] — so the band is **kept and validated**, not retired. Nothing moved: `RATE_MIN`/`RATE_MAX` unchanged, no assertion added or removed, the disposition is a ≈25-line in-comment measurement record. Probe `tests/validation/probe_straight_wire_convergence_npoints.py` (asserts nothing, 49 s at `-n 2`), anchor `test_straight_wire.py` **7 passed / 371 s / Status 0** from `main` with `E_Ω` fit 1.6854 and the h = 0.0025 record 1.0617170193e-01 (1.5e-9 of its 1e-4 band) untouched. **Two residual findings handed to the review, not acted on**: the sampler swing on this test's 0.4 R window is only 6–7% of the error (vs the 34% `MAG-19` measured on the 0.8 R window), but it still moves the *rate* by **37% of its own value**, and the n = 10 row clears `RATE_MIN` by **0.0246**)* | standard |

**`MAG-17` — the Coulomb-gauge multiplier does not vanish for a
divergence-free source: h-ladder discriminator** ✅ *(step 1 closed 2026-08-20, 07:30 implementer slot — and the chunk closes; audited COMPLIANT 2026-08-21 18:00 review)*.
- Verdict **DISCRETE-SOURCE**: multiplier spread 7.836781e+00 → 3.052022e+00 → 1.438617e+00 at h = 0.005 / 0.0035 / 0.0025 (29 190 / 82 819 / 208 049 cells), fitted log-log rate **2.4476** (pairwise 2.645 / 2.234) vs the pre-registered ≥ 0.7 band; the ASSEMBLY-DEFECT band (|rate| < 0.3) missed superlinearly. Base rung reproduces the `OPS-17` record to every printed digit. Negative control: the incompatible straight wire stays at 2.083064e+02, > 10× the loop's base-h spread (recorded 26.6×).
- Diagnosis: the anchor was wrong, not the constraint block — `p` absorbs the interpolated `J`'s O(h) discrete divergence, so "spread → 0 to solver tolerance" cannot hold on any single mesh. The strict xfail is retired; the claim now lives in `tests/solver/test_gauge_multiplier_convergence.py` (monotone decrease + rate ≥ the **unmoved** 0.7, deliberately not tightened; plus `test_multiplier_still_separates_an_incompatible_source`); `test_gauge_lagrange.py` keeps the wire-side scale gate.
- Standard tier, 97 s at `-n 2`, real build. Logs: `20260820T123307Z_MAG-17-step1-ladder.log` (ladder; exit 1 — a sign-convention fix in the fit sits between it and the record run, same spreads, band unmoved), `20260820T123823Z_MAG-17-step1-final2.log` (final, 6 passed), sizing probe `20260820T123124Z_MAG-17-step1-probe.log`.
- Carry-forwards: known-issues defect 2 retired in the same commit; none owed — the residual is benign; `MAG-15`'s open follow-ups are unaffected.

Full narrative: `docs/planning/plan-archive.md`, entry «§7 MAG-17 full narrative — archived 2026-08-23 (weekly review)».

**`MAG-16` — complex-build-safe magnetostatic energy** ✅ *(2026-08-05;
retires known-issues 8; full narrative in `docs/planning/plan-archive.md`)*.
`compute_magnetic_energy()` reduces with `np.real` and **raises** above
`ENERGY_IMAG_RTOL = 1e-8` (`abs()` rejected deliberately — it would absorb
a negative real part too). Cross-build pin 2.9e-07 (penalty; its operator
carries the gauge null space, hence `PIN_RTOL = 1e-5`) / 1.3e-13
(Lagrange); `Im W` exactly 0 by `ufl.inner` conjugation. The file joined
the `validation-complex` CI job, which is what stops the cast returning.

**`MAG-6` — the step record, compressed** *(steps 1–5; full plans,
journals and the step-4 √3 diagnosis archived in
`docs/planning/plan-archive.md`)*:
> * **Step 1** (2026-08-08): the boundary-mirror hypothesis died — the
>   symmetry metric is **rank-dependent 3.03×** through the CG1
>   interpolation (`curl A` is cell-wise constant; nodal averaging is
>   partition-dependent), while DG0 sampling is rank-stable; padding and
>   gauge both exonerated on the DG0 path. The negative control is not
>   directional because μ is uniform — the phantom is physically invisible.
> * **Step 2**: the DG0 metric converges at p ≈ 1.07 and meets the
>   untouched 0.350 at h = 0.010; CG1 does not converge on the same ladder.
>   The finest rung voided itself by its own rank-spread control (point
>   sampling under a fixed grid becomes a partition question again — the
>   ceiling on refining this estimator).
> * **Step 3** (closes the chunk, known-issues 4 retired): the test
>   samples DG0 at h = 0.010, both metrics re-pointed, tolerances
>   untouched; mirror `max_rel_diff` 0.324/0.303/0.308 at `-n 1/2/4`,
>   spread 7.00% vs the ≤ 10% gate. The centerline metric carries **no**
>   rank-stability claim.
> * **Step 4**: the centerline's "88% scatter" was two things — a **√3
>   probe bug** (`Function.eval` squeezes a single point to shape `(3,)`;
>   always `.reshape(-1, 3)` — `evaluate_vector_field_parallel` is immune
>   by construction) and **gauge contamination** at the sub-floor
>   `gauge_penalty=1e-3`; at the validated 1.0 the spread is 0.341%.
> * **Step 5**: the gate fixture solves at the validated floor
>   (`1e-3 → 1.0`); both metrics landed on their step-4 predictions to
>   < 0.01%. Of the other sub-floor call sites, only `examples/mri/01`
>   carried on-record numbers — handled by `EX-13`/`EX-16`.
>
> The metric gates **discretisation symmetry, not phantom physics**
> (uniform μ); the caveat lives in the test's module docstring.

**`MAG-13` step-2 lineage — the < 5% wire, measured three ways**
*(profile, re-gate, 2b, 2c, rungs 2–3, diag; full narratives archived in
`docs/planning/plan-archive.md`. `MAG-13` stays ✅ at its recorded numbers
throughout — this lineage extends measurements, it never reopened the
chunk.)*
> * **Rung 2** (2026-08-11): h = 0.00125, 1 097 873 cells, relL2 **5.6494%**
>   — on-rate (1.174) but over the 5% target. **Rung 3** (2026-08-13):
>   h = 0.001127, 1 520 152 cells, **3.7372%** in 420.3 s at `-n 8` — the
>   target reached by brute force; three-rung fitted rate 1.407, not a
>   converged constant. The per-radius "near-wire concentration" pattern
>   did not survive rung 3 (far radii worsened while the total fell —
>   mesh-realization noise, not a spatial map).
> * **Profile step** (2026-08-11/12, re-gated with biting exit codes): the
>   dense 45-radius map — error ∝ 1/r (slope −1.069) *and* a per-cell
>   staircase, both signatures of **cell-wise-constant B** (lowest-order
>   N1curl); local error is O(h/r), which is why the global rate reads ≈1.
> * **Step 2b** (2026-08-12): CG1-projected `curl A` reads **1.9557%**
>   where DG1 reads 4.7235% — on the *existing* mesh, for 2.71 s (1% of
>   the solve); the staircase breaks 8/8; the residual is band-flat ≈ 2%.
>   **Step 2c** (2026-08-13): the third rung confirms the recovery rate —
>   CG1 **p = 2.003** / DG1 p = 1.217 — with the honest caveat that
>   pairwise rates spread ±0.20 and the level came in 6.8% below the p = 2
>   extrapolation (floor-approach signature).
> * **Standing decision (weekly review's):** whether `compute_b_field`
>   moves to continuous CG1 recovery — cheap and second-order, but every
>   B-consuming test's recorded number would shift (a re-gating exercise).
>   Graded refinement stays the cheaper *mesh* route on cost alone. The
>   exit gate was bitten live 2026-08-15 (exit 1 at the 12.7485% smoke
>   rung, azimuthality PASS independently). Degree-2 `A` is on record
>   diverging on this fixture — not a free swap.
>   **Decided 2026-08-16 (weekly review): declined for now.** No §10
>   goal is currently limited by B-field accuracy — the wire gate stands
>   at 3.74% on the mesh route — so the re-gating tax (every recorded
>   B-consuming number shifts in one commit) buys nothing today.
>   Revisit when §10 subgoal-4 B1+ work opens: its gates are new, so
>   adopting CG1 recovery there carries no re-gating cost, and the
>   measured case (1.9557% vs 4.7235%, p = 2.003) is already on record.
>   Not an epitaph — the option stays live, the default stays DG1.
> * The two 2026-08-08 "mystery harness deaths" were the
>   **background-and-end-turn trap** (a headless session backgrounding a
>   harness run and ending its turn SIGKILLs the tree) — the reason every
>   scheduled-slot recipe runs harness commands foreground.

**`MAG-18` — sampler-independent straight-wire gate** ✅ *(closed 2026-08-23
03:00 review on the 2026-08-22 logs; commissioned 2026-08-22 18:00 review as
the disposal of `OPS-18` step 3 leg 2; re-gated on the 0.11 image
2026-08-23.)* Replaces the 10-point radial sampler — whose statistic swung
34% under its own `n_points` and whose 15% band already failed at 8 points on
0.7.2 — with an annulus-restricted domain L2, `E_Ω`. **Gated on 0.11, all
three anchors green, `7 passed` / exit 0 twice in-slot:** (i) `E_Ω` 25.2868 →
10.6172 → 6.6458% monotone at fitted rate 1.6854 ≥ 0.7; (ii) `-n 2` vs `-n 4`
agree to 4.86e-07 against the re-registered ≤ 1e-6 (the ≈1e-7 direct-LU
cross-width floor is why the original 1e-10 clause was not met); (iii)
natural-BC wall ratio 0.3285, strictly worse. `E_OMEGA_H0025_RECORD`
1.061717e-01 (v0.11.0) reproduces to 2.9e-09 of its 1e-4 band. Logs
`20260824T003059Z_MAG-18-regate-run1.log`, `…003650Z_…run2.log`,
`…003606Z_MAG-18-regate-n4.log`. **Full narrative archived 2026-09-09 → `docs/planning/plan-archive.md`.**

**`MAG-19` — dispose of the red straight-wire rate gate on 0.11** ✅
*(commissioned 2026-08-25 10:30 review, closed 2026-08-25; full 180-line
entry archived in `docs/planning/plan-archive.md`, 2026-09-06 weekly
review).* Result: both norms on the same four-rung ladder — the sampled
three-rung fit reproduces the red **1.9038** and the `E_Ω` fit reads
**1.6854** with the h = 0.0025 record at 2.094e-08 relative, so the `ANS-1`
import is right and the physics moved with the image; the pre-stated
decision rule selected neither branch (a second out-of-band sampled pair
at 0.5822; dropping h = 0.0018 alone returns 0.7309), and the ruling
re-registered the gate on the measured ladder without loosening. The MAG-19
audit narrative is the auditor agent's worked exemplar (`.claude/agents/
auditor.md`). Live carry-forward: none. Log:
`20260825T183555Z_MAG-19-step1-dualnorm-fits.log`.

**`MAG-20` — dispose of the residual two-sided sampled rate band in `test_straight_wire_convergence`** ✅
**2026-08-28** *(step 1 measured and disposed in the 00:00 slot; audited COMPLIANT 03:00 review; commissioned
2026-08-26 03:00 review from `MAG-19` step 2's residual. Full narrative archived in `docs/planning/plan-archive.md`, entry «§7 MAG-20 full narrative — archived 2026-08-30 (weekly review)».)*
> * **Measured (`20260828T050130Z_MAG-20-step1-npoints-probe.log`, 49 s, 0.11 / gmsh 4.15.2):** the two-rung
>   ([0.004, 0.0025], 38 740 / 147 235 cells) sampled fit at `n_points` 8 / 10 / 20 reads **0.7900 / 0.7246 / 0.9934**
>   (errors 21.5512 / 21.1826 / 22.6647% and 14.8669 / 15.0685 / 14.2097%, swing +7.00% / +6.04%). **No count crosses
>   either edge of [0.7, 1.5]** ⇒ under the pre-stated rule the band is kept and **validated**; `RATE_MIN`/`RATE_MAX`
>   unchanged, disposition is an in-comment measurement record; n = 8 reproduces `MAG-19`'s 0.7900 to four decimals.
> * **Anchor (§4):** `test_straight_wire.py` **7 passed / 369.95 s / Status 0** from `main`
>   (`20260828T050256Z_MAG-20-step1-anchor-module.log`), `E_Ω` fit **1.6854**, h = 0.0025 record 1.0617170193e-01 vs
>   tagged 1.0617170177e-01 (1.5e-09 relative, 1.5e-05 of its 1e-4 band); diff two pure-addition hunks, no `src/`.
> * **Findings 45–46, recorded not acted (03:00 ruling):** sampler swing on this window is 6–7% (vs `MAG-19`'s 34%)
>   yet moves the rate by 37% of its value; n = 10 clears `RATE_MIN` by only 0.0246. No third rung commissioned; if
>   the gate ever reads red at n = 8 the disposition is ruling (i)'s retire-with-basis; never tighten the record band below ≈1e-8.

**Open follow-ups in MAG:**

- `MAG-15` is a working option, not a finished subsystem: Dirichlet
  conditions on `A` are rejected (`bc_functions` raises — so
  `TimeHarmonicBoundaryCondition.PEC` would not work); the point-pin on `p`
  is not `H¹`-stable in 3D, so use `gauge_multiplier_spread()` rather than
  `max|p|`; it is not wired through `TimeHarmonicSolver` or the port entry
  points; and the degree-2 cost (≈7.5× the penalty) is unprofiled against
  MUMPS.
- `J·n ≠ 0` at the straight-wire end caps still stands unmeasured; capping
  the wire short of the end faces was never needed.

### GEO — Geometry & meshing

Independent of the §2.1 physics defect; meshes are meshes.

| ID | Title | Status | Tier |
|---|---|---|---|
| `GEO-1` | Parametric birdcage geometry generator *(2026-09-30 03:00 review: `GEO-1` is `birdcage_port_domain`'s default solid-conductor, box-port mode — the live generator, not a superseded one. Its module `tests/mesh/test_birdcage_port_tags.py` already gates the layout closed forms (1e-12), the tag set, the box tiling and the port-box volumes (1e-9); its one mesh test, `test_birdcage_volumes_partition_the_box`, is the open 🔴 of known-issues 2026-08-27 (`GEO-21`'s residual).)* **Done-when (amended 2026-10-02 18:00 review, rule (j); the 09-30 weekly's target; the review extended the 08-26 ruling (b) to this last call site, the weekly ratifies):** `test_birdcage_port_tags.py` green in one harness window, real, `-n 2` (plus `-n 1`), with `test_birdcage_volumes_partition_the_box` on the graded control h_c = 4.8e-3 (`GEO-21` ruling (b), version-tagged), its 1e-9 identities and region bands unchanged, the ungraded build asserted to raise on every rank, `test_birdcage_conductor_sizing.py` green in the same window, elapsed recorded ⇒ ✅; known-issues 2026-08-27 retires, 2026-08-25 stays open. *History: `docs/planning/chunks/GEO-1.md`.* | ✅ **Audit PASS 2026-10-03 03:00 review (chunk file).** 2026-10-03 (§9 item 186): `test_birdcage_port_tags.py` is green at `-n 2` (with `test_birdcage_conductor_sizing.py`, 14 passed) and `-n 1`. The partition test on h_c = 4.8e-3 reads 33 185 cells, identities 1.000000000000, conductor 0.8104, phantom 0.9705, and the ungraded build raises on every rank (`20261003T035654Z_GEO-1-W1.log`, `…035842Z_GEO-1-W2.log`). Known-issues 2026-08-27 retired; 2026-08-25 stays open. | standard |
| `GEO-2` | Port-face geometry robustness checks. **Done-when (2026-10-02 18:00 review, rule (j); §9 item 184):** the layout closed forms (`test_birdcage_port_tags.py:86`, 1e-12) and rejections (`:261`) green; `io.mesh_qa.cad_identity_report` reproduces the four port-box volumes 8.0e-7 m³ to 1e-9 on `birdcage_port_domain` at h_c = 4.8e-3 and flags exactly `exact_volume` when one P1 cell is re-tagged air; one harness window, real, `-n 2`, elapsed recorded ⇒ ✅. Not claimed: meshed port-face areas, sheet ports. *History: `docs/planning/chunks/GEO-7.md`.* | ✅ **Audit PASS 2026-10-03 03:00 review (chunk file `GEO-7.md`).** 2026-10-03 (§9 item 184): the layout closed forms and rejections are green. `cad_identity_report` reproduces the four port boxes at 8.0e-7 m³ (≤ 4.4e-16) and flags exactly `exact_volume` on the re-tag control (4.2e-2), at `-n 2` and `-n 1` (`20261003T005340Z_GEO-7-W1.log`, `…005520Z_GEO-7-W2.log`). | standard |
| `GEO-3` | Phantom placement presets (centered/off-center) — **Done-when (written 2026-09-29 03:00 review, lane 1, §9 item 73):** on `coil_phantom_domain` at h = 0.015, tag 3's volume-weighted first moment equals the placed cylinder's centre for `centered`, `off_center` and one explicit `phantom_offset_xy` override, each inside the in-run bound δ·r (in-plane) / δ·h/2 (axial) with δ = V_cyl/V_tag − 1 measured in the same run (M ⊂ C convex ⇒ m_M − c_C = δ(c_C − c_D)); premise 0 < δ < 0.1; the comparand-swap control misses by more than the bound; harness at `-n 2`, elapsed recorded ⇒ ✅ | ✅ **2026-09-29 (06:00 slot, §9 item 73):** δ = 0.016744 on all three meshes (V_C 5.026548e-4, V_M 4.943768e-4). Moment misses over their bounds are x 0.0284, y 0.0097, z 0.0038 on each placement: the same shape moved. The in-plane bound is asserted per component; the xy-norm reading, ≈ 0.030·δr, is computed from the logged components and is not gated. The comparand-swap control misses by 20.87× the bound (pred. ≈ 20.9×). 12 passed + 4 complex-skips, 50 s, `-n 2` (`20260929T112525Z_GEO-3-W1.log`). | standard |
| `GEO-4` | Air-box and boundary sizing heuristics. **Done-when (2026-10-02 18:00 review, rule (j); §9 item 185; the 09-30 weekly's target):** `MeshGenerator.two_torus_box_padding(τ)` inverts `PORT-1`'s free-exponent fit (reproduced bitwise from `port1_padding_fit.py`; D∞ +1.6934 pp, p 1.6574), rounded up to 1 cm, returning 0.11 m at τ = 6 pp and 0.14 m at τ = 4 pp; on the 3b-xi two-torus fixture the measured \|deficit(W) − D∞\| ≤ τ at both paddings; the 0.08 m rung reproduces 3b-xi within 0.05 pp and its box term exceeds both τ; `test_domain_sizing_heuristics.py` green in the same window; complex, `-n 2`, harness, elapsed recorded ⇒ ✅. Not claimed: a converged box-free deficit, any other fixture. *History: `docs/planning/chunks/GEO-4.md`.* | ✅ **Audit PASS 2026-10-03 03:00 review (chunk file).** 2026-10-03 (§9 item 185): `two_torus_box_padding` reproduces the probe's fit bitwise and returns 0.11 / 0.14 m at τ = 6 / 4 pp. The measured \|deficit − D∞\| is 5.719 / 3.901 pp (≤ τ), P0 is −8.0404 pp (0.008 pp from 3b-xi), and C1 / C2 fire (`20261003T035018Z_GEO-4-W1.log`, 29 passed; `…035410Z_GEO-4-W2.log`). The 0.14 m deficit −2.2071 pp is a fourth-rung record for `PORT-1` 3b-xi; the p = 3 fit would miss it by 0.46 pp. | heavy |
| `GEO-5` | Region-specific mesh resolution policy. Done-when (written 2026-09-30 03:00 review, rule (j); queued §9 lane 1). On `coil_phantom_domain` (the `GEO-17` … | ✅ Done 2026-09-30 (§9 item 103) with audit PASS: the phantom volume deficit falls at LSQ rate 1.9707 inside [1.7, 2.3], the equivalent-edge and partition gates hold on all four rungs, and the coil-slope control misses as predicted. Only the phantom leg is gated (coil and air stay on `GEO-17`'s sign gate), and a finer rung would need the a_eq/h_p band re-examined first. **Done-when:** every mesh's tags partition it to 1e-9; tag 3's volume never exceeds πr²H and its relative deficit falls strictly, with LSQ rate against the requested size in [1.7, 2.3]; tag 3's equivalent edge (6√2·V₃/N₃)^{1/3} lies in [0.8, 1.5]·h_p at every rung; the coil's deficit regressed on the same ladder misses the rate band; `test_mesh_tag_integrity.py` is green in the same window. Harness, real, `-n 2`, elapsed recorded ⇒ ✅. *History: `docs/planning/chunks/GEO-5.md`.* | standard |
| `GEO-6` | Geometry sanity report utility. **Done-when (2026-10-02 18:00 review, rule (j); §9 item 184):** `cad_identity_report` reproduces `coil_phantom_domain`'s `UNIFORM_VOLUMES_RECORD` (`test_mesh_tag_integrity.py:89–94`), box and partition identities at 1e-9, and `two_torus_domain`'s box and partition identities (`test_two_torus_conforming.py:97, :104`); its readings equal `tests/mesh/helpers.py`'s bitwise; an ungrouped cell flags exactly `partition`; harness, `-n 2` and `-n 1` ⇒ ✅. Not claimed: the legacy cell-ratio `coil_phantom_geometry_sanity_report` stays a print. *History: `docs/planning/chunks/GEO-7.md`.* | ✅ **Audit PASS 2026-10-03 03:00 review, tier relabelled standard (chunk file `GEO-7.md`).** 2026-10-03 (§9 item 184): `UNIFORM_VOLUMES_RECORD` is reproduced ≤ 1.7e-10. The box and partition identities hold ≤ 1.6e-14 on both generators, the readings are bitwise equal to `helpers.py`'s, and the ungrouped-cell control flags exactly `partition`, at `-n 2` and `-n 1` (`20261003T005340Z_GEO-7-W1.log`, `…005520Z_GEO-7-W2.log`). | standard (91 s at `-n 2`, 73 s at `-n 1`) |
| `GEO-7` | Mesh-tag QA diagnostic hardening. **Done-when (2026-10-02 18:00 review, rule (j); §9 item 184):** `io/mesh_qa.py` gains `cad_identity_report` (tag sets from owned cells only, reduced volumes, default band `VOLUME_PARTITION_BAND`), passing on all three live generators; a merged tag flags exactly `tag_set`; `tests/io/test_mesh_qa_diagnostics.py` and `test_mesh_tag_integrity.py` green in the same window; harness, `-n 2` ⇒ ✅. *History: `docs/planning/chunks/GEO-7.md`.* | ✅ **Audit PASS 2026-10-03 03:00 review (chunk file `GEO-7.md`).** 2026-10-03 (§9 item 184): `cad_identity_report` (owned-cell tag sets, reduced volumes, default band `VOLUME_PARTITION_BAND`) passes on all three live generators. The merged-tag control flags exactly `tag_set`, and the diagnostics and tag-integrity modules are green in the same window (`20261003T005340Z_GEO-7-W1.log`, 25 passed). | standard |
| `GEO-8` | **Make `two_torus_domain` a conforming mesh** | ✅ | standard |
| `GEO-9` | **`coil_phantom_domain` / birdcage meshes do not generate** | ✅ 2026-08-03 — step 1 (coil+phantom gated), step 2a (finalize + `bcast`: 180 s hang → 13 s), step 2b (`occ.fragment` rewrite; both identities 1.000000000000, whole `tests/mesh` green in CI). Retires known-issues 7 | standard |
| `GEO-10` | **`two_torus_domain` never emits its `outer_boundary` facet tag** (known-issues 10) | ✅ *(2026-08-06, 00:00 run; known-issues 10 retired)* | standard |
| `GEO-11` | **Boundary-classification margins under OCC bounding-box padding (CAD-only probe sweep)** | ✅ | smoke |
| `GEO-12` | **Widen the two `1e-9` wall tolerances and gate the `outer_boundary` group** (known-issues 12) | ✅ | standard |
| `GEO-13` | **Decouple `cylindrical_domain`'s wall tolerance from `resolution`** (known-issues 13) | ✅ | standard |
| `GEO-14` | **The shared ≈3% geometry floor: faceting vs resolution** (entry lives after `TH-11`, beside the fixtures it measures) | ✅ *(closed 2026-08-15 review on the refuted hypothesis: RESOLUTION, 3.643% → 1.781% at 55 251 cells, rate 1.77 in h — no faceting floor)* | standard |
| `GEO-15` | **Birdcage conductor sizing: is graded sizing a `PORT-9` prerequisite?** (the 0.7091 question; named prerequisite of `PORT-9` step 3) | ✅ 2026-08-16 (graded sizing recovers **0.9670** of the conductor's CAD mass at h_c = 1.6 mm vs **0.7403** baseline, gate cleared, `GEO-9` identities unmoved at < 1e-9; 41 s at `-n 2`; closed by the 10:30 review — the chunk was its one question, now answered by measurement. **2026-08-25: its gate `test_graded_conductor_sizing_recovers_the_cad_mass` is red on `main` on the 0.11 image — the ungraded baseline rung no longer meshes (conductor-sizing axis, not resolution; known-issues 2026-08-25). Disposition is `GEO-21`; the ✅ here is the 0.7.2 close and stands**. **2026-08-26: `GEO-21` ✅ — the gate is green again on a *coarse-graded* control (`BASELINE_CONTROL_RESOLUTION` = 4.8e-3, 0.846150 vs graded 0.966977), so what the live module measures is now **fine vs coarse grading**. The graded-vs-ungraded answer to the `PORT-9`-prerequisite question is the 0.7.2 close and must be cited as such — do not restate it off the module's present numbers**) | standard |
| `GEO-16` | **Emit the gap boxes' longitudinal port-sheet mid-plane in `two_torus_domain`** (the `PORT-9` step-1 mesh prerequisite; commissioned 2026-08-16 18:00 review) | ✅ **2026-08-25: `test_kwarg_off_reproduces_the_recorded_mesh` is red on `main` — the kwarg-off record 79 534 (0.7.2) vs 79 070 measured on 0.11, sheet exonerated by two independent no-sheet builds. Ruled 2026-08-25 18:00 review: re-record licensed, gate constant + both guide copies in one commit (§9 item 3, known-issues entry has the full scope); the ✅ stands. **LANDED 2026-08-25, 22:30 slot** — `NCELLS_UNGATED_RECORD` = **79 070** version-tagged to the 0.11 image with the 0.7.2 digit and both provenance logs in-comment; gate pair `5 passed / 55.84 s / Status 0` printing `cells=79070` and the meshed-band cross-check **0.974490841** inside 0.970–0.980 (`20260826T033222Z_GEO-16-rerecord-gate-pair.log`), `mesh:4` green with the sheeted build properly distinct at **79 940** (`…033350Z`, 31 s) and `mesh:1` at **79 070 / 14.1 s** (`…033431Z`, 16 s); the `mesh:1` docstring + guide and four `mesh:4` guide copies moved in the same commit. No band moved; known-issues entry retired** | standard |
| `GEO-17` | `coil_phantom_domain` region-resolution policy shrinks the coil volumes it refines (−21.68%/−22.62%; `OPS-17` step-2 defect 1, known-issues 2026-08-17; … | ✅ Done: step 1 (2026-08-20, audited compliant) made the region sizes actually apply, raising coil meshed/CAD to 0.8356. The `mesh:5` inverted control that lost separation on 0.11 was re-chosen on 2026-08-25 at `CONTROL_RESOLUTION` 0.018 with margins above the 0.05 separation, so the ✅ stands. *History: `docs/planning/chunks/GEO-17.md`.* | standard |
| `GEO-18` | Birdcage conductor gaps: cut the legs so the port boxes have terminals (`PORT-9` step-3 mesh prerequisite; commissioned 2026-08-20 03:00 review from step 3 legs (a)+(b) 🚫) | ✅ 2026-08-22 (*step 2 audited COMPLIANT 2026-08-22 03:00 review — every figure verbatim in `20260822T020113Z_GEO-18-step2.log:8569-8576`, no pre-existing test touched; one transparency note: the sheets-off control asserts the `110+i` **cell** tags absent, and the `210+i` facet absence the entry and commit claim is implied by that, not measured — `EX-28` is commissioned to assert it directly*; **step 1 ✅ 2026-08-20** — terminals exist: 2.236196e-04 m² per port, **0.988616** of the closed-form `2·π·r_leg²`, all four equal to the printed 7 digits, *audited COMPLIANT 2026-08-21 18:00 review — closed-form band, closure and mass identities verified against all three logs, the pre-derivation red disclosed*; **step 2 ✅ 2026-08-22** — the sheets exist and are exact: meshed sheet area `1.120000000e-04 m²` = the analytic `dx·g` at **1.000000000000** on all four ports, `w_eff = A/h` equal to the bbox extent to 1.000000000000, out-of-plane spread ≤ 2.512e-16 m, half-volumes 0.500000000000 each, C4 sheet spread **8.470e-16**, step 1's terminal band and closure re-asserted on the sheeted mesh) | standard |
| `GEO-19` | `birdcage_port_domain` at `leg_count = 16`, gapped + sheeted: identity family re-gated (C16), cost rung measured | ✅ Closed: the 16-leg gapped and sheeted birdcage's identity family is re-gated under C16, with the 4-leg control reducing to the old flat gate. Its cost rung is Phase 6's first measured rung, and `GEO-20` step 2 is unblocked. *History: `docs/planning/chunks/GEO-19.md`.* | heavy (probe first) |
| `GEO-20` | High-pass birdcage ring-gap port layout (`ring_gap_length`, `2·leg_count` ports, the `GEO-18` pattern on the end rings) | ✅ Closed: the 32-port ring-gap module landed on `main` with the 4-leg and 16-leg ring rungs reproduced, and the 32-port fixture's sheets are licensed to read. Nothing loosened and no band moved. *History: `docs/planning/chunks/GEO-20.md`.* | standard |
| `GEO-21` | Dispose of the red `GEO-15` graded-conductor gate on 0.11: the ungraded baseline (`conductor_resolution=None`) no longer meshes at any global resolution tried, so the whole gate is non-executing on `main` (found by `EX-30` leg (mesh); commissioned 2026-08-25 18:00 review; known-issues 2026-08-25) | ✅ **2026-08-26** *(status marker flipped by the 18:00 review audit — the prose entry and the `GEO-15` row already read ✅; history: step 1 measured 2026-08-26, blocked on a ruling — the candidate control `h_c = 3.2e-3` recovers **0.916742** of the CAD mass, which is **neither** pre-stated branch: not ≤ 0.90 "clearly below", not clearing 0.95, and **inside the module's own `CAD_MASS_GATE - 0.05` = 0.90 separation guard**, so branch (2) as named cannot produce a green gate without loosening that guard. Graded side green — 0.966977 at 98 666 cells. Coarse-ward ladder measured and handed over, nothing adopted; see the entry. **Ruled 2026-08-26 03:00 review: option (b), control = 4.8e-3** — separation 0.846150 vs 0.966977 with the 0.90 guard unmoved, demoted claim (fine-vs-coarse grading) stated; 6.4e-3 rejected for cliff adjacency. **Step 2 ✅ 2026-08-26 (04:30 slot) — landed as written and the chunk closes**: control `None` → 4.8e-3 version-tagged with the six-rung probe table in-comment, demoted claim (fine vs coarse grading) in the module docstring and the `mesh:3` guide, gate `1 passed in 41.11s` at 0.846150 / 0.916742 / 0.966977 with the 0.95 gate and the 0.90 separation guard unmoved, `mesh:3` green at separation 0.120826, docrefs `dead=0 guide=0 stale=13 exit=2` giving `meshing` 2 → 0; known-issues gate red retired, generator-continuum finding re-headed and open)* | standard |
| `GEO-22` | `straight_wire_domain` coarse-resolution floor on 0.11: bisect the `[0.008, 0.010)` threshold and land a measured guard so a too-coarse request raises legibly instead of aborting inside gmsh (owner for the `EX-30` leg (root) finding; commissioned 2026-08-26 18:00 review; known-issues 2026-08-25, re-headed 08-26 and again 08-28) | ✅ *(**demoted ✅ → 🧪 2026-08-29 10:30 review, §4 clause 3** — the closing evidence is a probe that asserts nothing plus a re-run of `test_geometry_failure_is_collective.py`, whose `total_caught == comm.size` is the did-raise property the 08-28 18:00 audit already classified as non-quantitative for `GEO-23`; the fix is step 2c, one asserted leg-D anchor — queued §9)* — **step 1 ✅ 2026-08-28 as a measured negative: there is no floor.** The 2.5e-4 sweep of `[0.008, 0.010]` reads **non-monotone on both geometries** — `h = 0.00875` fails while coarser rungs mesh — and reproduces **bit-identically** across two runs, so no `RESOLUTION_FLOOR` can be written and none was; `src/` untouched. **Step 2 ✅ 2026-08-29 (07:30 slot) ⇒ chunk ✅** — the size-field probe leg the 08-28 10:30 ruling commissioned reads **18/18 OK and 0/18 fallbacks** against leg C's 7 FAILs and 18/18 fallbacks, on a leg-C control that reproduces step 1 bit-identically in its own process; `src/` still untouched and the size-field *licence* is the weekly review's. **Step 2c ✅ 2026-08-29 (16:30 slot) ⇒ chunk back to ✅** — the probe's finding is now asserted: `tests/mesh/test_straight_wire_size_field_probe.py` reads **19 823 cells / 0 fallbacks** patched and **21 830 / 1** unpatched on the example geometry at `h = 0.008`, 0.00% against the pre-stated ±1% band, identical at `-n 1` (8 s) and `-n 2` (7 s); the perturbed-reference negative control reds on both rank streams and the restored constant re-runs green. `src/` untouched, no record moved, known-issues entry stays OPEN | smoke (probe `-n 1`) + standard (gate) |
| `GEO-23` | The 0.11 "Invalid boundary mesh (overlapping facets)" family: one owner for the three `OPS-26` census reds (coil+phantom generator, `birdcage_port_domain` … | ✅ Closed: the overlapping-facets family has one owner and its gate is load-bearing and rank-symmetric, with `src/` untouched. The `GEO-21` residual is unchanged and not this chunk's. *History: `docs/planning/chunks/GEO-23.md`.* | smoke (`-n 1` probes) + standard (`-n 2` ladders) |
| `GEO-24` | Give `birdcage_port_domain` the `shared_facet` ghost layer, and re-read every module that reconstructs a sheet on it at `-n 2` and `-n 12` | ✅ Closed: the `shared_facet` ghost layer is plumbed and changed nothing this family could see, so `GEO-20` step 1's all-ranks reading loses its width-conditional caveat. Nothing loosened, no band moved, and no record in `tests/` touched. *History: `docs/planning/chunks/GEO-24.md`.* | standard |
| `GEO-25` | F-human cost probe | ✅ Closed as the F-human cost probe: branch B's CAD identities and cell record are registered, not a solve, a port gate or a physics fixture. F-human physics stays Phase 6's and the weekly's. **Re-recorded 2026-10-11 (`GEO-35` 2b-ii part 2, (1\*) licence):** the rung module's default build drops the leg port boxes; the branch B / A cell records move to new `F_HUMAN_BRANCH_B_CELL_RECORD_UNBOXED` 500 425 (−0.8356 %) / `…_A_…_UNBOXED` 203 945 (−0.5035 %), reproduced exactly (`20261012T021903Z_GEO-35-s2bii-W3-rung.log:8772,17804`; 10 passed 4 skipped, 14 collected as at the census); the boxed 504 642 / 204 977 stay for their boxed consumers; partition, 32 terminal ratios and the CAD separation (0.965414 / 0.893028) unchanged. *History: `docs/planning/chunks/GEO-25.md`.* | heavy (probe first) |
| `GEO-27` | Phantom-resolution cost probe for the 1 g SAR rung | ✅ 2026-10-12: priced at 0.005 / 0.0025 m, degree 2, `-n 8` (130 427 / 201 072 cells, 836 890 / 1 284 856 dofs, body 132.86 / 319.07 s, 20.72 / 44.05 GiB); 1 g C4 pairs ≤ 0.0446 % at 0.0025 m (`20261012T034956Z_GEO-27-step2-W2-h0025.log`). Step 1 (2026-09-08, `mesh-probe`): the phantom-resolution cost table on the ragged F-small — `h_p` 0.0075 / 0.005 / 0.00375 / 0.0025 m → 120 499 / 129 505 / 144 212 / 199 920 cells, 4.96 cells across the 1 g ball at 0.0025; only the 0.0075 rung repeated (`20260908T033218Z_GEO-27.log`, 132 s). **Done-when (re-posed by the 2026-10-10 weekly; amended 2026-10-11 under rule (j)):** step 2 — the 1 g rung priced on the exact-rectangle degree-2 F-small fixture at `h_p` = 0.005 and 0.0025 m, `-n 8` (cells, degree-2 unknowns, wall time, summed peak memory in one table; a rung that does not fit the heavy tier is recorded as such), **and** `MAT-4` step 5's four cyclic 1 g C4 pairs asserted at the imported `C4_COVARIANCE_BAND` on the finest rung that fits, in one logged window ⇒ ✅. A measurement with one imported gate; no band moves, no convergence claim. *History: `docs/planning/chunks/GEO-27.md`.* | heavy |
| `GEO-28` | C4 census of the unloaded F-small birdcage mesh | 🚫 **Closed as overtaken 2026-10-05 (03:00 review, lane-1 disposition the 10-03 weekly asked for by 10-10):** not a ragged-cut census. Its only consumer, `WF-6`, closed ✅ on a re-scoped target, and nothing reads a quadrant census of the exact fixture. The record stays valid as a record. 🧪 Measured: the C4 census of the unloaded F-small birdcage mesh is delivered, and the near-field shell carries on the order of one cell at the global resolution, supporting `GEO-29`'s premise. It is measurement-only and never ✅. *History: `docs/planning/chunks/GEO-28.md`.* | standard (heavy by ceiling), `-n 2` — measured 29 / 30 s |
| `GEO-29` | Global-resolution cost ladder of the unloaded F-small birdcage port fixture | 🚫 **Closed as overtaken 2026-10-05 (03:00 review):** not a ragged-cut census. The exact fixture's cost ladder is `ANS-4` items 212 / 225: 117 206 / 282 198 / 259 524 cells, mesh 23.3 / 50.2 / 55.7 s (`docs/planning/chunks/ANS-4.md`). This record stays valid as a record. 🧪 Measured: the global-resolution cost ladder of the unloaded F-small fixture is delivered as a cost table. None of its finer rungs may be pinned as a version-tagged record without its own repeat. *History: `docs/planning/chunks/GEO-29.md`.* | heavy by ceiling, `-n 1` — measured **145 s** (probe 141.30 s) |
| `GEO-30` | Does the four-port birdcage fixture stay C4-symmetric under refinement? One geometric census across *both* refinement axes | 🚫 **Closed as overtaken 2026-10-05 (03:00 review), per the 10-02 ruling:** on exact rectangles, `ANS-4`'s (F) asserts every strip equals the SPEC rectangle at `STRIP_BAND` on the ×1, h = 0.0075 and conductor-×0.5 rungs. The printed strip facet counts are C4-equal (40 / 44 / 76 per port). This record stays valid as a record. 🧪 Measured 2026-09-09: the four-port fixture's geometric census across refinement is delivered and the mesh is excluded as the mechanism. The ruling on `ANS-4` step 2a and `WF-6` step 4f was left to the review. *History: `docs/planning/chunks/GEO-30.md`.* | standard by expectation, heavy by ceiling, `-n 2`, no solve — priced from `GEO-29`'s 21.6 / 26.4 / 32.9 / 45.6 s per mesh and `GEO-28`'s one-window census; **actual 133 s** for all four rungs in one window |
| `GEO-31` | The port gap-sheets' *triangulation* under the C4 rotation | 🚫 **Closed as overtaken 2026-10-05 (03:00 review), per the 10-02 ruling:** the in-sheet triangulation reached the ports only through the ragged midpoint filter (`_narrowed_radial`). An exact strip is a CAD-bounded rectangle. 🧪 Measured 2026-09-09: the port gap-sheets' triangulation differs under the C4 rotation, and refinement amplified an existing asymmetry rather than creating one. The follow-on is `PORT-18`. *History: `docs/planning/chunks/GEO-31.md`.* | standard by expectation, heavy by ceiling, `-n 2`, no solve — priced from `GEO-30`'s 133 s and **measured 144 s** (meshing 25.29 / 33.93 / 31.25 / 36.97 s) |
| `GEO-32` | Force a C4-congruent port-sheet cut by construction | 🚫 **Closed as overtaken 2026-10-05 (03:00 review), per the 10-02 ruling:** `c4_congruent_sheets` refuses `leg_sheet_strip_fraction` (`io/mesh.py:4266–4271`), so it cannot run on exact rectangles. Exact strips removed the C4 break it was built to fix, and the switch stays for the older records. 🧪 Measured and accepted by the 2026-09-10 18:00 review: a C4-congruent cut carries the whole read-back spread, and `ANS-4` step 2a′ showed it carried the `Z` break too. It is measurement-only. *History: `docs/planning/chunks/GEO-32.md`.* | standard (heavy by ceiling), `-n 2`, no solve — **measured 127 s + 81 s + 54 s regression** |
| `GEO-33` | F-human-32: the 32-leg birdcage rung | 🟡 Steps 0a and 0b are green (2026-10-02, §9 items 149 and 170): the two-sided Leontovich lip term on the cracked box reads Q_sep/Q₀ = 0.772322 against Pozar's 0.771914 (0.053 %), and the merged single-valued sheet control reads 1.00000000. Step 1 (the 32-leg cost probe) and a single priced degree-1 strip-sheet drive remain, sequenced after the F-human degree-order question. **Done-when:** 0a + 0b green ⇒ 🟡 (the lip term exists and is gated); ✅ additionally needs step 1 (the 32-leg cost probe) and the strip-sheet rung's single priced drive. *History: `docs/planning/chunks/GEO-33.md`.* | unmeasured (cost probe first) |
| `GEO-34` | STEP import through Gmsh/OpenCASCADE with named-group tagging and a case config file | 🟡 Steps 1a (STEP import with tag census), 1b (through the congruent-sheet config, 2026-09-26) and 2 (the preflight writer, 2026-09-30) are done, and §9 item 119 then wrote the import's own wall group so the wall re-enters step 2's identity. Steps 1c (independent kernel) and 3 remain open. Step 1d (§9 item 289) is BLOCKED at W0 (2026-10-07, `20261007T034333Z_GEO-34-step1d-W0-names.log`): Gmsh 4.15.2's OCC STEP writer keeps neither the face name nor the embedded rectangle, so a named interior face needs a STEP from another kernel (1c). U12a record (§9 item 333, 2026-10-08, 🧪 measurement, glyph unchanged; 1d family 2 of 4): no Gmsh-OCC export form (fragment / solid / compound / shell / two-file) × read (default / `highestDimOnly=False`) returns an **embedded** 0.48 m² face after `removeAllDuplicates`. The compound (r2 only) and the two-file merge (both reads) return it **free**. That is amber, read as red. The default read drops the compound's free face (`20261008T214725Z_GEO-34-step1d-probe.log`, repeat bit-identical). **Step 1d record (U12′, §9 item 985, 2026-11-14 12:00 slot; family 3 of 4, green):** a CadQuery 2.8.0-authored STEP (committed, sha256-pinned) with product `PortSheet` is mapped by name to `port` (tag 201) through two additive importer keys (`highest_dim_only`, `surface_groups`), embedded by `occ.fragment` and traced through its `out_map`; at `-n 1 / 2 / 4` the owned tag-201 area is 0.48 to ≤ 2.2e-16, `dS(201)` to ≤ 8.9e-16, 192 = 192 triangles, 0 one-sided (`20261114T181927Z_GEO-34-U12p-W1.log`); key-absent gates re-run green. CAD identity only, no port solve. **Step 1e** (§9 item 996, 2026-11-15, 🧪 measurement, glyph unchanged): Netgen's optimiser crashes with SIGSEGV on the embedded-face arm at every h over 0.04–0.16, and never on box-only or the face-spanning two-volume split (readout (a); `20261115T223203Z_GEO-34-1e-E-h0p08.log`). *History: `docs/planning/chunks/GEO-34.md`.* | standard (steps 1–2); heavy by prediction (step 3, probe first) |
| `GEO-35` | Drop the floating leg-port boxes from unbroken-leg birdcages | ✅ 2026-10-11 (amended Done-when, rule (j), met). Step 1 (item 275): `emit_leg_ports=False` leaves zero leg-port cells and the ring census intact, the `GEO-25` CAD identities hold (`20261007T004311Z_GEO-35-s1-W1-mesh.log`); item 357: transverse (N′) 0 (`20261010T112516Z_GEO-35-s2bi-p2-W2.log`). Step 2a (item 351): solid F-human drops the boxes, hole path stays boxed (369 387 cells), six `WF-7` constants re-recorded. **2b-ii:** `WF-7` step 3 + `EX-62` (item 385a: 532 935 cells, SAR10g(c\*) 2.079992000e-04 rel 9.864e-11, `20261012T020246Z_GEO-35-s2bii-W1.log`) and `GEO-25` transverse + `mesh:12` (item 385b: `F_HUMAN_BRANCH_B_CELL_RECORD_UNBOXED` 500 425 and provisional `…_A_…_UNBOXED` 203 945, exact, `20261012T021903Z_GEO-35-s2bii-W3-rung.log`, `20261012T022429Z_GEO-35-s2bii-W4-mesh12.log`) re-recorded under (1\*), each reproducing its 2b-i reading; gated anchors green at unchanged bands. `F_HUMAN_BRANCH_B_CELL_RECORD` 504 642 stays for its three boxed consumers (Q5), `mesh:9`–`11` keep their boxes (Q2), κ-era `WF-10` / `ANS-8` stay boxed with `emit_leg_ports=True` (Q3). Branch A's `…_UNBOXED` awaits the 10-12 weekly's ratification. *History: `docs/planning/chunks/GEO-35.md`.* | heavy (mesh test 300–366 s; re-record windows 230–424 s, `-n 2` / `-n 8`; audit 2026-10-12 03:00) |
| `GEO-36` | A rod in the phantom: the implant-in-phantom generator (Phase 7 (7b)) | ✅ Audit PASS 2026-10-02 (CAD identities at `-n 2`); the kept rod on F-small, F-human and in F-clinical's ASTM box is a family of records in the chunk file — entry 48 (I36), I57 / I67, C32 and entry 61 (I36-3T), entry 55 (I32, its tip gate re-posed on the consistent-mass projection as I32b and re-run as `xl` entry 63), and the stub proofs I32b, I85, I87′ (the `ROD_THERMAL` knob) and I85c. Open: entry 63's read and the steps serial on it (the probe readout, entry 67), the 11-14 weekly's. *History: `docs/planning/chunks/GEO-36.md`.* | standard |
| `GEO-37` | Does an interior 2D face reach DolfinX's facet tags on the 0.11 image? The CAD-port probe | ✅ GREEN as a probe 2026-10-04 (audit PASS 2026-10-05): on the 0.11 image an interior dim-2 group reaches `facet_tags` complete and unduplicated at `-n 1 / 2 / 4` (owned area − 0.48 ≤ 1.1e-16), so native facet tagging is the (9b) CAD-port route. Step 2 (2026-10-06, §10 U11) assembles `dS` over the tag inside `STRIP_BAND` on a `shared_facet`-ghosted mesh, with the unghosted `-n 4` hang reproduced as a record. **Done-when:** a box with one embedded interior rectangle tagged as a dim-2 physical group, read at `-n 1 / 2 / 4`: the tagged facet count equal across widths after reduction and the tagged area against the rectangle's closed form at 1e-9, no hang, in one logged window ⇒ ✅ as a probe, green or red. Green ⇒ native facet tagging for CAD-drawn ports; red ⇒ the importer inserts a sliver pair and reuses `_interface_facet_tags`. Measurement only, no `src/` change. *History: `docs/planning/chunks/GEO-37.md`.* | smoke by prediction |
| `GEO-38` | F-clinical-1.5T cost probe | ✅ as a probe (audit PASS); the clinical fixtures' windows, roots and matrices are records in the chunk file. Newest: **entry 56** (read first branch, 11-07 weekly) — C_in(F-clinical-1.5T, d2) = 185.7629 pF at 0.25 pF (`20261105T080008Z_GEO-38-C5b-prime.log:12814`); **entry 54** (read second branch) — the degree-2 1.5 T 64×64 a record, `C_mode`(d2) 185.4048 pF (C18), the asymmetric class not degree-1 discretisation (×1.44); **entry 57** (RUN 2026-11-08, first branch as the 03:00 review reads it) — C_in(F-clinical-3T, d1) = 44.0525 pF at 0.05 pF (`20261108T080008Z_GEO-38-C16r.log:13229`), entries 58 / 61 queued at it. C5e DONE 2026-11-08 (item 869: the 3 T degree-2 knob and save path on C11's module, U₂ 7 335 988 exact at build, `20261108T180244Z_GEO-38-C5e-W1-build.log`). C35 DONE 2026-11-08 (item 889: the r2 rung knob on C3's module, 912 804 cells / U₁ 1 184 518 exact at build, `20261108T193325Z_GEO-38-C35-W1-build.log`). C33c DONE 2026-11-14 (item 979): `c33.custody_json` formats c\*'s dz signed (`+.6f`), so (J58) reads 76 fields with 0 mismatches on entry 58's JSON, and the (N-sign) control fails on dz alone (`20261114T133421Z_GEO-38-C33c-W1.log`); P89′ reads 0 string-only (`20261114T133435Z_MAT-21-P89-W1c.log`). Open: C36 / C43 (`ANS-11`, items 871 / 880); C33 / C37 / C39 / C40 after their windows; entries 58, 59, 61 (`xl`) and 62 (`xxl`) queued. Narrative to 11-07: the chunk file. *History: `docs/planning/chunks/GEO-38.md`.* | heavy (probe first; one degree-1 column may need `xl` — price the mesh first) |
| `GEO-26` | Longitudinal ring-gap port sheets for `birdcage_port_domain` | ✅ Closed: the longitudinal ring-gap port sheets exist in `birdcage_port_domain`, and the 16-leg terminal-area states are measured and stated in the comment. The 4-leg band did not move. *History: `docs/planning/chunks/GEO-26.md`.* | standard (step 1 ≈ 50 s per width, measured; step 2 declared heavy, 158–160 s per width measured; step 3 heavy by measurement, 221–226 s per window) |


**`GEO-24` — plumb the birdcage ghost layer, then re-read the fixture's records at two widths** ✅ *(steps 1a
2026-08-28 21:00, 1b 22:30, 2a′ 2026-08-29 04:30, 2b 2026-08-29 06:00 ⇒ chunk closed; commissioned by the human
operator 2026-08-28; the table row carries every step's digits; the commissioning plan archived in
`docs/planning/plan-archive.md`, entry «§7 GEO-24 full narrative (commissioning plan; closure recorded in the table row) — archived 2026-08-30 (weekly review)».)*
> * **Cause and fix:** `birdcage_port_domain` called `_model_to_mesh` with no partitioner (`GhostMode.none`), so a
>   facet reconstruction across a partition boundary could miss its second cell. One kwarg at `io/mesh.py:3356`,
>   `create_cell_partitioner(GhostMode.shared_facet, 2)` (`470f410`); nothing else in `src/`.
> * **Before (1a, 14 windows / 668 s; 1b, 13 windows / 660 s):** every cell count identical across widths (116 085 /
>   98 666 / 128 111 / 114 655 / 116 475 / 307 296); two `-n 12` reds, both facet reconstruction — `ring_gaps` P8
>   closure **0.990103697427** (175 vs 176 facets) and `port_terminals` phantom↔air **245 / 0.935322** vs 255 / 0.979885;
>   the five `tests/validation/` consumers clean at both widths (σ_max 0.999992805, spreads 0.0553 / 0.0353 / 0.0214%).
> * **After (2a′ 6 windows / 246 s, `20260829T0930…–0934…_GEO-24-step2aP-*`; 2b 11 windows / 485 s,
>   `20260829T1100…–1108…Z_GEO-24-step2b-*`):** P8 **176 / 1.000000000000** at `-n 12`; phantom↔air **256 / 0.984183 at
>   `-n 1` / `-n 2` / `-n 12`** — ruled 08-29 03:00 a defect repair (the `-n 1` truth), not a re-baseline; every `-n 2`
>   validation digit reproduces 1b, the only `-n 12` movement `Z_11` at 4.1e-9. Known-issues entry RETIRED; `GEO-20` step 1's "1.000000000000 on all 12" loses its width caveat; `GEO-20` step 2 becomes an unblocked re-run. * **Carry-forward:** the two-torus `-n 12` gap-ratio drift (0.894141 → 0.894274 at 184 176 cells) is **`PORT-12`**.

> `GEO-4`'s substance is discharged for the two-torus fixture (`air_padding` +
> graded sizing), but it stays 🧪 until its own test executes. **Every other
> fixture in `io/mesh.py` still uses a single global `setSize` and tight padding,
> including coil+phantom** — expect the same boundary-mirror error that cost 20%
> on Helmholtz, and expect graded sizing to be equally necessary.

**`GEO-4` step 1 — off-centre domain sizing (known-issues 5)** ✅
*(2026-08-06; full diagnosis in `docs/planning/plan-archive.md`)*. The
oldest standing failure on `main` was the test's own assertion: the
`radial_clearance` guard means the coil always governs the box, so "an
offset phantom grows the box" is false for every meshable configuration.
The test now gates the containment identity with the clearance term
explicit plus the exact identity `clearance(centered) − clearance(shifted)
= 0.03 m`; the `OPS-11` `--deselect` is gone and `tests/mesh` runs
unexcluded in CI. Handed to a review: the overlap guard is z-blind. `GEO-4`
itself stays 🧪 (graded-sizing generalization is separate work).

**`GEO-10` — `outer_boundary` facet tag never emitted (known-issues 10)** ✅
*(2026-08-06; full narrative in `docs/planning/plan-archive.md`)*. Cause:
gmsh inflates OCC bounding boxes by the geometric tolerance (measured
1.000e-07) and the flat-against-wall test used `tol = 1e-9` — every wall
failed and the physical group was silently skipped. Fix: `1e-9 → 1e-6`.
Gate: tag sets exactly `{1}` / `{1, 201, 202}`, tagged area = analytic box
surface at ratio 1.000000000000000 (an identity — planar walls); neither
Helmholtz consumer depends on tag `1` (digit-identical regression).

**`GEO-8` — make `two_torus_domain` a conforming mesh** ✅ *(2026-08-01;
full narrative archived in `docs/planning/plan-archive.md`)*
> The fixture never fragmented — three disconnected components, `PORT-1`'s
> `Z₁₂ ≡ 0`. Fixed with `occ.fragment` + centroid/mass re-derivation of the
> physical groups (**fragment renumbers — never trust its returned tag
> order**). Mesh volume / analytic box 1.002633 → 1.000000000; Helmholtz
> centre-field error 1.731% → 0.728%. Measurement worth keeping: at uniform
> `resolution=0.01` a meshed torus retains only 0.598 of its analytic
> volume — **the wire needs `wire_resolution ≲ 0.4·minor_radius`** before
> any volume-based conformity statement means anything.

**`GEO-9` — `coil_phantom_domain` / birdcage meshes do not generate** ✅
*(steps 1, 2a, 2b ✅ 2026-08-03; known-issues 7 retired; full diagnosis in
`docs/planning/plan-archive.md`)*
> Two independent defects, neither the guessed one: `coil_phantom_domain`
> was innocent; the **birdcage** raised without reaching `gmsh.finalize()`,
> poisoning later meshes and hanging rank 1 at a skipped collective. Fixes:
> rank-0 body isolated with the failure `bcast` to every rank (180 s hang →
> 13 s prompt failure — the isolation gate's no-hang `allreduce` assertion
> degenerates at `-n 1`, never move it to a single-rank job); one
> `occ.fragment` against all tools with every group re-derived from the
> out-map. Both volume identities 1.000000000000, all four port boxes
> exact, conductor 0.7091 banded — that number is exactly what `GEO-4`'s
> open half (the birdcage's single global `setSize`) measures. The mesh is
> 8.95 s, not the known-issues "≈10 minutes" (that was the hang).

**`GEO-11` — boundary-classification margins under OCC padding** ✅
*(2026-08-06; full sweep narrative in `docs/planning/plan-archive.md`)*.
CAD-only two-sided margin gate over four fixtures
(`test_boundary_classification_margins.py`): found `GEO-10`'s exact defect
live in `loop_over_half_space_domain` (0/12 walls accepted) and
`sphere_in_box_domain` (0/7) — became known-issues 12 → `GEO-12`; and
`cylindrical_domain` under-separated at 4.50× — known-issues 13 →
`GEO-13`. Not covered, deliberately: `coil_phantom_domain` /
`birdcage_port_domain` (their CAD stages would have to be factored out of
the generators first — a review's call).

**`GEO-12` — widen the two `1e-9` wall tolerances** ✅ *(2026-08-06;
known-issues 12 retired; full narrative in
`docs/planning/plan-archive.md`)*. Both tolerances → `1e-6`; accepted
counts land exactly on each fixture's wall count; meshed facet-tag gate
(`test_wall_boundary_tag_areas.py`) asserts the analytic cube-surface area
at 1e-9. All six discarding callers re-run: no landed number moved a
digit.

**`GEO-13` — decouple `cylindrical_domain`'s wall tolerance from
`resolution`** ✅ *(2026-08-07; known-issues 13 retired; full narrative in
`docs/planning/plan-archive.md`)*. Tolerance is now
`0.01 × (outer_radius − inner_radius)` — chosen from a sweep window
`[1e-4, 0.05]` measured on all four repo call sites; classification
unchanged (3 of 6 surfaces), margins now 1.1e-04× / 1.0e+02×. All four
fixtures in the margins file assert live (the `GEO-11` sweep closes). New
precondition at the use site: a radial gap below ≈1e-4 m stops clearing
the OCC padding by 10× (smallest gap in the repo: 0.07 m).

**`GEO-15` — birdcage conductor sizing: is graded sizing a `PORT-9`
prerequisite?** ✅ 2026-08-16 *(entry written 2026-08-16, 03:00 daily
review; step 1 executed by the 07:30 implementer run; **closed by the
10:30 review** — the chunk was scoped as a single measured question and
step 1 answered it (audited §4-COMPLIANT: pre-stated ≥ 0.95 gate, logs,
elapsed recorded). The implementer left the flip to the review; the
remaining 3.3% is curvature faceting, not a mesh-size failure, and pinning
it down buys `PORT-9` nothing. One latent hazard the audit found — a
rank-local `perf_counter` budget break in the ladder test that could
desync collectives if it ever fires — is recorded in known-issues, not
here. Mesh-only: no solves.)* The birdcage mesh (`GEO-9`) keeps only **0.7091** of the
conductor's analytic volume under the single global `setSize = 0.015` —
part is the analytic sum double-counting the 8 leg∩ring junctions (CAD
masses give 0.9578), the rest is 0.015 against a 0.004 ring minor radius.
`GEO-8`'s measured rule (`wire_resolution ≲ 0.4·minor_radius`, i.e.
≲ 0.0016 here) says the conductor is ≈10× under-resolved, and a lumped
port on that surface inherits the coarse conductor boundary. This chunk
answers `PORT-9` step 3's open question by measurement.
> * **Step 1 ✅ 2026-08-16** — gate **0.967019 ≥ 0.95** of CAD mass at
>   h_c = 1.6e-3 (98 474 cells, 16.74 s; baseline global-`setSize`
>   0.740335 at 48 245 cells / 6.07 s; negative-control separation
>   0.2267); the junction double-count is **4.22%** (CAD 1.030097043e-04
>   vs analytic 1.075503356e-04 m³); `GEO-9` identities unmoved < 1e-9 on
>   every rung. Logs `20260816T123337Z_GEO-15-step1.log` (1 passed, 41 s,
>   `-n 2`) + `20260816T123433Z_GEO-15-step1-regression.log` (4 passed,
>   21 s). **Live carry-forwards:** new `birdcage_port_domain` kwargs
>   `conductor_resolution` / `conductor_refine_distance` /
>   `return_diagnostics`, all defaults unchanged; the working mechanism
>   is a Distance→Threshold field over the conductor surfaces — the
>   three `Mesh.MeshSizeFrom*` switches must be off or gmsh re-imposes
>   the coarse size; `PORT-9` may assume graded sizing and budgets from
>   98 k cells; the residual 3.3% is curvature faceting; `GEO-4` stays 🧪
>   (no solve ran). Full narrative + original plan:
>   `docs/planning/plan-archive.md`, archived 2026-08-16.

**`GEO-16` — emit the gap boxes' longitudinal port-sheet mid-plane in
`two_torus_domain`** ✅ *(closed 2026-08-17; commissioned 2026-08-16 18:00 review as the mesh prerequisite `PORT-9` step 1 named.)*
> **Result ✅ 2026-08-17** (`tests/mesh/test_two_torus_port_sheet.py`, `20260817T003627Z_GEO-16-regression.log`, 5 passed, 47.3 s at `-n 2`, standard). Opt-in `emit_port_sheet=False` (default) fragments each gap box on its mid-plane `z = ±separation/2` (dim-2 tool to the existing `occ.fragment`), halves carried as cell tags `101`/`111`, `102`/`112` by centroid z; sheet facet tags `211`/`212` rebuilt dolfinx-side via `_interface_facet_tags` (now accepts a sequence of cell-tag pairs; no dim-2 gmsh group, known-issues 9).
> **Anchor met:** MPI-reduced `dS` area per sheet **9.573030358733e-05 m²** vs CAD mid-plane **9.573030358733e-05 m²**, `meshed/CAD = 1.000000000000` (band 1e-9); 84 owned facets per sheet, asserted non-empty before the identity; out-of-plane spread 3.5e-18 m; sheets agree < 1e-12.
> **Measured extents (printed, never gated — what `PORT-9` step 1 must use instead of nominal dimensions):** `w = 1.200000000e-02 m`, `h = 7.977525299e-03 m`, `w/h = 1.504225878` squares, `area/(w·h) = 1.000000000` (CAD bbox `w/h = 1.504206917`, 5th-digit difference from gmsh's 1e-7 bbox inflation, `GEO-10`).
> **Negative controls held:** kwarg off — 79 534 cells, tag sets `{1,2,3,101,102}` / `{1,201,202}`, no `21x` group; 3b-iv gate reproduces `meshed/analytic = 0.974490841` bit-identical to the 2026-08-05 record. Fragmented-mesh port areas 1.563786482e-04 m² per port (same 0.9745 of analytic).
> **Carry-forward:** a caller selecting the gap volume by tag must take **both** halves of each box when the kwarg is on (`PORT-9` step 1).
> Full narrative: `docs/planning/plan-archive.md`, entry «§7 GEO-16 full narrative — archived 2026-08-23 (weekly review)».

**`GEO-17` — `coil_phantom_domain`'s region-resolution policy shrinks the coil volumes it refines** ✅ *(step 1
2026-08-20 06:00 slot; audited COMPLIANT 2026-08-21 18:00 review; commissioned 2026-08-17 10:30 review from
`OPS-17` step-2 defect 1 (−21.68% / −22.62% coil volume under a finer request). Full narrative archived in
`docs/planning/plan-archive.md`, entry «§7 GEO-17 full narrative — archived 2026-08-30 (weekly review)».)*
> * **Hypothesis refuted, mechanism measured:** there were no size fields at all — `gmsh.model.getBoundary` with the
>   default `combined=True` returned **0 points** for every region (`20260820T110127Z_GEO-17-step1-diag.log`), so the
>   policy meshed the coil at the air's 0.020 clamp. **Fix:** per-region sizes as a `Min` over per-volume `Constant`
>   fields (`VolumesList`, `IncludeBoundary=1`) set as the background mesh.
> * **Measured (`20260820T110549Z_GEO-17-step1-final.log`, `-n 2`, 13 s):** coil_1 **+10.7169%** (1.191750413e-04 →
>   1.319468693e-04 m³, meshed/CAD 0.754685 → **0.835563**), coil_2 **+10.7851%** (0.752565 → 0.833730), phantom +0.9374%,
>   air −0.2643%; both meshes partition their volume at **1.000000000000**; uniform column reproduces the `OPS-17`
>   record at 1e-9 (negative control). Probe forcing `MeshSizeExtendFrom*` off moved coil_1 −3.12% and was reverted
>   (`…110302Z_…-fieldfix.log` vs `…110407Z_…-probe-defaults.log`).
> * **Band replaced with its measurement:** the strict xfail's 5% "must not move" premise is false for a curved region;
>   now `test_region_resolution_policy_refines_the_tagged_volumes_toward_cad` (policy > uniform, meshed/CAD ≤ 1.0, coil recovery ≥ 0.755). Known-issues "Four defects" §1 resolved; `mesh:5` control re-chosen 2026-08-25 (row).

**`GEO-18` — birdcage conductor gaps: cut the legs so the port boxes have
terminals** ✅ *(both steps closed — step 1 2026-08-20, step 2 2026-08-22; commissioned 2026-08-20 03:00 review from `PORT-9` step 3 legs (a)+(b) 🚫. Review decision: cut the legs, not the end rings — drive direction `ẑ` for every port, square section `dx = dy` makes the four-port layout C4-invariant by construction; port azimuths move from 45° + k·90° to the leg positions k·90° in the gapped variant, a deliberate low-pass-birdcage physics change.)*
> * **Step 1 ✅ 2026-08-20** (`leg_gap_length` opt-in, default `None` bit-for-bit; `_birdcage_leg_gap_layout` helper). `g = 8 mm`, box `(1.400000e-02, 1.400000e-02, 8.000000e-03)` m, `h_c = 1.6e-3`. Per port: terminal area **2.236196e-04 m²** = **0.988616** of `2·π·r_leg² = 2.261947e-04 m²` (band [0.95, 1.0]), equal on all four ports to 7 digits; closure `(A_cond + A_air + A_phan)/A_box = 1.000000000000`; phantom-facing 0; gap volume meshed/analytic **1.000000000000**; `GEO-9` partition identities < 1e-9; gapped meshed/CAD conductor **0.970152** (≥ 0.95); **114 846** cells, mesh 22.61 s, rung 24.32 s. Mass identity re-derived off a cancellation (§4 `MAG-10`/`MAG-15` precedent): `CAD_gapped/(CAD_uncut − 4·π·r_leg²·g)` = **1.000000000192** (1e-9), the difference form read 0.999999994733. Negative control: kwarg off reproduces **98 474** cells, `EX-21` 0.967019, conductor-facing `0.000000e+00 m²` on all four ports. Logs `20260820T093433Z_GEO-18-step1.log`, `20260820T093603Z_GEO-18-step1-final.log` (8 passed, 136.61 s), `20260820T093830Z_GEO-18-step1-record.log` (1 passed, 45.16 s), all `-n 2`. `tests/mesh/test_birdcage_port_terminals.py` kept as the standing guard on the default geometry.
> * **Step 2 ✅ 2026-08-22** (`emit_port_sheets` opt-in, `ValueError` without `leg_gap_length`; attempt 1 parked on `attempt/GEO-18-step2-20260822T004500Z`, exit 124 from an `allreduce` inside `if comm.rank == 0`, log `20260822T003614Z_GEO-18-step2.log`; fix hoisted the collective). Per port (`20260822T020113Z_GEO-18-step2.log`, 2 passed, 53 s, `-n 2`): sheet **54 facets, 1.120000000e-04 m²**, meshed/analytic `dx·g` = **1.000000000000** (1e-9); `h = 8.000000000e-03 m`, `w_eff = A/h = 1.400000000e-02 m`, `w_eff/w_bbox = 1.000000000000`; out-of-plane spread **2.512e-16 / 9.714e-17 m** (band 1e-12); half-volumes **0.500000000000 / 0.500000000000**; step-1 gates survive (terminal ratio 0.988616, closure 1.000000000000, partition < 1e-9); C4 sheet spread **8.470e-16**. Sheeted mesh **116 416** cells (mesh 22.73 s, rung 24.77 s); cell tags `100+i`/`110+i`, sheet facet tags `210+i`. Negative control: sheets off reproduces step 1 exactly (114 846 cells, 0.988616, tags `[1, 2, 3, 101, 102, 103, 104]`, every `110+i`/`210+i` absent). Regression: birdcage mesh suite 10 passed, 186 s (`20260822T020224Z_GEO-18-step2-regression.log`).
> **Carry-forward:** `PORT-9` step 3's mesh prerequisite is discharged — the birdcage has terminals and a port sheet per port; step 3 re-runs unchanged (gates (i)–(iii) never moved). A gapped birdcage without lumped elements still cannot resonate — no port model, solve, impedance or resonance claim here. Ramp example `EX-28` ✅ 2026-08-23.
> Full narrative: `docs/planning/plan-archive.md`, entry «§7 GEO-18 full narrative — archived 2026-08-23 (weekly review)».

**`GEO-19` — `birdcage_port_domain` at `leg_count = 16`, gapped and sheeted: the identity family re-gated and
the cost rung measured** ✅ *(closed 2026-08-25 12:00 slot; step B 2026-08-25 under ruling (6\*), step C under the
10:30 construction-symmetry ruling; commissioned 2026-08-23 weekly review, 32-port directive item (a). Attempts,
rulings (4\*)/(6\*) and the class-gate derivation archived in `docs/planning/plan-archive.md`, entry «§7 GEO-19 full narrative — archived 2026-08-30 (weekly review)».)*
> * **Blockers cleared:** tag collision above 9 legs (upper halves → `200+i`); axis-aligned sheets → each box and
>   sheet built at azimuth 0 and taken to the leg by a snapped `affineTransform`. Step B invariance from `main`:
>   **116 085** sheeted / **114 655** gapped cells, C4 spread **6.050e-16**, terminals 0.988616 × 4, no-gap 98 666
>   (`20260825T003437Z_GEO-19-stepB-invariance-main.log`); the three `PORT-9` modules `19 passed` twice on the
>   mesh-tagged re-records (σ_max 0.999992805, separation 166.6766×, (d0) margin 2256.9707×). Open-limit (1e6 Ω)
>   column retired as a record-bearing fixture — no band widened.
> * **Step C at 16 legs (`20260825T170316Z_GEO-19-stepC-ruled.log` 2 passed / 117 s; `…T170523Z…-ruled-record.log`):**
>   partition / closure / halves / `dx·g` **1.000000000000**, C16 sheet spread **1.331e-15** vs 1e-12, conductor
>   meshed/CAD **0.981503** ≥ 0.95, separation 2.731265e-02 m vs floor 1.750000e-02 m (**1.560723×**); gate (ii) per
>   azimuth class — intra **1.923e-07 / 5.849e-08 / 6.144e-08** vs 1e-6, inter **8.431e-04** vs 5e-3; 4-leg control
>   one class at 3.184e-08, 116 085 cells delta 0. **Cost rung: 116 085 → 307 296 cells (2.6472×), mesh 22.99 → 74.37 s (3.2346×).** * **Carry-forward:** F-small only — **says nothing about F-human**; closed-form layout ceiling **N ≤ 25** on `ring_radius = 0.07` with 14 mm boxes (32 legs need ≥ 0.0876 m) flagged to §10 Phase 6; open-limit conditioning known-issues entry OPEN with (6\*)'s retire-when; thin adjacent/opposite separation flagged for Phase 6.

**`GEO-20` — high-pass birdcage: ring-gap port layout (`ring_gap_length`),
the `GEO-18` pattern on the end rings** 🟡 *(commissioned 2026-08-23 weekly
review — item (b) of the 32-port directive. First at `leg_count = 4`
(8 ring-gap ports), standard tier; the 16-leg / 32-port instantiation is a
second step serial on `GEO-19`. Mesh only, no solve.)* **Do:** opt-in
`ring_gap_length = g` removes the arc `|φ − φ_mid| ≤ g/(2·ring_radius)` at
the mid-azimuth between each adjacent leg pair on **both** end rings
(`2·leg_count` gaps), re-places a port box centred on the gap spanning it
exactly (the `GEO-18` step-1 construction rotated into the ring's local
frame: drive direction azimuthal `φ̂`, terminals = the two planar cut faces
of the ring tube, closed-form area `π·ring_minor_radius²` each), and — with
`emit_port_sheets` — the sheet on the plane through the gap centre normal
to `φ̂`. Leg gaps stay independent (`leg_gap_length` may be `None`: a
high-pass birdcage has uncut legs). **Gate (§4), pre-stated:** per ring
gap, terminal area in [0.95, 1.0] of `2·π·r_ring²` and equal across the
`2·leg_count` ports to 1e-5; closure `(A_cond + A_air + A_phan)/A_box =
1` to 1e-9; sheet meshed/analytic = 1 to 1e-9 with `C_N` spread < 1e-12
(N = `leg_count`) and top/bottom ring mirror symmetry < 1e-12; `GEO-9`
partition < 1e-9; conductor meshed/CAD ≥ 0.95. **Negative controls:**
kwarg off reproduces `EX-21`'s uncut record (98 474 cells, 0.967019)
bit-identically; leg gaps on + ring gaps on coexist with both identity
families holding (the 4-leg fixture becomes a 12-port mesh). **Done-when:**
the gates executed at 4 legs with elapsed recorded, and a `mesh:` example
owed by the daily review's ramp rule; step 2 (16 legs, 32 ports) re-runs
the same gates after `GEO-19` and is priced there.
> * **Step 1 (4 legs, 8 ring ports) ✅ 2026-08-24** *(`tests/mesh/test_birdcage_ring_gaps.py`, `20260824T124525Z_GEO-20-step1-ringgaps.log` 2 passed / 70.4 s and `20260824T124646Z_GEO-20-step1-ringgaps-pass2.log` 5 passed / 158.0 s, standard, `-n 2`, real)*. `ring_gap_length = 8.0e-03 m` ⇒ half-angle **5.714285714e-02 rad**, port width `w = 2·r_ring + 2·clearance` = **1.0e-02 m**; the ring-gapped rung meshes **110 786** cells in 20.9 s, the leg+ring 12-port rung **128 402** in 24.9 s *(corrected 10:30 review audit — 25.2 s was the red `124317Z` run's mesh time; the cited green log reads 24.93 s)*.
>   **The construction, and why the closed forms exist.** The `GEO-18` docstring's "the end-ring alternative gives oblique torus sections at 45 degrees and no closed form at all" is true of an *axis-aligned* box cutting the ring, and that is not what is built. Each ring is `leg_count` partial-torus arcs whose ends are the **radial** half-planes `phi = phi_c ± alpha`, so every cut face is an exact disk of area `pi·r_ring²`. The port solid spanning the gap is the `GEO-18` box **rotated into the gap's own frame**: the wedge `|phi − phi_c| ≤ alpha` ∩ `|z − z_ring| ≤ w/2` ∩ `|u − R| ≤ w/2` with `u = rho·cos(phi − phi_c)`. All six faces are planar, so `V = 2·R·w²·tan(alpha)` = 8.008718871e-07 m³, `A = 2·w²/cos(alpha) + 8·R·w·tan(alpha)` = 5.206757303e-04 m², and the mid-plane section `w²` = 1.0e-04 m² are all exact under a linear mesh — a constant-`rho` face would have turned all three into faceting bands. Corners are evaluated directly in global coordinates rather than built at `phi = 0` and rotated (`GEO-19` ruling (4\*)'s ulp lesson).
>   **Gates, all green twice in-slot.** Terminal **9.796288e-05 m²** = **0.974455** of the closed-form `2·pi·r_ring²` = 1.005309649e-04 m² (band [0.95, 1.0]), the 8 readings taking two values 4.1e-12 apart — spread **≈ 2.1e-08** against the 1e-5 equality gate *(corrected 10:30 review audit: the assertion is green but never prints the spread; the close wrote 2.1e-09, and `std/mean` derived from the printed per-port areas is 2.06e-12 / 9.796e-05 ≈ 2.1e-08 — still three decades inside the gate)*; closure `(A_cond + A_air + A_phan)/A` and port volume/analytic both **1.000000000000**; sheet meshed/analytic **1.000000000000** on 14 facets per port with out-of-plane spread **5.042e-18 m** measured along the sheet's own azimuthal normal; C4 spread and top/bottom mirror on volume and sheet below 1e-12; `GEO-9` partition < 1e-9; conductor meshed/CAD **0.969275** ≥ 0.95; phantom-facing area exactly 0 on all 8.
>   **Negative controls.** (i) Kwarg off reproduces the uncut birdcage — 4 port tags only, **98 666** cells (the 0.11 image's count, ratio 1.001950 against the module's 98 474 record, inside its own 1% band, so nothing is re-recorded here) and meshed/CAD **0.966977** vs `EX-21`'s 0.967019. (ii) Leg gaps + ring gaps together give the **12-port** mesh with *both* identity families exact: leg terminals reproduce `GEO-18` step 1's **0.988616** digit for digit and ring terminals 0.974455, closure and volume 1.000000000000 on all 12. (iii) `GEO-18`'s own two modules re-run green from the same tree.
>   **One finding, measured and not gated.** The union form of the mass identity — gapped CAD conductor = uncut CAD conductor − `2·leg_count·pi·r_ring²·g` — reads **0.999998939803**, 1.06e-06 off, well past the 1e-9 the leg cut achieved. It is *not* the arcs: Pappus on the ring primitives before any boolean reads **1.000000000000** on both the 8 arcs (4.099883683960e-05 m³) and the 2 uncut tori (4.421582771688e-05 m³), so the swept angles are exactly `2·pi/N − g/R`. The residual is OCC's quadrature on a union of 28 vs 20 curved pieces, differenced; the module therefore gates the primitive identity at 1e-9 and records the union ratio. `GEO-18` step 1 hit the same amplification (28× on its own difference) and moved its assertion for the same reason — there the primitive was a cylinder and needed no separate check.
> * **Step 2 attempt 1 (2026-08-28, 04:30 slot) — negative result, parked on `attempt/GEO-20-step2-20260828T094500Z`.** Summary in the table row; full measurement in `attempts.md` 2026-08-28T09:50Z and the known-issues entry of the same date. The branch is **kept**: its module (`tests/mesh/test_birdcage_ring_gaps_scaleup.py`, every band imported, both controls) is the fixture step 2a needs and is not reproduced anywhere on `main`.
> * **Ruled 2026-08-28 10:30 review — step 2a, the discriminator, is queued (§9 item 3); the `src/` fix is not.** The attempt's hypothesis (`_interface_facet_tags` matches a sheet facet by the tags of its two adjacent *owned* cells, so a port whose two half-boxes straddle a rank boundary loses its sheet) predicts one specific, cheap observable: **the set of broken ports moves with the rank count.** Measure that before any `src/` line moves. **Do (no `src/` change):** copy the parked module's test onto the working tree from the branch (`git show attempt/GEO-20-step2-20260828T094500Z:tests/mesh/test_birdcage_ring_gaps_scaleup.py > …` — do not merge the branch), add a per-port print of the rank owning the `PORT_LOWER+i` / `PORT_UPPER+i` cells (rank-local counts `allgather`ed — `cell_tags.values` is rank-local), and run the one test at `-n 4` and `-n 8` with `-s`, each its own footered command. **Anchor (§4), pre-stated:** the broken-port set at each width, printed as a table against the `-n 2` record {P30, P37, P45}; the prediction is a *different* set at each width **and** every broken port having its lower and upper half-box cells owned by different ranks with every intact port on one rank — that ownership-vs-reconstruction agreement across 32 × 3 ports is the quantitative assertion, and a single intact straddling port or a single broken same-rank port **refutes** the hypothesis. **Negative controls:** the 29 intact `1.000000000000` sheet readings at `-n 2` are the baseline the new widths are read against; all 32 port volumes at 1.000000000000 of `2·R·w²·tan α` must reproduce at both new widths (the volume identity does not route through the sheet and must not move). **Cost:** the `-n 2` record was 198 s; mesh time is rank-independent (gmsh on rank 0, ≈72 s), so expect ≈150–200 s per width, ≈ 400 s over two commands at `-k 30 400` each; standard tier, real build. **Traps:** `-s` is mandatory (the module's evidence lives in `print` — the `-n 1` probe carries no numbers for exactly that reason); a broken sheet gives `_sheet_azimuth_deg` a NaN centre — the `_report_safely` guard already converts that into a footered red, keep it; a fresh gmsh model per build. **Scope:** measurement only — no edit to `_interface_facet_tags`, no band, no record, the module lands on `main` only if the review later rules it in (leave it on the branch; commit the ownership print there or as a probe under `scripts/probes/`). **Negative result:** the same three ports broken at every width means the defect is *not* partition-driven and the hypothesis is wrong — record the table in the known-issues entry, chunk stays 🟡, stop; the hypothesis *confirmed* is also stop — the fix touches every sheet-reconstructing module (`GEO-16`, `GEO-18`, `GEO-19`, `PORT-9`'s three, `EX-31`/`EX-33`) and is commissioned as step 2b by a review with an explicit re-record licence, never in-slot.
> * **Step 2a ✅ executed 2026-08-28 (15:00 slot) — the ownership table is measured on the 32-port fixture, the broken set moves with the rank width, and it agrees with rank ownership port for port with no exception at either width.** Logs `20260828T200204Z_GEO-20-step2a-n4.log` and `20260828T200524Z_GEO-20-step2a-n8.log`, Status 1 / **189 s each**, standard tier, real, `-s`, the parked module plus an `allgather`ed per-port count of owned `PORT_LOWER+i` / `PORT_UPPER+i` cells; no `src/` change and nothing landed on `main` but logs and prose. Broken sheets: **{P25, P29, P37, P41, P45}** at `-n 4` and **{P17, P21, P26, P30, P37, P44, P48}** at `-n 8`, against the `-n 2` record {P30, P37, P45} — three different sets, and neither new set is a superset of the old (P45 breaks at 2 and 4, not 8; P30 at 2 and 8, not 4). The pre-stated quantitative assertion holds exactly: at both widths **broken ≡ straddling**, symmetric difference **empty in both directions**, 32 ports × 2 widths, so the refuting observation (an intact straddling port, or a broken same-rank port) did not occur. The 4-leg control is **0 broken / 0 straddling** at both widths. Every negative control reproduces digit for digit: 40/40 port volumes **1.000000000000** of `2·R·w²·tan α`, `GEO-9` partition and air-box closure **1.000000000000**, Pappus **1.000000000000**, kwarg-off 16-leg **307 296** cells / C16 spread **1.331e-15** (ratio 1.000000), 4-leg ring **110 786** cells (ratio 1.000000). Failure shape unchanged: the sheet is lost whole (0 facets) or as a fragment (5 facets / **0.315302109223** at `-n 4`, 6 facets / **0.449137697797** at `-n 8`) while the terminal stays at its intact value and the volume stays exact; closure drops (0.991120008826 / 0.991064589826) only on the one port per run that also loses an **air** facet — the 4-leg-at-`-n 12` signature. **This confirms the phenomenology the ruling asked for, and it is the ghost-layer cause, not the `_interface_facet_tags` one:** a straddling port is exactly a port whose interior sheet facets have a neighbour cell that `GhostMode.none` never materialises. Per the ruling, **confirmed ⇒ stop** — chunk stays 🟡, no band, no record, no `src/` line. The fix and the re-record sweep are **`GEO-24`**'s (already commissioned by the human operator); with this table `GEO-24` step 1 owes only the before/after readings of the modules on `main`, and `GEO-20` step 2 becomes a re-run once `GEO-24` step 2 is green. The instrumented module lands as `scripts/probes/geo20_step2a_ownership_scaleup.py` — the §7 entry's stated alternative to the branch, taken because a concurrent session held `examples/ports/05_*` modified mid-slot and that path differs between `main` and the branch, so a `git checkout` of the branch would have refused; the branch is unchanged and nothing under `tests/` moved. See the attempts.md entry of 2026-08-28T20:40Z for that anomaly.
> * **Step 2 attempt 2 ✅ executed 2026-08-29 (09:00 slot) — the re-run on the plumbed tree is green at `-n 2` and `-n 12`; `GEO-20` is ✅ and the 32-port ring-gap layout is built and gated.** Two windows, **372 s** total, standard tier, real, `-s`, `-k 30 570`, no `src/` change in the slot: `20260829T140037Z_GEO-20-step2-rerun-n2.log` (Status 0 / 188 s) and `20260829T140402Z_GEO-20-step2-rerun-n12.log` (Status 0 / 184 s), each `1 passed`. The module came back exactly as parked — `git checkout 31c08ed -- tests/mesh/test_birdcage_ring_gaps_scaleup.py`, the single file on the attempt branch, +569 lines; the birdcage generator's signature had not moved and no import failed. **The quantitative assertion (§4) is the 32-fold reconstruction identity family**, all read at both widths: 32/32 sheets meshed/analytic **1.000000000000** of `w²`, 32/32 boundary closures **1.000000000000** against 1e-9, 32/32 `volume/analytic` **1.000000000000** of the analytic wedge `2·R·w²·tan α`, plus `GEO-9` partition, air-box closure and the Pappus ring-arc check all at 1.000000000000 and the ruled per-azimuth-class terminal bands (intra 4.198e-08 / 4.498e-07 / 4.681e-07 / 8.997e-07 against 1e-6; inter 3.315e-07 against 5e-3). **Attempt 1's failure is repaired exactly where the `GhostMode.none` diagnosis said it would be:** P30 and P37 go 0 → **176 air facets**, P45 goes 5 facets / 0.315302109223 → **180 / 1.000000000000**, on the identical **265 621**-cell geometry — so `GEO-24`'s one-keyword plumb, not any change in this module, is what made the fixture readable in parallel. Both in-run controls reproduce (307 296 and 110 786 cells, ratio 1.000000 each); the single moving digit anywhere is the C16 sheet spread at the 1e-15 floor (1.331e-15 → 1.210e-15), the same pair `GEO-24` step 1a reported. The cost rung re-measures at 2.3976× cells and ≈3.15× mesh time, unchanged from attempt 1. Branch `attempt/GEO-20-step2-20260828T094500Z` **deleted** with this commit, per its disposition ("delete it with the commit that lands or retires the module"). The known-issues entry was already retired by `GEO-24` step 2b; it gains only a confirmation line, no re-opening.

**`GEO-21` — dispose of the red `GEO-15` graded-conductor gate on 0.11: the ungraded baseline no longer
meshes** ✅ *2026-08-26 (step 1 00:00 slot, ruled 03:00 review option (b), step 2 04:30 slot; commissioned
2026-08-25 18:00 review from `EX-30` leg (mesh). Full narrative archived in `docs/planning/plan-archive.md`, entry «§7 GEO-21 full narrative — archived 2026-08-30 (weekly review)».)*
`GEO-15`'s ✅ is the 0.7.2 close and stands; the `h_c = None` control aborted in gmsh ("overlapping facets") on 0.11.
> * **Step 1 (measure-first):** red reproduced `20260826T050100Z_GEO-21-step1-red-repro.log`; CAD-mass probe
>   `…050134Z_…-cad-mass-probe.log` — 3.2e-3 → **0.916742** (47 975 cells), 1.6e-3 → **0.966977** (98 666), CAD mass
>   1.030097043e-04 m³; control ladder `…050319Z_…-control-ladder.log` (`-n 1`) — 4.8e-3 → **0.846150** (33 185),
>   6.4e-3 → 0.767219 (27 912), 9.6e-3 → FAIL. Neither pre-stated branch fired (0.916742 fails the module's own
>   `CAD_MASS_GATE − 0.05` = 0.90 separation guard); ruling requested.
> * **Ruling (b):** control = `h_c = 4.8e-3` (buffer rung away from the meshability cliff); the claim is **demoted to
>   fine vs coarse grading** — "grading required" closed on 0.7.2 and stays there — stated in module, test, `mesh:3`
>   docstring and guide. **Step 2:** gate `20260826T093202Z_GEO-21-step2-gate.log` `1 passed in 41.11s`, 0.846150 /
>   0.916742 / 0.966977 asserted at the unmoved ≥ 0.95 gate and 0.90 guard (cleared by 0.0538); `mesh:3`
>   `…093403Z_…-mesh3.log` separation 0.120826; docrefs `…093552Z` `stale=13`, **`meshing` 2 → 0**; two latent bugs fixed.
> * **Carry-forward (OPEN):** known-issues entry re-headed on the generator-continuum finding — coarse conductor sizings, `None` included, cannot mesh on 0.11 — uncommissioned; `GEO-23`'s fourth site (`birdcage_port_domain`) lands here.

**`GEO-22` — `straight_wire_domain`: bisect the coarse-resolution floor and guard it** ✅ *(step 1 2026-08-28
07:30 slot; step 2 2026-08-29 07:30; demoted ✅ → 🧪 by the 08-29 10:30 review under §4 clause 3 and restored by
step 2c, 2026-08-29 16:30. Commissioned 2026-08-26 18:00 review as owner of `EX-30` leg (root)'s finding. Full
narrative archived in `docs/planning/plan-archive.md`, entry «§7 GEO-22 full narrative — archived 2026-08-30 (weekly review)».)*
> * **Step 1 — measured negative:** the 2.5e-4 sweep of [0.008, 0.010] is **non-monotone on both geometries**
>   (example FAIL at 0.00875 / 0.00925 / 0.00975 / 0.01000, gate FAIL at 0.00875 / 0.00900 / 0.01000) and bit-identical
>   twice (`20260828T123115Z_GEO-22-step1-bisect.log`, `…123205Z_…-repeat.log`); 21 830 cells at 0.008. **No
>   `RESOLUTION_FLOOR` exists and none was written; `src/` untouched.** Ruling 08-28 10:30: raise-path wrap via
>   `GEO-23` step 2a (gate `test_geometry_failure_is_collective.py`, `20260828T170254Z`), allowlist rejected.
> * **Step 2 — size-field probe (leg D, no `src/`):** gmsh `Distance`/`Threshold` field on the wire surface reads
>   **18/18 OK, 0/18 fallbacks** vs leg C's 7 FAIL / 18 fallbacks (`20260829T123331Z_GEO-22-step2-sizefield.log`, 33 s);
>   leg C control bit-identical in its own process (`…123413Z_…-legC-control.log`); gate re-run 1 passed (`…123459Z`).
> * **Step 2c — the assert:** `tests/mesh/test_straight_wire_size_field_probe.py` patched **19 823 cells / 0
>   fallbacks**, unpatched **21 830 / 1**, **0.00% vs ±1%** at `-n 1` (`20260829T213132Z_…-step2c-n1.log`, 8 s) and `-n 2`
>   (`…213148Z_…-n2.log`, 7 s); perturbed reference reds on both ranks (`…213211Z`), restore green (`…213227Z`). * **Carry-forward (OPEN):** landing the size field would move `mag:1`'s 21 830 and the three ladder records — **that licence is the weekly review's**, and the known-issues entry stays OPEN on exactly that decision.

**`GEO-23` — the 0.11 "overlapping facets" family: classify, ladder, own** ✅ *(step 1 2026-08-28 09:00 slot;
ruled 10:30 review; step 2a 12:00, step 2b 13:30 ⇒ chunk closed 2026-08-28; commissioned 2026-08-27 03:00 review
from `OPS-26` step 2's four sites. Full narrative and both step-2 rubrics archived in `docs/planning/plan-archive.md`, entry «§7 GEO-23 full narrative — archived 2026-08-30 (weekly review)».)*
> * **Step 1 (8 windows, 318 s, `20260828T140041Z…141217Z`):** **all four sites fail at `-n 1`** — nothing is
>   partition-dependent; the two "rank-divergent" observations were log-interleave artifacts (withdrawn by measurement).
>   The deadlock is a **raise-path** property (rows 1/3/4 Status 124 at 120–121 s; `birdcage_port_domain`, already
>   wrapped, footers Status 1 in 5 s). Ladders monotone: `cylindrical_domain` 0.032 → 1 213 cells, `coil_phantom_domain`
>   0.024 → 5 464. In-process ladders read gmsh contamination (`IndexError`), not geometry. Dead `test_cylindrical_domain.py`
>   → one asserting test at 1e-9; `tests/mesh` collects 58.
> * **Step 2a (raise path, 12 windows / 72 s):** `_raise_geometry_failure_on_every_rank` (`io/mesh.py:30–61`) wraps
>   `straight_wire` / `cylindrical` / `coil_phantom` (+76 lines, zero deletions); the three deadlocking rows footer
>   **Status 1 in 2–3 s** at `-n 2` (`20260828T170311Z/170331Z/170347Z`, ≈ 50× cheaper); `mag:1` 21 830 cells and
>   `B(3 mm) = 6.666667e-05 T` unmoved (`…170414Z`); `GEO-21` control 1 passed / 36.76 s.
> * **Step 2b (sizing, 8 windows / 40 s):** `0.04 → 0.032` and `0.03 → 0.024` at the three call sites; all green at `-n 1`
>   and `-n 2` (`20260828T183106Z…183242Z`) with counts **1213 / 5464 exact (0.00% vs ±1%)**; no physics assertion moved; the three known-issues entries RETIRED. * **Carry-forward:** the fourth site, `test_birdcage_volumes_partition_the_box` on `birdcage_port_domain`, is **still red** — it is `GEO-21`'s open generator-continuum entry; re-home or reopen as step 3 is a review's call.

### TH — Time-harmonic Maxwell (Phase 2)

| ID | Title | Status | Tier |
|---|---|---|---|
| `TH-1` | **Real complex time-harmonic formulation** | ✅ | standard |
| `TH-2` | *(history and §9 briefs: `docs/planning/chunks/TH-2.md`)* Time-harmonic API hardening. ✅ **Audit PASS 2026-10-02 10:30 review (chunk file).** **Closed 2026-10-02 07:30 slot (§9 item 153):** W1 complex `-n 2` 26 passed in 118 s, with `MAT-1` (P1)/(P2) at 0.0 on 74 020 dofs, (c1) 1.11e-16, (c2) 2.406e-14, (P3) 3.536 % and presets 3.523 / 3.548 / 3.557 % reproduced (`20261002T124713Z_TH-2-W1.log:499–512`). The real-build W2 gave 20 passed / 6 skipped (complex-only), raise tests green (`…124920Z_TH-2-W2-real.log:571`). **Done-when (written 2026-10-01 10:30 review, rule (j); §9 item 153):** one harness window, `-n 2`, complex (`FEM_EM_REQUIRE_COMPLEX=1`, `tests/environment` first), green on `test_material_map_rank_safety.py` (`OPS-13`: 1/162 and ∫σ = σ·V to 1e-12 on every rank; absent tag rejected on every rank), `test_phantom_material_model.py` (`MAT-1` (P1)/(P2) bitwise, (c1)/(c2) > 0, (P3) and three presets < 10 %) and the three raise-before-solve tests (`test_time_harmonic_smoke.py:625–696`); a real-build window recorded beside it ⇒ ⚠️ → ✅. Not claimed: the frequency-range and phantom-frequency guards (`core/time_harmonic.py:447–466`) are untested. *(Replaces the 2026-09-30 03:00 review's "retire or merge" note — the row's surface is gated by `OPS-13` / `MAT-1`, and this window executes those gates under this row.)* | ✅ | standard |
| `TH-3` | Boundary-condition option set. **Done-when (written 2026-09-29 18:00 review, rule (j); queued §9 lane 1).** On [0, 0.2]³ at 127.74 MHz, σ 0.7, εᵣ 78, degree 1: zero-data PEC constrains exactly 18n² owned dofs (2 592 / 10 368 at 12³ / 24³) and `NATURAL` 0. E_pec = (sin ky sin kz, sin kz sin kx, sin kx sin ky) under PEC and E_nat = (sin kx cos ky, −cos kx sin ky, 0) under `NATURAL` (k = π/L), each driven by −jωμ₀J = (2k² − k₀²ε_c)E with `project_source=False`, give relative L2 < 0.10 at 24³ and rate > 0.9 over 12³ → 24³ (`TH-1`'s bands, `test_time_harmonic_mms.py:187–197`). Each mode-swapped solve misses its field by > 0.10. `TH-1`'s MMS is green in the same window, and the placeholder PEC test is replaced. Harness `-n 2`, complex, elapsed recorded ⇒ ✅. No new enum member (⚠️ no-features rule). | ✅ *(**Closed 2026-09-30 (21:00 CDT slot take-next, §9 item 89):** zero-data PEC constrains exactly 18n² owned dofs (2 592 / 10 368 at 12³ / 24³) and `NATURAL` 0; E_pec under PEC reads relative L2 1.166e-1 → 5.852e-2 (rate 0.9942) and E_nat under `NATURAL` 7.528e-2 → 3.774e-2 (rate 0.9961), the swapped controls missing at 0.891 / 0.626; `TH-1`'s MMS premise green in the same window (16³ 5.659e-2, rate 0.9929 — its first 0.11-image run); 21 passed, 83 s, `-n 2` complex (`20260930T022051Z_TH-3-W1.log:90–99,280–297`). Tests only; the placeholder PEC test is replaced.)* | standard |
| `TH-4` | Convergence/conditioning diagnostics — **Done-when (written 2026-09-28 18:00 review, lane 1, §9 item 67):** on `tests/solver/test_convergence_diagnostics.py`'s fixture, (i) `residual_norm` equals the last `residual_history` entry (rel ≤ 1e-14), (ii) `converged_reason == 2` with h[-1] ≤ 1e-8·h[0] < h[-2], (iii) the explicitly recomputed preconditioned residual ‖B(b − Ax)‖/h[0] ≤ 1e-6, and the (1 + δ)x control ≥ 1e-5, all green through the harness, complex, `-n 2`, elapsed recorded ⇒ ✅. It relies on the `TH-2` ⚠️ API (`TH-3` ✅ 2026-09-30) and default boundary condition and does not revalidate them; "conditioning" has no landed code and is outside the ✅. **Closed 2026-09-28 (22:30 take-next):** PRECONDITIONED / LEFT as predicted; 1409 its; (i) rel 0.0; (ii) reason 2, h[-1]/h[0] 9.865e-9 < 1e-8 < h[-2]/h[0] 1.005e-8; (iii) explicit ‖B(b − Ax)‖/h[0] 9.865e-9 (= h[-1] to 1.000000); control 1.000e-4 (pred. 1.0e-4); 27 passed, 38 s, `20260929T034518Z_TH-4.log`. | ✅ | standard |
| `TH-5` | Absorbing boundary condition (ABC) | ✅ on step 2c and the §2.1 line (the first-order ABC's closed forms, 2a's (S) standing and (R) retired); step 3 (2026-10-22, §9 item 596) adds the PEC-backed graded PML on 1a's TEM line as a record on the optional step, at normal incidence; step 3b (2026-10-25, §9 item 628) adds the same PML at 60° oblique incidence, degree 2; step 3c (45°, the closed leg): item 644 BLOCKED 2026-10-26 at its 300 s ceiling, re-posed as §9 item 649 (attempt 2 of 4). **Step 3c green 2026-10-26 (item 649; registered record): the angular law at 45° on the closed leg.** At degree 2 and λ₀/40, \|Γ_h − Γ_c(45°)\| = 8.89e-8 against `PML_REFLECTION_BAR` 1e-4 (Jin 3e (9.133) over (9.165)), and the normal-incidence law misses by 6.87e-3 (68.7× the bar, asserted > 10×). Degree 1 3.8e-5 and (Cabc) ratio 22.7 printed (`20261026T094918Z_TH-5-step3c-W1b.log`, 109 s; 3b regression `20261026T095117Z_TH-5-step3c-Wreg.log`). *History: `docs/planning/chunks/TH-5.md`.* | heavy (relabelled 2026-10-08 18:00 review: step 2c's W1 measured 225 s) |
| `TH-6` | **Validation: plane wave in lossy half-space** | ✅ | standard |
| `TH-7` | **Validation: waveguide cutoff / coaxial line** | ✅ | standard |
| `TH-8` | **Validation: sphere in uniform field (quasi-static)** | ✅ | standard |
| `TH-9` | **Validation: PEC rectangular-cavity resonances** | ✅ *(2026-09-29, 06:00 slot, §9 item 74: `analytic_cavity_frequencies` now counts degenerate pairs through `analytic_box_cavity_modes`. The six solver modes match the closed form with the TE₁₁₁ / TM₁₁₁ pair counted twice, worst 0.0517 % against the imported 1 % band. The three original tests are unchanged and green (0.0436 % / 3.85 / 8 of 8). The frozen old helper is kept as the control. `20260929T112910Z_TH-9-KI-W2.log`; the W1 red on old `src/` is `…112807Z_TH-9-KI-W1.log`.)* | standard |
| `TH-10` | **Validation: lossy dielectric sphere in a full-wave field at 64/128 MHz (the first Larmor-regime gate)** | ✅ — **128 MHz re-recorded on the 0.11 image, `OPS-18` step 3 attempt 1 / step 3b (2026-08-23):** interior relL2 **1.826% → 1.769%**, *with* its mesh moving **55 251 → 55 241 cells** under the image's gmsh 4.15.2. The 64 MHz gate is **bit-identical** across the image change (3.643%, power 3.629%), so this is the same new-gmsh mesh drift as the `OPS-17` volume record, disposed the same way — a re-record by measurement under class ruling (1\*), no band and no gated physics claim moved (the §2.1 / §2.2 `TH-10` numbers and `EX-19`'s reproduction stand). | standard |
| `TH-11` | Coil-loading trend across the eddy→displacement transition (`MAT-6`'s ΔR machinery at rising f) | ✅ Closed on step 4's answer plus step 5's measured negative, on the `GEO-14` precedent: the degree-1 coil-loading ladder hits a superlinear memory wall. Any reopening needs an `xl` slot and a review re-pricing it. *History: `docs/planning/chunks/TH-11.md`.* | standard (steps 4–5 heavy; step 5d **`xl`**) |
| `TH-12` | Second-order elements (degree-2 N1curl): accuracy-per-DOF and cost, measured (operator directive 2026-08-18; decides the production element order for §10 Phase … | ✅ Closed at the 2026-09-02 weekly review, with the decision clause on production element order re-affirmed. Degree 2 buys accuracy per cell at a measured cost, and the coil rung meets the same memory wall. *History: `docs/planning/chunks/TH-12.md`.* | standard (step 2 heavy) |
| `TH-13` | The degree-2 gradient-subspace injector: feed model or any `W_m ≫ W_e` fixture? — the discriminator `TH-12` step 3 named (commissioned 2026-08-23 weekly review; cheap fixtures only) | ✅ *(closed 2026-08-31 on step 2 — (A) holds at round-off, the injector is the degree-1-only `H¹₀`-only source projection; audited COMPLIANT 03:00 review. Follow-on **step 3a** — the matched projection, opt-in, loop fixture only — scoped 03:00 review, §9 item 1, **executed 2026-08-31 04:30 slot and 🟡**: both anchors met with 6–14 orders of margin (residue 1.298386e-02 / 1.045186e-01 → **8.109635e-17 / 1.790460e-16** vs ≤ 1e-8; gradient share of `W_e` 99.98% / 99.9997% → **4.6e-23 / 3.1e-21** vs ≤ 1e-6; `W_e` to 0.018% / 2.6e-4 % of record vs ≤ 2% / ≤ 1%), default path bit-identical on control (b) at **0.000e+00**, but one of the two owed regression re-runs — `test_coil_loading_degree2.py` — **could not be executed** (exit 124 at 571 s, twice — step 3a″ on 2026-08-31 measured the cost as the degree-2 pair alone, ≥ 524 s, mesh 4.3 s; known-issues entry), so the coil identity tests' 1e-9 reds are unverified on this commit — the degree-1 half of the owed claim was re-observed by 3a″ (+0.00039 pp of record), and **step 3a‴ (module split, one σ-half per window) closed that gap on 2026-09-01**: the two reds are now **observed** at 3.8990e-09 / 3.7235e-09 against the unloosened 1e-9, one per half, each window footered inside its 600 s ceiling — see the step-3a‴ bullet; the lumped-sheet coil drive is `project_source=False` and out of 3a's reach, see entry)* | standard |
| `TH-15` | Internal perfect-electric-conductor bodies | ✅ **2026-10-19 weekly (row written 2026-10-20 03:00 review; audit in `TH-15.md`).** Steps 1, 2a, 2d, 3, 3a–3c green with asserted identities: the PEC sphere as a hole against its closed form (β 1.97 % vs 4.886 %, `20260905T170225Z_TH-15.log:290–296`); the birdcage 4×4 as a PEC hole through `PORT-9` / `PORT-11`'s imported gates at 10 / 64 / 128 MHz (`20260913T202954Z_TH-15.log:951–954, 2768`, `…202718Z_TH-15.log:951–954`, `…202828Z_TH-15.log:956–959`); the terminal excess attributed to the sheets' Cauchy–Schwarz deficit (`20260914T005452Z_TH-15-step3c.log:966–970`). Step 2 killed (2026-10-07 epitaph); step 3's default closed without a flip (2026-10-19 epitaph). **Done-when (as amended by the 2026-10-07 and 2026-10-19 weeklies):** steps 1, 2a, 2d, 3 and 3a–3c recorded green, the PEC-sphere closed form and the hole route's discrete power identity asserted, cell counts and elapsed times recorded, §2.1's conductor-model line stating what is gated (PEC interior bodies: sphere, two-torus mesh, birdcage) and what is not (no loss, no Q, no copper — `TH-14`'s) ⇒ ✅. No band moved; no new window. *History: `docs/planning/chunks/TH-15.md`.* | standard (step 3 heavy) |
| `TH-14` | Surface-impedance (Leontovich) boundary on conductor surfaces | ✅ Done: the Leontovich surface-impedance boundary passed the lossy-cavity Q_c anchor, the copper F-small birdcage port gates at 10 / 64 / 128 MHz, and on 2026-09-14 the Dodd–Deeds copper-slab floor step that re-closed the row after the audit demotion. The §2.1 conductor-model line is written. *History: `docs/planning/chunks/TH-14.md`.* | standard (step 3 heavy; slab step heavy, 340 s) |
| `TH-16` | Symmetry planes: per-face PEC / PMC on cut faces with port rescaling | ✅ Done 2026-09-30 (§9 item 93) with audit PASS: the quarter birdcage's three solves reproduce the 4×4 of its own mirror-built full mesh to 1.520e-12, with the PMC-for-PEC negative control asserted. Scope is one fixture at 10 MHz and degree 1, and the re-anchored gate was ratified by the 2026-09-30 weekly. *History: `docs/planning/chunks/TH-16.md`.* | standard |
| `TH-17` | Birdcage eigenmodes | ✅ Done 2026-10-01 (§9 item 131): gates (c) and (d′) passed, completing the amended Done-when together with (a), (b) and (e), and the 10:30 review audited it PASS. Gate (a) is closed as re-posed at `C_mode` rather than `C_tuned`, a reading the weekly confirms. **Done-when:** (a)–(e) through the harness at `-n 4`, the fixed-point shift recorded (< 0.1 % or the band re-read), elapsed recorded ⇒ **`TH-17` → ✅ and Phase 6's first physics target ticked.** *History: `docs/planning/chunks/TH-17.md`.* | heavy *(step 1 probe measured 47.8 s / 2.18 GiB at `-n 4`)* |
| `TH-18` | **Layered impedance boundary** — thin copper foil on a substrate (HFSS *Layered Impedance*), the construction of real coils; a refinement of `TH-14` — **feature ladder B4** (operator directive 2026-09-04; serial on `TH-14`) | ✅ **The σ = 10 S/m amendment is ratified by the 2026-09-30 weekly (§10 rulings).** *(**Closed 2026-09-30 (21:00 CDT slot take-next, §9 item 87; ladder B4):** on a driven TEM line (degree 2, 1440 cells) with a volumetric εr = 4 slab and the σ = 10 S/m foil as the only surface term, Γ at 127.74 MHz meets the closed-form stack for `t5` / `thin` / `pec` to d\|Γ\| ≤ 1.2e-6 and d∠Γ ≤ 1.42e-5 rad, inside B = (3e-6, 5e-5 rad) frozen from the `leo` premise 7.2e-7 / 1.35e-5 (`20260930T021219Z_TH-18-step2-W1.log:89,91`); negative control 0.092963 against 0.093 ± B and > 10B; copper 35 µm printed \|Γ − Γ_pec\| 5.24e-5; 21 passed, 54 s, `-n 2` complex (`20260930T021508Z_TH-18-step2-W2.log:100–116,189`). Tests only.)* *(step 1 green 2026-09-29, item 82: Done-when's step-1 limits and cavity ratios hold in one window with `TH-14` re-run)* | standard |
| `TH-20` | Symmetry planes at human scale: the quarter F-human (§10 chain H13) | ✅ Step 1's degree-1 quarter reassembly is asserted (2.106e-11 against 1e-9, item 233) and step 2a priced the degree-2 quarter (684 904 unknowns, ≤ 19.57 GiB at `-n 8`, item 350); step 2b-a (item 412) put the `FEM_EM_TH20_DEGREE` knob on step 1's module with (Deg) exact on the four sectors and the mirror-built mesh. ✅ **2026-10-20 (`xl` entry 37, row written by the 03:00 review; audit in `TH-20.md`):** step 2b's (A) **1.739e-11** against the imported 1e-9 at degree 2, 14 passed on all 16 ranks, Status 0, 936 s, 90.24 GiB (`20261020T070008Z_TH-20-step2b-H13b.log:10476`, `:11610–11615`); entry 37's rule fires its first branch — the quarter is licensed as the degree-2 human-scale route. C1 / C1x / C2 were printed, not asserted (the readout's parenthetical, ruled non-gating). **Done-when:** step 1's identity asserted in one logged window ⇒ 🟡; step 2a priced ⇒ record; step 2b's (A) — the degree-2 quarter reassembly against a direct degree-2 column — asserted ⇒ ✅. *History: `docs/planning/chunks/TH-20.md`.* | heavy |
| `TH-19` | Carry the matched source projection to the coil, and retest the production element order | ✅ **2026-10-17 weekly (row written 2026-10-18 03:00 review; audit below).** Steps 1–3 green with logged, asserted identities; step 3 holds both degree-2 power identities on the sheet-driven 4-leg birdcage at 10 / 128 MHz — (a) 4.6e-15 / 8.9e-15 ≤ `DISCRETE_IDENTITY_RTOL` 1e-6, (b) 6.8e-11 / 1.8e-12 ≤ `IDENTITY_TOLERANCE` 1e-9 (`20260913T183446Z_TH-19-step3-10MHz.log`, `…183723Z_…-128MHz.log`, 134 / 118 s at `-n 8`); outcome (a) accepted with no default change, and the production-order question is answered by the order rule (degree 2 default, 2026-10-01). **Done-when (2026-10-17 weekly, the banked draft with "⇒ ✅" and no new window, so the draft's "re-observed on the current image" clause is not carried, `TH-19.md` 2026-10-18 ruling):** steps 1–3 recorded green; step 3's (a) ≤ `DISCRETE_IDENTITY_RTOL` and (b) ≤ `IDENTITY_TOLERANCE` at 10 and 128 MHz; the production-order question answered by the order rule in force ⇒ ✅. The coil volume-current default-path red (known-issues) is not part of it. Audit PASS 2026-10-18 03:00 after the tier cell was corrected. *History: `docs/planning/chunks/TH-19.md`.* | heavy (steps 1–2 measured 475–484 s at `-n 8`; step 3 134 / 118 s) |

**`TH-10` — lossy dielectric sphere, full-wave, 64/128 MHz (Larmor gate)**
✅ *(steps 1–4 ✅ 2026-08-13, chunk closed by the 10:30 review; full step
narratives archived in `docs/planning/plan-archive.md`)*. The Mie-series
anchor lives in `utils/analytical.py` (`LossySphereSeries`, `e^{+jωt}` by
conjugating both `ε_c` and the field; six identity gates incl. a
conjugated-convention control at 2.1e+04× separation). Gated against the
series on an imposed total-field wall (the `TH-8` pattern): interior relL2
**3.643% at 64 MHz** (18.68× closer to the series than to quasi-statics)
and **1.826% at 128 MHz** (57.31×), both under 5% and decreasing; the
SAR-relevant ½∫σ|E|² to **3.629%** at 64 MHz with the quasi-static power
route missing by 58.14% (quasi-statics under-predicts absorbed power 2.4×
at 1.5 T). Monotonicity of the power error is asserted in-file since
2026-08-15 (8.387% → 3.629%, digits bit-identical to the record).
`GEO-14`'s discriminator later showed the ≈3% residual is mesh resolution
still converging (1.781% at 55 251 cells, rate 1.77 in h), not a faceting
floor — frequency-independent *and* mesh-limited. Standing caveats: the
step-4 negative-control margin is 1.16× (a fixture property); the series
reference vs quasi-statics separation is 55–69% in relL2 vs 102–155% in
max-norm — quote the norm with the number. Gates the volume integral only
— no mass averaging, no C95.3 wording, no coil. The coil-loading trend was
commissioned as `TH-11`.

**`TH-11` — coil-loading trend across the eddy→displacement transition** ✅ *(closed 2026-08-18, 10:30 review, by adjudication — step 5 a measured negative; commissioned 2026-08-13 10:30 review to validate `MAT-6`'s ΔR machinery at rising f.)*
> **⚠️ Memory-premise caveat (2026-08-24).** The container ceiling was raised **64 GiB → 128 GiB** by operator directive (§5.1; `/sys/fs/cgroup/memory.max` verified at 128.0 GiB). Step 5's negative — 2.81 M cells OOM, 0.99 M pegged at 64.00 GiB — was measured against the *old* wall, and both rungs sat inside 2× it. That makes the affordability claim **unmeasured, not false** — it is not evidence about this box any more, and must not be cited without re-measuring. **Nothing here is reopened by this note**: re-pricing a third rung is a **review** decision, and step 4's finding is the load-bearing one — the apparent trend was the resolution term, flat in f, and more RAM does not turn a resolution artefact into a trend. A revival needs a new argument *and* its finest rung priced first (§10 epitaph).
> * **Step 1 ✅ 2026-08-13** (`tests/validation/test_coil_loading_larmor_probe.py`, `20260814T003445Z_TH-11-step1-larmor-n2.log`, 138 619 cells, solves 30.5 + 27.0 s at `-n 2`): 64 MHz ΔR deviation **+10.2698%** vs Dodd–Deeds (1.5834% at 10 MHz) at 1.26 cells/δ; complex-power identities ≈1e-14, σ = 0 control exact.
> * **Step 2 ✅ 2026-08-15** (`test_coil_loading_larmor_resolution.py`, `20260816T003251Z_TH-11-step2-resolution-n2.log`, 390.9 s, 417 914 cells): `near` = 0.0025 reads **+2.8063%**, a −7.4635 pp move — RESOLUTION-DOMINATED.
> * **Step 3 ✅ 2026-08-16** (`test_coil_loading_transition_30mhz.py`, `20260816T183310Z_TH-11-step3-30mhz-n2.log`, 10 passed, 70.3 s): 30 MHz reads **+5.5912%**, ΔX ratio 0.9500; identity 2.7373e-14 / 1.6799e-14 (bound 1e-9); confound monotone (cells/δ 3.18 → 1.84 → 1.26).
> * **Step 4 ✅ 2026-08-17** (`test_coil_loading_richardson_ladder.py`; `20260817T033320Z_TH-11-step4-baseline.log` 138 s, `20260817T033547Z_TH-11-step4-fine-10mhz.log` 422 s, `20260817T034258Z_TH-11-step4-fine-30mhz.log` 383 s; heavy, `-n 2`): **fixed-f h-ladders are flat in f** — refining 0.005 → 0.0025 moves **+1.5834% → −0.2829%** (10 MHz) and **+5.5912% → +1.1119%** (30 MHz); h → 0 brackets **[−2.1492%, −0.9050%]** at 10 MHz, **[−3.3675%, −0.3812%]** at 30 MHz (p = 1 / p = 2), overlapping, both straddling ≈−1%. Identity ≤ 8.1597e-14 on six solves; baseline anchors reproduced to −0.00002 / −0.00000 pp (floor 0.01 pp). The apparent trend was the resolution term.
> * **Step 5 (64 MHz third rung, `near` = 0.00125) — measured negative.** Probe `20260817T123353Z_TH-11-step5-probe.log` (exit 124): 2 807 309 cells, 5.03 cells/δ, mesh 288.2 s. **5a ✅ 2026-08-17:** XDMF cache exact (`20260817T183751Z_TH-11-step5a-cache-third.log`, 2 807 309 cells, 126.4 s mesh, 14.8 s read-back; smoke `20260817T183709Z_TH-11-step5a-cache-smoke2.log`, `20260817T183248Z_TH-11-step5a-cache-smoke.log` exit 124); rank control `20260817T184026Z_TH-11-step5a-rank-control.log` (`-n 8`, 174 s) reproduces **+2.8063%** to **+0.00002 pp** against the 0.1 pp band; fine-rung solve 72–73 s at `-n 8`. **5b 🟡×2:** `-n 12` OOM-killed with the container (`20260818T003806Z_TH-11-step5b-third-loaded.log`, exit 137 at 518 s); `-n 8` `memory.peak` **68 719 480 832 B = 64.00 GiB vs `memory.max` 68 719 476 736 B**, timed out (`20260818T020143Z_TH-11-step5b-third-loaded-n8.log`, exit 137, 908 s); loaded/free split exact (`20260818T003418Z_TH-11-step5b-rehearsal.log`, 288 s, `DRIVE_SCALAR_BAND` 1e-12 measured 0.000e+00). **5c 🟡:** `near` = 0.0018 → **994 258 cells** (`20260818T093219Z_TH-11-step5c-cache.log`); loaded solve completed at `-n 8` in 320.5 s, ΔR +1.3628036e+00 Ω, but `memory.peak` 64.00 GiB = 100.0% of cap (`20260818T093314Z_TH-11-step5c-loaded-n8.log`); free solve exit 124 at 479.2 s (`20260818T093919Z_TH-11-step5c-free-ladder-n8.log`).
> **Adjudication (2026-08-18 10:30 review), relied on by §2.2:** the memory wall is superlinear in cells (0.42 M comfortable / 0.99 M pegged at 64.00 GiB / 2.81 M OOM — MUMPS fill-in); every degree-1 third rung is unaffordable or statistically useless (ratio ≈ 1.2 vs 0.01 pp floor); no 5d. The chunk's question is answered by step 4 (resolution, not physics; brackets flat in f at 10/30 MHz); **64 MHz has no h → 0 bracket** and the gated-bracket deliverable transfers to `TH-12` step 2. §2.2's extrapolation bullet carries the negative and moves only on a gated 64 MHz bracket. Branches `attempt/TH-11-step5b-20260818T024200Z`, `attempt/TH-11-step5c-20260818T101500Z` deleted; the non-uniform three-rung fit `(d_c − d_m)/(d_m − d_f) = (h_c^p − h_m^p)/(h_m^p − h_f^p)` (p by bisection) is recorded as a formula only, never exercised on data — re-implement, do not inherit.
> **Standing cautions:** ≈480 s is the safe container ceiling for a foreground slot (560 s exceeded the 660 s Bash wall under memory pressure); `memory.peak` is a lifetime high-water mark (see `TH-12` step 1).
> Full narrative: `docs/planning/plan-archive.md`, entry «§7 TH-11 full narrative — archived 2026-08-23 (weekly review)».

**`TH-12` — second-order elements (degree-2 N1curl): accuracy-per-DOF and
direct-solver cost, measured on gated fixtures** ⬜ *(commissioned
2026-08-18, operator directive — element order is to be evaluated as a
cross-phase lever, and the production element order for the §10 Phase-5/6
work is to be decided from this chunk's measurements, not assumed.)*
Motivation, all on record: every open accuracy question in the TH/MAT
lineage is a cells-per-skin-depth question (`TH-11` steps 2/4: the
apparent frequency trend was the resolution term), and `TH-11` step 5b
has now measured that the degree-1 route to the 64 MHz h → 0 bracket
**does not fit the box** — 2 807 309 cells OOM at every legal rank count
(64.00 GiB = `memory.max`; **⚠️ caveat 2026-08-24 — the ceiling is now 128
GiB by operator directive, see §5.1, so that OOM is a statement about the
old box and this chunk's step 2 must re-price rather than inherit it. The
degree-2 wall this chunk was scoped around, 61.94 GiB on the *coarse*
rung, was 96.8% of the old cap and is ≈48% of the new one**). Degree 2 is the other axis: ≈20 DOFs/tet vs
6, denser MUMPS blocks, but second-order field convergence ⇒ far fewer
cells at matched accuracy. The infrastructure exists —
`TimeHarmonicSolver(problem, degree=2)` builds `("N1curl", 2)` and the
DG output spaces follow `self.degree` — it has simply never been gated.
**Scope guard:** time-harmonic E-formulation only. The magnetostatic
A-formulation's degree-2 failure (penalty-gauge null-space contamination,
920% field error with a clean solver exit, `core/solvers.py`) is a
formulation property, stays barred, and is *not* evidence about this
chunk. Curved second-order **geometry** (gmsh mesh order 2 — the answer
to `GEO-15`'s 3.3% faceting residual) is a separate knob, out of scope
here; a `GEO` chunk may cite this entry.
> * **Step 1 (gate) — the sphere at degree 2** ✅ *(2026-08-18, `tests/validation/test_lossy_sphere_degree2.py`, `20260818T110442Z_TH-12-step1-sphere-degree2-rss.log`, 7 s at `-n 2`; identical accuracy digits in `20260818T110346Z`)*. Degree 2 on the coarse 5 866-cell rung reads **0.1405%** interior relL2 against the gate ≤ 3.643% (degree-1 fine-rung record at 17 670 cells) — 25.9× the accuracy at 3.01× fewer cells; ohmic-power error **8.3869% → 0.0058%**. Negative control: degree 1 reproduces its 8.387% power error to 0.0001 pp. Cost: 7 591 → **39 634 DOFs** (5.22×), solve 0.93 → **4.03 s** (4.32×), summed peak RSS 388 → **1 036 MiB** (2.67×). `|Im P|/Re P` = **0.000e+00** at both orders. Caution: cgroup `memory.peak` is a lifetime high-water mark — price memory by summed `ru_maxrss` or a freshly recreated container.
> * **Step 2 (reading) — the coil at degree 2** ✅ *(2026-08-18 attempt 2, `20260818T200059Z_TH-12-step2-full.log`, 546 s at `-n 8`; probe `20260818T183449Z_TH-12-step2-probe.log`, calibration `20260818T183730Z_TH-12-step2-calibrate.log` fitting RSS exponent p = 1.271 on degree-1 rungs; audited COMPLIANT 2026-08-18 18:00 review)*. Degree 2 on the 138 619-cell fixture at 10 MHz reads ΔR deviation **−0.8508%** (ΔR +3.1985142e-01 Ω, ΔX −5.6252149e-01 Ω, ΔX ratio 0.9134) vs the step-4 h → 0 bracket [−2.1492%, −0.9050%] — outside by 0.054 pp past the upper edge (5× the 0.01 pp floor, inside the bracket's 1.24 pp width); degree 1 on the same mesh +1.5834% (reproduced to −0.00002 pp), so order moved the deviation −2.434 pp. Cost: 162 710 → **882 296 DOFs** (5.423×), solve 12.4 + 12.2 s → **235.4 + 266.4 s** (≈20×, superlinear), summed peak RSS 6.66 → **61.94 GiB** (96.8% of `memory.max`; 29% above the 48.04 GiB projection — treat p = 1.271 as a floor for degree-2 pricing). Controls: cells 138 619, σ = 0 dissipation +0.0 at both orders, drive mismatch 9.2e-35 / 1.0e-34.
>   **Defect, not loosened:** the complex-power identity fails at degree 2 (4.5931e-09 loaded / 3.0030e-09 free vs 1e-9; degree 1 8.07e-15 / 8.71e-15) — `W_e` 2.03e-13 → **7.16e-06 J**, `Im Z` +9.02 → **−2 117 Ω**: the ungauged gradient null space at second order swamps the magnetic term; common-mode, cancels in ΔZ, so the reading stands. Module fails by default; known-issues entry open with three dispositions. **Adjudicated (18:00 review):** no affordable (order, h) route to the 64 MHz bracket exists on this box (recorded in §2.2); no rung swap scoped; dispositions (a)/(c) contingent on step 3.
> * **Step 3 (mechanism) — is the degree-2 `W_e` explosion generic to incompatible drives, or coil-feed-specific?** ✅ *(2026-08-19, `tests/validation/test_degree2_energy_mechanism.py`, `20260819T183425Z_TH-12-step3-warm.log`, 8 passed / 10 s at `-n 2`; cold-compile `20260819T183329Z_TH-12-step3-compile.log` identical ratios)*. Reading **`COIL-SPECIFIC`** (pre-registered ≤ 10×-on-both band): cross-order `W_e/W_m` moves **1.155×** on the incompatible-drive smoke fixture (2.164348 → 2.499688), **1.015×** on the sphere (1.068190 → 1.052552), vs the coil's **3.426e+07×** (6.677632e-06 → 2.287540e+02) — `J·n ≠ 0` is not sufficient; hypothesis refuted. Anchors: smoke degree-1 reproduces `POST-5` **1.199162e-06 W** at rtol 1e-6 on 1 405 cells; sphere pair reproduces step 1 (0.1405% / 0.0058%) on 5 866 cells, 7 591 / 39 634 DOFs; `|Im P|/Re P` < 1e-9 both orders. **Confound carried:** baseline `W_e/W_m` differs 2.16 / 1.07 / 6.7e-6, so a fixed absolute contamination moves the coil ≈1e6× more — the step does not separate "the feed model injects it" from "only a `W_m ≫ W_e` fixture displays it"; discriminating needs a magnetically-dominated compatible-drive fixture or the absolute gradient content of `E` (disposition (b) proper) — scoped as `TH-13` (2026-08-23 weekly review). Known-issues entry stays open, the two degree-2 coil identity tests stay failing, no coil number moved.
>   Full narrative (steps 1–3): `docs/planning/plan-archive.md`, entry «§7 TH-12 steps 1–3 closure narrative — archived 2026-08-23 (weekly review)».
> * **Decision clause:** results go to the weekly review, which sets
>   the production element order for the §10 Phase-5/6 breakdown (the
>   32-port directive's cost rung is then priced at the chosen order).
>   No recorded degree-1 number moves; every degree-2 number lands as a
>   new row beside its degree-1 sibling.

**`GEO-14` — the shared ≈3% geometry floor: discriminate faceting from
resolution** ✅ *(commissioned 2026-08-13, closed 2026-08-15 on a refuted
hypothesis; full entry archived in `docs/planning/plan-archive.md`)*.
Step 1's one-command discriminator read the pre-registered **RESOLUTION**
band: the 64 MHz interior residual falls 3.643% → **1.781%** at the priced
55 251-cell mesh (rate 1.77 in h — no floor), with the 128 MHz record
reproduced to 1.8e-05 as the exact negative control
(`20260813T213156Z_GEO-14-step1-discriminator.log`). There is no shared
faceting floor to grade against; surface-graded sizing (the never-scoped
step 2) stays available to any future chunk that finds an actual floor.
The re-aim at `MAG-13`'s wire was declined — its own rung ladder already
attributes that residual to resolution.

**`TH-1` — Real complex time-harmonic formulation** ✅ *(all five steps,
2026-07-30/31; full step journal archived in `docs/planning/plan-archive.md`)*
> Replaced the `E = −ωA` proxy with an actual frequency-domain solve:
>
> ```
> ∇×(μᵣ⁻¹∇×E) − k₀²ε_c E = −jωμ₀J,    ε_c = εᵣ − j·σ/(ωε₀)
> ```
>
> **Formulation notes (live constraints, not history):**
> - **The sign convention is part of the spec.** The equation assumes `e^{+jωt}`,
>   matching `ε_c = εᵣ − j·σ/(ωε₀)`. Every analytic gate must be derived in the
>   same convention or validation fails spuriously with conjugated fields.
>   `ufl.inner` conjugates its second argument in complex mode — `ufl.dot` for
>   the load silently flips the convention.
> - **Do not port the gauge penalty.** At ω > 0 the operator acts as `−k₀²ε_c`
>   on the gradient subspace — nonzero everywhere, dissipative wherever σ > 0.
>   The `MAG-10` disease is statics-only; a penalty here would *add* error.
> - **The silent-failure mode is near-resonance ill-conditioning.** With PEC
>   boundaries and lossless air the matrix is exactly singular at cavity
>   eigenfrequencies — and an MRI coil is deliberately operated near resonance.
>   MUMPS returns clean exit codes on near-singular systems. The `core/resonance.py`
>   energy-continuity guard (step 5) is the calibrated detector: threshold 50
>   fires at ≈4% fractional detuning, verified against the `TH-9` fixture
>   (pole-law 16.505× vs 16.0×, 3.16%).
> - `build_material_fields` returns σ and εᵣ only — per-tag `μᵣ` needs that
>   function extended (`POST-3` step 5's first code touch).
>
> Results: MMS gate `E_ex = (sin ky, sin kz, sin kx)` at relative L2
> 11.26% → 5.66% over a 2× refinement, **rate 0.9929** vs the O(h) expectation;
> assembled operator complex symmetric to 1e-10 and **not** Hermitian — the
> structural signature that the loss term survived assembly. `solve()` raises
> in real mode (`require_complex_mode`); `TimeHarmonicProblem.dirichlet_e_field`
> imposes an analytic total field on the exterior N1curl dofs, which is how
> every closed-form gate below drives its box. Logs
> `20260731T003553Z_TH-1-steps123-mms.log`, `20260731T021415Z_TH-1-step5.log`.

**`TH-9` — PEC cavity resonance gate** ✅ *(2026-07-30, `core/cavity.py` +
`tests/validation/test_cavity_resonances.py`, `20260730T154846Z_TH-9.log`)*
> First four modes of a 1.0 × 0.8 × 0.6 m cavity match the closed form to
> **0.0436%** (720 cells) / 0.0102% (2268); rate 3.85; the 8 gradient modes
> return as a machine-zero cluster (3.2e-15). **Traps that stand:** a
> 1.0 × 0.7 × 0.5 box is degenerate (two modes coincide); PEC rows need a large
> diagonal in `A` and **unit** diagonal in `B` or the GHEP orthogonalisation is
> invalid. This is the known-frequency fixture the `TH-1` resonance guard is
> verified against.

**`TH-6` — lossy plane wave vs closed form** ✅ *(2026-07-31,
`tests/validation/test_lossy_plane_wave.py`, `20260731T020427Z_TH-6-gate3.log`)*
> The exact source-free `E = ẑe^{−jkx}`, `k = k₀√ε_c` (`Im k < 0` branch),
> imposed as Dirichlet data (εᵣ = 78, σ = 0.7 S/m, 127.74 MHz); *interior*
> slopes measured: `α` to **0.019%**, `β` to **0.059%**, L2 rate 0.9998. Clears
> §10's < 5% MVP bar (the bar is on the field norm — 16³ landed at 5.41% and
> the fix was mesh, not tolerance). Fixed in passing: `post/evaluation.py`
> gathered into a `float64` buffer and had never been called under the complex
> build; it now follows the function's dtype.

**`TH-7` — Validation: waveguide cutoff** ✅ *(2026-07-31,
`tests/validation/test_waveguide_cutoff.py`, `20260731T123411Z_TH-7-gate-final.log`)*
> Evanescent TE₁₀ below cutoff — decay from the transverse geometry against the
> operator's *real* part, complementing `TH-6`'s `Im ε_c` decay. `γ` to
> **0.006%** at 2.4 GHz; L2 rate 1.0013; three-frequency sweep each within
> 0.066% with end-to-end ratio to 0.038% (a k₀-blind solver returns ratio 1);
> `|Im E_y|/|Re E_y|` exactly 0.0 (convention check); the `TH-1` energy guard
> asserted quiet in-band.

**`TH-8` — Validation: dielectric sphere in a uniform quasi-static field** ✅
*(2026-07-31, 15:00 run; `20260731T200457Z_TH-8-gate-final.log`)*
> `E_in = 3/(ε+2)·E₀` measured at **2.443%** at the finest of three meshes,
> fitted rate 1.9675 (superconvergence of the probe-averaged functional, not a
> better element); interior uniformity 0.080%. The load-bearing negative
> control: dropping the sphere from the `material_map` under the *same*
> Dirichlet data moves the interior to 2348% off — the gate cannot be passing
> by reading back its boundary data. New fixture
> `MeshGenerator.sphere_in_box_domain`; its sizing is a gmsh `Ball` field, not
> a `Distance` field (unsigned — would coarsen toward the centre, where the
> gate measures). The lossy-sphere extension it named became `MAT-4` step 1;
> the `k₀R → 0` low-frequency-breakdown regime is still unstressed.
> `TH-8` is a cheap closed-form gate in the same mould; it, `TH-7`, or `TH-6`
> would have caught the `E = −ωA` defect immediately.

> (`TH-4` ✅ 2026-09-28; the note below explains why it was 🧪 before closing.)
> `TH-4` is 🧪 rather than ⚠️ because PETSc residual/conditioning diagnostics are
> meaningful regardless of which weak form is assembled.

> `TH-5` demoted off the MVP path: a birdcage operates inside an RF shield, so a
> **PEC outer boundary is physically correct** for the Phase-5 deliverables.
> ABC/PML is needed only for unshielded free-space validation geometries.

**`TH-13` — the degree-2 gradient-subspace injector: feed model, or any
`W_m ≫ W_e` fixture?** ✅ *(commissioned 2026-08-23, closed 2026-08-31 on
step 2; step 3a 🟡 2026-08-31; full 577-line entry archived in
`docs/planning/plan-archive.md`, 2026-09-06 weekly review).* Result: (A)
holds at round-off — the injector is the degree-1-only `H¹₀` source
projection, so the degree-2 `W_e` explosion is a property of the drive
projection, not of the feed; audited COMPLIANT 03:00 review. Step 3a's
matched projection meets both anchors with 6–14 orders of margin (residue
8.1e-17 / 1.8e-16 vs ≤ 1e-8; gradient share of `W_e` 4.6e-23 / 3.1e-21 vs
≤ 1e-6), default path bit-identical on control (b). Live carry-forward: the
coil degree-2 identity reds stay open in known-issues as the reopening
condition of the `TH-12` production-order clause (§10). Logs:
`20260830T020301Z_TH-13-step1.log`, `20260830T200543Z_TH-13-step1prime-final.log`,
`20260831T021154Z_TH-13-step2.log` (16 in the archived entry).


**`TH-15` narrative** — moved byte for byte to `docs/planning/chunks/TH-15.md` (`OPS-47`).

**`TH-14` — surface-impedance (Leontovich) boundary on conductor surfaces** ✅
*(closed 2026-09-14 04:30 slot; result block written by the 2026-09-19 weekly,
full narrative in `docs/planning/plan-archive.md`, entry «§7 TH-14 full
narrative — archived 2026-09-19 (weekly review)».)* Three gates, digits in the
§7 table row: **step 1** lossy-wall `TH-9` cavity Q against Pozar's Q_c
(σ = 1e4 miss +0.010 % in a 5 % band, `20260914T004807Z_TH-14.log`);
**slab step** copper Dodd–Deeds floor under `MAT-6`'s loop (ΔR −0.299 % in a
2 % band, surface-loss identity 1.31e-8, `20260914T094224Z_TH-14-step2.log`);
**birdcage step** copper F-small hole passes the `PORT-9`/`PORT-11` gates at
10 / 64 / 128 MHz with the surface-loss identity ≤ 3.0e-13 and the σ-ladder
converging onto the PEC matrix (`20260914T021500Z_TH-14.log`). The Fresnel
plane-wave step is retired as written (no separating anchor at copper).
**Carry-forwards:** the solid-side bracket is asserted at 10 MHz only (no
stored σ = 800 4×4 at 64/128 MHz); no absolute S claim on the copper coil —
that is `ANS-6`.


**`TH-5` — absorbing boundary condition (HFSS *Radiation*)** 🟡 *(**feature
ladder B1**, operator directive 2026-09-04, §9 item 5; the first Tier-B item,
queued after the conductor lineage; step 1a DONE 2026-09-30, §9 item 94,
record in the §7 row.)* **Why.** A shielded birdcage's shield
is its physical boundary and the PEC box is right; everything else —
unshielded coils, leads near the phantom surface, any SAR comparison against
an open-boundary HFSS model — is not, and `PORT-1` already had to carry a
named PEC-box systematic (`D∞ = +0.0169` at p = 1.657, effective-range
extrapolation). Jin ch. 9.
> * **Step 1.** First-order ABC `n × (∇ × E) + jk₀ n × (n × E) = 0` on the
>   outer box as a third `TimeHarmonicBoundaryCondition` mode — a surface
>   term of the `PORT-9` hook's shape, Jin §9.2. Gate (closed form): an
>   electric dipole in free space, `|E|` on a sphere of sampling points vs
>   the exact field, box at 0.5 λ and 1 λ; pre-state the band from Jin's
>   first-order-ABC reflection-coefficient curve at the box's incidence
>   angles, not from a wish.
> * **Step 2 (the systematic collapses, or does not).** The two-torus
>   mutual with the ABC box at `PORT-1`'s three paddings: the box term
>   `D∞` must fall below its own PEC-box magnitude by ≥ 5× — a
>   pre-registered identity between two boundary models on one fixture.
>   *(**Superseded (rule (j), annotated 2026-10-08 by §9 item 320 per the
>   03:00 daily's ruling (5)):** the ≥ 5× bar is a banked negative (2a,
>   printed false at every rung); step 2 was re-posed as 2a (S) + 2b (M)
>   and closed by 2c under the 2026-10-07 weekly's ruling — (R) retired,
>   r(W) a printed record. See the §7 row and `docs/planning/chunks/TH-5.md`.)*
> * **Step 3 (optional, review's call).** PML only if step 1's band is
>   insufficient for the 64/128 MHz birdcage box (λ ≈ 2.3–4.7 m, so the
>   box is electrically small and a first-order ABC may already be enough).
> * **Done-when (§4).** Steps 1–2 executed and asserted; §2.1 gains the
>   boundary-model line. *(Rule (j) annotation 2026-10-08, §9 item 320: met
>   as re-posed — step 2 = 2a (P0)/(E)/(S) + 2b (M)/(c2), asserted in
>   `20261008T183105Z_TH-5-step2c-W1.log` / `20261008T183455Z_TH-5-step2c-W2.log`;
>   §2.1 line written. Done-when not moved.)* Step 3's "may already be
>   enough" remains contradicted by 1b-0's closed form (see below) and stays
>   open.
> * **Amended 2026-09-29 18:00 review (rule (j); §9 lane 2).** The step-1
>   phrase "a third `TimeHarmonicBoundaryCondition` mode" is superseded by
>   §9's standing rule (no new features in `TH-3`'s ⚠️ set): the ABC is a
>   surface term through the `PORT-9` hook (`core/boundary_terms.py`). The
>   dipole gate becomes step 1b. **🟡 when step 1a is green in one logged
>   window** *(extraction named by the 2026-09-30 03:00 review, rule (j),
>   after attempt 1's red — bands unchanged; journal
>   `docs/planning/reviews/2026-09-30-daily.md` §4)*: on a driven TEM line
>   the first-order ABC gives \|Γ\| ≤ 1e-3 at λ₀/80 and λ₀/160, strictly
>   decreasing over λ₀/20 … /160 with a degree-1 LSQ rate ≥ 1.5, Γ read by
>   the volume least-squares projection of E_z onto e^{∓jk₀x} over the
>   whole cell columns in x ∈ [0.15, 0.75] m (attempt 1's 41-point fit
>   printed beside it as the record, not asserted); on a TE10 guide at 60°
>   it gives Γ = −1/3 within 3e-3 by the 41-point fit (Jin §9.1.1 (9.4) /
>   §9.2 (9.48)), with the exact-angle condition's Γ = 0 inside the same
>   band. ✅ when step 1b (the dipole gate, its band from 1a's measured
>   R(θ)) and step 2 are green and §2.1 gains the boundary-model line.
>   **Attempt 1 (2026-09-29, §9 item 90), recorded:** by the 41-point fit
>   (N) read 4.916e-3 / 2.018e-3 / 1.531e-3 / 1.296e-3, rate 0.617
>   (`20260930T022458Z_TH-5-step1a-W1.log:95–102`); (O1), (O2) and the
>   numpy anchors green. Family count: attempt 2 of 4 is §9 item 94.
>   **Attempt 2 (2026-09-30 04:30 slot, §9 item 94), green ⇒ 🟡:** by
>   projection (N) reads 3.814e-3 / 9.500e-4 / 2.373e-4 / 5.932e-5, rate
>   2.0021; the 41-point fit on the same solves reproduces attempt 1
>   (`20260930T094044Z_TH-5-step1a2-W1.log:105–112`).
> * **Step 1b re-posed 2026-09-30 10:30 review (rule (j); for the weekly
>   to ratify; journal `docs/planning/reviews/2026-09-30-daily-1030.md`
>   §4).** "Band from 1a's measured R(θ) at the box's incidence angles"
>   cannot produce a number: R(θ) is a planar plane-wave result (its
>   direction-weighted rms over a cube is ≈ 0.113 whatever the box size)
>   and omits the curvature term 1/(2(k₀R)²), 5.1 % at 0.5 λ₀; a box has no
>   closed-form reflection, and a point-sampled \|E\| read is what floored
>   1a's attempt 1. The clause becomes two steps.
>   **1b-0 (§9 item 107, numpy only):** the Hertzian-dipole field and the
>   first-order condition's l = 1 reflection on a sphere,
>   Γ_sph(x) = −e^{−2jx}/(2x² − 1 + 2jx), x = k₀R, verified by
>   finite-difference residuals before any solve. **Green 2026-09-30
>   (§9 item 107; `docs/planning/chunks/TH-5.md`): 1b below is
>   effective.**
>   **1b (§9 item 114; effective only if 1b-0 lands green — if it is red
>   this clause is void and the anchor returns to the weekly):** on a
>   spherical truncation at k₀R = π (`ball_domain`, R = 0.5 λ₀), a
>   z-directed dipole field imposed on a core, Γ read by the two-function
>   volume projection onto {E_out, conj E_out} over whole cells in
>   0.22 λ₀ < r < 0.42 λ₀, degree 1 at λ₀/40: with Z_s = η₀f/g, \|Γ\| ≤ B;
>   with the first-order ABC, \|Γ − Γ_sph(π)\| ≤ B; the PEC wall misses
>   Γ_sph by > B. B is frozen from the first window as 3 × the
>   exact-impedance reading rounded up to one figure (`TH-18` step 2's
>   rule) and is **not frozen above 1.69e-2** (= \|Γ_sph(π)\|/3). The
>   0.5 λ₀ box and k₀R = 2π are printed beside it, never gated.
> * **Step 2 — not queued; for the weekly (2026-09-30 10:30 review).** At
>   10 MHz the two-torus box sits at k₀R ≈ 0.03, where \|Γ_sph\| → 1: the
>   first-order ABC is a reflector there, equal to the natural wall to
>   O(k₀R). On a sphere model the static reaction field is +2/R³ under PEC
>   and −1/R³ under a magnetic wall, so the box term would flip sign and
>   about halve — a predicted ceiling of ≈ 2× on a fall the step bars at
>   ≥ 5×. A prediction on a sphere model of a cubic box, unmeasured; it
>   rests on 1b-0's closed form. Rubric element 2 (ceiling before factor)
>   keeps the step out of the queue as written; the weekly re-poses it
>   (mixed average, second order, PML) or commissions the measurement as a
>   banked negative, and revisits step 3's "electrically small box … may
>   already be enough", which the same closed form contradicts.
>   *(Rule (j) annotation 2026-10-08, §9 item 320: the weekly re-posed it —
>   2a measured the banked negative (sign flip held; r(W) 0.4846 / 0.3567 /
>   0.2266, not a steady half), 2b measured the mixed average, 2c closed the
>   step; "not queued" is history. The step-3 question is untouched.)*

**`TH-16` — symmetry planes** ✅ *(**feature ladder B2**, operator directive
2026-09-04; **closed 2026-09-30 by step 3c, §9 item 93**, on the gate
re-anchored by the 2026-09-30 03:00 review.)* **Result.** Per-face PEC / PMC
on the eigen path (step 1: box spectra worst 0.0214 % / 0.0508 % / 0.0474 %,
flag-off reproduction 0.000e+00, `20260929T033718Z_TH-16-step1-W1.log`); the
driven path with port rescaling in `ports/symmetry.py` (step 2: impedance
halves across a PEC cut and doubles across a PMC cut; record in the §7 row);
the cutter and the quarter birdcage (step 3c): three quarter solves reproduce,
as an identity, the 4×4 of the quarter's own mirror-built full mesh — (A)
1.520e-12 against 1e-9, summed quarter memory 0.370 GiB against the 0.425 GiB
bar; record in the §7 row. **Carry-forwards:** the comparison with the
independently meshed `S_10MHZ_EPS0_RECORD` is a printed two-mesh record, not
a gate (`20260930T003844Z_TH-16-step3b-W1.log:945`); not claimed — absolute S,
a separately meshed birdcage, 64 / 128 MHz, degree 2; the F-human memory
lever (`GEO-25`) is undated. *Narrative: `docs/planning/plan-archive.md`.*

**`TH-17` — birdcage eigenmodes** ⬜ *(**feature ladder B3**, operator
directive 2026-09-04; serial on `TH-15` and `PORT-14`.)* The `TH-9`
eigensolver on `birdcage_port_domain` with the coil as a PEC hole and
`PORT-14`'s capacitor sheets in the ring gaps, phantom lossless first. Gate:
the N mode frequencies of the ladder network against the lumped-element
closed form (Phase 6's named first target; the ladder's `L`/`M` from the
`PORT-13` Z-matrix, its `C` from the sheets) to a pre-stated band the review
sizes from the closed form's own approximation. Lossy phantom (complex
eigenvalues, Q) is step 2 and is what `WF-5` measures.

**`TH-18` — layered impedance boundary** 🟡 *(**feature ladder B4**, operator
directive 2026-09-04; serial on `TH-14`.)* `Z_s` of a foil-on-substrate
stack from the transmission-line recursion; gate: normal-incidence
reflection from a thin copper foil on a dielectric slab vs the closed-form
stack. ~~Last in Tier B; do not queue before B1–B3 have landed.~~
*(Struck 2026-09-29 10:30 review: §9 rule (5) lane 2, operator directive
2026-09-28, releases Tier B now, in the order `TH-16` → `TH-18` → `TH-5`.)*
**Done-when (written 2026-09-29 10:30 review, rule (j)).** 🟡 when step 1
is green in one logged window, with `TH-14` step 1 re-run green in the same
window:
- `layered_surface_impedance` meets its limits: t → ∞ equals `TH-14`'s
  Z_c to ≤ 1e-12; t → 0 gives Z_b; a thin PEC-backed layer gives jωμ₀t; the
  two-layer cascade holds; the λ/4 dielectric inverter holds.
- On `TH-14`'s cavity: at t = 5δ the Q ratio matches the TL prediction to
  ≤ 1e-5; at δ/2 on PEC the Q ratio is 6.2427 within 0.5 % and the
  frequency-shift ratio is 0.9680 within 2 %.

~~✅ when step 2 runs this entry's gate in FEM: normal-incidence |Γ| and ∠Γ
of a copper foil (as the layered term) on a dielectric slab against
Γ = (Z_s − η₀)/(Z_s + η₀) from the closed-form stack, at two thicknesses
(t ≫ δ and t ≲ δ), with the band pre-registered from a Leontovich run of the
same fixture. Not claimed: birdcage use, any AED comparison.~~
*(Amended 2026-09-29 18:00 review, rule (j), **for the weekly to ratify**:
at 127.74 MHz a copper foil's |Γ| is within 3e-5 of 1 at every t ≥ δ/100 —
the arithmetic that retired `TH-14`'s Fresnel step — so the copper gate
cannot tell t ≫ δ from t ≲ δ. The rung moves to σ = 10 S/m; copper is
printed as the PEC-limit reading.)* ✅ when step 2 runs this entry's gate in
FEM: normal-incidence |Γ| and ∠Γ at the front of a volumetric dielectric
slab (εr = 4, 0.2 m), terminated by a free-space-backed conductor foil as
the layered term, against Γ = (Z_s − η₀)/(Z_s + η₀) from the closed-form
stack, at t = 5δ and t = δ/2, σ = 10 S/m, 127.74 MHz, with the band
pre-registered from a bulk-Leontovich run of the same fixture. Copper is
printed, not gated. Not claimed: copper, birdcage use, any AED comparison.

**Step 1 green 2026-09-29 (item 82), 🟡:** new `core/layered_impedance.py`;
`solve_impedance_wall_cavity_mode` gains `surface_impedance_fn` (unset ⇒
`TH-14`'s path). Numpy limits n1–n5 ≤ 1.8e-16 relative (n3 6.7e-9 ≤ 1e-7),
n6 min Re Z_s 4.5e-10 ≥ 0. Cavity: premise Q_leo 8.0176788e+02 (2.0e-8 rel);
(i) Q_t5/Q_leo 0.9998744405 (|Δ| 1.4e-8); (ii) Q_thin/Q_leo 6.2426070
(−0.0017 %), Q_thin 5005.12 (+0.0086 %); (iii) shift ratio 0.9674439
(−0.055 %); negative control 0.8398 > 0.5; `air` Q 399.23 printed. 32 passed,
67 s, `-n 2` complex, `TH-14` step 1 re-run green in the same window
(`20260929T200249Z_TH-18-step1-W1.log`).

### MAT — Materials & phantoms (Phase 3)

| ID | Title | Status | Tier |
|---|---|---|---|
| `MAT-1` | Gelled saline presets (low/mid/high σ). Done-when (written 2026-09-29 10:30 review, rule (j); queued §9 lane 1). On `MAT-4` step 1's lossy sphere (h = R/10, 64 … | ✅ Done 2026-09-29 (§9 item 83, attempt 2): all three presets match the lossy-sphere closed form at 64 MHz, and the phantom path assembles exactly the same material arrays and operator as the `material_map` path. The original solve-level 1e-12 premise was re-anchored upstream of the solve, so no solver-determinism claim is made. **Done-when:** On `MAT-4` step 1's lossy sphere (h = R/10, 64 MHz), each preset passed as `phantom_material` gives mean SAR and Im/Re E_z within 10 % of `uniform_sphere_sar_closed_form` *History: `docs/planning/chunks/MAT-1.md`.* | standard *(re-declared 2026-09-29 18:00 review from smoke: measured 118 s at `-n 2`)* |
| `MAT-2` | Materials demonstrably affect solved fields | ✅ | standard |
| `MAT-3` | Debye/Cole-Cole dispersion models | 🟡 Step 1 (Debye on `TH-6`'s plane wave, 2026-09-29) and step 1b (a two-pole Cole–Cole medium, α = 0.30 / 0.35, 2026-10-26, item 648) are green and consumed by `Study` (`WF-15` step 8, item 653); step 2 waits on the Gabriel 1996 tissue table landing in `docs/references/` (operator). **Done-when:** *Step 2:* one Gabriel 1996 tissue reproduces the paper's σ and εr at 64 and 128 MHz to a stated band, gated only after the parameter table and published values are committed to `docs/references/`. ✅ = steps 1 and 2 green. *History: `docs/planning/chunks/MAT-3.md`.* | standard (re-declared 2026-10-26 18:00 audit: step 1 100 s, 1b W1b 60 s) |
| `MAT-4` | SAR computation `σ|E|²/(2ρ)` | ✅ *(step 4 closed 2026-09-06: the C95.3 mass-averaging operator on the **coil-driven** F-small field — 10 g C4 identity 0.3303 / 0.0756 / 0.0574 / 0.3132% at the imported unmoved 5% band, whole-phantom identity 1.58e-14, `23 passed` / 130 s at `-n 4`, `20260907T003548Z_MAT-4-step4.log`. **Scope: 10 g at fixed `h` on one fixture at 10 MHz; the 1 g column is a printed record at 6.8383% worst, verdict (b); no absolute SAR, no C95.3 compliance claim, no Larmor, no convergence claim**. Step 5, 2026-09-08: on the 0.0025 phantom rung the 10 g pairs re-assert at 0.0309–0.0644% and the 1 g pairs read **0.0957–0.1305%**, printed. **Step 5b ✅ 2026-09-08: 1 g and 10 g are both C4-gated on the 0.0025 rung** at the same imported unmoved 5% band, and the rung's mesh (199 920 / 58 866) is a version-tagged record at the 1% `CELL_COUNT_BAND`; `26 passed` / Status 0 / 313 s at `-n 4`, `20260908T200629Z_MAT-4-step5b.log`)* **Externally anchored 2026-09-18 (`ANS-2` AGREE): the coil-driven phantom power and pointwise SAR this identity is built on agree with HFSS at 10 MHz on F-small (5 % integral, ±10 % pointwise, residual = feed); annotation only, no band moved, no AED value tracked.** | standard *(step 4 heavy by ceiling, measured 130 s; step 5 251 s at `-n 4`)* |
| `MAT-5` | Temperature-dependent conductivity. Done-when (written 2026-09-29 18:00 review, rule (j); step 1 queued §9 lane 3 | ✅ Done: step 1 (the Picard fixed point on the Robin ball) and step 2 (per-cell σ(T) inside the EM solve on `MAT-4`'s sphere) both passed 2026-09-30, audit PASS. The (D′) transfer gate's planted failures were measured the same day (§9 item 109), closing the audit caveat. **Step-3 record (2026-10-05, §9 item 239, no glyph):** σ(T) fed back on P1's coil-driven F-small phantom (10 MHz, degree 1, G = 1 W, α = 0.02 /K an input): ρ₄ +0.099100 against the re-solved-E power-iteration μ_pred +0.099081, \|ρ₄/μ_pred − 1\| 1.83e-4 (band 0.10; (K-conv) 2.6e-4, (M0) 9.9e-2 ≥ 5e-3); fixed point ‖δ₈‖/‖θ₈‖ 8.5e-8; energy identity at iterate 8 (E\*) 9.4e-15, (B′) 6.2e-15, (C′) 3.2e-15, (D′) 0; α = 0 = P1 ((R0a) 7.8e-13 / 3.2e-12, (R0b) P₁ 6.4e-14, mean θ(t\*) Δ 3.5e-7 K, (D0) 0, (R0c) 0.0 at the provisional `L0_BAND` 1e-9); frozen-\|E\| μ_E +0.1338 printed, 35.0 % from μ_pred (the δE share is visible). Records: α = 0 steady state θ̄ 5.407745 K, θ_s 4.825138 K, max θ 9.544040 K; coupled θ̄₈/θ̄₁ − 1 +11.15 %, P_diss,8/P₁ 1.112344. `20261005T094115Z_MAT-5-step3-r2-W1-alpha0.log:1874–1883`, `20261005T094412Z_MAT-5-step3-r2-W2-coupled.log:1874–1908`. **Done-when:** ✅ = steps 1 and 2. *History: `docs/planning/chunks/MAT-5.md`.* | standard *(step 2 heavy by ceiling: 178 s against 180)* |
| `MAT-6` | Dodd–Deeds coil-over-lossy-half-space impedance | ✅ Closed: the Dodd–Deeds coil-over-half-space impedance is gated in the eddy-current regime, with step 11 scoped by the 2026-09-02 weekly review to promote the slab-refined fixture to production. The coil-at-Larmor case remains an extrapolation. *History: `docs/planning/chunks/MAT-6.md`.* | heavy |
| `MAT-8` | Finite-wire correction to the Dodd–Deeds closed form | ✅ Closed: the finite-wire correction to the Dodd–Deeds closed form is implemented and gated by its filament and perfect-conductor limits, smoke tier with no solve. *History: `docs/planning/chunks/MAT-8.md`.* | smoke |
| `MAT-9` | Steady heat conduction with a volumetric source and a convective boundary | ✅ Done 2026-09-29 (§9 item 76): the steady and transient conduction steps are green and the held root anchor passes in its scale-relative form. The lane-3 opening was ratified by the 2026-09-30 weekly. **Done-when:** ✅ = steps 1 and 2 green. *History: `docs/planning/chunks/MAT-9.md`.* | standard |
| `MAT-10` | **SAR as a heat source: `MAT-4`'s FEM SAR drives `MAT-9`'s Robin-cooled ball on the sphere's own submesh — Phase 8 step 2, gel regime** (§9 rule (5) lane 3; opened by the 2026-09-29 18:00 daily review under lane 3's closed-form licence — **ratified by the 2026-09-30 weekly**, dashboard). **Done-when (rule (j)).** On `MAT-4` step 1's lossy sphere (R = 0.01 m, h = R/10 in / R/5 out, σ 0.57, εᵣ 78, 64 MHz, complex build, `-n 2`), with q_h the DG0 cell average of G·½σ\|E\|² transferred to a `create_submesh` of cell tag 1 and h = 100 W/m²/K (Bi = 5/3): (P) a uniform q_cf on the same submesh meets `robin_sphere_steady_temperature` to interior relative L2 ≤ 1 % (`MAT-9` (i), imported; also the first complex-build thermal solve); (A) θ_h[q_h] meets the closed form with q_cf = G·ρ·`uniform_sphere_sar_closed_form` to ≤ 11 % — a composed band, `MAT-4`'s 10 % mean-SAR band plus (P)'s 1 %, never tightened after reading; (B) ∫_sub q_h dV = G·`dissipated_power_w` on the parent over tag 1 to ≤ 1e-10; (C) the discrete balance ≤ 1e-9; the Re-E-only and missing-½ controls each miss (A)'s band; `MAT-9`'s modules re-run green in the real build if `thermal/steady.py` changes (rule (c)) ⇒ ✅. **Out of scope:** perfusion, `MAT-5`, coil-driven SAR, transient, mass-averaged SAR. | ✅ *(**Closed 2026-09-30 (21:00 CDT slot take-next, §9 item 88):** `MAT-4`'s FEM SAR (DG0 cell average of G·½σ\|E\|², 74 020 → 19 496-cell submesh) drives the Robin ball to (A) 3.25 % interior L2 against the closed form (band 11 %; SAR-offset leg 3.55 %, (P) leg 0.319 % — the legs partly cancel), with (B) transfer 2.8e-15 and (C) balance 1.2e-13; controls (a) Re-E-only 79.3 % / (b) missing-½ 106.5 %; 13 passed, 70 s, `-n 2` complex (`20260930T021802Z_MAT-10-step1-W1.log:157–167,233`). `steady.py` unchanged (item 85's spec reused), so W2 was skipped under rule (c); the lane-3 opening was ratified by the 2026-09-30 weekly.)* | standard |
| `MAT-11` | The Pennes perfusion sink on the Robin ball | ✅ Done 2026-09-30 (§9 item 97) with audit PASS: the perfused Robin ball meets its closed-form, rate and balance gates, both controls miss, and the rule (c) re-runs reproduce. A complex-build leg (§9 item 120) reproduced every digit; the row was ratified by the 2026-09-30 weekly. **Done-when:** (i) P1 interior relative L2 on θ ≤ 1 % at the 0.01 m rung (`MAT-9` (i), imported); (ii) P1 whole-ball LSQ rate ≥ 1.7 on the 0.02 / 0.01 / 0.005 m ladder; (iii) ∫q dV = ∮hθ dS + ∫wθ dV to ≤ 1e-9 on every rung; the perfusion-dropped and sign-flipped controls each miss (i)'s band; `MAT-9`'s, `MAT-5`'s and `MAT-10`'s modules re-run green at the default w = 0 (rule (c)) ⇒ ✅. *History: `docs/planning/chunks/MAT-11.md`.* | standard |
| `MAT-12` | Transient SAR heating of the gel ball: `MAT-10`'s FEM SAR through `MAT-9`'s θ-scheme | ✅ Done 2026-09-30 (§9 item 100, first try) with audit PASS: the complex-build transient driven by FEM SAR meets its closed-form, calorimetric-slope, balance and steady-limit gates, and both controls miss. The row was ratified by the 2026-09-30 weekly. **Done-when:** (E) per-step balance ≤ 1e-9 on every run; (∞) backward Euler 20 × 5 R²ρc/k equals `solve_steady_conduction` with q_h on the same submesh to ≤ 1e-8; the Re-E-only and quasi-static controls each miss (A_t)'s band; `MAT-9`'s modules re-run green in the real build (rule (c)) ⇒ ✅. *History: `docs/planning/chunks/MAT-12.md`.* | standard |
| `MAT-13` | The transient Pennes equation on the Robin ball: the perfusion sink in the θ-scheme | ✅ Done 2026-09-30 (§9 item 108), audited PASS at the 18:00 review: the transient Pennes θ-scheme meets its closed-form, CN rate, balance, steady-limit and w → 0 gates with all controls missing. The complex-build leg (§9 item 120) reproduced every digit, and the lane-3 row was ratified by the 2026-09-30 weekly. **Done-when:** the perfusion-dropped and sign-flipped controls each miss (P_w)'s band and the perfused-term-dropped balance misses 1e-3; `MAT-9`'s modules (real build) and `MAT-12`'s (complex build) re-run green on their recorded digits (rule (c)) ⇒ ✅. *History: `docs/planning/chunks/MAT-13.md`.* | standard |
| `MAT-14` | SAR-driven transient Pennes on the lossy sphere: `MAT-10`'s FEM SAR through `MAT-13`'s perfused θ-scheme | ✅ Done: the complex-build SAR-driven perfused transient was green on every gate with each control missing its band, audited PASS at the 2026-10-01 03:00 review. The row was ratified by the 2026-09-30 weekly. *History: `docs/planning/chunks/MAT-14.md`.* | standard |
| `MAT-15` | Coil-driven heating of the gel phantom on F-small | ✅ Audit PASS 2026-10-02 (step 2(b) against the finite Robin cylinder series, LSQ rate 2.1511). Steps 3–5 are records on tuned copper fixtures — step 5 (P17, 2026-10-10, provisional) heats H9′'s tuned copper F-human's own phantom, identities green, thermal mirror 1.16 % / 0.76 % (≤ 5 %). Step 6 (P22) = `xl` entry 38 at degree 2, PENDING its degree knob (P22a, §9 item 413). **Step 6a** (§9 item 413, P22a): `FEM_EM_MAT15_STEP5_DEGREE=2` (guarded with `FEM_EM_WF10_DEGREE`) on step 5's module; unknowns 2E + 2F = 2 570 024 exact on the hole mesh (build-only, `-n 2`, `20261013T110156Z_MAT-15-step6a-W1.log`); unset path re-gated (`20261013T110401Z_MAT-15-step6a-W2-regate.log`). **Step 6b** (§9 item 515, P22b): `FEM_EM_WF10_NODE=d2root` (with both degree knobs at 2) reaches step 5's Z_RC / Z_th, its (S_rep) path `h9p-node-115.5209-d2.json` and its prints, with Y_in and SAR10g(c\*_θ) printed; (Deg-2) 2 570 024 at `-n 2` with the knob (`20261018T094411Z_MAT-15-step6b-W1.log`); unset path re-gated (`20261018T094613Z_MAT-15-step6b-W2-regate.log`, `20261018T094815Z_MAT-15-step6b-W3-regate.log`). Dose at the degree-2 root: `MAT-20` step 4 (P40, §9 item 557, 2026-10-20, `20261020T110139Z_MAT-20-step4-W1.log`). **Done-when:** step 1 asserted green ⇒ 🟡; step 2 ⇒ ✅. *History: `docs/planning/chunks/MAT-15.md`.* | heavy (step 2(b) 252 s / 273 s at `-n 2`) |
| `MAT-17` | Duty-cycled RF heating of the Robin gel ball: the on/off source's history and periodic steady state — Phase 8 (8a) | ✅ Audit PASS 2026-10-03: `source_schedule` in `solve_transient_conduction` reproduces the duty-cycle closed form on `MAT-9`'s ball (interior 4.39e-4 / 1.07e-3, BE rate 0.996) with `MAT-9` … `MAT-14` re-gated. Step 3 (2026-10-06, §10 P12) is a record on `MAT-15` P1's coil-driven source: Duhamel superposition 3.4e-14, duty 1 = P1 at 2.1e-15, ⟨θ⟩(t\*)/(D·P1) = 0.986, the (d) control on the last on-step (ratified by the weekly). **Done-when:** *step 1* (numpy) `robin_sphere_duty_cycle_temperature` (from rest and periodic) matches its Fourier route (`MAT-11`'s sphere continued to w = iωρc) ≤ 1e-9·θ_ss(0) and the superposition of `robin_sphere_transient_temperature` ≤ 1e-10; the periodic cycle mean equals D·θ_ss ≤ 1e-12; duty 1 reproduces the step response ≤ 1e-13; Robin residuals ≤ 1e-10·qR/h; three controls miss ⇒ 🟡. *Step 2* (FEM): an additive per-step source schedule in `solve_transient_conduction` reproduces step 1 on `MAT-9`'s ball at `MAT-9` step 2's imported bands (interior ≤ 1 %, BE rate ∈ [0.9, 1.2], per-step balance ≤ 1e-9), with `MAT-9`…`MAT-14`'s transient gates re-run green ⇒ ✅. *History: `docs/planning/chunks/MAT-17.md`.* | standard |
| `MAT-18` | The gel ball in a container wall: conduction through a two-layer Robin sphere (Phase 8 (8a)) | ✅ Audit PASS 2026-10-03: steps 1, 2(a) and 2(b) are green — the two-layer Robin sphere closed form and the piecewise k and ρc FEM in the steady and transient solvers. The opening and the R0 re-pose were ratified by the 2026-10-03 weekly. **Step 3 record (2026-10-06, §9 item 270, §10 P11; `20261006T200649Z_MAT-18-step3-W1.log`, 266 s):** the wall as a series resistance (h_eff 6.5517 W/m²/K) on `MAT-15` P1's coil-driven phantom keeps (E′) ≤ 2.3e-14, (L0) ≤ 6.7e-16 and C4 ≤ 0.049 %, and raises ⟨θ⟩(t\*) by +4.63 % on all four drives — a record, no glyph. **Done-when:** 2(a) + 2(b) ⇒ ✅. Not claimed: contact resistance, a real container's properties, the cylinder (`MAT-15`), perfusion. *History: `docs/planning/chunks/MAT-18.md`.* | standard |
| `MAT-19` | The virtual temperature probe: a point time series with a stated probe volume | ✅ Done 2026-10-05 (§9 item 243): the owned-cell ball-average probe on `MAT-12`'s sphere reads (T) 1.7315e-03 against the 1 % band and is width-invariant to 6.661e-16. Step 2 (item 258, 2026-10-06) is a record on `MAT-15` P1's coil-driven gel phantom at `ANS-9`'s probe centres: (A) ≤ 2.19e-14 and (W) 3.714e-11 on a submesh partition at `-n 4`, ball/point at B −1.24 %; its identities green ⇒ a record in this row, no glyph. **Step 3 record (2026-10-08, §9 item 314, §10 P10′, attempt 2 of 4):** on `ANS-7`'s exact hole fixture with the phantom at 0.005 m (93 956 cells, 8 349 phantom cells), degree 2, Cu, retuned on C_in = 63.587076 pF (in-model zero f₀/64 MHz − 1 = −2.53e-6, printed), (AB′) 1.257e-5, (A) ≤ 5.3e-15, (W) 1.96e-11 at `-n 4` and 1.69e-11 at `-n 2`, (Σ) 1.215e-3 against `C4_COVARIANCE_BAND`; identities green ⇒ a record in this row, no glyph. **Done-when:** an additive probe post-processor in `thermal/` (ball-averaged temperature against time: an integral over owned cells with no point owner, width-invariant as `POST-8`. *Amended 2026-10-05 03:00 review, rule (j), item 243, because `MAT-12`'s thermal field lives on a `create_submesh` mesh, where the point-owner rule raises; the weekly ratifies*); on `MAT-12`'s SAR-heated sphere the series at three radii is within 1 % of the closed form at every recorded time, and the series is width-invariant to 1e-9 at `-n 1 / 2 / 4`; a displaced-probe control misses; elapsed recorded ⇒ ✅. *History: `docs/planning/chunks/MAT-19.md`.* | standard (step 1) / heavy (step 2, measured 248–290 s) |
| `MAT-20` | Thermal dose and the six-minute window: the regulatory post-processors against closed forms | ✅ 2026-10-11 (§9 item 370): steps 1–2 closed on the closed forms and on `MAT-17` step 3's FEM history; step 3 (§10 P20, a record) wrote dose and six-minute peak per accepted watt on `MAT-15` step 5's series. Step 4 (record, §10 P40, §9 item 557, 2026-10-20): on entry 42's degree-2 series (gitignored `theta-series-ccw-d2.json`), (Q) / (Win) green at `DOSE_BAR` / `AVG_BAR` on the 10 g ball-mean and volume-mean series; the d1 → d2 dose move per accepted watt printed beside θ's −1.27 % (volume mean: dose −3.28e-4, excess over the 20 °C baseline −1.313 %; 10 g ball: −9.2e-6 / −0.026 %) (`20261020T110139Z_MAT-20-step4-W1.log`, smoke 5 s, `-n 2`). R3: the 43 °C kink at ε = 1e-3, the small-slope limit and the constant segment at `DOSE_BAR` (§9 item 638; `20261025T200112Z_MAT-20-R3-W1.log`, smoke 14 s, `-n 1`). **Done-when:** *Step 2* — on `MAT-17`'s duty-cycled closed-form history the running mean over an integer number of periods equals D × the continuous-source value at 1e-12 and the window maximum sits at the pre-computable phase; on `MAT-17` step 3's FEM history — its fixture re-run at `-n 2` with every step's ⟨θ⟩ kept (ruling 319-3 option (a), ratified 2026-10-10; the log's end-time records reproduced, printed; amended 2026-10-11 under rule (j) because the log's printed digits cannot reach `D_BAND`) — the operator applied to the duty series reproduces Σ(−1)^m of it applied to the ten Heaviside series at `D_BAND`; a shifted-window control misses ⇒ ✅. *History: `docs/planning/chunks/MAT-20.md`.* | smoke (steps 1–2 closed form) / heavy (FEM clause 270 s, `-n 2`; step 3 W1 160 s, `-n 8`) |
| `MAT-21` | The rectangular gel box with Robin walls under a uniform source | ✅ on step 1 with step 0's numpy records; steps 2–50 (the walled, duty-cycled, pre-equilibration and capacitive-wall box, the wait sweep in closed form and FEM, the three closed-form SAR pictures at entry 59's readout points and their 3 T / 1.5 T lever (P83), P76's read of entry 59, P88 — the dose is not per-watt — and P89′, the lattice-ball superposition whose verdict line reads "lattice resolves the excess", and P90 — the same closed form on entry 62's degree-2 lattice: θ(c\*, t\*) moves +3.665e-03, so it predicts entry 69's 10 % clause does not fire at c\*) are records in the chunk file. Step 51 (P91′, §9 items 991 → 999, 2026-11-15, a 🟡 record, no band): the box volume mean ⟨θ⟩(t\*) per deposited watt in closed form through the s₀ time split — uniform 0.02521149 K/W, 18-ball 1.5 T degree-1 lattice 0.02713084 K/W, ratio 1.076130; (S-a) seam ≤ 3.3e-16, (Bq-N) ≤ 2.8e-16, (Bq-split) 4.07e-13, (N-h) at h = 1e-11 1.016e-13 / 1.27e-14, all against `B1_BAND` (`20261116T030255Z_MAT-21-P91-W1b.log`). Step 52 (P92, §9 item 1008, 2026-11-16, a 🟡 record, no band, provisional for the weekly): the same series on entry 58's 3 T degree-1 lattice ((K) entry 60's strings = entry 58's) — ⟨θ⟩_lat,3T(t\*) 0.02713068 K/W, 3 T / 1.5 T 0.9999941; (Bq-split) 4.803e-13, (N-h) 1.27e-14 against `B1_BAND`; control (N-lat) 5.948e-6 > `B1_BAND` (`20261116T104406Z_MAT-21-P92-W1.log`). Open: P71 (a `mat:` example), a 🟡 record. *History: `docs/planning/chunks/MAT-21.md`.* | step 0 smoke, no MPI; step 1 heavy, -n 2, complex |
| `IMP-2` | A PEC prolate spheroid in a lossy host under a uniform axial field | ✅ Both steps went green on 2026-10-01: step 0's numpy closed form, then step 1's FEM at aspect 2 and 5 inside the imported bands with the transverse control missing as predicted, audited PASS at the 18:00 review. Quasi-static, PEC, imposed trace only; no rod, no coil and no 1 g SAR claim. **Done-when:** step 0 green ⇒ 🟡; step 1 at both aspects ⇒ ✅. *History: `docs/planning/chunks/IMP-2.md`.* | heavy |
| `IMP-4` | *(history and §9 briefs: `docs/planning/chunks/IMP-4.md`)* **A PEC circular cylinder (infinite, 2D-invariant) in a lossy host at finite k — the lead-shaft canonical conductor (Phase 7 (7a))** (opened by the 2026-10-03 03:00 daily under rule (5) lane 3; ratified by the 2026-10-03 weekly). `IMP-1`'s host (σ 0.57, εᵣ 78, 64 MHz, `K_FIXTURE`), a = 1 mm; TM_z n = 0 and TE_z n = 1 drives. Jin 3e §11.2.1 eq (11.56) and Exercise 11.4 (after Harrington, not held in the repo); App. E (E.1), (E.19). *(The 10-02 18:00 review's "§11.3" was a mis-citation.)* | ✅ **Audit PASS 2026-10-03 18:00 review.** 2026-10-03 (§9 item 206; opening ratified by the 2026-10-03 weekly). Step 1 (degree 2, 25 568 cells, exact trace, one mesh — not a convergence claim): \|β − 1\| TM 0.036 % / TE 0.774 %, annulus Φ miss +1.01 % / −0.23 %, cavity dofs zero; natural-cavity controls miss (`20261003T200454Z_IMP-4-step1-W1.log`, `…200703Z_…-W2.log`). Step 0 ✅ (§9 item 202, `20261003T172120Z_IMP-4-step0-W1.log`). **Done-when:** *step 0* (numpy): Bessel residuals ≤ 1e-12 (scale before cancellation); Wronskian ≤ 1e-13; PEC traces ≤ 1e-14·E₀; 2πaH_φ(a) = 4E₀/(ωμ₀H₀⁽²⁾(ka)) ≤ 1e-12; TE static limit ≤ 1e-9; three controls miss ⇒ 🟡. *Step 1* (FEM, degree 2, `TH-15` hole route): additive `as_hole=True` on `cylindrical_domain`; exact trace on the outer wall; end planes PEC (TM) or natural/PMC (TE); per drive \|β − 1\| ≤ `TH15_BETA_BAND`, annulus ⟨\|E\|²⟩ on the same mask ≤ `MAT4_BAND`, cavity dofs zero, natural-cavity control misses ⇒ ✅. Not claimed: a finite rod or tip, a coil, heating, published data. | standard *(step 1 heavy)* |
| `IMP-5` | *(history and §9 briefs: `docs/planning/chunks/IMP-5.md`)* **A dielectric-coated PEC circular cylinder (infinite, 2D-invariant) in a lossy host at finite k — the insulated lead shaft (Phase 7 (7a))** (opened by the 2026-10-03 18:00 daily under rule (5) lane 3; ratified by the 2026-10-03 weekly). `IMP-4`'s host and a = 1 mm; lossless εᵣ 3 sleeve to b = 1.5 mm; TM_z n = 0 and TE_z n = 1; series citation in the chunk file. **Done-when:** *step 0* (numpy): the series' residual, interface, limit and quasi-static anchors at the bands listed in the chunk file; three controls miss ⇒ 🟡. *Step 1* (FEM, degree 2, `TH-15` hole route): additive `coat_radius=` on `cylindrical_domain(as_hole=True)` (sleeve cell tag 4, `None` bit-identical); sleeve via `material_map`; exact trace; ends as `IMP-4`; per drive \|β − 1\| ≤ `TH15_BETA_BAND`, sleeve and near-lead host annulus ⟨\|E\|²⟩ each ≤ `MAT4_BAND`, cavity dofs zero; TE coating-omitted and TM natural-cavity controls miss ⇒ ✅. Not claimed: k_z ≠ 0 lead propagation, lossy or thin insulation, a finite lead, a coil, heating. | ✅ 2026-10-03 22:30 slot (step 1, §9 item 211, `20261004T035027Z_IMP-5-step1-W1a.log` / `…035247Z_…-W1b.log`, 133 / 168 s, -n 2, complex; W0 `…034835Z_…-W0.log` re-gates the `coat_radius=None` branch at 25 568 cells): \|β − 1\| TM 0.031 % / TE 0.545 %; sleeve / host-annulus Φ misses +1.04 / +0.36 % (TM), −0.17 / −0.19 % (TE); controls c1 (\|β − 1\| 1.474) and c2 (0.427) miss (`…W2a` / `…W2b`). Audited PASS 2026-10-04 03:00 review. Step 0 ✅ (§9 item 210, `20261004T003556Z_IMP-5-step0-W1b.log`, 59 s): every closed-form anchor ≤ 4.3e-14 of its band's scale, TE t/t₀ = −2.158 + 0.166j (dipole reversed), three controls miss. One mesh, accuracy bands, not a convergence claim. TM n = 0 barely sees a non-magnetic sleeve (β moves ≈ 3e-6, predicted), so TE carries the coating test. | standard |
| `IMP-6` | A PEC oblate spheroid in a lossy host under a uniform quasi-static field | ✅ (audit PASS): steps 0 / 1 / 3 gate the PEC oblate spheroid's rim and pole enhancement at aspect 2 / 5; step 2's rim ladder is a frozen record (≈ +0.9 % at k → 0 on this comparand, not rim grain) and step 4 (I31) builds the kept aspect-5 disc in F-clinical's ASTM box with (D1)–(D10) green at `-n 4`. **Done-when:** *step 0* (numpy, §9 item 403): the depolarisation factors, both potentials and fields, and the rim / pole identities at the bands in the chunk file, with the prolate-formula control missing ⇒ 🟡. *Step 1* (FEM, degree 2, `IMP-2`'s hole route with an additive `oblate=True` on `spheroid_in_box_domain`, a = 6 mm, aspect 2 and 5, §9 item 404): per aspect \|β − 1\| ≤ `TH15_BETA_BAND`, rim-ball Φ ≤ `MAT4_BAND`, cavity dofs zero; the axial-drive control misses ⇒ ✅. *History: `docs/planning/chunks/IMP-6.md`.* | step 0 smoke, no MPI; step 1 heavy, -n 2, complex |
| `IMP-7` | SAR heating around a PEC-sphere implant in a Robin-cooled gel ball: the implant-heating canonical in closed form (Phase 7 subgoal (7e); opened by the … | ✅ on steps 0 and 1; the curved P2 rate gate, the hollow-sphere transient, its CN twin, the duty cycle and the dose are records in the chunk file. **The pole's dose gap (I48c … I54) is closed as a record by the 11-02 weekly, no further step** — four fifths survive h and Δt (I51), the Δt table has no rate (I54), the series' comparand moves ≤ 6.6e-08 (I53′), M_p falls with the cap, not the error (I52) — and I49's (Bd)(ii) stays a reading (three points, one fixture). *History: `docs/planning/chunks/IMP-7.md`.* | step 0 standard, no MPI, real; step 1 standard / heavy, `-n 2`, real |
| `IMP-8` | The heated rod in gel: the finite line source in closed form (Phase 7; opened 2026-11-02 10:30 daily review, rule (5) lane 3; the weekly ratifies) | ✅ on the closed form (ratified by the 11-02 weekly); the FEM half, the rod's own thermal properties and the F2182-shaped probe forms — steps 9–21 (I55 … I73, I75, I78, I80, I81, I86), every gate at its imported band — are 🟡 records in the chunk file. Open: nothing enumerated; I61 is sized as `EX-182` (`imp:23`). *History: `docs/planning/chunks/IMP-8.md`.* | standard, `-n 1`, no solve |
| `IMP-3` | A rod in the coil-driven phantom: with / without on the loaded F-small birdcage at 64 MHz (Phase 7 subgoal (7b)) | ✅ as a fixture (audit PASS): the rod's identities are green on the record mesh and the 1 g rung, the Palace twins agree at matched order, and the pointwise tip ratio η_rod(T1) at degree 2 is a **four-mesh record, 3.482328 / 3.556262 / 3.530461 / 3.539570 (±1.1 %), not a converged value** (`xl` entry 43 / I30 read by the 10-24 weekly, first branch; known-issues 2026-10-14 retired; caveats in `IMP-3.md`). **Step 11 (I39, §9 item 626, 2026-10-25): the first coil-driven implant-heating record.** Under the ccw drive with one 1 W normalisation and an adiabatic rod wall, `MAT-15`'s (E′) reads 3.1e-15 / 5.0e-15 and (S2) 5.2e-6 / 5.0e-6 on the without / with arms, and the without-rod (C4) is ≤ 2.07 %. The 1 g tip θ(t\*) ratio rod / without is 1.111 / 1.088 beside Q 4.03 / 3.64; the max θ ratio is 1.041 (`20261025T110300Z_IMP-3-I39-W1.log:1840–1851,3661–3690`). The with-rod (C4) control did not miss. **On a matched mesh pair (I41, §9 item 639, 2026-10-25)** — the kept-rod build with tag 403 mapped to gel, 66 693 thermal cells against the rod arm's 51 693, the 51 693 tag-3 cells asserted identical — the 1 g tip θ(t\*) ratio R / M is **1.1045 top / 1.1087 bottom** beside Q 4.09 / 4.11, the max θ ratio 1.014, P_phantom +2.27 %, every identity green and the rod-free (C4) ≤ 0.93 % (`20261025T201613Z_IMP-3-I41-W1b.log:3645,3683–3689,3776–3779`); one mesh, a record, not a converged ratio. No compliance claim. **Done-when:** the identities green with the rod in on `GEO-36`'s 167 371-cell record mesh ⇒ 🟡; the same identities green on the 1 g rung with the tip ratio printed there ⇒ ✅ as a fixture. The ratio's accuracy is `ANS-10`'s (commissioned by the weekly that reads this). No compliance claim; step 11's heating is a record (audit PASS 2026-10-25 10:30). *History: `docs/planning/chunks/IMP-3.md`.* | heavy |
| `IMP-1` | A PEC-sphere implant in a lossy host under an imposed quasi-static field | ✅ Closed 2026-10-01 on the bands-only closure ratified by the 2026-09-30 weekly and audited PASS; ρ is kept as a record with no convergence claim. Step 3 records follow: 3(a) done, 3(b) landed as a record 2026-10-02 (§9 item 165), read as "the floor is neither comparand model error nor box trace alone; discretisation not separated". 3(c) a/6 at degree 3 with the exact trace gives ρ³ −0.9410 % against −0.9522 % at degree 2, which is the pre-registered **geometric** reading: residual geometric h², C ≈ 34 % · (h/a)² (Δ_trace −0.8465 % as predicted; `20261002T142346Z_IMP-1-step3c-W2.log`). It is a record, with no glyph change. **Done-when:** ✅ on (i)–(iii) at three rungs, (M), (iv′) and the control when §9 item 123 lands them green on `main` *History: `docs/planning/chunks/IMP-1.md`.* | standard *(step 2 heavy)* |

> ~~`MAT-1` is `⚠️` not because the preset table is wrong but because nothing
> consumes it.~~ *(Retired 2026-09-29, §9 item 83: the presets are now
> consumed by the real solve and gated against the lossy-sphere closed form.)*

**`MAT-4` — SAR computation** ✅ *(closed 2026-09-06 on step 4, extended by steps 5 / 5b on 2026-09-08; the step-4 / 5 / 5b plans, execution journals and audits are archived verbatim in `docs/planning/plan-archive.md`, entry «§7 MAT-4 full narrative — archived 2026-09-13 (weekly review)»; steps 1–3's plans and control-ceiling arithmetic were archived there earlier.)*
> * **Steps 1–3** ✅ (2026-08-03 / 04 / 07): mean SAR vs the lossy-sphere closed form `σ|3E₀/(ε_c+2)|²/(2ρ)` **3.42% / 3.54%** at h = R/10 under a 10% bound; the averaging operator exact at 1 g and 10 g on an imposed field (R = 0.03 m, kernel mass 0.0120% / 0.0044%), quadrature degree 16 required. Standing traps: `ufl.real` around any UFL comparison off the origin; density via `build_density_field`.
> * **Step 4** ✅ (2026-09-06, `20260907T003548Z_MAT-4-step4.log`, 23 passed, **132 s** at `-n 4`; audited PASS 2026-09-07 03:00): the C95.3 mass-averaging operator on the **coil-driven** F-small field at 10 MHz on the 120 499 / 2 746-cell fine-phantom rung — the four cyclic **10 g** C4 pairs **0.3303 / 0.0756 / 0.0574 / 0.3132%** against the imported, unmoved 5% band; whole-phantom ball identity **1.576517e-14** vs 1e-10 (mass 4.096723e-14, fine-rung power record 5.587038273e-08 W to 5.398570e-11); mis-paired 180° control 86.0132 / 85.9249 / 85.9671 / 85.9582%; 1 g pairs printed 6.8383 / 2.9297 / 0.4116 / 4.7159% (verdict (b), the ≈10-cell ball). One `src/` change: `build_density_field` accepts `0.0` (regression `20260907T003812Z_MAT-4-step4-src-regression.log`, 6 passed).
> * **Step 5** ✅ (2026-09-08, `20260908T140457Z_MAT-4-step5.log` + rerun `20260908T141147Z_MAT-4-step5-rerun.log`, 18 passed, 251 / 245 s at `-n 4`; audited PASS 2026-09-08 10:30): the `GEO-27` rung (`phantom_resolution` 0.0025, **199 920 / 58 866** cells on both windows) — 10 g pairs 0.0309 / 0.0060 / 0.0394 / 0.0644%, ball identity 8.126833e-14; 1 g pairs **0.0957 / 0.1199 / 0.1305 / 0.1065%** printed under clause (A); 10 g peaks 6.178300937 / 6.176390316 / 6.176760066 / 6.174323671e-07 W/kg (−2.67 … −2.43% rung to rung), 1 g averages 5.515422939 / 5.510146420 / 5.516752827 / 5.509553429e-07 W/kg, phantom power +1.5005% — none a compliance figure.
> * **Step 5b** ✅ (2026-09-08, `20260908T200629Z_MAT-4-step5b.log`, 26 passed, **313 s** at `-n 4`; audited PASS 2026-09-09 03:00): the 1 g column **asserted** at the same unmoved 5% band (0.0957 / 0.1199 / 0.1305 / 0.1065%, factor 38 of headroom), ball identity 5.262457e-14, mis-paired 1 g control 87.0143 / 87.0546 / 87.0506 / 87.0592% printed (predicted ≈85%); `ONE_GRAM_RUNG_CELL_RECORD = 199_920` and `ONE_GRAM_RUNG_PHANTOM_CELL_RECORD = 58_866` version-tagged (1\*) at the imported `CELL_COUNT_BAND` 1%, never at equality.
> * **Scope of the ✅, verbatim:** 1 g and 10 g mass-averaged SAR on the coil-driven field — C4 covariance identities at fixed `h` on F-small at 10 MHz; **no absolute SAR, no C95.3 compliance or limit claim, no homogeneity, no Larmor, no convergence rate** (two rungs of one quantity), no quadrature drive. **Live carry-forwards:** the absolute / compliance claim waits on `ANS-2`'s AED half; the owed example landed as `EX-53` (`mri:3`); `QUADRATURE_DEGREE` 16 is imported from the step-3 module.

**`MAT-2` — conductivity demonstrably drives the solved field** ✅ *(2026-07-31,
`tests/validation/test_lossy_plane_wave.py::test_conductivity_measurably_changes_the_field`,
log `20260731T020427Z_TH-6-gate3.log`, 21 s at `-n 2`, complex build)*
> Stronger than the "differ by a stated threshold" originally planned: σ = 0.1
> and σ = 1.4 S/m are solved on the same 24³ box at 127.74 MHz and each interior
> decay constant is compared with *its own* closed form — 2.1193 vs 2.1243 Np/m
> (0.233%) and 21.8781 vs 21.9045 Np/m (0.121%) — and the **ratio** 10.3232 is
> compared with the closed-form 10.3116 (0.113%). A σ-independent solver (the
> retired proxy) returns ratio 1. `MAT-6` still owns the coil-loading claim.

> `MAT-6` is the quantitative teeth for `MAT-2`. Dodd & Deeds (1968) gives the
> closed-form impedance change of a circular coil above a layered conductive
> half-space — the project's headline physics, *"the phantom loads the coil"*, in
> closed form. Upgrades `MAT-2` from "fields differ by a threshold" to "the coil
> impedance change matches a published solution", and bridges `TH-1` to `PORT-1`.

**`MAT-6` — the step record, compressed** *(steps 1–10b; the 2026-08-16-form ten-bullet step record is archived verbatim in `docs/planning/plan-archive.md`, entry «§7 MAT-6 full narrative — archived 2026-09-13 (weekly review)», beside the full narratives and probes; the step-11 promotion is in `docs/planning/chunks/MAT-6.md`.)*
> * **Steps 1–3** ✅ (2026-07-31 / 08-04): Dodd–Deeds anchored on the perfect-conductor limit (two derivations to 0.0002%); gate at f = 10 MHz, σ = 100 S/m, W = 0.15, 138 619 cells — **ΔR +0.3276882 Ω vs +0.3225961 Ω, 1.58%** (< 5%), re-gated on the production projected drive at **1.5834%** (the projection a 5e-5 no-op on ΔR); ΔX ratio 0.8123 → 0.9200, sign + order only; null-tagging control 1.31e-08, σ-blind fails by 100%. Traps: `ufl.max_value` does not compile in the complex build; a killed run leaves a stale FFCx lock (`rm -rf ~/.cache/fenics`).
> * **Steps 4–5, 8–9** ✅ (08-05 → 08-12): the ΔX drive gap is finite-wire discretisation (collapses 215× under wire refinement; ΔR's wire term 1.5834% → 1.0562% at `resolution_wire` 0.001); the slab knob owns ΔR — `resolution_near` 0.005 → 0.0025 takes ΔR **1.5834% → 0.2829%** on 417 914 cells (promotion deferred until `ANS-1` was adjudicated; landed as step 11, 2026-09-06, 0.2747% on the 418 888-cell fixture); box truncation owns ΔX (0.9200 / 0.9849 / 0.9960 → r∞ = 1.0023 at recovered p = 3.045); ΔR is *not* box-converged at W = 0.25 (0.38 pp at W = 0.35, control refuted, band not widened); composed readings compose **signed ΔZ**, never percent errors.
> * **Steps 6, 7** 🚫 / ✅ (08-08, 08-11): the 697 401-cell combined fixture OOM'd at the then-16 G cap; the operator raised it to 64 G and the additivity reading landed — ΔX ratio **0.9835 vs 0.9843 predicted, −0.080 pp, ADDITIVE**; `-n 8` established on this family (2.08× over `-n 4`); headless sessions run harness commands foreground.
> * **Steps 10 / 10a / 10b** (08-12, 08-16): the 895 974-cell composed fixture overran ≥ 5.7× (one solve > ≈1 700 s at `-n 8`); fill-in exonerated (1.693× flops vs 1.28× cells), memory pressure (MUMPS in-core 69 894 MB vs the 65 536 MiB cap) and load balancing the surviving suspects; 10b commissioned 2026-08-16 (weekly review) as one `-n 12` solve with `ICNTL(14)` raised under `timeout -k 30 1200`, done-when ≤ 2× the 257 s prediction. The overrun proved plain `timeout` does not stop `mpiexec` (hence `-k 30`) and wedged the container (recovery `up -d --force-recreate`); the factor stays resident after `solve()` (`time_harmonic.py:453`).

### POST — Post-processing & field extraction

| ID | Title | Status | Tier |
|---|---|---|---|
| `POST-1` | Interface-aware field extraction reliability | ✅ 2026-10-05 (item 245, 7b′). Step 7a is done (2026-10-03, item 187); 7b measured red on the dev run (item 203 BLOCKED) and was ruled the estimator's mesh floor in fixed space, not a field C4 finding. **7b′ green** on `TH-16`'s mirror-built F-small mesh (un-narrowed whole-face sheets, 64 MHz, degree 1): (M0) / (M1) / (I3) / (X0a) F 296 = 2V − 4 / (X0b) / (X1) / (P) / (B) green; the D2 image pairs (G1′) ≤ 2.2e-13, (G2′) 40 readings ≤ 1.6e-12, (G3′) face-level ≤ 5.0e-11 against the 5 % band; (W0) bitwise at `-n 4`, N / D / R_n within 3.8e-12. Elapsed: W1 88 s (`-n 2`; fixture 23.1 s: quarter 8.1, mirror 1.0, sweep 11.2, extraction 2.7), W2 83 s (`-n 4`). R_n 3.31–3.95 per drive is a record (no band), and the dropped C4 contrast on this mesh reads up to 69.6 % (printed only). **Done-when:** 7a asserted — a partition-free interface set (one tag-1 and one tag-2 cell per point, F = 2V − 4), the N1curl tangential trace continuous to 1e-10, and the normal ε̃E-jump residual and both one-sided fields against `LossySphereSeries` decreasing over three rungs, with levels recorded ⇒ stays 🟡; 7b recorded on the loaded birdcage with its D2 image pairs gated at the band its brief imports, on `TH-16`'s mirror-built F-small mesh, elapsed recorded ⇒ ✅ *(the C4 statement was dropped by the 10-03 weekly's 7b′ re-pose, §10; amended 2026-10-05 03:00 review, rule (j), item 245)*. *History: `docs/planning/chunks/POST-1.md`.* | standard; 7a heavy (one ≈ 195 s re-gate) |
| `POST-2` | Energy/consistency diagnostics. Done-when (written 2026-09-29 18:00 review, rule (j); queued §9 lane 1). On `TH-6`'s plane wave (0.1 m, 127.74 MHz, σ 0.7, εᵣ … | ✅ Closed 2026-09-30 and audited PASS: the ratio and span diagnostics read the lossy plane-wave closed form on `TH-6` within 5 % with the σ-blind control missing, and the cross-unit balance metric was replaced the same day by a reference-speed form. No energy quantity is gated and detecting a wrong solve is not claimed. **Done-when:** On `TH-6`'s plane wave (0.1 m, 127.74 MHz, σ 0.7, εᵣ 78, 24³, degree 1), with B from `magnetic_flux_density_from_e` and both fields sampled over every owned cell by `compute_tagged_vector_magnitude_stats`: `e_to_b_mean_ratio` is within 5 % of ω/\|k\|, and `e_span_ratio` / `b_span_ratio` are within 5 % of 1 − e^{−αΔx_c}. *History: `docs/planning/chunks/POST-2.md`.* | standard |
| `POST-3` | Replace vacuous consistency metrics | ✅ 2026-09-30 (18:00 review's adjudication on the 0.11 closing windows, §9 item 116; record in the block below) | standard |
| `POST-4` | Centerline point evaluation is rank-count-dependent: attribute and fix the ownership tie-break in `evaluate_vector_field_parallel` | ✅ *(chunk closed 2026-08-12 — every step closed or dispositioned; note the title's premise was itself refuted, the tie-break was never the defect. Step 1 ✅ 2026-08-11 — ownership **refuted**, 0/120 multi-claims; locus is the Lagrange-P1 interpolation, 1.163e+04× separation. Step 2 🚫 skipped. Step 3 ✅ 2026-08-11 — the centerline samples the source fields: **23.5539% → 0.008613%**, a 2735× collapse; known-issues entry **retired**. Step 4 ✅ 2026-08-12 — the export-path P1 artifact is **bounded and attributed**: midpoint relative medians **51.17% / 52.47% / 20.18%** (`A`/`B`/`E`), vertex/midpoint separation **0.42–0.68×** so the step's vertex-localization hypothesis is **REFUTED**, and a DG1 target reproduces all three sources to round-off — 100% of it is the P1 continuity constraint. All four steps now closed or dispositioned)* *History: `docs/planning/chunks/POST-4.md` (§9 item 607, the export-site re-census, 2026-10-23).* | standard |
| `POST-5` | Real Poynting power balance: wrong-sign boundary flux on the time-harmonic smoke fixture + `poynting_power_balance` raises on scalar `sigma=0.0` (`OPS-17` … | ✅ Closed: `poynting_power_balance` no longer raises on scalar zero conductivity, and step 3 found the boundary leg sound against closed form, overturning step 2's verdict. *History: `docs/planning/chunks/POST-5.md`.* | standard |
| `POST-7` | **Is the ≈ 1 % degree-2 terminal accounting residual at human scale an energy leak? `PORT-16`'s exact discrete power identity at degree 2 on F-human** (opened 2026-09-26 by the weekly from `xxl` 09-26, `WF-7` step 0d: the terminal residual is 1.08–1.18e-2 on all 32 columns at degree 2 against 4.06–4.62e-3 at degree 1; §10 chain H7) | ✅ **Audit PASS 2026-10-06 18:00 review** (`…200124Z_POST-7-step3-W1.log:79–81`, `:225–226`; the step-2 identity `20261004T070008Z_POST-7-step2.log:10960`). **2026-10-06 (step 3, §9 item 269, `20261006T200124Z_POST-7-step3-W1.log`, 71 s): the degree-1 (4.061–4.624e-3) and degree-2 (1.081–1.179e-2) terminal residuals are registered as records in the `WF-7` module, string-equal to their log tokens; with step 2's degree-2 identity (1.207e-15, `xl` 2026-10-04) the Done-when is met under the 2026-10-05 weekly's cron-window ruling.** *History: `docs/planning/chunks/POST-7.md`.* **Done-when (§4):** the identity asserted at degree 2 by an executed window, elapsed and `memory.peak` recorded. Green ⇒ the residual is `PORT-16`'s κ-class terminal-form deficit and is registered as a record at both orders; red ⇒ known-issues entry and a degree-2 accounting diagnosis step — **no band moves either way**. Negative control (rule (e)): the identity evaluated with one sheet's dissipation omitted must miss by ≥ 10× the band. | heavy, then `xl` |
| `POST-8` | Pointwise samples of a degree-1 N1curl field depend on the MPI partition: `ans:2`'s Row 1 moves −19 % … +112 % in SAR between `-n 2` and `-n 4`, while every … | ✅ Done 2026-10-22 (§9 items 586, 587; audit PASS, 10:30 review). **3b-iii done** 2026-10-22 (§9 item 587): `ans:2` honours `FEM_EM_ANS_SKIP_PRIVATE`; tracked Row 1 re-issued as the owner-cell average (`pointwise_rule` `owner_cell`), within 1e-9 of 3a's `-n 2` record at `-n 2 / 1 / 4` (8.329e-11 / 1.838e-10 / 1.890e-10); CG1 row a labelled record; SPEC v2; consumer green; known-issues 2026-09-27 retired. History of the steps: step 1–3a done; step 3b (the owner-cell rule as the default for every point-sampling site) is landing site by site — 3b-i, 3b-ii t1 / t2 / e1 / e2a / e2d (`mag:2`, `mag:4`, `mag:6`, `ports:12`) and e3 (`th:4`, `mag:1`, `mri:1`) done, **3b-ii-e4 done** 2026-10-22 (§9 item 586): `ports:1` on the owner rule, four records re-based `(1*)` at `-n 2`, 1e-6 bands unmoved, green at `-n 1 / 2 / 4`; cross-width spread of the S-derived records printed (`-n 4` raw 5.5e-10 / ‖S‖₂ 5.1e-10, `-n 1` 2.5e-10 / 2.2e-10); e2 and 3b-ii closed; 3b-iii (`ans:2`'s Row 1 re-issued) done 2026-10-22. Known-issues 2026-10-19 (`th:4` at `-n 4`) retired 2026-10-21; the (S) form on ghosted meshes is the 2026-10-19 weekly's. **Done-when:** *(amended 2026-10-02 18:00 review, rule (j), as the 10:30 review directed; the weekly ratifies)* a sampling rule whose output is width-invariant to 1e-9 relative, asserted at `-n 1`, `-n 2` and `-n 4` on `ans:2`'s points; `ans:2`'s Row 1 quoted as the owner-cell average (3a gates both; the 2026-10-21 weekly picked the owner-cell one — Row 4 is already the 1 g average), width-invariant to 1e-9 relative at the same three widths; the rule adopted by every point-sampling site of the census as the weekly classifies it; and the tracked Row 1 re-issued as the averaged row ⇒ ✅. *History: `docs/planning/chunks/POST-8.md`.* | standard; steps 2–3b-iii heavy (≈ 150–390 s per window) |
| `POST-9` | Adjoint sensitivities of safety functionals | ⬜ Future work for Phase 9, not dated: the weekly places it, and it is not on §9 rule (5)'s lanes unless the operator releases it. Steps 0 to 4 are written but none has started. *History: `docs/planning/chunks/POST-9.md`.* | unmeasured (step 1 smoke / standard) |
| `POST-10` | The point-owner rule on a `create_submesh` mesh: key submesh cells by their parent's original cell index *(history and §9 briefs: … | ✅ Done 2026-10-05 (§9 item 254): `point_owner_cells` / `evaluate_vector_field_parallel` take `parent=(parent_mesh, parent_cells)`, owner integers equal the oracle at `-n 1 / 2 / 4`, (N) ≤ 2.220e-15 and every misuse raises on every rank. The known-issues entry is retired and the `POST-8` / `POST-1` modules re-ran green unedited. **Done-when:** on a 4³ cube cut at x < 0.75 (288 cells), at P0–P5 (interior, facet, vertex, boundary facet, boundary vertex, outside): (K) the submesh owner equals the geometric oracle's minimum over the kept cells, exactly; (I) values and validity are bitwise equal to the parent's `candidate_cells=keep` rule on a planted DG0; (N) N1curl `a + b×x` meets its closed form ≤ 1e-12; (W) all of these at `-n 1/2/4`. Negative controls: (C1) the planted field separates a wrong owner by ≥ 1 at facet and vertex points; (C2) a submesh call without `parent` raises on every rank; (C3) / (C4) misuse raises. *History: `docs/planning/chunks/POST-10.md`.* | standard |
| `POST-6` | Arbitrary multi-port drive superposition | ✅ 2026-09-13: the 32-port ccw quadrature drive on `PORT-13`'s 16-leg ring fixture, superposed through `superpose_drives` at 10 MHz, is C16-invariant — worst `|B₁⁺|` rotation spread 0.8102 %, mirror 0.6769 % (both ≤ the imported 5 %), exact power identity 3.961e-15 (`20260913T185043Z_POST-6-step3.log:11798–11800`), with step 1's four anchors and the 4-leg cw control by record green (`…185405Z_…:1913`). Re-scoped 2026-09-13 (weekly, `auditor` DEMOTE(scope) concurring) to close on step 3; the `WF-6`/`WF-7` re-pointing clause dropped; no homogeneity, absolute or Larmor claim. *History: `docs/planning/chunks/POST-6.md`.* | standard; step 3 heavy (191 s at `-n 8`) |

> *(Closed-step plans, execution journals and audits for `POST-1` and
> `POST-3` are archived verbatim in `docs/planning/plan-archive.md`.)*
>
> **`POST-1`** 🟡 — steps 1–6 all ✅; **what remains is the coil+phantom
> application, where the chunk earns its ✅.** Compressed record: three real
> defects found and fixed (the complex→float64 cast scoring phantom metrics
> on `Re(E)` — closed by `POST-3` step 4; the ghost-cell double-count in
> tagged-cell aggregation — step 1; the guardrail's rank-local fallback —
> step 2, decision now on the allreduced interior count). Steps 3/4/4b
> measured drop-set semantics on the `TH-8` sphere and a chordal-error-free
> planar fixture: the interface guardrail is **harmless for means**
> (0.01 pp) and **harmful for peaks** (dropping the interface layer costs
> 2.157× in peak error — the layer is 22% *more* accurate than the
> interior), and step-3's `Re E`-vs-`|E|` question was discharged as an
> exact equality on the lossless sphere. **Step 5 flipped the production
> default to `prefer_interior=False`** (all four sites in
> `post/phantom_fields.py`; parameter retained, `True` path pinned; no
> landed gate moved — `MAT-4`'s mean SAR is structurally insensitive).
> Step 6 gated CSV-export/stats parity bit-for-bit in both modes.
> Transferable lessons that outlive the steps: sample `e_complex`, never
> `Re` of a phasor, for magnitude anchors (the planar fixture scored
> 61.8232% on the substitution); thin tagged regions for guardrail tests
> must be hexahedra; the `_owned_cell_count` AttributeError escape hatch is
> pinned, not fixed.

> The old flagship metric `e_to_b_mean_ratio` is by construction a function
> of B = ∇×E/(−jω): it reads ω/|k| on a plane wave (`POST-2` ✅ gates that),
> but, B being derived from E, it cannot fail for a wrong solve;
> deprecated-as-a-gate in `post/consistency.py`. *(Wording corrected
> 2026-09-30 18:00 review; it read "a mesh length scale, not physics",
> which `POST-2`'s plane-wave gate contradicts. The same phrase survives in
> the comments at `post/power_balance.py:3–5` and `post/quicklook.py:99–101`
> — a comment-only rider for the next item that opens either file.)* **`POST-3` replaced it with identities that can
> fail for real reasons — all five of its own steps are ✅** *(full plans +
> journals in `docs/planning/plan-archive.md`)*: step 1 Poynting real-power
> balance (4.13% at 24³, rate 0.987, σ-blind control 95.2%); step 2 σ(x) as
> a DG0 field (4.49% two-slab, rate 0.9915); step 3 total-current
> divergence residual (rate 0.942 in h, CG2/CG1 vacuity separation 1.5e13;
> environment note: `pc_type hypre` SIGABRTs this image — use `gamg`);
> step 4 phasor-magnitude semantics (both `Re`-cast sites removed,
> identities exact); step 5 piecewise μᵣ through both legs (4.3284% at 32³,
> rate 0.9922, both vacuity controls fire at 3.69× / 5.10×).
> **Reciprocity is discharged by `PORT-1` step 2** (decision 2026-08-02:
> the reaction-route `‖Z−Zᵀ‖/‖Z‖` *is* the field-level reciprocity, at
> machine precision). **What remains:** the 🟡 → ✅ flip is a review's
> adjudication, nothing else is open.
> **Adjudicated ✅ 2026-09-30 18:00 review** (`auditor` PASS) on the 0.11
> closing windows, a pure re-run with no edit (§9 item 116):
> `20260930T200032Z_POST-3-closing-W1.log:88–150` (23 passed, 123 s) and
> `20260930T200244Z_POST-3-closing-W2.log:89–106` (15 passed, 105 s). The
> structured-box Poynting figures equal their 0.7.2 records to the digit
> (4.1307 % / 4.4932 % / 4.3284 %; divergence 9.316e-2 → 6.358e-2, rate
> 0.942 computed in the assert, not printed). Two records moved on the
> gmsh-meshed fixtures and are re-read here, no band touched: CG1 residual
> 6.85e-15 (0.7.2: 6.14e-15), so the CG2/CG1 separation is 1.36e13, not
> 1.5e13; the phasor fixture's Re-cast deficit 44.3870 % (0.7.2: 45.3999 %),
> spread 0.55754, 0.99 points above the (0.4340, 0.4740) band's lower edge.
> Reciprocity rests on `PORT-1` ✅ (0.7.2 log) and the 0.11 four-port gates
> (`PORT-9` / `PORT-11`), not on these windows. The phasor identities run on
> the `prefer_interior=True` sample set (`dc607021`), not the production
> default.


**`POST-6` narrative** — moved byte for byte to `docs/planning/chunks/POST-6.md` (`OPS-47`).

**`POST-5` — real Poynting power balance: wrong-sign flux + the scalar-σ
raise** ✅ *(closed 2026-08-19, step 4 at the 15:00 slot; steps 1–4 ✅ 2026-08-18/19, each audited COMPLIANT)*. **Premise refuted, helper fixed.** The "forbidden" flux sign was never a Maxwell violation: with an impressed `J` the identity is three-term, `−∮½Re(E×H̄)·n̂dS = ½∫σ|E|²dV + ½Re∫E·J̄dV`, and the omitted source term was the whole O(100%) imbalance. **Gated numbers:** step 1 h-ladder h ∈ {0.030, 0.020, 0.015} imbalance 116.7465% / 115.4059% / 114.4227%, fitted rate **0.0290** vs pre-registered ≥ 0.7, sign never corrects ⇒ SOURCE/ASSEMBLY; `∮x·n̂dS = 3|Ω|` ratio 1.000000000000 (1e-10 band) — `ds` is outward; σ-blind control exactly 0.000000e+00 W at all rungs (`== 0.0`). Step 2 closed azimuthal drive (div J = 0, J·n = 0): 105.9632%, still negative ⇒ source compatibility excluded. Step 3 on `TH-6` plane wave, legs vs closed forms (`2αβ = ωμ₀σ = 7.060162290693e+02` at rtol 1e-12, common value 1.241101e-04 W): flux leg 8.1205% (12³) → **4.1141%** (24³), rate 0.981, volume leg 0.0711% / 0.0174%, both inside the 10% `POST5_STEP3_LEG_BAND` — `H = ∇×E/(−jωμᵣμ₀)` and the facet assembly are correct; three-term residual axial **16.7465%**, azimuthal **5.9632%** vs 25% `SOURCE_TERM_RESIDUAL_MAX`. Step 4: `poynting_power_balance` gained `current_density` / `source_measure` and returns `source_power_w` plus always-present `two_term_power_scale_w` / `two_term_relative_imbalance`; `test_time_harmonic_smoke_solve_conserves_real_power` lost `xfail(strict=True)` and passes at 16.7465% vs 25% (two-term 116.7465% reproduced at rtol 1e-6); J = `fem.Constant` zero ⇒ source term exactly 0.0 W, 7 other keys bit-identical (imbalance 8.185716%); σ-blind separation re-derived — ceiling 5.97×, floor 3.0×, measured 83.2535% = 4.97×.
**Logs:** `20260818T215101Z_POST-5-step1-ladder2.log` (5 s), `20260818T215117Z_POST-5-step1-negcontrol.log` (8 passed, 129 s), `20260819T051150Z_POST-5-step2-closed-drive2.log` (4 s), `20260819T051210Z_POST-5-step2-smoke-full.log` (10 passed + 1 xfailed), `20260819T123438Z_POST-5-step3.log`, `20260819T124405Z_POST-5-step3-source.log` (5 passed), `20260819T201005Z_POST-5-step4-smoke-final.log` (12 passed, 8 s), `20260819T200934Z_POST-5-step4-smoke-diag.log`, `20260819T200651Z_POST-5-step4-negcontrol.log` (15 passed, 152 s).
**Carry-forwards:** (a) the smoke fixture's 16.7465% / 5.9632% residuals are curl-trace discretisation error at ≈9 cells/λ, quoted as gated, not explained; the two-term 116.7465% / 105.9632% readings stay computable by design. (b) Any `SpatialCoordinate`-in-a-facet-integral form on a gmsh mesh must pin `metadata={"quadrature_degree": …}` — an unpinned one killed two windows (`20260818T213256Z`, `20260818T214040Z`); 0-byte FFCx cache stubs are a live lock (`find /root/.cache/fenics -size 0` in any stalled-JIT preflight; see known-issues). (c) `test_time_harmonic_smoke.py` and `test_poynting_balance.py` no longer fit one 540 s window together. (d) Natural-BC driven fixtures satisfy `½∫σ|E|² + ½Re∫E·J̄ = 0` by construction, so the three-term residual there is exactly the flux over scale. Known-issues defect-3 entry corrected on the sign claim. Ramp example `EX-26` ✅.
Full narrative: `docs/planning/plan-archive.md`, entry «§7 POST-5 full narrative — archived 2026-08-23 (weekly review)».

**`POST-4`** ✅ *(closed 2026-08-12; full step plans + journals in
`docs/planning/plan-archive.md`)*. The chunk title's premise was refuted by
its own step 1: the `evaluate_vector_field_parallel` ownership tie-break
was never the defect (0/120 multi-claims) and that function was never
changed. What the chunk bought: **(step 1)** the 23% centerline rank
spread attributed to the Lagrange-P1 interpolation of non-conforming
fields (interpolant path 97.9755% vs source path 0.008426%, 1.163e+04×
separation; a 56× artifact present even at `-n 1`); **(step 3)** the mri
centerline printout now samples the source fields — spread 23.5539% →
**0.008613%**, a 2735× collapse; known-issues entry retired; **(step 4)**
the export-path P1 artifact bounded and attributed — midpoint relative
medians 51.17% / 52.47% / 20.18% (`A`/`B`/`E`), a DG1 target reproduces
all three sources to round-off, so 100% of it is the P1 continuity
constraint; **(step 5)** the DG1/VTX faithful-export route priced — exact
round-trip fidelity, 10.5× disk, no wall-clock cost; complex-build
`VTXWriter` emits two real arrays per field; VTX rows are dof coordinates
incl. ghosts. **Standing (2026-08-12 review call):** DG1/VTX is the
faithful-export direction, but adoption is **blocked on the operator's
one-click ParaView check** of a DG1 `.bp` (dashboard Waiting-on-you);
"P1 + caveat" is the standing answer and no example switches its export
until that check returns.

### PORT — Ports & S-parameters (Phase 4)

**All `⚠️` chunks below sit on the §2.2 placeholder.**

| ID | Title | Status | Tier |
|---|---|---|---|
| `PORT-0` | Quarantine the placeholder coupling model | ✅ | smoke |
| `PORT-1` | **Real port excitation from the solved field** | ✅ *(closed 2026-08-15 review: done-when met through the package — `‖S−Sᵀ‖/‖S‖ = 2.5494e-05` vs the 1e-3 gate on a field-derived S, two-torus fixture; step 4, `20260813T183606Z`)* The padding test's known-issues entry (2026-09-09) is banked as an OBSERVATION (R1, attempt 3 of 4, 10-31 weekly); it retires with the padding module or a new consumer. *History: `docs/planning/chunks/PORT-1.md`.* | standard |
| `PORT-2` | Port data model and tagging contract | ✅ **Closed 2026-10-05 (step 1, §9 item 247):** on the exact ×1 strip the `_four_port_rung` `PortDefinition`s round-trip at `-n 2` and `-n 1` — reduced facet set {211…214} ∋ positives, cell set ∋ negative 1, n_i = 40 = global union on every port, areas +8.9e-16 … +1.6e-15 of the 0.007 × 0.008 rectangle, quarters [0, 1, 2, 3], planted missing tag raises and planted P2↔P4 swap passes existence + area but fails the round trip naming both (`20261005T140943Z_PORT-2-step1-W1.log:1886–1910`, `20261005T141117Z_PORT-2-step1-W2.log:1856–1873`). At `-n 2` neither rank's local tag set holds all four sheet tags ((c3) [False, False], printed), so the rank-local read `OPS-14` named would raise on both. **Done-when (2026-10-05 03:00 review, rule (j); the 10-03 weekly's lane-1 target):** on the exact-strip ×1 birdcage, the `PortDefinition`s built as `_four_port_rung` builds them round-trip at `-n 2` and `-n 1`. That means: `validate_required_port_tags_exist` passes on the reduced tag set; every port's facet set is non-empty with an exact owned-count reduction; every area equals the CAD rectangle within `STRIP_BAND`; port IDs ↔ tags are bijective, P_i on leg i; a planted missing tag raises; a planted swap fails the round trip naming both ports ⇒ ✅. *History: `docs/planning/chunks/PORT-2.md`.* | smoke (unit module); standard (round trip, measured 87 s / 91 s) |
| `PORT-3` | *(history and §9 briefs: `docs/planning/chunks/PORT-3.md`)* Calibration checklist → executable checks. **Done-when (written 2026-10-01 10:30 review, rule (j); the 2026-09-30 weekly's lane-1 target; §9 item 154):** new `ports/network_gates.py`: `four_port_network_gates` with `PORT-9` / `PORT-11` bands as defaults (each equal to its test constant); on `LEG_D_S_MATRIX_10MHZ` and `PORT11_S_64MHZ_RECORD` its readings equal the in-test arithmetic bitwise, reproduce the logged spreads / σ_max at print precision (\|reading − record\| ≤ ½ unit in the last printed digit; the separation ratio is printed against its record and its `resolved` verdict asserted — 2026-10-02 10:30 review ruling on §9 item 154), all-pass, and its derived terminated Z reproduces the logged Z to 1e-6; four planted matrices each fail exactly their own gate; `tests/ports/test_port_definition.py` green in the same window; harness, real, `-n 1` ⇒ 🧪 → ✅. **Closed 2026-10-02 (§9 item 178, attempt 2 of 4):** W2 `20261002T201726Z_PORT-3-W2.log`, real, `-n 1`, Status 0, 26 passed / 4 skipped, elapsed 51 s (standard tier) — (A1) bitwise at 10 / 64 MHz; (A2) spreads / pooled / σ_max / column power within ½ ulp (tightest 4.147e-5 ≤ 5e-5); (A3) derived Z ≤ 1e-6; (N1) unitary phase plant reciprocity 1.878e-2 with Δσ_max 0.0, (N2)–(N4) each fail only their own gate; `test_port_definition.py` green in the same window. | ✅ **Audit PASS 2026-10-02 18:00 review.** | standard |
| `PORT-4` | Multi-port drive/termination consistency. **Done-when (written 2026-09-29 10:30 review, rule (j); queued §9 lane 1).** On `PORT-9` leg (d)'s 116 085-cell fixture at 10 MHz with P2 / P4 terminated in 200 / 25 Ω, Y = I_i(drive j)/V_src is symmetric to ≤ 1e-12 relative; the 50 Ω sweep reproduces `S_10MHZ_EPS0_RECORD` at `RECORD_RTOL`; the termination-blind control's worst pair is ≥ 0.5; the placeholder `test_single_port_excitation_returns_finite_estimates` is replaced; harness `-n 2` complex, elapsed recorded ⇒ ✅ (known-issues entry 6 retired). | ✅ *(2026-09-29, 12:00 slot, §9 item 77: premise 116 085 cells, 50 Ω S vs record 5.061e-15; unequal-termination ‖Y − Yᵀ‖_F/‖Y‖_F 8.513e-15 ≤ 1e-12 (50 Ω 5.967e-15 printed), every Y_ij/Y_ji = 1 to 12 digits; termination-blind control worst pair (P2, P4) 0.875000 ≥ 0.5 (pred. 0.875); records: sweep s_matrix asymmetry 1.604e-01, ‖Y − (s_to_z(S₅₀) + diag Z_p)⁻¹‖/‖Y‖ 4.154e-03 (5.7× `PORT-14`'s 200 Ω floor, printed); 16 passed, 91 s, `-n 2` complex, `20260929T170418Z_PORT-4-W1.log:1866,1900–1909,2002`. Tests only, known-issues entry 6 retired.)* | standard |
| `PORT-5` | S-matrix reciprocity/passivity metrics. Done-when (written 2026-09-30 03:00 review, rule (j); queued §9 lane 1). (a) Through `run_n_port_sparameter_sweep`, … | ✅ Closed 2026-09-30 and audited PASS at the 10:30 review: through the real assembler, synthesised port states reproduce known S matrices within 1e-13 and the sanity metrics return their closed forms, with the transposed-comparand control missing by exactly 1. The frequency-sweep ambition has no landed code and stays outside the closure. **Done-when:** (a) Through `run_n_port_sparameter_sweep`, port states synthesised from a known S₀ (a lossless reciprocal line, an ideal 3-port circulator, a seeded non-symmetric 4-port) assemble to S₀ within 1e-13, and the circulator misses its own transpose by exactly 1. *History: `docs/planning/chunks/PORT-5.md`.* | smoke *(measured 2 s; was standard)* |
| `PORT-6` | Frequency sweep orchestration. Done-when (written 2026-09-29 18:00 review, rule (j); queued §9 lane 1). `plan_frequency_sweep` equals the integer-lattice … | ✅ Closed 2026-09-30: the planner equals the integer-lattice closed form bitwise on the seeded 500-case family and the frozen `c1d167d` planner reproduces the `OPS-56` failures; the validity floor was fixed the same day, raising for steps under 16 ulp of the carrier. Planning only: orchestration of solves over a plan is outside the closure. **Done-when:** `plan_frequency_sweep` equals the integer-lattice closed form bitwise on a seeded 500-case family (integer-Hz carriers 1e6–3e8 Hz, steps 100 Hz–1 MHz, on- and off-grid endpoints, 0–3 refinement windows); `z0_ohm` / `port_order` pass through exactly; the frozen `c1d167d` planner reproduces the three `OPS-56` failures (`20260921T112323Z_OPS-56.log:34–36`); the existing 19 tests are green unchanged. *History: `docs/planning/chunks/PORT-6.md`.* | smoke |
| `PORT-7` | Touchstone metadata + parser cross-check — **Done-when (written 2026-09-28 18:00 review, lane 1, §9 item 64):** scikit-rf 2.1.0 (pinned, `OPS-60`) as the independent reader: its premise on spec-layout 2- and 3-port literals (1e-15); our writer → `skrf.Network` and `Network.write_touchstone` → `load_touchstone` identical to ≤ 1e-12 at N = 2–4 on the stored `PORT-15` records and on seeded non-reciprocal matrices; `nports` / `f` / `z0` metadata; the transpose control on the non-reciprocal fixture; through the harness, elapsed recorded ⇒ ✅. *(Drafting note: `touchstone.py:153–159` writes every N column-major on one line, while Touchstone 1.x lists N ≥ 3 row by row — predicted red at N ≥ 3; only a non-reciprocal fixture can see it.)* **✅ 2026-09-29 (22:30 CDT slot, §9 item 64):** premise green (skrf reads the spec-layout 2- and 3-port literals to 0.0); W1 matched the prediction exactly — N = 2 green (6.6e-14), N ≥ 3 red with skrf reading Sᵀ (max|S_read − Sᵀ| 4.7e-14 / 6.0e-14 at N = 3 / 4, max|ΔS| 0.439 / 0.566) and `load_touchstone` rejecting skrf's one-row-per-line N ≥ 3 files (ragged rows); the stored reciprocal records (a)/(b) passed writer → skrf *through* the bug (max|S − Sᵀ| 7.9e-15 / 1.4e-14). Branch fired: writer N ≥ 3 row-major, one matrix row per line, ≤ 4 pairs per line, N ≤ 2 untouched; loader collects tokens across lines, N ≥ 3 in C order. W2 19 passed at `-n 2`, 0 skipped: writer → skrf ≤ 6.3e-14, skrf → `load_touchstone` 0.0 on all five cases, `nports` / `f` (1e-3 Hz) / `z0` = 50 exact; asserted transpose control misses 0.337 / 0.439 / 0.566 (floor 1e-2, ceiling 1); existing `test_touchstone_export.py` green unedited (rule (c)). Smoke, 2 s (`20260929T033225Z_PORT-7-W2.log`; W1 `…033144Z_PORT-7-W1.log` red on a test-literal `repr` bug in the premise, W1b `…033200Z_PORT-7-W1b.log` the clean pre-fix reading). Known-issues entry opened and retired in the closing commit. | ✅ | smoke |
| `PORT-8` | Port-orientation sensitivity — **Done-when (written 2026-09-29 03:00 review, lane 1, §9 item 70):** on `PORT-9` leg (d)'s 116 085-cell fixture at 10 MHz, reversing P2 / P4's `LumpedSheetPortSpec.drive_direction` gives S′ = DSD and Z′ = DZD (D = diag(1, −1, 1, −1)) to ≤ 1e-12 relative against the unflipped sweep on the same mesh; the unflipped S reproduces `S_10MHZ_EPS0_RECORD` at `RECORD_RTOL`; the orientation-blind control ‖DSD − S‖/‖S‖ ≥ 0.5; the `PORT-0` placeholder tests in `tests/ports/test_port_orientation_sensitivity.py` replaced; harness at `-n 2`, complex, elapsed recorded ⇒ ✅ | ✅ **2026-09-29 (06:00 slot, §9 item 70):** Done-when met. The placeholders are replaced by `test_reversing_p2_p4_drive_conjugates_s_and_z_by_d`: 116 085 cells, and the unflipped S matches the record to 4.129e-15. ‖S′ − DSD‖/‖S‖ = ‖Z′ − DZD‖/‖Z‖ = 0.000e+00 (≤ 1e-12). Adjacent entries are negated, self and opposite unchanged, and the control reads 1.4887 ≥ 0.5 (pred. 1.489). 13 passed, 85 s, `-n 2` (`20260929T110818Z_PORT-8-W1.log:1866,1903–1904`). | standard |
| `PORT-9` | Lumped-element port boundary condition (the birdcage port model) | ✅ Closed 2026-08-25 by leg (d1′), as the note at the head of the history records; the two-torus class re-record that preceded it moved no band. *History: `docs/planning/chunks/PORT-9.md`.* | standard |
| `PORT-10` | The two `PORT-1` systematics: composition measured, not assumed | ✅ 2026-08-16 (cross-term **−0.0604 pp** inside the pre-stated ±0.5 pp) | heavy |
| `PORT-11` | Lumped-sheet ports on the gapped birdcage at 64 MHz (then 128): `PORT-9`'s three gates in the displacement-current regime | ✅ Closed: step 2's 4×4 at 64 MHz ran at standard tier after step 1's probe cleared the phantom resolution stop rule, and all three gates were green. *History: `docs/planning/chunks/PORT-11.md`.* *(2026-10-01, `PORT-21` step 5: with `FEM_EM_PORT11_WIDTH_CORRECTION=in_run` the three gates stay green at 64 MHz κ-corrected and the displaced control still breaks (iii′); `20261001T171613Z_PORT-21-step5-W2.log`.)* | heavy (probe first; **step 2 measured and ran standard**) |
| `PORT-12` | The two-torus gap-route record drifts with rank width on an already-plumbed fixture: `tests/validation/test_port_lumped_two_torus.py` reads gap ratio 0.894141 … | ✅ Closed: step 1's gap-ratio record reproduces at both rank widths inside the unmoved band, and the 1e-8 assert was probed load-bearing. No band widened, no record rewritten, and no root-cause claim is made. *History: `docs/planning/chunks/PORT-12.md`.* | standard (complex, 84 s per width) |
| `PORT-13` | Phase-6 ring-rung solve probe | ✅ Closed: steps 1, 2 and 3 are all closed and the full 32×32 on the 32-ring-port high-pass layout is reciprocal / passive / C16 per step 3. A step 4 is a review's to scope. *History: `docs/planning/chunks/PORT-13.md`.* | heavy (probe first) |
| `PORT-14` | Lumped RLC sheets | ✅ Closed 2026-09-13 at step 3: the κ-derived width route is registered at 64 MHz with corrected residuals well under 1e-3 and the `PORT-9` gate module green. κ is carried as the sheet's named systematic, and the 128 MHz pair stays a printed out-of-sample reading. *History: `docs/planning/chunks/PORT-14.md`.* | standard (steps 1–2b); heavy by ceiling for steps 1c–3 (`timeout -k 30 590`, windows ≤ 196 s — cell reconciled 2026-09-14 03:00 review on the auditor's caveat) |
| `PORT-15` | The circuit layer | ✅ 2026-09-14 02:07Z on gate (i) + the tuned `S₁₁` (one F-small 4-leg fixture, 64 MHz; gate (ii), the 32-port mode spectrum, is `TH-17`'s and not claimed). Step 1's ladder-network closed form and termination reduction hold as pure-numpy identities. Gate (i): discharged by `PORT-14` step 3, stored record here (`tests/validation/test_port_circuit_layer_field.py`, step 2 ✅ 2026-09-14 01:07Z — stored 64 MHz residuals C 5.36e-5 / L 2.0e-6 under 1e-3, records reproduced ≤ 1.8e-13, `20260914T010448Z_PORT-15.log`). Step 3 ✅: `C_tuned` = 15.570 pF (P2..P4, P1 driven) zeroes `Im Z_in` to 2.1e-15 relative; in-model vs circuit tuned `S₁₁` 8.26e-5, 2×2 6.12e-5 under 1e-3 (`20260914T020419Z_PORT-15.log`, 7 passed, 176 s). "Tuned" = series resonance, `R_in` 6.77 Ω, `|S₁₁|` 0.761 — not matched. *History: `docs/planning/chunks/PORT-15.md`.* | smoke (step 1, measured 4 s); standard for the field-side steps (step 3 window 176 s) |
| `PORT-16` | The ≈1 %-of-supplied accounting gap on the 4-leg fixture | ✅ Closed 2026-09-07: the exact discrete power identity closes on every single drive and on the superposed quadrature drive, and the accounting gap is attributed to the terminal form's Cauchy–Schwarz deficit. Step 3's gap h-rate is optional, unclaimed and the weekly's to commission. *History: `docs/planning/chunks/PORT-16.md`.* | standard (three windows ≤ 180 s at `-n 2`; the ×0.6 rung measured 136 s) |
| `PORT-17` | Wave ports / coax feeds — HFSS *Wave Port*; low priority: MRI coils are fed at lumped points through matching networks — **feature ladder C2** (operator directive 2026-09-04; commission only if a benchmark demands it) *(**renumbered from `PORT-16` by the 2026-09-09 03:00 daily review**, on the 02:15 weekly's finding 1: §7 carried two `PORT-16` rows, which breaks the stable-ID contract. The closed, audited, log-bearing chunk above keeps the ID; this unopened ladder entry moves. `PORT-17` was unused repo-wide. The §9 item-5 ladder table moved with it; no other file referenced this row)* | ⬜ | standard |
| `PORT-18` | Does the lumped-sheet port read its current off one side of a discontinuous component | ✅ **Closed 2026-10-05 (step 2, §9 item 246):** on the exact ×1 strip after a P1-driven degree-1 10 MHz solve, the `'+'` and `'−'` sheet currents agree to ≤ 5.0e-16 relative and `J[E]` ≤ 3.3e-15 on all four sheets at `-n 2` and `-n 1`, with the planted DG0 jump reading J = 2 exactly (`20261005T140253Z_PORT-18-step2-W1.log:1919–1946`, `20261005T140501Z_PORT-18-step2-W2.log:1897–1924`). The read-back's `'+'` restriction is not a source of asymmetry on exact-rectangle sheets; no `S₁₁` / `Z_in` / `f₀` figure is licensed. 🧪 Measured 2026-09-10: the premise was false — the sheet normal is azimuthal — and the '+'-side mechanism is excluded for both the read-back and the source term. The printed lead is taken up by `GEO-32`; `ANS-4` step 2a and `WF-6` step 4f stay blocked. **Done-when (2026-10-05 03:00 review, rule (j); the 10-03 weekly's lane-1 target):** on the exact-rectangle ×1 fixture, κ off, after a P1-driven degree-1 10 MHz solve at `-n 2` and `-n 1`, `sheet_terminal_current`'s `'+'` reading and its `'−'` mirror agree within `EXACTNESS_BAND` 1e-9 on all four sheets, and a planted DG0 discontinuity reads J = 2 ⇒ 🧪 → ✅. *(The 09-10 record holds no ragged one-sided discrepancy to collapse: it read the CAD face, where R⁺ = R⁻ to print precision.)* *History: `docs/planning/chunks/PORT-18.md`.* | standard (heavy by ceiling), `-n 2` + `-n 1`, no solve — **measured 154 s + 68 s**; step 2: standard by prediction, heavy by ceiling, `-n 2` + `-n 1`, one degree-1 solve — **measured 121 s + 117 s** |
| `PORT-19` | Reuse the factorisation across the drives of one sweep | ✅ Closed: one factorisation per sweep instead of one per port, with `S`, `Z` and the kept fields matching per-drive solves inside their asserted bounds at both rank widths. Step 6 recorded the two remaining example callers green on the reuse-on default. *History: `docs/planning/chunks/PORT-19.md`.* | smoke (step 1), standard (step 2), heavy (step 3) |
| `PORT-20` | The gap-voltage (current-drive) route must report the S-matrix | ✅ All three steps executed 2026-09-20: the current-drive route now reports S as `z_to_s(Z)`, the closed-form identity and both negative controls are asserted, and the S-derived records are re-recorded at unmoved bands. Left over: the `ans:3` generated outputs still carry the pre-correction S table, filed in known-issues with its unblock condition. **Done-when:** steps 1–3 executed by the implementer, the closed-form identity and both negative controls asserted, elapsed recorded, the known-issues entry retired, §2.1's `PORT-1` line gains "S on the current route is `z_to_s(Z)` since `PORT-20`". *History: `docs/planning/chunks/PORT-20.md`.* | smoke (step 1), standard (steps 2–3) |
| `PORT-21` | What the lumped-sheet port adds to `S₁₁`: the self-class feed systematic, measured on our side | ✅ **2026-10-05 (weekly): T6 decided on F18a / F18b and §2.2 rewritten.** Steps 1, 3, 4 and 5 landed (step 1 closed 2026-09-21; steps 4 and 5 on 2026-10-01 as records); step 2's private ruling is tonight's: on exact rectangles, κ off, the degree-2 self class sits within 1 % of AED First Order on the ×1, h = 0.0075 and conductor-halved meshes (0.73 / 0.97 / 0.77 % worst class), so the feed systematic this chunk was opened to measure was the ragged sheet at degree 2 and element order at degree 1 (§10 rulings 2026-10-05; AED figures private). **Done-when:** `PORT-21` → ✅ when T6 is decided either way (clause (ii) at ×1; clause (i′) on the h = 0.0075 exact rung, `xl` entry 27 / F18b; clause (i″) on the conductor-halved rung, F18a — as decided 10-05) and §2.2 is written. *History: `docs/planning/chunks/PORT-21.md`.* | heavy |
| `PORT-22` | Find the tuned birdcage's resonance by a driven sweep | ✅ Step 1 closed 2026-09-20: the driven sweep across 44–84 MHz matches the circuit-layer prediction at every grid point, finds one sign change of Im Z_in bracketing 64 MHz, and the half-C control misses by far more than the band. Re-run green 2026-10-01 through the in-run κ option under `PORT-21` step 5; still one fixture at degree 1 with no absolute input-match claim. **Done-when:** (i)–(iii) and the control asserted by the implementer, elapsed recorded; §2.1 gains "the F-small tuned resonance is located by a driven sweep consistent with the circuit layer across 44–84 MHz" — still one fixture, degree 1, series resonance, no absolute `S₁₁` (known-issues 2026-09-19), no mode-spectrum claim (`TH-17`'s). *History: `docs/planning/chunks/PORT-22.md`.* | heavy |
| `PORT-23` | The circuit layer through `scikit-rf`, gated against the raw (C2) reduction | ✅ **Audit PASS 2026-10-02 10:30 review; the N1 / N2 prediction misses are ruled in the chunk file.** **Closed 2026-10-02 (09:00 slot, §9 item 155):** step 2 green — on `ans:7`'s 51 stored grid 4×4 S, skrf `Z_bare` matches Route B to 5.3e-15 / 7.9e-15 (Cu / PEC), rows agree with the raw route to ≤ 4.3e-13 and `DENSE_ROWS` to ≤ 1.5e-11; 20 passed at `-n 2`, 56 s (`20261002T140150Z_PORT-23-step2-W1.log`). Step 1 done 2026-09-23: the `scikit-rf` termination reproduces the raw (C2) reduction to machine precision on the stored 64 and 10 MHz matrices, with the swapped-port control asserted. Step 2, re-posed 2026-10-01 as reproducing Route B on `ans:7`'s stored grid matrices (§9 item 155). **Done-when:** on `ans:7`'s 51 stored grid 4×4 S (`logs/ans7/B1.json` / `B2.json`, Cu and PEC), `terminate` (P2–P4 in `C_F`) reproduces Route B's `Z_bare` (`route_b_z` on `z_from_s`) to `S_VS_Z_ROUTE_RTOL` at every point; `parallel_rows` on the skrf-side `Z_in` gives f₀, R_in, Q_loaded within 1e-9 of the raw route and within `test_ans7_sparse_route_b.BANDS` of `metrics.json`'s rows; harness, complex, `-n 2` ⇒ `PORT-23` ✅. *History: `docs/planning/chunks/PORT-23.md`.* | standard (`-n 2`, stored matrices; seconds) |
| `PORT-24` | *(history and §9 briefs: `docs/planning/chunks/PORT-24.md`)* **Capacitor ESR: a 0.2 Ω series resistance on every lumped capacitor, in the circuit layer and in-model** (operator binding ruling 2026-10-06, §9; opened by the 2026-10-07 weekly; §10 chain T25 / T26, H17). *Step 1 (circuit layer, no solve):* an `esr_ohm` argument (default 0) on `ring_mode_rows` / `ring_mode_rows_v2` and on the circuit-layer terminations (`Z = r + 1/(jωC)`). *Step 2 (in-model):* a series-RC sheet carrying the same r and C through `series_rlc_impedance(R=r, L=0, C)`. **Done-when (§4):** *step 1* — with `esr_ohm` unset every stored record reproduces bitwise (the F-small 4×4 at 10 / 64 / 128 MHz, F-human's saved degree-1 32×32); at r = 0.2 Ω, Re λ_k(r) − Re λ_k(0) = r for every mode k at 1e-12 relative while Im λ_k and `C_mode` are unchanged within 1e-12 relative (re-stated from "bitwise" by the 2026-10-10 weekly, T25); a planted non-uniform r (one gap at 0.4 Ω) breaks the uniform shift on every mode ⇒ 🟡. *Step 2* — on the F-small exact four-port at 64 MHz, degree 1, `-n 2`: the sheet with (r, C) reproduces `series_rlc_impedance(0.2, 0, C)` at `PORT-14`'s imported band and with r = 0 reproduces the capacitor-only sheet's 4×4 at `LEG_D0_REPRODUCTION_BAND`; the exact power identity holds with the sheet's own dissipation counted; printed, Q(r)/Q_0 against 1/Q = 1/Q_0 + 1/Q_C (Q_0: the lossless-capacitor coil's loaded Q) ⇒ ✅. Lossless records are never re-read with r. Not claimed: a measured ESR, frequency dependence, any AED re-run. **Step 2 green 2026-10-08** (item 309, `91254c5a`; readings in `PORT-24.md`). **✅ 2026-10-11 (item 363, T25′): step 1 landed; (S)-Im ≤ 2.889e-15 / `C_mode` ≤ 4.571e-16 against 1e-12; `20261011T093040Z_PORT-24-step1-T25p.log`; audited PASS 10:30.** | ✅ | standard (step 1 smoke, no solve; step 2 measured heavy-ceiling, 91–140 s per arm at `-n 2`) |

**`PORT-16` (the accounting-gap row) — rulings, 2026-09-07 03:00 review: step 1 audited PASS; `POWER_BALANCE_BAND` ruled a record; step 2 rescoped; step 3 named.**
*Audit.* `auditor` PASS on all eight checks at `a2f8db0`: headers `fcfd101` = the closer's parent; 20 of 20 digits traced (`20260907T051231Z_PORT-16.log:1882–1893, 1895–1910, 1912–1915, 1917–1920, 1922–1924, 1997, 2065–2066`; window 1's sign signature `rel dev 2.000e+00` at `…050816Z:1882, 1885, 1888, 1891, 1917`); (i) and (iv) are executed asserts (`test_birdcage_power_identity.py::test_the_discrete_power_identity_closes_on_every_drive`, `::test_the_cauchy_schwarz_deficit_of_the_terminal_form_reproduces_the_gap`); `POWER_BALANCE_BAND` is imported nowhere in the module; `test_port_drive_superposition.py` is absent from the diff; no `src/`; 131 s inside standard. Re-traced by this review: `:1882` `rel dev 6.760e-15`; `:1917` `C - sheets_terminal 6.716202469e-05 W … rel dev 3.212e-13`. Both derivation repairs are ratified — the driven sheet's field-only form is not a theorem (its terminal current is driven by `E_t + E_src ĥ`), and a sign fixed in the test's own helper with the red window committed beside the green one is the MAG-10 / MAG-15 precedent. The control's 0.484102 outside the predicted [0.5, 2] window was printed, not asserted, exactly as rule (e) requires; the prediction was wrong by 3% and is not relied on anywhere.
*Ruling (1) — `POWER_BALANCE_BAND` is a record, not a conservation band.* The quantity it bounds, `|P_acc − P_vol|/P_acc` with `P_acc` built from terminal currents, is not an identity for the discrete solution: it misses by the terminal form's Cauchy–Schwarz deficit, now measured at 1.0106× the terminal sheet dissipation on every drive and reproduced to 1e-13. The three green single-drive asserts (`test_birdcage_b1_plus_map.py:432`, `test_port_birdcage_ring_matrix.py:273, 466`, `test_port_drive_superposition.py:682, 713`) **stay exactly as they are** — nothing green is touched — but are henceforth read as records of a ≈1% fixed offset at fixed `h`; a future red on any of them is attributed through step 1's split before anyone re-bands it. The **drive-level** assertion (`test_port_drive_superposition.py::test_the_drive_level_power_identity_closes`, the deliberate red at 11.648%) was asserting a fixed absolute deficit against a 5× smaller denominator; it is re-pointed at the identity that *is* exact — step 2. This is not a loosened bound: the replacement asserts two identities at 1e-6 and 1e-1 (measured 1e-14 / 1e-13 class on the singles) on the same solve where the old line asserted a non-identity at 1e-2, and the old residual stays printed with the band beside it. The `POST-6` known-issues entry stays open until step 2's commit retires it. The 18:00 review had listed this ruling for the 09-09 weekly; step 1's attribution is decisive enough to act on now, and the weekly may overrule it with step 3 in hand.
*Ruling (2) — step 2 rescoped* (§9 item 2, full item there): the exact identity (i) at the imported `DISCRETE_IDENTITY_RTOL` and the attribution (ii) at the imported `ATTRIBUTION_RTOL` on the **superposed** ccw / cw drives in the module that holds the red, with `w = e_k` reproducing step 1's single-drive integrals to 1e-12 as the control of the generalisation; the row closes ✅ when it lands, and its example is `ports:9`'s power table — no new one. *Step 3 (named, not queued):* the `h`-exponent of `C/sheets_terminal − 1` on `PORT-14` step 1b's ×0.75 / ×0.6 conductor rungs — the only physics question left in this row, worth a slot only if the weekly wants the offset's `h`-dependence on record before ruling on the band's future.
*Audit of the closure (step 2), 2026-09-07 10:30 review — PASS on all eight checks at `7a82688`.* Header `2f65cbe` = the closer's parent (`20260907T110826Z_PORT-16.log:1–16`); 10 of 10 digits traced (`:1942–1943`, `:1945–1946`, `:1948–1951`, `:2038`, `:2106–2107`); the disposed red's baseline `20260905T004202Z_POST-6.log` sits on `2d8ce09`, a strict ancestor; the anchors are executed asserts (`test_port_drive_superposition.py:959`, `:974`) against bounds imported from step 1's module; `POWER_BALANCE_BAND` (`test_birdcage_b1_plus_map.py:101`, 1e-2) is absent from the commit's diff and its single-drive assert at `test_port_drive_superposition.py:895–896` is unedited; no `src/`, `test_birdcage_power_identity.py` untouched; 104 s inside standard. Re-traced by this review: `:1942` `rel dev 5.877e-15`, `:1945` `P_acc − P_vol 3.511221378e-04 W … C − sheets_terminal 3.511221378e-04 W`. The auditor's one caveat — the `timeout -k 30 400` wrapper exceeds the 180 s standard figure — is the step-1 precedent (tier by measured elapsed) and changes nothing. Step 3 stays the weekly's; not queued.

**`PORT-12` step 2 — width-qualify the record and bound the parallel drift
(ruled 2026-08-30 02:15 weekly review, option (i) with a bounded envelope;
ruling text in the known-issues entry — that session died before writing
this paragraph, the 10:30 daily review wrote it from the ruling).** Smoke
tier, complex, `main`, **no `src/` change**: `tests/validation/
test_port_lumped_two_torus.py` only. **Build:** `REPRODUCTION_BAND = 1e-4`
stays and its comment states it is a **`-n 2` record** (every
`PORT-1`/`OPS-18` two-torus digit is quoted at that width); a new
pre-registered `PARALLEL_DRIFT_ENVELOPE = 3e-4` (provenance: step 1's
four-width table, max +2.06e-04 at `-n 8`, 1.46× headroom) is the band
`test_step_1_measurements_reproduce` uses when `comm.size > 2`, so the drift
is bounded on every width rather than hidden by a skip; and the lumped
route's width-flatness becomes the module's new negative control — assert
`Im Z12(lumped)` reproduces **1.029281338** at rtol **1e-8** (measured
flat to 2e-9 across all four widths), which is the reason the production
port model is unaffected. **Anchors, asserted:** at `-n 2` the gap ratio
**0.894141** inside 1e-4 (unchanged); at `-n 8` — the worst width — the gap
ratio inside 3e-4 of the record with the measured **+2.06e-04** printed;
`Im Z12(lumped)` at 1e-8 on both widths; the module's other four tests
unchanged. **Negative control:** run the `-n 8` window on `main` **first**
and footer its Status 1 (the existing red, `1 failed, 4 passed`), then the
patched tree green — the envelope must be visibly load-bearing; and the
1e-8 lumped assert must *fail* if pointed at the gap-route `Im Z12`
1.110303775 (a one-line probe, not committed). **Cost:** 84 s per width
(`20260829T1700…Z_PORT-12-step1-*`), three windows ≈ **260 s**, `-k 30
300` each, `tests/environment` first. **Traps:** complex build +
`FEM_EM_REQUIRE_COMPLEX=1`; `-n 8` is the hard case, not `-n 12`; `-s`;
do not touch `STEP1_GAP_RATIO_RECORD`. **Scope:** closes the chunk ✅ and
retires the known-issues entry; **not** a root-cause fix (option (iii)
declined — no birdcage or Larmor quantity reads a gap-route integral) and
**not** a widened record (option (ii) declined). **Negative result:** the
`-n 8` drift reading outside 3e-4, or the lumped route moving at 1e-8, is
a new fact about the fixture — record it on the known-issues entry, do not
widen the envelope, stop.


**`PORT-14` narrative** — moved byte for byte to `docs/planning/chunks/PORT-14.md` (`OPS-47`).

**`PORT-15` — the circuit layer (HFSS + Circuit)** ✅ *(closed 2026-09-14
02:07Z on gate (i) + the tuned `S₁₁`; result block written by the 2026-09-19
weekly, full narrative in `docs/planning/plan-archive.md`, entry «§7 PORT-15
full narrative — archived 2026-09-19 (weekly review)»; row history stays in
`docs/planning/chunks/PORT-15.md`.)* Digits in the §7 table row
(`20260914T010448Z_PORT-15.log`, `20260914T020419Z_PORT-15.log`). **Live
carry-forwards:** gate (ii), the 32-port mode spectrum, is `TH-17`'s and is
not claimed; "tuned" means series resonance on one F-small 4-leg fixture at
64 MHz, not matched.

**`PORT-1` — Real port excitation from the solved field** ✅ *(closed by the 2026-08-15 18:00 review. Full
plans, journals and adjudications for steps 1–4 and 3b(i–xviii) were archived on 2026-08-15; the 2026-08-15
compressed record itself is now in `docs/planning/plan-archive.md`, entry «§7 PORT-1 compressed record (2026-08-15 form) — archived 2026-08-30 (weekly review)» — grep there before re-deriving.)*
> * **Done-when, met:** `run_n_port_sparameter_sweep` reads the solved field — `‖S−Sᵀ‖/‖S‖ = 2.5494e-05` vs the 1e-3
>   gate, `‖S‖₂ = 0.861449`, `is_placeholder=False` (step 4, `20260813T183606Z_PORT-1-step4-packagegate.log`, 7 passed,
>   153.9 s); retiring heuristic differs by 3.078e-01. **Two-torus fixture through this entry point only** — not
>   birdcage ports, not S11/Z_in, not B1+. (Later re-records: 0.11 image 3.11213e-05 / 0.861356895 under `OPS-18` (1\*);
>   power-wave route 4.758625e-05 / 0.864809457 under `PORT-9` leg (d3).)
> * **Ladder, compressed:** f = 10 MHz, `ωM₁₂ = +1.241755 Ω`; reaction-route reciprocity 2.65e-13; electric-energy
>   excess removed by the solenoidal projection (ratio 0.999998, production drive since 2f); 3b found the factor 2 in
>   `_gap_arc_quadrature` (port voltage 0.894543 × ωM₁₂) and landed 3b-xvii (Faraday closure, 11× margin) and 3b-xviii
>   (`Im Z₁₂` vs the filamentary form at the unmoved 10%: corrected **0.939581, −6.04%**).
> * **Systematics (`ports/systematics.py`):** PEC box `D∞ = +0.0169` at `p = 1.657` (effective-range, never quoted
>   without its exponent); gap physics ÷(1 − 0.030224) (Jin 3e §10.4.2.1); composition additive (`PORT-10`, −0.0604 pp).
> * **Standing cautions:** no `Z_in`/`S₁₁` off the unprojected diagonal (`W_e/W_m = 6.524`); mutual is always the undriven port; `MUTUAL_TOLERANCE = 0.10` is measurement-justified — do not tighten; unconditional `create_entity_permutations()` before subdomain `dS` integrals (3b-iv hang); known-issues 3 (defect 1) and 11 stay open.

> **Two port tests are red and deliberately left red.** Both fakes set
> `current = voltage/z0` at the driven port, making it perfectly matched, so
> `b = (V − Z₀I)/2√Z₀ = 0` and the S-matrix diagonal is legitimately zero against
> an assertion demanding non-zero. Fixing them means tuning assertions to match a
> heuristic that `PORT-1` deletes; resolve them there.

> `PORT-3` was 🚫; its recorded failure was a docker preflight error, not a code
> failure, and `OPS-1` resolved that. Real status is unknown — hence 🧪, not ✅.

> `PORT-0` quarantine: `PlaceholderPortModelWarning` on every call, `is_placeholder`
> threaded through the result types, and `export_touchstone()` refuses flagged data
> unless `allow_placeholder=True`. Fabricated `.s2p` files can no longer leave the
> project looking authoritative. `PORT-6`/`PORT-7` are 🧪 rather than ⚠️ — sweep-grid
> generation and Touchstone *formatting* are correct independent of what fills the
> matrix.
> **Amended 2026-09-29 (`PORT-7`):** the formatting claim held for N ≤ 2 only — the
> writer listed N ≥ 3 column-major, which a Touchstone 1.x reader (scikit-rf 2.1.0)
> reads as Sᵀ; invisible on reciprocal data. Fixed and cross-checked by `PORT-7`.

**`PORT-5` step 1 — sweep-level sanity metrics on the field route** ✅
*(2026-08-16, `20260816T093556Z_PORT-5-step1-rerun.log`, 10 passed
149.1 s at `-n 2`, standard, wrap `timeout -k 30 500`.)* The metrics had
only ever seen placeholder or hand-built matrices — §10 target 3's
"`PORT-5`'s sweep-level path is untouched". Three cases now ride the
`PORT-1` step-4 module's own fixture (`tests/validation/test_port_package_sparameters.py`;
same module-scoped sweep pair, **no extra solves** — the metrics are pure
numpy). Measured on the report `run_n_port_sparameter_sweep` *returns*:
`passivity_max_sigma` = **0.861449197** against step 4's gated
`‖S‖₂ = 0.861449`, miss **1.97e-07** inside the pre-stated 1e-6, and
equal to `numpy.linalg.norm(S, 2)` to < 1e-12 (same quantity, two
implementations); `reciprocity_max_abs_delta` = 2.194793e-05, which for a
2×2 converts exactly (`‖S−Sᵀ‖_F = √2·max|Sᵢⱼ−Sⱼᵢ|`) to
**‖S−Sᵀ‖/‖S‖ = 2.549409e-05** against the gated 2.5494e-05, band 5e-7;
`passivity_max_column_power_sum` = 0.741345553 ≤ 1 (the second metric
step 4 never read); **no warnings** on the field route.
> **Negative controls, both executed.** The deprecated heuristic's S
> through the same metrics: `passivity_max_sigma` = **0.999985964171**,
> separation from the field route's **0.138537** > the pre-stated 0.13,
> `reciprocity_max_abs_delta` identically 0. An S with one off-diagonal
> perturbed by 2× the warning threshold: delta 9.999344e-02, both
> reciprocity warnings fire, and the untouched matrix still reports none
> — "no warnings" is evidence only because a warning can fire.
> **One §9 constant was wrong and is corrected with its measurement.**
> The item quoted the heuristic's `passivity_max_sigma` as exactly
> `1.000000000000`; that number is the *reaction-route* fixture's
> (`PORT-1` step 2 iv) and the hand-built unitary S of
> `test_port_reaction_impedance.py`, not this mesh's. Measured here:
> 0.999985964171, unitary to 1.4036e-05 (first run,
> `20260816T093226Z_PORT-5-step1.log`, 1 failed / 9 passed — the two
> anchor cases passed at their pre-stated bands in that same run). The
> premise assertion now reads "unitary to 5e-5" with the measurement in a
> code comment; the *discriminating* assertion, the 0.13 separation, was
> never moved.
>
> **Scope.** Metrics wiring only: no tolerance in `sparameters.py` moved,
> `PORT-5`'s frequency-sweep ambitions stay unscoped, and the claim is the
> two-torus fixture through this entry point — every `PORT-1` standing
> caution above applies unchanged.

**`PORT-9` — lumped-element port boundary condition (the birdcage port model)** ✅ *(closed 2026-08-25 06:00
slot, leg (d1′), **at 10 MHz on the gapped 4-leg birdcage**; scoped 2026-08-16 weekly review, Jin 3e ch. 11 port
hierarchy. Steps 1–2c (two-torus) archived 2026-08-23; step 3 legs (a)–(d1′) with rulings (2\*)/(4\*)/(5\*)/(6\*)
archived in `docs/planning/plan-archive.md`, entry «§7 PORT-9 full narrative (step 3 legs (a)–(d1′)) — archived 2026-08-30 (weekly review)».)*
> * **Two-torus (steps 1–2c):** resistive-sheet lumped port (`ports/lumped.py`), circuit identity to < 1e-12; narrowed
>   sheet f = 0.5 with `w = A/h` gates cross-route **1.8333%** vs 5%; sweep route `‖S−Sᵀ‖/‖S‖` 2.574249e-11.
> * **Step 3 on the birdcage** (mesh from `GEO-18`, then `GEO-19` step B's 116 085 cells; `Z_p = z0 = 50 Ω`, f = 0.5):
>   leg (c) one solve **7.55 s**, adjacent spread 0.0159%; (d0) termination margin **598.4002×** at 50 Ω (open 1e6 Ω
>   column near-degenerate — `Z_p` is the gap capacitor); (d) 4×4 reciprocity 2.495292352e-05, σ_max 0.862659137;
>   (d1) displaced rung lost reciprocity (5.57e-03) → (d2) readout *is* the source's adjoint (1.33e-10), the asymmetry
>   is the terminated-`Z` assembly → ruling (2\*): **power-wave S**; (d3) asymmetric two-torus 1.324004669e-16 vs old
>   1.143811489e-04 (9.525277e+11× separation); (d3b/d3c) birdcage fixed route 4.6e-15–1.2e-14, σ_max **0.999993391**,
>   spreads 0.0617 / 0.0359 / 0.0237%, cell record 116 416 → 116 368 (image) → 116 085 (step B), all re-records version-tagged under (1\*)/(5\*).
> * **Closing leg (d1′), `20260825T110438Z_PORT-9-step3d1.log` (13 passed 106.64 s) + `…110643Z_…-consumers.log`
>   (24 passed 222.15 s):** zero rung reproduces (d0)'s column ≤ 2.568e-10, σ_max **0.999992805**, spreads **0.0553 / 0.0353 / 0.0214%**, separation 166.6766×, margin 2256.9707×; leg 1 rotated 22.5° breaks (iii′) on all three classes — **6.2219 / 7.1142 / 2.8474%** vs 0.5% — while `‖S−Sᵀ‖/‖S‖` = 2.259e-14 (2.466e+11× from the pre-fix 5.57e-03). (iii′) 5% → 0.5% committed as a tightening, all consumers green. * **Carry-forward:** **no Larmor, resonance or tuning claim** (`PORT-11` owns 64/128 MHz); terminated `Z` is a diagnostic, never reciprocity-gated; power-wave residuals are order-of-magnitude only ((d3c) rule); open-limit conditioning known-issues entry OPEN; `GapVoltagePortSpec` cannot run on the `GEO-16`-fragmented mesh.

**`PORT-10` — the two `PORT-1` systematics: composition measured, not
assumed** ✅ *(2026-08-16; full 49-line entry archived in
`docs/planning/plan-archive.md`, 2026-09-06 weekly review).* Result:
cross-term **−0.0604 pp** against the pre-stated ±0.5 pp band — additive,
the sequential ladder in `ports/systematics.py` stands as measured; it
closed §10 Phase 5 subgoal 1's last open question. Live carry-forward:
none. Logs: `20260816T140457Z_PORT-10-costprobe.log`,
`20260816T140643Z_PORT-10.log`.

**`PORT-11` — lumped-sheet ports on the gapped birdcage at the Larmor frequencies** ✅ **2026-08-26** *(step 1
probe 08-25 19:30, step 2 64 MHz 08-26 06:00, step 3 128 MHz 08-26 16:30; both audited COMPLIANT by the 18:00
review; commissioned 2026-08-23 weekly review as §10 subgoal 2b, serial on `PORT-9`. Full narrative archived in
`docs/planning/plan-archive.md`, entry «§7 PORT-11 full narrative — archived 2026-08-30 (weekly review)».)* Same code path and fixture as `PORT-9` leg (d) — only the frequency moved.
> * **Step 1 (`20260826T003427Z/003559Z_PORT-11-step1*.log`, 14 passed twice):** phantom cells/δ **5.9213** (stop rule
>   ≥ 2), cells/λ 21.8936, loss tangent 1.8004; MUMPS frequency-flat (9.49 / 6.36 s vs 6.56 / 6.50 s at 10 MHz).
> * **Step 2, 64 MHz (`20260826T110434Z_PORT-11-step2.log`, 17 passed 177.48 s, 179 s):** (i) `‖S−Sᵀ‖/‖S‖`
>   **2.581325834e-14** vs 1e-3; (ii) σ_max **0.999721388** ≤ 1 + 1e-9, max column power 0.804704664; (iii′) class spreads
>   **0.0573 / 0.0599 / 0.0370%** vs 0.5%, anti-noise control 671.0527× vs 10×; 10 MHz control reproduces
>   `LEG_D_S_MATRIX_10MHZ` to **1.158e-10** vs 1e-6 and `LEG_D0_Z_COLUMN` to 2.568e-10; displaced control breaks (iii′) at
>   12.8947 / 27.7509 / 7.7239% with (i) at 1.25e-15; consumer `…110750Z` 5 passed 103.82 s.
> * **Step 3, 128 MHz (`20260826T213414Z_PORT-11-step3.log`, 18 passed 197.85 s, 201 s):** pre-gate rule by measurement —
>   loss tangent **0.9002** (displacement-dominated), **cells/λ 12.5024 ≥ 10**, cells/δ 5.1845; (i) **7.030990825e-15**;
>   (ii) σ_max **0.998974779**, column power 0.861668762; (iii′) **0.1012 / 0.0916 / 0.0654%**, control 576.9483×;
>   displaced 16.7006 / 34.6556 / 13.2091%; consumer `…213748Z` 16 passed 130.04 s. * **Carry-forward:** **self-consistency identities on one fixture — no resonance, tuning or absolute-accuracy claim**; feed systematics are the two-torus ones (`PORT-1` 3b-xviii, `PORT-10`); `|Im P|/Re P` (terminal power 1.755210 / 2.659902 at 64 / 128 MHz) printed never gated; no vessel-wall region; `ANS-4` handed to the 08-30 weekly review.

### WF — End-to-end workflow & MRI outputs (Phase 5)

| ID | Title | Status | Tier |
|---|---|---|---|
| `WF-1` | `femem run`: the MRI study CLI and YAML config | ✅ (item 256, U5): `femem run` builds Studies equal to the Python-built ones; every later YAML step, through `thermal.transient` on the EM sphere (step 29, U134) and `thermal.dose` / `DoseReadout` on the walled box (step 30, U139) — public keys for the 11-09 weekly to ratify — is a record in the chunk file. *History: `docs/planning/chunks/WF-1.md`.* | heavy by ceiling (four-port window 400 s; others standard / smoke) |
| `WF-2` | Reproducible output bundle manifest *(brief: `docs/planning/chunks/WF-2.md`)* | ✅ 2026-10-09 (00:00 slot, §9 item 346; re-pose ratified by the 2026-10-10 weekly; audit PASS 2026-10-11 03:00 review): `tests/package/test_wf2_bundle_reproducible.py` green at `-n 2` — W1 sphere 6 passed 48 s (`20261009T053015Z_WF-2-W1-sphere.log`), W2 four-port 5 passed 77 s (`20261009T053110Z_WF-2-W2-four-port.log`); (D) / (F) / (T) bitwise, (R) 1.453e-12 / 5.010e-16, (X) ≤ 6.222e-11 / 9.557e-11, all ≤ `L0_BAND`; (c1) / (c2) / (c3) / (c4) fire (one-ULP σ caught by `sha256` alone, SAR rel 7.566e-12). Was: ⚠️ Lane-1 re-pose queued as §9 item 346. **Done-when (§4) (written 2026-10-08 18:00 review, rule (j); ratified by the 2026-10-10 weekly, including reading "the manifest" as `run_study_file`'s bundle):** two `femem run`s of the same study file at `-n 2`, complex, degree 1 — `WF-1`'s lossy-sphere and F-small four-port files — write bundles whose deterministic fields (config `sha256`, family, degree, frequency, ranks, cells, Touchstone header bar the timestamp) are bitwise equal, and whose S and SAR agree with each other and with `MAT4_DEG1_RECORD` / `P1["64"]` within `L0_BAND` (the measured `-n 2` drift ≤ 4.45e-11; the old 1e-12 sat below it); a planted degree change, a P1↔P2 swap and a one-ULP σ config are each detected; harness, elapsed recorded ⇒ ✅. Not claimed: cross-width reproduction, degree 2, the `mri:1` example manifest, Z. | heavy (by ceiling) |
| `WF-3` | Quick-look phantom metrics report. Done-when (written 2026-09-30 03:00 review, rule (j); queued §9 lane 1). On `TH-6`'s plane wave at 24³ (`POST-2`'s fixture), … | ✅ Closed 2026-09-30 and audited PASS at the 10:30 review: the eight quick-look numbers parsed from the written markdown match the lossy plane-wave closed form within 5 % and the JSON summary, with the σ-blind and swap controls missing. The balance line, the status line on physical fields and SAR are not claimed. **Done-when:** the eight numbers parsed from the written markdown (\|E\| and \|B\| min / max / mean, the two shape ratios) are each within 5 % of the closed form on the owned centroids (e^{−αx_c}; × \|k\|/ω; ω/\|k\|); each equals its JSON `summary` value to 5e-7 relative; the JSON `summary` equals the metrics dict exactly; coverage reads 82 944/82 944/82 944 with no drops and no report warning flag; the σ-blind comparand misses by > 5 %; `POST-2`'s test is green in the same window. *History: `docs/planning/chunks/WF-3.md`.* | standard |
| `WF-4` | Scenario presets (debug / dev / benchmark-lite) | ✅ (§9 item 471b): three presets under `examples/presets/` each reproduce their named record through `femem run` at its imported band, and `femem run --estimate` prices each inside ×2 of its own run; only `benchmark-lite` tests the fit (the `debug` / `dev` checks are reproductions). **Done-when:** three YAML presets under `examples/presets/` that `femem run` accepts, each reproducing one named record at its imported band through the CLI — `debug` = `MAT4_DEG1_RECORD` (degree 1, `-n 2`, `L0_BAND`); `dev` = item 192's exact-rectangle four-port S `_record(1)` (degree 1, `-n 2`, `LEG_D0_REPRODUCTION_BAND`, `PORT-11`'s gates green); `benchmark-lite` = `ANS6_EXACT_STRIP_CLASS_MEANS_P2["64"]["pec"]` on the exact-rectangle PEC hole (degree 2, `-n 8`, `LEG_D0_REPRODUCTION_BAND`) — and `femem run --estimate` on each prints a `predicted_seconds` inside `BRACKET` (×2) of that preset's own `summary["elapsed_s"]` in the same window; asserted through the harness, elapsed recorded ⇒ ✅. *History: `docs/planning/chunks/WF-4.md`.* | heavy (`-n 8` window, ≈ 85 s measured on the Study route) |
| `WF-5` | Loaded birdcage: frequency shift & Q degradation | ✅ All three steps were asserted or recorded on 2026-10-01 and the Done-when is met, audited PASS at the 18:00 review with caveats (a thin A3 margin, the PEC control missing its tighter prediction). No h-convergence, absolute-Q, AED or tuning claim; the weekly still reconciles §10 T15 with the printed-only Δf₀. **Done-when:** step 1 asserted ⇒ 🟡; steps 1–2 asserted and 3 recorded, elapsed recorded ⇒ ✅. *History: `docs/planning/chunks/WF-5.md`.* | heavy |
| `WF-6` | B1+ field mapping and homogeneity (CV) | ✅ 2026-09-13 (step 5) on the re-scoped subgoal-4 target: on the unloaded F-small birdcage at 10 MHz, CG1, the worst-radius C4 four-copy spread of `\|B₁⁺\|` falls 5.2506 → 2.0719 % and the C4 covariance 3.6159 → 1.6815 % as resolution goes 0.015 → 0.012 m — both falls asserted, both spreads reproduced at rtol 1e-3 (`20260913T170259Z_WF-6-step5.log`, 24 passed, Status 0, 156 s at `-n 4`). A convergence statement only — no closed-form, homogeneity, absolute, C95.3 or Larmor claim; the flag-on ×0.0095 power residual stays banked OPEN in known-issues. *History: `docs/planning/chunks/WF-6.md`.* | heavy (step 1 standard, complex; steps 3h 194 s / 3f′ 117 s / 3i 96 s standard; step 4a 4 s smoke, `-n 1`, real; **step 4g 204 s at `-n 4`**; **step 5 156 s at `-n 4`**; **step 4h 199 s at `-n 4`** + 52 s flag-off control; **step 4i 198 s + 104 s at `-n 4`**; **step 4j 199 s + 101 s at `-n 4`**; step 4k 29 + 32 + 39 + 38 s smoke, `-n 1`, no solve) |
| `WF-7` | SAR10g hotspot identification | ✅ Closed at step 3 on 2026-09-25 (audited PASS 2026-09-26; *2026-10-10, `GEO-35` step 2a: steps 1–2 records now on the leg-box-free fixture, boxed values kept as `*_BOXED` and rebuildable with `FEM_EM_WF7_LEG_PORTS=boxed`*; *2026-10-10, `GEO-35` step 2b-i (§9 item 356): step 3's 10 g C16 (0.3373 %), mirror (0.0234 %), coverage (≤ 1.6e-13) and containment checks green on the leg-box-free refined mesh, 532 935 cells, `20261010T111003Z_GEO-35-s2bi-W2.log:10825–10834`; step 3's records still boxed until 2b-ii*; *2026-10-11, `GEO-35` 2b-ii part 1 (item 385a): step 3 re-recorded on the leg-box-free default, 532 935 / 22 822 cells, SAR10g(c\*) 2.079992000e-04, `20261012T020246Z_GEO-35-s2bii-W1.log`*): the 10 g SAR hotspot is located at human scale with the coverage, C16 and mirror identities green on one phantom at degree 1, 64 MHz. The degree-2 full 32×32 ran on 2026-09-26 as a record, and no compliance or absolute SAR claim is made. **Done-when:** (ii) the **10 g hotspot** — the argmax of `mass_averaged_sar` over a fixed candidate lattice — and its fifteen **C16 images** within the imported `C4_COVARIANCE_BAND` (5 %); (iii) the **mirror identity** — the cw drive's 10 g SAR at the mirrored hotspot vs the ccw drive's at the hotspot — within the same band; (iv) every ball's containment in the phantom. Knobs `FEM_EM_WF7_SAVE_NPZ` and `FEM_EM_WF7_ORIENTATION_CONTROL` (unset = inline) landed 2026-10-04 (§10 H15a, §9 item 231; records only). *History: `docs/planning/chunks/WF-7.md`.* | heavy; step 0b **`xxl`** (measured 673 s / 106.1 GiB, 2026-09-19); step 0c **`xl`** (cost probe); step 0d **`xxl`** |
| `WF-8` | Figures that equal their records — the visualization layer of Phase 9 (9d) | ✅ 2026-10-08 (§9 item 328): `post.figures` + `femem figure` — every plotted array bitwise equal to its source record (the S lines, the SAR slab, the probe series), the planted-wrong-array controls fire on every figure, `--check-links` dead = 0. **Done-when:** in one logged window at `-n 2`, each figure's plotted arrays are read back from the figure object and asserted **bitwise equal** to the source record they draw (the `.s4p` through `scikit-rf`, the XDMF array, the CSV column); a planted wrong-array control (one value perturbed at 1e-6) is detected on every figure; the files exist and are non-empty; `--check-links` over `docs/guide` stays at 0 dead ⇒ ✅. *History: `docs/planning/chunks/WF-8.md`.* | standard |
| `WF-10` | The tuned human-scale birdcage: tune F-human in the circuit layer, then B₁⁺ and the 10 g SAR hotspot | ✅ (audit PASS on steps 1, 1e and 2): the in-model Larmor root on copper with r = 0.2 Ω and the feed-bisector mirror gates at degree 1 — no tuning, absolute or compliance claim; steps 3a–5 and the degree-2 records (H18, H20a, H24, H26a, H27, H28: the degree-2 root C_in(d2) = 115.5209 pF and its node knob `FEM_EM_WF10_NODE=d2root`) are records. **Step 4 at the degree-2 root** (`xl` entry 41, H28′, 2026-10-22; degree 2, 2 570 024 unknowns, every gate green at C_in(d2) = 115.5209 pF): mean \|B₁⁺\|_ccw 2.792653e-07 T, SAR10g(c\*) 4.547175546e-03 W/kg at r = 45 mm, z_c, az₁₇ (5.402616e-02 W/kg/W of P_src), (M_B) 1.7093e-03, (iii-a) / (iii-b) 0.3052 / 0.3277 %; root-to-root −3.17 % / −6.30 % (per phantom watt +0.60 / +1.13 %) ⇒ the rule's second branch, **ruled by the 10-24 weekly** (H31's carry ratified; Phase 6 exited on it): the human-scale tuned figures are quoted at degree 2 only, at C_in(d2); the degree-1 figures are a record with the order move beside them (`20261022T070008Z_WF-10-step4-H28.log`). **Done-when:** steps 1, 1e and 2 asserted, elapsed recorded ⇒ ✅. *History: `docs/planning/chunks/WF-10.md`.* | heavy |
| `WF-11` | **The 3 T body coil: a 32-leg (64 ring-gap port) high-pass birdcage at human scale, tuned at 128 MHz, with `\|B₁⁺\|` and the 10 g SAR hotspot** (operator directive 2026-09-28, interactive session: on the list for the long horizon, **not yet** — the 16-leg / 32-port F-human at 64 MHz, `WF-10`, is the starting point for human-scale coils). The 128 MHz twin of `WF-7` + `WF-10` on `GEO-33`'s 32-leg rung, the geometry class of the 3 T research coil the project is ultimately compared against. **Serial on:** `GEO-33` (the 32-leg rung, strip conductor as two-sided Leontovich sheets — chain H6a / H6), `WF-10` ✅ (the tuning route and its mode selection proven at 64 MHz on 16 legs). **Parameters are generic by `GEO-33`'s rule** — leg count 32, a coil length and radius in the published clinical-birdcage range, the phantom's σ / εᵣ at 128 MHz from the literature (§11); nothing from any specific product. **Known physics to plan for, not a claim:** at 128 MHz the wavelength in a saline-class phantom is ≈ 26 cm, comparable to the body, so `\|B₁⁺\|` is inhomogeneous by physics (standing-wave / dielectric effects) — the gates are symmetry and circuit↔field identities as in `WF-7` / `WF-10`, never homogeneity; degree order is re-priced at 128 MHz (`WF-7` 0b measured 5.5 % at 64 MHz, expected larger at 128). Steps are written by the weekly when `GEO-33` step 1 has a price. **Not claimed:** compliance, a match to any product, a second phantom. | ⬜ Future work, undated (operator directive 2026-09-28). | unmeasured (serial on `GEO-33`'s cost probe) |
| `WF-12` | F-human on copper: the 32-port human-scale birdcage as a hole mesh with the Leontovich wall | 🟡 Steps 1–2 green (CAD identities, the Cu drive's exact power partition, the copper 32×32 at degree 1, 64 MHz); step 3 (H12c) is a pre-registered record — the identity at copper's `C_mode` 118.1186 pF reads 1.581167 inside [1.5, 1.7] — and step 3b (H12d) located the miss at the termination level (max \|δ̂\| 1.30e-2). Step 4 (H16, `xl` entry 32, 2026-10-13, a record): one degree-2 copper column, 2 570 024 unknowns, 514 s / 89.3 GiB at `-n 16`, (P) 1.468e-13 and NEG 119× asserted, S17,17 d1 → d2 move 3.550 % (`20261013T070008Z_WF-12-step4-H16.log:6711–6720`); step 5a (H16b, §9 item 415): `tests/validation/test_wf12_f_human_copper_s32_degree2.py` with `FEM_EM_WF12_STEP5` (`build` / `all`) and `FEM_EM_WF12_STEP5_SAVE_NPZ` (`-d2.npz` under `logs/wf12/`, degree-1 records refused); (A0), (Deg) 2 570 024 at `-n 2` (`20261013T111749Z_WF-12-step5a-W1.log`). **Done-when:** steps 1–3 ⇒ ✅, step 3 counting only with its circuit ↔ field identity at copper's `C_mode` green *History: `docs/planning/chunks/WF-12.md`.* | heavy (measured: step 1 W1 245 s `-n 1`, W2 199 s / 9.6 GiB `-n 8`; step 2 W2 311 s / 10.2 GiB `-n 8`, held drive ≤ 3.57 s) |
| `WF-9` | Implant transfer-function excitation (ISO 10974 Tier 3) | ✅ As a fixture 2026-10-04 (step 1, I8; audit PASS): the σ = 800 S/m rod in the coil-driven F-small phantom at 64 MHz — the adjoint fold equals the direct tip response to ≤ 3.7e-12, the N = 32 path fold to ≤ 1.47e-3 against 5 %, the Born fold misses. Step 2 (I16, 2026-10-06) repeats the identities at degree 2 on the same mesh with the TF order move 2.6–29.1 % per drive printed — two orders on one mesh, not a convergence claim. Step 3 (I18b, 2026-10-08, §9 item 317) records the insulated lead (rod 403 + εᵣ 3 sleeve 404 on `GEO-36` step 5's 190 747-cell build): with the operator-derived complex contrast the adjoint fold equals the direct tip response to ≤ 1.9e-12, the N = 32 path fold to ≤ 2.9e-3, the Born fold misses; tip ratio 1.87–1.95 coated vs 2.46–2.61 bare, shaft/tip \|TF\| 3.04 vs 8.67 — a record on one mesh, the Done-when is not extended (`20261008T142909Z_WF-9-step3-W1.log`, `…143309Z_WF-9-step3-W2.log`). **Done-when:** *step 1 (I8)* — the lead in the coil-driven F-small phantom at 64 MHz, on one mesh with the lead's cells kept and switched by material (σ = 800 S/m against gel), so the adjoint fold is an exact discrete identity; asserted: for each of the four port drives the adjoint TF folded with the lead-free field equals the direct with-lead tip response at `DISCRETE_IDENTITY_RTOL`, the 1-D path fold ∫TF(z)·E_tan(z) dz reproduces it within **5 %** on every drive, and the lead-free Born fold misses ⇒ ✅ as a fixture. Not claimed: a measured transfer function, heating, a real lead. *History: `docs/planning/chunks/WF-9.md`.* | heavy |
| `WF-13` | The package installs and runs outside the repository | ✅ (§9 item 220; audit PASS): installed with no `PYTHONPATH`, the package reproduces `TH-6` bit-identically and `MAG-14` at 0.671 %; steps 2–16 (`Study.estimate()` and `PRICE_TABLE`, rows FH1 / FH2 / FC1 / FC2 / FC3 / FC-ROD1, U118 / U119, the estimate pinned on every shipped file by `WF-17` step 30) are records in the chunk file. *History: `docs/planning/chunks/WF-13.md`.* | standard by prediction |
| `WF-14` | **Every public entry point declares its validation envelope — Phase 9 (9a)** (opened 2026-10-03 by the weekly; §10 chain U2) *(history and §9 briefs: `docs/planning/chunks/WF-14.md`)* | ✅ 2026-10-04 (§9 item 221; audit PASS 2026-10-05): `tests/package/validation_envelope.toml` maps every public name to collected gate nodes, snapshot (S), registry (R), collection (C), no-number (N) and both planted controls green in one window (`20261005T022839Z_WF-14-W1.log`, 4 s); steps 2–3 (the generated `validation.md`, item 375; the F27 / F26 h-statement nodes, item 433, 196 nodes) are records in the chunk file. **U120** (item 784, 2026-11-03): `WF-13`'s U47 rehearsal re-run unedited on the current tree — installed `femem run` SAR rel 6.125e-11 to `MAT4_DEG1_RECORD` with no `PYTHONPATH`, surface 46 names = snapshot (+`ThermalBoxGeometry`, `TransverseQuadraticHeatSource`, `ColeColeMaterial` since 10-14), hash equal, dropped-name control differs (`20261103T150106Z_WF-14-U120-rehearsal.log`, 31 s). **Done-when (§4):** a registry mapping each public name to the gate test nodes that cover it; asserted in one logged window: every name in the declared public surface has an entry, every entry names a test node `pytest --collect-only` finds, and a snapshot test of the surface fails when a planted unregistered name is added (the control) ⇒ ✅. The registry holds test node IDs only — no band, record or tolerance is restated in it. | smoke by prediction |
| `WF-15` | `Study` v0: one scriptable study object, on the lossy sphere | ✅ Closed 2026-10-05 (step 1): `Study` v0 on the lossy sphere is bitwise against the direct call up to the solve, (G) 0.041 % against the closed form at degree 2; steps 2–17 (thermal, probes, boxes, Cole–Cole, the implant rod and its heating through the four-port `Study` — step 14, U127′ — and `MAT-12`'s transient and duty cycle on the EM sphere at degrees 1 and 2 — U131 / U132 / U135; step 19, U147, the tuned two-port `duty` refused with a named `TypeError`, `20261110T121255Z_WF-15-s19-W1.log`) are records in the chunk file. Open: U122 on `WF-17`. *History: `docs/planning/chunks/WF-15.md`.* | heavy (W2a 428 s, W2b 314 s; W1 115 s standard; label corrected 2026-10-05 18:00 audit) |
| `WF-16` | `Study` on the four-port birdcage: S and Touchstone through the API | ✅ 2026-10-05 (item 244): the `Study` reproduces the exact-rectangle four-port records at degree 1 and 2 with `PORT-11`'s gates green and the envelope registered; steps 2–19 are records beyond the Done-when, the newest `Study.load_result` on the F-human (step 15), the F-clinical 64×64 (step 16, `C_mode` 195.0265 pF with no solve), a multi-frequency file (step 17, item 722, `tune(at_hz=…)` — its refusal moved from `load_result` to `tune`, four-port only, ratified by the 11-02 weekly) and **U115** (step 18, item 780, 2026-11-03): `Study.solve_sweep(…, fields_at=f)` keeps the port fields at one swept frequency — on `mri:38`'s rod four-port \|B₁⁺\| and SAR10g equal the private `_solve_port_fields` route at rel 0 (`IDENTITY_BAND`), one factorisation per frequency as with `fields_at=None`, an off-sweep f refused naming `fields_at` (`20261103T133533Z_WF-16-s18-W1.log`) and **U117** (step 19, item 781, 2026-11-03): the 3 T twin's saved 64×64 loads and tunes with no solve, `C_mode` 43.9719 pF at rel 4.2e-15 to C15 (`REDUCTION_BAND`), the 1.5 T file refused by frequency, the YAML route equal by `float.hex` (`20261103T134456Z_WF-16-s19-W1.log`) and **U129** (step 20, item 867, 2026-11-07): the degree-2 1.5 T saved 64×64 loads and tunes with no solve, `C_mode` 185.4048 pF at rel 1.5e-14 to C18 (`REDUCTION_BAND`), the degree-1 file missing by 5.19e-02 (`20261107T182729Z_WF-16-s20-W1.log`). **Done-when:** the same `Study` object builds the F-small four-port fixture on exact port rectangles, κ off, solves at 64 MHz and returns S and a Touchstone file equal to `ANS-4` step 5's registered exact-rectangle records at 1e-9, with `PORT-11`'s reciprocity / passivity / C4 gates imported and green; envelope registered ⇒ ✅. *History: `docs/planning/chunks/WF-16.md`.* | heavy by prediction |
| `WF-17` | The first tutorial is a test | ✅ (§9 item 257, U7): the first-study tutorial is a test; every later guide page, the guide-page re-run batches A–F with S24 (steps 39–42, 44, 45), and steps 43 and 48–52 (U152′, U158–U162: the implant page's one-step section, the degree-2 determinism pair, the clinical page's rod-in-box and closed-form gel-box sections, the estimates page's measured 3 T and 1.5 T sections) are records in the chunk file. Step 53 (U163, §9 items 992 → 998, 2026-11-15, a record, no band): the estimates page's generated "Measured degree-2 clinical prices beside the estimator" section — `xxl` entry 62's two calls and window beside FC2, ratios estimator / measured 0.9709 GiB, 0.2653 / 0.6723 / 0.4388 seconds (window / call 1 / call 2), custody on entry 62's log through the `.gz`-safe reader; the 3 T (Sec) and 1.5 T (Order) premises edited as the 11-14 weekly licensed (`20261116T014216Z_WF-17-s53b-W0.log`, 32 passed; S30 (Page) builds 5 / hits 23, diff 0, `20261116T014250Z_WF-17-s53b-W1-s30.log`). Ops, item 997 (2026-11-15, a record, no band): the five landed plain-text log readers (`_custody`, `_parse_c52`, the closed-form page's readers, the lattice `_c39`) read through `tests/log_records.read_log_text(errors="strict")`, and `tests/package/test_wf17_log_gz_readers.py` asserts gz-vs-plain text identity on all ten logs they read (`20261116T013304Z_WF-17-i997-W1.log`). Step 54 (U166, §9 item 1001, 2026-11-15, provisional, a record, no band): the estimates page's generated "Measured 3 T degree-2 price beside the estimator" section, between the 3 T and degree-2 sections — `xl` entry 64 (`16` / `1 057 895` / `377.12` / `8158`, custody `[]`, tag 403 absent) beside FC4 (degree 2, not priced: `rows_used` `()`, ratios `nan` / `nan`, CLI = API by printed string) and FC3 (degree 1, same mesh and width, ratios 0.1044 / 0.0166); (c-plant) shows the FC4 ratio is the estimator's (`20261116T015521Z_WF-17-U166-W0.log`, 45 passed; S30 (Page) builds 5 / hits 23, diff 0, `20261116T015532Z_WF-17-U166-W1-s30.log`). Open: U116 (`mri:38`'s guide on U115's public field route) and U122 (the implant page's capacitor-sheet section, after U121). *History: `docs/planning/chunks/WF-17.md`.* | heavy (one degree-2 solve, 282.4 s measured on `WF-15`) |
| `WF-18` | **Rename to PrecessEM — Phase 9 (9c)** (operator directive 2026-10-08, interactive session; for the weekly to place in chain U) *(history and §9 briefs: `docs/planning/chunks/WF-18.md`)* | ⬜ Not started. Names: brand / repo **PrecessEM**, import **`precessem`**, distribution `precessem`, CLI **`precess`** (replacing `fem_em_solver` / `fem-em-solver` / `femem`). The old import name keeps working for one release through a shim that warns; `femem` stays as an alias of `precess` for the same release. **Done-when (§4):** one atomic commit; in one logged window set, the full CI selection that passes on `main` before the rename passes after it with the same pass / skip / xfail counts, the examples census reads the same `ok` count with `missing=0 broken=0`, `WF-14`'s envelope snapshot and collection gates are green under the new name, `WF-13`'s installed-package run reproduces `TH-6` and `MAG-14` from `pip install .` as `precessem`, and a test asserts `import fem_em_solver` still resolves to the same objects with a `DeprecationWarning`; no band, record or tolerance moves (a moved record is a defect) ⇒ ✅. Not claimed: renaming the Docker services / containers, the repo directory or the GitHub repository (operator), or rewriting history in the plan, journals and logs. | standard by prediction (CI selection + census) |

> `WF-2`/`WF-3` produce structurally valid manifests and reports containing
> physically meaningless numbers. The plumbing is fine and becomes useful the
> moment `TH-1` lands.
>
> *(No longer true of either: `WF-3` ✅ 2026-09-30; `WF-2` ✅ 2026-10-09, §9 item 346 — `femem run`'s bundle carries S and SAR reproduced at `L0_BAND`, provisional until the 10-10 weekly ratifies.)*


**`WF-9` — implant transfer-function excitation (ISO 10974 Tier 3)** ⬜
*(**feature ladder C1**, operator directive 2026-09-04; Phase 7, unscoped
until the phase opens.)* Implant heating is assessed from the tangential
`E` along the lead trajectory in the phantom *without* the implant, folded
with a transfer function measured by stepping a local source along the lead.
Needs: a point-like impressed current source (the `project_source` machinery
at a chosen location/direction), a path-integral post-processor for
`E_tan(z)`, and the fold. Gate: the Tier-3 identity — heating from
`∫ TF(z)·E_tan(z) dz` against the direct lead-in-phantom solve on one
fixture — plus published measured lead-heating data where a geometry
matches. `WF-7`'s hotspot chunk is to be scoped as this pipeline's consumer,
not as a local-SAR peak alone.

**`WF-6` narrative** — moved byte for byte to `docs/planning/chunks/WF-6.md` (`OPS-47`).

### EX — Examples (§5.4 ramp)

Standalone example chunks enqueued by the daily review when a chunk closes a
quantitative gate (§5.4). Each is sized for one implementer run, executes via
`./run_examples.sh`, produces combined-XDMF that opens in ParaView, and
demonstrates a **gated** capability from an angle no existing example covers.

**Ramp accounting.** Phases 1–3 are at quota per §5.4's
`min(5, gating chunks closed ✅)` (backfill `EX-4`…`EX-12` discharged the
2026-08-09 audit's shortfalls); Phases 4/5 accrue as their gates close.
`mri:1` is the one ungated example, labelled as such in the file.

| ID | Title | Status | Tier |
|---|---|---|---|
| `EX-1` | Two-torus port fixture: conforming mesh, cell and facet tags in ParaView | ✅ | standard |
| `EX-2` | Cylindrical phantom domain: wall classification and tags in ParaView | ✅ | standard |
| `EX-3` | Mass-averaged SAR on the standard-masses sphere: point and 1 g/10 g fields in ParaView | ✅ | standard |
| `EX-4` | Lossy plane wave: decay and phase vs closed form (first time-harmonic example) | ✅ | standard |
| `EX-5` | PEC cavity resonances: eigenfrequencies vs closed form, mode field in ParaView | ✅ | standard |
| `EX-6` | Sphere in a uniform field: solved quasi-static response vs closed form | ✅ | standard |
| `EX-7` | Waveguide/coax: the `TH-7` gated quantity as a runnable example | ✅ | standard |
| `EX-8` | Resonance guard on a frequency sweep: the `TH-1` step-5 detector firing | ✅ | standard |
| `EX-9` | Measured h-convergence rate as an example output (Phase 1) | ✅ 2026-08-10 | heavy (reclassified from standard: 130 s of solve, §5.1 names convergence studies) |
| `EX-10` | Gauge cross-check: penalty vs Lagrange-multiplier Coulomb gauge (Phase 1) | ✅ | standard |
| `EX-11` | Dodd–Deeds coil loading: ΔR vs closed form, eddy currents in ParaView | ✅ | standard |
| `EX-12` | Examples hygiene: stale claims, dead references, the 2026-02 PNG | ✅ | smoke |
| `EX-13` | `examples/mri/01` at the validated gauge floor: rank-spread measured, on-record numbers refreshed | 🚫 | standard |
| `EX-14` | Straight-wire VTX export repair + the refcheck freshness branch exercised | ✅ (2026-08-10: round-trip max\|B\| identical to 12 digits, rel diff 0.000e+00 vs 1e-10; freshness branch fired, then green) | standard |
| `EX-15` | Every runnable example gets a step-by-step analysis guide (3 steps, operator directive) | ✅ (2026-08-11: all three steps landed; **16 of 16** runnable examples checked against 3 required headings, `PENDING_GUIDES` empty, negative controls fired in all three steps) | standard |
| `EX-16` | `examples/mri/01`: converge the frequency-domain solve, then re-measure the rank spread | 🚫 (2026-08-10: solve converges — `preonly`/LU, `reason=4` — and the spread does **not** move, 23.5539% vs the 23.5545% unconverged record; anchor FAIL, negative-result clause taken. Fix landed; the 23% is the centerline sampling path, 3215× the phantom path on the same fields) | standard |
| `EX-17` | Circular-loop VTX export repair: port the `EX-14` diff, same round-trip anchor | ✅ (2026-08-10: round-trip max\|B\| 7.756122914931e-05 T both ways, rel diff 0.000e+00 vs 1e-10; loop's analytic numbers unmoved, checker green) | standard |
| `EX-18` | Gap-voltage port pair → Z → S on the two-torus fixture (the 3b-xvii/xviii gated capability; first ports example) | ✅ (2026-08-13: raw 0.894543 × ωM₁₂ printed as the miss it is, corrected 0.939849 (−6.02%) inside the unmoved 10%; ‖S−Sᵀ‖/‖S‖ = 2.5494e-05, ‖S‖₂ = 0.861449 ≤ 1; blind-ladder negative control −98.26% asserted to fail; 134 s at `-n 2`) | standard |
| `EX-19` | Larmor lossy-sphere example (`TH-10`'s newly gated capability: first example solving at 64/128 MHz; rubric in the §9 item) | ✅ (2026-08-13: `th:6`, fixture imported from the `TH-10` test module; all four records reproduced through the example path — 3.643% / 1.826% interior relL2 at 18.68× / 57.31× separation, power 3.629% with the quasi-static route missing 58.140% — max drift **1.7e-04** vs a pre-stated 1% band; convergence and both negative controls executed in-run; 24 s at `-n 2`) | standard |
| `EX-20` | Package S-parameter sweep example (`PORT-1` step 4's newly gated capability: first example calling `run_n_port_sparameter_sweep` on the solved field — the entry-point angle `EX-18` does not cover; full rubric in the §9 item, commissioned 2026-08-15 review) | ✅ (2026-08-16: `ports:2`, one `run_n_port_sparameter_sweep(..., gap_voltage_ports=specs)` call → two solves → Z → S; **all four step-4 records reproduced inside the pre-stated 1% band, misses 3.33e-07 / 3.23e-07 / 3.67e-06 / 2.29e-07** — raw 0.894543 printed first and asserted to *fail* the 10% band, corrected 0.939849 (−6.02%) inside it, ‖S−Sᵀ‖/‖S‖ = 2.5494e-05, ‖S‖₂ = 0.861449 ≤ 1, `\|Z₁₂−Z₂₁\|/\|Z₂₁\|` = 5.8309e-04 printed; negative control executed in-run — the deprecated heuristic route on the same mesh/ports gives an identically-zero off-diagonal, max\|ΔS\| = 3.078e-01 with its `DeprecationWarning` shown; 178.2 s at `-n 2`, 178 055 cells; guide pass 19/19 green. **Named limitation on record:** the sweep returns no fields, so the combined XDMF costs one extra port-1 solve (23.0 s) — surfacing `TimeHarmonicFields` from `SParameterSweepResult` is unscoped. *Audit 2026-08-16, 03:00 review: COMPLIANT; tier reclassified standard → heavy per the `EX-9` precedent — 178.2 s sits at the 180 s standard boundary and the wrap was 500 s; the companion docrefs log exits 1 on 24 pre-existing stale artifacts from other examples, none EX-20's (journaled in attempts.md)*) | heavy |
| `EX-21` | Graded birdcage conductor mesh (`GEO-15`'s newly gated capability: first birdcage example of any kind — geometry angle no example covers; mesh-only, no solve; full rubric in the §9 item, commissioned 2026-08-16 10:30 review) | ✅ (2026-08-16: `mesh:3`, `examples/meshing/03_birdcage_graded_conductors.py` + same-stem guide; two rungs of the same fixture on the **CAD (occ) mass** denominator — graded `h_c` = 1.6 mm keeps **0.967019** ≥ the imported `CAD_MASS_GATE` = 0.95, baseline global `setSize` keeps **0.740335** and is asserted to *fail* the same gate (`EX-18` inverted-assertion pattern), **separation 0.226685**; `GEO-9` box-partition identities re-asserted on **both** rungs < 1e-9 and the conductor CAD mass identical across them < 1e-12; 48 245 → 98 474 cells, 26.0 s at `-n 2`, standard. Every constant imported from the `GEO-15`/`GEO-9` test modules, none restated (`ANS-1`). Logs `20260816T200348Z_EX-21-example-n2.log` and `20260816T200516Z_EX-21-example-n2-final.log` (ratios bit-identical across both); docrefs `20260816T200505Z_EX-21-docrefs-fix.log` — 24 dead references, all pre-existing staleness from other examples, **none EX-21's** (its one own violation was found and fixed by the first docrefs run). **Measured note for `PORT-9` step 3:** the graded birdcage is 98 474 cells, confirming that entry's 98 k budget) | standard |
| `EX-22` | Restore the absent example artifacts: refresh runs for `mag` 01/02/04/05/06 + `mri:1` (commissioned 2026-08-16 weekly review — see entry below) | ✅ (2026-08-19: the standing **`stale=24` backlog is 0** — `dead=0 guide=0 stale=0 stale_severity=report exit=0`, the first `exit=0` the checker has returned under the `OPS-19` contract, guide pass green at 24/24 runnable examples and 33 guides scanned. Three refresh runs, all `-n 2`, real build for `mag` / complex for `mri`, all **exit 0**: `-e 1,2,4` 230 s, `-e 5,6` 151 s, `-e mri:1` 6 s. Every recorded anchor reproduced by the examples' own asserts — `EX-14`/`EX-17` VTX round-trips **0.000e+00** against their 1e-10 tol on both `straight_wire_B.bp` (max\|B\| 4.463816061893e-05 T) and `circular_loop_B.bp` (7.756122914931e-05 T); `EX-10` gauge cross-check **0.0003% probe / 0.0033% volume** against the 5% `MAG-15` gate; `EX-9` fitted h-convergence rate **1.1009** inside the `MAG-13` (0.7, 1.5) band, byte-matching the record; `EX-19`-era helmholtz centre B_z 1.28% with the h-ladder 0.89% / 0.24% / 1.28% at 70 054 / 103 984 / 160 478 cells. `mri:1` is the labelled **ungated** example — its printed record reproduced digit-for-digit against the guide, not gated: 9 261 cells / 2 077 vertices, cell tags 385 / 350 / 493 / 8 033, phantom \|E\| 1.244231e+02 / 3.150176e+02 / 1.975909e+02, \|B\| 8.791014e-08 / 2.771692e-06 / 1.292004e-06, ratio 1.529336e+08, coverage 493/493/493 with 0 drops. Logs `20260820T003126Z_EX-22-mag-124.log`, `20260820T003532Z_EX-22-mag-56.log`, `20260820T003812Z_EX-22-mri1.log`, `20260820T003833Z_EX-22-docrefs.log`. The 2026-08-16 premise correction held throughout: nothing was absent, `dead=0` before and after — this was a pure freshness restore) | heavy |
| `EX-23` | Two-torus port-sheet mesh (`GEO-16`'s newly gated capability: first example with an interior sheet surface, facet tags rebuilt dolfinx-side — geometry angle no example covers; mesh-only, no solve; commissioned 2026-08-17 review) | ✅ (2026-08-17: `mesh:4`, `examples/meshing/04_two_torus_port_sheet.py` + same-stem guide; both sheets **84 facets**, meshed/CAD = **1.000000000000** inside the imported `AREA_IDENTITY_BAND` = 1e-9, 211/212 area symmetry bit-identical (< 1e-12), out-of-plane spread **3.469e-18** m; kwarg-off control reproduces **79 534** cells with cell tags `{1,2,3,101,102}`, facet tags `{1,201,202}` and sheet tags asserted *absent* (`EX-18` inverted-assertion pattern); extents printed not gated — w = 1.200000000e-02 m, h = 7.977525299e-03 m, **w/h = 1.504225878** against the generator's CAD-side 1.504206917; port areas 1.563786482e-04 m² on both 201/202, unmoved. 79 888 cells / 13.7 s sheet mesh, 79 534 / 12.2 s control, **26.0 s** in-script (30 s harness) at `-n 2`, standard. Every constant imported from `tests/mesh/test_two_torus_port_sheet.py` and the `PORT-1` facet module (`ANS-1`) except `SHEET_SYMMETRY_BAND = 1e-12`, which the test holds only as an inline literal and the example restates unloosened (10:30-review audit note). Logs `20260817T140233Z_EX-23-list.log`, `20260817T140242Z_EX-23-example-n2.log`, docrefs `20260817T140416Z_EX-23-docrefs.log` — `dead=0 guide=0 stale=24 exit=2`, staleness-only and **none EX-23's** (its own artifacts are fresh; the 24 are `EX-22`'s standing backlog), guide pass green: 22 runnable examples checked, 31 guide files scanned (the 03:00-era "31 guides green" conflated the two counts; corrected 10:30 review)) | standard |
| `EX-24` | Lumped-sheet port at interior width (`PORT-9` step 2b's newly gated capability: first example instantiating the lumped-element port BC | ✅ Closed: `ports:3` runs the lumped-sheet width ladder and the half-width gate on one mesh against the imported, unmoved band, with the inverted negative control asserted. *History: `docs/planning/chunks/EX-24.md`.* | standard |
| `EX-25` | Degree-2 Larmor sphere: accuracy-per-cost side by side (`TH-12` step 1's newly gated capability: first example at any element order other than 1 | ✅ Closed: `th:7` solves the Larmor sphere at degree 1 and degree 2 on one mesh, reproducing the `TH-12` step-1 records inside the imported band. *History: `docs/planning/chunks/EX-25.md`.* | standard |
| `EX-26` | Poynting power-balance audit (`POST-5`'s newly gated capability: `poynting_power_balance` with the impressed-source term | ✅ Closed as written: `th:8` audits the Poynting power balance on the driven cylinder and the `TH-6` plane wave, with the two-term reading asserted to miss the band. *History: `docs/planning/chunks/EX-26.md`.* | standard (measured smoke) |
| `EX-27` | Region-resolution policy on the coil+phantom mesh (`GEO-17`'s newly gated capability: first example whose subject is a mesh-*sizing policy* | ✅ Closed as written on the first run: `mesh:5` shows the region-resolution policy against clamps-only on the coil and phantom mesh, with the `GEO-17` records reproduced and no band moved. *History: `docs/planning/chunks/EX-27.md`.* | standard (measured smoke) |
| `EX-28` | Gapped birdcage with leg terminals and port sheets (`GEO-18`'s newly gated capability, both steps: first example with a discontinuous conductor | ✅ Closed as written on the first run: `mesh:6` shows the gapped birdcage with leg terminals and port sheets, meshed against analytic, with no band moved. *History: `docs/planning/chunks/EX-28.md`.* | standard |
| `EX-29` | Doc-reference checker freshness-gates every example's own `paraview_output/` (22 of 27 examples were never checked — known-issues 2026-08-23; commissioned 2026-08-23 weekly review) | ✅ | smoke |
| `EX-30` | Refresh the 13-example stale artifact set the checker could not see (10–17 d old on 2026-08-23; commissioned 2026-08-23 weekly review; re-scoped 2026-08-24 … | ✅ Closed: all four legs of the stale-artifact refresh landed, with the licensed in-class example re-records version-tagged. *History: `docs/planning/chunks/EX-30.md`.* | heavy (measured standard per leg) |
| `EX-31` | Ring-gapped birdcage with dual port families (`GEO-20` step 1's newly gated capability: first example with ring-gap terminals as exact disks from radial cut planes, and the first 12-port dual-family mesh — a geometry angle `EX-28` (leg gaps only) does not cover; mesh-only, no solve; commissioned 2026-08-24 10:30 review) | ✅ (2026-08-24: `mesh:7`, `examples/meshing/07_birdcage_ring_gap_ports.py` + same-stem guide; **closed as written on the first run**, every element of the rubric executed and no band moved — see the prose entry for the digits. *Audited COMPLIANT 2026-08-24 18:00 review: all three footers verified (75 / 72 / 1 s, statuses 0/0/2 with exit 2 = the `OPS-19` stale-only contract), the printed assertions match this entry to the last decimal, the gate-module strengthening is real (records asserted at their source, lines 403/458/538, imported by the example not restated), and `git show 7529fa4 -- tests/` is 24 insertions / 0 deletions — purely additive, nothing loosened*) | standard (measured standard) |
| `EX-32` | Birdcage 4-port power-wave S-matrix at 10 MHz (`PORT-9`'s newly gated capability: first example solving ports on the **birdcage** — every existing S-parameter example is two-torus (`EX-20`/`EX-24`/ports:1–3), and `EX-28`/`EX-31` are mesh-only; commissioned 2026-08-25 10:30 review, §5.4 ramp) | ✅ *(2026-08-26, green on the first run; every gate-module record reproduced exactly and the only reading that moved is the one (d3c) declares non-reproducible)* | standard (measured standard, 88 s) |
| `EX-33` | 16-leg gapped + sheeted birdcage mesh (`GEO-19`'s newly gated capability: first example above four legs — `EX-28`/`EX-31`/`mesh:3` are all 4-leg, so the geometry angle is new; mesh-only, no solve, no port claim; commissioned 2026-08-25 18:00 review, §5.4 ramp) | ✅ *(audited COMPLIANT 2026-08-26 18:00 review: three footers Status 0 at 131 / 126 / 1 s, all 16 claimed digits grep out of the run log, bands *and* `_assert_identity_family` imported from the gate module, the gate-module diff a single +6-line additive hunk)* | standard |
| `EX-34` | Birdcage 4-port S-matrix across the frequency ladder 10 / 64 / 128 MHz on **one** mesh (`PORT-11`'s newly gated capability: first example solving ports at a Larmor frequency — `EX-32` is 10 MHz only, `EX-19`/`EX-25` are Larmor on the sphere with no ports; the drive/output angle is the frequency ladder itself: loss tangent, cells/δ, cells/λ, `\|Im P\|/Re P` and the three gate readings side by side; commissioned 2026-08-26 18:00 review, §5.4 ramp) | ✅ 2026-08-28 (`ports:5`, `20260828T110615Z_EX-34-run2.log`, **139 s** Status 0 at `-n 2` complex; one 116 085-cell mesh, 12 driven solves, all three gates green on all three rungs; 128 MHz cells/λ 12.5024 ≥ 10) | standard |
| `EX-36` | Re-run the examples whose artifacts the 2026-08-28 rename (`67e4c1c`) orphaned | ✅ Closed: the examples the rename orphaned are re-run and the corpus census came back clean. The stale `ports:3` narrative ladder it found went to a review. *History: `docs/planning/chunks/EX-36.md`.* | standard per leg (105 / ≈ 500 / 447 / 935 s measured by `EX-30`) |
| `EX-35` | 16-leg **ring-gapped** birdcage mesh — the 32-ring-port high-pass layout (`GEO-20` step 2's newly gated capability, 2026-08-29: `EX-31`/`mesh:7` cuts the rings at 4 legs, `EX-33`/`mesh:8` cuts the **legs** at 16 — no example has the production high-pass topology at the production leg count; mesh-only, no solve, no port claim; commissioned 2026-08-29 10:30 review, §5.4 ramp) | ✅ | standard |
| `EX-37` | **`ANS-1`/`ANS-3` are unrunnable since the 2026-08-28 rename**: `67e4c1c` rewrote the two cases' `__import__("01_dodd_deeds_coil_loading")` / `("02_package_sparameter_sweep")` strings along with their artifact basenames, and no module of the prefixed name exists — the `ANS-4` slot's observation, verified by the 18:00 review from `git log -S`; restore the two strings, re-run `ans:1` and `ans:3`, census (commissioned 2026-08-30 18:00 review; known-issues entry) | ✅ *(2026-08-31, 19:30 slot — negative control Status 1/3 s, both cases green, 63 s / 128 s)* | standard (≈ 5 + 70 + 131 s + census) |
| `EX-38` | `\|B₁⁺\|` map of the loaded 4-leg birdcage at 10 MHz into ParaView (`WF-6` step 1's newly gated capability, 2026-08-30: the CG1-projected `\|B₁⁺\|` field under gates (i)/(ii) — no example writes a B₁⁺ field, `ports:4`/`ports:5` stop at S-matrices; the output-quantity angle; commissioned 2026-08-30 18:00 review, §5.4 ramp) | ✅ *(2026-08-31, 15:00 slot — `ports:6` green in 63 s: gate (i) 9.795751117e-03 vs record to 1.195e-08, CG1 C4 covariance 2.1870% vs the 5% band and the record to 1.643e-05, DG0 control 8.6516% to 3.227e-06, 51/51 valid, cell ratio 1.000000)* | standard (≈ 70 s) |
| `EX-39` | Quadrature-driven birdcage: co- and counter-rotating `\|B₁⁺\|` / `\|B₁⁻\|` maps at 10 MHz into ParaView (`WF-6` step 2's newly gated capability, 2026-08-31: exact superposition of the four single-drive fields with `e^{∓jkπ/2}` phases, C4-invariance and the mirror identity at the CG1 floor — the **drive** angle: no example superposes ports or shows the two rotating senses side by side; `EX-38` is the single-drive map; commissioned 2026-08-31 03:00 review, §5.4 ramp) | ✅ *(2026-08-31, 16:30 slot — `ports:7` green in 81 s: identity (a) 0.9818% and (b) 0.8087% vs the 5% band and step 2's records to 9.619e-06 / 3.585e-05, mis-paired control 95.1975% at 118× identity (b), gate (i) 9.795751117e-03 to 1.195e-08, 51/51 valid, cell ratio 1.000000)* | standard (≈ 100 s) |
| `EX-41` | **`mesh:6` and `mesh:7` get a footered run on the 0.11 image** — their newest log (`20260901T050408Z_EX-36-leg-mesh-b.log`) is the orphaned footerless window the 00:00 slot's executor left running; both printed "All identities hold" inside it (58.6 s / 94.8 s) but §4 wants a footer and an elapsed time, and the last footered run is `20260825T213323Z_EX-30-mesh-run-6to7.log`. One `example-runner` window, `-e mesh:6,mesh:7 -t 400`, anchors = each example's own identity asserts + `mesh:7`'s 12-port dual-family record reproduced to the digit; census `exit != 1`. Standard. Opened 2026-09-02 weekly review (examples-health pass) | ✅ *(2026-09-02, 15:00 slot — `mesh:6` footered 50.6 s (harness 54 s), `mesh:7` footered 83.5 s (harness 85 s); both "All identities hold"; `mesh:6` cells=116085 gate 0.970069/0.95, uncut control cells=98666 vs record 98474 ratio 1.001950; `mesh:7` GEO-20 ring-gapped rung cells=110786, leg+ring rung cells=128111, uncut control cells=98666 ratio 1.001950 — every digit matches `20260825T213323Z_EX-30-mesh-run-6to7.log`; census dead=0/exit=2 pre and post (stale 16→17, one item crossed the 48h threshold mid-run, unrelated to this chunk))* *Audited 2026-09-02 18:00 review — `auditor` PASS on all eight checks; the review re-traced `20260902T200352Z_EX-41-mesh7.log:6949` against `20260825T213323Z_EX-30-mesh-run-6to7.log:15030` itself (identical to every printed digit). The slot's two-windows-not-one deviation is ratified: two independent footers are stronger evidence at no extra cost, and multi-example items now say "one window per example".* | standard |
| `EX-42` | `mat:1` prints the finite-wire-corrected Dodd–Deeds beside the filament form and the FEM ΔR | ✅ Closed: `mat:1` prints the finite-wire-corrected Dodd–Deeds beside the filament form and the FEM ΔR, with the `MAT-8` record reproduced and no band moved. *History: `docs/planning/chunks/EX-42.md`.* | standard (≈ 60 s; measured 64 s) |
| `EX-43` | The first coil-driven SAR quantity in ParaView | ✅ Closed: the first coil-driven SAR quantity is in ParaView through `ports:9`, with every anchor met on the first run. *History: `docs/planning/chunks/EX-43.md`.* | standard (≈ 130 s; measured 77 s) |
| `EX-44` | The longitudinal ring-gap port sheet in ParaView | ✅ Closed: `mesh:10` shows the longitudinal ring-gap port sheet in ParaView, with both cell records and the closed-form halves reproduced and no solve. *History: `docs/planning/chunks/EX-44.md`.* | standard (≈ 70 s; measured 59 s) |
| `EX-45` | The 16-leg longitudinal ring-gap rung in ParaView, with its two-state terminal triangulation as a per-port cell field (`GEO-26` step 3's newly gated … | ✅ Closed: the 16-leg longitudinal ring-gap rung is in ParaView with its two-state terminal triangulation as a per-port cell field, and the census is clean after the full-filename fix. *History: `docs/planning/chunks/EX-45.md`.* | standard (82 s measured) |
| `EX-46` | The first field on the 32-ring-port birdcage in ParaView | ✅ Closed: the first field on the 32-ring-port birdcage is in ParaView, with one additive return key licensed and its gate module re-run green. *History: `docs/planning/chunks/EX-46.md`.* | heavy by ceiling; **measured 103 s** (example) + 143 s (gate re-run) |
| `EX-47` | The ring rung's mirror pair in ParaView | ✅ Closed: the ring rung's mirror pair is in ParaView as an example of two columns, not a gate, with no band moved. *History: `docs/planning/chunks/EX-47.md`.* | standard (host-runner, ≤ 600 s window) |
| `EX-48` | The ring rung's quarter turn in ParaView | ✅ Closed: the ring rung's quarter turn is in ParaView as two columns and one field-level reading, not a gate, with no record re-recorded. *History: `docs/planning/chunks/EX-48.md`.* | standard (host-runner, ≤ 600 s window; **measured 112 s**) |
| `EX-49` | An asymmetric drive through `ports.superpose_drives` in ParaView | ✅ Closed and audited PASS: an asymmetric drive through `ports.superpose_drives` is in ParaView. Its ruling gave §9 standing rule (e), the asserted-or-predicted label on every negative-control factor. *History: `docs/planning/chunks/EX-49.md`.* | standard (host-runner, ≤ 600 s window; **measured 73 s**) |
| `EX-50` | The conductor as a hole in ParaView | ✅ Closed: the conductor as a hole is in ParaView, with the census count matching its prediction and no `src/` or `tests/` change. *History: `docs/planning/chunks/EX-50.md`.* | standard (≈ 30 s at `-n 2`; **measured 4 s**) |
| `EX-51` | The F-human birdcage in ParaView | ✅ Closed: the F-human birdcage is in ParaView, with one additive keyword on the gate helper re-run green and no Phase 6 cost claim. **Re-recorded 2026-10-11 (`GEO-35` 2b-ii part 2):** `mesh:12` asserts the `…_UNBOXED` records, 500 425 / 203 945 exactly, all gates green, 230 s at `-n 2`, setup PNG redrawn without leg boxes (`20261012T022429Z_GEO-35-s2bii-W4-mesh12.log:17761,17767`). *History: `docs/planning/chunks/EX-51.md`.* | standard (host-runner, ≤ 600 s window; measured **196 s** at `-n 2`) |
| `EX-40` | `\|B₁⁺\|` maps at 64 and 128 MHz — the Larmor frequency ladder in ParaView (`WF-6` step 2b's newly gated capability, 2026-08-31: the five identities hold at 64/128 MHz on one mesh at 21.89 / 12.50 phantom cells/λ — the **frequency** angle: `EX-38`/`EX-39` are 10 MHz only, `EX-34` runs the ladder but stops at S-matrices; commissioned 2026-08-31 10:30 review, §5.4 ramp; **queued 18:00 review, §9 item 2**) | ✅ *(2026-08-31 21:00 slot, `20260901T020415Z_EX-40.log`, Status 0, **113 s** at `-n 2` — every anchor met on the first run: gate (i) 9.5231e-03 / 9.2445e-03 reproducing step 2b to 4.4e-08 / 5.4e-08; gate (ii) 2.2187% / 2.1315% inside the 5% band, records to 1.1e-05 / 1.8e-06; the mis-rotated P3@+90° control 24.7535% / 25.2589% outside it, 11.2× / 11.9× separation; cells/λ 21.8936 / 12.5024 above the imported floor of 10; 51/51 valid; one mesh, ratio 1.000000, `reused_mesh`)* | standard (≈ 130 s; measured 113 s) |
| `EX-52` | Imposed-field SAR on the lossy sphere against its closed form in ParaView | ✅ Closed: imposed-field SAR on the lossy sphere is shown against its closed form in ParaView, on an imposed uniform field with no coil claim. *History: `docs/planning/chunks/EX-52.md`.* | standard (host-runner, ≤ 500 s window) |
| `EX-53` | Mass-averaged 10 g SAR on the coil-driven field in ParaView | ✅ Closed: mass-averaged 10 g SAR on the coil-driven field is in ParaView, importing the gate module's fixture rather than re-implementing it. *History: `docs/planning/chunks/EX-53.md`.* | standard *(host-runner window ≤ 600 s; measured 91 s at `-n 4`)* |
| `EX-54` | **The birdcage as a PEC hole in ParaView — the 4×4 and the field on the conductor-free mesh** (`examples/ports/14_birdcage_pec_hole_ports.py` + same-stem guide; `example-runner`) | ✅ **DONE 2026-09-14 02:31Z.** `_hole_rung(FREQUENCY_HZ)` at 10 MHz, imported unmodified except five additive return keys (`mesh`, `cell_tags`, `wall_facet_tags`, `sheet_facet_tags`, `fields` — nothing renamed/removed; gate module re-run green in-slot, 6 passed / 4 skipped (solid control not requested), 32.70 s, `20260914T023104Z_EX-54-th15-rerun.log:961`). All three imported gates plus the power identity asserted and green in the example: reciprocity 1.508e-14 (band 1e-3), σ_max 0.999994234 (band 1 + 1e-9), class spreads self 0.0190 %/adjacent 0.0094 %/opposite 0.0059 % (band 0.5 %), power rel dev 1.327e-10 (band 1e-3) — `20260914T023145Z_EX-54.log:897`. Printed only (no pre-registered band): the hole's 4×4 beside `PORT-11`'s solid σ = 800 S/m 10 MHz record (`LEG_D_S_MATRIX_10MHZ`), `max\|ΔS\|` per C4 class self 9.769e-02 / adjacent 2.484e-02 / opposite 4.828e-02 (`:908`); the cavity wall's (tag 401) `\|n × E\|` reads 2.989e-16 over 4.052772e-02 m² — machine-zero, confirming the PEC enforcement (`:910`). 80 181 cells, 31.8 s total at `-n 2` (`:895,914`). Combined XDMF `examples/ports/paraview_output/ports_14_birdcage_pec_hole_ports_combined.xdmf` (`CellTags`, `E_magnitude` DG0, `facet_tags` grid with tag 401). Census: pre `dead=0 guide=0 stale=0 exit=0` (`20260914T022309Z_EX-54.log:39`), post `dead=0 guide=0 stale=0 exit=0` (`20260914T023230Z_EX-54-census-post.log:39`). No copper, Larmor-accuracy or `TH-14` claim. | standard (host-runner, measured 34 s example / 32.70 s gate rerun) |
| `EX-55` | The 32-port ccw quadrature drive on the 16-leg birdcage in ParaView (`examples/ports/15_birdcage_sixteen_leg_quadrature_b1.py` + guide; `example-runner`) | ✅ Closed 2026-09-14 and audited PASS: the 32-port ccw quadrature drive on the 16-leg birdcage asserts the C16 spread, the mirror identity and the exact power identity at their imported bands, with the census and the rule (a) gate re-run green. Nothing is open; no homogeneity, absolute or Larmor claim is made. **Done-when:** the three imported identities asserted, artifact named, census `exit != 1`, elapsed recorded. *History: `docs/planning/chunks/EX-55.md`.* | heavy (host-runner; measured 168 s at `-n 2`) |
| `EX-56` | **`\|B₁⁺\|` spread against resolution — the two-rung ladder in ParaView** (`examples/ports/16_birdcage_b1_resolution_ladder.py` + guide; `example-runner`) | ✅ **2026-09-14 06:00 slot:** fall asserted 5.2506 % → 2.0719 % at `-n 2` (records printed, rel 5.7e-06 / 6.9e-06, width disclosed), flagged 142 s / unflagged 136 s (`20260914T110437Z_EX-56.log:3618,3636`, `…110730Z_EX-56-control.log:3618,3629`), census 51/3/48/0 as predicted; **deviation:** two per-rung combined XDMFs, not one time-stepped file (the rungs' meshes differ). ~~⬜~~ *Opened 2026-09-13 18:00 review (step 5) for the gate `WF-6` step 5 closed 2026-09-13: `ports:8` (`EX-40`) ladders* frequency*; no example ladders resolution.* Angle: output quantity — the worst-radius C4 four-copy spread of CG1 `\|B₁⁺\|` on the unloaded F-small birdcage at 10 MHz at global sizing ×1 and ×0.012, imported from `tests/validation/test_birdcage_b1_plus_closed_form.py` (`LADDER`, `STEP5_RECORDED_SPREADS` 5.2506 % / 2.0719 % at `STEP5_SPREAD_RTOL` 1e-3 — imported records, asserted only as the module asserts them; the monotone-fall identity asserted), both rungs' `\|B₁⁺\|` written to one combined XDMF with the rung as a time step so ParaView steps through resolution, and the interior CV printed beside step 4a's filament closed-form CV (printed only). Trap: the ×0.0095 rung stays off (its power residual is a banked known-issues negative); `-n 4` is the record width and the runner's `-n 2` carries the 1e-4-class MUMPS drift (known-issues 09-09) — the example asserts the records at the module's rtol 1e-3 *only if* it runs at `-n 4`; at `-n 2` it prints them and asserts the fall. Done-when (§4): the fall asserted, records reproduced or their width disclosed, artifact named, census `exit != 1`, elapsed recorded. Predicted 204 s at `-n 4` (`20260909T200431Z_WF-6.log`); `-n 2` *predicted* ≈ 6–7 min, runner `-t 900`. A convergence statement only — no closed-form, homogeneity or Larmor claim. | heavy (host-runner) |
| `EX-57` | One setup figure per example guide | ✅ Closed 2026-09-28: the setup-figure census exits 0 on `main` with every example guide carrying a setup figure and docrefs clean. The re-render of the earlier figures with clipped legends stays open for a review to queue. **Done-when:** closes when the census exits 0 on `main` *History: `docs/planning/chunks/EX-57.md`.* | the example's own's own |
| `EX-58` | **The tuned birdcage — sweep, `C_tuned`, the in-model field** (`examples/ports/17_birdcage_tuned_circuit.py` + guide; `example-runner`) | ✅ **2026-09-14 06:00 slot (take-next):** `|Im Z_in|/|Z_in|` 2.113e-15 ≤ `TUNING_IM_Z_RTOL`, tuned `S₁₁` residual 8.255812e-05 ≤ `REDUCTION_BAND` (`20260914T112025Z_EX-58.log:43,1799`, flagged 58 s; control `…112312Z_EX-58-control.log:43,1799`, 57 s, same digits); rule (a) lift in `test_port_circuit_layer_field.py` (additive), gate re-run 7 passed 169 s (`…111515Z_EX-58-rule-a-gate.log:3753`); census 52/4/48/0; XDMF has two time steps, written with `XDMFFile` directly. ~~⬜~~ *Opened 2026-09-14 03:00 review (step 5) for the gates `PORT-15` steps 2–3 (closed 2026-09-14) and `PORT-14` step 3 (closed 2026-09-13): `ports:3` (`EX-24`) shows lumped* sheets *and nothing in the corpus terminates a port in a capacitor or tunes.* Angle: drive — `S_64MHZ_EPS0_RECORD`, `tuning_sweep`, `select_c_tuned`, `tuned_input` imported from `tests/validation/test_port_circuit_layer_field.py` (never copied); the `Im Z_in(C)` sweep printed, its zero asserted at the imported `TUNING_IM_Z_RTOL`; one in-model 64 MHz solve with `PORT-14`'s κ-corrected capacitor sheets at `C_tuned` (the module's terminated-network helper made importable additively under §9 rule (a), gate module re-run green in the same slot), the tuned `S₁₁` residual asserted at the imported `REDUCTION_BAND`; `\|E\|` on the phantom at `C_tuned` beside the 50 Ω baseline as two time steps of one combined XDMF. Done-when (§4): the two imported assertions green, artifact named, census `exit != 1`, setup figure per `EX-57`, elapsed recorded. Predicted ≈ 90 s at `-n 2` (`20260914T020419Z_PORT-15.log:3728`); runner `-t 600`. A series resonance on one fixture — no match, no mode frequency, no Larmor-accuracy claim. | standard (host-runner) |
| `EX-59` | **The copper birdcage — the surface loss density on the coil** (`examples/ports/18_birdcage_copper_leontovich.py` + guide; `example-runner`) | ✅ **Done 2026-09-14 07:30 slot** — imported port gates and the surface-loss identity asserted green, and the gate's 10 MHz copper digits reproduced (`20260914T124230Z_EX-59.log:904–918` against `20260914T021500Z_TH-14.log:1032–1055`; flagged 35 s, control 32 s). The facet field is **DG0 per adjacent cell, in W, not a density**; its owned sum equals `P_surf` to 1.5e-15, asserted and true by construction. The helpers are imported as they are, with no gate edit. Census 52/4/48/0 → 53/5/48/0, docrefs `exit=0`. ⬜ *Opened 2026-09-14 03:00 review (step 5) for the gated `TH-14` birdcage step (2026-09-14; the row is 🟡 on its slab step, the birdcage identities themselves are gated): `ports:14` (`EX-54`) is the PEC hole; no example carries a surface impedance.* Angle: output quantity — the 10 MHz copper configuration of `tests/validation/test_th14_birdcage_copper.py` (helpers module-private, §9 rule (a)); the three imported port gates and the surface-loss identity at `DISCRETE_IDENTITY_RTOL` asserted; the DG0 facet field `½Re(1/Z_s)\|n × E\|²` on tag 401 written through the facet grid (the `ports:14` pattern) beside `\|E\|` on the phantom; `P_coil/P_in` printed beside the σ = 800 solid's share. Done-when (§4): the imported assertions green, artifact named, census `exit != 1`, setup figure per `EX-57`, elapsed recorded. Predicted ≈ 60 s at `-n 2` (the gate fixture: 104.7 s for 3 f × 3 σ under reuse, `20260914T021500Z_TH-14.log:1348`); runner `-t 600`. No absolute copper-coil claim (`ANS-6`). | standard (host-runner) |
| `EX-60` | **The lossy-wall cavity — Q against σ beside Pozar** (`examples/time_harmonic/10_lossy_wall_cavity_q.py` + guide; `example-runner`) | ✅ **Done 2026-09-14 07:30 slot (take-next)** — both imported assertions green: `Q/Q_c − 1` +0.010 % ≤ 5 % and `Q(1e6)/Q(1e4)` 9.995026 ≤ 5 % from 10. The PEC control `|Im λ|/Re λ` 4.798e-19 and the exported mode's Rayleigh quotient (rel 7.3e-14) are also asserted, and the rung lines match `20260914T004807Z_TH-14.log:96–99` (`20260914T124948Z_EX-60.log:39–59`; flagged 11 s, control 7 s). `src/`: additive `return_mode` / `return_vectors` in `core/cavity.py`, with the gate re-run green (18 passed, 40 s). Census 53/5/48/0 → 54/6/48/0, docrefs `exit=0`. ⬜ *Opened 2026-09-14 03:00 review (step 5) for the gate `TH-14` step 1 (closed 2026-09-13): `th:2` (`EX-5`) is the PEC box; no example solves a lossy-wall eigenmode.* Angle: boundary model — `solve_impedance_wall_cavity_mode` (`core/cavity.py`) on the `TH-9` box at σ ∈ {1e4, 1e6, 5.8e7}, Pozar's `Q_c` imported from `tests/validation/test_cavity_leontovich_q.py` (never restated), the 1e4 rung asserted at the imported 5 % band and the `√σ` scaling identity asserted, copper printed as the PEC-limit reading; the TE₁₀₁ `\|E\|` to XDMF. Done-when (§4): the two imported assertions green, artifact named, census `exit != 1`, setup figure per `EX-57`, elapsed recorded. Predicted ≈ 1 min at `-n 2` (9.2 s of eigen-solves in a 45 s window, `20260914T004807Z_TH-14.log`); runner `-t 300`. One box, one mode; no coil. | standard (host-runner) |
| `EX-61` | The tuned birdcage's resonance curve | ✅ Done 2026-09-21 and audited PASS: `ports:19` shows the tuned birdcage's resonance curve with all three imported anchors green, f0, R_in and loaded Q printed but never gated, and the rule (a) lift re-gated green. Nothing is open; no absolute `S₁₁` or `Z_in` claim is made. **Done-when:** the three imported assertions green through `./run_examples.sh -e ports:19`, artifacts named in the guide by full filename, both censuses before / after, docrefs `exit != 1`, elapsed recorded. *History: `docs/planning/chunks/EX-61.md`.* | heavy |
| `EX-62` | The 10 g SAR hotspot at human scale | ✅ Done 2026-09-26 and audited PASS: `mri:4` reproduces the `WF-7` human-scale 10 g SAR hotspot gate (coverage identity, C16 and mirror spreads, contained averaging balls, negative control) at `-n 8`, with the rule (a) lift re-run green and censuses matched. Nothing is open; no absolute SAR or compliance claim is made. *Re-recorded 2026-10-11 (`GEO-35` 2b-ii part 1): leg port boxes dropped, 532 935 cells, SAR10g(c\*) 2.079992000e-04, worst C16 0.3373 %, mirror 0.0234 %, `20261012T021003Z_GEO-35-s2bii-W2-mri4.log`.* *History: `docs/planning/chunks/EX-62.md`.* | heavy (`-n 8`, record width) |
| `EX-64` | PEC wall against PMC wall | ✅ Done 2026-09-30 and audited PASS: `th:11` solves the PEC and `NATURAL` (PMC) walls on one box, every imported anchor and control green on the gate's digits, rule (a) re-gate green and censuses clean. Nothing is open. **Done-when:** asserts green through `./run_examples.sh -e th:11 -n 2`, artifacts named by full filename, both censuses before / after (`examples` +1), docrefs `exit != 1`, elapsed recorded. *History: `docs/planning/chunks/EX-64.md`.* | standard (host-runner, ≤ 300 s window) |
| `EX-65` | Which gel? The three gelled-saline presets on the lossy sphere, against the closed form (`examples/materials/03_gelled_saline_presets.py` + same-stem guide, … | ✅ Done 2026-09-30 and audited PASS: `mat:3` runs the three gelled-saline presets on the lossy sphere, SAR and Im/Re inside `MAT4_BAND` against the closed form and every cross-preset control missing, with the re-gate and censuses green. Nothing is open. **Done-when:** asserts green through `./run_examples.sh -e mat:3 -n 2`, artifacts named by full filename, both censuses before / after (`examples` +1), docrefs `exit != 1`, elapsed recorded. *History: `docs/planning/chunks/EX-65.md`.* | standard (host-runner, ≤ 300 s window) |
| `EX-66` | Quarter the box | ✅ Done 2026-09-30 and audited PASS: `th:12` cuts the driven box on PMC and PEC symmetry planes and reproduces the full model's port records within the imported band, with the pinned-cut control missing, the gate re-run green and censuses clean. Nothing is open. **Done-when:** asserts green through `./run_examples.sh -e th:12 -n 2`, artifacts named by full filename, both censuses before / after (`examples` +1), docrefs `exit != 1`, elapsed recorded. *History: `docs/planning/chunks/EX-66.md`.* | standard (host-runner, ≤ 180 s window) |
| `EX-67` | SAR heats the gel ball | ✅ Done 2026-09-30 and audited PASS: `mat:4` heats the gel ball from SAR, with the steady, transient and calorimetric-slope readings equal to the `MAT-10` and `MAT-12` records under their imported bands and both controls missing. Nothing is open; the demonstrated rows' lane-3 opening was ratified by the 2026-09-30 weekly. **Done-when:** asserts green through `./run_examples.sh -e mat:4 -n 2`, artifacts named by full filename, both censuses before / after (`examples` +1), docrefs `exit != 1`, elapsed recorded. *History: `docs/planning/chunks/EX-67.md`.* | standard (host-runner, ≤ 180 s window) |
| `EX-68` | The gel warms and conducts better | ✅ Done 2026-09-30 and audited PASS: `mat:5` runs its own σ(T) fixed-point loop over the `MAT-5` gate helpers, all eight iterates equal to the gate record and the frozen-σ control missing, with no gate edit. Nothing is open. **Done-when:** asserts green through `./run_examples.sh -e mat:5 -n 2 -t 600`, artifacts named by full filename, both censuses before / after (`examples` +1), docrefs `exit != 1`, elapsed recorded. *History: `docs/planning/chunks/EX-68.md`.* | heavy (host-runner, ≤ 600 s window) |
| `EX-69` | Perfusion carries half the heat | ✅ Done 2026-09-30: `mat:6` runs green with the steady and transient Pennes perfusion readings on the Robin ball inside their imported bands, the w = 0 controls missing, and census and docrefs clean. Nothing is open. **Done-when:** steady (i) `<= INTERIOR_BAND`, (iii) `<= BALANCE_BAND`, transient (P_w) `<= INTERIOR_BAND`, (E) `<= BALANCE_BAND`, all imported; controls w = 0 steady / transient each miss `INTERIOR_BAND`; asserts green through `./run_examples.sh -e mat:6 -n 2`, artifacts by full filename, both censuses before / after (`examples` +1), docrefs `exit != 1`, elapsed recorded. *History: `docs/planning/chunks/EX-69.md`.* | standard |
| `EX-72` | *(history and §9 briefs: `docs/planning/chunks/EX-72.md`)* **A metal sphere in the gel: the dipole it adds and the SAR it concentrates** (`imp:1`, a new `implants` example group; `example-runner`; §5.4 ramp for Phase 7 — opened 2026-09-30 by the weekly, serial on `IMP-1` ✅). Imports `IMP-1`'s gate records, never restates them; a ParaView file of the SAR enhancement around the sphere; the natural-cavity control. **Done-when (written 2026-10-01 03:00 review, rule (j)):** `examples/implants/01_pec_sphere_in_gel.py` + same-stem guide with `## Setup figure`, on a new `imp:` runner group (a disclosed `scripts/run_examples.sh` change, ratified by this review on the `ports:` / `th:` / `mat:` precedents `89fd522c`, `cdbc637c`, `486ddd04`); at `IMP-1`'s a/6 rung, degree 2: \|β − 1\| ≤ `TH15_BETA_BAND`, Φ miss ≤ `MAT4_BAND`, cavity max \|E\| < 1e-12, each printed beside the gate record; the natural-cavity control misses both bands; combined XDMF with SAR/SAR₀; census +1 ⇒ ✅ (§9 item 138). **✅ 2026-10-01 (12:00 slot):** β 0.985596, Φ_h / Φ_cf,h 1.606456 / 1.605264 (0.0743 %), cavity 5 670 dofs max 0 — every gate record digit for digit; control misses both (\|Δβ\| 1.4328, Φ 44.2975 %); census 66 → 67 `missing=0 broken=0`; 44 s at `-n 2` (`20261001T170201Z_EX-72.log`). *Audited PASS 2026-10-01 18:00 review (`…170201Z_EX-72.log:101–108,110,113,116–117`; control `…170250Z_EX-72-control.log:94–101,109–110`; records `20261001T034025Z_IMP-1-land-W3.log:372–379,383–384`).* | ✅ | standard |
| `EX-73` | *(history and §9 briefs: `docs/planning/chunks/EX-73.md`)* **SAR heats a perfused ball: the rise with the sink switched on** (`mat:7`; `example-runner`; §5.4 ramp for Phase 8 — opened 2026-09-30 by the weekly on `MAT-14` ✅). Imports `MAT-14`'s records; the w → 0 control beside it. | ✅ **Audit PASS 2026-10-02 03:00 review (`docs/planning/chunks/` file).** **2026-10-02 (00:00 CDT slot, §9 item 151):** every imported band green — (P_w′) 1.6089e-03, (A_w) 3.4979e-02, (E) ≤ 9.3e-13, (∞) 3.089e-13 — w → 0 strings equal both `MAT-12` records, w-dropped control 0.43826 misses `A_BAND`; census 68 → 69, docrefs exit 0; 37 s flagged at `-n 2` (`20261002T051055Z_EX-73.log`). *Was:* ⬜ **Done-when (written 2026-10-01 10:30 review, rule (j); §9 item 151):** `examples/materials/07_sar_heats_perfused_ball.py` + same-stem guide with `## Setup figure`, on `MAT-14`'s fixture (`MAT-12`'s submesh, R 0.01, h 100, w 24 000, complex, `-n 2`): (P_w′) ≤ `P_BAND`, (A_w) ≤ `A_BAND`, (E) ≤ `E_BAND` on every run and recomputed, (∞) ≤ `INF_BAND`, all imported and printed beside the gate record; the w = 0 legs reproduce `MAT12_P_RECORD` / `MAT12_A_RECORD`; the w-dropped q_h control misses `A_BAND`; combined XDMFs and the rise CSV named by full filename; census +1, docrefs `exit != 1`, elapsed recorded ⇒ ✅. | standard |
| `EX-74` | *(history and §9 briefs: `docs/planning/chunks/EX-74.md`)* **How long does a gel ball take to settle? Conduction to a Robin surface against the eigen-series** (`mat:8`; `example-runner`; §5.4 ramp for Phase 8 — opened 2026-09-30 by the weekly on `MAT-9` ✅). Imports `MAT-9`'s records; a cooling curve and the steady profile in ParaView. | ✅ **Audit PASS 2026-10-02 03:00 review (`docs/planning/chunks/` file).** **2026-10-02 (00:00 CDT slot, §9 item 152):** (i) 1.7085e-03 string-equal to `MAT9_RECORD`, balance ≤ 3.8e-14 over 13 runs, quasi-static control 0.94145 misses; printed BE / CN rates 0.984 / 2.001 on their records, t₉₉ FEM 87 864 s vs series 88 118 s (−2.9e-3), late slope −3.649 vs −λ₁² −3.637; settle and cooling curves written; census 69 → 70, docrefs exit 0; 31 s flagged at `-n 2` (`20261002T051536Z_EX-74.log`). *Was:* ⬜ *Re-scoped 2026-10-01 10:30 review (rule (j)): the eigen-series comparison and the steady profile are already shown by `mat:4` / `mat:6`; this example is the time axis.* **Done-when (§9 item 152):** `examples/materials/08_gel_ball_settling_time.py` + guide with `## Setup figure`; on `MAT-9`'s ball at 0.01 m, complex, `-n 2`: (i) ≤ `INTERIOR_BAND` on `MAT9_RECORD`; balance ≤ `BALANCE_BAND` on every run; the quasi-static control misses (i); the backward-Euler / CN rates, long-time limit, settle curve to 3 t_R, t₉₉ and the cooling curve printed and written (XDMF and CSV, full filenames); no gate edit; census +1, docrefs `exit != 1` ⇒ ✅. | standard |
| `EX-75` | *(history and §9 briefs: `docs/planning/chunks/EX-75.md`)* **Four-port quadrature on the tuned copper F-human — the k = 0-free drive, for completeness** (side branch; re-posed 2026-10-10 by the weekly on H9′'s fixture — §10 H19; the full row text is in the chunk file) | ✅ 2026-10-11 (§9 item 379, `mri:18`; 206 s flagged, 194 s control; figures census 105 → 106, `broken=0`; docrefs `dead=0`). Was ⬜ Queued as §9 item 379 (`mri:18`). **Done-when (written 2026-10-11 03:00 review, rule (j), from the weekly's re-posed row):** `examples/mri/18_f_human_four_port_quadrature.py` + same-stem guide (three required headings + `## Setup figure` + `## Related`, PNG ≤ 600 KiB); one `build_tuned_quadrature` build (`with_z50=False`), complex, degree 1, `-n 8`; every band imported, none restated, no gate or `src/` edit; H9′'s two-port gates re-asserted at their imported bands ((A0), (P_RC), (S_rep), (P_w), (Lin), (Σ_P), (Sense) > 1, (M_B) ≤ `C4_COVARIANCE_BAND`), record rows equal to item 358's W1 log at the printed digits; the four-port drive (slots 0 / 4 / 8 / 12, quadrature phases, same C) with its construction identities asserted and its k-spectrum, `\|B₁⁺\|`, C4 / C16 spreads and `Y_in` printed, not gated; both drives' `\|B₁⁺\|` on one colour scale; census +1, `broken=0`, no new dead docref; elapsed recorded ⇒ ✅. Not claimed: tuning, absolute field, homogeneity. | heavy (`-n 8`) |
| `EX-76` | *(history and §9 briefs: `docs/planning/chunks/EX-76.md`)* **The tuned birdcage's own modes: the eigen pencil puts the k = 1 pair at 64 MHz and the k = 2 singlet above it** (`ports:20`; `example-runner`; §5.4 ramp for `TH-17` ✅ `3ed8c4e9`; opened by the 2026-10-01 10:30 review, step 5 — no example calls `solve_general_mesh_modes`). **Done-when:** `examples/ports/20_birdcage_eigenmodes.py` + same-stem guide (three required headings + `## Setup figure`, `ports_20_birdcage_eigenmodes_setup.png`); imports `TH-17`'s fixture, helpers, κ and bands, never restates them; on the hole mesh at `C_mode`: (a) the k = 1 pair within κ of 64 MHz, (b) split ≤ `ADJACENT_SPREAD_BAND` and f_pair < f₂, (d′) the √2 law at 2·`C_mode` within κ; negative control: the 2·`C_mode` pair misses 64 MHz by > κ (record 27.01 κ); exported eigenvectors' residual ≤ 1e-6; `ports_20_birdcage_eigenmodes_combined.xdmf` (DG0 \|E\| of both k = 1 partners and k = 2, unit peak); additive `return_vectors` (`core/cavity.py`, `solve_modes`) disclosed with step 3 re-gated green; census +1, `missing=0 broken=0`. Not claimed: class labels, the circuit side, AED, B₁, Q. | ✅ **Audit PASS 2026-10-02 10:30 review (chunk file).** **2026-10-02 (09:00 slot, §9 item 156):** flagged run `20261002T140934Z_EX-76.log`, 71 s, `-n 4`: (a) 0.486 / 0.491 κ, (b) split 5.573e-05, f₂ 78.39 MHz, (d′) 0.356 κ, residuals ≤ 9.5e-12, 2C control 27.01 κ, all equal to W3's records at the printed digits; the unflagged control repeats them; re-gate 13 passed (`…140514Z_EX-76-th17-regate.log`, 198 s); census 73 → 74, `missing=0 broken=0`. | heavy |
| `EX-77` | *(history and §9 briefs: `docs/planning/chunks/EX-77.md`)* **A rod-tip stand-in in the gel: the slender spheroid's tip concentrates SAR ≈ 50× over the tip ball, and only when the field runs along it** (`imp:2`; `example-runner`; §5.4 ramp for `IMP-2` ✅ `cd438c8f`; opened by the 2026-10-01 18:00 review, step 5 — `imp:1` is the sphere; no example builds a spheroid or shows an orientation dependence). **Done-when:** `examples/implants/02_pec_spheroid_tip_in_gel.py` + same-stem guide (three required headings + `## Setup figure`) on the existing `imp:` group; it imports `IMP-2`'s fixture, helpers and bands and never restates them. At aspect 5, degree 2, `-n 2`: \|β − 1\| ≤ `TH15_BETA_BAND`, tip-ball \|Φ_tip,h/Φ_tip,cf,h − 1\| ≤ `MAT4_BAND`, cavity dofs > 0 at max \|E\| < 1e-12, each printed beside the W1 record. Negative control: the transverse drive on the same mesh misses both bands. Combined XDMF (DG0 SAR ratios axial / closed form / transverse, tip mask, CG1 E). The additive `return_fields` on `_solve_spheroid` is disclosed and the module re-gated green in-slot (rule (a)). Census +1 from the fresh pre-census, `broken=0`; docrefs `exit != 1`. Not claimed: aspect 2, a rod, a coil, finite k, 1 g SAR, heating, ASTM F2182, convergence. | ✅ **Audit PASS 2026-10-02 03:00 review (`docs/planning/chunks/` file).** **2026-10-02 (00:00 CDT slot, §9 item 162):** at aspect 5, (i) \|β − 1\| 0.9021 % ≤ `TH15_BETA_BAND`, (ii) Φ_tip,h 49.805748 vs cf 49.904262 (0.1974 %) ≤ `MAT4_BAND`, (iii) cavity 23 886 dofs, max 0; cells 44 304 and β 0.990979 equal the W1 record; transverse control β 4.4e-5, Φ miss 98.56 % — both outside their bands; rule (a) `return_fields` keyword on `_solve_spheroid`, re-gate 15 passed (`20261002T051831Z_EX-77-imp2-regate.log`); census 70 → 71, docrefs exit 0; 113 s flagged at `-n 2`, 6.06 GiB (`20261002T052321Z_EX-77.log`). *Was:* ⬜ Queued §9 item 162. | heavy |
| `EX-78` | *(history and §9 briefs: `docs/planning/chunks/EX-78.md`)* **Where the tuned birdcage's 1/Q goes: Q from the field's stored energy against the Y_in slope, loaded and unloaded, and the phantom's share** (`mri:5`; `example-runner`; §5.4 ramp for `WF-5` ✅ `76362ff5`; opened by the 2026-10-01 18:00 review, step 5 — `ans:7` and `ports:19` read Q off the reactance slope only, `ports:18` splits P_coil / P_phantom on the untuned 10 MHz coil; no example forms ω₀W/P, solves the unloaded coil, or shows Δ(1/Q)). **Done-when:** `examples/mri/05_birdcage_q_loaded_unloaded.py` + same-stem guide on the existing `mri:` group; it imports `WF-5`'s helpers, grids and bands and the tracked κ, never restating them. Cu column on `ANS-7`'s Route A hole mesh (80 181 cells), degree 1, `-n 4`: (A0) loaded secant f₀ / Q reproduce `rows_route_a.cu` ≤ `RECORD_RTOL`; loaded \|Q_field − Q_slope\|/Q_slope ≤ `B_Q_FROZEN`; (A1) P_phantom,U = 0.0; (A2) the unloaded discrete identity ≤ `DISCRETE_IDENTITY_RTOL`; (A3) unloaded Q two ways ≤ `B_Q_U_FROZEN`; (A4) \|Δ(1/Q)_slope − Π\|/Π ≤ B_Δ in-run. Each is printed beside the `WF-5` record. Negative control: the W_sheets-omitted budget misses all three Q bands. Combined XDMF whose map sums reproduce P_phantom / P_surf ≤ `DISCRETE_IDENTITY_RTOL`. No gate edit; census +1, `broken=0`; docrefs `exit != 1`. Not claimed: the PEC column, phantom resolution, absolute Q, tuning, AED. | ✅ **Audit PASS 2026-10-02 03:00 review (`docs/planning/chunks/` file).** **2026-10-02 (00:00 CDT slot, §9 item 163):** on `WF-5`'s Cu fixture (80 181 cells, -n 4): loaded Q_slope 919.0373774 on record (1.4e-13), Q two ways 1.0586e-3 ≤ `B_Q`, unloaded 1.6934e-3 ≤ `B_Q,U`, P_phantom,U 0, Δ(1/Q) vs Π 6.30e-4 ≤ B_Δ, identities ≤ 1.6e-11, map sums ≤ 2.8e-15; W_sheets-omitted control misses all three bands (0.2454 / 0.2449 / 0.2457); census 71 → 72, docrefs exit 0; 88 s flagged, 2.38 GiB (`20261002T053014Z_EX-78.log`). *Was:* ⬜ Queued §9 item 163. | standard |
| `EX-91` | *(history and §9 briefs: `docs/planning/chunks/EX-91.md`)* **The first study through the installed API: a `Study` on the lossy sphere reproduces `MAT-4`'s SAR** (`mri:7`; `example-runner`; opened 2026-10-05 by the weekly, §5.4 ramp — Phase 9 has seven gating chunks ✅ and no example; angle: the API itself) **Done-when:** `examples/mri/07_study_lossy_sphere.py` + same-stem guide (three required headings + `## Setup figure` + `## Related`, PNG ≤ 600 KiB); builds `WF-15`'s `Study` (geometry generator, `MAT-1` preset, imposed field, degree 2 by default — no `degree=` keyword), imports `WF-15`'s record and bands and `MAT-4`'s closed form, restates none; at `-n 2`, complex: mean SAR within `L0_BAND` of `WF-15`'s record and within `MAT4_BAND` of the closed form, both asserted in-script; control: a σ-shifted `Study` misses `MAT4_BAND`; combined XDMF (\|E\|, SAR) opens in ParaView; census +1, `broken=0`, docrefs green; elapsed recorded ⇒ ✅. Not claimed: anything `WF-15` does not. | ✅ **Audit PASS 2026-10-06 18:00 review** (`…202424Z_EX-91.log:185–198`; the 10th-digit flag/control drift sits 9× inside `L0_BAND` and the Done-when has no digit clause; the guide's "digit-equal" sentence corrected by the review). 2026-10-06 (§9 item 272, `20261006T202424Z_EX-91.log` 297 s, control `…202928Z_EX-91-control.log` 304 s): (Deg) 2, (R) 9.512e-11 ≤ `L0_BAND`, (G) 0.0413 %, (c1) 14.05 %, (c2) 3.49e-02, all asserted green; census 86 → 87, `broken=0`. The flagged and control SAR differ at the 10th digit (`-n 2` MUMPS, known-issues 2026-09-10) — disclosed for the review. | heavy (one degree-2 solve, 282 s measured on `WF-15`) |
| `EX-92` | **Four-port S from one `Study`: the F-small birdcage on exact rectangles, with its Touchstone file** (`mri:8`; `example-runner`; opened 2026-10-05 by the weekly, §5.4 ramp; angle: ports and S through the API) **Done-when:** `examples/mri/08_study_four_port_birdcage.py` + same-stem guide (as `EX-91`); builds `WF-16`'s `Study` at degree 1, 64 MHz, κ off, imports item 192's `P1` record, `LEG_D0_REPRODUCTION_BAND` and `PORT-11`'s gate bands, restates none; at `-n 2`, complex: S equals `P1` at the imported band, `PORT-11`'s reciprocity / passivity / class-spread gates green, the written `.s4p` read back through `scikit-rf` equals S at `PORT-7`'s imported band — all asserted in-script; combined XDMF (\|E\| for the P1 drive) opens in ParaView; census +1, docrefs green; elapsed recorded ⇒ ✅. Not claimed: degree 2, tuning, absolute accuracy beyond §2.2's licence. | ✅ *(**DONE 2026-10-07 (04:30 slot, §9 item 278), first try:** (W) width 2; (R) \|S − P1\|/\|P1\| 9.557e-11; (G) reciprocity 4.377e-16, σ_max 0.999760093664, spreads 0.0747 / 0.0691 / 0.0897 %, separation 451.71×; (T) `.s4p` via `skrf` 6.855e-14, file vs P1 9.556e-11; (c1) 4.9103e-02, (c2) 0.2492 (rotation 6.185e-04 printed); XDMF sweep S = Study S. `20261007T094132Z_EX-92.log:1885–1900` (51 s), control `…094301Z_EX-92-control.log` (45 s, same readings at the band); census 88 → 89 `broken=0`, docrefs 109 → 110 guides `dead=0`. The XDMF is written; it was not opened in ParaView by the slot. History: `docs/planning/chunks/EX-92.md`.)* | standard (51 s measured) |
| `EX-93` | *(history and §9 briefs: `docs/planning/chunks/EX-93.md`)* **The same studies from YAML: `femem run` on two files gives the Python route's outputs** (`mri:9`; `example-runner`; opened 2026-10-05 by the weekly, §5.4 ramp; angle: the YAML / CLI layer) **Done-when:** `examples/mri/09_study_from_yaml.py` + same-stem guide (as `EX-91`) + the two YAML files beside it; drives `fem_em_solver.cli:main` in-process on the sphere and four-port YAMLs (degree 1, `-n 2`), imports `WF-1`'s bands, restates none; asserted in-script: the YAML-built `Study` specs equal the Python-built ones field by field by `float.hex`, the SAR and S equal `MAT4_DEG1_RECORD` / `P1` within `L0_BAND` / `LEG_D0_REPRODUCTION_BAND`; control: a malformed YAML (unknown key) raises `StudyConfigError` naming the field; XDMF of the sphere's SAR opens in ParaView; census +1, docrefs green; elapsed recorded ⇒ ✅. Not claimed: a stable schema. | ✅ **Audit PASS 2026-10-06 18:00 review** (`…213328Z_EX-93.log:1997–2011`, census `…213628Z_EX-93-post-figures.log:124`). 2026-10-06 (GREEN; `20261006T213328Z_EX-93.log` 73 s, control 66 s) | standard (W1 100 s / W2 142 s measured on `WF-1`) |
| `EX-94` | **A temperature from one `Study`: steady conduction on the SAR-heated sphere** (`mri:10`; `example-runner`; opened 2026-10-05 by the weekly, §5.4 ramp; angle: the thermal result through the API; **serial on §10 U10 — `WF-15` step 2 — landing**) **Done-when:** `examples/mri/10_study_temperature.py` + same-stem guide (as `EX-91`); `Study.temperature(RobinCooling)` on `WF-15`'s sphere, imports `MAT-10`'s (A) / (P) records and (B) / (C) bands, restates none; at `-n 2`, complex: (A) reproduced at `.4e`, (B) / (C) inside their imported bands, asserted in-script; control: h → 10 h moves the surface rise by the closed form's factor; combined XDMF (SAR, θ) opens in ParaView; census +1, docrefs green; elapsed recorded ⇒ ✅. Not claimed: transient heating, a coil, tissue. | ✅ *(**DONE 2026-10-07 (04:30 slot, §9 item 279), first try, degree 1 per the brief:** (R) (A) 3.2464e-02 / (P) 3.1923e-03 equal `MAT10_A_PRINT` / `MAT10_P_PRINT`; (G) (A) 3.246e-02, (P) 3.192e-03, (B) 8.0e-16, (C) 9.5e-14 inside their imported bands; (Once) no second EM solve; (Sub) 19 496 cells; surface ratio 1.033687 (printed). Controls: (c-h) h → 10 h (A) comparator 6.0941e-01 > `A_BAND` (predicted [0.55, 0.70]), surface-rise ratio FEM 0.100000 vs the closed form's 0.1 to 6.2e-13 (printed); (c1) missing ½ 1.0649e+00 = its record. `20261007T094810Z_EX-94.log:116–130` (21 s), control `…094905Z_EX-94-control.log` (19 s, `.4e` prints identical); census 89 → 90 `broken=0`, docrefs runnable 89 → 90 `dead=0`. The XDMF is written; it was not opened in ParaView by the slot. History: `docs/planning/chunks/EX-94.md`.)* | standard (21 s measured, degree 1) |
| `EX-95` | **A frequency sweep from one `Study`: 10 / 64 / 128 MHz on one mesh, one multi-frequency Touchstone** (`mri:11`; `example-runner`; opened 2026-10-05 by the weekly, §5.4 ramp; angle: the sweep and the multi-frequency export; **serial on §10 U9 — `WF-16` step 2 — landing**) **Done-when:** `examples/mri/11_study_frequency_sweep.py` + same-stem guide (as `EX-91`); `Study.solve_sweep` on `WF-16`'s F-small at degree 1, `-n 2`, imports item 192's three `P1` records and `LEG_D0_REPRODUCTION_BAND`, restates none; asserted in-script: every frequency block equals its record at the imported band, `PORT-11`'s gates green per frequency, the `.s4p` carries exactly three frequencies and reads back through `scikit-rf`; XDMF of \|E\| at 64 MHz opens in ParaView; census +1, docrefs green; elapsed recorded ⇒ ✅. Not claimed: resonance, tuning, degree 2. | ✅ *(**DONE 2026-10-07 (04:30 slot, §9 item 280), first try:** (W) width 2; (M) one mesh, 117 206 cells; (R′) 6.601e-11 / 9.557e-11 / 9.447e-11 at 10 / 64 / 128 MHz; (T′) 3 frequencies, 4 ports, \|df\| 0, skrf ≤ 6.855e-14, file vs P1 ≤ 9.556e-11; (G′) green per frequency, min separation 113.78× at 10 MHz. Controls: (c4) min off-diagonal 0.5524 (64 → 128); (c5) 10 ↔ 128 swap 1.2763 / 1.2367; (c1′) degree-2 gaps 0.0102 / 0.0491 / 0.0555 (printed). `solve_sweep` 53.0 s (mesh 27.7 s + three ≈ 8.3 s sweeps). `20261007T095333Z_EX-95.log:1944–1964` (65 s), control `…095449Z_EX-95-control.log` (65 s, same at the band); census 90 → 91 `broken=0`, docrefs runnable 90 → 91 `dead=0`. The XDMF is written; it was not opened in ParaView by the slot. History: `docs/planning/chunks/EX-95.md`.)* | standard (65 s measured) |
| `EX-96` | Price a `Study` before you solve it: `Study.estimate()` on the lossy sphere, then the solve it priced (`mri:12`) | ✅ DONE 2026-10-07 (§9 item 296, first try): the estimate's 482 158 unknowns / 74 020 cells equal the solved space and `WF-13`'s records, the solve's mean SAR 1.527e-10 from `WF15_DEG2_RECORD` at `L0_BAND` and 0.0413 % from the closed form, leave-one-out P2 inside ×2; predicted 15.36 GiB / 246.5 s against measured 20.25 GiB / 263.8 s printed. Census 91 → 92; the XDMF was written and not opened in ParaView by the slot. *Re-pointed 2026-10-16 (item 490) after 471a: (c2) now asserts the refusal on a P1-only table; the full table prices degree 1, printed.* *History: `docs/planning/chunks/EX-96.md`.* | heavy (277 s measured) |
| `EX-97` | A temperature from a config file: `femem run` on a YAML with `thermal:` reproduces `MAT-10`'s records (`mri:13`) | ✅ DONE 2026-10-07 (§9 item 297, first try): `femem run` returns 0 with its `temperature` block equal to `run_study_file`'s, (A) / (P) equal `MAT10_A_PRINT` / `MAT10_P_PRINT` at `.4e`, (B) 6.8e-15 / (C) 1.3e-13, (IT) 4.874e-12; the 2G control exceeds `A_BAND` and `thermal:` on a port YAML raises. Base step U20 ratified by the 2026-10-07 weekly; census 92 → 93; the XDMF was written and not opened in ParaView by the slot. *History: `docs/planning/chunks/EX-97.md`.* | standard (50 s measured) |
| `EX-98` | **The container wall on the coil-driven gel phantom: `MAT-18` step 3's series-resistance wall under four birdcage drives** *(history and §9 briefs: `docs/planning/chunks/EX-98.md`)* (`mat:14`; `example-runner`; opened 2026-10-07 by the 10:30 daily review, §5.4 ramp — `MAT-18` step 3 recorded the wall on `MAT-15` P1's fixture and no example shows it; angle: the wall on the coil-driven phantom) **Done-when:** `examples/materials/14_coil_driven_phantom_in_container_wall.py` + same-stem guide (as `EX-79`); a rule (a) lift of the step-3 fixture body into `_run_walled(*, return_fields=False)`, with the module re-gated green in-slot and disclosed; `-n 2`, complex; every band imported, none restated; asserted in-script: (R) ≤ `R_BAND`, (E′) / (E) ≤ `E_BAND`, (L0) ≤ `B_BAND` (helper equivalence), C4 θ pairs ≤ `C4_COVARIANCE_BAND`, controls (d) > `D_CONTROL_FLOOR` and (W-ctl) > `W_CONTROL_FLOOR`; the wall effect printed beside the step-3 record; combined XDMF (θ with and without the wall, their difference) written; census +1, `broken=0`, docrefs `exit != 1`; elapsed recorded ⇒ ✅. Not claimed: a conforming wall, wall heat storage, contact resistance, a real container, absolute temperature. | ✅ 2026-10-07: rule (a) lift `_run_walled(*, return_fields=False)` re-gated green (31 passed, 265.55 s, `20261007T172528Z_EX-98-mat18s3-regate.log`); flagged `20261007T173058Z_EX-98.log` (202 s) and control `20261007T173430Z_EX-98-control.log` (205 s) both Status 0; (R) 4.985e-16, (E′) ≤ 1.145e-14, (E) 4.656e-14, (L0) helper equivalence ≤ 5.551e-16, C4 0.0402 / 0.0210 / 0.0292 / 0.0484 %, (d) 1.000000, (W-ctl) 4.63e-2 on all four; census 94 → 95, `broken=0`, docrefs exit 2 (stale only). | heavy (`-n 2`, ≈ 200–300 s per window) |
| `EX-99` | A duty-cycled coil drive heats the gel phantom: the on/off history, Duhamel's sum and the half-power stand-in (`mat:15`) | ✅ 2026-10-07 (audit PASS): the rule (a) lift re-gated green (19 passed), (Δ) 3.719e-14, (D′) 2.120e-15, (E) ≤ 3.819e-14, (Hist) ≤ 4.9e-15, (Δ-ctl) 0.298210; census 95 → 96. Item 306 (2026-10-08) printed the 10 g ball record — duty 1.251111962 K, P1 2.587441071 K, ratio 0.967065 — digit for digit with W1 `:1975`. *History: `docs/planning/chunks/EX-99.md`.* | heavy (`-n 2`, ≈ 250 s predicted, `-t 590`) |
| `EX-100` | **Probe readouts at the benchmark's centres: what a probe log of the coil-heated gel phantom reports** *(history and §9 briefs: `docs/planning/chunks/EX-100.md`)* (`mat:16`; `example-runner`; opened 2026-10-07 by the 10:30 daily review, §5.4 ramp from `MAT-19` step 2's record, §9 item 258; `mat:13` probes `MAT-12`'s sphere and `mat:9` heats this phantom without probes) **Done-when:** `examples/materials/16_probe_readouts_in_coil_phantom.py` + same-stem guide (three required headings + `## Setup figure` + `## Related`, PNG ≤ 600 KiB). The script composes `MAT-19` step 2's fixture from imported helpers (`_run_heated`, `_transient`, `_probe_set`) and imports every band and record, restating none, with no gate-module edit. Complex build. At `-n 2` it asserts (A0), (R0), (A) ≤ `B_BAND`, (U4) and (Σ) ≤ `C4_COVARIANCE_BAND`, plus the controls (C1″) > `WIDTH_BAND` and (f″) > `F_CONTROL_FLOOR`. A `-n 4` control window asserts (W) ≤ `WIDTH_BAND` against the flagged run's hash-guarded npz, and (U3c) iff G_w > 0. Printed only: ball/point at A–E, the 0.1 mm placement shift and an a-sweep at B. Outputs: a combined XDMF (41 steps: θ, q_h, probe masks) and a series CSV. Census +1, `broken=0`; docrefs `exit != 1`; elapsed recorded ⇒ ✅. Not claimed: a real probe's sensing volume, mass or lead; any `ANS-9` comparison; h-convergence; absolute temperature (1 W normalisation). | ✅ 2026-10-07: no lift, no gate-module edit; flagged `-n 2` `20261007T184931Z_EX-100.log` (202 s) and control `-n 4` `20261007T185259Z_EX-100-control-n4.log` (190 s) both Status 0; (A0) 58 866 cells, (R0) 2.547e-15, (A) 2.382e-14 / 8.664e-15, (U4) True, (Σ) 5.9190e-04, (W) `-n 4` vs flagged npz 8.526e-11 ≤ `WIDTH_BAND`; controls (C1″) 4.3978e-03, (f″) 0.8018, (U3c) 2.6047e-02 separate; probe log at t\* and (P4) equal to step 2's records; census 96 → 97, `broken=0`, docrefs exit 2 (stale only). | heavy (`-n 2` + `-n 4`; step-2 windows measured 248–290 s) |
| `EX-101` | **A sheet drawn inside a box becomes an integration surface: native interior facet tags and `dS` on a ghosted mesh** *(history and §9 briefs: `docs/planning/chunks/EX-101.md`)* (`mesh:13`; `example-runner`; opened 2026-10-07 by the 10:30 daily review, §5.4 ramp; angle: geometry/mesh input for the Phase 9 (9b) external-CAD user; base step §10 U11) **Done-when:** `examples/meshing/13_interior_sheet_surface_integral.py` + same-stem guide; real build, `shared_facet` ghosts, `-n 2` and `-n 4`; imports `GEO-37`'s fixture, closed forms and `STRIP_BAND`, restates none; asserted in-script: the reduced owned tag-201 count equals gmsh's triangle count, and ∫1, ∫g('+'), ∫(u⁺ + u⁻) over `dS(201)` and ∫1 `ds(101)` meet their closed forms at `STRIP_BAND`; controls: `dS(999)` assembles to 0.0 while `validate_required_port_tags_exist` raises on 999, and one planted dropped facet leaves the band; the combined XDMF (planted fields, rank, facet tags) is written; census +1, docrefs green; elapsed recorded ⇒ ✅. Not claimed: a port, a STEP import (U12 BLOCKED), an unghosted multi-rank run. | ✅ 2026-10-07: `mesh:13` runs at `-n 2` (4.0 s in-script, 6 s wall, `20261007T172140Z_EX-101.log`) and `-n 4` (`20261007T172154Z_EX-101-control-n4.log`); every anchor and control held at both widths; census 93 → 94, `broken=0`, docrefs `dead=0 guide=0`; outcome in `docs/planning/chunks/EX-101.md`. | standard (smoke by prediction) |
| `EX-102` | **Refine the mesh at an implant tip: the CAD does not move, the cap-ball cells grow as (r/4 ÷ h)³** *(history and §9 briefs: `docs/planning/chunks/EX-102.md`)* (`mesh:14`; `example-runner`; opened 2026-10-07 by the 18:00 daily review, §5.4 ramp; angle: the `implant_rod["tip_size"]` knob for the Phase 7 user who refines at a lead tip; base `GEO-36` step 4, §9 item 298) **Done-when:** `examples/meshing/14_implant_tip_refinement.py` + same-stem guide; real build, no solve, `-n 2` (flagged) and `-n 4` (control); imports `GEO-36` step 4's fixture, rungs, record and bands, restates none; asserted in-script per rung (unset, r/(4√2), r/8): (A0) CAD masses, wall area and centroid bitwise equal to the unset build's, `tip_size_m` = request; (A1) / (A3) at `CAD_ANALYTIC_BAND`, (A4) in `TERMINAL_AREA_BAND`, (A5) census, (A6) at `EXACT`·l; (D) N1curl d2 = 2E + 2F, d1 = E; (G) N_ball strictly increasing; unset cells against `ROD_1G_RUNG_CELL_RECORD` at `CELL_COUNT_BAND`; control (c1): two out-of-range `tip_size` values refused before gmsh; three combined XDMFs and the setup figure written; census +1, docrefs green; elapsed recorded ⇒ ✅. Not claimed: a solve, η_rod or tip SAR, an h-convergence statement, I15's solve windows, a ball-radius knob, the kept-cells build. | ✅ 2026-10-08: `mesh:14` runs at `-n 2` (flagged, 234.8 s in-script, 239 s wall, `20261008T003337Z_EX-102.log`) and `-n 4` (240.6 s, `20261008T003809Z_EX-102-control-n4.log`); every anchor and (c1) held at both widths, all integers width-invariant and equal to the `GEO-36` step 4 W2 records; N_ball ×2.7497 / ×7.8013 vs cubic 2.828 / 8.000 (printed); census 97 → 98, `broken=0`, docrefs `dead=0 guide=0`; outcome in `docs/planning/chunks/EX-102.md`. | heavy (by measurement) |
| `EX-103` | **The absorbing boundary on a small box: at 10 MHz the first-order ABC is a magnetic wall, not an open boundary** *(history and §9 briefs: `docs/planning/chunks/EX-103.md`)* (`th:16`; `example-runner`; opened 2026-10-08 by the 18:00 daily review, §5.4 ramp from `TH-5` step 2c, §9 item 320) **Done-when:** `examples/time_harmonic/16_abc_two_torus.py` + same-stem guide; complex build, `-n 2` (flagged and unflagged control); imports both `TH-5` gate modules' solver, constants and bands, restates none; one rule (a) lift (`return_fields=False` keyword on `_solve_z_column`) re-gated green at `-n 2`; asserted in-script: (P0) at `P0_BAND_PP`, (E) Re Z₁₁^ABC > 0 and PEC \|Re Z₂₁\| ≤ `IDENTITY_TOLERANCE`·\|Z₂₁\|, (S) B^ABC > 0 > B^PEC at 0.08 / 0.10 / 0.12 m, (c1) and (M) / (c2) at `MIX_BAR_DIVISOR`; r(W) printed as a record, never asserted; printed rows equal the step 2c logs to their printed digits; combined XDMFs and the setup figure written; census +1, no new dead docref; elapsed recorded ⇒ ✅. Not claimed: an open boundary at k₀R < π, the 64 / 128 MHz birdcage box, second-order ABC or PML. | ✅ 2026-10-08 (22:30 slot, §9 item 337): run of record `20261009T033559Z_EX-103.log` 241 s `-n 2`, every anchor green and every printed row equal to the step 2c records; re-gate W1 7 passed, control unflagged `-n 2` identical digits; census 98 → 99, docrefs dead 0 → 0. | heavy |
| `EX-104` | **Thermal dose and the six-minute average on closed-form histories: what CEM43 and a 360 s window report on a duty-cycled gel ball** *(history and §9 briefs: `docs/planning/chunks/EX-104.md`)* (`mat:17`; `example-runner`; opened 2026-10-08 by the 18:00 daily review, §5.4 ramp from `MAT-20` steps 1–2 closed-form half; no example computes a dose or a six-minute average) **Done-when:** `examples/materials/17_thermal_dose_duty_cycle.py` + same-stem guide; `-n 2`; imports `MAT-20`'s two gate modules' helpers and bands, restates none, no gate or `src/` edit; closed-form histories only (ramps, `MAT-9`'s step response at 37 °C, `MAT-17` case A periodic), never an FEM history; asserted: (A1), (A2), (B1), (B2) at the imported bars with recomputed rows equal to the gate logs' records; power six-minute max = 1 and whole-period mean = D; the centre duty dose exceeds the average-power dose under an asserted min T > 43 °C; controls (Rc), (Wc), (S1c) miss; closed-form P1 map XDMF + period series + centre CSV; census +1, no new dead docref; elapsed recorded ⇒ ✅. Not claimed: a dose limit or six-minute-SAR compliance, tissue / perfusion, `MAT-20`'s FEM clause, a coil source. | ✅ 2026-10-08 (22:30 slot, §9 item 338): run of record `20261009T035111Z_EX-104.log` 45 s `-n 2`, (A1) / (A2) / (B1) / (B2) / (Rc) / (Wc) equal the `MAT-20` step 1 / 2 logs at the printed digits, (S1) / (S2) / (S1c) / (J) green; unflagged control identical; census 99 → 100, docrefs dead 0. | standard |
| `EX-105` | **The capacitor ESR sheet on the copper birdcage: what 0.2 Ω does to Q** *(history and §9 briefs: `docs/planning/chunks/EX-105.md`)* (`ports:21`; `example-runner`; opened 2026-10-08 by the 18:00 daily review, §5.4 ramp from `PORT-24` step 2, §9 item 309) **Done-when:** `examples/ports/21_capacitor_esr_sheet.py` + same-stem guide; complex build, `-n 2`, all three arms (r = 0.2, r = 0, planted) in one run; imports `PORT-24` step 2's fixture body (rule (a) lift `_run_port24_arm`, additive return keys / keyword, gate re-run green at `-n 2` on all three arms in the same slot), helpers and bands, restates none; asserted in-script: every sheet's identity ≤ `REDUCTION_BAND` at r = 0.2 and r = 0, r = 0 vs the capacitor-only sheet ≤ `LEG_D0_REPRODUCTION_BAND`, the power identity ≤ `DISCRETE_IDENTITY_RTOL`; control: the planted 0.4 Ω P2 > `REDUCTION_BAND`; printed rows equal item 309's logs; Q_rest(r)/Q_rest(0) and the ESR loss share printed beside 1/Q = 1/Q_rest + 1/Q_C, never asserted; combined XDMF and setup figure; census +1, no new dead docref; elapsed recorded ⇒ ✅. Not claimed: an absolute Q, `f₀`, `S₁₁` or tuning figure; a measured ESR; degree 2; F-human; `PORT-24` step 1. | ✅ 2026-10-08 (22:30 slot, §9 item 339): run of record `20261009T040158Z_EX-105.log` 336 s `-n 2`, (I) / (Z) / (P) / (R0) green on every sheet and arm, planted 0.4 Ω misses at 10.95× band, every printed row equal to the `PORT-24` step 2 logs; three-arm re-gate green (same-width control, ruling (4)); census 100 → 101, docrefs dead 0. | heavy (≈ 335 s predicted) |
| `EX-107` | **Figures that equal their records** *(history and §9 briefs: `docs/planning/chunks/EX-107.md`)* (`mri:15`; `example-runner`; opened 2026-10-08 by the 18:00 daily review, §5.4 ramp from `WF-8` ✅, §9 item 328) **Done-when:** `examples/mri/15_figures_from_records.py` + same-stem guide; complex build, `-n 2`; the script writes its own three records through the gated fixtures (`WF-16` sweep `.s4p`, `WF-15` SAR XDMF, `MAT-19` probe CSV) and imports `post.figures`, `femem figure` and `WF-8`'s comparators, restating no band; asserted in-script: (S) 16 / 16 \|S\| / arg S lines `np.array_equal` (float64) to scikit-rf; (X) `real_SAR` slab colours and centroids equal to the h5py parse; (P) probe columns equal to a `csv` parse; (F) rank-0-only non-empty PNGs; (CLI) rc `[0, 0, 0]`; controls: a 1e-6 plant gives exactly one mismatching index per figure, column swaps are caught; `WF-8` W0 re-gated green after an additive path-keyword lift; census +1, no new dead docref; elapsed recorded ⇒ ✅. Not claimed: the sources' correctness, any S / SAR / temperature accuracy, style. | ✅ 2026-10-09 (00:00 slot, take-next, §9 item 341): run of record `20261009T051318Z_EX-107.log` 112 s `-n 2`, (S) 16 / 16, (X) 2 772 of 19 496, (P) 41 rows, (F) / (CLI) green, all bitwise; (c1) [1] / [0] / [5], (c2) 3 / 2 772 / 40; path-keyword lift re-gated 17 passed / 4 skipped (`…051030Z_EX-107-wf8-regate-W0.log`, `[wf8]` lines = W0); control 110 s; census 102 → 103, docrefs dead 0. | heavy (by ceiling, ≈ 125 s predicted) |
| `EX-106` | **10 g SAR and B₁⁺ from one `Study`: one field-keeping solve, two C4 identities on the four-port birdcage** *(history and §9 briefs: `docs/planning/chunks/EX-106.md`)* (`mri:14`; `example-runner`; opened 2026-10-08 by the 18:00 daily review, §5.4 ramp from `WF-16` steps 4–5, §9 items 321 / 322; `mri:8` shows S only) **Done-when:** `examples/mri/14_study_sar_and_b1.py` + same-stem guide; `WF-16`'s `Study` at degree 1, 64 MHz, κ off, `-n 2`, public methods only; imports every band, restates none; asserted in-script: (R) S vs `P1` ≤ `LEG_D0_REPRODUCTION_BAND`, (One) one field-keeping solve serves both, four 10 g C4 pairs and three \|B₁⁺\| images ≤ `C4_COVARIANCE_BAND`, every sample point valid, whole-phantom coverage ≤ `EXACT_IDENTITY_RTOL`; controls far-side ball and mis-rotation > the band (sign); printed rows equal the item 321 / 322 W1 logs at print precision; combined XDMF; control run at the band; census +1, no new dead docref; elapsed recorded ⇒ ✅. Not claimed: absolute SAR or B₁⁺, compliance, spatial peak, homogeneity, degree 2. | ✅ 2026-10-09 (00:00 slot, §9 item 340): run of record `20261009T050323Z_EX-106.log` 74 s `-n 2`, (R) / (One) / (ii) / (i) / (iv) / (B1) / (B3) green, worst 10 g pair 4.5024 % (0.90× band), controls 88 % / 35.0040 % above it, all 42 `REC` digits equal to the item 321 / 322 W1 logs; control run digit-identical; census 101 → 102, docrefs dead 0. | standard by prediction (≈ 95 s), heavy by ceiling |
| `EX-108` | **The `Study` on a copper birdcage: the same F-small four-port through the API on copper and PEC, and where the power goes** *(history and §9 briefs: `docs/planning/chunks/EX-108.md`)* (`mri:16`; `example-runner`; opened 2026-10-08 by the 18:00 daily review, §5.4 ramp from `WF-16` step 7, §9 item 324) **Done-when:** `examples/mri/16_study_copper_birdcage.py` + same-stem guide; complex, degree 1, 64 MHz, `-n 2`; imports `WF-16` step 7's gate module, records and bands, restates none; asserted in-script: both Studies 80 798 cells, Cu and PEC class means vs `ANS6_EXACT_STRIP_CLASS_MEANS_P1["64"]` at `LEG_D0_REPRODUCTION_BAND`, `p_surf` > 0 with the identity residual ≤ `DISCRETE_IDENTITY_RTOL`, `TH-14` gates on both columns, per-cell copper loss sums to `p_surf`; control (c1) PEC Study misses the Cu record > `NEGATIVE_CONTROL_FLOOR`; loss partition printed; combined XDMF + setup figure; control run at the band; census +1, no new dead docref; elapsed recorded ⇒ ✅. Not claimed: absolute S / Z_in, resonance, tuning or Q (`ANS-6` / `ANS-7` 🟡), degree 2. | ✅ 2026-10-09 (00:00 slot, take-next, §9 item 342): run of record `20261009T052039Z_EX-108.log` 76 s `-n 2`, (M) 80 798 cells both, (R) Cu ≤ 1.376e-10 / PEC ≤ 5.875e-11 vs the step 7 records, (b) `p_surf` residual 8.507e-13, (L) 3.211e-15; (c1) PEC misses Cu ≥ 2.811e-04 (= gate `:1826`); partition P_coil / P_in 0.4478617 printed (a record on this fixture at degree 1; `ANS-6` stays 🟡); control 75 s; census 103 → 104, docrefs dead 0. | standard (≈ 85 s predicted) |
| `EX-109` | *(history and §9 briefs: `docs/planning/chunks/EX-109.md`)* **The tuned copper F-human birdcage under the clinical two-port 90° drive: \|B₁⁺\| and the 10 g SAR hotspot at C_in,F-human = 118.2437 pF, r = 0.2 Ω** (`mri:17`; `example-runner`; opened 2026-10-10 18:00 daily review from H9′'s gates — `WF-10` 🟡, named for the weekly to ratify; no example shows a tuned coil at human scale) **Done-when:** `examples/mri/17_f_human_tuned_quadrature_b1_sar.py` + same-stem guide (three required headings + `## Setup figure` + `## Related`, PNG ≤ 600 KiB); one `build_tuned_quadrature` build, `with_z50=False`, complex, degree 1, `-n 8`; imports every band, restates none, no gate or `src/` edit. Asserted at imported bands, printed beside items 358 / 359's W1 logs: (A0), (P_RC), (S_rep), (P_w) and (Lin) per sense, (Σ_P), (Sense) > 1, the feed-bisector mirror (M_B) ≤ `C4_COVARIANCE_BAND`, coverage (i), the 10 g mirror at c\* and over the lattice, containment (iv); the DG0 phantom-loss map sums to P_phantom ≤ `EXACT_IDENTITY_RTOL`. Asserted controls: the mis-paired sense and the unswapped single feed both miss the band. C16, the azimuthal spread and the P21 partner are printed only. Outputs: `mri_17_f_human_tuned_quadrature_b1_sar_combined.xdmf` and a lattice CSV. Census +1, `broken=0`, no new dead docref; elapsed recorded ⇒ ✅. Not claimed: four-port (`EX-75`), f₀ / Q / S₁₁ / tuning accuracy, homogeneity, absolute B₁⁺ or SAR, compliance, C16, degree 2, AED. **✅ 2026-10-10 (19:30 slot, §9 item 362):** every gate green in the flagged run and in the unflagged control at `-n 8` (241 s / 234 s); record rows equal both W1 logs at the printed digits, except one round-off digit on Im Y17 under `cw`; (X1) 4.441e-16; census 104 → 105 with `missing=0 broken=0`; docrefs `dead=0` after one guide reword (`20261011T004348Z_EX-109.log`, `…T005311Z_EX-109-post-docrefs-r2.log`). Opened from H9′ for the weekly to ratify. | ✅ | heavy (`-n 8`) |
| `EX-110` | **The six-minute window on an FEM heating history: windowing a duty-cycled coil drive equals windowing its ten switch-on runs** *(history and §9 briefs: `docs/planning/chunks/EX-110.md`)* (`mat:18`; `example-runner`; opened 2026-10-11 by the 18:00 daily review, §5.4 ramp from `MAT-20` ✅, §9 item 370; `mat:17` windows closed forms only and `mat:15` shows this history unwindowed) **Done-when:** `examples/materials/18_six_minute_window_coil_heating.py` + same-stem guide (three required headings + `## Setup figure` + `## Related`, PNG ≤ 600 KiB); complex, `-n 2`; rule (a) lift in `tests/validation/test_mat20_fem_history.py` (`_fem_history(*, keep_duty_fields=False)`, `_identities(t, duty, hs)`; additive, prints unchanged), module re-gated green in-slot; every band imported, none restated, no `src/` edit. Asserted in-script through the imported helper: (H0) ≤ `L_BAND`, (FEM) and (FEM-2d) ≤ `D_BAND`, control (Wc-FEM) > `D_BAND`; example-side (Cm) the volume mean of the per-vertex six-minute field equals the windowed series ≤ `D_BAND`, its one-sample-short control (Cm-c) misses; records printed beside `20261011T172810Z_MAT-20-step2-fem-W1.log:1957–1974`. Outputs: θ and its six-minute mean as a 29-frame series XDMF, a t\* combined XDMF and a history CSV. Census +1, `broken=0`, no new dead docref; elapsed recorded ⇒ ✅. Not claimed: a limit or compliance, tissue, perfusion, absolute temperature, a periodic state, a CEM43 gate. | ✅ **Audit PASS 2026-10-12 18:00.** 2026-10-12 (12:00 slot, §9 item 399, attempt 2): owned-slice maxima guarded (`initial=-np.inf`); every anchor green and both controls miss in the flagged run and the unflagged control at `-n 2` (214 s / 206 s), records equal W1 at the printed digits; lift re-gated 13 passed (278 s); allgather `[(11906, 11906), (0, 0)]`; census 106 → 107, `missing=0 broken=0`, docrefs `dead=0` (`20261012T170911Z_EX-110.log`). | heavy (`-n 2`, ≈ 225 s per window, `-t 590`) |
| `EX-111` | **A disc-edge stand-in in the gel: a flat PEC body's rim concentrates the field only under an in-plane drive, and less than a rod tip of the same curvature** *(history and §9 briefs: `docs/planning/chunks/EX-111.md`)* (`imp:08`; `example-runner`; opened 2026-10-12 by the 18:00 daily review, §5.4 ramp from `IMP-6` ✅, §9 item 405; `imp:02` is the prolate tip, no example builds a flat body) **Done-when:** `examples/implants/08_pec_oblate_disc_edge_in_gel.py` + same-stem guide (three required headings + `## Setup figure` + `## Related`, PNG ≤ 600 KiB); complex, degree 2, `-n 2`, aspect 5, a = 6 mm; rule (a) lift `_solve_oblate(…, *, return_fields=False)` in `tests/validation/test_imp6_oblate_fem.py` (additive), module re-gated green in-slot; every band imported, none restated, no `src/` edit. Asserted through the imported helper: \|β − 1\| ≤ `TH15_BETA_BAND`, rim-ball \|Φ_rim,h/Φ_rim,cf,h − 1\| ≤ `MAT4_BAND`, cavity dofs > 0 at max \|E\| < 1e-12, each beside `20261012T184537Z_IMP-6-step1-W1b.log:189–192`; the axial drive on the same mesh misses both (> band; W2 `:190–191`). Printed only: oblate 1/L_x against `imp:02`'s prolate 1/L_z at aspect 5 and the same 6 mm long semi-axis. Outputs: combined XDMF (DG0 SAR ratios in-plane / closed form / difference / axial, rim mask, CG1 E) and a rim-line CSV of \|E\| against the closed form. Census +1, `broken=0`, docrefs `dead=0`; elapsed recorded ⇒ ✅. Not claimed: a square-edged plate, aspect 2, a coil, finite k, SAR averaging, heating, convergence. | ✅ **Audit PASS 2026-10-13 03:00.** 2026-10-12 (§9 item 405, 19:30 slot): (i) \|β − 1\| 0.6343 %, (ii) rim Φ miss 1.4181 %, (iii) cavity 19 952 dofs at 0, all equal to W1b; axial control misses both (1.0018, 97.8352 %); 86 s flagged (`20261013T003602Z_EX-111.log`), re-gate 16 passed in 237 s; census 108/108, `dead=0`. | heavy (`-n 2`, ≈ 85 s per run, `-t 400`; re-gate ≈ 237 s, `-t 590`) |
| `EX-112` | The F2182-sized gel box: a 90 mm slab to 3 % at its centre, cooled ≈ 80 % through its two large faces, and hours from steady (`mat:19`) | ✅ 2026-10-15 18:00 review: the two gating logs landed at `dcdbc42c` (operator, ruling (a)) carrying item 438's readings — (S1) 1.417055e-03, (T) 9.071528e-04, Status 0 (`20261014T033257Z_EX-112.log:56`, `:106`, `:136`; control `…033435Z_EX-112-control.log:49`, `:99`, `:129`), census 108 → 109; audit PASS (no hash had been recorded, so "unchanged" rests on the digits). *Was* 🧪 2026-10-14 03:00 on custody only. **Done-when:** `examples/materials/19_f2182_gel_box_slab_and_faces.py` + same-stem guide (three required headings + `## Setup figure` + `## Related`, PNG ≤ 600 KiB); complex, `-n 2`, the gate's m = 6 rung through its module-level helpers; no gate-module edit, no `src/` edit; every band imported, none restated. *History: `docs/planning/chunks/EX-112.md`.* | heavy (`-n 2`, ≈ 90 s per run, `-t 420`) |
| `EX-113` | `Study.tune`: the mode capacitance is read from Im z_k only, cancels mode 1 exactly, and the mode fields show which mode it tuned (`mri:19`) | ✅ 2026-10-14 (item 446, audit PASS 10:30): `mri:19` — `Study.tune` on the PEC hole and solid F-small Studies, every gate anchor equal to `20261014T004318Z_WF-16-s6-W1.log:2765–2781`, the example-side ESR / tuned-cancellation identities at `EXACT_IDENTITY_RTOL`, the mode fields printed; census 109 → 110, `broken=0`, docrefs `dead=0`. **Done-when:** `examples/mri/19_study_tune_mode_capacitance.py` + same-stem guide (three required headings + `## Setup figure` + `## Related`, PNG ≤ 600 KiB); complex, `-n 2`, the gate's PEC hole and solid Studies through the API; no `tests/` or `src/` edit; every band imported. *History: `docs/planning/chunks/EX-113.md`.* | heavy (`-n 2`, ≈ 125 s per run, `-t 420`) |
| `EX-114` | `femem run` thermal on the tuned copper F-human from YAML: one EM solve answers the file's thermal knobs, and the temperature field the API keeps to itself (`mri:20`) | ✅ 2026-10-15 (item 453, after items 447 / 449 BLOCKED on the record-check tolerance and a `snaps[-1]` `KeyError`, the one-line fix licensed by the 18:00 review): `mri:20` — `femem run` thermal on the tuned copper F-human from the tracked YAML at `-n 8`, one EM solve, every anchor and control green; census 110 → 111, `broken=0`, docrefs `dead=0`. **Done-when:** `examples/mri/20_f_human_thermal_from_yaml.py` + same-stem guide (three required headings + `## Setup figure` + `## Related`, PNG ≤ 600 KiB); complex, `-n 8`, `run_study_file` on the tracked YAML in place; no `tests/` or `src/` edit; every band imported. *History: `docs/planning/chunks/EX-114.md`.* | heavy (`-n 8`, ≈ 185 s per run, `-t 590`) |
| `EX-115` | `femem run` with a top-level `tune:` section (`mri:21`) | ✅ 2026-10-15 (items 461 / 477): `mri:21` — `femem run` with a `tune:` section in-process on the tracked file, every anchor at its imported band, the k = 2 and malformed-file controls miss; census and docrefs clean. **Done-when:** `examples/mri/21_femem_run_tune_section.py` + same-stem guide (three required headings + `## Setup figure` + `## Related`, PNG ≤ 600 KiB); complex, `-n 2`, `cli.main(["run", …])` in-process on the tracked `f_small_hole_pec_tune.yaml`, captured by a call-through wrapper; no `tests/` or `src/` edit; every band imported. *History: `docs/planning/chunks/EX-115.md`.* | heavy (`-n 2`, ≈ 85 s per run, `-t 420`) |
| `EX-116` | The three shipped presets through `femem`, priced then run (`mri:22`) | ✅ 2026-10-16 (§9 item 488): `mri:22` — the three presets priced at one width, the two degree-1 presets run with each price inside ×2 of its own run, `benchmark-lite` priced only; controls miss; census and docrefs clean. **Done-when:** `examples/mri/22_presets_priced_then_run.py` + same-stem guide (three required headings + `## Setup figure` + `## Related`, PNG ≤ 600 KiB); complex, `-n 2`; presets read in place; no `tests/`, `src/` or preset edit; every band imported from `test_wf4_presets.py`. *History: `docs/planning/chunks/EX-116.md`.* | heavy (`-n 2`, ≈ 170 s per run, `-t 590`) |
| `EX-117` | The clinical-scale body coil with the ASTM F2182 gel box in ParaView (`mesh:15`) | ✅ 2026-10-16 (§9 item 489): `mesh:15` — F-clinical-1.5T with the ASTM F2182 box in ParaView, `GEO-38` step 2's gates recomputed on this mesh, the box's surface and reach identities asserted; census and docrefs clean. **Done-when:** `examples/meshing/15_f_clinical_coil_astm_box.py` + same-stem guide (three required headings + `## Setup figure` + `## Related`, PNG ≤ 600 KiB); real build, `-n 4`, `_box_kwargs()` unchanged, mesh only; no `tests/` or `src/` edit; every band imported. *History: `docs/planning/chunks/EX-117.md`.* | heavy (`-n 4`, ≈ 450 s, `-t 600`, durable capture) |
| `EX-118` | `Study.tune` on the 32-port copper F-human and its k = 1 ring field (`mri:23`) | ✅ 2026-10-16 (§9 item 492): `mri:23` — `Study.tune` on the 32-port copper F-human in one solve returns the H12c `C_mode`, k = 2 misses, the k = 1 ring field exported; census and docrefs clean. **Done-when:** `examples/mri/23_f_human_study_tune_ring.py` + same-stem guide (three required headings + `## Setup figure` + `## Related`, PNG ≤ 600 KiB); complex, `-n 8`, `g._study()` unchanged (degree 1, Leontovich hole, κ off); no `tests/` or `src/` edit; every band imported. *History: `docs/planning/chunks/EX-118.md`.* | heavy (`-n 8`, ≈ 390 s per run, `-t 590`) |
| `EX-120` | Tune then load: `Study.tune(k=1)`'s `C_mode` on three gap sheets through `CapacitorSheets` | ✅ 2026-10-16 (audit PASS 18:00; §9 item 495): `mri:24` — tune then load through `CapacitorSheets`, the reduced network equal to the circuit reduction and the loaded field's power identity green; controls miss; census and docrefs clean. **Done-when:** `examples/mri/24_study_tune_then_capacitor_sheets.py` + same-stem guide (three required headings + `## Setup figure` + `## Related`, PNG ≤ 600 KiB); complex, `-n 2`, `g._study` unchanged (degree 1, PEC hole, κ off, lossless); no `tests/` or `src/` edit; every band imported. *History: `docs/planning/chunks/EX-120.md`.* | heavy by ceiling (`-n 2`, ≈ 83 s per run, `-t 420`) |
| `EX-121` | Point probes read straight off the gel-ball submesh | ✅ 2026-10-16 (audit PASS 18:00; §9 item 498): `mat:21` — `parent=` point probes on the gel-ball submesh, owner cells identical at `-n 2` / `-n 4`, every gate anchor at its band; census and docrefs clean. **Done-when:** `examples/materials/21_point_probes_on_gel_ball_submesh.py` + same-stem guide (three required headings + `## Setup figure` + `## Related`, PNG ≤ 600 KiB); complex, run of record `-n 2` plus a `-n 4` control; no `tests/` or `src/` edit; every band imported. *History: `docs/planning/chunks/EX-121.md`.* | heavy by ceiling (`-n 2` + `-n 4`, ≈ 70 s per run predicted, `-t 300`) |
| `EX-122` | A PEC disc in the coil-driven phantom: rim 1 g SAR with / without (`imp:9`) | ✅ 2026-10-16 (audit PASS 18:00; §9 item 499): `imp:9` — the PEC disc in the coil-driven phantom with / without, rim 1 g SAR and the PEC-versus-natural rim map, every anchor equal to W1; controls green; census and docrefs clean. **Done-when:** `examples/implants/09_pec_disc_in_birdcage_rim_sar.py` + same-stem guide (three required headings + `## Setup figure` + `## Related`, PNG ≤ 600 KiB); complex, `-n 2`, W1 record meshes through `test_imp3_disc_birdcage`'s `_case_disc` / `_report` / `_control_c1_disc` unchanged; no `tests/` or `src/` edit; every band imported. *History: `docs/planning/chunks/EX-122.md`.* | heavy (`-n 2`, ≈ 204 s, `-t 590`, durable capture) |
| `EX-123` | The `Study` transient API: a duty-cycled gel ball and a walled gel ball with no EM solve | ✅ 2026-10-17 (§9 item 501): `mri:25` — the `Study` transient API on the duty-cycled and walled gel balls with no EM solve, `MAT-17` A / B and `MAT-18` 2(b) reproduced at their records; controls miss; census and docrefs clean. **Done-when:** `examples/mri/25_study_transient_duty_and_wall.py` + same-stem guide (three required headings + `## Setup figure` + `## Related`, PNG ≤ 600 KiB); complex, `-n 2`, gate helpers from `test_wf15_study_transient` unchanged (complex-safe `m13._interior_error` substituted); no `tests/` or `src/` edit; every band imported. *History: `docs/planning/chunks/EX-123.md`.* | heavy by ceiling (`-n 2`, ≈ 250 s predicted, `-t 590`, durable capture) — measured 194 s |
| `EX-124` | `femem run` on a `capacitor_sheets` file | ✅ 2026-10-17 (§9 item 502): `mri:26` — `femem run` on a `capacitor_sheets` file: Route A equals the circuit reduction, the summary equals its Study, YAML equals Python, and the ESR file written from the tune summary carries r = 0.2 Ω; census and docrefs clean. **Done-when:** `examples/mri/26_femem_run_capacitor_sheets.py` + same-stem guide (three required headings + `## Setup figure` + `## Related`, PNG ≤ 600 KiB); complex, `-n 2`, four real builds (`cli.main(["run", …])` ×3, `g._py_caps()`), no build cache; no `tests/` or `src/` edit; every band imported. *History: `docs/planning/chunks/EX-124.md`.* | heavy by ceiling (`-n 2`, ≈ 165 s predicted per run, `-t 420`) |
| `EX-125` | The Pennes sink on the coil-driven gel phantom at four w (`mat:22`) | ✅ 2026-10-17 (§9 item 503): `mat:22` — the Pennes sink at four w on the coil-driven gel phantom through `MAT-15` step 4's own body, every identity at its band, the perfusion share printed per w; census and docrefs clean. **Done-when:** `examples/materials/22_pennes_sink_coil_driven_phantom.py` + same-stem guide (three required headings + `## Setup figure` + `## Related`, PNG ≤ 600 KiB); complex, `-n 2`; rule (a) lift — additive `return_fields=` on `_run_perfused` (keys `base`, `results`), re-gated green in-slot, disclosed; no other `tests/` or `src/` edit; every band imported. *History: `docs/planning/chunks/EX-125.md`.* | heavy (`-n 2`, ≈ 275 s predicted + ≈ 305 s re-gate, `-t 590`, durable capture) |
| `EX-126` | Two `femem run`s of one study file: the output bundle and its diff | ✅ 2026-10-17 (§9 item 504): `mri:27` — two `femem run`s of one file diffed, deterministic keys bitwise and S within `L0_BAND`, each plant caught by its layer; census and docrefs clean. **Done-when:** `examples/mri/27_femem_run_bundle_reproducible.py` + same-stem guide (three required headings + `## Setup figure` + `## Related`, PNG ≤ 600 KiB) + the one-ULP σ control YAML; complex, `-n 2`, `f_small_four_port.yaml` in place, three CLI runs into fresh gitignored dirs; no `tests/` or `src/` edit; comparator and bands imported from `test_wf2_bundle_reproducible`. *History: `docs/planning/chunks/EX-126.md`.* | heavy by ceiling (`-n 2`, ≈ 140 s predicted, `-t 420`) |
| `EX-127` | **The 10 g lattice under the four-port drive on the tuned copper F-human: the 12 z_c balls pass C4 at m = 1/2/3, the two-port misses, and the spread per radius reads drive anisotropy** (`mri:28`; `example-runner`; opened 2026-10-16 18:00 daily, §5.4 ramp from `WF-10` step 5 (H24, §9 item 422), §9 item 505; `mri:18` reads four-port fields only, `mri:17` the two-port lattice without C4) *(history: `docs/planning/chunks/EX-127.md`)* **Done-when:** `examples/mri/28_f_human_four_port_sar10g_c4.py` + same-stem guide (three required headings + `## Setup figure` + `## Related`, PNG ≤ 600 KiB); complex, `-n 8`, `build_tuned_quadrature(with_z50=False)` + the gate's `_four_port_drive` unchanged (degree 1, κ off, untuned four-port at the two-port root); no `tests/` or `src/` edit; every band imported. Asserted beside `20261013T141801Z_WF-10-step5-W1.log:6324–6405`: (C4₄) max_m ≤ `C4_COVARIANCE_BAND`; builder anchors, (S_rep), (P_w)/(Lin), (X4-K/P/L), (Σ_P)₄, (i)₄, (iv); records at 10^-(d−2). Example-side: (O) the z_c rows are C4 images ≤ `EXACT_IDENTITY_RTOL`; (X1)₄/(X1)₂ DG0 loss sums = P_phantom. Control: two-port max_m > band (asserted, backed 17.38 %); m = 2 blindness printed. Printed: spread per radius beside the two-port's and the verdict. Combined XDMF + z_c CSV. A leak-blocked log ⇒ BLOCKED. Census +1, `broken=0`, docrefs `dead=0`, elapsed recorded ⇒ ✅. Not claimed: tuning, f₀/Q/match, absolute SAR, compliance, h, degree 2, an external check. | ✅ 2026-10-17 (item 506, 04:30 slot; item 505 had been BLOCKED on a rank-0 tail crash, `20261017T033328Z_EX-127.log:6293`): `20261017T093909Z_EX-127.log` Status 0, `[capture] rc=0`, 376 s, (C4₄) 1.8127 / 0.7215 / 0.2645 % (`:6222`), control 17.3831 % > band (`:6223`), verdict `drive anisotropy` (`:6228`), 68 / 68 records max rel 7.962e-05 (`:6237`); figures 123 → 124 `broken=0`, docrefs `dead=0` exit 2; audit PASS 2026-10-17 10:30) | heavy (`-n 8`, ≈ 310 s, `-t 590`, durable capture) |
| `EX-128` | A σ = 800 S/m rod in the tuned copper F-human's gel phantom under the two-port drive: with / without on one mesh, the tip 1 g ratio, the adjoint fold at human … | ✅ 2026-10-17 (§9 item 512, attempt 2 of 4; audit PASS 18:00): `imp:10` — a σ = 800 S/m rod in the tuned copper F-human's phantom with / without on one mesh, every gate anchor equal to W1b, the tip 1 g ratio and map ratios printed; census and docrefs clean. **Done-when:** `examples/implants/10_sigma_rod_in_tuned_f_human_tip_sar.py` + same-stem guide (three required headings + `## Setup figure` + `## Related`, PNG ≤ 600 KiB); complex, `-n 8`, `_run_rod`'s body repeated with provenance (`0e7d3abd`, `:410–847`) through `_build_f_human_rod` / `_read_record` unchanged (kept-rod mesh, degree 1, κ off, r = 0.2 Ω, C_in 118.2437 pF); no `tests/` or `src/` edit; every band imported. *History: `docs/planning/chunks/EX-128.md`.* | heavy (`-n 8`, ≈ 316 s, `-t 590`, durable capture) |
| `EX-129` | **The F2182 gel box under the transverse eddy-SAR source 3Q(y² + z²)/(B² + C²): hot side walls, heat leaving through ±y, the thin-axis term in ParaView** (`mat:23`; `example-runner`; opened 2026-10-17 03:00 daily, §5.4 ramp from `MAT-21` step 9 (§10 P36, provisional), §9 item 509; `mat:19` is uniform, `mat:20` varies in x only) *(history: `docs/planning/chunks/EX-129.md`)* **Done-when:** `examples/materials/23_f2182_gel_box_transverse_eddy_source.py` + same-stem guide (three required headings + `## Setup figure` + `## Related`, PNG ≤ 600 KiB); complex, `-n 2`, m = 4 at P2 through the gate's helpers; no `tests/` or `src/` edit; every band imported. Asserted: (Pdeg); (S1-quad) / (T-quad) ≤ `INTERIOR_BAND` against `robin_box_quadratic_*` at `QUAD_TERMS`; (S2-quad) ≤ `S2_GAUSS_BAND`; (E-quad) ≤ `E_BAND`; (F0) face partition; (F1) Σ face powers = source power; backed controls (c-src), (c-src′), (c-z), (c-src-T) > band; (cF) +y dropped. Records printed beside `20261016T171904Z_MAT-21-step9-W2-quadratic-fem.log`. Combined XDMF + y-centre-line CSV. Census +1, `broken=0`, docrefs `dead=0`, elapsed recorded ⇒ ✅. Withdrawn if the weekly reverts P36. Not claimed: the box's EM SAR, a coil, other weights, a hand-off, an h-rate, compliance, an external check. **✅ 2026-10-17 06:00 slot (audit PASS 10:30; provisional on P36):** green on the first run of record, Status 0, 164 s (`20261017T110217Z_EX-129.log`); every anchor and backed control equals its gate record, (F0) 9 344, ±y share 0.2256 vs uniform 0.1262, census 124 → 125 `broken=0`, docrefs `dead=0`. | ✅ | heavy (`-n 2`, ≈ 145 s, `-t 590`) |
| `EX-130` | **The F2182 gel box through the `Study` API: `ThermalBoxGeometry` → `Study.temperature` → `Study.probe` at the six F2182 positions** (`mat:24`; `example-runner`; §5.4 ramp from `WF-15` steps 5–6, §9 items 521a / 523; queued 2026-10-18 10:30 as §9 item 534) | ✅ Done 2026-10-18 (§9 item 534, run of record `20261018T214035Z_EX-130.log`, Status 0, 21 s). *History: `docs/planning/chunks/EX-130.md`.* **Done-when:** `examples/materials/24_f2182_gel_box_study_api_probes.py` + same-stem guide (three required headings + `## Setup figure` + `## Related`, PNG ≤ 600 KiB); complex, `-n 2`, m = 4 P1; no `tests/` or `src/` edit; every band and record imported from `test_wf15_study_thermal_box`, `test_wf15_study_probe_thermal` and `test_mat21_box_probe_points_fem`. Asserted: (S1-m4) whole `== STEP1_M4_WHOLE`, interior ≤ `INTERIOR_BAND`, Gauss ≤ `S2_GAUSS_BAND`, cells / dofs equal to the record; (I-rep) ≤ `L0_BAND`; (A-box) ≤ `B_BAND`; (R-box) five balls ≤ `INTERIOR_BAND`; (V) six points valid; (M) inset balls inside the box; (c-wall) backed > band; (c-half) > band. Six point values and inset ball means printed beside the series. Combined XDMF + probe-table CSV. Census +1, `broken=0`, docrefs `dead=0`, elapsed recorded ⇒ ✅. | standard (`-n 2`, ≈ 40 s, `-t 180`) |
| `EX-131` | **The F2182 gel box under the duty-cycled transverse eddy-SAR source: temperature histories at the six F2182 probe points** (`mat:25`; `example-runner`; §5.4 ramp from `MAT-21` steps 10–11, P37 / P38, §9 items 520 / 524; queued 2026-10-18 10:30 as §9 item 535) | ✅ Done 2026-10-18 (§9 item 535, run of record `20261018T214758Z_EX-131.log`, Status 0, 189 s). *History: `docs/planning/chunks/EX-131.md`.* **Done-when:** `examples/materials/25_f2182_gel_box_duty_probe_histories.py` + same-stem guide (three required headings + `## Setup figure` + `## Related`, PNG ≤ 600 KiB); complex, `-n 2`, m = 4 at P2 through the gates' helpers; no `tests/` or `src/` edit; every band imported. Asserted: (Pdeg); (U); (V); (Pt) ≤ `INTERIOR_BAND`; (S0); (R36); (D′-quad) ≤ `D_BAND`; (T-qduty) ≤ `INTERIOR_BAND`; (E-qduty) ≤ `E_BAND`; (H40), (H≡), (D0-h) exact; backed controls (c-step), (c-u-pt), (c-z-pt) and predicted (cH) > band. History time-series XDMF + combined XDMF + probe-history CSV. Census +1, `broken=0`, docrefs `dead=0`, elapsed recorded ⇒ ✅. | heavy (`-n 2`, ≈ 215 s, `-t 590`) |
| `EX-132` | **`femem run` on the F2182 gel-box thermal files: three tracked `thermal_box` YAML files through the CLI, the JSON summary read back, θ_file = θ_Python** (`mri:29`; `example-runner`; §5.4 ramp from `WF-1` step 15, §9 item 521b; queued 2026-10-18 18:00 as §9 item 541) | ✅ Done 2026-10-19 (§9 item 541, run of record `20261019T050316Z_EX-132.log`, Status 0, 216 s). *History: `docs/planning/chunks/EX-132.md`.* **Done-when:** `examples/mri/29_femem_run_f2182_gel_box_files.py` + same-stem guide (three required headings + `## Setup figure` + `## Related`, PNG ≤ 600 KiB); complex, `-n 2`, steady m = 4 and duty m = 6, P1, through `cli.main` `check` / `run` on the four tracked `thermal_box` files; no `tests/` or `src/` edit; every band and record imported through `test_wf1_yaml_thermal_box`. Asserted: (K); (UT) cells and 41 multipliers; (I-Py) ≤ `L0_BAND` on both files; (R) whole `== STEP1_M4_WHOLE`, interior ≤ `INTERIOR_BAND`, Gauss ≤ `S2_GAUSS_BAND`, (T-duty) `== P26_T_DUTY`, balance ≤ `E_BAND`; (I=); from the on-disk JSON (P-json) ≤ `L0_BAND`, (S0-json); backed controls (c-wall), (c3), (cX). (G-json) printed only (ruling 18:00). Two combined XDMFs + summary-readback CSV. Census +1, `broken=0`, docrefs `dead=0`, elapsed recorded ⇒ ✅. | heavy (`-n 2`, ≈ 220 s, `-t 590`) |
| `EX-133` | **The F-clinical-3T twin in ParaView: 32 legs, 64 ring-gap sheets, closed PEC can Ø 0.743 × 0.70 m, ASTM F2182 gel box as the load** (`mesh:16`; `example-runner`; the 10-17 weekly's "`mesh:16` (the 3 T twin) after C8"; queued 2026-10-18 18:00 as §9 item 542) | ✅ 2026-10-21 (§9 item 542): `mesh:16`, the F-clinical-3T twin with the ASTM box, with C8 W3's gates recomputed on this mesh exactly. The can reads Ø 0.743 × 0.70 m off the vertices, the box sits 25 mm from each cap, and all four controls separate; census 134 → 135, `broken=0`, docrefs `dead=0` (`20261021T183238Z_EX-133.log`, 315 s). *History: `docs/planning/chunks/EX-133.md`.* **Done-when:** `examples/meshing/16_f_clinical_3t_coil_astm_box.py` + same-stem guide (three required headings + `## Setup figure` + `## Related`, PNG ≤ 600 KiB); real build, `-n 4`, the 3 T box kwargs by C8's own swap, mesh only; no `tests/` or `src/` edit; every band imported. Gate anchors recomputed beside C8 W3 (`20261018T171406Z_GEO-38-C8-W3-box-n4.log`): (T1) 64/64, (T2), (B1), (B2), (V1) ≤ `STRIP_BAND`, (V2) ≤ `CAD_ANALYTIC_BAND`, (K), (Tab) via `_assert_3t_table`; cells printed, not asserted. Example-side: (F), (R3), (Sc) can Ø 0.743 and length 0.70 from mesh vertices, (Zc) cap clearance 0.025 m, all ≤ `STRIP_BAND`, with backed controls. Combined XDMF with `cell_diameter` and the sheet + box facet tags. Census +1, `broken=0`, docrefs `dead=0`, elapsed recorded ⇒ ✅. Not claimed: any solve, port gate, tuning, 128 MHz, B₁/SAR, the C8b widths, the 3 T price. | heavy (`-n 4`, ≈ 335 s, `-t 600`, daytime) |
| `EX-134` | **The F-small four-port case through `Study.solve(backend="palace")`: the Study writes the case, a `fem-em-palace` window runs it, the Study reads S back** (`mri:30`; `example-runner`; §5.4 ramp from `ARCH-1` step 1 (U16, §9 items 517 / 517b / 533); queued 2026-10-18 18:00 as §9 item 544) | ✅ Done 2026-10-19 (§9 item 544, run of record `20261019T093732Z_EX-134.log`, Status 0, 60 s; Palace window `20261019T093301Z_EX-134-palace-d2.log`, 262 s). *History: `docs/planning/chunks/EX-134.md`.* **Done-when:** `examples/mri/30_study_palace_backend_f_small.py` + same-stem guide (three required headings + `## Setup figure` + `## Related`, PNG ≤ 600 KiB); complex, `-n 1`, degree 2, Palace `-np 8` order 2 into a fresh `paraview_output/` case directory; no `tests/` / `src/` / `scripts/` / `docker/` edit; every band imported from `test_arch1_palace_backend.py`. Asserted: (K); (W-case) manifest and config = U16's; (R) cells and order; (E2) ≤ `LEG_D0_REPRODUCTION_BAND`; (P); (Tie2) ≤ `TIE2_BAND` per class; backed (c-order) > `C_ORDER_FLOOR`; (c-freq) printed. Geometry/tag combined XDMF + S-class CSV. Census +1, `broken=0`, docrefs `dead=0`, elapsed recorded ⇒ ✅. Not claimed: CLI, order 1, fields from Palace, accuracy, any AED comparison. | heavy (`-n 1` ≈ 110 s, `-t 180`; Palace ≈ 265 s, `-k 30 590`) |
| `EX-135` | **`femem run` on the transverse-quadratic-source and box-probe `thermal_box` files: P36's record and U60's five probe readings from the JSON summary, θ_file = θ_Python** (`mri:31`; `example-runner`; §5.4 ramp from `WF-1` steps 16–17, §9 items 528 / 539; queued 2026-10-19 03:00 as §9 item 546) | ✅ Done 2026-10-19 (§9 item 546, run of record `20261019T123750Z_EX-135.log`, Status 0, 189 s). *History: `docs/planning/chunks/EX-135.md`.* **Done-when:** `examples/mri/31_femem_run_gel_box_source_and_probe_files.py` + same-stem guide (three required headings + `## Setup figure` + `## Related`, PNG ≤ 600 KiB); complex, `-n 2`, `cli.main` `check` / `run` on seven tracked `thermal_box` files (quadratic steady / CN + uniform P2 at P2; steady, wrong-wall, duty probe files + malformed at P1), none combining source and probes, no new YAML; no `tests/` or `src/` edit; every band and record imported. Asserted: (K); (F-diff), (F-none); (Pdeg), (UT); (I-Py) ≤ `L0_BAND`, (I-tr); (S1-quad) `== P36_S1`, (T-quad) `== P36_T`, Gauss ≤ `S2_GAUSS_BAND`, (U), balance ≤ `E_BAND`; (A-box) ≤ `B_BAND`, (R-box) ≤ `INTERIOR_BAND` (U60 digits printed, `EX-130` precedent); (I=), (U′); from the on-disk JSON (S-json), (P-json), (E-json), (R-json), (C-json); backed controls (c-src′), (c-src), (c-wall), (cX). Three combined XDMFs + two CSVs. Census +1, `broken=0`, docrefs `dead=0`, elapsed recorded ⇒ ✅. Not claimed: source + probes in one file, U66, a complex gate of the probe route, real build, compliance. | heavy (`-n 2`, ≈ 230 s, `-t 590`) |
| `EX-136` | **`femem run` on the duty-cycled transverse-quadratic `thermal_box` file and the two ball-probe `thermal_ball` files: P38's (T-qduty) and U60's three ball-probe readings read back from the CLI's JSON summary, θ_file = θ_Python** (`mri:32`; `example-runner`; §5.4 ramp from `WF-1` steps 18–19, §9 items 540 / 545; queued 2026-10-19 10:30 as §9 item 550) | ✅ Done 2026-10-19 (§9 item 550, run of record `20261019T171706Z_EX-136.log`, Status 0, 190 s). *History: `docs/planning/chunks/EX-136.md`.* **Done-when:** `examples/mri/32_femem_run_ball_probes_and_duty_quadratic_files.py` + same-stem guide (required headings + `## Setup figure` + `## Related`, PNG ≤ 600 KiB); complex, `-n 2`, `cli.main` `check` / `run` on four tracked files (quadratic duty box; m19 and duty ball-probe files + malformed), none combining source and probes, no new YAML, no `tests/` or `src/` edit, every band and record imported. Asserted: (T-qduty) `== P38_T_QDUTY` ≤ `INTERIOR_BAND`, schedule `== SCHEDULE`, balance ≤ `E_BAND`; (A-ball) ≤ `B_BAND`, (T*-ball) ≤ `P_BAND` (U60 digits printed); θ_file = θ_Python ≤ `L0_BAND`; the summaries and probe table read back from disk; controls (c-step), (c-disp), (cX), (cX-duty). Census +1, `broken=0`, docrefs `dead=0`, elapsed recorded ⇒ ✅. Not claimed: source + probes in one file, a complex gate of the ball-probe route. | heavy (`-n 2`, ≈ 210 s, `-t 590`) |
| `EX-137` | **An implant inside the standard's box in a clinical-scale coil: the kept 1 mm × 30 mm rod in F-clinical's ASTM F2182 gel box, in ParaView** (`mesh:17`; `example-runner`; §5.4 ramp from `GEO-36` step 8, item 530, `c07031df`; queued 2026-10-19 18:00 as §9 item 553) | ✅ 2026-10-21 (15:00 slot, §9 item 553): (G0)–(G8) and `REC` reproduce W1 at `-n 4`, (Cap) max d/r − 1 = 6.883e-15 and every control separates; build 453.0 s (gmsh 425.97 s, under the 450 s rule), harness 466 s, census 136 / `missing=0 broken=0`, docrefs `dead=0` (`20261021T200333Z_EX-137.log`). *History: `docs/planning/chunks/EX-137.md`.* **Done-when:** `examples/meshing/17_f_clinical_rod_in_astm_box.py` + same-stem guide (three required headings + `## Setup figure` + `## Related`, PNG ≤ 600 KiB); real build, `-n 4`, the gate's `_kwargs()` imported (rod at `_rod_centre()`, x = Lx/4), mesh only, 450 s gmsh stop rule armed; no `tests/` or `src/` edit; every band imported. Step 8's (G0)–(G8) recomputed beside W1 (`20261019T200702Z_GEO-36-rodbox-W1-n4.log`) at `STRIP_BAND` / `EXACT` / `TERMINAL_AREA_BAND`, `REC` integers exact, cells printed beside 1 217 025, not asserted. Example-side: (F), (R3), (Cap) every tag-403 vertex inside the capsule at `_rod_centre()` ≤ `STRIP_BAND`; controls (N1r), (N2r), (cB) mis-placed rod refused, (cR) step 6's centre, (cF), each labelled per rule (e). Combined XDMF with `cell_diameter`, cell tags and the sheet + box + 402 facet tags. Census +1, `broken=0`, docrefs `dead=0`, elapsed recorded ⇒ ✅. Not claimed: any solve, tip field, SAR, heating, ASTM placement or compliance, coat / hole rod, the `-n 2` width, the rod-solve price. | heavy (`-n 4`, ≈ 465 s, `-t 600`, daytime) |
| `EX-138` | `femem run` on the four-port F-small EM file with point probes: per-port point SAR and \|E\| read back from the JSON summary, file = Python Study, a 1 mm … | ✅ 2026-10-20 (13:30 slot, §9 item 569; audit PASS): one `-n 2` complex run reproduces item 564's records — 117 206 cells, S tie 9.557e-11, (I-Py) 6.005e-13 ≤ `L0_BAND`, the 1 mm displaced probe misses by ≥ 5.6184e-02 — census 132 → 133, docrefs `dead=0` (`20261020T184521Z_EX-138.log`). The first example to run `outputs.probes` on an EM file (its full Done-when is in the chunk file); withdrawn if U70 is reverted. *History: `docs/planning/chunks/EX-138.md`.* | standard by prediction (`-n 2`, ≈ 130 s), heavy by ceiling (`-t 300`) |
| `EX-139` | `femem run` on the four-port F-small sweep file: one multi-frequency `.s4p` and the JSON summary read back from disk, per-frequency S = the `P1` records, file … | ✅ 2026-10-20 (21:00 slot, §9 item 578; audit PASS 2026-10-21 03:00): one `-n 2` complex run reproduces item 572's per-frequency `P1` records at the sweep's three frequencies on one 117 206-cell mesh (S misses 6.601e-11 / 9.557e-11 / 9.449e-11, the `.s4p` read from disk equal, (I-Py) ≤ 6.328e-15, the 10 ↔ 128 MHz swap misses by 1.2367 / 1.2763) — census 133 → 134, docrefs `dead=0` (`20261021T020252Z_EX-139.log`). The first example to run `study.sweep_hz` through `femem run` (its full Done-when is in the chunk file); withdrawn if U71 is reverted. *History: `docs/planning/chunks/EX-139.md`.* | standard by prediction (`-n 2`, ≈ 115 s), heavy by ceiling (`-t 300`) |
| `EX-140` | **The gel heated by an implant's own SAR: `IMP-7`'s metal-core sphere in the two-layer gel ball, in ParaView** (`imp:11`; `example-runner`; §5.4 ramp from `IMP-7` step 1, item 594, `04885e88`; queued 2026-10-22 18:00 as §9 item 603) | ✅ 2026-10-23 (log `20261023T005150Z_EX-140.log`; (S2) 1.2080e-03, (S3) 5.039e-14, (c-iso) 2.7326e-02; census 137). **Done-when:** `examples/implants/11_pec_sphere_implant_heating_in_gel_ball.py` + same-stem guide; `./run_examples.sh -e imp:11 -n 2` (complex build) exits 0 with (S2) interior rel L2 ≤ `S2_BAND`, (S3) balance ≤ `S3_BAND`, (Src′) at `B1_BAND`, max \|Im T_h\| = 0 and the isotropic-source control missing, every band and the fixture imported from `IMP-7`'s modules; combined XDMF (T_h, θ_exact, error, q) and a setup figure written; census +1, `EX-57` `missing=0`, docrefs `dead=0`. *History and brief: `docs/planning/chunks/EX-140.md`.* | standard by prediction (`-n 2`, 30–70 s), `-t 180` |
| `EX-141` | **The same SAR source around a thermally insulating implant: `IMP-7`'s adiabatic hole on the gel-shell submesh, in ParaView** (`imp:12`; `example-runner`; §5.4 ramp from `IMP-7` step 1's hole gate, `20261022T202503Z_IMP-7-step1-W2.log`; opened and queued 2026-10-23 03:00 as §9 item 604) | ✅ 2026-10-23 (log `20261023T093235Z_EX-141.log`, `-n 2` complex, flagged, 22 s harness-wall, Status 0): (G) 150 334 submesh cells = parent tag-2, Robin area −5.0690e-04, rmin/rmax at round-off; (S2) interior rel L2 8.6093e-04 (whole 9.0332e-04); (S3) 2.856e-14; (c-iso) 3.6062e-02 = 3.6× the band; max \|Im T_h\| = 0; every gate record at the digit except (S3), round-off against the record's 2.571e-14 (audit PASS, 10:30 review). Printed: pole / equator 0.392616 / 0.328033 K (closed form 0.394437 / 0.327248); closed-form split hole 0.067189 K vs `imp:11` metal core 0.017949 K, 3.74×; inner flux / source 1.116e-02. One disclosed deviation: `_errors_c` local (the gate's `_errors` squares as `d**2`); no gate-module or `src/` edit. Census 137 → 138, `EX-57` `missing=0 broken=0`, docrefs `dead=0` exit 2 (`20261023T093322Z_EX-141-post-docrefs.log`). **Done-when:** `examples/implants/12_insulating_implant_heating_in_gel_shell.py` + same-stem guide; `./run_examples.sh -e imp:12 -n 2` (complex build) exits 0 with the hole's (G) submesh census, (S2) interior rel L2 ≤ `S2_BAND`, (S3) balance ≤ `S3_BAND`, max \|Im T_h\| = 0 and the isotropic-source control missing, every band, helper and the fixture imported from `IMP-7`'s modules; combined XDMF on the submesh (T_h, θ_exact, error, q) and a setup figure written; census +1, `EX-57` `missing=0`, docrefs `dead=0`. *History and brief: `docs/planning/chunks/EX-141.md`.* | standard by prediction (`-n 2`, 25–60 s), `-t 180` |
| `EX-142` | **Where the absorbing boundary works: a dipole in a sphere at k₀R = π reflects Γ_sph(π) ≈ 0.05, and the first-order ABC reads it** (`th:17`; `example-runner`; §5.4 ramp from `TH-5` step 1b, `20260930T185052Z_TH-5-step1b-W2.log`; opened and queued 2026-10-23 10:30 as §9 item 606 — `th:16` shows only the small-box failure). *History and brief: `docs/planning/chunks/EX-142.md`.* | ✅ 2026-10-23 (`20261023T170527Z_EX-142.log`: 76 s, exit 0, (B)/(E)/(A)/(c-PEC) asserted, every record equal to `W2`; census 138 → 139 as predicted, `missing=0 broken=0` (logged by the review, `20261023T230222Z_EX-57-census-review.log:35`), docrefs `dead=0` exit 2; audited PASS 10-23 18:00, chunk file) **Done-when:** `examples/time_harmonic/17_abc_dipole_sphere.py` + same-stem guide; `./run_examples.sh -e th:17 -n 2 -t 600` exits 0 (complex) with the 1b gate's four asserts — `BAND ≤ BAND_CEILING`, \|Γ_E\| ≤ `BAND` and \|Γ_A − Γ_sph(π)\| ≤ `BAND` at λ₀/40, PEC at λ₀/20 missing Γ_sph by > `BAND` — every band, helper and fixture imported by name from `tests/validation/test_th5_dipole_sphere.py`, records printed beside; combined XDMF + setup figure; census +1, `missing=0`, docrefs `dead=0`. | heavy |
| `EX-143` | **The first-order ABC at oblique incidence: a TE10 guide at 60° reflects R(θ) = −1/3, and the angle-matched condition absorbs** (`th:18`; `example-runner`; §5.4 ramp from `TH-5` step 1a (O), `20260930T094044Z_TH-5-step1a2-W1.log.gz`; opened and queued 2026-10-23 18:00 as §9 item 608 — `th:16` / `th:17` show normal-ish incidence only; one provisional rule (a) hoist for the weekly). *History and brief: `docs/planning/chunks/EX-143.md`.* | ✅ 2026-10-23 (`20261024T003217Z_EX-143.log`: 17 s, exit 0, (O1) 1.161362e-04 / (O2) 2.620831e-04 / (R) / (c-PEC) 0.6667 (x222) asserted; re-gate `20261024T003104Z_EX-143-regate.log` 8 passed; census 134 → 135 group-prefixed / 139 → 140 `--list` entries, figures `missing=0 broken=0`, docrefs `dead=0` exit 2; audited PASS 2026-10-24 03:00). **Done-when:** `examples/time_harmonic/18_abc_oblique_incidence.py` + same-stem guide; two disclosed rule (a) lifts in `test_th5_abc_reflection.py` (`O_BAND = 3e-3` hoisted by value into both (O) asserts; additive `return_fields` on `_te10`), module re-gated green at `-n 2`; `./run_examples.sh -e th:18 -n 2 -t 600` exits 0 (complex) asserting (O1) \|Γ_exact\| ≤ `O_BAND`, (O2) \|Γ_ABC − R\| ≤ `O_BAND`, and PEC missing R by > `O_BAND` (backed, `W1:129`), every band and helper imported by name; combined XDMF + setup figure; census +1, logged `missing=0`, docrefs `dead=0`. | heavy |
| `EX-144` | **A foil thinner than its skin depth on a cavity wall: Q rises 6.24× while the resonance barely moves** (`th:19`; `example-runner`; §5.4 ramp from `TH-18` step 1, `20260929T200249Z_TH-18-step1-W1.log.gz`; opened and queued 2026-10-23 18:00 as §9 item 609 — `th:10` shows the Leontovich cavity, `th:13` the foil on a line). *History and brief: `docs/planning/chunks/EX-144.md`.* | ✅ 2026-10-23 (log `20261024T003645Z_EX-144.log`; (P) 1.998e-08, (i) 1.446e-08, (ii) −0.0017 % / +0.0086 %, (iii) −0.0554 %, (c-bulk) 0.8398; census 135 → 136 group-prefixed / 140 → 141 `--list` entries; audited PASS 2026-10-24 03:00) **Done-when:** `examples/time_harmonic/19_layered_wall_cavity_q.py` + same-stem guide; `./run_examples.sh -e th:19 -n 2 -t 300` exits 0 (complex) with the step-1 gate's six asserts — (P) at `PREMISE_TOLERANCE`, (i) at `T5_TOLERANCE`, (ii) at `THIN_RATIO_TOLERANCE` and `Q_TOLERANCE`, (iii) at `SHIFT_TOLERANCE`, bulk-Leontovich control above `NEGATIVE_CONTROL_FLOOR` — every band, record and constant imported by name, no gate or `src/` edit; combined XDMF + setup figure; census +1, logged `missing=0`, docrefs `dead=0`. | standard |
| `EX-145` | **A dispersive (Debye) material: the plane wave decays as k₀√ε̂(f) at both Larmor frequencies, and ε̂ frozen at the wrong one misses ≈ 12 %** (`th:20`; `example-runner`; §5.4 ramp from `MAT-3` step 1 (gated; the row is 🟡 on step 2), `20260929T183252Z_MAT-3-step1.log.gz`; opened and queued 2026-10-23 18:00 as §9 item 610 — no example uses a frequency-dependent material). *History and brief: `docs/planning/chunks/EX-145.md`.* | ✅ 2026-10-23 (log `20261024T004035Z_EX-145.log`, `-n 2` complex, flagged, 33 s harness-wall / 31.0 s in Python, Status 0): (i) α 0.0262 %, β 0.0554 %, rel L2 3.29632e-02; (ii) α 0.0419 %, β 0.0475 %; (c-a) α miss 11.9990 % (2.4× `BAND`); (c-b) 21.9949 % (4.4× `BAND`); L2 rate 0.9998 printed; all equal to `W1:157–165` at the printed digits. Census 141 → 142 listed entries (predicted +1) (`20261024T003959Z_EX-145-census-pre.log`, `20261024T004126Z_EX-145-census-post.log`), figures `missing=0` (`20261024T004129Z_EX-145-figures.log`), docrefs `dead=0 guide=0` exit 2 (`20261024T004135Z_EX-145-post-docrefs.log`). `MAT-3` stays 🟡. Audited PASS 2026-10-24 03:00. **Done-when:** `examples/time_harmonic/20_debye_plane_wave.py` + same-stem guide; `./run_examples.sh -e th:20 -n 2 -t 600` exits 0 (complex) asserting the step-1 gate's `BAND` asserts — rel L2, α and β at 127.74 MHz (24³), α and β at 63.87 MHz, the wrong-frequency and dispersion-dropped controls missing α by > `BAND` — with `BAND`, `DEBYE` and the helpers imported by name, no gate or `src/` edit; combined XDMF + setup figure; census +1, logged `missing=0`, docrefs `dead=0`; `MAT-3` stays 🟡. | heavy |
| `EX-146` | **The wave impedance of a lossy plane wave: \|E\|/\|B\| from Faraday's law equals ω/\|k\|, and a σ-blind speed misses ≈ 21 %** (`th:21`; `example-runner`; §5.4 ramp from `POST-2` ✅, `20260930T094718Z_POST-2-balance-W1.log.gz`; opened and queued 2026-10-24 03:00 as §9 item 612 — no example checks a Faraday-derived B against a closed form). *History and brief: `docs/planning/chunks/EX-146.md`.* | ✅ 2026-10-24 (log `20261024T093646Z_EX-146.log`, `-n 2` complex, flagged, 31 s harness-wall / 29.2 s in Python, Status 0): (R) +0.0164 %, (S_E) +0.5378 %, (S_B) +0.0002 %, (Bal) 1.638957e-04, (c-σ) −21.1963 % (≈ 4.2× `POST2_BAND`); 12³ rung +0.0649 % / +1.1437 % / +0.0007 % / −21.1580 % and max-ratio miss +0.7185 % printed; all equal to `W1:108–130` at the printed digits. Census 142 → 143 listed entries (predicted +1) (`20261024T093605Z_EX-146-census-pre.log`, `20261024T093739Z_EX-146-census-post.log`), figures `examples=143 ok=143 missing=0 broken=0` (`20261024T093742Z_EX-146-figures.log`), docrefs `dead=0 guide=0` exit 2 (`20261024T093745Z_EX-146-post-docrefs.log`). `POST-2` unchanged. **Done-when:** `examples/time_harmonic/21_plane_wave_impedance.py` + same-stem guide; `./run_examples.sh -e th:21 -n 2 -t 600` exits 0 (complex) asserting `POST-2`'s 24³ asserts — mean-ratio, E- and B-span misses and the balance at ω/\|k\| below `POST2_BAND`, the σ-blind comparand missing by > `POST2_BAND` — with `POST2_BAND` and the helpers imported by name, no gate or `src/` edit; combined XDMF (E and B) + setup figure; census +1, logged `missing=0`, docrefs `dead=0`. | heavy |
| `EX-147` | A PEC disc face-on to the field: the pole enhancement 1/L_z of an oblate spheroid in gel, and the edge-on drive that misses it (`imp:13`; `example-runner`; … | ✅ 2026-10-24 (`imp:13`, audited PASS): the aspect-5 PEC disc face-on reproduces `IMP-6` step 3's axial gate and the in-plane control misses, census +1, docrefs dead=0; the picture solve's β print relabelled (item 615). *History: `docs/planning/chunks/EX-147.md`.* | heavy |
| `EX-148` | **A sampler-free convergence rate: the domain L2 error of \|B\| around a straight wire falls at rate ≈ 1.69, and the natural wall triples it** (magnetostatics 7; `example-runner`; §5.4 ramp from `MAG-18` ✅, `20260824T003650Z_MAG-18-regate-run2.log.gz`; opened and queued 2026-10-24 03:00 as §9 item 614 — `mag:6` names the rate gate but asserts only sampled monotonicity). *History and brief: `docs/planning/chunks/EX-148.md`.* | ✅ 2026-10-24 (log `20261024T094607Z_EX-148.log`, `-n 2` real, flagged, 178 s harness-wall / 176.0 s in Python, Status 0): (M) E_Ω 25.286826 % / 10.617170 % / 6.645807 % at 38 740 / 147 235 / 383 146 cells, (Rate) 1.6854 (2.41× the 0.7 floor), (c-BC) natural wall 32.315493 % against 10.617170 % (ratio 0.3285); h = 0.0025 reading 1.061717e-01 printed beside `E_OMEGA_H0025_RECORD` (not asserted); equal to `20260824T003650Z_MAG-18-regate-run2.log.gz:2681–2685`, `:3636` at the printed digits except rung 1 (2.5286825569e-01 against 2.5286824423e-01, relative 4.5e-8; 10:30 audit). Census 144 → 145 listed entries (predicted +1) (`20261024T094523Z_EX-148-census-pre.log`, `20261024T094910Z_EX-148-census-post.log`), figures `examples=145 ok=145 missing=0 broken=0` (`20261024T094939Z_EX-148-figures-2.log`; the first run `20261024T094912Z_EX-148-figures.log` read missing=1 because the guide's Setup figure section was not yet written), docrefs `dead=0 guide=0` exit 2 (`20261024T094941Z_EX-148-docrefs.log`). `MAG-18` unchanged. **Done-when:** `examples/magnetostatics/07_domain_l2_convergence_rate.py` + same-stem guide; `./run_examples.sh -e 7 -n 2 -t 900` exits 0 (real build) asserting `MAG-18`'s E_Ω ladder monotone, the fitted rate ≥ `E_OMEGA_RATE_MIN`, and the natural wall worse at h = 0.0025, every band and helper imported by name, the record printed not asserted, no gate or `src/` edit; combined XDMF + setup figure; census +1, logged `missing=0`, docrefs `dead=0`. | heavy |
| `EX-149` | **A thermal h-convergence rate: the P1 error on the perfused Robin ball falls at rate ≈ 2.06 against the Pennes closed form** (`mat:26`; `example-runner`; §5.4 ramp from `MAT-11` ✅ gate (ii), `20261001T004409Z_MAT-11-cx-W2.log.gz`; opened and queued 2026-10-24 10:30 as §9 item 616 — no thermal example asserts a spatial rate). *History and brief: `docs/planning/chunks/EX-149.md`.* | ✅ 2026-10-24 (§9 item 616; audit PASS 18:00 review; `20261024T170646Z_EX-149.log`: rate 2.060, interior 1.8285e-03, balance ≤ 3.9e-15, w = 0 control 1.1055; census 145 → 146, figures `missing=0`, docrefs `dead=0`). **Done-when:** `examples/materials/26_perfusion_ball_h_convergence_rate.py` + same-stem guide; `./run_examples.sh -e mat:26 -n 2 -t 420` exits 0 (complex) asserting `MAT-11` gate (ii) — balance ≤ `BALANCE_BAND` per rung, least-squares rate ≥ `RATE_FLOOR`, interior error at 0.01 m ≤ `INTERIOR_BAND` — and its control (a), w = 0 missing `INTERIOR_BAND`, every band and helper imported by name, no gate or `src/` edit; combined XDMF + setup figure; census +1, logged `missing=0`, docrefs `dead=0`. | heavy |
| `EX-150` | **The power-wave identity per drive: an S-column's accepted power equals the terminal accounting on the 4-leg birdcage, and deleting S_kk misses by exactly \|S_kk\|²/(1 − Σ\|S_ik\|²)** (`ports:22`; `example-runner`; §5.4 ramp from `POST-6` ✅ step 1b, `20260923T095500Z_POST-6-step3-rerun-WF-7-step2.log.gz`; opened and queued 2026-10-24 10:30 as §9 item 617 — `ports:13` prints only the drive-level ratio). *History and brief: `docs/planning/chunks/EX-150.md`.* | ✅ 2026-10-24 (§9 item 617; audit PASS 18:00 review; `20261024T171129Z_EX-150.log`, Status 0, 51 s, `-n 2`). **Done-when:** `examples/ports/22_power_wave_identity_per_drive.py` + same-stem guide; `./run_examples.sh -e ports:22 -n 2 -t 600` exits 0 (complex) asserting step 1b per drive — the identity within `STEP1B_POWER_WAVE_RTOL`, the S-diagonal control ≥ `STEP1B_CONTROL_MIN_MISS` and equal to its ceiling at 2 × `STEP1B_POWER_WAVE_RTOL`, the single-drive residuals ≤ `POWER_BALANCE_BAND` — every band and helper imported by name, no gate or `src/` edit; combined XDMF + setup figure; census +1, logged `missing=0`, docrefs `dead=0`. | heavy |
| `EX-151` | **A time-step convergence rate under a switched source: backward Euler on the duty-cycled gel ball converges at order ≈ 1.00, and Crank–Nicolson's order 2 falls outside that band** (`mat:27`; `example-runner`; §5.4 ramp from `MAT-17` ✅ step 2 gate (ii), `20261003T171113Z_MAT-17-step2-W1b.log.gz`; opened and queued 2026-10-24 10:30 as §9 item 618 — `mat:8` prints constant-source rates, unasserted). *History and brief: `docs/planning/chunks/EX-151.md`.* | ✅ 2026-10-24 (§9 item 618; audit PASS 18:00 review; `20261024T171609Z_EX-151.log`, Status 0, 33 s; BE rate 0.996, CN rate 2.104; census 147 → 148, `missing=0`, `dead=0`). **Done-when:** `examples/materials/27_duty_cycle_time_step_rate.py` + same-stem guide; `./run_examples.sh -e mat:27 -n 2 -t 420` exits 0 (complex) asserting step 2 (ii) — (S0) on every run, balance ≤ `BALANCE_BAND`, the BE least-squares rate inside `BE_RATE_BAND` — and the CN + Rannacher rate above `BE_RATE_BAND`'s top, every band and helper imported by name (complex-safe norms from `MAT-13`), no gate or `src/` edit; combined XDMF + ladder CSV + setup figure; census +1, logged `missing=0`, docrefs `dead=0`. | heavy |
| `EX-152` | **The lumped-sheet power identity per sheet, and the sheet excess falling with the sheet width, on the F-small exact-rectangle four-port** (`ports:23`; `example-runner`; §5.4 ramp from `ANS-4` H33, `20261025T095010Z_ANS-4-H33a-W1.log`; opened and queued 2026-10-25 10:30 as §9 item 640 — no `ports:` example computes a per-sheet identity). *History and brief: `docs/planning/chunks/EX-152.md`.* | ✅ 2026-10-25 (§9 item 640; `20261025T200450Z_EX-152.log`, Status 0, 98 s, `-n 2`; census `20261025T200655Z_EX-152-census-post.log` 149, figures `missing=0`, docrefs `dead=0`). **Done-when:** `examples/ports/23_lumped_sheet_excess_vs_width.py` + same-stem guide; `./run_examples.sh -e ports:23 -n 2 -t 600` exits 0 (complex) asserting H33a's gates at f = 0.25 and 0.5, degree 1 — (Wd) ≤ `STRIP_BAND`, (CS) ≤ `DISCRETE_IDENTITY_RTOL` with the planted-V control missing by ≥ `BLIND_FACTOR`×, (P_w) / (P_acc-1) with the conductor-blind control missing, (R0) at the record width ≤ `LEG_D0_REPRODUCTION_BAND` — every band and helper imported by name, no gate or `src/` edit; the excess per width printed; combined XDMF + CSV + setup figure; census +1, logged `missing=0`, docrefs `dead=0`. | heavy |
| `EX-153` | **A graded PEC-backed PML at 60° incidence: the TE10 guide's reflection follows R(0)^cos θ, and the ABC at the same angle reflects ≈ 333× more** (`th:22`; `example-runner`; §5.4 ramp from `TH-5` step 3b, `20261025T111709Z_TH-5-step3b-W1.log`; opened and queued 2026-10-25 10:30 as §9 item 641 — no example shows a PML). *History and brief: `docs/planning/chunks/EX-153.md`.* | ✅ 2026-10-25 (§9 item 641; `20261025T202452Z_EX-153.log`, Status 0, 113 s, `-n 2`; census `20261025T202713Z_EX-153-census-post.log` th:22 added, figures `missing=0`, docrefs `dead=0`). **Done-when:** `examples/time_harmonic/22_pml_oblique_incidence.py` + same-stem guide; `./run_examples.sh -e th:22 -n 2 -t 420` exits 0 (complex) asserting step 3b's gates — Γ_c and the angular identity ≤ `IDENTITY_TOLERANCE`, (B) degree 2 ≤ `PML_REFLECTION_BAR`, (Cabc) > `CABC_RATIO_BAR`, (NC-0) > `NC0_BAR` — every band and helper imported by name, no gate or `src/` edit; combined XDMF + setup figure; census +1, logged `missing=0`, docrefs `dead=0`. | heavy |
| `EX-154` | **The implant-heated gel ball at P2: `imp:11`'s rung, the error cut by element order, the l = 2 pole moment beside the series** (`imp:14`; `example-runner`; §5.4 ramp from `IMP-7` step 2 / I40, `20261025T123436Z_IMP-7-step2-W1.log`; opened and queued 2026-10-25 10:30 as §9 item 642 — no thermal example runs at degree 2). *History and brief: `docs/planning/chunks/EX-154.md`.* | ✅ 2026-10-25 (§9 item 642; `20261025T203007Z_EX-154.log`, Status 0, 157 s, `-n 2`; P2 interior 8.4803e-04, (c-beta0) 5.9510e-02; census 150 → 151, figures `missing=0`, docrefs `dead=0`). **Done-when:** `examples/implants/14_implant_heating_p2_metal_core.py` + same-stem guide; `./run_examples.sh -e imp:14 -n 2 -t 420` exits 0 (complex) asserting on the metal core at `LADDER[1]` balance ≤ `S3_BAND` at P1 and P2, interior rel L2 ≤ `S2_BAND` at P1 and P2, and the β = 0 control at P2 missing `S2_BAND` — every band and helper imported by name (a local complex-safe norm), no gate or `src/` edit; combined XDMF + setup figure; census +1, logged `missing=0`, docrefs `dead=0`. | heavy |
| `EX-155` | **The rod's own SAR heating the coil-driven F-small phantom, read like for like on a matched mesh pair** (`imp:15`; `example-runner`; §5.4 ramp from `IMP-3` I39 / I41, `20261025T201613Z_IMP-3-I41-W1b.log`; opened and queued 2026-10-25 18:00 as §9 item 643 — no example heats a phantom with a rod in it). *History and brief: `docs/planning/chunks/EX-155.md`.* | ✅ 2026-10-26 (§9 item 643; `20261026T033157Z_EX-155.log`, Status 0, 233 s, `-n 2`; (W0) 0 / 4 958, (C4) arm M worst 0.9318 %; tip 1 g θ R/M 1.104515 / 1.108684 printed; census 151 → 152, figures `missing=0`, docrefs `dead=0`; one XDMF per arm, since the arms' cell counts differ). **Done-when:** `examples/implants/15_rod_heating_matched_phantom.py` + same-stem guide; `./run_examples.sh -e imp:15 -n 2 -t 420` exits 0 (complex) asserting (A0m) the matched census exact before any solve, (W0) the rod-wall facet census 0 / parent count, (B) / (E) / (E′) / (S2) on both arms at `B_BAND` / `E_BAND` / `S2_BAND`, and (C4) on the rod-free arm at `C4_COVARIANCE_BAND` — every band and helper imported by name, no gate or `src/` edit; the tip 1 g θ ratio printed as a record; combined XDMF + setup figure; census +1, logged `missing=0`, docrefs `dead=0`. | heavy |
| `EX-156` | **The implant-heated gel ball on curved (P2-geometry) cells: the chord floor removed, the error cut ≈ 14× at the same rung** (`imp:16`; `example-runner`; §5.4 ramp from `IMP-7` step 3 / I42, `20261026T003501Z_IMP-7-step3-W1.log`; opened and queued 2026-10-26 03:00 as §9 item 650 — no example runs on second-order geometry). *History and brief: `docs/planning/chunks/EX-156.md`.* | ✅ 2026-10-26 (§9 item 650; `20261026T095740Z_EX-156.log`, Status 0, 180 s, `-n 4` complex; (cmap) 2 / 1, (S3c) 4.3e-13, (S2c) 6.2631e-05, (c-β0c) 5.9317e-02 — each equal to I42's real-build record; curved / affine whole 0.0722; census 152 → 153, figures `missing=0`, docrefs `dead=0`; setup figure rendered from the affine mesh of the same rung, disclosed). **Done-when:** `examples/implants/16_implant_heating_curved_p2.py` + same-stem guide; `./run_examples.sh -e imp:16 -n 4 -t 540` exits 0 (complex) asserting (cmap) the curved / affine coordinate degree 2 / 1, (S3c) balance ≤ `S3_BAND` and (S2c) interior ≤ `S2_BAND` on curved `LADDER[1]` at P2, and (c-β0c) the β = 0 control missing `S2_BAND` — every band and helper imported by name, no gate or `src/` edit; the volume errors and the curved / affine ratio printed; combined XDMF + setup figure; census +1, logged `missing=0`, docrefs `dead=0`. | heavy |
| `EX-157` | **A broadened two-pole (Cole–Cole) dielectric: the plane wave decays as k₀√ε̂(f), and both a Debye-reduced and a frozen-frequency medium miss** (`th:23`; `example-runner`; §5.4 ramp from `MAT-3` step 1b, §9 item 648; opened and queued 2026-10-26 03:00 as §9 item 652, serial on item 648 — `th:20` shows one Debye pole only). *History and brief: `docs/planning/chunks/EX-157.md`.* | ✅ 2026-10-26 (§9 item 652; `20261026T110247Z_EX-157.log`, Status 0, 35 s, `-n 2` complex; (Ev) 1.687e-16, (i) α / β 0.0415 % / 0.0947 % at 127.74 MHz, rate 0.9994, (ii) 0.0334 % / 0.0520 % at 63.87 MHz, (a′) 20.0464 %, (b′) 30.1096 % — each equal to `MAT-3` step 1b's W1b record; census 153 → 154, figures `missing=0`, docrefs `dead=0`). **Done-when:** `examples/time_harmonic/23_cole_cole_plane_wave.py` + same-stem guide; `./run_examples.sh -e th:23 -n 2 -t 420` exits 0 (complex) asserting (Ev), (i), (ii), (a′) and (b′) as `MAT-3` step 1b's module does, `BAND` and every helper imported by name, no gate or `src/` edit; combined XDMF + setup figure; census +1, logged `missing=0`, docrefs `dead=0`; `MAT-3` stays 🟡. | standard |
| `EX-158` | **The implant-heated gel shell switched on, continuous and duty-cycled: the CN transient on curved cells against the hollow-sphere series** (`imp:17`; `example-runner`; §5.4 ramp from `IMP-7` steps 5–6 / I46′, I47; opened and queued 2026-10-27 18:00 as §9 item 678). *History and brief: `docs/planning/chunks/EX-158.md`.* | ✅ 2026-10-28 (§9 item 678; `20261028T022737Z_EX-158.log`, Status 0, 169 s, `-n 4` complex; (T) 1.8004e-04 / 1.0444e-04, (E) 5.912e-12 / 5.867e-12, (S0) 0, (c-cont) 0.53873 — (T), (S0), (c-cont) equal to I46′'s / I47's real-build records, (E) within round-off of them (5.81e-12; audit 2026-10-28 03:00); census 154 → 155, figures `missing=0`, docrefs `dead=0`; setup figure rendered from the affine mesh of the same rung, disclosed). **Done-when:** `examples/implants/17_implant_heating_transient_duty.py` + same-stem guide; `./run_examples.sh -e imp:17 -n 4 -t 480` exits 0 (complex) asserting (T), (E), (S0) and (c-cont) at I46′'s / I47's imported bands, helpers imported by name, no gate or `src/` edit; combined XDMF + setup figure; census +1, logged `missing=0`, docrefs `dead=0`; `IMP-7` does not move. | heavy |
| `EX-159` | **The gel box with a meshed acrylic wall: order-2 convergence to the two-layer slab, and the bare box missing** (`mat:28`; `example-runner`; §5.4 ramp from `MAT-21` step 20 / P48b; opened and queued 2026-10-27 18:00 as §9 item 679). *History and brief: `docs/planning/chunks/EX-159.md`.* | ✅ 2026-10-27 22:30 slot (`20261028T033215Z_EX-159.log`). **Done-when (met):** `examples/materials/28_gel_box_meshed_acrylic_wall.py` + same-stem guide; `./run_examples.sh -e mat:28 -n 2 -t 300` exits 0 asserting (V), (S1), (E) and (c-h) at P48b's imported bands, helpers imported by name, no gate or `src/` edit; combined XDMF + setup figure; census +1, logged `missing=0`, docrefs `dead=0`; `MAT-21` does not move. | standard |
| `EX-160` | **Tune a saved four-port with no field solve: `Study.load_result` and `femem run` on a Touchstone file, and a swapped file missing** (`mri:35`; `example-runner`; §5.4 ramp from `WF-16` step 14 / U94 and `WF-1` step 24 / U99; opened and queued 2026-10-27 18:00 as §9 item 684). *History and brief: `docs/planning/chunks/EX-160.md`.* | ✅ 2026-10-27 22:30 slot (`20261028T033751Z_EX-160.log`). **Done-when (met):** `examples/mri/35_tune_a_saved_result.py` + same-stem guide; `./run_examples.sh -e mri:35 -n 2 -t 300` exits 0 (complex) asserting (CM), (S0), (=) and (c-sw) at U94's / U99's imported bands and records, no gate or `src/` edit; combined XDMF + setup figure; census +1, logged `missing=0`, docrefs `dead=0`. | standard |
| `EX-161` | **The gel box relaxing from a pre-equilibration gradient: the drift an F2182 measurement subtracts, against the triple Robin series** (`mat:29`; `example-runner`; §5.4 ramp from `MAT-21` step 19 / P47b′; opened 2026-10-27 18:00 as §9 item 683, **serial on item 676**). *History and brief: `docs/planning/chunks/EX-161.md`.* | ✅ 2026-10-27 22:30 slot (`20261028T034336Z_EX-161.log`). **Done-when (met):** `examples/materials/29_gel_box_initial_gradient_drift.py` + same-stem guide; `./run_examples.sh -e mat:29 -n 2 -t 300` exits 0 asserting (I), (E), (Sym′) and (c-Dz) at item 676's imported bands, helpers imported by name, no gate or `src/` edit; combined XDMF + setup figure; census +1, logged `missing=0`, docrefs `dead=0`. | standard |
| `EX-162` | **The walled gel box heating up: the acrylic wall's heat capacity in the transient, against the two-layer slab series** (`mat:30`; `example-runner`; §5.4 ramp from `MAT-21` step 20c / P48c; opened 2026-10-28 03:00 as §9 item 686). *History and brief: `docs/planning/chunks/EX-162.md`.* | ✅ 2026-10-28 04:30 CDT slot (`20261028T093706Z_EX-162.log`; audited PASS 10:30 review). **Done-when (met):** `examples/materials/30_gel_box_wall_heat_capacity_transient.py` + same-stem guide; `./run_examples.sh -e mat:30 -n 2 -t 480` exits 0 asserting (V), (E), (T) at t\* and t₂ and the massless-wall control (c-m) at item 682's imported bands, helpers imported by name, no gate or `src/` edit; combined XDMF + setup figure; census +1, logged `missing=0`, docrefs `dead=0`. | heavy (`-n 2`) |
| `EX-163` | **Tune a saved human-scale birdcage with no second solve: the 32-port copper F-human through `Study.load_result`** (`mri:36`; `example-runner`; §5.4 ramp from `WF-16` step 15 / U100′; opened 2026-10-28 10:30 as §9 item 690). *History and brief: `docs/planning/chunks/EX-163.md`.* | ✅ 2026-10-28 12:00 slot (`20261028T170623Z_EX-163.log`, `-n 8` complex, 238 s harness; census 159 → 160, `missing=0`, `dead=0`). **Done-when (met):** `examples/mri/36_tune_a_saved_human_scale_result.py` + same-stem guide; `./run_examples.sh -e mri:36 -n 8 -t 540` exits 0 (complex) asserting (CM), (T=), (S0) and the swapped-file control (c-sw) at U100's imported bands, helpers imported by name, no gate or `src/` edit; combined XDMF (cell tags + σ) + setup figure; census +1, logged `missing=0`, docrefs `dead=0`. | heavy (`-n 8`) |
| `EX-164` | **The walled gel box under a duty-cycled source: the acrylic wall's heat capacity when the heating switches, against the duty-cycled two-layer slab series** (`mat:31`; `example-runner`; §5.4 ramp from `MAT-21` step 21b / P49b; opened 2026-10-28 18:00 as §9 item 692). *History and brief: `docs/planning/chunks/EX-164.md`.* | ✅ 2026-10-30 15:00 CDT slot (`20261030T201705Z_EX-164.log`). **Done-when (met):** `examples/materials/31_gel_box_wall_duty_cycled_transient.py` + same-stem guide; `./run_examples.sh -e mat:31 -n 2 -t 480` exits 0 asserting (Sch), (V), (E), (S0), (T) at t\* and t₂ against `two_layer_slab_duty_cycle_temperature` and the continuous-source control (c-c) at item 688's imported bands, helpers imported by name, no gate or `src/` edit; combined XDMF (θ at t\* / t₂ and the error) + setup figure; census +1, logged `missing=0`, docrefs `dead=0`. | heavy (`-n 2`) |
| `EX-165` | **The thermal dose around the implant hole: CEM43 and the six-minute window on FEM point histories, against the hollow-sphere series** (`imp:18`; `example-runner`; §5.4 ramp from `IMP-7` steps 7–8 / I48a, I48b; opened 2026-10-31 03:00 as §9 item 704). *History and brief: `docs/planning/chunks/EX-165.md`.* | ✅ 2026-10-31 04:30 CDT slot (`20261031T094231Z_EX-165.log`). **Done-when (met):** `examples/implants/18_implant_thermal_dose.py` + same-stem guide; `./run_examples.sh -e imp:18 -n 4 -t 480` exits 0 asserting (E), (S0), (T) at 4 t\*, (Q-hd), (Win-hd), the R swap and the continuous-series control at item 700's imported bands, helpers imported by name, no gate or `src/` edit; combined XDMF (θ at 4 t\*, the series, the error) + setup figure; census +1, logged `missing=0`, docrefs `dead=0`. | heavy (`-n 4`) |
| `EX-166` | **The thermal dose in the walled gel box: CEM43 and the six-minute window on FEM point histories, against the two-layer series** (`mat:32`; `example-runner`; §5.4 ramp from `MAT-21` step 23 / P50b; opened 2026-11-01 10:30 as §9 item 732). *History and brief: `docs/planning/chunks/EX-166.md`.* | ✅ 2026-11-01 (`20261101T211827Z_EX-166.log`, Status 0, 117 s, 44 gates green; closing window `20261102T013146Z_EX-166-close.log`, `./scripts/run_examples.sh`, Status 0, 112 s, `OPS-74`) **Done-when:** `examples/materials/32_gel_box_wall_thermal_dose.py` + same-stem guide; `./run_examples.sh -e mat:32 -n 2 -t 480` exits 0 asserting (E), (S0), (T) at 4 t\*, (Q-wf), the window check, the massless-wall and continuous-source controls at item 708's imported bands, helpers imported by name, no gate or `src/` edit; combined XDMF + setup figure; census +1, logged `missing=0`, docrefs `dead=0`. | heavy (`-n 2`) |
| `EX-167` | **A rod in the coil-driven phantom through a `Study`, against the same study with the rod removed** (`mri:37`; `example-runner`; §5.4 ramp from `WF-15` step 10 / U97′ and `WF-17` step 25 / U98; opened 2026-11-01 10:30 as §9 item 733). *History and brief: `docs/planning/chunks/EX-167.md`.* | ✅ 2026-11-01 (`20261101T212457Z_EX-167.log`, Status 0, 211 s, 12 gates green; closing window `20261102T013649Z_EX-167-close.log`, `./scripts/run_examples.sh`, Status 0, 203 s, `OPS-74`) **Done-when:** `examples/mri/37_implant_rod_through_a_study.py` + same-stem guide; `./run_examples.sh -e mri:37 -n 2 -t 560` exits 0 asserting cells, S vs direct at `IDENTITY_BAND`, (i), (ii), (E) and the off `Study`'s miss at item 709's imported bands, helpers imported by name, no gate or `src/` edit; combined XDMF + setup figure; census +1, logged `missing=0`, docrefs `dead=0`. | heavy (`-n 2`) |
| `EX-168` | **The dose around the implant hole against the duty cycle: the closed-form matrix and the FEM at D = 0.25** (`imp:19`; `example-runner`; §5.4 ramp from `IMP-7` steps 10–11 / I49, I50; opened 2026-11-01 10:30 as §9 item 734). *History and brief: `docs/planning/chunks/EX-168.md`.* | ✅ 2026-11-01 (`20261101T223205Z_EX-168.log`, `-n 4` complex, flagged, Status 0, 190 s harness-wall / 186.2 s in Python; (T) 7.3744e-05, (E) 5.910e-12, (S0) 0, (c-cont) 7.7871e-01, (c-50) 5.1989e-01, R swap ≥ 6.8545e+06; 12 matrix cells and 4 FEM cells within `DOSE_BAR` / `AVG_BAR`; census 164 → 165 as predicted, figures `ok=165 missing=0 broken=0`, docrefs `dead=0 guide=0` exit 2; closing window `20261102T013342Z_EX-168-close.log`, `./scripts/run_examples.sh`, Status 0, 183 s, `OPS-74`) **Done-when:** `examples/implants/19_implant_dose_against_duty.py` + same-stem guide; `./run_examples.sh -e imp:19 -n 4 -t 480` exits 0 asserting (E), (S0), (T) at 4 t\*, (Q-hd), (Win-hd), (Q-m), (Win-m) and the R-swap, continuous and D = 0.5 controls at items 711 / 719's imported bands, helpers imported by name, no gate or `src/` edit; combined XDMF + setup figure; census +1, logged `missing=0`, docrefs `dead=0`. | heavy (`-n 4`) |
| `EX-169` | **What skipping the F2182 wait costs: the gel box's initial gradient biasing the sourced rise, against the triple Robin series** (`mat:33`; `example-runner`; §5.4 ramp from `MAT-21` steps 24 and 26 / P51, P52; opened 2026-11-01 10:30 as §9 item 735). *History and brief: `docs/planning/chunks/EX-169.md`.* | ✅ 2026-11-01 (`20261101T223957Z_EX-169.log`, `-n 2` complex, flagged, Status 0, 74 s harness / 72.0 s in-script; (I) 1.125097e-03 / 9.076702e-04, (E) ≤ 2.809e-12, (Sym′) ≤ 6.316e-14 g C, (c-ref) |Δ| 4.353e-11, (c-Dz) 6.7596e-02; census 165 → 166 as predicted, figures `ok=166 missing=0 broken=0`, docrefs `dead=0 guide=0` exit 2; closing window `20261102T013021Z_EX-169-close.log`, `./scripts/run_examples.sh`, Status 0, 74 s, `OPS-74`) **Done-when:** `examples/materials/33_gel_box_pre_equilibration_bias.py` + same-stem guide; `./run_examples.sh -e mat:33 -n 2 -t 300` exits 0 asserting (I), (E), (Sym′), the reference reproduction and the Dirichlet-z control at item 720's imported bands, helpers imported by name, no gate or `src/` edit; combined XDMF + setup figure; census +1, logged `missing=0`, docrefs `dead=0`. | standard (`-n 2`) |
| `EX-170` | **A rod in the coil-driven phantom on a two-frequency sweep, against its point solves** (`mri:38`; `example-runner`; §5.4 ramp from `WF-15` step 11 / U109; opened 2026-11-01 18:00 as §9 item 748). *History and brief: `docs/planning/chunks/EX-170.md`.* | ✅ 2026-11-02 (closing window `20261102T060640Z_EX-170.log`, `./scripts/run_examples.sh -e mri:38 -n 2 -t 560`, complex, Status 0, 256 s harness / 253.6 s in-script; cells 168 507 on all four builds; (S) 0.000e+00 at both frequencies, (i) 1.767e-15 / 6.133e-16, (ii) σ_max 0.999767876817 / 0.999152329797, column power 0.8205266 / 0.8702433, (c-f) 5.5237e-01; 13 gates, none red; census 166 → 167 as predicted, figures `ok=167 missing=0 broken=0`, docrefs `dead=0 guide=0` exit 2; flagged figure window `20261102T060200Z_EX-170-figure.log`, closes nothing) **Done-when:** `examples/mri/38_implant_rod_on_a_sweep.py` + same-stem guide; the closing window exits 0 asserting cells, sweep = point solves, (i), (ii) and the 64-vs-128 MHz control at U109's imported bands, helpers imported by name, no gate or `src/` edit; combined XDMF + setup figure; census +1, logged `missing=0`, docrefs `dead=0`. | heavy (`-n 2`) |
| `EX-171` | **What skipping the F2182 wait costs in thermal dose: FEM histories against the series** (`mat:34`; `example-runner`; §5.4 ramp from `MAT-21` steps 28–29 / P55, P56; opened 2026-11-01 18:00 as §9 item 749). *History and brief: `docs/planning/chunks/EX-171.md`.* | ✅ 2026-11-02 (closing window `20261102T061623Z_EX-171.log`, `./scripts/run_examples.sh -e mat:34 -n 2 -t 300`, complex, Status 0, 105 s harness / 102.2 s in-script; (I) 1.125097e-03 / 9.076702e-04, (E) ≤ 2.809e-12, (Sym′) ≤ 7.895e-14 g C, (Q-wf) worst 5.386e-16 on nine histories, (c-dose) 3.4217e-01, (c-ref) |Δ| 4.353e-11; census 167 → 168 as predicted, figures `ok=168 missing=0 broken=0`, docrefs `dead=0 guide=0` exit 2; flagged figure window `20261102T061416Z_EX-171-figure.log`, closes nothing) **Done-when:** `examples/materials/34_gel_box_wait_in_dose.py` + same-stem guide; the closing window exits 0 asserting (I), (E), (Sym′), (Q-wf) and the dose and reference controls at item 736's imported bands, helpers imported by name, no gate or `src/` edit; combined XDMF + setup figure; census +1, logged `missing=0`, docrefs `dead=0`. | standard (`-n 2`) |
| `EX-173` | **A dose read between two mesh planes is the interpolant's** (`mat:35`; `example-runner`; §5.4 ramp from `MAT-21` step 30 / P57; opened 2026-11-02 03:00 as §9 item 756). *History and brief: `docs/planning/chunks/EX-173.md`.* | ✅ 2026-11-02 (closing window `20261102T121038Z_EX-173.log`, `./scripts/run_examples.sh -e mat:35 -n 2 -t 300`, complex, Status 0, 71 s harness / 68.2 s in-script; (I) 1.125097e-03, (E) ≤ 2.809e-12, (Q-wf) worst 3.844e-16, (Grid) 1.542e-16 C, (c-rep) rel ≤ 3.208e-07, (c-plane) 3.5830e-01; census 169 → 170 as predicted, figures `ok=170 missing=0 broken=0`, docrefs `dead=0 guide=0` exit 2; flagged figure window `20261102T120845Z_EX-173-figure.log`, closes nothing) **Done-when:** `examples/materials/35_gel_box_dose_between_mesh_planes.py` + same-stem guide; the closing window exits 0 asserting (I), (E), (Q-wf), (Grid), (c-rep) and (c-plane) at item 742's imported bands; helpers imported by name, no gate or `src/` edit; combined XDMF + setup figure; census +1, logged `missing=0`, docrefs `dead=0`. | standard (`-n 2`) |
| `EX-174` | **The dose gap around the implant hole on a Δt ladder** (`imp:20`; `example-runner`; §5.4 ramp from `IMP-7` step 15 / I54, serial on step 14; opened 2026-11-02 03:00 as §9 item 757). *History and brief: `docs/planning/chunks/EX-174.md`.* | ✅ 2026-11-02 (closing window `20261102T123354Z_EX-174.log`, `./scripts/run_examples.sh -e imp:20 -n 4 -t 560`, complex, Status 0, 475 s harness / 472.6 s in-script; (T) 9.0012e-05 / 8.4187e-05 / 8.3973e-05, (E) ≤ 2.399e-11, (S0) 0, (Grid) 0, (c-cont) 5.3906e-01 on all three; census 170 → 171 as predicted, figures `ok=171 missing=0 broken=0`, docrefs `dead=0 guide=0` exit 2; flagged figure window `20261102T122428Z_EX-174-figure.log`, closes nothing) **Done-when:** `examples/implants/20_implant_dose_time_step_ladder.py` + same-stem guide; the closing window exits 0 asserting (T), (E), (S0), (Grid) on each of three Δt and (c-cont) at item 745's imported bands; helpers imported by name, no gate or `src/` edit; combined XDMF + setup figure; census +1, logged `missing=0`, docrefs `dead=0`. | heavy (`-n 4`) |
| `EX-175` | **Refresh the stale example artifacts, in batches, on `main`** (opened 2026-11-02 10:30 daily review; the weekly ratifies; `EX-22` and `EX-30`'s kind of work). *History and brief: `docs/planning/chunks/EX-175.md`.* | ✅ 2026-11-02, **ratified by the 11-02 weekly** — all four batches landed (§9 items 762, 763, 764, 766): 21 examples green on `main`, no red, no kill, nothing edited; `stale=168` after it (`20261102T210817Z_EX-175-b4-docrefs.log:208`). A re-run, no `EX-*` chunk re-closed; the weekly made the refresh a standing lane-1 supply (`EX-177`, then `EX-179`). | each example's own, ≤ 600 s a window |
| `EX-176` | **The same skipped wait under two schedules** (`mat:36`; `example-runner`; §5.4 ramp from `MAT-21` steps 32–33 / P59, P60; opened 2026-11-02 10:30 as §9 item 765). *History and brief: `docs/planning/chunks/EX-176.md`.* | ✅ 2026-11-02, `20261102T194429Z_EX-176.log`, Status 0, 242 s, `-n 2`; census 171 → 172. **Done-when:** `examples/materials/36_gel_box_wait_across_schedules.py` + same-stem guide; the closing window exits 0 asserting (I), (E), (S0), (Q-wf), (c-dose) and (c-cont) at items 747's and 736's imported bands; helpers imported by name, no gate or `src/` edit; combined XDMF + setup figure; census +1, logged `missing=0`, docrefs `dead=0`. | heavy (`-n 2`) |
| `EX-177` | **Refresh the stale example artifacts, round two, the `mri:` and `ports:` groups first** (opened 2026-11-02 18:00 daily review; the weekly ratifies or withdraws; `EX-175`'s kind of work). *History and brief: `docs/planning/chunks/EX-177.md`.* | ✅ 2026-11-02 22:30 (`6f9df9fb`), **ratified by the 11-02 weekly** — all six batches landed (§9 items 769–774): 35 examples green on `main`, no red, no kill, nothing edited; stale 168 → 108 (`20261103T050726Z_EX-177-b10-docrefs.log:148`). A re-run, no `EX-*` chunk re-closed; its Done-when (no script, guide, gate or band edited; a red or kill filed, nothing fixed in passing) is the standing refresh lane's, continued by `EX-179`. | each example's own, ≤ 480 s a window |
| `EX-178` | **The same watts on a rod as a line and as two tips** (`imp:21`; `example-runner`; §5.4 ramp from `IMP-8` steps 0–1; opened 2026-11-02 18:00 as §9 item 768). *History and brief: `docs/planning/chunks/EX-178.md`.* | ✅ 2026-11-02, §9 item 768: `imp:21` green on `./scripts/run_examples.sh -e imp:21 -n 2 -t 180` (`20261103T014422Z_EX-178.log`, `Status: 0`, 4 s), (Rate) 2.0047, (Far2) 1.871e-06, (c-one) 4.994e-01 at `IMP-8`'s imported bands; census 172 → 173, `missing=0`, docrefs `dead=0`. | standard (`-n 2`) |
| `EX-179` | **Refresh the stale example artifacts, round three — the standing refresh lane** | 🟡 Opened 2026-11-03 03:00 review under the 11-02 weekly's refresh-lane rule: the daily keeps one batch of ≈ 6 examples queued while `docrefs` reads stale > 0, oldest first, `mri:` / `ports:` / `mesh:` / `ans:` before the rest, ≤ 480 s each; `mag:` and the `-n 4` / `-n 8` examples outside until priced. `stale=108` at `20261103T050726Z_EX-177-b10-docrefs.log:148`; history `docs/planning/chunks/EX-179.md`. **Done-when:** per batch as `EX-177`'s — windows exit 0 with no script, guide, gate or band edited; named artifacts leave the stale list by name; a red or kill filed in known-issues, nothing fixed in passing; the chunk stays 🟡 while the lane runs and is closed or succeeded by the weekly. | standard / heavy per example `-t` |
| `EX-180` | **The skipped wait inside the walled box, and the wall's heat capacity in the drift** (`mat:37`; `example-runner`; §10 P66, §5.4 ramp from `MAT-21` step 36 / P63; opened 2026-11-04 03:00 as §9 item 813). *History and brief: `docs/planning/chunks/EX-180.md`.* | ✅ 2026-11-04, `20261104T104403Z_EX-180.log`, Status 0, 58 s, `-n 2`; census 173 → 174. **Done-when:** `examples/materials/37_gel_box_wall_pre_equilibration_bias.py` + same-stem guide; the closing window exits 0 asserting (I) on R1, (I-mw) on R3, (E), R0's item-682 control and the wall-free control at item 805's imported bands, max\|Im θ\| 0; helpers imported by name, no gate or `src/` edit; combined XDMF + setup figure; census +1, logged `missing=0`, docrefs `dead=0`. | heavy (`-n 2`) |
| `EX-181` | **A rod under a duty cycle, as a line and as two tips: the six-minute window and the dose** (`imp:22`; `example-runner`; §5.4 ramp from `IMP-8` steps 4 and 7 / I59, I62; opened 2026-11-04 10:30 as §9 item 816; closed form only, not §10 I61). *History and brief: `docs/planning/chunks/EX-181.md`.* | ✅ 2026-11-04, §9 item 816: `imp:22` green on `./scripts/run_examples.sh -e imp:22 -n 2 -t 180` (`20261104T180706Z_EX-181.log`, `Status: 0`, 5 s). (Duh) worst 2.952e-16, (Win-d) ≤ 5.720e-16, max T 35.00473 °C (tips, 0.5 W), reproduction control 3.413e-10 ≤ 1e-9, (c-cont) 6.9531e-01, all at items 787 / 811's imported bands. Census 174 → 175, `missing=0`, docrefs `dead=0`. | standard (`-n 2`) |
| `EX-182` | **The rod's heat in FEM: the line source and its two-tip twin on a coarse box** (`imp:23`; `example-runner`; §10 I61, chain I; opened 2026-11-05 03:00 as §9 item 829). *History and brief: `docs/planning/chunks/EX-182.md`.* | ✅ 2026-11-05, audit PASS (18:00 review; Done-when's `dead=0` read `dead=2` from the filed `ANS-11` SPEC entry, unchanged by this example), §9 item 829, `20261105T194909Z_EX-182.log`, Status 0, 190 s, `-n 2`; 86 726 cells (scale 1.8); (E) ≤ 2.176e-14, (S0) ≤ 4.441e-16, (I) line ≤ 8.730e-03 and tips ≤ 6.810e-03 against 1e-2; census 175 → 176, `missing=0`. **Done-when:** `examples/implants/23_rod_line_source_fem_twin.py` + same-stem guide; the closing window `./scripts/run_examples.sh -e imp:23 -n 2 -t 480` exits 0 on I55's construction at ≤ 120 000 cells (the W0 probe's count asserted), degree 2, CN N = 40, complex, asserting (E) ≤ `E_BAND` and (S0) at `B_BAND` on the line and two-tip runs and max\|Im θ\| ≤ 1e-12; (I) at I55's four points printed beside `INTERIOR_BAND`, asserted only where the flagged probe run read it inside the band; Φ₂ / Φ at (2, 15) mm printed beside 3.179482; bands and helpers imported by name, no gate or `src/` edit; combined XDMF + setup figure; census +1, logged `missing=0`, docrefs `dead=0` ⇒ ✅. Not claimed: a rate, a wall, the rod's heat capacity, entry 48's rod. | heavy (`-n 2`, `-t 480`) |
| `EX-183` | **The drift through a wall with heat capacity: the composite comparand, closed form and FEM** (`mat:38`; `example-runner`; §10 P71, chain P; opened 2026-11-05 18:00 as §9 item 848). *History and brief: `docs/planning/chunks/EX-183.md`.* | ✅ 2026-11-05 (19:30 slot; `20261106T015820Z_EX-183.log`, Status 0, 38 s; audit PASS 2026-11-06 03:00 review, guide line 55 corrected) **Done-when:** `examples/materials/38_gel_box_wall_capacity_drift.py` + same-stem guide; the closing window `./scripts/run_examples.sh -e mat:38 -n 2 -t 480` exits 0 on item 836's fixture asserting (I) on R1 and R2 against P67's composite drift, (E), (S0) on R2, the unique-z precondition, R0's item-682 control and the series-only control (c-ser) at the imported bands, max\|Im θ\| 0; the massless comparand printed beside; helpers imported by name, no gate or `src/` edit; combined XDMF + setup figure; census +1, docrefs `dead=0` (or the filed `ANS-11` entry's `dead=2` alone) ⇒ ✅. Not claimed: the six-wall box, a rate, a recommended wait, a cited PMMA ρc. | heavy (`-n 2`, `-t 480`) |
| `EX-184` | **The rod's own heat capacity and conductivity: the composite cylinder, closed form and FEM** (`imp:24`; `example-runner`; §10 I72, chain I, daily provisional under rule (6) — the 11-07 weekly ratifies; opened 2026-11-05 18:00 as §9 item 850). *History and brief: `docs/planning/chunks/EX-184.md`.* | ✅ 2026-11-05 (21:00 slot; `20261106T031102Z_EX-184.log`, Status 0, 44 s; census 178 → 179; audit PASS 2026-11-06 03:00 review) **Done-when:** `examples/implants/24_rod_composite_cylinder_heat.py` + same-stem guide; the closing window `./scripts/run_examples.sh -e imp:24 -n 2 -t 480` exits 0 on item 837's curved-P2 cylinder asserting (I) / (I0), (E), (S0), the k_r = `K_GEL` control, and item 843's (St) / (Lim) at their imported bands, max\|Im θ\| 0; the Ti / gel / hole table printed; helpers imported by name, no gate or `src/` edit; combined XDMF + setup figure; census +1, docrefs `dead=0` (or the filed `ANS-11` entry's `dead=2` alone) ⇒ ✅. Not claimed: a finite rod in gel, a Ti property, a rate. | heavy (`-n 2`, `-t 480`) |
| `EX-185` | **The implant as the heat source: the core-source composite cylinder, closed form and FEM** (`imp:25`; `example-runner`; §10 I74, chain I, daily provisional under rule (6) — the 11-07 weekly ratifies; opened 2026-11-06 03:00 as §9 item 857). *History and brief: `docs/planning/chunks/EX-185.md`.* | ✅ 2026-11-06 (04:30 slot; `20261106T105319Z_EX-185.log`, Status 0, 43 s; census 179 → 180; audit PASS, 2026-11-06 10:30 review) **Done-when:** `examples/implants/25_rod_core_source_heat.py` + same-stem guide; the closing window `./scripts/run_examples.sh -e imp:25 -n 2 -t 480` exits 0 on item 852's curved-P2 cylinder with the source in the core, asserting (I) / (I0), (E), (S0), the k_r = `K_GEL` steady-axis control, and item 847's (St) / (Lim) at their imported bands, max\|Im θ\| 0; helpers imported by name, no gate or `src/` edit; combined XDMF + setup figure; census +1, docrefs `dead=0` ⇒ ✅. Not claimed: a finite rod, a Ti property, an implant SAR shape. | heavy (`-n 2`, `-t 480`) |
| `EX-186` | **The dose under a wall with heat capacity, D = 0.25: closed form and FEM** (`mat:39`; `example-runner`; §10 P72, chain P, daily provisional under rule (6) — the 11-07 weekly ratifies; opened 2026-11-06 03:00 as §9 item 858). *History and brief: `docs/planning/chunks/EX-186.md`.* | ✅ 2026-11-06 (04:30 slot; `20261106T105958Z_EX-186.log`, Status 0, 30 s; audit PASS, 10:30 review) **Done-when:** `examples/materials/39_gel_box_wall_capacity_dose.py` + same-stem guide; the closing window `./scripts/run_examples.sh -e mat:39 -n 2 -t 480` exits 0 asserting P70's (Q-d) / (Win-d) / (D50-cf) with their controls and the massless-wall control, and P68's FEM (I) tw0, (E), (Q-wf), at the imported bands, max\|Im θ\| 0; the composite / six-face / z-only dose table printed; helpers imported by name, no gate or `src/` edit; combined XDMF + setup figure; census +1, docrefs `dead=0` ⇒ ✅. Not claimed: the six-wall box, a dose band, a recommended wait, a cited PMMA ρc. | heavy (`-n 2`, `-t 480`) |
| `EX-187` | **The SAR-heated gel sphere over time, through one `Study`** (`mri:40`; `example-runner`; §5.4 ramp from `WF-15` step 15 / U131; opened 2026-11-07 03:00 as §9 item 864). *History and brief: `docs/planning/chunks/EX-187.md`.* | ✅ 2026-11-07 (§9 item 864; closing log `20261107T181301Z_EX-187.log`). `mri:40` is green: (A_t) 3.4708e-02 and (S2) 4.9994e-05 equal `MAT-12`'s prints, (c-qs) 1.0108e+00 misses the band, census 181 → 182. **Done-when:** `examples/mri/40_study_sar_driven_transient_sphere.py` + same-stem guide; the closing window `./scripts/run_examples.sh -e mri:40 -n 2 -t 180` exits 0 asserting U131's (UT), (I-dir), (A_t) with `MAT-12`'s `.4e` record, (S2) with its record, (E) and the (c-qs) control at the imported bands, through `Study.temperature(…, transient=)` at degree 1; helpers imported by name, no gate or `src/` edit; combined XDMF + setup figure; census +1, docrefs `dead=0` ⇒ ✅. Not claimed: degree 2, duty, probes or a YAML key on the sphere. | standard (`-n 2`, `-t 180`) |
| `EX-188` | **The gel / Ti / hole bracket on one page: the rod as gel, as titanium and as a hole, heated at its tips and through its own body** (`imp:26`; `example-runner`; §5.4 ramp from `IMP-8` I66 / I69 / I73; §10 I77; opened 2026-11-08 03:00 as §9 item 882). *History and brief: `docs/planning/chunks/EX-188.md`.* | ✅ 2026-11-08 (§9 item 882; closing log `20261108T211230Z_EX-188.log`; audit PASS 18:00). `imp:26` is green on 64 140 cells: all five (E) ≤ 1.71e-14, (S0) ≤ 4.441e-16, and the three controls move by +0.50244 / −0.078325 / −0.030925 (tips / Ti, tips / hole, rod / Ti); census 182 → 183. **Done-when:** `examples/implants/26_rod_bracket_gel_ti_hole.py` + same-stem guide; the closing window `./scripts/run_examples.sh -e imp:26 -n 2 -t 480` exits 0, on one coarse box of ≤ 120 000 cells (the W0 figure), asserting (E) ≤ `E_BAND` and (S0) within `B_BAND` on all five CN runs (tip source on gel / Ti / hole rods, rod source on gel / Ti), max\|Im θ\| ≤ 1e-12, and the control that each non-gel arm's θ(2, 0) mm is off its gel run's by more than `CONTROL_FLOOR`; the table of θ(2, 0), θ(2, 15) mm and the tip-point ratio per arm printed beside the 485 291-cell records, (I) on the gel runs printed with `INTERIOR_BAND` beside; helpers imported by name, no gate or `src/` edit; combined XDMF + setup figure; census +1, docrefs `dead=0` ⇒ ✅. Not claimed: the gates' mesh, a rate, a Ti property, the F2182 probe (I75), an implant SAR shape. | heavy label (`-n 2`, `-t 480`) |
| `EX-189` | **The dose vs the wall's heat capacity: the D = 0.25 dose across the wait, closed form and FEM** (`mat:40`; `example-runner`; §10 P77, chain P; opened 2026-11-08 03:00 as §9 item 886, serial on P74). *History and brief: `docs/planning/chunks/EX-189.md`.* | ✅ 2026-11-08 (§9 item 886; closing log `20261109T030424Z_EX-189.log`, 85 gates, 0 red, 31 s; audit PASS 2026-11-09 03:00). (I) tw1 4.712673e-03 ≤ `INTERIOR_BAND`, (E) ≤ 1.939e-13 on all three calls, hand-off 0 K; composite / massless − 1 at P5 falls with the wait, +3.7292e-02 → +1.0069e-02; census 183 → 184. **Done-when:** `examples/materials/40_gel_box_wall_capacity_dose_vs_wait.py` + same-stem guide; the closing window `./scripts/run_examples.sh -e mat:40 -n 2 -t 480` exits 0 asserting P73's (Lim) against P70's construction, (Mono), (Q-d) / (Win-d) and the massless-wall control at the four waits, and P74's FEM (I) on the t_w = 1/λ₁ run, (E) on both of its calls, the hand-off precondition and (Q-wf), at the imported bands, max\|Im θ\| 0; the composite / massless dose-by-wait table and the FEM dose at P5 beside P73's printed; helpers imported by name, no gate or `src/` edit; combined XDMF + setup figure; census +1, docrefs `dead=0` ⇒ ✅. Not claimed: the six-wall box, a dose band, a recommended wait or thickness, a cited PMMA ρc. | heavy (`-n 2`, `-t 480`) |
| `EX-190` | **The F2182-shaped probe on the rod bracket: gel, titanium and hole rods read 1 and 2 mm beyond the tip, continuous and at D = 0.5** (`imp:27`; `example-runner`; §5.4 ramp from `IMP-8` I75 / I78; §10 I82, chain I; opened 2026-11-08 18:00 as §9 item 896, rule (6) provisional). *History and brief: `docs/planning/chunks/EX-190.md`.* | ✅ 2026-11-09 (04:30), §9 item 896: `20261109T105242Z_EX-190.log` exit 0, 152 s; audit PASS 2026-11-09 10:30. **Done-when:** `examples/implants/27_rod_f2182_probe_bracket.py` + same-stem guide; the closing window `./scripts/run_examples.sh -e imp:27 -n 2 -t 480` exits 0 on `imp:26`'s 64 140-cell box, asserting (E) ≤ `E_BAND` and (S0) within `B_BAND` on three rod-source CN runs (gel / Ti / hole), max\|Im θ\| ≤ 1e-12, `imp:26`'s rod / gel and rod / Ti θ(2, 0), θ(2, 15) reproduced within one sixth-digit unit, each non-gel arm off gel by more than `CONTROL_FLOOR`, θ_probe(900 s) > θ_probe(675 s) at (0, ±16) mm, and 1 − θ(17)/θ(16) > `INTERIOR_BAND` on every arm; the probe table at 16 / 17 mm, continuous and D = 0.5 by superposition, printed beside I75's records and I78's closed form; helpers imported by name, no gate or `src/` edit; combined XDMF + setup figure; census +1, docrefs `dead=0` ⇒ ✅. Not claimed: the gates' mesh, the clinical fixture, a switched duty solve, any F2182 or published value. | heavy label (`-n 2`, `-t 480`) |
| `EX-191` | **The rod's finite radius at the probe: gel and Ti beside the uniform-cylinder closed form** (`imp:28`; `example-runner`; §5.4 ramp from `IMP-8` I83; §10 I84, chain I; opened 2026-11-09 10:30 as §9 item 920, rule (6) provisional). *History and brief: `docs/planning/chunks/EX-191.md`.* | ✅ 2026-11-09 (12:00 slot, fourth item), §9 item 920: run of record `20261109T183047Z_EX-191-r2.log` exit 0, 104 s harness (102.0 s in-script), re-run after a figure-window digit was struck from the guide (first closing log `20261109T182711Z_EX-191.log`, same gate digits); census 185 → 186; audit PASS 2026-11-09 18:00. **Done-when:** `examples/implants/28_rod_finite_radius_probe.py` + same-stem guide; the closing window `./scripts/run_examples.sh -e imp:28 -n 2 -t 480` exits 0 on `imp:27`'s 64 140-cell box, asserting (E) ≤ `E_BAND` and (S0) within `B_BAND` on the gel and Ti rod-source CN runs, max\|Im θ\| ≤ 1e-12, `imp:27`'s gel / Ti θ(900 s) at ±16 / ±17 mm and D = 0.5 means reproduced within one sixth-digit unit, (Lim) at `B1_BAND` at every rung, 1 − θ̄_gel/(SCALE·Φ_line) and \|θ̄_Ti/(SCALE·Φ_a) − 1\| > `INTERIOR_BAND` at 16 / 17 mm, and θ(16) > θ(17) > θ(18) > θ(20); I83's (G) printed beside `INTERIOR_BAND`, not asserted; helpers imported by name, no gate or `src/` edit; combined XDMF + setup figure; census +1, docrefs `dead=0` ⇒ ✅. Not claimed: an h-statement or rate, the conducting rod's closed form, the clinical fixture, any F2182 or published value. | heavy label (`-n 2`, `-t 480`) |
| `EX-192` | **Entry 59's readout points under three closed-form SAR pictures: uniform, eddy and end-weighted per deposited watt, beside entry 52's SAR10g shape** (`mat:41`; `example-runner`; §10 P82, chain P; opened 2026-11-09 10:30 as §9 item 922, rule (6) provisional, serial on item 921). *History and brief: `docs/planning/chunks/EX-192.md`.* | ✅ 2026-11-09 (§9 item 922; closing log `20261109T193426Z_EX-192.log`, 29 gates, 0 red, 4.9 s in-script; audit PASS 2026-11-09 18:00; census +1, 186 → 187 on the docrefs basis). P80's and P81's gates re-asserted on the page, (V-pt) worst 1.226e-14 (eddy; U3 1.124e-15, end 1.297e-15) ≤ `B1_BAND`, (cV) least miss 4.141e-03 > `B1_BAND`; max\|Im θ\| 0; census 180 → 181, docrefs `dead=0 guide=0`. **Done-when (as opened):** `examples/materials/41_gel_box_readout_points_three_sar_pictures.py` + same-stem guide; the closing window `./scripts/run_examples.sh -e mat:41 -n 2 -t 480` exits 0, asserting P80's and P81's (P-eq), (Q2-p) / (Q2-s), (Rep-s), (Lim), (Lin), (Ball), (K-lat) and the controls (c-u), (c-e), (c-ball) at the imported bands, plus (V-pt) at the 24 vertices and the (cV) separation; max\|Im θ\| 0; helpers imported by name, no gate or `src/` edit; combined XDMF + setup figure; census +1, docrefs `dead=0` ⇒ ✅. Not claimed: which picture entry 59 follows, a fitted source, any entry-59 / entry-60 figure, a dose band. | heavy (`-n 2`, `-t 480`) |
| `EX-193` | **The field-strength lever on one page: F-clinical at 1.5 T and 3 T per accepted watt, factored into the box's power share and SAR per box watt** (`mri:41`; `example-runner`; §5.4 ramp from `GEO-38` C33 / `xl` entries 52, 58; §10 C47, chain C; opened 2026-11-09 10:30 as §9 item 924, rule (6) provisional, serial on the 11-09 weekly's reading of entry 58). *History and brief: `docs/planning/chunks/EX-193.md`.* | ✅ 2026-11-09 (§9 item 924, 22:30 CST slot; audit PASS 2026-11-10 03:00; closing log `20261110T043404Z_EX-193.log`, Status 0, 4 s, 35 of 35 gates green; census +1, 187 → 188 on the docrefs `runnable` basis, `dead=0 guide=0`). (R33) four ratio strings and both c\* labels equal C33's; (Dec-S) rel 9.601e-08 ≤ 8.450e-07, (Dec-B) rel 3.249e-07 ≤ 1.177e-06; controls (N-norm), (N-plant), (N-swap) red as posed. **Done-when (as opened):** `examples/mri/41_f_clinical_field_strength_lever.py` + same-stem guide; closing window `./scripts/run_examples.sh -e mri:41 -n 2 -t 120` exits 0 with (X), (D), (R33) C33's four ratios and c\* labels reproduced by string, (Dec-S) / (Dec-B) the factorisations at half a unit, (Lat), (Arg), and controls (N-norm), (N-plant), (N-swap) asserted; imported helpers only, no gate-module / page / `src/` edit; setup figure a byte-identical copy of `mesh:16`'s render; census +1, docrefs `dead=0 guide=0` ⇒ ✅. | standard |
| `EX-194` | **The field-strength lever on the heating: P83's closed-form prediction beside C39's reading of `xl` entries 59 / 60, per accepted watt** (`mat:42`; `example-runner`; §5.4 ramp from `GEO-38` C39 / `MAT-21` step 46; §10 P84, chain P; opened 2026-11-11 03:00 as §9 item 954, rule (4), serial on item 953). *History and brief: `docs/planning/chunks/EX-194.md`.* | ✅ 2026-11-11 (§9 item 954, 04:30 slot; closing log `20261111T104011Z_EX-194.log`, Status 0, 2 s, 78 of 78 gates green; census +1, 189 → 190 on the docrefs `runnable` basis, `dead=0 guide=0`, figures `broken=0`). (R39) all 22 ratio strings equal C39's; (Share) rel 1.048e-09 ≤ 1.389e-08; the ccw hotspot ratio 4.042515 is 1.094495× P83's power share and 0.985633× the SAR lever. | smoke |
| `EX-195` | **The finite-radius two-tip twin at the probe: I86's closed form beside the point twin, the line and Φ_a on the 16 / 17 / 18 / 20 mm ladder** (`imp:29`; `example-runner`; §5.4 ramp from `IMP-8` I86; §10 I89, chain I; opened 2026-11-12 03:00 as §9 item 964, rule (4)). *History and brief: `docs/planning/chunks/EX-195.md`.* | ✅ *(**DONE 2026-11-12 (07:30 slot, §9 item 964), first try:** `20261112T133229Z_EX-195.log` Status 0, 3 s; every gate green — see item 964; census 191 → 192, docrefs `dead=0 guide=0`, figures `broken=0`.)* Opened 2026-11-12 03:00 as §9 item 964: a numpy record page, no mesh or solve, reproducing I86's logged strings (`20261110T120247Z_IMP-8-item934-W1.log`). **Done-when:** `examples/implants/29_rod_tip_twin_finite_radius_probe.py` + same-stem guide; the closing window `./scripts/run_examples.sh -e imp:29 -n 1 -t 120` exits 0 asserting (X) custody of I86's tracked log, (K) θ(900 s) per watt 60.34935 at 16 mm, the point twin 63.14, the radius factor 0.9558023 and the tip weighting 6.561211 by string, (K-rec) the rest of I86's record block by string, (Lim) / (Lim-a) at `B1_BAND`, the control 1 − Φ₂,a/Φ₂ > `INTERIOR_BAND` at 16 mm; helpers imported by name, comparands parsed never typed, no gate or `src/` edit; census +1 on the docrefs `runnable` basis, `dead=0 guide=0` ⇒ ✅. Not claimed: FEM, box walls, rod conduction, any F2182 or published value. | standard (`-n 1`, `-t 120`) |
| `EX-172` | **Sweep once, tune at any of its frequencies, with no second solve** (`mri:39`; `example-runner`; §5.4 ramp from `WF-16` step 17 / U108 and `WF-1` step 27 / U111; opened 2026-11-01 18:00 as §9 item 750). *History and brief: `docs/planning/chunks/EX-172.md`.* | ✅ 2026-11-02 (closing window `20261102T062422Z_EX-172.log`, `./scripts/run_examples.sh -e mri:39 -n 2 -t 300`, complex, Status 0, 55 s harness / 52.5 s in-script; (C=) rel 4.985e-14 / 5.218e-15 / 1.912e-13 at `REDUCTION_BAND`, (I=) API = YAML by `float.hex` at all three, (f=) exact, (S0) 0 calls, (c-next) 0.975063 / 0.754562, (Ref) `ValueError` at 1.92e8 Hz; census 168 → 169 as predicted, figures `ok=169 missing=0 broken=0`, docrefs `dead=0 guide=0` exit 2; flagged figure window `20261102T062150Z_EX-172-figure.log`, closes nothing) **Done-when:** `examples/mri/39_tune_any_frequency_of_a_sweep.py` + same-stem guide; the closing window exits 0 asserting (C=), (I=), (f=), (S0), the next-frequency control and the out-of-span refusal at items 722 / 730's imported bands, on this run's own file; helpers imported by name, no gate or `src/` edit; combined XDMF + setup figure; census +1, logged `missing=0`, docrefs `dead=0`. | standard (`-n 2`) |
| `EX-119` | **The F2182 gel box heated by an end-weighted source — UFL, DG0 and P2 hand-offs in ParaView** (`mat:20`; `example-runner`; opened 2026-10-16 03:00 daily, §5.4 ramp from `MAT-21` P32 / P34 / P35, §9 item 493; `mat:19` is uniform-q only) *(history: `docs/planning/chunks/EX-119.md`)* **Done-when:** `examples/materials/20_f2182_gel_box_end_weighted_source.py` + same-stem guide (three required headings + `## Setup figure` + `## Related`, PNG ≤ 600 KiB); complex, `-n 2`, m = 4 at P2 through the gate modules' helpers; no `tests/` or `src/` edit; every band imported. Asserted: (Pdeg); (S1-sep) / (S1-dg0) / (S1-p2s) / (T-p2s) ≤ `INTERIOR_BAND`; (S2-p2s) ≤ `S2_GAUSS_BAND`; (E-p2s) ≤ `E_BAND`; (F0) face partition; (F1) Σ face powers = source power ≤ `S2_GAUSS_BAND`; controls (c-p2s), (cF). Records printed beside `20261016T033407Z_MAT-21-step8-W2-p2-source.log`. Combined XDMF (P1 q / θ fields, DG0 q, face tags) + x-centre-line CSV. Census +1, `broken=0`, docrefs `dead=0`, elapsed recorded ⇒ ✅. Not claimed: coil SAR, the DG0 staircase series, an h-rate, compliance, an external check. | ✅ 2026-10-16 (§9 item 493, 06:00 slot): `examples/materials/20_f2182_gel_box_end_weighted_source.py` + guide + setup PNG; flagged `20261016T112833Z_EX-119.log` (Status 0, 107 s) and control `20261016T113027Z_EX-119-control.log` (Status 0, 103 s), `-n 2`: (Pdeg) 116 025 ×4; (S1-p2s) 3.486949e-05, (T-p2s) 1.388286e-04, (S2-p2s) 1.677e-14, (E-p2s) 1.079e-13 at the gate's digits; (F0) 9 344 exact, (F1) ≤ 1.831e-13; (c-p2s) 0.51239, (cF) 0.06973; census figures 115 → 116 `broken=0`, docrefs `dead=0` exit 2. | heavy (`-n 2`, ≈ 190 s per run, `-t 590`) |
| `EX-79` | *(history and §9 briefs: `docs/planning/chunks/EX-79.md`)* **Coil-driven SAR heats the gel phantom: four drives, four hot spots, one C4 rotation** (`mat:9`; `example-runner`; §5.4 ramp for Phase 8 — opened 2026-10-02 03:00 review, §10 P3, on `MAT-15` step 1 (§9 item 130); `mat:4` / `mat:7` heat a ball under an imposed field, `mri:3` shows the drives' SAR without heat). **Done-when:** `examples/materials/09_coil_driven_heating_phantom.py` + same-stem guide (three required headings + `## Setup figure`, `materials_09_coil_driven_heating_phantom_setup.png`); it calls `_run_heated(return_fields=True)` — rule (a) lift of step 1's fixture body, additive, module re-gated green in-slot, disclosed — and imports its bands (`is` checks), never restating them. On `ANS-2`'s 0.0025 m rung, gel constants, G = 1 W/P₁, complex, `-n 2`: (B) ≤ `B_BAND`, V_sub ≤ `V_BAND`, (E)/(E′) ≤ `E_BAND`, (S2) ≤ `S2_BAND`, the four θ 10 g C4 pairs ≤ `C4_COVARIANCE_BAND`, each printed beside `20261001T123309Z_MAT-15-step1-W1.log:1924–1954`. Controls (asserted, backed `:1949–1953`): (d) doubled source > `D_CONTROL_FLOOR`; (f) far-side ball > `F_CONTROL_FLOOR`. `materials_09_coil_driven_heating_phantom_combined.xdmf` (P1 θ at t\*, DG0 q_h, four drives). Census +1 from the fresh pre-census, `broken=0`; docrefs `exit != 1`; elapsed recorded. Not claimed: absolute temperature, quadrature drive, convergence, perfusion, tissue, `ANS-9`. | ✅ **Audit PASS 2026-10-02 10:30 review (chunk file).** **2026-10-02 (§9 item 172).** Flagged `-n 2` run green in 162 s (`20261002T114023Z_EX-79.log`): (B) 2.2e-14, V_sub 3.4e-13, (E) 8.8e-12, (E′) 2.3e-14, (S2) 5.75e-6, C4 0.020–0.049 %; controls (d) 1.000000 and (f) 75.25–75.28 %; the five `is` checks True. Re-gate 35 passed with W1's digits (`…112729Z_EX-79-mat15-regate.log`); census 72 → 73, `broken=0`. | heavy (`-n 2`, ≈ 220 s, `-t 480`) |
| `EX-80` | *(history and §9 briefs: `docs/planning/chunks/EX-80.md`)* **The human-scale coil on copper: the F-human 32×32 through the Leontovich wall, and where each drive's absorbed power goes** (`mri:6`; `example-runner`; opened 2026-10-02 10:30 daily review, §5.4 ramp from `WF-12` step 2 (§9 item 169), deferred to H12b by the 03:00 review; `ports:18` / `mri:5` map copper loss on F-small and `ans:8` / `mri:4` solve F-human at σ = 800, but no example solves F-human on copper) | ✅ **Audit PASS 2026-10-02 18:00 review. DONE 2026-10-02 (15:00 CDT slot, §9 item 179):** re-gate `20261002T202528Z_EX-80-wf12-regate.log` 13 passed, 311.9 s, W2 digits; flagged `20261002T203221Z_EX-80.log` 278 s (`-n 8`), unflagged `20261002T203707Z_EX-80-control.log` 270 s, records unchanged to round-off; census figures 75 → 76 / ok 76 / broken 0, docrefs 75 → 76 runnable, exit 2, stale 6 unchanged. Total ≈ 17 min of compute. **Done-when:** `examples/mri/06_f_human_copper_loss_partition.py` + same-stem guide (three required headings + `## Setup figure`, PNG ≤ 600 KiB). It calls `_run_step2(32, True, field_hook=…)`, the rule (a) lift of the gate test's body (additive, module re-gated green in-slot with W2's digits, disclosed), and imports every band and record. At `-n 8`, degree 1, complex: (A0), (R), (P) ×32, (ii), (iii), (iv′) and (PEC) at the gate's imported bands, each printed beside `20261002T094751Z_WF-12-step2-W2.log:6352–6414`. Asserted control: the P_surf-blind miss ≥ `BLIND_FACTOR`× on every drive. `mri_06_f_human_copper_loss_partition_combined.xdmf` (DG0 \|E\|, phantom loss density, Cu coil loss per cell; Cu and PEC P17), whose map sums reproduce P_phantom / P_surf ≤ `DISCRETE_IDENTITY_RTOL`. Census +1, `broken=0`; docrefs `exit != 1`. Not claimed: \|B₁⁺\| (H9), `C_mode` (H12c), degree 2, AED, absolute loss. | heavy (`-n 8`, ≈ 300 s, `-t 480`) |
| `EX-81` | *(history and §9 briefs: `docs/planning/chunks/EX-81.md`)* **The phantom's own cylinder against its closed form: uniform heating with Robin walls beside the double eigen-series, steady and at t\*** (`mat:10`; `example-runner`; opened 2026-10-02 10:30 daily review, §5.4 ramp from `MAT-15` step 2(b) (§9 item 168); `mat:4` / `6` / `7` / `8` compare balls and `mat:9` heats this cylinder with no comparand) | ✅ **Audit PASS 2026-10-02 18:00 review.** 2026-10-02 (16:30 CDT slot; §9 item 180; outcome in `docs/planning/chunks/EX-81.md`). Done-when met: `examples/materials/10_finite_cylinder_robin_series.py` + same-stem guide (three required headings + `## Setup figure`). It calls `_run_transient(return_fields=True)`, the rule (a) lift of the `transient` fixture body (additive, transient nodes re-gated green in-slot, disclosed), and imports every band (`is` checks). On the 0.0025 m rung, `-n 2`, complex: cells = `FINEST_CELLS`; (T), (E) and (∞) at `INTERIOR_BAND` / `E_BAND` / `INF_BAND`; steady interior ≤ `INTERIOR_BAND` and Gauss ≤ `S2_GAUSS_BAND`, composed from the gate's `_exact` / `_rel`. Each is printed beside `…093610Z_…W5-transient.log:1875–1883` / `…093152Z_…W4-steady.log:7174–7176`. Asserted control: (c2) Dirichlet roots miss (record 0.77946). `materials_10_finite_cylinder_robin_series_combined.xdmf` (P1 θ_h, series, error, steady and t\*; DG0 interior mask). Census +1, `broken=0`; docrefs `exit != 1`. Not claimed: the LSQ rate, EM, perfusion, `ANS-9`. | heavy (`-n 2`, ≈ 270 s, `-t 480`) |
| `EX-82` | *(history and §9 briefs: `docs/planning/chunks/EX-82.md`)* **A cracked cavity: each lip of the split wall is its own Leontovich wall, Q drops to Pozar's 0.772 × Q₀, and a single-valued sheet in the same plane does nothing** (`th:15`; `example-runner`; opened 2026-10-02 10:30 daily review, §5.4 ramp from `GEO-33` step 0b (§9 item 170); `th:10` shows the uncracked box only) | ✅ **Audit PASS 2026-10-02 18:00 review. DONE 2026-10-02 (15:00 CDT slot, §9 item 181):** every anchor and the control green first run; flagged `20261002T202252Z_EX-82.log` 15 s (`-n 2`), unflagged `20261002T202332Z_EX-82-control.log` 8 s, records unchanged; census figures 74 → 75 / ok 75 / broken 0, docrefs 74 → 75 runnable, exit 2 (stale 6, unchanged). *Was:* ⬜ Queued §9 item 181. **Done-when:** `examples/time_harmonic/15_cracked_box_two_sided_lips.py` + same-stem guide (three required headings + `## Setup figure`). It uses `cracked_box_mesh` and `solve_impedance_mesh_mode(…, return_vectors=True)`, a rule (c) additive keyword with the `TH-14` and `GEO-33` modules re-gated green in-slot and disclosed, and imports every band and literal from `test_geo33_two_sided_leontovich`. At `-n 2`, degree 2: (Q0), (Qs) both pair members and (R) within `Q_TOLERANCE`; (F) and (D); each printed beside `20261002T095638Z_GEO-33-step0b-W2.log:91–101`; each exported vector's Rayleigh quotient = λ to 1e-6. Asserted control: the septum sheet reads Q_sheet/Q₀ within `Q_TOLERANCE` of 1, with its miss at the 29.548 % ceiling. `time_harmonic_15_cracked_box_two_sided_lips_combined.xdmf` (DG0 \|E\| of both pair members, \|curl E\|, half index; lip facet tags). Census +1, `broken=0`; docrefs `exit != 1`. Not claimed: a partial crack, strip sheets, a birdcage, copper, `GEO-33` step 1. | standard (`-n 2`, ≈ 40 s, `-t 180`) |
| `EX-83` | *(history and §9 briefs: `docs/planning/chunks/EX-83.md`)* **A metal rod in the coil-driven phantom: the 64 MHz birdcage's 1 g tip SAR with and without it, beside `IMP-2`'s closed-form estimate, and the \|E\| / SAR maps at the tip** (`imp:3`; `example-runner`; opened 2026-10-02 18:00 daily review, §5.4 ramp from `IMP-3` (§9 item 174), carrying the `GEO-36` example the 10:30 review deferred; `imp:1` / `imp:2` impose a quasi-static field with no coil or rod, `mri:3` has no implant) | ✅ **Audit PASS 2026-10-03 03:00 review (chunk file).** 2026-10-03 (§9 item 191): every anchor holds at its imported band on the 1 g rung at `-n 4` (223 s), the Q > 1 control holds on all 12 tip × drive pairs, and the census went 77 → 78 with `broken=0` (`20261003T051229Z_EX-83.log`); the re-gates are digit-identical to W1 / W2. **Done-when:** `examples/implants/03_rod_in_birdcage_tip_sar.py` + same-stem guide (three required headings + `## Setup figure`, PNG ≤ 600 KiB). It calls `_run_one_gram_rung(field_hook=…)`, the rule (a) lift of `test_imp3_one_gram_rung`'s body with an additive `field_hook=None` through `_case` / `_arm` (both nodes re-gated green in-slot, disclosed), and imports every band and record. On the 1 g rung at `-n 4`, degree 1, complex: (A0), (i), (ii), (P), (M), (E), (S) at the gate's imported bands, each printed beside `20261002T170843Z_IMP-3-step1-W2.log:3798–3877`; Q and Q / estimate printed, never asserted. Asserted control: the rod-blind arm under-reads the 1 g tip SAR on all six drives at both tips (Q > 1; record min at W2:3867). `implants_03_rod_in_birdcage_tip_sar_combined.xdmf` on the rod mesh (DG0 \|E\| and local SAR for P1 and ccw, 1 g ball masks, tag 402), whose P1 loss-map sum reproduces P_phantom ≤ `DISCRETE_IDENTITY_RTOL`. Census +1, `broken=0`; docrefs `exit != 1`. Not claimed: the ratio's accuracy (`ANS-10`), absolute SAR, heating, C1 on the rung, F2182 / 10974. | heavy (`-n 4`, ≈ 280 s, `-t 480`) |
| `EX-84` | `imp:4` — one adjoint solve gives an implant lead's transfer function | ✅ Audit PASS 2026-10-03: `imp:4` holds every anchor and (X0)–(X3) at its imported band at `-n 2`, census 78 → 79. **Done-when:** `examples/implants/04_lead_transfer_function_adjoint.py` + same-stem guide (three required headings + `## Setup figure`, PNG ≤ 600 KiB). *History: `docs/planning/chunks/EX-84.md`.* | heavy (`-n 2`, ≈ 150 s, `-t 420`) |
| `EX-85` | *(history and §9 briefs: `docs/planning/chunks/EX-85.md`)* **Duty-cycled RF heating of the gel ball: the on/off history against the closed form, the periodic peak against the time-averaged-power steady state, and what a continuous source at the same average power hides** (`mat:11`; `example-runner`; opened 2026-10-03 18:00 daily review, §5.4 ramp from `MAT-17` ✅ (§9 item 200, `f320a1f6`, audited PASS); `mat:4` / `6` / `7` / `8` / `10` heat with an always-on source and no example imports `source_schedule`) | ✅ 2026-10-04 (§9 item 207, `20261004T004009Z_EX-85.log`, 317 s, complex build): the gate's (i) reproduced digit-identical (4.3870e-04 / 1.0683e-03), (H) / (AVG) inside `INTERIOR_BAND`, both controls miss; census 79 → 80, `broken=0`. Audited PASS 2026-10-04 03:00 review. **Done-when:** `examples/materials/11_duty_cycled_gel_ball.py` + same-stem guide (three required headings + `## Setup figure`, PNG ≤ 600 KiB), importing every band and helper from `test_mat17_duty_cycle_fem.py`, its step-1 module and `MAT-13`'s complex-safe twins; no gate or `src/` edit. On `MAT-9`'s 0.005 m ball, CN + Rannacher, `-n 2`, complex: (S0) on every scheduled run; (iii) ≤ `BALANCE_BAND`; the gate's (i) for cases A / B ≤ `INTERIOR_BAND`, printed beside `20261003T171113Z_MAT-17-step2-W1b.log:98, :101`; two mid-phase instants (H) and the average-power anchor (AVG) ≤ `INTERIOR_BAND`. Asserted controls: schedule ignored misses; the case-B average-power stand-in misses (factor printed). Combined XDMF (history + t_end). Census +1, `broken=0`; docrefs `exit != 1`. Not claimed: the FEM periodic state, SAR / coil drive, perfusion, any regulatory statement. | heavy (`-n 2`, ≈ 80–240 s, `-t 420`) |
| `EX-86` | *(history and §9 briefs: `docs/planning/chunks/EX-86.md`)* **The gel ball in a container wall: the temperature kinks at r = a while the flux stays continuous, and the wall raises the core and slows the settling** (`mat:12`; `example-runner`; opened 2026-10-03 18:00 daily review, §5.4 ramp from `MAT-18` ✅ `c371c5a4`, audited PASS; every thermal example solves one homogeneous medium) | ✅ 2026-10-04 (§9 item 208, `20261004T010436Z_EX-86-retry2.log`, 98 s, complex build; two earlier windows red on example-script bugs, disclosed in the chunk file): (S) 1.1206e-03, (T1) 8.2066e-04 and both controls digit-identical to `MAT-18`'s records; census 80 → 81, `broken=0`. Audited PASS 2026-10-04 03:00 review. **Done-when:** `examples/materials/12_gel_ball_in_container_wall.py` + same-stem guide (three required headings + `## Setup figure`, PNG ≤ 600 KiB); imports `MAT-18`'s fixture, helpers and comparands and `MAT-13`'s `INTERIOR_BAND` / `BALANCE_BAND`, restates none; no gate or `src/` edit. At `-n 2`, degree 1, complex, 0.005 rung: steady interior ≤ `INTERIOR_BAND` and balance ≤ `BALANCE_BAND`; CN + Rannacher N 40 at t\* interior ≤ `INTERIOR_BAND`; every run ≤ `BALANCE_BAND`; \|Im T\| = 0; each printed beside `…172548Z_MAT-18-step2a-W1.log` / `…184850Z_MAT-18-step2b-W1b.log`. Asserted controls: uniform k and uniform ρc miss `INTERIOR_BAND`. Printed: the profile across r = a, the core / wall / film budget, t₉₉ against the bare ball. `materials_12_gel_ball_in_container_wall_combined.xdmf`. Census +1, `broken=0`; docrefs `exit != 1`. Not claimed: contact resistance, a real container, degree 2, SAR / EM, perfusion. | heavy (`-n 2`, ≈ 150 s, `-t 420`) |
| `EX-87` | *(history and §9 briefs: `docs/planning/chunks/EX-87.md`)* **A lead shaft in the gel: the wire shorts the field that runs along it and doubles the field that crosses it** (`imp:5`; `example-runner`; opened 2026-10-03 18:00 daily review, §5.4 ramp from `IMP-4` ✅ `9b8fefc7`, audited PASS; no example builds a cylinder, evaluates a Hankel-function field or drives one conductor in both polarisations) | ✅ 2026-10-04 (§9 item 209, `20261004T033913Z_EX-87.log`, 146 s, complex build; one earlier window red on an example-script XDMF bug, disclosed in the chunk file): TM / TE β, Φ_h / Φ_cf,h and cavity readings digit-identical to `IMP-4` W1 `:133–142`; rule (a) re-gate 15 passed (`20261004T033046Z_EX-87-imp4-regate.log`); C1 / C2 miss; census 81 → 82, `broken=0`. Audited PASS 2026-10-04 03:00 review. **Done-when:** `examples/implants/05_lead_shaft_polarisation_in_gel.py` + same-stem guide (three required headings + `## Setup figure`, PNG ≤ 600 KiB); imports `IMP-4`'s fixture, helpers and bands, restates none; rule (a) keyword-only `return_fields` / `trace` on `_solve_cylinder`, module re-gated green in-slot, disclosed. At `-n 2`, degree 2, complex, one mesh, per drive (TM_z n = 0, TE_z n = 1): \|β − 1\| ≤ `TH15_BETA_BAND`, annulus \|Φ_h/Φ_cf,h − 1\| ≤ `MAT4_BAND`, cavity dofs > 0 at max \|E\| < 1e-12, printed beside `20261003T200454Z_IMP-4-step1-W1.log:133–142`; DG0 FEM and closed-form SAR-ratio maps sum over the annulus to Φ_h / Φ_cf,h ≤ 1e-12. Asserted controls (wired), the swapped end condition; factors printed. `implants_05_lead_shaft_polarisation_in_gel_combined.xdmf`. Census +1, `broken=0`; docrefs `exit != 1`. Not claimed: a finite, insulated or coated lead, a tip, a coil, heating, convergence. | heavy (`-n 2`, ≈ 200 s, `-t 420`) |
| `EX-88` | *(history and §9 briefs: `docs/planning/chunks/EX-88.md`)* **An insulated lead shaft in the gel: the sleeve takes the field that crosses the wire and moves its hot spot 90°; the field along the wire does not notice it** (`imp:6`; `example-runner`; opened 2026-10-04 03:00 daily review, §5.4 ramp from `IMP-5` ✅ `d2bc2ca9`, audited PASS; no example puts a dielectric layer on a conductor or compares coated and bare on one mesh) **Done-when:** `examples/implants/06_insulated_lead_shaft_in_gel.py` + same-stem guide (three required headings + `## Setup figure` + `## Related`, PNG ≤ 600 KiB); imports `IMP-5`'s fixture, helpers and bands, restates none; rule (a) keyword-only `return_fields` / `trace` on `_solve_coated`, module re-gated green in-slot, disclosed. At `-n 2`, degree 2, complex, one mesh, per drive: \|β − 1\| ≤ `TH15_BETA_BAND`, sleeve and host-annulus \|Φ_h/Φ_cf,h − 1\| ≤ `MAT4_BAND`, cavity dofs > 0 at max \|E\| < 1e-12, printed beside `20261004T035027Z_IMP-5-step1-W1a.log` / `…035247Z_…-W1b.log`; DG0 maps sum to Φ ≤ 1e-12 (region by cell tag). Asserted controls: TE coating omitted (W2a record) and the bare wire on the same mesh both miss. TE normal-E jump (≈ 59×), rotated hot spot and +30 % host SAR printed; TM coating effect (≈ 3e-6) stated below FEM resolution. `implants_06_insulated_lead_shaft_in_gel_combined.xdmf`. Census +1, `broken=0`; docrefs `exit != 1`. Not claimed: k_z ≠ 0 lead propagation, lossy or real insulation, a finite lead, a tip, a coil, heating, convergence. | ✅ 2026-10-04 (§9 item 223; audit PASS 10-05): `imp:6` ran at `-n 2` in 154 s (`20261005T033928Z_EX-88.log`, Status 0; unflagged control `20261005T034212Z_EX-88-control.log` readings unchanged); `IMP-5` re-gated 17 passed, digit-identical (`20261005T033404Z_EX-88-imp5-regate.log`); asserted C1 and bare-wire B0 both miss; census 83 → 84, `ok` 82 → 83, `broken=0`, docrefs 83 → 84, `exit=2` (stale only). | heavy (`-n 2`, ≈ 300 s, `-t 600`) |
| `EX-89` | `imp:7` — a lead's transfer function in the coil-driven gel phantom: one adjoint solve predicts the tip response to every drive (§5.4 ramp from `WF-9`) | ✅ 2026-10-05 (§9 item 238): `imp:7` W1 / W2 14 passed (184 / 167 s), (T3) ≤ 4.3e-12, C4 24.53 %, census 84 → 85, `broken=0`. Coated-lead rows added 2026-10-08 (§9 item 332): the `WF-9` step 3 lift re-gated at `-n 2` / `-n 4`, rows equal `20261008T142909Z_WF-9-step3-W1.log` digit for digit (280 s). **Done-when:** `examples/implants/07_lead_tf_in_the_coil.py` + same-stem guide (three required headings + `## Setup figure` + `## Related`, PNG ≤ 600 KiB); imports `WF-9` step 1's fixture, helpers and bands, restates none; rule (a) lift of `test_wf9_step1_identities`'s body with an additive return, module re-gated green at `-n 2` / `-n 4` in-slot, disclosed. *History: `docs/planning/chunks/EX-89.md`.* | heavy (`-n 2`, ≈ 280 s, `-t 540`) |
| `EX-90` | *(history and §9 briefs: `docs/planning/chunks/EX-90.md`)* **What a temperature probe in the heated gel ball reports: a ball-averaged series with a stated volume, beside the exact ball mean and the point value** (`mat:13`; `example-runner`; opened 2026-10-05 10:30 daily review, §5.4 ramp from `MAT-19` ✅ `0a0f85b1`; no example uses `thermal.probe` or states a probe volume) | ✅ 2026-10-05 15:00 CDT slot (§9 item 255): `mat:13` ran at `-n 2`, complex, 81 s (`20261005T201522Z_EX-90.log`; unflagged control 79 s, `20261005T201654Z_EX-90-control.log`, readings md5-identical); lift re-gated W1 `-n 2` 14 passed 106 s, W2 `-n 4` 14 passed 103 s, digits equal to `MAT-19` records, (W) 8.882e-16; census `examples` 85 → 86, `ok` 85 → 86, `broken=0`, docrefs runnable 85 → 86, `dead=0 guide=0`, exit 2 (stale 45 → 45). **Done-when:** `examples/materials/13_virtual_probe_in_gel_ball.py` + same-stem guide (three required headings + `## Setup figure` + `## Related`, PNG ≤ 600 KiB); imports `MAT-19` step 1's fixture, helpers and bands, restates none; rule (a) lift of the `run` fixture body and `test_mat19_probe_series`'s body with additive returns, module re-gated green at `-n 2` / `-n 4` in-slot, disclosed. At `-n 2`, P1, complex, 19 496 cells: (A0), (T) ≤ `P_BAND` at three centres × 40 times, (C1) misses, printed beside `20261005T095808Z_MAT-19-step1-W1.log`; an a = 1 / 4 / 8 mm sweep at the centre printed against its exact ball mean, and (C3) the point comparand at a = 8 mm printed beside its predicted miss (rule (e)). `materials_13_virtual_probe_in_gel_ball_combined.xdmf` (41 steps, θ_h / θ_exact / error, probe masks) + series CSV / plot. Census +1, `broken=0`; docrefs `exit != 1`. Not claimed: a real probe's mass, lead or sensing volume, radii other than 2 mm as gated, convergence, phantom / coil, compliance. | heavy (`-n 2`, ≈ 120 s, `-t 300`) |
| `EX-70` | A foil thinner than its skin depth | ✅ Done 2026-10-01: `th:13` reproduces the `TH-18` layered-foil reflection record for every termination to the record's digit, with the bulk-term control missing as predicted, the rule (a) lift re-gated green and censuses clean. Nothing is open; no copper gate or oblique-incidence claim is made. **Done-when:** asserts green through `./run_examples.sh -e th:13 -n 2`, artifacts named by full filename, both censuses before / after (`examples` +1), docrefs `exit != 1`, elapsed recorded. *History: `docs/planning/chunks/EX-70.md`.* | standard (host-runner, ≤ 180 s window) |
| `EX-71` | Does the total current close? ∇·J_tot at CG2 against its vacuous CG1 twin, on `POST-3`'s two-slab box … | ✅ Done 2026-09-30: `th:14` shows the total-current divergence on the two-slab box at CG2 against its vacuous CG1 twin, matching the `POST-3` record to the digit with the control and the gate re-run green and censuses clean. Nothing is open. **Done-when:** asserts green through `./run_examples.sh -e th:14 -n 2`, artifacts named by full filename, both censuses before / after (`examples` +1), docrefs `exit != 1`, elapsed recorded. *History: `docs/planning/chunks/EX-71.md`.* | standard (host-runner, ≤ 180 s window) |

**`EX-57` — one setup figure per example guide (standing, recurring)** 🟡 *(opened 2026-09-13 by operator directive: "at least one good visualization of the problem for each example in the examples folder, inserted into the .md files describing each example so people can see what the basic setup is"; step 0 executed interactively the same evening)*. Every script `./run_examples.sh --list` enumerates owes its same-stem guide a **`## Setup figure`** section that embeds one committed PNG from the group's `figures/` directory (`figures/<basename>_setup.png`, ≤ 600 KiB) and a caption naming the regions, the tags, the slice plane and what the reader should notice. **The infrastructure (step 0):** `src/fem_em_solver/post/setup_figure.py::write_setup_figure` — rank-safe (owned cells gathered to rank 0), headless PyVista, opt-in through `FEM_EM_SETUP_FIGURES=1` so the committed PNG never moves on a corpus run, two panels (tagged regions clipped in 3-D beside a mid-plane slice with mesh edges), one colour per region *class* so copper is copper in every figure; `scripts/testing/check_example_setup_figures.py` — the census (`ok` / `missing` / `broken`; exit 0 / 2 / 1 on the docrefs contract; `--next` prints the next example owing one in the runner's `--list` order); the docrefs checker searches `figures/` beside `paraview_output/`; exemplar `mesh:3` (`examples/meshing/03_birdcage_graded_conductors.py` + guide).
**Recurrence (the operator's "recurring background task"):** the daily review's step 6 appends **one `EX-57` figure item** to §9 On deck on every active day, naming the script `check_example_setup_figures.py --next` prints, until the census reads `missing=0`; the drained-queue fallback in §9 is the same item, so a slot that empties the queue draws the next figure instead of stopping. **Batched since 2026-09-25 (operator directive):** an item — queued or fallback — names up to three scripts, `--next 3` cut at the first that leaves the leading script's group directory, and **one** `example-runner` run draws them in order, one commit per figure, so a slot that runs out of clock keeps what it committed; the setup-figure and docrefs pre-censuses run once before the first figure and the post-censuses once after the last, predicting `missing` down by the number committed. Everything below that says "each item" applies to each figure. Reason (interactive transcript analysis, 2026-09-20 → 09-25): `example-runner` took ≈ 31 % of scheduled-session context, ≈ 65 turns per single figure, much of it the same start-up reads (this row, `setup_figure.py`, the harness log tails) repeated for every figure. Each figure is one example: add the `write_setup_figure` call right after the example builds its mesh (region names from the fixture's own tag map, air hidden, phantom translucent, the slice through the plane that shows the physics — the port plane, the coil mid-plane, the wave's propagation axis), run the example **once with the flag** through the harness at its recorded width, `git add` the PNG, write the `## Setup figure` section with its caption, then run both censuses (docrefs `exit != 1`, setup-figure `broken=0`) and re-run the example **without** the flag to show the default path is unchanged (same printed records). **Done-when per item (§4):** the census delta predicted before and read after (`missing` down by exactly one, `broken=0`), the example's imported assertions still green in the flagged run, elapsed recorded, PNG ≤ 600 KiB. **The committed PNG needs no pin edit** (`OPS-59`, 2026-09-20): `tests/unit/test_doc_reference_exit_codes.py` admits `examples/**/figures/*_setup.png` by *rule* — the census's own script↔figure pairing — so a figure lands without touching the artifact record. **A non-setup artifact is the opposite:** any other committed file under `examples/` (an `ans:` `metrics.json`, a validation PNG) still has to be named by hand in `LISTED_EXAMPLE_ARTIFACTS`, in the commit that commits it, or that test goes red on `main`. **Done-when for the chunk:** `check_example_setup_figures.py` exits 0 on `main`. **Traps:** the flagged run must be **the harness** (a bare shell render is unlogged — `EX-43` precedent); a region wholly inside the clipped half-space is skipped, not drawn — choose `clip_normal` so the ports survive; examples that build several meshes draw the *gated* rung; magnetostatics examples with a wire source and no conductor tag name the wire region so the copper colour applies; Ansys-benchmark cases (`ans:`) draw the geometry only — never a number from `aed_results/`. **Cost per item:** the example's recorded window plus the render (measured step 0: **+2 s** on `mesh:3`) plus censuses, ≈ 20 slot-min; **≈50 examples ⇒ ≈50 items**, one per active day, ≈ 10 weeks at the wind-down cadence, faster whenever a slot drains the queue. Tier: the example's own.
**Step 0 ✅ 2026-09-13 (interactive, on a worktree of `1abd55f`):** `mesh:3` at `-n 2` with the flag — Status 0, **31 s** (record 27.7 s ⇒ render ≈ +2–3 s), PNG **401 KiB**, gate records **0.966977 / 0.846150** reproduced digit for digit in both the flagged and the unflagged run; unflagged control Status 0, 29 s, no figure line. Census after: **`examples=48 ok=1 missing=47 broken=0`** (exit 2). Two earlier flagged attempts were red and fixed in place: an empty clipped port box (PyVista refuses an empty mesh — now skipped, and the half-space clip is opt-in) and a stray indent. Docrefs census on the worktree read `dead=88 guide=0`, all 88 the fresh checkout's absent `paraview_output/` artifacts and none the new PNG — not a corpus reading; `main`'s own census stands.
**Logs:** `20260914T024424Z_EX-57-step0.log` (flagged, Status 0, 31 s), `20260914T024455Z_EX-57-step0-control.log` (unflagged, Status 0, 29 s), `20260914T024524Z_EX-57-step0-census.log` (both censuses); reds `20260914T022747Z_EX-57-step0.log` (Status 142, the empty clipped box), `20260914T024312Z_EX-57-step0.log` (Status 1, the indent).

**`EX-26` — Poynting power-balance audit** ✅ *(2026-08-20, 12:00 slot; commissioned 2026-08-20 03:00 review, §5.4 ramp on `POST-5` step 4; audited COMPLIANT 2026-08-21 18:00 review)*. `examples/time_harmonic/08_poynting_power_balance.py` + same-stem guide, `th:8`. **Closed as written, both fixtures on one run, no band moved.** Driven cylinder three-term **16.7465%** inside the imported `POYNTING_IMBALANCE_MAX` = 25%, two-term 116.7465% asserted to *miss* (inverted control); `TH-6` plane wave source-free **8.185716%**, legs 8.1205% / 0.0711% inside `POST5_STEP3_LEG_BAND` = 10%; J = 0 source term `== 0.0` W with 7 other keys bit-identical; σ-blind residual 83.2535% = **4.97×** (floor 3.0×, ceiling 5.97×); impressed-source term = 100.0% of the largest term. All 8 records inside a 1% band, worst drift 3.00e-04. Restated with provenance: `TH6_RECORD_IMBALANCE` = 0.08185716, `TH6_RECORD_FLUX_ERROR` = 0.081205, `TH6_RECORD_DISSIPATED_ERROR` = 0.000711, `TH6_CELLS` = 10368. Two combined XDMFs (`E` CG1; `B` and `½Re(E×H̄)` as DG0 — honest resolution of a degree-1 `curl E`, faceted in ParaView by choice). 1 405 / 10 368 cells, 4.7 s in-script, 8 s harness at `-n 2`.
**Logs:** `20260820T170422Z_EX-26-example-n2.log` (exit 0), `20260820T170540Z_EX-26-docrefs.log` (`dead=0 guide=0 stale=0 exit=0`).
**Carry-forwards:** commissioned standard, **measured smoke** — the 152 s estimate belonged to the `TH-6` file's other tests, not the 12³ rung; demonstration only, no SAR or coil-loading claim.
Full narrative: `docs/planning/plan-archive.md`, entry «§7 EX-26 full narrative — archived 2026-08-23 (weekly review)».

**`EX-25` — degree-2 Larmor sphere: accuracy-per-cost side by side** ✅ *(2026-08-19, 09:00 slot; commissioned 2026-08-18 10:30 review, §5.4 ramp on `TH-12` step 1)*. New `examples/time_harmonic/` example on the `TH-10` lossy saline sphere (`EX-19` fixture, imported): the same coarse 5 866-cell mesh solved at degree 1 and degree 2 in one run, side-by-side table + combined XDMF of both solutions (one mesh, one CG1 export space — the only variable is the element). **Closed as written, no band moved:** all four records inside the imported 1% `REPRODUCTION_BAND` (worst drift 1.48e-03, degree-2 power error) — degree 2 **0.1405%** relL2 / **0.0058%** power, degree 1 **8.1541%** / **8.3869%**; `|Im P|/Re P` exactly 0.0 at both orders; inverted control asserted both ways (degree 1 misses the 3.643% fine-rung record, degree 2 beats it); DOFs asserted exactly **7 591 / 39 634**. Cost ratios 5.22× DOFs → 2.02× wall, 2.74× RSS — sublinear, the opposite shape from `TH-12` step 2's coil (≈20× wall for 5.42× DOFs); the 2.02× includes mesh + assembly (step 1's solve-only 4.32×), neither gated. Constants from `tests/validation/test_lossy_sphere_degree2.py` and the `TH-10` module (`ANS-1`); complex-mode XDMF `real_*`/`imag_*` split (`OPS-21`); summed `ru_maxrss`, not `memory.peak`.
**Carry-forwards:** sphere fixture only — no production-order claim (that is the weekly review's per the `TH-12` decision clause; decided 2026-08-23, §10); the example prints both cost shapes side by side precisely so this fixture's cheap second order is not generalised.
Full narrative: `docs/planning/plan-archive.md`, entry «§7 EX-25 full narrative — archived 2026-08-23 (weekly review)».

**`EX-24` — lumped-sheet port at interior width** ✅ *(2026-08-18, 19:30 slot; commissioned 2026-08-17 18:00 review + 2026-08-18 03:00 addendum, §5.4 ramp on `PORT-9` steps 2b/2c; audited COMPLIANT 2026-08-19 03:00 review)*. Example on the `GEO-16` port-sheet solve fixture: the f ∈ {1.0, 0.735, 0.5} width ladder as three lumped-BC solves plus, at f = 0.5, the two-port sweep through `LumpedSheetPortSpec`; combined XDMF with sheet tags; same-stem guide. Constants imported from `tests/validation/test_port_lumped_narrowed_sheet.py`, `tests/validation/test_port_lumped_sheet_sweep.py` and the two-torus module (`ANS-1`). **Closed as written, both legs on one run, no band moved:** ladder 7.7095% / 3.6730% / **1.8333%** (gate at f = 0.5 against 5%), f = 1.0 reproducing both step-1 records (7.7095%, gap ratio 0.894310) inside 1e-4 *and* asserted to miss the band, open-limit identity ≤ **1.8e-15** per width, sweep reciprocity **2.574296e-11** against 1e-3, gap route asserted flat across the ladder (0.894310 → 0.894349, drift 3.9e-5 against `REPRODUCTION_BAND`). Both legs share one mesh (non-mutating midpoint filter): 237.5 s in-script against the ≈260 s estimate, standard tier at `-n 2`.
**Logs:** `20260819T003401Z_EX-24-example-n2.log` (every plan number matched digit-for-digit; footers 4 s / 239 s / 1 s); docrefs exit 2 staleness-only, all 24 stale artifacts pre-existing `EX-22` backlog.
**Carry-forwards:** the sweep's cross-route sits **≈0.23 pp below** the ladder's at the same width (1.6079 / 1.5950% vs 1.8333%) — the impressed-sheet drive reads closer to the centreline than the impressed-gap drive, a systematic `PORT-9` step 3 should expect rather than debug; `w = A/h`, never the bbox extent; two-torus only, no birdcage claim.
Full narrative: `docs/planning/plan-archive.md`, entry «§7 EX-24 full narrative — archived 2026-08-23 (weekly review)».

**`EX-23` — two-torus port-sheet mesh example** ✅ *(2026-08-17, 09:00
slot; commissioned 2026-08-17 review, §5.4 ramp on `GEO-16`'s newly gated
capability)*. **Closed as written** — every element of the rubric below
executed, no band moved and none needed to move: both area identities
landed at `1.000000000000`, the kwarg-off control reproduced 79 534 cells
with the sheet tags absent, and the docrefs checker exits 2
(staleness-only, none of it this example's). Two things worth carrying
forward: the sheet mesh costs **+354 cells** over the control
(79 888 vs 79 534) — the fragment is nearly free, which is the number
`PORT-9` step 3 should budget from; and the measured `w/h = 1.504225878`
sits **1.3e-05 relative** above the generator's own CAD-side
`squares_w_over_h = 1.504206917`, the arc-chord difference between the
CAD surface and its triangulation, printed side by side in the example.
*(Original plan below.)*
`examples/meshing/01` shows the gapped two-torus and `03` the graded
birdcage; no example shows an **interior sheet surface** — the
known-issues-9 pattern (facet tag rebuilt from cell tags on the dolfinx
side, never a dim-2 gmsh group) that `GEO-16` gated and `PORT-9` now
builds on. **Do:** `examples/meshing/04_two_torus_port_sheet.py` + the
same-stem guide (`EX-15` rule: same commit), dispatched through
`./run_examples.sh`, producing combined-XDMF (cells + facets, tags
`211`/`212` visible) that opens in ParaView. Import every constant from
`tests/mesh/test_two_torus_port_sheet.py`, none restated (`ANS-1`
pattern): assert both sheets' meshed/CAD area to `AREA_IDENTITY_BAND`
(= 1e-9; `GEO-16`'s gate, `1.000000000000` on record), the 211/212
symmetry identity < 1e-12, and the kwarg-off negative control (recorded
cell count `79 534`, no `21x` tags — the `EX-18`/`EX-21`
inverted-assertion pattern); print the measured extents
(`w/h = 1.504225878` on this fixture's parameterisation). **Tier/cost:**
standard, mesh-only, no solve; `GEO-16`'s runs were 31/49 s at `-n 2` —
budget ≈60 s, `timeout -k 30 480`. **Traps:** the sheet kwarg splits
each gap volume into **two** cell tags (`101`+`111`, `102`+`112`) —
select both; docrefs checker gates on **`exit != 1`** (`OPS-19`
contract), staleness is information; facet counts asserted `> 0` before
any area identity (vacuous-pass guard). **Scope:** mesh-only — no port,
no solve, no `Z` claim; the solve fixture's differently-parameterised
`w/h = 0.745249896` belongs to `PORT-9`, not here. **Negative result:**
the area identity failing on the example path is a regression against
`GEO-16`'s gate — known-issues entry, report, stop.

**`EX-22` — restore the absent example artifacts** ✅ *(closed 2026-08-19,
19:30 slot — see the §7 table row for the full closure numbers. Done-when
met in three runner commands rather than the entry's two: `mag` was split
`-e 1,2,4` / `-e 5,6` so neither container window could approach the
foreground ceiling, since the only prior all-mag timing on record (204 s)
predates examples 05 and 06 and would have left the group unpriced at
≈390 s. Measured: 230 s + 151 s + 6 s. The checker reads `exit=0` for the
first time under the `OPS-19` contract; nothing was ever dead, so the
2026-08-16 premise correction is confirmed a second time and the chunk was
a freshness restore, never a recovery)* *(commissioned
2026-08-16, weekly review — examples-health audit)*. Six examples' gated
`paraview_output/` outputs are **absent on disk**, not merely stale:
every `straight_wire_*`, `helmholtz_*`, `gauge_cross_check_*`,
`h_convergence_rate_*` and `mri_coil_phantom_*` artifact plus
`circular_loop_B.bp` is gone (gitignored, so unrecoverable without
reruns) — the source of the doc-reference checker's 24 standing
stale-reference violations (`OPS-19` owned the *policy* split and landed
2026-08-16 — those 24 now score `exit=2`, staleness-only, and no longer
mask a real defect; this chunk restores the *artifacts*).
**Premise correction, measured 2026-08-16 16:30 slot
(`20260816T213312Z_OPS-19-step1-rerun.log` lines 44–68):** the 24 are
`dead=0 stale=24` — every one of them **exists** in `paraview_output/`,
aged 145.5–151.4 h, including `circular_loop_B.bp`. "Absent on disk" does
not hold at this commit; the chunk's refresh work is unchanged, but its
done-when should be read as 24 → 0 *stale*, and whether any artifact is
genuinely missing needs re-auditing before the runs are sized. **Do:**
runner refresh runs, `mag` 01/02/04/05/06 and `mri:1`, through the
harness. **Done-when:** each run exits 0 with its recorded anchors
reproduced by the examples' own asserts (`EX-14`/`EX-17` round-trip
identities included), the six examples' artifacts exist on disk, and the
doc-reference checker's stale/dead count drops 24 → 0 (guide pass stays
green — 21/21 since `EX-21`; the checker under the `OPS-19` contract then
reads `exit=0`). **Tier:** heavy — `EX-9`'s convergence example alone is ≈130 s;
run as two runner commands (`mag` group, then `mri:1`), each wrapped
`timeout -k 30 500`. **Trap:** `mri:1` is the labelled *ungated* example
— reproduce its printed record, do not invent a gate for it. **Negative
result:** an example that no longer reproduces its record is a real
regression — known-issues entry, report, stop; never refresh past a
failure.

**`EX-4`…`EX-12` — the 2026-08-09 backfill, all ✅ by 2026-08-10** *(full
plans + closure narratives in `docs/planning/plan-archive.md`)*. Common
pattern, held across all nine: fixture and constants **imported** from the
gate test, never restated; the gate record reproduced digit for digit
through the example path; the exported array itself gated, not merely
written; runner-dispatched harness logs. Per-example anchors on record:
`EX-4` α/β to 0.0185% / 0.0593% (`th:1`); `EX-5` four cavity modes to
≤ 0.0436%, Rayleigh-quotient export identity 3.48e-15 (`th:2`); `EX-6`
interior E_z 2.443%, volume-integral corroboration 0.014% (`th:3`); `EX-7`
γ = 37.650399 Np/m at 0.006%, TE₁₀ profile RMS 0.200% (`th:4`); `EX-8`
resonance guard 137.554 vs threshold 50, pole law 3.156% (`th:5`); `EX-9`
h-convergence rate 1.1009 in the gate's (0.7, 1.5) band, heavy tier, plus
the CG1-export finding (smoothing costs 7.89 pp on a 1/r field); `EX-10`
gauge cross-check 0.0004% probe / 0.0033% volume with an 11-order |A|
separation (`-e 5`); `EX-11` Dodd–Deeds ΔR 1.5834% byte-matching `MAT-6`
step 3, σ = 0 control exactly zero (`mat:1`, feeds `ANS-1`); `EX-12`
examples hygiene + the doc-reference checker (its freshness/`.bp` findings
spawned `EX-14`/`EX-17`).

**`EX-15` ✅ 2026-08-11** *(operator directive 2026-08-10; full plan +
three step closures in `docs/planning/plan-archive.md`)*. Every runnable
example ships a same-stem analysis guide with three required headings
(policy stated in §5.4), enforced by a guide pass in
`scripts/testing/check_example_doc_references.py` that reads the example
set from `./run_examples.sh --list`. Closed at **16 of 16** guides,
`PENDING_GUIDES` empty (`20260811T110627Z_EX-15-step3-refcheck-final.log`);
negative controls (missing guide, missing heading) fired in all three
steps. Standing rule: a new example must ship its guide in the same
commit; on-record numbers in guides are copied from §7/gate records and
cited by log name, never re-measured.

**`EX-13` 🚫 2026-08-10, closed negative** *(full plan + result in
`docs/planning/plan-archive.md`)*. The `MAG-6` gate-fixture rank-spread
reading (0.024%) does not transfer to `examples/mri/01`: floor spread
23.5545% vs the < 5% anchor, sub-floor 23.3010% (no discrimination —
`TimeHarmonicSolver.solve` ignores `gauge_penalty`, so the E leg is inert
by construction). The 23% measured the unconverged GMRES iterate; salvage
scoped as `EX-16`.

**`EX-16` 🚫 2026-08-10, closed negative** *(full plan + result in
`docs/planning/plan-archive.md`)*. Converging `examples/mri/01`'s
time-harmonic solve (direct `preonly`/LU, `reason=4`) did **not** move the
23% `-n 2`/`-n 4` centerline spread (23.5539% vs the 23.5545% unconverged
record) — the convergence hypothesis is refuted. The decisive positive
control: the 493-point phantom-region sampler agrees to 0.007326% on the
same fields, **3215×** tighter, so the defect is the centerline
point-evaluation path (on-axis points on shared mesh edges, the `MAG-6`
step-4 mechanism). Known-issues entry stays open, re-pointed at
`evaluate_vector_field_parallel` (assigned `POST-4`). The code fix landed
on its merits; `WF-1` stays 🧪.

**`EX-17` ✅ 2026-08-10** *(full plan + closure in
`docs/planning/plan-archive.md`)*. The `EX-14` diff ported to
`02_circular_loop.py`: round-trip max |B| 7.756122914931e-05 T both ways,
rel diff 0.000e+00 vs 1e-10 on a mesh 30× larger
(`20260810T200154Z_EX-17-gate-mag2.log`, 124 s); the loop's analytic
numbers (6.3046% / 13.5037%) unmoved; known-issues entry retired.

**`EX-18` ✅ 2026-08-13** *(doc repairs 2026-08-16; full plan + closure in
`docs/planning/plan-archive.md`)*. `ports:1`, the first ports example:
gap-voltage pair → Z → S on the two-torus fixture, reproducing the
3b-xviii digits — raw 0.894543 × ωM₁₂ printed first and labelled the miss
it is, corrected 0.939849 (−6.02%) inside the unmoved 10%,
‖S−Sᵀ‖/‖S‖ = 2.5494e-05, ‖S‖₂ = 0.861449; blind-fixture control −98.26%
asserted to fail (`20260813T110940Z_EX-18-example-n2-v3.log`, 135 s). The
systematics ladder was lifted into `src/fem_em_solver/ports/systematics.py`
with a bit-identity test, so example and gate share one definition. Doc
repairs 2026-08-16: guide pass 3 → 0 violations, the 400× band-margin
comment corrected to 7.7× with the band value untouched
(`20260816T033121Z_EX-18-docrefs-fix.log`).

**`EX-14` ✅ 2026-08-10** *(full plan + closure in
`docs/planning/plan-archive.md`)*. Straight-wire VTX export repaired
(writers handed the Lagrange interpolants, split `try`); round-trip
identity exact — read-back max |B| 4.463805898300e-05 T, rel diff
0.000e+00 vs 1e-10 (`20260810T140337Z_EX-14-gate-mag1-v2.log`). The
freshness negctl caught a second real defect: a `.bp` is a directory whose
own mtime never updates — fixed with `artifact_mtime()`. Filed, not fixed:
the identical defect in `02_circular_loop.py` (became `EX-17`).


**`EX-1` ✅ 2026-08-07** *(landed 2026-08-06, demoted for missing runner
logs, restored by `20260807T003044Z_EX-1-runner-mesh1.log`; full history in
`docs/planning/plan-archive.md`)*. `mesh:1` builds the gapped two-torus
fixture (79 534 cells) and asserts the `GEO-8`/`GEO-10` closed-form
identities allreduced — area and volume ratios 1.000000000000, wire ratios
0.963633 / 0.963756 — through the runner. The lesson that outlived it:
§5.4 examples close only on a logged `./run_examples.sh --list` + `-e`
dispatch, not a byte-equivalent direct invocation.

**`EX-2` ✅ 2026-08-07** *(full plan + closure in
`docs/planning/plan-archive.md`)*. `mesh:2` builds `cylindrical_domain()`
at the `GEO-13` defaults and reproduces its classification record exactly
(3 of 6 surfaces accepted, wall 1.111111e-04 × tol, interior
9.999989e+01 × tol); partition identity 1.000000000000000; the plan's
per-tag inscription-band premise was refuted by measurement — at
`resolution = 2 × inner_radius` the inner cylinder is a gmsh heptagonal
prism, gated in closed form at 1.11e-16 (asserted 1e-12). Caller audit:
the 28% inner-volume deficit is latent in every repo caller, armed only if
a test ever gates an inner-region quantity. Log `20260807T140554Z_EX-2.log`.

**`EX-3` ✅ 2026-08-08** *(full plan + closure in
`docs/planning/plan-archive.md`)*. `mri:2` reproduces the `MAT-4` step-3
record digit for digit through the example path: `SAR_avg/SAR_point` =
1.00000000 at 1 g and 10 g, kernel mass 0.0120% / 0.0044%, pointwise vs
closed form 4.96e-16; the exported DG0 array re-averages to the closed form
to 1.32e-15; surface-ball negative control 2.1894 vs lens ceiling 2.1681.
Log `20260808T020414Z_EX-3-gate.log` (14 s, `-n 2`). Imposed field only —
no SAR-on-a-coil claim.

**`EX-29` — the doc-reference checker must freshness-gate every example's own `paraview_output/`** ✅ *(closed
2026-08-24 06:00 slot; commissioned 2026-08-23 weekly review; smoke, `-n 1`, no solve. Full narrative archived in
`docs/planning/plan-archive.md`, entry «§7 EX-29 full narrative — archived 2026-08-30 (weekly review)».)* Defect: `check_example_doc_references.py` resolved only the repo-root directory
and exempted any basename found under `examples/`, so **22 of 27** runnable examples were never freshness-checked.
> * **Done:** artifacts resolved at their example-relative path with `--max-age-s` applied to all; exemption
>   restricted to `git ls-files`-tracked paths, **pinned by path** in `COMMITTED_EXAMPLE_ARTIFACTS` (three tracked
>   artifacts exist — the "assert empty" instruction would have been false); orphaned
>   `examples/magnetostatics/paraview_output/` deleted.
> * **Gated:** `tests/unit/test_doc_reference_exit_codes.py` **15 passed in 3.71 s** (`20260824T110512Z_EX-29-unit.log`,
>   `-n 1`), whole `tests/unit` 22 passed (`…110540Z_…-run2.log`); negative control on the same tree — pre-fix
>   `stale=24` (`20260824T110150Z_EX-29-prefix-control.log`) vs post-fix **`stale=55`** (`…110531Z_EX-29-census.log`),
>   gate (d)'s independent walk reproducing `stale=55 checked=58 hidden_pre_fix=32` exactly.
> * **Findings:** in-container `git` needs `-c safe.directory=` (dubious-ownership on the bind mount, silently
>   turning tracked artifacts into `dead=1`); known-issues entry closed. `EX-30` re-sized from 55, not 24.

**`EX-30` — refresh the stale example-artifact set the checker could not see** ✅ *(2026-08-26, all four legs;
commissioned 2026-08-23 weekly review, re-scoped 2026-08-24 10:30 review to legs (th) / (root) / (mesh) / (ports)
from the honest `stale=55` census. Every leg's journal, red and census derivation archived in `docs/planning/plan-archive.md`, entry «§7 EX-30 full narrative — archived 2026-08-30 (weekly review)».)*
> * **Leg (th) ✅ 08-25:** all eight `th:` examples green, 105 s; licensed 128 MHz alignment 0.01826 → **0.01769** /
>   57.31 → **59.16** from `20260822T123746Z_OPS-18-step3-th10-rerun.log`, `th:6` reproducing 1.769% (2.02e-04) and
>   59.16× (5.45e-05), 64 MHz 3.643% unmoved; `time_harmonic` census 4 → 0. Reds spawned `OPS-24`, `OPS-25`.
> * **Leg (mesh) 08-25 (4/7 green, licence unused):** `mesh:1` 79 070, `mesh:2` 5 717 cells; reds spawned `GEO-21`
>   (graded-conductor gate), the `GEO-16` 79 534 → 79 070 record class (→ `OPS-27` precedent), `mesh:5`'s control (→ `GEO-17` row).
> * **Leg (ports) ✅ 08-26 (935 s):** `ports:2` restated `‖S−Sᵀ‖/‖S‖` 2.5494e-05 → **4.758625e-05** and `‖S‖₂` 0.861449 →
>   **0.864809457** (the gate's power-wave digits, `PORT-9` leg (d3)), `ports:1` → 3.11213e-05 / 0.861356895; `ans:1`
>   ΔR 1.5838% vs the 1.5834% record; `ports` 4 → 0, `ans` 2 → 0 (`20260826T12…Z_EX-30-ports-*.log`).
> * **Leg (root) ✅ 08-26 (447 s):** `mag:1` `resolution` 0.01 → **0.008**, 21 830 cells, `B(3 mm) = 6.666667e-05 T`
>   (`20260826T170155Z_EX-30-root2-run-mag1.log`); `mag:6` green through `MAG-19` on 21.8417 / 15.3848 / 4.4605%; (1\*)
>   guide re-records `mag:2` 409 596 / 6.2134%, `mag:4` 69 918 / 103 950 / 160 677, `mri:1` 1.979842e+02 (commit `c466143`);
>   earlier `mag:5` 14 055 cells all 7 assertions, `mri:2` SAR 3.31e-16 / 2.81e-15, `mat:1` ΔR 1.5838%. * **Census:** **`dead=0 guide=0 stale=0 exit=0`** (`20260826T171345Z_EX-30-root2-postcensus.log`), the first clean corpus-wide reading since `EX-29`. **Carry-forward:** the `straight_wire_domain` floor entry → `GEO-22` (OPEN).

**`EX-31` — ring-gapped birdcage with dual port families** ✅ *(2026-08-24; `mesh:7`,
`examples/meshing/07_birdcage_ring_gap_ports.py` + same-stem guide, closed as written on the first run; commissioned
2026-08-24 10:30 review as `GEO-20` step 1's ramp obligation. Result and commissioning text archived in
`docs/planning/plan-archive.md`, entry «§7 EX-31 full narrative (result block + commissioning block) — archived 2026-08-30 (weekly review)».)*
> * **Run:** `20260824T200613Z_EX-31-example-n2.log`, exit 0, **70.6 s in-script / 75 s harness at `-n 2`**, real build;
>   every gated figure reproduces `GEO-20` step 1's log (`20260824T124525Z`) to the printed digit.
> * **Ring-gapped rung (110 786 cells = record; 4 leg + 8 ring ports):** terminals **0.974454791 / 0.974454832** of
>   `1.005309649e-04 m²` inside [0.95, 1.0] and the `0.974455 ± 1e-5` record band, spread 2.099e-08 vs 1e-5; closure and
>   volume/analytic **1.000000000000** on all eight; sheets 14 facets, `1.000000000e-04 m²`, meshed/analytic
>   1.000000000000; C4+mirror spreads 1.666e-15 / 2.443e-16 vs 1e-12; Pappus 1.000000000000; conductor 0.969275 ≥ 0.95.
> * **12-port rung (128 402 cells = record):** both identity families exact; leg terminals 0.988615809–0.988615855.
>   **Negative control (inverted):** `ring_gap_length=None` gives tags `[1, 2, 3, 101-104]`, 98 666 cells (ratio
>   1.001950 of the 98 474 record, inside 1%), `_global_facet_count = 0` on every ring sheet group 215–222.
> * **Hoisted to the gate (`ANS-1`):** `RING_TERMINAL_RATIO`, `RING_GAP_CELL_RECORD`, `LEG_RING_CELL_RECORD` now asserted
>   in `tests/mesh/test_birdcage_ring_gaps.py` (`20260824T200739Z_EX-31-gate-module.log`, 2 passed / 70.4 s). Docrefs
>   `…200857Z` `stale=55 exit=2`, 28 runnable examples scanned.

**`EX-32` — birdcage 4-port power-wave S-matrix at 10 MHz** ✅ 2026-08-26 *(15:00 slot, green on the first run;
commissioned 2026-08-25 10:30 review as `PORT-9`'s §5.4 ramp. Full narrative archived in `docs/planning/plan-archive.md`, entry «§7 EX-32 full narrative — archived 2026-08-30 (weekly review)».)*
`examples/ports/04_birdcage_four_port_sparameters.py` + guide, `ports:4` — the first example solving a port on the coil;
it calls leg (d)'s own `build_four_port_sweep()` (lifted to module level, additive), never re-implementing the sweep.
> * **Run:** `20260826T200545Z_EX-32-run1.log`, **88 s** / 85.0 s in-script, Status 0, `-n 2` complex. **Every gate record
>   reproduced:** 116 085 cells at ratio **1.000000**, σ_max **0.999992805**, max column power **0.793823974**, class
>   means 2.338160261e+01 / 1.700854304e+01 / 1.606048044e+01 Ω with spreads **0.0553 / 0.0353 / 0.0214%** vs the
>   imported 0.5%, separation **166.6766×** vs 10×; P1 column vs `LEG_D0_Z_COLUMN` **1.071e-10 … 2.568e-10** vs 1e-9.
> * **Negative control asserted:** the retired `PORT-0` heuristic on the same mesh keeps `is_placeholder=True`, prints an
>   identically zero off-diagonal and separates from the field-derived S by **6.446452e-01** vs the `EX-20` 2e-3 floor.
> * **(d3c) reading:** `‖S − Sᵀ‖/‖S‖` 4.183068067e-13 here vs the module's ≈2.152e-14 — both ≈11 decades under 1e-3; first
>   evidence the "order of magnitude only" rule is ≈1.3 decades wide. Nothing re-recorded.
> * Gate module green after the refactor (`20260826T200746Z_EX-32-gate.log`, 16 passed in 71.98s); census
>   **`dead=0 guide=0 stale=0 exit=0`**, 30/30 guided (`…200908Z_EX-32-census.log`). Scope: 10 MHz only, `PORT-9` caveat verbatim.

**`EX-34` — birdcage 4-port S-matrix across the Larmor frequency ladder, one mesh** ✅ 2026-08-28 *(06:00 slot,
green on the second run — the first hit the runner trap: `run_examples.sh` drives `docker` and must run on the
host; commissioned 2026-08-26 18:00 review as `PORT-11`'s ramp. Full narrative archived in `docs/planning/plan-archive.md`, entry «§7 EX-34 full narrative — archived 2026-08-30 (weekly review)».)*
`examples/ports/05_birdcage_larmor_frequency_ladder.py` + guide, `ports:5`: **one** 116 085-cell mesh (ratio 1.000000,
24.0 s) reused by all three rungs (asserted `reused_mesh`), twelve driven solves, `_four_port_rung` imported.
> * **Run:** `20260828T110615Z_EX-34-run2.log`, **139 s** / 136.8 s in-script, Status 0, `-n 2` complex. Gates imported
>   and asserted on every rung (1e-3 / 1 + 1e-9 / 0.5% / 10×): **10 MHz** 1.657e-14, σ_max 0.999992805, 0.793823974,
>   0.0553 / 0.0353 / 0.0214%, 166.6766×; **64 MHz** 1.179e-15, 0.999721388, 0.804704664, 0.0573 / 0.0599 / 0.0370%,
>   671.0527×; **128 MHz** 5.457e-15, 0.998974779, 0.861668762, 0.1012 / 0.0916 / 0.0654%, 576.9483×. Pre-gate stop
>   rule first: cells/λ **12.5024** ≥ 10 (cells/δ 5.1845), loss tangent 11.5225 → 1.8004 → 0.9002 up the ladder.
> * **Anchors:** 10 MHz reproduces leg (d)'s 4×4 to **1.158e-10** vs 1e-6 and (d0)'s column to 2.568e-10; Larmor rungs
>   reproduce `PORT-11` step 2/3 inside 1% (σ_max 2.814e-10 / 4.374e-11; worst spread misses 1.075e-03 / 6.755e-04 =
>   print precision). **Negative control at 128 MHz:** `PORT-0` heuristic zero off-diagonal, separation **1.585461e+00** vs 2e-3.
> * Gate module green after the additive `reuse=` parameter (`20260828T111019Z_EX-34-gate.log`, 5 passed in 103.07s);
>   census `dead=0 guide=0`, 31/31 guided (`…111008Z_EX-34-census2.log`). **Scope verbatim in the guide: self-consistency
>   identities on one fixture — no resonance, tuning, B1+/SAR or absolute-accuracy claim.**

**`EX-33` — 16-leg gapped + sheeted birdcage mesh** ✅ 2026-08-26 *(13:30 slot, green on the first run;
commissioned 2026-08-25 18:00 review as `GEO-19`'s §5.4 ramp; mesh-only, no solve, no F-human claim. Full narrative
archived in `docs/planning/plan-archive.md`, entry «§7 EX-33 full narrative — archived 2026-08-30 (weekly review)».)* `examples/meshing/08_birdcage_sixteen_legs.py` + guide, `mesh:8`.
> * **Run:** `20260826T183240Z_EX-33-run1.log`, **131 s** / 127.7 s in-script, Status 0, `-n 2` real. The `GEO-19` identity
>   family asserted by the gate module's own `_assert_identity_family` on this run's mesh: partition / air box
>   1.000000000000, 32 halves 0.500000000000, 16 sheets `dx·g` 1.000000000000, C16 sheet spread **1.331e-15** vs 1e-12,
>   closure 1.000000000000, conductor meshed/CAD **0.981503**, separation margin **1.560723×**; azimuth classes
>   `aligned` (8) 0.988615772, 22.5° (4) 0.989367514, 67.5° (4) 0.989449735 — intra **1.923e-07 / 5.849e-08 / 6.144e-08**
>   vs 1e-6, inter **8.431e-04** vs 5e-3.
> * **Negative control asserted:** the in-run 4-leg build reports **one** azimuth class (spread 3.184e-08, inter-class
>   0.000e+00) at **116 085** cells, relative 0.000e+00 against the step-B record.
> * **Cost rung (printed):** 4 → 16 legs, cells 116 085 → 307 296 (**2.6472×**, record reproduced exactly), mesh 26.51 →
>   84.25 s (**3.1777×**) — mesh seconds grow faster than cells.
> * Gate module green after the additive `_measure` return (`20260826T183618Z_EX-33-gate.log`, 2 passed in 124.56s);
>   census **`dead=0 guide=0 stale=0 exit=0`**, 29/29 guided (`…183831Z_EX-33-census.log`). No band or record moved.

**`EX-35` — 16-leg ring-gapped birdcage mesh: the 32-ring-port high-pass layout** ✅ *(2026-08-29 15:00 slot,
green first attempt; commissioned 2026-08-29 10:30 review as `GEO-20`'s §5.4 ramp; mesh-only, legs uncut — 32 sheets,
not 48; no solve, no F-human claim. Full narrative archived in `docs/planning/plan-archive.md`, entry «§7 EX-35 full narrative — archived 2026-08-30 (weekly review)».)*
`examples/meshing/09_birdcage_sixteen_ring_gaps.py` + guide, `mesh:9`, no runner edit needed.
> * **Run:** `20260829T200308Z_EX-35-run1.log`, **104 s** at `-n 2` / 101.1 s in-script, Status 0. Every commissioning
>   record reproduced to the digit: **265 621 cells** (relative 0.000e+00), terminals **0.974454791–0.974455668**, C32
>   sheet spread **4.985e-16** vs 1e-12, conductor meshed/CAD **0.976465** ≥ 0.95, Pappus on the 32 arcs
>   3.134786420778e-05 / 3.134786420778e-05 = 1.000000000000, partition / closure / wedge volume / `w²` sheets
>   1.000000000000 on all 32. Anchor: `_assert_ring_identity_family` imported from the gate module, both rungs.
> * **Azimuth classes 4 at 16 legs / 1 at 4:** intra 4.198e-08 / 4.498e-07 / 4.681e-07 / 8.997e-07 vs 1e-6, inter-class
>   **3.315e-07** vs 5e-3 — four orders inside, against the leg family's 8.431e-04 (`EX-33`): ring-gap cut faces are
>   exact planar disks whose triangulation barely notices azimuth.
> * **Negative control:** in-run 4-leg ring rung **110 786** cells (relative 0.000e+00 vs `RING_GAP_CELL_RECORD`), one
>   class, all eight terminals asserted against 0.974455. **Cost rung (printed):** 110 786 → 265 621 cells (2.3976×),
>   mesh 22.29 → 66.95 s (3.0042×); leg-arc clearance 3.744468e-03 m at 16 legs (Phase 6 count study).
> * Gate module green after the additive `_measure_ring` return (`20260829T200504Z_EX-35-gate-rerun.log`, 1 passed in 183.33 s).

**`EX-36` — re-run the examples the 2026-08-28 rename orphaned until the
census reads `dead=0`** ✅ *(commissioned 2026-08-30 weekly review; four
legs by runner group 2026-08-31 … 09-01; full 166-line entry archived in
`docs/planning/plan-archive.md`, 2026-09-06 weekly review).* Result: legs
(th), (mesh), (root+mri+mat) and (ports+ans) each re-ran through the
harness with pre/post census, `dead=53 exit=1` → `dead=0`; the mesh leg's
first window lost its footer (`…050408Z`, superseded by the remainder
windows). All four legs were `example-runner` executions, the Sonnet-tier
experiment's first sample (§10 agent-value ledger). Live carry-forward:
none — the artifact-naming rule it enforced is `example-runner.md`
non-negotiable 6. Logs: `20260831T123115Z_EX-36-leg-th-a.log`,
`20260901T050142Z_EX-36-legs-precensus.log`,
`20260901T093312Z_EX-36-leg-mesh-remainder-precensus.log` (9 in the
archived entry).

**`EX-37` — restore the `ANS-1`/`ANS-3` example imports the rename broke**
✅ *(2026-08-31 19:30 slot; full 69-line entry archived in
`docs/planning/plan-archive.md`, 2026-09-06 weekly review).* Result: the
`__import__` strings re-pointed to the prefixed module names; negative
control Status 1 / 3 s, both cases green (63 s / 128 s); records reproduce
the 08-16 figures with AED columns still blank. Live carry-forward: the
run-to-run Z/S scatter it measured (≤ 5e-8) is the basis of the 1e-6
reproduction controls (ruled 2026-09-02 weekly). Logs:
`20260831T003025Z_EX-37.log`, `20260831T003409Z_EX-37-docrefs.log`.

**`EX-38` — the first `|B₁⁺|` field in ParaView: loaded 4-leg birdcage at
10 MHz** ✅ *(2026-08-31, 15:00 slot; commissioned 2026-08-30 18:00 review,
§5.4 ramp, from `WF-6` step 1 ✅.)*
`examples/ports/06_birdcage_b1_plus_map.py` + same-stem guide (`ports:6`):
two lumped-sheet solves on the 116 085-cell mesh and one combined XDMF
carrying the CG1 `B` phasor (real/imag), `|B₁⁺|` and the phantom cell tags.
Every constant, band and record imported from
`tests/validation/test_birdcage_b1_plus_map.py` (`ANS-1`'s rule), none
restated. Landed `20260831T200401Z_EX-38.log`, Status 0; docrefs
`20260831T200629Z_EX-38-docrefs2.log`, `guide=0`. No homogeneity, absolute
or tuning claim — the map is `WF-6` step 1's gated capability and nothing
more. **Full narrative archived 2026-09-09 → `docs/planning/plan-archive.md`.**

**`EX-39` — the quadrature drive in ParaView: `|B₁⁺|` and `|B₁⁻|` of the
loaded 4-leg birdcage at 10 MHz** ✅ *(commissioned 2026-08-31 03:00 review,
§5.4 ramp, from `WF-6` step 2 ✅; landed 2026-08-31, 16:30 slot.)*
`examples/ports/07_birdcage_b1_quadrature_map.py` + same-stem guide
(`ports:7`): four lumped-sheet solves on the 116 085-cell mesh, the four DG0
curls superposed with `e^{∓jkπ/2}` on the fixture's azimuth index for the co-
and counter-rotating senses. Imports `WF-6` step 2's gated identities
(C4-invariance 0.9818%, mirror 0.8087% against the 5% band, centre
polarisation purity 127.9) and restates none of them. Logs
`20260831T213402Z_EX-39.log` (Status 0),
`20260831T213623Z_EX-39-docrefs.log` (`guide=0`). No homogeneity, absolute or
tuning claim. **Full narrative archived 2026-09-09 → `docs/planning/plan-archive.md`.**

**`EX-40` — the Larmor `|B₁⁺|` ladder in ParaView: 64 and 128 MHz maps of the
loaded 4-leg birdcage** ✅ *(commissioned 2026-08-31 10:30 review, §5.4 ramp,
from `WF-6` step 2b ✅; queued 2026-08-31 18:00 review as §9 item 2, executor
`example-runner`; closed 2026-08-31 21:00 slot, first run green.)*
`20260901T020415Z_EX-40.log` (Status 0, **113 s** at `-n 2`; mesh 24.8 s, six
kept solves ≈ 34 s) plus `20260901T020734Z_EX-40-docrefs2.log`. Both Larmor
rungs on one 116 085-cell mesh (ratio 1.000000, `reused_mesh = True` and
`sweeps[label]["mesh"] is base["mesh"]` asserted), 51 of 51 sample points
valid at each. Gates (i)/(ii) and the P3 control are imported from `WF-6`
step 2b, none restated. No homogeneity, absolute, tuning or convergence claim
at either Larmor frequency. **Full narrative archived 2026-09-09 → `docs/planning/plan-archive.md`.**

### ANS — Ansys benchmark cases (§5.4)

Commissioned by the weekly planning review only, on gated physics only; the
human operator replicates each case in Ansys Electronics Desktop and the
next weekly review adjudicates the returned numbers.

| ID | Title | Status | Tier |
|---|---|---|---|
| `ANS-1` | Loop over a lossy slab at 10 MHz: runnable half of the first AED benchmark | ✅ Closed: no disagreement with the independent solver to diagnose, ΔX stays reported and never gated, and the numeric ruling lives only in the gitignored docs/private adjudication. Any `ans:` case that commits its own `metrics.json` pins that path in the same commit. *History: `docs/planning/chunks/ANS-1.md`.* | standard |
| `ANS-2` | Coil-driven SAR in the loaded four-leg birdcage at 10 MHz: pointwise SAR, whole-phantom power, mass-averaged SAR | ✅ Adjudicated AGREE 2026-09-18 on coil-driven whole-phantom power and pointwise SAR, and the step 4 halved-phantom rung run 2026-09-25 was ruled branch (a): the residual is the feed, and the pointwise band tightens to the measured spread for the weekly to carry into the private file. Step 2 remains an optional record step; no compliance, homogeneity, B₁⁺, Larmor or quadrature claim is licensed. *History: `docs/planning/chunks/ANS-2.md`.* | heavy *(step 1 198 s at `-n 4`; step 4 heavy by prediction, cost probe first)* |
| `ANS-3` | Two coaxial gapped loops at 10 MHz: runnable half of the second AED benchmark (2-port Z/S; `ANS-2` reserved by §10 for the future B1+/SAR case) | ✅ Adjudicated 2026-09-19 by the weekly: AGREE on the primary row, the mutual coupling, while the S rows DISAGREE through a defect on our side, the current-drive route's S table, handed to `PORT-20`. The self-impedance rows are informative only, as pre-declared; nothing else is open. *History: `docs/planning/chunks/ANS-3.md`.* | heavy |
| `ANS-4` | Gapped four-leg birdcage, phantom-loaded, four lumped ports: 4×4 S-matrix at 10 / 64 / 128 MHz | ✅ The self-consistency gates are green on every rung; on the exact SPEC rectangles (117 206 cells), κ off, `PORT-11`'s gates are green at 10 / 64 / 128 MHz at degree 1 and 2 and the degree-1 4×4s are registered (step 5, item 192). Step 6a (items 212 / 225) lands the resolution / conductor-factor / κ knobs and the conductor-halved exact rung at degree 2 with class moves of 0.30–0.60 % against item 192 printed. **The h-statement (F27, 2026-10-12 weekly, `xl` entry 30, `20261011T070008Z_ANS-4-step6c.log`):** the degree-2 exact four-port S is h-converged to ≤ 0.1 % on three rungs (self-class 0.0075 → 0.005 change 0.014 / 0.057 / 0.066 % at 10 / 64 / 128 MHz) and stays within 1 % of AED First Order on every class. Next: F33 / F35 (§9 items 416 / 417, the public verdict sections and h-table). **Palace at matched order** (`OPS-61`, v0.18.1 order 2 on the shared meshes) ties our degree-2 class means to ≤ 2.6e-10 per class on every strip rung bought — ×1 and h = 0.0075 at 10 / 64 / 128 MHz, h = 0.005 at 10 and 64 MHz — and its own h-moves equal ours to the printed digit (`ans:4` `COMPARISON.md` § Palace at matched order). *History: `docs/planning/chunks/ANS-4.md`.* | heavy (≈ 160 s at `-n 2`, one command; measured 125 / 128 s); step 2 **`xl`** |
| `ANS-5` | Pin the element-order correspondence in every ANS `SPEC.md`/`COMPARISON.md` | ✅ Done 2026-08-31: every ANS spec and comparison table now pins our degree 1 to HFSS Zero Order with a First Order sensitivity column, per the 2026-08-30 weekly ruling, and both generators were re-run green. Nothing is open; step 3 was ruled not taken. *History: `docs/planning/chunks/ANS-5.md`.* | smoke (no compute; step 1b ≈ 200 s, measured 61 + 133 s + 125 s control + 1 s census) |
| `ANS-7` | The tuned copper F-small birdcage: `Z_in(f)`, f₀ and Q against AED | ✅ **2026-10-17 weekly (T41).** `f₀` AGREE with its h-statement (T13: within 0.11 % of AED First Order at ×1, 0.03 % at h = 0.005, first order in h). Q and `R_in` converge in h — Cu rows at `C_F` on h = 0.015 / 0.01 / 0.005 m: Q 977.26 / 970.44 / 965.55, `R_in` 18 441.75 / 18 313.60 / 18 223.23 Ω (successive-change ratios 0.717 / 0.705 on steps where first order gives 1.00 and second order 0.60; T30's printed rates p = 1.645 / 1.679; Richardson limits Q∞ ≈ 963.3, `R_in`∞ ≈ 18 182 Ω; brackets in the 10-11 ruling's form [960.7, 965.6] / [18 133, 18 223] Ω; `20261011T094910Z_ANS-7-item364-T30.log`) to a level ≈ 3 % below AED First Order (T13's relative reading, re-used; Palace ties ours to six digits, T43 / T44) — **the absolute Q / `R_in` is INCONCLUSIVE, not an agreement**, pending AED's own pass history (T48); T30's "≤ ½" test is retired (§10 epitaph 2026-10-17). F31 (`xl` entry 34) reads the copper-loss partition the same way (`ANS-6` row). The ESR steps 11–14 (T26 / T27 / T34–T38) and Palace's matched-order driven reads at ×1 and h = 0.01 at both Larmor frequencies (T39–T44) are records, ratified. **Done-when:** ✅ when T30 / F31 write Q's h-statement — met 2026-10-17 (T41; audit PASS 2026-10-18 03:00); until then the row held 🟡 on the Q clause. *History: `docs/planning/chunks/ANS-7.md`.* | heavy |
| `ANS-6` | Copper birdcage, phantom-loaded, four lumped ports: the `ANS-4` fixture with a σ = 5.8e7 S/m coil | ✅ **2026-10-12 weekly (F26 / F19), audit PASS 2026-10-13 03:00.** On exact SPEC rectangles through the hole route `TH-14`'s gates and the copper surface-loss identity (≤ 2.43e-12) are green at 10 / 64 / 128 MHz at degrees 1 and 2, the σ-ladder strictly monotone, class means registered; at h = 0.0075 (226 904 cells, 64 MHz) the Cu − PEC class gap is within 3.2 % of ×1 (`20261012T070009Z_ANS-6-step5b.log:4595`), and the external reading is AGREE within 0.33 % (×1, three frequencies) / 0.21 % (h-rung) of AED First Order. Degree-2 absolute S on this fixture only — no tuning, Q or loss-h claim (F31, `xl` entry 34, a record on the row). **Done-when:** the weekly's adjudication of the exact-rectangle re-read ⇒ ✅ (met). **Palace at matched order** (`OPS-61`, v0.18.1 order 2) ties our degree-2 Cu and PEC class means to ≤ 2.6e-10 per class on the ×1 and h = 0.0075 hole meshes at 64 MHz, its ×1 → 0.0075 h-moves equal to ours to the printed digit (`ans:6` `COMPARISON.md` § Palace at matched order). *History: `docs/planning/chunks/ANS-6.md`.* | heavy |
| `ANS-8` | Human scale: the tuned F-human birdcage, two-port 90° quadrature | ✅ **2026-10-12 weekly on step 1 (entry 29 / H15c), audit PASS 2026-10-13 03:00.** Our side reproduces (v2 rows `20261011T171000Z_ANS-8-v3b-R.log:41–80`; the saved degree-2 σ = 800 32×32 to 4.3e-10, `20261011T123443Z_ANS-8-step1c-W1.log:94–95`); mode content AGREES with AED within 0.6 % (10-03) and degree-2 `C_mode` / `f₀` AGREE with AED First Order within 0.7 % / 0.4 % (one mesh, σ = 800; human-scale tuning quantities quoted at degree 2 only). **Step 2** (the tuned copper fixture with ESR, SPEC v3; our side landed 2026-10-11, `20261011T170627Z_ANS-8-v3b-V3.log`, accepted-power normaliser ruled 10-12) stays open on the row — H22 (§9 item 411, DONE 2026-10-13: `COMPARISON.md` re-stated on SPEC v3's accepted power, SPEC.md "Operator session brief" added, `20261013T094942Z_ANS-8-item411-W1.log`) has taken it to the operator; the AED v3 column is the operator's. SPEC v4 (2026-10-18, §9 item 513): C = 115.5209 pF, R_p ‖ C = 2 317 Ω ‖ 115.52 pF; degree-1 values a labelled record. The lumped-sheet finding (known-issues 2026-10-11) is chain H's (H33 / H34, 10-24 weekly); H35 (§9 item 623) carries row 2. H25 (2026-10-18, §9 item 516, no solve): on entry 36's saved degree-2 copper 32×32, `C_mode`(d2, Cu) = 114.5221 pF beside 118.1186 (degree 1) and the 115.5209 pF ESR root (−0.865 %); circuit-route (Norton feeds, not the in-model row) 2b ccw feed-ring share 0.968088, k = 0 leak 0.025825 (`20261018T095559Z_ANS-8-H25-W1.log`). **Done-when:** the AED half lands and the weekly adjudicates share and leak at a stated relative level ⇒ ✅ (met on step 1). *History (step 1c / v3 readings moved verbatim 2026-10-13): `docs/planning/chunks/ANS-8.md`.* | heavy (our side: H3's 32×32 exists) |
| `ANS-9` | EM–thermal coupling against AEDT HFSS + Ansys thermal (Icepak FEA): SAR from the coil heats the gelled-saline phantom (operator directive 2026-09-30, … | ✅ **2026-10-14 weekly — re-adjudicated on step 1d's column (audit PASS 2026-10-15 03:00 daily review; numbers in the gitignored `docs/private/ans9-readjudication-2026-10-14.md`).** Every row AGREES with AEDT HFSS First Order + IcepakFEA on the production convention (degree 2, exact rectangles, h_p 0.0025 m): ⟨θ⟩ within 0.8 % at every report time, the initial slope within 0.5 %, the surface heat within 1.3 % on the Simpson figure (bands 5 %); probes A–D within 2.5 % with B > D > C held (band 10 %); rows 5 / 7 records (the maximum within 0.5 %, 1.3 mm apart; incident power within 4.4 %, `ANS-2`'s feed question). The 10-10 DISAGREE at D (11 %) was 13 % phantom resolution + 8 % element order on our side (steps 1b / 1c, the 10-12 weekly); known-issues 2026-10-10 retired. Step 2 (P27, two-way σ(T)) is **not posed** as a slot step (§10 epitaph 2026-10-14) — an `xl` probe if the operator requests a two-way comparison. *Earlier:* step 1 adjudicated 2026-10-10 (AGREE on the hand-off, DISAGREE at D on the degree-1 ragged column). **Step 1d landed 2026-10-13 (§9 item 410): `ans:9`'s primary column is degree 2 on exact rectangles, h_p 0.0025 m (route B), rows 1 / 3 / 6 asserted, degree-1 column kept and labelled** — (R1) 1.1e-14, (R3) 5.6446e-6, (R6) 2.2e-15; θ_D 0.756873 K against degree-1 0.697444 (×1.0852), ⟨θ⟩₄₀ ×1.0001; the runnable-half clause is met again on the new column, glyph stays 🟡 (`20261013T093826Z_ANS-9-step1d-W2-h0025.log`). **Done-when:** the runnable half — `ans:9`, SPEC rows 1–7, rows 1 / 3 / 6 asserted at `MAT-15`'s imported bands — ⇒ 🟡; the weekly's adjudication of the operator's run ⇒ ✅. *History: `docs/planning/chunks/ANS-9.md`.* | heavy |
| `ANS-11` | The ASTM-box-loaded F-clinical-1.5T under the two-port quadrature drive at its Larmor root, 64 MHz — `\|B₁⁺\|`, the located 10 g hotspot, the feed quantities (chain C's C7; commissioned 2026-11-04 weekly on entry 52's record) | 🟡 **Done-when (1)'s degree-1 half landed 2026-11-05 (§9 item 833):** the `ans:11` record page parses entry 52's tracked log, imports C6's `_lattice` / `_p_acc_spec`, and its (X) extraction, (D) derived-digit and (R) byte-regeneration anchors are green (`20261105T223815Z_ANS-11-page.log`, 30 PASS, 3 s; W1 `…223748Z_ANS-11-page-W1d.log`). **Done-when (1) met 2026-11-14 (C37, §9 item 983, 12:00 CST slot):** the degree-2 column is parsed from entry 62's tracked log (`20261114T080009Z_GEO-38-C19.log`), (X2) / (D2) / (J62, 53 fields, 0 mismatches) / (R) green, closing window `20261114T181245Z_ANS-11-C37-page-r.log`; the 2 % order rule does not fire on the four per-watt rows (largest |move| 1.26 %, gel fraction; moves also carry the 193.5996 → 185.7629 pF node change). Glyph stays 🟡 for Done-when (2). SPEC v1 written (`examples/ansys_benchmarks/ans11_f_clinical_1p5t_box_loaded_64MHz/SPEC.md`): geometry, materials, 64 ESR capacitor sheets, the two feeds and the pre-registered rows are final; **SPEC v2 (2026-11-07 weekly): the node is C = 185.7629 pF, entry 56's degree-2 root at 0.25 pF** (R_p = 896.1 Ω ‖ C at 64 MHz, Q_C = 66.9) — build and enter C (dashboard W10); the comparison script and the page's node section landed as C36 / C43 (items 871 / 880); **SPEC v2.1 (2026-11-09 weekly): row 3 gains the tie clause** (the location half against the tie set of lattice balls within 1 % of our maximum — the degree-1 column's top two balls are mirror balls 210 mm apart within ≈ 6.3e-04), nothing the operator builds changes. Our degree-1 column is entry 52's record (`20261103T080009Z_GEO-38-C6.log`); **the adjudication column is C19 at degree 2 (`xxl` entry 62)** by the 10-12 epitaph. **Done-when:** (1) the runnable half lands — an `ans:11` page that imports C6's module and never re-implements it, `COMPARISON.md` with our degree-1 and degree-2 rows and the AED columns blank, its closing window `./run_examples.sh`; (2) the weekly adjudicates the operator's First Order column against C19 by the SPEC's rows (rows 2 / 3 / 5 within 10 % / 15 % + one ball radius / 5–10 %) — AGREE / DISAGREE / INCONCLUSIVE public with a relative level, every AED number in the gitignored `aed_results/` and `COMPARISON_private.md` only. No tuning, compliance or homogeneity claim at any outcome. *History and brief: `docs/planning/chunks/ANS-11.md` (the `ans:11` record page, §9 item 833, 2026-11-05 03:00 review).* | heavy (`xxl` column); smoke (the page) |
| `ANS-12` | The ASTM-box-loaded F-clinical-3T under the two-port quadrature drive at its Larmor root, 128 MHz — `\|B₁⁺\|`, the located 10 g hotspot, the feed quantities (chain C; commissioned 2026-11-09 weekly on `xl` entry 58's record, by entry 58's decision rule; the 3 T twin of `ANS-11`) | 🟡 **Done-when (1) met 2026-11-10** (§9 items 930 + 931: the `ans:12` page, and `COMPARISON.md` through the `ans:11` page's `"ans12"` case — degree-1 rows from entry 58, degree-2 / AED blank, row 3's tie set; closing window `20261110T105035Z_ANS-12-C49-ans12-page.log`, (R-cmp) live, 0 [FAIL]); (2) open. **SPEC v1 written 2026-11-09** (`examples/ansys_benchmarks/ans12_f_clinical_3t_box_loaded_128MHz/SPEC.md`): geometry (rung Ø 0.704 m, shield Ø 0.743 × 0.70 m, 32 legs, 64 ring-gap sheets), materials, the two feeds P33 / P41 and the pre-registered rows are final; **the node is TBD** — C_in(F-clinical-3T, d2) at 0.25 pF from `xl` entry 64 (Sunday 11-15, the three-node grid) and its refinement (entry 56's pattern), written as SPEC v2 by the weekly that reads the refinement; until then the operator builds the geometry and does **not** solve. Our degree-1 column is entry 58's record (`20261109T080008Z_GEO-38-C6-3T.log`, at the degree-1 root 44.0525 pF — a record, never adjudicated); **the adjudication column is the 3 T degree-2 set (`xxl` entry 66, PENDING)** by the 10-12 epitaph. **Done-when:** (1) the runnable half lands — an `ans:12` record page that imports C21's module (`tests/validation/test_geo38_f_clinical_3t_box_loaded.py`) and never re-implements it, parses entry 58's tracked log through C33's parser with (X) extraction, (D) derived-digit and (R) byte-regeneration anchors green, and `COMPARISON.md` with our degree-1 rows, the degree-2 and AED columns blank, the tie set of row 3 listed; closing window `./run_examples.sh` (§10 C48 / C49); (2) the weekly adjudicates the operator's First Order column against entry 66 by the SPEC's rows (rows 2 / 3 / 5 within 10 % / 15 % + one ball radius or the tie set / 5–10 %) — AGREE / DISAGREE / INCONCLUSIVE public with a relative level, every AED number in the gitignored `aed_results/` and `COMPARISON_private.md` only. No tuning, compliance or homogeneity claim at any outcome; the gel's 128 MHz σ / ε_r stay the provisional W3 ruling until the F2182 clause is in the repo. *History and brief: `docs/planning/chunks/ANS-12.md`.* | heavy (`xxl` column); smoke (the page) |
| `ANS-13` | The ASTM-box-loaded F-clinical birdcages' gel heating at 64 and 128 MHz against HFSS + IcepakFEA, one-way, per watt deposited in the box — θ(t\*) at c\* and six probes, the calorimetric slope, the surface heat at t\* (chains P / C; commissioned 2026-11-11 weekly on `xl` entries 59 / 60's records, by their decision rules; `ANS-9`'s pattern at clinical scale) | 🟡 **SPEC v1 written 2026-11-11** (`examples/ansys_benchmarks/ans13_f_clinical_box_heating_64_128MHz/SPEC.md`): the EM half is `ANS-11` / `ANS-12` unchanged (the 3 T node TBD with `ANS-12` v2), the thermal half is `MAT-15`'s route on the box (k 0.6, ρc 4.18e6, h 10 to 20 °C on six faces, t\* 2821.5 s), seven probes, the source fixed at 1 W deposited; our degree-1 columns are entries 59 / 60's records (`20261110T080008Z_GEO-38-P33.log`, `20261111T080008Z_GEO-38-P33-3T.log`); **the adjudication columns are the degree-2 twins (`xxl` entries 69 / 70, PENDING on P87)** by the 10-12 epitaph. **The `ans:13` record page landed 2026-11-12** (04:30 slot, §9 item 958 / §10 C53; closing `20261112T103734Z_ANS-13-page.log`, 165 [PASS] / 0 [FAIL]: (X) 52 of 52 SPEC record-table tokens, (D-dep) 24, (D-c\*), (D-ratio) 22 C39 strings, (D-lever), (D-cal), (R) `SPEC.md:140–148` byte-identical; (R-cmp) absent until C54). **Done-when (1) met 2026-11-14** (13:30 slot, §9 item 984 / §10 C54): `COMPARISON.md` with the two degree-1 columns, degree-2 and AED blank, (K) `ANS-11` / `ANS-12` byte-identical; closing `20261114T193413Z_ANS-13-C54b-ans13-page.log`, (R-cmp) live, 166 [PASS]. Glyph stays 🟡 on (2). **Done-when:** (1) the runnable half lands — an `ans:13` record page that imports P61's and C31's modules and C46's reader and never re-implements them, parses entries 59 / 60's tracked logs with (X) extraction, (D) derived-digit (per-deposited-watt quotients, C39's ratios by string) and (R) byte-regeneration anchors green, and `COMPARISON.md` with our degree-1 rows, the degree-2 and AED columns blank; closing window `./run_examples.sh` (§10 C53 / C54); (2) the weekly adjudicates the operator's First Order + IcepakFEA columns against entries 69 / 70 by the SPEC's rows (rows 3 / 4 within 5 %; row 2 within 10 % at seven probes with the wall-probe ordering preserved) — AGREE / DISAGREE / INCONCLUSIVE public with a relative level, every AED number in the gitignored `aed_results/` and `COMPARISON_private.md` only. No compliance, tissue, dose-per-watt or container-wall claim at any outcome. *History and brief: `docs/planning/chunks/ANS-13.md`.* | heavy (`xxl` columns); smoke (the page) |

**`ANS-5` — pin the element-order correspondence in the benchmark specs** ✅
*(ruled 2026-08-30 weekly review, steps 1–2 executed 2026-08-31; full
245-line entry archived in `docs/planning/plan-archive.md`, 2026-09-06
weekly review).* Result: our production `degree 1` (N1curl, 6 unknowns/tet)
is HFSS **Zero Order**, the adjudication column; AED's default **First
Order** (20/tet) is the order-sensitivity column; Mixed Order forbidden;
`ANS-1` in scope for the formulation/order line. README correspondence
table and a *Basis / element order* paragraph in all three SPECs; the
`ANS-1`/`ANS-3` `COMPARISON.md` rows regenerated from their generators
(step 1b). The 6 / 20 correspondence was **confirmed by HFSS's own matrix
sizes on 2026-09-04** (`ANS-4` row). Live carry-forward: none. Logs:
`20260831T110240Z_ANS-5-step1b-ans1.log`,
`20260831T110908Z_ANS-5-step1b-ans3-final.log`,
`20260831T111136Z_ANS-5-step1b-census.log`.

**Operator ruling 2026-10-01 (interactive session) — AED benchmark default:
First Order, copper / PEC coils; for the weekly to carry out.** It amends the
rule above for **new** AED runs; columns already run are unaffected.
(1) **First Order is the adjudication column** for new cases and new columns.
**Zero Order is run where it makes sense, by HFSS's own guidance** (operator
clarification, same day): on the AED side, First Order by default and Zero
Order allowed where the help page says it fits ("useful when a model requires
a mesh with more than 100,000 tetrahedra, but the model size is small compared
to wavelength"). A SPEC that runs Zero Order states why. Until `T6` passes
(degree 2 with the κ-corrected width becomes the order for lumped-sheet S,
§10 chain F), a case that adjudicates a figure our degree-1 production path
produces also qualifies. *Note for the weekly:* every coil case to date meets
HFSS's criterion (64 MHz λ ≈ 4.7 m in air against a ≈ 0.3 m coil, ≥ 0.5 M
tets set by wires, gaps and skin depth). So it is a cost argument, and the
SPEC must say what the Zero Order column is compared against. `ANS-8` v2's
AED Zero vs First Order columns on one fixture are the first measurement of
how order-sensitive the adjudicated rows are. Mixed Order stays forbidden
until the weekly rules otherwise.
(1b) **Our solver, long-term direction (operator, same day): degree 2 by
default; degree 1 allowed where a convergence ladder licenses it for the
reported quantity, the licence stated per case.** This is HFSS's rule
translated, not copied. HFSS's Zero Order is safe because adaptive refinement
drives ΔS to convergence; we solve on fixed meshes, so electrical size alone
cannot license degree 1 here. Only a ladder on the quantity reported can.
The criterion also fails inside tissue: εᵣ ≈ 78 gives λ ≈ 0.5 m at 64 MHz and
≈ 0.26 m at 128 MHz, against a 0.3–0.4 m phantom. The flip to degree 2 as the
production default stays gated by `T6`; this direction does not move that
gate. The costs it carries forward: degree 2 at human scale is `xxl`-tier
today (F-human 32×32 ≈ 3.3 M unknowns, ≈ 103 GiB) unless a degree-2
iterative solver is made to work; there is no DOLFINx analogue of Mixed
Order (no variable polynomial degree across a mesh); and the lumped-sheet κ
is order-dependent (`PORT-21` step 5, T11b).
**Long-term goal (Phase 9 era): adaptive meshing** — an a-posteriori error
estimator converging on a stated quantity (S, `C_mode`, or a SAR / B₁⁺
functional). It is what would let degree 1 on an adapted mesh be licensed the
way HFSS licenses Zero Order, and so let this solver follow HFSS's rule
literally. It matters most when Phase 9's external-CAD (STEP) path brings
geometries nobody will ladder by hand. **Route (operator, same day): adaptive
remeshing through gmsh, not DOLFINx `refine`.** `dolfinx.mesh.refine` (0.11)
bisects edges with Plaza closure but places new boundary vertices on the old
straight edges. Curved coil tubes, tori and phantom walls therefore stay at
their first-pass faceting however far it refines, and those are the regions
an estimator refines most. Remeshing regenerates from the OCC geometry each
pass, so curved surfaces are re-resolved. Physical groups live on geometry
entities, so coil / phantom / port-sheet tags survive with no parent-cell
bookkeeping. The gmsh side exists in the image (4.15.2, checked 2026-10-01):
per-element target sizes as a view (`gmsh.view.addModelData`) feeding a
`PostView` field, or `gmsh.model.mesh.setSizeCallback`, merged into the
generators' existing `Min` of `Distance` / `Threshold` / `Box` / `Ball`
fields so a pass can only refine past a known-good start. **To build:**
(1) the estimator in DOLFINx (a residual estimator for curl-curl; none is
provided — the real cost); (2) the bridge: gather η per cell to rank 0 via
`original_cell_index`, target sizes `h_new = h_old·(η_target/η)^(1/p)`,
remesh, redistribute through `model_to_mesh`; (3) gradation and floor /
ceiling clamps, isotropic first; (4) a record form for adapted meshes
(converged quantity and final-step change, not cell counts). Per pass: one
remesh (F-human build ≈ 110 s) plus one solve. A cheaper first step that is
useful on its own: **the estimator without refinement**, reporting where η
concentrates to check that a ladder refines the right region. `refine`
keeps a role for quick uniform checks on polyhedral geometry. The weekly
places the goal on the §6 phase map; it is not a near-term chunk.
(2) **Coils are copper (Finite
Conductivity σ = 5.8e7, Solve Inside off) with PEC as the control column.**
**No new σ = 800 S/m coil columns.** Existing σ = 800 AED data (`ANS-4`, the
`ANS-8` v1 column) stays in use for the chains already reading it (chain F /
`PORT-21`, H8b). *Why:* order-matched agreement at Zero Order has proved
fragile. `ANS-4`'s 10 MHz degree-1 agreement was re-read as partly
coincidental. The clean agreements are our degree 2 against AED First Order
(`ANS-6` step 2c, `ANS-4` at 128 MHz, the κ hypothesis). AED's adaptive Zero
Order mesh never matched our discretisation anyway. σ = 800 is an artifact of
our skin-depth mesh rule, on which the 2026-09-04 directive already lets no
parity claim stand. *What it costs until our side catches up* (the reasons
the operator weighed): our production order is still degree 1 (`T6` not
met); our degree-2 lumped-port S carries the κ confound until the in-run κ
(`PORT-21` step 5) is checked by T11b; and copper at human scale has no
our-side counterpart until `WF-12` unblocks or `GEO-33` lands. New
First Order / copper columns are therefore **reference data until those
land**. The ruling raises `T6`, T11b and `WF-12` in priority; it does not
change their gates. *For the weekly:* amend this entry's rule and the ANS
README correspondence table and SPEC template; set ANS-9's SPEC (and
`ANS-8` step 2's) to the new default; and state which open SPECs keep a Zero
Order or σ = 800 column under the exception. `ANS-8` step 1's v2 run
(started 2026-10-01 with copper Zero Order, σ = 800 First Order and copper
First Order) predates the ruling and stands as commissioned.

**`ANS-1` ✅ 2026-08-09** *(scoped 2026-08-09, weekly review; full plan and
closure narrative in `docs/planning/plan-archive.md`)*. Runnable half of
the first AED benchmark, dispatched through the runner's `ans:` group
(`./run_examples.sh -e ans:1 -n 2 -t 180`, log
`20260809T183731Z_ANS-1.log`, 70 s): ΔR = +3.2770406e-01 Ω, **1.5834%**
from Dodd–Deeds against the 2% ceiling and 1.387e-08 relative from the
`MAT-6` pin against 1e-3; σ = 0 control exact (0.0 W, 0.0 A/m², asserted
with no tolerance); energy identity ratio 1.0000. `metrics.json`,
`COMPARISON.md` (AED columns blank per SPEC), and the |J| XDMF landed in
the case directory; every constant, mesh, and drive is imported from the
`MAT-6`/`EX-11` modules, so the benchmark cannot drift from the gate. The
AED half is the operator's (§5.4 Waiting-on-you).

**`ANS-3` ✅ 2026-08-16** — two coaxial gapped loops at 10 MHz: runnable half *(the execution record and the original commissioning plan are archived verbatim in `docs/planning/plan-archive.md`, entry «§7 ANS-3 full narrative — archived 2026-09-13 (weekly review)».)*
> * **Done-when, met** (`20260816T110354Z_ANS-3-runnable-half-n2.log`, **131 s** wall at `-n 2`, 178 055 cells — mesh 35.9 s, sweep 46.3 s, export solve 21.4 s; `./run_examples.sh -e ans:3 -n 2 -t 500`): every `EX-20` anchor reproduced inside the pre-stated 1% band, misses ≤ **3.67e-06** — raw mutual 0.894543, corrected 0.939849, ‖S−Sᵀ‖/‖S‖ = 2.5494e-05, ‖S‖₂ = 0.861449; Im Z₂₁ = +1.110803269e+00 Ω vs ωM₁₂ = 1.241755 Ω, |Z₁₂−Z₂₁|/|Z₂₁| = 5.8309e-04 (reported, not gated). **Negative control printed first:** the raw rung −10.55% asserted to *fail* the unmoved 10% band; the corrected rung −6.02% inside.
> * **Artifacts:** `metrics.json` (complex 2×2 Z and S, ladder, identities, mesh/timing), `COMPARISON.md` (our columns filled, AED columns blank; Zero / First Order columns added by `ANS-5` step 1b), combined XDMF; every constant imported from `examples/ports/02_package_sparameter_sweep.py` (`EX-20`) and `ports/systematics.py` — `ANS-1`'s rule. Incidental fix landed with it: `scripts/run_examples.sh` issues `timeout -k 30` inside the container.
> * **Carry-forward:** the AED half is the operator's (dashboard Waiting-on-you, behind `ANS-2`) and is also `PORT-10`'s independent adjudication input; re-run green 2026-08-31 by `EX-37` (128 s), the 08-16 / 08-26 records standing.

**`ANS-4` — loaded birdcage 4-port S-matrix at 10 / 64 / 128 MHz: runnable half** ✅ **2026-08-30** *(the commissioning text and the 16:30 execution journal are archived verbatim in `docs/planning/plan-archive.md`, entry «§7 ANS-4 full narrative — archived 2026-09-13 (weekly review)»; the step-2 ladder history and the adjudications are in `docs/planning/chunks/ANS-4.md`.)*
> * **Done-when, met** (`20260830T213415Z_ANS-4-run1.log`, Status 0, **125 s** at `-n 2`, `./run_examples.sh -e ans:4 -n 2 -t 500`, heavy): every gate on every rung of one 116 085-cell mesh (ratio 1.000000, `reused_mesh` asserted) — reciprocity 1.469e-14 / 1.126e-15 / 8.763e-16, σ_max 0.999992805 / 0.999721388 / 0.998974779, C4 spreads 0.0553 / 0.0353 / 0.0214%, 0.0573 / 0.0599 / 0.0370%, 0.1012 / 0.0916 / 0.0654% (pooled separations 166.7× / 671.1× / 576.9× vs the 10× floor); 128 MHz cells/λ **12.5024** ≥ 10.
> * **Records reproduced:** 10 MHz leg (d)'s 4×4 to **1.158e-10** (band 1e-6) and leg (d0)'s column to 2.568e-10 (1e-9); 64 / 128 MHz worst 1.075e-03 / 6.755e-04 against `PORT-11` steps 2/3 (band 1e-2). **Negative control, printed first:** the retired `PORT-0` heuristic at 128 MHz, max|off-diagonal| 0.000000e+00, separation **1.585460** vs the 2e-3 floor.
> * **Artifacts:** `04_birdcage_four_port_10_64_128MHz.py` + same-stem guide, `metrics.json`, `COMPARISON.md` (two blank AED columns, Zero and First Order, per `ANS-5`), `paraview_output/ans4_birdcage_four_port_128mhz_combined.xdmf` (P1-driven at 128 MHz); census `20260830T213718Z_ANS-4-docrefs2.log` `guide=0`. Everything is imported from the `ports:5` / `PORT-9` / `PORT-11` modules, so the benchmark cannot drift from the gate.
> * **Carry-forward:** self-consistency identities on one fixture; the AED adjudication lives on the row, in `docs/planning/chunks/ANS-4.md` and in gitignored `docs/private/` — AGREE at 10 MHz (2026-09-06); the Larmor verdict is the 2026-09-13 weekly's (§2.2, §6, §10).

---


**`ANS-6` — copper birdcage, the `ANS-4` fixture on a realistic conductor**
⬜ *(commissioned 2026-09-04 by operator directive, interactive session — an
exception to "commissioned by the weekly planning review only", recorded
here as such; the weekly review owns its spec text and adjudication as for
every other case. Serial on `TH-15` ✅ and `TH-14` ✅.)* **Why.** `ANS-4`
is the absolute check of the port model at the Larmor frequencies, but on
an 800 S/m coil that no real birdcage resembles; its coil-loss share is
wrong by construction and every tuning-workflow parity claim (§1 step 2)
needs the conductor model checked against AED *before* Phase 6 leans on it.
One case, same geometry, one thing moved: the coil conductivity.
> * **Step 1 (spec).** `examples/ansys_benchmarks/birdcage_copper_four_port/SPEC.md`
>   by copying `ANS-4`'s spec with these changes and no others: coil
>   σ = **5.8e7 S/m** (copper, μᵣ = 1), *Solve Inside off*, HFSS **Finite
>   Conductivity** boundary on every coil face; a **second AED column with
>   the coil as PEC** (HFSS *Perfect E* on the same faces) matching `TH-15`;
>   the `ANS-5` order rule unchanged (Zero + First Order, Mixed forbidden);
>   quantities as `ANS-4` plus the **coil-loss share `P_coil/P_accepted`**
>   per frequency (HFSS: surface loss on the finite-conductivity boundary
>   from the field calculator) — the row that discriminates the conductor
>   model. Our side: `TH-14` (copper) and `TH-15` (PEC) columns.
> * **Step 2 (runnable half).** As `ANS-4`'s script, importing every band
>   from the `TH-14`/`TH-15` gate modules; `metrics.json`, `COMPARISON.md`
>   (AED columns blank by construction) and the private-mode writer
>   (`OPS-32`'s pattern, `aed_results/ans6_aed_results.json`). Heavy tier.
> * **Step 3 (operator replication).** The `ANS-4` PyAEDT script (untracked,
>   `aed/`) adapted: material copper, `solve_inside = False`,
>   `assign_finite_conductivity` on the coil faces, a second design with
>   `assign_perfecte_to_sheets` on them. Goes to the dashboard's
>   Waiting-on-you when step 2 closes. Numbers private, verdict public, as
>   always.
> * **Done-when (§4).** Step 2's harness log with the imported gates
>   asserted; adjudication by the first weekly review after the AED numbers
>   land — AGREE on the S classes **and** on the coil-loss share is what
>   licenses §2.1 to say "conductor model checked externally".

## 8. Legacy ID mapping

Commit messages, `docs/testing/logs/*.log`, and the retired
`docs/testing/pending-tests.md` (removed 2026-08-04; in git history) use older
IDs. Two generations collided — `E1`–`E4` refer to *different chunks* in
the ROADMAP than in `pending-tests.md`. Resolve via this table.

| Legacy (ROADMAP gen-2) | New ID | | Legacy (gen-1, in `pending-tests.md`) | New ID |
|---|---|---|---|---|
| A1 | `MAG-6` | | C1 (coil+phantom B-solve) | `MAG-1` |
| A2 | `OPS-3` | | C2 (sanity validation metrics) | `MAG-6` |
| A3 | `OPS-4` | | D1 (freq-domain scaffold) | `TH-2` |
| A4 | `GEO-7` | | D2 (phantom material MVP) | `MAT-1` |
| A5 | `OPS-5` | | D3 (E/B extraction) | `POST-1` |
| B1–B6 | `GEO-1`…`GEO-6` | | E1 (port data model) | `PORT-2` |
| C1 | `TH-2` | | E2 (birdcage port tags) | `GEO-1` |
| C2 | `MAT-1` | | E3 (port excitation hook) | `PORT-4` |
| C3 | `TH-3` | | E4 (N-port S-assembly) | `PORT-5` |
| C4 | `POST-1` | | E5 (Touchstone export) | `PORT-7` |
| C5 | `POST-2` | | E6 (calibration checklist) | `PORT-3` |
| C6 | `TH-4` | | F1 (run-and-log metadata) | `OPS-6` |
| D1–D6 | `PORT-3`…`PORT-8` | | F2 (manual checklist doc) | `OPS-7` |
| E1–E4 | `WF-1`…`WF-4` | | | |
| F1–F3 | `OPS-6`…`OPS-8` | | | |

---

## 9. Immediate sequencing

Phase-1/2 analytic gates are closed (§6); the working front is ports
beyond the two-torus fixture and the Larmor-regime validation gate.

*(§9 head items 1–6 — the history pointers, the operator's Ansys **feature ladder** (item 5, ranks A1–C4) and the XL-tier pointer — moved verbatim 2026-10-24 to `docs/automation/daily-review.md` § "Standing §9 rules" (operator ruling, interactive session: §9 holds the queue; the standing blocks live with the protocol that applies them). Citations such as "§9 item 5" and "ladder B1" keep their meaning.)*

**Standing rules.** Do not add new features to `⚠️` subsystems. Do not
trust a chunk's status without a log — any §7 status that is not `✅`
reads "unknown", not "probably fine". Since `OPS-19` (2026-08-16), the
docrefs checker exits 0/1/2 (clean / hard violation / staleness-only) and
prints a machine-readable `RESULT:` line — chunks that run examples gate
on **`exit != 1`**, and read staleness (`exit 2`) as information, not
failure.

### On deck — maintained by the scheduled daily review

The next scheduled implementer run takes the **first** item below that is not
marked done or blocked (see `docs/automation/implementer-run.md`). **Restock
floor (operator directive 2026-10-01; was ≥ 240 *predicted* minutes from
2026-09-13): ≥ 160 measured slot-minutes of open, unblocked work and ≥ 5
items** — each item's class median from the slots' own commit times, running
total written beside the list, shortfall stated in minutes when the floor is
not met (daily-review.md step 6). **Item briefs live in
`docs/planning/chunks/<CHUNK>.md` (`### §9 item N`); this section carries
two-line entries** (operator 2026-10-01; open items migrate at the first review
after the directive). A brief may pre-register `PREDICTED RED: <assertion>
range [a, b]`; a red inside its range with all else green parks the item and
the slot continues (implementer-run.md step 2). Items are ordered, each sized
for one run: ≤ 1 h wall clock, ≤ 20 min per compute command. Prefer items
that do not depend on each other; where the critical path is genuinely serial,
say so in the item. If every item is done or blocked, the drain instruction at the
end of this section applies: **stop and journal**.

**Standing review rules — moved 2026-10-24 (operator ruling, interactive session).** Rules (1)–(4) of 2026-09-12 (take-next, the four-attempt family cap, the status an item can move, the enumerated-chain licence), rule (5) (the backlog lanes), rule (6) (the chain-extension licence, **six** provisional steps per chain since 2026-10-24) and the two binding rulings live verbatim in `docs/automation/daily-review.md` § "Standing §9 rules"; they keep their numbers, so "§9 rule (N)" still cites them. **The binding rulings in one line each, for every slot:** port / lumped-element sheets are exact CAD rectangles and κ is not applied to them (2026-10-02); every lumped capacitor carries a 0.2 Ω series ESR, entered in HFSS as the parallel R_p = r(1 + Q²) at a named frequency (2026-10-06).

Last reviewed **2026-11-16, 03:00 review** (Monday; `claude-opus-5-5`, the standing daily choice). Full journal: `docs/planning/reviews/2026-11-16-daily.md` (previous: `…/2026-11-15-daily-1800.md`). These rulings bind the next slots:
- **A test that reads a tracked log reads it through `tests/log_records.read_log_text(…, errors="strict")` (or `resolve_log_path`), never `open()` / `read_text()` on `docs/testing/logs/<name>.log`** (18:00 11-15, provisional until the 11-16 weekly). `housekeeping.py` gzips every log at 7 days, cited or not (`:138–143`); citation exempts a log from deletion only. A plain reader goes red on the first sweep after its log turns 7 days old. Item 997 (19:30 11-15) routed the known readers through it and retired known-issues 2026-11-15. Every new brief that reads a log states this.
- **The 11-14 weekly ruled entry 63's branch (`PROJECT_PLAN.md` §10, "Rulings, 2026-11-14").** (S2-tip) takes its pre-registered branch: on the clinical body the tip read is printed, and `S2_BAND` never moves. The (E) red at Δt₁/10 is a floor reading. The asserted (E) returns to entry 55's definition (CN and the Δt₁ BE step), and the Δt₁/10 residual is printed: that is **I32c (item 986)**. Entry 63 is re-run as **`xl` entry 71** once I32c lands, and entry 67 rides behind 71. Known-issues 2026-11-06 stays open until entry 71 reads green. I76 is re-posed onto entry 63's log as **I76′ (987)**. U153 and 70′ stay serial on entry 71 green. No slot re-runs entry 63 or widens `E_BAND` / `S2_BAND`.
- **The after-window licence, extended by the 11-14 weekly.** A window read RED mechanically opens nothing until the weekly rules its branch. The weekly may then re-pose a step onto the red log where that step's inputs are untouched by the red gates (I76′).
- **A generator constant that a landed test pins is never extended by a brief; new keys live in their own dict (`TAG403_15T`). Before extending a generator constant, a brief greps `gen.<NAME>` across `tests/`** (03:00 11-14, ratified by the 11-14 weekly, now standing).
- **A brief's "≤ 1e-12" parenthetical on a one-step thermal identity is overruled by the 11-10 rulings** (10:30 11-13, ratified, standing). Item 969's K10-es (E) moved 1.521e-10 → 1.522e-10, a reading inside its band.
- **A reader requires exactly one distinct match, with the count printed, and never takes "the first"** (10:30 11-12, ratified, standing; `MAT-21.md` and `WF-17.md` § "Review ruling 2026-11-12 10:30"). A `-n 2` log prints its tables once per rank, and C45 prints a self-copy control table before the real one. Custody compares values (`float(a) == float(b)`), and a print-form-only difference is printed, not red. A reader premise that is wrong about a print count is the review's to rule, never the slot's.
- **The after-window licence (11-11 weekly).** A step written "serial on entry N green" is takeable from the first slot after a review has read the window **mechanically** (Status 0, `[capture] rc=0`, orphans 0 / 0, every asserted gate green by the `log-pathologist` reading in the ledger row) and written the row. A step whose *branch* depends on the decision rule still waits for the weekly. **Entry 64 (C5f) ran 02:00 11-15** (Status 0, `[capture] rc=0`, 8 158 s, orphans 0 / 0, `20261115T080008Z_GEO-38-C5f.log`); the ledger row carries the mechanical read (`xl-ledger.md:84`), so U166 is queued (item 1001, serial on 998).
- **The module hold is the window's import closure, read forward** (03:00 11-12, ratified by the 11-14 weekly, standing). A queued window's own module, every module that window imports, and its own module's importers are held until it runs. A module is **not** held merely because it imports *another* window's module that the window does not import back. Importing a held module is always allowed.
- **A solve output's print is gated against its named constant (18:00 11-10; ratified by the 11-11 weekly, standing; `WF-17.md` § "Review ruling (2026-11-10 18:00 daily — a solve output's print").** Where a node names a module constant for a solve output and asserts it only through a band, the constant is the record and the print is a reading of it. The print is shown beside both the constant and the last run's print. The 1-unit / 2-unit rule binds only (b) prints that have no named constant. The degree-2 sphere SAR's `.10e` print scatters 4–10 units between same-commit runs, and the degree-1 print scatters 1.99e-12. Item 943 closed on it.
- **A one-step slope identity at `IDENTITY_BAND` is below the absolute-T floor (18:00 11-10; ratified 11-11).** (E) is 1.2e-10 on one BE step at Δt₁/10. No brief poses a 1e-12 identity on a one-step thermal increment; U152′ (965) prints (Y-S2) / (Y-Π) as readings.
- **"Every record at its prior" is ruled (10:30 11-10; ratified by the 11-11 weekly, standing; `WF-17.md` § "Review ruling (2026-11-10 10:30 daily").** A **record** is a value the node asserts by equality, or a named print at ≤ 11 significant digits. A record off its prior by one unit in its last place is printed beside it and the slot continues; by two or more units it is BLOCKED. A **reading** is printed beside its prior and gated by its band only; a reading moved inside its band never blocks. Readings are identity residuals, distances-to-record (`rel …`) and prints finer than 11 digits. Precedent: item 928. Items 938–940 stopped on readings, so their known-issues entries are retired and the steps continue as 943–945. The ruling binds every chain U re-run brief (941–948).
- **A queued window's module is not touched before that window runs** (11-09 weekly ruling (c)), and a default-off knob does not lift it. Measured by import closure at this review (11-11 ruling (d), ratified):
  - **held:** entry 64's `test_geo38_f_clinical_3t_esr_root.py` (running 02:00 11-15), entry 65's `test_geo38_f_clinical_matrix.py` (11-16), and `xxl` entry 69's `test_geo38_f_clinical_box_heating.py` (11-28), each with its import closure;
  - **free:** entry 63's `test_geo36_rod_f_clinical_heating.py` (ran 11-13; I32c edits it, and entry 71 is PENDING, not queued), and entry 62's modules (ran 11-14);
  - each brief below states its own import-closure check.
  A slot confirms the window's `chore(xl)` commit is on `main` and its queue file is gone, or it skips the item.
- **Never delete a log, red or not.** Item 907's executor discarded a red W0 and its index row (`GEO-38.md` § item 907, review note). The record stands because the committed windows ran on the committed text, but the next deletion parks the item.
- **The `xl` readers' "(Pend)" premise tests are tree-independent** since item 918 (04:30 11-09; known-issues 2026-11-09 retired). Entries 59–61 landing no longer turns a node red.
- **Every new step's reader of a pending window parses its tracked log, not the host JSON** (`/logs/` is gitignored; P61 and C31 write no JSON). Items 893, 898 and 905–907 apply it; the host JSON, where present, is only a custody check. This is a route reading for the weekly, not a target change.
- **Never `--no-verify` past the leak hook, whatever the match looks like.** A record the hook refuses has its log moved unchanged to `logs/leak-held/` and the rest committed through the hook. **Since the 11-07 weekly a slot applies this itself on the first encounter of a launcher's refused record** (`implementer-run.md` step 1): the journal names the hook's line by number only, never the token.
- **The tip calorimetric-slope gate (S2-tip) is ruled (11-07 weekly, option (b)):** it is asserted against the consistent-mass projection ⟨Πq⟩_B at the unmoved `S2_BAND`, on one BE step at Δt₁/10, with ⟨q⟩_B and the DG0 ratio printed. I32b (item 890, re-posed from 870) lands it and U127′ (item 884) carries it into the `Study`. Known-issues 2026-11-06 resolves on 870 + 884 + `xl` entry 63. **Entry 63 is red (03:00 11-13): the rod's R_Π at Δt₁/10 reads 1.092e-04 / 1.104e-04, scaling O(Δt) from Δt₁ — the entry stays open for the 11-14 weekly.** No slot widens `S2_BAND`.
- **A §9 DONE label and a chunk-file Outcome name the local (CST) slot, never the UTC hour** (items 850, 855 and 859 carried the UTC hour; corrected).
- **Rule (6) since the 11-14 weekly (every chain reset to 0; U159′ / U160 / U161 / U162 ratified): I, C, P and U each 0 of 6.** **03:00 11-16: P92 (item 1008) is chain P's first provisional step — P 1 of 6; I, C, U 0 of 6.** Items 986–996 are the 11-14 weekly's own enumerated steps (rule (4)) or lane targets, not provisional. **U163's layout is ruled** (two one-line premise edits in landed nodes, licensed by the weekly; item 992). **U13 is re-posed** as U13a (item 993) then U13b (serial on U13a; not yet queued). **18:00 11-15 re-poses (rule (6)'s re-pose clause; family counts kept; no provisional step):** 993 → **1000** (ruling (a), whole port boxes; option (b), a `src/` lineage for a face shared by two volumes, goes to the weekly); 991 → **999** (the s₀ time split; (N-h) on N = (40, 40, 40); no `src/` edit); 992 re-runs as **998** after 997.
- **An example guide quotes closing-log digits only, and is not edited after its closing window** (`EX-183`'s figure-window digit struck at this review; `EX-184`'s post-closing guide edit accepted with the caveat). If a guide must change, re-run the closing window. A rounded bound on a closing-log digit ("≤ 2e-14" over 1.499e-14) meets the rule; a digit from another window does not (`EX-191` audit, 18:00). **A census count names its basis:** the docrefs `runnable` line, not the `ns:N ->` count, which reads 6 lower (`EX-192` audit, 18:00).
- **`OPS-82` has landed:** `_rod_source` / `_tips_source` build in both builds. A complex-build refusal must never be asserted through `fem.form` inside `pytest.raises`: `ArityMismatch` is a `BaseException` raised on rank 0 alone, and the peer hangs to the ceiling (item 849's first W2).
- **The rule (e) reminder applies to every brief:** a control factor that no log backs is printed as *predicted*, or asserted at the gate's own band.
- **A refresh window that goes red is a finding: stop, file it, fix nothing.** `EX-179`'s run-line-less examples run at `-t 480`; a kill there is filed, never raised.
- **A chunk file cites every log by its full basename.** `housekeeping.py`'s `CITATION` regex does not match `…193432Z_EX-175-mat4`, so an abbreviated log is non-gating and is swept (18:00 journal §3).
- **Examples close on `./scripts/run_examples.sh -e <id> -n <n> -t <s>` through the harness** (10-31 weekly). **The form is now observed to run, no `-T` needed** (`OPS-74`, and `EX-170` … `EX-172` closed on it, audited clean). A closure whose last window is anything else is audited with the caveat. The flagged `FEM_EM_SETUP_FIGURES=1` window is a separate earlier window that closes nothing. **A guide quotes a round-off term as a bound, or names the window its digits came from** (03:00 ruling, `EX-171`'s disclosure).
- **A root grid's reading is the bracketing pair's interpolation,** never a multi-point fit (10-31 weekly, entry 51). Each node's printed forms go through the leak guard before a queue file.
- **A crash in a module's own non-gate code before any gate value exists may be fixed and re-run once at the same `-t`,** first log kept and disclosed (items 678 and 700, ratified narrowly by the 10-31 weekly). "A kill parks the item" governs size and timeout kills.
- **`main` has no known red test node in an ordinary window since item 997 (19:30 11-15, `464e1207`, `20261116T013304Z_WF-17-i997-W1.log`)**, which retired known-issues 2026-11-15 (the 11-15 01:45 sweep's gzip reds: two measured 1.5 T nodes, four inferred `MAT-21` nodes). *Before that sweep:* `main` had no known red test node in an ordinary window since item 983 (12:00 11-14) retired `test_ans11_record_page.py::test_c19_absent`** (I85's `test_i85_k1_12_knob_parsing` was retired by item 978, 07:30 11-14) (C44's two premise nodes, red since `3b5b7e1f`, retired by item 918 at 04:30 11-09) (the knob-gated body `test_geo36_i32_rod_heating` is red in its `xl` window, known-issues 2026-11-06) since item 845 (19:30 slot, `20261106T013115Z_ANS-11-spec-W1.log`: the two `test_doc_reference_exit_codes.py` nodes green, docrefs `dead=0`). Item 841 retired the `test_wf14_validation_page.py` entry (`20261105T180811Z_WF-14-item841-W2.log`). **Never `gunzip` a record to get green.**
- **For each brief, the newest dated "Review ruling" paragraph in its chunk file binds over the brief.** A brief that adds a tagged block to a page with an exact fence pin says so in its own text (item 737).
- **Path checkouts take code only.** Record rows (test-results) are appended by hand, and `git diff --stat` on record files shows insertions only.
- **The environment-first ruling stays lifted** (10:30, after `OPS-72`). A window may lead with the full `tests/environment` (≈ 45 s), or keep `test_complex_mode.py` first plus one `W-env`.
- **Count pytest verdicts per rank, by token, never by line adjacency** (`log-pathologist`, 2026-10-26, item 644). In a 2-rank stream a trailing `PASSED` can be the other rank's verdict.
- **No window may need more than ≈ 600 s container-side** (`timeout -k 30 600` at most). The Bash foreground cap is 660 s in every session, slot or spawned, and no harness call is backgrounded; the heavy tier's 1200 s is unreachable in one call. A brief that sizes a window above it is unrunnable (item 635) — split it or change its width, never background it.
- **A reproduction control names its numeric bound in the brief, or the print whose last digit sets it** (18:00 ruling on I41, `IMP-3.md`). An executor that sets a tighter bound and relaxes it after a red records both values and the measurement.
- **Daytime rule (standing for clinical-scale builds):** an item marked *daytime* runs only in the 12:00, 13:30 and 15:00 slots, one per slot, as the slot's first item, and only with `uptime` loadavg (1 min) < 8. "Slot-first" means first, ahead of every list position. **No daytime item is open** (981 landed 12:00 11-14). C52's build read gmsh 408 s of the 450 s rule (C30's re-run read 449.40 s; the weekly's flag).
- **(S) on ghosted meshes is ruled** (10-19 weekly). No item asserts a cross-batch byte equality on any mesh. **`.eval` bytes depend on batch composition** (item 537), so every bitwise comparison is call by call on the same batch.
- **Item 526's W5b log, `EX-142`'s first log, `xl` entry 51's log (since 2026-10-28) and item 851's W3 log (since 2026-11-05, `20261106T031528Z_IMP-8-item851-W3.log`) and `xxl` entry 54's log (since 2026-11-07, `20261107T080008Z_GEO-38-C5d.log`, one MUMPS timing line) stay in the gitignored repo-root `logs/leak-held/`** until the operator rules. Any window that hits the leak guard moves its log there unchanged (never to `docs/private/`) and reports. Entry 51's match was the window's own zero-padded node capacitance: **a review writing a node grid into an `xl` command checks each node's printed `.12e` form against the guard first** (entry 49's note).
- **`xl` and `xxl`:**
  - `xxl`: 50, 54, 62 RUN and READ (62 by the 11-14 weekly, green branch). **69 (P33-d2) QUEUED for 11-28** by the 11-14 weekly (`xxl-queue.d/20261128-GEO-38-P33-d2.env`; the launcher takes it on 11-21 if 66 is not READY by the 11-20 18:00 review). **Floor met: 1 ahead.** 66 (11-21) waits on entry 68's read; 70 on 68 / 66.
  - `xl`: 45–63 RUN and READ (61 first branch, 63 red with its branch ruled; the 11-14 weekly).
    - **64 (C5f) RUN 2026-11-15** (Status 0, `[capture] rc=0`, 8 158 s, `memory.peak` 377.1 GiB at `-n 16`, orphans 0 / 0; `20261115T080008Z_GEO-38-C5f.log`); ledger row filled by the 10:30 review from a `log-pathologist` reading (`xl-ledger.md:84`); the 11-16 weekly reads the decision rule.
    - **65 (C35w) RUN 2026-11-16** (Status 0, `[capture] rc=0`, 1 244 s, orphans 0 / 0; `20261116T080008Z_GEO-38-C35w.log`); ledger row filled by the 03:00 review from a `log-pathologist` reading (`xl-ledger.md:85`), so **C51 is queued (item 1003)**; the branch (mesh's / fixture's) is the 11-16 weekly's.
    - **PENDING:** 71 (I32″, Tuesday 11-17 at the earliest; on I32c = item 986, whose landing commit marks it READY), 67 (I32-3T, Thursday 11-19 at the earliest; on I32c and entry 71 green), 68 (C5g; on entry 64's read, the 11-16 weekly — 11-18), 70′ (Ti rod; the daily's variant licence after entry 71 green).
    - **Floor NOT met: 0 ahead of 4** (`xl-queue.d` empty). 6 of 6 windows used in the trailing 7 days (11-10 … 11-16); 11-17 has one window free (cutoff 11-11: five charged). 71 and 68 fill 11-17 / 11-18 only if I32c lands and the 11-16 weekly reads 64; 11-19 is 67's, 11-20 70′'s, each on 71 green. No variant written (journal §6b).
  - No slot writes an `xl` command.
- **Two heavy windows never run at once.** Every heavy item's pre-check is `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3`, which must read 0.







**Standing execution rules — moved 2026-10-24 (operator ruling).** The compose-allow constraint, the docker-socket runner-trap diagnosis, the per-item rules (a)–(j) and the durable-capture idiom live verbatim in `docs/automation/implementer-run.md` § "Standing execution rules (from §9)" — **every slot and every review applies them**; they keep their letters, so "rule (i)" still cites them.

*(Items 19–61 — collapsed one-line DONE records, 2026-09-21 … 09-28 — are
archived verbatim in `docs/planning/plan-archive.md` § "§9 DONE items
19–61" (moved by the 2026-09-29 03:00 review, §9 size).)*


*(Items 106–255 — the archive pointers and one-line records of the 2026-09-30 … 10-05 reviews — are moved verbatim to `docs/planning/plan-archive.md` § "§9 archive pointers and one-line records, items 106–255 (2026-09-30 … 10-05)" (2026-10-06 03:00 review, §9 size). Still live from them: item 133 (H8b) is banked; item 203 waits on the weekly's re-pose.)*

*(The 2026-10-17 18:00 review's digest and restock, and the one-line records of items 511 and 512, are archived verbatim in `docs/planning/plan-archive.md` § "§9 digest and restock of the 2026-10-17 18:00 review, …" (moved by the 2026-10-18 03:00 review).)*







**Every item:** the checklist (bare harness calls, `timeout: 660000`, `-s`, BLOCKED in the same commit, one heavy window at a time, the latest Review ruling binds, no `docs/private/`) lives verbatim in the same `implementer-run.md` section, moved 2026-10-24.

*(The pointer paragraphs to the queue blocks of 2026-10-25 18:00 … 2026-11-07 18:00, each with its one-line DONE records, are moved verbatim to `docs/planning/plan-archive.md` § "§9 archive-pointer paragraphs of 2026-11-04 … 11-07 18:00" (2026-11-08 10:30 review, §9 size). Still live from them: 831 (superseded by 884), 839 and 869, below.)*

*(The archive-pointer paragraphs of the 2026-11-09 03:00 … 2026-11-11 18:00 queue blocks are moved verbatim to `docs/planning/plan-archive.md` § "§9 archive-pointer paragraphs of 2026-11-09 03:00 … 11-11 18:00" (2026-11-12 03:00 review, §9 size). Still live from them: 929 and 957, below; 941 / 933 re-posed as 965 / 968.)*

*(The 2026-11-12 10:30 queue block — the one-line records of items 958 … 957 and its running total — is moved verbatim to `docs/planning/plan-archive.md` § "§9 queue block of the 2026-11-12 10:30 review" (2026-11-12 18:00 review, §9 size). Its open entries 968 / 969 / 970 / 959 are carried below unchanged.)*

*(The 2026-11-12 18:00 queue block — item 975's DONE record, the open entries as written at 18:00, the not-queued list and the running total — is moved verbatim to `docs/planning/plan-archive.md` § "§9 queue block of the 2026-11-12 18:00 review" (2026-11-13 03:00 review, §9 size). Its open entries are carried below, with 968 / 969's dependency marked MET and 970 BLOCKED.)*

*(The 2026-11-13 03:00 queue block — the open entries as written at 03:00, the not-queued list and the running total — is moved verbatim to `docs/planning/plan-archive.md` § "§9 queue block of the 2026-11-13 03:00 review" (2026-11-13 10:30 review, §9 size). Its open entries are carried below; 968 / 969 collapse to DONE lines.)*

*(The 2026-11-13 10:30 queue block — the one-line DONE records of 968 / 969 / 977 / 976, the open entries as written at 10:30, the not-queued list and the running total — is moved verbatim to `docs/planning/plan-archive.md` § "§9 queue block of the 2026-11-13 10:30 review" (2026-11-13 18:00 review, §9 size). Its open entries are carried below.)*

*(The 2026-11-14 03:00 queue block — as edited in place from 11-13 18:00, with the 06:00 / 07:30 slots' markings: 982 / 978 / 979 DONE, 959 BLOCKED, its not-queued list and running total — is moved verbatim to `docs/planning/plan-archive.md` § "§9 queue block of the 2026-11-14 03:00 review" (2026-11-14 10:30 review, §9 size). Its open entries are carried below.)*

*(The 2026-11-14 10:30 queue block — the one-line DONE records of 982 / 978 / 979 / 981 / 983 / 985 / 984, the BLOCKED lines of 980 / 959, the open entries' markings as the 12:00 / 13:30 slots left them, the not-queued list and the running total — is moved verbatim to `docs/planning/plan-archive.md` § "§9 queue block of the 2026-11-14 10:30 review" (2026-11-14 18:00 review, §9 size). Its BLOCKED entry 970 is carried below.)*

*(The 2026-11-14 18:00 queue block — 0 open, 970 BLOCKED, its not-queued list and running total — is moved verbatim to `docs/planning/plan-archive.md` § "§9 queue block of the 2026-11-14 18:00 review" (drafted by the 2026-11-15 03:00 review, which timed out uncommitted; landed from `recovered/20261115T1200Z` by the 10:30 review). 970 is superseded by 987 (I76′).)*

*(The 2026-11-15 10:30 queue block — the DONE lines of 987 / 988 / 989 / 990 / 996, the BLOCKED lines of 991 / 992 / 993, 986 HELD, the U164 / U165 return note, the not-queued list and the running total — is moved verbatim to `docs/planning/plan-archive.md` § "§9 queue block of the 2026-11-15 10:30 review" (2026-11-15 18:00 review, §9 size). 986 is carried below; 991 / 992 / 993 are re-posed or re-run as 999 / 998 / 1000.)*

*(The 2026-11-15 18:00 queue block — the DONE lines of 997 / 998 / 1001 / 999 / 1000 / 1002, 986 HELD, the not-queued list and the running total — is moved verbatim to `docs/planning/plan-archive.md` § "§9 queue block of the 2026-11-15 18:00 review" (2026-11-16 03:00 review, §9 size). 986 is carried below.)*

**On-deck queue, 2026-11-16 03:00 review: C51 on entry 65's green window, refresh batches 18–20, and one provisional chain-P step (P92).** Since 18:00, the 19:30 and 21:00 slots landed all six queued items (997 / 998 / 1001 / 999 / 1000 / 1002); 21:00 (after its items), 22:30 and 00:00 drained. **Entry 65 (C35w) ran green 02:00 11-16; its ledger row carries the mechanical read (`xl-ledger.md:85`).**
- Take the first open item **in list order**. An item whose stated dependency or date has not arrived is skipped at no cost. A daytime item is slot-first in its slots.
- **No item starts a window after 01:45** (the `xl` window is 02:00).
- Every item below is independent of the others.

986. HELD (2026-11-15 10:30 review; module hold, the 11-16 weekly's) — `GEO-36` I32c (chain step (no solve), rule (4); implementer; standard, -n 1 / 2 / 3; ≈ 10.4 measured / 45 predicted min; depends on the weekly lifting the hold). **No slot takes it.**
   The (E) / (S2-tip) gate repair on `test_geo36_rod_f_clinical_heating.py`; I32b's stub literals reproduced by string; landing ⇒ entry 71 READY → docs/planning/chunks/GEO-36.md#§9-item-986
1003. **DONE (04:30 slot 2026-11-16)** — W1 6 passed (`20261116T103321Z_GEO-38-C51-W1.log`, 3 s harness); (P1) 34 / 34 + RR strings, (RR) / (CT) 34 / 34 both rungs by value, (Cus) 912 804 / 16 ranks / P33…P96 / `SAVE_PATH_R2`, (C13r) 0.8883, (CN) 34 / 34 differ; record: asymmetric class r2 / r1 1.1268, median 1.5010, min 1.0597 (other ring, 2), max 2.3375 (same ring, 16), 3 classes ≥ 2×; W-env 12 passed, 57 s. C51 DONE; the branch is the 11-16 weekly's. — `GEO-38` C51 (chain C, rule (4); implementer; standard, -n 1, no solve; ≈ 10.4 measured / 20 predicted min; depends on nothing — entry 65 ran green and its ledger row is written)
   The asymmetric class on two rungs: a new module parses entry 65's and entry 44's tracked class tables (exactly one match each, via `read_log_text`), re-derives them from each npz ((RR) / (CT) / (Cus) custody, (P1) against C13's records), and prints the per-class r2 / r1 ratio; no branch asserted → docs/planning/chunks/GEO-38.md#§9-item-1003
1004. DONE (04:30 slot 2026-11-16; no red, stale 124 → 122) — `EX-179` batch 18, `mat:34` (`-n 2 -t 300`) then `mri:39` (`-n 2 -t 300`) refresh (lane 1, refresh lane; implementer; standard, -n 2; ≈ 3.6 measured / 8 predicted min; **takeable from the 04:30 slot 11-16** — the pre docrefs window must name `materials_34_` and `mri_39_`, else stop at no cost)
   Each example's own gates (closing line, Status 0); `materials_34_` / `mri_39_` stale 2 → 0 by name, post `stale=` = pre − 2 → docs/planning/chunks/EX-179.md#§9-item-1004
1005. OPEN — `EX-179` batch 19, `mat:35` (`-n 2 -t 300`) refresh (lane 1, refresh lane; implementer; standard, -n 2; ≈ 3.6 measured / 6 predicted min; **not before 06:13 CST 11-16** — the pre docrefs window must name `materials_35_`, else stop at no cost)
   The example's own gates ("All imported gates hold.", Status 0); `materials_35_` stale 1 → 0 by name, post `stale=` = pre − 1 → docs/planning/chunks/EX-179.md#§9-item-1005
1006. OPEN — `EX-179` batch 20, `imp:19` (`-n 4 -t 480`) refresh (lane 1, refresh lane; implementer; standard, -n 4; ≈ 3.6 measured / 8 predicted min; **daytime: 12:00 / 13:30 / 15:00 slots 11-16, loadavg < 8, slot-first**, under the 11-05 widened-ceiling ruling — the pre docrefs window must name `implants_19_`, else stop at no cost)
   The example's own gates ("All imported gates hold.", Status 0); `implants_19_` stale 1 → 0 by name, post `stale=` = pre − 1 → docs/planning/chunks/EX-179.md#§9-item-1006
1008. **DONE (04:30 slot 2026-11-16)** — W1 9 passed (`20261116T104406Z_MAT-21-P92-W1.log`, 5 s harness); (K) 18 / 18 labels + sar10g strings equal, (K-91) 0.02521149 / 0.02713084 one match each, (Bq-N) 6.4e-16, (Bq-split) 4.803e-13, (Sym) 0 / 0, (N-eq-m) 1.3e-16, (Lim-∞) 0, (N-h) 1.27e-14, (N-lat) 5.948e-6 > `B1_BAND`, (Env) holds (envelope 2.2515e-4); record ⟨θ⟩_lat,3T(t\*) 0.02713068 K/W, 3 T / 1.5 T 0.9999941, lattice / uniform 1.076123, surface fraction 0.012443. P92 DONE (provisional, for the weekly). — `MAT-21` P92, step 52 (chain P / lane 3 closed form; rule (6) provisional; implementer; smoke, -n 1, complex; depends on nothing) → docs/planning/chunks/MAT-21.md#§9-item-1008
**Not queued:**
- **After their windows or the weekly:** U153, 70′ and entry 67 (entry 71 green, which waits on 986's hold); U13b (U13a′ landed, but its port reading is the weekly's: `ports/_sheet_geometry.py:93–97`).
- **Refresh:** `imp:20` (`-n 4 -t 560`; last close 475 s, over the 11-05 ruling's 400 s, so the weekly prices it); `mri:37` (`-n 2 -t 560`, ages in 06:20:42 11-16; last window 201 s; whether `-t 560` is inside the lane is the 10:30 review's call).
- **The weeklies':** 986's hold; entry 65's branch (with C51) and entry 64's read; entry 68's nodes and `ANS-12` SPEC v2; U164 / U165's anchors; option (b) for U13; whether `housekeeping.py` exempts ledger-cited logs from gzip.
- **The operator's items:** T48, W2, U31, `MAT-3` step 2, `OPS-54`.
- **Lanes:** lane 1 — 1004–1006 (refresh); lane 2 — no target; lane 3 — P92 (as chain P); no `EX-57` filler (census `missing=0`).

**Running total (03:00 11-16): 5 open items, ≈ 25.0 measured / 62 predicted minutes; 986 HELD.**
- **Floors:** **5 items met (5).** **160 measured minutes NOT met: short ≈ 135.0.** Every enumerated, re-posable or after-window step is queued. U167 is written but handed to the weekly (it needs a premise token in U166's landed node, `WF-17.md` § item 1007). The rest wait on windows or tonight's weekly.
- **Expect:** 04:30 takes 1003 → 1004 → 1008 by take-next; 1005 after 06:13 (06:00 by take-next or 07:30); 1006 at 12:00 (daytime). The 06:00 … 09:00 slots likely drain.
- **Class medians** (newest five landed per class; BLOCKED not counted): chain solve **9.9**; chain no-solve **10.4**; lane 1 **3.6**; lane 3 **3.8**; example **5.0**; ops **5.1**; daytime build **17.4**.

*(The per-review journal — slot recap, completion audits, plan-work notes,
§10 assessment — lives in the review commits and
`docs/planning/plan-archive.md`, not here.)*

**If every item is done or blocked: stop and journal.** The `EX-57` drained-slot fallback is retired (2026-10-10 weekly epitaph; `missing=0` since 09-29); its paragraph is archived verbatim in `docs/planning/plan-archive.md` § "§9 drained-queue fallback paragraph". `implementer-run.md`'s matching clause is the operator's to drop; until then a slot that finds `missing=0` journals one line and stops.

Every frequency-domain command needs `source /usr/local/bin/dolfinx-complex-mode`
**and** `FEM_EM_REQUIRE_COMPLEX=1`, with `tests/environment` first in the pytest
path list (*2026-10-26 … 10-27: `tests/environment/test_complex_mode.py` only, under
the environment-first ruling; lifted 2026-10-27 10:30 after `OPS-72` — the full
`tests/environment` (≈ 45 s) may lead again, or run as its own `W-env`*), so an environment regression fails before the formulation tests get
blamed.

---

## 10. Success criteria and long-horizon roadmap

### MVP (end of Phase 2)
- [x] Time-harmonic solver reproduces the analytic lossy plane-wave solution to < 5% *(3.61% in L2; decay constant 0.019%, `TH-6`)*
- [x] Helmholtz coil magnetostatic result matches analytic to < 5% *(0.04%)*
- [x] Phantom σ and εᵣ measurably change the solved field *(σ: interior decay
  constants each match their own closed form and their ratio is 10.3232 vs the
  closed-form 10.3116, `MAT-2`; εᵣ: at εᵣ = 78 the measured β = 27.03 rad/m
  matches the εᵣ-dependent closed form 27.02 to 0.059% where vacuum would give
  2.68, `TH-6`. The loaded-**coil** claim landed 2026-07-31: the FEM ΔR of a
  loop over a conductive half-space matches Dodd–Deeds to 1.58%, `MAT-6` step
  2b — in the eddy-current regime, not yet at saline/Larmor; step 3 extended it
  to the production projected drive at 1.5834%, 2026-08-04.)*

### Target (end of Phase 4)
- [x] Loaded birdcage + phantom simulation runs end to end **at 10 MHz —
  and, since 2026-08-26, at 64 and 128 MHz under the same three gates**
  *(ticked 2026-08-25 with `PORT-9` ✅ — the box's own condition below was
  "this box ticks at 10 MHz when leg (d1′) lands", and it landed: the
  displaced 4×4 breaks gate (iii′) on all three classes by two orders
  while reciprocity holds at 2.259e-14. **The Larmor-frequency half
  landed 2026-08-26 — `PORT-11` ✅ at 64 and 128 MHz** under the same
  three gates, frequency the only knob (10 MHz control to 1.158e-10):
  the loaded birdcage runs end to end **at the Larmor frequencies** as a
  self-consistency identity set. Still nothing here licenses an
  absolute-accuracy, resonance or tuning claim — 18:00 review. Original
  text, kept:)* *(the mesh half is
  done: both fixtures generate and are identity-gated in CI as of 2026-08-03,
  `GEO-9` steps 1 + 2b, and graded conductor sizing is gated as of
  2026-08-16, `GEO-15`. The two-torus excitation lineage closed —
  `PORT-1` ✅ 2026-08-15. What remains is ports on the birdcage itself:
  `PORT-9` 🟡 — as of 2026-08-23 the gapped, loaded birdcage has a solved
  lumped-sheet port column at 10 MHz with two pre-stated gates passed
  (C4 spread 0.0152–0.0159% vs 5%; 50 Ω termination margin 598× vs 10×),
  and as of the same day the **4×4 with reciprocity/passivity/C4 is
  solved and green** (leg (d): 2.495292352e-05 vs 1e-3, σ_max
  0.862659137, class spreads 0.0199 / 0.0180 / 0.0108% vs 5%). What is
  left is leg (d1), the geometric negative control of the C4 gate; this
  box ticks at 10 MHz when it lands, and the Larmor-frequency claim is
  `PORT-11`. **Leg (d1′) landed 2026-08-25 and `PORT-9` is ✅** — on
  `GEO-19` step B's mesh and the (d3) power-wave assembly the records are
  σ_max 0.999992805 and spreads 0.0553 / 0.0353 / 0.0214% against the
  tightened (iii′) 0.5%.)*
- [ ] S-parameters derived from the solved field, not a coupling heuristic
  *(the **route** is done — `PORT-1` ✅ 2026-08-15:
  `run_n_port_sparameter_sweep` reads the solved field end to end,
  `‖S−Sᵀ‖/‖S‖ = 2.5494e-05` vs the 1e-3 gate, heuristic retired behind a
  `DeprecationWarning`; reproduced through `EX-20` and `ANS-3` to
  ≤ 3.67e-06. Corrected 2026-08-16, weekly review — the previous
  parenthetical's "still calls the heuristic" predated step 4. The box
  stays unticked because the only fixture with ports is the two-torus
  pair: it ticks when the same machinery gates on the birdcage,
  `PORT-9`.)*
- [ ] S-matrix satisfies reciprocity and passivity within stated tolerance
  *(demonstrated on that same fixture, and as of `PORT-1` step 3a, 2026-08-03,
  through **`PORT-5`'s own metrics** rather than the test's arithmetic:
  `passivity_max_sigma = 1.000000000000` and unit column power sums to `1e-9`,
  `reciprocity_max_abs_delta = 3.4981e-13`. Left open because the matrix is
  still a two-loop air fixture's; what step 3a removed was the "placeholder
  matrices only" objection. **The sweep-level clause is now discharged** —
  `PORT-5` step 1, 2026-08-16: the report `run_n_port_sparameter_sweep`
  returns on the field route reproduces the gated `‖S‖₂` to 1.97e-07 and
  `‖S−Sᵀ‖/‖S‖` to 9e-11, warning-free, with both negative controls
  executed. What keeps the box unticked is the *fixture*, not the route.)*
- [ ] B1+ field matches literature/measured data qualitatively *(routes through
  the coil+phantom fixture, which `GEO-9` step 1 gated on 2026-08-03; nothing
  has yet computed B1+ on it.)*

### Long-horizon roadmap — owned by the weekly planning review

The weekly review (docs/automation/weekly-review.md) maintains this section
with brutal realism: phase goals, subgoals within phases, and dated
assessments extrapolated from **measured pace**, never from hope. Rules of
engagement: a subgoal that has not moved in a month is rescoped or killed, not
carried; parity claims are per-workflow, never per-product (§1); a phase goal
without a named validation target (closed form, literature value, or AED
comparison) is not a goal; every phase milestone lands an `examples/` case
and, where gated physics supports it, an Ansys benchmark case (§5.4).

First weekly review 2026-08-09: re-derived from measured pace per the
seeding note; owned here thereafter. "Phase 5 (current)" is the organizing
goal, not a claim that §6's phases 2–4 are closed — their open gates
(`PORT-1`, the Larmor-regime validation, coil-driven SAR) *are* the content
of the loaded-birdcage goal, and §6 stays authoritative for gate status.

*(Superseded weekly blocks, 2026-08-09 … 09-06 — five pace ledgers, two agent-value measurements, the 09-02 rulings, the 09-02 and 09-06 Phase-5 exit assessments, examples health 09-06 — moved verbatim to `docs/planning/plan-archive.md` 2026-09-19. The current block of each kind is the 2026-09-13 one below; the trend is in `docs/status/weekly-progress.md` from 2026-09-19.)*

*(Subgoal 4's B₁⁺ closed-form validation target, written 2026-09-06 and killed 2026-09-13 — epitaph below: `docs/planning/plan-archive.md`, moved verbatim 2026-09-26.)*

*(Phase-5 exit assessment 2026-09-09, including the dated scoping rule the 2026-09-13 exit decision applied: `docs/planning/plan-archive.md`, moved verbatim 2026-09-19.)*

**Phase 5 — loaded birdcage RF: exited on F-small 2026-09-13; first
human-scale safety quantity 2026-09-25.** How each subgoal ended (the
subgoal texts and their dated assessments 08-09 … 09-02:
`docs/planning/plan-archive.md`, moved verbatim 2026-09-26):

1. *Port-estimator adjudication* — closed 2026-08-16 (`PORT-1` ✅,
   `PORT-10` ✅: the two systematics are additive).
2. *Honest S-parameters from the package, on the birdcage* — `PORT-9` ✅
   2026-08-25, `PORT-11` ✅ 2026-08-26: identities at 10 / 64 / 128 MHz.
   Absolute `S₁₁` is chain F's, below.
3. *Larmor-regime validation* — the sphere closed (`TH-10` ✅); the
   coil-loading trend was killed 2026-08-23 (epitaph) and later measured
   flat in f on one `xl` window (`TH-11` 5d, a record).
4. *B₁⁺ and SAR maps* — `WF-6` ✅ 2026-09-13 on the two-rung convergence
   statement, `MAT-4` ✅, `ANS-2` AGREE 2026-09-18; at human scale `WF-7` ✅
   2026-09-25 (degree 1, one phantom, untuned coil, symmetry-gated).

*(Phase-5 exit assessments 2026-08-09 / 08-16 / 08-23: `docs/planning/plan-archive.md`, moved verbatim 2026-09-19.)*

*(Production element order as decided 2026-08-23 — degree 1 for coil-fed solves, degree 2 for imposed-field phantom solves; the identity miss it rested on was retired by `TH-19` steps 1–3: `docs/planning/plan-archive.md`, moved verbatim 2026-09-26. The decision in force is in the 2026-09-26 rulings below.)*

**Operator directive 2026-08-25 — two fixtures (in force; the full text and
the fixture-scale finding behind it: `docs/planning/plan-archive.md`, moved
verbatim 2026-09-30).** **F-small** is the gate fixture (`ring_radius`
0.07 m; its records do not move; cheap enough for an ordinary slot).
**F-human** is the production fixture (≈ 0.15 m radius, ≈ 0.30 m long, the
high-pass ring-gap topology, 32 ports); the Phase 5/6 deliverables run on
it. Its cost question was answered by `GEO-25` (504 642 cells, §6).

**Phase 6 — tuning.** Mode spectrum of the birdcage (the `TH-9` eigensolver
machinery on the birdcage mesh), lumped capacitors at the gap/port level,
and a circuit co-simulation loop: S-parameters from the EM solve, tuning
and matching in a circuit layer — the HFSS + Circuit split. Validation
targets: birdcage mode frequencies against the lumped-element
ladder-network closed form; tuned S11/capacitor values against an AED
HFSS + Circuit benchmark case (`ANS-7`, adjudicated 2026-09-30, below).
Near-resonance solves are §2.1's ill-conditioning trap by construction.
*(Assessments 2026-08-09 / 08-16 / 08-23: `docs/planning/plan-archive.md`,
moved verbatim 2026-09-30. The current assessment is the 2026-09-30 block
below.)*

*(Phase 6 assessment 2026-09-06 — the r³ cost-model epitaph, the `PORT-14` `REDUCTION_FLOOR_F_SMALL` ruling, the first measured Phase-6 pace: `docs/planning/plan-archive.md`, moved verbatim 2026-09-26.)*

**Operator directive 2026-08-17 (binding on this phase's scoping):** the
production target for real MRI-safety work at 1.5 T is a **high-pass
birdcage** (capacitors in the end-ring segments, not the rungs) with
**32 ports** — i.e. 16 rungs × 2 end rings, one lumped port per ring
gap. Everything gated so far is the 4-leg fixture with one port box per
leg-pair; a high-pass topology needs ring-gap port surfaces the current
`birdcage_port_domain` layout does not emit, and nothing above
`leg_count = 4` has ever been meshed, identity-gated, or costed. When
this phase is broken down, the breakdown must therefore include, before
any tuning claim: (a) parametric leg count re-gated at 16 (the `GEO-9`
identity family + graded sizing at the larger conductor count, with a
measured cost rung — cell count scales with legs); (b) the high-pass
ring-gap port layout as a mesh chunk (the `GEO-16` pattern, on the end
rings); (c) `PORT-9` step 3's circulant-symmetry gate generalized C4 →
C_N (the 32-port S-matrix is block-circulant in the 16 ring-gap pairs);
(d) the AED HFSS + Circuit benchmark commissioned at the production
rung count, not at 4. The 4-leg fixture stays the cheap validation
vehicle — first gates land there — but Phase 6 does not close on it.

**Operator directive 2026-09-04 (binding on this phase's scoping) — the
conductor model.** Every gated coil is solved inside at σ = 800 S/m, the
ceiling the skin-depth-≥-wire-radius mesh rule allows (§7 `TH-15` entry);
copper is δ ≈ 21 µm at 10 MHz and cannot be a volume. Tuning, matching and
per-unit-accepted-power B₁⁺/SAR are loss-partition quantities, so no
parity claim in this phase stands on the 800 S/m coil. Added as a
prerequisite subgoal, before any tuning claim: (e) **realistic conductors
by both routes** — `TH-15` (interior PEC, the lossless limit) then `TH-14`
(Leontovich surface impedance, copper), each gated on a closed form and on
the F-small birdcage's identities, and checked externally by **`ANS-6`**
(the `ANS-4` fixture with a copper coil, AED Finite-Conductivity and PEC
columns). Item (d)'s HFSS + Circuit benchmark at 16 legs is then
commissioned on the copper coil, not on 800 S/m. No date: the first step is
a mesh-probe whose outcome (hole meshing inside or outside the `GEO-23`
family) sets the pace; the review dates it after `TH-15` step 0 reports.
**Same directive, the feature ladder (§9 item 5):** Tier A — `POST-6`
(multi-port superposition), `PORT-14` (RLC sheets), `PORT-15` (the circuit
layer) — is queued *ahead* of the conductor lineage because it is linear
algebra on gated solves and, with `PORT-13`'s 32×32, it is the tuning
workflow itself; `TH-17` (eigenmodes, B3) is Phase 6's first physics target
and is serial on `TH-15` + `PORT-14`. The Phase 6 breakdown's items
(a)–(d) stand; (e) conductors and the ladder are added, not substituted.
Circuit-layer detail is deliberately left to the reviews when the phase
opens (operator: the area is well-covered by literature; ladder-network
closed forms are the entry point).

**The element-order lever (operator directive 2026-08-18, cross-phase).**
Second-order (degree-2 N1curl) elements are to be *evaluated by
measurement* (`TH-12`) and, if they win on accuracy-per-DOF, adopted as
the production order for the Phase-5/6 solves — every open TH/MAT
accuracy question is a cells-per-δ question, `TH-11` step 5b has measured
that the degree-1 route to the 64 MHz bracket does not fit the box, and
the production 32-port case only compounds that. Honest scope of the
lever: it applies to the time-harmonic E-formulation lineage (Phases 2–6,
including SAR fields and port integrals); it does **not** apply to the
closed Phase-1 magnetostatics (degree-2 A is on record diverging under
the penalty gauge — a formulation property, and Phase 1 is not worth
re-gating). Curved second-order *geometry* is the separable second half
of the same idea (it is the answer class for `GEO-15`'s 3.3% faceting
residual) and awaits its own `GEO` chunk. No dated estimate until
`TH-12` step 1 lands a number. *(2026-08-23: steps 1–3 landed; the
decision is recorded above — degree 1 in production for coil-fed solves,
degree 2 for imposed-field phantom solves, `TH-13` the reopening
condition.)*

**Clinical-scale coil fixtures — F-clinical is the production fixture from
here on (operator directive 2026-10-14, interactive session; enumerated as chain C
by the 2026-10-14 weekly, below).** F-human (0.15 m radius, 0.30 m long) is head-coil
sized; clinical body coils sit just outside a **0.70 m patient bore** (wide-bore
1.5 T / 3 T systems). The goal is **two viable clinical body-coil models, one
at 1.5 T and one at 3 T — the 1.5 T one first.** Nominal geometry (from
published clinical-scale models; the weekly may refine a value only with a
cited source, and records which):

| Fixture | f | Legs | Rung (conductor) diameter | RF shield diameter | Rung length | Source |
|---|---|---|---|---|---|---|
| **F-clinical-1.5T** (first) | 64 MHz | 32, high-pass | ≈ 0.71 m | ≈ 0.745 m | ≈ 0.45 m | 1.5 T body-coil rung length (US patent 11 442 125: 45 cm rungs on a 65 cm carrier); diameter set just outside the 0.70 m bore |
| **F-clinical-3T** (second) | 128 MHz | 32, high-pass | 0.704 m | 0.743 m (shield length 0.70 m) | 0.58 m | the published 3 T 70 cm-bore body-coil model (MRM, "Efficient high-frequency body coil for high-field MRI"); Skyra 70 cm bore modelled as 32-rung high-pass |

Shield within ≈ 2 cm of the end-ring diameter (same patent). Conventions are
the standing ones: copper via Leontovich on the hole route, **exact CAD
rectangle** capacitor / port sheets with κ off (ruling 2026-10-02), the 0.2 Ω
capacitor ESR (ruling 2026-10-06), degree 2 as the default order with degree 1
only ladder-licensed (order rule 2026-10-01). **Load:** the ASTM F2182 gel box
(650 × 420 × 90 mm, `MAT-21`'s fixture) placed as the standard places it —
this fixture is what makes the real ASTM geometry fit, and it **answers W5**
(the ASTM phantom did not fit F-human's bore). **Placement:** the **Phase 6
exit stays on F-human** (its records stay valid and are frozen as an
intermediate fixture; no Phase 6 re-scope); F-clinical is the fixture for all
**new** human-scale work — Phase 7 implant heating, Phase 8 heating, future
`ANS` cases and the Phase 9 examples that claim clinical relevance. **Ladder
(1.5 T first, then the 3 T twin):** (1) `GEO-38` cost probe — CAD identities
asserted, cell counts at two mesh-size rungs, degree-1 / degree-2 unknown
counts and a measured price; (2) mesh gates (port census, tags, exact
sheets); (3) the untuned **64×64** at degree 1 (32 legs ⇒ 64 ring-gap sheets,
operator ruling 2026-10-15) with `PORT-11`'s reciprocity / passivity /
symmetry gates; (4) the in-model Larmor tune with ESR (the `WF-10`
route); (5) degree 2 on `xl` / `xxl`, priced by step 1; (6) the ASTM box
loaded: `|B₁⁺|`, 1 g / 10 g SAR, located hotspot; (7) an `ANS` case for the
operator's AED replication. Not claimed by any step until it lands: tuning,
compliance, or a scanner-specific model.

**Phase 7 — implants.** Parametric implant geometry first (wires, rods,
plates in the phantom; CAD import later), mesh grading around thin
conductors (the `MAG-13` 1/r lesson, made worse by skin depth), local SAR
and near-implant hot spots. Validation targets: published measured implant-
heating data, and AED comparisons — this is where they matter most. No
dated estimate: no measured pace for this work type exists. **Opened
2026-09-28 in parallel with Phase 6** (operator directive; §9 rule (5)
lane 3): it consumes `WF-7`'s SAR and the F-small phantom, not the tuned
coil. **Subgoals (2026-09-30 weekly; numbered steps in chain I below):**
**(7a) canonical conductors in a lossy host under an imposed field, against
closed forms** — the sphere (`IMP-1`), then the same sphere against an
*exact finite-k* l = 1 field so the comparand has no model error, then the
PEC prolate spheroid, whose tip enhancement 1/L_z is the closed-form
stand-in for a rod tip (`IMP-2`); **(7b) a rod in the coil-driven phantom**
— a generator with CAD identities (`GEO-36`), then with / without the rod on
the F-small birdcage at 64 MHz, gated on identities, the tip 1 g SAR ratio
printed and compared with (7a)'s closed-form estimate from the local
background field (`IMP-3`); its number is adjudicated against AED
(`ANS-10`, commissioned by the weekly that reads `IMP-3`); **(7c) the
published target** — the ASTM F2182 reference rod in gelled saline: the
reference document is the operator's to supply (dashboard) and writing the
comparison is the subgoal; serial on H9, undated; **(7d)** `WF-9`, the
ISO 10974 Tier 3 transfer function, after (7b). Gating chunks: `IMP-1`,
`IMP-2`, `GEO-36`, `IMP-3`, `WF-9`. *(2026-10-03 weekly: (7a) extended
by the lead-shaft canonicals `IMP-4` — a PEC cylinder at finite k — and
`IMP-5`, the same with a dielectric coat; (7b) done on identities; (7d)'s
✅ gate is `WF-9` step 1, chain I8; `ANS-10` re-homed to Palace, I9.)*
*(2026-10-14 weekly: (7a) closed by `IMP-6` — the disc edge — with its
step-2 family frozen after I29; (7b)'s tip ratio is a record, not resolved
by the cap ball (entry 33, I17 / I30); all **new** human-scale implant work
moves to F-clinical (chain C, I31) — the F-human rod records stay as
records.)* *(2026-10-21 weekly: **(7e) implant-driven heating in closed
form, then FEM** — `IMP-7`: `IMP-1`'s PEC sphere's dipole SAR heating a
Robin-cooled gel shell around a conducting core, the steady temperature an
exact l = 0 + l = 2 Legendre series (Carslaw & Jaeger ch. IX; it reduces to
`MAT-18` at β = 0 and to `IMP-1`'s Φ), then the FEM on the two-layer ball
against it at `MAT-18`'s gates — the closed-form half of the
implant-heating question that (7c) will adjudicate against the F2182
reference; opened from the 10-20 18:00 daily's proposal with (Rk) printed,
not asserted. Steps I37 / I38 in chain I.)*

**Phase 8 — thermal.** Pennes bioheat with SAR as the source term;
phantom-regime validation first (gel: no perfusion, so the equation reduces
to heat conduction + source, which analytic solutions cover). Mathematically
the easiest phase; the risk is validation data and the EM–thermal
interface, not the solver. No dated estimate, same rule. **Opened
2026-09-28 in parallel with Phase 6** (operator directive; §9 rule (5)
lane 3); the gel regime needs a prescribed SAR source, not a tuned coil.
**Subgoals (2026-09-30 weekly; numbered steps in chain P below). Gel
first, perfusion later (operator).** **(8a) the closed-form lane on the
sphere — done:** conduction, FEM SAR as the source, the Pennes sink steady
and transient, the SAR-driven transient, σ(T) fed back (`MAT-9` … `MAT-14`,
`MAT-5` ✅ 2026-09-29 / 30; imposed uniform field, a small sphere, no
external check). **(8b) coil-driven heating of the gel phantom on F-small**
(`MAT-15`): the coil's SAR as the source on `ANS-2`'s 10 MHz fixture,
gated on the scheme's discrete energy identity, the calorimetric initial
slope `dT/dt|₀ = SAR/c` and the C4 identity of the four drives; and the
thermal solve on the phantom's own geometry against the finite-cylinder
Robin eigen-series. **(8c) the external check** — `ANS-9` (HFSS + Icepak),
on `ANS-2`'s fixture because its EM half is already adjudicated (integral
within 5 %); commissioned by the first weekly after `MAT-15` step 1 lands.
**(8d)** two-way σ(T) coupling on the coil-driven fixture (`ANS-9` step 2)
and **(8e)** an ASTM-style gel phantom at human scale, serial on H9 — both
unenumerated. Gating chunks: `MAT-9` … `MAT-14`, `MAT-5`, `MAT-15`,
`ANS-9`. *(2026-10-03 weekly: (8b) done — `MAT-15` ✅; (8a) extended by
`MAT-17` (duty cycle), `MAT-18` (container wall) and `MAT-19` (the virtual
probe) — what a gel-phantom measurement has; (8c) commissioned, SPEC in
`examples/ansys_benchmarks/ans9_coil_driven_gel_heating_10MHz/`; (8d) is
chain P6.)* *(2026-10-14 weekly: (8c) done — `ANS-9` ✅, every row AGREE
on the degree-2 column; (8e)'s closed-form half is `MAT-21` ✅ and its
fixture question is answered by F-clinical — the ASTM box heated by the
clinical coil is C6 / P33; (8d)'s external two-way comparison (P27) is not
posed until the operator requests it.)*

**Phase 9 — a usable package for MR-safety research (redefined 2026-09-30, operator directive, interactive session).** The phase is mostly about making the repository something the operator can use as an **external package** in their own MR-safety research: write a simulation script against a usable API, and bring geometry from a **separate CAD tool** as a STEP file. **Start: opened 2026-10-03 in parallel with Phase 6** (operator directive 2026-10-03, interactive session — supersedes the 2026-09-30 "after Phase 6 closes"): **(9a) the study API + YAML, (9c) packaging and (9d) documentation open now**; (9b) continues as `GEO-34`; **`POST-9`, AMR and `ARCH-1` still wait for Phase 6's exit** (`ARCH-1`'s own embargo, 10-17, is unchanged). **When a Phase 6 chain step and a Phase 9 step compete for a slot, the Phase 6 step goes first.** The weekly enumerates the first 9a / 9c / 9d steps at its next session (the 2026-10-03 one) under the usual rules — validation envelope per entry point, gated results reproduced through the API within their unmoved bands. Purpose (operator): a usable v0 package by the end of November 2026, ahead of the repository going public before 2027-01-10. Subgoals (the weekly writes them at step grain when the phase opens): **(9a) Study-level public API — Python first, YAML as a thin layer** (operator choice): one scriptable study object over what exists today — geometry (STEP + group map, or the built-in birdcage / phantom / implant generators), materials (a library: gel presets `MAT-1`, Cole–Cole tissues `MAT-3`, copper), excitation (ports, lumped RLC, quadrature drive, k = 1 tuning), solve (frequency / sweep, degree, ranks) and results (S / Z / Touchstone, `|B₁⁺|`, 1 g / 10 g SAR, loss partition, temperature via `thermal/`, exports to ParaView / CSV / JSON); `femem run study.yaml` drives the same API. Every public entry point **declares its validation envelope** — the gates that cover it — so the verification record travels with the API. A public-surface snapshot test and semver guard the surface. **(9b) External-CAD path:** model in any CAD tool → STEP with named bodies / faces → a config mapping names to roles (conductor, gap / port sheet, capacitor, phantom, implant, air box) → a preflight report and figure (tagging, sheet congruence, meshability) → solve; generalises `GEO-34` and closes its step 1c (a STEP from an independent kernel). **Conductor sheets** (operator 2026-09-30): the CAD path accepts a coil conductor drawn as **2D sheets** (the operator's AEDT practice) as well as solids — PEC sheets as interior n × E = 0 facets, copper sheets through the two-sided (decoupled-lip) model of H6a / `GEO-33` step 0 — with the same interior-surface tagging question as the ports below. **Ports as literal 2D sheets** (operator 2026-09-30): the port's electrical definition should be chain F's **F14 voltage-functional element** (one V, one I — HFSS's convention) once it lands, so the CAD sheet's integration line is used, not only drawn; a port, lumped RLC or impedance sheet is a named **2D face** in the STEP plus an **integration line** / direction for its orientation — the HFSS convention — not the built-in fixtures' split half-boxes (`PORT_LOWER/UPPER + i`), which exist only because interior dim-2 physical groups did not reach DolfinX's facet tags on the 0.7 image (known-issues 9/10, 2026-08-05; `_interface_facet_tags` rebuilds each sheet from a cell-tag interface). **First step, a probe on the 0.11 image:** a box with one embedded interior rectangle tagged as a dim-2 group, read at `-n 1 / 2 / 4` — facet count and area against the closed form, no hang. Green ⇒ native facet tagging for 2D ports; red ⇒ the importer inserts a thin sliver pair behind the user's 2D face and reuses `_interface_facet_tags`, so the user still draws only a sheet. The built-in fixtures keep their boxes until switched and re-recorded (as `GEO-35`). **(9c) Packaging:** pip-installable and importable outside the repo (no `PYTHONPATH=src`), a published container image carrying real and complex builds, MPI launch helpers and a resource estimate before a solve. **(9d) Documentation:** user guide, API reference, tutorials drawn from `examples/`. **Clean-up queued for the phase:** `GEO-35` (drop the inert floating leg-port boxes from unbroken-leg birdcages; operator 2026-09-30). **(9e) Included, not leading:** `POST-9` adjoint sensitivities (below), `ARCH-1` (field exports from HFSS / Palace through one interface — **framed by the operator 2026-10-05 as a pluggable EM backend, below**), AMR (ladder C4), MPI scaling. **Exit gate:** the operator writes a new study **from the docs alone**, in a separate directory and environment, importing the *installed* package, from a STEP file made in a separate CAD tool, and reproduces a gated result within its band (e.g. `ANS-2`'s coil-driven SAR). *Previous text (2026-09-28), kept:* Besides AMR (ladder C4), MPI scaling and sweeps, the phase carries **`POST-9`, adjoint sensitivities of safety functionals** (operator directive 2026-09-28): on a reciprocal, complex-symmetric system the adjoint re-uses the forward factorization, so gradients with respect to every tissue parameter, lumped capacitor or port drive cost one back-substitution — the capability that separates this platform from a closed solver for worst-case implant-heating search and first-order uncertainty on safety margins. Validation targets: closed-form derivatives on the `TH-1` slab and the `TH-10` sphere, Taylor tests, a measured cost ratio. Steps 1–3 depend on nothing in Phases 6–8. No dated estimate.

**Epitaphs.**
- *2026-11-14 — entry 63's decision rule, with no branch for a gate the
  window asserted.* The rule named the (S2-tip) and EM classes; the window
  went red on (E) as well, and for three days nothing serial on it could
  open while 21 slots drained. A window's rule names a branch for every
  gate class its module asserts — the brief lists them — or the window is
  not READY.
- *2026-11-14 — I32b's extra step inside a shared maximum.* The Δt₁/10 BE
  step was added for the tip rows and folded into the arm's asserted (E),
  never printed alone; the stub, at ≈ 25× the clinical source density,
  read it at 9e-11 and could not see the clinical floor (≈ 1.6e-08). A
  quantity added to an asserted maximum is a new assertion: it gets its
  own band and its own print, or it stays printed.
- *2026-11-14 — U13 as posed on 10-05: a 1e-9 reproduction of a
  primitive-mesh record through a STEP import.* Three readings (1b's
  "mesh differs", the OCC writer's dropped names, the named-sheet route's
  skipped optimiser) say no import route meshes identically to the
  primitive generator. A route change is gated on identities the new route
  can hold (CAD at `STRIP_BAND`, the port identities), with the old record
  printed as a distance.
- *2026-11-11 — P83's "(Sym) = 1", the field-strength heating lever
  predicted from pictures that share a shape by construction.* Three
  closed-form SAR pictures, each scaled by its side's power share, could
  only return the power share (3.69); the FEM read 4.04 at c\* because the
  3 T source is not the 1.5 T one scaled. A lever's shape part is computed
  from each side's own data (P89), or the prediction is the power share
  restated and is labelled so.
- *2026-11-11 — U152 as posed: a 1e-12 identity on a one-step slope.* The
  route carries absolute temperature, so its own exact balance reads 1e-10
  on that step; the identity could not pass on either route. An identity
  band is set against the route's own arithmetic floor, printed first
  (`(E)` was in the log before the band was chosen), and a slope read on a
  derived increment is gated by its physics band, never by the identity
  band of a finished field.
- *2026-11-07 — "a ≈ 40 % rod-model lever", read on two meshes.* The
  11-04 reading set the kept-cell rod's η_rod(T1) on I41's build against
  the PEC-hole rod's on entry 43's tip-refined ladder and called the
  difference the model's. On one mesh and one definition (I67) the model
  moves it 3 %; the rest was the mesh. A lever is read with everything
  but the lever held fixed, or it is not a lever yet.
- *2026-11-07 — the tip calorimetric gate, red on every field-derived
  source it met.* (S2-tip) compared a 1 g ball mean of the scheme's P1
  temperature with the ball mean of the raw DG0 source; the scheme sees
  the source through its consistent-mass projection, and on a source
  with cell-scale jumps the two ball means differ by 1e-3 to 1e-2 with no
  Δt dependence. A gate written on a derived quantity names the discrete
  operator the quantity is read through; "the same identity the global
  gate passes" is a hypothesis about the operator, not a band.
- *2026-11-02 — "a 3 T C6 body before entry 52 reads would be written
  blind", two days of chain C at zero.* The 1.5 T body's lattice, sense
  and tags were all ruled on 10-31 and its glue ran on stubs; the 3 T
  twin changes a radius, a frequency and a node. Six reviews declined
  to write it and chain C, whose two windows were the only physics read
  this interval, supplied no slot step. A successor body that composes
  already-ruled routes is written on the window's *command*, not on its
  result; what the result decides is the window's rule, not the body's.
- *2026-11-02 — the pole's dose gap, five steps to a record.* A
  1e-4-relative gap between a FEM dose and its series on a 537-cell
  fixture bought three refinement axes (h, Δt, modes) and a duty trend
  to show it was none of them. A gap on a record is bounded before it is
  chased: the first step prints the comparand's own floor (truncation,
  interpolation, round-off) beside the gap, and a gap inside that floor
  is a record on its first reading.
- *2026-10-31 — the R1 family (`PORT-1` gap-voltage padding), banked at
  three attempts.* Three windows measured the same 3.1 % against a 2.5 %
  band under two owner rules; the route it gates has had no consumer
  since the lumped-sheet port landed. A revalidation item on a route
  nothing consumes is a record the day its first window agrees with its
  known-issues entry, not a family.
- *2026-10-31 — "every step names its successor" against a lost weekly.*
  The successors carried the dailies through one lost review and no
  further: the cap of six per chain is counted per weekly, so the weekly
  that did not run froze three chains for 27 hours while their next
  steps were written and declined. A supply rule that assumes the
  supplier shows up needs a clause for when it does not — asked of the
  operator (W9), and tonight every chain is enumerated two intervals deep.
- *2026-10-31 — the first `xl` body that ran unstubbed (entry 52's C6a),
  four slot steps and one BLOCKED to make safe.* The weekly's step named
  the readouts and left the lattice, the sense, the coverage ball and
  the density map to execution; each was wrong once. A window body's
  geometry, sense and tags are named in the step or marked the daily's
  to fill before it is queued, and its glue runs on stubs (items 701 /
  702's pattern) as part of the step, not as a rescue.
- *2026-10-26 — H34, Palace's lumped-port excess, withdrawn.* A
  cross-code identity was posed on a field export Palace does not write
  and on a mesh that was not the fixture under test. A twin step names
  the file the other code writes before it is enumerated.
- *2026-10-26 — P46, the CN Δt-rate against the series, closed.* The
  comparand carried a Δt-independent spatial error of the same size as
  the time error at the finest step, so the slope read the floor. A rate
  is posed against a comparand whose error shrinks with the quantity
  being varied, or against the scheme itself — and the second is already
  `mat:27`'s record.
- *2026-10-26 — "minutes per item will rise with heavier classes",
  falsified.* The solve-bearing class medians 10.5 minutes; no mix of
  ten-minute steps reaches a forty-minute floor. The floor is met by
  depth of enumeration, not by the weight of each step.
- *2026-10-24 — entry 43's decision interval [3.48, 3.56], read as weak.*
  The interval was the spread of the three rungs it was meant to
  adjudicate, so a fourth reading inside it could only confirm that the
  spread is the spread. A decision interval is derived from something the
  new window does not already contain — a closed form, a finer comparand,
  an independent route — or the window is a record and is called one.
- *2026-10-24 — a disclosed band substitution, read at the readout's
  band.* Entry 44's module asserted 5 % where its readout named 0.5 %;
  the swap was disclosed in a docstring before the window and never put
  to a review, and the window ran. A disclosed substitution is a question
  for the next review, not a ruling; until it is ruled, the band the
  window was bought against is the one it is read at.
- *2026-10-21 — H30 as posed on rows 3–6 and CS_s, narrowed.* A read-only
  step was written on a window's record before the record's keys were
  listed; two of its rows and its diagnostic were in no file the window
  writes. A step posed on a record names the keys the record carries.
- *2026-10-21 — eight questions returned to the weekly, twelve idle
  periods.* Each question was slot-sized and each waited 24–48 h for this
  session while the chain it blocked sat at its provisional cap. A
  returned question is answered in the interval it is asked, or its chain
  is given a step that does not depend on it.
- *2026-10-19 — `TH-15` "step 3's default", closed without a flip.* Six
  weeks on 🟡 waiting for a default change that every consumer had already
  made by keyword; a kwarg every fixture names is the convention, and a
  default flip that only re-records is not a step.
- *2026-10-19 — (S) as a cross-batch byte test, retired.* It asserted that
  two evaluation paths agree to the byte without checking that they
  evaluated the same copy of the cell; on ghosted meshes they never do. A
  bitwise claim names the batch it is bitwise in.
- *2026-10-14 — the `IMP-3` cap-ball ladder, stopped at r/8.* Three rungs
  and an `xl` night refined the ball around the tip while the near-T1 edge
  sat outside it; the third rung moved the ratio back. A ladder on a
  pointwise quantity names the cell the point sits in before its first rung
  is bought.
- *2026-10-14 — P27 as a slot step, withdrawn.* Sixteen re-solves with no
  reusable factorisation never fit a 590 s window; a two-way thermal
  coupling is priced as an `xl` probe when its external comparand is
  requested, not posed as a slot step to be returned twice.
- *2026-10-12 — P21's D/C decision rule, killed.* The far-side finding was
  posed as the ratio of probe D to probe C, and C is a near-null (0.007–0.3
  W/kg) that moved ×270 between meshes while D moved 3 %; the ratio read 2.1
  and decided nothing. A decision rule is written on the quantity the
  comparand reports (θ_D), never on a ratio against a probe the fixture does
  not resolve.
- *2026-10-12 — the degree-1 `ans:9` column as AED's comparand, retired.*
  The thermal benchmark's one DISAGREE was 13 % phantom resolution and 8 %
  element order on our side. A benchmark column is recorded at the production
  convention (degree 2, exact rectangles, a resolved phantom) before it is
  adjudicated; the degree-1 column stays as the record of what a coarse column
  costs.
- *2026-10-10 — the `EX-57` drained-slot fallback, retired.* Its 09-26
  condition ("when its last 15 land, a drained slot stops") has held since
  10-07; four slots this interval opened it to read `missing=0` and stop. A
  fallback that can only confirm it has nothing to do is journal noise; the
  implementer protocol's owner drops the clause (flag 7).
- *2026-10-10 — `WF-10` step 1d, the circuit-layer tune of F-human, retired.*
  The circuit layer predicted the copper coil's capacitor to O(1) (H12c) and
  the in-model root found it in two heavy windows (H17). A tuning construction
  that needs the field solve anyway is posed on the field solve.
- *2026-10-07 — `TH-15` step 2, the two-torus PEC-hole unitarity gate, killed.*
  Seven sub-steps and a month located the 2 % `Z` asymmetry in the
  point-sampled path voltage of one fixture; meanwhile every hole-route
  birdcage (`ANS-6`, `ANS-7`, `WF-12`) passed reciprocity at 1e-14 on the
  lumped-sheet route, which samples no path. A gate whose only consumer has
  been superseded is closed, not carried; the known-issues entry stays as the
  two-torus record.
- *2026-10-07 — "H9 serial on H12c green", re-posed.* The circuit layer's
  prediction of a copper coil's tuning capacitor misses the field solve by
  O(1) at Q ≈ 67 while the 50 Ω closure is exact; the field solve with
  capacitor sheets is itself the full-wave answer. A human-scale tuning claim
  is posed on the in-model root (H17), and the circuit layer is a prediction
  to be explained (H12d), not a gate to wait on.
- *2026-10-03 — H8b, the tuned σ = 800 F-human, banked at three attempts.*
  A coil with Q ≈ 11 was tuned because it was the matrix on disk; its 0.79
  share and its identity miss were its own loss and a retired correction.
  A tuning step is posed on the conductor the coil is built from.
- *2026-10-03 — the κ-corrected `xl` / `xxl` rungs (entries 22, 23), closed.*
  A night and an `xxl` window were commissioned to test a correction; three
  ordinary slots on exact rectangles removed the thing being corrected.
  Before a ladder is bought to characterise an artefact, one slot asks
  whether the artefact can be built out.
- *2026-10-03 — "gap-voltage E-field loss in the phantom", dead.* The 8–11 %
  excess phantom loss it explained is 0.3 % at degree 2. A mechanism named
  for a degree-1 residual waits for the order move before it is written
  down.
- *2026-09-30 — the uncorrected global-resolution `xl` ladders, closed.*
  Six windows refined everything except the region the quantity lives in,
  and three of them reproduced a known answer. A ladder is commissioned on
  the axis its gap lives on, or not at all.
- *2026-09-30 — T7 (`PORT-22` step 2 on SPEC v1), retired as a record.* The
  step was written on a three-capacitor series resonance no birdcage is
  built with. A benchmark configuration is checked against how the device
  is built before a step is posed on it.
- *2026-09-30 — the 09-26 "drop a slot" recommendation, withdrawn.* Idle
  slots were read as spare capacity; they were unsupplied capacity
  (operator directive 2026-09-28). A drained queue is a planning failure.
- *2026-09-30 — `WF-1`, `WF-4`, `WF-8` as Phase 5 chunks, folded into
  Phase 9 (9a / 9d).* None has a validation target and none moved in over a
  month; they return as package and documentation work when Phase 9 opens.
- *2026-09-26 — the "successive-change ratio ≳ 1.5" convergence criterion,
  killed.* On rungs of unequal step a first-order method gives 3.0 and no
  convergence at all gives 1.7; the ladder it was applied to never refined
  the port neighbourhood. A convergence rule names the step ratios and the
  region refined, or it is not a rule.
- *2026-09-26 — `TH-17` gate (a), "mode 1 at 64 MHz with `C_tuned`",
  killed as posed.* The capacitor was chosen for a three-capacitor series
  resonance and the eigenproblem carries four; the two numbers answered
  different equations. A consistency gate is written from the identity, not
  from the frequency one hopes to see.
- *2026-09-26 — `EX-57` as the drained-slot sink, not renewed.* Two
  intervals running, 8 sessions in ≈20 drew only figures. When its last 15
  land, a drained slot stops.
- *2026-09-19 — "AGREE by mechanism" as a verdict class, killed.* The 09-13
  64 MHz `ANS-4` verdict was banked on a mechanism instead of a rung; the
  rung ran (`xl` 09-16) and its pre-registered rule fired on the self class.
  A verdict is read off a rung or it is INCONCLUSIVE.
- *2026-09-19 — the 09-13 "2–3 days per numbered step" forecast, killed.*
  Seven of ten chain steps landed in under a day; a closing step written
  after its diagnosis is priced by its own window, not its family's history.
- *2026-09-13 — the closed-form B₁⁺ band gate on F-small, killed.* A 2 % CG1
  floor cannot host a 5 % band with a 1.1 % comparand bar and a 3–8 % miss;
  subgoal 4 exits on the two-rung convergence statement (`WF-6` step 5 ✅).
  Revive only on a degree-2 or curved-geometry rung. Full text: archive.
- *2026-09-13 — `TH-14`'s Fresnel plane-wave step, retired as written.* At
  copper the Fresnel Γ differs from PEC by ≈ 4e-6, which no FEM band
  separates; the lossy-wall cavity Q replaced it and gated at +0.010 %.
- *2026-09-02 — a `MAG-20` third h-rung, killed.* `MAG-20` closed with a
  fitted rate inside its band on two-and-a-half rungs; a third rung would
  sharpen a Phase-1 number that gates nothing on the mission's path this
  quarter. Revive only if a Phase-5/6 magnetostatic control needs it.
- *2026-08-23 — Phase 5 subgoal 3, second target: "the coil-loading trend
  vs frequency crossing out of the eddy-current regime", killed.* Eight
  gated items and two weeks showed the three-point trend was the
  resolution term (`TH-11` step 4, flat in f), and that the 64 MHz h → 0
  bracket that would have licensed a trend claim does not fit a 64 GiB
  box at degree 1 (2.81 M cells OOM, 0.99 M pegged) or degree 2 (61.94
  GiB on the *coarse* rung). The lesson: a trend target needs its finest
  rung priced before it is named.
  **Caveat added 2026-08-24: the box is no longer 64 GiB.** The operator
  raised the container ceiling to **128 GiB**, so the affordability half
  of this epitaph is now unmeasured rather than false — 2.81 M cells at
  degree 1 and the degree-2 coarse rung both sat *inside* 2× their old
  wall. The *physics* half is untouched and is the load-bearing one:
  `TH-11` step 4 showed the three-point trend was the resolution term,
  flat in f, and more RAM does not make a resolution artefact a trend.
  Reviving this target therefore needs a new argument, not just a bigger
  rung — the weekly review owns that call. The Larmor coil question lives on as
  `PORT-11` (identity gates on the loaded birdcage) and subgoal 4, not as
  a ΔR(f) fit.

*(Examples-and-benchmarks ramp checks 2026-08-16 / 08-23 / 09-02 and the 2026-08-16 ratification of the operator-executed weekly: `docs/planning/plan-archive.md`, moved verbatim 2026-09-26. The current examples-health and benchmarks blocks are the 2026-09-26 ones below.)*

*(Superseded weekly blocks, 2026-09-09 and 2026-09-13 — assessment, pace ledger, third agent-value measurement, the Phase-5 exit decision, the three frozen families, the `ANS-4` Larmor verdict of 09-13, the `TH-19` / `POST-6` / `TH-11` 5d rulings, the tuned-birdcage chain as first written, XL/XXL, examples health, benchmarks, hygiene, flagged-for-the-daily — moved verbatim to `docs/planning/plan-archive.md` 2026-09-19. The current block of each kind is the 2026-09-19 one below.)*

*(Superseded weekly blocks, 2026-09-19 — pace ledger, the 09-13 forecast audit, the Phase-5 exit block, the `ANS-3` / `ANS-4` / `WF-7` 0b / `GEO-33` / production-order rulings, chains T / F / H as first written, the Phase 6 date arithmetic, fourth agent-value measurement, XL/XXL, examples health, benchmarks, hygiene, flagged-for-the-daily — moved verbatim to `docs/planning/plan-archive.md` 2026-09-26. The current block of each kind is the 2026-09-26 one below.)*

*(Superseded weekly blocks, 2026-09-26 — pace ledger, the 09-19 forecast audit, the `ANS-6` / chain F2 / production-order / `PORT-21` / `TH-17` (a) / `WF-7` 0d / `ANS-2` / `GEO-33` / B1 / `OPS-61` / `OPS-64` rulings, chains T / F / H as written on 09-26 with the operator's 09-28 … 09-30 additions (H8′ and its ruling, H-Cu, the strip-sheet end state — all restated in the chains below), the Phase 6 date arithmetic, fifth agent-value measurement, XL/XXL, examples health, benchmarks, hygiene, flagged-for-the-daily — moved verbatim to `docs/planning/plan-archive.md` 2026-09-30. The current block of each kind is the 2026-09-30 one below.)*

*(Superseded weekly blocks, 2026-09-30 — pace ledger, the 09-26 forecast audit, the `ANS-7` / `ANS-6` 2c / F2′ / κ-reading / production-order / `TH-5` / `TH-17` (c) / T7 / chain H / `OPS-61` / `POST-8` / `ANS-9` rulings and ratifications, chains T / F / H / P / I as written on 09-30 with the dailies' provisional steps of 10-02 / 10-03 and their slot outcomes, the supply arithmetic — moved verbatim to `docs/planning/plan-archive.md` 2026-10-03. The current block of each kind is the 2026-10-03 one below.)*

**Operator directive 2026-10-03 (interactive session) — the exact-rectangle
re-records lead chains F and T this week; Palace carries routine external
checks.** (a) **Priority.** The top slot-executable work on chains F and T
until the next weekly is rebuilding the F-small `ANS-4` / `ANS-6` / `ANS-7`
fixtures' port and lumped-RLC sheets as exact CAD rectangles (the
2026-10-02 sheets-and-κ ruling, §9, item (iii)), re-recording each under
the (1*) licence with no band moved, and solving each fixture's Palace twin
on the same mesh (`OPS-61`, §5.4 directive 2026-10-02). Then the `S₁₁`
(10 / 64 MHz) and `f₀` gaps are re-read against the private columns. The
hypothesis under test: the ragged sheet geometry is the common cause of
both open external findings (known-issues 2026-09-19 and 2026-09-30). It
is a hypothesis, not a result, and its re-reading follows the
pre-registered rules. (b) **No new AED run is requested for this.** The
tracked SPECs of `ANS-3` / `ANS-4` / `ANS-6` / `ANS-7` already define every
port and RLC sheet as a planar rectangle (checked 2026-10-03), so the
rebuild changes our side only; the existing private columns are the
comparand. (c) **AED is reserved for adjudication.** Routine external
checks go to Palace, which runs headless in the `fem-em-palace` service and
needs no operator time. New AED requests are for physics Palace cannot
model (thermal coupling, `ANS-9`) or for a verdict that Palace cannot settle,
so the operator's queue stays off the critical path.

*(Superseded weekly blocks, 2026-09-30 — the Phase 6 exit arithmetic, sixth agent-value measurement, XL/XXL, examples health, benchmarks, hygiene, flagged-for-the-daily: `docs/planning/plan-archive.md`, moved verbatim 2026-10-03.)*

*(Superseded weekly blocks, 2026-10-03 — pace ledger, the 09-30 forecast audit, the F9 / F16 / T13 / T6 / `ANS-8` / H8b / H12b′ / H13a / F14 / windows / `POST-1` 7b / `WF-9` rulings and ratifications, chains F / T / H / P / I / U as written on 10-03 with the dailies' provisional steps of 10-04 / 10-05 and their slot outcomes, the supply arithmetic — moved verbatim to `docs/planning/plan-archive.md` 2026-10-05. The current block of each kind is the 2026-10-05 one below.)*

*(Superseded weekly blocks, 2026-10-03 — the Phase 6 exit arithmetic, the Phase 9 no-date note, seventh agent-value measurement, XL/XXL, examples health, benchmarks, hygiene, flagged-for-the-daily: `docs/planning/plan-archive.md`, moved verbatim 2026-10-05.)*

*(Superseded weekly blocks, 2026-10-05 — pace ledger, the 10-03 forecast audit, the T6 / `OPS-61` / `POST-7` / T21′–T22′ / H14 / I13 / P8 / U8 rulings and ratifications, chains F / T / H / P / I / U as written on 10-05 with the dailies' provisional steps of 10-06 and their slot outcomes, the supply arithmetic — moved verbatim to `docs/planning/plan-archive.md` 2026-10-07. The current block of each kind is the 2026-10-07 one below.)*

*(Superseded weekly blocks, 2026-10-05 — the Phase 6 exit arithmetic, the Phase 9 first-rate note, eighth agent-value measurement, XL/XXL, examples health, benchmarks, hygiene, flagged-for-the-daily: `docs/planning/plan-archive.md`, moved verbatim 2026-10-07.)*

**Binding ruling — capacitor ESR (operator 2026-10-06; full text §9).** Every
lumped capacitor in a **new** step or fixture carries a frequency-independent
0.2 Ω series resistance, in the circuit layer (`Z = r + 1/(jωC)`), in-model (a
series-RC sheet) and in every new SPEC (as HFSS's parallel R ‖ C at a named
frequency). Lossless records stay valid as records, labelled lossless; r never
re-reads a red. **Adoption decided tonight (below, rulings):** `PORT-24` carries
the identities; H17 / H9′ and `ANS-7` step 11 are posed with r; H12c is not
re-posed with r; `ANS-8` step 2's SPEC carries r, its existing columns do not.

**Operator directive 2026-11-14 (interactive session) — Phase 9 v0 items.**
(a) **Docs hosting: GitHub Pages.** The 9d site is built with the existing
`mkdocs` / `mkdocs-material` / `mkdocstrings` toolchain (`pyproject.toml`
docs extra) and published to GitHub Pages. This settles one of the three
operator items bounding v0. The weekly poses the chain-U steps (a build
that is a test, the Pages config); a slot never pushes, and the
publish itself waits on the manual push to origin. No AED value or private
file ever enters the built site. (b) **CadQuery:** the operator ran
`scripts/maintenance/install_cadquery_image.sh` on 11-14 ≈ 03:18; `GEO-34`
step 1c reads its log (`logs/maintenance/cadquery_image_*.log`) before
claiming the kernel present. (c) **W9 is deferred, not declined** — the
operator will rule on rule (6)'s cap unit later; until then W8 stands as
written (six per chain per weekly held), and W9 stays on the dashboard
(no need to re-ask each review).

*(Superseded weekly blocks, 2026-10-07 — pace ledger, the 10-05 forecast audit, the `TH-5` 2a / H12c / H16 / I15 / P10 / `GEO-34` 1d / ESR-adoption / `TH-15` 2 / `GEO-35` / `WF-8` rulings and ratifications, chains F / T / H / P / I / U as written on 10-07 with the dailies' provisional steps of 10-08 … 10-10 and their slot outcomes, lanes 1–2, the supply arithmetic — moved verbatim to `docs/planning/plan-archive.md` 2026-10-10. The current block of each kind is the 2026-10-10 one below.)*

*(Superseded weekly blocks, 2026-10-07 — the Phase 6 exit arithmetic, the Phase 9 rate note, ninth agent-value measurement, XL/XXL, examples health, benchmarks, hygiene, flagged-for-the-daily: `docs/planning/plan-archive.md`, moved verbatim 2026-10-10.)*

*(Superseded weekly blocks, 2026-10-10 — pace ledger, the 10-07 forecast audit, the T13 / T29 / T25 / H9′ / entry 29 / H13b / U21 / `ANS-9` / `OPS-70` rulings and ratifications, chains F / T / H / P / I / U as written on 10-10 with the dailies' provisional steps of 10-11 / 10-12 and their slot outcomes, lanes 1–2, the supply arithmetic — moved verbatim to `docs/planning/plan-archive.md` 2026-10-12. The current block of each kind is the 2026-10-12 one below.)*

*(Superseded weekly blocks, 2026-10-10 — the Phase 6 exit arithmetic, the Phase 9 rate note, tenth agent-value measurement, XL/XXL, examples health, benchmarks, hygiene, flagged-for-the-daily: `docs/planning/plan-archive.md`, moved verbatim 2026-10-12.)*

*(Superseded weekly blocks, 2026-10-12 — pace ledger, the 10-10 forecast audit, the F27 / F26 / entry 29 / accepted-power / `ANS-9` D / T30 / `MAT-21` (P1) / H18 / chain-U rulings and ratifications, chains F / T / H / P / I / U as written on 10-12 with the dailies' provisional steps of 10-13 / 10-14 and their slot outcomes, lanes 1–2, the supply arithmetic — moved verbatim to `docs/planning/plan-archive.md` 2026-10-14. The current block of each kind is the 2026-10-14 one below.)*

*(Superseded weekly blocks, 2026-10-12 — the Phase 6 exit arithmetic, the Phase 9 rate note, eleventh agent-value measurement, XL/XXL, examples health, benchmarks, hygiene, flagged-for-the-daily: `docs/planning/plan-archive.md`, moved verbatim 2026-10-14.)*

*(Superseded weekly blocks, 2026-10-14 — pace ledger, the 10-12 forecast audit, the `ANS-9` / entry 33 / I28 / T42–T44 / F37 / entries 36–39 / P27 / U25(b) / U42 rulings and ratifications, chains F / T / H / P / I / U / C as written on 10-14 with the dailies' provisional steps of 10-15 … 10-17 and their slot outcomes, lanes 1–2, the supply arithmetic — moved verbatim to `docs/planning/plan-archive.md` 2026-10-17. The current block of each kind is the 2026-10-17 one below.)*

*(Superseded weekly blocks, 2026-10-14 — the Phase 6 exit arithmetic, the Phase 9 rate note, twelfth agent-value measurement, XL/XXL, examples health, benchmarks, hygiene, flagged-for-the-daily: `docs/planning/plan-archive.md`, moved verbatim 2026-10-17.)*

*(Superseded weekly blocks, 2026-10-17 — pace ledger, the 10-14 forecast audit, the T41 / entry 35 / entry 39 / entry 38 → 42 / entry 40 / I30–I31 / H28 / `TH-19` / `ARCH-1` / `OPS-54` / `POST-8` 3b rulings and ratifications, chains T / H / P / I / U / C as written on 10-17 with the dailies' provisional steps of 10-18 / 10-19 and their slot outcomes, lanes 1–2, the supply arithmetic — moved verbatim to `docs/planning/plan-archive.md` 2026-10-19. The current block of each kind is the 2026-10-19 one below.)*

*(Superseded weekly blocks, 2026-10-17 — the Phase 6 exit arithmetic, the Phase 9 rate note, thirteenth agent-value measurement, XL/XXL, examples health, benchmarks, hygiene, flagged-for-the-daily: `docs/planning/plan-archive.md`, moved verbatim 2026-10-19.)*

*(Superseded weekly blocks, 2026-10-19 — pace ledger, the 10-17 forecast audit, the entry 42 / entry 36 / U65–U67 / execution / (S)-on-ghosted-meshes / `th:4` / `TH-15` / P39 / U68 / H21 / `OPS-54` / known-issues / declined rulings and ratifications, chains T / H / P / I / U / C as written on 10-19 with the dailies' provisional steps of 10-20 / 10-21 and their slot outcomes, lanes 1–2, the supply arithmetic — moved verbatim to `docs/planning/plan-archive.md` 2026-10-21. The current block of each kind is the 2026-10-21 one below.)*

*(Superseded weekly blocks, 2026-10-19 — the Phase 6 exit arithmetic, the Phase 9 rate note, fourteenth agent-value measurement, XL/XXL, examples health, benchmarks, hygiene, flagged-for-the-daily: `docs/planning/plan-archive.md`, moved verbatim 2026-10-21.)*

*(Superseded weekly blocks, 2026-10-21 — pace ledger, the 10-19 forecast audit, the entry 37 / entry 40 / provisional-steps / `ports:1` / `POST-8` 3b-iii / H30 / `ring_mode_rows` / entry 46 / C11 gel / `IMP-7` / U75–U76 / P45 / `TH-5` step 3 / known-issues / declined rulings and ratifications — moved verbatim to `docs/planning/plan-archive.md` 2026-10-24. The current block of each kind is the 2026-10-24 one below.)*

*(Superseded weekly blocks, 2026-10-21 — chains T / H / P / I / U / C as written on 10-21 with the dailies' provisional steps of 10-22 … 10-24 and their slot outcomes, lanes 1–3, the supply arithmetic — moved verbatim to `docs/planning/plan-archive.md` 2026-10-24. The current chains are the 2026-10-24 ones below.)*

*(Superseded weekly blocks, 2026-10-21 — the Phase 6 exit arithmetic, the Phase 9 rate note, fifteenth agent-value measurement, XL/XXL, examples health, benchmarks, hygiene, flagged-for-the-daily: `docs/planning/plan-archive.md`, moved verbatim 2026-10-24.)*

*(Superseded weekly blocks, 2026-10-24 — pace ledger, the 10-21 forecast audit, the entry 41 / Phase 6 exit / entry 43 / entry 44 / lumped-sheet finding / U79–U84 / known-issues 2026-09-07 / R1–R2 / `MAT-16` / W8 / known-issues / declined rulings and ratifications — moved verbatim to `docs/planning/plan-archive.md` 2026-10-26. The current block of each kind is the 2026-10-26 one below.)*

*(Superseded weekly blocks, 2026-10-24 — chains T / H / P / I / U / C as written on 10-24 with the dailies' provisional steps of 10-25 / 10-26 and their slot outcomes, lanes 1–3, the supply arithmetic — moved verbatim to `docs/planning/plan-archive.md` 2026-10-26. The current chains are the 2026-10-26 ones below.)*

*(Superseded weekly blocks, 2026-10-24 — the Phase 6 exit arithmetic, the Phase 9 rate note, sixteenth agent-value measurement, XL/XXL, examples health, benchmarks, hygiene, flagged-for-the-daily: `docs/planning/plan-archive.md`, moved verbatim 2026-10-26.)*

*(Superseded weekly blocks, 2026-10-26 — pace ledger, the 10-24 forecast audit, the entry 45 / entry 46 / lumped-sheet finding / class-spread finding / P46 / R1′ / environment-first / provisional-steps / known-issues / declined rulings and ratifications — moved verbatim to `docs/planning/plan-archive.md` 2026-10-31. The current block of each kind is the 2026-10-31 one below.)*

*(Superseded weekly blocks, 2026-10-26 — chains T / H / I / C / P / U as written on 10-26 with the dailies' provisional steps of 10-27 … 10-31 and their slot outcomes, lanes 1–3, the supply arithmetic — moved verbatim to `docs/planning/plan-archive.md` 2026-10-31. The current chains are the 2026-10-31 ones below.)*

*(Superseded weekly blocks, 2026-10-26 — the Phase 6 exit note, the Phase 9 rate note, seventeenth agent-value measurement, XL/XXL, examples health, benchmarks, hygiene, flagged-for-the-daily: `docs/planning/plan-archive.md`, moved verbatim 2026-10-31.)*

*(Superseded weekly blocks, 2026-10-31 — pace ledger, the 10-26 forecast audit, the entry 47 / entry 51 / entry 50 / entry 52 / C17 / provisional-steps / U97′ / R1 / I48c / U102 / `OPS-54` / example-closure / known-issues / declined rulings and ratifications — moved verbatim to `docs/planning/plan-archive.md` 2026-11-02. The current block of each kind is the 2026-11-02 one below.)*

*(Superseded weekly blocks, 2026-10-31 — chains T / H / I / C / P / U as written on 10-31 with the dailies' provisional steps of 11-01 / 11-02 and their slot outcomes, lanes 1–3, the supply arithmetic — moved verbatim to `docs/planning/plan-archive.md` 2026-11-02. The current chains are the 2026-11-02 ones below.)*

*(Superseded weekly blocks, 2026-10-31 — the Phase 6 exit note, the Phase 9 rate note, eighteenth agent-value measurement, XL/XXL, examples health, benchmarks, hygiene, flagged-for-the-daily: `docs/planning/plan-archive.md`, moved verbatim 2026-11-02.)*

*(Superseded weekly blocks, 2026-11-02 — pace ledger, the 10-31 forecast audit, the entry 48 / entry 49 / provisional-steps / `IMP-8` / `EX-175`–`EX-177` / closing-window / `OPS-77`–`OPS-78` / dose-gap / entries 52–56 / `OPS-54` / known-issues / declined rulings and ratifications, and the two epitaphs added that night (now in the standing list above) — moved verbatim to `docs/planning/plan-archive.md` 2026-11-04. The current block of each kind is the 2026-11-04 one below.)*

*(Superseded weekly blocks, 2026-11-02 — chains T / H / I / C / P / U as written on 11-02 with the dailies' provisional steps of 11-03 / 11-04 and their slot outcomes, lanes 1–3, the supply arithmetic — moved verbatim to `docs/planning/plan-archive.md` 2026-11-04. The current chains are the 2026-11-04 ones below.)*

*(Superseded weekly blocks, 2026-11-02 — the Phase 6 exit note, the Phase 7 / 8 / 9 rate notes, nineteenth agent-value measurement, XL/XXL, examples health, benchmarks, hygiene, flagged-for-the-daily: `docs/planning/plan-archive.md`, moved verbatim 2026-11-04.)*

*(Superseded weekly blocks, 2026-11-04 — pace ledger, the 11-02 forecast audit, the entry 52 / entry 53 / provisional-steps / I57 / C26 / P65 / two-tag helper / I61 / `mri:38` / C30 / refresh-lane / entry 59 sense / known-issues / declined rulings and ratifications, and the two epitaphs added that night (now in the standing list above) — moved verbatim to `docs/planning/plan-archive.md` 2026-11-07. The current block of each kind is the 2026-11-07 one below.)*

*(Superseded weekly blocks, 2026-11-04 — chains T / H / I / C / P / U as written on 11-04 with the dailies' provisional steps of 11-05 … 11-07 and their slot outcomes, lanes 1–3, the supply arithmetic — moved verbatim to `docs/planning/plan-archive.md` 2026-11-07. The current chains are the 2026-11-07 ones below.)*

*(Superseded weekly blocks, 2026-11-04 — the Phase 6 exit note, the Phase 7 / 8 / 9 rate notes, twentieth agent-value measurement, XL/XXL, examples health, benchmarks, hygiene, flagged-for-the-daily: `docs/planning/plan-archive.md`, moved verbatim 2026-11-07.)*

*(Superseded weekly blocks, 2026-11-07 — pace ledger, the 11-04 forecast audit, the entry 54 / entry 56 / entry 55 / provisional-steps / re-poses / U124 / C18 / leak-held / refresh-lane / declined / operator rulings / known-issues / branches rulings, and the two epitaphs added that night (now in the standing list above) — moved verbatim to `docs/planning/plan-archive.md` 2026-11-09. The current block of each kind is the 2026-11-09 one below.)*

*(Superseded weekly blocks, 2026-11-07 — chains T / H / I / C / P / U as written on 11-07 with the dailies' provisional steps of 11-08 / 11-09 and their slot outcomes, lanes 1–3, the supply arithmetic — moved verbatim to `docs/planning/plan-archive.md` 2026-11-09. The current chains are the 2026-11-09 ones below.)*

*(Superseded weekly blocks, 2026-11-07 — the Phase 6 exit note, the Phase 7 / 8 / 9 rate notes, twenty-first agent-value measurement, XL/XXL, examples health, benchmarks, hygiene, flagged-for-the-daily: `docs/planning/plan-archive.md`, moved verbatim 2026-11-09.)*

<!-- 2026-11-09 weekly blocks follow -->

*(Superseded weekly blocks, 2026-11-09 — pace ledger, the 11-07 forecast audit, the entry 57 / entry 58 / provisional-steps / re-poses / `OPS-85` / `ANS-11` v2.1 / `duty` / guide-page / entry 64 / declined / operator rulings / known-issues / branches rulings — moved verbatim to `docs/planning/plan-archive.md` 2026-11-11. The current block of each kind is the 2026-11-11 one below.)*

*(Superseded weekly blocks, 2026-11-09 — chains T / H / I / C / P / U as written on 11-09 with the dailies' provisional steps of 11-10 / 11-11 and their slot outcomes, lanes 1–3, the supply arithmetic — moved verbatim to `docs/planning/plan-archive.md` 2026-11-11. The current chains are the 2026-11-11 ones below.)*

*(Superseded weekly blocks, 2026-11-09 — the Phase 6 exit note, the Phase 7 / 8 / 9 rate notes, twenty-second agent-value measurement, XL/XXL, examples health, benchmarks, hygiene, flagged-for-the-daily: `docs/planning/plan-archive.md`, moved verbatim 2026-11-11.)*

<!-- 2026-11-11 weekly blocks follow -->

*(Superseded weekly blocks, 2026-11-11 — pace ledger, the 11-09 forecast audit, the entry 59 / entry 60 / after-window-licence / provisional-steps / U152′ / I87′ / `ANS-13` / declined / operator rulings / known-issues / branches rulings — moved verbatim to `docs/planning/plan-archive.md` 2026-11-14. The current block of each kind is the 2026-11-14 one below.)*

*(Superseded weekly blocks, 2026-11-11 — chains T / H / I / C / P / U as written on 11-11 with the dailies' provisional steps of 11-12 … 11-14 and their slot outcomes, lanes 1–3, the supply arithmetic — moved verbatim to `docs/planning/plan-archive.md` 2026-11-14. The current chains are the 2026-11-14 ones below.)*

*(Superseded weekly blocks, 2026-11-11 — the Phase 6 exit note, the Phase 7 / 8 / 9 rate notes, twenty-third agent-value measurement, XL/XXL, examples health, benchmarks, hygiene, flagged-for-the-daily: `docs/planning/plan-archive.md`, moved verbatim 2026-11-14.)*

<!-- 2026-11-14 weekly blocks follow -->

**Pace ledger, interval 2026-11-11 21:22 → 11-14 21:00 CST (2.99 days,
measured from `28372b5b`; sources: 81 commits, the attempts journal,
`20261115T030109Z_WEEKLY-COUNTS.log`, the nine daily journals' slot tables,
the §7 glyphs).** 24 §9 items DONE (929, 957, 958, 960–962, 964–969,
972, 973, 975–979, 981–985) in 35 `ok` implementer sessions, of which
**21 whole slots found nothing admissible** (11-11 22:30; 11-12 00:00,
09:00, 15:00, 16:30, 21:00, 22:30; 11-13 00:00, 06:00, 07:30, 09:00,
13:30, 16:30, 21:00, 22:30; 11-14 00:00, 04:30, 09:00, 15:00, 16:30,
19:30), one stopped on a BLOCKED item alone (11-13 19:30, item 980) and
**13 working slots did all 24**; 0.7 per session, 1.8 per working slot.
**179.9 measured minutes on the 24 = ≈ 7.5 measured minutes per item**
(11-11: 6.3; 11-09: 8.7; 11-07: 7.1), plus 31.4 on four BLOCKED attempts
(963, 971, 980, 959 — each a brief or premise defect, each re-posed and
landed as 972 / 973 / 982 / 984 within a day); 257 wall-minutes (11-11:
205). **One row reached ✅ — `EX-195` (`imp:29`) — 1 / 2.99 = 0.3/day; no
physics-gate row** (11-11: 1.0 / 0). Steps: I 4 (C40, I87′, I85c, I89);
C 3 (C52, C37, C33c) + 3 windows (entries 61, 62 green; 63 red); P 5 (P88,
P89′, P87, C53, C54); U 9 (U156–U162, U152′, U12′); refresh 3 batches
(`EX-179` 14, 15, 16) = **24 slot steps + 3 windows in 2.99 days = 8.0
slot steps/day** (11-11: 11.1; 11-09: 22.7). **Four of the 24 were the
dailies' provisional steps under rule (6)** (U159′, U160, U161, U162 —
all chain U); the rest were enumerated chain steps, re-poses, lane items
or examples. By phase: Phase 6 none (exited); Phase 7 four; Phase 8 five;
the clinical fixture three steps and three windows; Phase 9 nine;
revalidation none; examples three batches + one new. Known-issues: 1
opened (2026-11-14, the Netgen SEGV), 4 retired (all 2026-11-14).
Progress against cost: `docs/status/weekly-progress.md` § 2026-11-14.

**The 11-11 forecast, audited.** "≈ 76 now + ≈ 31 after windows + ≈ 120
provisional ≈ 230 of 320–1 440 — ≈ 90 short at the measured wall rate,
≈ 1 210 short at the floor": the eleven now-items and the after-window
steps landed (≈ 110 measured minutes against the ≈ 107 forecast — the
shape held), the dailies wrote **four** provisional steps (≈ 29 minutes,
a quarter of the ≈ 120 forecast, all in chain U), and **21 slots
drained**. The measured total (≈ 180, ≈ 211 with the blocked attempts)
fell ≈ 77 short of the measured wall rate (257) and ≈ 1 260 short of the
floor. The binding constraint was again the chains' comparands behind
windows (61–65) — and, new this interval, **a window that went red on a
gate its pre-registered rule did not name**: entry 63's (E) at Δt₁/10
held four chain-I steps and entry 67 behind this session for three days
(epitaphs below). Rule (6) did not bind: I / P / C sat at 0 of 6
throughout; U reached 4 of 6.

**Rulings, 2026-11-14 (no AED figure was opened tonight; `aed_results/`
carries nothing newer than 10-06 and the private adjudications are 10-14
or earlier).**

* **Entry 61 (I36-3T, `xl` 11-12) — read: the first branch; I36-3T lands
  as a record.** Status 0, 848 s, 34.23 GiB at `-n 8` (under the [35, 70]
  GiB bracket — a price reading, as entry 60's was); every asserted gate
  green ((A0) 1 128 430 rel 0, (U1), census 64, (D-rod), (M2) / (Z),
  (R1) 1.887e-14 / 5.699e-15, (P) ≤ 1.810e-12 with NEG ≈ 2.0e-02 on both
  arms, (S0), arm M vs arm R S to ≈ 1e-4;
  `20261112T080008Z_GEO-36-step11-I36-3T.log:12441–12450`). **The rod in
  the clinical box at 128 MHz, degree 1, at the 3 T degree-1 root
  44.0525 pF:** tip 1 g Q 4.748393 / 4.884059 (top / bottom), 4.775105 /
  4.748393 (apices), η_rod(T1) 3.407173, the rod's share of the volume
  power 7.434370e-05; Im Y_in P33 +1.880e-03 / P41 −1.974e-03 S (signs
  opposite to entry 48's — a record, never a red). **Against entry 48
  (C40, item 967): Q per ball ≈ 1.08–1.10×, η_rod(T1) 0.935×, rod share
  1.99×** — the frequency lever on one rod model and one definition, at the
  gel's 1.5 T σ / ε_r (provisional, W3). By the rule the 3 T implant-heating
  window is posed — entry 67, pre-registered 11-09 and **re-dated tonight
  behind entry 71** (below). An internal statement: no published comparand
  (W3), no tissue, no compliance claim. No band moves.
* **Entry 62 (C19, `xxl` 11-14) — read: the green branch; the degree-2
  clinical column exists.** Status 0, `rc1` / `rc2` 0, 10 902 s (call 1
  — C6 at degree 2 with the box-wide grid — 4 301.9 s: build 480.2 s,
  MUMPS 166.6 + 2 691.8 s, grid 626.7 s, readouts 794.9 s; call 2 — I36
  at degree 2 — 6 590.9 s), `memory.peak` 486.57 GiB at `-n 16`, inside
  both brackets and under the 12 000–14 000 s / 500–540 GiB prediction;
  every asserted gate green on both calls ((Deg) U₂ = 2E + 2F exact, (R1)
  ≤ 2.4e-13, (P) ≤ 4.3e-12 with NEG ≈ 2e-2, (P_RC), (A0m), (C6b) 90 of 91
  grid points, coverage (i); `20261114T080009Z_GEO-38-C19.log:12147–12410`,
  `:25185–25195`). **`ANS-11`'s `COMPARISON.md` degree-2 rows are filled
  from the tracked log (C37, item 983, `20261114T181245Z_ANS-11-C37-page-r.log`;
  (J62) 53 fields, 0 mismatches against the host JSONs, custody only) and
  the degree-1 → degree-2 moves are records: the per-watt figures move
  −0.46 / −0.74 / −0.49 / −1.26 % with the node change (193.5996 →
  185.7629 pF) included — under the 10-24 exit's 2 % rule, so the rule
  does not fire:** the clinical degree-1 `|B₁⁺|` / SAR10g figures (entry
  52) stand as records with the move beside them, and the adjudication
  column is the degree-2 one, as SPEC v2 already says. The box-wide 10 g
  maximum g\* (x+157.5 mm, 2.340 × c\*'s) is a record of the C6b reader at
  degree 2; the rod's degree-2 reading (call 2: tip 1 g Q 4.145 / 4.133 /
  4.116 / 4.089, η_rod(T1) 3.443491, share 3.1332e-05) is **C56**'s to read
  against entry 48 under the same rule (below — by the ledger rows the rod's
  Q moves ≈ −5 … −8 %, so the rule is predicted to fire there). Entry 62's
  measured window prices entries 66 and 69 tonight. No band moves.
* **Entry 63 (I32′, `xl` 11-13) — read: RED, and the branch is ruled.**
  Status 1, 995 s, 43.97 GiB at `-n 8`, orphans 0 / 0; one failing node
  with three reds (`20261113T080008Z_GEO-36-step10-I32-prime.log:12532`):
  `gel twin (E) 9.862e-09`, `rod (E) 9.220e-09` against `E_BAND` 1e-09,
  `rod (S2-tip) 1.104e-04` against `S2_BAND` 1e-04. Every EM gate, (B),
  (E′), the global (S2), (S0) and the printed CN / Δt₁ (E) are green, and
  every θ record (tip 1 g θ(t\*) 3.3567e-02 / 3.3684e-02 K rod, ratio
  1.0890 / 1.0898 to the gel twin; the six probes; the rod share
  3.740924e-05) equals entry 55's at every printed digit. **(a) The
  (S2-tip) red takes its pre-registered branch:** the reading is recorded
  — R_Π 1.0920e-04 / 1.1039e-04 at Δt₁/10 against 1.1012e-03 / 1.0839e-03
  at Δt₁, ratios 10.08 / 9.82, the O(Δt) step term, ≈ 25–50× the gel
  twin's (2.7e-06 / 4.5e-06) because the rod's tip balls touch its
  adiabatic wall; the DG0 form gap is gone (⟨Πq⟩/⟨q⟩ − 1 reproduces entry
  55's R_DG0) — and **the tip read drops to printed on the clinical body**;
  the comparand stays the question, `S2_BAND` never moves. A two-point
  intercept (−1.0e-06 / +2.2e-06) cannot show a zero Δt-independent part;
  a third Δt can, and I32c adds it. **(b) The (E) red is not a chosen
  gate:** I32b folded the extra BE step at Δt₁/10 into the arm's asserted
  maximum (module `:354`, `:418–420`), never printed it alone, and the stub
  (≈ 25× the clinical source density) read it at 8.804e-11 and could not
  see it. By the 11-10 18:00 ruling's class — the route carries T from
  20 °C through direct LU and forms θ by subtraction, so a one-step
  identity's floor is ≈ ε·T∞·ρc·V/(Δt·∫q) ≈ 1.6e-08 at Δt₁/10 at
  ≈ 40.7 W/m³ per watt — **the Δt₁/10 residual is a reading, and the
  asserted (E) returns to entry 55's definition (CN and the Δt₁ BE step)
  with the Δt₁/10 residual printed beside its floor estimate (I32c)**. The
  23× / 28× growth from Δt₁ (not 10×) is unexplained and is I32c's
  separating measurement (the same step with T carried from 0). Nothing
  is loosened: `E_BAND` stays 1e-09 on the gates entry 55 passed, and the
  step I32b added by accident becomes what it should have been — a print.
  **(c) The EM record and every θ record stand as records (they are entry
  55's).** The window is re-run as **entry 71** (I32″) at entry 63's
  measured price once I32c lands; **known-issues 2026-11-06 stays open
  until entry 71 reads green**; I76 (item 970) is re-posed as **I76′** on
  entry 63's log tonight (the F2182 prints are per-watt θ(900 s) records
  the two red gates do not touch; the red is disclosed in its table), U153
  stays serial on entry 71 green, 70′ (the Ti rod) on entry 71 green. **(d)
  Entry 67 rides I32c's repaired gates**, re-dated behind entry 71
  (Thursday 11-19 at the earliest). **Declined:** a `src/` θ-carrying
  thermal solve to lower the floor (no consumer needs a 1e-9 identity on a
  Δt₁/10 step; declined on 11-11 as scaffolding, declined again); widening
  either band; re-running entry 63 unchanged; buying entry 67 on I32b's
  gates.
* **P89′'s verdict — read.** "Lattice resolves the excess": the shape
  part SP = 1.1041889 from each side's own 18-point lattice (entries 52 /
  58, C33's JSON strings re-read equal), the predicted lever 3.6935 × SP =
  4.078321 against C39's measured 4.0425 (0.9 % apart), where the pictures'
  3.6935 missed by ≈ 9 %
  (`20261112T182020Z_MAT-21-P89-W1b.log`). **The 3 T source's shape at
  c\* lives at the lattice's 100 mm scale, not below it** — a record; the
  11-11 epitaph's lesson stands and the question it posed is closed. **P90
  (below)** carries the superposition to the degree-2 lattice (entry 62)
  to predict entry 69's order move at c\*.
* **The after-window licence — applied as ruled.** C40 (967) on entry 61
  and C37 (983) on entry 62 landed under it from the 03:00 / 10:30
  mechanical reads; nothing serial on entry 63 opened (red). The licence
  now reads one line further: **a window read RED mechanically opens
  nothing until the weekly rules its branch, and the weekly's ruling may
  re-pose a step onto the red log where the step's inputs are untouched by
  the red gates** (I76′ tonight).
* **The dailies' provisional steps — U159′ (973), U160 (975), U161 (977),
  U162 (980 → 982): all ratified** (tags dropped); each landed as a record
  page section or a generated section with no band set or moved. **The
  dailies' rulings — ratified as written, now standing:** (a) 10:30 11-12
  — a reader that meets a doubled print (a `-n 2` log, a self-copy control
  table) requires exactly one distinct match with the count printed and
  never takes "the first"; custody compares values (`float(a) ==
  float(b)`) and a print-form-only difference is printed, not red; (b)
  03:00 11-12 — the module hold is the window's import closure read
  forward: its own module, every module it imports and its own module's
  importers, never a module that merely imports *another* window's module;
  (c) 10:30 11-13 — a brief's "≤ 1e-12" parenthetical on a one-step
  thermal identity is overruled by the 11-10 rulings (969's (E) moved one
  unit at 1.52e-10, a reading); (d) 10:30 11-14 — a reader of a pending
  window parses the tracked log and the host JSON is custody only (C37's
  route; the premise node `test_c19_absent` retired); (e) 03:00 11-14 — a
  generator constant a landed test pins is never extended by a brief;
  new keys live in their own dict (`TAG403_15T`), and a brief greps
  `gen.<NAME>` across `tests/` before extending one; (f) 10:30 11-13 —
  I85's guard count is a record re-recorded at 20 with the new member
  asserted, not derived from the list it counts. **Re-poses and readings,
  each ratified:** C54 → 984 (serial on 983, attempt 2 of 4); U162 → 982;
  P89 → 972 and U159 → 973 (the reader premise); `ANS-13` SPEC v1.0.1 (the
  3 T cw factor `1.785856659e-01` from `…P33-3T.log:12293` replacing a
  token found in no log, two citations corrected; no value, no band); the
  Netgen workaround in `src/fem_em_solver/io/step_import.py` (the
  named-sheet route skips the optimiser pass) **accepted as disclosed under
  rule (c), with its diagnosis posed as `GEO-34` 1e** (lane 1, below) —
  the importer's default path is unchanged; I90 **declined for cause,
  ratified** (the 3 T probe figure per accepted watt is entry 67's number;
  I86's kernel cancels in the ratio); the 18:00 11-14 "no U163" and "U13
  returned" — both ruled below. **Reverted: nothing.**
* **U163 — the degree-2 clinical price on the estimates page, layout
  ruled (chain U, now).** The two landed slicing premises are re-posed so
  that a third measured section fits: (Sec) slices the 3 T heading →
  **the next `## ` heading or EOF**, and (Order) counts "| `" rows between
  the estimate table's own heading and **the first `## ` heading after
  it** (`tests/package/test_wf17_estimate_measured_3t.py`,
  `…_15t.py` — one-line premise edits in landed nodes, licensed tonight,
  each disclosed with its before / after line in the brief; the records
  they pin are untouched). U163 then adds a generated **"Measured degree-2
  clinical prices beside the estimator"** section after the 3 T one: `xxl`
  entry 62's two calls (C6 d2 + grid 4 301.9 s; I36 d2 6 590.9 s; window
  10 902 s, 486.57 GiB at `-n 16`, from the ledger row by string) beside
  FC2 (U106 / U113's degree-2 clinical estimator row), the estimator /
  measured ratios printed, a tag-403 column, custody on entry 62's log;
  (Sec-d2) the section equals a fresh render; (Page) the S30 node;
  `docrefs` 0 dead. A record, no band set. Options (a) (rows that dodge
  the "| `" count) and (c) (a degree-2 row inside the 1.5 T section, which
  (K15)'s `rows_used == ("FC1",)` forbids) are declined.
* **U13 — the anchor rewritten (chain U, 9b).** The 10-05 target —
  item 192's `P1` at `LEG_D0_REPRODUCTION_BAND` 1e-9 through a STEP import
  — is withdrawn (epitaph below): no import route reproduces a
  primitive-mesh record (1b's own OCC round trip prints "mesh differs",
  `tests/mesh/test_step_import_four_port.py:26–28`; the OCC writer drops
  named sheets, known-issues 2026-10-07; the named-sheet route skips the
  Netgen pass, known-issues 2026-11-14). **Re-posed as (c) then (b):**
  **U13a** — a CadQuery-authored F-small four-port (the 1d authoring
  script's pattern, four named sheets at the generator's exact-rectangle
  positions) imported through U12′'s route: (CAD) the four sheets tagged
  with fragment lineage and each area at `STRIP_BAND` at `-n 1 / 2 / 4`;
  (Q) a mesh-quality reading (min dihedral angle, max aspect) on the
  un-optimised mesh printed beside the primitive mesh's; no solve.
  Validation target: the CAD identities at `STRIP_BAND`. **U13b** (serial
  on U13a) — a `Study` on that imported mesh at 10 MHz gated on `PORT-11`'s
  identities at their imported bands (`RECIPROCITY_BAND`,
  `PASSIVITY_SIGMA_TOLERANCE`, the C4 class-spread band), item 192's `P1`
  printed as a distance-to-record reading, Touchstone written. Validation
  target: the three identities. That is the external-CAD STEP → solve path
  (9b) stated on gates the route can hold.
* **Entry 69 (P33-d2) — READY and QUEUED tonight** (step 3b): P87 MET
  (969), entry 62 read above; the command in entry 59's shape with the
  degree knob at the 1.5 T degree-2 root 185.7629 pF (node forms cleared,
  `20261115T031801Z_LEAK-CHECK-P33-D2.log`), `-n 16`, priced from entry
  62's call 1 — predicted ≈ 4 000–5 000 s, 460–500 GiB. Saturday 11-28;
  **lexical order gives entry 66's `20261121-…` file precedence if both
  are queued, and if 66 is not READY by the 11-20 18:00 review the launcher
  takes 69 on 11-21 instead** (no dependency is lost; 66 then runs 11-28).
  **Entry 66 re-priced** from entry 62's measured window × 0.93 (cells and
  unknowns): ≈ 9 500–11 000 s, 440–480 GiB; its (iii) is MET (entry 61
  read); (i) waits on entry 68 → the 11-18 weekly.
* **Entry 71 (I32″) — pre-registered PENDING on I32c; entry 67 re-dated.**
  Entry 63's command byte for byte except the chunk and raw-log names,
  entry 63's measured price (995 s, 43.97 GiB), Tuesday 11-17 at the
  earliest; its rule names a branch for every asserted gate class (below,
  `xl-pending.md` entry 71). Entry 67 follows it (Thursday 11-19 at the
  earliest, after 68 on 11-18) once 71 reads green mechanically. 70′ (the
  Ti rod) stays the daily's variant licence after 71 green — Friday
  11-20's candidate.
* **Declined:** a `src/` θ-carrying thermal solve (above); widening
  `E_BAND` / `S2_BAND`; re-running entry 63 unchanged; buying 67 on I32b's
  gates; a 3 T degree-2 thermal window before the 3 T degree-2 root (entry
  70 stays PENDING on 68 / 66); I90; folding §9's standing bullets into
  `daily-review.md` (the 10-24 move was the operator's; the daily may
  propose the list on the dashboard, the operator moves it); editing
  `CLAUDE.md`, `docker-compose.yml`, `check_private_leak.py`,
  `housekeeping.py` or the dashboard's numbered items (operator-owned /
  the daily's).
* **Operator rulings in force:** W8 (rule (6)'s cap of six per chain per
  weekly held); **W9 deferred by the operator on 11-14** — it stays on the
  dashboard without being re-asked and did not bind this interval (I / P /
  C at 0 of 6, U at 4); the **2026-11-14 Phase 9 v0 directive** — docs on
  GitHub Pages with the existing `mkdocs` toolchain (U164 / U165 below; a
  slot never pushes, the publish waits on the manual push to origin, no
  AED value or private file enters the site), CadQuery installed 11-14
  (U12′ read its log and landed); the 2026-10-02 sheets ruling; the
  2026-10-06 ESR; the 2026-10-14 F-clinical directive. Rule (6) counts
  from tonight: every chain at 0 of 6.
* **Known-issues:** 1 opened (2026-11-14, the Netgen SEGV — resolves with
  `GEO-34` 1e's diagnosis or U13a's quality reading, the daily's row); 4
  retired and archived; **2026-11-06 ruled tonight** ("Entry 63's half"
  — the branch above; open until entry 71 reads green); 2026-10-27 waits
  on entry 65 / C51 (the 11-16 weekly); 2026-10-17 (`ANS-7` Q) on T48;
  `GEO-21` and the `OPS-43` observation unchanged.
* **Branches:** the re-poses landed and are ratified, so
  `attempt/WF-17-s43-20261110T2111Z` (U152′'s base, 965 landed),
  `attempt/MAT-21-P89-20261112T1210Z` (972), `attempt/WF-17-s49-20261112T1402Z`
  (973) and `attempt/WF-17-s52-20261114T0135Z` (982) **may be deleted by
  the daily once it confirms every log they cite is on `main`**;
  `attempt/MAT-21-step18-20261025T112325Z` (P46's record),
  `attempt/GEO-36-I32b-20261108T103328Z` and
  `attempt/WF-15-s14-20261105T212843Z` kept for their cited logs; no
  `recovered/*`.

**Epitaphs added tonight:** three, in the standing list above (entry 63's
rule without a branch for its own gates; I32b's gate folded into a shared
maximum; U13 as posed).

**Live chains and lanes — open steps numbered, each with its validation
target and serial prerequisites (§9 rules (4), (5), (6)). Landed and
removed since 11-11: C52, C53, C54, C40, C37, P87 / C55, P88, P89′, I87′,
I89, I85c, C33c, U152′, U156–U162, U12′; entries 61, 62 read, 63 read red
(archive, 11-11 block, with every slot outcome). Every chain ends in a
named comparand; slot steps that do not wait on a window are marked (now);
steps serial on a window are takeable under the after-window licence.**

*Chain T — `ANS-7` ✅ (T41); what remains is the external reading of Q.*

* **T48 (operator, dashboard W7).** As written 10-17: the existing
  `ans:7` AED project's adaptive-pass history of `f₀`, `R_in` and Q at
  64 MHz, read privately; public verdict toward / away. No new AED
  design.
* **T46 (heavy, `-n 8`; serial on T48 reading "flat or away").** As
  written 10-14.

*Chain H — Phase 6 exited; no open slot step; W2 is the operator's.*

* **W2 (operator):** `ANS-8` step 2 at the root is READY since H29;
  adjudicated when the column lands. **H12d** — a record from H25; not
  posed. **H6, H11** — as before (archive, 09-26 block).

*Chain I — Phase 7. Entry 63 is read red and its route is repaired by a
slot step (I32c, now); the clinical figure is bought by entry 71 (Tuesday
11-17 at the earliest) and the 3 T twin by entry 67 after it. I76′ and C56
are window-independent readers (now). The published comparand (7c) is
still W3.*

* **I32c (now; chain step, no solve, standard, `-n 1 / 2`; serial on
  nothing — entry 63 has run, the module is free).** The gate repair in
  `tests/validation/test_geo36_rod_f_clinical_heating.py`, by edit of the
  lines the ruling names and addition elsewhere: (i) `out["E"]` is the
  maximum over the CN transient and the Δt₁ BE step, as entry 55 asserted
  it; the Δt₁/10 step's residual is printed as `(E10)` beside its floor
  estimate ε·T∞·ρc·V_sub / (Δt · ∫q) (the pathologist's form), a reading
  with no band; (ii) a third one-step BE at Δt₁/100 joins the tip rows —
  R_Π at Δt₁, Δt₁/10 and Δt₁/100 printed with the two-point and
  three-point Richardson intercepts (readings); on the clinical body
  (S2-tip) is **printed**, not asserted (the pre-registered branch); on the
  stub, where I32b's proof passes it, (S2-tip) stays asserted at `S2_BAND`
  on the Δt₁/10 row so the comparand's form remains gated where it is
  gateable; (iii) **the separating measurement on the stub:** the same
  Δt₁/10 step with T carried from 0 (`t_inf = t_init = 0` through the
  helper's kwargs, or `T_INF` monkeypatched in this module), its (E)
  printed beside the 20 °C one and their ratio (the floor class predicts
  ≳ 1e3; printed, no band). I85's module (which composes `_run_heat` under
  C21's patch set) inherits the change: its stub log is re-run and the
  literal lines carrying the changed prints re-recorded under the (1\*)
  licence, every other line byte-identical. **Validation target:** I32b's
  stub literals reproduced by string except the lines this step changes
  (listed in the brief), `20261108T181211Z_GEO-36-I32b-L-W0b.log:103–130`;
  W0b / W0 / W0n3 green at `-n 2 / 1 / 3`; Wreg green; `_tip_slope_rows`'s
  existing rows unchanged (U127′ imports it). No `src/` edit; no band
  moves. Landing ⇒ entry 71 READY.
* **I76′ (now; chain step, no solve, standard, `-n 1`; §9 item 970
  re-posed; serial on nothing).** I76 on entry 63's log as it is: the
  clinical F2182-shaped readout — entry 63's probe prints (θ(900 s) 1 mm
  beyond each tip per watt deposited in the rod, both arms) beside I75's
  bracket arms, I81's placement ladder and I83's Φ_a, through I79's
  reader with `ENTRY63_LOG` = `20261113T080008Z_GEO-36-step10-I32-prime.log`;
  the table header states the window's status (EM green, (E) at Δt₁/10
  and (S2-tip) above band as readings — I32c) and the brief's "every
  asserted gate green" pre-check is replaced by "the ledger row is filled
  and the weekly's 11-14 branch ruling is on `main`". A record. Validation
  target: the log's printed strings reproduced by custody; I83's Φ_a at its
  imported band. **DONE (12:00 slot, 2026-11-15; §9 item 987):** 24 passed
  in 4 s at `-n 1` (`20261115T180304Z_GEO-36-I76p-W1b.log:248`, after a
  print-label-only fix to W1's equally green `…180232Z_GEO-36-I76p-W1.log`);
  (K-63) failures none `:89`, (K-c63) non-empty, (P63) |diff| 0 on both
  arms, (L17) 0, (N-row) miss 2.499734 K/W. Record (no band, not a matched
  source): clinical F2182 θ(900 s) per rod watt rod 397.8928 / 401.3429,
  gel twin 250.369 / 250.8551 K/W; rod / twin 1.589225 / 1.599899; rod tip
  mean / SCALE·Φ_a(17 mm) 59.80377. §9 item 970 superseded (DONE-by-987).
* **Entry 71 = I32″ (`xl`, Tuesday 11-17 at the earliest; PENDING on
  I32c).** Entry 63's command and price with I32c's module. Its rule
  (every asserted gate class named): every asserted gate green ((B), (E)
  on CN + Δt₁, (E′), the global (S2), (S0), the EM gates) ⇒ **the first
  clinical-scale implant-heating figure lands as a record** (no band on
  the θ ratio; no tissue, no compliance claim, no external check) and
  known-issues 2026-11-06 retires; the tip rows and `(E10)` are readings;
  (E) on CN / Δt₁ or (E′) or (S2) red ⇒ the step returns to the weekly, no
  second window; an EM gate red ⇒ as entry 48's rule; killed or OOM ⇒ the
  last flushed phase line is the lower bound.
* **Entry 67 = I32-3T (`xl`, Thursday 11-19 at the earliest; PENDING on
  I32c and entry 71 read green mechanically).** Entry 71's twin at
  128 MHz on entry 61's mesh, drive and node through I85's module; its
  rule is entry 71's; the 3 T / 1.5 T ratio of the tip θ per accepted watt
  through I88's reader is the reading (a record).
* **C56 (now; chain step, no solve, standard, `-n 1`; shared with chain
  C; serial on nothing).** The rod at degree 2: C45's `parse_i36_log` on
  entry 62's call-2 block beside entry 48's — tip 1 g Q on the four balls,
  η_rod(T1), the rod's power share, Im Y_in at the root — the degree-1 →
  degree-2 move per figure printed and **the 10-24 exit's order rule
  applied to the clinical rod figures** (a move > 2 % on a per-watt figure
  ⇒ the clinical rod figures are quoted at degree 2 only, with the degree-1
  values as records beside them). PREDICTED by the ledger rows: the rule
  fires (Q ≈ −5 … −8 %). Validation target: custody of both logs' strings;
  the verdict line is the record. Landing ⇒ `IMP-3` / `GEO-36`'s clinical
  rod figures carry their order statement. **DONE (12:00 slot, 2026-11-15;
  §9 item 988):** the rule **fires** on all 10 rule rows — Q −7.068 /
  −6.654 / −6.740 / −5.818 %, η_rod(T1) −5.459 %, arm M share −16.24 %,
  SAR_M / P_vol −7.378 / −6.934 / −7.089 / −6.187 % (node change included);
  the gel twin's SAR_R / P_vol moves < 0.4 %; custody 29 / 29 fields, 0
  value mismatches on both host JSONs
  (`20261115T180742Z_GEO-36-C56-W1-n1.log:79–100`; GEO-38.md § item 988
  Outcome).
* **70′ (the Ti-rod clinical variant — entry 71's command with
  `FEM_EM_GEO36_S8_ROD_THERMAL=ti`; `xl`, the daily's priced-family
  licence after entry 71 reads green; Friday 11-20's candidate).** Not
  numbered until then.
* **(7c)** — the operator's document (W3); undated. **I17 / I30** —
  records (10-14 epitaph).

*Chain C — F-clinical, the production fixture for new human-scale work.
It outranks chain U for a slot. C56 (above) and C58 are slot steps (now);
C51 waits on entry 65 under the after-window licence; entries 64 / 65
QUEUED, 69 QUEUED (`xxl`), 66 / 68 / 70 PENDING.*

* **Entries this week:** 64 = C5f (`xl` Sunday 11-15, the 3 T degree-2
  root grid; read by the 11-16 weekly); 65 = C35w (`xl` Monday 11-16, the
  r2 h-reading of the asymmetric class; the 11-16 weekly); 68 = C5g (`xl`
  Wednesday 11-18 at the earliest, on entry 64's read — the 11-16 weekly
  writes its nodes); 66 = C19-3T (`xxl` 11-21, on entry 68's read — the
  11-18 weekly writes its command; re-priced tonight); 69 = P33-d2 (`xxl`
  11-28, **READY and QUEUED tonight**); 70 = P33-3T-d2 (`xxl` 12-05, on
  68 and 66).
* **C58 (now; chain step, no solve, standard, `-n 1`; serial on
  nothing).** The entry-66 reader — C37's 3 T twin written ahead of its
  window: the `ans:12` page's degree-2 column parses entry 66's tracked log
  when the ledger names it, with a tree-independent "(Pend)" premise (the
  item 918 pattern) while it does not; stub-tested on entry 62's log shape
  — (K) the same parser reproduces C37's 53 degree-2 fields from
  `20261114T080009Z_GEO-38-C19.log` by string; `docrefs` 0 dead.
  Validation target: (K) exact; no solve. Landing ⇒ entry 66's reading is
  a slot step on its mechanical green, not a weekly's. **DONE (12:00 slot,
  2026-11-15; §9 item 989):** (K-n) 53 fields, (K) 0 mismatches against
  `ans:11`'s tracked `metrics.json`, (K-src) node 185.7629; the page prints
  "(Pend)" (ledger rows 0, distinct 0); 65 `[PASS]` / 0 `[FAIL]`, (R-cmp)
  5051 B == 5051 B (`20261115T181338Z_ANS-12-C58-page.log`, 2 s); W1 93
  passed per rank (`20261115T181205Z_ANS-12-C58-W1.log`, 27 s); docrefs
  `dead=0 guide=0`.
* **C51 (lane 1 → chain step, standard, `-n 1`, no solve; serial on
  entry 65 green — the after-window licence).** As written 11-09: C13's
  class-anatomy reader on entry 65's r2 matrix beside entry 44's r1; the
  branch is the 11-16 weekly's (known-issues 2026-10-27). **DONE
  2026-11-16 (04:30 slot, §9 item 1003, attempt 1 of 4;
  `tests/validation/test_geo38_c51_class_two_rungs.py`):** W1 6 passed
  (`20261116T103321Z_GEO-38-C51-W1.log`); custody green on both rungs;
  the asymmetric class (other ring, 1 step) reads r2 / r1 = 1.1268
  (1.0436 % / 0.9262 %), median class 1.5010, 3 / 34 classes ≥ 2×; the
  C13 statistic 0.8883 (r1) → 0.9360 (r2). A record; no verdict.
* **`ANS-12` SPEC v2 (the node)** — the 11-18 weekly's, after entry 68.
* **C56** — see chain I.

*Chain P — Phase 8, gel first. P90 and P91 are slot steps (now, closed
form, numpy); entry 69 is queued; entry 70 waits on the 3 T root.*

* **P90 (now; lane 3 → chain step, smoke, numpy; serial on nothing).**
  `MAT-21` step 50 — **the lattice superposition's degree-1 → degree-2
  move at 1.5 T:** P89′'s machinery on entry 62's degree-2 lattice (the 18
  SAR10g values per accepted watt from the tracked log, through C37's
  parser — (K) equal to C37's cells by string) beside entry 52's degree-1
  lattice; θ(c\*, t\*) per deposited watt from each, the move printed —
  **the closed-form expectation for entry 69's row-2 move** (its rule's
  10 % clause); (Lim) / (Bq) at the imported bands; control: the two
  lattices differ by more than `INTERIOR_BAND` at c\*. A record. **DONE
  (12:00 slot, 2026-11-15; §9 item 990):** (K) 18 / 18 / 2 strings equal
  C37's block; (Lim) 0.000e+00, (Bq)(i) / (ii) ≤ 6.058e-14 on θ₂ and the
  move; control on the raw strings −2.3493e-02 (normalised w(c\*) move
  +3.6659e-03, printed); θ(c\*, t\*) per deposited watt 0.9208001 →
  0.9241751 K/W, move +3.665e-03 — "the 10 % clause does not fire at c\*"
  (`20261115T181714Z_MAT-21-P90-W1.log`, 8 passed, 6 s).
* **P91 (now; lane 3 → chain step, smoke, numpy; serial on nothing).**
  `MAT-21` step 51 — **`ANS-13` row 5 in closed form:** the volume-mean
  rise ⟨θ⟩(t\*) per deposited watt of the Robin box under (a) the uniform
  source (step 1's triple series integrated over the box) and (b) P89′'s
  lattice superposition at 1.5 T; both printed with their ratio, (Lim) at
  `B1_BAND` against the series' own antiderivative; control (N-h): h → 0
  gives P_dep·t\*/(ρc·V_sub) to 1e-12. A record — the number P87's print
  (entry 69's row 5) is read against.
* **P92 — DONE (provisional, for the weekly; 04:30 slot 2026-11-16, `20261116T104406Z_MAT-21-P92-W1.log`, 9 passed, 5 s): ⟨θ⟩_lat,3T(t\*) 0.02713068 K/W, 3 T / 1.5 T 0.9999941.** (lane 3 → chain step, smoke, numpy; serial on nothing) *(daily, 2026-11-16 — provisional)*.
  `MAT-21` step 52 — **`ANS-13` row 5's 3 T half in closed form:** P91′'s
  box-mean series on entry 58's degree-1 3 T lattice ((K) entry 60's strings
  equal), beside the uniform re-print, ratios printed; identities at
  `B1_BAND`; control: the 1.5 T lattice separates by > `B1_BAND` (predicted
  5.95e-6, ceiling 2.25e-4). A record: ⟨θ⟩'s shape part ≈ 1 on the lattice.
  §9 item 1008.
* **P27** — not posed (10-14 epitaph). Perfusion — later (operator).

*Chain U — Phase 9, the usable package. Chain C outranks chain U for a
slot. U163, U13a, U164, U165 are independent (now); U13b follows U13a;
U166 waits on entry 64 and U153 on entry 71 under the after-window
licence.*

* **U163 (now; 9d; heavy by measurement, complex, `-n 2`, no solve;
  serial on nothing).** As ruled above — the two premise edits and the
  degree-2 measured section. Brief: the daily's, from the ruling.
  **DONE 2026-11-15 (19:30 slot, §9 item 998, attempt 2 of 4; `WF-17`
  step 53):** section landed through Wgen, ratios 0.9709 / 0.2653 /
  0.6723 / 0.4388 (`20261116T014216Z_WF-17-s53b-W0.log`,
  `20261116T014250Z_WF-17-s53b-W1-s30.log`).
* **U13a (now; 9b; standard, `-n 1 / 2 / 4`, no solve; serial on
  nothing).** As ruled above — the CadQuery F-small four-port's CAD
  identities at `STRIP_BAND` and the mesh-quality reading. Price from 1d's
  W1 (11 s) and 1b's build (≈ 2–3 min at `-n 2`).
* **U13b (9b; standard or heavy by measurement, complex, `-n 2`; serial
  on U13a).** As ruled above — the `Study` on the imported F-small gated on
  `PORT-11`'s identities, `P1` printed as a distance.
* **U164 (now; 9d; ops, `-n 1`, no solve; serial on nothing).** **The
  docs build is a test:** `python3 -m mkdocs build --strict` into a
  temporary site directory inside the container (the toolchain is in the
  image, `docker/Dockerfile:22–24`), asserting exit 0, zero `WARNING`
  lines, and that the built site carries no `aed_results/`,
  `COMPARISON_private` or `docs/private` path and passes
  `check_private_leak.py`; the build time printed (a price). Validation
  target: the three assertions. If `mkdocs` is not importable in the
  container, the finding is the step's result (operator item: the docs
  extra in the image), not a red.
* **U165 (now; 9d; ops, no solve; serial on nothing).** **The Pages
  configuration, unpushed:** `.github/workflows/docs.yml` building the
  site with the `docs` extra on push to `main` and deploying with the
  official Pages actions, `site_url` set in `mkdocs.yml`; a test that the
  workflow YAML parses, triggers on `main` only, installs only the `docs`
  extra and never references a private path. A slot never pushes; the
  publish waits on the operator's manual push (origin is stale).
  Validation target: the four assertions.
* **U166 (9d; standard, `-n 1`, no solve; serial on entry 64 green — the
  after-window licence).** U161's pattern: the estimates page's measured
  3 T section gains entry 64's row (the degree-2 root grid at `-n 16`,
  from the ledger row by string) beside the 3 T degree-2 estimator row
  (U118's branch at degree 2 if `Study.estimate` carries it, else a new
  pinned row — the brief says which), ratios printed, custody on the log.
  A record, no band.
  **DONE 2026-11-15 (19:30 slot, §9 item 1001, attempt 1 of 4; `WF-17`
  step 54):** own section between the 3 T and degree-2 sections; FC4 not
  priced (`nan`), FC3 beside it 0.1044 / 0.0166 — the brief's
  provisional "not priced" choice stands for the weekly
  (`20261116T015521Z_WF-17-U166-W0.log`,
  `20261116T015532Z_WF-17-U166-W1-s30.log`).
* **U153 (9d; standard; serial on entry 71 green — the after-window
  licence).** The clinical fixtures page quotes the clinical
  implant-heating figure as a record (entry 71's tip θ per accepted watt,
  I76′'s table) with its disclaimers; `docrefs` 0 dead.
* **U116, U122** (`WF-17`, open since 11-02) — as written; low priority.
  **U31 (`WF-18`)** — the operator's interactive commit (W6).

**Lane 1 and lane 2, targets for the daily.** Lane 1: **the refresh lane**
— batch 17 (`materials_33_`, `materials_32_`, `implants_19_` age in ≈ 19:31
11-15; `materials_34_`, `mri_39_` ≈ 00:20 11-16), one batch per review
while a ≤ 480 s or priced example reads stale, the standing floor (≈ 100
heavy lines) excluded; **`GEO-34` 1e** (now; `mesh-probe`, measurement
only, `-n 1`): the Netgen SEGV on the 1d fixture across a mesh-size ladder
and with / without the embedded face — the illegal-tet count and whether
the crash reproduces, a table and a probe script, no fix; it resolves
known-issues 2026-11-14 or names its scope. C51 / U166 / U153 after their
windows. Lane 2: Tier B is closed; no target — the daily invents none.
Lane 3: P90, P91 (closed forms, numpy); `MAT-3` step 2 waits on the Gabriel
table (operator).

**Supply to the 11-16 weekly (Monday 21:00), with the arithmetic.** 23
slots fire before it (22:30 tonight, twelve Sunday, ten Monday — no 21:00
on the weekly's evening). Demand: at the dailies' floor (160 measured
minutes per four slots) **≈ 920 minutes**; at this interval's measured 7.3
wall-minutes per `ok` session (257 / 35) **≈ 170**; at its 7.5 measured
minutes per item and 8.0 items per day, ≈ 16 items ≈ 120 per two days.
Supply enumerated above at the 18:00 daily's class medians (chain solve
9.9, chain no-solve 10.8, lane 1 3.6, lane 3 3.9, example 5.0, ops 5.1,
daytime build 17.4): now — I32c, I76′, C56, C58, U163 = 5 × 10.8 = 54;
U13a 9.9; U164, U165 = 2 × 5.1 = 10.2; P90, P91 = 2 × 3.9 = 7.8; `GEO-34`
1e 3.6; batch 17 (dated 19:30 11-15) 3.6 ⇒ **12 items, ≈ 89 minutes
now**; after windows under the after-window licence — U166 (64, Sunday
03:00), C51 (65, Monday 03:00), U13b (after U13a), U153 (71, 11-17) =
≈ 10 × 4 ≈ 42; rule (6) at six per chain across I / C / P / U ≈ 24 × 5 ≈
120 if the dailies find targets — every chain again has a
window-independent pattern to copy (a reader on a landed log, a closed
form, a page section as a test). **≈ 89 now + ≈ 42 after windows + ≈ 120
provisional ≈ 250 of 120–920 — meets the measured item and wall rates,
≈ 670 short of the floor; stated.** The shortfall is structural (the 10-26
epitaph: 4–11 minutes a step against a 40-minute floor). This interval's
new lesson is in the epitaphs: a window's red that its rule cannot dispose
of costs more slots than the window did. Idle slots are never evidence for
fewer slots.

**Phase 6 — exited 2026-10-24 (archive, 10-24 rulings).** Nothing moved
this interval; W2's column is the external reading, undated. On
F-clinical the 1.5 T degree-2 column exists (entry 62) and its order
move is ≤ 1.26 % on every per-watt figure — an internal statement outside
the exit.

**Phase 7 — rate and what the date depends on.** 4 slot steps in 2.99
days (1.3/day; 11-11: 1.5) — the 3 T rod record read (C40 on entry 61),
I87′'s knob and I85c's re-record stub-proved, `imp:29`. **The first
clinical implant-heating window ran red on a gate nobody chose (entry
63)**; the route is repaired by I32c (a slot step, now) and the figure is
bought by entry 71 (Tuesday 11-17 at the earliest), the 3 T twin by entry
67 (Thursday 11-19). The rod's degree-2 reading (C56) is a now-step. **(7c)
is still the operator's reference (W3)**; no published comparand, no
tissue, no external implant check; no date.

**Phase 8 — rate and what the date depends on.** 5 slot steps in 2.99
days (1.7/day; 11-11: 2.0) — `ANS-13`'s runnable half landed (C53 / C54:
Done-when (1) met), its degree-2 body (P87) stub-proved, and two closed
forms closed two questions (P88: a dose is not per-watt; P89′: the lattice
resolves the lever's excess). **`ANS-13`'s 1.5 T adjudication column is
queued — `xxl` entry 69, Saturday 11-28 (or 11-21 if entry 66 is not
READY)**; the 3 T column (entry 70) waits on the 3 T degree-2 root (68 →
66 → 70, 12-05 at the earliest). The external date is the operator's (W10
/ W11 / W12 build sequence) and the `xxl` calendar's; no compliance claim.

**Phase 9 — rate and what the date depends on.** 9 U steps in 2.99 days
(3.0/day; 11-11: 4.0) — the implant page's one-step section (U152′), the
determinism reading (U158), the clinical page's rod-in-box and
closed-form sections (U159′, U160), the estimates page's measured 3 T and
1.5 T sections (U157, U161, U162), the clinical heating records quoted
(U156), **and the CadQuery STEP import with a named port sheet (U12′) —
9b's kernel is in**; 192 runnable examples (+2: `imp:29`, `ans:13`). **Two
of the three operator items bounding v0 were settled on 11-14** (CadQuery
installed; docs hosting = GitHub Pages); the third is the rename (W6), and
the publish itself waits on the manual push. The enumerated 9b / 9d steps
(U13a / U13b, U163–U166, U153 — nine with their follow-ons) are ≈ 3 days
of slots at the measured 3.0/day, ≈ 11-18; **v0 is not dated beyond that:
it waits on W6 and the push.** Not claimed: a release, a stable schema, a
hosted site, a solve on an external-kernel STEP (U13b).

**Agent value, twenty-fourth measurement.** (a) Auditor: 1 glyph closure
audited (`EX-195`, first try, PASS); `ANS-11` / `ANS-13` Done-when (1)
read on their closing windows; no demotion — the twentieth interval
without a physics demotion. (b) Pathologist: three window readings
(entries 61, 63, 62), all COUNTABLE — entry 61 claims 1, 2, 4, 5 CONFIRMED
and claim 3's memory bracket REFUTED (a price reading, outside the rule);
entry 63 claims 1–3, 6, 7 CONFIRMED and 4 / 5 CONFIRMED as to class, not
as to a zero intercept or the growth factor (the reading this session
ruled on); entry 62 claims 1–3, 5 CONFIRMED and claim 4's moves
UNCOUNTABLE from the log alone (C37 computed them). (c) Reviews 9 / 9
`ok`, implementer 35 / 35, windows 3 / 3 ran (one red by its own gate).
(d) Navigator: not measured this week. (e) Example-runner: one new example
and three refresh batches, audited PASS, zero re-poses. **Footers: seven
agent files at 62–75 days (2026-08-31 / 09-13) — overdue; the refresh needs
the operator's worktree session (W4).**

**XL and XXL slots, 2026-11-14 (step 3b).** Three windows ran and each
decided something: entry 61 (I36-3T, `xl` 11-12) — the 3 T rod-in-box
record and the frequency lever on the rod figures (C40); entry 63 (I32′,
`xl` 11-13) — **red: the two one-step thermal gates conflict at one Δt at
clinical source density, and the (E) red was a gate nobody chose** — the
route is repaired (I32c) and re-bought (entry 71); entry 62 (C19, `xxl`
11-14) — the degree-2 clinical column exists and the order rule does not
fire for the box-loaded figures (`ANS-11`'s column filled). **Kept:** 64
(C5f, `xl` Sunday 11-15, QUEUED), 65 (C35w, Monday 11-16, QUEUED).
**READY and QUEUED tonight:** **69** (P33-d2, `xxl` Saturday 11-28;
`xxl-queue.d/20261128-GEO-38-P33-d2.env`; priced from entry 62's call 1 ≈
4 000–5 000 s, 460–500 GiB; lexical order lets 66 go first on 11-21 if it
is READY by then, else 69 runs 11-21 and 66 on 11-28). **Pre-registered
PENDING tonight:** **71** (I32″ — entry 63 re-run on I32c's module, `xl`
Tuesday 11-17 at the earliest, entry 63's measured price; READY when I32c
lands, the clerk queues `xl-queue.d/20261117-GEO-36-step10-I32-dprime.env`).
**Re-dated:** 67 (I32-3T, Thursday 11-19 at the earliest, on I32c and 71's
mechanical green); **re-priced:** 66 (≈ 9 500–11 000 s, 440–480 GiB from
entry 62; (iii) MET). 68 (11-18, on the 11-16 weekly's read of 64) and 70
(12-05, on 68 / 66) unchanged; 70′ the daily's variant licence after 71
green (Friday 11-20's candidate). **Not spent:** nothing — every `xl`
night 11-15 … 11-20 has a named candidate (64, 65, 71, 68, 67, 70′), the
last two conditional. **Backlog floor:** `xl` — two QUEUED (64, 65) plus
71 READY-on-I32c and 67 / 68 PENDING: the floor of four is met once I32c
lands (expected Sunday); `xxl` — 69 QUEUED (met), 66 / 70 PENDING. Budget:
six `xl` per trailing seven days holds on every night (11-15 … 11-20 are
six nights with 11-14 the `xxl` night); `xxl` one per seven — 11-14, 11-21,
11-28, 12-05.

**Examples health, 2026-11-14 (step 4).** **192 runnable examples by the
docrefs `runnable` line** (`20261114T193433Z_ANS-13-C54b-docrefs.log:159`;
+2 this interval: `imp:29`, `ans:13`). Ramp (gating ✅ / owed / present):
**Phase 9: 7 / 5 / 35 — met; Phase 7: 6 / 5 / 30 — met; Phase 8: 9 / 5 /
43 — met;** Phase 6: ≥ 22 — met, exited with its full five; Phases 1–5
unchanged; no gate closed this interval, so no ramp count moved. `EX-195`
audited PASS; the last `docrefs` window reads 0 dead. **The staleness
picture:** 121 stale artifact lines (11-11: 122) — batches 14 / 15 / 16
refreshed three lines while the fortnight rule aged in ≈ 3 a day; the
standing stale floor (the 11-07 ruling) is unchanged and not pursued;
batch 17's five lines age in 11-15 / 11-16. The examples the chains name
next: none tonight (I76′, C56, P90, P91 are readers and closed forms whose
`mat:` / `imp:` pages the daily may size as examples once landed). The
XDMF-opens-in-ParaView clause stays the operator's standing question.

**Benchmarks, 2026-11-14 (step 5).** *Adjudicated:* no case gained AED
numbers since 11-11 (`aed_results/` nothing newer than 10-06; the private
adjudications 10-14 or earlier). *Commissioned:* **none** — the only new
gated physics since `ANS-13` (11-11) is entry 62's degree-2 column, which
is `ANS-11`'s adjudication column, not a new case; entry 61's 3 T rod
record waits on (7c) for a comparand. *Columns:* `ANS-11`'s degree-2
column is **filled** (C37; SPEC v2 unchanged); `ANS-13`'s `COMPARISON.md`
exists with its degree-1 columns and blank degree-2 / AED columns (C54);
its 1.5 T degree-2 column is `xxl` 69 (queued), the 3 T one entry 70;
`ANS-12` SPEC v2 (the node) is the 11-18 weekly's after entry 68. W10 /
W11 / W12 lines refreshed on the dashboard; W2 (`ANS-8` step 2) and W7
(T48) are the operator's, unchanged. The implant `ANS` case waits on W3.
*Palace:* none commissioned (no thermal solver; the clinical EM twins are
unpriced).

**Plan hygiene, 2026-11-14.** Rotation first (`5332c22e`): 70 attempts
entries (562 lines), 4 retired known-issues entries (the 2026-11-14 set),
3 §7 rows over 2 KB to their chunk files (`GEO-36` 2 174 B, `MAT-21`
2 267 B, `WF-17` 4 343 B), the 11-11 §10 blocks to the archive (three
parts, 689 lines); PROJECT_PLAN.md 5 510 → 4 823 lines before tonight's
block — **the 4 000-line guide is unmet** by §7's ≈ 440 ID rows plus six
chains; declined again (rows over 2 KB after the rotation: none). **The
§7 ANS table was broken** — a blank line after the `ANS-9` row left
`ANS-11` … `ANS-13` rendering as plain text (the census read them
UNPARSEABLE); removed in the rotation commit. `CLAUDE.md`'s exit
sentence landed (operator, `7c2b7bc2`). §9 is the daily's (27.9 KB
against its ≈ 25 KB target; the fold of its standing bullets into
`daily-review.md` is declined above). No protocol file changes tonight.

**Flagged for the 2026-11-15 03:00 daily review — not weekly-owned.**

1. **Queue, the now-steps first, then dated, then after windows:** I32c
   (first — entry 71 waits on it), I76′ (970 re-posed), C56, C58, U163,
   U13a, U164, U165, P90, P91, `GEO-34` 1e; batch 17 at 19:30 11-15; then
   U166 after entry 64's mechanical read (Sunday 03:00), C51 after entry
   65's (Monday 03:00), U13b after U13a, U153 after entry 71's. **Size
   each restock to two intervals (eight slots).** Rule (6) counts from
   tonight: every chain at 0 of 6, and every chain has a
   window-independent pattern to copy.
2. **Rows to write at full rubric grain and the edits ruled tonight:**
   `GEO-36` (entries 61 / 63 read; I32c, I76′; entries 71, 67 re-dated,
   70′), `GEO-38` (entry 62 read; C56, C58; entry 66 re-priced, 69 queued),
   `MAT-21` (P89′ read; P90, P91), `WF-17` (U163 with the two premise
   edits, U164, U165, U166, U153; the four provisional tags dropped),
   `WF-16` / `GEO-34` (U13a, U13b, 1e; U13's epitaph), `IMP-3` (C56's
   order statement when it lands), `ANS-11` (the degree-2 column read),
   `ANS-13` (SPEC v1.0.1 ratified), known-issues 2026-11-06 (ruled —
   written tonight) and 2026-11-14 (the daily's "Resolves with" row:
   `GEO-34` 1e or U13a).
3. **`xl` clerk:** 69's queue file is written tonight (diff it once
   against the entry); 71 PENDING on I32c — the first review after I32c
   lands checks the new printed forms (`(E10)`, the Δt₁/100 row, the
   intercepts, the T-from-0 line) against the leak guard, marks READY and
   writes `20261117-GEO-36-step10-I32-dprime.env`; 67 READY on I32c + 71's
   mechanical green (11-19); 68 after the 11-16 weekly; 66 re-priced, its
   command the 11-18 weekly's; 70 unchanged; 70′ after 71 green. The `xl`
   floor is met once I32c lands; the `xxl` floor is met.
4. **§2 and the dashboard:** the headline gains entry 62's order
   statement (the clinical degree-1 `|B₁⁺|` / SAR10g figures stand, moves
   ≤ 1.26 % at the degree-2 root), entry 61's 3 T rod record (Q 4.75–4.88,
   η 3.41, share 7.4e-05; the lever on entry 48) and entry 63's finding
   (red on an unchosen gate; the figure is entry 71's); item 4's list
   clears (all ratified); item 4c closes (ruled); item 1 reads 69 QUEUED
   and 71 pending I32c; item 1a the twelve now-items; W10 / W12 refreshed
   tonight (weekly-owned lines); W4's ages re-counted 62–75 days; W9
   deferred, not re-asked.
5. **Audits owed:** none newly (I32c, I76′, C56, C58, P90, P91 are
   records on ✅ rows; U13a / U13b change `WF-16`'s step 3 — audit when
   U13b lands).
6. `docs/testing/known-issues.md`: 2026-11-06 carries tonight's ruling;
   2026-11-14 gets its "Resolves with" (the daily's row); nothing else.

---

## 11. Key technical reference

**Function spaces** — H(curl)/Nédélec for `E` and `A`; H(div) for `B`; L2 for
scalar potential in A-V formulations.

**Boundary conditions** — PEC `n×E = 0`; PMC `n·B = 0`; ABC for radiation; PML via
complex coordinate stretching; waveguide ports for S-parameter extraction.

**Materials** — `ε = ε₀(ε' − jε'')`, `μ = μ₀(μ' − jμ'')`; anisotropic tensors;
frequency-dependent dispersion. Gelled saline at 128 MHz (3T): `σ ≈ 0.6–0.9 S/m`,
`εᵣ ≈ 78–80`; phantom diameter 16–20 cm (head), 30–40 cm (body); ≈1% agarose.

**Linear solvers** — MUMPS direct for small/medium; GMRES + ILU iterative; complex
systems required for time-harmonic (`TH-1`).

**Post-processing** — `SAR = σ|E|²/(2ρ)` [W/kg]; `J = σE`; Poynting `S = ½Re(E×H*)`.

### Resources
- [FEniCSX Tutorial](https://jorgensd.github.io/dolfinx-tutorial/) ·
  [API docs](https://docs.fenicsproject.org/dolfinx/v0.7.0/)
- Similar work: [Elmer FEM](https://www.elmerfem.org/) (EM module),
  [OpenEMS](https://openems.de/) (FDTD), [scikit-rf](https://scikit-rf.readthedocs.io/)
