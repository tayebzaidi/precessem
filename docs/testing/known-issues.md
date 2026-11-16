# Known Issues — failing tests and open defects

**Purpose: tell you whether a failure is yours.**

Every entry below was verified failing *before* the work recorded against it, by
checking out the stated commit and re-running. If you hit one of these, you did not
break it. If you hit something not listed here, you probably did.

That distinction is the whole point of this file. Establishing it after the fact
costs a `git stash` / re-run cycle per failure; it was done five times over on
2026-07-27 and is not worth repeating.

Last audited: 2026-07-27 against `ce92e8c`.

## How to check a failure against this baseline

```bash
# Is this failure mine, or was it already there?
git stash                              # or: git checkout <base> -- <paths>
mpiexec -n 8 python3 -m pytest <the failing test> -q
git stash pop
```

If it fails at the base commit too, add it here rather than fixing it in passing —
unless fixing it is the task.

---

## Failing tests

### 🟡 OPEN 2026-11-06 (`xl` entry 55 / `GEO-36` I32, 02:00 window; filed by the 03:00 daily review): the tip calorimetric-slope gate (S2-tip) misses `S2_BAND` by 6–64× on every field-derived source it has met, while the global (S2) on the same run is green

| | |
| --- | --- |
| **Symptom** | `test_geo36_i32_rod_heating` (`-n 8`, `FEM_EM_GEO36_S8_HEAT=1`) is red on one assertion: `rod (S2-tip) 6.420e-03` against `S2_BAND` 1e-04 (`20261106T080008Z_GEO-36-step10-I32.log:12523`). Both rod-arm tip balls miss: top 2.4589e-03 and bottom 6.4201e-03 (`:12512–12513`). The gel twin, with no rod and no wall, prints 3.1051e-03 / 1.3389e-03 (`:12501–12502`); it is not asserted there (`_thermal_failures` gates the tip on the rod arm only). Every other gate is green: (A0), U₁, census 64, (D-rod), (R1), (P) / NEG on both arms, (B), (E), (E′), and the global (S2) 1.2537e-05 / 1.2541e-05 (`:12509`, `:12520`; `log-pathologist`, 2026-11-06 03:00 review). The same gate failed on `WF-15` step 14 / U127 (§9 item 831): top 4.3368e-03 and bottom 5.9425e-04, nearly unchanged at Δt₁/10, while (S2) scaled 10× (`20261105T212316Z_WF-15-s14-W1-r3.log` on `attempt/WF-15-s14-20261105T212843Z`). |
| **Verified at** | `main` at `e0865297` by execution, in the `xl` window: 1015 s, `memory.peak` 44.12 GiB. Status 1, `[capture] rc=1`, no orphans. |
| **Cause (measured on the stub: form CONFIRMED, `OPS-83` 2026-11-06, `20261106T120232Z_OPS-83-W1.log:78–84`)** | On I32a's 12³ stub with a cube-scale checkerboard q_r = q_u·(1 + ½·s), one BE step: R_DG0 = \|ρc⟨θ¹⟩_B/(Δt⟨q⟩_B) − 1\| is 1.5793e-03 / 2.1076e-03 (tip / interior ball) at Δt₁ and 1.5732e-03 / 2.1146e-03 at Δt₁/10 (Δt-ratio 1.004 / 0.997, Δt-independent), while against the projection, R_Π = \|ρc⟨θ¹⟩_B/(Δt⟨Πq⟩_B) − 1\|, it is 6.8172e-06 / 7.9282e-06 → 6.8552e-07 / 8.5219e-07 (ratio 9.944 / 9.303, O(Δt)); ⟨Πq⟩_B/⟨q⟩_B − 1 is 1.573e-03 / 2.115e-03, the whole of the gap. The lumped Π_L does not close it (R_L 5.15e-05 / 1.54e-05, ratio 1.135 / 1.846). On the uniform source all three coincide (K10-S tip 8.0400e-06 reproduced to 1.2e-11). Global (S2) 2.0247e-05 / 2.0246e-06 on q_r. The pre-stated reading printed "form CONFIRMED". *Open clause as filed (now answered below):* "that the clinical-scale reds (entry 55, U127) are wholly this form gap; the stub shows the mechanism at 1e-3, not its size on a degree-1 \|E\|²." **Measured on a field-derived source (`OPS-84`, 2026-11-06, `20261106T180958Z_OPS-84-W1.log:3747–3762`):** on I41's matched rod pair (degree-1 \|E\|², arm M rod-as-gel / arm R PEC rod, ccw and P1 drives, both 1 g tip balls), R_DG0 at Δt₁ is 1.645e-04 … 9.435e-03 (Δt-ratio 0.972–1.074 on seven of eight rows; 2.245 on the eighth, M ccw bottom, whose R_DG0 is 1.645e-04) and tracks \|⟨Πq⟩_B/⟨q⟩_B − 1\| (6.3e-05 … 9.59e-03) — the same 1e-4 to 1e-2 range as U127 and entry 55. Against the consistent-mass projection R_Π is 6.22e-05 … 1.520e-04 at Δt₁ with Δt-ratio 9.94–10.06 (O(Δt)), and 6.26e-06 … 1.53e-05 at Δt₁/10, every row under `S2_BAND`. The lumped Π_L does not close it (R_L 2.7e-02 … 1.0e-01, Δt-independent). Readings printed: "form ACCOUNTS" on 3 rows, "form + O(Δt)" on 4 (R_Π 1.05e-04 … 1.52e-04 at Δt₁ — ruling 5's conduction term, under the band at Δt₁/10), "INCONCLUSIVE" on 1 (M ccw bottom), "RESIDUAL" on none. So on this mesh pair the tip gap is wholly the comparand's form plus an O(Δt) step term; nothing Δt-independent is left beyond the projection. **Still not measured:** the same split on U127's exact-rectangle mesh or on F-clinical (entry 55) — this is the ragged-sheet I41 pair only. Hypothesis as filed, kept for the record: the gate is \|ρc·⟨θ¹⟩_B / (Δt₁·⟨q⟩_B) − 1\|, with ⟨·⟩_B a 1 g ball mean and q the DG0 cell average of ½σ\|E\|². In the Δt → 0 limit, one backward-Euler step gives ρc·θ¹/Δt₁ → M⁻¹b: the consistent-mass L2 projection of q onto the P1 temperature space. That projection preserves the whole-domain mean, so (S2) is green, but it does not preserve a ball mean. On a DG0 source with cell-scale jumps (a degree-1 N1curl \|E\|²), the gap is Δt-independent: this is U127's control. The uniform-source stub cannot see it (I32a's K10-S tip 8.04e-06, `20261101T103929Z_GEO-36-I32a-W0b.log:76`). The gel twin's miss says the gap is not about the rod. Both readings are open: the form of the comparand, or physics the gate is meant to catch. `OPS-83` (§9 item 859) separates them on the stub. |
| **Scope** | Every (S2-tip)-style gate on a field-derived DG0 source: `test_geo36_rod_f_clinical_heating.py` (entry 55) and `WF-15` step 14's `Study` route; a grep of `tests/` at this review finds the tip form in the first module only (U127's is on its attempt branch). **Not affected:** P61's modules (`test_geo38_f_clinical_box_heating.py`, its 3 T twin — entries 59 and 60) assert only the global (S2) (`:403`), and entry 61 is EM-only. Entry 55's EM record and θ figures stand under its pre-registered decision rule ("an energy or slope identity red with the EM gates green ⇒ a known-issues entry on the thermal route at clinical scale, the EM record stands"). |
| **Resolves with** | The 11-07 weekly's re-pose of the tip gate, for U127 and I32 together: (a) a band, or (b) the mass-consistent comparand ⟨Πq⟩_B at `S2_BAND`, with ⟨q⟩_B printed. `OPS-83`'s measurement is the evidence for choosing between them. Never loosen `S2_BAND` in passing. |
| **Ruling (2026-11-07 weekly; row written by the 2026-11-08 03:00 daily)** | **Option (b).** (S2-tip) is asserted against the consistent-mass projection ⟨Πq⟩_B at the unmoved `S2_BAND`, on one BE step at Δt₁/10, with ⟨q⟩_B, the ratio ⟨Πq⟩_B / ⟨q⟩_B − 1 and the Δt₁ read printed beside it (PROJECT_PLAN §10, rulings 2026-11-07, entry 55). **This entry resolves when all three land green in that form:** §10 **I32b** (the module and the shared thermal helper, stub-run; §9 item 870), **U127′** (`WF-15` step 14, attempt 2 of 4; §9 item 884, serial on 870), and **`xl` entry 63** (I32′, Friday 11-13, PENDING on I32b). `S2_BAND` does not move. Entry 55's EM record and θ figures stand as records meanwhile. |
| **U127′'s half (green, 2026-11-08 16:30 slot, §9 item 884)** | On U127's exact-rectangle rod mesh through the `Study` (drive P1 at 1 V, degree 1), R_Π at Δt₁/10 is **1.3394e-05 / 9.6264e-06** (top / bottom) ≤ `S2_BAND` 1e-04, asserted; R_Π at Δt₁ 1.3407e-04 / 9.6450e-05 (Δt-ratio 10.010 / 10.019, O(Δt)); the retired form R_DG0 4.2166e-03 / 6.8113e-04 at Δt₁/10 (attempt 1's reading reproduced) and ⟨Πq⟩_B/⟨q⟩_B − 1 = −4.2033e-03 / +6.9077e-04 — the whole of the DG0 gap; (B-vol) helper ball volume = the `Study`'s to 0 relative, quadrature 16 = 16 (`20261108T223537Z_WF-15-s14p-W1.log:5677–5679`). So the exact-rectangle mesh's tip gap is the comparand's form plus an O(Δt) step term, as on `OPS-84`'s ragged pair. **Still open:** `xl` entry 63 (I32′, F-clinical); this entry retires when it is green in this form. |
| **Entry 63's half (RED, 2026-11-13 `xl` window; filed by the 03:00 daily review from a `log-pathologist` reading)** | `20261113T080008Z_GEO-36-step10-I32-prime.log`: Status 1, 995 s, 43.97 GiB, orphans 0 / 0. Every EM gate, (B), (E′), the global (S2) and (S0) green; θ(t\*), probes and the rod share equal entry 55's at every printed digit. Three reds (`:12532`): **(S2-tip)** in the projected form — the rod's R_Π at Δt₁/10 is 1.0920e-04 / 1.1039e-04 against `S2_BAND` 1e-04, falling 10.08 / 9.82× from Δt₁ (1.1012e-03 / 1.0839e-03, `:12521–12522`): the O(Δt) step term, ≈ 25–50× larger than the gel twin's (2.7e-06 / 4.5e-06 at Δt₁/10, printed) — the rod's tip balls touch its adiabatic wall. The DG0 form gap is removed (⟨Πq⟩/⟨q⟩ − 1 −1.36e-03 / −5.34e-03 accounts for entry 55's R_DG0). A two-point intercept is −1.0e-06 / +2.2e-06 under a linear model only; a zero Δt-independent part is **not shown**. **And a new red, (E) 9.862e-09 / 9.220e-09 on both arms:** the residual of I32b's added BE step at Δt₁/10, which `out["E"]` now includes (module `:354`, `:418–420`); the printed CN / Δt₁ (E) are green, so this is not the assertion entry 55 passed. Class (pathologist, medium): finite-precision floor — the one-step balance is an exact identity under direct LU, T carried absolute from 20 °C, floor ≈ ε·T∞·ρc·V/(Δt·∫q) ≈ 1.6e-08 at Δt₁/10 at the clinical box's ≈ 40.7 W/m³ per watt (the 11-10 18:00 ruling's floor, `WF-17.md:3465–3478`; I32b's stub, at ≈ 25× the source density, read 8.804e-11). The 23× / 28× growth over Δt₁ is unexplained. **So on this route at clinical source density the two one-step gates conflict:** a smaller Δt clears the tip term and raises the (E) floor. Separating measurement named (not run): the same one-step BE on this submesh with `t_inf = t_init = 0`. **Owner: the 11-14 weekly** — the pre-registered rule has no branch for this pair. `S2_BAND` and `E_BAND` do not move. **Ruled 2026-11-14 (weekly, PROJECT_PLAN §10):** (S2-tip) takes its pre-registered branch — the reading is recorded and the tip read drops to printed on the clinical body; the (E) red is the Δt₁/10 step I32b folded into the asserted maximum, a reading of the absolute-T floor by the 11-10 18:00 ruling's class, so the asserted (E) returns to entry 55's definition (CN + the Δt₁ BE step) with the Δt₁/10 residual printed beside its floor estimate (§10 I32c, a slot step; its stub also runs the T-from-0 separating measurement); EM and θ records stand as entry 55's. **Resolves with:** I32c landed and `xl` entry 71 (I32″, 11-17 at the earliest) read green — this half stays open until then. No band moves. |

### 🟡 OPEN 2026-10-07 (`GEO-34` step 1d W0, 03:30 implementer slot, §9 item 289): Gmsh 4.15.2's OCC STEP writer drops entity names and does not keep an embedded face

| | |
| --- | --- |
| **Symptom** | `tests/mesh/probe_geo34_step_names.py` (a probe, not collected): a unit box fragmented with a 0.8 × 0.6 rectangle at z = 0.47. The face is named `PortSheet` and the volume `AirBox`, each by `setEntityName` and by a named physical group. `gmsh.write(….step)` writes 396 entities, and neither name string occurs in the file. `occ.importShapes` returns 1 volume named `Shapes/Open CASCADE STEP translator 7.6 1/…`, 6 unnamed faces of area 1.0, no physical groups, and nothing embedded. The 0.48 m² face is absent before and after `removeAllDuplicates` (`20261007T034333Z_GEO-34-step1d-W0-names.log:44–74`). |
| **Verified at** | `e4b2f445` plus the probe, by execution: one rank, 1 s, Status 0. |
| **Cause** | Not diagnosed beyond the reading. The writer exports the solid without the fragment's internal face and without names. Whether a different export option or a compound/shell export would keep them was not measured. |
| **Scope** | No test on `main` is red. It blocks `GEO-34` step 1d (§9 item 289, §10 U12): Gmsh's own writer cannot author a named-interior-face STEP. **Addendum (U12a probe, §9 item 333, 2026-10-08):** the form was measured, and it does not help. Five export forms were each read back two ways. The fragment and the solid alone carry no 0.48 m² face. A compound file (solid plus a free rectangle) keeps the face only under `highestDimOnly=False`, and two files merged in one read keep it under both reads. In every form the face is **free, never embedded** after `removeAllDuplicates` (`getEmbedded` is `[]`). No form writes a name string. Evidence: `tests/mesh/probe_geo34_step1d_export_forms.py`, `20261008T214725Z_GEO-34-step1d-probe.log`, and a bit-identical repeat at `…214738Z_GEO-34-step1d-probe-repeat.log`. The entry stays open because names are still dropped. |
| **Resolves with** | A STEP from a kernel that writes named, embedded faces (`cq-python`, dashboard item 10, `GEO-34` step 1c), then the weekly re-posing U12. Never a geometric-selector substitute. |

### 🟡 OPEN 2026-09-28 (`WF-10` step 1b W1, 21:00 implementer slot, §9 item 62) — at the k = 1 block's `C_mode` = 113.108 pF, the F-human in-model vs circuit 2×2 identity reads **3.4865e-03 against `REDUCTION_BAND` 1e-3**, where step 1 read 6.99e-04 at 19.7327 pF on the same fixture

| | |
| --- | --- |
| **Symptom** | `test_wf10_step1b_mode_tuning` (on `attempt/WF-10-step1b-20260929T020339Z`, not on `main`), `FEM_EM_WF10_STEP1B=1`, `-n 8`, complex, degree 1: `identity residual ||MEAS - PRED||/||PRED|| = 3.486500e-03 (band 1e-03, ASSERTED)` (`20260929T020339Z_WF-10-step1b.log:10648`). Same fixture (507 266 cells), κ 4.590991e-3, 32×32 reciprocity 8.53e-15 and σ_max 0.999923 green. The feed-ring k = ±1 share also misses its 0.9 floor (0.7898 in-model, 0.7882 circuit, `:10649`). |
| **Verified at** | `33ad054` plus the parked test, by execution: one window, 451 s, rc 1. |
| **Cause** | **Not diagnosed.** The residual is 5× step 1's at a capacitance 5.7× larger; whether the κ width correction, the sheet model at larger C, or proximity to a resonance of the terminated network drives it was not measured. **Step 1c (2026-10-01 06:00 slot, §9 item 128) measured it — reading "mixed", not classified:** over six points 19.73 → 113.11 pF the residual grows ≈ ∝ C (log-log slope +1.03; +0.95 over anchors and off-points) with a ≈ 1.7× bump at the upper cond peak (105.02 pF, 6.106e-03) but **none** at the lower one (75.70 pF, next to the k = 2 anti-mode, 2.017e-03); residual / cond spans 8.7e-05 … 2.9e-04, so the residual does not track cond. Rule-sensitivity flagged: the upper off-point lands 0.09 % below `C_mode`, so the monotonicity test turns on a 3.51e-03 vs 3.49e-03 difference (`20261001T110240Z_WF-10-step1c-W1.log`, `20261001T110957Z_WF-10-step1c-W2.log`). Entry stays OPEN; a band argument is the review's or the operator's (ruling (b)). **Step 1d (2026-10-01 09:00 slot, §9 item 133) adds the parallel-drive reading:** with C on all 32 gaps and P17 / P21 at Z_th = 50·Z_C/(50 + Z_C) at the same `C_mode` (reproduced to 4.3e-13), the identity reads **2.399e-03** (`20261001T140319Z_WF-10-step1d-W1.log:10658`, `test_wf10_step1d_parallel_drive` on `attempt/WF-10-step1d-20261001T140319Z`) — 0.69× the H8-configuration 3.49e-03 at the same C, still 2.4× the band; the two circuit routes agree to 9.1e-15 (`:10654`), so the gap is circuit-vs-field, not circuit algebra. Not diagnosed further. **Copper reading (`WF-12` step 2′, H12b′, §9 item 175, 2026-10-02 13:30 slot) — pre-registered classification printed "MIXED", not re-classified:** on H12b's copper hole fixture (369 387 cells, Leontovich Cu wall) in H8a's configuration at **uncorrected** widths (κ off per the operator ruling 2026-10-02; the σ = 800 comparands were κ-corrected, so ρ compares across that difference), the identity residual at the six shared C reads r_Cu = 1.726e-02 / 6.262e-02 / 6.587e-02 / 2.238e-01 / 2.798e-01 / 2.747e-01 at 19.73 / 75.70 / 88.35 / 105.02 / 113.01 / 113.11 pF, **ρ = r_Cu / r_800 = 24.7 / 31.0 / 18.4 / 36.6 / 79.6 / 78.8** (predicted [0.5, 2] — out), log-log slope +1.496 (σ = 800 +1.033; predicted 1.03 ± 0.3 — out), r_Cu(`C_mode`) 2.747e-01 (predicted [1.7e-3, 7e-3] — out). residual / cond 1.37e-3 … 6.17e-3. Copper's own spectrum moved: `select_c_tuned` 95.91 pF (σ = 800 19.73), cond maxima at 35.29 / 37.68 pF (σ = 800 75.70 / 105.02), where r_Cu = 1.495e-01 / 1.680e-01. Every 2×2 gate stays green (reciprocity ≤ 6.8e-14, column Σ\|S\|² ≤ 0.979); the two circuit routes agree to ≤ 1.8e-14, so the gap is again circuit-vs-field (`20261002T183158Z_WF-12-step2p-W1.log:6492–6503`, `20261002T183901Z_WF-12-step2p-W2.log:6403–6416`). Not diagnosed; whether κ-off vs κ-on or the conductor drives the ×18–80 is not separable from this pair of windows. **Copper, H8b's configuration (`WF-12` step 3, H12c, §9 item 261, 2026-10-06 04:30 slot; parked `attempt/WF-12-step3-20261006T0933Z`):** at copper's own `C_mode` = 118.1186 pF, with C on all 32 gaps and P17 / P21 at Z_th, the identity reads **1.581167** (`20261006T093326Z_WF-12-step3-W1.log`; every other anchor green, the two circuit routes 1.613e-13). A 1e-4 relative C change moves the circuit 2×2 by 9.3e-3, and no uniform circuit-side C shift within ±10 % brings it below 1.432 (`20261006T093934Z_WF-12-step3-diag.log:35–43`). **`log-pathologist` ruling (2026-10-06 10:30 review):** no defect found in the module (port order, Z_C, Z_th, held factor and lifted helpers all checked). But H12b′'s curve above is **H8a's** configuration (feeds at 50 Ω), so it is not the same quantity, and the item's `PREDICTED RED` range drawn from it did not apply. The in-model route has never been tied back to the 32×32: its closure at 50 Ω is the discriminator (§9 item 285). **(Z50) closes (`WF-12` step 3 attempt 2, §9 item 285, 2026-10-06 12:00 slot; parked `attempt/WF-12-step3-a2-20261006T1702Z`):** the in-model pair with all 32 sheets at 50 Ω reproduces s32's 2×2 to **6.761e-14**, and the identity at copper's `C_mode` re-reads **1.581167** (3.5e-13 relative to attempt 1), so the O(1) residual is a property of the copper network at `C_mode` in H8b's configuration, not of the in-model route. Shared-C point: H8b's configuration at 113.1081 pF reads **3.519e-01** (printed; H12b′'s 2.747e-1 there is H8a's configuration). Cause still not diagnosed (`20261006T170222Z_WF-12-step3-a2-W1.log`). **Termination fit (`WF-12` step 3b / H12d, §9 item 313, 2026-10-08 04:30 slot; stored npz only, no solve; green):** perturbing every port termination Z_i → Z_i(1 + δ_i) closes the `C_mode` miss — fit ratio **8.0e-14** of ‖CIRC(0) − MEAS‖_F = 0.7167, **max \|δ̂\| 1.3006e-02** (all 32 within 1.12e-2 … 1.30e-2, Re ≈ +0.1 … +0.9 %, Im ≈ −0.95 … −1.21 %), k-spectrum **k = 0**-dominated (0.9897 other ring, 0.9618 feed ring; k = ±1 ≤ 0.0042); cross-check δ̂ on the 113.1081 pF H8b pair **0.3519 → 0.0180** (printed); a uniform complex δ on the 30 gaps alone leaves 1.60e-2 of the miss (predicted control). Reading: the termination hypothesis survives — the O(1) miss is a common-mode ≈ 1.2 % per-sheet impedance offset between the circuit model and the field, not route noise; *which* physical term it is (copper-wall current beside the sheet, a sheet self-term, a convention) is still not measured (`20261008T093959Z_WF-12-step3b-W1.log`). |
| **Scope** | `WF-10` steps 1b and 1d (both parked). Nothing on `main` runs the test; no band moved. |
| **Resolves with** | The weekly's re-scope of H8′ (§9 item 62, BLOCKED). **Operator ruling 2026-09-30:** both bars kept; the share floor is read on C-on-every-gap with parallel drives (PROJECT_PLAN §10 H8′); the band stays 1e-3 until a diagnostic window sweeping C across 19.7 → 113 pF, the identity printed at each point, separates a C trend from a resonance. |

### 🟡 OPEN 2026-09-11 (`WF-6` step 4h, 09:00 implementer slot) — the unloaded F-small birdcage's ×0.0095 global-resolution rung **misses power accounting by ≈ 1.97e-02 on both ports even with the C4-congruent sheet cut**: the cut removes the P1/P2 split, not the residual

| | |
| --- | --- |
| **Test** | `tests/validation/test_birdcage_b1_plus_closed_form.py::test_power_accounting_closes_on_every_rung[x0.0095]`. **Reachable only opt-in** (`FEM_EM_WF6_C4_CONGRUENT=1`); the default `LADDER` (`x1`, `x0.012`) does not contain the rung, so nothing on `main` is red by default. |
| **Verified at** | `bf6ae76` + the step-4h test-side keywords. `20260911T140329Z_WF-6-step4h.log` — `-n 4`, complex, `FEM_EM_REQUIRE_COMPLEX=1`, keys `x1 x0.0095`, `timeout -k 30 590`, `-v -s`: **1 failed / 18 passed / 4 skipped in 196.74 s**, elapsed 199 s, `[capture] rc=1` (`:4149, :4372`). Default-path control `20260911T140714Z_WF-6-step4h-flagoff-x1.log`: 6 passed, rc 0. |
| **Symptom** | `rung x0.0095 [P1]: power accounting misses by 1.968410e-02 of the supplied 6.789015267e-03 W (phantom 0.000000000e+00, conductor 4.371025298e-04, sheets 6.218277077e-03); imported band 1e-02` (`:4147`). P2 reads **1.968464e-02** (`:4035`). The flag-off step-4f record on the same rung was P1 **1.853642e-02** / P2 **1.419812e-02** (`20260909T123716Z_WF-6.log:5777–5778, 5789–5795`), so on/off = 1.0619 / 1.3864 (`:4036–4037`). |
| **What is green in the same window** | The flag-on ×1 rung: 116 118 cells, residuals 9.795942e-03 / 9.796517e-03, covariance 3.5271%, cw/ccw 15.03× (`:2032–2044`). On ×0.0095 itself: covariance 1.3181%, cw/ccw 53.09×, 21 of 21 points (`:4031–4033`). |
| **Cause — not diagnosed** | The congruent cut **was** the carrier of the port-to-port split: P1 and P2 now agree to 2.7e-05 relative, as they do on both coarser rungs. The residual itself is **common-mode** and did not fall, so it is not a sheet-cut asymmetry. Unseparated candidates include the accounting's own terms at finer `h` (the conductor-less residual rose from 7.52e-02 at ×1 to 8.41e-02 at ×0.0095, `:2043, :4034`) and the sheet-power estimator's discretisation. Not measured in this slot. |
| **Consequence** | ×0.0095 stays out of the default `LADDER`; no `WF-6` record, band or rung is re-registered; `WF-6` stays 🟡. |
| **Resolves with** | A review scoping a diagnosis of the common-mode ≈ 2e-2 residual at ×0.0095 (for example, a per-term accounting ladder across ×1 / ×0.012 / ×0.0095 with the flag on). `POWER_BALANCE_BAND` is not to be widened. |
| **Ruled 2026-09-11 10:30 review — accepted as a negative result; the terminal-form deficit is the first candidate** | Re-read at `:4034–4035` (×0.0095) and `:2043–2044` (×1). The rung stays out of `LADDER`. `PORT-16` attributed the loaded ×1 fixture's ≈1 % gap to the Cauchy–Schwarz deficit of `½\|I\|²Re Z_p` below the sheet field's dissipation (`tests/validation/test_birdcage_power_identity.py` docstring (iv)), a term a finer rim can grow, so §9 item 3 (`WF-6` step 4i) prints the exact shares beside the terminal residual on the flag-on rungs, asserting the discrete identity at 1e-6. **Status caveat:** this window's harness footer reads `- Status: 0` (`:4375`), and so does its test-results row, over `[capture] rc=1` (`:4372`), because the capture command dropped `; exit $rc`. The rc line is authoritative, and `OPS-45` makes the harness say so. |
| **Ruled 2026-09-11 18:00 review — cause located, mechanism not; entry stays OPEN** | `WF-6` step 4i (`20260911T200251Z_WF-6-step4i.log:4057, :4061`, re-read here) shows the terminal residual is exactly `(C − terminal)/supplied`, the Cauchy–Schwarz deficit of the terminal sheet form: `C/terminal` = 1.010592 at ×1, 1.008756 at ×0.012 (`20260911T200634Z_WF-6-step4i-x0.012.log:2018`) and **1.021491 at ×0.0095**, uniform across the four sheets to ±1e-6. **The accounting is not switched to the exact form:** `supplied_terminal = P_vol + C` closes to ≈1e-15 on every rung because it is an algebraic identity of the discrete solve, so that gate could never fail. The terminal residual stays the gate, ×0.0095 stays out of `LADDER`, and `POWER_BALANCE_BAND` does not move. Still undiagnosed: why the sheet field's non-uniformity doubles on ×0.0095 while ×0.012 lowers it (strip-edge field resolved by more facets across the width, or a rim-cell artefact of that mesh). `WF-6` step 4j (§9) localises it. Resolves with that reading plus a review ruling on what the rung's power gate should be. |
| **Measured 2026-09-11 22:30 slot (`WF-6` step 4j) — located on the sheet; entry stays OPEN** | `20260912T033518Z_WF-6-step4j.log` (keys `x1 x0.0095`, `-n 4`, 1 failed [this entry's test, by design] / 24 passed / 4 skipped in 198.02 s, elapsed 199 s, `[capture] rc=1`) and `20260912T033918Z_WF-6-step4j-x0.012.log` (17 passed / 3 skipped in 99.33 s, elapsed 101 s, rc 0). C/terminal reproduces 4i to ≤ 3.5e-7 on every rung and drive. The sampled `R_s,tan` equals C/terminal to ≤ 3e-8 and splits C/terminal − 1 into two parts. The drive-component variance is 8.145e-3 / 6.691e-3 / **1.2045e-2** (×1 / ×0.012 / ×0.0095). The in-plane transverse term is 2.447e-3 / 2.065e-3 / **9.446e-3**, so it carries 64 % of the ×1 → ×0.0095 growth. On ×0.0095, one facet (mid-gap, v 0.42; width edge, u 0.16; 6.7 % of the area) carries 42 % of the drive variance and has the largest in-plane field. The gap-end bins carry 17 % on 41 % of the area. Every rung has only 26–29 facets per sheet and ≈ 3.5–3.8 across the width. **The evidence favours a rim-facet artefact over a resolved strip edge or the terminals.** The mechanism for the in-plane field on that facet is not measured. |
| **Ruled 2026-09-12 03:00 review — 4j accepted; a no-solve census queued; entry stays OPEN** | Re-read `20260912T033518Z_WF-6-step4j.log:4105, :4108`. The facet at (u 0.164, v 0.422) is 6.73 % of P1's sheet area and carries 42.25 % of the drive variance with `<\|E_w\|>` 20.9 V/m against 3–8 elsewhere. `WF-6` step 4k (§9 item 5, `mesh-probe`) measures that facet's aspect ratio, adjacent-tet quality and rim/terminal adjacency on the three flag-on rungs (116 118 / 148 988 / 197 284 cells) without a solve. The rung's power gate is not ruled until that reads; `POWER_BALANCE_BAND` does not move. |
| **Measured 2026-09-12 16:30 slot (`WF-6` step 4k) — the sliver hypothesis is refuted on shape; entry stays OPEN** | `tests/mesh/probe_wf6_sheet_facet_census.py`, no solve, `-n 1` complex, one process per rung: `20260912T213317Z_WF-6-step4k-x1.log`, `20260912T213407Z_WF-6-step4k-x0.012.log`, `20260912T213439Z_WF-6-step4k-x0.0095.log`, `20260912T213518Z_WF-6-step4k-x0.0095-repeat.log` (all Status 0). Cell and P1 facet counts reproduce 4j. The (0.164, 0.422) facet reads area share 6.726 % (`…-x0.0095.log:1928`), and the repeat build is hash-identical (`:2007`). On ×0.0095 the facet is among the best-shaped on its sheet: edge ratio 1.195 and circumradius/inradius 2.047 rank 5/29, below Q1. Its worse adjacent tet reads 3r/R 0.7864 and min dihedral 46.08°, rank 23/29, above Q3. It is the largest facet on the sheet (2.174× median) and has a lateral-rim edge, as do 8 of 29 facets, but no terminal edge (`:1928–1933`). The top-1 facet on ×1 and ×0.012 is the same kind of large mid-gap lateral-rim facet. **Not a sliver; the item's "not a rim facet" clause is not literally met** — the in-plane field on that facet needs a field-side reading (the N1curl trace on its adjacent cells), which a review scopes. |
| **Ruled 2026-09-12 18:00 review — 4k accepted; the step-4 family is frozen; entry stays OPEN as a banked negative, no further slot goes to it** | Re-read `20260912T213439Z_WF-6-step4k-x0.0095.log:1928–1933`: the facet ranks 5/29 on edge ratio and R/r, 23/29 on adjacent-tet 3r/R and dihedral, 29/29 on area, lateral-rim edge, no terminal edge. Not a sliver. Under the four-attempt family cap (operator directive 2026-09-12, enacted by this review in §9) the `WF-6` step-4 family — 4h, 4i, 4j, 4k since 4g's green — is frozen: the four attempts all diagnosed a rung already removed from `LADDER`, so none could move a status. The field-side reading 4k proposes is **declined**, not queued. What retires this entry is one of: the 2026-09-13 weekly putting ×0.0095 back on the ladder with its residual re-registered by measurement under rule (f) (a numbered `WF-6` step, not a 4l), or a decision that the two-rung ladder on `main` is the chunk's convergence statement and the ×0.0095 reading stands as provenance only. The mechanism stays **not explained**; `LADDER`, bands and accounting unchanged. |

### 🟡 OBSERVATION (was OPEN) 2026-09-10 (`OPS-43` (d) gate, 06:00 implementer slot) — the complex MUMPS smoke solve is **not bit-reproducible across processes at `-n 2`**: 2 of 22 identical child solves drifted by 1 ULP. **Re-headed 2026-09-10 13:30 slot: no gate depends on this any more** — the (d) gate landed on `main` re-anchored at child `-n 1` (`20260910T183305Z_OPS-43d.log`, `Status: 0` `:245`, 54 s). Not retired: the drift is still undiagnosed.

| | |
|---|---|
| **Gate that no longer depends on it (2026-09-10 13:30 slot)** | `tests/solver/test_solver_progress_inert.py::test_solver_progress_variable_is_bit_inert_at_n1_and_visible_only_when_set`, on `main`. It asserts one `float.hex` value across 4 unset + 4 set interleaved child `-n 1` solves: norm2 `{'0x1.2bbe0e158fdb3p-5': 8}` (`:90`), three-term imbalance `{'0x1.56f8034472e19p-3': 8}` (`:91`). The `-n 2` children are printed only. In this window 12 of 12 read norm2 `0x1.2bbe0e158fda0p-5` (6 unset, 6 set) and three-term `0x1.56f8034472e29p-3` (`:108–111`), so there was no drift this time and the running tally is **2 of 34**. |
| **Test (original, 06:00)** | `tests/solver/test_solver_progress_inert.py::test_solver_progress_variable_is_inert_on_result_and_visible_only_when_set` — the parked `-n 2` single-pair design, superseded by the test above. Fixture = `tests/solver/test_time_harmonic_smoke.py`'s cylinder at h = 0.03 (1 405 cells, 2 004 dofs), axial drive, degree 1, default MUMPS options; each solve a separate child `mpiexec -n 2` launched from rank 0 with a PMI-scrubbed env. |
| **Symptom** | `20260910T110332Z_OPS-43d.log` (footer `Status: 1` `:388`, elapsed 36 s `:389`): `AssertionError: FEM_EM_SOLVER_PROGRESS is not inert: norm2 unset 0.03658964873147431 (0x1.2bbe0e158fda1p-5) vs set 0.0365896487314743 (0x1.2bbe0e158fda0p-5), difference -6.938893903907228e-18` (`:244`); 6 of 8 compared scalars differ at 1 ULP-scale (`:230–237`), e.g. three-term imbalance `0x1.56f8034472e2cp-3` vs `…e29p-3`. |
| **Why it is not the variable (measured)** | Same module, three further windows, 4 then 4 then 20 children: `20260910T110452Z_OPS-43d.log` (`Status: 0` `:383`, 37 s) all four -n 2 children (unset, set, unset, set) = `…fda0p-5`, i.e. window 1's **unset** child was the outlier; `20260910T110629Z_OPS-43d.log` (`Status: 0` `:383`, 37 s) the **set_repeat** child drifted to `…fd9fp-5` (`:243`) while mesh-geometry sum `-0x1.12f6eb68be2c0p-1`, `‖A‖_F` `0x1.35c8f773e9117p+14`, `‖b‖₂` `0x1.afb410f1f7f48p-3` and cells/rank `[692, 713]` were bit-identical to the other three — the drift enters in the factor/solve, not the mesh or assembly, and in both settings; `20260910T110803Z_OPS-43d.log` (`Status: 0` `:403`, 61 s, `FEM_EM_OPS43D_PROBE=8`) 12/12 `-n 2` children `…fda0p-5` and 8/8 `-n 1` children `…fdb3p-5` (`:269–270`). Tally at `-n 2`: **2 drifted of 22** (one unset, one set); at `-n 1`: 0 of 8 (small sample). |
| **Cause** | **Not diagnosed.** Localised to after assembly (A/b fingerprints identical) on MUMPS 5.8.2, `#MPI = 2, without OMP`, `OPENBLAS_NUM_THREADS=1`. Candidates unseparated: a timing- or alignment-dependent path in parallel MUMPS / its BLAS kernels. |
| **Consequence** | A single-pair bit-identity on this solve cannot gate "inert" — it would fail roughly one pair in six in CI for reasons unrelated to the variable. **Not relaxed to a tolerance** (§9 item 2 forbids it). The (d) gate's other legs were green in every window: `[solve]` lines 0 unset / 2 set; options dict unset == pre-`81861d0` baseline literal, set == baseline + `mat_mumps_icntl_4: 2`; MUMPS `ICNTL(4)` read back `[0, 0]` / `[2, 2]`. (d1)/(d2) stay in place. |
| **Resolves with** | A review ruling on the anchor's design, e.g. bit-identity at child `-n 1` (8/8 reproducible so far), or majority/mode over N repeats per setting, or a discriminator that compares the unset-vs-set spread against the same-setting repeat spread. Not an in-slot decision. **Ruled 2026-09-10 10:30 review:** bit-identity at child `-n 1`, asserted as one `float.hex` value across 4 unset + 4 set interleaved children; `-n 2` children printed only (PROJECT_PLAN §7 `OPS-43`, §9 item 1). When that gate lands green this entry stops being a *failing-test* entry. The drift itself stays open as an observation until its cause is diagnosed. |
| **Exposure census (2026-09-10 10:30 review)** | An `Explore` sweep of the 167 modules under `tests/` finds **no** assertion comparing a MUMPS solve result against a hard-coded record tighter than **1e-9** relative. Tightest found: `LEG_D0_REPRODUCTION_BAND = 1.0e-9`, `tests/validation/test_port_birdcage_four_port.py:165`, spot-checked by the review. The only bit-exact comparisons on solve output (`np.array_equal`, `tests/validation/test_port_drive_superposition.py:736, 825`) are between quantities from the same solve in one process, so run-to-run drift cannot reach them. A 1-ULP (≈1e-16) drift therefore exposes no existing record. Any future record at `-n 2` tighter than ≈1e-14 relative on a MUMPS solve would be exposed. |
| **Second fixture (2026-09-29, `MAT-1` §9 item 81; copied by the 18:00 review so the data outlive the `MAT-1` entry)** | On `MAT-4` step 1's lossy sphere (74 020 cells, `-n 2` complex, one process, identical inputs), mean SAR / ⟨E_z⟩ drift between repeated solves by up to **4.45e-11 / 1.93e-11** relative threaded and **4.24e-11 / 2.14e-11** at OMP / OpenBLAS threads = 1 (`20260929T183957Z_MAT-1-probe-noise.log:105–106`, `20260929T184102Z_MAT-1-probe-omp1.log:106–107`) — within one process, five orders above this entry's cross-process 1 ULP. **The census's exposure threshold tightens accordingly: on lossy fixtures, any solve-output record tighter than ≈ 5e-11 at `-n 2` is exposed.** Still undiagnosed. |
| **A gate inside the exposure line (2026-11-02 03:00 review, from `OPS-74` and `EX-172`)** | `IDENTITY_BAND` = 1e-12 compares two separately built and factorised solves in the `WF-15` / `EX-167` (S) gate, on the σ-loaded rod phantom at `-n 2`: inside the ≈ 5e-11 line above. Read so far: 1.295e-15 (`20261101T212457Z_EX-167.log:5571`), then 6.598e-14 (`20261102T013649Z_EX-167-close.log:5568`), ≈ 51 ×, with no `src` or script change. `EX-172`'s 10 MHz `C_mode` also moved between windows of one slot: `0x1.53d78952a75e9p-29` (`20261102T062150Z_EX-172-figure.log:1907`) against `0x1.53d78952a7944p-29` (`20261102T062422Z_EX-172.log:1911`), rel ≈ 1.4e-13. **No gate has failed on it, and no band moves.** `OPS-77` (§9 item 755) measures the distribution; `(I=)` hex comparisons within one run are unaffected, and a hex record across runs is safe only on a saved file (`ITEM730_HEX`). |
| **Third fixture: the (S) distribution (2026-11-02, `OPS-77`, §9 item 755)** | On the same σ-loaded rod phantom at `-n 2`, five samples, each its own process, read by `scripts/testing/read_s_identity_record.py` (`20261102T122110Z_OPS-77-reader.log:36–44`): 1.295e-15, 1.357e-15, 6.598e-14, 9.801e-14, **1.147e-13** (`20261102T121721Z_OPS-77-mri37-r2.log:5568`). The two newest windows are identical and read ≈ 85 × apart, 1.357e-15 then 1.147e-13. The largest is 0.115 × `IDENTITY_BAND`, headroom 8.7 ×, so **the tail branch fired**: one read in (1e-13, 1e-12], none above. All five windows were green. Five samples set no failure rate, and the cause stays undiagnosed. `mri:38`'s (S) reads 0.000e+00 at 64 and 128 MHz (`20261102T060640Z_EX-170.log:5643`, `:5646`). **Corrected 2026-11-02 18:00 review: this row called that "the floor when both sides share a route", but those are two factorisations (`20261102T180848Z_OPS-78-census.log:85–86`), so it is a zero read inside the exposed class and explains nothing about `mri:37`.** A neighbour on another fixture, `WF-16`'s four-port (I), reads 1.004e-13 (`20261101T105341Z_WF-15-s10-Wreg.log:3716`). **No band moves;** whether the gate's design changes is the weekly's. |
| **The exposed set (2026-11-02, `OPS-78`, §9 item 767)** | Every use of a `*_BAND` name ≤ 1e-11 classified by reading the call (`20261103T013749Z_OPS-78-census.log`, `Status: 0`, counts equal `grep -rnw` per name, three controls as required, `:74–117`): **369 lines, 92 separate solve, 51 shared route, 217 not a solve output, 9 out of scope, 0 unread.** 21 separate-solve gates carry a row in `docs/planning/chunks/OPS-78.md`, and **none reads above its band.** Top three by headroom: `MAT-21` (H0) 3.801e-13 = 0.38 × 1e-12, 2.6 × (`20261014T003128Z_MAT-21-step2-a2-W2-identities.log:123`); `MAT-20` / `mat:18` (FEM-2d) 3.476e-14 = 0.35 × 1e-13, 2.9 × (`20261012T170911Z_EX-110.log:1920`); (FEM) 3.434e-14, 2.9 × (`20261011T172810Z_MAT-20-step2-fem-W1.log:1969`). The tightest port gates are `WF-17 s25` (P=Y) 1.180e-13 and (S) 1.172e-13, 8.5 ×. The tight rows are thermal, and whether their solver is the factorisation this entry localised was not read. No cause, failure rate or design is claimed, and **no band moves.** |
| **Weekly ruling 2026-11-02 — the gates' design stands** | `IDENTITY_BAND` 1e-12 across two factorisations sits ≈ 9× above the worst measured (S) sample (1.147e-13, `OPS-77`) and every one of `OPS-78`'s 21 separate-solve gates is ≥ 2.6× inside its band; no band moves and no gate changes design. **Rule:** a future red on a separate-solve (S) gate is read first against `OPS-77`'s five samples and `OPS-78`'s table — a red inside the measured drift is this observation firing and is filed here; a red outside it is a defect. `mri:37` non-zero 5 of 5 against `mri:38` zero 6 of 6 is a recorded difference between two builders with no claim attached (PROJECT_PLAN §10, rulings 2026-11-02). |
| **The exposed set, refreshed at HEAD (2026-11-04, `OPS-78`, §9 item 815)** | The census re-run on `59784023` (`20261104T180148Z_OPS-78-census-r2.log`, `Status: 0`, `:75–145`): **402 lines, 100 separate solve, 51 shared route, 242 not a solve output, 9 out of scope, 0 unread**, the three controls as required. Three new separate-solve gates, none above band: `WF-17 s28` (Id) 6.377e-14 (15.7 ×), `WF-15 s13` (Id) 4.722e-16, `WF-16 s18` (I) 0.000e+00 on 5 reads. The top three by headroom are unchanged. **`mri:38`'s (S) is no longer zero on every read:** `EX-170`'s 2026-11-03 re-run reads 9.342e-14 and 1.076e-13 (`20261103T225531Z_EX-170.log:5643`, `:5646`), so it is now 8 reads, 6 exactly 0. The ruling row above says "zero 6 of 6", which was the record then. This is recorded for the weekly with no claim, and **no band moves.** |
| **A pointwise read at the 1e-9 line (2026-11-08, `WF-1` step 29 / U134, §9 item 873, parked; filed by the 10:30 daily review)** | On the `MAT-4` / `MAT-10` sphere at EM degree 1, `-n 2` complex, the YAML route's q_h and an independently solved Python Study's q_h differ by an **owned-max relative 3.47e-09 (A) / 9.52e-10 (S)** (`20261108T122244Z_WF-1-s29-W1b.log:310–311`, on `attempt/WF-1-s29-20261108T1224Z`). The integrated source power agrees to 2.9e-12 / 7.4e-13. One BE step carries that into θ, and (IT) S reads 4.500e-09 / 1.122e-09 against `L0_BAND` 1e-9 (W1 / W1b). `SOLVE_FLOOR_RECORD` (5.962e-13, 4.450e-11) records *integrated* quantities, so it does not cover this. The pointwise figure is ≈ 78× its upper end, and above `SOLVE_RECORD_FLAG` 1e-9 (`tests/materials/test_phantom_material_model.py:240`). **Not established** (`log-pathologist`, 10:30 review): that the pointwise spread is a solver floor rather than a route difference. **Resolves with** §9 item 891's (Q-rep) print, which compares two same-route solves pointwise in one process: if they also differ by ≥ 1e-9 the floor is real and gets its own record; then a weekly ruling on whether a pointwise comparand across two solves may gate at `L0_BAND`. No band moves. |

### 🟡 OBSERVATION (was OPEN; banked 2026-10-31 by the weekly review — the R1 family closed at attempt 3 of 4, the strict `xfail` is the record) 2026-09-09 (filed by `TH-15` step 3b, 22:30 implementer slot; the defect itself dates from `PORT-1` step 3b-xii, 2026-08-07) — `tests/validation/test_port_gap_voltage_padding.py::test_box_enlargement_discriminates_between_the_two_routes` is **red on `main` and has been since it landed**: the estimator/control deviation reads **3.111508e-02** at `air_padding = 0.10` against the pre-decided 2.5% threshold, which is `PORT-1`'s **disposition (ii)** — a real estimator bias, not a truncation artefact

> This module was **never observed in a footered run** before this slot.
> `OPS-26` step 2 findings 21 and 33 deferred it structurally — `Status 124`
> at `-n 2` / 400 s and again at `-n 2` / 590 s, `collected 2 items` and
> then the first test's name and nothing else, reason recorded as "module
> fixture alone > 590 s at `-n 2`", with the note that its `scope="module"`
> `gap_ports_padded` fixture makes a by-name split useless. `TH-15` step 3b
> ran it at **`-n 4`** and it footered in 537.36 s, so the deferral is
> retired by width and the red behind it is visible for the first time.

| | |
|---|---|
| **Test** | `tests/validation/test_port_gap_voltage_padding.py::test_box_enlargement_discriminates_between_the_two_routes`. The module's other **12** names pass in the same run. |
| **Symptom** | `AssertionError: the estimator/control deviation at air_padding = 0.1 is 3.111508e-02, above the pre-decided 2.5%: enlarging the PEC box from 0.08 did not pull the two routes together (padding 0.08 read 3.022400e-02)` → `assert 0.03111507676007086 <= 0.025` |
| **Verified at** | `20260910T033937Z_TH-15.log` — `-n 4`, complex build (`FEM_EM_REQUIRE_COMPLEX=1`), `tests/environment` first, `-s -v --tb=short`, `timeout -k 30 590`: **1 failed, 12 passed in 537.36 s**, `Status: 1`, elapsed **540 s**. Tree = `main` + `TH-15` step 3's additive `gap_cell_tags` field. |
| **Cause** | Not diagnosed further here, and **not this slot's change**: the module names no `gap_cell_tag` and the step-3 default selection is byte-identical, while the two printed deviations reproduce `PORT-1` step 3b-xii's 2026-08-07 record (−3.0188e-02 / −3.1267e-02 at padding 0.10, −3.0224e-02 at 0.08) measured a month before the field existed. The test was landed *carrying* this red by `a755afb` (`PORT-1` step 3b-xvii) as the pre-decided discriminator between disposition (i) (truncation) and (ii) (estimator bias); it lands on (ii). |
| **Fix** | **Deliberately not fixed, and `REACTION_CONSISTENCY_TOLERANCE` (3%) is not to be moved** — the test's own message says so. Owner is `PORT-1`'s successor to step 3b-xii. What this entry adds is that the red is now *measured on the 0.11 image at `-n 4`* rather than inferred, and that the module is runnable in a 590 s window at that width. |
| **Owner since 2026-10-24; shrink ruled 2026-10-26 (weekly, §10 R1′)** | Lane 1, R1: re-run under the owner rule (`POST-8`'s point-sampled-voltage fix); under 2.5 % ⇒ green by consequence and this entry retires, not ⇒ strict `xfail` with this entry as reason, never loosened. Items 635 / 637 (2026-10-25) were killed before the assert: the legacy module takes 537 s at `-n 4` and the owner rule plus load exceeds a 600 s foreground call (`20261025T170055Z_PORT-1-R1-W1.log:844`, killed at 9 of 10 solves). **R1′:** the padding assertion reads only 3 of the fixture's 10 solves (the two port-driven ones and the σ = 0 control, `test_port_gap_voltage_padding.py:107–125`), and the σ ladders do not use the owner rule — so an additive keyword on `gap_ports_padded` computes only those three per box, `-n 4`, predicted ≈ 150 s; the `owner=None` threading already on `attempt/PORT-1-R1-20261025T171120Z` (`9b8ed86f`) is the starting point. R1 family attempt 3 of 4. |
| **Re-measured under the owner rule 2026-10-27 (R1′, §9 item 666, 07:30 slot) — still red, now strict `xfail`** | `20261027T125439Z_PORT-1-R1p-W1.log` (`-n 4`, complex, `padding_only=True` + `owner="lowest_original_cell"`, 3 solves, **1 failed / 5 passed in 218.41 s**, elapsed 220 s; FAILED ×4 ranks): worst deviation **3.101935e-02** (`:728`; columns −3.0515e-02 / −3.1019e-02, `:726–727`) against the unmoved 2.5 % — the owner rule moves it by ≈ 1e-4 (3.111508e-02 without it). Negative control unchanged: wedge-only ratio 0.5181, deviation −0.4819 > 0.20. 195 837 cells; control Re Z21 = 0 exactly. **Rank ownership is excluded as the cause**; the red is disposition (ii), an estimator bias. The test now carries `xfail(strict=True)` naming this entry; `20261027T125835Z_PORT-1-R1p-W1x.log` shows XFAIL ×4 ranks, **5 passed, 1 xfailed in 207.73 s**, Status 0, elapsed 209 s. `CONVERGENCE_THRESHOLD` and every assert unmoved. Entry stays open: a strict XPASS will flag any future fix. |
| **Banked 2026-10-31 (weekly review, §10 rulings and epitaph)** | The owner rule moved the deviation by ≈ 1e-4 (3.101935e-02 against 3.111508e-02), so the 3.1 % is the padding estimator's own bias — disposition (ii) as filed — and the gap-voltage padding route has had no consumer since the lumped-sheet port landed (`PORT-9`, 2026-08-25). A fourth attempt would re-measure the same number; the family is closed at 3 of 4 and this entry is an **observation**, its record the strict `xfail`. **Retires** when `test_port_gap_voltage_padding.py`'s padding route is deleted (an `OPS` clean-up) or a consumer of the estimator reappears and re-opens it as a defect. No band moves; `attempt/PORT-1-R1-20261025T171120Z`'s threading is not needed. |

### 🟡 OPEN 2026-09-07 (`WF-6` step 4, 15:00 implementer slot) — the **unloaded** F-small birdcage's FEM `|B₁⁺|` reads **26.62% below** the filament closed form at the centre plane, with every convention verified: the FEM's end rings contribute **10.07%** of the legs' centre field where Kirchhoff's filament rings contribute exactly **50%**

| | |
| --- | --- |
| **Where this fires** | `tests/validation/test_birdcage_b1_plus_closed_form.py::test_the_fem_b1_plus_matches_the_filament_closed_form[0…10]` — all eleven points. **Parked, not on `main`**: branch `attempt/WF-6-step4-20260907T205600Z` (`499c527`). `main` carries only this entry and the three logs. |
| **Verified at** | `b328185` + the parked worktree. `20260907T200630Z_WF-6.log` — **11 failed / 16 passed in 90.55 s** at `-n 4` complex, elapsed 93 s, `Status: 1` (`:1918–1938` the table, `:2499–2520` the asserts, `:2573–2577` the footer). Rule-(c) control: the additive `phantom_material` keyword is a no-op on the gate — `20260907T200858Z_WF-6.log` **5 passed / 44.01 s**, `Status: 0`, elapsed 46 s (`:84, :96–97`). |
| **Literal symptom** (`…200630Z:1924`) | `[ 0]            centre  FEM 8.097478100e-08 T  closed form 1.103500413e-07 T  dev  26.6201%`; the eleven points span **16.7641% … 27.3560%** against the pre-registered `CLOSED_FORM_BAND = 5.0e-2`, which was **not widened**. |
| **Cause — the ring current, not diagnosed further** | The three conventions the 10:30 review pre-paid are all *verified in the same run*, so the item's exit-1 (convention error, > 50%) does not apply: the fixture's leg azimuths are 360 / 90 / 180 / 270 deg, one common offset of **0.000000 deg** with **0.000e+00 deg** spread over the four legs (`:1919–1920`); every sheet's `drive_direction` is asserted `(0, 0, 1)`, which is the closed form's `+ẑ` leg convention; the mode-2 negative control (β) is **green** at **7.743627e-03** of the mode-1 centre against 5.00e-02 (`:1938`), so the FEM does track the closed form's exact zero. The discriminator is control (α): the **legs-only** closed form reads **7.356669419e-08 T**, i.e. the filament's Kirchhoff rings add exactly **50%** (step 4a's `R²/ρ² = 0.5`, `20260907T123926Z_WF-6.log:51`) while the FEM's coil adds only **10.07%** — FEM/legs-only **1.1007** against the filament's **1.500**. The four superposed leg currents are equal to 4 digits (1.820730e-02 / 1.820756e-02 / 1.820518e-02 / 1.820726e-02 A, `:1921`), so the *leg* current the sheets report is not the issue; the ring current the FEM's coil actually carries is ≈5× below `cumsum(I) − mean(cumsum(I))` on those same currents. Not diagnosed further in-slot. |
| **Consequence** | **No absolute `\|B₁⁺\|` claim exists.** §2's B₁⁺ row does **not** gain the "closed-form-gated" clause the item authorized, and `WF-6` stays 🟡 with its symmetry-only gates untouched — every one of them is invariant to a uniform 27% scale, which is exactly why this comparison was scoped. |
| **Resolves with** | A review scoping the ring-current measurement: read the FEM's actual end-ring current (a surface integral of `J` over a ring cross-section, the `WF-6` step 3g integral machinery) and compare it to the Kirchhoff `cumsum` on the measured leg currents. The two candidate mechanisms it separates are (a) the leg current is not uniform along `z` — the sheet reads it at the gap, `z = 0` — and (b) the return path is not the rings at all, the σ = 800 S/m conductor's ≈36 Ω leg resistance being comparable to the 50 Ω port so that displacement current to the PEC box competes with the ring path. |
| **Ruled 2026-09-07 18:00 review — the comparand, twice, not the coil** | (1) The generator places the rings at `z = ±leg_spacing/2 = ±55 mm` (`io/mesh.py:3772, 3821`) and runs the legs to `±70 mm` with dead-end stubs beyond (`:3858`); the item passed `coil_length=0.14`, so the filament's rings sat at ±70 mm on legs 27% longer than the fixture's current path — the free-space form's 1.50 L is 1.414 L at the true planes. (2) The air is closed by a **PEC** box (`:4092–4101`) at ≈ ±0.11 / ±0.11 / ±0.10 m, 4.5 cm from the rings and 4 cm from the legs; `B_n = 0` images at first order cut the rings to ≈ 0.30 L and the legs to ≈ 0.74 L, total ≈ 1.04 L against the measured 1.1007 L (`…200630Z:1937`), the corner images adding the rest. The "rings add 10%" reading compared a boxed FEM with an unboxed, mis-positioned legs-only form and says nothing about the ring current. Disposition (rule (f)): the comparand is re-registered as the same filament form at `coil_length = LEG_SPACING` summed over the box's image lattice, `CLOSED_FORM_BAND = 5.0e-2` unmoved and asserted against it — **`WF-6` step 4b (§9 item 1)**; the ring/leg currents are still measured, as a Kirchhoff check for the ≲ 6% displacement-to-wall term the band's margin must absorb — **step 4c (§9 item 4)**, printed only. Entry retires with 4b's landing commit. |
| **Step 4b executed 2026-09-07 19:30 slot — the ruling is confirmed by measurement and the miss collapses from 26.62% to 2.33% at the centre, but the pre-registered `r = 0.5R` exit fires: ten of eleven points are inside 5%, the eleventh is 7.09%** | `20260908T004020Z_WF-6.log`, **1 failed / 15 passed in 68.81 s** at `-n 4` complex, elapsed **71 s**, `Status: 1` (`:1881–1909` the table, `:1998` the footer); the same module with `tests/environment` and no `-s` at `20260908T003808Z_WF-6.log`, **1 failed / 26 passed in 100.06 s**, elapsed **102 s**, `Status: 1` (`:176–177` the single assert, `:473`). **The re-registered comparand is right and the review's arithmetic lands on the nose:** centre comparands in units of `L` = the free-space legs-only form at `COIL_LENGTH` (7.356669419e-08 T) read **1.5000 L** free space at `COIL_LENGTH` (review 1.500), **1.4140 L** free space at `LEG_SPACING` (review 1.414), **0.7321 L** boxed legs-only (review ≈ 0.74) and **1.0756 L** boxed legs + rings against the FEM's **1.1007 L** (review ≈ 1.04 + corners) — `:1902–1907`. **The eleven points** (`:1889–1899`): centre **2.3345%**, then `+x̂` 2.5175 / 1.3374 / 2.3495 / 3.8483 / 3.5983%, `+ŷ` 1.6429 / 0.1095 / 1.7611 / 4.4223 / **7.0875%**, median 2.3495%, worst point [10] = `(0, 0.035, 0)` m = `r = 0.5R` along `+ŷ`. Control (α′) is **green and asserted**: the step-4 free-space comparand still misses the centre by **26.6201%** (`:1901`), reproducing `…200630Z:1924` to the digit, so the gate discriminates. Control (β) green at **7.743627e-03** (`:1909`). **Two readings the item did not predict, both worth a review's attention.** (a) The PEC box measured off the mesh is **±0.120 / ±0.120 / ±0.100 m**, not the predicted ±0.11 / ±0.11 / ±0.10 (`:1882–1883`) — the review's radial arithmetic was 1 cm short, which is why the first-order 1.04 L under-reads the 1.0756 L the exact lattice gives. (b) **The lattice truncation drift is 4.642e-02 at `N = 3` and 3.297e-02 at `N = 4`** against the predicted ≤ 1e-2 (`:1908`) — i.e. the comparand's own error bar is the same size as `CLOSED_FORM_BAND`, which the item's clause "if it exceeds the band's margin at `N = 4`, report it as the finding" pre-registers as a finding. A 7.09% point against a comparand carrying a ≈3–5% truncation bar is **not** evidence that the FEM is wrong. (c) The FEM's own `+x̂` / `+ŷ` pair at `r = 0.5R` differs by **3.4%** (9.805561792e-08 vs 1.013569652e-07 T), a C4 asymmetry of the same order as the miss. **Anchor (ii), the numpy identities on the new function** (`20260908T003713Z_WF-6.log`, `-n 1`, **1 failed / 9 passed in 5.44 s**, elapsed 7 s): (ii-a) the `N = 0` term equals `birdcage_filament_field` bit for bit — **green**, `np.array_equal`; (ii-b) `P_11 I` equals the currents rotated two legs and each mirror is an involution with zero sum — **green**; (ii-c) the wall-normal identity **converges hard but lands at 1.077e-02 against the pre-registered 1e-02** — `|B_n|/|B_free|` at the six wall centres falls from max **9.990e-02** at `N = 2` to max **1.077e-02** at `N = 3` (`:60–61`), a factor 9.3, which is what a *correct* image map does; the residual is the unpaired outermost shell of the cube truncation and is the same truncation the drift in (b) measures. The band was **not** widened. One in-slot definition correction, disclosed: the item wrote the ratio as `|B_n|/|B|`, but `|B_total|` is the cancelled field — at the `+ŷ` wall of a mode-1 drive it is purely normal by symmetry, so the literal ratio reads **1.000e+00** at both `N = 2` and `N = 3` whatever the lattice does (`20260908T003604Z_WF-6.log:60–61`, the first run). The denominator is therefore the source coil's own `N = 0` field at the same point, and the drive was changed to `(1, 2, −3, 0)` so no wall reads zero by symmetry rather than by construction. **Nothing landed on `main`.** Parked on `attempt/WF-6-step4b-20260908T004458Z`; §9 item 1 marked 🚫 in the same commit (rule (d)). **Resolves with** a review ruling on the *comparand's* truncation, which is now the binding uncertainty and not the coil: either a lattice acceleration / larger `N` with the drift asserted below the band, or a re-registration of the gate as "band + measured drift". `CLOSED_FORM_BAND` stays 5.0e-2, §2's B₁⁺ row gains nothing, `WF-6` stays 🟡. |
| **Ruled 2026-09-08 03:00 review — acceleration, not a wider band** | The plain cube-shell drift falls as **1/N** (5.249e-02 / 4.642e-02 / 3.297e-02 at `N = 2, 3, 4`, ratios 0.88 / 0.71 against 2/3 / 3/4) — the signature of a dipole lattice whose shell sums alternate in sign with `N`, a Leibniz series whose partial sums straddle the limit. The mean of consecutive partial sums `S̄_N = (S_N + S_{N−1})/2` cancels the `1/N` term and converges as `1/N²`; predicted drift of `S̄` at `N = 6` below 3e-3, against ≈ 1.4e-2 for the plain sum at the same order. **`WF-6` step 4d (§9 item 1)** gates the eleven points against `S̄_6` at the unmoved 5% with the comparand's own drift asserted ≤ 1e-2 *first* (stop there if it fails — then cube-shell averaging is not the acceleration and a point-dipole tail is the next route), the wall identity re-asserted on `S̄_6`, and the FEM's C4 spread at the four rotated copies of each radial point printed (the 3.4% at `r = 0.5R` is the mesh-side candidate for the last point; `GEO-28`, §9 item 5, censuses the mesh's quadrants). "Band + measured drift" is **not** adopted — it would license a 3–5% comparand bar against a 5% band. Entry retires with 4d's landing commit. |
| **Step 4d executed 2026-09-08 04:30 slot — the ruling's premise is confirmed at the *interior* points and anchor (i) would pass at 2.242e-03, but anchor (ii) FAILS: shell averaging **destroys** the PEC wall identity, because the wall residual is not an alternating tail — it has a hard even/odd parity in `N` and averaging pairs a good odd order with a bad even one** | Nothing landed on `main` beyond this record and the two logs; the code is parked on `attempt/WF-6-step4d-20260908T094500Z`. **(1) The shell terms alternate, as ruled** — `20260908T093332Z_WF-6.log:40–46`, the signed `(S_N − S_{N−1})·ŷ` at the box centre in units of the free-space `S_0` centre field, mode-1 unit drive, box ±0.120 / ±0.120 / ±0.100 m: **+1.737937e-01, −2.954646e-02, +2.496714e-02, −1.833887e-02, +1.468487e-02, −1.222862e-02** for `N = 1 … 6`. Strict alternation with magnitudes falling as ≈ 1/N; the 03:00 ruling's premise is measurement now, not inference. **(2) At the 22 interior sample points the acceleration works and works well** (`:38–39`): plain drift `max|S_N − S_{N−1}|/|S_N|` = 3.259e-01 / 5.249e-02 / 4.642e-02 / 3.297e-02 / 2.712e-02 / **2.208e-02** at `N = 1 … 6` (reproducing 4b's 4.642e-02 and 3.297e-02 to the digit), against the shell-averaged `max|S̄_N − S̄_{N−1}|/|S̄_N|` = 1.316e-01 / 4.160e-03 / 6.058e-03 / 3.329e-03 / **2.242e-03** at `N = 2 … 6` — a **9.8×** reduction at `N = 6`, inside the item's ≤ 1e-2 anchor (i) and inside its ≤ 3e-3 prediction. **(3) Anchor (ii) is red at 3.292e-02 against the unmoved 1e-02** (`20260908T093523Z_WF-6.log:78`, `1 failed / 10 passed in 10.08 s` at `-n 1` complex, elapsed 12 s, `Status: 1`). `|B_n|/|B_(N=0)|` at the six wall centres, drive `(1, 2, −3, 0)`, box ±0.11 / ±0.11 / ±0.10: the **plain** sums read max **7.508e-02 / 9.990e-02 / 1.077e-02 / 7.082e-02 / 6.046e-03 / 5.980e-02** at `N = 1 … 6` (`:66–71`) — odd orders converge (7.5e-02 → 1.077e-02 → **6.046e-03**, the last *inside* the band), even orders barely move (9.99e-02 → 7.08e-02 → 5.98e-02). The **shell-averaged** sums therefore read max 4.625e-01 / 1.241e-02 / 4.457e-02 / 3.002e-02 / 3.843e-02 / **3.292e-02** (`:72–77`): `S̄_6` is essentially half of `S_6`'s residual, three times *worse* than the plain `S_3` the step-4b construction used. The band was **not** widened and the module stops there; anchors (iii) and (iv) were not measured, so no FEM window was spent. The three supporting identities are green in the same run: (ii-a) `N = 0` equals `birdcage_filament_field` bit for bit, (ii-a′) `ladder[m] == birdcage_filament_field_in_pec_box(image_order=m)` **bit for bit** at `m = 0, 1, 2` and `S̄` 1-based, (ii-b) `P_11` is the two-leg rotation. **Cause — measured, one mechanism named, not proved** | The interior field and the wall-normal field are *different functionals of the same truncation*. Interior: the outermost cube shell contributes a sign-alternating, ≈ 1/N term, so `(S_N + S_{N−1})/2` cancels it — confirmed. On a wall at `x = +X_b`, the identity `B·n̂ = 0` comes from the pairwise cancellation of images `i` and `1 − i` about that wall; a cube truncated symmetrically at `|i| ≤ N` breaks that pairing at the ends, and the residual it leaves flips character with the parity of `N` rather than its sign. Averaging an odd order against an even one therefore inherits the even order's residual at half strength instead of cancelling anything. **Not diagnosed further in-slot.** | **Consequence** | Step 4 does not land; §2's B₁⁺ row gains nothing and `WF-6` stays 🟡; `CLOSED_FORM_BAND` is unmoved at 5.0e-2 and `IMAGE_ORDER` is unmoved at 3 on `main`. §9 item 1 marked 🚫 in this commit (rule (d)). **Resolves with** a review ruling on a comparand whose *both* diagnostics are green at once. The measurement above makes one candidate cheap and dominant: **the odd-order cube sums alone** — `S_5` is 6.046e-03 on the wall (inside 1e-2) and the odd-order pair drift `|S_5 − S_3|/|S_5|` is computable from the same ladder; an "odd-order shell average" `(S_{N} + S_{N−2})/2` would test the same idea one parity at a time. The alternative the 03:00 review already named — a point-dipole tail for `N > 3` — is untouched by this result. Whatever is chosen must be asserted on **both** the interior drift and the wall identity, since this slot's finding is precisely that the two disagree about which truncation is good. |
| **Ruled 2026-09-08 10:30 review — the truncation's sign is now measured and it widens the miss: no comparand fix lands step 4** | From the 4d shell terms, `S_∞ − S_3 ≈ t_4 + t_5 + t_6/2 = −0.00977` on `S_3 = 1.16929` (units of the free-space `S_0` centre field) — the converged lattice is **0.84% below** the `N = 3` sum 4b gated against, so under any converged comparand the centre miss is ≈ 3.2% and the `+ŷ` `r = 0.4R` / `0.5R` points ≈ 5.3% / 7.9% (up to ≈ 6.1% / 8.8% where the off-centre shell terms are twice the centre's). Step 4c's measured currents cannot absorb it (26 / 41 mm stations +0.3% / +0.2% above the terminal currents; ring arcs 1.0–1.3% *below* Kirchhoff). The comparand is finished with the **odd-order sum `S_11`** (wall identity green at odd `N` by measurement — 6.046e-03 at `N = 5`; Leibniz bound `\|t_12\|` ≈ 0.5% at the centre, ≈ 1.1% worst); odd-parity averaging is rejected (both odd sums lie on the same side of the limit). **`WF-6` step 4e (§9 item 1)** asserts the wall identity and the Leibniz bound on `S_11`, then measures the eleven points against it with anchor (iii) *predicted red* at [9] and [10]; on that outcome the `src/` comparand and the unit identities land on `main`, the gate module parks with its table, and this entry gains the table. The FEM side is then the open question — a degree-1 solve on a 0.015 m air mesh never `h`-refined on this fixture — priced by `GEO-29` (§9 item 5, the global-`resolution` ladder) and solved by step 4f (the 18:00 review's to scope). `CLOSED_FORM_BAND` unmoved; the eleven points are never pruned. Entry retires only with a landing that gates all eleven. |
| **Step 4e executed 2026-09-08 12:00 slot — anchor (i) is RED at the first stop, and the comparand route itself is now refuted: the cube-truncated image lattice's PEC wall-normal residual does not converge to zero in *either* parity, it converges to ≈ 3.4e-02** | `20260908T170345Z_WF-6.log`, **1 failed / 21 passed in 128.05 s** at `-n 1` complex with `tests/environment`, elapsed **130 s**, `Status: 1`; no FEM window was spent and nothing landed on `main` beyond this record and the log (code parked on `attempt/WF-6-step4e-20260908T170345Z`). Max over the six wall centres of `\|B_n\|/\|B_(N=0)\|`, drive `(1, 2, −3, 0)`, one `N = 12` ladder call (56 953 filament evaluations, 113.9 s, `:85–99`): **7.508e-02, 9.990e-02, 1.077e-02, 7.082e-02, 6.046e-03, 5.980e-02, 1.384e-02, 5.399e-02, 1.834e-02, 5.040e-02, 2.127e-02, 4.796e-02** at `N = 1 … 12` (`N ≤ 6` reproduce 4d's `…093523Z:66–71` to the digit). The **odd** sub-sequence has a *minimum at `N = 5`* and rises: 6.046e-03 → 1.384e-02 → 1.834e-02 → **2.127e-02**, increments +7.79e-03 / +4.50e-03 / +2.93e-03; the **even** sub-sequence falls with decrements −5.81e-03 / −3.59e-03 / −2.44e-03. Both parities approach a **common non-zero limit ≈ 3.2–3.5e-02**, which is where the shell averages sit and stay (3.682e-02 / 3.392e-02 / 3.617e-02 / 3.437e-02 / 3.584e-02 / 3.462e-02 at `N = 7 … 12`, `:98`). **`S_5`'s 6.046e-03 was a crossing, not a convergence**, and the 10:30 review's "the odd sequence extrapolates below 3e-3" is refuted. **Mechanism (named, not proved):** the sum is conditionally convergent, so its value depends on the summation order; the wall identity `B·n̂ = 0` needs image `i` paired with image `1 − i` about that wall, and a symmetric cube truncated at `\|i\| ≤ N` never pairs the outermost shell — that unpaired shell subtends a *fixed* solid angle at the wall for every `N`, so its contribution tends to a constant, exactly the `1/N²` two-sided approach to a constant measured above. **Consequence:** no cube order in any parity satisfies the identity inside the unmoved 1e-02 except the accidental crossing at `N = 5`, and gating on a crossing is fitting. `CLOSED_FORM_BAND` 5.0e-2, `WALL_NORMAL_BAND` 1.0e-2 and `IMAGE_ORDER` 3 all unmoved on `main`; §9 item 1 marked 🚫 (rule (d)); `WF-6` stays 🟡. **Resolves with** a review scoping a summation whose order is not a free choice — the point-dipole tail (near images exact, far lattice by its multipole limit, absolutely convergent) or an Ewald split — or, cheaper, abandoning the closed-form comparand for this fixture and gating `\|B₁⁺\|` by `h`-convergence instead (`GEO-29` prices the ladder). Not a band widening and not `N = 5`. |
| **Ruled 2026-09-09 03:00 review — the closed-form comparand route is finished for this fixture, and step 4 stops chasing it** | 4e did not merely miss a band; it **refuted the route**. The cube-truncated image lattice's wall-normal residual has a minimum at `N = 5` and *rises* thereafter in the odd parity while the even parity falls, both approaching a common non-zero ≈ 3.4e-02 — the sum is conditionally convergent and the unpaired outermost shell subtends a fixed solid angle at each wall, so **no cube order in any parity satisfies the identity** and `S_5`'s 6.046e-03 was a crossing. Neither of the two remaining routes (point-dipole tail, Ewald split) is worth a slot: the 10:30 review's own arithmetic already showed the converged comparand sits 0.84% *below* the `N = 3` sum, i.e. every miss grows, and 4c's measured currents cannot absorb it (26 / 41 mm stations +0.3% / +0.2%; ring arcs 1.0–1.3% below Kirchhoff). The truncation is not the binding error. **Disposition: the analytic comparand is set aside — not widened, not re-registered, not used.** `CLOSED_FORM_BAND` stays 5.0e-2, `WALL_NORMAL_BAND` stays 1.0e-2, `IMAGE_ORDER` stays 3 on `main`, and the eleven points are never pruned. What remains is exactly the residue the 10:30 ruling named, and two probes have now narrowed it to one suspect: `GEO-28` **excluded the mesh's symmetry** (every quadrant mass spread ≲ 0.1% against a 3.4% field effect, `20260908T183317Z_GEO-28.log:1800–1830`) and `GEO-29` **measured its resolution** — interior mean circumradius 2.19e-02 m at ×1, **1.46× the nominal 0.015 m and larger than the 1.4e-02 m shell the gate samples in**, falling to 1.71e-02 / 1.37e-02 on the next two rungs (`20260909T003221Z_GEO-29.log:7038–7042`). **`WF-6` step 4f (§9 item 3)** therefore measures the FEM against **itself and a symmetry identity** rather than a closed form — the field's C4 four-copy spread per radius across the 0.015 / 0.012 / 0.0095 rungs, with `GEO-28`'s ≈ 0.1% mesh spread as the floor, the ×1 rung's recorded 3.4% as the asserted reproduction control, and `POST-6` step 1's cw drive (95.1975%) as the asserted negative control at ≥ 10× separation (ceiling 28×). A symmetry identity is §4-compliant without a closed form. Under rule (e) the *fall itself* is **predicted, not asserted** — no prior `h`-ladder on this fixture backs it — so 4f closes step 4 only if the fall lands, and a spread that does not fall is the more informative outcome: it moves the 3.4% off discretisation and onto the degree-1 formulation or the CG1 `curl E` estimator. §2's B₁⁺ clause gains nothing either way and `WF-6` stays 🟡; the 2026-09-13 weekly holds the dated Phase-5 exit decision (gate, or a measured convergence statement) against 4f's result. Entry retires only with a landing that gates all eleven points, or with the weekly converting the target. |
| **Step 4f executed 2026-09-09 07:30 slot — the `h`-ladder runs to completion and the answer is a *partial* yes: the C4 four-copy spread falls **5.2506% → 2.0719% → 1.9514%** and then stalls, an order above `GEO-28`'s ≈ 0.1% mesh floor** | `20260909T123716Z_WF-6.log` (collect-only smoke `20260909T123702Z_WF-6.log`, 40 items / 4 s), `-n 4`, complex build, `timeout -k 30 560`, **4 failed / 32 passed / 4 skipped in 451.24 s**, harness elapsed **453 s**, `Status: 1`. **Parked, not on `main`**: branch `attempt/WF-6-step4f-20260909T124552Z` (`ab2a2cf`); `main` carries this row, the two logs and the §7 annotation. **The ladder** (`…123716Z:2080–2086, 3901–3907, 5789–5795`), `resolution` 0.015 / 0.012 / 0.0095 m at cell counts 116 085 / 149 049 / 197 393, **all three reproduced to ratio 1.000000**: worst-radius (0.5R) C4 four-copy spread of the ccw-quadrature `\|B₁⁺\|` **5.2506% / 2.0719% / 1.9514%** (ratios to ×1 1.0000 / **0.3946** / **0.3717**); the gated C4 covariance falls with it, **3.6159% / 1.6815% / 1.5029%** against the imported, unmoved 5% band (green on all three rungs); the eleven-point drift is median 0.36% / 0.91%, max 2.46% / 1.84%. `GEO-29`'s interior mean `h` over the same rungs is 2.19e-2 / 1.71e-2 / 1.37e-2 m, so **the ×0.012 → ×0.0095 step buys 0.12 pp on a 2 pp spread** — refinement owns ≈ 60% of the effect and then stops, and the residual ≈ 2% is the degree-1 N1curl solve's or the CG1 `curl E` estimator's, not `h`'s. **Literal symptoms, the three reds.** (1) `test_the_x1_rung_reproduces_the_recorded_c4_spread[x1]` — `the x1 rung's worst-radius C4 four-copy spread reads 5.2506% against the recorded 3.3106% — 58.60% relative, outside the 10% reproduction control`. (2) `test_the_cw_drive_spread_dwarfs_the_ccw_spread_on_every_rung[x1]` — `the cw drive's worst-radius C4 spread 50.0268% is only 9.53x the ccw quadrature's 5.2506%`, against a 10× bar; the same test passes at **19.52×** and **58.38×** on the finer rungs. (3) `test_power_accounting_closes_on_every_rung[x0.0095]` — `power accounting misses by 1.853642e-02 of the supplied 6.796871053e-03 W (phantom 0.000000000e+00, conductor 4.378232190e-04, sheets 6.233058154e-03); imported band 1e-02` (P2 1.419812e-02). **Cause of (1) and (2) — diagnosed, and it is the pre-registration, not the fixture:** the ×1 rung reproduces step 4b's eleven-point table **character-for-character** (9.805561792e-08 T at `+x̂`, 1.013569652e-07 T at `+ŷ`, `20260908T004020Z_WF-6.log:1894, 1899`), so the reproduction control is met in the strongest possible sense; the anchor equated a **four**-copy spread with a **two**-copy record, and the `−x̂` copy **1.024222080e-07 T** — a quantity no prior run measured — is the maximum that record never contained. (2) inherits the same mis-specification through its 95.2/3.4 = 28× sizing. **Nothing was re-tuned in-slot**; re-registering the anchor (on the two-copy statistic, or on the four-copy one at its own measured value) is a review's ruling. **Cause of (3) — not diagnosed:** the power residual is 9.796e-03 at ×1 and 8.114e-03 at ×0.012 and **degrades to 1.85e-02 under refinement**; on the same rung the sweep prints `Z` class spreads **1.3475 / 0.4544 / 0.9165%** against `PORT-9`'s 0.5% (ungated in this module). This is the *global*-resolution twin of `ANS-4` step 2a's finding that refining `conductor_resolution` breaks the same fixture's C4 symmetry — two independent refinements, the same class of symptom. **Consequence** | No absolute `\|B₁⁺\|` claim exists and none is created; §2's B₁⁺ clause keeps its wording, `CLOSED_FORM_BAND` (5.0e-2), `POWER_BALANCE_BAND` (1e-2), `C4_COVARIANCE_BAND` (5e-2) and `CELL_COUNT_BAND` (1e-2) are all unmoved on `main` and on the branch, the closed-form comparand path is **deleted** from the gate (ruling (1)) rather than widened, and `WF-6` stays 🟡. The 2026-09-13 weekly's dated Phase-5 exit decision now has a **measured convergence statement** in place of a gate. **Resolves with** a review doing three things: re-registering the ×1 anchor against the four-copy statistic's own measured 5.2506% (or the two-copy 3.3106% it was written for), re-sizing the cw separation bar off that number, and scoping the ×0.0095 power residual — which must be diagnosed before any finer rung of this fixture is trusted by `ANS-4` or `WF-6` alike. |
| **Ruled 2026-09-09 10:30 review, ruling (3) — all three are done, and the third is scoped *out* of the ladder rather than absorbed into it** | (1) **The ×1 anchor is re-registered on the statistic that exists**, the four-copy spread's own measured **5.2506%** inside 10% relative (`20260909T123716Z_WF-6.log:2080–2086`). The fixture is exonerated, not accommodated: the ×1 rung reproduced step 4b's eleven-point table *character-for-character* (`20260908T004020Z_WF-6.log:1894, 1899`), so the reproduction control was met in the strongest available sense and the `−x̂` copy 1.024222080e-07 T is simply a quantity the two-copy record never contained. (2) **The cw bar is re-sized to ≥ 5×, off a ceiling the original never computed.** At ×1 the cw spread is 50.0268% against the ccw 5.2506%, so **9.53× is arithmetically the maximum this rung can show** and the pre-registered 10× bar was *unreachable*, not merely unmet; 4f measured 9.53× and 19.52× on the two rungs kept, so 5× holds with margin and would still hold if the ccw spread halved again. Neither number was ever an imported band or a green record, so nothing is loosened — this is the rubric's own ceiling-first rule applied after the fact. (3) **The ×0.0095 power residual is not re-registered, not widened, and not carried.** `POWER_BALANCE_BAND` stays 1e-2 and the rung is **dropped from `WF-6` step 4g's ladder** because this row itself forbids trusting a finer rung of this fixture until the residual is diagnosed, and the 02:15 weekly's Phase-5 exit clause names a two-rung monotone fall as sufficient on its own terms. Its readings stay on record above and are not deleted. **The residual's suspected mechanism is now a measurement, not a hypothesis:** the P1/P2 split at ×0.0095 is a *symmetry* break on the global-resolution axis, the twin of `ANS-4` step 2a's break on the conductor axis, and `GEO-30` (§7, §9 item 1) censuses both with no solve. If `GEO-30` finds the mesh symmetric at ×0.0095, the rung returns to the ladder — that is the next review's call, not step 4g's. `CLOSED_FORM_BAND` (5e-2), `C4_COVARIANCE_BAND` (5e-2) and `CELL_COUNT_BAND` (1e-2) are unmoved, the eleven points are never pruned, and the deleted closed-form comparand is not re-introduced. Entry stays open until the residual is diagnosed. |
| **Ruled 2026-09-09 18:00 review, ruling (2) — the ×0.0095 rung stays OUT, and the residual ≈ 2% is now attributed away from `h` by two independent measurements** | Step 4g's item said the dropped rung comes back "if item 1 finds the mesh symmetric at ×0.0095", and `GEO-30` did **not** find that. It found the mesh symmetric in *mass* (coil-volume spread 0.17% at that rung) and in gap-sheet *area* (6.05e-16) — but the same rung's gap-sheet **facet counts** read **70/76/70/76**, C2 rather than C4, which is exactly the port-pair split this row's own power residual shows (P1 1.853642e-02 / P2 1.419812e-02 against five-figure P1/P2 agreement on both coarser rungs). The rung's symmetry defect is therefore *reproduced by a second, independent, no-solve measurement on a different quantity*, not excluded. **The rung is not rehabilitated and is not re-introduced to the ladder; no band is widened and none of 4f's readings are deleted** — they stand in the §7 `WF-6` row. **What is settled, and it is the substantive result of the interval:** the two-rung ladder landed green on `main` (`49432f9`, `20260909T200431Z_WF-6.log`, Status 0, 204 s, 28 passed / 2 skipped) with the four-copy C4 spread falling **5.2506% → 2.0719%** (ratio 0.3946) at unmoved bands and a cw-vs-ccw separation of 9.53× / 19.52× against the re-sized 5× bar, so refinement owns ≈ 60% of the effect and **the surviving ≈ 2% is not `h`'s** — it belongs to the degree-1 N1curl solve, the CG1 `curl E` estimator, or the sheet reconstruction `GEO-31` now measures. That is the measured convergence statement the 2026-09-13 weekly's Phase-5 exit decision was waiting for, and it is handed to that review rather than converted into a gate here. **Disposition: `GEO-31`** (§9 item 1) is the shared next measurement for this row and the `ANS-4` step-2a row above — a single probe, as `GEO-30` was. `WF-6` stays 🟡, §2's B₁⁺ clause does not move, and this entry stays open until the residual is diagnosed. |
| **Step 4g executed 2026-09-09 15:00 slot — the two re-registered anchors and the re-sized negative control are GREEN on both kept rungs, and the landing is on `main` rather than a branch** | `20260909T200431Z_WF-6.log` (collect-only smoke `20260909T200247Z_WF-6.log`, 12 items / 5 s), `-n 4`, complex build, `timeout -k 30 500`, **`28 passed, 2 skipped` in 202.58 s** (`:3917`), `Status: 0`, harness elapsed **204 s** (`:4116-4117`). Every 4f number on the two kept rungs reproduces to the printed digit (`:1982-1994, 3803-3815, 3828`): spread **5.2506% -> 2.0719%** (ratio 0.3946), cells 116 085 / 149 049 at ratio 1.000000, covariance 3.6159% / 1.6815%, power residuals 9.795836e-03 / 9.796294e-03 and 8.113516e-03 / 8.111819e-03, cw separations **9.53x / 19.52x** against the new 5x bar, drift median 0.9482% / max 2.9092% (printed). The eleven-point table reproduces step 4b character-for-character (`:2000, 2005`). **This entry does not retire.** It stays open on the clause it was ruled open on — the ×0.0095 power residual is still undiagnosed and its rung is still out of the ladder — and nothing here rehabilitates it. No band moved: `CLOSED_FORM_BAND` (5e-2), `POWER_BALANCE_BAND` (1e-2), `C4_COVARIANCE_BAND` (5e-2), `CELL_COUNT_BAND` (1e-2) all imported and unmoved, eleven points never pruned, closed-form comparand still deleted. **One correction to the item's premise, recorded not absorbed:** `main` never carried 4f's additive `resolution` / `phantom_material` keywords, so the slot's first run died on `TypeError: build_four_port_sweep() got an unexpected keyword argument 'phantom_material'` (`20260909T200300Z_WF-6.log:716, 737`); `tests/mesh/test_birdcage_port_sheets.py` and `tests/validation/test_port_birdcage_four_port.py` were taken from the same branch, both `None`-defaulted and bit-for-bit inert for every gate. `WF-6` stays 🟡; no closed-form B₁⁺ gate exists. |

### 🔴 OPEN 2026-08-25, re-headed 2026-08-26 (`GEO-21` step 2) — `birdcage_port_domain` **cannot mesh a coarse conductor sizing on the 0.11 image**: `conductor_resolution=None` and everything coarser than ≈4.8 mm abort in gmsh with "Invalid boundary mesh (overlapping facets)"

> **Retire-when gains (2026-10-03, §9 item 186):** "the item-186 control goes
> red" — `test_ungraded_build_still_hits_the_documented_floor`
> (`tests/mesh/test_birdcage_port_tags.py`) asserts the ungraded build still
> raises on every rank; green at `-n 2` / `-n 1` on 2026-10-03.

> **RULED 2026-08-30, weekly planning review — the floor stays a documented
> limitation, deliberately uncommissioned; it becomes a stated trap, not a
> chunk.** `GEO-23` has since classified the whole "overlapping facets"
> family as geometry-deterministic and closed with "land no fix" by
> commission; every gate that reads this generator is green on a graded
> control (`h_c = 4.8e-3`, ruling (b) of 08-26), and the production
> fixtures — F-small at 116 085 cells and the 16-leg rungs — all mesh with
> `conductor_resolution = 0.4 × ring_minor_radius = 1.6 mm`, three times
> finer than the floor. A root-cause hunt inside gmsh's boolean fragment
> buys nothing the mission needs. **Where it bites next, stated now:** the
> F-human cost probe (`GEO-25`, §7) scales `ring_radius` 0.07 → 0.15 m and
> may be tempted to coarsen the conductor to hold the cell count — it may
> not go coarser than 4.8 mm on the conductors without first re-measuring
> this floor at that scale, and a probe rung that aborts with this string
> is to be recorded against this entry, not diagnosed in-slot. Re-opens as
> a chunk only if that probe finds no affordable rung *above* the floor.

> **✅ The gate-red portion of this entry RETIRED 2026-08-26** (`GEO-21` step 2,
> 04:30 implementer slot), exactly as the retire-when below specifies.
> `test_graded_conductor_sizing_recovers_the_cad_mass` and `mesh:3` are **green
> on `main`**: the negative control moved `None` → `BASELINE_CONTROL_RESOLUTION`
> = 4.8e-3 (the 03:00 review's ruling (b)), version-tagged with the six-rung
> probe table in-comment, and the demoted claim — **fine vs coarse grading**,
> no longer "grading required" — stated in the module docstring and the `mesh:3`
> guide. Gate `20260826T093202Z_GEO-21-step2-gate.log`, `1 passed in 41.11s`,
> Status 0, 43 s, `-n 2`: control **0.846150** at 33 185 cells, 3.2e-3
> **0.916742** at 47 975, graded **0.966977** at 98 666 — every step-1 probe
> figure reproduced exactly, now through the gate's own assertions. Consumer
> check `20260826T093403Z_GEO-21-step2-mesh3.log`, Status 0, 29 s, separation
> 0.120826. `CAD_MASS_GATE`, the `- 0.05` separation guard and `CONDUCTOR_RUNGS`
> are all unmoved; nothing was loosened.
>
> **What stays open is the generator limitation in the heading**, which
> `GEO-21` step 1 measured to be *wider* than "the ungraded path": 9.6e-3 fails
> the same way at a fourth distinct surface pair, so `conductor_resolution=None`
> is the coarsest point of a continuum whose coarse end stopped meshing at the
> 0.11 merge, not a special broken path. Hardening `birdcage_port_domain`
> against it would have to cover coarse *graded* sizings too. Still
> **deliberately not commissioned** — no production path uses a coarse
> conductor sizing now that the control sits at 4.8e-3, and that decision is
> recorded here rather than silently made. **Retire-when:** a chunk that
> hardens the generator, or a gmsh/image change that makes the coarse end mesh
> again — re-measure the ladder before retiring, do not infer it.
>
> **Verified at** the `GEO-21` step 2 landing commit.
>
> ---
>
> *Everything below is the history of the retired gate red, kept because the
> measurements in it are what the disposition rests on.*
>
> **Where this fired.**
> `tests/mesh/test_birdcage_conductor_sizing.py::test_graded_conductor_sizing_recovers_the_cad_mass`,
> and through it `examples/meshing/03_birdcage_graded_conductors.py` (`mesh:3`),
> which imports `CONDUCTOR_RUNGS` / `CAD_MASS_GATE` / `_check_geo9_identities`
> from that module per `ANS-1`. Both build `baseline = _mesh(conductor_resolution=None)`
> **first**, so both abort before the graded rung — the rung that actually
> carries the gate — ever runs.
>
> **Literal symptom** (`20260825T213821Z_EX-30-mesh-birdcage-gate-probe.log`,
> `-n 2`, real, **`1 failed in 2.51s`**, `Status: 1`; the example's own abort is
> `20260825T213142Z_EX-30-mesh-run-1to5.log`, `Status: 1`):
>
> ```
>     baseline = _mesh(conductor_resolution=None)
> E   RuntimeError: birdcage_port_domain geometry generation failed on rank 0
> E   Exception: Invalid boundary mesh (overlapping facets) on surface 59 surface 79
> ```
>
> **Localised, not inferred** — `tests/mesh/probe_birdcage_conductor_resolution.py`
> (new with this entry; measurement only, asserts nothing, imported by nothing),
> `20260825T213926Z_EX-30-mesh-birdcage-resolution-probe.log`, `Status: 0`,
> 39 s, `-n 1`:
>
> ```
> Leg A -- the fixture's global resolution (0.015), conductor sizing swept:
>   h_c = None    (baseline) FAIL  Invalid boundary mesh (overlapping facets) on surface 59 surface 79 (1.8 s)
>   h_c = 3.2000e-03         OK       47975 cells (10.4 s)
>   h_c = 1.6000e-03         OK       98666 cells (20.7 s)
>
> Leg B -- the baseline's conductor sizing (h_c = None), global resolution stepped finer:
>   h = 0.0150              FAIL  ... on surface 59 surface 79 (1.7 s)
>   h = 0.0130              FAIL  ... on surface 48 surface 48 (1.5 s)
>   h = 0.0110              FAIL  ... on surface 65 surface 65 (1.3 s)
> ```
>
> **It is the conductor sizing, not the resolution.** Both `GEO-15` rungs mesh
> at the *same* global 0.015 that the baseline fails at, and refining the global
> size does not walk out of the failure — three finer steps fail on three
> *different* surface pairs. This is the **opposite** reading from the
> `straight_wire_domain` entry below, where `resolution` alone explained
> everything and every geometry failed at h = 0.01; the two findings are the
> same *family* (0.11 gmsh meeting a parameter set no green gate exercises) but
> not the same axis, and a single ruling will not cover both.
>
> **Consequence.** The `GEO-15` graded-conductor CAD-mass gate — the one
> `EX-21` and `mesh:3` rest on — has been **non-executing on `main` since the
> 0.11 merge**, unobserved, in the same class as `OPS-24`'s cavity gate and leg
> (root)'s `MAG-13` convergence gate. What it *would* have measured is now
> partly known: the graded rung meshes at **98 666 cells**, against the
> 2026-08-16 record of 98 474 (0.7.2) — but that number comes from this probe,
> which does not run the gate's assertions, so it is a bracket, not a re-record.
>
> **Not diagnosed further, and deliberately not fixed.** Whether the baseline
> control moves to a sizing that meshes, `birdcage_port_domain` is hardened
> against the ungraded path, or the inverted control is re-chosen is a
> `GEO-15`/`EX-21` ruling, not `EX-30`'s — this chunk files reds, it does not
> repair them, and no band was touched.
>
> **Verified at** `9b679d8`.
>
> **RULED 2026-08-25, 18:00 review — the disposition is commissioned as
> `GEO-21` (§7), following the `MAG-13`→`MAG-19` precedent: `GEO-15`'s ✅ is
> the 0.7.2 close and stands; a new chunk disposes of the 0.11 red.**
> Measure-first, decision rule pre-stated in the `GEO-21` entry: probe
> `h_c = 3.2e-3`'s CAD-mass recovery (10 s, already priced) — if it sits
> clearly below the gate the way `h_c = None`'s 0.7403 did, the baseline
> control moves there version-tagged (old `None` in-comment citing the
> resolution probe log); if it *clears* the gate, the inverted premise has no
> meshable carrier and the finding is reported, never manufactured around.
> The generator finding — `birdcage_port_domain(conductor_resolution=None)`
> cannot mesh on 0.11 at any global resolution tried — **stays in this entry
> even after the gate goes green**: hardening the ungraded path is
> deliberately not commissioned while no production path uses it, and that
> decision is recorded here rather than silently made. **Retire-when
> (narrowed):** the gate-red portion retires with the commit that lands
> `GEO-21` green; the ungraded-path generator limitation then re-heads this
> entry and stays open.
>
> **MEASURED 2026-08-26 (`GEO-21` step 1, 00:00 implementer slot) — the
> ruling's own branch is excluded, and the generator finding gets wider.**
> The red reproduced unchanged first (`20260826T050100Z_GEO-21-step1-red-repro.log`,
> `1 failed in 4.80s`, same surfaces 59/79). Then the number the ruling turned
> on, measured through the gate module's own `_mesh`
> (`20260826T050134Z_GEO-21-step1-cad-mass-probe.log`, `-n 2`, Status 0, 35 s):
>
> ```
>   h_c = 3.2000e-03  cells=  47975  meshed/CAD=0.916742
>   h_c = 1.6000e-03  cells=  98666  meshed/CAD=0.966977
> ```
>
> **0.916742 is neither branch** — not ≤ 0.90, not clearing 0.95 — and it sits
> inside the gate module's pre-registered guard
> `baseline_ratio < CAD_MASS_GATE - 0.05` = 0.90, failing it by 0.016742. So
> moving the baseline control to `h_c = 3.2e-3` would relocate this red to the
> separation guard rather than clear it, and no licence permits loosening a
> guard that says the premise needs re-examining. The graded rung is *not* in
> question: 0.966977 ≥ 0.95, 98 666 cells matching this entry's bracket exactly.
>
> **The generator limitation is broader than "the ungraded path".** A
> coarse-ward graded ladder (`20260826T050319Z_GEO-21-step1-control-ladder.log`,
> `-n 1`, Status 0, 30 s; measurement only, no control adopted) reads
> 3.2e-3 → 0.916742 (width control exact vs `-n 2`), 4.8e-3 → 0.846150,
> 6.4e-3 → 0.767219, and **9.6e-3 → FAIL, "Invalid boundary mesh (overlapping
> facets) on surface 54 surface 86"** — the same failure family as
> `h_c = None`, at a fourth distinct surface pair. `conductor_resolution=None`
> is therefore not a special path that broke: it is the **coarsest point of a
> continuum whose coarse end stopped meshing at the 0.11 merge**. That widens
> the finding this entry keeps open after the gate red retires — hardening
> `birdcage_port_domain` would have to cover coarse graded sizings too, not
> just the `None` default. Still deliberately not commissioned; still no
> production path uses it.
>
> **Gate red stays OPEN**: disposition needs a review ruling between (b) a
> coarse-graded control (4.8e-3 or 6.4e-3, both separating) and (c) retiring
> the baseline comparison — options and their cost in claim-strength are in the
> §7 `GEO-21` entry. Retire-when is unchanged.
>
> **RULED 2026-08-26 03:00 review: option (b), control = `h_c = 4.8e-3`**
> (0.846150 — clears the 0.90 separation guard by 0.0538 with the guard
> unmoved; `6.4e-3` rejected for cliff adjacency, the cliff having moved
> once already at the 0.11 merge; (c) rejected because a graded-side-only
> assertion cannot distinguish broken grading from common-mode drift). The
> demoted claim — **fine-vs-coarse grading**, no longer "grading required" —
> lands in-comment and in the guide with the landing, `GEO-21` step 2,
> §9 item 1. On that landing this gate red retires; the
> generator-continuum finding above re-heads and stays open, still
> deliberately uncommissioned.
>
> **Verified at** `ab55ff1`. **Landed 2026-08-26 — see the retirement note at
> the head of this entry.**

### 🚫 OPEN — the gapped birdcage's **open-limit (1e6 Ω) driven self-impedance `Z₁₁` is not mesh-converged**: it moves ≈40% under a 0.24% cell-count change, while the terminated column moves 1.9e-02 (`GEO-19` step B attempt 2, 2026-08-24)

**Test id:** `tests/validation/test_port_birdcage_termination_probe.py::test_the_open_control_reproduces_leg_c_before_the_knob_turns`
(and, on the same fixture,
`tests/validation/test_port_birdcage_lumped_column.py::test_adjacent_ports_of_the_driven_leg_agree_and_the_opposite_one_does_not`).
**Neither is red on `main`, and neither can be: step B landed 2026-08-25
under ruling (6\*) and both open-limit record assertions are gone with it**
(the Z-column reproduction retired, the degeneracy ordering assertion
retired; both quantities are still solved and printed as diagnostics, with
their digits kept in-comment as mesh-tagged history). `19 passed` twice from
`main` at 116 085 cells — `20260825T003622Z_GEO-19-stepB-port9-run1.log`,
`20260825T003832Z_...-run2.log`. This entry is no longer about a red; it is
the standing record of the **conditioning finding itself**, which is
unmeasured and unfixed. The attempt branch it once pointed at is deleted;
the content is on `main`.

**Literal symptom**, `20260824T183519Z_GEO-19-stepB-port9-measure.log`
(`3 failed, 16 passed` / 117.80 s / Status 1), the open (1e6 Ω) column at
116 085 cells against its record at 116 368:

```
Z_11 +9.201557829e+02-4.718342449e+03j  vs record +7.111692404e+02-3.351665665e+03j
Z_21 +1.390012417e+01-1.872224592e+03j  vs record +1.224919287e+01-1.878346946e+03j
Z_31 +1.322525314e+01-1.872769896e+03j  vs record +1.193721196e+01-1.878700877e+03j
Z_41 +1.465032447e+01-1.872096207e+03j  vs record +1.231338434e+01-1.878312313e+03j
```

`|Z₁₁|` goes 3.42e+03 → 4.81e+03 Ω, a **40.6%** move; the three mutuals move
0.3%; the driven current moves 1.381e-03. On the *terminated* (50 Ω) fixture
the same mesh change moves `Z₁₁` by 1.852e-02 and `Z₄₁` by 5.9e-04.

**Cause — measured, and it is not a defect in step B.** The mesh moved for the
reason ruling (4\*) already adjudicated (the local-frame construction is
exact-onto at the four axis azimuths but not bit-identical to the old one;
gmsh tie-breaking turns ≈5 ulps into 116 368 → 116 085 cells). What is new is
the *sensitivity*: at `Z_p = 1e6 Ω` the port is very nearly open, so `I₁` is a
near-cancellation residual (≈1e-9 A) and `Z₁₁ = V₁/I₁` inherits its
conditioning. The same fixture's degeneracy margin flips in step with it —
leg (c)'s magnitude-only reading 5.0594× → **0.7906×**, and the complex form
6.9398× → 1.5951×, both already below leg (d0)'s 10× floor *before* step B.
The terminated fixture shows the opposite behaviour and gets **better**:
margin 253.2002× → 2256.9707×, class separation 150.3584× → 166.6766×, every
intra-class spread down (0.0617/0.0359/0.0237% → 0.0553/0.0353/0.0214%).

**Why this is not disposed of by a re-record.** §9 item 2 licensed a
mesh-tagged re-record of the moved records and a pre-registered disposition for
the degeneracy gate. It did not anticipate a 40% move, and pinning `Z₁₁` at a
1e-9 print band on a quantity with no demonstrated mesh stability would record
noise as a fact. The reading this entry offers the review: **the open-limit
column is a diagnostic, not a record-bearing fixture**, and the anti-degeneracy
role it was carrying is already carried, with two decades more margin, by leg
(d0)'s terminated discrimination gate and leg (d)'s 4×4 class separation —
both gated on `main`, both green, both improved by step B.

**Flagged to the weekly review** (§10 Phase 6): if the open-limit column is
retired as a record, `PORT-9` leg (c)'s reproduction anchor needs a replacement
on the terminated fixture, and leg (d1′) should be re-scoped to match.

**Ruled 2026-08-24, 18:00 review — (6\*), option (A): the open-limit column
is retired as a record-bearing fixture; step B lands with the retirement
(§9 item 1), leg (c)'s anchor re-sites on its driven `I₁` + the terminated
fixture, (d1′) re-scoped.** The entry stays OPEN as the record of the
conditioning finding itself. **Retire when:** an h-refinement rung measures
the open column's conditioning, or Phase 6 adjudicates that no open-limit
quantity is record-bearing for the tuning workflow — whichever a review
commissions first. The Phase 6 flag above stands.

**Executed 2026-08-25 (§9 item 1, 19:30 slot).** The retirement is in code
on `main` and the three modules are green twice in-slot. Nothing about the
*finding* changed: `|Z₁₁|` at `Z_p = 1e6 Ω` is still unconverged, still
undiagnosed beyond the near-cancellation argument above, and no h-refinement
rung has been run. What changed is that no record now rests on it. The two
gates that absorbed the anti-degeneracy duty read **2256.9707×** (leg (d0),
floor 10×) and **166.6766×** (leg (d) class separation, floor 10×) on step
B's mesh, both improved. Retire-when is unchanged.

**Verified at:** `cc4ab78` (`main`) + `6c1f54e`
(`attempt/GEO-19-stepB-20260824T183000Z`), 2026-08-24; re-verified on `main`
2026-08-25 with the retirement landed.

### 11. `two_torus_domain` gap-box terminal facet tags include lateral strips below `gap_overhang ≈ 6e-4`

| | |
|---|---|
| **Test** | `tests/validation/test_port_gap_voltage_impedance.py` disc-area band (on `attempt/PORT-1-step3bv-20260806T004500Z`, not on `main`) — and any future test that gates on the `201`/`202` facet areas at small overhang |
| **Symptom** | Measured facet-group area `1.643447371e-04 m²` per port at `gap_overhang = 2e-4` — `1.0241 ×` the exact oblique cut `1.604721580e-04 m²`, i.e. **above** a value an inscribed linear-tet section must sit below |
| **Cause** | Measured, not guessed (`PORT-1` step 3b-v, `20260806T003559Z_PORT-1-step3bv-gate.log`): at overhang 2e-4 the tube protrudes 0.2018 mm through the gap box's `−x` face over `2.821 mm < \|y\| < 3.989 mm` (box `min x` = 1.480000e-02, tube `min x` at `y = half_y` = 1.459821e-02), so the fragment-boundary intersection that defines tags `201`/`202` picks up the arc-end disc pair **plus two lateral strips**. The "gap box contains the arc ends" invariant fails below overhang ≈ 6e-4. 3b-iv's disc-area band was measured at overhang 1e-3, where the tube clears the face by 0.598 mm — it does not transfer to small overhang, and the mesh is not what is wrong. |
| **Verified at** | `main` fixture geometry as of `7747999`; the failing band assertion lives only on the parked branch |
| **Fix options, neither taken in-slot** | Raise `GAP_OVERHANG` back above ≈6e-4 (changes the comparison geometry all 3b measurements share), or make the band overhang-aware (the strip area is computable from the same arithmetic above). A per-geometry decision for whoever next gates on these tags. |
| **Owned by** | `PORT-1` steps 3b-vi/3b-vii note it as a trap (do not gate on the 2xx areas at overhang 2e-4); entry leaves with the commit that restores an exact terminal-area anchor at the geometry it is asserted on |
| **Second measurement, 2026-08-07 (`PORT-1` step 3b-x)** | The strips bias *any* facet average, not just the area. Reading the terminal angle as `arcsin(⟨y⟩/a)` over tags `201`/`202` gives **0.173852206 rad against the exact 0.175335123** — 1.48e-3 short, because the strips sit at `\|y\| < half_y` (`20260807T093604Z_PORT-1-step3bx-gate-n2.log`). The workaround that step took, and the one to reuse: gate the interface's **extreme** reach (`max \|y\|` over the tagged facets' nodes, exact to 5.6e-17 rad), since every strip point lies inside the box while the box face is a plane its nodes sit on exactly. Print the mean beside it — it is this entry's magnitude on the live fixture. |

## Non-test issues

### The container-side `timeout` in the standard harness recipe does not reliably stop an `mpiexec` job, and an overrun can wedge the container (`MAT-6` step 10, 2026-08-12)

**Verified at `648b216`, 00:00 implementer slot.** The recipe every heavy
chunk uses — `... bash -lc 'source ... && timeout <s> mpiexec -n <N> python3
...'` — was run with `timeout 590`. It should have fired at 05:11:23Z. The
ranks were still burning 8–12 cores at 05:31Z, ≈1 700 s into a solve, and the
harness never wrote a footer
(`20260812T050133Z_MAT-6-step10-probe.log`). **Treat the container-side
`timeout` as best-effort, not a guaranteed stop**; the recipe had no
`-k`/kill-after, so a job that ignores or outlives SIGTERM keeps the cores.
Repair landed 2026-08-12 (03:00 daily review): the canonical recipes in
CLAUDE.md, PROJECT_PLAN §5, and docs/automation/ now read
`timeout -k 30 <s>`. This entry stays for the container-wedge behavior and
its recovery below, which the `-k` reduces but does not provably eliminate.

**The overrun then wedged the container**, and the usual levers failed in
order: `docker compose exec` hung twice with no output (>2 min each);
`docker compose restart` and `docker compose kill` both returned
`Error response from daemon: ... tried to kill container, but did not receive
an exit event`; a later exec failed with `OCI runtime exec failed: ... error
executing setns process: exit status 1`. **What recovered it:**

```bash
docker compose -f docker/docker-compose.yml up -d --force-recreate
```

Verified clean afterwards — exec responds, `/sys/fs/cgroup/memory.max` still
`68719476736`, zero stray `python3`, host load 12.2 → 8.9 as the orphaned
ranks died. Reach for the force-recreate rather than repeating the
restart/kill pair. Not diagnosed: whether the wedge is specific to a job that
survived its `timeout`, or to any long `docker compose exec` under load.

### A killed `mpiexec` job leaks its `/dev/shm/mpich_shm_*` segment, and once the container's 64 MB `/dev/shm` fills, **every** later run dies with `EXIT CODE: 135` (interactive session, 2026-08-28)

**Verified at `15e596f`, interactive session.** An `examples/run_examples`
invocation was interrupted by the human operator mid-run. Every subsequent
`mpiexec` — including ones that had been green minutes earlier — then aborted
immediately with

```
=   BAD TERMINATION OF ONE OF YOUR APPLICATION PROCESSES
=   PID <pid> RUNNING AT <container-id>
=   EXIT CODE: 135
=   CLEANING UP REMAINING PROCESSES
```

**135 = 128 + 7 = SIGBUS**, raised when MPICH touches a shared-memory page it
cannot back. `df -h /dev/shm` inside the container read **64M used, 0 avail,
100%**, holding **28** orphaned `mpich_shm_*` segments dated 2026-08-27 through
2026-08-28 — three of them 16.9 MB, from the interrupted run. No MPI process
was alive; nothing owned them. **Recovery is cheap** — no force-recreate, no
container restart:

```bash
docker compose -f docker/docker-compose.yml exec -T fem-em-solver \
  bash -lc 'rm -f /dev/shm/mpich_shm_*; df -h /dev/shm'
```

Verified clean afterwards: shm 0%, and `mpiexec -n 2 python3 -c "from mpi4py
import MPI; c=MPI.COMM_WORLD; print(c.rank, c.allreduce(1))"` returns `2` on
both ranks, exit 0.

**The segment size is a pure function of local rank count**, ≈1.06 MB/rank,
independent of problem size or mesh — measured live, and the 16-rank row is not
extrapolation (it is the size of the orphans the interrupted run left):

| ranks | segment | orphaned runs before `/dev/shm` is full |
|---|---|---|
| 1 | 1.06 MB | 60 |
| 2 | 2.12 MB | 30 |
| 8 | 8.46 MB | 7 |
| 12 | 12.69 MB (measured live at `15e596f`) | 5 |
| 16 | 16.92 MB | 3 |

So a *single* run of any legal width fits comfortably; what kills the container
is **accumulation**. A clean exit releases its segment, so this is purely a
kill/timeout aftermath — and each kill permanently burns a slot until someone
clears it. At the project's `-n 12` ceiling the **sixth** killed run poisons the
box; the standard `-n 2` recipes give 30. Note that `timeout -k 30` firing
(exit 124) is itself a kill and leaks a segment like any other.

**Do not read a post-kill `EXIT CODE: 135` as a regression** — same discipline
as the FFCx JIT-cache entry below, and the two compound: a killed run can leave
*both* a poisoned JIT cache and a leaked shm segment, so clear both before
trusting any failure that follows a kill. Suspect this whenever a previously
green command fails instantly with 135 and no Python traceback.

**Not fixed, and one-line to fix:** `docker/docker-compose.yml` sets no
`shm_size`, so the container is on Docker's 64 MB default. `shm_size: 2gb`
would raise the orphan headroom to ≈120 sixteen-rank runs and costs nothing
unused (`/dev/shm` is tmpfs — RAM is consumed only by pages actually touched).
Deliberately left unapplied here; it needs a container re-create to take effect.
Related: the **12-rank ceiling is not enforced against an interactive session** —
`scripts/automation/hooks/bash_guard.py:36-43` is a PreToolUse hook, so it
denies `mpiexec -n >12` only from agent sessions. The 16.9 MB orphans above are
what an unguarded human-typed `n` looks like.

### A killed harness run poisons the FFCx JIT cache, and the *next* run fails unrelated forms with "JIT compilation timed out" (`OPS-17` step 3 attempt 3, 2026-08-18)

**Verified at `2f97048`, 00:00 implementer slot.** Two complex-mode legs were
killed at their ceilings by `timeout -k 30 570` (exit 124,
`20260818T050123Z_OPS-17-step3c-complex-portgap.log` and
`20260818T051115Z_OPS-17-step3c-complex-remainder.log`). In the *second* of
those, `tests/solver/test_coil_phantom_magnetostatics.py::test_coil_phantom_magnetostatics_matches_the_two_loop_closed_form`
FAILED at 67% — a test that is green in real mode and whose gated quantity
(17.1233% L2 against a 30% band, `OPS-17` step 2) has nothing to do with the
build mode. Re-run alone it fails in **14.09 s** with

```
RuntimeError: Failed just-in-time compilation of form: JIT compilation timed
out, probably due to a failed previous compile.
Try cleaning cache (e.g. remove /root/.cache/fenics/libffcx_forms_<hash>.c)
or increase timeout option.
```

(`20260818T052132Z_OPS-17-step3c-coilphantom-complex.log`, exit 1, 15 s.) The
form's compile was interrupted mid-flight when the first leg was killed; the
lock file it left in `/root/.cache/fenics/` makes every later process wait out
the JIT timeout and raise. **The failure is an artifact of the kill, not a
regression** — do not open a chunk against the test on this evidence, and do
not trust *any* failure in a run that follows a killed one until the cache is
cleared. Suspect this whenever a fast, previously-green test fails in seconds
with a `dolfinx/jit.py` `RuntimeError` rather than an assertion.

**⚠️ PARTLY SUPERSEDED 2026-08-18 (`OPS-17` step 3 leg (b1), 07:30 slot).** The
cache-poisoning mechanism above is real and confirmed: `rm -rf
/root/.cache/fenics` is sufficient (no force-recreate needed), and the JIT
`RuntimeError` message does change with cache state. But the **conclusion drawn
about this particular test was wrong**. On a genuinely cold cache the test does
*not* return to green — it fails in **5.58 s** with a real complex-mode defect
(`ComplexComparisonError`, next entry). The poisoned cache was masking a
pre-existing failure, not creating one. Read this entry as "a killed run makes
the *message* untrustworthy", not "a killed run makes the *failure* spurious":
after clearing the cache you must still re-run and read the new message.

Measured cache-state → message map for this one test, all at `-n 2` complex:

| Cache state | Result | Message |
|---|---|---|
| poisoned by a kill mid-compile | FAILED 14.09 s | `JIT compilation timed out, probably due to a failed previous compile` |
| warm from a *completed* prior leg | FAILED 13.92 s | `Failed just-in-time compilation of form: Compilation failed on root node.` |
| cold (`rm -rf /root/.cache/fenics`) | FAILED **5.58 s** | `ComplexComparisonError: You can't compare complex numbers with max.` |

**Sizing corollary, measured 2026-08-18 (leg (b1) attempt 2, 09:00 slot).**
Clearing the cache is correct but it is not free, and it is the dominant cost
of a complex leg — not the solves. Same commit, same command shape, complex
`-n 2`: `tests/solver` on a **cold** cache exhausted a 480 s window at 61%
(`20260818T140137Z_OPS-17-step3e-complex-solver-tail.log`, exit 124), while on
a **warm** cache the same directory (minus the coil-phantom file) is
**46 passed / 2 xfailed in 111.22 s**
(`20260818T141104Z_OPS-17-step3e-complex-solver-warm.log`, exit 0) — ≈2.7× the
41 s real-mode number, i.e. the recorded 2.6× rule holds. The file the cold leg
died in, `test_gauge_penalty.py`, is 8 passed in 20.33 s standalone warm
(`20260818T141020Z_OPS-17-step3e-complex-gaugepenalty.log`), so a cold-leg
death location says nothing about which test is expensive. **Never infer a
per-test cost from a cold-cache run, and never let compilation and measurement
share one window** — size the first post-clear command as a throwaway warm-up.

### The complex-power identity reads 3–5e-9 at degree-2 N1curl on the coil fixture, against a 1e-9 family bound — and the quantity it gates is 99.6% spurious electric energy (`TH-12` step 2, 2026-08-18)

**Verified at `92fc3e7`, 15:00 implementer slot,
`20260818T200059Z_TH-12-step2-full.log`** (exit 1, 546 s, `-n 8`, complex
build, `TH12_STEP2_MODE=full`, 138 619 cells at 10 MHz).

| | |
|---|---|
| **Tests** | **Re-pointed 2026-09-01 (`TH-13` step 3a‴, the module split):** `tests/validation/test_coil_loading_degree2_pair.py::test_complex_power_identity_holds_at_degree2`, once per σ-half — `TH12_DEGREE2_HALF=loaded` and `=free`, one red in each of the two windows. The former names `tests/validation/test_coil_loading_degree2.py::test_complex_power_identity_holds_at_this_order[loaded-2]` / `[free-2]` no longer exist: that module is now degree-1-only and green. |
| **Symptom** | `AssertionError: complex-power identity broken on the loaded solve at degree 2: reaction -2.117210e+03 Ohm vs energy -2.117210e+03 Ohm, relative 4.5931e-09` (free: 3.0030e-09). The degree-1 rows of the same run, same mesh, same process are green at **8.0743e-15 / 8.7088e-15** — three orders of magnitude *inside* the bound, so this is an order effect, not a fixture or reduction defect. |
| **Cause** | Diagnosed by inspection of the printed energies, not fixed. `Im Z = 4ω(W_m − W_e)/I′²` is exact for the discrete solution, so the residual is arithmetic cancellation — and at degree 2 the cancellation is catastrophic because `W_e` explodes. Degree 1: `W_m` 3.04e-08 J, `W_e` **2.03e-13 J** ⇒ `Im Z` = **+9.02 Ω**, the physical loop reactance. Degree 2: `W_m` 3.13e-08 J (unmoved, 3%), `W_e` **7.16e-06 J** — 3.5e7× larger — ⇒ `Im Z` = **−2.117e+03 Ω**. The ungauged curl-curl operator's gradient null space is vastly richer at second order (882 296 DOFs vs 162 710) and the penalty-free formulation lets irrotational content sit in `E` at an amplitude that swamps the magnetic term. |
| **What it does and does not invalidate** | The spurious term is **common-mode**: it cancels in the loaded−free difference, so `ΔX` moves only −0.5666 → −0.5625 Ω (0.7%) and `ΔR` is unaffected. `TH-12` step 2's `ΔR` reading is therefore not contaminated by this. What is destroyed is the **identity's discriminating power at degree 2 on this fixture**: it now gates a number that is 99.6% non-physical, so passing it would mean little and failing it at 5e-9 means only that 2 117 Ω minus 2 117 Ω leaves 1e-5 Ω of round-off. |
| **Not** | **Not** to be fixed by widening `IDENTITY_TOLERANCE` — the bound is the `TH-11` step-2f family's, it is met at 1e-14 at degree 1 in the very same process, and widening it would hide the `W_e` explosion that is the actual finding. Not the magnetostatic degree-2 null-space failure either (different formulation, barred by `TH-12`'s scope guard) — though the two rhyme, and that rhyme is the hypothesis worth testing. |
| **Fix** | Not attempted; `TH-12` step 2 is a reading with no gate. The candidate dispositions, cheapest first: (a) re-anchor the identity at degree ≥ 2 on the **difference** `Im ΔZ` rather than the absolute `Im Z`, which is the quantity every downstream claim actually uses; (b) measure the gradient content of `E` directly (Helmholtz split, or `‖∇·(εE)‖` against a null-space projector) and report it as the degree-2 cost line; (c) treat it as a formulation defect and price a gauged/regularised second-order path. (a) is a test change, (b) is a measurement, (c) is a chunk. |
| **Resolves with** | `TH-12` **step 3** (commissioned 2026-08-18 18:00 review): the affordable form of disposition (b) — measure whether the `W_e` explosion is generic to incompatible drives (smoke fixture, `J·n ≠ 0`) or coil-feed-specific, at degrees 1/2 on the smoke + sphere fixtures at smoke cost. Dispositions (a)/(c) stay contingent on its reading — both need the 62 GiB coil solve to verify at degree 2. Until then this file **fails by default** (`TH12_STEP2_MODE` defaults to `full`), and `probe`/`calibrate` modes remain green. |
| **Step-3 reading (2026-08-19, `20260819T183425Z_TH-12-step3-warm.log`, `tests/validation/test_degree2_energy_mechanism.py`, 8 passed / exit 0 / 10 s at `-n 2`)** | **COIL-SPECIFIC** at the pre-registered ≤ 10×-on-both band. Cross-order move in `W_e/W_m`: smoke fixture (incompatible axial drive, `J·n ≠ 0` on the end caps, 1 405 cells) **1.155×** (2.164348 → 2.499688); lossy sphere (imposed field, 5 866 cells) **1.015×** (1.068190 → 1.052552); this coil **3.426e+07×** (6.677632e-06 → 2.287540e+02). So **`J·n ≠ 0` is not sufficient** — the incompatible-drive hypothesis in the "Cause" row above is refuted as *stated*, and the ungauged-gradient explanation now has to say why only this fixture displays it. **Confound the step could not separate, measured:** the fixtures' baseline `W_e/W_m` spans **2.16 / 1.07 / 6.7e-6**, so a contamination of fixed *absolute* size moves the quasi-static coil's ratio ≈1e6× more than either cheap fixture's; "the coil's feed model injects it" and "only a `W_m ≫ W_e` fixture can display it" both survive the reading. Splitting them needs a magnetically-dominated fixture with a compatible drive, or the absolute gradient content of `E` (disposition (b) proper, unscoped). **This entry stays open**; the two degree-2 identity tests stay failing at the unloosened 1e-9 bound, and no coil number moved. |
| **`TH-13` step-1 reading (2026-08-30, `20260830T020301Z_TH-13-step1.log`, `tests/validation/test_degree2_gradient_discriminator.py`, 1 failed / 12 passed / 1 skipped, exit 1, 36 s at `-n 2`)** | **The discriminator did not discriminate — a pre-registered negative result on both of its clauses, and the failing precondition assert is a deliberate red on `main`.** The fixture was `POST-5` step 2's closed azimuthal loop (`div J = 0`, `J·n = 0`) on the smoke box's own 1 405-cell mesh at 10 MHz, the missing magnetically-dominated-plus-compatible-drive cell of the step-3 table. (i) **Precondition failed:** degree-1 `W_e/W_m` = **1.952350e-02** against the pre-registered ≤ 1e-2 — a factor of 1.95 miss, so the fixture is *not* magnetically dominated and is not that cell. The number is unsurprising in hindsight: `W_e/W_m ~ ω²` at fixed impressed current, and the smoke box's own 2.164348 at 127.74 MHz scaled by (10/127.74)² predicts 1.33e-2. (ii) **Verdict IN-BETWEEN:** the loop's cross-order move is **5.156e+01×**, between the pre-registered 10× (FEED) and 1e3× (CLASS), so the reading is recorded and no band was invented in-slot (test skips, per §7). What the run *does* establish: **both step-3 controls reproduced to the digit on this code path** — smoke 1.155× (2.164348 → 2.499688) and sphere 1.015× (1.068190 → 1.052552) at the 1% band — so nothing has drifted since 2026-08-19, and `|Im P|/Re P` = 0.000e+00 at both orders on the loop. The suggestive part of the in-between number: degree 2 lifts the loop's `W_e` 63.7× (5.621559e-19 → 3.579741e-17 J) while `W_m` moves 1.23×, landing `W_e/W_m` at **1.006682** — i.e. at the same O(1) equipartition the smoke (2.50) and sphere (1.05) fixtures already sit at, and nowhere near the coil's degree-2 **229**. That is consistent with an absolute gradient contamination that saturates at equipartition, which neither band was written to catch. **This entry stays open**; nothing is fixed, no coil number moved, the two degree-2 identity tests stay failing at the unloosened 1e-9 bound. **Rescope hypothesis for the review** (not executed in-slot — §7 pins 10 MHz): the ω² scaling puts the precondition in reach at ≤ 7 MHz and puts it at ≈2e-4 at 1 MHz, which is also far enough below equipartition that a 1e3× move is *representable* — at 1.95e-2 a CLASS reading was arithmetically capped at ≈50× before the run started, so the fixture could not have returned CLASS whatever the physics. |
| **`TH-13` step-1′ reading (2026-08-30, `20260830T200543Z_TH-13-step1prime-final.log`, same module, 1 failed / 13 passed / 1 skipped, exit 1, 32 s at `-n 2`)** | **The ω² rescope is refuted by measurement: `W_e/W_m` on this fixture is frequency-INDEPENDENT, so frequency is not the knob that makes it magnetically dominated, and the precondition fails again at 1 MHz.** The step ran the same fixture, mesh, drive and forms at `LOOP_FREQUENCY_HZ = 1.0e6` alongside the 10 MHz row. Pre-registered in §7: degree-1 `W_e/W_m` ≈ **1.95e-4**, inside the unmoved ≤ 1e-2 band by 50×, with ≈5e3× of headroom so a 1e3× CLASS move would finally be representable. **Measured: 1.926692e-02** — **98.7×** the prediction and **0.9869×** the 10 MHz reading across a full decade of ω. The individual energies barely move either: `W_e` 5.621559e-19 → 5.544787e-19 J, `W_m` 2.879380e-17 → 2.877879e-17 J, dissipated power 1.139571e-09 → 1.124008e-09 W. Both `E` and `H` are therefore essentially unchanged by the frequency, and the quasi-static `E ~ ωA`-at-fixed-impressed-current argument that both the step-1 hindsight paragraph and this rescope rest on **does not describe this solve** — that assumption should not be reused without measuring it. **Consequences:** (a) the precondition assert stays red on `main`, now on the 1 MHz row instead of the 10 MHz one — the red moved, it did not multiply, and it was not loosened **[RETIRED 2026-09-01 by `TH-13` step 4: the assert was re-pointed at the matched-path solve at the *unchanged* 1e-2 band, where it reads 3.424858e-06, and this 1.926692e-02 became its lower-bounded control — see the step-4 row below]**; (b) the ω²-prediction clause was pre-registered as an assertion and is **demoted to a printed record** with the measurement in-comment (`OMEGA_SQUARED_PREDICTED_RATIO`), because a second red on a premise measurement has already killed gates nothing the precondition does not; (c) the **saturation trap §7 named fired**: the 1 MHz row's degree-2 `W_e/W_m` is **1.010649**, the same O(1) equipartition as 10 MHz's 1.006682, with a cross-order move of **5.246e+01×** against 10 MHz's 5.156e+01× — the contamination saturates at equipartition regardless of baseline, so **"cross-order move in `W_e/W_m`" is the wrong discriminant** and no choice of frequency will make it the right one. **Negative controls all green:** the 10 MHz row reproduced step 1 to rtol 1e-3 on both its degree-1 ratio (1.952350e-02) and its move (5.156e+01×); step 3's smoke 1.155× / sphere 1.015× at the 1% band; the `POST-5` anchor 1.199162e-06 W at rtol 1e-6; `|Im P|/Re P` = 0.000e+00 on all four loop solves. **This entry stays open**; nothing fixed, no coil number moved, the two degree-2 identity tests stay failing at the unloosened 1e-9 bound. **Reading for the review:** step 1′ has exhausted the frequency route, and the saturation says the *ratio* is the wrong observable — which points squarely at **step 2 (the absolute gradient content of `E`, `‖∇φ‖/‖E‖`)** as the next move, not a third fixture rescope. |
| **Review reading 2026-08-30 18:00 — a mechanism nobody had read off the code, pre-registered as `TH-13` step 2 (rescoped; §9 item 2)** | `core/source_projection.py::remove_gradient_content` (the `project_source=True` default, `PORT-1` step 2f) projects the drive against **`("Lagrange", 1)` ∩ H¹₀ only** (`:108`, `:134–136`), at every solve degree and under every boundary mode — while a degree-2 N1curl test space contains `∇Lagrange₂ ⊋ ∇CG1`, and the default `NATURAL` (PMC) boundary admits `∇q` for every `q ∈ CG`, not only `H¹₀`. Testing the assembled weak form with `v = ∇ψ` gives, exactly, `P_∇ₚ E_h = c · P_∇ₚ J′` with `c` fixed by the form's coefficients — on a homogeneous box `\|c\| = 1/\|σ + iωε\|`, frequency-flat for `σ ≫ ωε` (0.7 vs 0.0434 S/m at 10 MHz), which is the step-1′ reading. If the identity holds at both degrees, the degree-2 "explosion" is the drive's *unremoved* gradient residue answered by Ohm's law — common to the CG1-projected coil drive (residue in `∇Lagrange₂`) and to the unprojected port sheets (residue = their full discrete divergence, closed by the phantom's σ, which is physics) — and the fix is a degree-/boundary-matched projection (step 3), not a gauged formulation. The `‖∇φ‖/‖E‖` framing is retired before running: that ratio is bounded by 1 and a quasi-static conduction field is itself a gradient, so it cannot tell physical from spurious. **Nothing measured yet**; this row is a derivation the run is meant to refute. |
| **`TH-13` step-2 reading (2026-08-31, `20260831T021154Z_TH-13-step2.log`, same module, 1 failed / 15 passed / 1 skipped, exit 1, 32 s at `-n 2`) — DISPOSITION: (A) HOLDS. The injector is the degree-1-only, `H¹₀`-only source projection.** | The pre-registered identity `‖∇χ − c∇φ‖/‖∇χ‖ ≤ 1e-6` reads **2.970e-12 / 2.640e-11** (1 MHz, degrees 1/2) and **3.697e-13 / 2.586e-12** (10 MHz) — round-off at every row, four orders inside the bar, with the load-bearing probe firing exactly as pre-registered (`c` mistuned by 1.1× moves the residual to **1.000e-01** on all four rows, so the reading is not insensitive to `c`). `c` was built from the form's own coefficients and agrees with `−1/(σ + jωε₀εᵣ)` to `|c|·|σ + jωε| = 1.000000000` at both frequencies; its magnitude moves only 1.425834 → 1.428544 across a decade of ω, which **is** step 1′'s frequency-flatness, now explained rather than observed. **(B), recorded not gated — the mechanism carries the whole effect:** the gradient part of `E` accounts for **99.98%** of the measured `W_e` at degree 1 and **99.9997%** at degree 2 (1 MHz; 98.24% / 99.97% at 10 MHz), and the residue grows across order by `‖P_∇₂J′‖/‖P_∇₁J′‖` = **8.049884** at *both* frequencies — whose square, **64.8**, is the 63.7× degree-2 `W_e` lift step 1 measured. So the degree-2 electric energy this entry calls "spurious" is, to four digits, the drive's unremoved gradient residue divided by `(σ + jωε)`: there is no free "ungauged null-space mode" (`k² ≠ 0` pins it), and the ungauged-second-order-gradient-space story in the "Cause" row above is **superseded**. **Projection control (`PORT-1` step 2d/2e precedent):** an extra `project_source=False` degree-1 solve reads `‖P_∇₁J‖/‖J‖` = **7.589863e-02** against the projected drive's **1.298386e-02** — `remove_gradient_content` removes 5.8× of the CG1 content and leaves the rest, which is exactly the non-`H¹₀` part a PMC box still tests. **Disposition of this entry:** the diagnosis is closed, the *fix* is not — a degree-/boundary-matched projection is `TH-13` step 3, a `src/` change a review prices and then re-measures on the coil (where `TH-12` step 2 read 229×). Until it lands this entry **stays open**, the two degree-2 identity tests stay failing at the unloosened 1e-9 bound, and no coil number moved. The step-1′ precondition assert also stays red on the 1 MHz row — a deliberate red held per §7, not this step's to retire. **[Retired 2026-09-01 by step 4, on this step's own identity: the 1.926692e-02 it was asserting is the injector's residue, not the fixture.]** |
| **`TH-13` step-3a reading (2026-08-31, `20260831T094852Z_TH-13-step3a-final.log`, same module, 1 failed / 18 passed / 1 skipped, exit 1, 37 s at `-n 2`) — THE FIX LANDS, AND IT DOES NOT REACH THIS ENTRY'S FIXTURE.** | The matched projection is now in `src/`: `remove_gradient_content` takes `degree` (default 1) and `pin_exterior` (default `True`), and `TimeHarmonicSolver.solve(project_source="matched")` sets them from the solve's own degree and boundary mode. On the loop fixture it does exactly what step 2's identity predicted. **(i)** `‖P_∇ₚJ′‖/‖J′‖` falls from **1.298386e-02** (degree 1) and **1.045186e-01** (degree 2) to **8.109635e-17** and **1.790460e-16** — a separation of **1.601e+14×** / **5.838e+14×** against the pre-registered ≤ 1e-8 and ≥ 1e3× bands. **(ii)** the gradient share of `W_e` falls from 99.98% / 99.9997% to **4.618447e-23** / **3.109722e-21** (1 MHz) and **4.390171e-25** / **2.738225e-23** (10 MHz), against ≤ 1e-6; and `W_e` itself collapses from 5.544787e-19 → **9.856327e-23 J** (degree 1) and 3.592428e-17 → **9.349492e-23 J** (degree 2) at 1 MHz, i.e. to 0.018% and 2.6e-4 % of record against the pre-registered ≤ 2% / ≤ 1%. **Recorded, not gated:** the matched degree-2/degree-1 `W_e` ratio is **9.485777e-01×** where the default path's was 63.7× — the 63.7 *was* the residue, and what is left is the solenoidal electric energy of an unconverged 1 405-cell loop, cross-order-flat. `W_e/W_m` under `"matched"` is 3.424858e-06 / 2.630270e-06 at 1 MHz (the step-1′ precondition is **not** re-asserted there); `\|Im P\|/Re P` = 0.000e+00 on every matched row. **Default path unmoved:** control (b) — PEC at degree 1, where `"matched"` *is* `True` — reads `W_e` = 5.995936714066138e-23 J on both paths, **0.000e+00 relative**; and every step-1/1′/2 assertion above reproduces to the digit (residuals 2.970e-12 / 2.640e-11 / 3.697e-13 / 2.586e-12, mistuned 1.000e-01, `‖P_∇₂J′‖/‖P_∇₁J′‖` = 8.049884, unprojected control 7.589863e-02, `POST-5` 1.199162e-06 W). `test_dodd_deeds_projected_drive.py` re-run green, **15 passed / exit 0 / 79 s** (`20260831T093558Z_TH-13-step3a-dodd-regression.log`), matching its 2026-08-27 record. **THIS ENTRY STAYS OPEN, and 3a cannot close it:** the lumped-sheet birdcage drive is a *surface* term with `current_density=None` and `project_source=False` (`ports/lumped.py:429`), so the matched projection never touches the coil's 229× — that residue is the sheets' own discrete divergence closed by the phantom's σ, and whether it should be removed at all is a formulation ruling for the weekly review (a possible step 3b). The two degree-2 coil identity tests stay failing at the unloosened 1e-9 bound, **unverified on this commit** — see the entry below; the step-1′ precondition red also stays on the default path **[retired 2026-09-01 by step 4 — the row below]**. |
| **`TH-13` step-4 reading (2026-09-01, `20260901T110318Z_TH-13-step4.log`, same module, **19 passed / 1 skipped, Status 0, 38 s at `-n 2`, complex, standard tier**) — THE STEP-1′ PRECONDITION RED IS RETIRED BY RE-POINTING, NOT BY LOOSENING; THE MODULE GOES EXIT 0 FOR THE FIRST TIME.** | The 2026-08-31 10:30 review ruled (PROJECT_PLAN §7 `TH-13` step-4 bullet) that the precondition was measuring the **injector**, not the fixture: step 2's identity says the default path's degree-1 `W_e/W_m` = **1.926692e-02** is, to four digits, the CG1∩H¹₀ projection's unremoved gradient residue answered through Ohm's law. Executed: `test_the_loop_fixture_is_magnetically_dominated` now asserts the **matched-path** degree-1 reading **3.424858e-06** against the **unchanged** ≤ 1e-2 band — **2.920e+03× inside**, reproducing `20260831T094852Z_TH-13-step3a-final.log:3683` to every digit — and the default-path 1.926692e-02 is asserted in the same test as a **lower-bounded control** (`> 1e-2`, plus rtol 1e-3 against the step-1′ record, `STEP1PRIME_DEFAULT_DEGREE1_RATIO`), so the residue stays *measured* on `main` rather than vanishing: a change that silently removes it on the default path now turns this module red. **No band, tolerance or record moved, no new test function** (20 collected, as before), **`git diff -- src/` empty**. Every other reading reproduces to the digit: cross-order moves 5.246e+01 / 5.156e+01× (loop 1 / 10 MHz) and 1.155 / 1.015× (smoke / sphere), step-2 residuals 2.970e-12 / 2.640e-11 / 3.697e-13 / 2.586e-12 with the mistuned-`c` probe at 1.000e-01, `‖P_∇₂J′‖/‖P_∇₁J′‖` = 8.049884 at both frequencies, projection control 7.589863e-02 vs 1.298386e-02, matched degree-2/degree-1 `W_e` 9.485777e-01×, `\|Im P\|/Re P` = 0.000e+00 on every row. **THIS ENTRY STAYS OPEN:** step 4 retires only the step-1′ precondition red (`main`'s residual deliberate reds at `-n 2` go 9 → 8). The two degree-2 **coil** identity tests still fail at the unloosened 1e-9 bound — observed 2026-09-01 by step 3a‴ at 3.8990e-09 / 3.7235e-09, one per σ-half — no coil number moved, and the sheet-drive residue (step 3b) remains a formulation ruling for the weekly review. |
| **`TH-19` steps 1–2 reading (2026-09-10, 22:30 implementer slot, `-n 8`, complex, `timeout -k 30 600`) — ON THE COIL, THE MATCHED PROJECTION REMOVES THE EXPLOSION: OUTCOME (a) ON BOTH σ-HALVES.** | **Step 1, default path, loaded** (`20260911T033210Z_TH-19-step1.log`, 1 failed / 17 passed / 2 skipped, Status 1, 475 s): the red reproduces — relative **2.3898e-09** vs 1e-9 (`:538`), `W_e` 7.859344e-06 J, `W_m` 3.135703e-08 J, `W_e/W_m` 2.506405e+02 (`:437`); energies and `Im Z` equal the 2026-09-01 record to every printed digit, the residual digits moved (3.8990e-09 → 2.3898e-09, cancellation round-off). **Step 2, `TH19_PROJECT_SOURCE=matched`** (degree-2 half only; same module, same unmoved 1e-9): loaded `20260911T034033Z_TH-19-step2-loaded.log`, 18 passed / 2 skipped, Status 0, 484 s — residual **2.0421e-14**, `Im Z` **+9.307547 Ω**, `W_e` **2.043905e-13 J**, `W_m` 3.135703e-08 J, `W_e/W_m` **6.518171e-06** (`:398–399`); free `20260911T034858Z_TH-19-step2-free.log`, 18 passed / 2 skipped, Status 0, 483 s — residual **5.9380e-15**, `Im Z` **+9.871903 Ω**, `W_e` **3.405566e-13 J**, `W_m` 3.325847e-08 J, `W_e/W_m` **1.023970e-05** (`:420–421`). `W_e` falls 3.845e7× (loaded) / 2.308e7× (free, against the 2026-09-01 free record), `W_m` is unmoved, and `Im Z` returns to the physical sign and scale (degree-1 default +9.017628 / +9.584200 Ω). Degree-1 ΔR control +1.5838% (+0.00039 pp) and degree-1 residuals ≤ 1.3e-14 in all three windows. So the coil's degree-2 `W_e` explosion **was** the CG1∩H¹₀ projection's residue, as `TH-13` step 2 derived on the loop. **THIS ENTRY STAYS OPEN:** the test only goes green under the env; on the default path (`project_source=True`, what `main` runs) the two identity tests still fail at 1e-9, and the birdcage's lumped-sheet drive bypasses the projection altogether (`ports/lumped.py:480–483`). Changing the default or retiring this entry is the weekly's. |
| **`TH-19` step 3 reading (2026-09-13, 13:30 implementer slot, `-n 8`, complex, `timeout -k 30 590`) — ON THE SHEET-DRIVEN BIRDCAGE, DEGREE 2 IS IDENTITY-CLEAN WITHOUT ANY PROJECTION, AND THE 3a EXPECTATION ABOVE IS NOT BORNE OUT.** *(Row added by the 2026-09-13 18:00 review.)* | `tests/validation/test_birdcage_power_identity_degree2.py`, the 116 085-cell 4-leg fixture, single P1 sheet drive, bypass untouched (`20260913T183446Z_TH-19-step3-10MHz.log`, 15 passed, 134 s; `…183723Z_…-128MHz.log`, 15 passed, 118 s; readings `:1946–1955` in each): degree-2 exact power identity (a) **4.575e-15 / 8.928e-15** vs 1e-6 and reactive identity (b) **6.789e-11 / 1.816e-12** vs the imported 1e-9 at 10 / 128 MHz; degree 1 (a) 9.656e-16 / 3.527e-16, (b) 2.975e-12 / 4.734e-14, reproducing `PORT-16` step 1's records to ≤ 1e-10. `W_e/W_m` 2.343550e-02 → 2.228745e-02 (10 MHz) and 1.209660e-01 → 1.174350e-01 (128 MHz) from degree 1 to 2 — **`W_e` falls ≈ 7 %**, it does not explode. So the 3a clause above ("that residue is the sheets' own discrete divergence") is not what the sheets do: the surface drive carries no coil-class `W_e` excess at degree 2. This entry's coil-fixture reds are untouched by it; the production-order decision is the 09-16 weekly's. Caveat: (b)'s 10 MHz degree-2 margin is 15× against a band whose round-off has moved 3.9e-09 → 2.4e-09 between runs on the coil (row above) — one observation. |

### The gauge penalty is a workaround, not a gauge — closed by decision (2026-07-28)

`DEFAULT_GAUGE_PENALTY = 1.0` prices the curl-curl operator's gradient null space
rather than removing it (`PROJECT_PLAN.md` §7 `MAG-10`). This was previously the
highest-risk open item; the decision is recorded here so it is not reopened:

- `MAG-15` landed `GaugeMethod.LAGRANGE` — an (A, p) saddle point that removes the
  null space with no parameter — as a cross-check and diagnostic
  (`tests/solver/test_gauge_lagrange.py`). The penalty at 1.0 stays the production
  default on cost grounds (≈2× at degree 1, ≈7.5× at degree 2).
- Tree-cotree gauging is rejected: `TH-1`'s E-field formulation has no static null
  space at ω > 0 (the operator acts as −k₀²ε_c on the gradient subspace), so deeper
  magnetostatic gauge machinery has no Phase-2 payoff.
- The risk does **not** transfer to `TH-1` as gauge cancellation. The Phase-2
  silent-failure analog is *near-resonance ill-conditioning*, tracked in the
  `TH-1` formulation notes in `PROJECT_PLAN.md` §7.

### Air-box sizing is not generalised

Only `two_torus_domain` has `air_padding` and graded sizing. Every other fixture in
`io/mesh.py` — including the coil+phantom geometry all MRI work depends on — still
uses a single global `setSize` with padding tied to feature size. On Helmholtz that
pattern produced a **20.4% error that did not improve across a 7× refinement**, and it
was invisible until an analytic comparison existed. Coil+phantom has no analytic
reference yet. See `PROJECT_PLAN.md` §9.


### Every XDMF/VTX export of a Nédélec or DG field ships a Lagrange-P1 interpolant that disagrees with the solved field at O(20–52%) — everywhere in the cell, not only at vertices (`POST-4` step 4, 2026-08-12)

**What this affects.** Eleven interpolation sites in `examples/` build a
`("Lagrange", 1, (3,))` (or scalar `("Lagrange", 1)`) function for export;
**ten** of them are fed by a source that is not continuous in that space —
N1curl (`A`, `E`) or DG (`B`):

| file | line | source |
| --- | --- | --- |
| `examples/mri/01_coil_phantom_fields.py` | 413 | `A` N1curl, `B` DG1, `E` N1curl |
| `examples/magnetostatics/01_straight_wire.py` | 185 | `B` DG1 — **evaluated, not only exported** (see below) |
| `examples/magnetostatics/02_circular_loop.py` | 259 | `A` N1curl, `B` DG1 |
| `examples/magnetostatics/04_helmholtz_analytic_comparison.py` | 136 | `A` N1curl, `B` DG1 |
| `examples/magnetostatics/05_gauge_cross_check.py` | 198 | `A` N1curl, `B` DG1 (both gauges) |
| `examples/time_harmonic/01_lossy_plane_wave.py` | 109 | `E` N1curl |
| `examples/time_harmonic/02_pec_cavity_resonances.py` | 145 | eigenmode, N1curl |
| `examples/time_harmonic/03_dielectric_sphere_in_uniform_field.py` | 258 | `E` N1curl |
| `examples/time_harmonic/04_evanescent_waveguide_decay.py` | 246 | `E` N1curl |
| `examples/time_harmonic/05_resonance_guard_sweep.py` | 127 | `E` N1curl |

`examples/magnetostatics/06_h_convergence_rate.py:164` is the exception that is
*safe by construction*: it exports the CG1 function it also asserts on, so the
asserted field and the exported field are the same object.

**Measured magnitude** (`examples/mri/01` debug preset, 9261 cells, `-n 2`,
400 cell midpoints + 400 vertices, `20260812T003454Z_POST-4-step4-anchor-n2.log`,
exit 0, 4 s). Pointwise relative median of |P1 − source| / |source|:

| field | midpoints | vertices | scaled median (mid / vtx) |
| --- | --- | --- | --- |
| `A` (N1curl) | **51.17%** | 27.33% | 0.1032 / 0.0432 |
| `B` (DG1) | **52.47%** | 38.39% | 0.1590 / 0.0766 |
| `E` (N1curl) | **20.18%** | 15.79% | 0.1633 / 0.1116 |

**Two things this measurement settled, both against the prior expectation.**
(1) The artifact is *not* localized at shared vertices — the midpoint
disagreement is **larger** than the vertex disagreement in all three fields
(separation 0.4185× / 0.4818× / 0.6835×, where the step-4 entry predicted
O(50×) the other way). The wrong vertex dofs define the P1 interpolant over
the whole cell, so the interior inherits them; a vertex sample can also, by
chance, draw the same cell trace on both paths, which is why vertices read
*quieter*. (2) None of it is interpolation error. Interpolating the same three
sources onto a **DG1** target — same degree, no dofs shared between cells —
reproduces them to round-off (scaled median 3.25e-17 / 0.0 / 0.0), because all
three are degree-1 discontinuous polynomials and are represented exactly. The
entire disagreement is the P1 continuity constraint. Negative control: a
conforming P1 source round-tripped through the same machinery agrees to
**0.000e+00** against a 1e-10 bound.

**Consequences.** ParaView renderings of `A`, `B`, and `E` from these exports
are a continuous *approximation* of a field that is not continuous, off by tens
of percent pointwise at debug resolution; they are qualitative pictures, not
data. No gate cites them. The one case that goes beyond visualization is
`examples/magnetostatics/01_straight_wire.py:185`, which interpolates `B` to P1
and then **evaluates the radial profile from the interpolant** — that printout
carries this artifact and should be read as indicative only; the `MAG-13`
convergence numbers do not come from it.
*(Superseded 2026-10-21, noted by the 2026-10-23 10:30 review: since
`a49b4e16` (`POST-8` step 3b-ii-e3b, §9 item 580) `mag:1` samples the DG
source `B` under the owner rule through `evaluate_vector_field_parallel`
(`examples/magnetostatics/01_straight_wire.py:274–276`); `B_lag` is kept only
as the export, so the printout no longer carries this artifact. The site
table above dates from 2026-08-12 (11 sites) and is stale — a grep of
`examples/` now finds ≈ 82 Lagrange-1 constructions; the re-census is §9 item
607. The entry itself — the P1 export's disagreement and the DG1-vs-P1
decision — stays open.)*

#### Re-census 2026-10-23 (`POST-4` §9 item 607)

Probe `scripts/probes/post4_p1_site_census.py` (stdlib `ast`, no dolfinx),
log `docs/testing/logs/20261023T171007Z_POST-4-census.log` (smoke, host
`python3`, exit 0, 2 s, at `0c47b4e4`). It walks 151 `examples/**/*.py` files
and finds **82** `functionspace(…, ("Lagrange"|"CG"|"P", 1 …))` constructions;
the independent line-regex control in the same window also counts **82**
(assert (ii) PASS). The brief's literal grep counts 83 lines — the extra one is
prose in the `examples/mri/01_coil_phantom_fields.py:45` docstring, not a
construction. All ten August file paths still exist (assert (i) PASS); their
sites have moved to the current lines `mri/01:437`, `magnetostatics/01:251`,
`/02:315`, `/04:157`, `/05:218`, `time_harmonic/01:111,121`, `/02:147,160`,
`/03:264,274`, `/04:268,278`, `/05:129,137` (the second line in each
`time_harmonic` file is the scalar `|E|` built from the vector interpolant).
The full per-site row list (space, `Function`, `.interpolate` source) is in
the log. Classification by reading the code (executor, 2026-10-23):

| source | (x) export only | (e) evaluated / asserted | total |
| --- | --- | --- | --- |
| (N) N1curl / `ufl` expression of one | 52 | 2 | 54 |
| (D) DG present (with or without N1curl) | 2 | 4 | 6 |
| (C) conforming (Lagrange, CG1 projection, analytic) | 21 | 1 | 22 |
| (U) unknown | 0 | 0 | 0 |
| total | 75 | 7 | 82 |

Where one space carries several interpolants the row is classed by its worst
source (N/D over C) and by any use (e over x). The (C) rows are the thermal
exports (P2 temperatures, closed-form series), the `|B₁⁺|` maps built from
CG1 L2 projections (`ports/06`, `/07`, `/08`, `/16`, the `mri/19`–`/26`
tuning maps), and `meshing/13:225` (an analytic planted field, the one (C)(e)).
The 54 (N)(x) rows are the `_paraview_fields` / `_export` boilerplate of the
time-harmonic, port, implant and Study examples (`E` or `|E|` of an N1curl
phasor). `mri/14:367` carries both a CG1-projected `|B₁⁺|` (C) and `|E|`
(N) — counted (N)(x).

**Every (N)- or (D)-sourced site in class (e)** — six; none is fixed here:

| site | source | what is done with the P1 interpolant | what it prints / asserts |
| --- | --- | --- | --- |
| `examples/ansys_benchmarks/ans2_birdcage_coil_driven_sar_10MHz/02_birdcage_coil_driven_sar_10MHz.py:336` | `E` N1curl | `_pointwise_sar` samples the CG1 `e_cg` at the four named points through `evaluate_vector_field_parallel` ("legacy CG1" SAR row) | the legacy CG1 SAR row, its driven-diagonal spread, and `cg1_separation` against the owner-cell row (the gated row 1 is the owner-cell average, not this one); asserts `cg1_separation.min() > PLANT`, i.e. the CG1 row is asserted to **differ** from the gated one |
| `examples/time_harmonic/05_resonance_guard_sweep.py:129` | `E` N1curl | `_magnitude_field` reduces the vertex `|E|` of the CG1 interpolant to a global peak | `[paraview] \|E\| peaks at … V/m` near / quiet and their factor; asserts `peak_ratio > 1.0` (a sign test, not a magnitude) |
| `examples/magnetostatics/06_h_convergence_rate.py:192` | `B` DG | the CG1 `b_num_cg` is sampled under the owner rule against the analytic wire field | `exported fld : x% at the same points` (export error and its offset from the finest-mesh error); asserts `export_error < errors[0]` |
| `examples/magnetostatics/05_gauge_cross_check.py:218` | `B` DG (and `A` N1curl), both gauges | volume L2 of the difference of the two gauges' CG1 `B` interpolants | `volume L2 rel diff (export): x%`; asserts `< VOLUME_AGREEMENT_RTOL` |
| `examples/magnetostatics/01_straight_wire.py:251` | `B` DG (and `A` N1curl) | `_check_vtx_roundtrip` compares the global max `\|B\|` of `B_lag` with the `.bp` read back | `in-memory max\|B\| = …` and the read-back; asserts round-off agreement — a file-fidelity identity, P1 against its own file, so it carries no artifact into a physics number |
| `examples/magnetostatics/02_circular_loop.py:315` | `B` DG (and `A` N1curl) | same `_check_vtx_roundtrip` as `mag:1` | same as `mag:1` |

The first four print a number derived from a P1 interpolant of a
discontinuous source; of those, `ans2:336` and `magnetostatics/06:192` sample it
pointwise, which is exactly the regime the measurement above puts at tens of
percent at debug resolution. No gate in the list compares such a number to a
closed form within a tight band (the `ans2` and `th:5` asserts are separation /
sign tests; `mag:6` asserts the export is no worse than the coarsest-mesh
error). Whether a guide quotes any of these printouts as a number was not
audited in this slot. The DG1-vs-P1 decision stays with the review.

**Per-site disposition (`POST-4` R2, §9 item 636, 2026-10-25).** Guides audited
against the six class-(e) sites: `ans2:336` — none (the guide already labels the
CG1 row a record); `mag:6:192` — none (the page quotes the export's own error,
labelled "exported", owner-sampled since `POST-8`); `mag:1:251` / `mag:2:315` —
none (file-fidelity identities). **`th:5:129`** — re-quoted: the guide's peak
`|E|` (6.0531e+03 / 6.1951e+02 V/m, factor 9.77) was the P1 export's; the native
N1curl field at the mesh vertices through `owner="lowest_original_cell"` gives
5.8778e+03 / 5.9984e+02 V/m, factor 9.80 (`20261025T183230Z_POST-4.log`); the
export now prints 6.0350e+03 / 6.1719e+02, factor 9.78. **`mag:5:218`** —
re-quoted: the guide's 0.0033% volume L2 difference was labelled the export's;
on the native DG `B` fields it is 0.0033%, and the CG1 export now reads 0.0040%
(`20261025T183305Z_POST-4.log`; native volume integral assembled on the DG
functions, since an integral is not a point evaluation). **Control missed at
both sites:** the export line does not reproduce the old guide digits (th:5
6.0350e+03 vs 6.0531e+03; mag:5 0.0040% vs 0.0033%) — the guides' export
figures are from older logs and the P1 interpolant has moved with the image; no
assert depends on either number. The sites' existing asserts are green and
unedited. *Review ruling 2026-10-25 18:00:* the control holds on the current image —
both P1 lines reproduce the 2026-08-25 post-image logs digit for digit, and the guide
digits are 2026-08-10 records from before the 08-22 image change (`POST-4.md`, item 636).

**Scope of this entry.** Measurement only — no export code was changed, and no
ParaView claim elsewhere in the repo is withdrawn by it. Whether to export DG1
instead (faithful, larger files, discontinuous rendering) or to keep P1 with
this caveat is a review decision, not this step's. Probe:
`scripts/probes/post4_step4_probe.py`.

## Recording a new entry

Add an entry when you find a failure you are **not** fixing. Include: the test id, the
literal symptom, the commit you verified it against, the cause (or an explicit "not
diagnosed" — an honest gap is more useful than a guess), and which chunk resolves it.

Remove the entry in the same commit that fixes the test.

If the retirement itself needs recording (digits, log lines, what was ruled),
rewrite the entry's **heading** to begin `RETIRED YYYY-MM-DD` and leave the
body: the weekly review moves every such entry, verbatim, to
`known-issues-archive.md` and leaves one index line under "Retired entries"
below (`rotate_plan_archive.py known-issues`, since 2026-09-19). This file
holds OPEN entries only — that is what makes "is this failure mine?" quick to
answer.

## 2026-09-20 — **A *piped* or subshell-wrapped harness call is denied the docker socket in ≈1 s; a plain `cd <repo> && run_and_log.sh …` is NOT. Never an outage — never park a slot on it** — `permission denied while trying to connect to the docker API at unix:///var/run/docker.sock`

**The rule, in one line:** invoke `scripts/testing/run_and_log.sh …` as the
whole command with **no pipe and no `bash -c`** wrapper, and read the log with
the Read tool instead of `| tail`. A `cd <repo> &&` prefix is fine.

**Why this entry is long:** two executors in the same slot reached two *wrong*
conclusions from this denial — the first that it was a project-wide outage, the
second that any compound command triggers it. Both are refuted below by
measurement. The history is kept because the first reading cost a whole
attempt.

**Seen** by the 2026-09-20 04:30 implementer slot's delegated executor
(`TH-17` step 1b) at `1ba7302`: `20260920T093410Z_TH-17.log` — Status 1,
Elapsed 1 s, no test ran, `## Output` is the one denial line. The 03:00 review
met the same denial the same morning
(`20260920T080639Z_OPS-53-callsite-pin.log`, Status 1, 0 s).

**Disproof of the "outage" reading — read this before parking a slot.** The
executor concluded from the two denials that the harness could reach no socket
at all and that no scheduled slot could execute any compute; it parked its item
on that basis. The slot owner re-ran the *identical* command minutes later in
the same session and it ran **green**: `20260920T093731Z_TH-17.log` (a trivial
`echo` probe through the harness, Status 0) and then
`20260920T093810Z_TH-17.log` — the byte-identical window the executor had
declared impossible — **Status 0, 232 s, 3 passed** at `-n 4`. So the denial is
**transient and per-invocation**, not a configuration boundary, and the
"sandbox binds the socket per allowlisted command, so harness children get
none" cause read is **wrong**: harness children reach the socket routinely.

**Mechanism (the `PORT-20` executor's diagnosis, corrected).**
`.claude/settings.json` exempts the harness from the OS sandbox by prefix
(`sandbox.excludedCommands: ["docker *", "scripts/testing/run_and_log.sh *",
…]`), and the sandbox's write allowlist does **not** include
`/var/run/docker.sock`. A harness call that loses that exemption runs sandboxed
and the socket write is refused in ≈1 s, before any test runs. What loses it,
measured in this slot:

- **A pipe loses it.** `20260920T095034Z/095043Z/095130Z_PORT-20.log` and
  `…095149Z_PORT-20.log` — Status 1, ≤ 1 s, all piped (the last one with a bare
  `scripts/testing/…` prefix and still denied). The same window text with the
  pipe removed: `…095207Z_PORT-20.log`, Status 0, 2 s, 3 passed.
- **An explicit `bash -c` subshell loses it.** `docker compose … exec -T
  fem-em-solver bash -lc 'echo hi'` succeeds directly, while
  `bash -c '<the same docker command>'` is denied.
- **A `cd <repo> &&` prefix does NOT lose it** — this is where the executor's
  write-up was wrong, and it named that shape explicitly. Three green windows in
  this slot were invoked as `cd /home/…/fem-em-solver && scripts/testing/run_and_log.sh …`:
  `20260920T093731Z_TH-17.log` (Status 0), `20260920T093810Z_TH-17.log`
  (Status 0, 232 s, 3 passed) and a probe run specifically to settle it,
  `20260920T100019Z_PORT-20.log` (Status 0). So compounding is not the trigger;
  piping and subshell-wrapping are.

The original `TH-17` denial (`…093410Z`) and the 03:00 review's
(`…080639Z_OPS-53-callsite-pin.log`) have no surviving record of their outer
wrapper, so they are consistent with the pipe mechanism but not proof of it.
Whether this is also what PROJECT_PLAN §9 records for `./run_examples.sh`
(2026-08-29, 08-30, 09-01) is **not** established.

**What to do.** Write the call bare and unpiped; a `cd` prefix is safe. If a
window is denied anyway, **retry it once, unpiped** — do not park an item, do
not write an infrastructure blocker, and do not propose an allowlist change on a
≈1 s Status-1 window. That misreading cost the 04:30 slot its first executor's
entire attempt; the item's (i) and (ii) then ran green, unchanged, in one
window. Bypassing the harness stays correctly refused by
`scripts/automation/hooks/bash_guard.py`, and no operator action is needed.

## 2026-10-17 — 🟡 OPEN (`ANS-7` T41, weekly review; §10 chain T) — the tuned F-small birdcage's **Q and `R_in` converge in h to a level ≈ 3 % below AED First Order** (Q), the direction T13 read as "away"; Palace at matched order reproduces ours to six digits, so the level is shared by both FEM codes — INCONCLUSIVE as an absolute Q claim until AED's own pass history is read

| | |
| --- | --- |
| **Symptom** | None in CI — every `ANS-7` gate is green. The gap is against the external code only, and it does not close with h: Cu rows at `C_F` on h = 0.015 / 0.01 / 0.005 m read Q 977.26 / 970.44 / 965.55 and `R_in` 18 441.75 / 18 313.60 / 18 223.23 Ω (`20261011T094910Z_ANS-7-item364-T30.log`), successive-change ratios 0.717 / 0.705 on steps where first order gives 1.00 and second order 0.60 — converging at T30's printed p = 1.645 / 1.679, Richardson limits Q∞ ≈ 963.3 and `R_in`∞ ≈ 18 182 Ω (0.24 % / 0.23 % below the finest rung). |
| **What agrees** | `f₀` (within 0.11 % at ×1, 0.03 % at h = 0.005 — T13); the mode identity; the copper contribution to 1/Q (within 5 %); Palace's `f₀` / Q / `R_in` at matched order on ×1 and h = 0.01 at both Larmor frequencies (T39–T44), including the ×1 → 0.01 Q move to five / six significant digits. The copper-loss partition on the `ANS-6` hole fixture converges the same way (F31: 0.4760 → 0.4611 → 0.4576 on h = 0.015 / 0.0075 / 0.005, ratio 0.23 against first-order 0.33 / second-order 0.19). |
| **What does not** | Q and `R_in` against AED First Order: ≈ 3 % low at h = 0.005 (T13's relative reading) and the h-limit sits at the same level. The T13 split put the gap in the phantom's part of 1/Q at degree 1; at degree 2 the phantom part was within 0.3 % and the conductor part 4.4 % high — the degree-2 attribution is not re-read tonight. |
| **What is and is not known** | Ours converges; AED's First Order is one discretisation whose pass-to-pass trend on this design is unread. A 3 % Q gap is 14–19 linewidths' worth of nothing in `f₀` (which agrees) and is the size of the conductor part's own 4.4 % excess. Not a sheet effect (exact rectangles, κ off), not an order effect (both orders read the same direction), not a Palace-vs-ours effect (tied). |
| **Consequence** | **No Q, `Z_in` or match figure of ours is quoted as agreeing with AED.** `ANS-7` is ✅ on its Done-when (the h-statement is written) with the absolute Q INCONCLUSIVE — not DISAGREE (one external column, its own convergence unread; §10 epitaph 2026-09-19 on verdict classes). |
| **Next** | §10 T48 (operator, dashboard W7): the existing `ans:7` AED project's adaptive-pass history of `f₀` / `R_in` / Q at 64 MHz, read privately by the weekly; AED moving toward ours ⇒ the level is AED's discretisation and this entry retires; AED flat or moving away ⇒ T46 (the strip-edge loss band on our side) runs. Numbers: gitignored `docs/private/t13-t29-ans9-adjudication-2026-10-10.md` (T13's reading, re-used tonight — no AED figure was opened). |

## 2026-10-27 — 🟡 OPEN (`GEO-38` C15, §9 item 655, 04:30 CDT implementer slot; observed, not fixed) — F-clinical-3T, class "other ring, 1 step": absolute spread **1.1235×** the adjacent class's ("same ring, 1 step")

- **Test id:** none red — `tests/validation/test_geo38_c15_3t_c_mode.py::test_geo38_c15_3t_c_mode` (`FEM_EM_GEO38_C15=1`, `-n 1`) prints the ratio and asserts nothing on it (the brief's rule); every anchor in the window is green.
- **Literal symptom:** `20261027T093207Z_GEO-38-C15-W1.log:87` — `anatomy other ring,  1 steps  n =  128  mean |S| 9.942989475e-02  abs spread 5.806345e-04  rel   0.5840%  abs / (same ring, 1 step) 1.1235`, against the reference `:80` ("same ring, 1 step" abs spread 5.168109e-04, rel 0.1383 %). SUMMARY `:93`: 12 / 34 classes above `ADJACENT_SPREAD_BAND` 0.005, largest ratio 1.1235, one class with ratio > 1. The next largest is "other ring, 2 steps" at 0.8650 (`:88`; the brief's arithmetic predicted every other qualifying class ≤ ≈ 0.6 — a reading outside that prediction, not a red); the other ten are ≤ 0.5280.
- **Matrix custody (green in the same window):** entry 47's saved `logs/geo38/c11-cu-s64-d1-128mhz.npz` (`20261027T070008Z_GEO-38-C11.log`) — (RR3) reciprocity 8.858135e-15 and σ_max 0.999805773843 string-equal to entry 47's prints (`:56`); (CT3) 34 / 34 class lines equal `C11.log:14750–14783` (`:57`). The circuit layer's k = 1 block: anti 1.000000, off-block 1.168811e-03 (`:76`).
- **Rule that fired:** entry 47's decision rule as annotated by the 2026-10-26 weekly (§10 rulings; the retired 2026-10-24 entry below): "a known-issues row only if the absolute spread exceeds the adjacent class's". On the 1.5 T twin the same statistic read 0.8883 (C13).
- **Cause:** not diagnosed. Whether this is the mesh or the fixture (leg-box corner, shield can) is not measured; C15 claims neither.
- **Consequence:** no band moves; C16 / `xl` entry 49 still follows (the brief: neither outcome is a red). **Owner: the weekly review** (the 10-28 session never ran; written 2026-10-31).
- **Retirement condition (weekly 2026-10-31, §10 rulings):** the first reading is `xxl` entry 54 (the 1.5 T twin's 64×64 at degree 2 on the same mesh family, Saturday 11-07), read by C18 beside entry 44's degree-1 classes. **If the 1.5 T classes' relative spreads fall by ≥ 2× from degree 1 to degree 2 on the same mesh**, the 3 T reading is discretisation at 128 MHz (the same mesh is coarser in λ) and this entry **re-files as a record**. **If they do not**, the asymmetry is the fixture's (leg-box corner, shield can) and a 3 T degree-2 matrix (§10 C5e, then a priced window) is the direct reading before any band or geometry moves. Reciprocity (8.9e-15) and the copper identity hold either way, so no gate is touched and entry 49 runs as queued.
- **h-reading commissioned (weekly 2026-11-09):** the degree-1 64×64 on the r2 rung (`xl` entry 65, Monday 11-16, C35's knob on C3's module; C51 reads the class table beside r1's). The asymmetric class's spread on r2 ≥ 2× r1's 0.9262 % ⇒ the mesh's and this entry retires with the h-statement; within 2× ⇒ the fixture's and a geometry step (the leg-box corner first) is posed by the 11-16 weekly. No band moves.
- **h-reading (C51, §9 item 1003, 04:30 CST implementer slot 2026-11-16; `20261116T103321Z_GEO-38-C51-W1.log:94–166`):** entry 65's r2 matrix (912 804 cells; custody green — (RR) 1.573980e-14 / 0.999976063007 and (CT) 34 / 34 by value against `20261116T080008Z_GEO-38-C35w.log`) beside entry 44's r1. r2 is the **coarser** rung, so "≥ 2× with h" reads r2 / r1 ≥ 2. The asymmetric class (other ring, 1 step) reads **1.0436 % on r2 vs 0.9262 % on r1, r2 / r1 = 1.1268** (`:113`, `:129`); over 34 classes the median r2 / r1 is 1.5010, min 1.0597 (other ring, 2 steps), max 2.3375 (same ring, 16 steps), **3 classes at ≥ 2×** (same ring 2 / 15 / 16; `:130`). The C13 abs-spread statistic (other ring 1 step / same ring 1 step) reads 0.8883 (r1) → 0.9360 (r2) (`:166`). Printed only; the branch (retire with the h-statement, or pose the leg-box-corner step) is the 11-16 weekly's. Entry stays OPEN.
- **Order reading (weekly 2026-11-07, on `xxl` entry 54 read by C18, `20261107T182319Z_GEO-38-C18-W1.log:81–116`): the clause's second branch.** On the 1.5 T twin the relative spreads do **not** fall by ≥ 2× from degree 1 to degree 2 on the same mesh: the asymmetric class (other ring, 1 step) falls 0.9262 % → 0.6424 % (×1.44), the median class ×1.36 (min 0.57, max 2.07 over 34 classes; 1 class at ≥ 2×), and the C13 statistic (other ring 1 step / same ring 1 step) *rises* 0.8883 → 1.1578 at degree 2. So the asymmetry is not degree-1 discretisation; it is the fixture's (leg-box corner, shield can) or the mesh's. **Not re-filed; stays OPEN.** The 3 T degree-2 matrix is *not* the next reading (it would read the same class at the same order on the same mesh family); the next reading is **§10 C35** — the degree-1 64×64 on the r2 mesh rung (912 804 cells, `GEO-38` step 1a's second rung), where a fall with h says the mesh and no fall says the fixture. No band moves; reciprocity and the copper identity hold at both orders.
- **Verified at:** `94e22033` + the C15 module (the commit that lands this entry).

## 2026-11-14 — 🟡 OPEN (`GEO-34` step 1d / U12′, §9 item 985, 12:00 CST implementer slot; worked around in `src/`, not fixed; filed by the 18:00 daily review) — gmsh's Netgen optimiser reports "456 illegal tets" and the rank dies with SEGV on a mesh with an `occ.fragment`-embedded named face

- **Test id:** none red. The first W1 window of `tests/mesh/test_step_import_cq_port_sheet.py` at `-n 1` died: `20261114T181903Z_GEO-34-U12p-W1.log:111–112` "Remove Illegal Elements" / "456 illegal tets", `:116` "Caught signal number 11 SEGV", `:125` Status 15. It ran on the CadQuery `PortSheet` box after `generate(3)`, inside `gmsh.model.mesh.optimize("Netgen")`.
- **Workaround on `main`:** `src/fem_em_solver/io/step_import.py:320–330`. With the `surface_groups` key set, the Netgen pass is skipped and only `generate(3)`'s default tet optimisation runs. With the key absent the call is unchanged. The second W1 (`20261114T181927Z_GEO-34-U12p-W1.log`) is green at `-n 1 / 2 / 4`, and W2a / W2b re-ran the key-absent gates green.
- **Cause:** not diagnosed. The slot's comment ("Netgen's optimiser does not respect an embedded face") is a hypothesis; nobody has run the optimiser on an embedded face that came from another route.
- **Scope measured by `GEO-34` step 1e (§9 item 996, 16:30 CST 11-15, `mesh-probe`, gmsh only, no `import_step`; added by the 18:00 daily review):** on the 1d fixture the crash is not specific to the STEP route. A face embedded in **one** volume (`getEmbedded(3,1)` non-empty) dies with Status 139 at every h of 0.16 / 0.12 / 0.08 / 0.06 / 0.04 (`20261115T223140Z_GEO-34-1e-E-h0p16.log` … `20261115T223220Z_GEO-34-1e-E-h0p04.log`); (E, 0.08) repeats bit-identical three times (10316 tets, 456 illegal tets; `20261115T223203Z_GEO-34-1e-E-h0p08.log`, `…223323Z_…-r2.log`, `…223331Z_…-r3.log`). The box alone (N) and a face that is the **shared boundary of two volumes** (S) never crash (`20261115T223231Z_GEO-34-1e-N-h0p16.log` … `20261115T223434Z_GEO-34-1e-S-h0p04.log`). Pre-optimise minSICN > 0 everywhere, and every arm prints "illegal tets", so the count is not the signature; only E shows a first-pass `badmax` of 1e+24 or more. Readings for the Cause cell, not a diagnosis. The workaround stays as it is.
- **Scope:** every named-sheet STEP import, which is the route U13 (`WF-16` step 3, a `Study` on a STEP-imported F-small) would take, gets a mesh without the Netgen pass. Records measured on an optimised mesh (item 192's `P1`, `GEO-34` steps 1a / 1b) are therefore not reproducible through this route by construction. The mesh's quality without the pass has not been measured.
- **Resolves with:** the weekly's re-pose of U13 (2026-11-14 18:00 journal §5). Either a measured quality reading on the un-optimised named-sheet mesh with the step's gates green, or an optimiser route that keeps the embedded face. Never by re-enabling the pass on the key-absent path.
- **Verified at:** `b8214906` (U12′'s commit), from the two W1 logs; the review re-ran nothing.

## Retired entries — full text in `known-issues-archive.md`

- ✅ RETIRED 2026-11-15 (`WF-17` §9 item 997, 19:30 CST implementer slot) — was 🟡 OPEN 2026-11-15 (`WF-17` U163, item 992) — the 1.5 T estimate-page custody nodes raised `FileNotFoundError` on two logs the weekly sweep had gzipped
  - **Resolution:** the five plain-text readers (`_custody`, `_parse_c52`, the closed-form page's `read_one` / `is_file` / two `splitlines` reads, the lattice `_c39`) read through `tests/log_records.read_log_text(errors="strict")`. The pre-edit W0 measured 6 red per rank: the two 15t nodes **and the four p80 / p81 `_parse_c52` nodes, red since the 11-15 sweep with no entry until the 18:00 review inferred them; found and fixed here** (`20261116T013050Z_WF-17-i997-W0.log`). W1 55 passed per rank, 16 s, with the new gz-vs-plain identity module `tests/package/test_wf17_log_gz_readers.py` (`20261116T013304Z_WF-17-i997-W1.log`). The eight nodes due at the 11-22 sweep are covered by that module's opposite-form copies. Exempting cited logs from gzip remains the weekly's call. Full text in `known-issues-archive.md`.

- ✅ RETIRED 2026-10-21 (`POST-8` step 3b-ii-e3b, §9 item 580, 06:00 CDT implementer slot) — was 🟡 OPEN 2026-10-19 (`POST-8` step 3b-ii-e1, item 547) — `th:4` red at `-n 4`, its CG1-export profile reading 2.25 % RMS against `PROFILE_RMS_MAX`
  - **Resolution:** the profile now samples the DG split and the cross-check the raw N1curl `fields.e_complex`, both under the owner rule; bands unmoved. `th:4` green at `-n 1 / 2 / 4` with profile 1.591 % at every width, cross-check ≤ 1.887e-16, (S-sb) 0 at every width, (O) call count / line / `n_inc` / owner indices exact (`20261021T110150Z_POST-8-step3b-ii-e3b-th4-n2.log`, `…110229Z_…th4-n4.log`, `…110303Z_…th4-n1.log`). The CG1 export route remains partition-dependent and is printed, labelled, not asserted: 2.249 % at `-n 4` (`cells_ghost` [559, 652, 520, 670]), 0.124 % at `-n 2` / `-n 1`. Full text in `known-issues-archive.md`.

- ✅ RETIRED 2026-10-12 (`WF-1` step 9, §9 item 401, 12:00 CDT implementer slot) — was 🟡 OPEN 2026-10-12 (found by `WF-1` step 7a, item 390a) — `femem check` printed `family: four_port_birdcage` for an `f_human_birdcage` YAML file
  - **Fix (`cli.py` `check` branch only):** it prints the parsed `kind` from `study_yaml._load` (the private name, imported inside the package; no wrapper added). The other five lines, return codes and the `StudyConfigError` path are unchanged. `tests/package/test_wf1_check_family.py` (`20261012T172853Z_WF-1-s9-W1.log`, 6 passed, `-n 2`): (F) on 11 tracked files the printed family equals `_load(f)[1]`; (F-h) F-human `f_human_birdcage`; (U) the other five lines equal the pinned `b01dc790` `cli.py`; (N) the logged pre-change line differs. Re-gate 13 passed (`20261012T172909Z_WF-1-s9-RG.log`); installed `femem` 1 passed (`20261012T172921Z_WF-1-s9-W9-installed-femem.log`).

- ✅ RETIRED 2026-10-12 (`ANS-7` lane 1, §9 item 400, 12:00 CDT implementer slot) — was 🟡 OPEN 2026-10-12 (found by `ANS-7` step 13b, item 392) — a later `XT64R` / `X4` / `XA` run of `ans:7` would drop the 127.74 MHz ESR section from `COMPARISON.md`, because `_write_comparison` rewrote the file wholesale and only `XT128R`'s `_write_comparison_esr128` inserted that section
  - **Fix (render path only):** `_write_comparison` ends by calling `_write_comparison_esr128(m)` when `tuned_larmor_esr` holds `"127.74"`; nothing else in the script changed. `tests/unit/test_ans7_comparison_writer.py`, render-only on tmp copies: (R) the committed metrics render the committed `COMPARISON.md` byte for byte (3543 B, heading count 1); (N-old) the pinned `b01dc790` writer gives 2886 B, count 0; (N-key) without the key the new writer equals (N-old) byte for byte; (R0) repo files untouched. `20261012T170114Z_ANS-7-item400-W1.log`, 5 passed. No `ans:7` stage re-run.

- ✅ RETIRED 2026-10-12 (`OPS-66` step 3, §9 item 393, 09:00 CDT implementer slot) — was 🟡 OPEN 2026-10-12 (found by `OPS-66` step 2, item 384) — `test_rotate_known_issues_entries.py::test_a2_old_subset_of_new_and_difference_is_exactly_the_four` was red on clean `main`: its `archive ⊆ FOUR` clause asserted a snapshot of the archive, and later rotations moved three more closed undated RETIRED headings there (`20261012T035840Z_OPS-66-s2-W0-main.log`)
  - **Fix (test only):** the snapshot clause is replaced by "the four ⊆ live ∪ archive new − old". The closed-only check on both files and the live `== FOUR` clause are unchanged. Controls, both asserted: the archive's new − old beyond the four is non-empty (3, so the old clause would still fail), and a planted `## RETIRED 2026-10-12 PARTIAL planted` heading (admitted, not closed) is caught by the closed check. The brief's `RETIRED-NOT` form is not admitted by `KI_ENTRY_RE` (printed). `20261012T141002Z_OPS-66-s3-W1.log`, 14 passed.

- ✅ RETIRED 2026-10-11 (`WF-16` step 11, §9 item 388, 22:30 CDT implementer slot) — was 🟡 OPEN 2026-10-11 — `Study.sar` / `b1_plus` on a `FourPortBirdcageGeometry(as_hole=True)` Study raised `material_map … [1]`, and on `as_hole=False` with `LeontovichWall()` silently returned the σ = 800 fields
  - **Faces:** (F1) `_study_solve_port_fields` mapped conductor tag 1 unconditionally, so on the hole mesh (no tag-1 cells) `sar` / `b1_plus` / `b1_plus_at` raised `ValueError: material_map references tags that do not exist in problem.cell_tags: [1]. Known tags: [2, 3, 101, 102, 103, 104, 201, 202, 203, 204]` (`20261012T033148Z_WF-16-s11-W0r.log:942`, unmodified `src/` at `e122cead`); (F2) it never called `_conductor_setup`, so a solid-route `LeontovichWall` Study was not refused — on a no-mesh stub it reached `built["cell_tags"]` and raised `KeyError: 'cell_tags'` (same log :1190) instead of `solve()`'s `ValueError` ("as_hole").
  - **Fix:** `_study_solve_ports`' conductor lines copied, in order, into `_study_solve_port_fields` (no shared helper). W1 (`20261012T033238Z_WF-16-s11-W1.log`, 12 passed): (S=) field route vs `solve()` on the same Cu Study 7.820e-16 (band 1e-9); (P) kept Cu P1 field residual 8.505e-13 (rtol 1e-6), P_surf 1.516873149e-06 W; (N2) P_surf-blind residual 3.1559e-04 = 315.6× the rtol. Default route re-gated: RG1 / RG2 (R) 9.557e-11.

- ✅ RETIRED 2026-10-08 (`WF-17` API pages, §9 item 331, 12:00 implementer slot) — was 🟡 OPEN 2026-10-08 (`PORT-24` step 1, item 308) — the committed API reference was stale for `load_study` / `run_study_file`; regenerated, `test_wf17_api_reference.py` 12 passed, guide links 0 dead
- ✅ RETIRED 2026-10-06 (`ANS-7` step 7′, §9 item 262, 06:00 implementer slot) — was 🟡 OPEN 2026-10-04 (`ANS-7` step 7, item 217) — `test_t21_t2_inmodel_zero` red by pre-registration at C_mode; the in-model zero at C_in reads −2.536e-6, inside `T21P_ZERO_RESOLUTION` 1.287e-5
- RETIRED 2026-09-19 — Fourteen committed setup figures carried legends clipped mid-word — a false-artefact mode the setup-figure census cannot see (2026-09-19, 0…
- ✅ RETIRED 2026-09-13 (`TH-15` step 3c, 19:30 implementer slot) — was OPEN 2026-09-13 (`TH-15` step 3, 15:00 implementer slot; filed by the 18:00 review on a `lo…
- ✅ RETIRED 2026-09-12 (`PORT-19` step 5, 06:00 implementer slot) — `tests/validation/test_port19_factor_reuse.py::test_a_kept_fields_are_per_drive_and_match` was…
- ✅ RETIRED 2026-09-11 (`ANS-4` step 2a″, 04:30 implementer slot — the guard landed and its control passed) — ~~🔴 OPEN 2026-09-09 (`ANS-4` step 2a, 04:30 implemen…
- ✅ RETIRED 2026-09-08 (`OPS-41`, 07:30 implementer slot) — ~~🟡 OPEN 2026-09-07 (`TH-15` step 2d, 12:00 implementer slot) — `tests/validation/test_port_package_sp…
- ✅ RETIRED 2026-09-13 (`ANS-4` Larmor verdict, weekly review — banked by the operator's interactive session after the scheduled session died on the account limit…
- ✅ RETIRED 2026-09-06 (`OPS-39`, 15:00 implementer slot) — ~~🟡 OPEN 2026-09-06 (`TH-15` step 2a, 00:00 implementer slot; filed by the 10:30 review) — the step-0 …
- ✅ RETIRED 2026-09-11 (`PORT-14` step 1e, 07:30 implementer slot) — systematic registered as a record — ~~🟡 OPEN 2026-09-05 (`PORT-14` step 1, 21:00 implementer …
- ✅ RETIRED 2026-09-07 (`PORT-16` step 2, 06:00 implementer slot) — ~~🟡 OPEN 2026-09-04 (`POST-6` step 1, 19:30 implementer slot) — the **drive-level** power iden…
- ✅ RETIRED 2026-09-03 (`GEO-26` step 3, 22:30 implementer slot) — ~~🔴 OPEN 2026-09-03 (`GEO-26` step 2, 12:00 implementer slot) — the **longitudinal** ring-gap s…
- ✅ RETIRED 2026-09-03 (`OPS-33`, 00:00 implementer slot) — ~~🟡 OPEN 2026-09-02 (`OPS-32` slot) — `ans:3`'s four reproduction controls cannot be re-registered at …
- ✅ RETIRED 2026-09-02 (`WF-6` step 3h, 19:30 implementer slot) — ~~🔴 OPEN 2026-08-31 (`WF-6` step 3, 13:30 implementer slot) — the coil-driven **point-SAR** map …
- ✅ RETIRED 2026-09-03 (`WF-6` step 3f′, 21:00 implementer slot) — ~~🔴 OPEN 2026-09-02 (`WF-6` step 3f, 15:00 implementer slot) — the gated `|B₁⁺|` **C4 identitie…
- ✅ RETIRED 2026-08-27 (`OPS-27` step 2, 21:00 implementer slot) — a **tenth** stale exact cell-count record, and it is a **fifth** mesh with **no sibling**: `tes…
- ✅ RETIRED 2026-08-27 (`OPS-27` steps 1 and 2) — **five more stale exact cell-count records across five modules, and the class collapses to THREE shared meshes, …
- ✅ RETIRED 2026-08-28 (`OPS-27` step 3, 16:30 implementer slot) — `test_coil_loading_larmor_third_rung.py` asserted an **exact** fine-rung cell count recorded on…
- ✅ RETIRED 2026-08-27 (`OPS-27` step 1, 19:30 implementer slot) — `test_coil_loading_larmor_mesh_cache.py` asserted an **exact** third-rung cell count recorded o…
- ✅ RETIRED 2026-08-27 (`OPS-27` step 1, 19:30 implementer slot) — `test_geometry_floor_discriminator.py` asserted the **pre-`OPS-18`** 128 MHz record (1.8260%) a…
- ✅ RETIRED 2026-08-28 (`OPS-28`, 22:30 implementer slot) — the whole of `tests/ports/test_port_orientation_sensitivity.py` dies on `'_DummyComm' object has no at…
- ✅ RETIRED 2026-08-28 (`GEO-23` step 2b, 13:30 implementer slot) — ~~🔴 OPEN 2026-08-27, re-headed 2026-08-28 (`GEO-23` step 1) — a **fifth** site of "Invalid bou…
- ✅ RETIRED 2026-08-28 (`GEO-23` step 2b, 13:30 implementer slot) — ~~🔴 OPEN 2026-08-27 — `test_boundary_condition_selection.py` **deadlocks the whole command** a…
- ✅ RETIRED 2026-08-28 (`GEO-23` step 2b, 13:30 implementer slot) — ~~🔴 OPEN 2026-08-27 — `test_phantom_field_metrics_and_exports_are_finite` aborts in gmsh with …
- ✅ RETIRED 2026-08-28 (`GEO-23` step 1 (d), 09:00 slot) — ~~DEAD MODULE, filed 2026-08-27 (`OPS-26` step 2 leg (a)) — `tests/mesh/test_cylindrical_domain.py` col…
- ✅ RETIRED 2026-08-25 (`EX-30` item 3 half A, 22:30 implementer slot) — `test_kwarg_off_reproduces_the_recorded_mesh` was **red on `main`**: the `GEO-16` kwarg-o…
- ✅ RETIRED 2026-08-25 (`EX-30` item 3 half B, 22:30 implementer slot) — `mesh:5`'s **inverted control lost its separation**: the clamps-only mesh *cleared* the 0…
- ✅ CLOSED 2026-09-01 (`OPS-30`, 21:00 implementer slot) — two `scripts/probes/` scripts were **never migrated to dolfinx 0.11**: they construct `fem.petsc.Linear…
- ✅ CLOSED 2026-08-25 (`MAG-19` step 2, 21:00 implementer slot) — `tests/validation/test_convergence.py::TestConvergence::test_h_refinement_straight_wire` was **r…
- ✅ RETIRED 2026-08-30 (weekly planning review — size-field licence **denied**, the `GEO-23` step-2a wrap ruled sufficient) — ~~`MeshGenerator.straight_wire_domai…
- ✅ RETIRED 2026-08-25 by `OPS-24` — `core/cavity.py` was **never migrated to dolfinx 0.11**: `assemble_matrix(..., diagonal=)` no longer exists, so the whole `TH…
- ✅ RETIRED 2026-08-25 by `OPS-25` — `th:7` calls `Function.interpolate(cells=)`, removed in 0.11 — the **only** such site in the repo, so the example has diverge…
- ✅ RETIRED 2026-08-25 by `EX-30` leg (th) — `th:6`'s **128 MHz** interior relL2 does not reproduce the `TH-10` record on the 0.11 image (1.76864% vs 1.826%, 3.14…
- ✅ RETIRED 2026-08-25 (`GEO-19` step C, 07:30 implementer slot) — `birdcage_port_domain(emit_port_sheets=True)` **cannot build any birdcage with more than four l…
- ✅ RETIRED 2026-08-25 (`GEO-19` step C, 12:00 implementer slot) — `GEO-18` step 1's **1e-5 terminal-equality band is a C4 band, not a C_N one**: at 16 legs the t…
- ✅ RETIRED 2026-08-24 (`PORT-9` leg (d3c), 12:00 implementer slot) — two birdcage **reproduction controls** were red on `main`: the 0.11 image meshes the gapped …
- ✅ RETIRED 2026-08-23 (`OPS-18` step 3b merge) — the two-torus port fixtures **SIGABRT'd in `gmsh.model.mesh.generate` only in the 0.11 image**: numpy 2 renders …
- ✅ RETIRED 2026-08-23 (`OPS-18` step 3b merge) — two two-torus **reproduction records** moved in the 0.11 image at 1e-4 while every physics identity held (`OPS-1…
- ✅ RETIRED 2026-08-23 (`MAG-18` re-gate, 19:30 implementer slot) — `test_straight_wire_b_field` failed **only in the 0.11 image**: the discretization error moved…
- ✅ RESOLVED 2026-08-23 (`OPS-18` step 3b) — `test_region_resolution_policy_refines_the_tagged_volumes_toward_cad` failed **only in the 0.11 image**: the uniform-…
- ✅ FIXED 2026-08-19 (`OPS-22`) — the magnetostatic loop-drive fixtures were complex-hostile: `ufl.max_value` / `<=` geometry predicates in the current-density ca…
- 1. ✅ RETIRED 2026-07-31 — stale test double, `DummyMagnetostaticSolver`
- 2. ~~Residual-trend classifier disagrees with its test~~ — RESOLVED 2026-08-08 (`OPS-12`)
- 4. ✅ RETIRED 2026-08-08 — coil+phantom B-field symmetry exceeded tolerance (`MAG-6` step 3)
- 5. ✅ RETIRED 2026-08-06 — domain sizing heuristic, off-centre phantom
- 8. ✅ RETIRED 2026-08-05 — magnetostatic energy raised `TypeError` in the complex build
- 9. ✅ RETIRED 2026-08-05 — `-n 2` hang on `two_torus_domain`'s port facets
- 10. ✅ RETIRED 2026-08-06 — `two_torus_domain`'s outer-boundary facet group never reached the dolfinx facet tags (`GEO-10`)
- 12. ✅ RETIRED 2026-08-06 — `loop_over_half_space_domain` and `sphere_in_box_domain` never declared their `outer_boundary` group (`GEO-12`)
- 13. ✅ RETIRED 2026-08-07 — `cylindrical_domain`'s classification margin was 4.50× its tolerance (`GEO-13`)
- 7. ✅ RETIRED 2026-08-03 — birdcage mesh fails to generate (`GEO-9`, steps 1 + 2a + 2b)
- ~~`check_example_doc_references.py` freshness-gates only 5 of 27 examples — every `stale=24, none of them mine` line is not an all-clear~~ — RESOLVED 2026-08-24…
- ~~`test_coil_phantom_magnetostatics` fails in the complex build on a cold FFCx cache: `ComplexComparisonError`~~ **RESOLVED 2026-08-19 (`OPS-20`, 06:00 implemen…
- ✅ RETIRED 2026-09-01 — `test_coil_loading_degree2.py` no longer returns inside its 570 s ceiling (exit 124 at 571 s, **twice** on 2026-08-31)
- ✅ RETIRED 2026-08-11 — "unexplained" mid-command termination of the logging harness was the background-and-end-turn trap (2026-08-08, 15:00 and 19:30 implemente…
- ✅ RETIRED 2026-08-09 — "rank-dependent DG0 centerline sample" was the probe's own bug (`MAG-6` step 4, second pass)
- ✅ RETIRED 2026-08-04 — reaction Z-matrix diagonal is negative where it must be inductive
- ✅ RESOLVED 2026-08-03 — "birdcage suite is over the compute budget" was the hang, not meshing cost
- ✅ RESOLVED 2026-07-30 — truncation-wall modeling floor on the wire/loop fixtures (MAG-13)
- ✅ RESOLVED 2026-08-01 — `two_torus_domain` was not a conforming mesh (`GEO-8`)
- ✅ RETIRED 2026-08-10 (`EX-17`) — `02_circular_loop.py` never wrote its VTX/`.bp` output
- ✅ RETIRED 2026-08-11 (`POST-4` step 3) — `examples/mri/01` centerline samples were rank-dependent at ≈23%, at and below the gauge floor
- ✅ RESOLVED 2026-08-23 (`OPS-18` step 3a `5df1e39`, landed on `main` by 3b) — the two-torus and straight-wire solves are run-to-run non-deterministic at ≈1e-10 r…
- ✅ RETIRED 2026-08-25 — Gate (iii) is blind to a broken C4 on the *opposite* class, and the lumped-sheet 4-port sweep loses reciprocity by 223× on an asymmetric …
- ✅ RETIRED 2026-08-29 (`GEO-24` step 2b, 06:00 implementer slot) — ~~the **32-port** ring-gap sheet reconstruction is rank-width dependent~~ ~~**`birdcage_port_d…
- ✅ RETIRED 2026-08-30 (`PORT-12` step 2, 13:30 implementer slot) — ~~the two-torus gap-route reproduction record drifts with rank width on a fixture that already…
- ✅ RETIRED 2026-08-30 (`WF-6` step 1d, 12:00 implementer slot) — ~~the first `|B₁⁺|` map is C4-covariant to only ≈9%, against a 5% pre-registered discretisation …
- ✅ RETIRED 2026-09-07 (`TH-15` step 2d, 12:00 implementer slot) — ~~The gap-displacement current and the conduction current are different quantities on an **undr…
- ℹ️ RESOLVED-IN-CODE 2026-09-22 (`OPS-62`) — the `ANS-4` step 3e XL window exited 1 on a spurious guard; its four rungs' data are intact
- ✅ RETIRED 2026-09-30 (§9 item 119, 19:30 CDT slot; opened 2026-09-30 by item 118) — `import_step`'s `outer_boundary_tag` created no facet group on the F-small S…
- ✅ RETIRED 2026-09-29 (`MAT-1` §9 item 83, 19:30 CDT slot; opened 2026-09-29 by item 81) — `tests/materials/test_phantom_material_model.py::test_gelled_saline_pr…
- ✅ RETIRED 2026-09-29 (opened and retired in the same commit, `PORT-7`, 22:30 CDT slot, §9 item 64) — `export_touchstone` wrote N ≥ 3 S-matrices column-major on …
- ✅ RETIRED 2026-09-29 (06:00 slot, §9 item 74) — was OPEN 2026-09-29 (filed by `TH-16` step 1, 22:30 CDT slot, §9 item 65; found while drafting by the 2026-09-28…
- ✅ RETIRED 2026-09-30 (18:00 daily review, on §9 item 112's six logs and a `log-pathologist` ruling: located and mitigated, not fixed — the drift needs ≥ 2 ranks…
- ✅ RETIRED 2026-09-30 (`PORT-5` ✅, 04:30 slot, §9 item 98) — 3. Port tests assert a non-zero power wave where the placeholder's fake makes it exactly zero (re-da…
- ✅ RETIRED 2026-09-29 (`PORT-4` ✅, 12:00 slot, §9 item 77) — 6. Rank-dependent: single-port excitation
- ✅ RETIRED 2026-09-30 (04:30 slot, §9 item 96: the balance is now dimensionless against an optional `reference_speed_m_per_s`, `None` with no warning when none i…
- ✅ RETIRED 2026-09-30 (04:30 slot, §9 item 95: tolerance floored at 4 ulp of the carrier, last point pinned to `stop`, steps under 16 ulp raise; probe `duplicate…
- ✅ RETIRED 2026-09-30 (13:30 slot, §9 item 111: the break decision goes through `_ladder_budget_allows`, an `MPI.MIN` allreduce of `remaining − 2·last_cost`; `te…
- ✅ RETIRED 2026-10-03 06:00 slot (§9 item 197 green re-run on the rule (f) band) — was OPEN 2026-10-02 (`TH-20` step 0a, §9 item 176 W1, landed by §9 item 182) —…
- ✅ RETIRED 2026-10-02 13:30 slot (§9 item 175 green re-run, κ off per the operator ruling) — was OPEN 2026-10-02 (`WF-12` step 2′ W1, 12:00 implementer slot, §9 …
- ✅ RETIRED 2026-10-01 21:00 slot (§9 item 158, on the 18:00 review's rule (f) ruling) — was OPEN 2026-10-01 (`WF-12` step 1 W1, 06:00 implementer slot, §9 item 1…
- ✅ RETIRED 2026-10-03 (22:30 CDT slot, §9 item 186; `GEO-1`) — was: 🔴 OPEN 2026-08-27 (`OPS-26` step 2 leg (a), second slot) — `test_birdcage_volumes_partition_t…
- RETIRED 2026-09-20 (`OPS-59` (b), 07:30 slot) — `tests/unit/test_setup_figure_title.py::` the call-site count pin (`EXPECTED_EXAMPLE_CALL_SITES = 18`) was stale…
- ✅ RETIRED 2026-09-21 (`OPS-59` (a), 06:00 implementer slot) — ~~2026-09-19 — `tests/unit/test_doc_reference_exit_codes.py::test_the_in_tree_exemption_cannot_sil…
- ✅ RETIRED 2026-09-20 (`PORT-20` steps 1–3, 04:30 + 09:00 implementer slots) — ~~🟡 OPEN 2026-09-19 (`ANS-3` adjudication, weekly review) — on the gap-voltage / i…
- ✅ RETIRED 2026-09-27 (03:00 daily review; fixed by `4719e82`, 2026-09-26 09:00 slot) — was: 2026-09-20 — 🟡 OPEN (`PORT-20` step 3, 09:00 implementer slot) — `an…
- ✅ RETIRED 2026-10-05 (`POST-10` step 1, 15:00 CDT slot, §9 item 254) — was 🟡 OPEN 2026-10-04 (`ANS-9` step 1a) — `evaluate_vector_field_parallel(owner="lowest_original_cell")` raised `Mixed topology unsupported.` on a `create_submesh` mesh
- ✅ RETIRED 2026-10-04 as moot (03:00 daily review, on the 2026-10-03 weekly's "retire it as moot or fix the test, your call") — was 🟡 OPEN 2026-10-03 (`xxl` wind…
- ✅ RETIRED 2026-10-07 (`OPS-68`) — was `MAG-18` anchor (ii) is unreachable as pre-registered: the magnetostatic solve's own cross-width floor is ≈1e-7, not 1e-10…
- ✅ RETIRED 2026-10-07 (`OPS-67`) — was 🟡 OPEN 2026-09-03 (`PORT-13` step 1, 22:30 implementer slot) — `GEO-20`'s ring-gap port sheets are **transverse** sections…
- ✅ RETIRED 2026-10-07 (`OPS-68`) — was `check_example_doc_references.py` exit codes: 0 clean, 1 real defect, 2 staleness only (`OPS-19`, 2026-08-16 — contract, n…
- ✅ RETIRED 2026-10-07 (`OPS-67`; resolved 2026-09-01 by `TH-13` step 3a‴) — `test_coil_loading_degree2.py` no longer returns inside its 570 s ceiling (exit 124 a…
- ✅ RETIRED 2026-10-07 (`OPS-68`) — was Four defects surfaced by replacing finiteness-only tests with real anchors (`OPS-17` step 2, 2026-08-17)
- ✅ RETIRED 2026-10-07 (`OPS-67`; repaired 2026-08-28 by `OPS-28`) — `PORT-1` step 4's `allgather` reduction broke `_DummyComm`, and one orientation test regresse…
- ✅ RETIRED 2026-10-07 (`OPS-68`) — was `EX-49`'s negative-control ceiling-first floor (≥2× the band) was a plan-time estimate that the measured fixture falls 1.2…
- 2026-09-19 — ✅ RETIRED 2026-10-05 (weekly review, T6 decided on F18a / F18b) — was 🟡 OPEN (`ANS-4` step 3 / 3c private readouts, weekly review) — the order-matc…
- ✅ RETIRED 2026-10-09 (`OPS-70` fix, §9 item 347) — was 🔴 OPEN 2026-10-08 (`OPS-70` step 1 W7, 21:00 implementer slot, §9 item 335): `OPS-70: tests/package/test_…
- ✅ RETIRED 2026-10-08 (OPS-69, §9 item 336) — was 🟡 OPEN 2026-10-08 (found by `EX-89`, §9 item 332, 12:00 slot; not fixed in passing) — the examples docrefs cens…
- 2026-10-10 — ✅ RETIRED 2026-10-14 (weekly review; `ANS-9` re-adjudicated on the degree-2 exact-rectangle column, every row AGREE, probe D within 2.5 %) — was 🟡 …
- 2026-10-17 — ✅ RETIRED 2026-10-17 (opened and retired in the same commit, `WF-13` step 9, §9 item 507, 04:30 CDT implementer slot) — `femem run --estimate` on a…
- 2026-09-30 — ✅ RETIRED 2026-10-17 (weekly review, T41: the `f₀` clause was degree-1 discretisation — T13, 2026-10-10, within 0.11 % / 0.03 % at degree 2 on exac…
- 2026-10-13 — ✅ RETIRED 2026-10-17 (weekly review: `IMP-6`'s step-2 family is frozen on the row after I29 — the rim truncation below ρ_r/3 stands as a documented…
- ✅ RETIRED 2026-10-22 (`POST-8` 3b-iii, §9 item 587; retirement row below) — was 🟡 OPEN 2026-09-27 (03:00 daily review, from the 2026-09-26 09:00 slot's disclosu…
- ✅ RETIRED 2026-10-24 as a record (weekly review: `TH-15` step 2 was killed 2026-10-07 and `TH-15` closed ✅ 2026-10-19, so "stays open until step 2 closes" can n…
- ✅ RETIRED 2026-10-24 (weekly review, on `xl` entry 43's first branch; was 🟡 OPEN 2026-10-14 — `IMP-3` step 7 / I15 part 3, `xl` entry 33, 2026-10-14 02:00; rule…
- ✅ RETIRED 2026-10-27 (`OPS-72`, §9 item 665, 07:30 CDT implementer slot — full `tests/environment` 44.90 s at `-n 1`, 44.73 s at `-n 2`; the environment-first r…
- ✅ RETIRED 2026-10-26 (weekly review, on the 10-24 condition's first clause — H33's width law at both orders; re-filed in PROJECT_PLAN §2.2 (4) as a record; H34 …
- ✅ RETIRED 2026-10-26 (weekly review, by the entry's own first branch — C13 read the class's absolute spread at 0.888× the adjacent class's: a relative artefact …
- ✅ RETIRED 2026-11-05 (`OPS-81`, §9 item 828, 12:00 implementer slot — the scope docstrings and the guard landed; the positive control raised and the callers' ga…
- ✅ RETIRED 2026-11-09 (§9 item 918, 04:30 implementer slot — the premise tests read empty tmp log roots; three W0s green, the control red as registered) — was 🟡 …
- RETIRED 2026-11-10 (18:00 daily review's ruling, `WF-17.md` § "Review ruling (2026-11-10 18:00 daily": the print equals the node's own named record `WF15_DEG2_R…
- RETIRED 2026-11-10 (10:30 daily review's ruling, `WF-17.md` § "Review ruling (2026-11-10 10:30 daily": the moved value is a band-gated residual, not a record — …
- RETIRED 2026-11-10 (10:30 daily review's ruling, `WF-17.md` § "Review ruling (2026-11-10 10:30 daily": the moved value is a band-gated residual, not a record — …
- RETIRED 2026-11-10 (10:30 daily review's ruling, `WF-17.md` § "Review ruling (2026-11-10 10:30 daily": the moved value is a band-gated residual, not a record — …
- RETIRED 2026-11-14 (`ANS-11` C37, §9 item 983, 12:00 CST implementer slot — `test_c19_absent` retired in favour of `test_c37_custody`; the default render reads …
- RETIRED 2026-11-14 (`GEO-36` I85c, §9 item 978, 07:30 CST implementer slot — the literal re-recorded at 20 with `FEM_EM_GEO38_C21_DEGREE` asserted as a member; …
- RETIRED 2026-11-14 (`GEO-38` C33c, §9 item 979, 07:30 CST implementer slot — `:265` formats dz `+.6f`; (J58) 0 mismatches, P89′ 0 string-only) — was 🟡 OPEN 2026…
- RETIRED 2026-11-14 (`GEO-34` step 1d / U12′, §9 item 985, 12:00 CST implementer slot — the "Resolves with" clause's additive read option landed as `"highest_dim…
