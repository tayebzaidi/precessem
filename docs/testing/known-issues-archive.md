# Known-issues archive — retired entries moved out of known-issues.md

Verbatim RETIRED entries moved out of `docs/testing/known-issues.md` by
`scripts/maintenance/rotate_plan_archive.py known-issues` (first batch
2026-09-19), in their original order — never summarized, never edited.
known-issues.md answers "is this failure mine?", which only OPEN entries can;
a retired entry is history. Grep here when a retirement's reasoning, digits or
log lines are needed; the live file keeps one index line per entry.

---

## Batch moved 2026-09-19

### RETIRED 2026-09-19 — Fourteen committed setup figures carried legends clipped mid-word — a false-artefact mode the setup-figure census cannot see (2026-09-19, 03:00 review, from the 2026-09-18 04:30 slot's finding; narrowed to leg B's nine 2026-09-19, retired when leg B landed)

**Leg A (the five magnetostatics figures) is re-rendered and this entry
narrows to leg B's nine (`mesh:3`, `th:10`, `mri:1/2/3`, `ports:15–18`),
2026-09-19:** `mag:1/2/4/5/6` re-rendered flagged (`FEM_EM_SETUP_FIGURES=1`,
`-n 2`, real build, recorded width) at `1c12233` and read with the Read
tool; every entry's legend now shows whole (`mag:2`'s previously read
"wire (con" / "air (hidde" now reads "wire (conductor)" / "air (hidden)"
in full — the negative control). Every printed identity reproduced to the
digit against its original `EX-57` flagged window: `mag:1` mesh
21830/4662, relL2 51.9781%, max rel 76.7331%→76.7332%
(`20260914T125223Z_EX-57-straight-wire.log:271-272` vs
`20260919T111506Z_EX-57.log:271-272`) — the max-rel last-printed-digit
move is the same mesh-partition round-off class the original item already
noted against the guide's own table (76.7330%), not a regression, and
`mag:1` carries no `assert`; `mag:2` mesh 409596 cells, relL2 6.2134%, max
rel 11.6541%, energy 2.466102e-08 J, byte-for-byte
(`20260914T125551Z_EX-57-circular-loop.log:264,289-290,294` vs
`20260919T111554Z_EX-57.log:264,289-290,294`); `mag:4` cells
69918/103950/160677, centre B_z 3.563601e-09/3.519075e-09/3.483786e-09,
rel err 0.92/0.34/1.34%, mean/max 2.15/7.92, 1.03/4.64, 1.56/5.33%, CV
0.075/0.028/0.056% (`20260914T140147Z_EX-57-helmholtz.log:248-250,427-429,626-628`
vs `20260919T111830Z_EX-57.log:248-250,427-429,626-628`); `mag:5` probe
vector L2 0.0003%, volume L2 0.0040%, max|A| ratio 2.773e-11, multiplier
spread nan/2.083e+02, "All assertions hold"
(`20260916T095415Z_EX-57-mag5-flagged.log` vs
`20260919T112024Z_EX-57.log:219,220,235-237,255`); `mag:6` errors
21.8417%/15.3848%/4.4605%, rate 1.9038, exported fld 16.8915%
(`20260916T095753Z_EX-57-mag6-flagged.log` vs
`20260919T112054Z_EX-57.log:185,341,572,585-587,589,592`). Every caption's
existing prose (region names, slice plane, hidden air) checked against the
new image and still true; no caption or `.py` comment needed correction, so
no re-run beyond the one flagged window per example was triggered. Census
unchanged, `examples=54 ok=18 missing=36 broken=0`
(`20260919T112334Z_EX-57.log:89`); docrefs `dead=0 guide=0 stale=25 exit=2`
(`20260919T112335Z_EX-57.log:64`, ≠ 1). Leg B (the other nine) is untouched
and this entry stays open until it lands.

**Leg B (`mesh:3`, `th:10`, `mri:1/2/3`, `ports:15–18`) is re-rendered and
this entry is retired, 2026-09-19:** all nine re-rendered flagged
(`FEM_EM_SETUP_FIGURES=1`, `mesh:3` real / the rest complex, `-n 2` except
`mri:3` at its recorded `-n 4`) at `f8018fc` and read with the Read tool;
every legend now shows whole. Negative control, `ports:15`: `git show
f8018fc:examples/ports/figures/ports_15_birdcage_sixteen_leg_quadrature_b1_setup.png`
shows the legend truncated mid-word — `"ring port 17 … ring p"` — where the
re-rendered image reads `"ring port 17 … ring port 38"` in full. Every
printed identity reproduced to the digit (or to the same solver-round-off
class already on record) against each example's original `EX-57` flagged
window — `mri:1` `|E|` 1.979842e+02 / `|B|` 1.294602e-06
(`20260916T110133Z_EX-57.log` vs `20260919T112931Z_EX-57.log:438-439`);
`mri:2` closed-form/point/DG0/1g/10g digits byte-identical
(`20260916T123256Z_EX-57.log` vs `20260919T112951Z_EX-57.log:122-137`);
`mri:3` all four C4 pair spreads and mis-paired controls byte-identical,
only the ≤1e-10-budget identity residual moved 6.217e-15→5.107e-15 and the
step-3f comparison 5.300e-11→4.453460e-11, both round-off
(`20260916T140657Z_EX-57-mri3-flagged-r2.log` vs
`20260919T113410Z_EX-57.log:1858,1865-1868,1885`); `mesh:3` `GEO-15` graded
0.966977 / control 0.846150 / separation 0.120826 unmoved
(`20260919T113114Z_EX-57.log:2869-2875`); `th:10` every gate/scaling/copper
`Z_s`/`lambda`/`f`/`Q` line and `[a]`/`[b]`/`[c]` assertions byte-identical,
only the PEC-pencil control (4.798e-19→8.063e-20) and Rayleigh residual
(7.304e-14→7.473e-14) moved at round-off, both ≪ their 1e-10/1e-6 budgets
(`20260914T124948Z_EX-60.log` vs `20260919T113011Z_EX-57.log:51,54-55,59`);
`ports:15` (i)/(ii)/(iii)/n_valid byte-identical
(`20260914T095534Z_EX-55.log` vs `20260919T113830Z_EX-57.log:11641-11647`);
`ports:16` x1/x0.012 spreads and the monotone-fall assertion byte-identical
(`20260914T110437Z_EX-56.log` vs `20260919T113547Z_EX-57.log:1840,3615,3618`);
`ports:17` `C_tuned`, `Z_in`, `|S11|`, both residuals byte-identical
(`20260914T112025Z_EX-58.log` vs `20260919T113256Z_EX-57.log:41,43,1799-1800`);
`ports:18` gates/reciprocity/P_src/P_sheet/P_phantom/P_surf/ratios
byte-identical, only the ≤1e-3-budget reciprocity and ≤1e-6-budget residual
moved at round-off (`20260914T124230Z_EX-59.log` vs
`20260919T113208Z_EX-57.log:910,913-914,918`). Two guide captions' exact
byte-size parentheticals (`ports:16`, `ports:17`) drifted by under 3 KiB
from the longer un-clipped legend text and were corrected in the same
commit; every other caption sentence checked against the new image was
still true, so no other correction and no `.py`-comment re-run was
triggered. Census unchanged, `examples=54 ok=18 missing=36 broken=0`
(`20260919T114427Z_EX-57.log:89`); docrefs `dead=0 guide=0 exit=2`
(`20260919T114436Z_EX-57.log:62`, ≠ 1). Both legs have now landed; this
entry is retired.

**Half (b) of this entry — `write_setup_figure` having no title-length guard —
is fixed and retired by `OPS-53` (2026-09-19):** the helper now carries a
measured `MAX_TITLE_CHARS = 130` (bracket [126 ok, 150 clipped] at font size 11
on its fixed canvas) and raises `ValueError` from `check_title_length` as its
first statement, on every rank, before the opt-in gate and before any plotting
import. `tests/unit/test_setup_figure_title.py` asserts 131 raises with both
integers in the message and 130 does not, `ast`-scans all 18 `examples/` call
sites plus the docstring exemplar (13 constant titles inside the limit, 5
f-string/`%` titles listed and left to the runtime guard), and asserts the
absent-guard control at the pinned pre-change sha `8ef79a1`: 7 passed at `-n 1`
(`20260919T110950Z_OPS-53.log:51`, 1 s) and at `-n 2` on both ranks
(`20260919T111004Z_OPS-53.log:68,71`, 2 s). The census is unchanged —
`examples=54 ok=18 missing=36 broken=0` (`20260919T111013Z_OPS-53.log:89`, exit 2 on the 36 missing, as before) —
and the flagged `mesh:5` re-run is green with every `GEO-17` digit identical to
`20260918T140647Z_EX-57.log:952–963` and a byte-identical PNG
(`20260919T111033Z_OPS-53.log:952–963,976`, 12 s). What remains open is half
(a), below.

| | |
| --- | --- |
| **Symptom** | The PNGs for `mesh:3`, `th:10`, `mag:1/2/4/5/6`, `mri:1/2/3` and `ports:15–18` were rendered before `c9cc369` and show legends whose right edge overruns the subplot viewport — entries truncated mid-word. |
| **Verified at** | `e64cae5`, by the 04:30 slot reading PyVista 0.48.4's `renderer.py::map_loc_to_pos` (the `upper right` family anchors `x = 1 − size[1] − border`, off the box's *height*) and viewing the images. This review read the journal entries; it viewed no image. |
| **Cause** | A PyVista anchoring quirk, fixed for new renders by `c9cc369` (legend anchored `center left`, sized to the longest grouped label); committed PNGs do not change until re-rendered. |
| **Not caused by** | The examples: every imported identity is green in every figure's flagged window, and the figures are opt-in (`FEM_EM_SETUP_FIGURES`). The census is not wrong either — `broken` checks that the file parses, which it does. |
| **Scope** | Reader-facing artefacts only; no test, gate or number. **Not defects, recorded here so nobody files them again:** legends elide past two same-class labels (`first … last`) and one colour per region class cannot show a split between two same-class regions (`mesh:4`) — corpus design conventions, the operator's to change. |
| **Resolved by** | §9 items 6–7 (2026-09-19): both legs re-rendered — A (magnetostatics, five, `997a5fe`) and B (`mesh:3`, `th:10`, `mri:1/2/3`, `ports:15–18`, this commit) — each PNG read and each caption sentence checked against the new image. Retired. |

### ✅ RETIRED 2026-09-13 (`TH-15` step 3c, 19:30 implementer slot) — was OPEN 2026-09-13 (`TH-15` step 3, 15:00 implementer slot; filed by the 18:00 review on a `log-pathologist` ruling) — on the birdcage-as-PEC-hole route the **four-port terminal power sum exceeds the field accounting by an unattributed 7.6937e-05 W**: `Re Σ½V I*` = 7.700077682e-05 W against a phantom loss of 6.376395218e-08 W, 1 208×, on a mesh that has no other lossy volume

**Retired by step 3c: the excess equals a printed term. It is the terminal sheet form's Cauchy–Schwarz deficit, and it is not a hole-route readout systematic.**
The window was `20260914T005452Z_TH-15-step3c.log`: `-n 2`, complex build, `-s`, 21 passed, `[capture] rc=0`, 102 s.
- **Hole.** The excess over the volume loss reads **7.693701287e-05 W**. `C_total − Σ½|I|²Re Z_p` reads **7.693701287e-05 W**, identical to every printed digit (`:969–970`).
- **Algebra.** `V = V_s − I Z_p` per port, so `Σ½Re(V I*) − P_vol = (½Re(V_s Ī₁) − P_src) + (ΣP_sheet,field − Σ½|I|²Re Z_p)`. On the hole that is +4.685842647e-03 − 4.608905634e-03 W.
- **Solid.** The same split closes on the solid, where the excess is 6.716202469e-05 W, again equal to `C − terminal` (`:2790–2791`). That is `PORT-16`'s 10 MHz gap.
- **Ratio.** The driven-port pooled `C/terminal − 1` is **1.058874954e-02** on the hole and **1.059204217e-02** on the solid, so hole/solid = **0.9997**. The predicted 5× filing threshold is not approached (`:2796`).
- **Why 1 208×.** The hole's only lossy volume is the phantom, 6.376e-08 W. The sheets' ≈1.06 % deficit, 7.69e-05 W (2.8955e-02 of `P_src`), therefore dominates the terminal sum. On the solid, the 4.48e-04 W conductor loss masks the same deficit.

The row's power sentence is for the next review to re-register on this term (rule (h)). The body is kept below for the record.

**Where:** `tests/validation/test_th15_birdcage_pec_hole.py`, 10 MHz, `-n 2`,
`20260913T202311Z_TH-15.log:958–963` (Status 1 — the item's registered anchor
`Re P_in = ½∫_phantom σ|E|²` at 1e-3 went red on the terminal form) and
`20260913T202954Z_TH-15.log:954` (the same terminal and phantom digits,
re-printed "record only").

**What the logs establish:** `P_src` 2.657078677e-03 W and `Σ P_sheet,field`
2.657014913e-03 W differ by exactly the phantom loss 6.376395215e-08 W
(3.634e-10 relative) — `PORT-16`'s exact discrete identity holds on the hole.
The three undriven terminal terms (≈ −6.089e-04 VA each) equal `½|I|²·50` from
the printed currents to four figures, so the excess sits at the **driven
port**: its terminal term 1.897496397e-03 VA against a driven-sheet
dissipation of ≈ 8.365e-04 W by subtraction (the per-sheet split is computed
and discarded unprinted, `test_th15_birdcage_pec_hole.py:298–307`). The
module's comment (`:218–223`) attributes the gap to the driven terminal's
voltage carrying the source emf; **no printed quantity establishes that.**
The in-slot substitute anchor `P_src − ΣP_sheet = P_phantom` is not this
entry's fix: with the conductor solids removed (`…202954Z:80–81`, 10 fragment
volumes) it is the mesh's construction, not a finding, and it was ruled a
record by the 18:00 review.

**Why it matters:** the hole route's S-matrix is read from the same terminal
V and I. The identity gates (reciprocity, passivity, C4) are blind to a
common-mode readout error, and `TH-15` step 2 already measured a 2.07 % `Z`
asymmetry on the two-torus hole attributed to the point-sampled
`_path_voltage` (entry below, 2026-09-07). Whether this is the same mechanism
at the birdcage's surface-only terminals is the question.

**Discriminator (queued as `TH-15` step 3c, 18:00 review):** print the
per-sheet dissipation, the driven port's `C/terminal − 1` (`PORT-16`'s
`_exact_shares`, or `ports/shares.py` once `PORT-14` step 3 lands) on the
hole *and* the solid in one window, and a `½∫_{Ω∖phantom} σ|E|²` term. The
solid's 10 MHz `C/terminal − 1` is a record (2d: 1.059204e-02 on P1,
`20260912T123429Z_PORT-14-step2d-10mhz.log`); the hole's is the number this
entry waits on. Retires when the excess is attributed to a printed term and
the row's power sentence is re-registered on it, or is filed as the hole
route's terminal-readout systematic with its size.

### ✅ RETIRED 2026-09-12 (`PORT-19` step 5, 06:00 implementer slot) — `tests/validation/test_port19_factor_reuse.py::test_a_kept_fields_are_per_drive_and_match` was **red at `-n 2` in 2 of 2 windows** and **bit-identical at `-n 1`**. The kept `E` phasors differed by ≈4e-11 relative against the 1e-12 band, while `S`/`Z` agreed to ≤ 1.3e-14.

**Retired by the fixing commit, per the 03:00 ruling below.** The field half now asserts at `-n 1` only; at `-n 2` it prints and skips. `REUSE_REPRODUCTION_RTOL` stays 1e-12.
- `-n 1`: 18 passed, `Status: 0`, 111 s. `S`, `Z` and all four kept `E` read 0.000e+00 (`20260912T110109Z_PORT-19-step5-n1.log:1872–1873, :1950–1955`).
- `-n 2`: 17 passed, 1 skipped, `Status: 0`, 90 s. `S` 1.013e-14, `Z` 1.124e-14. Kept `E` P1 3.455e-11, P2 4.455e-11, P3 4.013e-11, P4 3.537e-11 were printed, not asserted (`20260912T110311Z_PORT-19-step5-n2.log:1891–1892, :1907, :2046–2051`).
- The `-n 2` MUMPS drift itself stays undiagnosed under the 2026-09-10 observation.

The body is kept below for the record.

| | |
|---|---|
| **Test** | `tests/validation/test_port19_factor_reuse.py::test_a_kept_fields_are_per_drive_and_match` (anchor (A), field half; `REUSE_REPRODUCTION_RTOL = 1.0e-12`). Verified at `42aad21`. |
| **Symptom** | `20260912T050110Z_PORT-19-step4-w1.log` (after the leg-offset module, same process; 1 failed / 22 passed, `Status: 1`, 217 s): `AssertionError: kept E phasors differ between paths beyond 1e-12: {'P1': 3.4554005026312407e-11, 'P2': 4.602345577977818e-11, 'P3': 4.013259366604542e-11, …}`; the printed line reads P1 3.455e-11, P2 4.602e-11, P3 4.013e-11, P4 3.537e-11, while `S` 5.015e-15 and `Z` 6.429e-15 (both in band). `20260912T050930Z_PORT-19-step4-isolation.log` (the module alone, step 2's exact shape; 1 failed / 17 passed, `Status: 1`, 119 s): P1 3.495e-11, P2 4.190e-11, P3 4.025e-11, **P4 0.000e+00**; `S` 1.130e-14, `Z` 1.297e-14. The other anchors and the stale-factor control (1.073e-02 at `S₂₂`, identical to step 2) are green in both. |
| **Discriminator (measured)** | `20260912T051457Z_PORT-19-step4-n1.log`, the same module at **`-n 1`**: 18 passed, `Status: 0`, 137 s; `S`, `Z` and all four kept `E` phasors **0.000e+00** (`:1872–1873`). Step 2's only prior window, `20260911T110544Z_PORT-19-step2.log:1891–1892` (`-n 2`), read 0.000e+00 everywhere. No `src/` commit between step 2's log and today touches the reuse path (`a20a8d6` adds a default-off test keyword to `build_four_port_sweep`, forwarded and ignored under `build_only`). |
| **Cause** | **Not diagnosed; consistent with the 2026-09-10 `-n 2` MUMPS non-reproducibility observation below**, here at ≈4e-11 on a field norm of the 139 140-dof birdcage rather than 1 ULP on a smoke scalar: at `-n 1` both paths are bit-identical, at `-n 2` the per-drive path's four separate factorisations are not bit-identical to the one held factor, and the deviation varies run to run (P4 3.5e-11 then 0.0). The 2026-09-10 exposure census predicted exactly this: "any future record at `-n 2` tighter than ≈1e-14 relative on a MUMPS solve would be exposed" — this band was added after the census. |
| **Consequence** | `main` carries one intermittent red at `-n 2` in this module. The reuse default is **not** implicated: `S`/`Z` agree to ≤ 1.3e-14 at `-n 2` and bit-for-bit at `-n 1`, and the stale-factor control separates at 1.073e-02. **Not relaxed** (CLAUDE.md hard rule); the default was not flipped. |
| **Resolves with** | A review ruling on anchor (A)'s field half, e.g. bit-identity asserted at `-n 1` (the `OPS-43` (d) precedent) with `-n 2` printed only, or a band registered by measurement of the `-n 2` repeat spread. Not an in-slot decision. |
| **Ruled 2026-09-12 03:00 review** | Anchor (A)'s field half is **re-anchored at `-n 1`** (bit-identity: S, Z and all four kept `E` read 0.000e+00, `20260912T051457Z_PORT-19-step4-n1.log:1872–1873`) and **printed only at `-n 2`**. `REUSE_REPRODUCTION_RTOL` stays 1e-12, and the S/Z halves stay asserted at every width (≤ 1.297e-14 measured at `-n 2`, `20260912T050930Z_PORT-19-step4-isolation.log:1891`). Rationale: the per-drive comparand — four separate factorisations — is itself not bit-reproducible at `-n 2` (the 2026-09-10 observation below, whose census predicted exactly this exposure), so the `-n 2` field comparison measures MUMPS run-to-run drift, not the reuse path. This is a width re-anchoring on the `OPS-43` (d) precedent, not a loosened bound; a `-n 2` band, if ever wanted, is registered from a measured repeat spread (≥ 5 windows). `PORT-19` step 5 (§9 item 2) lands it; this entry retires with that commit. |

### ✅ RETIRED 2026-09-11 (`ANS-4` step 2a″, 04:30 implementer slot — the guard landed and its control passed) — ~~🔴 OPEN 2026-09-09 (`ANS-4` step 2a, 04:30 implementer slot)~~ — the ×0.75 rung of the `ANS-4` convergence ladder **breaks C4 symmetry**: refining `conductor_resolution` moves the `Z` class spreads from **0.1012 / 0.0916 / 0.0654%** at ×1 to **0.5390 / 0.4591 / 1.6886%**, so the finer rung fails the imported, unmoved `ADJACENT_SPREAD_BAND` = 0.5% and is not a usable ladder point

| | |
| --- | --- |
| **Test** | `tests/validation/test_ans4_resolution_ladder.py::test_every_rung_passes_the_imported_port11_gates`, whenever the ladder contains a rung other than ×1. On `main`. |
| **Verified at** | `d6cd0fb` + this slot's one-knob edit. `20260909T093534Z_ANS-4-step2a.log` — `-n 4`, complex build, `tests/environment` first, `-v -s --tb=short`, `timeout -k 30 560`, `FEM_EM_ANS4_STEP2_DEGREE2=0 FEM_EM_ANS4_STEP2_RUNGS="1.0 0.75"`: **1 failed / 14 passed in 130.88 s**, `Status: 1`, elapsed **132 s** (`:3716–3717`, footer `:4008–4009`). Collect-only smoke first, 4 tests, `Status: 0`, 3 s (`20260909T093524Z_ANS-4-step2a.log`). |
| **Symptom** | `[ANS-4 step2] gates x0.75 degree 1` (`:3692`): `Z` class spreads **self 0.5390%**, adjacent 0.4591%, **opposite 1.6886%** against the imported 0.5% — the assert fires on `self` first (dict order), `opposite` is 3.4× the band. The same rung's **S** class spreads (`:3704–3706`) read self 0.4433%, adjacent 0.0549%, **opposite 2.0731%**. The ×1 rung on the identical route is green at **0.1012 / 0.0916 / 0.0654%** (`:3691`). |
| **What is green in the same window** | ×1 rung **116 085 cells, ratio 1.000000** against `GEO-19` step B's `STEP2_CELL_COUNT` inside `STEP2_CELL_COUNT_BAND` 2e-2 (`:3680`); the refinement negative control — ×0.75 meshes **161 695** cells against ×1's 116 085, +39.3%, reproducing `PORT-14` step 1b's measurement to the digit (`:3685–3686`), so `conductor_resolution` demonstrably refines; reciprocity `‖S − Sᵀ‖/‖S‖` = **9.194053e-16** (×1) and **2.271834e-15** (×0.75) against 1e-3; `σ_max` = **0.998974779044** / **0.998823907170** against 1 + 1e-9. Mesh 21.6 / 29.9 s, four drives 18.5 / 28.6 s at `-n 4`; ladder built in 104.2 s. |
| **Cause — not diagnosed** | Two candidates, unseparated in this slot: (a) **gmsh congruence** — the four legs are meshed independently and only the ×1 conductor size happens to land four near-congruent leg meshes, so a refinement that is not C4-covariant injects a geometric asymmetry (the `GEO-26` step 2 class of defect: a terminal triangulation that is not `C_N`-covariant, retired there by giving the rung its own measured record); (b) a genuine per-leg physical asymmetry the coarse ×1 mesh under-resolves. (a) is the more likely on the evidence that the *reciprocity* and *passivity* readings are untouched (1e-15, 0.9988) — a broken solve would move those, a broken symmetry of the discretisation would not. **Not measured**: the four legs' individual cell counts / sheet areas at each rung, which is the discriminator. |
| **Consequence** | `ANS-4` step 2a **does not deliver its ladder**: only the ×1 rung is a usable point, so no Richardson h → 0 estimate exists, no "move from ×1" trend is trustworthy beyond the single ×0.75 reading, and **`ANS-4`'s Larmor verdict stays INCONCLUSIVE**. The 128 MHz S entries were printed anyway for the record (see the §7 `ANS-4` row). No band was moved and no rung was kept by loosening one. The second pre-registered window (`"0.6 0.45"`) was **not run** — the item's own negative-result clause says stop. `main` now carries a red test in this module whenever the ladder is not ×1-only. |
| **Resolves with** | A review ruling on whether the ladder needs a **C4-covariant refinement** of the conductor mesh (the `GEO-26` step 2 remedy: measure the per-leg congruence at each rung and either enforce it in the generator or give each rung its own measured spread record, rule (f) class) before the convergence measurement can be made at all. **`ADJACENT_SPREAD_BAND` is not to be widened** — it is `PORT-11`'s gate and the ×1 rung meets it at 0.1%. Cheapest next measurement: per-leg cell counts and gap-sheet areas on the ×1 and ×0.75 meshes, no solve. |
| **Ruled 2026-09-09 10:30 review, ruling (1) — this entry and the `WF-6` step-4f entry below are the same open question, and one no-solve census answers both** | The discriminator this entry names as "not measured" is also what the `WF-6` step-4f entry needs, and 4f supplies the sharper evidence: refining the *global* `resolution` (a different axis from this entry's `conductor_resolution`) split the power residual **between the two ports** at ×0.0095 — P1 1.853642e-02 / P2 1.419812e-02, against P1/P2 agreement to five figures on both coarser rungs (`20260909T123716Z_WF-6.log:2089–2090, 3910–3911, 5777–5778`) — on the rung that also printed ungated `Z` class spreads 1.3475 / 0.4544 / 0.9165%. Two refinement axes, one fixture, one symptom class. **Disposition: `GEO-30`** (§7, §9 item 1) measures the per-leg cell counts, quadrant volumes and gap-sheet areas across **both** axes with no solve, exactly the measurement this row asks for. Nothing is ruled about candidates (a) vs (b) here, no band is widened or re-registered, and no rung is rehabilitated — `ANS-4`'s Larmor verdict stays INCONCLUSIVE. This entry stays open, and the ruling on the C4-covariant refinement is written from `GEO-30`'s table by the next review. |
| **Ruled 2026-09-09 18:00 review, ruling (1), from `GEO-30`'s table — the mesh is excluded in mass and in sheet *area*, and the one surviving mesh-side candidate is the sheet *triangulation*** | `GEO-30` ran green and answered with a negative (`20260909T171106Z_GEO-30.log`, `-n 2`, no solve, 133 s, all three anchors green): across all four rungs no **conductor** mass quantity leaves 0.35% — coil volume ≤ 0.173%, leg volume ≤ 0.328% — and the four gap-sheet **areas** agree to 1e-15 at every rung (6.05e-16 / 3.63e-16 / 8.47e-16 / 6.05e-16). On the very ×0.75 mesh that produced this row's 1.6886% `Z` spread, the coil-volume spread is **0.13%**: an order under `ADJACENT_SPREAD_BAND` and two under the symptom. **So candidate (b) — a genuine per-leg physical asymmetry the coarse mesh under-resolves — is not supported, and candidate (a) is narrowed rather than confirmed:** it is not gmsh congruence *in mass*. **The one quantity that tracks the failures exactly is the gap-sheets' facet count**: C4-equal on the two rungs that behaved (58/58/58/58 at ×1, 62/62/62/62 at `resolution` 0.012, spread 0) and **C2, not C4**, on precisely the two that broke (**80/74/80/74** at ×0.75, spread 7.79e-2; **70/76/70/76** at 0.0095, spread 8.22e-2) — opposite ports equal, adjacent ports differing, which is this row's own adjacent-vs-opposite class split. Identical area, different cut. **Disposition: `GEO-31`** (§7, §9 item 1, `mesh-probe`, no solve) reads the four sheets' triangulation directly — sorted per-facet areas, facet-centroid and boundary-vertex sets under the C4 rotation map — and returns the discriminator *same support, different interior cut* (a reconstruction/quadrature question) versus *different support or boundary* (a generator defect). **Nothing is ruled about the mechanism here, no rung is rehabilitated, `ADJACENT_SPREAD_BAND` is not widened, and `ANS-4`'s Larmor verdict stays INCONCLUSIVE.** This entry stays open. |
| **Update 2026-09-10 19:30 slot (`ANS-4` step 2a′) — the cause is the port-sheet cut, and a flag-on build clears the red** | With `GEO-32`'s opt-in `c4_congruent_sheets` threaded through `FEM_EM_ANS4_STEP2_C4_CONGRUENT=1` and nothing else changed, the same window (`-n 4`, `RUNGS="1.0 0.75"`, `DEGREE2=0`, `timeout -k 30 560`, `-s`) is **15 passed, Status 0, 147 s** (`20260911T003204Z_ANS-4-step2a-prime.log`). The ×0.75 `Z` spreads read **0.0481 / 0.0649 / 0.0555 %** (`:3913`) and the S31 class spread 0.0267 %, against the unmoved 0.5 %. The flag-off control, run in the same slot, reproduces this entry's symptom to the printed digit (`20260911T003442Z_ANS-4-step2a-prime.log:3694, :3718`; `S` entries and σ_max identical to `20260909T093534Z`). **So candidate (a) is confirmed in its narrowed form, the sheet triangulation, and (b) is excluded.** The test stays red on the **default** (flag-off) path whenever the ladder is not ×1-only, because the default was deliberately not flipped. This entry stays open for the review to retire or re-head: flip the default, or require the flag for any ladder beyond ×1. No band moved. |
| **Ruled 2026-09-11 03:00 review — require the flag beyond ×1; do not flip the mesher default** | Flipping `c4_congruent_sheets`' default in `MeshGenerator.birdcage_port_domain` would move the `GEO-19` 116 085 record, and every gate that imports it, to 116 118. That is a record sweep bought for a test-side problem. **Instead**, `test_ans4_resolution_ladder.py` refuses to build a `conductor_resolution` ladder with any factor ≠ 1.0 while `FEM_EM_ANS4_STEP2_C4_CONGRUENT` is **unset**, raising before any mesh on every rank and naming this entry. An **explicit** `=0` stays allowed, because it is the documented flag-off negative control (`20260911T003442Z_ANS-4-step2a-prime.log`). §9 item 1 (03:00 queue) implements the guard with its own negative control. **This entry retires in that item's commit**, since the default path can then no longer reach the red silently. No band moves. |
| **Retired 2026-09-11 04:30 slot (`ANS-4` step 2a″) — guard landed, control green, flag-on ladder green to ×0.45** | `_require_explicit_c4_flag` in `tests/validation/test_ans4_resolution_ladder.py` runs at the top of the `ladder` fixture, before any mesh. It raises `ValueError` naming this entry and both remedies when `FEM_EM_ANS4_STEP2_C4_CONGRUENT` is unset and any `conductor_resolution` factor ≠ 1.0. `RUNGSPEC`/`RESOLUTION` paths untouched; explicit `=0` still allowed. **Guard control (flag unset, `RUNGS="1.0 0.75"`, `-n 2`, `timeout -k 30 120`):** `20260911T093156Z_ANS-4-step2a-dprime-w0.log` — all four ladder tests **ERROR at setup** with the message (`:81–104`), no rung/cells line (nothing meshed), `11 passed, 4 errors` (the 11 are `tests/environment`, `:172`), Status 1, **31 s**. With the flag on, ×0.6 and ×0.45 also pass the unmoved 0.5 % spread band (`20260911T093524Z_ANS-4-step2a-dprime-w2.log:5966–5968`; §7 `ANS-4` step 2a″). The default path can no longer reach the red silently. No band moved and the mesher default stays off. |

### ✅ RETIRED 2026-09-08 (`OPS-41`, 07:30 implementer slot) — ~~🟡 OPEN 2026-09-07 (`TH-15` step 2d, 12:00 implementer slot) — `tests/validation/test_port_package_sparameters.py` fails **three digit-reproduction records at `-n 4`** and passes all 17 at `-n 2`: the records are `-n 2` records and the module is rank-width-sensitive at 1e-4, four decades above its 1e-6 reproduction bands~~ — **the records are now declared `-n 2` records in the module (`RECORD_RANK_WIDTH = 2`), the `-n 4` miss is a pinned strict `xfail` behind `FEM_EM_RECORD_WIDTH_OVERRIDE=1`, and the 1e-4 is attributed by measurement to the point-sampled gap voltage `V`, not to the current `I`**

| | |
| --- | --- |
| **Where this fires** | `tests/validation/test_port_package_sparameters.py::test_package_sweep_reproduces_the_gated_mutual`, `::test_sanity_report_reproduces_the_gated_metrics_on_the_field_route`, `::test_reciprocity_warning_fires_on_an_asymmetrised_field_smatrix` — **at `-n 4` only**. Not a `TH-15` failure: the step-2d diff (`ports/gap_voltage.py`, `ports/excitation.py`) is additive and leaves the `current_route="conduction"` arithmetic byte-identical, which is what the `-n 2` window proves. |
| **Verified at** | `8ef690a` + the step-2d worktree. `-n 4`: `20260907T170528Z_TH-15.log` — **3 failed / 3 passed in 145.46 s**, elapsed 147 s, `Status: 1` (`:743–749, :755, :882–883`). `-n 2`, same tree, `tests/environment` first: `20260907T170838Z_TH-15.log` — **17 passed in 187.37 s**, elapsed 188 s, `Status: 0` (`:146, :214–215`). |
| **Literal symptom** (`…170528Z:743–749`) | `raw mutual 0.8942257288323762 does not reproduce the re-based 0.11 record 0.8945163786446685 within 1.0e-06 relative (missed by 3.249e-04)`; `passivity_max_sigma 0.864692569 does not reproduce the step-4 record 0.864809 within 1.0e-06`; `the perturbation did not land in the reciprocity metric: 9.988682e-02 vs 1.000000e-01`. |
| **Cause — the records' width, not diagnosed further** | Every one of the three asserts reproduction of a *recorded digit string* at 1e-6 relative. The records were established at `-n 2` (`20260904T110501Z_OPS-37.log:12`, `mpiexec -n 2`, 17 passed / 171.42 s). Changing the rank count repartitions the 184 176-cell two-torus, which moves the solved S-matrix at the 1e-4 level — the size of all three misses. No band or physics gate fails: the mutual, reciprocity and passivity *bands* are green in the same `-n 4` run; only the digit records miss. Whether the module can be made width-invariant, or its records should be declared `-n 2`-only in the module docstring, is a review's call. |
| **Consequence** | Run this module at `-n 2`, the width its records were set at and the width CI uses. `TH-15` step 2d's rule-(c) re-run was executed there and is green. |
| **Resolves with** | A review scoping either (a) a docstring/`pytest` marker declaring the records `-n 2`-only, or (b) a width-invariance measurement on this fixture (the `OPS-18` "same command, different digits" entry below is the adjacent precedent — that one is same-width run-to-run at 1e-10; this is 1e-4 across widths, six decades larger, so they are probably not the same defect). |
| **Ruled 2026-09-07 18:00 review — option (a), plus one attribution print** | The solve is a MUMPS direct factorisation (`core/time_harmonic.py:546–548`), which moves a well-conditioned answer at ≈1e-12 with partition, not 1e-4; the likelier carrier is the point-sampled path voltage — `_path_voltage` evaluates the N1curl `E` at points along the gap path, and a point on a partition boundary can be owned by a different cell at a different width, where only the tangential component is continuous. **`OPS-41`** declares `RECORD_RANK_WIDTH = 2`, skips the three record tests at other widths *after printing the reading beside the record*, pins the `-n 4` miss as a strict `xfail` behind an override variable, and prints `V` and `I` per drive at both widths to attribute the 1e-4 (predicted: `V` moves, `I` does not). Until it lands, every rule-(c) gate on this module is written `-n 2`. Entry retires with `OPS-41`'s landing commit. |
| **Retired by (`OPS-41`, 2026-09-08 07:30 slot) — the prediction is confirmed by measurement: `V` carries the 1e-4, `I` does not** | Test module only, no `src/` change, no band moved (`REPRODUCTION_BAND_RELATIVE`, `PASSIVITY_REPRODUCTION_BAND`, `SYMMETRY_RATIO_BAND`, `MUTUAL_TOLERANCE`, `HEURISTIC_SEPARATION_FLOOR` all untouched). Three windows, all `Status: 0`. **`-n 2`** (`20260908T123716Z_OPS-41.log`, `tests/environment` first, 17 passed in 179.30 s, elapsed 181 s, `:813`, `:881–882`): the records reproduce as before — raw miss 4.269e-10, corrected 4.190e-10 (`:716`), `passivity_max_sigma` miss 2.616e-10, symmetry-ratio miss 4.889e-11 (`:735`), perturbation recovery 4.070e-04 relative (`:743`). **`-n 4`** (`20260908T124026Z_OPS-41.log`, **3 passed, 3 skipped** in 134.84 s, elapsed 136 s, `:759`, `:771–772`): the three skips fire *after* the readings print, and every width-independent assertion still runs and is green at that width — the mutual band (raw −10.58%, corrected −6.05% against the unmoved 10%), reciprocity, passivity, and the heuristic separation 3.031654e-01 against the 2.0e-3 floor (`:714`, `:737`). **`-n 4` + `FEM_EM_RECORD_WIDTH_OVERRIDE=1`** (`20260908T124317Z_OPS-41.log`, **3 passed, 3 xfailed** in 140.38 s, elapsed 142 s, `:751`, `:763–764`): the three records miss at the 1e-4 class as a *strict* `xfail` — raw 3.249e-04 / corrected 3.189e-04 (`:702`), `passivity_max_sigma` 1.169e-04 and symmetry ratio 8.634e-05 (`:738`), perturbation recovery 1.132e-03 (`:748`) — reproducing `…170528Z:743–749` to the digit. **Attribution** (repr-precision `V_i^(k)` / `I_k` prints, `…123716Z:709–715` vs `…124026Z:719–725`): the driven currents move `|ΔI|/|I|` = **5.8e-11** (`I_P1`) and **1.5e-10** (`I_P2`), the undriven receiver currents **7.6e-10** / **3.6e-09**; the gap voltages move `|ΔV|/|V|` = **3.249e-04** (`V_P2^(P1)`, i.e. exactly the raw-mutual miss), **1.176e-04** (`V_P1^(P2)`), **5.617e-04** (`V_P1^(P1)`) and **2.287e-02** (`V_P2^(P2)`, the ungated `Z22` diagonal). So the point-sampled `_path_voltage` carries the width sensitivity — six to eight decades above the facet-integrated current — corroborating `TH-15` step 2f's independent finding (`d348db0`) that the point-sampled voltage, not the field, is that fixture's 2% `Z` asymmetry. **A fix (edge-integrated or face-averaged `_path_voltage`) is a separate `src/` chunk this measurement licenses**; this entry retires because the domain of the records is now declared and the off-width miss is pinned rather than forgotten. |

### ✅ RETIRED 2026-09-13 (`ANS-4` Larmor verdict, weekly review — banked by the operator's interactive session after the scheduled session died on the account limit) — ~~🟡 OPEN 2026-09-06 (`ANS-4` adjudication, weekly review run interactively) — the 4-leg birdcage's 4×4 S-matrix **agrees with HFSS at 10 MHz and is inconclusive at 64 / 128 MHz**: a miss that grows with frequency, several times larger than HFSS's own Zero-vs-First-Order shift, on a fixed 116 085-cell mesh that has never had a convergence ladder on this fixture~~ — **the cause is identified and the verdict is banked: at degree 2 (HFSS First Order by `ANS-5`) on a mesh of comparable unknown count, step 2d's order-matched rung reads AGREE at 128 MHz, and the mechanism — our degree-1 gate fixture unconverged while AED's run had converged — gives AGREE by mechanism at 64 MHz pending `ANS-4` step 3 (`xl`, pre-registered in §10). The 09-06 'element order is excluded' inference was the error: AED was order-insensitive because it had converged. Numbers private; §2.2 / §6 carry the verdict word.**

| | |
|---|---|
| **Test** | None red: every gate on the fixture is a self-consistency identity (reciprocity / passivity / C4) and all pass at all three frequencies (`PORT-9`, `PORT-11`, `ANS-4` runnable half). This entry records that the **absolute** comparison is open at the Larmor frequencies, so no absolute-accuracy claim may cite the 64 / 128 MHz S-matrix. |
| **Log** | Ours: `docs/testing/logs/20260905T041254Z_ANS-4-private.log` (Status 0, 128 s at `-n 2`). AED side: private (`aed_results/`, `COMPARISON_private.md`, `docs/private/ans4-adjudication-2026-09-06.md` — licence terms; no number here). |
| **Symptom** | Qualitative only: at 10 MHz couplings, self term and the S-derived terminated-drive input impedance agree to the few-percent level the two AED orders bracket; at 64 and 128 MHz the miss grows with frequency on every C4 class. AED's two orders agree with each other to well under the gap everywhere. |
| **Cause** | Not diagnosed by measurement. Element order is **excluded** (the AED order shift is several times smaller than the gap). First suspect: our mesh — ≈ 12 cells/λ and ≈ 5 cells/δ in the 800 S/m legs at 128 MHz, no h-ladder on this fixture, and the sign of the miss matches more conductor loss being resolved inside the legs. Second: a Larmor-frequency feed-model systematic of the lumped sheet. Only moving *our* discretisation separates them. |
| **What is unmoved** | No band, record or gate moved. The 10 MHz 4×4 gains "externally checked, AGREE" (§2.2, §6). `ANS-1`'s AGREE verdict is unaffected. |
| **Disposition** | **`ANS-4` step 2 (`xl`)** — the first XL-slot run (§5.1): four degree-1 conductor-resolution rungs plus one degree-2 solve at 128 MHz in one window, Richardson h→0 printed; decision rule pre-registered in the §7 row (≥ 50 % of the private gap closed → discretisation; < 25 % → feed model → a birdcage-sheet Larmor systematics chunk; between → one more rung). Retire this entry when the 09-09 or a later weekly rules. |
| **Update 2026-09-10 — the cause is identified and the entry is ready to retire** | `ANS-4` step 2d ran the order-matched comparison (`20260910T180911Z_ANS-4-step2d-capture.log`, 14 passed / 1 failed, 7225 s, `-n 16`, peak 290.2 GiB). At **degree 2** — which is HFSS First Order by `ANS-5` — on a mesh of comparable unknown count to AED's own converged run, the disagreement falls to a small fraction of what the degree-1 fixture showed *(the percentages were redacted 2026-09-10 by the 18:00 daily review: an ours-vs-AED gap is a benchmark number, CLAUDE.md hard rules; they live only in the private note below)*, and the degree-2 sequence converges (successive change 1.19 % → 0.62 %, ratio 1.90). **The 2026-09-06 reasoning was wrong**: it excluded element order because AED's own order shift was small, but AED's answer was order-insensitive because *its* run had converged and ours had not — the asymmetry was the signal, not an exclusion. Numeric ruling in gitignored `docs/private/ans4-step2d-2026-09-10.md`. **Retire this entry when the 2026-09-13 weekly banks the verdict** — the reading the evidence supports is AGREE at 64/128 MHz, and that is a review's to bank, not a run's. |

### ✅ RETIRED 2026-09-06 (`OPS-39`, 15:00 implementer slot) — ~~🟡 OPEN 2026-09-06 (`TH-15` step 2a, 00:00 implementer slot; filed by the 10:30 review) — the step-0 mesh probe's `_census` **double-counts ghost entities at any rank count above 1**: under `GhostMode.shared_facet` the solid two-torus tag census overshot the owned cell count by **2972** at `-n 2`, reading tag 2 at **9448 against 9348** (1.0697% > the 1% band)~~ — **the probe's `_census` now masks on `index_map(dim).size_local`, gated on the exact identity `Σ_tag census[tag] == size_global` at `-n 1` and `-n 2` with the naive overshoot asserted equal to `Σ_ranks num_ghosts` (measured 256 at `-n 2`, 0 at `-n 1`); the probe's solid variant at `-n 2` reproduces step 0's `-n 1` digits exactly**

| | |
|---|---|
| **Test** | None red on `main`: the gate module `tests/mesh/test_two_torus_conductor_hole.py` masks its own census on `size_local` (`_owned_census`, `:119–134`) and asserts the per-tag census sums to the global owned count (`:219`), so it is green. The defect lives in `tests/mesh/probe_two_torus_conductor_hole.py::_census` (`:431–440`), a measurement-only probe that asserts nothing — and in anything else that imports it at width > 1. |
| **Log** | `docs/testing/logs/20260906T050521Z_TH-15.log` — the raw solid census at `:1464` (9556 / 9448 / 113 483 / 13 661 / 13 648 / 13 658 / 13 694, sum 187 148 against `n_cells` 184 176), the failing assert at `:1483–1485`, Status 1 / Elapsed 130 s at `:1632–1634`. The green re-run with the masked census, same tree: `20260906T050829Z_TH-15.log:1464` reads 9471 / 9348 / 110 696, `12 passed, 4 skipped` at `:1543`. |
| **Symptom** | `_census` calls `np.unique(values)` on the rank-local `tags.values`, which on a `shared_facet`-ghosted mesh includes every ghost cell the rank holds, then `allreduce`s the counts — so a cell shared by *k* ranks is counted *k* times. The total cell count, the sheet areas and the cavity area were already exact on the failed run; only the per-tag counts drifted, by exactly the ghost count. |
| **Cause** | Rank-safety: `cell_tags.values` is rank-local *and* ghost-inclusive (CLAUDE.md hard rule); the probe was written for and read at `-n 1` (step 0, 2026-09-05), where owned = local and the two censuses agree. |
| **What is unmoved** | No record moved: step 0's `-n 1` digits (161 461 / 184 176; 9471 / 9348; areas to 3.7e-7) are the ones the gate module now holds at `-n 2`. No band moved — the measurement was wrong, not the record. |
| **Disposition** | **CLOSED by `OPS-39`, 2026-09-06 15:00 slot.** `_census` gained two optional keyword arguments (`indices=`, `n_owned=`) that mask on `index_map(dim).size_local` before counting, plus a probe-local `_owned_census(msh, tags, dim, comm)` wrapper now used at all three probe call sites; the positional `(values, comm)` signature and the `{tag: count}` return type are unchanged, so the gate module's import stays green (`20260906T200326Z_OPS-39.log`, `5 passed`, 57 s at `-n 2`). New gate `tests/mesh/test_probe_census_rank_safety.py`: on a 3072-cell `create_box` with `shared_facet` ghosting the masked census sums to `size_global = 3072` exactly and the naive sum overshoots by exactly `Σ_ranks num_ghosts` — **0 at `-n 1`** (`20260906T200217Z_OPS-39.log:43–45`) and **256 at `-n 2`**, naive `{1: 774, 2: 770, 3: 1784}` sum 3328 (`20260906T200227Z_OPS-39.log:51–53`). The probe's solid variant re-run at `-n 2` prints step 0's `-n 1` digits to the integer: `{1: 9471, 2: 9348, 3: 110696, 101: 13661, 102: 13648, 111: 13658, 112: 13694}`, 184 176 cells (`20260906T200239Z_OPS-39.log:950–952`) — so there is no rank-count mesh difference, only the ghost double-count. |

### ✅ RETIRED 2026-09-11 (`PORT-14` step 1e, 07:30 implementer slot) — systematic registered as a record — ~~🟡 OPEN 2026-09-05 (`PORT-14` step 1, 21:00 implementer slot) — the lumped sheet is **single-mode to 3.4e-03, not 1e-3**: the field-solved terminated 3×3 misses the circuit reduction of the 50 Ω 4×4 by **1.596e-03 (C = 100 pF)**, **3.371e-03 (L = 1 µH)** and **7.250e-04 (R = 200 Ω)** against the pre-stated 1e-3 band~~ — **the three residuals are now the F-small (1\*) record `REDUCTION_FLOOR_F_SMALL`, asserted at rtol 1e-3 at `-n 2` (|ratio − 1| ≤ 2.4e-07); `REDUCTION_BAND` = 1e-3 is unmoved and printed as a named systematic of the lumped sheet**

| | |
|---|---|
| **Test** | `tests/validation/test_port_lumped_rlc_termination.py::test_the_terminated_solve_matches_the_circuit_reduction` — the other two tests in the module pass. **Green since 2026-09-11 (step 1e, record-asserted).** Was a **deliberate red on `main`**, the `POST-6` step-1 precedent one slot earlier: §9 item 2's negative-result protocol for a residual in (1e-3, 1e-2] is "known-issues entry with the three residuals, row 🟡, band not widened", and no band is widened in-slot. |
| **Log** | `docs/testing/logs/20260905T020428Z_PORT-14.log` — residuals and the meas/pred 3×3s at `:1857–1878`, Γ = 0 control at `:1881–1884`, 50 Ω baseline at `:1854`, terminations and per-case wall times at `:1850–1853`, failure at `:1891–1892`, `1 failed, 2 passed … 102.72s` at `:1896`, Status 1 / Elapsed 105 s at `:1911–1912`. Heavy tier by ceiling (`timeout -k 30 560`), `-n 2`, complex build, 13 solves on one 116 085-cell mesh. Verified at the closing commit of the same slot. |
| **Symptom** | `‖S'_meas − S'_pred‖_F/‖S'_pred‖_F` with port P1 terminated and P2/P3/P4 driven: **C = 100 pF (Z_p = −j1.591549e+02 Ω, \|Γ\| = 1.000000) → 1.595580e-03**; **L = 1 µH (+j6.283185e+01 Ω, \|Γ\| = 1.000000) → 3.370512e-03**; **R = 200 Ω (Γ = +0.600000, \|Γ\| = 0.6) → 7.249519e-04**, the last one *inside* the band. The identity `S' = S_aa + S_ab Γ (I − S_bb Γ)⁻¹ S_ba` is exact for a single-mode port, so the residual is a measurement of the sheet's port model, not a discretisation error in the usual sense. |
| **Cause** | **Not diagnosed by measurement, but the ordering is the hypothesis: the residual tracks \|Γ\|, not the element type.** The two lossless terminations (\|Γ\| = 1, all power returned to the sheet) read 1.6e-03 and 3.4e-03; the lossy one (\|Γ\| = 0.6) reads 7.2e-04, ≈2–5× smaller. A single-mode port's higher-order sheet content is re-excited by whatever returns, so a total reflection stresses the assumption hardest — consistent with "the sheet is not exactly one mode" and inconsistent with a complex-`Z_p` arithmetic bug (which would not care about \|Γ\|, and which would not leave R = 200 Ω inside the band while C = −j159 Ω is outside). Not measured here: whether refining the sheet, or widening the interior-width fraction, moves it. |
| **What is unmoved** | The 50 Ω baseline reproduces `PORT-9` on this run's mesh — `‖S − Sᵀ‖/‖S‖ = 1.464324816e-14` against the imported 1e-3 `RECIPROCITY_BAND`, `σ_max = 0.999992805` against `1 + PASSIVITY_SIGMA_TOLERANCE` (`:1854`) — so this is not a broken 4×4. The ceiling-first Γ = 0 negative control **passes with room**: the 4×4's own coupling term is `Δ = 3.218888e-01 / 3.254627e-01 / 2.112830e-01`, all far above the 5e-3 floor, and `S_aa` misses the measured 3×3 by `3.219520e-01 / 3.267853e-01 / 2.120063e-01` — 200–330× the band, so the gate above is resolving the termination and not passing on a weakly-coupled network (`:1881–1884`). |
| **Disposition** | `REDUCTION_BAND` stays **1e-3** and is not widened; `PORT-14` step 1 is 🟡 — complex `Z_p` on a lumped sheet exists and solves, the reduction identity is measured, the band is missed by ≤ 3.4×. `PORT-15`'s gate (i) should compare against the measured 3×3 with these three numbers in hand rather than re-deriving them. Expected on `main`; **not yours**. **Review 2026-09-05 03:00:** the |Γ| ordering is accepted as the working diagnosis (not an arithmetic bug — one that left R = 200 Ω inside the band while C = −j159 Ω sits outside would be a strange one); **`PORT-14` step 1b** (§9 item 4) re-measures the two lossless terminations on two conductor-resolution rungs with the residuals printed and the band unmoved. This entry leaves with a gated residual ≤ 1e-3 or with a review that re-registers the band on measured evidence — never in-slot. **Step 1b (2026-09-05, 09:00 slot) refuted the resolution hypothesis** — the residual is non-monotone in `conductor_resolution` (×2.63 at ×0.75, ×0.93 at ×0.6) — and named `sheet_width_m` as the suspect. **Step 1c (2026-09-05, 21:30 slot) confirms it on the fixed 116 085-cell gate mesh** (`20260905T213322Z_PORT-14-step1c.log`, `9 passed, 6 deselected in 230.46s`, Status 0, 258 s): perturbing only the width the sheet law is *told* moves the residual by ×5.80 / ×5.83 at +5%, ×3.75 / ×3.72 at −5% and ×5.63 / ×5.68 at alternating ±5.3% (C / L), with every perturbed 4×4 still reciprocal to ≤ 1.9e-14 and passive to σ_max ≤ 0.999997273 and the mesh bitwise unchanged. Both directions *raise* the residual, so the nominal `A/h` is already near the optimum; a three-point fit (arithmetic on the log, not a measurement) puts the common zero-crossing at ε\* ≈ −0.011 for both elements with a fitted minimum ≈ 0. So the cause hypothesis above is superseded: the sheet's single-mode residual is dominated by the effective-width constant, not by \|Γ\| or by mesh density. The band is still not widened and this test is still deliberately red on `main`. **Review 2026-09-05 18:00:** reading (1) licenses **step 1d** (§7 `PORT-14` "Step 1d", §9 item 1): the sheet geometry `build_four_port_sweep` already measures (`area`, `h` = bounding-box extent along the drive, `w_bbox`, `facets`) is read per port and compared against the fitted ε\* to name a *geometric* origin for the 1.1%, and the fit's prediction is tested on a fourth point — the told width scaled by (1 + ε\*) — with the residual printed. No constant is tuned to fit: an effective-width correction enters the law (a step 1e, `src/`) only if a measured geometric ratio predicts ε\* on its own. This entry still leaves only with a gated residual ≤ 1e-3 or a review that re-registers the band on measured evidence. **Step 1d (2026-09-05, 19:30 slot) selects pre-registered reading (2)** (`20260906T003627Z_PORT-14-step1d.log`, `4 passed, 15 deselected in 110.42s`, Status 0, 137 s): ε\* recomputed in code from step 1c's residuals is **−0.010735 (C) / −0.010970 (L)**, and configuration D — all four told widths × 0.989147732, same mesh, same `reuse` route — reads **5.756561e-05 (C)** and **1.196780e-04 (L)**, i.e. ×0.036 of step 1's records and 0.058× / 0.120× the band (`:2005–2007`, `:2094–2105`), so the `\|r₀ + kε\|` model holds on a fourth point (reading (3) excluded) — but the residual is a **fitted** number and closes nothing. **No measured geometric ratio predicts ε\***: on all four identical sheets the bounding-box fill / ragged-edge ratio is **−0.204336** (inverse +0.256812) and the told `leg_gap_length` is reproduced **exactly** by `h_bbox` (8.000000000e-03 m, ε_geom = +0.000000), against the pre-registered ±0.003256 window — nearest miss 0.010852 = 3.3× the window (`:2001–2004`, `:2008–2040`). So the 1.09% is a **field effect of the sheet** (edge fringing — the single-mode residual proper), not a mis-measured terminal separation: **no constant enters the law and there is no step 1e**. D's own 4×4 is a valid passive reciprocal network (cells 116 085 bitwise, reciprocity 9.998294990e-15, σ_max 0.999993774, `:2099`) with the Γ = 0 control asserted at Δ = 0.321025 / 0.324791 and a 321× / 325× miss (`:2110–2111`). The band is still not widened and this test is still deliberately red on `main`; re-registering it on this measured evidence is the weekly review's call. **Weekly review 2026-09-06 ruled:** the band is not re-registered, and the measured floor becomes a (1\*) record on F-small. **RETIRED by `PORT-14` step 1e (2026-09-11, 07:30 slot):** `test_the_terminated_solve_matches_the_circuit_reduction` now asserts each residual against `REDUCTION_FLOOR_F_SMALL` (cited to `20260905T020428Z_PORT-14.log:1858, 1865, 1872`) at `REDUCTION_FLOOR_RTOL` = 1e-3. The assert runs only at `RECORD_RANK_WIDTH` = 2, on the `OPS-41` precedent. It prints each residual against the unmoved 1e-3 band as "systematic, not gated". `20260911T123334Z_PORT-14-step1e.log`: C 1.595580e-03 (ratio 0.999999762), L 3.370512e-03 (1.000000106), R 7.249519e-04 (0.999999966) (`:1888, :1895, :1902`). Cross-element negative control `|r_C/rec_L − 1|` = 0.526606 and `|r_R/rec_C − 1|` = 0.545650 (`:1913, :1917`). 50 Ω baseline and Γ = 0 control unmoved (`:1884`, `:1922–1924`). `15 passed, 16 skipped … 110.80s` (`:2050`), Status 0 / 113 s (`:2118–2119`), verified at `5e8f42a` + the step 1e diff. §9's residual-red tally goes 4 → 3. |

### ✅ RETIRED 2026-09-07 (`PORT-16` step 2, 06:00 implementer slot) — ~~🟡 OPEN 2026-09-04 (`POST-6` step 1, 19:30 implementer slot) — the **drive-level** power identity `½aᴴ(I − SᴴS)a = ½∫σ|E_w|²` misses by **11.648%** against the imported 1% `POWER_BALANCE_BAND`, and the **single drives on the same fixture miss the same identity by ≈13%**: the band was pre-registered against the wrong denominator~~ — **the drive-level test now asserts the exact discrete identity on the superposed field (rel dev 5.877e-15 / 1.435e-14 against 1e-6) and its attribution (9.881e-14 / 1.581e-13 against 1e-1); the 11.648% is printed as a record beside the untouched band, which keeps its value and its single-drive asserts**

| | |
|---|---|
| **Test** | `tests/validation/test_port_drive_superposition.py::test_the_drive_level_power_identity_closes` — the other nine tests in the module pass (six from step 1, three added by step 1b). A **deliberate red on `main`**: the §9 item's negative-result protocol for anchor (iii) is "known-issues entry, row 🟡", and no band is widened in-slot. |
| **Log** | `docs/testing/logs/20260905T004202Z_POST-6.log` — readings at `:1915–1920`, failure at `:3772–3773`, `1 failed, 23 passed … 144.48s` at `:3839`, Status 1 / Elapsed 146 s at `:3917`. Heavy tier by ceiling (`timeout -k 30 590`), `-n 2`, complex build. Verified at the closing commit of the same slot. |
| **Symptom** | ccw quadrature: `P_acc = ½aᴴ(I − SᴴS)a` = **3.014424803e-03 W** (available 1.000000000e-02 W) against a superposed-field volume loss `½∫σ|E_w|²` = **2.663302665e-03 W** (phantom 3.796523707e-07, conductor 2.662923013e-03) — residual **1.164806e-01** vs the 1e-2 band. cw is the same to seven digits (3.014424670e-03 / 2.663302637e-03, `:1919–1920`). Cross-term share 40.4652%, blind sum `Σ|w_k|²P_k` = 1.794631250e-03 W (`:1918`). |
| **Cause** | **Not diagnosed by measurement, but arithmetic on the same log makes one hypothesis dominant: the band's denominator, not the superposition.** `WF-6` step 1's gate (i) reads its 9.795751e-03 residual against the **supplied** power (6.856240413e-03 W, `20260831T033704Z_WF-6-step2.log:4684`), of which the sheets take 92.48%. The same fixture's single-drive form of *this* identity is `P_acc,1 = supplied − sheets = 6.856240413e-03 − 6.340800348e-03 = 5.154401e-04 W` against a volume loss of `4.482780e-04 W` — a **13.0%** miss, i.e. the four-port drive reproduces the single-drive miss rather than adding to it. The superposition itself is exact to round-off on the same solves: the package path and the fixture's DG0 path agree to **1.2e-15 / 1.9e-15 / 1.2e-16** on the three identity readings (`:1911–1913`), and all four single-drive accountings reproduce step 1's 9.795751e-03 (`:1914`). The two candidate mechanisms left are (a) the 1% band was pre-registered against `supplied` and re-used against a quantity ≈13× smaller, or (b) a real un-accounted loss channel of ≈3.5e-04 W common to every drive. The single-drive arithmetic above is *computed from a log, not asserted in-run*; the step that closes this must measure it. |
| **What is unmoved** | Anchors (i), (ii) and (iv) all pass: C4 **0.9818%** and mirror **0.8087%** through `ports.superpose_drives` reproduce `WF-6` step 2's records at rel dev 9.619e-06 / 3.585e-05, the cw negative control reproduces 95.1975% at 5.248e-07 (`:1911–1913`); `w = e_k` returns drive k's dof array bit for bit; `S(w₁)+S(w₂) = S(w₁+w₂)` and the terminal weighted sums hold at 1e-12; the four single-drive residuals read 9.795751e-03 / 9.796209e-03 / 9.794985e-03 / 9.795283e-03 (`:1914`). `tests/validation/test_port_birdcage_four_port.py` is green in the same window (the `run_lumped_sheet_port_case` gate re-run for the additive `return_fields` keyword). |
| **Disposition** | `POWER_BALANCE_BAND` stays **1e-2** and is not re-pointed at a different denominator in-slot — that is a review's call, and it needs the single-drive `P_acc` measured in-run beside the superposed one, not read out of a log by hand. `POST-6` step 1 is 🟡: the entry point exists and its field/terminal/identity anchors are gated; the drive-level power statement is not. Expected on `main`; **not yours**. **Review 2026-09-05 03:00:** the review's own arithmetic on the same log — `WF-6` step 1's absolute accounting gap is 9.795751e-03 × 6.856240413e-03 = **6.716e-05 W**, which is 13.0% of `supplied − sheets` = 5.154401e-04 W, and the superposed drive's 3.511e-04 W on 3.014e-03 W is the same 11.6% fraction — makes (a) the working diagnosis; what separates (a) from (b) is whether the S-derived `P_acc,k` equals `supplied_k − Σ sheets` (the power-wave identity), which no run has asserted. **`POST-6` step 1b** (§9 item 2) asserts it at rtol 1e-6 and prints the four single-drive residuals; the review after it re-points or keeps the band. This entry leaves with that decision's commit, not before. **Step 1b executed 2026-09-05, 06:00 implementer slot — the identity holds at machine precision and hypothesis (a) is now measurement, not arithmetic on a log.** See the step-1b row below. **Review 2026-09-05 10:30:** the band decision is held for the 2026-09-06 weekly, which the 03:00 review had already named as the actor; the recommendation (re-scope the drive-level test to superposed-vs-mean-single-drive residual, and give the 6.7e-05 W systematic its own chunk) is in the `POST-6` §7 entry. Band still 1e-2, test still red, entry still open. |
| **Step 1b, measured** | `docs/testing/logs/20260905T110305Z_POST-6.log`, `1 failed, 21 passed … 102.66s` at `:2034`, Elapsed 104 s at `:2112`, heavy by ceiling (`timeout -k 30 600`), `-n 2`, complex. **`P_acc,k = ½\|a_k\|²(1 − Σ_i\|S_ik\|²)` equals `supplied_k − Σ_i sheets_ik` to 1.683e-15 / 0.000e+00 / 1.679e-15 / 0.000e+00** for P1–P4 (`:1924–1939`) against the pre-registered 1e-6 — the S-matrix's power-wave normalisation and the sheet accounting are the **same number**, so mechanism "the two accountings disagree" is **excluded**. The four single-drive `\|P_acc,k − P_vol,k\|/P_acc,k` read **1.303004e-01 / 1.302723e-01 / 1.300031e-01 / 1.302052e-01**, mean **1.301952e-01**; the superposed ccw residual 1.164806e-01 is **0.8947×** that mean (`:1940`) — the four-port drive misses *slightly less* than its own single drives, so nothing about the superposition adds error. Negative control (drive k's `S_kk` deleted from its column) misses by 7.659732e-01 / 7.656687e-01 / 7.634896e-01 / 7.652279e-01, each equal to its closed-form ceiling `\|S_kk\|²/(1 − Σ_i\|S_ik\|²)` computed from the assembled 4×4 before the assertion, against a claimed floor of 1e-3 (`:1927`, `:1931`, `:1935`, `:1939`). The superposed `P_acc` reproduces step 1's 3.014424803e-03 W at rtol 1e-9. **What is still open:** the ≈6.7e-05 W absolute gap common to every drive is unexplained — 13.0% of accepted power, 0.98% of supplied. Step 1b settles *which denominator the 1% band was registered against*, not *what the gap is*; (b) "a real un-accounted loss channel" remains live at the same absolute size. |
| **`PORT-16` step 1 — mechanism (b) is excluded and the gap is attributed** | `docs/testing/logs/20260907T051231Z_PORT-16.log`, `16 passed … 129.04s` (`:1997`), Status 0 / Elapsed 131 s (`:2065–2066`), standard tier (`timeout -k 30 300`), `-n 2`, complex. Nothing leaks: the **field-level** identity `P_src,exact = P_vol + P_sheet,exact` (the real part of `a(E,E) = L(E)`, exact for the discrete solution) closes at **6.760e-15 / 9.656e-15 / 1.048e-14 / 4.414e-15** on P1–P4 against 1e-6 (`:1882–1893`), and at 1.369e-15 with every sheet's `Z_p` doubled (`:1922`). The ≈6.7e-05 W is therefore **not** an un-accounted loss channel: it is the **Cauchy–Schwarz deficit of the terminal-current sheet form**. `½\|I\|²Re Z_p` is the uniform-field *lower bound* on `½Re(Y_s)∫\|E_t + E_src ĥ\|² dA`, the ratio reads **1.010592–1.010593 on all sixteen** sheet readings (`:1895–1910`), and `C − sheets_terminal` reproduces the gap to **3.212e-13 / 4.770e-13 / 5.086e-13 / 2.276e-13** with `P_vol + C = supplied_terminal` to all ten printed digits (`:1917–1920`). The §9 item's own attribution (iii) is arithmetic and splits the gap into **−54.28×** and **+55.28×** itself (`:1912–1915`) — it attributes nothing. **This entry still does not close**: `POWER_BALANCE_BAND`'s drive-level use is `POST-6` step 2's to retire, and the deliberate red is untouched. |
| **Ruling, 2026-09-07 03:00 review** | `POWER_BALANCE_BAND` is a **record, not a conservation band**: the terminal form `P_acc = ½aᴴ(I − SᴴS)a` misses the volume loss by its own Cauchy–Schwarz deficit (1.0106× the terminal sheet dissipation, reproduced to 1e-13 on every single drive), so this entry's 11.648% is a fixed ≈6.7e-05 W offset divided by a 5× smaller denominator, not a defect in the superposition. **Disposition:** `PORT-16` step 2 (§9 item 2) re-points `test_the_drive_level_power_identity_closes` at the exact identity `P_src,exact(w) = P_vol(w) + P_sheet,exact(w)` (imported `DISCRETE_IDENTITY_RTOL` 1e-6) and its attribution `P_acc(w) − P_vol(w) = C(w) − sheets_terminal(w)` (imported `ATTRIBUTION_RTOL`), with the terminal residual printed beside the band as a record; the band's value and its four green single-drive asserts are not touched. **This entry is retired by step 2's commit, not by the ruling** — the test is red on `main` until then. |
| **CLOSED by `PORT-16` step 2, 2026-09-07 06:00 slot** | `docs/testing/logs/20260907T110826Z_PORT-16.log`, `23 passed … 101.91s` (`:2038`), Status 0 / Elapsed 104 s (`:2106–2107`), standard by measurement, `timeout -k 30 400`, `-n 2`, complex. `test_the_drive_level_power_identity_closes` no longer asserts the terminal form on the ccw / cw quadrature drives. It asserts, on the same superposed field: **(i)** `P_src,exact(w) = P_vol(w) + P_sheet,exact(w)` at the imported `DISCRETE_IDENTITY_RTOL` 1e-6 — ccw `P_src,exact 3.837471142e-03 W`, `P_vol 2.663302665e-03 W`, `P_sheet,exact 1.174168477e-03 W`, **rel dev 5.877e-15**; cw **1.435e-14** (`:1942–1943`); **(ii)** the attribution `P_acc(w) − P_vol(w) = C(w) − sheets_terminal(w)` at the imported `ATTRIBUTION_RTOL` 1e-1 — ccw `P_acc − P_vol 3.511221378e-04 W` against `C − sheets_terminal 3.511221378e-04 W`, **rel dev 9.881e-14**; cw **1.581e-13** (`:1945–1946`), i.e. the *whole* of the drive-level residual is step 1's Cauchy–Schwarz deficit, on a drive none of step 1's solves is. **Control of the generalisation (asserted, 1e-12):** `w = e_k` reproduces step 1's own single-drive `_exact_shares` — `p_src`, `sheet_field_total`, `sheet_ceiling_total` — at **rel dev 0.000e+00 on all twelve readings** (`:1948–1951`). **Nothing widened:** `POWER_BALANCE_BAND` keeps its value 1e-2, its import and its green single-drive asserts (residuals 9.795751e-03 / 9.796209e-03 / 9.794985e-03 / 9.795283e-03, `:1914`); the 1.164806e-01 terminal residual is printed on both senses with the band beside it and asserted nowhere (`:1917`, `:1921`); the blind-sum cross-term control is unchanged. No `src/` change; `test_birdcage_power_identity.py` unedited (no helper signature changed, so rule (c)'s re-run was not triggered — the step-1 import is lazy to avoid a genuine collection-order cycle). |

### ✅ RETIRED 2026-09-03 (`GEO-26` step 3, 22:30 implementer slot) — ~~🔴 OPEN 2026-09-03 (`GEO-26` step 2, 12:00 implementer slot) — the **longitudinal** ring-gap sheet's terminal triangulation is **not C16-covariant** at 16 legs: the terminal disks take exactly **two** discrete areas, 5 of the 16 gap azimuths landing on the low one, so two azimuth classes read an intra-class covariance of **9.990e-05** against the mode's measured 2.0e-5 band~~ — **the 18:00 ruling stood on re-measurement: the 16-leg rung got its own record band 2.0e-4, the reading reproduced to the digit at both widths, and the two-state census (10 of 32 low) is now printed in the log**

| | |
|---|---|
| **Test** | `tests/mesh/test_birdcage_ring_sheet_orientation.py::test_the_longitudinal_ring_sheets_span_the_gap_chord_and_split_the_box[16]` — the `[4]` case (step 1) stays green. A **deliberate red on `main`**, kept so the reading is executable; the band is not widened. |
| **Log** | `-n 2`: `docs/testing/logs/20260903T170351Z_GEO-26.log:53436–53439` (the four classes), failure at `:53448`, `1 failed, 1 passed` / Status 1 / 160 s at `:53453,53470–53471`. `-n 12`: `…170701Z_GEO-26.log:53526–53529`, `:53543`, Status 1 / 158 s at `:53700–53701`. Heavy tier by declaration, 158–160 s measured. Verified at the `GEO-26` step-2 commit. |
| **Symptom** | The four azimuth classes' intra-class terminal-area covariance, **identical to four digits at both widths** (`-n 2` / `-n 12`): `'11.250 deg'` **9.989957e-05 / 9.989957e-05**; `'33.750 deg'` **3.792060e-11 / 3.792088e-11**; `'56.250 deg'` **3.792129e-11 / 3.792129e-11**; `'78.750 deg'` **9.990206e-05 / 9.990206e-05**. The assert that fires is `_assert_ring_identity_family`'s per-class one (`tests/mesh/test_birdcage_ring_gaps_scaleup.py:513`), `9.989956525036291e-05 < 2e-05`. |
| **Cause** | Measured, not inferred: across the 32 ring terminals the meshed area takes exactly **two** values — `9.791961125e-05 m²` and `9.792939386e-05 m²`, differing by 9.99e-05 relative — and the low value is taken by the 5 gap azimuths 11.25 / 78.75 / 101.25 / 191.25 / 281.25 deg (10 of the 32 ports, both rings) and the high value by the other 11. Five of sixteen is not a subgroup of C16, so this is not the azimuth-class effect `GEO-19` step C ruled and `_azimuth_class` folds: it is the **bistable** inscribed triangulation of a disk that step 1 measured as constrained — the longitudinal sheet's two `φ` edges are *diameters* of the two terminal disks (`RING_LONGITUDINAL_CELL_RECORD`'s comment), and gmsh resolves that constrained disk one of two ways depending on where the diameter falls against the surrounding air mesh. At four legs all four gaps happen to land the same way, which is why step 1 read a clean 1.605e-05 and this rung does not. |
| **What is unmoved** | Everything else on the rung, at both widths: 32 sheets at `φ̂`-extent/chord, `ẑ`-extent/`w` and area/(chord·`w`) = 1.000000000000, out-of-plane ≤ 1e-12 m; both halves at their closed forms with the sum at `ring_port_volume_m3`; **C32 sheet spread 6.035e-16 (`-n 12`: 5.998e-16) and top/bottom mirror spread 5.551e-16** against 1e-12; the `GEO-9` partition, the air-box closure, the Pappus arcs and the terminal ratio band all exact; **inter**-class spread inside its 5e-3 ceiling; and the negative control — the default orientation at 16 legs — reproducing `EX-35`'s **265 621 cells at ratio 1.000000** with all 32 transverse sheets at `φ̂`-extent ≤ 1.741094e-17 m (the largest of the 32) (`…170351Z:26462–26494`, `…170701Z:26542`). The 16-leg longitudinal cell count is **270 728**, identical at both widths. |
| **Disposition** | `LONGITUDINAL_TERMINAL_INTRA_BAND` stays **2.0e-5** and `TERMINAL_INTRA_CLASS_BAND` stays 1e-6 — the 10:30 review pre-registered a class above the band as a **stop, never a widening**. `PORT-13` stays 🚫: its step 1 needed this rung's record, and while the record exists (270 728) the rung it was measured on carries an unadjudicated generator finding. A review owns the next move — either rule the bistable triangulation acceptable (it is a terminal-*area* reading, every terminal is inside [0.95, 1.0], and `PORT-13` integrates over the **sheet**, which is exact to 1e-16 here) and re-register the band with this measurement, or classify by the actual two-valued partition instead of `_azimuth_class`. Expected on `main`; **not yours**. **Ruled 2026-09-03 18:00 review — the bistable triangulation is acceptable and the 16-leg rung gets its own record band; the `_azimuth_class` fold is kept.** Grounds: the quantity is a terminal-*area* covariance, a mesh-reproducibility record, not a physics band or a closed-form comparison; every one of the 32 terminals is inside its [0.95, 1.0] closed-form band; the sheet `PORT-13` integrates over is exact to 1e-16 on all 32; the 9.99e-05 two-state amplitude is a decade under the 1e-3 reciprocity band and three under the 5% C4 class band, and the mechanism is measured (constrained-diameter triangulation, two discrete areas, 5 of 16 azimuths) rather than inferred. Classifying by the two-valued partition was rejected: it would fit the assertion to the artifact and hide a third state if one ever appeared. The 4-leg band stays **2.0e-5** (step 1's own measurement, unmoved); the fix is `GEO-26` **step 3** (§9 item 3): a version-tagged `LONGITUDINAL_TERMINAL_INTRA_BAND_16 = 2.0e-4` for the `[16]` case only, with all four class readings and the two areas in the constant's comment (the `RING_LONGITUDINAL_CELL_RECORD` pattern), asserted at both widths, plus the 5-of-16 low-state count printed so a change of state is visible. This entry retires with that commit. `PORT-13` step 1 is unblocked by the ruling alone (§9 item 4) — it reads the record and the sheet, not this test's colour. |
| **Retired by (`GEO-26` step 3, 2026-09-03 22:30 slot)** | The prescription ran exactly as written. `LONGITUDINAL_TERMINAL_INTRA_BAND_16 = 2.0e-4` registered version-tagged (0.11 image) in `tests/mesh/test_birdcage_ring_sheet_orientation.py` with all four class readings, both terminal areas and the five low-state azimuths in its comment, passed for the `[16]` case only; `LONGITUDINAL_TERMINAL_INTRA_BAND` (2.0e-5) and `TERMINAL_INTRA_CLASS_BAND` (1e-6) unmoved, no `src/` change, no helper renamed. **Negative control first**, footered: the `[16]` case pointed back at 2.0e-5 reproduces this entry's failure exactly — `assert np.float64(9.989956525036291e-05) < 2e-05`, `1 failed, 3 passed`, Status 1, 226 s (`20260904T033238Z_GEO-26.log:26282,26286`) — then reverted. **Green at both widths**, `4 passed` / Status 0, heavy by measurement: `-n 2` `20260904T034052Z_GEO-26.log` (222 s), `-n 12` `20260904T034441Z_GEO-26.log` (224 s). The four classes read 9.989957e-05 / 3.792060e-11 (`-n 12` 3.792088e-11) / 3.792129e-11 / 9.990206e-05 against 2.0e-4, and the newly printed state census reads **2 states, 9.791961125e-05 m² ×10 and 9.792939386e-05 m² ×22 — 10 of 32 on the low area** at both widths, the pre-registered expectation to the count. Every step-2 anchor unmoved (32 sheets at 1.000000000000 on the three ratios, both halves at closed form, C32 6.035e-16 / 5.998e-16, `RING_LONGITUDINAL_SCALED_CELL_RECORD` 270 728 at ratio 1.000000 — reachable for the first time — the 265 621 and 110 786 controls at 1.000000, and `[4]` at step 1's digits). **One correction to this entry's Cause row:** the state census shows the **4-leg** rung is two-state as well (6 at 9.793917647e-05 m², 2 at 9.794074883e-05 m², `…034052Z:15663`), so "at four legs all four gaps happen to land the same way" is wrong; the two states there are 1.605e-05 apart rather than 9.99e-05, which is why that rung reads inside 2.0e-5. Nothing about the disposition changes — the comment in the code now says what was measured. |

### ✅ RETIRED 2026-09-03 (`OPS-33`, 00:00 implementer slot) — ~~🟡 OPEN 2026-09-02 (`OPS-32` slot) — `ans:3`'s four reproduction controls cannot be re-registered at **1e-6**: three of the four sit at 1.9e-06 … 3.0e-05 against the `PORT-1` step-4 **log record**, which is not a run-to-run comparison~~ — **the diagnosis was right and the fix was the record, not the band: all four `RECORDED_*` re-based onto the 0.11 image, then 1e-6 registered (relative on raw / corrected / ‖S‖₂, absolute on the symmetry residual) with the superseded digits asserted to fail it**

| | |
|---|---|
| **Log** | `docs/testing/logs/20260902T183603Z_OPS-32.log:1405` — `[ANS-3] reproduction of the PORT-1 step-4 record (band 1% relative): raw 2.98e-05, corrected 2.92e-05, symmetry 1.91e-06, ‖S‖_2 3.07e-10`, Status 0, 171 s. Verified at the `OPS-32` commit. |
| **Symptom** | `OPS-32` was scoped (2026-09-02 weekly, §9 item 2) to re-register both benchmark scripts' in-script reproduction controls at 1e-6 against `EX-37`'s measured ≤ 5e-8 run-to-run Z/S scatter. `ans:4` already carries 1e-6 / 1e-9 (`FREQUENCY_CONTROL_BAND`, `LEG_D0_REPRODUCTION_BAND`; measured 1.157e-10 and 2.568e-10 in `…183909Z_OPS-32.log:4801`) and needed nothing. `ans:3` at 1e-6 would fail on three of four entries. |
| **Cause** | Not scatter. `ans:3`'s misses are taken against `RECORDED_*` constants transcribed from `20260813T183606Z_PORT-1-step4-packagegate.log`, a different code version and a different image; the 2.9e-05 raw/corrected offset is a standing record-vs-today gap, and the symmetry entry moved 1.71e-06 → 1.91e-06 between two runs of the same code (the residual is ≈4.8e-05 absolute, so its *relative* miss amplifies scatter). Re-basing those constants is an `OPS-27`-class record job, not a band tightening. |
| **Disposition** | The band stays at `EX-20`'s 1%, with the measurement recorded in a code comment at `examples/ansys_benchmarks/ans3_two_torus_gap_ports_10MHz/03_two_torus_gap_ports_10MHz.py:113`. The writer half of `OPS-32` landed; the control half is this finding. A follow-up would re-record `RECORDED_RAW_RATIO` / `RECORDED_CORRECTED_RATIO` / `RECORDED_S_SYMMETRY_RESIDUAL` from a current run under the (1*) licence **and only then** register a 1e-6 control. Expected on `main`; **not yours**. |
| **Retired by (`OPS-33`, 2026-09-03)** | The follow-up ran exactly as prescribed. All four `RECORDED_*` in `examples/ports/02_package_sparameter_sweep.py` re-based onto the 0.11 image under the in-class (1\*) licence (0.8945163788281 / 0.9398215452105 / 4.7586341120262e-05 / 0.8648094567341), the superseded v0.7.2 digits kept as `SUPERSEDED_V072_RECORD`, and the control registered at **1e-6 relative** on raw / corrected / ‖S‖₂ plus **1e-6 absolute** on the symmetry residual — absolute because that entry is itself ≈4.8e-05, so `EX-37`'s 5e-8 S-entry scatter is ≈1e-3 of it *relatively* and a relative 1e-6 is arithmetically unreachable there (this entry's own 1.71e-06 → 1.91e-06 wobble). Measured green: `ans:3` `20260903T050551Z_OPS-33.log:1409` raw 0.00e+00, corrected 4.63e-14, symmetry 5.99e-14 abs, ‖S‖₂ 7.57e-15 (Status 0, 165 s); `ports:2` `20260903T050845Z_OPS-33.log:1425` raw 6.51e-10, corrected 6.39e-10, symmetry 1.87e-10 abs, ‖S‖₂ 6.64e-11 (Status 0, 206 s) — the second run is a genuinely independent solve and still sits ≈3 decades inside the band, so the 1e-6 figure is a scatter control, not a tautology. Negative control asserted in both: the superseded digits miss the same run by 2.98e-05 / 2.92e-05, three decades outside. No physics band or gate moved. |

### ✅ RETIRED 2026-09-02 (`WF-6` step 3h, 19:30 implementer slot) — ~~🔴 OPEN 2026-08-31 (`WF-6` step 3, 13:30 implementer slot) — the coil-driven **point-SAR** map misses every C4 / mirror identity by **25–40%** against the same 5% band the `|B₁⁺|` map meets at ≈2%: the pointwise `|E|` estimator has its own, much larger floor, and nobody had measured it~~ — **the construction was the mechanism; the identity is gated as cell integrals and the five pointwise asserts are records**

> **✅ RETIRED 2026-09-02, 19:30 implementer slot (`WF-6` step 3h) — no band
> was ever moved and none was widened to retire this.** The five asserts in
> `tests/validation/test_birdcage_sar_map.py` now assert their step-3
> `STEP3_PRIMAL_IDENTITY_RECORDS` (25.1096 / 40.5462 / 30.0142 / 38.6120 /
> 28.1459%) at rtol 1e-3 **and** that each reading still *exceeds*
> `C4_COVARIANCE_BAND` — the pointwise construction retired as a gate with its
> measurement kept and its miss asserted, so a change that silently made it
> pass is visible rather than absorbed. The identity itself is gated where six
> rungs said it belonged: as **cell integrals of the primal `E`** in
> `tests/validation/test_birdcage_sar_integral.py`, twelve C4 pairs and the
> quadrature four-quadrant spread asserted `≤ 5%` (imported, unmoved) and
> against step 3g's table at rtol 1e-3 — 0.7149 … 1.5200%, spread 0.4641%,
> partition identity exact, mis-paired control 89–159× larger. `main` is
> **`0 failed`** on both modules: `20260903T003309Z_WF-6-step3h.log`,
> **`109 passed` / Status 0 / 194 s** at `-n 2` complex with
> `tests/environment` (`test_birdcage_sar_integral.py` 22 items,
> `test_birdcage_sar_map.py` 76 items, no other count moved). **Retired-by
> row: `WF-6` step 3h, this commit; ruled 2026-09-02 18:00 review, §7 `WF-6`
> rulings blockquote, last paragraph.** The gate this registers is exactly
> one thing — a C4 symmetry identity of quadrant powers on one fixture at
> 10 MHz at fixed `h` — with no mirror identity, no absolute SAR, no
> homogeneity, no C95.3, no Larmor and no convergence claim, and `WF-6`
> stays 🟡.
>
> The historical record below is kept verbatim: it is the only place the six
> rungs that excluded the projector, the estimator's degree and (at fixed `h`)
> the mesh are written down.

> **Five deliberate reds on `main`** (retired above), all in the new
> `tests/validation/test_birdcage_sar_map.py`:
> `test_single_drive_sar_map_is_c4_covariant` at all three parametrised angles,
> `test_quadrature_sar_map_is_c4_invariant`, and
> `test_reversing_the_rotation_sense_equals_reflecting_the_sar_map`. Kept red on
> purpose, per the step's own pre-registered negative-result clause ("an identity
> missing is an estimator finding exactly like step 1's DG0 floor — known-issues
> with every reading, keep the assert, stop"). **No band was moved and none may
> be moved in-slot.**

| | |
|---|---|
| **Verified at** | `75f60bc` + this commit, 0.11 image, complex build, `-n 2`, 2026-08-31 — `20260831T183526Z_WF-6-step3.log`, `5 failed, 16 passed` / Status 1 / **96 s**, with `tests/environment`. |
| **Literal symptom** | Relative ℓ² of `σ\|E\|²/(2ρ)` over the 51 phantom centroids, band **5.0%** imported unmoved from step 1d: `SAR_P2(Rx)` vs `SAR_P1(x)` **25.1096%**; `SAR_P4(−Rx)` **40.5462%**; `SAR_P3(180°)` **30.0142%**; quadrature `SAR_ccw(Rx)` vs `SAR_ccw(x)` **38.6120%**; mirror `SAR_cw(Mx)` vs `SAR_ccw(x)` **28.1459%**. |
| **What is *not* wrong** | Every control and cross-check in the same run passes. The two negative controls hold with room to spare — mis-rotated `SAR_P3(Rx)` vs `SAR_P1(x)` **129.8187%**, quadrature vs single-drive **334.5786%**, both asserted `> 5%`. `mean_sar`'s phantom `dissipated_power_w` on P1 reads **5.637745667e-08 W**, reproducing step 1's gate-(i) record *to every printed digit* — the same `σ\|E\|²` integrated instead of sampled, so this module is demonstrably on step 1's solve. The phantom σ premise is asserted and holds: **0.5 S/m** flat over 537 tag-3 cells. `point_sar` raises rather than zero-filling, so all 51 points of every rotated and mirrored image set were evaluated. |
| **Cause (diagnosed, on this run's own evidence)** | The **estimator**, not the field and not the new arithmetic. The three single-drive readings (i) use *only* step 1's four solved fields and step 1d's own image sets — no superposition, no mirror, no code this step introduced — and they already miss by 25–40%. What changed against the B₁⁺ legs is what is read: `|B₁⁺|` comes from an L²-projected CG1 `B` (step 1d's ruling, 2.19 / 2.11 / 1.89%), while SAR is read pointwise off the **primal N1curl `E`** through its DG interpolant, with no projection. A lowest-order Whitney `E` is tangentially continuous and normally discontinuous, so its *pointwise* value at a cell centroid carries an O(h) per-cell error that the C4 rotation does not share; squaring it doubles the relative error. 25–40% is ≈ 2× a ≈13–20% pointwise `|E|` floor — the same shape of finding as step 1's 8.65% DG0 curl, one estimator further out. |
| **A degeneracy worth recording** | The mirror-omitted comparison `SAR_ccw(Mx)` vs `SAR_ccw(x)` reads **28.1445%** against identity (iii)'s **28.1459%** — the mirror moves the reading by 1.4e-3 pp. The 10:30 scoping had listed the mis-paired sense as a *control*; it is not one for SAR (a magnitude-squared has no ± senses to mis-pair), and the module therefore prints it ungated and asserts the two controls that do separate instead. Documented in the module docstring; the disposition is the review's. |
| **Recorded, ungated** | Point SAR at 1 V per port: P1 single drive peak **7.630679e-07** W/kg, mean **1.453536e-07** W/kg, peak/mean **5.2497**; quadrature (ccw) peak **2.065442e-06** W/kg, mean **7.706353e-07** W/kg, peak/mean **2.6802**. Not a safety figure — no converged mesh, no real drive, no mass averaging. |
| **Consequence** | `WF-6` stays **🟡**: steps 1, 2 and 2b remain ✅ and **no SAR claim of any kind exists** — not homogeneity, not absolute, not C95.3. The `|B₁⁺|` gates are untouched (step 2's module re-run green and unmoved this slot, `20260831T183734Z_WF-6-step3-step2-regression.log`, `6 passed` / Status 0 / 80 s, identities 0.9818 / 0.8087%, control 95.1975%). |
| **Resolves with** | A `WF-6` step-3b estimator leg on the step-1b pattern, which a review must scope: read the *same* identities off an L²-projected CG1 `E` (`post.project_to_cg1` already takes any vector field) beside the primal column, on the same 51 points, and decide between (a) a pointwise-`E` estimator floor — CG1 lands inside 5% while the two controls survive the projection — and (b) something the rotation does not share in the field itself, which the 180° column would show. Cheap: no new solve, ≈ the same 96 s. Until then the five asserts stay exactly as written. **Scoped 2026-08-31 18:00 review as `WF-6` step 3b, §9 item 1** — the primal readings above become asserted records, the CG1 readings are printed beside them with a pre-registered three-way verdict ((a) estimator floor / (b) rotation-specific field feature, the 180° column deciding / (c) the fixture's cells do not resolve a quadratic-in-`E` map), and no band moves in-slot; whichever verdict prints is journaled here. The mis-paired-sense comparison is struck from SAR scopings. |
| **Step 3b executed 2026-08-31 19:30 slot — verdict (c) prints, and a diagnostic says (c)'s *reason* is not the one pre-registered: `project_to_cg1` does not fit the N1curl `E` at all** | `20260901T003300Z_WF-6-step3b.log` and, with the added diagnostic, `20260901T003548Z_WF-6-step3b-diagnostic.log`, both `5 failed, 25 passed` / Status 1 / **105 s** and **100 s** at `-n 2` complex with `tests/environment`; the 5 failures are this entry's five primal asserts, unmoved. **Anchor held to every digit:** the primal column reproduced step 3's five identities (25.1096 / 40.5462 / 30.0142 / 38.6120 / 28.1459%) and both controls (129.8187 / 334.5786%) at rtol 1e-3 — nothing about the fixture moved. **The CG1 column (printed, not gated) is *worse* than the primal one at every identity:** 152.0459 / 109.7797 / 169.5050 / 53.1869 / 40.8440% against primal 25.11 / 40.55 / 30.01 / 38.61 / 28.15%; CG1 controls 163.6144 / 75.9135%, both asserted `> 5%` and both holding. So the pre-registered verdict evaluates to **(c)**, all five outside. **But (c)'s pre-registered reading — "the ≈1 cm phantom cells do not resolve a quadratic-in-`E` map" — is contradicted by the two projection diagnostics in the same run:** the CG1 phantom power `½∫σ\|E_cg1\|²` reads **1.990062891e-05 W** against the primal record 5.637745667e-08 W (**+35 198.9%**, where step 1d's `B` projection moved its mean by 0.38%), and `‖E_cg1 − E‖/‖E‖` over the phantom reads **1876.1871%**. A projection that fits its argument cannot be 19× larger than it in the very region being read. The measured finding is therefore about **`post.project_to_cg1` applied to an N1curl `E`**, not about the phantom's resolution: the estimator step 1d registered for the DG0 `B` phasor is not, as-landed, a usable `E` estimator on this fixture, and no SAR conclusion — floor, resolution or otherwise — can be drawn until that is diagnosed. **No band moved, no assert loosened, no SAR claim exists**, `WF-6` stays 🟡. The mis-paired-sense comparison stayed struck; this module carries no step-2 numeric literals, so the scoped `STEP2_*` collateral was a no-op. |
| **Resolves with (step 3b's own finding)** | A diagnosis of `project_to_cg1` on N1curl input, which a review must scope — the candidates the run does not separate: (1) the global L² fit is dominated by the sheet/conductor-edge `E` singularities, so the phantom sees a mass-matrix tail rather than the local field, in which case an `E` estimator has to be restricted to the phantom subdomain or projected cellwise; (2) the CG (Jacobi) mass solve at `ksp_rtol` 1e-12 does not converge on a vector CG1 space over this 116 085-cell mesh, which the helper does not check and no caller has ever asked it; (3) an element-side mismatch that value-shape `(3,)` does not catch. The three are cheaply separable in one slot: print the KSP converged reason and iteration count, and re-read `‖E_cg1 − E‖/‖E‖` over the *whole* mesh beside the phantom-restricted figure. Until then the CG1-`E` column is a measurement of the projector, not of SAR, and the `\|B₁⁺\|` gates are untouched by it (they project `B`, whose 0.38% figure is unchanged). |
| **Step 3c executed 2026-09-01 07:30 slot — the projector *is* a projector: candidates 2 and 3 are refuted, candidate 1 is measured** | `20260901T123421Z_WF-6-step3c.log`, `5 failed, 38 passed` / Status 1 / **103 s** at `-n 2` complex with `tests/environment`; the 5 failures are this entry's five primal asserts, unmoved to the digit. Every step-3b record reproduced at rtol 1e-3 and is now asserted (primal and CG1 identity columns, both control sets, CG1 phantom power 1.990062891e-05 W, the 1876.1871% phantom residual). **(1) The mass solve converges:** PETSc `converged_reason` **2** (`KSP_CONVERGED_RTOL`) in **26** iterations on **64 191** CG1 dofs at `ksp_rtol` 1e-12 — candidate 2 refuted; the solver is now readable through an opt-in `return_diagnostics` kwarg on `post.project_to_cg1` (default off, `B` callers untouched). **(2) The projector reproduces exactly what `CG1³` contains:** `f = a + b × x` (in `N1curl₁ ∩ CG1³`), interpolated complex into the solve's own N1curl space and projected, leaves `‖P f − f‖/‖f‖` = **1.326607e-13** against the 1e-10 anchor; the control's control `x² ê_x` leaves **9.882703e-02** (> 1e-3) — candidate 3 refuted. **(3) The domain table names candidate 1:** `‖E_cg1 − E‖/‖E‖` = **32.7802%** whole mesh, **1876.1871%** phantom, **838.8978%** phantom core (33 of 537 owned tag-3 cells with no vertex on the phantom boundary). The whole-mesh figure is 57× *below* the phantom's, so the global L² fit is a fit of the sheet / conductor-edge `E` that dominates `‖E‖`, and the phantom — orders of magnitude lower `\|E\|` — gets that fit's tail; excluding the σ-interface layer only halves the phantom figure, so the interface smear is secondary, not the mechanism. **No band moved, no assert loosened, nothing re-registered, no SAR claim exists**, `WF-6` stays 🟡. |
| **Resolves with (step 3c's finding)** | Nothing is wrong with `post.project_to_cg1`; the **use** is wrong. A *global* L² projection is not an `E` estimator inside a low-field subdomain of a fixture that carries a huge-field region, so `WF-6` **step 3d** — a review's to scope — must build the `E` estimator either on a **phantom submesh** (project on the subdomain the SAR map is read over, where the fit is not dominated by fields outside it) or **cellwise** (a per-cell L² projection, no global coupling), then re-read the same five identities and two controls beside the primal column with the same pre-registration discipline. The `\|B₁⁺\|` gates stay untouched under either: they project `B` on the whole mesh where no comparable field-magnitude contrast exists, and their 0.38% figure is unchanged. Until step 3d lands, the five primal asserts stay exactly as written and no SAR claim of any kind exists. **Scoped 2026-09-01 10:30 review as `WF-6` step 3d, §9 item 9** — the phantom-restricted route, on the parent mesh (mass matrix and load integrated over tag 3 only, CG1 dofs with no phantom support pinned to zero), with the best-approximation inequality `‖P_Ω E − E‖_Ω ≤ ‖P E − E‖_Ω` (the 1876.1871% record) as its asserted identity, the exact-reproduction control re-run under the restriction, and the five identities read off it beside the primal column under 3b's pre-registered verdict; the cellwise route is struck by derivation (a per-cell L² projection onto DG1³ reproduces a degree-1 N1curl field exactly, so it *is* the primal column). |
| **Step 3d executed 2026-09-01 13:30 slot — the phantom-restricted estimator is honest and the identities still miss: verdict (c), now uncontradicted** | `20260901T183416Z_WF-6-step3d.log`, `5 failed, 52 passed` / Status 1 / **123 s** at `-n 2` complex with `tests/environment`; the 5 failures are this entry's five primal asserts, unmoved to the digit (25.1096 / 40.5462 / 30.0142 / 38.6120 / 28.1459%). Every 3b/3c record reproduced at rtol 1e-3, and step 3c's other two domain figures (whole mesh **32.7802%**, phantom core **838.8978%**) are now asserted too. The estimator is a test-local `_project_to_cg1_restricted` on the **same** `("Lagrange", 1, (3,))` space and the **parent** mesh (no submesh, no cross-mesh N1curl interpolation): mass and load over `dx(3)`, every CG1 dof with no phantom-cell support pinned to zero through a `dirichletbc` from a zero `Function` over owned **and ghost** blocks — **170** free of 21 397 owned blocks, pinned max \|value\| **exactly 0.000e+00** at `-n 2`, all six restricted solves `converged_reason` **2** in 21–25 its on 64 191 dofs. **(i) The best-approximation inequality holds with a 100.20× margin:** `‖P_Ω E − E‖_Ω/‖E‖_Ω` = **18.7238%** against the global fit's 1876.1871% on the same phantom — 3c's diagnosis confirmed constructively, the 1876% *was* the global fit's tail. **The restricted phantom power corroborates independently:** `½∫_Ω σ\|E_Ω\|²` = **5.440097168e-08 W**, **−3.51%** from the primal record 5.637745667e-08 W, against the global projection's **+35 198.9%**. **(ii) The restriction is a projection:** `a + b × x` restricted-projects to **4.385695e-13** (anchor 1e-10) and `x² ê_x` to **3.741459e-01**, above the arithmetic `(h/D)² ≳ 1e-4` floor. **The five identities (printed, NOT gated) improve 3–6× to 8.2868 / 9.4743 / 7.3477 / 6.8146 / 6.1185%** (primal 25.11–40.55%, global CG1 40.84–169.51%) with both restricted controls surviving and asserted (**123.6255%**, **333.0778%**) — **and still miss the 5% band at all five**, so the pre-registered verdict evaluates to **(c)**. Unlike step 3b, no diagnostic in the run contradicts (c)'s pre-registered reading: with the projector exonerated (3c) and the global-fit tail removed (this step), the residual 6–9.5% is what an *honest* estimator leaves, i.e. the fixture's ≈1 cm phantom cells reading a quadratic-in-`E` map. **No band moved, no assert loosened, no SAR gate registered, nothing under `src/`, no SAR claim exists**, `WF-6` stays 🟡. |
| **Resolves with (step 3d's finding)** | A **finer phantom rung** — the same five identities on a mesh whose phantom resolves `\|E\|` better, read off the restricted estimator this step landed, to see whether 6–9.5% falls toward the `\|B₁⁺\|` map's ≈2.2% floor. That is a **weekly-review** question, not a daily one: it is a new mesh and a new four-drive sweep, so it must be costed against `TH-11`'s memory wall before anyone queues it, and the outcome decides whether a SAR gate can ever be registered on this fixture family. Two cheaper things a daily review may scope first, in either order: (1) **promote `_project_to_cg1_restricted` into `post/`** now that its best-approximation anchor, its exact-reproduction control and its ghost-safe pinning are all measured and asserted — it is the only `E` estimator this repo has that reproduces the phantom power to 3.5%; (2) **re-register nothing** — the five primal asserts stay exactly as written and red until a rung actually lands, because the restricted column is printed, not gated, by design. The `\|B₁⁺\|` gates remain untouched throughout (they project `B` on the whole mesh, 0.38%). **Both of (1)'s halves scoped 2026-09-01 18:00 review as `WF-6` step 3e (§9 item 1): the promotion into `post/` as `project_to_cg1_restricted`, its anchors being step 3d's records reproduced through the packaged path with the global projector's 1876.1871% as the asserted negative control, plus the docstring warning on `post.project_to_cg1` that a global L² fit is not an `E` estimator inside a low-field subdomain — 3c's finding, which until then lives only here.** The same review also scoped **step 3e′ (§9 item 4)**, which is *not* the finer-mesh rung: nothing on this fixture has separated **estimator degree** from **mesh h**, and a CG2³ restricted projection of the same four solved fields on the same mesh separates them for six mass solves and no curl-curl solve. Its anchors are theorems — the restricted best-approximation residual cannot increase with degree (`≤ 18.7238%`), and `x² ê_x` lies in `CG2³` exactly, so where CG1 left 3.741459e-01 the CG2 fit must reproduce it to ≤ 1e-10, a nine-decade pre-registered separation. Its expected outcome is the *null* one (the five identities move by < ≈1 pp), which corroborates (c) and leaves the finer-mesh rung — still the **2026-09-02 weekly review's**, still costed against `TH-11`'s memory wall — as the remaining candidate. **Step 3e′ executed 2026-09-02 06:00 slot and printed (γ), with the pre-registered cause of (γ) excluded by its own anchors — the result nobody scoped for.** All six anchors green on the first run: the CG2³ restricted fit is *strictly better* in the norm it minimises (`‖P²_Ω E − E‖_Ω/‖E‖_Ω` **14.4724%** vs CG1's 18.7238%, the degree-monotonicity theorem asserted), both fields it contains come back to ≈1e-12 (`a + b × x` **1.363313e-12**, `x² ê_x` **1.505524e-12** — the pre-registered flip landed, the *same* exact source reading **6.659346e-02** at degree 1, a separation of **10.6 decades**), the pin is exactly **0.000e+00** over owned and ghost blocks (1 004 free of 160 537 owned CG2 blocks, 481 611 dofs, vs CG1's 170 / 21 397 / 64 191), all six CG2 solves reason **2** in 39–48 its, both controls survive (**123.2927% / 327.6543%**), and every CG1 record of steps 3d/3e reproduced unmoved in the same window. And yet the five identities got **worse**: **19.3491 / 17.2097 / 16.0699 / 14.4087 / 11.3230%** against CG1's 8.2868 / 9.4743 / 7.3477 / 6.8146 / 6.1185% (+11.06 / +7.74 / +8.72 / +7.59 / +5.20 pp), phantom power 5.519662942e-08 W (−2.09% from the primal record, better than CG1's −3.51%). So a **globally better L² fit of `E` over the phantom is pointwise worse for these C4 identities at these 51 points** — the CG1 fit's extra smoothing was flattering the identities, not resolving them, and neither "estimator degree" nor "the projector" is the mechanism verdict (c) named. Nothing was rescoped in-slot: no band moved, no gate registered, the five primal asserts stay red, `WF-6` stays 🟡. **What now resolves this** is a review ruling on the *construction*, not another rung: an identity read from a fitted field inherits the fit's pointwise error, so either (a) the identities must be read as **integrals** over the phantom (a C4-covariance of `∫σ|E|²` over rotated subregions, which the power figures suggest is stable to ≈2%) rather than at 51 sampled points, or (b) the 51-point sample set must be shown to be the mechanism (step 1c did exactly this for the `|B₁⁺|` DG0 column and found it was *not*), or (c) the finer-mesh rung (step 3f) runs anyway and is read knowing that a better fit may *raise* these numbers. Logs `20260902T110503Z_WF-6-step3e-prime.log` (`5 failed, 82 passed` / Status 1 / 125 s) and `20260902T111000Z_WF-6-step3e-prime-verdict.log` (the printed table, 121 s). **Step 3g executed 2026-09-02 16:30 slot and took option (a) of that list — read the identities as INTEGRALS — and it prints clause (a) at *fixed* `h`.** `tests/validation/test_birdcage_sar_integral.py` on the **default** 116 085-cell mesh (this entry's own mesh, no `phantom_resolution`) forms `P_j^{(k)} = ½∫_{tag 3} σ w_j |E^{(k)}|²` with the smooth azimuthal partition `w_j = ((c_j + √(c_j²+ε²))/2)²`, `c_j = (x cos θ_j + y sin θ_j)/√(x²+y²+ε²)`, `ε = 1e-9 m` (a regularised `sqrt`, because `ufl.max_value` is the `OPS-22` complex-build trap), off the **primal** N1curl `E` with no projection, no estimator and no sample set. The twelve C4 pairs `\|P_{j+1}^{(k+1)} − P_j^{(k)}\|/P_j^{(k)}` read **0.7149 / 1.1908 / 1.4417 / 0.9703 \| 1.5200 / 0.2132 / 0.3377 / 0.9086 \| 1.0569 / 0.8355 / 0.2780 / 0.2302%** — worst **1.5200%**, all inside the unmoved 5% band, against **this entry's own pointwise primal 25.11–40.55%** and the pointwise restricted-estimator column's 6.1–9.5% on the *same four solves and the same mesh*; the quadrature drive's four-quadrant spread is **0.4641%**. Asserted and green: the partition identity `Σ_j P_j^{(k)} = P_phantom^{(k)}` to every printed digit for all five drives at rtol 1e-10, the gate-(i) drive's total **5.637745667e-08 W** (step 1's record, digit-for-digit), and the mis-paired 180° control strictly larger at **96.1655 / 97.4944 / 95.5869%**, ratios **89.09 / 130.89 / 159.27×**. `20260902T213441Z_WF-6-step3g.log`, `20 passed` / Status 0 / 106 s at `-n 2` complex. **So the mechanism this entry has been hunting since step 3 is the CONSTRUCTION, not the code, the projector, the estimator degree or (at fixed `h`) the mesh: a quadratic identity read from a field at sampled points inherits that field's pointwise error squared, and the same solves that miss by 25–41% pointwise satisfy the C4 statement to ≤ 1.52% as integrals.** This entry stays **OPEN** and the five primal asserts stay exactly as written and red: nothing was re-registered, no band moved, **no SAR gate exists**, and `WF-6` stays 🟡. **What now resolves it is a review ruling on which construction the repo gates** — the integral one at 1.52% on the coarse mesh (step 3g) or the pointwise restricted-estimator one at 2.5–3.5% on the finer phantom (step 3f) — plus, either way, a decision on what to do with five pointwise asserts the repo has now measured to be the wrong statement of the identity. Neither is a slot's call. |
| **Step 3e executed 2026-09-02 00:30 slot — the restricted estimator is packaged, and the move moved nothing** | `20260902T003813Z_WF-6-step3e-table.log`, `5 failed, 69 passed` / Status 1 / **122 s** at `-n 2` complex with `tests/environment` in the same window (module-only re-run `20260902T003518Z_WF-6-step3e.log`, `5 failed, 58 passed` / Status 1 / 98 s; `tests/environment` alone `20260902T003443Z_WF-6-step3e-env.log`, `11 passed` / Status 0 / 29 s). The 5 failures are this entry's five primal asserts, **unmoved to the digit** (25.1096 / 40.5462 / 30.0142 / 38.6120 / 28.1459%). Step 3d's test-local `_project_to_cg1_restricted` moved **verbatim** into `src/fem_em_solver/post/faraday.py` as **`post.project_to_cg1_restricted`**, signature `(field, cell_tags, *, name, tag, ksp_rtol=1e-12, return_diagnostics=False)` — `return_diagnostics` default flipped **off** to match `project_to_cg1`, PETSc prefix renamed off the step-scoped name to `fem_em_restricted_cg1_mass_`, exported from `post/__init__.py` and both `__all__`s, with the test module its first caller at `return_diagnostics=True`. **17 new test items, all green, every one a step-3d record re-measured through the packaged path:** `‖P_Ω E − E‖_Ω/‖E‖_Ω` **18.7238%**; the **negative control** — the *global* `post.project_to_cg1` on the same field over the same phantom, read in the same run at **1876.1871%** — separating by **100.20×** against the pre-registered **50×** floor (a promotion that kept integrating over `dx` would read 1.0×); `a + b × x` **4.385695e-13** and `x² ê_x` **3.741459e-01**, bounds 1e-10 / 1e-4 unmoved; pinned max \|value\| **0.000e+00** over owned **and** ghost blocks, census asserted at **170** free of **21 397** owned blocks on **64 191** dofs (globally reduced, rank-count independent); six solves `converged_reason` **2** in **25 / 25 / 25 / 25 / 21 / 25** its, now asserted `== 2` and `21 ≤ its ≤ 25` beside the pre-existing `> 0`; restricted phantom power **5.440097168e-08 W** (−3.5058% from the primal record); the five identities **8.2868 / 9.4743 / 7.3477 / 6.8146 / 6.1185%** and both controls **123.6255%** / **333.0778%**, all at `CG1_RECORD_RTOL`. **Collateral, closing the second half of the step-3d ask:** `post.project_to_cg1`'s docstring now carries 3c's finding as a `.. warning::` — a *global* L² fit is not a field estimator inside a low-field subdomain — with the 32.7802 / 1876.1871 / 838.8978% domain table and a pointer at the restricted sibling, so it no longer lives only in this file. **This packages an estimator; it does not register a gate.** The five primal asserts stay exactly as written and red, the module still exits 1, no band, tolerance or record moved, nothing was re-registered, no SAR claim exists, verdict (c) is unaffected, and `WF-6` stays **🟡**. `B` callers are untouched by construction. This entry stays **OPEN**: the finer-phantom rung above is still the resolution. |
| **Step 3f executed 2026-09-02 15:00 slot — the finer-phantom rung prints (a): all five identities land INSIDE the 5% band, and verdict (c)'s attribution to the phantom's `h` is confirmed** | `20260902T170559Z_WF-6-step3f.log`, `3 failed, 27 passed` / Status 1 / **175 s** (pytest 173.01 s) at `-n 2` complex with `tests/environment`, standard tier against the item's 150–200 s estimate and its 600 s ceiling. New module `tests/validation/test_birdcage_sar_fine_phantom.py`: the whole four-drive sweep rebuilt through `build_four_port_sweep(phantom_resolution=0.0075)` — **no `reuse=`**, which hands back the coarse mesh — and every column of step 3e re-read on it. **(iii) The knob reached the constructor:** **120 499** cells and **2 746** tag-3 cells, asserted at exact equality against step 3f₀'s measurement (coarse 116 085 / 537, 5.1136×). **The five identities off `post.project_to_cg1_restricted` (degree 1), printed not gated: 3.3600 / 3.4442 / 3.4525 / 3.0332 / 2.5465%** against the coarse mesh's 8.2868 / 9.4743 / 7.3477 / 6.8146 / 6.1185% — **−4.93 / −6.03 / −3.90 / −3.78 / −3.57 pp, ratios 0.4055 / 0.3635 / 0.4699 / 0.4451 / 0.4162** (mean 0.42, against the 0.25 a purely second-order residual would give under one halving of `h`), every one inside the imported and unmoved 5% band, so the pre-registered clause evaluates to **(a)**. Both negative controls survive and are asserted: mis-rotated **121.0800%** (coarse 123.6255%), quadrature-vs-single-drive **384.1297%** (coarse 333.0778%). **(ii) The estimator is still an estimator on the new mesh**, every anchor green: the **same-mesh** best-approximation inequality **12.5225% ≤ 1626.2098%** (a 129.86× separation, both measured in this run — the 10:30 review's sharpening; the coarse 18.7238% is printed, not asserted), `a + b × x` **9.947634e-13** (bound 1e-10) with `x² ê_x` at **2.142147e-01** (floor 1e-4), pinned max \|value\| **0.000e+00** over owned and ghost blocks with **722** free of 22 147 owned CG1 blocks on **66 441** dofs, all six restricted solves reason **2** in 20–25 its. Restricted phantom power **5.499426495e-08 W**, **−1.5681%** from this mesh's primal 5.587038273e-08 W (the coarse mesh read −3.5058%). **(i) Gate (i)'s power accounting is unmoved**: residual **9.795780e-03** (P1) / **9.796465e-03** (P2) inside the unmoved 1e-2 band, conductor-drop control 7.52e-02. **No SAR gate was registered and no band moved** — registering the first coil-driven SAR gate on this estimator at this rung is explicitly the *next review's* ruling, per the pre-registration; `WF-6` stays **🟡** and the five primal asserts in `test_birdcage_sar_map.py` stay red and untouched. **Read with the caveat the module prints:** the sample set is the mesh's own tag-3 centroids, so it grew 51 → **373** points with the refinement; step 1c measured that the sample set is not the mechanism for the `\|B₁⁺\|` column (±2 pp centroid vs ring set), which is why this is read as an `h` rung, but a review wanting the gate may want that control re-run on this column. The rung also turned up a **coil-side** finding, its own entry below: the `\|B₁⁺\|` C4 identities are not mesh-converged either. |
| **Ruled 2026-09-02 18:00 review — the integral construction is the gate; these five asserts become records** | Two constructions printed (a) on the same day. The review picked step 3g's **integral** one (§7 `WF-6` rulings blockquote, last paragraph): it is a statement about the primal field with no estimator, projector or sample set in the path, its partition identity is an exact anchor, and it is cheaper by a whole mesh; 3f's pointwise-restricted 2.5–3.5% is corroboration under `h`. **`WF-6` step 3h** (§9 item 1) registers the twelve integral C4 pairs at the unmoved imported 5% band and re-states this entry's five pointwise asserts as `STEP3_PRIMAL_*` record reproductions at rtol 1e-3 — the pointwise construction retired as a gate with its measurement kept and asserted to *still* exceed the band, not a band widened. This entry **retires with 3h's commit**, which must point its retired-by row here. Until then the five stay red exactly as written. |

### ✅ RETIRED 2026-09-03 (`WF-6` step 3f′, 21:00 implementer slot) — ~~🔴 OPEN 2026-09-02 (`WF-6` step 3f, 15:00 implementer slot) — the gated `|B₁⁺|` **C4 identities are not mesh-converged**: halving the *phantom's* cell size alone moves them **2.19 / 2.11 / 1.89% → 0.62 / 0.60 / 0.56%**, a 1.3–1.6 pp *improvement*, against a pre-registered 0.5 pp no-move anchor~~ — **the sample set is not the mechanism (±1.06 pp ring vs centroid on all eight identities), the anchor is one-sided, and the fall is recorded**

> **✅ RETIRED 2026-09-03, 21:00 implementer slot (`WF-6` step 3f′) — nothing
> was widened to retire this and no `|B₁⁺|` figure was replaced.** The three
> asserts in `test_the_b1_plus_c4_identities_do_not_move_with_the_phantom_h`
> are now **one-sided** (`reading ≤ coarse record + 0.5 pp`: an identity may
> not get *worse* under refinement; a fall is the convergence measurement the
> rung exists to make), the 0.5 pp ceiling itself is unchanged, the 5% band is
> unchanged and still asserted in the same test, and the fine readings are
> recorded **beside** the coarse ones as `STEP3F_B1_PLUS_FINE_RECORDS` (0.6177
> / 0.5966 / 0.5647%, asserted at rtol 1e-3) so the improvement is a number the
> suite pins rather than a slack the one-sidedness hides. The confound this
> entry named — the phantom's `h` and the centroid sample set moving together —
> is **measured and excluded**: step 1c's 96-point rotation-invariant ring set,
> read on this same mesh and the same four solves, agrees with the 373-centroid
> set within **+0.11 / +0.09 / +0.05 pp** on the `|B₁⁺|` identities and
> **+0.24 / −1.05 / −0.71 / −1.01 / −0.41 pp** on the five restricted-CG1 SAR
> identities — all eight inside step 1c's measured ±2 pp bar, asserted, with
> all three negative controls surviving on the ring set (24.1868% / 123.3351% /
> 375.0478%). `20260903T020607Z_WF-6-step3f-prime.log`, **`54 passed` / Status
> 0 / 117 s** at `-n 2` complex with `tests/environment`; the module went `3
> failed, 27 passed` → `0 failed`, 43 items.

| | |
|---|---|
| **Verified at** | `dae3987` + this commit, 0.11 image, complex build, `-n 2`, 2026-09-02 — `20260902T170559Z_WF-6-step3f.log`, `3 failed, 27 passed` / Status 1 / **175 s**, with `tests/environment`. |
| **Literal symptom** | On the `phantom_resolution = 0.0075` mesh (120 499 cells, 2 746 tag-3, `WF-6` step 3f₀) the CG1 `\|B₁⁺\|` C4 covariance identities read **P2@+90° 0.6177%**, **P4@−90° 0.5966%**, **P3@180° 0.5647%** on 373 of 373 phantom centroids, against the coarse mesh's recorded 2.1870 / 2.1146 / 1.8911% — moves of **−1.5693 / −1.5180 / −1.3264 pp**, above the 0.5 pp ceiling the 2026-09-02 weekly review pre-registered for this rung's anchor (i). |
| **What is *not* wrong** | Everything else on the coil side reproduces. Gate (i)'s three-way power accounting closes at **9.795780e-03** (P1) and **9.796465e-03** (P2) inside the unmoved 1e-2 band — the P1 figure agreeing with the coarse mesh's record 9.795751e-03 to 3e-6 relative — with its conductor-drop negative control at 7.52e-02. The mis-rotated `\|B₁⁺\|` control survives at **25.4563%** (coarse 23.2642%), so the finer map has not lost its azimuthal structure. All three identities are far *inside* the 5% band; the band itself is untouched and is asserted in the same test. |
| **Cause (measured, not diagnosed away)** | The `\|B₁⁺\|` C4 identity is **not** at a mesh-converged ≈2% floor, as the step-1d ruling's band provenance assumed. Halving the phantom's `h` alone — nothing else in the fixture changes, the cells outside tag 3 moved 1.9% (step 3f₀) — cuts the identity by **3.5×** at all three angles simultaneously. That is the signature of a *discretisation* floor still falling with `h`, not of a converged physical asymmetry; it is the same mechanism, one estimator in, as the SAR identities' 8.3–9.5% → 2.5–3.4% fall in the same run. Two things move together with the refinement and this run does not separate them: the phantom's `h`, and the centroid sample set that is built from it (51 → 373 points). Step 1c separated them for the *coarse* DG0 column (±2 pp between a centroid set and a rotation-invariant ring set) but nobody has done it here. |
| **Consequence** | Gate (ii) is **not** invalidated — it asserts `≤ 5%` and the new reading is 0.6%, four times *further* inside the band. What is invalidated is the claim that its 2.19 / 2.11 / 1.89% record is a converged floor with "2.3× of headroom": the headroom is at least 8× on a phantom meshed at 7.5 mm. No `\|B₁⁺\|` claim changes (they were always symmetry identities and the §2 language already says so), no absolute or homogeneity claim is created, and `WF-6` stays **🟡**. The three reds are this module's only failures; `test_birdcage_b1_plus_map.py` and the rest of the `WF-6` family are untouched and unmoved on the coarse mesh. |
| **Resolves with** | A **review ruling**, not another implementer slot, and it is the deferred *absolute-convergence rung* arriving early — the 2026-09-02 weekly parked that rung "behind step 3f, whose anchor (i) is its first data point for free", and this is that data point. The options, none of which a slot may take: (1) re-read anchor (i) as a **one-sided** anchor (the identities must not get *worse* by more than 0.5 pp) and record the improvement as the convergence measurement it is — the direction is favourable and the same rung's SAR column depends on it; (2) commission the ring-set control on this mesh (step 1c's construction, already in `test_birdcage_b1_plus_map._ring_points`) so the `h` half and the sample-set half are separated before either number is quoted; (3) run a third rung (`phantom_resolution = 0.00375`, ≈ 8× the tag-3 cells again — cost it against step 3f₀'s +4 414 cells per halving before queueing) to see whether the identities are heading to zero or to a floor. Until a review rules, the three asserts stay exactly as written and red, and **no `\|B₁⁺\|` figure is re-recorded**. |
| **Ruled 2026-09-02 18:00 review — option (2), with (1) folded in** | **`WF-6` step 3f′** (§9 item 2) runs step 1c's 96-point ring set (`test_birdcage_b1_plus_map._ring_points`) on the 0.0075 mesh for both the `\|B₁⁺\|` and the SAR columns (anchor: ring vs centroid within ±2 pp, step 1c's measured separation), re-reads anchor (i) **one-sided** (an identity may not get *worse* than its coarse record by more than 0.5 pp), and records the fine-rung figures 0.6177 / 0.5966 / 0.5647% as `STEP3F_B1_PLUS_FINE_RECORDS` **beside** the coarse ones, never replacing them. Option (3), the 0.00375 rung, is not queued: 3f₀ priced one halving at 5.1× the phantom cells, so the next needs a cost probe (the weekly's). This entry **retires with 3f′'s commit**; the three asserts stay red until then. |
| **Retired by** | `WF-6` step 3f′, 2026-09-03 21:00 implementer slot, this commit — `20260903T020607Z_WF-6-step3f-prime.log:4732–4748` (the ring-vs-centroid table and the one-sided anchor), `54 passed` / Status 0 / **117 s**, `-n 2` complex, standard tier under the 600 s ceiling. All three additions the ruling asked for landed in one window on the same four solves: **(A)** the ring set on both columns, ±2 pp asserted at all eight identities (worst |Δ| **1.0548 pp**, `SAR_P4(−Rx)`); **(B)** the one-sided anchor green at all three angles with the fine records reproduced; **(C)** the twelve C4 **integral** pairs of step 3g/3h's construction on this mesh, printed not gated — **0.1550 / 0.2786 / 0.0220 / 0.1466 \| 0.0276 / 0.2370 / 0.1754 / 0.1177 \| 0.0457 / 0.0621 / 0.1952 / 0.2998%**, worst **0.2998%** against the coarse mesh's 1.5200%, so the pre-registered clause reads **(a)**: *the integral gate's headroom does not shrink with `h`* — 5.1× more headroom on the finer phantom, with the partition identity exact at rtol 1e-10 on all four drives and the P1 total 5.587038273e-08 W agreeing with `mean_sar`'s independent assembly. Every step-3f record reproduced unmoved (120 499 / 2 746 cells, 12.5225% ≤ 1626.2098%, 9.947634e-13, the five SAR identities 3.3600–2.5465%, the primal phantom power). **No band moved, no gate was registered here** (step 3h owns the SAR gate, at fixed `h`), no `\|B₁⁺\|` figure was replaced, and `WF-6` stays 🟡. Option (3) — the 0.00375 rung — is still unqueued and still wants a cost probe. |

### ✅ RETIRED 2026-08-27 (`OPS-27` step 2, 21:00 implementer slot) — a **tenth** stale exact cell-count record, and it is a **fifth** mesh with **no sibling**: `test_dodd_deeds_reactance_wire_resolution.py` recorded 366 207 and 0.11 meshes 365 970

> **✅ RETIRED 2026-08-27, 21:00 implementer slot (`OPS-27` step 2).** The
> literal at `:268` is re-recorded 366_207 → **365_970**, exact equality,
> version-tagged with the 0.7.2 digit and the census log in-comment
> (`GEO-16` precedent, no band). Anchor: the projected half re-run from
> `main` at the same four node ids as the census —
> `20260828T022006Z_OPS-27-step2-wire-projected.log`, **4 passed in
> 459.44s**, Status 0, elapsed 462 s, `-n 2` complex, both rank streams
> identical, against the census's `1 failed, 3 passed`. Same collected
> count (4), so exactly this one name flipped. The module's `growth`
> denominator `138_619` at `:263` is **not** a record (it feeds the printed
> ratio and the unmoved `2.0 < growth < 3.5` band) and was left at its
> 0.7.2 digit — step 1 finding 40's prose class.

> **CANDIDATE OWNER: `OPS-27`** (§9 item 2), extending the leg (f) entry
> below. This one refines that entry's shape rather than repeating it: the
> class is now **ten red names over nine modules carrying five distinct
> meshes**, and this fifth mesh is recorded in **exactly one** place, so
> "one re-record retires up to four names" is an upper bound, not the
> typical case — the sweep must still be by value, but the value→module
> map is ragged (4, 3, 1, 1, 1).

| | |
|---|---|
| **Test (1 new name, footered Status 1, at `4586a13`)** | `tests/validation/test_dodd_deeds_reactance_wire_resolution.py::test_the_refinement_landed_on_the_wire_and_not_on_the_far_field` (1 of the module's 6 validation names; the other 5 are green — 2 in the same run, 2 `pinned` in `20260827T201823Z`, the 6th selected in the same node-id list and green) |
| **Log** | `docs/testing/logs/20260827T202222Z_OPS-26-step2g-dodd-wireres-projected.log` — `1 failed, 3 passed in 434.40s`, **Status 1**, elapsed 436 s, `-n 2`, **complex build**, heavy tier, four names selected by node id. Both rank streams identical (434.36 / 434.40 s). |
| **Symptom** | `AssertionError: the mesh is deterministic and only resolution_wire moved, so the cell count must be the probe's 366207; got 365970`, at `test_dodd_deeds_reactance_wire_resolution.py:268`. |
| **Cause — not diagnosed in code** | Same class: an exact-equality cell count recorded on 0.7.2, never swept when `OPS-18` moved the image. Drift **−237 on 366 207 = −0.0647%**, a *fifth* distinct per-mesh figure and the second negative one (−0.093%, −0.065%, +0.032%, +0.075%, +0.233%) — confirming again that no unmeasured record can be predicted from a measured one. The test's stated premise ("the mesh is deterministic and only `resolution_wire` moved") is not what failed; the *record* is stale. |
| **Sibling sweep (zero compute)** | `grep -rn '366207\|366_207' tests/` returns **only** `:268` and its message at `:270` — no sibling module carries this mesh, unlike 138 619 (three modules) and 417 914 (three). |
| **Disposition** | Filed, **not fixed** — `OPS-26` step 2 is a census and lands no fix. Expected on `main`; **not yours**. |

### ✅ RETIRED 2026-08-27 (`OPS-27` steps 1 and 2) — **five more stale exact cell-count records across five modules, and the class collapses to THREE shared meshes, not nine independent constants**

> **✅ THE 417 914 AND 697 401 FAMILIES RETIRED 2026-08-27, 21:00
> implementer slot (`OPS-27` step 2) — this entry is now fully retired.**
> Two edits closed all three remaining names:
> `test_dodd_deeds_resistance_slab_resolution.NCELLS_FINE` and
> `test_coil_loading_larmor_resolution.NCELLS_FINE`, both 417_914 →
> **418_888**, and
> `test_dodd_deeds_reactance_combined_knobs.NCELLS_COMBINED` 697_401 →
> **697_926** — exact equalities, version-tagged, no band, `git diff --
> src/` empty. Anchors from `main`, both rank-stream-identical:
> `20260828T020157Z_OPS-27-step2-slab.log` — **16 passed in 479.37s**,
> Status 0, elapsed 482 s, `-n 2` complex (census: `1 failed, 15 passed`);
> `20260828T021014Z_OPS-27-step2-knobs.log` — **15 passed in 577.00s**,
> Status 0, elapsed 579 s, `-n 8` complex (census: `1 failed, 14 passed`).
> Collected counts identical to the census runs (16 / 15), so exactly the
> two stale-record names flipped. `larmor_resolution` was **edited but not
> re-run in this slot** — it asserts the same 418 888 the slab run
> measured on the same mesh; its re-run is owed to the next census (see
> the `third_rung` entry below, which shares the constant by import).
> **The owed re-run landed 2026-08-28, 16:30 slot (`OPS-27` step 3):**
> `20260828T213049Z_OPS-27-step3-larmor-resolution.log` — **17 passed in
> 424.32s**, Status 0, elapsed 426 s, `-n 2` complex, both rank streams
> identical, against the census's `1 failed, 16 passed` on the same
> collected 17. Nothing is owed on this entry any more.

> **✅ THE 138 619 FAMILY RETIRED 2026-08-27, 19:30 implementer slot
> (`OPS-27` step 1); the 417 914 and 697 401 families closed by step 2
> above.** Four of this entry's six names are green
> from `main` — `richardson_ladder::test_the_rung_has_its_recorded_cell_count`
> `[10MHz]` and `[30MHz]`, `larmor_probe::test_the_mesh_is_the_mat6_step3_baseline`,
> `transition_30mhz::test_the_mesh_is_the_step1_baseline`. Finding 30's
> "the unit of repair is the mesh, not the file" is **confirmed to its
> strongest form here**: all four names read one constant,
> `test_coil_loading_larmor_probe.NCELLS_BASELINE`, imported by the other
> two modules (and by `test_coil_loading_degree2.py`), so **one edit**
> 138_619 → **138_490** retired four reds. Logs, both `-n 2` complex,
> standard tier, both rank streams identical:
> `20260828T003400Z_OPS-27-step1-richardson.log` — **25 passed in 147.00s**,
> Status 0, elapsed 149 s (census: `2 failed, 23 passed`);
> `20260828T003636Z_OPS-27-step1-probe-30mhz.log` — **23 passed in 149.09s**,
> Status 0, elapsed 150 s (census: `2 failed, 21 passed`). Same collected
> counts as the census runs (25 / 23), so only the two stale-record names
> per run flipped. **Still open in this entry:**
> `slab_resolution::test_the_refinement_landed_in_the_slab_and_not_on_the_wire_or_far_field`
> and `larmor_resolution::test_the_mesh_is_the_mat6_step8_fine_rung`
> (417 914 → 418 888) and
> `combined_knobs::test_the_combined_mesh_is_the_probes_and_both_knobs_moved`
> (697 401 → 697 926).

> **CANDIDATE OWNER: `OPS-27`** (§9 item 2), whose rubric as written names
> **two** sites. The census has now found **nine red names over eight
> modules**, and this entry is the one that changes the shape of the job:
> the nine reds carry only **four distinct** 0.7.2-era cell counts, because
> the same mesh is recorded independently in several modules. Re-recording
> is therefore four measurements and \~nine edits, not nine measurements —
> but a per-file fix that misses a sibling leaves a red behind.

| | |
|---|---|
| **Tests (5 new names, all footered Status 1, all at `cf03754`)** | `test_coil_loading_richardson_ladder.py::test_the_rung_has_its_recorded_cell_count[10MHz]` and `[30MHz]`; `test_coil_loading_larmor_probe.py::test_the_mesh_is_the_mat6_step3_baseline`; `test_coil_loading_transition_30mhz.py::test_the_mesh_is_the_step1_baseline`; `test_dodd_deeds_resistance_slab_resolution.py::test_the_refinement_landed_in_the_slab_and_not_on_the_wire_or_far_field`; `test_dodd_deeds_reactance_combined_knobs.py::test_the_combined_mesh_is_the_probes_and_both_knobs_moved`; `test_coil_loading_larmor_resolution.py::test_the_mesh_is_the_mat6_step8_fine_rung` |
| **Logs** | `20260827T183121Z_OPS-26-step2f-richardson.log` (`2 failed, 23 passed in 141.15s`, `-n 2` complex, 143 s); `20260827T185143Z_OPS-26-step2f-probe-30mhz.log` (`2 failed, 21 passed in 137.35s`, `-n 2` complex, 139 s); `20260827T183401Z_OPS-26-step2f-dodd-slab-resolution.log` (`1 failed, 15 passed in 429.91s`, `-n 2` complex, 431 s); `20260827T184138Z_OPS-26-step2f-dodd-combined-knobs.log` (`1 failed, 14 passed in 568.26s`, `-n 8` complex, 570 s); `20260827T185422Z_OPS-26-step2f-larmor-resolution.log` (`1 failed, 16 passed in 428.37s`, `-n 2` complex, 430 s). Rank streams identical on all five. |
| **Symptom** | Each is `AssertionError: … the count must be … <recorded>; got <measured>`, with only three distinct pairs across seven of the nine census names: **138 619 → 138 490** (−129, **−0.093%**) in `richardson_ladder` ×2, `larmor_probe`, `transition_30mhz`; **417 914 → 418 888** (+974, **+0.233%**) in `slab_resolution`, `larmor_resolution` and (leg (e)) `third_rung`; **697 401 → 697 926** (+525, **+0.075%**) in `combined_knobs`. The fourth quantity is leg (d)'s **2 807 309 → 2 808 204** (+0.032%) in `mesh_cache`. |
| **Cause — not diagnosed in code** | 0.11's gmsh meshes the *same* geometry to a slightly different cell count than 0.7.2 did; the records are exact equalities never swept when `OPS-18` re-recorded. **Drift is not one signed constant** — it is −0.093%, +0.032%, +0.075%, +0.233% on four meshes, so it is per-mesh, not a global offset, and cannot be predicted for an unmeasured record. Every affected module's *other* names are green in the same run: the physics reproduces, the cell count does not. |
| **What this adds to `OPS-27`'s rubric** | (1) The sweep must be **by cell-count value across modules**, not per file — 138 619 alone appears in three modules, 417 914 in three. (2) `grep -rn '0\.7\.2' tests/` is insufficient (leg (e)'s note, now confirmed twice): none of these five was reachable by version tag, only by reading a red's assertion message. (3) One re-record of a shared mesh retires up to four red names at once, so the job is **four measurements**, not nine. |
| **Disposition** | Filed, **not fixed** — `OPS-26` step 2 is a census and lands no fix. Expected on `main`; **not yours**. |

### ✅ RETIRED 2026-08-28 (`OPS-27` step 3, 16:30 implementer slot) — `test_coil_loading_larmor_third_rung.py` asserted an **exact** fine-rung cell count recorded on the 0.7.2 image (417 914) and 0.11 meshes 418 888: the **fourth** site of the stale-record class, and the largest drift yet at 0.233%

> **✅ RETIRED 2026-08-28, 16:30 implementer slot (`OPS-27` step 3) — the
> owed re-run executed and the module is green.**
> `20260828T213807Z_OPS-27-step3-thirdrung.log` — **18 passed in 291.03s**,
> Status 0, elapsed 293 s, `-n 8` complex, `TH11_STEP5_RUNG=fine`, heavy
> tier; all eight rank streams identical (290.98–291.05 s). The census read
> `1 failed, 17 passed` on the same collected 18, so exactly the one
> stale-record name flipped and no other name's status moved. No code
> changed in this slot — the import of `NCELLS_FINE` from
> `larmor_resolution` (finding 41) is what re-recorded it in step 2; this
> entry stayed open only for the missing execution. **Cost note (a
> correction to finding 25's inference):** the module's *cold* price was
> never measured — the 900 s window ruled for it returned at 291 s, i.e.
> the 304 s warm figure was essentially the whole price, and this slot's
> `larmor_resolution` window ran first as ruled, so the figure is again
> warm-fixture. The ≥ 500 s cold estimate is **not** confirmed; nothing on
> record measures this module cold.

> **🟡 RE-RECORDED but NOT RE-RUN in the 21:00 slot.** This module holds no
> constant of its own: `:443` asserts the `expected` drawn from the `fine`
> entry of its rung table, which is `NCELLS_FINE` **imported from**
> `tests/validation/test_coil_loading_larmor_resolution` (step 1
> finding 38's import-graph rule, repeating). So the single 417_914 →
> **418_888** edit in `larmor_resolution` re-records this name too, and
> `git diff -- tests/validation/test_coil_loading_larmor_third_rung.py` is
> empty. The same mesh value was **measured green** in the slab run
> (`20260828T020157Z_OPS-27-step2-slab.log`, 16 passed), but this module
> was not itself re-executed: it is a warm-cache-only ≈304 s at `-n 8`
> (`OPS-26` finding 25) and did not fit beside the slot's three windows.
> **Re-run owed to the next census.** Entry stays filed until then.

> **CANDIDATE OWNER: `OPS-27`** (§9 item 2) — same class, same remedy
> (re-record exact, version-tagged, `GEO-16` precedent, no band). `OPS-27`'s
> §7 rubric as written names **two** sites; this is a third and the review
> should decide whether to fold it in or leave it for the sweep clause.
> Filed here so it is not lost either way.

| | |
|---|---|
| **Test** | `tests/validation/test_coil_loading_larmor_third_rung.py::test_the_rung_is_inside_the_priced_ceiling` (1 of the module's 7 collected names; the other 6 are green in the same run) |
| **Log** | `docs/testing/logs/20260827T171110Z_OPS-26-step2e-thirdrung-destubbed.log` — `1 failed, 17 passed in 302.31s`, **Status 1**, elapsed 304 s, `-n 8`, **complex build**, `TH11_STEP5_RUNG=fine`, heavy tier. All eight rank streams identical (302.29–302.35 s). Footered, so it counts as a red under the census's fail-closed control. Verified at `b9b8cc7`. |
| **Symptom** | `AssertionError: frequency does not reach the mesh generator, so the count must be the rung's recorded 417914; got 418888`, at `test_coil_loading_larmor_third_rung.py:443`. |
| **Cause — not diagnosed in code** | An **exact-equality** cell-count record made on the **0.7.2** image, compared against a mesh generated by 0.11's gmsh. Drift is **+974 cells on 417 914 = +0.233%** — an order of magnitude larger in relative terms than the `mesh_cache` site's 0.032%, but still far below any geometry change, and the same sign. The test's own premise (*"frequency does not reach the mesh generator"*) is **not** what failed: the count is frequency-independent as claimed, it is the *record* that is stale. Corroborated by the other 6 names in the module — including both complex-power identities, the free-solve dissipation identity and `test_the_fine_rung_reproduces_step2s_recorded_deviation` — passing in the same run, i.e. the physics on this rung reproduces its `TH-11` step-2 record while only the cell count does not. |
| **Class — fourth instance this census** | Same class as leg (d)'s finding 23 (`mesh_cache`, 0.032%), leg (c)'s finding 19 (`test_geometry_floor_discriminator.py`) and `GEO-16`'s two-torus red: **exact-equality records made on 0.7.2, never swept after the 0.11 bump**. Four sites, three distinct quantities (residual, cell count ×2, separation ratio). Note for `OPS-27`: its planned `grep -rn '0\.7\.2' tests/` sweep would **not** necessarily have found this one — the constant at `:443` must be reached by *reading the assertion messages of the census reds*, not only by the version tag, so the sweep clause should be widened to exact-equality mesh counts regardless of tag. |
| **Disposition** | Filed, **not fixed** — `OPS-26` step 2 is a census and lands no fix. Expected on `main`; **not yours**. |

### ✅ RETIRED 2026-08-27 (`OPS-27` step 1, 19:30 implementer slot) — `test_coil_loading_larmor_mesh_cache.py` asserted an **exact** third-rung cell count recorded on the 0.7.2 image (2 807 309) and 0.11 meshes 2 808 204: a 0.032% gmsh drift against an equality record

> **✅ RETIRED 2026-08-27, 19:30 implementer slot.** `NCELLS_THIRD`
> 2_807_309 → **2_808_204**, version-tagged with the 0.7.2 digit and
> `20260827T141059Z_OPS-26-step2d-meshcache-real.log` in-comment, docstring
> copies moved. The disposition line below argues an equality on a mesher
> count is "arguably the wrong shape"; **that call was not taken** — the
> review ruled exact-and-version-tagged (`GEO-16` precedent, no band), and
> `NCELLS_THIRD_CEILING` is unchanged. Green from `main`:
> `20260828T003915Z_OPS-27-step1-meshcache-real.log` — **12 passed, 4
> skipped in 254.75s**, Status 0, elapsed 256 s, `-n 2` **real** build,
> standard tier, both rank streams identical; against the census's `1
> failed, 11 passed, 4 skipped`, i.e. the same 4 `tests/environment`
> complex-only skips and the one red flipped.

> **OWNER ASSIGNED 2026-08-27, 10:30 review: `OPS-27`** (§9 item 2) — the
> stale-record class, one chunk for both open sites; re-record exact and
> version-tagged on the `GEO-16` precedent, no band. Retires with it.

| | |
|---|---|
| **Test** | `tests/validation/test_coil_loading_larmor_mesh_cache.py::test_the_cached_rung_is_the_priced_mesh` (1 of the module's 5 collected names; the other 4 are green in the same run) |
| **Log** | `docs/testing/logs/20260827T141059Z_OPS-26-step2d-meshcache-real.log` — `1 failed, 11 passed, 4 skipped in 217.70s`, **Status 1**, elapsed 219 s, `-n 2`, **real build**, standard tier. Both rank streams identical (217.68 / 217.70 s). Footered, so it counts as a red under the census's fail-closed control. Verified at `5590b81`. |
| **Symptom** | `AssertionError: the third rung meshed to 2808204 cells, not the probe's recorded 2807309: the fixture changed rather than being re-meshed`, `assert 2808204 == 2807309`. The 4 skips are all in `tests/environment` (complex-only tests correctly skipped in the real build), **not** in the census roots. |
| **Cause — not diagnosed in code** | An **exact-equality** cell-count record made by the `TH-11` step-5 probe on the **0.7.2** image, compared against a mesh generated by 0.11's gmsh. The drift is **+895 cells on 2.81 M = +0.032%** — far too small to be a geometry change and exactly the size of a mesher-version tetrahedralisation difference. No `src/` or fixture edit is implicated: the assertion message's own alternative ("the fixture changed rather than being re-meshed") is the *other* branch, and nothing in `git log` touches this fixture since `OPS-18`. |
| **Class — third instance this census** | This is the same defect class as leg (c)'s finding 19 (`test_geometry_floor_discriminator.py`, a pre-`OPS-18` 128 MHz constant) and `GEO-16`'s two-torus red: **records not swept after the 0.11 re-record**. Three sites now. The census's value here is that the class, not the individual constant, is the finding — a review should commission one sweep chunk over all exact-equality records made on 0.7.2 rather than three one-constant fixes. |
| **Disposition** | Filed, **not fixed** — `OPS-26` step 2 is a census and lands no fix (`OPS-26` §7 scope). The fix is a one-constant re-record on the 0.11 image with the basis stated in-comment, and an equality assertion on a mesher cell count is arguably the wrong shape in the first place (a ±0.1% band would survive an image bump); both calls belong to a chunk, not to this slot. Until then this red is expected on `main` and is **not yours**. |

### ✅ RETIRED 2026-08-27 (`OPS-27` step 1, 19:30 implementer slot) — `test_geometry_floor_discriminator.py` asserted the **pre-`OPS-18`** 128 MHz record (1.8260%) and measured `OPS-18`'s re-recorded 1.7686%: a stale constant, not a physics regression

> **✅ RETIRED 2026-08-27, 19:30 implementer slot — landed exactly as ruled.**
> `RECORD_128_RELL2` 0.01826 → **0.017686** and `RECORD_128_SEPARATION`
> 57.31 → **59.16**, version-tagged with the 0.7.2 digits and both census
> logs in-comment (`GEO-16` precedent), every docstring copy moved and the
> dated 2026-08-13 result block annotated rather than rewritten. **No band
> touched** (`REPRODUCTION_BAND` 1%, `CELL_COUNT_BAND` 1% unchanged) and the
> 64 MHz leg untouched. Green from `main`:
> `20260828T003300Z_OPS-27-step1-geomfloor.log` — **12 passed in 46.45s**
> (11 `tests/environment` + this module's one test), Status 0, elapsed 49 s,
> `-n 2` complex, smoke tier, both rank streams identical.

> **OWNER ASSIGNED 2026-08-27, 10:30 review: `OPS-27`** (§9 item 2) —
> `RECORD_128_RELL2` / `RECORD_128_SEPARATION` → `OPS-18`'s 0.017686 /
> 59.16, every docstring copy moved in the same commit. Retires with it.

| | |
|---|---|
| **Test** | `tests/validation/test_geometry_floor_discriminator.py::test_larmor_sphere_residual_at_the_priced_fine_mesh` (the module's only test) |
| **Log** | `docs/testing/logs/20260827T125507Z_OPS-26-step2c-v42-geomfloor.log` — `1 failed in 22.15s`, **Status 1**, elapsed 23 s, `-n 2`, **complex build** (`FEM_EM_REQUIRE_COMPLEX=1`), smoke tier. Both rank streams identical. Footered, so it counts as a red under the census's fail-closed control. Verified at `58c77d9`. |
| **Symptom** | `AssertionError: the 128 MHz relL2 moved to 1.7686% from the recorded 1.8260% (3.14% > 1%) at the mesh it was recorded on. That is a regression in the fixture or this file, not a geometry finding — the 64 MHz reading above must not be interpreted until it is explained`, `assert 0.03141140883816601 < 0.01`. The 64 MHz leg passes; only the 128 MHz record fires. |
| **Cause — not diagnosed in code, but the number is already explained on `main`** | The measured 1.7686% is **`OPS-18`'s value**. PROJECT_PLAN §2 and CLAUDE.md both record that `TH-10`'s 128 MHz figure "is 1.769% on the 0.11 image `main` boots — re-recorded with its mesh by `OPS-18`, 2026-08-22", against the 1.826% originally recorded on the 0.7.2 image at `TH-10` closure (2026-08-13). This file's constant is therefore the **pre-`OPS-18`** one: the assertion is comparing a 0.11 measurement against a 0.7.2 record, and its own message ("a regression in the fixture or this file") names the right disposition — it is *this file*. The test is doing its job; nobody updated it when `OPS-18` re-recorded. |
| **What this does *not* mean** | It is **not** evidence against `TH-10`, and not a new 0.11 break. The 3.14% is the distance between two *recorded* numbers, not a drift in the solve; the solve reproduces `OPS-18` to the digits printed. Do not re-open `TH-10` on this. |
| **Disposition** | Filed, not fixed — `OPS-26` step 2 is a census and lands no fix (`OPS-26` §7, leg (c) scope). The fix is a one-constant re-record with the `OPS-18` mesh cross-referenced in-comment, and it belongs to a chunk that can re-run the priced fine mesh and state the basis; a review should commission it. Until then this red is expected on `main` and is **not yours**. |
| **Census accounting** | Counted as the census's 4th red repo-wide and as leg (b)'s owed complex conversion: the name carried by leg (b) as `deferred — complex-only, SKIPPED in the real build` is now **observed**, and it resolved to a red rather than the green the conversion pattern produced six times in leg (d). |

### ✅ RETIRED 2026-08-28 (`OPS-28`, 22:30 implementer slot) — the whole of `tests/ports/test_port_orientation_sensitivity.py` dies on `'_DummyComm' object has no attribute 'allgather'`: an **`OPS-14` rank-safety reduction broke its test double**, and nothing scheduled has run the module since

> **RETIRED 2026-08-28 by `OPS-28`.** `_DummyComm` gained
> `allgather(value) -> [value]` beside its `allreduce`
> (`tests/ports/test_port_orientation_sensitivity.py:23-31`); `src/` is
> untouched, the reduction stays, the deprecated route stays runnable.
> Bracketed by measurement: `20260828T033037Z_OPS-28-red-baseline.log`
> reproduces the `AttributeError` on both names (`3 failed, 14 passed in
> 1.50s`, Status 1, 3 s) and `20260828T033055Z_OPS-28-gate.log` is
> `2 failed, 15 passed in 0.79s` (Status 1, 2 s) on the identical command.
> **Disposition of the two names.**
> `…::test_port_orientation_flip_changes_induced_voltage_sign` is **green**
> — it asserts the sign-flip identity on the placeholder route
> (`V(P2) = +5.000000e-02 V` aligned, `−5.000000e-02 V` flipped, equal in
> magnitude to `rel=1e-12`, coupling factor `+1.0e-01 → −1.0e-01`), and
> that name is now retired outright.
> `…::test_port_orientation_flip_changes_off_diagonal_sparameter_sign`
> **now reaches its S-matrix assertion and is red there** — it belongs to
> **entry 3** and only to entry 3, which is re-dated below with the
> correction this run measured (the vanishing wave on *this* test is the
> **off-diagonal**, not the diagonal). The `AttributeError` symptom this
> entry filed is gone from both names, so the entry retires whole.

| | |
|---|---|
| **Tests** | `tests/ports/test_port_orientation_sensitivity.py::test_port_orientation_flip_changes_induced_voltage_sign` (**new** — not previously filed anywhere)<br>`tests/ports/test_port_orientation_sensitivity.py::test_port_orientation_flip_changes_off_diagonal_sparameter_sign` (already filed under **known-issues 3**, but for a *different* symptom — see "Relationship to entry 3") |
| **Log** | `docs/testing/logs/20260827T093747Z_OPS-26-step2b-p03-orient.log` — `2 failed, 3 warnings in 0.69s`, **Status 1**, elapsed 2 s, `-n 2`, real build, smoke tier. Footered, so these count as reds under the census's fail-closed control. |
| **Symptom** | Both tests raise, identically and immediately: `src/fem_em_solver/ports/excitation.py:265: in run_placeholder_port_coupling_case` → `{int(v) for values in problem.mesh.comm.allgather(` → `AttributeError: '_DummyComm' object has no attribute 'allgather'`. The second test reaches it one frame deeper, through `sparameters.py:368: in run_n_port_sparameter_sweep`. |
| **Cause — diagnosed, one line, not a 0.11 break** | `OPS-14` hardened `run_placeholder_port_coupling_case` to reduce the rank-local `cell_tags.values` before validating port tags (`excitation.py:262-268`, whose own comment says it is "fixed here so the deprecated route stays *runnable* as `PORT-1` step 4's negative control"). The module's stub comm — `tests/ports/test_port_orientation_sensitivity.py:16-21` — defines `rank` and a `staticmethod allreduce` and **nothing else**, so the new `allgather` call finds no attribute. The reduction is correct and must not be reverted; the double is what is stale. Note the irony recorded verbatim: the change made to keep the deprecated route runnable is what stopped it running, because the module that exercises it is not in any scheduled command. |
| **Relationship to entry 3** | Entry 3 ("Port tests assert a non-zero S-matrix diagonal on a matched port") lists `…::test_port_orientation_flip_changes_off_diagonal_sparameter_sign` as red for the *zero-diagonal* reason. That is no longer the symptom: the test now dies earlier, in the tag reduction, and never reaches its S-matrix assertion. Entry 3's diagnosis is unfalsified but currently **unreachable** on that test. `…::test_port_orientation_flip_changes_induced_voltage_sign` was never in entry 3 at all — it does not touch an S-matrix diagonal — so it is a genuinely new red. Both entries stay open; whichever chunk fixes the double must re-read entry 3 afterwards to see what the module actually asserts. |
| **Why the census counts this as a hit** | This is precisely the class `OPS-26` was commissioned to find (§7 `OPS-26`, "a module that no scheduled command runs can be red or dead indefinitely") — but it is a **new sub-class**: not an un-migrated 0.11 call site (step 1 swept those and `src/`+`tests/` are clean at 434 sites), and not a gmsh regression. It is **test-double drift behind a rank-safety fix**, invisible to `check_dolfinx_api_migration.py` by construction, since `comm.allgather` is a valid mpi4py API and `_DummyComm` is not a DolfinX type. Step 1's static sweep structurally cannot see this; only execution can. |
| **Verified at** | `b39799e`, real build, `-n 2`. |
| **Fix** | **Deliberately not fixed** — `OPS-26` is a census and files rather than repairs, per the item's own rule. The repair is one line (give `_DummyComm` an `allgather` returning `[value]`), but it belongs with a `PORT-0`/`PORT-1` owner who can also dispose of entry 3's assertions on the same module, and it needs a decision this census may not take: whether the deprecated placeholder route is kept runnable at all. |

### ✅ RETIRED 2026-08-28 (`GEO-23` step 2b, 13:30 implementer slot) — ~~🔴 OPEN 2026-08-27, re-headed 2026-08-28 (`GEO-23` step 1) — a **fifth** site of "Invalid boundary mesh (overlapping facets)"~~ — **`tests/materials/test_phantom_material_model.py:110` moved from `resolution=0.03` (does not mesh) to step 1's coarsest meshing rung `0.024`; the module is `4 passed` at `-n 1` and `-n 2` complex, `5464` cells to the digit**

> **GEOMETRY HALF CLOSED 2026-08-28 (`GEO-23` step 2b, 13:30 slot) — the
> sizing lever.** The call site is the only thing that moved: `0.03 → 0.024`,
> the coarsest rung step 1's monotone `-n 1` ladder measured as meshing on
> `coil_phantom_domain` (5 464 cells; the finer rungs are 0.0192 / 0.01536 /
> 0.012288). Nothing in `src/` and no band, tolerance or record was touched —
> the module pins no cell count, so the assertions that had to survive the
> re-mesh are its physics ones, and they all pass. Complex build,
> `FEM_EM_REQUIRE_COMPLEX=1`: `4 passed in 2.66s` at `-n 1`
> (`20260828T183204Z_GEO-23-step2b-phantommaterial-n1.log`, Status 0, 4 s) and
> `4 passed in 1.63s` at `-n 2`
> (`20260828T183214Z_GEO-23-step2b-phantommaterial-n2.log`, Status 0, 3 s),
> against the census red's `1 failed, 3 passed`. The printed global cell count
> is **5464 at both widths**, reproducing the step-1 ladder's 5 464 exactly
> (0.00% against the ±1% band). Both halves of this entry are now closed.

> **DEADLOCK HALF CLOSED 2026-08-28 (`GEO-23` step 2a, 12:00 slot).** The
> `-n 2` command no longer hangs: `coil_phantom_domain`'s rank-0 gmsh build is
> now wrapped, so rank 1 raises
> `RuntimeError: coil_phantom_domain geometry generation failed on rank 0
> (resolution=0.03); this is rank 1` instead of blocking in `_model_to_mesh`.
> **Status 1 in 2 s** where step 1 recorded **Status 124 at 120 s**, summary
> unchanged at `1 failed, 3 passed`
> (`20260828T170347Z_GEO-23-step2a-phantommaterial-n2.log`; `-n 1` control
> unchanged, `20260828T170340Z_…-n1.log`). **The entry stays OPEN for the
> geometry red itself** — the sizing lever (`GEO-23` step 2b) is what retires
> it; only the cost of observing it changed.
>
> **CORRECTED 2026-08-28 by measurement (`GEO-23` step 1, 09:00 slot).** This
> entry's "PASSED on one rank, FAILED on the other" claim does not survive an
> `-n 1` run. At `-n 1` the test is **Status 1 in 2 s** with the same
> `overlapping facets` string (`20260828T140613Z_GEO-23-step1b-phantommaterial-n1.log`,
> `1 failed, 3 passed in 1.18s`), so the failure is a deterministic property of
> the geometry. Re-reading the `-n 2` log with that in hand
> (`20260828T140622Z_…-n2.log`, lines 48–56): the second rank prints its three
> passes and then **nothing** — it never reaches the fourth name. **Absence of
> a verdict is not a PASS.** The `-n 2` Status 124 is the surviving rank
> blocking in the collective after the other raises — a **teardown/raise-path**
> effect, not a partition-dependent trigger. Consequently the row below headed
> "What it adds to `GEO-23`" is withdrawn: there is **no** second
> partition-dependent site, and the resolution-floor reading is *supported*,
> not weakened — `GEO-23` step 1's ladder shows this generator failing only at
> its own `resolution=0.03` and meshing at 0.024 / 0.0192 / 0.01536 / 0.012288.
> The entry stays OPEN (it retires with `GEO-23` step 2); only its *diagnosis*
> is corrected.

| | |
|---|---|
| **Test** | `tests/materials/test_phantom_material_model.py::test_phantom_material_assignment_and_time_harmonic_pipeline_wiring` |
| **Log** | `docs/testing/logs/20260827T093043Z_OPS-26-step2b-materials-complex.log` — complex build, `FEM_EM_REQUIRE_COMPLEX=1`, `-n 2`, `timeout -k 30 180`. **Status 124, elapsed 181 s.** The pytest summary printed (`1 failed, 14 passed, 30 warnings in 20.79s`) and then teardown never returned — the same eat-the-window mechanism as leg (a)'s findings 2 and 11. |
| **Symptom** | Interleaved rank streams show the *same* test resolving differently on the two ranks: `PASSED [ 66%]` on one, `FAILED [100%]` on the other, with `Exception: Invalid boundary mesh (overlapping facets) on surface 1 surface 1`, then `MPI_Abort` on the kill. |
| **Why it was run** | Leg (a) left this name as the one remaining `tests/materials` runtime skip — complex-only, skipped in the real build — and §9 item 1 owed a ≈30 s complex command to convert it to green. **The conversion failed**: it converts to a rank-divergent abort instead, so the name stays **`deferred — rank-divergent gmsh abort, no Status-0/1 footer`** in the leg (a) table rather than becoming green. Leg (a)'s totals are unchanged at 184 / 189 (182 green, 2 red, 5 deferred); only this deferral's *reason* is upgraded from "runtime skip" to the above. |
| **Not counted as a red** | Same ruling as leg (a) finding 11, applied for consistency: the census's fail-closed control admits reds **only** from Status-0/1 footered runs, and this run has neither. Filed as a mechanism finding. |
| **What it adds to `GEO-23`** | `GEO-23` was commissioned (03:00 review) to own the "overlapping facets" family across four generators. This is a **fifth** call site and — with `test_boundary_condition_selection.py` — the **second** where the trigger is demonstrably partition-dependent rather than a property of the geometry. Two independent rank-dependent sites materially strengthens leg (a)'s reading that the resolution-floor explanation the three earlier entries share is incomplete, and `GEO-23` step 1's 2 × 3 Status-by-rank-width table should include this module. |
| **Verified at** | `b39799e`, complex build, `-n 2`. |
| **Fix** | Not fixed, not diagnosed further — `OPS-26` files, `GEO-23` owns. |

### ✅ RETIRED 2026-08-28 (`GEO-23` step 2b, 13:30 implementer slot) — ~~🔴 OPEN 2026-08-27 — `test_boundary_condition_selection.py` **deadlocks the whole command** at `-n 2`~~ — **`tests/solver/test_boundary_condition_selection.py:26` moved from `resolution=0.04` (does not mesh) to step 1's coarsest meshing rung `0.032`; the module is `3 passed, 1 skipped` at `-n 1` and `-n 2`, `1213` cells to the digit**

> **GEOMETRY HALF CLOSED 2026-08-28 (`GEO-23` step 2b, 13:30 slot) — the
> sizing lever.** `_make_problem`'s `cylindrical_domain` call moved
> `0.04 → 0.032`, the coarsest rung step 1's monotone `-n 1` ladder measured
> as meshing (1 213 cells; the finer rungs are 0.0256 / 0.02048 / 0.016384).
> Nothing in `src/`, no band and no record: the module pins no cell count, so
> the assertions that had to survive the re-mesh are the boundary-condition
> physics ones, and they all pass. Real build: `3 passed, 1 skipped in 0.94s`
> at `-n 1` (`20260828T183106Z_GEO-23-step2b-bcsel-n1.log`, Status 0, 2 s) and
> `3 passed, 1 skipped in 0.80s` at `-n 2`
> (`20260828T183116Z_GEO-23-step2b-bcsel-n2.log`, Status 0, 2 s), against the
> census red's `1 failed, 2 passed, 1 skipped`; the one skip is the
> `complex_only` name, unchanged. Printed global cell count **1213 at both
> widths**, reproducing the step-1 ladder exactly (0.00% against ±1%). Both
> halves of this entry are now closed.

> **DEADLOCK HALF CLOSED 2026-08-28 (`GEO-23` step 2a, 12:00 slot) — the
> headline "deadlocks the whole command at `-n 2`" no longer holds.**
> `cylindrical_domain`'s rank-0 gmsh build is now wrapped, so rank 1 raises
> `RuntimeError: cylindrical_domain geometry generation failed on rank 0
> (resolution=0.04); this is rank 1` instead of blocking in `_model_to_mesh`.
> **Status 1 in 2 s** where step 1 recorded **Status 124 at 120 s**, summary
> unchanged at `1 failed, 2 passed, 1 skipped`
> (`20260828T170311Z_GEO-23-step2a-bcsel-n2.log`; `-n 1` control unchanged,
> `20260828T170303Z_…-n1.log`). **The entry stays OPEN for the geometry red
> itself** — `GEO-23` step 2b's sizing move is what retires it.
>
> **CORRECTED 2026-08-28 by measurement (`GEO-23` step 1, 09:00 slot) — this
> entry's own requested `-n 1` command has now been spent, and it resolves the
> pre-stated reading against rank-dependence.** At `-n 1` the test is **Status
> 1 in 3 s** with the same `overlapping facets on surface 1 surface 1`
> (`20260828T140041Z_GEO-23-step1a-bcsel-n1.log`,
> `1 failed, 2 passed, 1 skipped in 1.85s`). The failure is therefore a
> deterministic property of the geometry, and `GEO-23` step 1's ladder confirms
> it: `cylindrical_domain` fails **only** at this fixture's own
> `resolution=0.04` and meshes at 0.032 (1 213 cells) / 0.0256 / 0.02048 /
> 0.016384.
>
> **The "one rank raises while the other returns" observation is a
> log-interleave artifact.** In the `-n 2` log
> (`20260828T140055Z_GEO-23-step1a-bcsel-n2.log`, lines 48–54) the apparent
> `PASSED [ 25%]` on the failing name is the *other* rank's verdict for its own
> **first** test, appended mid-line to this rank's name line — the percentages
> settle it (`[ 25%]` cannot be the third of four tests), and the failing
> name's only verdict anywhere in the log is `FAILED [ 75%]`. This is exactly
> the interleave trap `GEO-23`'s §7 entry warns about. What *is* rank-dependent
> is only the **teardown**: one rank raises, the survivor blocks in the
> collective, Status 124 at 121 s. `birdcage_port_domain` is the control —
> it re-raises the rank-0 throw as a `RuntimeError` on every rank and footers
> at Status 1 in 5 s, so wrapping the raise is a step-2 lever that would retire
> the deadlock without touching any mesh. Entry stays OPEN (retires with
> `GEO-23` step 2); only its diagnosis is corrected.
>
> **OWNER ASSIGNED 2026-08-27, 03:00 review: `GEO-23`** step 1 (a) — the
> `-n 1` command named at the foot of this entry is its first move.
>
> **Not a census red.** `OPS-26`'s fail-closed control admits a red only from
> a run with a Status-0-or-1 footer of its own. Both runs below ended
> **Status 124**, so the four tests in this module are counted
> `deferred — module-scoped command deadlocked; no footer`, and the failures
> here are recorded as a **mechanism finding**, not as a counted red. The next
> chunk that owns this must re-establish it in a footered run before treating
> any name below as failing.
>
> **Where this fired.** `tests/solver/test_boundary_condition_selection.py`,
> `mpiexec -n 2`, both builds, `OPS-26` step 2 leg (d) (2026-08-27 00:00
> implementer slot), each as a command containing **only this module**:
>
> | log | build | timeout | Status / elapsed |
> |---|---|---|---|
> | `20260827T050123Z_OPS-26-step2a-legd-m01-bcsel.log` | real | `-k 30 240` | **124 / 241 s** |
> | `20260827T051201Z_OPS-26-step2a-legd-m01-bcsel-complex.log` | complex, `FEM_EM_REQUIRE_COMPLEX=1` | `-k 30 180` | **124 / 180 s** |
>
> **The finding — the same test PASSES on one rank and FAILS on the other, in
> the same run.** The complex log is the clean read of it (the real log's two
> streams interleave mid-line and its per-test attribution cannot be trusted):
>
> ```
> ...::test_time_harmonic_solver_boundary_natural_selects_empty_dirichlet_set PASSED [ 46%]
> ...::test_time_harmonic_solver_boundary_natural_selects_empty_dirichlet_set FAILED [ 93%]
> ```
>
> (`20260827T051201Z_..._bcsel-complex.log:69,76`.) The failing rank raises
>
> ```
> Exception: Invalid boundary mesh (overlapping facets) on surface 1 surface 1
> ```
>
> **That is the mechanism of the hang.** One rank raises inside gmsh and
> leaves the collective; the other completes and blocks in the next
> collective; nothing progresses until `timeout -k 30` sends KILL, and the
> trailer is a PETSc `Caught signal number 15` followed by
> `Abort(59) ... MPI_Abort(MPI_COMM_WORLD, 59)`. This is the `mag:1`-class
> divergence the leg-(c) slot inferred for `tests/solver` as a root
> (`PROJECT_PLAN.md` §7 `OPS-26`, finding 8) — leg (d) localizes it to **this
> one module**: the other **twelve** `tests/solver` modules each returned a
> Status-0 footer when run as their own command.
>
> **Two consequences worth keeping.**
> 1. **"Overlapping facets" is not deterministic.** Every prior entry for that
>    string (`GEO-21` open birdcage; the coil+phantom generator; and
>    `birdcage_port_domain`, the entry below) reads as a property of a
>    geometry. Here the *same* geometry, in the *same* run, meshes on one rank
>    and fails on the other — so at least on this call path the trigger is
>    rank-partition-dependent, which a resolution-floor reading does not
>    explain. Stated as an observation from two runs, not a diagnosis.
> 2. **The second test's failure is the leg-(c) candidate signature.**
>    `test_time_harmonic_solver_boundary_pec_is_applied_to_solve_path` fails
>    with `IndexError: index 0 is out of bounds for axis 0 with size 0` —
>    the identical string the leg-(c) slot recorded on 21 of 23 discarded
>    names. It appears here on a module run in isolation on a swept cache,
>    which makes "one shared cause cascading" more plausible than it was; it
>    is still not established, because this run has no footer either.
>
> **Cache state is exonerated.** The slot opened with
> `find /root/.cache/fenics -name '*.c' -size 0 -print -delete`, which printed
> **nothing** (`20260827T050052Z_OPS-26-step2a-legd-collect.log:34`), and the
> real-build run above was the *first* command after it. The poisoned-stub
> entry lower in this file does not apply.
>
> **Next step for the owning chunk (a `solver`/`mesh` chunk, not the census):**
> run this module at `-n 1`. If it returns a Status-1 footer, the failure is
> real and the deadlock is purely the rank asymmetry; if it goes green, the
> failure itself is partition-dependent. One smoke-tier command settles it.
> `OPS-26` deliberately did not spend it — `-n 1` is not the census's recorded
> width and an observation at it would not count.

### ✅ RETIRED 2026-08-28 (`GEO-23` step 2b, 13:30 implementer slot) — ~~🔴 OPEN 2026-08-27 — `test_phantom_field_metrics_and_exports_are_finite` aborts in gmsh with "Invalid boundary mesh (overlapping facets)" on the **coil+phantom** geometry~~ — **`tests/post/test_phantom_field_metrics.py:35` moved from `resolution=0.03` (does not mesh) to step 1's coarsest meshing rung `0.024`; the module is `2 passed` at `-n 1` and `-n 2` complex, `5464` cells to the digit**

> **GEOMETRY HALF CLOSED 2026-08-28 (`GEO-23` step 2b, 13:30 slot) — the
> sizing lever.** Same generator and same rung as the
> `test_phantom_material_model.py` entry above (step 1 found the two modules
> call `coil_phantom_domain` with byte-identical kwargs, so the unit of repair
> is the generator call): `0.03 → 0.024`, 5 464 cells. Nothing in `src/`, no
> band, no record; the module pins no cell count, so the assertions that had
> to survive the re-mesh are its phantom |E|/|B| metric and export ones, and
> they all pass. Complex build, `FEM_EM_REQUIRE_COMPLEX=1`:
> `2 passed in 1.71s` at `-n 1`
> (`20260828T183223Z_GEO-23-step2b-phantommetrics-n1.log`, Status 0, 3 s) and
> `2 passed in 1.67s` at `-n 2`
> (`20260828T183231Z_GEO-23-step2b-phantommetrics-n2.log`, Status 0, 3 s),
> against the census red's `1 failed, 1 passed`. Printed global cell count
> **5464 at both widths**, reproducing the step-1 ladder exactly (0.00%
> against ±1%). Both halves of this entry are now closed.

> **DEADLOCK HALF CLOSED 2026-08-28 (`GEO-23` step 2a, 12:00 slot).** The
> predicted removal below is now measured: with `coil_phantom_domain`'s rank-0
> build wrapped, rank 1 raises
> `RuntimeError: coil_phantom_domain geometry generation failed on rank 0
> (resolution=0.03); this is rank 1` and the `-n 2` command footers at
> **Status 1 in 3 s** where step 1 recorded **Status 124 at 120 s**, summary
> unchanged at `1 failed, 1 passed`
> (`20260828T170331Z_GEO-23-step2a-phantommetrics-n2.log`; `-n 1` control
> unchanged, `20260828T170323Z_…-n1.log`). **The entry stays OPEN for the
> geometry red itself** — `GEO-23` step 2b's sizing move is what retires it.
>
> **MEASURED 2026-08-28 (`GEO-23` step 1, 09:00 slot) — geometry-deterministic,
> with a measured floor one step away.** `-n 1` **Status 1, 3 s**
> (`1 failed, 1 passed in 1.17s`,
> `20260828T140352Z_GEO-23-step1b-phantommetrics-n1.log`); `-n 2` **Status 124,
> 120 s** with a complete summary then `MPI_Abort`
> (`20260828T140401Z_…-n2.log`) — so the abort is a property of the geometry
> and only the *hang* is rank-dependent (see the `birdcage_port_domain` entry
> above: wrapping the raise removes it). The `GEO-23` ladder, one process per
> rung at `-n 1`, shows `coil_phantom_domain` failing **only** at this
> fixture's own `resolution=0.03` and meshing at **0.024 (5 464 cells) /
> 0.0192 (9 330) / 0.01536 (16 177) / 0.012288 (28 485)** — monotone, one
> 0.8-step from green. **This site and
> `tests/materials/test_phantom_material_model.py` are the SAME call**:
> byte-identical `coil_phantom_domain` kwargs, so one sizing change would
> retire both reds. Entry stays OPEN (retires with step 2).
>
> **OWNER ASSIGNED 2026-08-27, 03:00 review: `GEO-23`** step 1 (b)/(c).
>
> **Where this fired.**
> `tests/post/test_phantom_field_metrics.py::test_phantom_field_metrics_and_exports_are_finite`,
> complex build, `FEM_EM_REQUIRE_COMPLEX=1`, `mpiexec -n 2`. Found by the
> `OPS-26` step 2 execution census (2026-08-26 19:30 implementer slot) — i.e.
> by the instrument built to find exactly this class, not by an example.
>
> **Literal symptom** (`20260827T004755Z_OPS-26-step2a-red-tb.log`, isolated,
> `--tb=short`):
>
> ```
> raise Exception(logger.getLastError())
> E   Exception: Invalid boundary mesh (overlapping facets) on surface 1 surface 1
> ```
>
> with the fragment census printed immediately before the abort:
> `[coil-phantom-mesh] fragment volumes=4 masses[m^3]: 1:1.381745e-04,
> 2:1.381745e-04, 3:2.261947e-04, 4:9.865456e-03 | air=4 coil_1=1 coil_2=2
> phantom=3`. **`1 failed, 1 passed in 1.24s`** — the failure is immediate and
> cheap, not a cost or JIT problem. The sibling test in the module,
> `test_evaluate_on_cells_fallback_skips_invalid_cell_point_pairs`, is green.
>
> **Second-order damage — this red burns the whole command.** After the
> failure the ranks diverge and teardown never completes: the isolated run
> printed its full pytest summary at 1.24 s and then sat until `timeout -k 30
> 200` killed it (**Status 124, elapsed 201 s**, PETSc error trailer after the
> summary line). In the census batch the same divergence hung the *next*
> module — `tests/post/test_phantom_phasor_semantics.py` never produced a
> result and the 900 s window was consumed
> (`20260827T003201Z_OPS-26-step2a-complex.log`, Status 124, 901 s). Anyone
> batching `tests/post` must expect to lose the tail of the command, not just
> this test. Budget it isolated.
>
> **Cause.** Not diagnosed. The symptom string is **identical** to the
> `GEO-21` entry immediately below, which measured the coarse end of
> `birdcage_port_domain`'s conductor sizing to have stopped meshing at the
> 0.11 merge. This occurrence is on a **different generator** (the
> coil+phantom fragment path, `air`/`coil_1`/`coil_2`/`phantom`), which is
> evidence that the 0.11 gmsh regression is **not birdcage-specific** — but
> that is a hypothesis stated from one shared error string, not a measurement.
> Whoever takes it should first check whether the fixture's resolution is at
> the coarse end of a continuum, as `GEO-21` step 1 found for the birdcage.
>
> **Not.** Not a tolerance or physics failure — no assertion is reached. Not
> the retired `DummyMagnetostaticSolver` red under §1 (that was an
> `AttributeError` and its code path was deleted by `TH-1`). Not diagnosed as
> real-mode-affected: the census only reached this module in the complex build,
> so the real-mode disposition of this test is **unmeasured**.
>
> **Filed, not fixed** — `OPS-26` step 2 is a census and files reds by
> construction (§9 item 1: "a red found here is filed, never fixed or
> re-recorded in-slot").
>
> **Retire-when:** the test is green through the harness on `main` in both
> builds, with the fixture's meshing resolution recorded.
>
> **Verified at** `18bb604` (tree clean at slot start; no source edited this
> slot).

### ✅ RETIRED 2026-08-28 (`GEO-23` step 1 (d), 09:00 slot) — ~~DEAD MODULE, filed 2026-08-27 (`OPS-26` step 2 leg (a)) — `tests/mesh/test_cylindrical_domain.py` collects **zero tests**~~

> **RETIRED by conversion, not deletion.** The module is now one asserting
> test, `test_cylindrical_domain_tag_volumes_partition_the_mesh`: the identity
> the old script only `print`ed, in quantitative form — the inner and outer tag
> volumes sum to the mesh volume at the shared helper's **1e-9** band (and the
> helper *reduces*, where the old `(ct.values == 1).sum()` counts were
> rank-local and never were) — plus an outer > inner ordering assertion that
> catches a tag swap the sum alone cannot see. `1 passed in 1.38s` at `-n 2`,
> and `tests/mesh --collect-only` now reports **58** where this module
> contributed **0** before
> (`20260828T141217Z_GEO-23-step1d-cylindrical-module.log`, Status 0, 5 s). The
> collection-time mesh build is gone with it — the build now happens inside the
> test. Its `resolution=0.02` was deliberately left unmoved: `GEO-23` step 1's
> ladder puts this generator's floor at 0.04-fails / 0.032-meshes, so 0.02 is
> comfortably inside the meshing range and is not that chunk's to change.
>
> ~~**OWNER ASSIGNED 2026-08-27, 03:00 review: `GEO-23`** step 1 (d) — convert
> to one asserting test or delete; either retires this entry.~~
>
> The file is a module-level *script*, not a test module: it calls
> `MeshGenerator.cylindrical_domain(...)` and `print`s at import time and
> defines no `test_*` function. `--collect-only` on the seven leg-(a)
> directories lists every other module in `tests/mesh` and **omits this one
> entirely** (`20260827T003050Z_OPS-26.log`, 189 collected, Status 0, 5 s) —
> so it contributes 0 to the census denominator while looking, by filename,
> like coverage.
>
> Worse than dead: pytest **imports** it during every collection of
> `tests/mesh`, so its mesh build runs — as collection-time work no
> disposition covers — and any exception it raises would surface as a
> collection error rather than a test failure.
>
> This is the `OPS-26` class in its purest form: a module no scheduled command
> can report on. **Filed, not fixed** — turning it into a real test (or
> deleting it) is a `mesh`-owning chunk's call, not the census's.
>
> **Retire-when:** the file either defines an asserting `test_*` function
> observed green through the harness, or is removed.
>
> **Verified at** `18bb604`.

### ✅ RETIRED 2026-08-25 (`EX-30` item 3 half A, 22:30 implementer slot) — `test_kwarg_off_reproduces_the_recorded_mesh` was **red on `main`**: the `GEO-16` kwarg-off cell record read 79 534, the 0.11 image meshes **79 070** (`EX-30` leg (mesh), 2026-08-25)

> **Where this fires.**
> `tests/mesh/test_two_torus_port_sheet.py::test_kwarg_off_reproduces_the_recorded_mesh`
> (`NCELLS_UNGATED_RECORD = 79_534`, line 78), and through the `ANS-1` import
> also `examples/meshing/04_two_torus_port_sheet.py` (`mesh:4`), which fails on
> the same constant.
>
> **Literal symptom** (`20260825T213632Z_EX-30-mesh-gate-probe.log`, `-n 2`,
> real, **`1 failed, 5 passed, 4 warnings in 42.06s`**, `Status: 1`):
>
> ```
> E   AssertionError: the default path meshed 79070 cells against the recorded 79534:
>     the opt-in sheet perturbed the mesh every gated PORT-1 / PORT-10 number was measured on
> E   assert 79070 == 79534
> ```
>
> **The assertion's own premise is *not* what broke.** The message blames the
> opt-in port sheet, and the sheet is innocent: two independent no-sheet builds
> in this leg agree exactly at **79 070** — `mesh:1`
> (`01_two_torus_ports.py`, which does not assert a cell count and ran **green**,
> `[mesh] 79070 cells built in 14.2 s`) and `mesh:4`'s own kwarg-off control
> (`[control] emit_port_sheet=False: 79070 cells in 13.6 s`) — while the
> *sheeted* build is a properly distinct 79 940. So the default path is
> self-consistent on 0.11 and it is the 79 534 record, measured on 0.7.2 in
> `20260817T003524Z_GEO-16.log`, that is stale. The five other assertions in the
> module pass, including the CAD port-interface area and the 0.970–0.980
> meshed-band cross-check the constant's own comment names as its guard.
>
> **Deliberately not re-recorded**, though `EX-30` leg (mesh) holds an in-class
> (1\*) example-record licence for moved cell counts: this constant lives in a
> **gate module**, not an example, and the licence does not reach it. Aligning
> the example's guide to 79 070 while the gate still asserts 79 534 would create
> exactly the example/gate divergence `ANS-1` exists to prevent, so the
> `mesh:1` guide's "79 534 cells / 12.9 s" was left standing too. The review
> owns the call: re-record the gate constant to 79 070 (and the `mesh:1` guide
> with it), or treat the 464-cell move as a regression to diagnose.
>
> **Verified at** `9b679d8`.
>
> **RULED 2026-08-25, 18:00 review — re-record, licensed and scoped to this
> constant.** `NCELLS_UNGATED_RECORD` 79 534 → **79 070**, version-tagged to
> the 0.11 image (gmsh 4.15.2), old digit in-comment citing both
> `20260817T003524Z_GEO-16.log` and this entry's probe log; the `mesh:1` and
> `mesh:4` guide copies move in the **same commit** so example and gate
> cannot diverge. Basis: two independent no-sheet builds agree exactly at
> 79 070 with the sheeted build properly distinct at 79 940 (the sheet is
> exonerated), the module's five other assertions are green including the
> 0.970–0.980 meshed-band guard the constant's comment names, and the −0.58%
> move is in family with the measured 0.11-gmsh mesh motion already ruled
> re-recordable under (1\*) (two-torus solve fixture −0.40%, `TH-10` −0.02%).
> Not a regression to diagnose: the record is a change-detector, and what it
> detected is the documented image change. No band moves. Landing is §9
> item 3 (18:00 queue). **Retire-when:** the commit that lands the re-record
> with the gate module and `mesh:4` green.
>
> **✅ RETIRED 2026-08-25, 22:30 implementer slot — landed exactly as ruled.**
> `NCELLS_UNGATED_RECORD` 79 534 → **79 070**, version-tagged to the 0.11 image
> with the old digit and both provenance logs in-comment; the `mesh:1`
> docstring + guide and the four `mesh:4` guide copies moved in the same commit.
> Green on all three anchors: the gate pair
> (`20260826T033222Z_GEO-16-rerecord-gate-pair.log`, `-n 2`, real, **5 passed
> in 55.84 s**, Status 0, elapsed 57 s) printing `[GEO-16 control] cells=79070`
> and the 0.970–0.980 meshed-band cross-check at **0.974490841**; `mesh:4`
> (`20260826T033350Z_GEO-16-rerecord-mesh4.log`, Status 0, 31 s) —
> `[mesh] 79940` sheeted, `[control] emit_port_sheet=False: 79070 cells in
> 13.9 s (record 79070)`, sheet tags absent; `mesh:1`
> (`20260826T033431Z_GEO-16-rerecord-mesh1.log`, Status 0, 16 s) —
> `[mesh] 79070 cells built in 14.1 s`. No band moved. Also re-recorded under
> the same (1\*) class: the `mesh:4` guide's sheet-facet count 84 → **82** and
> its wall-time/cell-count row, both un-asserted guide figures.

### ✅ RETIRED 2026-08-25 (`EX-30` item 3 half B, 22:30 implementer slot) — `mesh:5`'s **inverted control lost its separation**: the clamps-only mesh *cleared* the 0.755 floor it is asserted to fail, by 6e-6 (`EX-30` leg (mesh), 2026-08-25)

> **Where this fires.** `examples/meshing/05_region_resolution_policy.py`
> (`mesh:5`), line 255 — the `EX-18`/`EX-20` inverted-assertion pattern:
> `assert recovery["clamps_only"][tag] < POLICY_MIN_CAD_RECOVERY`.
>
> **Literal symptom** (`20260825T213601Z_EX-30-mesh-run-5.log`, `-n 2`, real,
> `Status: 1`, 7 s):
>
> ```
>   clamps_only  cells=   19618  mesh=  2.40 s
>   policy       cells=   20745  mesh=  2.65 s
> AssertionError: clamps-only mesh recovers 0.755006 of tag 1 (coil_1)'s CAD volume,
>   clearing the 0.755 floor the policy is supposed to be needed for
>   (on record: 0.754685 / 0.752565). The control no longer separates —
>   the premise needs re-examining, not the floor.
> ```
>
> **This is example-side only; the gate module is green.**
> `tests/mesh/test_mesh_tag_integrity.py` passed all four of its tests in the
> same probe run (`20260825T213632Z_EX-30-mesh-gate-probe.log`) because it
> asserts the floor **one-sidedly on the policy mesh** (`policy_volumes[tag] /
> cad_volume >= POLICY_MIN_CAD_RECOVERY`, line 248) and never asserts that the
> clamps-only control fails it. The stricter inverted claim exists only in the
> example, which is why only the example is red.
>
> **The margin is the whole finding: 0.755006 against a 0.755 floor** — the
> control clears it by 6.0e-6 relative, having sat at 0.754685 on 0.7.2, a
> 3.2e-4 move. The example's own comment anticipates this shape ("a baseline
> sitting at 0.949 would clear 'fails the gate' while saying nothing"), and
> `CONTROL_SEPARATION = 0.05` exists for the `mesh:3` fixture but has no
> counterpart here.
>
> **Deliberately not re-recorded and not widened.** The (1\*) licence covers
> moved cell counts and CAD masses, not a control's separation premise, and
> moving `POLICY_MIN_CAD_RECOVERY` or the record to recover the assertion would
> be loosening a gate — the assertion message itself says the premise is what
> needs re-examining. A `GEO-17` ruling: re-choose the clamps-only control so it
> fails by a stated margin, or retire the inverted claim.
>
> **Verified at** `9b679d8`.
>
> **RULED 2026-08-25, 18:00 review — re-choose the control, measure-first;
> demote to report only if no separating control exists.** The example may
> not carry an inverted assertion whose separation is 6e-6 on a record that
> moved 3.2e-4 under an image change. The landing (§9 item 3, second half)
> probes the clamps-only configuration at 2–3 coarser sizings (≈2.4 s per
> mesh) and adopts the coarsest sensible control that **fails the 0.755
> floor with ≥ 0.05 relative separation** — the `CONTROL_SEPARATION`
> precedent the example itself cites for the `mesh:3` fixture — old sizing
> in-comment, version-tagged. If no clamps-only sizing both fails the floor
> and remains a legitimate control, the inverted assertion is demoted to a
> printed report and the guide states the finding honestly: on 0.11 the
> clamps-only mesh already recovers the floor, so the policy's necessity is
> not demonstrable on this fixture at this floor. `POLICY_MIN_CAD_RECOVERY`,
> the one-sided gate-module assertion, and the records never move.
> **Retire-when:** the commit that lands `mesh:5` green under whichever
> branch the measurement selects.
>
> **✅ RETIRED 2026-08-25, 22:30 implementer slot — the re-choose branch, not
> the demotion.** A separating sizing exists, so the inverted assertion stays
> an assertion. The probe
> (`20260826T033622Z_GEO-17-mesh5-sizing-probe.log`, `-n 1`, 8 s) measured coil
> meshed/CAD at four uniform sizings — h = 0.015: 0.755006 / 0.750454 (margin
> **−0.000006**, the red); **h = 0.018: 0.649812 / 0.648431 (margin +0.105188 /
> +0.106569, SEPARATES)**; h = 0.020: 0.595547 / 0.579713; h = 0.025:
> 0.471986 / 0.510423. `CONTROL_RESOLUTION = 0.018` adopted as the first
> candidate that separates, with the whole table in-comment; the probe stopped
> there rather than hunting a margin.
>
> The control is a **third build**, not a re-pointed one: `UNIFORM_VOLUMES_RECORD`
> is a `GEO-17` gate constant measured at h = 0.015, so moving the clamps-only
> sizing would have broken negative control (a) — the `OPS-17` table
> reproduction at 1e-9 — which is not licensed to move. The clamps-only mesh
> therefore stays as that reproduction and as the baseline the refine/coarsen
> sign identities read against; the new `coarse_control` build carries only the
> inverted assertion, now gated at `CONTROL_SEPARATION = 0.05` rather than on
> the bare `<` that went red. `SIZING_SEPARATION` is asserted against **both**
> baselines (+0.078411 / +0.085109 against clamps-only, +0.183605 / +0.187132
> against the coarse control) so keeping the tighter of the two is not quietly
> dropped. `POLICY_MIN_CAD_RECOVERY`, the one-sided gate-module assertion and
> every record are untouched; `tests/mesh/test_mesh_tag_integrity.py` was not
> edited at all.
>
> Green: `20260826T033758Z_GEO-17-mesh5-control-rechoice.log` and the
> post-doc-edit confirm `20260826T033959Z_GEO-17-mesh5-confirm.log` (both
> `-n 2`, real, Status 0, 8 / 9 s), printing the control's failing margin.

### ✅ CLOSED 2026-09-01 (`OPS-30`, 21:00 implementer slot) — two `scripts/probes/` scripts were **never migrated to dolfinx 0.11**: they construct `fem.petsc.LinearProblem` without 0.11's required `petsc_options_prefix`

> **Where this fires.** `scripts/probes/mag13_step2b_recovery.py:180` and
> `scripts/probes/post3_step3_debug.py:55`. Neither is a test, neither is an
> example, and **nothing scheduled runs either of them** — which is exactly the
> shape `OPS-26` exists to enumerate. They are one-off diagnostic probes kept
> for their write-ups.
>
> **Literal symptom** (static, from the sweep — these scripts were not
> executed; `20260825T200918Z_OPS-26.log`, `Status: 1`, 4 s):
>
> ```
> scripts/probes/mag13_step2b_recovery.py:180: [missing-required]
>   dolfinx.fem.petsc.LinearProblem — required parameter
>   'petsc_options_prefix' not supplied
> scripts/probes/post3_step3_debug.py:55: [missing-required] (same)
> ```
>
> 0.11 made `petsc_options_prefix` a keyword-only parameter with no default, so
> each line raises `TypeError` on the first call. `src/` and `tests/` carry
> **zero** such sites (log `20260825T200851Z_OPS-26.log`, 434 resolved call
> sites over 29 APIs, `violations=0`) — the migration is complete everywhere a
> gate can reach; these two are outside that reach.
>
> **Not diagnosed further, and deliberately not fixed** — `OPS-26` files
> survivors rather than repairing them in-slot, and a probe script's value is
> its recorded output, not its ability to re-run. Fixing them is a two-line
> change whenever someone wants one of them again.
>
> **Verified at** `e26c128` (tree state of the 15:00 slot).
> **Scoped 2026-09-01 18:00 review as `OPS-30`, §9 item 2** — a
> `petsc_options_prefix` at each of the two sites and nothing else, landed in
> the same commit as the survivor-set pin below (which goes red in *either*
> direction, so the two cannot move separately). Anchor: the migration sweep,
> which reads `violations=0` over 434 resolved call sites in `src/` + `tests/`
> and exactly **2** at the two named `file:line` sites today, must read **0**
> at those two with the `src/` + `tests/` count unchanged.
>
> **Retire when** both sites take `petsc_options_prefix`, or the scripts are
> deleted; `tests/environment/test_dolfinx_api_migration.py::test_filed_survivors_outside_the_gated_roots_are_unchanged`
> pins the survivor set at exactly these two and goes red in **either**
> direction, so this entry cannot rot silently.
>
> **Closed 2026-09-01 21:00 slot by `OPS-30`.** Both sites now pass
> `petsc_options_prefix` (`fem_em_probe_mag13_step2b_` /
> `fem_em_probe_post3_step3_`) and nothing else changed — the probes were not
> modernised or re-run, their value is still their recorded output. Count
> identity in both directions, smoke, `-n 1`: the `examples` + `scripts` sweep
> read **2** violations at the two named sites pre-fix and **0** after, with
> its census unchanged at **82 files / 320 resolved call sites / 22 APIs**,
> while `src` + `tests` stayed **177 / 484 / 30** at `violations=0`
> (pre-fix `20260902T020122Z_OPS-30.log`, `SURVIVOR_STATUS=1`, Status 0, 12 s;
> post-fix `20260902T020238Z_OPS-30.log`, `SURVIVOR_STATUS=0`, `3 passed`,
> Status 0, 37 s). The negative control was measured live pre-fix rather than
> read from the August log, so the sweep is confirmed to reach `scripts/`.
> The pin moved in the same commit and was **strengthened rather than
> deleted**: `FILED_SURVIVORS` is now the empty set, guarded by census floors
> over `examples` + `scripts` and a reachability floor of ≥ 2 resolved
> `dolfinx.fem.petsc.LinearProblem` sites under `scripts/probes`, so "no
> violations" can no longer be produced by a sweep that misses the two files.
> The `src` + `tests` figures here (434 sites / 29 APIs, measured 2026-08-25)
> are the August numbers; the tree has since grown to 484 / 30 — the gate
> asserts floors, not equalities.

### ✅ CLOSED 2026-08-25 (`MAG-19` step 2, 21:00 implementer slot) — `tests/validation/test_convergence.py::TestConvergence::test_h_refinement_straight_wire` was **red on `main`**: the fitted rate is **1.90** against the `MAG-13` band `[0.7, 1.5]`, because the finest rung's error collapsed 9.26% → **4.4605%** on the 0.11 image

> **Where this fires.** `tests/validation/test_convergence.py::TestConvergence::test_h_refinement_straight_wire`,
> real build, on **`main`** at `878fa3e`. Not a worksite, not an example path:
> the gate itself. It was found from the example side — `examples/magnetostatics/06_h_convergence_rate.py`
> imports `RATE_MIN`, `RATE_MAX`, `RESOLUTIONS`, `solve_h_refinement` and
> `fit_convergence_rate` from this very module (the `ANS-1` rule, already
> applied), so `-e 6` runs the gate's computation — and the gate was then
> probed directly to confirm it, rather than inferred.
>
> **Literal symptom** (`20260825T141636Z_EX-30-root-mag6-gate-probe.log`,
> `1 failed in 143.11s`, `Status: 1`, `-n 2`, both rank footers identical):
>
> ```
> Convergence rate: 1.90
> Expected rate for linear elements: ≈1.0
> AssertionError: Convergence rate 1.90 outside [0.7, 1.5] (expected ≈1.0 for
>   N1curl degree 1); errors [0.21841667267163878, 0.15384842035994292,
>   0.04460534278989355] at h [0.004, 0.0025, 0.0018]
> ```
>
> **What moved, rung by rung.** The band is unchanged; the error ladder is not.
> Against the `MAG-13` record (`20260730T125522Z_MAG-13.log`, 0.7.2):
>
> | h (m) | cells (0.11) | rel L2 (0.11) | rel L2 on record | move |
> | --- | --- | --- | --- | --- |
> | 0.0040 | 38 740 | 21.8417% | 22.19% | −0.35 pp |
> | 0.0025 | 147 235 | 15.3848% | 12.75% | +2.63 pp |
> | 0.0018 | 383 146 | 4.4605% | 9.26% | **−4.80 pp** |
>
> The middle rung is exactly the 147 235 cells / 15.3848% the retired
> `MAG-18` entry below already names, so two of the three rungs are accounted
> for by the documented 0.11 gmsh mesh motion. The **finest rung is not**: a
> better-than-halved error is what levers the fitted slope from 1.10 to 1.90.
> The sequence is still monotone, so the example's own negative control
> (monotone decay) passes — it is the *rate* that breaks, on the upper edge
> the gate's docstring says has teeth precisely because "a rate well above 1.5
> means one resolution in the sequence is anomalous".
>
> **Update 2026-08-25 (`MAG-19` step 1, 13:30 implementer slot) — the
> discriminating measurement is done, and it selects *neither* pre-stated
> reading.** Both norms were run on the *same* four solves (added rung
> h = 0.0030) via `tests/validation/probe_straight_wire_dual_norm.py`;
> log `20260825T183555Z_MAG-19-step1-dualnorm-fits.log`, 160 s at `-n 2`,
> Status 0. The red reproduces digit for digit (three original rungs within
> 1.321e-06 relative of the row above; sampled three-rung fit **1.9038**) and
> so does the `E_Ω` negative control through the imported machinery (fit
> **1.6854**; the h = 0.0025 `E_Ω` record within **2.094e-08**) — so the
> `ANS-1` import is sound and what moved is the measurement.
>
> | h (m) | cells | sampled 10-pt | `E_Ω` |
> | --- | --- | --- | --- |
> | 0.0040 | 38 740 | 21.841675% | 25.286827% |
> | 0.0030 | 88 018 | 18.473177% | 14.288381% |
> | 0.0025 | 147 235 | 15.384843% | 10.617170% |
> | 0.0018 | 383 146 | 4.460528% | 6.645807% |
>
> Sampled pairwise rates 0.5822 / 0.7456 / 1.9894 / 1.0034 / 2.7819 / 3.7690;
> `E_Ω` 1.9843 / 1.8464 / 1.6735 / 1.6288 / 1.4985 / 1.4261. Reading (a) fails
> because the sampled ladder carries a **second** out-of-band pair that avoids
> h = 0.0018 (0.004→0.003 at 0.5822, on the rung (a) would promote); reading
> (b) fails because the sampled norm is *not* unstable everywhere — dropping
> h = 0.0018 alone returns the fit to **0.7309**, inside the band. The new
> constraint: `E_Ω` is the stable instrument but sits **above 1.5 on every
> subset** (fits 1.6661–1.8588), so transferring the two-sided [0.7, 1.5] onto
> it would be red on arrival — its live gate is one-sided ≥ 0.7, which it
> meets on 6/6 pairs. **Nothing was re-recorded and no band was moved**;
> `MAG-19` stays 🟡 with three options for the review in its §7 entry. This
> entry stays open: the gate is still red on `main` for the reason above.
>
> **Cause: not diagnosed.** Two readings fit and this slot did not
> discriminate them: (a) the h = 0.0018 rung's mesh moved enough on 0.11 that
> its sampled 10-point error is anomalous, in which case the sequence needs
> re-choosing the way `MAG-13` excluded h = 0.0035; (b) the sampled 10-point
> norm is the wrong instrument on 0.11 and the `MAG-18` `E_Ω` norm — which
> **is** green on 0.11 at rate 1.6854 (see the retired entry below) — is the
> one to gate. Deciding between them is a `MAG-13`/`MAG-18` question, not an
> `EX-30` one; **nothing was re-recorded and no band was moved.**
>
> **Owning chunk:** `MAG-13` (the band and the resolution sequence), with
> `OPS-18` as the image that moved the ladder under it. **Consequence for
> `EX-30`:** `examples/magnetostatics/06_h_convergence_rate.py` exits 1 for
> this reason and nothing else — its artifact
> (`h_convergence_rate_combined.xdmf`) is written *before* the assertion and
> did refresh, so it is not in the stale census.
>
> **RULED 2026-08-25, 18:00 review — option (i): the rate duty transfers to
> the `E_Ω` ladder under `E_Ω`'s own one-sided ≥ 0.7 criterion, which
> `MAG-18` already gates live and green on 0.11.**
> `test_h_refinement_straight_wire` keeps its monotone-decay assertion and
> prints the error table as a report; the two-sided [0.7, 1.5] on the sampled
> 10-point statistic is **retired with its basis stated**, not widened: that
> statistic swings 34% of its own value under its sampler (measured on both
> images, `OPS-18` step 3 attempt 5) and the band already failed on 0.7.2 at
> `n_points = 8` — the gate was passing on a sampler choice, and `MAG-18`
> built `E_Ω` precisely to replace it. On option (iii)'s question: the upper
> edge is **not** re-imposed on `E_Ω` — no two-sided band has ever been
> validated on 0.11 for either statistic, under-convergence (the ≥ 0.7 side)
> is the failure mode a rate gate exists to catch, and a superconvergence
> guard, if ever wanted, must be commissioned with its own measured basis
> (a weekly-review question, deliberately not opened here). Option (ii)
> rejected: a 0.03 margin on a statistic with a 34% sampler swing is not a
> gate. Landing is `MAG-19` step 2 (§9 item 2, 18:00 queue); `mag:6` imports
> `RATE_MIN`/`RATE_MAX`/the fit from this module, so the landing reconciles
> the example to the transferred duty in the same commit (old text
> in-comment) and runs `-e 6` as the consumer check. **Retire-when:** the
> commit that lands `MAG-19` step 2 green including `mag:6`.
>
> **RETIRED 2026-08-25, 21:00 implementer slot — the retire-when is met, in
> one commit, on four logged runs at `-n 2`, real build, on `main`.**
> (1) *The red reproduced first*, before anything was edited:
> `20260826T020124Z_MAG-19-step2-red.log`, Status 1, 145.27 s — 21.8417% /
> 15.3848% / 4.4605% at 38 740 / 147 235 / 383 146 cells, rate **1.90**, the
> same digits `20260825T141636Z` and `MAG-19` step 1 recorded. (2) *The
> disposition's own green*: `20260826T020508Z_MAG-19-step2-green.log`,
> `1 passed` / Status 0 / 142.36 s on **bit-identical** errors (the fit still
> prints, at 1.9038, as a report beside the retired band and the duty owner).
> (3) *Negative control — `MAG-18`'s gate module green **untouched***, zero
> edits to `test_straight_wire.py`: `20260826T020739Z_MAG-19-step2-mag18.log`,
> `7 passed` / Status 0 / 362.68 s, `E_Ω` 25.2868 → 10.6172 → 6.6458% at fitted
> rate **1.6854 ≥ 0.7**, the h = 0.0025 record 1.0617170177e-01 and the
> natural-BC ratio 0.3285 — all three reproducing the 2026-08-23 re-gate, so
> the duty was transferred to a gate that is executing and green, not to a
> claim. (4) *Consumer check*: `20260826T021403Z_MAG-19-step2-e6.log`,
> `-e 6 -n 2`, Status 0 / 148 s, "All assertions hold", printing the same three
> errors, the report-only rate and the duty owner.
>
> **No band moved anywhere.** `RATE_MIN` / `RATE_MAX` keep their values and
> their names; what changed is that nothing on this ladder asserts on them.
> **One residual, filed rather than fixed:** `test_straight_wire.py::TestStraightWire::test_straight_wire_convergence`
> still gates a *two-rung, 8-point sampled* fit on the same two-sided band, and
> it is green (fitted **0.7900** in run (3)). It was outside `MAG-19`'s scope
> and inside the module this landing had to leave untouched, so it was left
> alone and is named in-comment at the constants. Whether the ruling's "no
> upper edge on a sampled statistic" extends there is a review question, not an
> implementer's.
>
> **RULED 2026-08-26 03:00 review: commissioned as `MAG-20`** (§7,
> measure-first, not queued this interval). Ruling (i) is *not* inherited by
> fiat onto a currently-green test — the residual gets its own `n_points`
> sweep on its own two rungs with a pre-stated decision rule: retire the
> two-sided band under the ruling-(i) pattern if the fit crosses either edge
> under the sampler; keep it, recorded as *validated*, if it is stable at
> every count. Either outcome is a measurement.
>
> **RESIDUAL CLOSED 2026-08-28 (`MAG-20` step 1, 00:00 implementer slot) — the
> band is KEPT, validated by measurement.** The sweep on this test's own two
> rungs and its own 0.4 R window fits **0.7900 / 0.7246 / 0.9934** at
> `n_points` 8 / 10 / 20 — **no crossing of either edge of [0.7, 1.5]** — so the
> pre-stated rule's *keep* branch fires and nothing moved: `RATE_MIN`/`RATE_MAX`
> unchanged, no assertion added or removed, the disposition is an in-comment
> measurement record at the assertion. Probe
> `tests/validation/probe_straight_wire_convergence_npoints.py`, log
> `20260828T050130Z_MAG-20-step1-npoints-probe.log` (49 s, `-n 2`); anchor
> `test_straight_wire.py` `7 passed / 369.95 s / Status 0`
> (`20260828T050256Z_MAG-20-step1-anchor-module.log`) with `E_Ω` 1.6854 and the
> h = 0.0025 record 1.0617170193e-01 untouched. **Two findings left open for the
> review, not defects on `main`:** the sampler swing on this window is 6–7% of
> the error (vs 34% on the 0.8 R window), but it still moves the rate by 37% of
> its own value, and the n = 10 row clears `RATE_MIN` by only **0.0246**. The
> test is green at all three counts; the thin margin is recorded in the
> `MAG-20` §7 entry as findings 45–46.

### ✅ RETIRED 2026-08-30 (weekly planning review — size-field licence **denied**, the `GEO-23` step-2a wrap ruled sufficient) — ~~`MeshGenerator.straight_wire_domain` has a **coarse-resolution floor on the dolfinx 0.11 image**: `resolution = 0.01` aborts inside gmsh with duplicated facets for every geometry tried, `0.008` and finer mesh, and the threshold between them is unbisected — no guard exists, so a too-coarse request still fails illegibly~~ (originally 🔴 OPEN 2026-08-25, re-headed 2026-08-26 by `EX-30` leg (root))

> **RULED 2026-08-30, weekly planning review — the size field is NOT
> licensed into `src/`, and this entry retires on the wrap.** The 08-28
> ruling left one question here: whether the gmsh Distance→Threshold size
> field that `GEO-22` step 2/2c measured (18/18 OK, 0 fallbacks, 19 823 vs
> 21 830 cells at `h = 0.008`, asserted at `-n 1` and `-n 2`) should land in
> `straight_wire_domain`. It should not, for three reasons that are each
> sufficient: (1) it moves four **Phase-1** records for no physics gain —
> `mag:1`'s 21 830 and the three `MAG-13`/`MAG-18` ladder records — and §10's
> element-order note already rules that Phase 1 is complete and not worth
> re-gating; (2) the defect the entry names — a too-coarse request failing
> *illegibly* — is fixed: `GEO-23` step 2a's raise path lands the failure
> on every rank in seconds, and `GEO-22` step 1 proved there is no floor to
> guard (the `[0.008, 0.010]` sweep is non-monotone and bit-reproducible),
> so a guard would be a fiction; (3) the example that fired this moved to
> `resolution = 0.008` on 2026-08-25 and has been green since. The
> size-field probe and its asserted gate (`tests/mesh/test_straight_wire_size_field_probe.py`)
> stay as the on-record measurement of *why* the coarse rung fails (a
> fallback triangulation under the unfielded mesher), available to any
> future generator that needs it — and the first place to look if a
> **new** geometry meets the same string is `GEO-23`'s classification, not
> this entry. Nothing re-recorded, no band moved, `src/` untouched.

> **Where this fires.** `./run_examples.sh -e 1`, real build, on **`main`** at
> `878fa3e`. The crash is in the mesh generator, before any solve, so the
> example produces nothing: its **seven** repo-root artifacts
> (`straight_wire_{A,B}.bp`, `straight_wire_{A,B,B_analytical}.xdmf`,
> `straight_wire_combined.{h5,xdmf}`) are the entire remainder of leg
> (root)'s stale census.
>
> **Literal symptom** (`20260825T140159Z_EX-30-root-run-mag1to2.log`,
> `Status: 124` — the 124 is a post-`MPI_Abort` teardown hang against the
> runner's `-t 300`, not a compute overrun; the failure itself is immediate):
>
> ```
> Info    : Reconstructing mesh...
> Info    :  - Creating surface mesh
> Info    : Found two duplicated facets.
> Info    :   1st: [145, 20, 216] #1
> Info    :   2nd: [145, 20, 216] #1
> Error   : Invalid boundary mesh (overlapping facets) on surface 1 surface 1
> ...
>   File "/workspace/src/fem_em_solver/io/mesh.py", line 304, in straight_wire_domain
>     gmsh.model.mesh.generate(3)
> Exception: Invalid boundary mesh (overlapping facets) on surface 1 surface 1
> ```
>
> **`MeshGenerator.straight_wire_domain` is not broken in general** — the
> *same* generator meshed three times in the same slot, in the run immediately
> after, at 38 740 / 147 235 / 383 146 cells
> (`20260825T141141Z_EX-30-root-run-mag6.log`). The parameter sets differ, and
> the example's is exercised by **no** gate:
>
> | | wire_length | domain_radius | resolution |
> | --- | --- | --- | --- |
> | `01_straight_wire.py:118-121` | 0.3 | 0.04 | **0.01** |
> | `test_straight_wire.py:62-65` | 0.20 | 0.03 | 0.0025 |
> | `test_convergence.py:41-61` | 0.2 | 0.03 | 0.004 / 0.0025 / 0.0018 |
>
> **Cause, localised — it is `resolution` alone, and the geometry is
> irrelevant.** `tests/validation/probe_straight_wire_mesh_resolution.py`
> walks the two axes separately
> (`20260825T142512Z_EX-30-root-mag1-mesh-probe.log`, `Status: 0`, **29 s**,
> `-n 1`; the probe documents why `-n 1` — a gmsh throw on rank 0 deadlocks
> the collective `_model_to_mesh` on the others):
>
> ```
> Leg A -- the example's geometry (L = 0.3, R = 0.04), resolution swept:
>   h = 0.0100  FAIL  Invalid boundary mesh (overlapping facets) on surface 1 surface 1 (0.3 s)
>   h = 0.0080  OK       21830 cells (2.6 s)
>   h = 0.0060  OK       34250 cells (4.2 s)
>   h = 0.0050  OK       55306 cells (7.0 s)
>   h = 0.0040  OK       98778 cells (13.0 s)
>
> Leg B -- the example's resolution (h = 0.01), geometry stepped to the gate's:
>   L = 0.30  R = 0.040 FAIL  ...
>   L = 0.30  R = 0.030 FAIL  ...
>   L = 0.20  R = 0.040 FAIL  ...
>   L = 0.20  R = 0.030 FAIL  ...
> ```
>
> So: **`h = 0.01` is the only failing rung, it fails for every geometry
> tried including the gate's own `L = 0.20 / R = 0.030`, and everything from
> `h = 0.008` down meshes.** This is therefore *not* "the example's box is
> unusual" — it is a coarse-resolution floor in `straight_wire_domain` on the
> 0.11 image, and the gate's geometry would hit it too if any gate ever ran
> that coarse. None does, which is why only the example sees it. The failure
> is instant (0.2–0.3 s) and perfectly reproducible across all five cases.
>
> **Still not diagnosed:** *why* 0.01 specifically. The obvious suspect is the
> wire cylinder — `resolution = 0.01` is 1.67× the wire *diameter*
> (2·0.003 m), so the wire surface cannot carry a well-formed facet loop —
> but `h = 0.008` is still 1.33× that diameter and meshes fine, so the
> threshold is not simply "coarser than the wire". Bisecting 0.008–0.010 was
> not run. `straight_wire_domain` is untouched since the 0.11 API migration
> (`d176bc1`, `OPS-18` step 2), so this is image behaviour, not a repo
> regression.
>
> **Owning chunk:** unassigned — a review call, and the localisation above
> makes it a cheap one. Two fixes are available: move the example off
> `resolution = 0.01` (0.008 costs 2.6 s of meshing and is the nearest
> working rung), or give `straight_wire_domain` a guard/clamp so a
> too-coarse request fails legibly instead of inside gmsh. `EX-30` leg (root)
> has no licence for either. **Nothing was re-recorded.**
>
> **RULED 2026-08-25, 10:30 review: the example moves to
> `resolution = 0.008`** (old 0.01 kept in-comment citing the probe log
> `20260825T142512Z_EX-30-root-mag1-mesh-probe.log`), executed inside
> `EX-30` leg (root)'s completion slot — owning chunk is now `EX-30`.
> **No guard is written in-slot**: the threshold is unbisected
> (0.008 works, 0.010 fails, the boundary between them is unmeasured), and
> a guard constant without a measured boundary would encode a guess. The
> coarse-resolution floor itself stays on record **in this entry** as
> documented 0.11-image behaviour after the example fix lands — the entry
> then re-heads as the floor finding (the example symptom retired), and
> retires fully only when a measured-threshold guard lands in
> `straight_wire_domain` or the upstream image moves the floor.
>
> **EXAMPLE SYMPTOM RETIRED 2026-08-26 (`EX-30` leg (root) completion, 12:00
> implementer slot); the floor finding above stays OPEN.** The ruling was
> executed exactly as written: `01_straight_wire.py:120` moved
> `resolution` `0.01` → **`0.008`**, with the old value and the probe's full
> reasoning in-comment at the constant. `-e 1 -n 2` is now **green from
> `main`** (`20260826T170155Z_EX-30-root2-run-mag1.log`, Status 0, **9 s**,
> real build) and meshes at **21 830 cells / 4 662 vertices** — the probe's
> `h = 0.0080 OK 21830 cells` reproduced exactly, which is the confirmation
> that the localisation was right and not a coincidence of that probe's
> geometry. Closed forms unmoved and reproduced: `B(3 mm)` analytic
> `6.666667e-05 T` = `μ₀I/2πr`, analytic decay ratio `B(3 mm)/B(38 mm)` =
> **12.67** = 38/3. The example's derived figures moved with the mesh
> (relL2 65.8739% → **51.9781%**, max rel 85.2498% → **76.7330%**, numerical
> decay 29.83 → **20.31**, energy 2.307201e-08 → **2.630243e-08 J**) and are
> re-recorded version-tagged in `01_straight_wire.md` under the leg's (1\*)
> guide-table licence, old digits in-comment. The seven `straight_wire_*`
> artifacts cleared: census `stale=7` → **`stale=0`**
> (`20260826T170118Z_…-precensus.log` → `20260826T171345Z_…-postcensus.log`,
> `dead=0 guide=0 exit=0`).
>
> **What is still open, and it is the whole reason this entry survives:** no
> guard was written, so `straight_wire_domain(resolution=0.01)` still aborts
> inside gmsh with `Invalid boundary mesh (overlapping facets)` rather than
> raising something a caller can read. The threshold in `[0.008, 0.010)`
> remains unbisected and *why* 0.01 specifically remains undiagnosed (the
> wire-diameter hypothesis is contradicted by 0.008 working at 1.33× the
> diameter). **Retire-when:** a measured-threshold guard lands in
> `straight_wire_domain`, or the upstream image moves the floor and that is
> measured. **Owning chunk:** unassigned — `EX-30` owned only the example
> fix, and it is done.
>
> **OWNER ASSIGNED 2026-08-26, 18:00 review: `GEO-22`** — bisect
> `[0.008, 0.010)` on both geometries with the existing `-n 1` probe, land a
> `ValueError` guard at the *measured* `h_ok`, gate it with the exception
> type, the probe cell count and `mag:1`'s unmoved 21 830; the straight-wire
> gate ladders are the negative control. Full rubric in the §7 entry;
> queued as a §9 spare. This entry retires when that guard lands.
>
> **RE-HEADED 2026-08-28 (`GEO-22` step 1, 07:30 implementer slot) — there is
> no floor, and no guard is possible: the failure is NON-MONOTONE in
> `resolution` on both geometries, deterministically so.** The bisection ran
> as commissioned and returned a measured negative. `probe_…_resolution.py`
> gained a leg C that sweeps the whole open interval `[0.008, 0.010]` on a
> uniform **2.5e-4** grid — nine rungs — on the example's geometry *and* the
> gate's, at `-n 1` (`20260828T123115Z_GEO-22-step1-bisect.log`, `Status 0`,
> **23 s**; repeat `20260828T123205Z_GEO-22-step1-bisect-repeat.log`,
> `Status 0`, **22 s**). A uniform sweep was run rather than a true bisection
> deliberately: it costs about what three bisection steps would and it is the
> only form of the measurement that can *see* non-monotonicity. It did.
>
> | `resolution` | example (L = 0.3, R = 0.04) | gate (L = 0.20, R = 0.030) |
> | --- | --- | --- |
> | 0.00800 | OK 21 830 | OK 8 262 |
> | 0.00825 | OK 18 745 | OK 8 004 |
> | 0.00850 | OK 17 644 | OK 7 755 |
> | 0.00875 | **FAIL** | **FAIL** |
> | 0.00900 | OK 14 709 | **FAIL** |
> | 0.00925 | **FAIL** | OK 6 894 |
> | 0.00950 | OK 17 683 | OK 6 768 |
> | 0.00975 | **FAIL** | OK 12 200 |
> | 0.01000 | **FAIL** | **FAIL** |
>
> Every failing cell is the same literal `Invalid boundary mesh (overlapping
> facets) on surface 1 surface 1`, in 0.2–0.3 s, as before.
>
> **Three findings, in order of consequence.**
>
> 1. **The failing set is interleaved, not a floor.** `h = 0.00875` fails on
>    both geometries while the *coarser* 0.00900 (example) and 0.00925 /
>    0.00950 / 0.00975 (gate) mesh. So "everything from 0.008 down works,
>    everything above fails" — the reading this entry has carried since
>    2026-08-25 — is **false**; it was an artefact of the old leg-A ladder
>    sampling only 0.010 and 0.008 and nothing between. There is no threshold
>    to encode, so `GEO-22` step 1's pre-registered stop condition fires and
>    **no guard was written.** A `resolution > RESOLUTION_FLOOR` guard at any
>    constant would either reject meshing rungs or admit failing ones.
> 2. **It is deterministic.** The two runs above are independent invocations
>    and reproduce **bit-identically** — same OK/FAIL at all 18 cells, same
>    cell count to the digit. So this is not run-to-run instability in gmsh's
>    randomised insertion; it is a reproducible function of `(geometry,
>    resolution)`. That matters for step 2: a retry-with-jitter fix would have
>    to *perturb* the request, not merely repeat it.
> 3. **The cell count is non-monotone in `h` too, and by a lot.** Example:
>    0.00900 → 14 709 but the *coarser* 0.00950 → 17 683 (+20%). Gate:
>    0.00950 → 6 768 but 0.00975 → **12 200**, a **1.80×** jump for a coarser
>    request. So the mesher's whole response to `resolution` is discontinuous
>    in this band, and the failures are the visible part of that.
>
> **Mechanism, localised but not diagnosed.** Every rung in the sweep — the
> meshing ones included — prints `[ 0%] NNN triangles are equivalent` on
> surface 1 (the wire cylinder) and then falls back
> `Frontal-Delaunay` → **`MeshAdapt`** for that surface alone. So
> Frontal-Delaunay is producing coincident triangles on the wire at every
> size in this band; whether gmsh's fallback then yields a boundary the 3D
> reconstruction accepts is what varies rung to rung. That is consistent with
> all three findings and with the wire-diameter suspicion being *near* the
> mark without being the mechanism (0.008 is 1.33× the 0.006 m diameter and
> also triggers the fallback — it just survives it). **Not chased in-slot:**
> `GEO-22` step 1's scope explicitly excludes diagnosing gmsh.
>
> **Not run, deliberately:** `mag:1` and the straight-wire gate ladders were
> the negative controls *for a guard*. No guard landed, no `src/` line
> changed, so there is nothing for them to control and spending 370 s on them
> would prove only that an unmodified generator is unmodified. The example's
> own 0.008 rung is reproduced inside this very sweep at **21 830 cells**,
> which is the `EX-30`/`mag:1` record to the digit — that is the control that
> was worth having, and it is free.
>
> **Retire-when — restated, because the old condition is now unreachable:** a
> *measured-threshold* guard cannot land. This entry retires when either (a)
> `straight_wire_domain` stops emitting the coincident-triangle surface mesh
> (an upstream image move, or a generator change that meshes the wire surface
> differently — e.g. an explicit size field on the cylinder rather than a
> global `resolution`), or (b) a `GEO-22` step 2 lands a guard of a *different
> shape* — a post-mesh validity check, or a documented allowlist of verified
> rungs — and a review rules that shape sufficient. **Owning chunk:**
> `GEO-22`, step 1 done as a measured negative, step 2 a review's call.
>
> **RULED 2026-08-28, 10:30 review — guard shape.** The allowlist is
> **rejected** (nine rungs on two geometries is a sample, not a truth). The
> post-mesh wrap is **adopted** and lands as part of `GEO-23` step 2a (§9
> item 1): `straight_wire_domain` gets the `birdcage_port_domain` raise path
> (catch on the building rank, `bcast`, raise on every rank) so a failing
> rung footers in seconds at `-n 2` instead of deadlocking, with a gate at
> `h = 0.00875` asserting the raise on every rank. The size field (branch
> (a) above) is the only candidate that could *fix* the fallback and would
> move `mag:1`'s 21 830 and the three ladder records — its *licence* is the
> 2026-08-30 weekly review's; its *measurement* is §9 item 5, a no-`src/`
> probe leg predicting 0 fallbacks / 18 of 18 OK. This entry retires with
> the step-2a wrap **only if** the review then rules the wrap sufficient;
> otherwise it stays open pointing at the size-field decision.
>
> **WRAP LANDED 2026-08-28 (`GEO-23` step 2a, 12:00 slot) — the adopted half is
> done; the entry stays OPEN.** `straight_wire_domain` now carries the shared
> `_raise_geometry_failure_on_every_rank` raise path, and the gate is
> `tests/mesh/test_geometry_failure_is_collective.py`: the `mag:1` example
> geometry at `h = 0.00875` raises on **every** rank (the caught flag is
> `allreduce`d, and rank 1's message is asserted to name both the generator and
> the resolution), `1 passed in 0.91s` at `-n 2` inside a `-k 30 60` window
> (`20260828T170254Z_GEO-23-step2a-gate-n2.log`, Status 0, 3 s; `-n 1`
> `20260828T170247Z_…-n1.log`). The wrap moved no mesh — `mag:1` re-reads
> **21 830 cells** and `B(3 mm) = 6.666667e-05 T`
> (`20260828T170414Z_GEO-23-step2a-control-mag1.log`). What is **not** done is
> the other clause of the restated done-when: the size-field probe table (§9
> item 5) is unrun, and no fallback line count changed, so the coincident-
> triangle emission this entry is about is untouched. Retirement still needs a
> review's ruling that the wrap is sufficient, or the size field.
>
> **SIZE-FIELD PROBE RUN 2026-08-29 (`GEO-22` step 2, 07:30 implementer slot) —
> the hypothesis is CONFIRMED: a wire-surface size field removes the fallback
> entirely and every rung meshes. The entry stays OPEN; nothing landed in
> `src/`.** `probe_straight_wire_mesh_resolution.py` gained a **leg D** that
> re-runs leg C's nine rungs on both geometries with a gmsh
> `Distance`/`Threshold` field anchored on the generator's own `wire_surface`
> physical group — `SizeMin = wire_radius = 0.003`, `SizeMax` = the rung's own
> `h`, `DistMin = 0.003`, `DistMax = 0.006` — installed by patching
> `gmsh.model.mesh.generate` for the duration of one call (so geometry,
> fragment, physical groups, the raise path and `_model_to_mesh` are all the
> shipped code) with `Mesh.MeshSizeFromPoints` /
> `…ExtendFromBoundary` / `…FromCurvature` off so the generator's
> `setSize(points, resolution)` cannot override the field.
> `20260829T123331Z_GEO-22-step2-sizefield.log`, `Status 0`, **33 s**, `-n 1`.
>
> | `resolution` | example, leg C | example, leg D | gate, leg C | gate, leg D |
> | --- | --- | --- | --- | --- |
> | 0.00800 | OK 21 830 | OK 19 823 | OK 8 262 | OK 10 196 |
> | 0.00825 | OK 18 745 | OK 18 807 | OK 8 004 | OK 9 596 |
> | 0.00850 | OK 17 644 | OK 17 563 | OK 7 755 | OK 9 248 |
> | 0.00875 | **FAIL** | OK 16 655 | **FAIL** | OK 8 892 |
> | 0.00900 | OK 14 709 | OK 15 909 | **FAIL** | OK 8 579 |
> | 0.00925 | **FAIL** | OK 15 464 | OK 6 894 | OK 8 144 |
> | 0.00950 | OK 17 683 | OK 14 980 | OK 6 768 | OK 7 918 |
> | 0.00975 | **FAIL** | OK 14 331 | OK 12 200 | OK 7 757 |
> | 0.01000 | **FAIL** | OK 13 837 | **FAIL** | OK 7 407 |
>
> **The two numbers the ruling asked for.** Leg D reads **18/18 OK** and
> **0/18 rungs with a `triangles are equivalent` line**; the whole log contains
> **zero** occurrences of that string and **zero** of `MeshAdapt`, against
> **18** occurrences (exactly one per rung) in both leg C runs. So
> Frontal-Delaunay now completes the wire surface unaided at every size in the
> band, and the seven step-1 failures are gone with it — including
> `h = 0.01000`, the rung that opened this entry on 2026-08-25.
>
> **Negative control, executed in its own process and reproduced
> bit-identically:** leg C re-run
> (`20260829T123413Z_GEO-22-step2-legC-control.log`, `Status 0`, **20 s**)
> returns step 1's table cell for cell — same OK/FAIL in all 18 cells, same
> cell counts to the digit (21 830 at the example's 0.008; the gate's 6 768 at
> 0.00950 and the 12 200 at the coarser 0.00975), same two `NON-MONOTONE`
> verdicts. Leg C and leg D are separate command-line modes and were run as
> two commands, so the control saw exactly the process history step 1 gave it
> (`GEO-23` finding F). The change is therefore the size field's and not the
> day's or the process's.
>
> **A third reading, free:** leg D's cell count is **monotone decreasing in
> `h`** on both geometries (example 19 823 → 13 837, gate 10 196 → 7 407),
> where leg C's jumps around (gate 6 768 at 0.00950 → 12 200 at the coarser
> 0.00975, 1.80×). So the discontinuous response to `resolution` recorded as
> step 1 finding 3 is *also* the wire surface, not the volume mesher.
>
> **Still not licensed.** The field is in the probe only; `src/` is untouched
> and no record moved. Landing it in `straight_wire_domain` would move
> `mag:1`'s **21 830** (leg D reads 19 823 at the same `h`) and the three
> straight-wire ladder records, so the re-record call is the 2026-08-30 weekly
> review's, per the 08-28 10:30 ruling. **Retire-when, unchanged in substance:**
> this entry retires when either a review rules the `GEO-23` step-2a wrap
> sufficient, or the size field lands in `straight_wire_domain` under that
> licence with the moved records re-recorded. **Owning chunk:** `GEO-22`, now
> **✅** — both clauses of its restated done-when are met (wrap + gate landed
> 08-28; this table recorded 08-29), and the size-field decision is the weekly
> review's, not a `GEO-22` step.
>
> **NOW GATED, ENTRY STILL OPEN (`GEO-22` step 2c, 2026-08-29 16:30 slot).**
> The two numbers above are no longer prose only:
> `tests/mesh/test_straight_wire_size_field_probe.py` builds the example
> geometry at `h = 0.008` twice in one process — under the probe's imported
> `_SizeFieldPatch` and without — and asserts **19 823 cells with 0
> fallbacks** patched against **21 830 cells with ≥ 1 fallback** unpatched,
> each count within ±1%. Both readings reproduce **exactly** at `-n 1`
> (`20260829T213132Z_GEO-22-step2c-n1.log`, 8 s) and `-n 2`
> (`…213148Z_…-n2.log`, 7 s). This changes nothing about the licence: the
> field is still probe-only, `src/` is still untouched, and `mag:1`'s 21 830
> is still the shipped record — the gate now *pins* both sides of the
> comparison, so a future landing of the field in `straight_wire_domain`
> must move this module's references deliberately rather than silently.

### ✅ RETIRED 2026-08-25 by `OPS-24` — `core/cavity.py` was **never migrated to dolfinx 0.11**: `assemble_matrix(..., diagonal=)` no longer exists, so the whole `TH-9` cavity + resonance-guard family was **non-executing on `main`** (`EX-30` leg (th), 2026-08-24)

**Retired.** The keyword was renamed `diagonal=` → `diag=` in 0.11 with unchanged
semantics (0.11 docstring: "Rows/columns that are constrained by a Dirichlet
boundary condition are zeroed, with the diagonal to set to `diag`"), verified by
introspecting the installed `dolfinx.fem.petsc.assemble_matrix` signature rather
than assumed. Both sites migrated; all four tests green at `-n 2` complex, and
every recorded figure reproduces the pre-0.11 record **to the printed digit**:
per-mode errors 0.0123 / 0.0153 / 0.0201 / **0.0436%** worst-mode against the
closed form on 720 cells / 5330 dofs, refinement 0.0436% → 0.0102% at fitted
rate 3.85, null cluster 8/8 with max |λ| = 5.560e-14, guard 137.554 (near-
resonant) vs 21.951 (clear) against the 50.0 threshold. `13 passed in 29.71s`
(with `tests/environment`), Status 0, 31 s harness —
`20260825T020157Z_OPS-24-green-quoted.log`; the red baseline was reproduced
in-slot first (`4 failed, 9 passed in 1.83s`,
`20260825T020052Z_OPS-24-red-baseline.log`). No band, tolerance or recorded
eigenfrequency was touched. Original entry below, for the audit trail.

**Test ids — four, all red on `main` with no local change:**

```
tests/validation/test_cavity_resonances.py::test_pec_cavity_resonances_match_closed_form
tests/validation/test_cavity_resonances.py::test_pec_cavity_resonances_improve_under_refinement
tests/validation/test_cavity_resonances.py::test_n1curl_gradient_modes_form_a_clean_zero_cluster
tests/validation/test_resonance_guard.py::test_energy_continuity_guard_fires_near_a_cavity_mode
```

and the two examples on the same code path, `th:2`
(`examples/time_harmonic/02_pec_cavity_resonances.py`) and `th:5`
(`examples/time_harmonic/05_resonance_guard_sweep.py`), which crash outright.

**Literal symptom**, identical at every site:

```
TypeError: assemble_matrix() got an unexpected keyword argument 'diagonal'
```

from `src/fem_em_solver/core/cavity.py:129`
(`A = assemble_matrix(stiffness, bcs=[bc], diagonal=bc_diagonal)`) and its
sibling `cavity.py:131` (`B = assemble_matrix(mass, bcs=[bc], diagonal=1.0)`),
reached through `_cavity_forms` from `cavity.py:229` (`solve_pec_cavity_modes`)
and `cavity.py:324`.

Gate probe, complex build, `FEM_EM_REQUIRE_COMPLEX=1`, `tests/environment`
first, `-n 2`: **`4 failed, 9 passed in 2.11s`**
(`docs/testing/logs/20260824T213908Z_EX-30-th-cavity-gate-probe.log`, Status 1,
4 s harness). Example crashes:
`20260824T213123Z_EX-30-th-run-1to4.log` (`th:2`) and
`20260824T213228Z_EX-30-th-run-5to8.log` (`th:5`).

**Verified at:** `main` @ `7529fa4` (2026-08-24), 0.11.0.post0 image.

**Cause — diagnosed, and it is ours.** `dolfinx.fem.petsc.assemble_matrix`
dropped the `diagonal=` keyword between 0.7.2 and 0.11. `OPS-18` step 2 migrated
the codebase to the 0.11 API but missed this module, and the miss was invisible
because **nothing scheduled runs these two test modules** — they are not in any
chunk's verification path, so the 0.11 merge (2026-08-23) landed on a green-
looking tree with a dead subsystem inside it. **All 9 `tests/environment` tests
pass in the same probe run**, so this is not an environment or complex-mode
regression.

**What this costs, stated plainly.** `TH-9`'s closed-form eigenfrequency
comparison and the energy-continuity resonance guard have produced **no number
since 2026-08-23**. Any reading of either as "gated" is unsupported until this
is fixed — per §9's standing rule, a status without a log reads "unknown".

**Retire when:** the two calls are migrated to the 0.11 signature and all four
tests above run green, with the closed-form eigenfrequency comparison quoted.
Scoped by `EX-30` leg (th) as an `OPS-18` follow-on for a review to queue; it
also unblocks 2 of the 6 `time_harmonic` artifacts that leg could not refresh.
**Commissioned as `OPS-24`** (2026-08-24 18:00 review, §9 item 2); §2.1 now
carries the non-executing caveat on the cavity figure until this retires.

### ✅ RETIRED 2026-08-25 by `OPS-25` — `th:7` calls `Function.interpolate(cells=)`, removed in 0.11 — the **only** such site in the repo, so the example has diverged from the gate it claims to import (`EX-30` leg (th), 2026-08-24)

**Retired by hoist, per the ruling — the divergence is gone, not just the
`TypeError`.** The five lines the example had re-derived (CG2 vector space →
`Function` → sphere-cell index array → restricted `interpolate`) now live once,
in the gate module, as `series_interior_function(series, msh, cell_tags)`
(`test_lossy_sphere_fullwave.py:367`); the gate's `_power_rung` and the example
both call it, so there is no second copy left to rot. The example's private
`cells=` line is deleted, not repaired; the surviving call site is the gate's
already-migrated `cells0=`.

**Evidence the hoist is behaviour-preserving.** The gate's own power figures
reproduce the pre-refactor green log
(`20260822T123746Z_OPS-18-step3-th10-rerun.log`) **bit-identically to all ten
printed digits** — `P_series(meshed)` = 1.048951142e-07 W at the coarse rung and
1.066439173e-07 W at the fine rung, errors 8.387% / 3.629%, quadrature-16
recheck 1.24e-16 — and that is the *only* quantity the moved code produces.
`13 passed in 25.28s`, Status 0, 27 s harness, `-n 2` complex with
`tests/environment` first, covering `test_lossy_sphere_fullwave.py` **and**
`test_lossy_sphere_degree2.py` (`20260825T033221Z_OPS-25-gate-green.log`).

**`th:7` green end-to-end**, `./scripts/run_examples.sh -e th:7 -n 2 -t 300`,
Status 0, 14 s (`20260825T033152Z_OPS-25-th7-green.log`), asserting both
element-order records against their own 1% band: degree 1 relL2 8.1541%
(drift 4.00e-06) / power 8.3869% (1.18e-05), degree 2 relL2 0.1405%
(5.50e-05) / power 0.0058% (1.48e-03). The red was reproduced in-slot first
(`20260825T033114Z_OPS-25-red-baseline.log`, Status 1, `TypeError` at line 198).
No record, band or assertion moved anywhere — this was a hoist. Original entry
below, for the audit trail.

**Test id:** none — no test asserts this. The failing artifact is the example
`examples/time_harmonic/07_element_order_lossy_sphere.py` (`th:7`), which cannot
run at all. `tests/validation/test_lossy_sphere_degree2.py`, the gate it builds
on, does **not** exercise this call and is not implicated by this entry.

**Literal symptom**
(`docs/testing/logs/20260824T213804Z_EX-30-th-run-7to8.log`, Status 1, 2 s,
`-n 2`, complex build):

```
File "/workspace/examples/time_harmonic/07_element_order_lossy_sphere.py", line 198, in _row_and_fields
  e_series_fn.interpolate(_series_interior_interpolant(series), cells=sphere_cells)
TypeError: Function.interpolate() got an unexpected keyword argument 'cells'
```

**Verified at:** `main` @ `7529fa4` (2026-08-24), 0.11.0.post0 image.

**Cause — diagnosed.** 0.11 renamed the cell-restriction argument of
`Function.interpolate`. A repo-wide grep for `interpolate(...cells=)` across
`src/`, `tests/` and `examples/` returns **exactly one hit**: that line. So this
is not a migration class, it is a single site — and the reason it is a *lone*
site is the interesting part. The example's own banner says its fixture is
"imported wholesale from tests/validation/test_lossy_sphere_degree2.py and the
TH-10 module it builds on", but this interpolation step is the example's own
code, not imported. The gate is green and the example is dead, which is exactly
the drift the `ANS-1` rule (records and machinery live in the gate, examples
import them) exists to prevent.

**Retire when:** the call is migrated *or*, preferably, the interpolation is
hoisted into the gate module and imported — a review should decide which, since
repairing it in place preserves the divergence. Unblocks 2 of the 6
`time_harmonic` artifacts `EX-30` leg (th) could not refresh.
**Ruled and commissioned as `OPS-25`** (2026-08-24 18:00 review, §9 item 3):
hoist and import — in-place repair is rejected because it preserves the
divergence this entry documents.

### ✅ RETIRED 2026-08-25 by `EX-30` leg (th) — `th:6`'s **128 MHz** interior relL2 does not reproduce the `TH-10` record on the 0.11 image (1.76864% vs 1.826%, 3.14% drift) while **64 MHz reproduces to 4.04e-05 on the same mesh** (`EX-30` leg (th), 2026-08-24)

**Retired by the licensed version-tagged re-record — the diagnosis below was
right and needed no new measurement.** The example's two restated 128 MHz
constants now carry the 0.11 digits from `TH-10`'s own gate re-run
(`20260822T123746Z_OPS-18-step3-th10-rerun.log`): `RECORD_INTERIOR_L2[128 MHz]`
0.01826 → **0.01769** and `RECORD_SEPARATION[128 MHz]` 57.31 → **59.16**, with
the 0.7.2 digits and their 55 251-cell mesh kept beside them in-comment. The
1% reproduction band did **not** move, and the 64 MHz constants were not
touched.

`th:6` runs green end-to-end, `-n 2`, complex, `20260825T050232Z_EX-30-th-run-5to6.log`
(Status 0, 55 s for `th:5` + `th:6`): 128 MHz fine rung **relL2 1.769% against
1.769% (drift 2.02e-04)**, **separation 59.16× against 59.16× (drift
5.45e-05)** at 55 241 cells, and 64 MHz unmoved at 3.643% / 18.67× (drifts
4.04e-05 / 2.96e-04) — the licensed pair now reproduces three decades inside
its own band, and the frequency asymmetry that made this entry is gone because
it was never physics. Original entry below, for the audit trail.

**Test id:** the assertion is in the example,
`examples/time_harmonic/06_larmor_lossy_sphere.py:186` (`_check_record`), which
gates against a 1% reproduction band. The owning gate module for the record
itself is `TH-10`'s; **it has not been re-run on 0.11**, which is the open
question below.

**Literal symptom**
(`docs/testing/logs/20260824T213236Z_EX-30-th-run-6to8.log`):

```
[64 MHz]  h_sphere = 0.00833 ( 17667 cells): relL2 = 3.643%, separation 18.67x
[64 MHz]  vs the TH-10 record: relL2 3.643% against 3.643% (drift 4.04e-05),
          separation 18.67x against 18.68x (drift 2.96e-04) — band 1%
[128 MHz] h_sphere = 0.00833 ( 17667 cells): relL2 = 3.302%, separation 31.75x
[128 MHz] h_sphere = 0.00556 ( 55241 cells): relL2 = 1.769%, separation 59.16x
AssertionError: [128 MHz] fine-rung interior relL2: this run measured 0.0176864
against the `TH-10` record 0.01826, a drift of 3.14% outside the 1% reproduction band
```

**Verified at:** `main` @ `7529fa4` (2026-08-24), 0.11.0.post0 image, `-n 2`,
complex build.

**Cause — NOT diagnosed. One thing is excluded by measurement: the mesh.** Every
0.11 record motion recorded in this file so far has been a moved mesh underneath
a converged solve. This one is not: the 128 MHz fine rung meshes to **55 241
cells, which is `TH-10`'s own re-recorded count**, and the coarse rung to 17 667,
also on record. Same code, same mesh, same rung — 3.14% different answer at
128 MHz and 4.04e-05 at 64 MHz. The frequency asymmetry is the whole finding;
128 MHz is the harder-conditioned rung (`|m|k₀a` = 1.374 vs 0.850, series N = 8
vs 7, last-term bound 7.2e-16 vs 8.1e-16, so the series itself is not the
suspect).

**Why this matters more than an example red.** 3.643% / 1.826% is the pair
CLAUDE.md and PROJECT_PLAN §2 quote as `TH-10`'s close — the Larmor-frequency
validation gate. The 64 MHz half reproduces; **the 128 MHz half is currently
unreproduced on the image `main` boots.**

**Nothing was re-recorded and no band was moved** — the assertion's own message
and §9 item 4's negative-result clause both forbid it, and the drift is 3× the
band on a quantity §2 depends on.

**The next measurement, and it is one standard-tier command:** re-run `TH-10`'s
own gate module on 0.11. If the gate reproduces 1.826%, this is an example/gate
path divergence (cf. the `th:7` entry above, same family, same slot). If the gate
also measures ≈1.7686%, it is a real 0.11 motion in a §2 figure and the §2
sentence needs revising. Until one of those is done, this is undiagnosed.

**DIAGNOSED 2026-08-24, 18:00 review — from documentation; the measurement
already existed.** `TH-10`'s gate module *was* re-run on 0.11, on 2026-08-22,
by `OPS-18` step 3 attempt 1: green log
`20260822T123746Z_OPS-18-step3-th10-rerun.log` (11 passed, exit 0) prints the
128 MHz fine rung at **relL2 1.769%, separation 59.16×, 55 241 cells** — the
re-record 1.826% → 1.769% is explicit in §7 `OPS-18` ("*with its mesh*,
55 251 → 55 241"). `th:6` measured 1.76864% / 59.16× on 55 241 cells: the
gate's own 0.11 digits to 2e-4, against its never-updated restated constants
(`RECORD_INTERIOR_L2[128 MHz]` 0.01826, `RECORD_SEPARATION[128 MHz]` 57.31).
So this is the **example/gate-divergence branch** — same class as the `th:7`
entry — not a physics motion; the "mesh did not move" observation above was
correct and is exactly why: the mesh had already moved *at the re-record*,
and the record moved with it while the example's copy did not.

**Retire when:** `th:6`'s two 128 MHz constants are re-recorded
version-tagged from the `20260822T123746Z` log (0.7.2 digits kept beside)
and the example runs green — licensed to `EX-30` leg (th)'s re-run
(§9 item 4, 2026-08-24 18:00 review). §2 and CLAUDE.md now quote the
version-tagged pair.

**Also on that log — the `Status: 124` is a teardown hang, not a compute
overrun.** After the assertion fires (≈40 s of real work), MPI deadlocks in
`mpi4py.MPI.commlock_free_cb` during interpreter shutdown (`SystemError: …
returned a result with an exception set`); the container-side `timeout -k 30 300`
fired and PETSc reported signal 15. The `-k 30` worked as designed —
`docker compose ps` read Up and `pgrep -c python3` read **0** immediately after,
no wedge and no force-recreate. Worth knowing: an assertion failure in this
example does not exit cleanly under `mpiexec`, so budget for the full timeout
when re-running it red.

### ✅ RETIRED 2026-08-25 (`GEO-19` step C, 07:30 implementer slot) — `birdcage_port_domain(emit_port_sheets=True)` **cannot build any birdcage with more than four legs**: the mid-plane sheet is an axis-aligned rectangle (`GEO-19` attempt 1, 2026-08-23)

**Retired on the run this entry named as its retire-when.** Step B's
local-frame construction landed on `main` 2026-08-25, and step C built the
16-leg fixture from `main`: **307 296 cells / 74.18 s of mesh time**, all 16
sheets at `dx·g` = 1.000000000000 with a C16 area spread of **1.331e-15** and
out-of-plane extents ≤ 1.736e-17 m in each port's own frame, 32 half-boxes at
0.500000000000, `GEO-9` partition 1.000000000000
(`20260825T124357Z_GEO-19-stepC-run1.log`, `-n 2`, 114 s). The
`NotImplementedError` is unreachable: ports P2/P4/… sit at 22.5° + 45k and
build. What step C found instead is a *band-domain* question on gate (ii)'s
equality half, filed as its own entry below — the capability limit this entry
described is gone.

**Test id:** no test asserts this on `main` — the module that hits it,
`tests/mesh/test_birdcage_port_scaleup.py`, is parked on
`attempt/GEO-19-20260823T214500Z` rather than landed red.

**Literal symptom**, at `leg_count = 16` with `leg_gap_length = 8e-3`,
`emit_port_sheets = True`:

```
NotImplementedError: emit_port_sheets builds axis-aligned rectangles, so
every leg must sit on a coordinate axis; port P2 is at 22.500 degrees
```

raised from `src/fem_em_solver/io/mesh.py:3189`, `-n 2`, 1.4 s,
`docs/testing/logs/20260823T213546Z_GEO-19-step1.log`.

**Verified at:** `main` @ the `GEO-19` blocker-A commit (2026-08-23), on the
0.11.0.post0 image.

**Cause — diagnosed, not a mystery.** The sheet is entered into the OCC
fragment as an axis-aligned dim-2 tool, and the half-assignment that follows
tests a *single Cartesian centroid coordinate* against the plane's offset
(`mesh.py:3270-3278`). Both steps assume the leg's radial direction is `x̂`
or `ŷ`, which holds only for `leg_count <= 4`. Any other count puts legs at
intermediate azimuths and the construction has no rectangle to place. This is
a genuine capability limit, not a regression: nothing above four legs had ever
been sheeted, so it had never been reached.

**Not to be confused with** the `100+i`/`110+i` tag collision in the same code
path, which was the *outer* guard (`leg_count <= 9`) and **is fixed** — the
upper base is now `200+i`, verified inert against `GEO-18` step 1/2 at
`3 passed` / 116 368 cells / C4 spread 6.050e-16, digit-identical to the
pre-change run (`20260823T213647Z_GEO-19-tagfix-regression.log`).

**Retire when:** the sheet is built in the leg's local `(r̂, ẑ)` frame with the
half test taken along the leg's radial normal, and `GEO-19`'s gates (i)–(v)
run at 16 legs. **Fix scoped 2026-08-23 18:00 review as `GEO-19` step B
(§9 item 3, rotated construction + 4-leg invariance control); the 16-leg
gates are step C (§9 item 5). This entry retires with step B's commit** —
the gates run is a separate deliverable and does not hold it open.

> **Step B attempt 1, 2026-08-24 03:30Z implementer slot — the rewrite exists,
> is green, and is parked rather than landed.** `attempt/GEO-19-stepB-`
> `20260824T034500Z` (`12737a8`) carries it: box and sheet both built at
> azimuth 0 and taken to the leg's azimuth by one snapped rotation about `ẑ`,
> half-assignment by signed projection on the plane's own normal, the raise
> deleted. The `GEO-18` step-1 + step-2 invariance control is `3 passed` /
> Status 0 **twice in-slot** (90.08 s / 88.97 s,
> `20260824T033811Z_GEO-19-stepB-snapped-run1.log`,
> `20260824T033956Z_…-run2.log`), reproducing the record's terminal ratios
> (0.988616 × 4) and C4 sheet spread (**6.050e-16**, exactly) with every
> analytic identity exact and the CAD digit-identical.
> **What holds it back is not the geometry but its consumers:** the mesh cell
> count moves (116 368 → **116 085** sheeted, 114 855 → **114 655** gapped),
> and on that moved fixture three `PORT-9` birdcage assertions go red
> (`20260824T034214Z_GEO-19-stepB-port9-regression.log`, `3 failed, 16 passed`
> / 124.68 s): leg (c)'s driven current deviates 1.376e-03 from record, leg
> (d0)'s `Z_11` 1.840e-02 against a 1e-9 print band, and leg (c)'s
> class-degeneracy gate **flips** — the opposite port sits 0.0321% from the
> adjacent pair's mean against the pair's own 0.0407% spread, i.e. inside it.
> Re-recording those digits is §9 item 4's licence, and the flipped gate is a
> ruling, so this entry stays open until a review disposes of it.

> **Ruled 2026-08-24 03:00 review — ruling (4\*), full text in §9.** The
> rewrite is adjudicated correct; the cell-count digit-for-digit expectation
> was unsatisfiable (the old construction sat ≈5 ulps off exact positions and
> gmsh tie-breaking amplifies that to ≈1e-3 in cell count; CAD and every
> analytic identity reproduce digit for digit). Sequencing: `PORT-9` leg
> (d3b) re-records the birdcage class on the **unmoved** mesh first (§9
> item 1), then step B lands with a mesh-tagged re-record (§9 item 2), whose
> licence includes the pre-registered degeneracy-gate disposition. **This
> entry retires with §9 item 2's commit.**

Blocked `GEO-19` until step C ran; never blocked `GEO-20` (ring-gap ports
sit at different azimuths and are scoped to their own local frame from the
start). Related and separate: the layout clearance floor independently caps
this geometry at `N <= 25` legs — see the `GEO-19` §7 entry.

### ✅ RETIRED 2026-08-25 (`GEO-19` step C, 12:00 implementer slot) — `GEO-18` step 1's **1e-5 terminal-equality band is a C4 band, not a C_N one**: at 16 legs the terminal disks spread **8.434e-04**, two decades wide of it, while every ratio stays inside the closed form's [0.95, 1.0] (`GEO-19` step C, 2026-08-25)

**Test id:** `tests/mesh/test_birdcage_port_scaleup.py::`
`test_sixteen_leg_identity_family_and_cost_rung`, parked on
`attempt/GEO-19-stepC-20260825T125000Z` (`e7a3926`) rather than landed red.

**Literal symptom**, at `leg_count = 16`, `-n 2`
(`docs/testing/logs/20260825T124357Z_GEO-19-stepC-run1.log`, 114 s):

```
AssertionError: 16 legs: the 16 terminal areas differ by 8.434e-04 relative
against the pre-stated 1e-05; the ports are not the same disk under a rotation
assert np.float64(0.0008433598402112139) < 1e-05
```

**Verified at:** `main` @ `d74c7a5`, 0.11.0.post0 image, step B's local-frame
construction.

**Cause — diagnosed, and it is arithmetic rather than geometry.** The 16
meshed terminal ratios take **three** values, and they sort by azimuth:

| azimuth | ratio to `2·π·r_leg²` | ports |
|---|---|---|
| 0/45/90/…/315° | 0.988615667 … 0.988615857 | P1 P3 P5 P7 P9 P11 P13 P15 |
| 22.5/157.5/202.5/337.5° | 0.989367491 … 0.989367549 | P2 P8 P10 P16 |
| 67.5/112.5/247.5/292.5° | 0.989449699 … 0.989449760 | P4 P6 P12 P14 |

Within each class the spread is ≤ 2e-7 — *tighter* than the 1e-5 band. The
band was measured at C4, where the four ports are related by exact 90°
coordinate permutations, so their inscribed triangulations are the same
arithmetic and agree to 3.2e-08 (this run's control reproduces that:
**3.184e-08**). At 22.5° there is no coordinate permutation; the disk's
inscribed polygon lands on different nodes, and an inscribed triangulation
under-reads a circle by ≈1.1% to begin with. The observed 8.4e-04 spread is
**thirteen times smaller than that under-read** — it is the discretization
error's own azimuthal variation, not a broken port.

**Not a widening candidate in-slot.** The `GEO-19` §7 entry's negative-result
clause is explicit ("a gate red at 16 legs is a generator finding at the new
count — known-issues + §7 annotation, stop; never widen a `GEO-18` band"), so
the band is untouched and the module is parked. Everything else step C
pre-stated is green and is recorded in the §7 entry.

**Retire when:** a review rules on which of these the equality gate is asserting
— exact C_N symmetry of the *construction* (in which case the reading is
per-azimuth-class and the gate becomes an intra-class one at ≈1e-6, tighter
than today's) or agreement of the *discretization* across azimuths (in which
case the band is h-dependent and wants a refinement rung, not a constant). Both
readings are consistent with the measurement; choosing between them is a
ruling, not an implementer judgement.

**RULED 2026-08-25, 10:30 review: construction symmetry.** The ≤ 2e-7
intra-class tightness says the construction is exactly C_N-covariant, and the
inter-class 8.4e-04 is the inscribed triangulation's azimuthal variation
(13× below its own ≈1.1% under-read scale) — gating it at 1e-5 gates the
mesh, not the generator. The reading becomes per-azimuth-class: intra-class
equality asserted at **1e-6** (tighter than today; basis ≤ 2e-7), inter-class
spread asserted under a coarse **5e-3** discretization ceiling (basis
8.434e-04; half the under-read scale) so a broken port cannot hide, the
absolute [0.95, 1.0] band and every C4 module's 1e-5 unmoved. Full landing
instructions in the `GEO-19` §7 entry. **This entry now retires with the
commit that lands the ruled module green from `main`** (§9 item 1,
2026-08-25 10:30 queue).

**RETIRED 2026-08-25, 12:00 slot — the ruled module is green from `main`,
and the class table above reproduces port for port.**
`20260825T170316Z_GEO-19-stepC-ruled.log` (`2 passed` / **117 s**, Status 0)
and the record run `20260825T170523Z_GEO-19-stepC-ruled-record.log`
(`2 passed` / **115 s**, Status 0, `-s`). The partition is taken from the
mesh's own coordinate mirrors (`_azimuth_class`), not from the measured
areas, and it lands on exactly the three classes tabulated above:

| class | ports | meshed/analytic | intra-class spread (band 1e-6) |
|---|---|---|---|
| aligned (0/45/…/315°) | 8 | 0.988615772 | **1.923e-07** |
| 22.5/157.5/202.5/337.5° | 4 | 0.989367514 | **5.849e-08** |
| 67.5/112.5/247.5/292.5° | 4 | 0.989449735 | **6.144e-08** |

Inter-class spread **8.431e-04** against the 5e-3 ceiling. The four-leg
control returns **one** class at **3.184e-08** — the back-compat identity the
ruling asked for, i.e. the reading reduces to the old flat gate exactly — with
116 085 cells (delta **0**) and C4 sheet spread 6.050e-16. No band outside
`tests/mesh/test_birdcage_port_scaleup.py` moved; the C4 modules keep their
1e-5 and `TERMINAL_AREA_BAND` keeps [0.95, 1.0].

### ✅ RETIRED 2026-08-24 (`PORT-9` leg (d3c), 12:00 implementer slot) — two birdcage **reproduction controls** were red on `main`: the 0.11 image meshes the gapped birdcage at **116 368** cells, the records were taken at **116 416** (`PORT-9` leg (d3b), 2026-08-24)

**Test ids** (both `@complex_only`, both failed on `main` with no local change):

```
tests/validation/test_port_birdcage_termination_probe.py::test_the_open_control_reproduces_leg_c_before_the_knob_turns
tests/validation/test_port_birdcage_four_port.py::test_the_driven_column_reproduces_leg_d0
```

**Literal symptom**, `-n 2`, complex build, `2 failed, 17 passed`:

```
AssertionError: the open solve's driven current +9.992781266e-07+3.346865998e-09j A
deviates 6.829e-06 from leg (c)'s recorded +9.992734880e-07+3.351870842e-09j A
AssertionError: Z_11 = +2.172952668e+01+7.461413742e+00j Ohm deviates 1.449e-04 from
leg (d0)'s recorded +2.173224483e+01+7.459491479e+00j Ohm against the 1e-09
print-precision band
```

`docs/testing/logs/20260824T093133Z_PORT-9-step3d3b-run1.log` and
`…T093526Z_…-run2.log` — both runs identical in every Z digit.

**Verified at:** `main` @ `082e30f` (2026-08-24), 0.11.0.post0 image.

**Cause — measured, not a mystery: the mesh moved, and it is not the route.**
All three modules printed `116368 cells (record 116416, ratio 0.999588)`. The
`PORT-9` leg (c)/(d0)/(d) records were taken 2026-08-22/23 **before** the
`OPS-18` step 3b merge put main on dolfinx/gmsh 0.11, and 0f8ea96 already
measured 116 368 on 2026-08-23 *both before and after* its own tag change, so
the tag encoding is excluded. This is the same 1e-4 record motion the retired
0.11 entry below recorded for the two-torus family, recurring for the birdcage
family, which nothing re-gated on 0.11. The `PORT-9` leg (d3) route change
cannot be the cause: it touched only `_assemble_sparameter_matrix`, and these
two modules never call the sweep's S assembly at all. **Every physics identity
still held** on the moved mesh — sheet area `5.930614898e-05 m²`, `h` exactly
`8.000000000e-03 m`, out-of-plane `8.882e-19 m`, all four sheets identical, C4
class spreads 0.0617 / 0.0359 / 0.0237% and the (d0) discrimination margin
253.2002× against a 10× floor.

**Retired by:** ruling (5\*) (2026-08-24 10:30 review) granted the image-caused
re-record; `PORT-9` leg (d3c) executed it in the 12:00 slot. The two controls
are now recorded image-tagged at 116 368 — leg (c) `I_1` =
+9.992781266e-07 + 3.346865998e-09j A, leg (d0) `Z_11` =
+2.172952668e+01 + 7.461413742e+00j Ω, both beside the pre-0.11 digits in the
modules' own comments — and both runs of the re-record read **19 passed** with
every edited constant reproducing to ≤ 2.4e-10 against its 1e-9 band
(`20260824T170332Z_PORT-9-step3d3c-run1.log`,
`20260824T170544Z_PORT-9-step3d3c-run2.log`). No band was moved.

### ✅ RETIRED 2026-08-23 (`OPS-18` step 3b merge) — the two-torus port fixtures **SIGABRT'd in `gmsh.model.mesh.generate` only in the 0.11 image**: numpy 2 renders `!r` of a numpy scalar as `np.float64(…)` inside a gmsh `MathEval` string (`OPS-18` step 3 attempt 2, 2026-08-22; entered by the 10:30 review)

> **RETIRED 2026-08-23, 15:00 implementer slot.** The entry's own retirement
> condition is met: the fixing commit `445a3ea` (`float()` coercion at the
> four `MathEval` sites, attempt 3) reached `main` with the `OPS-18` step 3b
> merge, and `main` now boots the 0.11 image the defect was scoped to. The
> two-torus fixtures the SIGABRT killed are green on that image at
> `19 passed` / exit 0 twice (step 3a attempt 8,
> `20260823T170403Z_OPS-18-step3a-leg1-run1.log`,
> `20260823T170821Z_…-run2.log`). Original entry follows.
>
> **Where this fires.** `tests/validation/test_port_package_sparameters.py`
> and `tests/validation/test_port_lumped_two_torus.py` (every test that
> builds `two_torus_domain`), on the `attempt/OPS-18` worksite only. `main`
> boots 0.7.2 with numpy 1.x and is unaffected.
>
> **Literal symptom** (`20260822T140912Z_OPS-18-step3-port1-rerun.log`,
> `Status: 134` at 12 s, `-n 2`, complex build):
>
> ```
> Error   : Error [mathex::parseatom()]: invalid token on expression
> terminate called after throwing an instance of 'std::runtime_error'
> ```
>
> **Cause, measured by two probes.** gmsh in the 0.11 image is
> 4.15.2-git-657c8e9 and parses the exact gap-arc expression when its numbers
> are plain floats (`20260822T141005Z_…-gmsh-mathex-probe.log`); the image's
> numpy 2.4.6 renders a numpy scalar's `repr` as `np.float64(0.00591…)`
> (`20260822T141027Z_…-numpy-repr-probe.log`). `two_torus_domain` builds its
> `MathEval` size field by f-string `!r` interpolation (`io/mesh.py:1550`–
> `1557`) of `arc_half_y`, a numpy scalar, so the string literally contains
> `np.float64(`. Image debt of ours that numpy 1.x masked — not an upstream
> regression, not a gated number moving; no band, assertion or record touched.
>
> **Fix owner:** §9 item 3a — `float()` coercion at the sites plus a `src/`
> sweep for `!r` in any string handed to gmsh, PETSc options or a shell (the
> birdcage fixtures are predicted to carry siblings). This entry leaves with
> that commit.
>
> **FIXED on the branch 2026-08-22, 12:00 implementer slot (`445a3ea`,
> `attempt/OPS-18`) — this entry retires when that branch merges (§9 item
> 3b (iii)).** `arc_half_y`, `major_radius` and `z_c` are coerced to plain
> `float` before the f-strings. Negative control discharged against the red
> baseline quoted above: the same batch now meshes and runs to a footer,
> `17 passed / 2 failed` in 260.93 s
> (`20260822T170346Z_OPS-18-step3-port1-coerced.log`) — the SIGABRT is gone
> and the two remaining failures are the moved reproduction records in the
> next entry, a different defect. **Sweep result, measured:** `src/` carries
> **53** `!r` interpolations; the **4** two-torus `MathEval` sites are the
> only ones handed to a foreign parser (the other 49 are Python exception
> messages, which is what `!r` is for), and `MathEval` has exactly **one**
> call site in all of `src/` (`io/mesh.py:1568`). The prediction that the
> birdcage fixtures carry siblings is **wrong, by measurement** — they build
> no `MathEval` field. The class is closed at one instance.

### ✅ RETIRED 2026-08-23 (`OPS-18` step 3b merge) — two two-torus **reproduction records** moved in the 0.11 image at 1e-4 while every physics identity held (`OPS-18` step 3 attempt 3, 2026-08-22)

> **RETIRED 2026-08-23, 15:00 implementer slot.** The entry's stated closing
> condition — "this entry closes with step 3b's merge commit, which is what
> carries these branch-only writes to `main`" — is met: `5df1e39` (step 3a
> attempt 8) is on `main` with this merge, and the four re-records plus the
> two the same assertion loop unmasked are all written version-tagged, no
> band moved. See the attempt-8 resolution block at the end of this entry
> for the digits. Original entry follows.
>
> **Where this fires.** `tests/validation/test_port_package_sparameters.py::test_sanity_report_reproduces_the_gated_metrics_on_the_field_route`
> and `tests/validation/test_port_lumped_two_torus.py::test_step_1_measurements_reproduce`,
> on the `attempt/OPS-18` worksite only. `main` boots 0.7.2 and is unaffected.
> These are the two failures left in the log that discharged the SIGABRT
> above; they are a **different defect** from it, and they only became
> visible once the mesh built again.
>
> **Literal symptom** (`20260822T170346Z_OPS-18-step3-port1-coerced.log`,
> `2 failed, 17 passed` in 260.93 s, `Status: 1`, `-n 2`, complex build,
> both rank footers identical):
>
> ```
> passivity_max_sigma 0.861356895 does not reproduce the step-4 record
>   0.861449 within 1.0e-06   (moved 9.210529e-05)
> gap ratio: 0.894141 against step 1's record 0.894310 — moved by 1.69e-04,
>   above 1e-04
> ```
>
> and, printed but not asserted in the same test's stdout:
> `||S-S^T||/||S|| = 3.112128e-05` against the record `2.5494e-05`
> (band 5.0e-07).
>
> **What did *not* move.** Every gate that scores physics rather than a
> recorded digit string passes, in the same run: reciprocity at
> 3.112128e-05 is 32× inside the `PORT-1` band of **1e-3**; passivity holds
> (σ_max 0.8614 < 1, column power sum 0.7411 < 1); the open-limit reduction
> to the sheet average and the cross-route transverse-average identity both
> PASS. 17 of 19 tests pass.
>
> **Cause, not diagnosed — one hypothesis, and it is testable.** gmsh moved
> 4.11 → 4.15.2 with the image, and the meshes it emits move with it: the
> straight-wire fixture in the entry below went 145 900 → **147 235** cells
> on the same requested `h`. A two-torus mesh perturbed at that level moves
> a solved S-matrix at 1e-4, which is exactly the size of both misses and
> ≈7× the tightest band (1e-6). The alternative — that 0.11 changed an
> assembly or interpolation the port route depends on — is not excluded by
> anything measured here. Deciding between them needs the two-torus cell
> count printed on both images; nobody has printed it.
>
> **Update, 2026-08-22 (`OPS-18` step 3 attempt 4): the mesh did move, and
> by 24× more than the records did.** `tests/mesh/probe_two_torus_cell_count.py`
> builds this exact fixture — every argument imported from `_build` in the
> test that owns the record — and prints the reduced global counts:
>
> | | dolfinx 0.7.2 / gmsh 4.11.1 | dolfinx 0.11.0.post0 / gmsh 4.15.2 | Δ |
> |---|---|---|---|
> | cells | 184 919 | 184 176 | **−743, −4.017e-03** |
> | vertices | 31 676 | 31 550 | −126, −3.978e-03 |
> | wire tags `1`/`2` | 9 463 / 9 380 | 9 556 / 9 448 | +0.98% / +0.72% |
> | gap halves `101`/`102`/`111`/`112` | 13 627 / 13 599 / 13 763 / 13 771 | 13 661 / 13 648 / 13 658 / 13 694 | ≤ 0.8% |
>
> Logs `20260822T183313Z_OPS-18-step3-twotorus-cells-072.log` (Status 0, 33 s,
> `main` source on the 0.7.2 image) and
> `20260822T183626Z_OPS-18-step3-twotorus-cells-011.log` (Status 0, 34 s,
> `attempt/OPS-18` source on the 0.11 image), both `mpiexec -n 2`. The
> 0.7.2 run had to be taken on `main`'s source: the branch's `io/mesh.py`
> imports `dolfinx.io.gmsh`, which 0.7.2 does not have, so the two runs
> differ by the step-2 migration as well as the image — a caveat on
> attribution, not on the counts.
>
> **So the mesh hypothesis is *consistent* here**: a 4.0e-03 mesh
> perturbation moving two solved records by 9.2e-05 and 1.7e-04 is a 24-40×
> attenuation, which is what a converged quantity on a perturbed mesh does.
> It is still not a proof — the same slot **refuted** the mesh explanation
> for the straight-wire entry below, where the counts moved by < 0.16% —
> so the two failures should not be assumed to share a cause. **The
> re-record decision remains the review's**; what this update supplies is
> the evidence it asked for.
>
> **No band, assertion or record was touched**, and none should be until
> the cause is decided: re-recording a solved S-matrix to make a version
> bump land is how a version bump hides a physics change. §9 item 3a's
> negative-result clause fired on this and stopped the leg.
>
> **Ruling, 18:00 review 2026-08-22 — re-record licensed, narrowly.** The
> three numbers (`passivity_max_sigma` at 1e-6, the gap ratio at 1e-4, the
> reciprocity record at 5e-7) are reproduction records of a solved field on
> a named mesh, every physics gate in the same run holds, and the mesh
> moved 4.017e-03 — the same grounds on which `TH-10`'s 55 251 → 55 241
> re-record was granted. Conditions: version-tagged (0.7.2 value + 184 919
> cells stay beside the new value + 184 176 cells + `0.11.0.post0 / gmsh
> 4.15.2`), bit-identical across two runs in the slot before the digit
> string is written, no band moves, branch-only until 3b merges. Queued as
> §9 item 4; this entry closes with that commit. Full text in the §7
> `OPS-18` entry.
>
> **Update, `OPS-18` step 3a attempt 7 (2026-08-23, 06:00 slot) — the three
> licensed records are written and green, and writing them surfaced a
> fourth and a fifth of the same kind.** On `attempt/OPS-18` at `44b5600`,
> image `v0.11.0`, complex, `-n 2`, two same-slot runs
> (`20260823T110726Z_OPS-18-step3a-leg1-confirm.log`, `1 failed / 18
> passed` / 234.88 s, and `20260823T112102Z_…-leg1-confirm-rerun.log`,
> same counts / 226.36 s). `passivity_max_sigma` 0.861356895, `‖S−Sᵀ‖/‖S‖`
> 3.11213e-05 (six digits per (b′)) and the gap ratio 0.894141 all
> reproduce as written in both runs; physics green in both (reciprocity
> 2.679e-05 inside 1e-3, σ_max 0.861357 < 1).
>
> The remaining failure is **not** a regression: `test_step_1_measurements_
> reproduce` checks three records in one loop, and the gap ratio was the
> first, so fixing it unmasked the other two.
>
> | record | 0.7.2 | 0.11.0 | move | band | run-to-run move |
> |---|---|---|---|---|---|
> | `STEP1_LUMPED_RATIO_RECORD` | 0.829782 | 0.828893 | 8.89e-04 | 1e-4 | 6.6e-10 (6.6e-06 of band) |
> | `STEP1_CROSS_ROUTE_RECORD` | 0.077095 | 0.077431 | 3.36e-04 | 1e-4 | identical to 6 printed digits |
>
> Both are step-1 reproduction records of the *same* solved field on the
> *same* fixture whose mesh moved 184 919 → 184 176 cells — the class
> ruling (1) licensed — and both satisfy (b′)'s reproduction condition on
> this slot's two runs. But ruling (1) enumerates **three** numbers, and an
> implementer does not extend a review ruling in-slot, so **neither was
> written**. What is owed is one review decision: extend ruling (1) to
> these two, or rule them differently. Nothing else in leg 1 is open.
>
> **Resolved, `OPS-18` step 3a attempt 8 (2026-08-23, 12:00 slot) — the
> 10:30 review's class ruling (1\*) admitted both, they are written, and
> leg 1 is green.** On `attempt/OPS-18` (`main` merged at `070b1b5`),
> image `v0.11.0`, complex, `-n 2`, two same-slot runs both **`19 passed`**
> / exit 0 (`20260823T170403Z_OPS-18-step3a-leg1-run1.log`, 238.64 s;
> `20260823T170821Z_OPS-18-step3a-leg1-run2.log`, 238.73 s), where attempt
> 7 read `1 failed / 18 passed`. `STEP1_LUMPED_RATIO_RECORD` **0.828893**
> and `STEP1_CROSS_ROUTE_RECORD` **0.077431** reproduce identically to six
> printed digits in both runs — (b′) arithmetic: the lumped route's
> `Im Z12` moves 1.029281339 → 1.029281338 Ω (1e-9 absolute, 1e-5 of the
> 1e-4 band once divided through by ωM₁₂), the cross-route print does not
> move at all. Physics green in both: reciprocity 2.679e-05 inside 1e-3,
> σ_max 0.861357 < 1, open-limit and cross-route decomposition PASS, the
> 5% cross-route MISS unchanged and still pre-stated. The fixture meshes
> to 184 176 cells in both runs. **The loop unmasked no further record**,
> and the only other consumer of these constants,
> `test_port_lumped_narrowed_sheet.py`'s `f = 1.0` negative control, is
> green on the new values in the same slot (`12 passed` / 142.72 s / exit
> 0, `20260823T171239Z_OPS-18-step3a-narrowed-sheet.log`; that rung prints
> 7.7431%, gap 0.894141, lumped 0.828893 — the identical digits). No band
> moved. **This entry closes with step 3b's merge commit**, which is what
> carries these branch-only writes to `main`.

### ✅ RETIRED 2026-08-23 (`MAG-18` re-gate, 19:30 implementer slot) — `test_straight_wire_b_field` failed **only in the 0.11 image**: the discretization error moved 12.75% → 15.3848% against a 15% band, on a mesh that grew 145 900 → 147 235 cells (`OPS-18` step 3 attempt 3, 2026-08-22)

> **Closing note (2026-08-23, 19:30 slot).** This entry stayed open after
> step 3b only because the `MAG-18` `E_Ω` ladder was unobserved on 0.11
> (see the scope update below). It has now been observed, on `main`, on
> the image the project boots, twice in the same slot: `7 passed` /
> exit 0 at `-n 2` (`20260824T003059Z_MAG-18-regate-run1.log` and
> `20260824T003650Z_MAG-18-regate-run2.log`, 296 s each), plus a `-n 4`
> record probe (`20260824T003606Z_MAG-18-regate-n4.log`, 32 s). The mesh
> growth this entry is about is *real and unchanged* — the h = 0.0025
> rung still meshes to **147 235** cells on 0.11 — but in the `MAG-18`
> norm it costs the gate nothing: `E_Ω` = 25.2868 → 10.6172 → 6.6458%,
> monotone, fitted rate **1.6854** against 0.7.2's 1.6842 (7e-04 of
> move), cross-width agreement 4.86e-07 ≤ 1e-6, natural-BC ratio 0.3285.
> The 10-point number this entry names is still 15.3848% on 0.11 — it is
> printed, reproduced under the version-keyed control row to 2.7e-06,
> and gates nothing. Nothing was re-recorded and no band moved. Retired,
> not merely superseded.

> **Where this fires.** `tests/validation/test_straight_wire.py::TestStraightWire::test_straight_wire_b_field`,
> real build, on the `attempt/OPS-18` worksite only. `main` boots 0.7.2 and
> is unaffected.
>
> **Scope update, `OPS-18` step 3b (2026-08-23, 15:00 slot) — read this
> first.** `main` now boots the 0.11 image, so the "worksite only" clause
> above is spent. The failure it names, however, **can no longer fire**: the
> 15% band is retired by `MAG-18` (`d494d81`, on `main` since 2026-08-22) —
> `test_straight_wire_b_field` now *prints* the 10-point number as
> "[reported, not gated]" and gates nothing, and the replacement gates are
> `TestStraightWire::test_domain_l2_*` on the sampler-free `E_Ω` and its
> convergence rate. This entry is therefore **superseded, not resolved**: it
> stays open because the thing it is really about — that the 0.11 image's
> gmsh moves the wire mesh (145 900 → 147 235 cells) and the solved error
> with it — has **not** been re-measured in the `MAG-18` norm on 0.11.
> Step 3b did not run the real-mode `MAG` leg (its anchor is the environment
> + mesh-tag family and the two collects), so **the `E_Ω` ladder and its
> rate band on 0.11 are owed and unobserved**; the numbers on record
> (25.3787 / 10.7288 / 6.6708%, rate 1.6842) are 0.7.2 numbers. A review
> should queue that leg before reading any `MAG` figure as re-gated.
>
> **Literal symptom** (`20260822T171401Z_OPS-18-step3-real-mag2.log`,
> `1 failed, 17 passed, 4 skipped` in 272.43 s, `Status: 1`, `-n 2`, both
> rank footers identical):
>
> ```
> Cells: 147235
> Relative L2 error: 15.3848%
> AssertionError: Relative error 15.3848% exceeds 15%
> ```
>
> **What did *not* move.** The rest of the `MAG` family passes on 0.11 in
> the same run: `test_straight_wire_convergence`,
> `test_analytic_bc_improves_on_natural_bc` (the `MAG-13` claim itself),
> both `test_circular_loop` gates, and all 7 `test_mutual_inductance_reference`
> tests. The failure is one number, on the one test whose band was set from
> a measured h-error ladder.
>
> **Cause, measured in part.** The 15% band is not a physics tolerance: the
> comment at `tests/validation/test_straight_wire.py:174`–`186` records it
> as the measured error of *this mesh* — `h=0.0025`, **145.9k cells**,
> **12.75%** on an O(h^1.2) ladder that has no plateau. The 0.11 image's
> gmsh (4.15.2-git-657c8e9 vs 4.11) meshes the same requested `h` to
> **147 235** cells, and this test's error is what a still-converging
> discretization does when its mesh changes — the number was always
> mesh-specific, and 15% was 1.18× a measurement rather than a bound with
> headroom.
>
> **Not diagnosed:** whether 15.3848% is *only* the mesh change. A 1.9%
> cell-count change producing a 21% error change is steeper than the
> recorded ladder (38.8k→145.9k cells for 22.19%→12.75%), so the mesh
> hypothesis is **not sufficient on its own** and something else may be
> contributing. Re-running the recorded ladder's other two rungs on 0.11
> would settle it; that is a slot's work and was not done here.
>
> **Update, 2026-08-22 (`OPS-18` step 3 attempt 4): the ladder was re-run on
> both images, and the mesh explanation is refuted.**
> `tests/validation/probe_straight_wire_ladder.py` drives the recorded rungs
> through the same `_solve_straight_wire` / `_sample_radial` / `ErrorMetrics`
> the gated test uses (imported, not restated), `mpiexec -n 2`, real build:
>
> | h | cells 0.7.2 | error 0.7.2 | cells 0.11 | error 0.11 | Δcells | Δerror |
> |---|---|---|---|---|---|---|
> | 0.0040 | 38 750 | **22.1925%** | 38 740 | **21.8417%** | −0.13% | −1.6% |
> | 0.0025 | **145 884** | **12.7485%** | 147 235 | **15.3848%** | **+0.93%** | **+20.7%** |
> | 0.0018 | 383 248 | **9.2568%** | 383 146 | **4.4605%** | −0.03% | **−51.8%** |
>
> Logs `20260822T184158Z_OPS-18-step3-wire-ladder-072.log` (Status 0, 98 s),
> `20260822T183710Z_OPS-18-step3-wire-ladder-011.log` (Status 0, 105 s) and,
> for the gated rung's 0.7.2 control,
> `20260822T185944Z_OPS-18-step3-wire-h0025-072.log` (Status 0, 27 s); the
> 0.11 `h=0.0025` row is quoted from the two runs above, which agree.
>
> **The 0.7.2 column is a clean control on all three rungs**: the July
> record reproduces to **+0.011%**, **−0.012%** and **−0.035%**, so the
> ladder is not stale and the deltas are the image.
>
> **And the outlier is not rank-dependent.** The gated rung re-run on 0.11
> at `-n 4` is **bit-identical** to the `-n 2` result — 147 235 cells,
> 15.3848% (`20260822T184951Z_OPS-18-step3-wire-h0025-n4.log`, Status 0,
> 28 s). Serial (`-n 1`) is a **sizing** finding, not a result: exit 124 at
> the 400 s ceiling (`20260822T185030Z_…-n1.log`), not retried. So the
> mesh count 147 235 and the 15.3848% are both stable and reproducible
> across rank widths, and partitioning is excluded as the cause.
>
> **What it shows.**
> 1. **Not the mesh.** Both probed rungs mesh to within 0.13% of their
>    recorded cell counts on 0.11, and their errors move by −1.6% and
>    −51.8%. A mesh that does not move cannot explain an error that halves.
> 2. **The convergence *rate* moved.** Fitted over the same two endpoints,
>    0.7.2 gives ln(22.1925/9.2568)/ln(0.004/0.0018) = **1.10** — the
>    recorded O(h^1.2) — and 0.11 gives ln(21.8417/4.4605)/ln(2.2222) =
>    **1.99**. At the fine end 0.11 is **2.1× more accurate**, not less.
> 3. **The gated rung is an outlier on its own ladder.** The 0.11 fit
>    predicts 21.8417% · (0.0025/0.004)^1.99 = **8.6%** at h = 0.0025;
>    it measures 15.3848%, **1.8×** that. It is also the only rung whose
>    cell count moved appreciably (+0.92%, 6-30× the others).
>
> **Consequence for the band.** 15% was set as 1.18× a 12.75% measurement on
> a still-converging ladder. On 0.11 that ladder converges at \~O(h²) and the
> h = 0.0018 rung already reaches 4.46% — inside the 5% target `MAG-13`'s
> comment calls out of reach at ≈1.1M cells. So the honest disposal is
> probably *not* "loosen 15%": it is to find why h = 0.0025 sits off its own
> ladder on 0.11, since two rungs either side of it behave better than the
> record. **Nothing was touched here** — no band, no assertion, no record.
>
> **What is left to test.** Rank dependence and mesh instability are both
> excluded by the `-n 4` re-run above, so the remaining candidates are
> (a) a genuine non-monotonicity of this discretization near h = 0.0025
> that 0.11's slightly different mesh happens to land on, and (b) something
> in the sampling / point-location path that is sensitive to where a
> 1e-3-perturbed mesh puts cell boundaries relative to the 10 sample radii.
> (b) is the cheaper to test: re-run the gated rung on 0.11 with a
> different `n_points` (8 and 20) and see whether 15.3848% moves toward the
> ≈8.6% the 0.11 fit predicts. If it does, the band was measuring the
> sampler as much as the solve.
>
> **No band was touched.** Loosening 15% to accommodate a version bump
> would erase the only instrument that shows the ladder moved. §9 item 3a's
> "a moved gated physics number is a known-issues entry and a stop" clause
> fired on this.
>
> ---
>
> **Update, `OPS-18` step 3 attempt 5 (2026-08-22, 15:00 slot) — candidate
> (b) is answered, and it answers more than it was asked.** The probe now
> solves each rung once and samples the *same* field at every requested
> `n_points`, so a spread within a row is the sampler's alone. Errors at
> n_points **8 / 10 / 20**, `-n 2`, both images
> (`20260822T200411Z_…-wire-h0025-npoints.log` Status 0 31 s,
> `20260822T200503Z_…-wire-ladder-npoints-011.log` Status 0 106 s,
> `20260822T201014Z_…-wire-ladder-npoints-072.log` Status 0 126 s):
>
> | h | 0.7.2: 8 / 10 / 20 | 0.11: 8 / 10 / 20 |
> |---|---|---|
> | 0.0040 | 18.6850% / **22.1925%** / 20.9923% | 18.5328% / **21.8417%** / 22.0704% |
> | 0.0025 | **15.8028%** / **12.7485%** / 11.4984% | 16.6033% / **15.3848%** / 13.6986% |
> | 0.0018 | 11.5626% / **9.2568%** / 7.5722% | 4.9201% / **4.4605%** / 4.8086% |
>
> Every n_points = 10 column still reproduces its record (bold, ≤ 0.035%
> on 0.7.2), so the sweep is anchored rather than a new measurement.
>
> 1. **The metric is sampler-fragile on *both* images, and worse on the
>    old one.** The 10-point radial L2 spans 34% of its own value on
>    0.7.2's gated rung and 43% on its fine rung, versus 21% and 10% on
>    0.11. This is a property of a 10-sample radial estimator, not of the
>    upgrade.
> 2. **The 15% band already fails on 0.7.2 at n_points = 8** — 15.8028%,
>    on the image the record was taken on. The gate's 1.18× headroom is
>    *inside* the statistic's own sampler spread, so it was passing on a
>    sampler choice rather than on a margin. That is the more important
>    finding here than anything version-specific.
> 3. **The outlier survives the control.** At fixed n_points the 0.11
>    gated rung is worse at every count (8: 15.80 → 16.60; 10: 12.75 →
>    15.38; 20: 11.50 → 13.70), and no n_points brings it near the 0.11
>    fit's 8.6%. Candidate (b) is therefore **excluded** as the *cause* of
>    the outlier, alongside partitioning and mesh instability; candidate
>    (a) — a real non-monotonicity near h = 0.0025 — is what remains.
>
> **Still nothing touched** — no band, no assertion, no record. But the
> disposal question has changed shape: it is no longer only "may this
> record move", it is "is a 10-point radial L2 a gateable statistic at
> all". A band 1.18× a number that swings 34% under its own sampler is
> not measuring the discretization. The review owns that call; a
> defensible fix (raise `n_points`, or gate a sampler-independent norm)
> is a `MAG` chunk, not an `OPS-18` clause.
>
> **Ruling, 18:00 review 2026-08-22 — the gate is replaced, not
> re-banded.** "Loosen 15%" refused (standing rule; and 0.11 is the more
> accurate solver on this ladder, so a loosened band would record the
> wrong fact). "Chase the 0.11 non-monotonicity" refused as an `OPS-18`
> clause: the statistic it would chase swings 34% under its own sampler on
> the image that recorded it, so it cannot adjudicate a version bump either
> way. **`MAG-18`** is commissioned (§7; §9 item 1): an annulus-restricted
> domain L2 of `|B_h| − |B_ana|`, assembled not sampled, gated on a
> pre-registered rate ≥ 0.7 plus rank-independence and the natural-BC
> control, measured on 0.7.2 / `main` first; the 10-point number becomes
> reported-not-gated. `OPS-18` then re-measures leg 2 in that norm on 0.11
> (§9 item 4). **This entry stays open**: the h = 0.0025 non-monotonicity
> on 0.11 is observed and unexplained, and `MAG-18` does not claim to
> resolve it — it leaves only when someone explains or retires it with a
> commit.
>
> ---
>
> **Update, `MAG-18` executed 2026-08-22 (log
> `20260823T003518Z_MAG-18-full.log`, `7 passed` / 270.64 s / `-n 2`,
> 0.7.2 / `main`) — the ruling is implemented, and the symptom at the top
> of this entry can no longer fire.** `rel_error < 0.15` is gone from
> `test_straight_wire_b_field`; the number is printed, and reproduced
> under assertion at `n_points` 8 / 10 / 20 in
> `test_domain_l2_record`. So on 0.11 that test will now *pass* — the
> 15.3848% it used to fail on is still there, still unexplained, and this
> entry stays open for it, but a run hitting it will see a printed line
> rather than a red assertion. **Do not read the green as an explanation.**
> The gate is now `E_Ω` = 25.3787 / 10.7288 / **6.6708%** on the recorded
> ladder, rate **1.6842**, monotone, with the natural-BC control at
> 32.3117% vs 10.7288%.
>
> **Update, `OPS-18` step 3a attempt 7 (2026-08-23, 06:00 slot) — leg 2 is
> green on 0.11 and the entry still stays open.** On `attempt/OPS-18` at
> `44b5600`, real, `-n 2`, `20260823T111216Z_OPS-18-step3a-leg2-confirm.log`:
> **`11 passed, 4 skipped`** / 293.59 s / exit 0, against attempt 6's `1
> failed / 10 passed / 4 skipped`. `E_Ω` = 25.2868 / 10.6172 / 6.6458% on
> the recorded ladder, rate **1.6854**, monotone, natural-BC 32.315493% vs
> 10.617170% (ratio 0.3285). The `E_Ω` h = 0.0025 record is now written
> version-tagged as 1.061717e-01 at 147 235 cells and measures
> 1.0617170177e-01; the retired sampler control is keyed by image and
> reproduces the 0.11 triplet 16.603276 / 15.384842 / 13.698645% to
> ≤ 3.3e-06. **The 15.3848% non-monotonicity is exactly that middle
> number** — it is now recorded and asserted rather than failing, which is
> what "reported-not-gated" means. Still unexplained; this entry stays
> open until someone explains or retires it.

### ✅ RESOLVED 2026-08-23 (`OPS-18` step 3b) — `test_region_resolution_policy_refines_the_tagged_volumes_toward_cad` failed **only in the 0.11 image**: the uniform-sizing meshed volume moved 4.251e-04 relative against its `OPS-17` record (`OPS-18` step 2, 2026-08-22)

> **RESOLVED 2026-08-23, 15:00 implementer slot (`OPS-18` step 3b), by
> re-record under the class ruling (1\*)** — the disposition the entry
> below asked step 3 to make, made by comparison rather than argument.
> The image's gmsh moves three of the four uniform-sizing volumes;
> `UNIFORM_VOLUMES_RECORD` in `tests/mesh/test_mesh_tag_integrity.py` now
> carries the v0.11.0 values version-tagged beside the v0.7.2 ones —
> tag 1 1.191750413e-04 → **1.192257046e-04** (+4.251e-04 relative),
> tag 2 1.188402981e-04 → **1.185069486e-04** (−2.805e-03), tag 3
> 4.943767949e-04 **unmoved**, tag 4 1.143560787e-02 → **1.143589055e-02**
> (+2.472e-05). **The 1e-9 band is untouched**, as are the `GEO-17` sign
> and CAD-recovery gates, which stay green on their own digits (policy /
> uniform recovery 0.833417 / 0.755006 and 0.835563 / 0.750454, phantom
> 0.992751 / 0.983531). The identity the drift could have broken still
> closes exactly: the tagged-volume partition ratio is
> **1.000000000000** on both sizings and on the plain integrity mesh.
> Condition (b′) is met with room to spare — **every printed digit is
> identical across the two `-s` confirming runs of the slot**, run-to-run
> move 0.0 (`20260823T200533Z_OPS-18-step3b-confirm-run2.log`,
> `20260823T200550Z_OPS-18-step3b-confirm-run3.log`, both `7 passed,
> 4 skipped` / exit 0 / `-n 2`, both rank footers identical), against the
> red baseline `1 failed, 6 passed, 4 skipped` reproduced this same slot
> on the rebuilt image (`20260823T200356Z_OPS-18-step3b-confirm.log`,
> Status 1). Original entry follows.
>
> **Where this fires.** `tests/mesh/test_mesh_tag_integrity.py::test_region_resolution_policy_refines_the_tagged_volumes_toward_cad`,
> on the `attempt/OPS-18` worksite branch only. `main` boots 0.7.2 and is
> unaffected — this is not a failure on `main`, and no `main` run should ever
> see it.
>
> **Literal symptom** (`20260822T110624Z_OPS-18-step2-shim-runtime.log`,
> `1 failed, 6 passed, 4 skipped` in 15.85 s, `-n 2`, real build, both rank
> footers identical):
>
> ```
> AssertionError: uniform sizing moved tag 1 (coil_1) by 4.251e-04 against its
> OPS-17 record: 1.191750413e-04 -> 1.192257046e-04 m^3
> ```
>
> **Verified at** `cc431c9` + the step-2 shim, i.e. the first commit at which
> `io/mesh.py` imports under 0.11 at all — the `gmshio` collect error masked
> every runtime number before it, so this could not have been observed earlier.
>
> **Cause: the image's gmsh, not the migration.** The failing quantity is a
> *meshed volume* of a fixed CAD region under fixed sizing — a mesh-generator
> output, not a solved one. The 0.7.2 → 0.11 image change carries a new gmsh
> along with it, and the §9 item-3 trap clause anticipates exactly this class
> ("the new image carries a new gmsh — a moved *cell count* with identities
> intact is re-recorded with a note; a moved gated *physics* number is a
> finding"). The three identity-shaped assertions in the same file — tag
> integrity and the region-resolution policy pair — **pass unchanged**, so the
> tagging and the sizing policy are intact; only the record's last digits moved.
> The drift is 4.3e-04 relative, four orders below any gated physics band.
>
> **Not diagnosed further, and deliberately not fixed here.** `OPS-18` step 2's
> done-when is collect-level, and the re-record/finding disposition belongs to
> **step 3**, which owns §5.3's environment table and every gated band. No
> assertion, band or record was touched by step 2 — the record moving is
> evidence about the image, and re-recording it before the re-gate leg has run
> would spend the trap clause's discrimination for nothing.
>
> **For step 3:** this is the first observed instance of the predicted gmsh
> drift and it is in the *re-record* category on its face. Expect siblings in
> other volume/cell-count records; a gated *physics* number moving is the
> different animal the clause stops on.

### ✅ FIXED 2026-08-19 (`OPS-22`) — the magnetostatic loop-drive fixtures were complex-hostile: `ufl.max_value` / `<=` geometry predicates in the current-density callable (`OPS-17` leg (b2), 2026-08-19)

> **FIXED 2026-08-19, 04:30 implementer slot (`OPS-22` step 1).** All three
> files repaired; nothing marked `@real_only`. The complex build now runs
> `test_circular_loop.py`, `test_helmholtz_magnitude.py` and
> `test_helmholtz_v2.py` as one command — **5 passed in 412.12 s, exit 0**,
> both ranks identical (`20260819T094710Z_OPS-22-step1-complex-all.log`) —
> and the printed digits equal the real-mode record to the last figure
> (loop relL2 7.0658% / max 13.8212%; Helmholtz centre 0.728%, mean 0.644%,
> CV 0.1602%). Real mode is unmoved: identical digits in
> `20260819T093105Z` (before), `20260819T093529Z` (after the predicate fix)
> and `20260819T095414Z` (after the real-part fix, 5 passed / 199.91 s).
> The 5 tests are unblocked for `OPS-17` leg (b2).
>
> **The fix had two layers, and the second was not predicted here.**
> (1) As diagnosed: `ufl.max_value(rho, 1e-12)` → `ufl.sqrt(x²+y²+1e-24)`,
> and `(...) <= a²` → `ufl.le(ufl.real(...), a²)`. (2) With the form
> compiling, the complex run reached the *assertions* and failed at
> `ValueError: Unknown format code '%' for object of type 'complex'`:
> `evaluate_vector_field_parallel` returns the complex scalar type even for
> a real-valued magnetostatic solution. Both comparing tests now assert
> `max|Im B_z| ≤ 1e-12·max|B_z|` and compare on `np.real` (a no-op in real
> mode). **`OPS-20` should expect the same second layer** once its
> predicate is fixed, as should
> `examples/magnetostatics/02_circular_loop.py:173` and
> `04_helmholtz_analytic_comparison.py:79`, which still carry the idiom and
> are unexercised in complex mode (follow-up, not done in this slot).
>
> The **poisoned-stub trap** documented below stands unchanged and is
> independently useful; the sweep was clean before and after this work.
>
> *(Diagnosis, retained as the record of how this was found:)*
> **CAUSE DIAGNOSED 2026-08-19, 22:30 slot (attempt 2).** The heading below
> originally read "cannot JIT-compile one form … cause not diagnosed". Two
> `--tb=long` runs on a **verified-stub-free** cache localized it, and the
> repo already half-knew the answer. Both failure modes come from the same
> place — the **load form `L`** built at `src/fem_em_solver/core/solvers.py:385`
> (`LinearProblem(a, L, …)`) from the test's `current_density` callable, which
> encodes "inside the wire" with ordering comparisons that UFL forbids on
> complex-typed operands:
>
> ```python
> in_wire_1 = ((rho - MAJOR_RADIUS) ** 2 + (x[2] - z1) ** 2) <= MINOR_RADIUS**2   # ordering
> rho_safe  = ufl.max_value(rho, 1e-12)                                            # max on complex
> ```
> (`tests/validation/test_helmholtz_magnitude.py:83–87`, same shape at
> `test_circular_loop.py:54` and `test_helmholtz_v2.py:46`.)
>
> **Two distinct symptoms, one cause.** Which one you see depends on whether
> UFL's comparison checker catches it before FFCx starts:
>
> | file | symptom | cost |
> |---|---|---|
> | `test_helmholtz_magnitude.py` | `ComplexComparisonError: Ordering undefined for complex values.` (`ufl/algorithms/comparison_checker.py:49`) — raised in UFL, **before** compilation; the form repr shows `Conditional(OrCondition(LE(Sum(Power(…SpatialCoordinate…)))))` | **13.10 s**, then the ≈300 s non-collective exit hang (3b-xiii family) |
> | `test_circular_loop.py` | passes the checker, then `RuntimeError: Failed just-in-time compilation of form: Compilation failed on root node.` — the compiler's own message is swallowed by FFCx | **113.38 s**, nearly all in the doomed compile |
>
> **The codebase already carries the workaround** and its rationale, in
> comments written by earlier chunks: regularise *inside* the `sqrt` instead of
> with `ufl.max_value` — see `tests/validation/test_dodd_deeds_impedance.py:237–239`,
> `test_port_reaction_impedance.py:200–202`, `tests/mesh/test_two_torus_conforming.py:164`,
> all of which say in so many words that "the magnetostatic loop fixture's
> `max_value` form does not compile here". Those files were made complex-safe;
> these three (plus `examples/magnetostatics/02_circular_loop.py:173` and
> `04_helmholtz_analytic_comparison.py:79`, unexercised in complex mode) were
> not. So this is **fixture debt, not a solver defect** — no `src/` magnetostatic
> path is implicated, and real mode is green throughout.
>
> **Same family as the `test_coil_phantom_magnetostatics` entry below** (`You
> can't compare complex numbers with max.`), whose "cause not diagnosed / the
> comparison most likely enters through a DolfinX/UFL helper" line is now
> **superseded**: it enters through the drive callable, exactly as here.
> `OPS-20`'s step-1 diagnostic should start from that hypothesis.

**Verified at `c612920`, 22:30 implementer slot** — `20260819T033938Z_OPS-17-step3g-helmholtz-magnitude-isolated.log`
(exit 124 at the 300 s ceiling, but the full `--tb=long` traceback and the
`1 failed … in 13.10s` footer print before the hang) and
`20260819T034936Z_OPS-17-step3g-circularloop-onaxis-clean.log` (**exit 1**,
115 s, `1 failed, 2 deselected in 113.38 s`, both ranks). The cache was swept
with `find /root/.cache/fenics -name '*.c' -size 0` and the one stub deleted
**before** these runs, so neither reading is a stub artifact.

*(Prior text of this entry, from attempt 1 on 2026-08-19 at `e2295bf`, three
independent runs; the poisoned-stub trap below stands unchanged and is
independently useful.)*
Real mode is unaffected — these files are green in the real-mode leg (a) sweep.

| | |
|---|---|
| **Tests** | `tests/validation/test_circular_loop.py::TestCircularLoop::test_circular_loop_on_axis` (FAILED) and `::test_circular_loop_field_symmetry` (never returns; takes the window to `exit 124`); `tests/validation/test_helmholtz_magnitude.py::test_helmholtz_centre_field_magnitude` (FAILED, `ComplexComparisonError`); `tests/validation/test_helmholtz_v2.py::test_helmholtz_field_uniformity_two_torus` (never returns — it carries the same `max_value` idiom and was the test hanging batch A). **5 tests**, complex build only |
| **Symptom** | `RuntimeError: Failed just-in-time compilation of form: Compilation failed on root node.` on rank 0 (`dolfinx/jit.py:91`); rank 1 raises the same `RuntimeError` with `JIT compilation timed out, probably due to a failed previous compile. Try cleaning cache (e.g. remove /root/.cache/fenics/libffcx_forms_3b01242391fa699f45d97f502c916e1a1c96c1e6.c)`. Teardown then emits `AttributeError: 'LinearProblem' object has no attribute '_solver'` (a secondary, not the fault). The call is **109.07 s** — nearly all of it in the failing compile. |
| **Evidence it is not a cache artifact** | Three runs, same commit. (i) In a mixed batch: FAILED at 31%, then the next test in the file hung to `exit 124` (`20260819T021242Z_OPS-17-step3f-complex-validation-subset2.log`, 481 s). (ii) Isolated, `--tb=long`: `1 failed, 2 deselected in 109.58 s` (`20260819T022120Z_OPS-17-step3f-complex-circularloop-onaxis.log`, exit 1). (iii) **After deleting every 0-byte stub in `/root/.cache/fenics`**, the same test FAILED again and **re-created the identical hash `3b01242…` at 0 bytes** (`20260819T022356Z_OPS-17-step3f-complex-circularloop-repaired.log`, exit 124, 421 s). A cache artifact does not survive its own repair, and does not regenerate the same hash at zero length. |
| **Cause** | ~~Not diagnosed.~~ **Diagnosed 2026-08-19 attempt 2 — see the block at the top of this entry.** The complex-hostile geometry predicate (`ufl.max_value` / `<=`) in the fixture's `current_density` callable. FFCx still swallows the root-node compiler message for `test_circular_loop`, so the *compiler's* words remain unrecovered; the offending construct does not. |
| **Second-order damage — the poisoned-stub trap** | A 0-byte stub left by any killed compile makes **every later run that needs that form** fail this way rather than recompiling, and the message blames the cache, not the form. A stale stub from **2026-08-18 14:02** (leg (b1) attempt 2's era) was still sitting in the cache when this slot started — i.e. this has been silently mis-attributed before. Sweep with `find /root/.cache/fenics -name '*.c' -size 0` before trusting any "JIT compilation timed out" message; that check is cheap and should precede any cache-clear. Note this cuts **against** reflexively clearing `~/.cache/fenics`: the targeted delete is the diagnostic, and here it *exonerated* the cache. |
| **Not** | Not the `>12× real` cost rule (withdrawn 2026-08-18) and not a solve cost — 109 s of compile, 0 s of setup. Not a physics or tolerance failure: no assertion is ever reached. Not a real-mode issue. |
| **Fix** | Not attempted — `OPS-17` is a bookkeeping leg and does not touch `tests/`. The fix is mechanical and already precedented in this repo: regularise inside the `sqrt` (`test_dodd_deeds_impedance.py:237`, `test_port_reaction_impedance.py:200`) instead of `ufl.max_value`, and express the wire predicate without an ordering comparison on complex-typed operands. Whoever takes it should decide per file whether the complex build needs the magnetostatic path at all — `@real_only` is a legitimate, cheaper disposition (same call `OPS-20` faces). Do **not** conclude the forms are un-compilable: the sibling files prove the same physics compiles once the predicate is complex-safe. |
| **Resolves with** | **`OPS-22`** (commissioned 2026-08-19 03:00 review; §7 entry with the full rubric): per-file fix-or-mark of the three fixtures, real-mode records asserted unmoved, complex run completing with a footer. `OPS-20` stays a separate chunk (same family, different file) and its step 1 is re-pointed at the drive callable. Blocks 5 tests in `OPS-17` leg (b2) from being observed in a completed leg until it lands. |

### 1. ✅ RETIRED 2026-07-31 — stale test double, `DummyMagnetostaticSolver`

The two phantom tests
(`tests/materials/test_phantom_material_model.py::test_phantom_material_assignment_and_time_harmonic_pipeline_wiring`,
`tests/post/test_phantom_field_metrics.py::test_phantom_field_metrics_and_exports_are_finite`)
failed with `AttributeError: 'DummyMagnetostaticSolver' object has no attribute
'last_solve_diagnostics'`. No fix ever landed: `TH-1` steps 1–3 deleted the
`last_solve_diagnostics` read from the time-harmonic path, so the failing code
path went away under the entry.

Re-verified through the harness at `424faed`:
`20260731T170152Z_KI-1-retire-gate.log` — complex build, `FEM_EM_REQUIRE_COMPLEX=1`,
`-n 2`, **10 passed in 4.6 s**; and `20260731T170140Z_KI-1-real-mode-iomatpost.log`
— real build, the CI `validation` job's exact command with the `--deselect`s
removed, **15 passed, 2 skipped in 0.5 s**. Both `--deselect`s are gone from the
`validation` job and both files are now listed in `validation-complex`, which is
where the two `@complex_only` tests actually execute. The test-side
`ComplexWarning` casts were fixed in the same commit; the remaining one at
`post/phantom_fields.py:88` was recorded under `POST-1` in `PROJECT_PLAN.md` §7
and is **fixed as of 2026-08-04** (`POST-3` step 4 — both cast sites removed,
statistics taken on the phasor magnitude, gated by
`tests/post/test_phantom_phasor_semantics.py`). `POST-1` stays ⚠️ for the
interface-guardrail machinery, which is unrelated to the cast.

The heading is kept (rather than deleted with the entries renumbered) so the
numbering of entries 2–6 stays stable — several commits and CI comments refer to
them by number.

### 2. ~~Residual-trend classifier disagrees with its test~~ — RESOLVED 2026-08-08 (`OPS-12`)

Both failures were on the **code** side, and they were two different defects,
not one:

1. `classify_residual_trend()` carried undocumented thresholds
   (`f >= 0.75` ⇒ `mostly-decreasing`, `f >= 0.5` ⇒ `mixed`, where `f` is the
   fraction of non-increasing steps). Nothing — not the docstring, not the
   label names — licensed that asymmetry: it gave increases a band of width
   0.5 and decreases 0.25, so a history that decreased on two of every three
   iterations was reported `mixed`, and **no history of four or fewer samples
   could reach `mostly-decreasing` at all**. The labels now partition by the
   sign of `f - 0.5`, which is what their names mean, and the docstring states
   the table. The test's `[1.0, 0.4, 0.45, 0.1] -> mostly-decreasing`
   (`f = 2/3`) was right and is unchanged.
2. The second failure was never `assert diagnostics is not None` — the
   recorded symptom was wrong. It was `assert diagnostics.converged`, with
   `converged_reason = -3` (`KSP_DIVERGED_ITS`): the fixture asked for
   gmres+jacobi at `ksp_rtol = 1e-8` with `ksp_max_it = 300`, and that solve
   needs **1409** iterations on the 1405-cell fixture. The cap was
   under-resourced; the assertion was not wrong and was not touched.

A third defect fell out of the diagnosis: the time-harmonic path never called
`ksp.setConvergenceHistory()` (the magnetostatic path always has), so
`residual_history` came back **empty** and `residual_trend` was permanently
`unavailable` — the classifier was unreachable in production and the test's
membership assertion passed vacuously. Armed, and the test now gates
`len(history) == iterations + 1` and `trend == classify_residual_trend(history)`.

| | |
|---|---|
| **Fixed in** | `OPS-12`, 2026-08-08 — `src/fem_em_solver/core/solvers.py`, `src/fem_em_solver/core/time_harmonic.py`, `tests/solver/test_convergence_diagnostics.py` |
| **Gate** | 18 passed at `-n 2` under the complex build in 0.93 s (`20260808T050622Z_OPS-12-gate-final.log`); the classifier identity is asserted with `==` on an 11-row parameterized family spanning both sides of `f = 0.5` and of the retired `f = 0.75`, with negative controls |
| **CI** | `tests/solver/test_convergence_diagnostics.py` is now in the `validation-complex` job, which is what this entry's old status line named as its exit condition |

### 4. ✅ RETIRED 2026-08-08 — coil+phantom B-field symmetry exceeded tolerance (`MAG-6` step 3)

**Resolution.** The estimator, not the physics, owned every number below. Both
sampled metrics in the test now go through **DG0** instead of CG1, and the
fixture is refined one rung to `resolution = 0.010` m; **no tolerance moved**
(`PHANTOM_SYMMETRY_REL_TOL = 0.35` and `PHANTOM_CENTERLINE_JUMP_RATIO_MAX =
0.60` are untouched, and the symmetry assertion was *tightened* — the
permissive `or max_abs_diff < …` escape is gone). Red first at `-n 1`
(`max_rel_diff=0.728`, `20260808T170126Z_MAG-6-step3-redbaseline-n1.log`);
green at `-n 1/2/4` with `max_rel_diff` = **0.323844 / 0.302661 / 0.308407**,
three-way spread **7.00%** against the pre-registered ≤ 10%, the `-n 1`
value byte-reproducing the step-2 record
(`20260808T170549Z_MAG-6-step3-gatefinal-n1.log`,
`…170529Z_…-n2.log`, `…170515Z_…-n4.log`).

**What this test now gates: discretisation symmetry, not phantom physics.**
`mu` is uniform, so the phantom is invisible to the solve and the exact
mismatch is 0 by construction. The fixture caveat below is preserved in the
test's module docstring; it was not "fixed", it was made explicit.

**Carried forward, not closed:** the centerline smoothness metric is *also*
CG1-owned and was re-pointed at DG0 in the same commit (CG1 read 0.705 and
0.732 on two identical `-n 4` runs against a 0.60 bound — rank-dependent and
not run-to-run reproducible; DG0 reads 0.227869 on the same solve). DG0 leaves
it passing at all three rank counts but still **88%** rank-scattered
(0.473300 / 0.268765 / 0.251746 at `-n 1/2/4`) — no rank-stability claim is
made for it, and sizing that second estimator is unscoped work for a review.
Related: the DG0 symmetry metric moved 6.8% between two identical `-n 2` runs
(meshing is not bit-reproducible), inside the spread gate but relevant to any
future tightening.

<details>
<summary>Historical record (the failing entry as it stood)</summary>


| | |
|---|---|
| **Test** | `tests/validation/test_coil_phantom_bfield_metrics.py::test_coil_phantom_bfield_metrics_are_finite_smooth_and_symmetric` |
| **Symptom** | **Rank-dependent** (`MAG-6` step 1, 2026-08-08). `max_rel_diff` reads **0.727907 at `-n 1` (fails)**, 0.240541 at `-n 2` (passes), 0.321468 at `-n 4` (passes), against the same `0.350` tolerance on the same 19 792-cell mesh. Historical record: 0.557 against 0.350, with `max_abs_diff=7.090e-07` against `6.523e-08`. |
| **Cause** | **The test's estimator, partially diagnosed** (`MAG-6` step 1). `curl A` for N1curl degree 1 is cell-wise constant, and the test interpolates it into **CG1** before sampling — a nodal value where the field jumps, supplied by whichever cell the partition happens to give. Sampling the same field at the same points through **DG0** is rank-stable (0.513648 / 0.534746 / 0.538472 at `-n 1/2/4`, 4.8% spread) where CG1 swings **3.03×**; the assembled `‖B‖_L2` moves 1.84% through CG1 but only **0.09%** through DG0, so the solve is not the suspect. What the ≈0.53 rank-stable residual itself *is* was answered by step 2 (row below): **coarse-mesh discretisation**, falling at p ≈ 1.07 under refinement. Tracked as `MAG-6` (legacy A1). |
| **Do not read the green `-n 2` run as a fix** | The test passes in CI (which runs `-n 2`) for a non-physical reason. The rank-stable estimate 0.51–0.54 is close to the historical 0.557 and **above** the 0.350 tolerance; the record was not wrong, the estimator drifted into a lucky partition. Retiring this entry on the strength of a green `-n 2` run would bury a real finding. |
| **Ruled out** (`MAG-6` step 1) | **Boundary-mirror artifact:** on the rank-stable DG0 path, growing `air_padding` 1.5× moves the metric 0.534746 → 0.534772, a **0.005%** move against the "≥ 2× drop" the hypothesis predicted. **Gauge penalty:** the test solves at `gauge_penalty=1e-3`, 1000× below `DEFAULT_GAUGE_PENALTY = 1.0` and the source of the 9 `GaugeContaminationWarning`s per run, but re-solving at 1.0 moves `max_rel_diff` 0.240541 → 0.241846 (`-n 2`) and `‖B‖_L2` by 0.016% — `MAG-10`'s catastrophe needs degree 2, this fixture is degree 1. |
| **Fixture caveat** | `MagnetostaticProblem` is built with a **uniform** `mu = MU_0`, so the phantom is physically invisible — the test's off-centre-phantom asymmetry control moves only the mesh, and on the rank-stable path it *decreases* the metric (0.476684 vs 0.534746) instead of increasing it. Any symmetry claim on this fixture is a claim about discretisation, not about the phantom. |
| **Note** | The coil+phantom fixture uses a single global `setSize` and tight air padding — the same pattern that cost **20% error** on Helmholtz until `air_padding` was decoupled (see `docs/validation/helmholtz.md`). Kept for the air-box generalisation follow-up; it is *not* the owner of this metric. |
| **Verified pre-existing at** | `ce92e8c` and earlier (0.559 at `HEAD`, 0.557 after the gauge change — the gauge default is not the cause). Reproduced 2026-08-08 at `ed2d7e5` in the rank-dependent form above: `20260808T033401Z_MAG-6-step1-rankcheck.log`, `20260808T034013Z_MAG-6-step1-sampling.log`, `20260808T033802Z_MAG-6-step1-gauge.log`. Probe: `scripts/probes/mag6_step1_probe.py`. |
| **What the ≈0.53 is** (`MAG-6` step 2, 2026-08-08) | **Coarse-mesh discretisation, measured.** On the fixture's own resolution knob with the probe grid frozen, the DG0 metric falls **monotonically**: 0.534746 (h = 0.015 m, 19 792 cells) → 0.312197 (h = 0.010, 55 784) → 0.255165 (h = 0.0075, 124 179) at `-n 2`, total ratio **2.0957**, observed rate **p = 1.067** — the O(h) that N1curl degree 1 allows. Independently at `-n 4`: 0.537750 → 0.304356 → 0.292706, ratio 1.8372. Rung 1 byte-reproduces step 1 (0.534746 at `-n 2`, 0.513648 at `-n 1`). **Caveat, load-bearing:** the finest rung is **void by its own rank control** — 14.71% spread between `-n 2` and `-n 4` against the ≤ 10% band — while its assembled `‖B_dg0‖_L2` moves only **0.0079%** across the same ranks, so the sampling, not the solve, destabilises as cells shrink under a fixed grid. The band holds on the two controlled rungs alone (4.69% and 6.40% spread; 1.713× drop across a 1.5× refinement). |
| **CG1 does not converge either** (`MAG-6` step 2) | On the identical solves, the test's own CG1 path reads 0.240541 → 0.760519 → 0.723637 at `-n 2` — **non-monotone and mostly rising** under refinement where DG0 falls. Step 1 showed CG1 owns the rank-dependence; refinement does not rescue it. |
| **The number the tolerance question was waiting on** | The DG0 metric meets the **unmodified 0.350** at **h = 0.010 m** (0.312197 at `-n 2`, 0.304356 at `-n 4`) for 55 784 cells, 6.4 s mesh + 2.0 s solve at `-n 2` — standard tier. So candidate (i) (re-point the estimator at DG0) plus one refinement rung gates green **without raising 0.350**. The choice is still a review's; nothing was re-pointed and the tolerance is untouched. |
| **Verified at** (step 2) | `20260808T123206Z_MAG-6-step2-meshprobe.log` (cell counts), `20260808T123245Z_MAG-6-step2-hconv-n2.log` (35 s), `20260808T124355Z_MAG-6-step2-hconv-n4.log` (32 s), `20260808T123335Z_MAG-6-step2-hconv-n1.log` (**exit 124** — sequential LU 13.0 s → 132.4 s over rungs 1–2, rung 3 past the 600 s ceiling; that rung's `-n 1` control was dropped on cost per §7, not retried longer). |
| **Estimator adjudicated 2026-08-08 (10:30 review)** | Candidate (i) is taken on step 2's licensed numbers: `MAG-6` step 3 (§7, queued §9 item 1) re-points the test's sampling at DG0 and refines the fixture to h = 0.010 m, gating against the **untouched** 0.350 with a ≤ 10% three-way rank-spread gate (6.40% on record). Candidate (ii) — real material contrast — is rejected for this chunk: it changes the physics under the metric; the uniform-μ caveat transfers into the test docstring instead. **This entry retires with step 3's landing commit.** |
| **Superseded next step scoped 2026-08-08 (03:00 review)** | `MAG-6` step 2 (§7, queued §9 item 3): measure the **DG0 metric's `h`-convergence** on a three-rung refinement ladder — the one measurement that decides whether ≈0.53 is coarse-mesh discretisation (it falls with `h`) or a mesh-independent defect (it plateaus). Candidates (i) re-point the estimator at DG0 and (ii) give the phantom real material contrast are deliberately deferred behind it: (i) would flip the test hard-red today with no licensed number to gate against, and (ii) changes the physics under the metric. The 0.350 tolerance stays untouched until step 2's bands read out. |

</details>

### 5. ✅ RETIRED 2026-08-06 — domain sizing heuristic, off-centre phantom

| | |
|---|---|
| **Test** | `tests/mesh/test_domain_sizing_heuristics.py::test_coil_phantom_domain_sizing_accounts_for_off_center_phantom_extent` |
| **Symptom** | `assert 0.09 > 0.09` |
| **Was pre-existing at** | `794d2f1` (pre-session); reproduced one last time 2026-08-06 at `d4e278d` (`20260806T033155Z_GEO-4-step1-precontrol.log`, 1 failed 3 passed in 1.31 s) |
| **Diagnosis** (`GEO-4` step 1) | **The test's assertion, not the arithmetic.** The air box is centred on the origin, so its half-width is `max(coil_major + coil_minor, |offset| + phantom_radius) + padding` — the off-centre phantom enters through the second term of the max. `coil_phantom_domain` rejects any placement with `|offset| + phantom_radius >= coil_major - coil_minor` (the `radial_clearance <= 0` guard), so the phantom's outer radius is **always** strictly below the coil's outer radius and the max is always won by the coil. The property "an offset phantom grows the box" is therefore unattainable for every meshable configuration, not merely unexercised by the chosen 0.03 m offset (phantom reaches 0.07 m vs the coil's 0.09 m). Test and code landed together in `2c52f05`; the test never passed. |
| **Resolution** | The strict `>` was **not** relaxed. The test now gates the containment identity with the clearance term explicit — `half_width == max(coil_outer, |offset| + r_phantom) + padding` for both presets, and `clearance(centered) - clearance(shifted) == 0.03` exactly (the whole offset is spent out of the phantom's wall clearance). Two tests added: the phantom-governed branch of the max keeps a strict `>` alive where it exists (arithmetic only, outside the meshable envelope), and a zero-padding negative control still detects an undersized domain. `coil_phantom_domain_sizing_diagnostics` gained `phantom_outer_radial_extent_m`, `phantom_boundary_clearance_m`, `phantom_offset_radius_m`, `phantom_governs_radial_extent`; **no sizing number changed**, so no meshed fixture moved. |
| **Verified at** | `20260806T033316Z_GEO-4-step1-gate.log` (6 passed, 1.36 s, `-n 2`) and `20260806T033327Z_GEO-4-step1-mesh-regression.log` (whole `tests/mesh`, no `--deselect`, **27 passed 1 skipped in 85.3 s**) |
| **Back in CI** | The `OPS-11` `--deselect` was removed in the same commit; the `Mesh generation suite` step now excludes nothing. |
| **Cited wrongly as "known-issues 6"** | Commit `3ac025c` and `docs/testing/attempts.md:1903,1907` call this entry 6. **It was entry 5.** Entry 6 is the rank-dependent single-port excitation test in `tests/solver`, unrelated to `tests/mesh`. |

**Open follow-up for a review (not a failing test):** because the overlap guard
is z-blind, a *short* phantom that would clear the torus tubes in z is rejected
just the same. If off-centre placements that radially govern the box are ever
wanted, that guard — not the sizing heuristic — is what must change.

### 8. ✅ RETIRED 2026-08-05 — magnetostatic energy raised `TypeError` in the complex build

**Fixed 2026-08-05 (16:30 run) by `MAG-16`.** Both tests pass at `-n 2` under
`dolfinx-complex-mode` with their identity assertions untouched
(`20260805T213601Z_MAG-16-gate-complex-final.log`, 10 passed in 4.9 s with
`tests/environment` first), and the complex-mode `tests/solver` sweep went from
4 standing failures to 2 — the remaining two are **entry 2**, unrelated
(`20260805T213408Z_MAG-16-regress-complex.log`, 2 failed 34 passed in 28.3 s).

`compute_magnetic_energy` now takes `np.real` of the reduced scalar and raises
if `abs(Im W)/abs(Re W)` exceeds `ENERGY_IMAG_RTOL = 1e-8`; `abs()` was
rejected on purpose, since it would swallow a genuine imaginary part *and* a
negative real one. The suspicion recorded below — that the imaginary part is
round-off — was **too pessimistic**: it is **exactly 0.0** in both gauges,
because the magnetostatic load is real and `ufl.inner` conjugates its second
argument, so the assembled integrand is `μ⁻¹|curl A|²/2` and the complex build
merely stores a real number in a complex slot.

The value had indeed never been compared across builds, so it was pinned before
being trusted: the real-build energies were captured *pre-fix*
(`20260805T213144Z_MAG-16-probe-real.log` — `1.121469318858e-08 J` penalty,
`1.121466766900e-08 J` Lagrange) and the complex build reproduces them to
`2.9e-07` and `1.3e-13` respectively. The penalty gauge's run-to-run wander
(`1.9e-08…2.9e-07`) is its κ ≈ 1e10 operator, not the reduction.

The file is now listed in the `validation-complex` CI job, which is what stops
this from recurring — nothing had ever run it under the complex build until a
`POST-3` step-5 regression sweep did by hand. Original entry follows.

| | |
|---|---|
| **Tests** | `tests/solver/test_energy_and_point_evaluation.py::test_energy_matches_explicitly_reduced_assembly`<br>`tests/solver/test_energy_and_point_evaluation.py::test_energy_satisfies_discrete_work_energy_identity` |
| **Symptom** | `TypeError: float() argument must be a string or a real number, not 'complex'` at `src/fem_em_solver/core/solvers.py:661` (`MagnetostaticSolver.compute_magnetic_energy`) |
| **Cause** | Not diagnosed beyond the mechanism: in the complex build every `fem.Function` is complex, so the assembled energy scalar is complex-typed (with a round-off imaginary part) and the unconditional `float(...)` refuses it. The magnetostatic energy is real by construction, so the fix is presumably `np.real(...)` before the cast — but that is `MAG` work and the value has never been checked against the real-build number, so it is recorded rather than patched in passing. Both tests pass in the real build. |
| **Verified pre-existing at** | `aabb0a7` — reproduced with the `POST-3` step-5 diff stashed: `2 failed, 2 passed in 4.46 s` (`20260805T003945Z_POST-3-step5-preexisting.log`, `-n 2`, complex build). Found by that step's regression sweep, which is the first time this file was run under `dolfinx-complex-mode`. |

Owned by chunk `MAG-16` (§7, written by the 2026-08-05 10:30 review); this
entry leaves with `MAG-16`'s fixing commit — which is the commit carrying the
retirement header above.

### 9. ✅ RETIRED 2026-08-05 — `-n 2` hang on `two_torus_domain`'s port facets

**Diagnosed and fixed 2026-08-05 (12:00 run); `tests/mesh/test_two_torus_port_facets.py`
is on `main` and green at `-n 2` in 20 s
(`20260805T171107Z_PORT-1-step3biv-parallel-gate-fixed.log`, 2 passed).**

**Cause: a lazy collective, reached on only one rank.** An interior-facet
(`dS`) assembly requires `Topology::create_entity_permutations()`, and the
dolfinx assembler calls it *lazily* — only once a rank finds integration
entities for the form's subdomain id. Under this fixture's partition each rank
owns the facets of exactly one port (rank 0 → tag 201, rank 1 → tag 202), so
assembling tag 201 put rank 0 inside that collective while rank 1, with no
201 facets, sailed past it into the next one. Hence the two ranks' *different*
SIGTERM stacks (`create_entity_permutations` vs mpi4py `MPI_Comm_dup`) — a
mismatched collective, not a slow one. The fix is one hoisted line in
`_facet_group_area`: call `create_entity_permutations()` unconditionally on
every rank before building the form.

**What the discriminating run was.** The same computation, marker-instrumented,
ran to completion as a *script* at `-n 2` (exit 0, 12 s,
`20260805T170545Z_PORT-1-step3biv-dS-localise.log`) while the pytest gate hung
— and the script's only extra call was the explicit
`create_entity_permutations()`. Markers inside the gate then pinned the hang to
`_facet_group_area` at tag 201
(`20260805T170743Z_PORT-1-step3biv-pytest-localise.log`, exit 124).

**The ghost-mode hypothesis was necessary but not sufficient.** The
`shared_facet` partitioner now plumbed into `two_torus_domain` does what the
entry below predicted — `cells_ghost` 0 → 239/231 per rank
(`20260805T170109Z_PORT-1-step3biv-ghostprobe.log`, 14 s) — but the gate still
hung with it alone (`20260805T170140Z_…-parallel-gate.log`, exit 124, 181 s).
Both changes are kept: an interior-facet integral does need both cells of every
facet, and the fixture's docstring records the requirement.

**Generalisation, untested elsewhere:** any `dS` integral over a subdomain that
some rank does not touch is exposed to this. Only this fixture is fixed.

#### Superseded diagnosis, 2026-08-05 (22:30 run)

**Retitled and half-refuted 2026-08-05 (22:30 run).** The mesh is innocent: with
the gmsh dim-2 physical groups removed entirely and the identical facet set
rebuilt on the dolfinx side from the distributed cell tags, `model_to_mesh`,
`create_entities(fdim)` and `create_connectivity(fdim, tdim)` all return at
`-n 2` on the gate's own mesh in **14 s** (marker probe
`20260805T034007Z_PORT-1-step3biv-hang-localise-fine.log`, exit 0; 39578/39956
cells per rank, 116 interface facets found per port). The `-n 2` gate still
times out, so *something* on this path hangs — but it is downstream of the tags,
in the `dS` facet-area assembly, and the paragraph below misattributes it.

**Leading hypothesis for the next attempt, measured not guessed:**
`gmshio.model_to_mesh` passes no partitioner, so the mesh is built with the
default ghost mode and the probe measures `cells_ghost=0` on **both** ranks. An
interior-facet (`dS`) integral needs both cells behind every facet, which a mesh
with no ghost cells cannot supply on a partition-boundary facet. First move:
hand `model_to_mesh` a `shared_facet` partitioner and re-measure. Second
observation from the same probe, relevant either way: with the current
partition each rank sees exactly **one** port (rank 0 tag 201, rank 1 tag 202),
so any per-port assertion is rank-local until it is reduced.

The original entry follows, kept because its serial measurements stand.

#### Original entry, 2026-08-05 (21:00 run) — the `model_to_mesh` attribution is superseded above

| | |
|---|---|
| **Tests** | `tests/mesh/test_two_torus_port_facets.py` (both tests) — **not on `main`**; the file and the mesh change were parked *because* of this, so nothing on `main` is red. *(The branch this row named, `…021000Z`, was deleted by the 2026-08-05 10:30 review; the current code is on `attempt/PORT-1-step3biv-20260805T034500Z`.)* |
| **Symptom** | With `PORT-1` step 3b-iv's dim-2 physical groups `201`/`202` (the gap↔conductor shared surfaces, which are **interior** facets) added to the gapped fixture, `mpiexec -n 2` hangs inside `gmshio.model_to_mesh` and is killed by `timeout` at the 180 s ceiling. Both ranks' loguru stacks are identical and spinning in `MPI_Testall` ← `MPI::compute_graph_edges_nbx` ← `IndexMap::index_to_dest_ranks` ← `Topology::create_entity_permutations` ← `create_entities`. gmsh finishes (`Done optimizing mesh (Wall 7.14s)`) ≈10 s in; the remaining ≈168 s is the hang. No test code runs — the hang is before `model_to_mesh` returns. Log: `20260805T020301Z_PORT-1-step3biv-costprobe.log`, exit 124. |
| **Cause** | Not diagnosed. Bounded from two sides by measurement: `-n 1` completes the identical case in 22.5 s with correct areas (`20260805T020843Z_PORT-1-step3biv-serial-gate.log`, 2 passed), and `-n 2` on the same fixture **without** the new facet groups is green today (`tests/mesh/test_two_torus_gapped.py`). So it is neither cost nor the gapped geometry: it is the distribution of facet tags whose facets are interior to the partitioned mesh. The `2xx` groups are the fixture's first interior dim-2 groups — the only pre-existing one is the outer boundary. |
| **Verified pre-existing at** | Not pre-existing — introduced by the parked branch, which is why it is parked. Recorded here so the next attempt starts from the stack trace instead of re-deriving it. |

### 10. ✅ RETIRED 2026-08-06 — `two_torus_domain`'s outer-boundary facet group never reached the dolfinx facet tags (`GEO-10`)

The group was never *declared*, so nothing downstream could lose it. gmsh
inflates an OCC entity's bounding box by its geometric tolerance — measured at
exactly **`1.000e-07`** on all six walls of this box
(`20260806T050143Z_GEO-10-probe.log`) — and the fixture's flat-against-wall
test used `tol = 1e-9`. All six walls failed it, `boundary_surfaces` came out
empty, and the `if boundary_surfaces:` guard silently skipped
`addPhysicalGroup`. The chunk's prime suspect, fragment renumbering, is
**refuted**: the group is re-derived from bounding boxes after `fragment` +
`synchronize`, so renumbering never reaches it.

Fixed by widening that one tolerance to `1e-6` — 10× above the measured
padding and four orders below the nearest interior face's `2.000e-02`
residual, so the interior-face protection the tight test existed for is
intact. ~~No other fixture is affected: the rest of `io/mesh.py` uses a
`< resolution` wall test (loose by ≈4 orders), and only `two_torus_domain` had
tightened it.~~ **Corrected 2026-08-06 by `GEO-11` measurement:** two other
fixtures had also tightened it to `1e-9` and have the identical defect —
`loop_over_half_space_domain` and `sphere_in_box_domain`, see entry 12. This
retirement stands (its own gate is unaffected); only the generality claim was
wrong.

Gated by `tests/mesh/test_two_torus_outer_boundary.py`: tag sets exactly `{1}`
ungapped and `{1, 201, 202}` gapped, and the assembled `ds` area over tag `1`
equals the analytic box surface `2(LW+LH+WH) = 3.220000000000e-02 m²` at
ratio **`1.000000000000000`** (`-n 2`,
`20260806T050313Z_GEO-10-gate-n2.log`, 25 s) and `1.000000000000001` (`-n 1`,
`…050350Z_…-n1.log`, 24 s) — planar walls, so this is an identity at `1e-9`,
not a band.

The open question the entry recorded is now answered: **neither Helmholtz
consumer depends on tag `1`.** Both were re-run with the group present and
their gated numbers are digit-identical — `MAG-14`'s centre-field error is
still `0.728%` (`…050656Z_GEO-10-helmholtz-regression.log`, 2 passed, 11 s).
The port-facet gate likewise reproduces `A_201 = A_202 = 1.563786482e-04 m²`
at `0.974490841` of analytic (`…050620Z_GEO-10-portfacet-digits.log`), so
adding a boundary group moved no interface tag. Full `tests/mesh` at `-n 2`:
**29 passed, 1 skipped, 107.64 s** (`…050421Z_GEO-10-mesh-regression.log`).

### 12. ✅ RETIRED 2026-08-06 — `loop_over_half_space_domain` and `sphere_in_box_domain` never declared their `outer_boundary` group (`GEO-12`)

Fixed by widening both `tol` from `1e-9` to `1e-6` (`io/mesh.py`, the two sites
below) — exactly `GEO-10`'s fix, for the identical measured cause. Post-fix
(`20260806T183203Z_GEO-12-probe.log`): the loop fixture accepts **10 of 12**
dim-2 entities (the cube's four `z = 0`-split sides plus top and bottom; the
two rejected are the torus surface at `9.000010e-02` and the air/slab interface
at `1.000001e-01`) and the sphere fixture **6 of 7** (the sphere surface is
rejected at `1.500001e-01`). Both land on wall ratio
`1.0000000000287557e-01` — the same `1e-7/1e-6 = 0.1` `GEO-10` designed — so
`tests/mesh/test_boundary_classification_margins.py` now *asserts* the
two-sided margin for them instead of pinning a defect; the two `pytest.skip`s
are gone.

The group was invisible because nothing gated it, so the tolerance landed with
a meshed gate: `tests/mesh/test_wall_boundary_tag_areas.py` asserts per fixture
that facet tag `1` exists, carries an allreduced facet count > 0, and that its
assembled `ds` area equals the analytic cube surface `6(2W)² = 2.4e-01 m²` —
an identity on planar walls, gated at `1e-9` relative. Measured at `-n 2`
(`20260806T183328Z_GEO-12-gate.log`, 3.2 s): loop **1958 facets, ratio
1.000000000000000**; sphere **988 facets, ratio 0.999999999999999**.

The latency claim is now measured, not assumed. All six downstream callers plus
`tests/post/test_drop_set_semantics_sphere.py` were re-run and **no landed
number moved a digit**: `MAT-6` step 3 `dR` rel. error `1.5834%` / `dX` ratio
`0.9200`; step 4 `1.5763%` / `0.9849` (projected) and `1.5713%` / `0.8740`
(pinned), identical to `20260805T200455Z`/`20260805T200938Z`; `TH-8`/`MAT-4`
mass-averaged SAR ratio `0.999846` and the `POST-1` sphere table `4.2530%`
(`20260806T183745Z_GEO-12-callers-A.log`, 24 passed 210 s;
`20260806T184151Z_GEO-12-callers-B.log`, 8 passed 574 s). Whole `tests/mesh` at
`-n 2`: **35 passed, 2 skipped, 118.29 s**
(`20260806T183404Z_GEO-12-mesh-regression.log`).

Entry 13 is **not** covered by this fix and stays open — `cylindrical_domain`'s
margin is a different mechanism (tolerance coupled to `resolution`).

| | |
|---|---|
| **Test** | None fails. Measured and pinned by `tests/mesh/test_boundary_classification_margins.py` (`GEO-11`), which **skipped** the margin assertion for these two after pinning their numbers |
| **Symptom** | Both fixtures' wall test classifies **zero** surfaces — 0 of 12 and 0 of 7 dim-2 entities — so `boundary_surfaces` is empty and the `if boundary_surfaces:` guard silently skips `addPhysicalGroup`. Facet tag `1` does not exist in the returned `facet_tags` |
| **Cause** | Measured, not guessed (`20260806T140325Z_GEO-11-probe.log`): identical to retired entry 10. Both use `tol = 1e-9` (`io/mesh.py` \~lines 1384 and 1532) against gmsh's OCC bounding-box padding, which is the same **`1.000e-07`** `GEO-10` measured on `two_torus_domain` — the tolerance sits **100× below** the padding it must clear. Retired entry 10's closing claim that "no other fixture is affected" was wrong: it checked the `< resolution` fixtures and missed these two, which had also tightened the test |
| **Live impact** | **None — latent.** Every caller of both generators discards the facet tags (`msh, cell_tags, _ = MeshGenerator...`) in `test_dodd_deeds_impedance.py`, `test_dodd_deeds_projected_drive.py`, `test_dodd_deeds_reactance_box_size.py`, `test_dielectric_sphere.py`, `test_lossy_sphere_sar.py`, `test_mass_averaged_sar.py`, and imposes its wall condition geometrically instead. No landed `MAT-6`, `TH-8` or `MAT-4` number reads the missing group, so none of them is wrong |
| **Verified at** | `main` as of `2cad984` |
| **Fix, not taken in-slot** | Widen both to `1e-6`, exactly as `GEO-10` did — the nearest interior faces sit at `9.000e-02` and `1.500e-01`, so `1e-6` keeps 5 orders of interior-face protection. `GEO-11`'s plan reserves any tolerance change for a review with the numbers in hand, which these are. Whoever takes it must add a facet-tag assertion at the same time: the defect was invisible precisely because nothing gates the group |
| **Owned by** | `GEO-12` (commissioned by the 2026-08-06, 10:30 review — §9 item 2, with the tolerance decision taken); **discharged 2026-08-06, 13:30 implementer slot** — see the retirement note above |

### 13. ✅ RETIRED 2026-08-07 — `cylindrical_domain`'s classification margin was 4.50× its tolerance (`GEO-13`)

The tolerance was the *mesh size*: `abs(r_max - outer_radius) < resolution`. It
is now `0.01 × (outer_radius - inner_radius)` — a fraction of the radial gap,
so the margin is a ratio of geometry to geometry and no longer moves when a
caller coarsens the mesh. Measured across all four argument sets the repo calls
the generator with (`20260807T033127Z_GEO-13-probe.log`): the fraction window
where **both** sides of the `GEO-11` identity hold is `[1e-4, 0.05]`, and `0.01`
sits in the middle of it. At defaults the interior margin goes
**`4.499995×` → `99.99989×`** (floor `10×`) with the accepted side at
`1.111111e-04×` (ceiling `0.1×`), and the classification itself is **unchanged**
— still 3 of 6 surfaces, on every one of the four geometries.

The failure mode the entry named is reproduced in the probe before the fix: at
`resolution = 0.09` (the gap) the old predicate accepts **6 of 6** surfaces, the
inner cylinder swept whole into `outer_boundary`.

`tests/mesh/test_boundary_classification_margins.py` now **asserts** the
two-sided margin for this fixture instead of pinning it — the pin and its skip
are gone, and the file reads the fraction from the generator so the two cannot
drift apart. All four fixtures in that file are now live: **5 passed in 1.05 s**
(`20260807T033236Z_GEO-13-margins.log`, `-n 1`). Regression: whole `tests/mesh`
**36 passed, 1 skipped in 110.34 s** (`20260807T033250Z_GEO-13-mesh-regression.log`,
`-n 2`) — one skip fewer than the 35/2 on record, which is this fixture. Callers
green at `-n 2`: `4 passed, 1 skipped in 0.97 s`
(`20260807T033454Z_GEO-13-callers.log`; the skip is the complex-mode PEC test).

**Live precondition, new:** the tolerance scales with the gap, so a gap below
\~`1e-4` m stops clearing the `1.000e-07` gmsh OCC bounding-box padding by 10×.
Recorded at the use site in `io/mesh.py`; no caller is near it (smallest gap in
the repo is `0.07` m).

<details>
<summary>Original entry</summary>

| | |
|---|---|
| **Test** | None fails. Measured and pinned by `tests/mesh/test_boundary_classification_margins.py` (`GEO-11`), which skips the margin assertion after pinning the ratio |
| **Symptom** | The nearest surface the wall test *rejects* sits at residual `8.999990e-02` = **`4.499995 ×`** `tol`, below the `10×` floor `GEO-11` asserts. The accepted side is fine (`5.000e-06 × tol`) |
| **Cause** | Measured (`20260806T140325Z_GEO-11-probe.log`): the test is `abs(r_max - outer_radius) < resolution` with `tol = resolution = 0.02` at defaults. The inner cylinder's surface and end caps sit at `r_max = inner_radius = 0.01`, a residual of `0.09`. The margin is a ratio of geometry to *mesh size*, so it shrinks as either radius gap narrows or `resolution` coarsens — at `resolution ≥ 0.09` the inner cylinder would be swept into `outer_boundary` |
| **Live impact** | **None at defaults.** `tests/mesh/test_cylindrical_domain.py` passes; the classification is correct today, only under-separated |
| **Verified at** | `main` as of `2cad984` |
| **Fix, not taken in-slot** | Either decouple the tolerance from `resolution` (a geometric fraction of `outer_radius - inner_radius` is the natural choice) or document the sizing precondition. A per-fixture decision for a review, per the `GEO-11` plan |
| **Owned by** | `GEO-13` (filed 2026-08-06, 18:00 review, with the geometric-fraction fix and the un-skip of the margin assertion; was `GEO-11`, whose sweep found it); **discharged 2026-08-07, 22:30 implementer slot** — see the retirement note above |

</details>

### 7. ✅ RETIRED 2026-08-03 — birdcage mesh fails to generate (`GEO-9`, steps 1 + 2a + 2b)

All three tests are green and the whole of `tests/mesh` less known-issues 5 is
`20 passed, 1 skipped, 1 deselected in 42.15 s`, exit 0
(`20260803T200504Z_GEO-9-step2b-gate.log`, `-n 2`) — the CI command verbatim,
with the birdcage `--ignore` removed in the same commit. **Step 2b** replaced the
`occ.cut(..., removeTool=False)` with a single `occ.fragment` of the air box
against every tool, so the legs and rings that pierce each other by construction
are booleaned into conforming pieces instead of being meshed twice; the physical
groups are re-derived from the fragment out-map (26 volumes, 20 of them
conductor) and the port boxes get the 3-D groups they never had. Measured:
`V_mesh/V_box = 1.000000000000` and `Σ(tagged)/V_mesh = 1.000000000000`, both
gated to `1e-9`, and every port box exact to `1e-9` of `dx·dy·dz`. The
rank-local `set(np.unique(cell_tags.values))` was fixed to `global_cell_tag_set()`
in the same commit — it was real and it fired: at `-n 2` rank 0 reported P2/P3
missing and rank 1 reported P1/P4 missing on an otherwise correct mesh
(`20260803T200151Z_GEO-9-step2b-probe.log`).

The history below is kept because several commits and the `GEO-9` §7 entry refer
to it, and because the diagnosis — one poisoned generator hanging every later
mesh in the process — is the reusable part.

<details>
<summary>Original entry (steps 1 and 2a)</summary>

#### Birdcage mesh fails to generate *(coil+phantom half resolved 2026-08-03, `GEO-9` step 2a)*

| | |
|---|---|
| **Tests** | ~~`tests/mesh/test_birdcage_port_tags.py::test_birdcage_like_mesh_has_core_and_port_tags`~~ — the "still red" note here is **historical**; the overlapping-facets geometry was fixed by `GEO-9` step 2b and the test itself no longer exists: `OPS-17` step 2 (2026-08-17) removed it as finiteness-only, its mesh-side content having been subsumed by `test_birdcage_volumes_partition_the_box` in the same file<br>~~`tests/mesh/test_coil_phantom_mesh.py::test_coil_phantom_mesh_generates_required_tags_centered_preset`~~ **passes since 2026-08-03**<br>~~`tests/mesh/test_coil_phantom_mesh.py::test_coil_phantom_mesh_off_center_preset_moves_phantom_without_overlap`~~ **passes since 2026-08-03** |
| **Symptom** | `gmsh.py:2006: Exception: Invalid boundary mesh (overlapping facets) on surface 3 surface 49` (birdcage) and `dolfinx/io/gmshio.py:118: AssertionError` ×2 (coil+phantom) — the meshes never reach dolfinx |
| **Cause** | **Diagnosed 2026-08-03 by `GEO-9` step 1, and it is a single cause, not two.** `birdcage_port_domain` raises inside its `comm.rank == rank` block (the overlapping-facets error) and therefore never reaches its `gmsh.finalize()`. The process is left with gmsh initialised and mid-command, so the *next* generator in the same pytest process gets `Warning : Gmsh has aleady been initialized` and `Info : I'm busy! Ask me that later...`, every subsequent `occ` call is silently refused, and `model_to_mesh` reads the stale birdcage model — whose mixed element types are what `gmshio.py:118` asserts on. **The coil+phantom generator is innocent:** in a fresh process all three of its tests pass in 4.8 s (`20260803T033050Z_GEO-9-before.log`), and they fail again the moment the birdcage file runs first (`20260803T033119Z_GEO-9-order-probe.log`, 3 failed 2 passed in 3.47 s — the same two failures). Found 2026-08-01 by the `GEO-8` run, which ran `tests/mesh` as a regression sweep; the whole directory is in no CI job, which is why these were invisible. Both fixtures are untouched by `GEO-8` (it changed `two_torus_domain` only). |
| **Verified pre-existing at** | `63c94f2` — the `GEO-8` diff touches neither generator; log `20260801T004839Z_GEO-8-unrelated-failures.log`, 3 failed 2 passed in 3.5 s at `-n 2` |
| **Note** | ~~Likely the same overlapping-geometry family `GEO-8` just fixed for `two_torus_domain`: `coil_phantom_domain` and the birdcage fixture should be audited for missing `occ.fragment`.~~ **Half wrong — corrected 2026-08-02, 18:00 review, by code reading (not by execution).** `coil_phantom_domain` **already fragments** (`io/mesh.py:1616`), so its cause is downstream of fragment; the visible fragility is the group re-derivation at `io/mesh.py:1622-1634`, which assumes fragment returns exactly four volumes and leaves any extra piece with no physical group — which is what `gmshio.py:118` asserts on. The **birdcage** does match the family, differently: it uses `occ.cut(..., removeTool=False)` (`io/mesh.py:1970`), so the conductors, phantom and port boxes are never booleaned against *each other* and the port boxes overlap the legs by construction — hence "overlapping facets". The birdcage's port boxes also receive no 3-D physical group (`io/mesh.py:1985-1988`). |
| **Owned by** | `GEO-9` (created 2026-08-02): step 1 coil+phantom, step 2 birdcage. Two of the four §10 Target criteria route through these two fixtures. **Step 1 landed 2026-08-03** — the coil+phantom half is now gated by `tests/mesh/test_coil_phantom_conforming.py` (volume-partition identity, `20260803T033659Z_GEO-9-step1-gate.log`, 8 passed 1 skipped in 22.25 s) and the generator raises with the volume count and per-volume masses if fragment ever does return other than four grouped volumes. **All three tests above still fail in a shared process** and will until step 2 fixes the birdcage: `20260803T033733Z_GEO-9-order-probe-after.log` shows the new guards do *not* fire (gmsh is already busy, so they never execute), still `gmshio.py:118`. **Step 2's first action is the cheap half of this — wrap the birdcage's rank-0 block in `try/finally: gmsh.finalize()`**, which stops one broken generator from poisoning every later mesh in the process, independently of fixing the geometry. Split out as **`GEO-9` step 2a** at the 2026-08-03 03:00 review, with step 2b holding the `occ.fragment` geometry rewrite. |
| **A poisoned process hangs — added 2026-08-03, 03:00 review** | Both `GEO-9` order probes report pytest finishing (3.47 s and 3.29 s) but the **harness** exits **124 at the 180 s ceiling**, with `Loguru caught a signal: SIGTERM` (`docs/testing/test-results.md:136,139`; logs `20260803T033119Z_GEO-9-order-probe.log`, `…033733Z_GEO-9-order-probe-after.log`). The step-1 prose quotes only the pytest wall time and does not say this. Two consequences: `tests/mesh` cannot enter CI with the birdcage in it (`OPS-11`) because the job would burn its whole `timeout-minutes` rather than fail fast; and `GEO-9` step 2a gets a sharper anchor than "the tests pass" — harness exit **124 at 180 s → 0 in seconds**. If a `try/finally` does not fix it, suspect an MPI collective the raising rank never reaches, not gmsh state. |
| **Excluded from CI, and the reason has changed — `OPS-11`, 2026-08-03** | The `validation` job's `Mesh generation suite` step `--ignore`s `test_birdcage_port_tags.py`. The exclusion no longer rests on the hang or the budget: post-2a the file fails **promptly** (the unexcluded-directory control `20260803T170132Z_OPS-11-fullsweep.log` is 2 failed 18 passed 1 skipped in 31.85 s, harness exit 1 in 33 s, where the pre-2a probes exited 124 at the 180 s ceiling). It is ignored because it is deliberately red until `GEO-9` step 2b, and a permanently-red test in CI hides regressions behind an expected failure. **Remove the `--ignore` in the commit that fixes the geometry.** The same control also shows the coil+phantom tests now passing *with the birdcage in the same process* — the step-2a poisoning fix holding under the condition that used to break it. |
| **✅ Two-thirds resolved 2026-08-03 by `GEO-9` step 2a — and the hang was *both* causes, not one** |<!-- retired: see the summary above --> The review's warning was right: `try/finally` alone would not have fixed the hang. Rank 0 raised and skipped the collective `gmshio.model_to_mesh`, so rank 1 blocked in it forever — that is the exit 124, and gmsh contamination is a second, independent defect. `birdcage_port_domain` now builds its model in `_build_birdcage_port_model`, calls `gmsh.finalize()` (guarded by `gmsh.isInitialized()`) when that raises, and `comm.bcast`es the failure so **every** rank raises before any enters `model_to_mesh`. The birdcage still fails loudly with the original message. **Before** (`20260803T123116Z_GEO-9-step2a-before.log`, re-run at the working commit, not quoted): 5 failed 2 passed in 3.16 s of pytest, harness **exit 124 at 180 s**. **After** (`20260803T123549Z_GEO-9-step2a-after.log`, same command): **1 failed 6 passed in 12.10 s, harness exit 1 at 13 s** — only the birdcage test itself. Gated by `tests/mesh/test_birdcage_finalize_isolation.py` (`20260803T123657Z_GEO-9-step2a-gate.log`, 1 passed in 5.30 s, **exit 0**), which runs the two generators in the poisoning order in one process and asserts `V_mesh/V_box = 1.000000000000` and `Σ(tagged)/V_mesh = 1.000000000000` to `1e-9` afterwards. `tests/mesh` less the birdcage file is now 17 passed 1 skipped, 1 failed in 28.46 s (`20260803T123714Z_GEO-9-step2a-sweep.log`) — the single failure is **entry 5**, unrelated. |

</details>

### ~~`check_example_doc_references.py` freshness-gates only 5 of 27 examples — every `stale=24, none of them mine` line is not an all-clear~~ — RESOLVED 2026-08-24 (`EX-29`; found by the weekly review 2026-08-23)

**Found by the 2026-08-23 weekly examples audit (read-only, not a test
failure).** The checker's `--output-dir` defaults to the single repo-root
`paraview_output/` (`scripts/testing/check_example_doc_references.py:241-242`),
and its `in_tree_artifacts` exemption (`:276-298`) treats any referenced
artifact whose basename also appears anywhere under `examples/` as
"committed next to its own case — existence is enough". Only
`mag` 1/2/4/5/6 and `mri:1` write to the repo-root directory; the other
**22 of 27** runnable examples write `Path(__file__).parent/"paraview_output"`,
so they get an existence-only check and are **never** freshness-gated. The
exemption's premise is false: `.gitignore` ignores `paraview_output/` at
every depth and `git ls-files examples/ | grep paraview_output` returns
nothing — no artifact under `examples/` is committed evidence. That is why
the stale count has read exactly 24 on every docrefs log since `OPS-19`
and always names the same set.

**Consequence.** On 2026-08-23 the artifacts the checker cannot see are
10–17 days old for 13 examples (`mesh:1` ≈17 d, `mesh:2` ≈16 d, `mri:2`
≈15 d, `mat:1` / `th:1`–`th:4` / `ans:1` ≈14 d, `th:5` ≈13 d, `ports:1` /
`th:6` ≈10 d), all predating `OPS-17`'s test replacement; no EX chunk's
"none of them mine" reading ever covered them. A second, smaller finding:
`examples/magnetostatics/paraview_output/` (2026-08-03/04 `circular_loop_*`
files) is an orphan — `02_circular_loop.py` has written to the repo root
since `EX-17`, and nothing regenerates that directory.

**RESOLVED 2026-08-24 by `EX-29`** (`tests/unit/test_doc_reference_exit_codes.py`,
15 tests, green twice in-slot; logs `20260824T110512Z_EX-29-unit.log` /
`20260824T110540Z_EX-29-unit-run2.log`). The checker now resolves each
referenced artifact in the citing guide's own `paraview_output/` first and the
shared `--output-dir` second, and exempts only paths `git ls-files` reports as
tracked. Measured on the same tree, same slot: **pre-fix `stale=24`**
(`20260824T110150Z_EX-29-prefix-control.log`) → **post-fix `stale=55`**
(`20260824T110531Z_EX-29-census.log`, `dead=0 guide=0 exit=2`), with
**32 of 58** resolved artifact references sitting outside the repo-root
directory — the set the old exemption hid. The orphaned
`examples/magnetostatics/paraview_output/` is deleted. `EX-30` refreshes the
now-visible stale set and is re-sized from 55 by a review.

Two corrections to this entry's own text, both measured: (a) the tracked set
is **not** empty — `git ls-files` reports three committed artifacts under
`examples/` (`ansys_benchmarks/*/metrics.json` ×2,
`magnetostatics/straight_wire_validation.png`), which are exactly the
"committed next to its own case" exemption the rule was written for and are
pinned by path in the test; (b) `git` inside the `fem-em-solver` container
fails with `fatal: detected dubious ownership in repository at '/workspace'`
(root over a host-owned bind mount), which would have silently emptied the
exemption and reported the two tracked `metrics.json` as **dead references** —
the checker passes `-c safe.directory=` for both the repo root and the docs
root. Any future in-container `git` call needs the same.

### ~~`test_coil_phantom_magnetostatics` fails in the complex build on a cold FFCx cache: `ComplexComparisonError`~~ **RESOLVED 2026-08-19 (`OPS-20`, 06:00 implementer slot)** (`OPS-17` step 3 leg (b1), 2026-08-18)

> **Resolved 2026-08-19, `OPS-20` step 1 — fixed, not marked.** The Cause row
> below is wrong and the 03:00 review's re-pointing was right: the `max` never
> entered through a DolfinX/UFL helper. It came from the test's **imported**
> drive callable — this file uses
> `tests/validation/test_circular_loop.azimuthal_current_density`, whose
> `ufl.max_value` `OPS-22` had already replaced with a regularise-inside-the-
> `sqrt` form at the 04:30 slot. So the commissioned defect was **already dead
> on arrival** and a free grep of the import line localized it; no cold-cache
> window was spent. What remained was the *second layer* `OPS-22` warned
> `OPS-20` to expect: with the form compiling, the complex run reached the
> print block and died at
> `ValueError: Unknown format code '%' for object of type 'complex'`
> (`test_coil_phantom_magnetostatics.py:145`) — `evaluate_vector_field_parallel`
> returns the complex scalar type although this magnetostatic solution is
> real-valued. Note the **rank split** that message produces: only rank 0
> executes the print block, so the diagnosis run read `1 failed` on rank 0 and
> `1 passed` on rank 1 in the same command. Fixed with the `OPS-22` idiom —
> assert `max|Im B_z| ≤ 1e-12·max|B_z|`, then compare on `np.real` (a new
> complex-mode assertion; exactly zero and a no-op in real mode). The
> non-collective ≈300 s exit hang is gone with the raise: every run in this
> slot returned a footer in ≤ 8 s.
>
> **Numbers** (`-n 2`, all four runs footered, exit 0 except the diagnosis):
> real-mode control before any edit `20260819T110051Z` — 1 passed / 5.81 s,
> L2 **17.1233%**, the `OPS-17` step-2 record to the digit; complex diagnosis
> `20260819T110111Z` — 1 failed (rank 0) / 6.19 s, the `ValueError` above with
> a user frame under `--tb=long`; complex after the fix `20260819T110144Z` —
> **1 passed / 5.11 s, L2 17.1233%**, both ranks identical, i.e. the complex
> build passes the *same* 30% gate at the *same* digits; real-mode re-run
> `20260819T110156Z` — 1 passed / 3.36 s, **17.1233%** unmoved, so the fix
> moves no real-mode digit. No `@real_only` marker anywhere, so the complex
> collect count is **unchanged at 49** and `OPS-17`'s bookkeeping does not
> move. Stub sweep `find /root/.cache/fenics -name '*.c' -size 0` clean before
> and after.
>
> **For `OPS-17` leg (b1):** the coil-phantom exclusion recorded when that leg
> closed is now discharged in principle — the file is green in the complex
> build. A whole-`tests/solver` complex batch attempted as confirmation in this
> slot **timed out at 89%** (`20260819T110220Z`, exit 124, 481 s) and is
> therefore **uncounted** — no footer, no count claim. The coil-phantom test
> itself is visible PASSED on both ranks at 10% in that log. Separate finding
> for the review: complex `tests/solver` fit 111.22 s warm on 2026-08-18 and no
> longer fits a 480 s window, which points at cold forms added since (`POST-5`
> step 2 is the candidate), not at this fix.

**Verified at `93fc531`, 07:30 implementer slot,
`20260818T124742Z_OPS-17-step3d-coilphantom-complex-cleancache.log`**, `-n 2`,
complex build + `FEM_EM_REQUIRE_COMPLEX=1`, `/root/.cache/fenics` removed
immediately before the run.

| | |
|---|---|
| **Test** | `tests/solver/test_coil_phantom_magnetostatics.py::test_coil_phantom_magnetostatics_matches_the_two_loop_closed_form` |
| **Symptom** | `ufl/algorithms/comparison_checker.py:66: ufl.algorithms.comparison_checker.ComplexComparisonError: You can't compare complex numbers with max.` — **1 failed in 5.58 s**, raised during form compilation, before any assertion runs. |
| **Mode** | Complex build **only**. Real mode is green and gates 17.1233% L2 against a 30% band (`OPS-17` step 2). Nothing about the gated quantity is complex-valued; this is a UFL-level rejection of a `max`-style comparison in a form that is assembled with a complex scalar type. |
| **Cause** | **Not diagnosed.** The run used `--tb=line`, which printed only the UFL frame and no user frame, so the offending expression is not localized. `grep` for `max_value`/`min_value`/`conditional(` across `src/` finds exactly one hit (`src/fem_em_solver/post/sar.py:286`), which this test does not exercise — so the comparison most likely enters through a DolfinX/UFL helper (a cell-size or clamp expression) rather than a literal `ufl.max_value` call. **One command settles it:** re-run this file alone on a cleared cache with `--tb=long`. |
| **Not** | Not the cache artifact the previous entry claimed, and not a regression from any recent chunk — no completed complex leg had ever reached this file before today, so it has no known-green complex history. |
| **Scope** | `OPS-17` is test-hygiene bookkeeping and deliberately did not fix it. Whoever fixes it should record whether the complex build ever needs this magnetostatic path at all — if not, an explicit `@real_only` marker is a legitimate disposition and is cheaper than making the form complex-safe. |
| **Resolves with** | `OPS-20` (commissioned 2026-08-18 10:30 review): one `--tb=long` cold-cache diagnosis, then fix-or-mark, real-mode record 17.1233% re-asserted unmoved either way. **Re-pointed 2026-08-19 03:00 review:** the Cause row's "DolfinX/UFL helper" hypothesis is disfavoured — every other instance of this error class was a fixture-side `max`/ordering predicate in the test's own drive callable (see the `OPS-22` entry above); start the diagnosis by grepping this test's drive construction. |

### ✅ RETIRED 2026-09-01 — `test_coil_loading_degree2.py` no longer returns inside its 570 s ceiling (exit 124 at 571 s, **twice** on 2026-08-31)

**Retired by `TH-13` step 3a‴, 22:30 implementer slot 2026-09-01 — the case
was shrunk, the ceiling was not raised, and all three windows returned
footers.** The degree-2 pair now lives in
`tests/validation/test_coil_loading_degree2_pair.py` and runs one σ-half per
process (`TH12_DEGREE2_HALF=loaded|free`); the original module defaults to
`probe`. Three runs, all `-n 8`, complex build:
`20260901T033335Z_TH-13-step3a3-original-probe.log` (**8 passed / 1 skipped,
exit 0, 49 s**, `timeout -k 30 180`);
`20260901T033434Z_TH-13-step3a3-degree2-loaded.log` (**1 failed / 6 passed /
2 skipped, exit 1, 374 s**) and
`20260901T034059Z_TH-13-step3a3-degree2-free.log` (**1 failed / 6 passed /
2 skipped, exit 1, 405 s**), both `timeout -k 30 600` — 1.5–1.6× margin inside
the ceiling and inside the 660 s Bash window. The pre-registered success
condition (a footered `exit 1` with exactly one degree-2 identity red per half)
is met exactly, so the owed observation is discharged: see the step-3a‴ row on
the entry above. The 138 490-cell mesh anchor, the +1.5838% degree-1 ΔR control
(+0.00039 pp off record, floor 0.01 pp) and the cost probe's under-cap verdict
(162 558 → 881 476 DOFs, 5.42×; 52.35 GiB against the 102.40 GiB threshold) were
re-observed in **all three** windows. `full` mode is retained for interactive
use only and is documented as such; nothing else about it was changed. The
original entry follows unedited.

**`TH-13` step-3a‴ reading (2026-09-01, the two half-mode windows above) — THE
REDS ARE OBSERVED AGAIN, FOR THE FIRST TIME SINCE 2026-08-18, AND THEY ARE
STILL RED.** `loaded`: `Im Z` reaction −2.323123e+03 Ω vs energy −2.323123e+03 Ω,
relative **3.8990e-09** against the unloosened 1e-9 (`W_m` 3.1357e-08 J, `W_e`
7.8593e-06 J); `P_loss` +1.3543068e-01 W. `free`: −2.322561e+03 Ω vs
−2.322561e+03 Ω, relative **3.7235e-09** (`W_m` 3.3258e-08 J, `W_e`
7.8594e-06 J); `P_loss` **+0.0000000e+00 W** exactly, the σ = 0 control green at
second order. Both sit inside the 08-18 record's band (4.5931e-09 / 3.0030e-09)
on a mesh whose cell count, degree-1 ΔR control and degree-1 identity residuals
were re-observed in the same process, so the residuals are an order effect and
not a fixture or reduction drift — the 08-18 diagnosis is unmoved and the `W_e`
explosion (7.86e-06 J against the degree-1 rows' green 1e-14-level residuals) is
reproduced. Nothing was loosened, no band moved, no coil number moved. **This
entry stays open**; what closed is only the "cannot be run at all" entry below.
The three cross-half degree-2 checks (the drive control, the ΔZ signs, and the
ΔR reading against step 4's bracket) are unobservable one half at a time and
stay as 2026-08-18 recorded them — the two tests `pytest.skip` with that reason.

### ✅ RETIRED 2026-08-11 — "unexplained" mid-command termination of the logging harness was the background-and-end-turn trap (2026-08-08, 15:00 and 19:30 implementer slots)

**Retired by the 2026-08-11 10:30 review — cause named, with wrapper-log
evidence; nothing host-side ever killed anything.** Both slots' automation
wrapper logs end with the session announcing it is waiting on a
*backgrounded* harness run, then `exit=0` at the exact minute the harness
log stops: `logs/automation/20260808T200001Z_implementer.log` ends "Stage 2
is running … I'll report when the monitor fires", `exit=0` 20:16:37Z
(death ≈20:15Z); `logs/automation/20260809T003001Z_implementer.log` ends
"Waiting on the background solve … I'll report when it lands", `exit=0`
00:33:04Z (the log's last flushed write, to the minute). A headless
`claude -p` session that ends its turn exits the CLI and SIGKILLs its
process tree — harness included: footerless log, no journal, dirty tree.
This is the same mechanism, established the same way, as the three `MAT-6`
step 7 deaths root-caused by the 2026-08-11 03:00 review (attempts.md
2026-08-11T08:00Z). It explains every anomaly this entry catalogued: the
6.7× spread in time-to-death (it is when the session chose to end its turn,
not a resource ceiling), the never-restarted container, the missing
attempts.md entries (the session died with its journal unwritten), and the
stage-uncorrelated kills (MESH_ONLY ran clean *foreground*). The fix is the
foreground recipe already landed in implementer-run.md and the §9 rubric
trap list (03:00 review commit `bc86367`); this commit names the cause,
which is the condition this entry set for leaving. **The standing
instruction below is lifted**: `MAG-13` step 2 is retryable under the
foreground recipe (re-queued §9 item 3, 2026-08-11 10:30), and the
operator's host-observables ask (`dmesg`, `journalctl`, WSL2 reclaim) is
withdrawn from Waiting-on-you. Original record retained below.

**Observed twice — the standing instruction below has fired; `MAG-13` step 2
is escalated, not retryable.** The second occurrence is recorded after the
first-occurrence text; the two are compared at the end of this entry.

Observed first during `MAG-13` step 2's stage-2 solve
(`timeout 1200 mpiexec -n 4 python3 scripts/probes/mag13_step2_probe.py`,
started 20:04:51Z), the harness died ≈660 s in — the log
(`20260808T200451Z_MAG-13-step2-solve-n4.log`) ends mid-Netgen with **no
`## Exit` block** and no `test-results.md` row. That is not the command
exiting: `run_and_log.sh` writes an Exit block even on non-zero status (it
did so for the same slot's stage-1 run forty minutes earlier), so the
harness process itself was killed from outside. The kill fell well inside
the command's own `timeout 1200` and well before the slot's 65-minute hard
kill (21:05Z). Not a container-cgroup OOM signature either — those
manifest as the *command* dying with signal 9 / exit 137 and an Exit block
written (see `MAT-6` step 6's two on-record examples). The session died
with its harness: no attempts.md entry, nothing committed or parked
(journaled by the 16:30 slot; artifacts landed by the 18:00 review,
`8b8a706`).

**Second occurrence — 2026-08-08, 19:30 slot.** The same command on the same
probe (`timeout 1200 mpiexec -n 4 python3 scripts/probes/mag13_step2_probe.py`,
started 00:31:25Z) died again, log
`20260809T003125Z_MAG-13-step2-solve-n4-cap16G.log`: 43 437 B, 627 lines, no
`## Exit` block, no `test-results.md` row, stopping mid-Netgen volume
optimisation (`Total badness = 1.36536e+06`). Last flushed write 00:33:04Z —
**≈ 99 s after start**, versus ≈660 s for the first occurrence. That slot also
left no attempts.md entry; its tree was journaled by the 21:00 slot and parked
by the 22:30 slot on `recovered/20260809T033023Z`. The log's one durable
measurement: `CGROUP_MEMORY_MAX=17179869184` (16 G) confirmed at the kernel
before any solve, so `MAT-6` step 6's cap no longer rests on a compose-file
read.

**Common factors, and the ones ruled out (22:30 slot, 2026-08-09):**

| | 15:00 slot | 19:30 slot |
|---|---|---|
| log | `…200451Z_…-n4.log` | `…003125Z_…-cap16G.log` |
| died after | ≈660 s | ≈99 s (flushed-output bound) |
| stage reached | `Done optimizing mesh (Wall 149.77s)` | mid-volume-optimisation |

Shared: same probe, same `-n 4`, same `run_and_log.sh` → `docker compose exec`
path, both truncated with no Exit block, neither at its `timeout`, neither at
its session hard kill, neither with a kernel OOM signature. **Ruled out — the
container did not restart.** `docker inspect` at 03:30Z reports
`StartedAt = 2026-08-08T20:00:21Z`, `RestartCount = 0`; the container has been
continuously Up across **both** deaths (20:15Z and 00:33Z). So the cause is not
a container/cgroup restart, and it is not inside Docker's lifecycle — the
host-side process tree is being killed. The 6.7× spread in time-to-death also
argues against a deterministic per-run resource ceiling, which would land at a
repeatable point on a fixed fixture.

**MESH_ONLY discriminator executed — 2026-08-09, 07:30 slot; the physics is
exonerated and the kill is *not* stage-correlated.** The queued next step ran
(`20260809T123053Z_MAG-13-step2-meshonly-diag.log`, exit 0, 188 s harness-wall,
668 lines, `## Exit` block present, `test-results.md` row written): the same
probe, same `-n 4`, same harness path, `MAG13_STEP2_MESH_ONLY=1`, FFCx cache
cleared first. It **reproduced the mesh rung exactly** — **1 097 873 cells**,
equal to the 2026-08-08 record digit for digit, in 185.7 s (record 192.7 s,
−3.6%); the log is structurally identical to the record's, both with their two
`Done optimizing mesh` lines at the same line numbers (486 / 663; fine
volume-optimisation wall 142.4 s here vs 147.8 s on record). Container state
before and after: `StartedAt = 2026-08-08T20:00:21Z`, `RestartCount = 0`, Up
17 h — unchanged across this run and both deaths.

**What that rules out, and the part that contradicts the pre-decided reading.**
Branch (a) of the §7 plan fired literally — MESH_ONLY completes — but its
inference ("the kill is specific to the longer/heavier solve stage") does **not
survive this run's own comparison**: the 19:30 death stopped mid-Netgen *volume
optimisation* of the fine mesh (`Total badness = 1.36536e+06`, before any
`Done optimizing mesh (Wall 14x s)` line and before any solve), i.e. inside the
very phase MESH_ONLY has now completed **twice** at the same rank count and
resolution. So no stage owns the kill: one death is in the mesh phase, one past
it in the solve, and the mesh phase runs clean on demand. Combined with the
6.7× spread in time-to-death and the never-restarted container, the surviving
hypothesis is a **non-deterministic host-side kill of the process tree,
uncorrelated with the computation** — WSL2 memory reclaim / a host session
supervisor are the live candidates, and none of it is observable from inside
the container.

**Standing instruction (unchanged in force):** a run that finds its log
truncated this way should treat the measurement as *unobserved*, not failed.
Do **not** spend a slot on the `MAG-13` step 2 solve — the stage-2 attempt
stays blocked pending a review, and the diagnostic budget for it is spent
(three data points, physics exonerated). What is left needs the **human
operator**: host-side observables at the timestamps of the two deaths
(20:15Z and 00:33Z on 2026-08-08/09) — `dmesg -T`, `journalctl -k`, WSL2
`vmmem` reclaim, any host cron/session supervisor. On the dashboard's
Waiting-on-you as of this slot. Cause unknown; this entry leaves only with a
commit that names and fixes it.

### ✅ RETIRED 2026-08-09 — "rank-dependent DG0 centerline sample" was the probe's own bug (`MAG-6` step 4, second pass)

**There is no rank-safety defect.** The 00:00 slot completed step 4's missing
rungs and the claimed defect is **refuted with its own signature measured**:
the probe's `instrumented_eval` inflated `|B|` by exactly **√3** at any point
whose claiming rank held *only that one point*. `Function.eval` squeezes its
return to shape `(3,)` for a single point, so `rank_vals[k]` was the scalar
x-component and `values[i] = rank_vals[k]` broadcast it across all three
components. Measured ratio at the two affected points
(`20260809T050838Z_MAG-6.log`, `WRITECHECK` lines):
`4.852607687905e-07 / 2.801654354883e-07 = 1.7320508` and
`2.853753669222e-07 / 1.647615449126e-07 = 1.7320508` — √3 to 8 digits. The
production path `post/evaluation.py::evaluate_vector_field_parallel` is
**immune by construction**: it assigns `values[rank_indices] = rank_values`, and
numpy broadcasts a `(3,)` row into a `(1, 3)` slice correctly. Nothing under
`src/` was ever wrong; the one-line fix is in the probe
(`scripts/probes/mag6_step4_probe.py`, `.reshape(-1, 3)`).

**The rank-invariance identity holds.** On the gate's own evaluation path at
the validated `gauge_penalty=1.0`, the centerline jump ratio reads
**0.251272 / 0.250416 / 0.250453** at `-n 1 / -n 2 / -n 4` — a three-way spread
of **0.341%**, against the ≤ 10% band, with the mirror-symmetry control on the
same solves at 0.311226 / 0.311166 / 0.311157 (**0.022%**). Logs
`20260809T050259Z` (`-n 1`, 152 s), `20260809T050621Z` (`-n 2`, 10 s),
`20260809T050202Z` (`-n 4`, 12 s), all `_MAG-6.log`.

**So step 3's 88% scatter is gauge contamination after all** — the mechanism
the first pass believed it had refuted. It "refuted" it by comparing 0.250406
at `-n 2` against 0.328496 at `-n 4`, but that `-n 4` number was a probe
call-1 value carrying the √3 bug; the same run's library-path value is
0.250417. At the sub-floor `gauge_penalty=1e-3` the gate fixture uses, the
scatter is real; at 1.0 it is 0.341%. **A fix chunk on the DG0 evaluation path
is not needed and must not be scoped.** Whether the *gate fixture* should stop
solving below the validated gauge floor is a separate, live question for a
review — `MAG-6` stays ✅ either way, passing its untouched 0.60 bound at every
rank count measured. Confirming run after the probe fix:
`20260809T050930Z_MAG-6.log` — all four evaluations in one process bitwise
identical, zero `WRITECHECK` DIFF, metric 0.250457 at `-n 4`.

*Original entry, retained for the record:*

**One cell's `B` value depends on the rank count, on a mesh and an owning cell
that are provably identical.** Diagnosing `MAG-6` step 3's 88% rank scatter in
the centerline jump-ratio metric, step 4 instrumented the point evaluation
(`scripts/probes/mag6_step4_probe.py`; no `src/` change) and eliminated both
mechanisms step 3 had proposed:

- **Not mesh noise.** `MESH_FINGERPRINT` (global cell count plus reduced
  midpoint moments) is `cells=55784 m1=-4.9768680987…e+00 m2=7.977798997317e+02`
  in **every** run — `-n 2`, `-n 4`, and a fixed-rank repeat — agreeing to 12
  significant digits (the last digit of `m1` moves with reduction order alone).
  The mesh is reproducible run to run; step 3's 6.8% "mesh drift" attribution
  does not survive.
- **Not partition-owned sampling.** `CENTERLINE_MULTICLAIM = 0/9` and
  `CENTERLINE_MULTICELL = 0/9` at both rank counts: every centerline point is
  claimed by exactly one rank and collides with exactly one cell, so the
  `links[0]` / rank-order-overwrite ambiguity in
  `post/evaluation.py::evaluate_vector_field_parallel` never fires here. The
  **chosen cell midpoints are identical across `-n 2` and `-n 4` for all nine
  points**, printed to 9 decimals.
- **Not gauge contamination.** The fixture solves at `gauge_penalty=1e-3`,
  below the validated floor of 1, and raises `GaugeContaminationWarning`. Re-run
  at `gauge_penalty=1.0` the scatter persists: jump ratio 0.250406 at `-n 2`
  vs 0.328496 at `-n 4` (31%).

**What is left.** At `gauge_penalty=1.0`, eight of the nine centerline points
are rank-invariant to ≈5 significant digits between `-n 2` and `-n 4`. The
entire metric spread is set by **one point, i=1 at z = -0.0225 m**, whose
sampled `|B|` reads `2.813455e-07` at `-n 2` and `4.852531e-07` at `-n 4` —
**72% apart, in the same cell** (midpoint
`(-1.204260909e-03, +4.174143551e-03, -2.041163735e-02)` in both). Same mesh,
same cell, same field definition, different rank count. That is a rank-safety
defect in the solve/interpolation path, not in the metric's reduction and not
in the sampling.

**A second, possibly related signal, not yet diagnosed:** the probe evaluates
the same `Function` at the same points twice in one run — once instrumented,
once through the library — and compares them exactly. They agree at `-n 2`
(`LIB_AGREES_WITH_INSTRUMENTED = True`) and **disagree at `-n 4`** (`False`)
in the same process, on the same unchanged `b_dg0`. Two identical evaluations
inside one run should be bitwise equal. Whether this is the same defect or an
independent one is unmeasured — the probe prints only the boolean, not the
magnitude of the difference.

**Impact and scope.** No gate is affected today: `MAG-6` stays ✅, the
centerline jump ratio passes its untouched 0.60 bound at every rank count
measured (0.250406–0.328496), and the ≤ 10% rank-stability claim belongs to the
mirror-symmetry metric alone — which on these same solves reads
0.306591 / 0.309126 / 0.310501 / 0.311161 / 0.311162, a **0.15% spread**, and is
therefore the in-fixture control showing the defect is localised to the
centerline sample rather than global to the solve. No fix is attempted here:
step 4 is diagnosis-only and a fix is a review-scoped chunk, not an in-slot
improvisation. Logs: `20260809T033322Z`, `…033350Z`, `…033403Z`, `…033514Z`,
`…033555Z`, `…033608Z`, all `_MAG-6.log`.

### ✅ RETIRED 2026-08-04 — reaction Z-matrix diagonal is negative where it must be inductive

**Fixed by `PORT-1` step 2f**: `TimeHarmonicSolver.solve()` now drives with the
CG1-weakly-solenoidal part of the prescribed current by default
(`project_source=True`, helper
`src/fem_em_solver/core/source_projection.py`), and the diagonal of
`test_port_reaction_impedance.py` is **gated**, not printed. Gate
`20260804T111102Z_PORT-1-step2f-gate.log`, 12 passed 1 deselected in 58.9 s at
`-n 2`, on the same fixture every measurement below was taken on:

| quantity | production path, projected | this entry's unprojected number |
|---|---|---|
| `Im Z₁₁`, reaction / energy routes | **`+7.437243e+00 Ω`** (both) | `−4.108550e+01 Ω` |
| `Im Z₂₂`, reaction / energy routes | **`+7.436633e+00 Ω`** (both) | `−4.092413e+01 Ω` |
| ratio to Grover `ωL = 6.818343 Ω` | **1.090770 / 1.090680** | −6.03 |
| complex-power identity residual | `4.0412e-11` / `9.1813e-11` (gated `< 1e-9`) | `1.8128e-10` |
| driven current | `I′ = 0.969001 A` | `I = 0.969009 A` |

The three gates are `test_projected_port_diagonal_is_inductive` (sign, a
priori, both ports, both routes), `..._satisfies_the_complex_power_identity`
(bookkeeping, `< 1e-9`) and `..._matches_grover` (the independent physics
anchor, band `(1.042, 1.140)` carried over from step 2e's measurement, which
the production path reproduced to 2e-5). The number reproduced step 2e's
hand-rolled `+7.437243e+00 Ω` to all seven printed figures — same physics,
now on the path callers actually use.

Nothing was widened to retire this: the diagonal moved from *ungated* to
gated, and the three files that pin the unprojected numbers (steps 2b, 2d, 2e)
now pass `project_source=False` explicitly and reproduce them unchanged
(`20260804T111221Z_PORT-1-step2f-regress-diagnosis.log`,
`20260804T111607Z_PORT-1-step2f-regress-remainder.log`). The original entry
follows, unedited, because the diagnosis chain in it is the reason the fix is
believable.

---

Found 2026-08-02 by `PORT-1` step 1; **diagnosed 2026-08-03 by step 2b to the
electric energy — see the update at the end of this entry — and still not
fixed.** The diagonal remains ungated. Recorded here because the number is wrong
in a way a later run would otherwise re-discover from scratch.

On the two-torus air fixture (a = 0.04 m, r_wire = 0.005 m, d = 0.04 m,
f = 10 MHz), the reaction integral `Z_i1 = −(1/(I₁Iᵢ))∫E₁·Jᵢ dV` returns
`Im Z₁₁ ≈ −40.9 Ω` (`-4.069329e+01j` at
`20260802T183226Z_PORT-1-step1-solve008.log:442`, and −41.09 / −40.97 at the two
boxsens configurations). A lossless loop must be inductive, `+ωL`; a Grover
estimate `ωL ≈ μ₀ωa(ln(8a/r_wire) − 2) ≈ 6.8 Ω` — **hand-evaluated, in no log**.
So the diagonal is wrong in sign and ≈6× in magnitude, while the *off*-diagonal
on the same solve is right in sign and within 5–9% of its closed form. That
contrast points at the self-term (the source's own singular field inside the
driven wire entering `∫E·J` over the source region) rather than at a global
convention error, but nothing has been measured to confirm it.

**Consequence:** no input impedance and no `S₁₁` derived from this path means
anything yet. `PORT-1` step 2 therefore leaves the diagonal **ungated** —
deliberately, and *not* by widening a bound. `PORT-1` step 2b (§7) owns the
diagnosis, and its anchor is the complex-power identity
`Im Z₁₁ = 4ω(W_m − W_e)/|I₁|²` as an independent second derivation.

**Update 2026-08-03 — step 2b executed; the reaction integral is exonerated and
the anomaly is localised to the electric energy.**
`tests/validation/test_port_self_impedance_energy.py`, 3 passed in 43.5 s at
`-n 2` (log `20260803T050252Z_PORT-1-step2b-gate.log`), one mesh at padding
0.08 / h_far 0.03 and one solve:

| quantity | value |
|---|---|
| `Im Z₁₁`, reaction integral | `−4.108550e+01 Ω` |
| `Im Z₁₁`, complex-power route `4ω(W_m−W_e)/I²` | `−4.108550e+01 Ω` |
| relative disagreement | **1.8128e-10** (gated `< 1e-9`) |
| `4ωW_m/I²` (inductive part) | `+7.437 Ω` vs Grover `ωL = 6.818 Ω`, **ratio 1.0908** |
| `4ωW_e/I²` (capacitive part) | `+48.52 Ω` |
| `W_e/W_m` | 6.524 |

So the guess above — that the self-term of `∫E·J` is the bug — is **wrong**. The
two derivations agree to 1.8e-10, and the magnetic half is the physical loop
inductance to 9.1% of Grover, which is within what the PEC box at this padding
can plausibly account for. The whole of `−40.9 = 7.44 − 48.52` is an *electric*
energy excess in the solved field.

**Leading hypothesis, not yet measured:** low-frequency breakdown of the
curl-curl formulation. At ω → 0 the operator acts on the gradient subspace as
`−k₀²ε_c`, so any residual non-solenoidal component of the discretised impressed
current — the analytic azimuthal `J` is exactly divergence-free and tangent to
the torus surface, but the *faceted* meshed boundary is only approximately so —
is amplified into a spurious electrostatic field that contributes to `W_e` and to
nothing else.

**Correction 2026-08-03 (03:00 review): the ω-sweep named here as "the
discriminating measurement" does not discriminate.** Restricting the solved
equation to the gradient subspace gives `E_g = jωμ₀J_g/k₀² = jJ_g/(ωε₀) ∝ 1/ω`,
hence `W_e ∝ ω⁻²` and `4ωW_e/I² ∝ ω⁻¹` — **the same `1/(ωC)` a physical
capacitance gives.** Gradient-space contamination *is* a spurious electrostatic
response, so it necessarily scales like one, and the sweep separates only the
capacitive family (`ω⁻¹`) from an induction-driven `E = −jωA` (`4ωW_e/I² ∝ ω³`).
It is a cheap sanity check, not the discriminator; do not spend a run on it
expecting an answer.

What does settle it, and why this fixture makes it decisive: the two-torus
fixture has **no conductors** — the tori are tagged *air* subdomains carrying an
impressed `J`, the only metal is the outer PEC wall, and `∇·J = 0` analytically
means no charge. There is therefore **no physical capacitance available to
find**, so the open question is quantitative: does the gradient content of the
load account for all 48.52 Ω? Because the N1curl/CG1 discrete sequence is exact,
that is answerable by two assemblies and no extra solve —
`∫E_h·∇q dV = (j/(ωε₀))∫J·∇q dV` must hold for every `q ∈ CG1 ∩ H¹₀` — with the
energy share following from one cheap scalar Poisson solve. `PORT-1` step 2d in
`PROJECT_PLAN.md` §7 carries the plan; `tests/validation/test_current_divergence.py`
(`POST-3` step 3) is the second, structural route to the same question.

**Update 2026-08-03 — step 2d executed; the answer is "all of it", and the cause
is now measured rather than hypothesised.**
`tests/validation/test_port_gradient_load.py`, 7 passed in 41.5 s at `-n 2` (log
`20260803T183556Z_PORT-1-step2d-gate.log`), same mesh, one curl-curl solve and
one CG1 Poisson solve:

| quantity | value |
|---|---|
| identity `∫E_h·∇q = (j/ωε₀)∫J·∇q`, relative residual | **4.4916e-09** (gated `< 1e-7`) |
| blind control, `j` dropped | 1.4142e+00 = `√2`, as the identity implies |
| `‖P_G J‖²` | `2.534713e-02` (two routes agree to 7.9e-15) |
| `4ωW_e^spur/I² = ‖P_G J‖²/(ωε₀I²)` | **`4.852262e+01 Ω`** |
| measured `4ωW_e/I²` | `4.852271e+01 Ω` |
| **ratio** | **0.999998** |

So the "leading hypothesis" above is confirmed *quantitatively*: the gradient
content of the **discretised** impressed current, amplified by `1/(ωε₀)` on the
subspace where the operator acts as `−k₀²ε_c`, is two-parts-in-a-million the
entire electric-energy excess. Nothing is left for a second mechanism. The
negative diagonal is an artifact of the current representation — not the
reaction integral (step 2b), not physical capacitance (none exists here).

The identity's bound was raised 1e-9 → 1e-7 after a first run measured
4.4916e-09 and failed the plan's house guess; the gate reproduced 4.4916e-09
bit-for-bit. It is a solve-accuracy number, not a physics one — a rank-count,
mesh or solver change that moves it is information, so re-measure rather than
widening again. PROJECT_PLAN §7 `PORT-1` step 2d carries the full reasoning.

**Step 2e executed 2026-08-04 — the fix works, and this entry still stays
open.** `tests/validation/test_port_solenoidal_drive.py` drives the same mesh
with `J′ = J − P_G J` and gets the prediction to three figures
(`20260804T050616Z_PORT-1-step2e-gate.log`, 9 passed in 41.8 s at `-n 2`):

| quantity | projected drive | unprojected (this entry) |
|---|---|---|
| `Im Z₁₁`, both routes | **`+7.437243e+00 Ω`** | `−4.108550e+01 Ω` |
| ratio to Grover's `ωL = 6.818343 Ω` | **1.090770** | −6.03 |
| `4ωW_e/I′²` | `8.761041e-05 Ω` | `4.852271e+01 Ω` |
| `‖P_G J′‖²/‖J′‖²` | `4.5758e-33` | `8.175e-06` |
| complex-power identity residual | `1.6242e-14` | `1.8128e-10` |
| meshed current | `I′ = 0.969001 A` | `I = 0.969009 A` |

`4ωW_m/I′² = 7.4373 Ω` is step 2b's number unchanged, as it must be — the
projection moves `W_e`, not `W_m` — so the fixture's inductance was physical
throughout and the sign is now explained *and* demonstrated, not just
diagnosed.

The entry stays open because **the production driver still builds the
unprojected load**: `TimeHarmonicSolver.solve()` assembles
`−jωμ₀∫J·v̄` from the caller's `current_density` with no projection, so the
diagonal in `test_port_reaction_impedance.py` is still negative and still
ungated. Making the projection the port-excitation default is its own step.

**Consequence, unchanged:** no input impedance and no `S₁₁` off this path means
anything. The off-diagonal is unaffected — `PORT-1` step 2 gates `Im Z₁₂` to
9.35% of `ωM₁₂` and that number does not go through `W_e`.

Remove this entry with the commit that explains the sign.

*(That step is 2f and it landed 2026-08-04 — see the retirement header at the
top of this entry. The off-diagonal did move slightly, as the projection
changes the field and not only its gradient part: `Im Z₁₂` went from
`+1.125614e+00 Ω` (−9.35% of `ωM₁₂`) to `+1.142011e+00 Ω` (**−8.03%**), toward
the closed form, under the unchanged 10% gate.)*

### ✅ RESOLVED 2026-08-03 — "birdcage suite is over the compute budget" was the hang, not meshing cost

The 10:30 review's reinterpretation was right and `GEO-9` step 2b's cost probe
settles the figure. The claim was that
`tests/mesh/test_birdcage_port_tags.py` takes **≈10 minutes** on its own and
that a full `tests/mesh` run exceeded a 700 s bound. **Measured at the fixed
geometry, default parameters, no coarsening:** the file is **8.95 s** of pytest,
harness 10 s (`20260803T200151Z_GEO-9-step2b-probe.log`), and the whole
directory at `-n 2` is **42.15 s**, exit 0
(`20260803T200504Z_GEO-9-step2b-gate.log`). The old number was a poisoned
process burning the harness `timeout` — pytest reported in ≈3 s while the run
exited 124 at 180 s — so `resolution` never had to be coarsened from 0.015 and
the file needs no exclusion from routine runs. It is in CI as of this commit.

The **latent rank-local tag bug** recorded here is also fixed, and it was not
latent: `set(np.unique(cell_tags.values))` is per-rank, and at `-n 2` on the
newly-working mesh rank 0 reported P2/P3 missing while rank 1 reported P1/P4
missing. `GEO-9` step 2b switched the test to
`tests/mesh/helpers.py::global_cell_tag_set()` — the fix that was written and
reverted on request pending exactly this rework. The assertion content (core +
per-port tags all present) is unchanged.

### ✅ RESOLVED 2026-07-30 — truncation-wall modeling floor on the wire/loop fixtures (MAG-13)

The natural BC `n×H = 0` on a truncation wall contradicts Ampère's law for any net
enclosed current (`∮H·dl = I` vs `H_φ(R) = 0` forced at the wall), so it puts a
floor under these fixtures that no refinement removes.

**Straight wire: fixed 2026-07-30** (`MAG-13` steps 1–3). `test_straight_wire.py`
now imposes the analytic `A_z` on the exterior via
`core.solvers.exterior_dirichlet_bc`; measured 35.13% → 22.19% at h=0.004 on the
same mesh, and 22.19% → 12.75% → 9.26% across h = 0.004/0.0025/0.0018, i.e. still
converging at \~O(h^1.2) with no plateau. `J·n ≠ 0` at the end caps remains (the
`MAG-15` multiplier spread measures it) but is not what was dominating.

**Loop: fixed 2026-07-30** (`MAG-13` steps 4–5). `test_circular_loop.py` now
imposes the Jackson 5.37 off-axis `A_φ` on the outer sphere through the same
helper, and `test_convergence.py::test_h_refinement_straight_wire` fits the rate
over three resolutions (**1.10**, bound `[0.7, 1.5]`) instead of two.

One measurement worth keeping: on the *loop*, the analytic wall is ≈20% **worse**
at fixed h than the natural one (16.23% vs 14.98% at h=0.0035; 10.37% vs 8.86% at
h=0.0025, on-axis `B_z` over `|z| ≤ 0.4 R`). Unlike the wire's Ampère-law
contradiction, the loop's natural-BC bias is a PMC image term of order
`(a/R)³ ≈ 3.7%`, which is *smaller* than the O(h) error that degree-1
interpolation of `A_φ` injects through the boundary data. What the Dirichlet wall
buys is the limit: 16.23% → 10.37% → 7.07% at h = 0.0035/0.0025/0.002 converges
monotonically to the analytic field (fitted ≈1.4), where the natural wall
converges to a different field. The loop tolerance is therefore tightened
10% → 8% *at h = 0.002* (411k cells) rather than at the old h = 0.0025 — no
assertion was loosened to accommodate the better boundary condition.

Remaining, not blocking: `J·n ≠ 0` at the wire end caps, and the < 5% wire target
needs h ≈ 0.00125 (≈1.1M cells, > 5 min at `-n 2`) — graded refinement (`MAG-9`),
not more uniform h.

### ✅ RESOLVED 2026-08-01 — `two_torus_domain` was not a conforming mesh (`GEO-8`)

The fixture added two tori and `occ.addBox` over them and never fragmented, so
gmsh meshed a solid box plus two torus islands: total mesh volume exceeded the
analytic box by exactly the two torus volumes (ratio `1.002633`), and driving
torus 1 with a tag-restricted source gave `∫|E|² dV` over tags (1, 2, 3) =
`2.0537e-04, 0, 0` — the field could not leave the driven island, which is why
`PORT-1` measured `Z₁₂ ≡ 0`.

Fixed by `occ.fragment` plus centroid/mass re-derivation of the physical
groups. Both signatures are now gated by
`tests/mesh/test_two_torus_conforming.py`: volume ratio `1.000000000`
(log `20260801T003528Z_GEO-8-after.log`) and air/driven `∫|E|²` ratio
`1.4118`, undriven/driven `5.2088e-08` (log
`20260801T003600Z_GEO-8-field-gate-numbers.log`). The Helmholtz users improved
rather than regressed: centre-field error `1.731% → 0.728%`.

### ✅ RETIRED 2026-08-10 (`EX-17`) — `02_circular_loop.py` never wrote its VTX/`.bp` output

Fixed by the one-file port of the `EX-14` diff: the writers now take the
`A_lag`/`B_lag` Lagrange interpolants the example already builds, each under
its own `try`, and `_check_vtx_roundtrip()` reads `circular_loop_B.bp` back
through ADIOS2 — in-memory and read-back max |B| both **7.756122914931e-05 T**,
relative difference **0.000e+00** vs the 1e-10 tolerance
(`20260810T200154Z_EX-17-gate-mag2.log`, exit 0, 124 s, `-n 2`). The pre-fix
state on record: `⚠ VTX output failed (ADIOS2 may not be available): Only
(discontinuous) Lagrange functions are supported` once per rank, and
`paraview_output/circular_loop_A.bp` opening with zero ADIOS2 variables.

### ✅ RETIRED 2026-08-11 (`POST-4` step 3) — `examples/mri/01` centerline samples were rank-dependent at ≈23%, at and below the gauge floor

**Closed by sampling the source fields.** The example's centerline table now
evaluates `E`/`B` **as solved** through `evaluate_vector_field_parallel`,
instead of the `("Lagrange", 1, (3,))` interpolants it builds for the XDMF
export — the locus the step-1 revision below attributed the spread to. Measured
across `-n 1/2/4` on the post-fix example:

```
                            max rank spread over the -n 1/2/4 pairs
  |E| centerline                     0.000000%   <- every printed digit
  |B| centerline                     0.008613%   <- magnetostatic solve noise
  collapse from the 23.5539% record  2735x       (anchor demanded >= 235x)
```

plus faithfulness — the printed values now equal step 1's measured source
values (|E| to 3.1e-7, |B| to 7.6e-5, the |B| leg's own floor); `|E|` at
z = −0.045 m reads **1.368268e+02** where the interpolant printed
**7.670127e+03**, the 56× artifact — and non-regression: the phantom-region
aggregates reproduce their `EX-16` record to 0.005745% (`-n 2`) and 0.002218%
(`-n 4`), inside the phantom path's own 0.007326% floor.

**Verified at:** `a34c3e6` + the `POST-4` step-3 diff, logs
`20260811T183229Z_POST-4-step3-n1.log`,
`20260811T183211Z_POST-4-step3-n2.log`,
`20260811T183222Z_POST-4-step3-n4.log`, anchor
`20260811T183503Z_POST-4-step3-anchor.log` (PASS, 1 s),
`scripts/probes/post4_step3_spread.py`.

**What is NOT closed by this.** The exported XDMF/VTX fields are still P1
interpolants of Nédélec/DG sources and still carry the vertex-convention
artifact — the fix moved the *printout* off that path, it did not repair the
path. `POST-4` step 4 bounds that residue on the export paths and is open. And
rank-invariance is not physics: `examples/mri/01` is ungated by design
(`WF-1` 🧪), so a faithful printout of an ungated proxy field is still ungated.

*The original entry and its two cause revisions are retained below for the
audit trail.*

#### Original entry (opened `EX-13`, 2026-08-10)

**Test id:** none — `./run_examples.sh -e mri:1 -n 2` vs `-n 4` (the example
is ungated by design, `WF-1` 🧪).
**Symptom:** the five printed centerline `(z, |E|, |B|)` pairs differ between
rank counts far beyond sampling noise. Worst pair at the validated gauge floor
(`gauge_penalty=1.0`), z = +0.0225 m:

```
-n 2:  z=+0.0225 m -> |E|=2.708874e+02, |B|=4.055231e-07
-n 4:  z=+0.0225 m -> |E|=2.592948e+02, |B|=5.304733e-07
```

— **23.5545%** max relative spread across the five pairs (|E| alone: 15.6832%
at z = +0.0450 m). The phantom-region aggregates are stable to ≈0.1%, so this
is specific to the centerline point samples.

**Verified at:** `3c9c0bf`, logs `20260810T050150Z_EX-13-floor-n2.log` /
`…-floor-n4.log` (floor) and `…050120Z_EX-13-subfloor-n2.log` /
`…050133Z_EX-13-subfloor-n4.log` (sub-floor `1e-3`), spread computation
`20260810T050319Z_EX-13-spread.log`.

**Cause:** partly diagnosed, and *not* the gauge. Sub-floor spread is
**23.3010%** — 0.9892× the floor value, i.e. indistinguishable — and the |E|
legs are bit-identical between the two settings because
`TimeHarmonicSolver.solve` ignores `gauge_penalty`
(`src/fem_em_solver/core/time_harmonic.py:351`). The leading suspect is the
unconverged frequency-domain solve: GMRES stops at `ksp_max_it=180` with
`converged=False (reason=-3)` and `residual_norm=1.684628e+00` (identical at
both rank counts), so the returned iterate — not a converged solution — is
partition-dependent. The magnetostatic |B| leg moves < 0.6% with the gauge, so
it cannot account for a 23% spread either. Not confirmed: no run yet with a
converged KSP. *(That suspect was tested and refuted on 2026-08-10 — see the
revision below.)*

**Cause, revised 2026-08-10 (`EX-16`) — the convergence suspect is refuted
and the sampling path is confirmed as the owner.** The demo now solves
direct (`ksp=preonly, pc=lu, converged=True (reason=4)` at both rank
counts) at `gauge_penalty=1.0`, and the spread **does not move**:

```
                       unconverged (EX-13)   converged (EX-16)
max spread, all pairs        23.5545%             23.5539%
  |B| leg (magnetostatic)    23.5545%             23.5539%
  |E| leg (time-harmonic)    15.6832%             13.4499%
```

The anchor's own max is carried by the **magnetostatic |B|** leg, which the
frequency-domain fix cannot touch and which reproduces the unconverged
record to 1.0000×. Converging the KSP bought only 15.68% → 13.45% on the
|E| leg. **Positive control, same two runs, same fields:** the 493-point
phantom-region sampling path agrees across rank counts to **0.007326%**
(|B| mean and all three |E| stats bit-identical) — **3215×** tighter than
the centerline path. Same solve, same field, two samplers: the defect is in
the **centerline point-evaluation path**, not in the solve, the gauge, or
the KSP. Leading (undiagnosed) mechanism: the centerline points sit at
x = y = 0, on mesh edges of the axis, so ownership in
`evaluate_vector_field_parallel` is partition-dependent — the mechanism
`MAG-6` step 4 already characterised for its own centerline metric.

**Verified at (revision):** `34f18de` + the `EX-16` diff, logs
`20260810T170234Z_EX-16-direct-n2.log`, `20260810T170309Z_EX-16-direct-n4.log`,
spread computation `20260810T170457Z_EX-16-spread-v2.log`.

**Not fixed here:** both `EX-13` and `EX-16` were scoped as example/hygiene
chunks with a < 5% rank-stability anchor and an explicit report-and-stop
clause for a converged solve that still spreads. Repairing on-axis point
evaluation is solver-side work on `post.evaluation`, not an example edit.
**Resolved by:** assigned to `POST-4` (2026-08-10, 18:00 review) — step 1
diagnoses the ownership mechanism on the `EX-16` fixture (claim multiplicity,
per-claiming-cell disagreement, `valid_mask`, and an off-axis ε-nudge
discriminator); step 2, conditional on step 1 confirming, replaces
last-writer-wins with a minimum-global-cell-index tie-break in
`evaluate_vector_field_parallel` (`MAG-6` step 4 is the precedent read
first — it measured 0/9 multi-claims on *its* fixture, so the mechanism is
unproven here until step 1 measures it). Entry leaves only with a `POST-4`
step 2 commit whose collapse anchor lands (23.5539% → ≤ 0.1% across
`-n 1/2/4`); stays open on any other outcome.

**Cause, revised 2026-08-11 (`POST-4` step 1) — the ownership tie-break is
refuted; the locus is the interpolation into Lagrange P1, upstream of
`evaluate_vector_field_parallel`.** The probe rebuilt this exact fixture and
solved it at `-n 1/2/4` on a **byte-identical mesh** (9261 cells, coordinate
moments identical to 12 digits at all three rank counts), then instrumented the
point evaluation: it evaluated at *every* colliding cell, not just `links[0]`,
and reduced the claiming `(rank, global cell)` sets across ranks. The
mechanism's necessary condition is simply absent —

```
rows (5 points x 4 fields x 2 point sets x 3 rank counts): 120
MULTI_RANK_CLAIMS   = 0/120
MULTI_CELL_CLAIMS   = 0/120   (links[0] can only bite here)
MASK_INVALID        = 0/120   (silent zero-fill candidate)
CROSS_CELL_DISAGREE = 0/120
```

— every centerline point is claimed by exactly one cell on exactly one rank,
every mask is full, and the ε-nudge to x = y = 1e-6 m does **not** collapse the
spread (97.9755% on axis → 97.9754% nudged, 1.00×) where the chunk's anchor
demanded ≥ 235×. Both candidate mechanisms — partition-dependent cell choice and
silent zero-fill — are therefore dead on this fixture, matching `MAG-6` step 4's
0/9 rather than contradicting it.

**Where the spread actually lives.** The probe sampled four fields at the same
points on the same solves: the Lagrange-P1 interpolants the example prints
(`E_lag`, `B_lag`) and the fields they were interpolated *from* (`E_src`,
`B_src`).

```
                     max rank spread over -n 1/2/4 pairs, on axis
  E_lag (P1 interp)                 97.975464%
  B_lag (P1 interp)                 49.126566%
  E_src (source field)               0.000000%   <- bit-identical
  B_src (source field)               0.008426%
  separation, same points/solves     1.163e+04x
```

The solve is rank-invariant to round-off — `E_src` is bit-identical across all
three rank counts and `B_src`'s 0.008426% sits at the phantom-path control's
0.007326% scale. The 23% enters at `fem.Function.interpolate` into
`("Lagrange", 1, (3,))`: the P1 vertex dof of a field that is not continuous
there is written from whichever adjacent cell writes last locally, which is a
property of the partition. The `-n 2` vs `-n 4` per-point table reproduces
`EX-16`'s record exactly (`B_lag` 23.5539% at z = +0.0225 m), confirming the
probe re-created the fixture rather than a neighbour of it. The `-n 1` leg,
unmeasured before now, is the worst: `E_lag` at z = −0.045 m reads
**7.670127e+03** against 1.5646e+02 / 1.5528e+02 at `-n 2` / `-n 4` and against
the source field's 1.368268e+02 — a 56× interpolation artifact at a single
vertex, present at *every* rank count and merely varying with it.

**Verified at (second revision):** `1c1dc13` + the `POST-4` step-1 diff, logs
`20260811T140414Z_POST-4-step1-n1.log`,
`20260811T140345Z_POST-4-step1-n2.log`,
`20260811T140402Z_POST-4-step1-n4.log`, attribution
`20260811T140549Z_POST-4-step1-attribution.log` (anchor PASS).

**Resolved by (re-pointed):** `POST-4` step 2 as scoped — a min-global-cell
tie-break in `evaluate_vector_field_parallel` — **cannot fix this** and is
skipped per its own conditional clause. The owner is the P1 interpolation of a
non-P1-conforming field (or the example's decision to sample the interpolant
rather than the source field at all); re-scoping that is the next review's.
This entry stays open, with its exit condition unchanged in kind: it leaves with
a commit whose fix collapses the on-axis spread across `-n 1/2/4`, wherever that
fix turns out to live. *(That exit condition was met by `POST-4` step 3 on
2026-08-11 — the second option, the example sampling the source field; see the
retirement block at the top of this entry.)*

### ✅ RESOLVED 2026-08-23 (`OPS-18` step 3a `5df1e39`, landed on `main` by 3b) — the two-torus and straight-wire solves are run-to-run non-deterministic at ≈1e-10 relative, so "bit-identical reproduction" is not an achievable criterion (`OPS-18` step 3a attempt 6, 2026-08-23)

**Found:** 2026-08-23, on `attempt/OPS-18` at `9b3c9e2`, image
`0.11.0.post0` / gmsh 4.15.2 / numpy 2.4.6, `-n 2`, in two back-to-back
runs of the *same* command on an *unchanged* tree
(`20260823T050426Z_OPS-18-step3a-leg1-run1.log`,
`20260823T050903Z_OPS-18-step3a-leg1-run2.log`).

| | |
|---|---|
| **Tests** | `tests/validation/test_port_package_sparameters.py::test_sanity_report_reproduces_the_gated_metrics_on_the_field_route`, `tests/validation/test_port_lumped_two_torus.py::test_step_1_measurements_reproduce` (both already red on 0.11 for the *moved-record* reason); the same effect is visible in real mode in `test_straight_wire.py`. |
| **Symptom** | Repeating one command reproduces every *outcome* and every band, but not every digit: `passivity_max_sigma` 0.8613568946068969 → 0.86135689450373 (1.2e-10 rel), two-torus gap ratio 0.8941410489050936 → 0.8941410492011536 (3.3e-10), `‖S−Sᵀ‖/‖S‖` 3.112128e-05 → 3.112131e-05 (1.0e-06 rel, 3e-11 abs). In one real-mode run `E_Ω` at h = 0.0025 printed 1.0617170184e-01 (ladder) and 1.0617170177e-01 (record test) — 7e-10 apart, same field, same run. |
| **Cause** | Not separated further in-slot; the signature (same mesh, same tree, same image, ≈1e-10) is floating-point summation/factorisation order inside MPI-parallel assembly and the direct solve, i.e. the same class as `MAG-18`'s 7.28e-08 `-n 2` vs `-n 4` floor. It is **not** the version bump: both runs are on one image. |
| **Consequence** | Any acceptance criterion phrased as "reproduces bit-identically" — including `OPS-18` ruling (1) condition (b), 2026-08-22 18:00 review — cannot be met by these fixtures. Records whose written precision is coarser than ≈1e-9 relative are unaffected in practice (0.861356895 and 0.894141 both reproduce as written); the `‖S−Sᵀ‖/‖S‖` record at 7 significant digits does **not**, and would have to be written as 3.11213e-05. |
| **Not a band question** | Every physics band in the same runs is orders of magnitude above the wobble (reciprocity 2.679e-05 vs 1e-3; the symmetry record's own band is 5e-7 absolute vs a 3e-11 move). Nothing here licenses loosening anything. |
| **Resolves with** | A review restating the reproduction criterion at a stated tolerance (proposal in attempts.md 2026-08-23T05:25Z: agreement to ≤ 1e-9 relative across two runs, record written only to digits both runs share). No code fix is implied. |
| **Ruling, 2026-08-23 03:00 review** | Restated as **(b′)**, per-record rather than one relative number — the proposed "≤ 1e-9 relative" would itself reject `‖S−Sᵀ‖/‖S‖` (1.0e-06 relative, cancellation-amplified): *the move across two same-slot runs must be ≤ 1% of the record's own unmoved band, and the value is written only to the digits both runs share, never fewer than the band resolves.* All four records pass (1.2e-4, 3.3e-6, 6e-5 and 7e-6 of their bands); the symmetry record is written as 3.11213e-05. Full text in PROJECT_PLAN §7 `OPS-18`. **This entry closes with the `OPS-18` 3a commit that writes the records** (§9 item 2). |
| **Closed, `OPS-18` step 3b (2026-08-23, 15:00 slot)** | The stated condition is met: the records-writing commit is **`5df1e39`** (3a attempt 8), on `main` with this step's merge. The criterion that replaces "bit-identical" is (b′), quoted verbatim from the ruling above: *"the move across two same-slot runs must be ≤ 1% of the record's own unmoved band, and the value is written only to the digits both runs share, never fewer than the band resolves."* Every record written under it since — the four of attempt 7, the two of attempt 8 (1.2e-4 / 3.3e-6 / 6e-5 / 7e-6, then 6.6e-6 and below the printed digits, all as fractions of their own bands), and step 3b's four `UNIFORM_VOLUMES_RECORD` volumes (move 0.0 across two runs) — satisfies it with no band moved anywhere. **The underlying wobble is not fixed and is not a defect**: it is the ≈1e-10 assembly/factorisation-order floor named in the Cause row, the same class as `MAG-18`'s 7.28e-08 cross-width floor, and it stays the reason no criterion in this project may say "bit-identical" of a *solved* number. |

### ✅ RETIRED 2026-08-25 — Gate (iii) is blind to a broken C4 on the *opposite* class, and the lumped-sheet 4-port sweep loses reciprocity by 223× on an asymmetric layout (`PORT-9` step 3 leg (d1), 2026-08-23)

**Retirement evidence (leg (d1′), 06:00 slot, `20260825T110438Z_PORT-9-step3d1.log`,
`13 passed 106.64 s`; confirmed by `20260825T110643Z_PORT-9-step3d1-consumers.log`,
`24 passed 222.15 s`).** Both findings are disposed on `main`, and the module that
carried them is now on `main` and green rather than parked:

* **Finding 1 (gate blindness) — disposed by tightening, not by widening.** Gate
  (iii)'s 5% became **(iii′) 0.5%** at its single source with all three consumers
  re-run green under it (leg (c) 0.0407%, leg (d0) 0.0040%, leg (d) 0.0553 /
  0.0353 / 0.0214%). Under (iii′) the 22.5° rotation breaks **all three** classes,
  not two: self **6.2219%**, adjacent **7.1142%**, and the opposite class — the
  one this entry named as blind at 1.6476% against 5% — **2.8474%**, an amplifi­
  cation of 133.11× over the symmetric rung. The 08-23 10:30 review's open
  question about whether the opposite class belongs in the geometric control is
  therefore answered affirmatively by measurement.
* **Finding 2 (reciprocity) — disposed by the (d3) power-wave assembly.** On the
  same displaced fixture the route now reads `‖S−Sᵀ‖/‖S‖` = **2.259e-14** against
  the unmoved 1e-3, versus this entry's **5.57e-03**: a **2.466e+11×** separation
  against the (d3) ruling's ≥ 100× bar, and the first test of that fix on a
  fixture that is both 3D and asymmetric. Per (d3c) the reading is an order of
  magnitude only (the confirm run read 6.846e-14).

`PORT-9` ✅ 2026-08-25 at 10 MHz; the branch
`attempt/PORT-9-d1-20260823T124500Z` is deleted. Original text kept verbatim
below.

---

Two findings from one run, both **measured, neither disposed** — the leg's own
negative-result clause (§7 `PORT-9` step 3 leg (d1)) sends both to the review.
The tests are **not on `main`**: they are parked on
`attempt/PORT-9-d1-20260823T124500Z` at `bbe657f`, so nothing here is red in CI.

| | |
|---|---|
| **Tests** | `tests/validation/test_port_birdcage_leg_offset_sweep.py::test_gate_iii_detects_the_broken_c4`, `::test_the_displaced_rung_stays_reciprocal` (parked branch only) |
| **Log** | `docs/testing/logs/20260823T140422Z_PORT-9-step3d1.log` — `2 failed, 7 passed` / 119 s / `-n 2`, complex build, standard tier |
| **Fixture** | `GEO-18`'s gapped, sheeted birdcage, two rungs of the same code path: `leg_azimuth_offsets_rad` all zero (116 416 cells) and leg 1 alone at **+22.5°** (116 944 cells); four driven lumped-sheet solves per rung at `Z_p = z0 = 50 Ω`, 10 MHz, `f = 0.5`, `w = A/h`. |
| **The comparison is controlled** | The zero rung reproduces leg (d)'s recorded 4×4 **entry by entry** to ≤ **2.969e-10** relative against the 1e-9 print-precision band (worst of sixteen), with `‖S−Sᵀ‖/‖S‖` = 2.495292352e-05 and `σ_max` = 0.862659137 identical to nine digits. The offset knob and the frame-aware sheet narrowing this leg added do not move the solve, so every difference below belongs to the displacement. |
| **Finding 1 — symptom** | Displaced, gate (iii)'s three class spreads read **self 5.1819%**, **adjacent 7.1147%**, **opposite 1.6476%** against the unmoved 5% band (symmetric rung: 0.0199 / 0.0180 / 0.0108%; amplifications 260.89× / 395.76× / 152.49×). The adjacent class detects the broken C4 by 1.42×; the **opposite class does not** — it is inside the band it passes on a symmetric layout. The leg's anchor required **both** off-diagonal classes to exceed the band. |
| **Finding 1 — cause** | Geometric and expected in direction, not in size: rotating leg 1 by 22.5° moves the P1–P3 separation 180° → 157.5° while P2–P4 stays 180°, so the opposite class mixes two separations 22.5° apart, where the adjacent class mixes 67.5° / 90° / 112.5°. Leg (d0) measured only 5.9% between 90° and 180°, so a 22.5° perturbation of the *opposite* pair is a second-order effect on an already-flat part of the coupling curve. Not diagnosed further in-slot. |
| **Finding 2 — symptom** | On the displaced rung `‖S−Sᵀ‖/‖S‖` = **5.570640234e-03** against step 2c's unmoved **1e-3** band (`‖Z−Zᵀ‖/‖Z‖` = 7.440778193e-03) — a **223×** rise from the same code path's 2.495292352e-05 on the symmetric rung. `σ_max` = 0.865743230, still passive. Reciprocity is a property of the materials, so this is a systematic of the route or the discretisation, never of the physics. |
| **Finding 2 — cause** | Not diagnosed. The one measured asymmetry that tracks it: the midpoint interior-width filter keeps **26** facets on the rotated port's sheet against **27** on the other three, so P1's `w = A/h` is **7.272128105e-03 m** against 7.413268623e-03 m elsewhere — a 1.9% width difference entering `LumpedSheetPortSpec.sheet_width_m`, hence the V/I estimate, asymmetrically between driven and undriven readings. On the symmetric rung all four sheets are identical and the systematic cancels exactly. This is a hypothesis with a measurement attached, not a diagnosis. |
| **Not a mesh defect** | The negative control of the control is green on **both** rungs: every sheet is a full rectangle of the closed-form `dx·g` = 1.120000000e-04 m², meshed/analytic **1.000000000000** to the 1e-9 band, planar to ≤ 1.7e-17 m in its own port frame, and narrower than the full sheet after filtering. Both rungs' meshes are conforming and every `GEO-18` identity holds. |
| **Consequence** | `PORT-9` stays **🟡**: step 3's gate (iii) is validated as a symmetry gate on the adjacent class only, and the displaced spreads cannot be read as pure geometry while finding 2 stands. §2.2's "no coil has ports" sentence is **unmoved**. Nothing licenses widening (i)–(iii) or the 5% band. |
| **Resolves with** | A review ruling on both: whether gate (iii) is re-specified (a tighter band, an adjacent-class-only statement, or a different invariant), and how the reciprocity systematic is disposed — a per-port equal-facet-count narrowing rule is the obvious first probe, and is code work, not a band question. |
| **Ruling, 2026-08-23 10:30 review** | **The width hypothesis is refuted from the log's own Z** (`…step3d1.log:9327-9330`): a 1.9% readout-width asymmetry on P1 would put a common factor 0.981 on every `Z₁ⱼ/Zⱼ₁` and leave the other pairs at 1; measured `\|Z_ij/Z_ji\|` = 0.99589 / 1.00109 / 1.00625 on row 1 and 1.00523 / **1.01041** / 1.00515 on the pairs that do not involve P1 — the worst pair is **P2–P4, neither port moved**. The asymmetry is global, 0.2–1.6% per pair, the order of the discretisation on a mesh gmsh regenerated whole. Reading: the route's `V` readout is not the impressed source's adjoint, so `Z − Zᵀ` is a local-discretisation residual that cancels only when every port sees the same local mesh — which every fixture this route has been measured on provides (two identical tori, C4 birdcage). **Step 2c's 2.6e-11 is evidence of a symmetric fixture, not of a reciprocal discretisation.** Disposal: `PORT-9` leg (d2), an asymmetric two-torus (`f` = 0.5 / 0.735) with pre-registered predictions (A: O(1e-2); B: ≤ 1e-9), §9 item 2. **Gate (iii) re-specified as (iii′) ≤ 0.5%** (25× above the measured symmetric floor, 3.3× below the weakest displaced class; a tightening; leg (d) stays ✅ under it). The (d1) re-run is serial on (d2). This entry closes with the (d1′) commit. |
| **Finding 2 disposed — leg (d2), 2026-08-23** | **Hypothesis A is refuted and A′ stands: the readout *is* the source's adjoint; the asymmetry is the terminated-`Z` assembly.** `tests/validation/test_port_lumped_sheet_asymmetric.py`, `9 passed` / 198 s and 191 s at `-n 2`, complex, standard (`20260823T183434Z_PORT-9-step3d2.log`, `20260823T183823Z_PORT-9-step3d2-repeat.log`). Two sweeps on **one** 184 919-cell two-torus mesh — control `f` = 0.5/0.5 and asymmetric `f` = 0.5/0.735, `w₂/w₁` = 1.472822047, `Z_p` = 1e6 Ω. Control reproduces step 2c: `‖S−Sᵀ‖/‖S‖` = 2.574356760e-11, **1.078e-15** from the record (band 1e-9). Asymmetric: **8.255602536e-09** — 320.7× the control but **5 orders inside** the unmoved 1e-3, i.e. prediction **B** at the Frobenius grain, so A's O(1e-2) does not happen. Mechanism, asserted at a pre-stated 1e-6: (i) `I₁(drive 2)` = `I₂(drive 1)` to **1.33e-10** — the transadmittance is discretely symmetric, so the current readout **is** the impressed source's adjoint (same facet set `S_i`, same weighting `ĥ_i/(R_i h_i)`, same vector the source is built from, on a complex-symmetric operator); (ii) `Z₁₂/Z₂₁` = `I₁(d1)/I₂(d2)` to **1.33e-10** — `_assemble_impedance_matrix` divides column *j* by the **driven** port's own current, so what it calls `Z` is a *terminated* transimpedance, not the open-circuit matrix reciprocity makes symmetric, and `Z_ij/Z_ji` collapses exactly to the ratio of the two driven-port self-currents (1 for equivalent ports, nothing in particular otherwise). Both runs identical to 8–10 digits. **Read the per-pair number, not the Frobenius one:** here `\|Z₁₂/Z₂₁\|` = **0.997537168** (phase −0.020146017°), a **0.25%** per-pair asymmetry — the same order as (d1)'s 0.2–1.6% table — which the Frobenius ratio hides because at `Z_p` = 1e6 Ω the kΩ diagonal (6.21 − 2.93j / 3.73 − 3.28j) drowns the ≈1.13 Ω mutuals, whereas the birdcage's 50 Ω termination puts `Z₁₁` ≈ 21.7 Ω beside 17 Ω mutuals and lets the same per-pair asymmetry surface as 5.57e-03. **So (d1)'s reciprocity miss is not a discretisation residual and not birdcage-specific: it is the assembly's per-column normalisation, made visible by a matched termination.** Not fixed in-slot per the leg's scope — the fix is an assembly change (an open-circuit `Z`, or `S` from power waves as `_assemble_sparameter_matrix` already does) and it **moves the 2b/2c/(c)/(d0)/(d) records**, which is a review's ruling. |
| **Ruling (2\*), 2026-08-23 18:00 review — fix scoped** | The fix is the **power-wave S assembly** on the gated routes (`S_ij = b_i/a_j`, `a_j` = `V_src/(2√z0)` at a matched drive — symmetric by mechanism identity (i)); the open-circuit-`Z` alternative is rejected on leg (c)'s near-degenerate 1e6 Ω column. The terminated `Z` stays as a documented diagnostic, never reciprocity-gated. Scoped as `PORT-9` legs **(d3)** (two-torus + class re-record under the (1\*) pattern, §9 item 2) and **(d3b)** (birdcage re-record, §9 item 4); **(d1′) serial on (d3b)**. This entry still closes with the (d1′) commit. |

### ✅ RETIRED 2026-08-29 (`GEO-24` step 2b, 06:00 implementer slot) — ~~the **32-port** ring-gap sheet reconstruction is rank-width dependent~~ ~~**`birdcage_port_domain` is built with no ghost layer (`GhostMode.none`), so interior port facets on a partition boundary are unclassifiable — at *every* leg count; it reaches the 4-leg fixture at `-n 12`**~~ — **fixed by the `shared_facet` plumb landed in `470f410` (`GEO-24` step 2a′) and cleared by the two-family re-read; the two-torus `-n 12` drift found inside this entry lives on as `PORT-12` below**

| | |
| --- | --- |
| **Test id** | `tests/mesh/test_birdcage_ring_gaps_scaleup.py::test_thirty_two_ring_ports_at_sixteen_legs` — **not on `main`**; parked on `attempt/GEO-20-step2-20260828T094500Z`. Nothing on `main` is red because of this; the entry exists so the next attempt does not re-derive the measurement. |
| **Verified at** | `e4510a5`, 0.11 image, real build, 2026-08-28. |
| **Symptom** | `birdcage_port_domain(leg_count=16, ring_gap_length=8e-3, emit_port_sheets=True)` — 265 621 cells, 48 ports — is **green at `-n 1`** (`20260828T093352Z_GEO-20-step2-probe1.log`, `1 passed` / 275 s) and **red at `-n 2` on the identical geometry** (`20260828T093839Z_GEO-20-step2-record.log`, Status 1 / 198 s). At `-n 2`, three of the 32 ring sheets do not reconstruct: **P30 and P37 return 0 facets** (meshed/analytic 0.000000000000 of `w²`) and **P45 returns 5 facets** (0.315302109223); the other 29 read **1.000000000000**. The gate that fires is P30's boundary closure, **0.981164653445** against 1e-9. A missing sheet also makes `_sheet_azimuth_deg` return NaN (`inf + -inf` bbox centre), which raises `ValueError: cannot convert float NaN to integer` inside `_azimuth_class` — caught by the module's `_report_safely` guard, so the run fails on the gate rather than hanging. |
| **Not the CAD, the cut, or the volume tagging** | At **both** widths: all 32 port volumes are 1.000000000000 of the analytic wedge, the `GEO-9` partition and air-box closure are 1.000000000000, the ring arcs satisfy Pappus at 1.000000000000, the conductor keeps 0.976465 of its CAD mass, no ring port touches the phantom, and all 32 terminals read 0.974454791–0.974455668 (spread 2.572e-07). The defect is confined to the facet-set reconstruction. |
| **Not the port count alone, and not the ring construction alone** | In the *same* `-n 2` run, the 4-leg ring-gapped fixture (8 ring ports, 110 786 cells) and the 16-leg **leg**-gapped fixture (16 ports, 307 296 cells, C16 sheet spread 1.331e-15) both reproduce their records digit for digit and are green. It takes 32 small port solids for the failure to appear. |
| **Cause** | **Not diagnosed.** Hypothesis with the signature attached: `_interface_facet_tags` matches a facet by the cell tags of its two adjacent **owned** cells, so a port whose `PORT_LOWER+i` and `PORT_UPPER+i` half-regions land on different ranks contributes no facet on either — which is exactly the 0-facet / partial-facet pattern. Discriminator that costs no `src/` change: re-run at `-n 4` / `-n 8` and see whether the broken port *set* moves with the rank count. |
| **Consequence** | `GEO-20` stays **🟡**; step 2 landed **no band and no record**. The 16-leg ring cost rung (110 786 → 265 621 cells, 2.3976×; mesh 23.30 → 72.23 s, 3.1003×) is measured and is safe to cite — it is a cell count and a wall time, not a reconstruction reading. Nothing licenses reading the 32-port fixture's sheets, or building a port model on it, until this is disposed. `GEO-19`'s 16-leg fixture passing at `-n 2` is luck of the partition, not immunity. |
| **Resolves with** | `GEO-20` step 2 attempt 2: confirm or refute the partition reading at `-n 4`/`-n 8` first. If confirmed, the fix is in `_interface_facet_tags` (ghost-layer-aware matching), which is `src/` work touching every module that reconstructs a sheet and could move existing records — a **review's ruling**, not an in-slot fix. |

> **✅ CONFIRMED CLOSED 2026-08-29, 09:00 implementer slot (`GEO-20` step 2
> attempt 2) — the fixture this entry is about now reads green in parallel and
> the module is ON `main`.** The "Test id" row above is superseded:
> `tests/mesh/test_birdcage_ring_gaps_scaleup.py` was restored from `31c08ed`
> and landed with this commit; `attempt/GEO-20-step2-20260828T094500Z` is
> deleted. Re-run on the plumbed tree (`470f410`) at both widths:
> `20260829T140037Z_GEO-20-step2-rerun-n2.log` (Status 0 / 188 s) and
> `…140402Z_…-rerun-n12.log` (Status 0 / 184 s), `1 passed` each. The
> **Symptom** row's three broken sheets are repaired at the identical 265 621
> cells — **P30 and P37 0 → 176 air facets**, **P45 5 facets /
> 0.315302109223 → 180 / 1.000000000000** — and all 32 sheets, all 32
> closures and all 32 `volume/analytic` read **1.000000000000** at `-n 2` and
> `-n 12` alike. The **Cause** row's "not diagnosed" hypothesis
> (`_interface_facet_tags`) is superseded by the `GhostMode.none` diagnosis
> and its one-keyword fix; `_interface_facet_tags` was never touched. The
> **Consequence** row's embargo is lifted: the 32-port fixture's sheets are
> licensed to read.

> **RULED 2026-08-28, 10:30 review — `GEO-20` step 2a queued (§9 item 3):**
> the `-n 4`/`-n 8` discriminator with a per-port print of which rank owns
> each `PORT_LOWER+i`/`PORT_UPPER+i` cell set, no `src/` change; the
> hypothesis predicts a *moving* broken-port set and exact agreement between
> "half-boxes on different ranks" and "sheet broken" across 32 × 3 ports. The
> `_interface_facet_tags` fix is withheld until a review holds that table.
> The attempt branch is kept as the fixture. Entry retires with the fix
> commit, whichever chunk lands it.

> **🔬 DIAGNOSED 2026-08-28, interactive session (at `883c52e`) — the
> hypothesis above is WRONG in its location, and the cause is now measured.
> The defect is not in `_interface_facet_tags`; it is that
> `birdcage_port_domain` never asks for a ghost layer.**
>
> Found by the human operator running
> `examples/meshing/07_birdcage_ring_gap_ports.py` at **`-n 12`**: `port P8
> closure 0.990103697427 on the doubly-gapped mesh`, on the **4-leg,
> 12-port** rung — the fixture `GEO-20` step 1 closed as exact. So the class
> is *not* confined to 16 legs / 32 ports.
>
> **The `-n 4`/`-n 8` discriminator the ruling asked for, answered on the
> small fixture** (probe: build the doubly-gapped rung, print every port's
> closure; logs `20260828T1727*`–`20260828T1730*Z_GEO-20-widthprobe-n{1,2,4,8,12}.log`,
> ≈32 s each, **128 111 cells at every width**):
>
> | ranks | 1 | 2 | 4 | 8 | 12 |
> |---|---|---|---|---|---|
> | red ports | — | — | — | — | **P8** |
>
> Same geometry, same cell count; only the partition moves. The defect is a
> **partition-density** property, not a leg-count or port-count property —
> which is why 16 legs breaks at `-n 2` and 4 legs survives to `-n 8`.
>
> **The lost quantity is one facet.** P8 reads **175** air facets against 176
> on every other ring port; its conductor count (102) and terminal ratio
> (0.974454791) are identical to P7/P11/P12, i.e. the *sheet* reconstructs
> correctly and only the boundary partition loses area.
>
> **Mechanism, measured** (`20260828T173142Z`/`20260828T173221Z_GEO-20-facetprobe-n{1,12}.log`)
> — facets the rank **owns** whose second cell is not present locally, not
> even as a ghost:
>
> ```
> WIDTH  1: owned facets touching P8 = 1417, locally one-sided (counts==1) =  0
> WIDTH 12: owned facets touching P8 = 1417, locally one-sided (counts==1) = 19  (ranks 2 and 5)
> ```
>
> `_interface_facet_tags` selects `interior = counts == 2`, so all 19 drop out
> of every interface group. Most are sheet-plane facets that belong to none of
> the three groups anyway; exactly one faced air, and that is the missing
> 0.99%.
>
> **Cause — diagnosed, one line.** `birdcage_port_domain` calls
> `_model_to_mesh(gmsh.model, comm, rank, gdim=3)` (`io/mesh.py:3356`) with
> **no partitioner**, i.e. gmshio's default `GhostMode.none` — no ghost layer
> at all. `two_torus_domain` was given
> `create_cell_partitioner(GhostMode.shared_facet, 2)` in `PORT-1` step
> 3b-iv for exactly this reason, and that site's comment says so verbatim:
> the port facets are *interior*, classifying one needs the tag of the cell
> on both sides, and "`shared_facet` is what makes that cell present as a
> ghost … Plumbed here only, so no other fixture changes partition." The
> birdcage reconstructs interior port facets the same way and was never given
> the ghost layer that reconstruction requires.
>
> **Confirmed by construction, then reverted** (`20260828T173412Z_GEO-20-ghostfix-probe-n12.log`):
> with the same partitioner plumbed into `birdcage_port_domain`, `-n 12` reads
> **`1.000000000000` on all 12 ports**, P8 back to **176** air facets,
> **128 111** cells unchanged, mesh 26.08 s unchanged. The patch was **not
> landed** — it changes the partition of the fixture `GEO-19`, `GEO-20`,
> `PORT-9` and `PORT-11` all share, which is a re-record question. `GEO-24`
> owns it.
>
> **What this revises in the rows above.** "Resolves with" is superseded: the
> fix is **not** ghost-layer-aware matching inside `_interface_facet_tags`
> (no reconstruction logic changes), and it does **not** touch every module
> that reconstructs a sheet — it is one keyword at one call site. The
> ruling's predicted "moving broken-port set" is confirmed as *width*-moving;
> the per-rank ownership table it asked for is no longer needed to locate the
> cause, though `GEO-24` step 1 still owes the before/after readings.
>
> **A record this puts a caveat on** — ***caveat DROPPED 2026-08-29 by
> `GEO-24` step 2a′/2b; the plumb is landed and the reading is measured true
> at 12 ranks. Kept struck-through for the audit trail:*** ~~`GEO-20` step
> 1's "closure and volume/analytic `1.000000000000` on all 12" is
> **width-conditional** — true at ≤ 8 ranks, false at 12. Every reading taken
> through `_interface_facet_tags` on this fixture inherits that caveat until
> the plumb lands.~~ `GEO-19`'s own remark that passing at `-n 2` is "luck of the
> partition, not immunity" is now measured rather than suspected.

> **📐 `GEO-20` STEP 2a EXECUTED 2026-08-28 (15:00 slot, at `61e97f1`) — the
> ownership table, on the 16-leg / 32-port fixture itself. The set moves with
> the width and agrees with rank ownership port for port, 32 × 2 with no
> exception.** Logs `20260828T200204Z_GEO-20-step2a-n4.log` (Status 1 / 189 s)
> and `20260828T200524Z_GEO-20-step2a-n8.log` (Status 1 / 189 s), standard
> tier, real build, the parked module plus a per-port `allgather`ed count of
> owned `PORT_LOWER+i` / `PORT_UPPER+i` cells. No `src/` change.
>
> | ranks | broken sheets (of 32) | ports whose two half-boxes are not on one rank | sets agree |
> |---|---|---|---|
> | 2 (recorded 04:30) | P30, P37, P45 | *(not instrumented)* | — |
> | 4 | **P25, P29, P37, P41, P45** (5) | P25, P29, P37, P41, P45 | **yes**, ∅ either way |
> | 8 | **P17, P21, P26, P30, P37, P44, P48** (7) | P17, P21, P26, P30, P37, P44, P48 | **yes**, ∅ either way |
>
> The 4-leg / 8-ring-port control in the same runs is **0 broken, 0 straddling
> at both widths** — consistent with the width probe above, where the small
> fixture only breaks at `-n 12`. The set is not nested and not monotone in
> membership (P45 breaks at 2 and 4 but not 8; P30 at 2 and 8 but not 4), which
> is what a partition boundary sweeping through a fixed geometry looks like and
> what no geometry-deterministic defect can produce.
>
> **Failure shape, unchanged from `-n 2`:** a broken port loses its sheet
> **entirely** (0 facets) or keeps a fragment — P29/P45 at `-n 4` **5 facets /
> 0.315302109223** of `w²`, P26 at `-n 8` **6 facets / 0.449137697797** — while
> its *terminal* stays exact (0.974454791–0.974455668, the intact value) and
> its **volume/analytic is 1.000000000000**. Boundary closure drops only on the
> one port per run that also loses an **air** facet (P29 at `-n 4`,
> 0.991120008826, 179 vs 180; P44 at `-n 8`, 0.991064589826, 175 vs 176) —
> the same one-facet loss the 4-leg fixture shows at `-n 12`.
>
> **Negative controls, both widths, digit for digit:** all 40 port
> volume/analytic readings **1.000000000000**, `GEO-9` partition and air-box
> closure **1.000000000000**, ring arcs against Pappus **1.000000000000**,
> conductor 0.976465 / 0.969275 of CAD, kwarg-off at 16 legs **307 296** cells
> (ratio 1.000000) with C16 sheet spread **1.331e-15**, 4-leg ring rung
> **110 786** cells (ratio 1.000000). The volume identity does not route
> through `_interface_facet_tags` and did not move at any width — so the defect
> is confined to facet reconstruction, exactly as the ghost-layer diagnosis
> says.
>
> **Reading.** This is confirmation of the *phenomenology* the 10:30 ruling
> asked for, on the large fixture, and it is consistent with the ghost-layer
> cause rather than with the `_interface_facet_tags` location — ownership
> predicts breakage perfectly because a straddling port is precisely one whose
> interior sheet facets have a neighbour cell that `GhostMode.none` does not
> materialise. Per the ruling, **confirmed is stop**: no `src/` line moved, no
> band, no record. The instrumented module is
> `scripts/probes/geo20_step2a_ownership_scaleup.py` (the parked module verbatim
> plus the ownership print — the attempt branch still holds the uninstrumented
> original; see attempts.md 2026-08-28T20:40Z for why it did not go to the
> branch). The fix and the re-record sweep stay
> `GEO-24`'s, and `GEO-24` step 1 now owes only the before/after readings on
> the modules that already live on `main`.

> **📐 `GEO-24` STEP 1a EXECUTED 2026-08-28 (21:00 slot, at `deef8c5`) — the
> `main`-side "before" table for the seven `tests/mesh/` consumers, at `-n 2`
> and `-n 12`, no `src/` change. Every cell count is identical across widths;
> `-n 2` is green everywhere; two modules are red at `-n 12`, both by facet
> reconstruction.** Fourteen windows, one module per window, `-s`, standard
> tier, real build, `-k 30 400` (`-k 30 570` for `port_scaleup`);
> logs `20260829T0200*`–`20260829T0213*Z_GEO-24-step1a-*-n{2,12}.log`,
> **668 s of compute** in total.
>
> | module | cells (`-n 2` / `-n 12`) | `-n 2` | `-n 12` |
> |---|---|---|---|
> | `test_birdcage_port_sheets` | 116 085 / 116 085 (control 114 655) | ✅ 2 passed, 52 s | ✅ 2 passed, 50 s |
> | `test_birdcage_port_terminals` | 98 666 / 98 666 | ✅ 1 passed, 22 s | ❌ **1 failed**, 22 s |
> | `test_birdcage_ring_gaps` | 128 111 / 128 111 (control 98 666) | ✅ 2 passed, 74 s | ❌ **1 failed, 1 passed**, 74 s |
> | `test_birdcage_leg_gaps` | 114 655 / 114 655 (control 98 666) | ✅ 1 passed, 44 s | ✅ 1 passed, 44 s |
> | `test_birdcage_leg_offset` | 116 085 / 116 475 / 116 085, both widths | ✅ 6 passed, 76 s | ✅ 6 passed, 75 s |
> | `test_birdcage_port_sheet_prerequisite` | 98 666 / 98 666 | ✅ 1 passed, 22 s | ✅ 1 passed, 21 s |
> | `test_birdcage_port_scaleup` | 307 296 / 307 296 (control 116 085) | ✅ 2 passed, 109 s | ✅ 2 passed, 108 s |
>
> **Consumer list re-derived by construction, no difference:**
> `grep -rln birdcage_port_domain tests/ examples/` ∩ the
> `_interface_facet_tags` / `port_sheet` users is exactly these seven under
> `tests/mesh/` — the 18:00 review's list is correct as written.
>
> **The prediction from the diagnosis is confirmed digit for digit.**
> `test_birdcage_ring_gaps` fails at `-n 12` on
> `port P8 closure 0.990103697427` (assert `0.009896302572964588 < 1e-09`) —
> the exact digit the width probe recorded, on the exact port, at an
> unchanged **128 111** cells. Every other reading in that module is
> identical to its `-n 2` value: all 12 ports `volume/analytic
> 1.000000000000`, all 8 ring sheets `1.000000000000`, ring terminals
> 0.974454791 / 0.974454832, leg terminals 0.988615826–0.988615858, Pappus
> 1.000000000000 both gapped and uncut, kwarg-off control **98 666** cells
> at ratio 1.001950.
>
> **New information — the defect reaches a *second* surface, not only port
> sheets.** `test_birdcage_port_terminals` is red at `-n 12` on its
> **positive control**, the phantom↔air interface:
> `phantom surface measures 1.939344e-02 m^2, 0.935322 of the closed-form
> 2.073451e-02 m^2` against the `[0.95, 1.0]` inscribed band — **245 facets
> at `-n 12` against 255 at `-n 2`** (0.979885), i.e. **10** interface facets
> lost, the same one-cell-per-partition-boundary shortfall an order larger.
> All four port boxes in that module stay exact at both widths (air 24 facets
> / 5.200000e-04 m², closure 1.000000000000, conductor 0 facets). So the
> `GhostMode.none` reconstruction gap is **not** specific to port sheets:
> *any* interior material interface on this fixture inherits it. **Step 2a's
> gate must include this reading** — the plumb should return 255 facets /
> 0.979885 at `-n 12`, and if it does not, that is a finding.
> **⚠️ The 255 / 0.979885 in this paragraph is DEFECTIVE (annotated
> 2026-08-29, step 2a′): it is one facet short. `-n 1` on `main`, no plumb,
> reads 256 / 0.984183 (`…050535Z_…-terminals-n1-main.log`) — the serial
> truth, which needs no ghost layer — so step 1a's `-n 2` reading inherited
> the same `GhostMode.none` gap it was meant to be the reference for. The
> post-plumb reading is 256 / 0.984183 at `-n 1`, `-n 2` and `-n 12` alike.
> No test carried the 255 (the band is `[0.95, 1.0]`), so nothing in
> `tests/` moved. Ruled a defect repair, not a re-baseline, by the
> 2026-08-29 03:00 review.**
>
> **Negative control as pre-stated, both widths, digit for digit:** the
> terminal ratios and port-volume identities — neither of which routes through
> a facet reconstruction — are identical at `-n 2` and `-n 12` in every
> module, including the 16-leg scale-up's three azimuth classes
> (0.988615772 / 0.989367514 / 0.989449735, intra-class spreads
> 1.923e-07 / 5.849e-08 / 6.144e-08, inter-class 8.431e-04) and its C16 sheet
> spread (1.331e-15 at `-n 2`, 1.210e-15 at `-n 12` — the only digit that
> moves anywhere in the table, at the 1e-15 floor).
>
> **Cost finding: nothing was unmeasured.** `test_birdcage_port_scaleup` — the
> module the review flagged as most likely to overrun at `-n 12` — completed
> in **108 s**, well inside `-k 30 570`; `GEO-19` step C's exit 124 at 561 s
> was a *bundled* window, not this module's own price. `-n 12` costs the same
> wall clock as `-n 2` throughout (±2 s), since the mesh is built on rank 0
> either way.
>
> **Nothing on `main` is red at `-n 2`,** so CI is unaffected; the two reds
> above are `-n 12`-only and are the measurement `GEO-24` step 1a was
> commissioned to take.
>
> **📐 `GEO-24` STEP 1b EXECUTED 2026-08-29 (22:30 slot, at `d5b4586`) — the
> `main`-side "before" table for the five `tests/validation/` consumers, at
> `-n 2` and `-n 12`, complex build, no `src/` change. Every cell count is
> identical across widths, every module is green at *both* widths, and every
> printed S/Z identity is digit-identical — the validation family shows **no**
> `-n 12` red at all. The one red in the slot is the pre-stated negative
> control, and it is `-n 12`-only.** Thirteen windows, one module per window,
> `-s`, standard tier, `-k 30 480` (`-k 30 400` for the two cheapest,
> `-k 30 300` for the last control); logs
> `20260829T0331*`–`20260829T0342*Z_GEO-24-step1b-*.log`, **660 s of compute**
> in total (env gate + 10 module windows + 2 control windows).
>
> | module | cells (`-n 2` / `-n 12`) | `-n 2` | `-n 12` |
> |---|---|---|---|
> | `test_port_birdcage_lumped_column` | 116 085 / 116 085 | ✅ 2 passed, 33 s | ✅ 2 passed, 31 s |
> | `test_port_birdcage_four_port` | 116 085 / 116 085 | ✅ 5 passed, 51 s | ✅ 5 passed, 40 s |
> | `test_port_birdcage_larmor_probe` | 116 085 / 116 085 | ✅ 3 passed, 40 s | ✅ 3 passed, 32 s |
> | `test_port_birdcage_termination_probe` | 116 085 / 116 085 | ✅ 4 passed, 39 s | ✅ 4 passed, 32 s |
> | `test_port_birdcage_leg_offset_sweep` | 116 085 + 116 475, both widths | ✅ 5 passed, 96 s | ✅ 5 passed, 77 s |
>
> **Consumer list re-derived by construction, no difference:**
> `grep -rln birdcage_port_domain tests/ examples/` ∩ the
> `_interface_facet_tags` / `port_sheet` users is exactly these five under
> `tests/validation/` — the 18:00 review's list is correct as written. Note
> that only `_lumped_column` calls `birdcage_port_domain` directly; the other
> four reach the fixture through its `_build` helper and appear in the
> intersection through their own `_interface_facet_tags` / sheet imports.
> `test_port_birdcage_larmor_gate.py` and `_larmor_gate_128.py` are **not** in
> the intersection (they reconstruct nothing themselves) and were not read.
>
> **Digits, both widths, identical unless stated.** All four
> `Z_{11,21,31,41}` reproduce their `PORT-9` records at rel. deviation
> 1.07e-10 – 2.57e-10 in every module that gates them; 4×4 gates in
> `_four_port`: `||S−S^T||/||S||` **8.141422487e-15** at `-n 2` and
> **1.116856988e-13** at `-n 12` (band 1e-3, both PASS; `||Z−Z^T||/||Z||`
> 8.814400605e-05 vs …604e-05 reported), `sigma_max(S)` **0.999992805** and
> max column power sum **0.793823974** at both, C4 class spreads
> **0.0553 / 0.0353 / 0.0214 %** at both (band 0.5%), pooled off-diagonal
> 9.2115% and separation 166.6766×. `_lumped_column`: all four sheets 26
> facets / 5.835298880e-05 m² / `w = A/h` 7.294123600e-03 m, out-of-plane
> 0.000e+00 m, at both widths. `_larmor_probe`: `Z_11` +2.215494591e+01
> +7.460189773e+00j (10 MHz) and +2.647082952e+01+4.646185233e+01j (Larmor
> rung), identical at both. `_termination_probe`: margin **2256.9707×**,
> spread **0.0040%**, open control 1.5951× / 0.0407%, `I_1` reproducing at
> 5.9e-12 / 1.1e-11. `_leg_offset_sweep`: displaced rung 116 475 cells (ratio
> 1.003360) with class spreads **6.2219 / 7.1142 / 2.8474 %** vs the zero
> rung's 0.0553 / 0.0353 / 0.0214 %, `sigma_max` 0.999992337 vs 0.999992805,
> `||S−S^T||/||S||` zero 1.044255156e-14 → 6.958642293e-14 and displaced
> 2.009039801e-14 → 4.532499019e-13 across widths (band 1e-3).
>
> **⚠️ The pre-stated negative control did not hold, and it is the step's new
> information.** `tests/validation/test_port_lumped_two_torus.py` — the
> fixture that *already* has `create_cell_partitioner(GhostMode.shared_facet,
> 2)` (`PORT-1` step 3b-iv) — is **green at `-n 2`** and **red at `-n 12`**:
>
> ```
> gap ratio: 0.894274 against step 1's record 0.894141 — moved by 1.33e-04,
>   above 1e-04; step 2's reads changed step 1's solve
> ```
>
> (`20260829T034112Z_…-twotorus-control-n12.log`, `1 failed, 4 passed`,
> Status 1, 84 s; `20260829T034253Z_…-twotorus-control-n2.log`, `5 passed`,
> Status 0, 84 s, gap ratio **0.894141** = the record exactly.) The mesh does
> **not** move — **184 176 cells at both widths** — and the four other tests in
> the module pass at `-n 12`, so this is not the `GhostMode` reconstruction
> defect: it is a *solved* quantity (the gap-route `V` line integral,
> `Im Z12 = 1.110469250` at `-n 12` vs `1.110303775` at `-n 2`, 1.5e-4
> relative) drifting with partition count on an already-plumbed fixture, at
> the same order as the `OPS-18` step-3 re-record band it is gated against
> (1e-4). **Consequence for step 2b:** a `-n 12` *solved* digit that moves at
> 1e-4 after the birdcage plumb is not by itself evidence the plumb failed —
> the gate must distinguish facet-reconstruction readings (which must return
> exactly 1.000000000000) from solve-derived digits at 1e-4. Whether the
> two-torus 1e-4 band should be width-qualified is a **review's call**, not
> this chunk's; nothing was loosened here.
>
> **Nothing on `main` moved at `-n 2` in this family** — every `-n 2` reading
> matches its `PORT-9` / `PORT-11` record, so CI is unaffected and no record
> is owed a re-write from step 1b. Step 2b (the plumb re-read of this family)
> is unblocked.
>
> **📐 `GEO-24` STEP 2a EXECUTED 2026-08-29 (00:00 slot, at `169c28c`) — the
> plumb repairs both `-n 12` reds, every cell count and every other `-n 2`
> digit holds, and it moves *one* `-n 2` reading, which is §9 item 4's
> pre-stated stop. The plumb is REVERTED on `main` and parked on
> `attempt/GEO-24-step2a-20260829T052300Z` (`e1dede8`); nothing under `src/`
> landed.** Sixteen windows, one module per window, `-s`, standard tier, real
> build, `-k 30 400` (`-k 30 570` for `port_scaleup`); logs
> `20260829T0501*`–`20260829T0519*Z_GEO-24-step2a-*.log`, **≈ 870 s** of
> compute (14 table windows + 2 serial diagnostics + 1 control window).
>
> The patch is one keyword at `io/mesh.py:3356` —
> `partitioner=create_cell_partitioner(GhostMode.shared_facet, 2)`, the
> `two_torus_domain` site's kwarg and comment, nothing else in `src/`.
>
> | module | cells (`-n 2` / `-n 12`) | `-n 2` before → after | `-n 12` before → after |
> |---|---|---|---|
> | `test_birdcage_port_sheets` | 116 085 / 116 085 (control 114 655) | ✅ → ✅ 2 passed, 66 s | ✅ → ✅ 2 passed, 62 s |
> | `test_birdcage_port_terminals` | 98 666 / 98 666 | ✅ → ✅ 1 passed, 26 s **(digit moved, below)** | ❌ → ✅ 1 passed, 26 s |
> | `test_birdcage_ring_gaps` | 128 111 / 128 111 (control 98 666) | ✅ → ✅ 2 passed, 85 s | ❌ → ✅ 2 passed, 85 s |
> | `test_birdcage_leg_gaps` | 114 655 / 114 655 (control 98 666) | ✅ → ✅ 1 passed, 51 s | ✅ → ✅ 1 passed, 51 s |
> | `test_birdcage_leg_offset` | 116 085 / 116 475 / 116 085, both widths | ✅ → ✅ 6 passed, 84 s | ✅ → ✅ 6 passed, 83 s |
> | `test_birdcage_port_sheet_prerequisite` | 98 666 / 98 666 | ✅ → ✅ 1 passed, 23 s | ✅ → ✅ 1 passed, 23 s |
> | `test_birdcage_port_scaleup` | 307 296 / 307 296 (control 116 085) | ✅ → ✅ 2 passed, 118 s | ✅ → ✅ 2 passed, 117 s |
>
> **Gate clause 1 — every cell count identical: ✅.** All seven modules read
> the same counts as step 1a at both widths, and every kwarg-off control
> reproduces (98 666 at ratio 1.001950, 114 655, 116 085, the scale-up's
> `cells 116085 vs 116085 (delta 0, relative 0.000e+00)`). The plumb changes
> partitioning, not meshing, exactly as predicted.
>
> **Gate clause 3 — both previously-red `-n 12` readings now green: ✅.**
> `test_birdcage_ring_gaps` port P8 returns to **176 air facets, closure
> 1.000000000000** from 175 / **0.990103697427**, at an unchanged 128 111
> cells, with all 12 ports' `volume/analytic` and all 8 ring sheets still
> 1.000000000000. `test_birdcage_port_terminals`' phantom↔air control returns
> to **256 facets** from 245, and all four port boxes stay exact (air 24
> facets / 5.200000e-04 m², closure 1.000000000000, conductor 0).
>
> **Gate clause 2 — every `-n 2` digit identical: ❌ in one cell, and that is
> the stop.** `test_birdcage_port_terminals`' positive control reads
> **256 facets / 2.040655e-02 m² / 0.984183** of the closed-form
> 2.073451e-02 m² after the plumb, against step 1a's recorded **255 facets /
> 0.979885** at `-n 2`. The test passes either way (the band is [0.95, 1.0]),
> but the digit moved, and §9 item 4's negative-result clause is explicit:
> a moving `-n 2` digit stops the chunk for a review. Every *other* `-n 2`
> digit in all seven modules is identical to step 1a's — the C4 sheet spread
> 6.050e-16, leg terminals 0.988615825–0.988615858, ring terminals
> 0.974454791 / 0.974454832, Pappus 1.000000000000, and the scale-up's three
> azimuth classes 0.989367514 / 0.989449735 / 0.988615772 with intra spreads
> 5.849e-08 / 6.144e-08 / 1.923e-07, inter-class 8.431e-04, C16 sheet spread
> 1.331e-15 at `-n 2` and 1.210e-15 at `-n 12` (the 1e-15 floor, as before).
>
> **What the moved digit actually is — measured, not argued.** Two extra
> serial windows settle it. On the **plumbed** tree at `-n 1`:
> **256 facets / 0.984183**. On **`main`** at `-n 1`, plumb reverted:
> **256 facets / 0.984183** — identical (`…050500Z_…-terminals-n1-plumbed.log`
> and `…050535Z_…-terminals-n1-main.log`, Status 0, 26 s / 25 s, 98 666 cells
> both). A single rank has no partition boundary, so it needs no ghost layer
> and reads the true interface either way. **So 256 is the truth, and the
> recorded `-n 2` value of 255 / 0.979885 was itself one facet short from the
> same `GhostMode.none` gap** — the record was *defective at every parallel
> width*, not partition-dependent physics. After the plumb the reading is
> 256 / 0.984183 at `-n 1`, `-n 2` and `-n 12` alike, i.e. width-independent
> and equal to the serial truth. This is the fix working on a surface where
> step 1a had mistaken a short reading for the baseline.
>
> **Negative controls, pre-stated, with the plumb applied: green and
> unmoved.** `tests/mesh/test_two_torus_port_sheet.py` and
> `tests/mesh/test_cylindrical_domain.py` — untouched fixtures — are
> `4 passed` / Status 0 / 33 s with the `GEO-16` control at its usual
> **79 070** cells and tags `[1, 2, 3, 101, 102] / [1, 201, 202]`
> (`…051958Z_…-controls-n2.log`). No other fixture's partitioner was touched.
>
> **Cost: nothing unmeasured.** `-n 12` again costs the same wall clock as
> `-n 2` throughout (±2 s), the mesh being built on rank 0 either way;
> `port_scaleup` took 117 s at `-n 12` inside `-k 30 570`.
>
> **What a review owes this chunk.** The measurement says the plumb is safe
> and correct on the mesh family and that step 1a's 255 is a defective
> record to be re-written to 256 / 0.984183 with the `-n 1` provenance —
> but re-writing records is **step 3's**, and item 4's scope is explicit that
> no record moves in step 2a. So: rule on the 255 → 256 re-record, then the
> parked branch lands as-is. Step 2b (the validation family) is unaffected by
> this stop — it reads a family that showed no reconstruction red at all —
> but it should not run before the plumb's disposition is ruled, since it
> measures the same patch.
>
> **RULED 2026-08-29, 03:00 review — defect repair, not re-baseline; land
> the parked commit.** The `-n 1` control on `main` (no plumb,
> `…050535Z_…-terminals-n1-main.log`, 256 facets / 0.984183, 98 666 cells)
> is the serial truth: one rank has no partition boundary and reads the
> interface whole. Step 1a's `-n 2` **255 / 0.979885** was therefore one
> facet short from the same `GhostMode.none` gap — a *defective record*,
> not partition-dependent physics — and no test carries it
> (`test_birdcage_port_terminals` gates `[0.95, 1.0]`). Disposition:
> §9 item 1 (`GEO-24` step 2a′) cherry-picks `e1dede8`, re-reads
> `port_terminals` at `-n 1/2/12` and `ring_gaps` at `-n 2/12` on the
> landed tree, and annotates the 255 in this table as defective; §9 item 2
> (step 2b, the validation family) retires this entry on green; the
> two-torus `-n 12` finding above is **split out** as its own entry
> (`PORT-12`, below) so this entry can retire without losing it; `GEO-20`
> step 2 is re-run as §9 item 4 after the plumb lands. No band, tolerance
> or record in `tests/` moves under this ruling.
>
> **Step 2a′ ✅ 2026-08-29, 04:30 slot — the plumb is LANDED on `main` and
> the defective digit is disposed of.** `e1dede8` cherry-picked onto
> `31a4e0b` (`470f410`; `git diff HEAD~1 -- src/` is the one
> `io/mesh.py:3356` kwarg + comment hunk, working tree clean after). Re-read
> on the landed tree in six windows / **246 s** of container time:
> `test_birdcage_port_terminals` at `-n 1` / `-n 2` / `-n 12`
> (`…093031Z`, `…093103Z`, `…093130Z`, Status 0, 23 / 21 / 22 s) and
> `test_birdcage_ring_gaps` at `-n 2` / `-n 12` (`…093201Z`, `…093326Z`,
> Status 0, 75 s each), all prefixed `GEO-24-step2aP-`.
>
> **Every anchor met, digit for digit.** Phantom↔air positive control
> **256 facets / 2.040655e-02 m² / 0.984183 at all three widths** on
> **98 666** cells — the serial truth now reproduced in parallel, which is
> the repair. `ring_gaps` port P8 at `-n 12`: **176 air facets / closure
> 1.000000000000** (was 175 / 0.990103697427) on **128 111** cells, with all
> 12 ports' `volume/analytic` 1.000000000000 and all 8 ring sheets
> `meshed/analytic` 1.000000000000; the ring-gapped rung is **110 786**
> cells and Pappus reads 1.000000000000 gapped and uncut. Every other
> printed digit is identical to the step-2a table: leg terminals
> 0.988615826 / 0.988615832 / 0.988615854 / 0.988615858, ring terminals
> 0.974454791 / 0.974454832, the four leg port boxes air 24 facets /
> 5.200000e-04 m² / closure 1.000000000000 / conductor 0, and the kwarg-off
> control 98 666 cells / 0.966977. `-n 2` and `-n 12` agree everywhere.
>
> **Negative control, pre-stated: green and unmoved.**
> `test_two_torus_port_sheet` + `test_cylindrical_domain` `4 passed`,
> Status 0, 30 s, `GEO-16` control at **79 070** cells with tags
> `[1, 2, 3, 101, 102] / [1, 201, 202]` (`…093452Z_…-controls-n2.log`). No
> cell count moved anywhere, so the landed commit is the one that was
> measured. The `attempt/GEO-24-step2a-20260829T052300Z` branch is deleted.
>
> **This entry stays open**: step 2b (the validation family, §9 item 2)
> retires it. Nothing loosened, no band moved, no record in `tests/`
> touched.
>
> **✅ Step 2b 2026-08-29, 06:00 slot — the validation family re-reads clean
> on the plumbed tree at both widths, and THIS ENTRY IS RETIRED.** Eleven
> windows / **485 s** of container time at `470f410`, complex build,
> `FEM_EM_REQUIRE_COMPLEX=1`, `-s`, `-k 30 480`, one module per width per
> window, no `src/` change in the slot; logs
> `20260829T1100*`–`20260829T1108*Z_GEO-24-step2b-*.log`, all **Status 0**.
> Environment gate `11 passed` / 21 s first.
>
> | module | cells (`-n 2` / `-n 12`) | `-n 2` | `-n 12` |
> |---|---|---|---|
> | `test_port_birdcage_lumped_column` | 116 085 / 116 085 | ✅ 2 passed, 33 s | ✅ 2 passed, 30 s |
> | `test_port_birdcage_four_port` | 116 085 / 116 085 | ✅ 5 passed, 49 s | ✅ 5 passed, 39 s |
> | `test_port_birdcage_larmor_probe` | 116 085 / 116 085 | ✅ 3 passed, 38 s | ✅ 3 passed, 33 s |
> | `test_port_birdcage_termination_probe` | 116 085 / 116 085 | ✅ 4 passed, 38 s | ✅ 4 passed, 33 s |
> | `test_port_birdcage_leg_offset_sweep` | 116 085 + 116 475, both widths | ✅ 5 passed, 97 s | ✅ 5 passed, 74 s |
>
> **The gate was pre-stated in two classes, and both are met** (the 03:00
> ruling's lesson, and step 1b's).
>
> *(i) Reconstruction readings — required identical to the digit, and they
> are.* All four `_lumped_column` sheets **26 facets / 5.835298880e-05 m² /
> `w = A/h` 7.294123600e-03 m / out-of-plane 0.000e+00 m** at `-n 2` and
> `-n 12` alike (same in `_four_port`), full-sheet bbox 1.400000000e-02 m,
> filtered bbox 9.167340025e-03 m; every cell count identical across widths
> and equal to step 1b's records at ratio **1.000000** (the displaced rung
> 116 475, ratio 1.003360 against the shared record, exactly as before).
>
> *(ii) Solved digits — required inside each module's own in-file band, and
> every one passes at both widths.* `_larmor_probe`: `Z_{11,21,31,41}`
> reproduce their `PORT-9` records at rel. deviation **1.071e-10 –
> 2.568e-10** (`-n 2`) and **1.071e-10 – 2.566e-10** (`-n 12`), phantom
> `cells/delta` 5.9213 (10 MHz control 12.0002). `_four_port`:
> `sigma_max(S)` **0.999992805** and max column power sum **0.793823974** at
> both widths, C4 class spreads **0.0553 / 0.0353 / 0.0214 %** (band 0.5%) at
> both, pooled off-diagonal 9.2115% and separation 166.6766×,
> `||S−S^T||/||S||` 1.044255156e-14 (`-n 2`) / 1.897457072e-14 (`-n 12`),
> band 1e-3. `_termination_probe`: margin **2256.9707×** and spread
> **0.0040%** at both, open control 1.5951× / 0.0407%. `_leg_offset_sweep`:
> displaced rung breaking (iii′) at **6.2219 / 7.1142 / 2.8474 %** against
> the zero rung's 0.0553 / 0.0353 / 0.0214 %, amplifications 112.58× /
> 201.52× / 133.11×, at both widths. The only `-n 12` movement anywhere is
> in the last displayed digits of `Z_11` in `_lumped_column`
> (+9.201557829e+02−4.718342449e+03j → +9.201557791e+02−4.718342444e+03j,
> **4.1e-9** relative) and in the two Frobenius asymmetry residuals at the
> 1e-14 floor — reported, orders inside every band, and far below the 1e-4
> the `PORT-12` precedent set as the threshold worth reporting at all.
>
> **Negative control** (pre-stated: this family has no kwarg-off control, so
> the control is step 1b's own table): every `-n 2` digit reproduces step
> 1b's `main`-side reading, so the plumb changed nothing this family could
> see — which is the expected result, step 1b having shown the family took no
> damage from `GhostMode.none` in the first place. The two-torus control was
> **deliberately not re-run**: its `-n 12` red is `PORT-12`'s and is not
> moved by this patch.
>
> **Disposition.** `GEO-24` **✅** — the ghost layer is plumbed and both
> consumer families are re-read at `-n 2` and `-n 12`. `GEO-20` step 1's
> "closure and volume/analytic `1.000000000000` on all 12" loses its
> **width-conditional caveat** (above): it is now measured true at 12 ranks
> on the landed tree (`GEO-24` step 2a′, port P8 176 air facets / closure
> 1.000000000000). `GEO-20` step 2 remains 🟡 as a **re-run** of the parked
> 32-port module, now unblocked. Step 1a's `-n 2` **255 / 0.979885** stands
> annotated *defective* (one facet short; `-n 1` truth 256 / 0.984183).
> Nothing loosened, no band moved, no record in `tests/` touched.

### ✅ RETIRED 2026-08-30 (`PORT-12` step 2, 13:30 implementer slot) — ~~the two-torus gap-route reproduction record drifts with rank width on a fixture that already has the `shared_facet` ghost layer~~ — **the record is now stated as the `-n 2` statement it always was, and the drift is *bounded* on every wider width by a pre-registered `PARALLEL_DRIFT_ENVELOPE = 3e-4` rather than left red: `-n 8` (the worst width) reads +2.06e-04 and passes, `-n 2` reads the record exactly, and the lumped route's width-flatness is asserted at 1e-8 as the module's negative control** (found 2026-08-29 by `GEO-24` step 1b)

| | |
| --- | --- |
| **Test id** | `tests/validation/test_port_lumped_two_torus.py` — the `STEP1_GAP_RATIO_RECORD` reproduction assertion (`REPRODUCTION_BAND` 1e-4). Green at `-n 2`, which is what CI runs; **red at every parallel width above 2** (`-n 4`, `-n 8`, `-n 12` measured). |
| **Verified at** | `d5b4586`, 0.11 image, complex build, 2026-08-29 (`20260829T034253Z_GEO-24-step1b-twotorus-control-n2.log`, `5 passed` / 84 s; `20260829T034112Z_GEO-24-step1b-twotorus-control-n12.log`, `1 failed, 4 passed` / Status 1 / 84 s). Widths 4 and 8 added by `PORT-12` step 1 at `c4630ed` (`20260829T170059Z_PORT-12-step1-twotorus-n4.log`, `1 failed, 4 passed` / Status 1 / 87 s; `20260829T170240Z_PORT-12-step1-twotorus-n8.log`, `1 failed, 4 passed` / Status 1 / 81 s; env gate `20260829T170032Z_PORT-12-step1-env.log`, `11 passed` / 21 s). |
| **Symptom** | `gap ratio: 0.894274 against step 1's record 0.894141 — moved by 1.33e-04, above 1e-04` at `-n 12`; `-n 2` reads **0.894141**, the record exactly. `Im Z12` 1.110303775 (`-n 2`) → 1.110469250 (`-n 12`), 1.5e-4 relative. The mesh does **not** move — **184 176** cells at both widths — and the other four tests in the module pass at `-n 12`. |
| **Four-width table** (`PORT-12` step 1, 2026-08-29) | All four widths on **184 176** cells. Gap route / lumped route / cross-route and `Im Z12`:<br>`-n 2` — gap **0.894141** (= record), `Im Z12(gap)` 1.110303775, `V_gap` +1.365256733e-02+1.079044036e+00j; lumped **0.828893**, `Im Z12(lumped)` 1.029281338; cross-route 7.743060e-02.<br>`-n 4` — gap **0.894274** (+1.33e-04), `Im Z12(gap)` 1.110469342, `V_gap` +1.368962224e-02+1.079204774e+00j; lumped **0.828893**, `Im Z12(lumped)` 1.029281337; cross-route 7.754834e-02.<br>`-n 8` — gap **0.894347** (+2.06e-04), `Im Z12(gap)` 1.110559796, `V_gap` +1.370291038e-02+1.079292623e+00j; lumped **0.828893**, `Im Z12(lumped)` 1.029281336; cross-route 7.761484e-02.<br>`-n 12` — gap **0.894274** (+1.33e-04), `Im Z12(gap)` 1.110469250, `V_gap` +1.373904726e-02+1.079204448e+00j; lumped **0.828893**, `Im Z12(lumped)` 1.029281338; cross-route 7.753298e-02. |
| **Negative control (held)** | Every reconstruction digit is **identical at all four widths**: 184 176 cells, sheet 212 **1583 owned facets**, meshed/CAD area **1.000000000000**, `w` 1.040000000e-02 m, `h` 1.395505060e-02 m, `w/h` 0.745249896 squares, out-of-plane spread **0.0e+00 m**, meshed/analytic gap volume **1.000000000000**. This is **not** the `GEO-24` class of defect on a plumbed fixture. |
| **Not the `GEO-24` defect** | `two_torus_domain` has carried `create_cell_partitioner(GhostMode.shared_facet, 2)` since `PORT-1` step 3b-iv; the quantity that moves is a *solved* line integral (`V = −∫E·dl` on the gap route), not a facet reconstruction. It was run as `GEO-24` step 1b's pre-stated negative control ("an already-plumbed fixture is width-independent") and the control **failed**. |
| **Cause** | **Classified, not yet root-caused. The drift is an evaluation-path effect confined to the gap route, and it is *non-monotone* in width.** The two pre-stated candidates separate cleanly: (a) a *solve-side* drift would move all three routes together — it does not. The **lumped route reads the same solved field through the sheet's own law and is flat to 2e-09** across all four widths (`Im Z12(lumped)` 1.029281338 / …337 / …336 / …338; `I_sheet` −4.122422e−08−1.000166e−06j at every width), and the step-2 *surface* read of the same field, `mean E.yhat over the sheet`, is **bit-identical to every printed digit at all four widths** (shadow −2.958541e+00−7.177866e+01j, fringe +8.607682e-03−1.009219e-02j, ratio 0.000185). So the solved field itself is width-independent to ≈1e-9, five orders below the gap route's 1.3e-04–2.1e-04 motion. (b) The **gap route alone moves**, and not monotonically: +1.33e-04 at `-n 4`, +2.06e-04 at `-n 8`, back to +1.33e-04 at `-n 12` — `-n 8`, not `-n 12`, is the worst width, so this is not a "more partitions ⇒ more drift" law. The cross-route figure tracks the gap route exactly (it is derived from it), and the step-2 path/projection residual is likewise non-monotone: 0.0689 / 0.0632 / 0.0662 / 0.0836 pp at 2 / 4 / 8 / 12. One sub-shape worth the root-cause hunt: `Re V_gap` **is** monotone in width (1.365256733e-02 → 1.368962224e-02 → 1.370291038e-02 → 1.373904726e-02, 6.5e-03 relative from `-n 2` to `-n 12`) while `Im V_gap` is not — consistent with the `V = −∫E·dl` path picking up partition-dependent contributions where it crosses a partition boundary, rather than with a converged field being integrated correctly. |
| **Consequence** | The 1e-4 band is a `-n 2` statement until this is classified; every `PORT-1`/`OPS-18` two-torus record is quoted at `-n 2`. `GEO-24` step 2b (the birdcage validation family after the plumb) must not read a `-n 12` solved digit moving at ≤ 1e-4 as evidence the plumb failed. Nothing loosened, no record re-written. **Step 1 adds:** because the drift is on the evaluation path and not the solve, any *other* `V = −∫E·dl` gap-route reading in the package is suspect at parallel width in the same way — but no lumped-sheet port reading is, the lumped route being flat to 2e-09 here. |
| **Resolves with** | ~~`PORT-12` step 1~~ **✅ 2026-08-29** (the four-width table and the classification above). **Step 2 is the 2026-08-30 weekly review's call** with this table in hand: the drift is *not* monotone and *not* shared by all three routes, so the "solver-side fix" option in the original framing is off — the choices are (i) width-qualify `REPRODUCTION_BAND` as a `-n 2` statement, (ii) widen it to a pre-registered parallel band ≥ 2.1e-04 with the non-monotone table as the justification, or (iii) commission a root-cause step on the gap-route line integral's partition crossing. |
| **RULED 2026-08-30, weekly planning review — option (i) with a bounded envelope; option (iii) declined with an epitaph** | `REPRODUCTION_BAND` = 1e-4 stays what it is and is **stated as a `-n 2` record** (every `PORT-1`/`OPS-18` two-torus digit is quoted at that width; the constant's comment says so). At `comm.size > 2` the same assertion runs against a **separate pre-registered `PARALLEL_DRIFT_ENVELOPE = 3e-4`** whose provenance is the four-width table above (max observed +2.06e-04 at `-n 8`, 1.46× headroom), so the drift is *bounded on every width CI might run* rather than hidden by a skip; and the lumped route's width-flatness — `Im Z12(lumped)` within **1e-8** relative across widths, measured 2e-9 — is asserted as the module's new negative control, because that flatness is the reason the production port model is unaffected. **Not (ii)** as framed: widening the *record* band would let the `-n 2` record itself drift. **Not (iii):** the gap route is the two-torus `V = −∫E·dl` estimator, which no birdcage or Larmor quantity reads — `PORT-9`/`PORT-11` are lumped-sheet throughout and flat here to 2e-9 — so a root-cause hunt on a partition-crossing line integral is off the mission's shortest path; it re-opens the day a production quantity is read through a gap-route integral at parallel width. This is `PORT-12` **step 2**, scoped in §7 for the daily review to queue (smoke, real+complex, `-n 2` and one `-n > 2` window). This entry **retires with step 2's landing**. |
| **RETIRED 2026-08-30 by `PORT-12` step 2 (13:30 implementer slot)** | The ruling above is implemented in `tests/validation/test_port_lumped_two_torus.py`, tests only, no `src/` change. Three windows / **310 s**, complex, `FEM_EM_REQUIRE_COMPLEX=1`, `-s`, `-k 30 300`, `tests/environment` first in every window. **Negative control first, on unpatched `main`:** `-n 8` footered **Status 1**, `1 failed, 15 passed` / 105 s (`20260830T183101Z_PORT-12-step2-control-n8-main.log`), `gap ratio: 0.894347 … moved by 2.06e-04, above 1e-04` — the envelope is visibly load-bearing. **Patched:** `-n 8` **17 passed** / Status 0 / 103 s (`20260830T183340Z_…-patched-n8.log`) printing `gap ratio at -n 8: 0.894347 … drift +2.06e-04 against the 3e-04 envelope`, lumped ratio drift **−2.99e-07**, cross-route **+1.84e-04**, `Im Z12(lumped) 1.029281338` at relative **2.344e-10**; `-n 2` **17 passed** / Status 0 / 102 s (`20260830T183533Z_…-patched-n2.log`) with the gap ratio **0.894141** — the record exactly, inside the unmoved 1e-4 — and `Im Z12(lumped)` at **4.649e-10**. **The 1e-8 assert was probed load-bearing** as pre-stated: pointed at the gap route's `Im Z12` 1.110303775 it fails at relative **7.297e-02** (`20260830T183730Z_…-probe-n2.log`, Status 1, one-line edit reverted, not committed). `STEP1_GAP_RATIO_RECORD` and `REPRODUCTION_BAND` untouched; no record re-written, no band widened. |

### ✅ RETIRED 2026-08-30 (`WF-6` step 1d, 12:00 implementer slot) — ~~the first `|B₁⁺|` map is C4-covariant to only ≈9%, against a 5% pre-registered discretisation band~~ — **the ≈9% was the DG0 cell-scatter floor, and gate (ii) now reads the L²-projected CG1 field: 2.1870 / 2.1146 / 1.8911% at +90° / −90° / 180° against the unmoved 5% band, `main` green** (found 2026-08-29 by `WF-6` step 1)

| | |
| --- | --- |
| **Test id** | `tests/validation/test_birdcage_b1_plus_map.py::test_b1_plus_map_is_c4_covariant_under_the_drive_rotation` — gate (ii). Red at `-n 2`, which is what CI runs. The module's other two tests, including **gate (i)**, the three-way power accounting, pass. |
| **Verified at** | `bea89f3` + this commit, 0.11 image, complex build, `-n 2`, 2026-08-29 (`20260829T183450Z_WF-6-step1.log`, `1 failed, 13 passed` / Status 1 / 89 s, with `tests/environment`; re-read with the ungated diagnostics at `20260829T183728Z_WF-6-step1-diagnostic.log`, `1 failed, 2 passed` / Status 1 / 87 s). |
| **Symptom** | `|B1+| from the P2 drive at the 90deg-rotated point disagrees with the P1 drive by 8.6516% in relative l2 over 51 phantom centroids, outside the pre-registered 5.0% discretisation band`. The sample set is the 51 tag-3 cell centroids with `r ≤ 0.02 m`, `|z| ≤ 0.02 m` on the 116 085-cell `GEO-19` step-B fixture; `|B₁⁺|` there reads mean 2.077398e-08 T, max 2.834980e-08 T, min 1.457925e-08 T at `V_src = 1 V`. |
| **The band is not passing on noise, and the estimator resolves azimuth** | The **180° negative control holds with room**: the P3 drive against the same 90°-rotated points reads **27.3161%**, 3.2× the failing reading and 5.5× the band. So the comparison is measuring the drive's azimuth, not returning scatter for everything. |
| **The miss is systematic, not a few outlier cells** | Printed and never gated: the pointwise `|B₁⁺|` deviation over the 51 points is median **6.7395%**, p90 **15.0357%**, max **17.5662%** — a broad distribution, so no handful of cells can be blamed. The **second instance of the same identity** — the P4 drive at −90°, which the module solves for exactly this purpose — reads **9.5808%**, alike to P2's 8.6516%. Both 90° instances agreeing rules out anything peculiar to P2 and points at the shared mechanism: `B` is DG0 on a gmsh mesh that is not itself C4-symmetric, so a sample point and its rotated image sit in *different* cells, and the DG0 cell-scatter of a curl is the floor of this comparison. |
| **Cause** | **Not root-caused; two candidates left open for a review.** (a) The pre-registered 5% band simply underestimated DG0 cell-to-cell scatter for a curl at this mesh resolution — the reading would then be an honest measurement of the estimator's floor, and the fix is a better estimator (a CG1 projection of `B`, or cell-volume-weighted comparison, or sampling on a rotation-invariant point set rather than on centroids) rather than a looser band. (b) A real C4 asymmetry in the solved field beyond what the terminal quantities show — but the same fixture's `Z` classes spread ≤ 0.5% (`PORT-9` gate (iii′)), which bounds the *terminal* asymmetry an order of magnitude below this. Nothing was widened, nothing was refitted, and both gates stand as pre-registered. |
| **Not a power-accounting problem** | Gate (i) closes at **9.795751e-03** of the supplied power at the P1 drive and **9.796209e-03** at P2, inside its 1e-2 band, with the shares 0.0008% phantom / 6.5374% conductor / 92.4822% sheets; its negative control (drop the conductor term) misses by 7.517001e-02, 7.7× the band. So the field is energetically accounted for at both drives that gate (ii) compares. |
| **Consequence** | `WF-6` stays **🧪** — step 1's `post/` helpers (`magnetic_flux_density_from_e`, `b1_plus`) are landed and exercised, gate (i) is closed, and **no B₁⁺ homogeneity, CV or absolute-accuracy claim exists**. `main` carries this one red deliberately, per the step's own pre-registered negative-result clause ("record the mismatch, keep the asserts, never widen either band in-slot"). |
| **Resolves with** | `WF-6` steps **1b** (CG1-projected `B` on the same 51 points, the 180° identity read for the first time) and **1c** (a 96-point rotation-invariant ring sample at DG0) — both scoped by the 2026-08-29 18:00 review and queued §9 items 1 and 3, independent of each other, measurement only. The 5% band must not be moved without their tables in hand; the review then either re-registers gate (ii) on the better estimator (candidate (a)) or commissions a field-side hunt (candidate (b)). Noted at scoping: 51 centroids in a 0.02 × 0.04 m cylinder means ≈ 1 cm phantom cells — the DG0 curl scatter candidate is the fixture's resolution, not a mystery. |
| **Step 1b executed 2026-08-29 19:30 — the estimator leg, verdict (a)** | `20260830T003238Z_WF-6-step1b.log`, `1 failed, 15 passed` / Status 1 / 98 s (the one failure is *this* entry's gate (ii), unchanged and deliberately red; both new tests pass). The DG0 `B_phasor` was L²-projected onto `("Lagrange", 1, (3,))` through a Hermitian mass-matrix `LinearProblem` (CG/Jacobi, `ksp_rtol` 1e-12 — never `interpolate`, which is ill-defined at vertices for a DG0 field), and `\|B_x + jB_y\|/2` formed from the projected vector at the *same* 51 points. Anchors all reproduced: DG0 P2-at-+90° **8.6516%**, gate (i)'s P1 residual **9.795751e-03**, both at rtol 1e-4; `valid` all-true, 51 of 51, on every rotated image. **The three-angle × two-estimator table:** <br>`P2 @ +90°` — DG0 **8.6516%** (med 6.7395, p90 15.0357) │ CG1 **2.1870%** (med 1.5240, p90 3.3040) <br>`P4 @ −90°` — DG0 **9.5808%** (med 5.1948, p90 13.4830) │ CG1 **2.1146%** (med 1.5757, p90 3.4706) <br>`P3 @ 180°` — DG0 **8.5970%** (med 5.1290, p90 13.6265) │ CG1 **1.8911%** (med 1.3170, p90 2.8471) <br>`P3 @ +90°` (mis-rotated control) — DG0 **27.3161%** │ CG1 **23.2642%**, both outside the 5% band, asserted. `\|B₁⁺\|` mean over the set 2.077398e-08 T (DG0) vs 2.069556e-08 T (CG1) — the projection moves the magnitude by 0.38%, not the map. |
| **What the table decides** | The pre-registered verdict is **(a) — the estimator floor**, unambiguously: CG1 is inside 5% at *all three* covariance angles, a factor 4–5 below DG0, while the mis-rotated control survives the projection at 23%, so the smoothing has not smoothed the map away. Read the 180° column, which step 1 never had: DG0 gives **8.5970%** there against 8.6516% at +90° — the *same* miss at both angles, which is the signature of a scatter floor and the opposite of what a C2-preserving, C4-breaking field asymmetry (candidate (b)) would produce. Candidate (b) is not supported by any reading on this fixture. **No band was moved:** re-registering gate (ii) on the CG1 estimator, with this table as the new band's provenance, is a **review's call** — the slot recorded and stopped, per the step's scope. `WF-6` stays 🧪 and the gate stays red until then. |
| **Step 1c executed 2026-08-29 22:30 — the sample-set leg, the set is not the mechanism** | `20260830T033147Z_WF-6-step1c.log`, `1 failed, 18 passed` / Status 1 / **97 s** (the one failure is again *this* entry's gate (ii), untouched). Estimator held at DG0, sample set replaced by one closed under the C4 rotation: rings at `r ∈ {0.005, 0.010, 0.015, 0.020}` m × `z ∈ {−0.015, 0, +0.015}` m × 8 azimuths in 45° steps, azimuth start jittered 3.7° off the coordinate planes — 96 points, every ±90° and 180° image a member of the set. Anchors: `valid` **96 of 96** on all four drives and every rotated image; centroid-set DG0 P2-at-+90° reproduced **8.6516%** and gate (i)'s P1 residual **9.795751e-03**, both rtol 1e-4; the mis-rotated control P3-at-+90° **25.8213%**, asserted outside the band. **Ring-set table (centroid-set figure, delta):** <br>`P2 @ +90°` **9.9271%** (med 6.9433, p90 16.2927) — centroid 8.6516%, **+1.28 pp** <br>`P4 @ −90°` **9.9519%** (med 7.3968, p90 16.4548) — centroid 9.5808%, **+0.37 pp** <br>`P3 @ 180°` **8.4706%** (med 5.7448, p90 13.5804) — centroid 8.5970%, **−0.13 pp**. `\|B₁⁺\|` over the ring set, P1 driven: mean 2.023327e-08 T, max 3.263326e-08, min 1.419703e-08. |
| **What the ring set decides** | The pre-registered verdict is **"sample set is not the mechanism"** — all three angles land within ±2 pp of the centroid set (max \|Δ\| 1.28 pp), so the centroid set's lack of closure under the rotation was not manufacturing the miss and the ≈9% floor is the **DG0 scatter itself**. This is the one thing step 1b could not distinguish, and it corroborates 1b's verdict (a) from the opposite direction: 1b changed the estimator and the miss fell 4–5×; 1c changed the sample set and the miss did not move. The 180° column agrees with +90° on the ring set too (8.47 vs 9.93%), so candidate (b) remains unsupported. **Per-ring structure, for the review to read against the coil geometry:** no monotone radial trend — 6.33…11.65% at `r = 0.010`, 4.61…12.63% at `r = 0.020`, the single lowest ring being the outermost top one (`r = 0.020, z = +0.015`, 4.61 / 6.21 / 3.96%) and the highest an inner one (`r = 0.005, z = −0.015`, 11.25 / 12.27%); ring-to-ring spread of the same order as the overall figure is what a per-cell scatter looks like. **No band was moved**; gate (ii) stays red, `WF-6` stays 🧪, and re-registering the gate on the CG1 estimator remains the review's call. |
| **RULED 2026-08-30, weekly planning review — gate (ii) is re-registered on the CG1-projected estimator; the DG0 5% assertion is replaced, not loosened** | The two legs read as a pair: change the estimator and the miss falls 4–5× (1b); change the sample set and it does not move (1c). That is the pre-registered candidate (a), and the 180° column (DG0 8.60% ≈ +90° 8.65%) rules out candidate (b) on this fixture. **Ruling:** (1) the production `|B₁⁺|` map estimator for gates and examples is `b1_plus` of the **L²-projected CG1 `B`** (mass-matrix `LinearProblem`, never `interpolate`), with the DG0 field kept as the raw curl; (2) gate (ii) becomes the CG1 covariance identity at **all three angles** (+90°, −90°, 180°) on the 51 centroids, band **5%** unchanged in value but now with a *measured* provenance — CG1 reads 2.19 / 2.11 / 1.89%, p90 ≤ 3.47%, so the band carries 2.3× headroom and is a discretisation floor someone has measured; (3) the mis-rotated 180°-vs-90° control stays asserted **> 5%** under CG1 (23.26%); (4) the DG0 readings are **printed and recorded**, not gated — the old DG0 assert is removed with its record (8.6516 / 9.5808 / 8.5970%) cited in-comment, because gating a curl at DG0 on a non-symmetric mesh gates the mesh, not the map (the `GEO-19` step-C precedent). No CV, homogeneity or absolute-accuracy claim follows from (ii) closing — it is a symmetry identity. This is **`WF-6` step 1d** (§7), one smoke/standard slot on the existing fixture; steps 2–3 stay serial on it. This entry **retires with step 1d's landing**; `main`'s deliberate red retires with it. |
| **RETIRED 2026-08-30 — step 1d executed, the ruling is in the code** | `20260830T170242Z_WF-6-step1d.log`, **19 passed / Status 0 / 97 s** at `-n 2` with `tests/environment`, complex build. The projector moved into the package as `fem_em_solver.post.project_to_cg1` (mass-matrix `LinearProblem`, CG/Jacobi, `ksp_rtol` 1e-12) and `test_b1_plus_map_is_c4_covariant_under_the_drive_rotation` is now the CG1 covariance identity at all three angles: **2.1870% / 2.1146% / 1.8911%** at +90° / −90° / 180°, each also asserted against step 1b's record at rtol 1e-3, all against the **unchanged** `C4_COVARIANCE_BAND = 5e-2`. Negative controls both hold: the mis-rotated P3-at-+90° reads **23.2642%** under CG1, outside the band; the ungated DG0 column reads **8.6516 / 9.5808 / 8.5970%**, unmoved (asserted at rtol 1e-4 — a moved DG0 column would mean the projection changed the field, not the estimator), `valid` 51 of 51 on every image. The `src/` change owed two example re-runs, both green through the §9 runner substitution (the docker-socket denial recurred): `ports:4` 76 s Status 0 (`20260830T170431Z`), `ports:5` 127 s Status 0 (`20260830T170559Z`); doc-reference census `exit=1` with `dead=53 stale=2`, none of them a `ports_04_*` or `ports_05_*` artifact — the 53 are `EX-36`'s, unchanged. `WF-6` is **🟡 with step 1 ✅** — a symmetry identity only, and **no B₁⁺ homogeneity, CV or absolute-accuracy claim follows**. |

*(The `ANS-1` / `ANS-3` `ModuleNotFoundError` entry — the 2026-08-28 artifact
rename that also rewrote the two `__import__` strings — was **removed
2026-08-31 by `EX-37`**, which restored both strings and re-ran both cases
green. The negative control that finally observed the symptom in a harness
log, and the two green re-runs, are recorded in the `EX-37` §7 entry.)*

### ✅ RETIRED 2026-09-07 (`TH-15` step 2d, 12:00 implementer slot) — ~~The gap-displacement current and the conduction current are different quantities on an **undriven** port, and no `h` reconciles them — `TH-15` step 2c's anchor (i) misses by ≈100% and its negative control fails (2026-09-07, 04:30 implementer slot)~~ — **retired as ruled: the two deleted anchors compared two different quantities to a third that neither is, so nothing was loosened and nothing was refined. The route landed on `main` at step 2d under the 10:30 review's ruling (1), gated on the open-circuit condition instead: `\|I_disp,undriven\|/I_drive` = **1.541112e-06** against the pre-registered `OPEN_CIRCUIT_BAND` = 1e-4 on both drives, with the conduction volume average asserted ≥ 100× above it and measured **782.9×** / **1096.8×** (`20260907T170302Z_TH-15.log:1050–1062`). (ii) 2.243038e-02 / 2.236565e-02, (iii) 5.2613e-04, (iv) 0.909618 reproduce the branch digits exactly (`:1071, :1076, :1078–1079`); 14 passed in 134.38 s at `-n 4`, elapsed 136 s (`:1150, :1349`). The deleted readings 9.998674e-01 / 9.993525e-01 and the 0.99× separation are kept in the module docstring with their log lines.**

| | |
| --- | --- |
| **Chunk** | `TH-15` step 2c (§7 `TH-15`, §9 item 1 of the 2026-09-07 03:00 review), the `current_route="gap_displacement"` readout on `GapVoltagePortSpec`. Code is parked, not landed — `attempt/TH-15-step2c-20260907T094500Z` (`f942dc2`). Nothing on `main` is red because of this. |
| **Verified at** | `20260907T093440Z_TH-15.log` (`-n 4`, complex, `tests/environment` first, **2 failed / 13 passed in 129.69 s**, elapsed 131 s, `:1140`) and `20260907T093741Z_TH-15.log` (`-s`, **2 failed / 2 passed in 102.27 s**, elapsed 103 s, `:1069`), on the solid two-torus, 184 176 cells, one 2×2 sweep at 10 MHz. Both logs are on the parked branch. Preceded by a 4 s collect-only smoke, `20260907T093428Z_TH-15.log`. |
| **The four readings the §7 item's stop text asks for** | (i) **undriven** port continuity `\|I_disp/I_cond − 1\|` = **9.998674e-01** (drive P1, port P2, `…093741Z:1000`) and **9.993525e-01** (drive P2, port P1, `:1005`) against the pre-registered **0.10** — *miss*. (ii) **driven** port continuity `\|(I_drive + I_disp)/I_cond − 1\|` = **2.243038e-02** / **2.236565e-02** (`:1013`, `:1018`) — *pass*. (iii) reciprocity `‖S − Sᵀ‖/‖S‖` = **5.2613e-04** inside the imported `S_SYMMETRY_BAND` 1e-3 (`:1030`) — *pass*. (iv) the corrected mutual `Im Z₁₂/ωM₁₂` = **0.909618** (raw 0.865226) inside the imported `MUTUAL_TOLERANCE` (`:1029`) — *pass*, beside the conduction route's 0.939822. Gap material read off the problem: `σ_gap = 0.000000e+00 S/m`, `ε_gap/ε₀ = 1.000000` (`:982–986`). |
| **The asserted negative control also fails** | The item asserted (rule (e) label: *asserted*, backed by `20260907T020614Z_TH-15.log:1052`) that the gap route beats the Ampère-loop route's undriven record `0.032 + 0.190j`, predicted ≈50× separation. Measured on the identical port of the identical solve: gap route **9.998674e-01** vs loop route **9.860947e-01**, separation **0.99×** (`…093741Z:1022–1025`). The new route is *not* better conditioned on the undriven port's current than the one it was scoped to replace. |
| **Cause — a definition mismatch, not resolution** | The undriven gap current **1.5407e-06 A** is *correct* as what it is: the gap capacitor's own current `jω(ε₀A_gap/g)V₂` with V₂ = 1.067 V (`:983`), 1e-6 class. The conduction route's `I_cond = σ/L ∫_conductor E·φ̂ dV` is a **volume average of the current density around a ring that is open at the gap**. On an open-circuited port the wire current vanishes at the gap faces and peaks opposite them, so the volume average and the gap-face current are simply not the same quantity and continuity between them does not hold at any `h`. On the *driven* port the impressed 1 A dominates both terms, which is exactly why (ii) reads 2.2% and (i) reads 100%. The 03:00 review's ruling (1) predicted the gap field is "the port's own by two decades" — that part is true; what does not follow is that the conduction route measures the same current. |
| **Caveat on the passing anchors** | (iii) and (iv) pass on off-diagonal currents ≈1e-3 of the conduction route's. That is a scale-invariance of the two-port `Z → S` reduction (both routes' `Z₁₂ = V₁/I₂` consume only the *driven* current), not independent evidence for the route. Do not cite (iii)/(iv) as validating the undriven readout. |
| **Nothing was loosened** | The 0.10 bound, `S_SYMMETRY_BAND`, `MUTUAL_TOLERANCE` and the loop-route record are all as pre-registered; the failing asserts were left failing on the branch and the negative-result exit was taken. Rule (c)'s `test_port_package_sparameters.py` re-run was **not** executed — nothing landed on `main` to regress. |
| **Resolves with** | A **review ruling**, not an implementer's call (§9 items 1 and 6 are 🚫 pending it). Two candidates: (a) assemble the full 2×2 from **driven**-port gap currents alone — one solve per port, each read at its own driven gap, which is precisely what anchors (ii)/(iv) support and which needs no conductor cells (so it survives on the PEC hole); (b) the `n × H` facet form over the conductor surface tag, the last candidate the step-2b ruling named. `TH-15` stays 🟡; both `attempt/TH-15-step2b-…` and `attempt/TH-15-step2c-…` are kept until the ruling. |
| **Ruling, 2026-09-07 10:30 review — (a), and it is the package's existing definition** | `_assemble_impedance_matrix` (`src/fem_em_solver/ports/sparameters.py:234–253`) builds `Z[i, k] = V_i / I_k` from the **driven** port's current in solve `k` and the path voltages of every port in that solve; the undriven current has never entered `Z` on any route, so anchor (i) gated a quantity no consumer reads. The premise it rested on — that the gap-face current and the conduction volume average are the same current on an open ring — is refuted by the readings above and by the definitions: the volume average is the ring's induced current closing through its own stray capacitance, the gap face carries the gap capacitor's current, and an open-circuit `Z` needs the *terminal* current to vanish — which is what the gap route reads (1.5407e-06 A against a 1 A drive, `…093741Z:983, 1000, 1005`). **Disposition:** `TH-15` step 2d (§9 item 1) lands `f942dc2`'s code by path checkout, **deletes** the undriven-continuity test and the loop-record control (each compared two different quantities to a third that neither is; their readings stay in the module docstring with these log lines), and adds the open-circuit anchor `\|I_disp,undriven\| / I_drive ≤ 1e-4` with the conduction volume average asserted ≥ 100× above it as the separation control (783× / 1097× measured). (ii)–(iv) are untouched. **This entry is retired by step 2d's commit, not by the ruling.** The `n × H` facet form is not scoped. |

---

## Batch moved 2026-09-26

### ℹ️ RESOLVED-IN-CODE 2026-09-22 (`OPS-62`) — the `ANS-4` step 3e XL window exited 1 on a spurious guard; its four rungs' data are intact

| | |
| --- | --- |
| **Symptom** | `20260922T070007Z_ANS-4-step3e.log` has Status 1 and one failure, `test_every_finer_rung_actually_refines`: `resolution 0.015 gave 116085 cells, not more than 0.015's 116085` (`:8802–8805, 8812`). |
| **Cause** | The guard asserted cells, not unknowns, on a ladder whose first two rungs are the **same mesh** at degree 1 then degree 2 by construction (`FEM_EM_ANS4_STEP2_RUNGSPEC="0.015:1 0.015:2 0.0075:2 0.005:2"`); the discretisation refined ≈139 302 → ≈742 944 unknowns (`:2452, 4542`). Not a physics or mesh defect. |
| **Read the log as** | A **usable measurement with one spurious red.** All four rungs built and solved; the record-rung check, the imported `PORT-11` reciprocity / passivity / C4 gates on every rung, and the full readout all passed (`:8814–8817` and following). |
| **Fixed by** | `OPS-62` (2026-09-22): the guard now counts unknowns via `_estimated_dofs`, and a control asserts that exact ladder shape with that window's own cell counts passes. Whether step 3e is re-run is a review's call, not this entry's. |
| **Correction + ruling (2026-09-23 03:00 review, `log-pathologist` reading)** | The window had **two** reds, not one — `2 failed, 15 passed` on all 16 ranks (`:8969–8970`): the guard above and `test_pinned_scikit_rf_version_is_exact` (`ModuleNotFoundError`, `:8892–8895`), the deliberate `OPS-60` red (entry retired 2026-09-23 when `OPS-60` closed ✅, `20260923T093017Z_OPS-60.log:76–77`), picked up because entry 7's command runs `tests/environment` first. Neither is physics. The four rungs' true sizes are 139 140 / 745 938 / 1 819 094 / 3 828 494 unknowns (ZMUMPS `N`, `:2114`, `:4204`, `:6281`, `:8450`). **Step 3e is not re-run**: its readout is complete and usable. `xl-ledger.md` row filled; `xl-pending.md` entry 7 marked RUN. The next `ANS-4` window on the same module (3f, 09-23) ran 18 / 18 green under the fixed guard (`20260923T070007Z_ANS-4-step3f.log:4923–4924`, `:5046`). |

---

## Batch moved 2026-09-30

### ✅ RETIRED 2026-09-30 (§9 item 119, 19:30 CDT slot; opened 2026-09-30 by item 118) — `import_step`'s `outer_boundary_tag` created no facet group on the F-small STEP

> **Retired by the exterior-skin predicate** (`20261001T003201Z_GEO-34-wall-W1.log`,
> 5 passed, 89 s, `-n 2` real). `import_step` now takes the wall as the
> combined boundary of all volumes (no tolerance) and raises if it is empty;
> the preflight's dolfinx-side fallback is removed. **Cause, measured:** every
> wall's OCC bbox is 2.000e-7 thick on its normal axis against the old
> predicate's tol 2.400e-10 (`:3815–3821`), so the thickness term alone
> rejected all six; the pinned `4857675d` predicate accepts 0 of the 6 new
> walls (`:3822`). (A1) owned tag-1 facets 3310 = owned exterior facets 3310,
> 0 not on exactly one cell (`:3813`); (A2) ∫ds(1) = 0.3072 m², rel
> −3.3e-16 (`:3814`); CAD wall area rel −2.2e-16 (`:3812`). Congruent
> config: 3310 = 3310, 0 off, and C_file = C_solve now includes the import's
> own wall (`20261001T003338Z_GEO-34-wall-W2.log:5718–5722`, 4 passed,
> 105 s); (n1) unchanged (`20261001T003527Z_GEO-34-wall-W3.log`, 1 passed,
> 109 s); 1b′ re-gate Z/S digits identical to `20260926T093715Z`
> (`20261001T003722Z_GEO-34-wall-W4.log`, 19 passed, 217 s complex). Cells
> 116 662 / 116 685 unchanged. Original entry follows.
>
> **Original heading:** 🟡 OPEN 2026-09-30 (`GEO-34` step 2, 15:00 CDT slot, §9 item 118) — `import_step`'s `outer_boundary_tag` creates no facet group on the F-small STEP: the import returns `facet_tags = None`

| | |
| --- | --- |
| **Tests** | `tests/mesh/test_step_preflight.py::test_census`, `::test_no_facet_grid` in their first form (both errored in fixture `solve_census`); no committed test fails now. |
| **Symptom** | `20260930T201313Z_GEO-34-step2-W1.log` (commit `9fc225d0` + the item's new files): `AttributeError: 'NoneType' object has no attribute 'indices'` at `_global_facet_count(msh, facet_tags, 1, comm)`; the preflight printed `facets per tag {1: 0, 211: 66, …}`. The green re-run prints `import's own wall-tag facets 0; facet_tags is None: True` (`20260930T201605Z_GEO-34-step2-W1.log:5717`). Congruent and 1a configs both hit it (W2 `20260930T201805Z_GEO-34-step2-W2.log:3762`, `:5588`). |
| **Verified at** | `fcffa701`. |
| **Cause** | Not diagnosed. Leading hypothesis: `io/step_import.py:196–211` matches wall surfaces with `tol = 1e-9 · extent` against `gmsh.model.getBoundingBox`, and OCC bounding boxes are padded. The model bbox read back gives a wall surface of 3.072005440e-01 m² against the meshed exterior's 0.3072 m², rel −1.77e-6 (`20260930T201605Z_GEO-34-step2-W1.log:5730`), which is far above that tol, so `walls` stays empty and no group is added. Steps 1a/1b never read the wall tag. |
| **Workaround** | `io/step_preflight.build_preflight_mesh` tags the wall dolfinx-side as the mesh's exterior facets (`exterior_facet_indices`) only when the import produced none, and prints that it did. The test counts C_solve's wall independently (owned facets with one attached cell): 3310 = 3310. |
| **Resolves with** | A chunk that fixes the wall predicate in `import_step` (a tolerance tied to the OCC bbox padding, or a comparison against the imported air box's faces), with 1a re-run green and a wall-count assert. The preflight fallback then goes dead and can be removed. |

### ✅ RETIRED 2026-09-29 (`MAT-1` §9 item 83, 19:30 CDT slot; opened 2026-09-29 by item 81) — `tests/materials/test_phantom_material_model.py::test_gelled_saline_presets_match_lossy_sphere_closed_form` **premise was red intermittently**

> **Retired by the re-anchored premise** (`20260930T003229Z_MAT-1-reanchor-W1.log`,
> 17 passed, 123 s, `-n 2` complex). The solve-output 1e-12 gate is withdrawn
> (not widened); the assembly identity it stood for reads exactly:
> **(P1)** σ / εᵣ / μᵣ max \|Δ\| over owned dofs **0.0**, 0 differing, owned
> 74 020 = global cells (`:485–487`); (P1′) 19 496 owned σ == 0.57 dofs on
> both paths (`:488`); **(P2)** same CSR pattern, ‖A_φ − A_map‖_∞ = **0.0**,
> ‖A_map‖_∞ = 4.894e4 (`:489–490`); **(c1)** 1-ULP σ reads max \|Δσ\| =
> 1.1102230246251565e-16 = 2⁻⁵³ on 19 496 dofs, Δεᵣ = Δμᵣ = 0 (`:491`);
> **(c2)** σ·(1 + 1e-9) reads ‖A″ − A_map‖_∞ = 1.178e-9 (2.41e-14
> relative) (`:493`). The printed solve-level record this run: 3.392e-11 /
> 1.281e-11 (SAR / ⟨E_z⟩, `:494`), inside the measured floor. The floor
> itself stays documented under the `OPS-43` observation below. Original
> entry text follows.
>
> **Original heading:** 🔴 OPEN 2026-09-29 (`MAT-1`, 18:30 implementer slot, §9 item 81) — premise red intermittently: `phantom_material` vs `material_map` on the same mesh reads **2.544e-11 / 1.200e-11** (mean SAR / ⟨E_z⟩) against the item's 1e-12 bound; the solve's own run-to-run floor is above that bound

> **Literal symptom** (`20260929T184225Z_MAT-1-W1-final.log:485,492`, Status 1,
> 78 s, `-n 2` complex, 0.11 image, tree at `a097bed2` + the uncommitted test):
> `AssertionError: phantom_material and material_map assemble different fields:
> SAR rel diff 2.544e-11, <E_z> rel diff 1.200e-11`. The identical code in the
> window two runs earlier read **9.047e-13 / 5.962e-13** and passed
> (`20260929T183740Z_MAT-1-W1.log:485`, 17 passed, 118 s).
>
> **Cause — the bound, not the path (measured, not fully diagnosed).** A
> scratch probe solving the `material_map` case twice and the phantom case
> once on the same R/10 mesh (74 020 cells) read map-vs-map **9.797e-13 /
> 2.052e-12** threaded (`20260929T183957Z_MAT-1-probe-noise.log:105`) and
> **9.047e-13 / 5.962e-13** at OMP / OpenBLAS threads = 1
> (`20260929T184102Z_MAT-1-probe-omp1.log:106`), and the third solve
> **4.45e-11 / 1.93e-11** (threaded, `…probe-noise.log:106`) and
> **4.24e-11 / 2.14e-11** (threads = 1, `…probe-omp1.log:107`) off the first,
> through the phantom path. Identical inputs through the same assignment path do not
> reproduce to 1e-12, so the premise's "predicted 0.0" assumed a bit-exactness
> the `-n 2` complex MUMPS solve does not have — the same family as the
> 2026-09-10 `OPS-43` observation below (1-ULP cross-process drift, re-anchored
> at `-n 1`), here amplified to ≈1e-11 by the lossy-sphere conditioning. The
> drift's source inside the solve is not diagnosed; setting threads to 1 did
> not remove it.
>
> **Everything else in the item is green on W1** (`…183740Z_MAT-1-W1.log:485–489`):
> `material_map` vs closed form 3.536 % (pred ≈3.5 %); presets low / mid / high
> mean-SAR 3.523 / 3.548 / 3.557 %, Im/Re 0.404 / 0.835 / 1.499 % at |k_in|R
> 0.1640 / 0.2010 / 0.2341 (band 10 %); σ-blind separation 4.529 (ceiling 4.531,
> printed only). The wiring test with its dead `DummyMagnetostaticSolver`
> patch stripped passes (asserts unchanged).
>
> **Disposition belongs to the review:** the 1e-12 bound was pre-registered and
> is not moved in-slot. Candidates the measurement supports — a bitwise check
> on the assembled DG0 σ / εᵣ arrays (the path's actual claim) plus a solve
> bound above the measured ≈5e-11 floor, or the `OPS-43` `-n 1` re-anchor. The
> row stays ⚠️.
>
> **Ruled 2026-09-29 18:00 review.** The premise's physics claim is that the
> phantom branch *assembles the same material* as `material_map`; that is
> checkable exactly, upstream of MUMPS (`build_material_fields`,
> `core/time_harmonic.py:243–289`). It is re-anchored as bitwise DG0 σ / εᵣ /
> μᵣ equality on every owned dof plus an operator identity
> (‖A_φ − A_map‖_∞ == 0.0, same pattern), with a 1-ULP σ control and a 1e-9
> σ control each reading > 0 (§9 item 83; `MAT-1` Done-when amended). The
> solve-level difference becomes a printed record; **the 1e-12 bound is not
> moved** — it is withdrawn as a gate on a solve output, not widened. `-n 1`
> is declined: at `-n 2` the bitwise compare exercises the rank-local tag
> loop, ghosts and `scatter_forward`, and `OPS-43`'s `-n 1` precedent fitted
> a solve-output claim. This entry stays 🔴 OPEN until item 83 lands; the
> landing commit retires it with the (P1)/(P2) readings. The floor data are
> copied to the `OPS-43` observation below so they outlive this entry.

### ✅ RETIRED 2026-09-29 (opened and retired in the same commit, `PORT-7`, 22:30 CDT slot, §9 item 64) — `export_touchstone` wrote N ≥ 3 S-matrices column-major on one line, so a Touchstone 1.x reader got Sᵀ; `load_touchstone` rejected spec-layout N ≥ 3 files

| | |
|---|---|
| **Tests** | `tests/io/test_touchstone_skrf_crosscheck.py::test_our_writer_read_by_skrf[c_nonrecip_3port]` / `[c_nonrecip_4port]`, `::test_skrf_writer_read_by_load_touchstone[*]` except `c_nonrecip_2port`, `::test_transpose_negative_control` |
| **Symptom** | Pre-fix W1b (`20260929T033200Z_PORT-7-W1b.log`, commit `8ac4279d`): `assert 0.43924608386804553 <= 1e-12` (N = 3), `assert 0.5658913216122953 <= 1e-12` (N = 4), with skrf's read matching Sᵀ to 4.7e-14 / 6.0e-14; `load_touchstone` on skrf-written N ≥ 3 files: `ValueError: setting an array element with a sequence. The requested array has an inhomogeneous shape`. |
| **Cause** | `ports/touchstone.py` flattened every N with `order="F"` on one line and read back the same way, so the in-house round trip cancelled; Touchstone 1.x lists N ≥ 3 row by row, one matrix row per line. Invisible on reciprocal data — the stored `PORT-15` records passed writer → skrf through the bug (max\|S − Sᵀ\| 7.9e-15). |
| **Fix** | Writer N ≥ 3 row-major, one matrix row per line, ≤ 4 pairs per line (N ≤ 2 untouched); loader collects tokens across lines, one record per 1 + 2N² values, N ≥ 3 in C order. `20260929T033225Z_PORT-7-W2.log`: 19 passed, `test_touchstone_export.py` green unedited. No repo caller of `export_touchstone` and no tracked `.sNp`, so no record moved. |

### ✅ RETIRED 2026-09-29 (06:00 slot, §9 item 74) — was OPEN 2026-09-29 (filed by `TH-16` step 1, 22:30 CDT slot, §9 item 65; found while drafting by the 2026-09-28 18:00 review) — `analytic_cavity_frequencies` counts every (m, n, p) once, so degenerate TE / TM pairs appear once

> **Retired.** The helper now delegates to `analytic_box_cavity_modes`, and 346.40 MHz is counted twice. The docstrings are corrected (353.53 MHz is the sixth mode).
> - W1 on the old `src/` was red as predicted: (i) failed at rtol 1e-14, and (ii) read 9.36 % worst (`20260929T112807Z_TH-9-KI-W1.log:208,221`).
> - W2 after the fix: 28 passed, and the six modes read ≤ 0.0517 %. The old helper, frozen as the control, counts the pair once and has a 2.06 % mode-5 gap. `TH-9`'s original three tests and `TH-16`'s tests are green (`20260929T112910Z_TH-9-KI-W2.log`).

| | |
| --- | --- |
| **Symptom** | `core/cavity.py::analytic_cavity_frequencies` appends each index triple once. On `TH-9`'s 1.0 × 0.8 × 0.6 box TE₁₁₁ / TM₁₁₁ (346.40 MHz) is a double eigenvalue, returned once, so its 5th value 353.53 MHz is really the 6th mode, and `test_cavity_resonances.py:13–14`'s docstring "fifth mode at 353.53 MHz" is wrong. The solver returns both 346.40 MHz modes (346.5470 / 346.5627 MHz on the 1.0 × 0.8 × 1.2 box, `20260929T033718Z_TH-16-step1-W1.log:109–110`). |
| **Verified at** | `d3073242` plus `TH-16` step 1: `analytic_box_cavity_modes(edges, ())` returns 346.40 MHz twice (the step's numpy multiset test). |
| **Cause** | Multiplicity rule missing: two independent polarisations when all three indices are non-zero, one when exactly one is zero. |
| **Scope** | `TH-9`'s four-mode gate is unaffected (the first four values are simple). `solve_pec_cavity_modes` takes its shift-invert target from the helper, which shifts the target only once `n_modes` ≥ 5. Do not use the helper past mode 4; `analytic_box_cavity_modes` (`TH-16` step 1) counts multiplicity. |
| **Resolves with** | A chunk that points `analytic_cavity_frequencies` at the multiplicity-aware helper and corrects the docstring, with `TH-9` re-run green. Not done in `TH-16` step 1 (scope: additive, `solve_pec_cavity_modes` untouched). |

### ✅ RETIRED 2026-09-30 (18:00 daily review, on §9 item 112's six logs and a `log-pathologist` ruling: located and mitigated, not fixed — the drift needs ≥ 2 ranks; see the last two rows) — was: The `B_DG1` midpoint relative median is not run-to-run reproducible at `-n 2` — both `POST-4` probes reprint it with a ≈1e-6–4e-6 relative spread on a bit-identical mesh, so the `OPS-50` step-2 1e-6 cross-instrument anchor fails on B and holds exactly on A and E (2026-09-19, `OPS-50` step 2, §9 item 4)

| | |
| --- | --- |
| **Symptom** | On the 0.11 image, complex build, `-n 2`, same mesh (`cells=9291`, fingerprint `m1=6.071004883645e+01 m2=1.297311089112e+02` identical in every window), the `B_DG1` midpoint relative median is reprinted differently by each run, while `A_N1curl` and `E_N1curl` are bit-identical across every window. Rows, as printed: **step-4 probe run 1** `B_DG1 mid_rel_med=4.717160e-01` (`20260919T110204Z_OPS-50-step2-probe4.log:399`); **step-4 probe run 2**, same command, 34 s later, `B_DG1 mid_rel_med=4.717154e-01` (`20260919T110258Z_OPS-50-step2-probe4-repeat.log:399`) — self-spread **1.27e-06**; **step-5 probe 2026-09-18** `REPRO B p1_mid_rel_med=4.717157e-01` (`20260918T093644Z_OPS-50-probe.log:400`); **step-5 probe 2026-09-19** `REPRO B p1_mid_rel_med=4.717141e-01` (`20260919T110225Z_OPS-50-step2-probe5.log:400`) — self-spread **3.39e-06**. Cross-instrument (step-4 run 1 vs step-5 today) **4.03e-06**, against the `OPS-50` step-2 asserted anchor of **1e-6**. A: `4.532338e-01` in all four. E: `2.175825e-01` in all four. The step-5 probe's `RT_DOF B field_max` moves with it (`1.680772e-05` → `1.680775e-05`); `RT_DOF A` and `E` do not. |
| **Verified at** | `77a761f`, by execution: the four windows named above, 4–5 s each, `-n 2`, complex build, standard tier. |
| **Cause** | **Not diagnosed.** The scale is right for a non-deterministic reduction or point-ownership/ordering path that only `B` takes — `B` is the one field of the three that is assembled as a `DG1` projection of `curl A` rather than read from an `N1curl` solve vector, so its dof values are the output of a rank-partitioned assembly whose summation order is not pinned at `-n 2`; a 400-point *median* then lands on a different sample element when two neighbouring points swap. A/E reproducing bit-exactly is consistent with that and rules out the evaluation and `.bp` round-trip paths (both fields go through them unchanged, and `RT_DOF A/E max_abs_diff=0.000000e+00` in every window). Not established by measurement — no single-rank or repeated-seed window was run. |
| **Not caused by** | The step-5 probe's P1 path diverging from the step-4 path: the two instruments agree to **0.00e+00** on A and E on the same mesh, and their B difference (4.03e-06) is inside the 1.27e-06–3.39e-06 self-spread each instrument shows on B alone. Not the fixture drift of the entry above either — that is 7.8–11.4 %, five orders larger, and both instruments read it identically. |
| **Scope** | The printed `B_DG1` / `p1_mid_rel_med` figure in the two `POST-4` probes. No test and no gate reads it. It bounds how tightly any future re-record of `STEP4_MID_REL_MED[B]` can be stated: **≈4e-6 relative is the floor**, so a version-tagged re-record must carry that as its reproducibility, not the printed 7 digits. |
| **Resolves with** | A one-window measurement, not yet scheduled: run `post4_step4_probe.py` three times at `-n 1` and three at `-n 2` and tabulate the `B_DG1 mid_rel_med` spread at each width. If the spread collapses at `-n 1`, it is the rank-partitioned assembly and the fix is a pinned reduction (or recording B's median to 5 digits); if it survives at `-n 1`, the hypothesis above is refuted and the locus is elsewhere. `OPS-50` stays 🟡 until a review rules on the re-record. **2026-09-20 03:00 review:** ruled — the re-record is licensed with B at five significant digits (`OPS-50` step 3, §9 item 15); **this entry stays open** after it. One further reading from the same logs: the *source* field's own `rms_src` for `B_DG1` MID differs between the two step-4 windows (`3.739430e-06` at `…probe4.log:385`, `3.739427e-06` at `…probe4-repeat.log:385`) while A's and E's are identical to the printed digit — so the spread is already in the projected `B` before either probe samples it, which is consistent with the hypothesis above and narrows the locus to the `curl A → DG1` projection solve. The `-n 1` / `-n 2` window is still unscheduled: no status rides on it (rubric element 7). **2026-09-21 (`OPS-50` step 3b):** `OPS-50` closed ✅ and the sibling probe-pin entry retired; B's reprinted median this slot, `4.717148e-01` / `4.717147e-01` (`20260921T111404Z_OPS-50-step3b-prechange.log:400`, `20260921T111436Z_OPS-50-step3b-probe.log:400`), sits inside the spread above; **this entry stays open**. |
| **Measured 2026-09-30 13:30 slot (§9 item 112) — the spread collapses at `-n 1`; retirement held for the review** | Six windows of `post4_step4_probe.py`, the recorded command with only the rank count changed. The mesh fingerprint (`cells=9291`, `m1=6.071004883645e+01 m2=1.297311089112e+02`) and `CTRL_P1` 0 are identical in all six. **`-n 1`** (`20260930T183717Z_…-n1-r1.log`, `…183728Z_…-n1-r2.log`, `…183741Z_…-n1-r3.log`, Status 0): `B_DG1 mid_rel_med` reads 4.576852e-01 ×3, and the whole `DISAGREE B_DG1 MID` line is identical, `rms_src` 3.739380e-06 ×3. **Spread 0.** **`-n 2`** (`…183752Z_…-n2-r1.log`, `…183758Z_…-n2-r2.log`, `…183803Z_…-n2-r3.log`): 4.717148e-01 / 4.717157e-01 / 4.717137e-01, spread **4.2e-6** relative, and `rms_src` 3.739429e-06 / 3.739428e-06 / 3.739428e-06. A and E are bit-constant **within** each width: 4.290218e-01 / 2.046582e-01 at `-n 1`, and 4.532338e-01 / 2.175825e-01 at `-n 2` (this entry's values). The `-n 2` footers read Status 1 from the probe's own `FAIL PIN E_N1curl` (vertex 2.126638e-01 above midpoint 1.731650e-01). That pin already read FAIL with the same figures on 2026-09-20 (`20260920T125627Z_OPS-50-step3-probe.log:418`), so it is a pre-existing probe pin, not a run failure. **Reading:** the item's first branch holds on the numbers (`-n 1` spread 0, `-n 2` > 0), which locates the drift in the rank-partitioned `curl A → DG1` projection, and B's five-digit record (`OPS-50` step 3b) is the mitigation. **Why not retired in-slot:** the item's instrument anchor asks for A = 4.532338e-01 and E = 2.175825e-01 "in all six" windows. At `-n 1` they read otherwise, because the P1 comparator's vertex pick is partition-dependent (`POST-4`, the item's own trap). They are constant within each width, so the instrument did not change run to run. The anchor's literal text is still missed, and the review rules whether that was a drafting slip. **If so, retire this entry on these six logs; no re-run is needed.** |
| **Ruled and retired 2026-09-30 18:00 daily review** | **The anchor was a drafting slip** (the 10:30 review's): its only source was a `-n 2` log, and the item's own trap said the `-n 1` P1 median need not equal the `-n 2` one. A `log-pathologist` ruling on the six logs: claims 1–4 CONFIRMED (`-n 1` B line identical ×3, `…-n1-r*.log:382,396`; `-n 2` spread 4.24e-6, `…-n2-r*.log:386,400`; fingerprint and `CTRL_P1` identical ×6; A / E at `-n 2` equal the 09-19 record to every printed digit, `rms_src` A included) and nothing in the logs shows the instrument changed; `-n 1` constancy is shown only within these three runs (no earlier `-n 1` log of this probe exists). **What the retirement says, precisely:** the run-to-run drift needs ≥ 2 ranks; the `curl A → DG1` projection is the best-supported locus, not a proven one (A's seven-digit stability bounds an A-solve drift only to ≈ 1e-7, the same order as B's `rms_src` movement). It still occurs at `-n 2`: **located and mitigated, not fixed.** Two corrections to the Scope row: pooling every `-n 2` print on record (4.717160e-01 … 4.717137e-01) gives **4.9e-6**, not ≈ 4e-6; and B does not reproduce at five significant digits literally (today's three round to 4.7171 / 4.7172 / 4.7171) — the record 4.7172e-01 is guarded by `PIN_REPRO_RTOL = 0.02` (`post4_step5_probe.py:123,137`), which is what the mitigation actually rests on. The Status-1 footers' exact precedent is this probe's own `20260919T110204Z_OPS-50-step2-probe4.log.gz:408`; the `FAIL PIN E_N1curl` pin does not fire at `-n 1`, so the pin itself is rank-count-dependent. A pinned reduction is a design note, not an open defect. |

### ✅ RETIRED 2026-09-30 (`PORT-5` ✅, 04:30 slot, §9 item 98) — 3. Port tests assert a non-zero power wave where the placeholder's fake makes it exactly zero (re-dated 2026-08-28 by `OPS-28`)

> **Retired 2026-09-30.** The last red name,
> `test_sparameter_assembly.py::test_n_port_sweep_assembles_finite_matrix_with_expected_shape`,
> is replaced by `…::test_n_port_sweep_assembles_the_known_s_matrix` (three
> typed S₀, among them an ideal circulator) and the closed-form metric
> tests: 14 passed, 2 s, `-n 1`, real build
> (`20260930T095825Z_PORT-5-W1.log:53–137,158`). The matched-port fake whose
> diagonal was legitimately zero is gone; nothing was tuned to the
> heuristic. Everything below is history.

> **Re-symptomed 2026-09-29 (`PORT-8` ✅, 06:00 slot).** The orientation
> name is gone: `PORT-8` replaced both placeholder tests in
> `tests/ports/test_port_orientation_sensitivity.py` with one real-solve
> test, reversing the lumped-sheet drive on P2/P4 ⇒ S′ = DSD
> (`20260929T110818Z_PORT-8-W1.log`, 13 passed). The entry stays open
> for `test_sparameter_assembly.py` alone, the 3-port zero-diagonal fake.
> The orientation rows below are history.

> **Re-dated 2026-08-28 (`OPS-28`, 22:30 slot) — still open, and the
> *reason* is confirmed while the entry's old one-line statement was
> imprecise for one of the two names.** With the `_DummyComm` double
> repaired, `…::test_port_orientation_flip_changes_off_diagonal_sparameter_sign`
> reaches its S-matrix assertion for the first time since `OPS-14` and dies
> at `tests/ports/test_port_orientation_sensitivity.py:115`,
> `assert aligned_s21.real > 0.0` → `assert np.float64(0.0) > 0.0`. On that
> fixture the **diagonal is not zero at all** — the log prints
> `S11 = 9.047e-01 − 1.289e-02j, S22 = 9.047e-01 − 1.289e-02j` — it is the
> **off-diagonal** that vanishes, because the placeholder gives the
> *undriven* port `V = 5.000000e-02 V`, `I = 1.000000e-03 A` at
> `Z₀ = 50 Ω`, i.e. `V = Z₀I` exactly, so `b = (V − Z₀I)/(2√Z₀) = 0` and
> `S21 = S12 = 0` identically. Same mechanism as the 3-port diagonal case
> below, different matrix entry, so the entry's title is corrected above
> and both names stay filed. `OPS-28`'s scope was the double only — no
> assertion moved, no `sparameters.py` edit, `git show -- src/` empty.
> Log: `20260828T033055Z_OPS-28-gate.log` (`2 failed, 15 passed in 0.79s`,
> Status 1, elapsed 2 s, `-n 2`, real build, smoke). Verified at the
> `OPS-28` commit. Disposition unchanged: `PORT-0`/`PORT-1` own it.

| | |
|---|---|
| **Tests** | `tests/ports/test_sparameter_assembly.py::test_n_port_sweep_assembles_finite_matrix_with_expected_shape` (zero **diagonal**, 3-port fake)<br>`tests/ports/test_port_orientation_sensitivity.py::test_port_orientation_flip_changes_off_diagonal_sparameter_sign` (zero **off-diagonal**, 2-port fake — measured 2026-08-28, see the note above) |
| **Symptom** | `assert np.all(np.abs(diagonal) > 0.0)` fails on `array([0.+0.j, 0.+0.j, 0.+0.j])`; on the orientation test, `assert aligned_s21.real > 0.0` fails on `np.float64(0.0)` |
| **Cause** | The fakes set `current = voltage / port.z0_ohm` — a perfectly matched port. The power wave `b = (V − Z₀I)/(2√Z₀)` is then `0` exactly, so the corresponding S entry is *legitimately* zero and the assertion cannot hold: on the 3-port fake that is every diagonal term, on the 2-port orientation fake it is the off-diagonal terms (the undriven port is the matched one there). |
| **Fix** | **Deliberately not fixed.** These exercise the placeholder coupling model's arithmetic (see `PORT-0`). Repairing them means tuning assertions to match a heuristic that `PORT-1` deletes. Resolve them there. |
| **Verified pre-existing at** | `53f6428` and earlier |
| **Progress 2026-08-04** | `PORT-1` step 3b-ii, the step that would replace these fakes with a driven gap port, was **attempted and parked** (`attempt/PORT-1-step3bii-20260804T141200Z`). The drive works — reciprocity `2.2840e-04`, undriven port open at `2.32e-03` — but `Im Z₁₂` is +72.12% off `ωM₁₂`, traced to the gap box's 1.83×-oversized cross-section rather than to the solve. These two tests stay red and unchanged; see PROJECT_PLAN §7 `PORT-1` step 3b-ii for the measurement and the ranked successor. |
| **Progress 2026-08-04 (2)** | Step 3b-iii tested that diagnosis and **refuted it** (parked on `attempt/PORT-1-step3biii-20260804T173000Z`; logs `20260804T170301Z_PORT-1-step3biii-costprobe.log`, `20260804T170439Z_PORT-1-step3biii-sweep-o5e4.log`). Shrinking only the transverse overhang gives `Im Z₁₂/ωM₁₂` = `+1.7210` (fringe 0.4546), `−0.2391` (0.3509), `+0.3317` (0.2739) — non-monotone and **sign-changing**, so a box-volume average is not a port voltage at any overhang. The shadow-restricted average is stable at `0.687–0.814 ×` across the same three geometries. The replacement route is now the facet-integral voltage (§7 step 3b-v on 3b-iv's tags); these two tests stay red and unchanged. |
| **Progress 2026-08-06** | Step 3b-v ran that facet-integral route and it is **also excluded** (parked on `attempt/PORT-1-step3bv-20260806T004500Z`; log `20260806T003559Z_PORT-1-step3bv-gate.log`, 3 failed / 7 passed, 67.6 s at `-n 2` — no hang, 3b-iv's hoisted `create_entity_permutations()` held). `V = −⟨E·ŷ⟩_{terminal discs, gap side}·L` gives `\|Im Z₁₂\|/ωM₁₂` = **4.845** (+384.54%) against the full-box `0.332` and tube-shadow `0.763 / 0.814` **on the same solve**, with reciprocity degrading `1.15e-4 → 1.79e-2`. Cause: `E·ŷ` on a terminal is the surface-charge-dominated *normal* component (measured gap/wire jump ratio 2.9e-5–4.6e-5), so a two-endpoint trapezoid samples exactly where the integrand peaks. Both the box family and the terminal-facet route are now excluded by measurement; a successor must integrate the *tangential* component along the whole gap path. `MUTUAL_TOLERANCE` unmoved. These two tests stay red and unchanged. |
| **Progress 2026-08-06 (2)** | Step 3b-vi ran the tangential path integral — `V = −∫E·t̂ dl` along the torus centreline arc through the gap, Gauss–Legendre nodes strictly interior to the arc, sampled via `post.evaluation.evaluate_vector_field_parallel` — and it is **unresolved, not refuted** (parked on `attempt/PORT-1-step3bvi-20260806T094500Z`, `ee5f0cb`; logs `20260806T093603Z_PORT-1-step3bvi-gate-n2.log` and `20260806T093808Z_PORT-1-step3bvi-quadrature-sweep-n2.log`, 4 failed / 4 passed, 136.13 s at `-n 2`). Value, off one solve: path **0.468933 / 0.499728** × ωM₁₂ against facet `4.802 / 4.889`, full-box `0.332`, tube-shadow `0.763 / 0.814` — a *third* distinct answer, below the shadow family, at −51.6% with reciprocity `6.3e-2`. But the plan's own quadrature precondition fails first: the proposed `(33, 65)` node pair disagrees by `1.07e-1`, and the sequence extended to 4097 nodes converges only at `O(1/n)` to a `~1e-3` plateau — because N1curl makes only the *facet*-tangential component continuous, the arc's tangent is not facet-tangential, and with `h_wire = 2.5e-3` against arc length `a·g = 1.2e-2` only ≈5 cells span the path. Four estimator families now give four answers spanning a factor 15 on one solved field. The successor is **arc refinement, not more quadrature nodes**; if `~0.48` survives it, the family question is closed negatively and the suspects become finite-σ terminal penetration or the `ωM₁₂` reference. `MUTUAL_TOLERANCE` unmoved. These two tests stay red and unchanged. |
| **Progress 2026-08-06 (3)** | Step 3b-vii ran that arc refinement and **the family question closes negatively** (parked on `attempt/PORT-1-step3bvii-20260806T170000Z`, `bc8c04e`; probe `20260806T170559Z_PORT-1-step3bvii-probe.log`, gate `20260806T170835Z_PORT-1-step3bvii-gate-n2.log`, 2 failed / 10 passed, 165 s at `-n 2`). A coordinate-defined `MathEval`+`Threshold` size field on the fragmented model put `h_gap = 3e-4` — 40 cells across `a·g` — in a tube around each gap arc: 124 753 → 178 055 cells (1.427×), gap-tagged cells 1569 → 24 430 per port. It fixed the discretization and not the answer. Reciprocity **6.3e-2 → 3.8823e-3** (inside the 1e-2 band for the first time) and the fixed-order quadrature residual improved ≈3×, but (129, 257) still disagree by `1.1444e-2` and the converged path value is **0.493653 / 0.491744** × ωM₁₂ (0.4808 at 4097 nodes) against 3b-vi's 0.468933 / 0.499728 — a discretization-level shift, `Im Z₁₂` at −50.73%. The built-in control holds: facet `5.165/5.169`, box `0.3496/0.3492`, shadow `0.8566/0.8386` all moved only a few % off the refined solve, so the solve did not change underneath the estimator. Four sampling geometries, four answers, and the ≈0.48 survives refinement: **the deficit is not the sampling geometry**. Remaining suspects are finite-σ terminal penetration at `δ = 1.125 r_wire` and the filamentary `ωM₁₂` reference itself — a review's adjudication. `MUTUAL_TOLERANCE` unmoved. These two tests stay red and unchanged. |
| **Progress 2026-08-07** | Step 3b-viii audited the **`ωM₁₂` reference** — the first of the two remaining suspects — in closed form, no solve and no mesh (`tests/validation/test_mutual_inductance_reference.py`, on `main` and green; log `20260807T020314Z_PORT-1-step3bviii-gate.log`, 7 passed in 0.43 s at `-n 1`, smoke). **The suspect is retired.** (i) The vector-potential route the gates use and an independent elliptic-integral reimplementation (`M = μ₀√(ab)[(2/k − k)K − (2/k)E]`, `m = k²`) agree to **1.5e-15** at the fixture's `d`, and to ≤ 7.5e-14 across `d/4 … 4d`; `ω·M` reproduces the logged `1.241755 Ω` to 3.1e-7. A vacuity control confirms the identity has teeth — feeding SciPy the modulus where the parameter belongs moves `M` by **140%**, eleven orders above the 1e-9 gate. (ii) The finite-cross-section correction at `r/a = 0.125`, the filament kernel averaged over both minor discs at uniform current density (Gauss–Legendre × periodic trapezoid, converged to 6.7e-16 between orders), is **`M_tube/M_fil = 1.004809992`, i.e. +0.481%** — and it has the *wrong sign* to help, since a larger reference makes the deficit slightly worse (0.4937 → 0.4914 × ωM). A 0.5% correction cannot produce a factor 2, consistent with step 2's field-level agreement at −9.35%. **The reference is exonerated; finite-σ terminal penetration (step 3b-ix) is now the only named suspect.** `MUTUAL_TOLERANCE` unmoved. These two tests stay red and unchanged. |
| **Progress 2026-08-07 (2)** | Step 3b-ix closed the loop and **found the cause of the factor 2 — it is not physics, it is the estimator's integration limits.** Parked on `attempt/PORT-1-step3bix-20260807T050000Z` (`6caec85`); log `20260807T050654Z_PORT-1-step3bix-gate-n2.log` on that branch, 178 055 cells, `-n 2`, 227 s. Faraday on the closed centreline circle, tiled into four segments each verified against the DG0 material indicator before any solve (0 misassigned of 5392 nodes), gives on the undriven port: `V_gap = 0.493653 / 0.491744`, **`V_buried = 0.399972 / 0.402239`**, `V_wire = 0.002394 / 0.002316`, **sum = 0.896019 / 0.896299 × ωM₁₂** — inside the pre-set `1 ± 0.15` band and matching step 2's independent reaction route (−9.35%, of which −9.36% is the PEC box at padding 0.08). **Finite-σ terminal penetration is retired**: at `σ × {1, 2, 4}` (δ/r_wire 1.125 → 0.796 → 0.563) `V_wire/ωM` = 0.002394 → 0.001856 → 0.000727 — the signature is real and 200× too small to matter — while `V_gap/ωM` = 0.493653 → 0.490837 → 0.485059, i.e. it *falls*. The missing half is the **buried** dielectric: `GAP_BURIAL` puts the gap region out to `±arcsin(half_y/a) = ±0.175335` rad while `_gap_arc_quadrature` integrates only the nominal `±GAP_ANGLE/2 = ±0.15`, so 1.013 mm of arc per side — 0.8% of the loop's length, and exactly where the terminal fields live — carries 45% of its EMF. **Terminal to terminal the gap-port voltage is 0.8936 × ωM₁₂, not 0.4937.** Both suspects the 2026-08-06 review named are now dead. Successor **step 3b-x** (correct the limits off the port facet tags, then a box-padding sweep — the residual −10.6% is the PEC box's, and a `MUTUAL_TOLERANCE` edit is *not* the remedy) is written into PROJECT_PLAN §7 for a review to scope. `MUTUAL_TOLERANCE` unmoved at 0.10; nothing under `src/` changed. These two tests stay red and unchanged. |
| **Progress 2026-08-07 (3)** | Step 3b-xii ran the pre-decided box discriminator and it lands on **disposition (ii): the residual is a real estimator bias, not the PEC truncation box.** Parked on `attempt/PORT-1-step3bxii-20260807T170000Z` (`87bf35d`, carrying the full 3b-ix → 3b-x-b lineage); probe `20260807T170143Z_PORT-1-step3bxii-probe.log` (59 s), gate `20260807T170430Z_PORT-1-step3bxii-disc-n2.log` (`-n 2`, standard, **353 s**, 5 passed + the discriminator red). The gapped fixture was rebuilt at `air_padding = 0.10` — nothing else moved, and the same `_solve_gap_ports` the 0.08 gates use — at 194 985 cells (1.0951× the 178 055 at 0.08, under the plan's 230 000 stop rule; padding 0.08 re-meshed at **exactly** 178 055, so the fixture identity holds). All four route values: at padding 0.08 estimator `0.894543 / 0.894022`, control `0.922423`, deviation `−3.0224e-02`; at padding 0.10 estimator **`0.924103 / 0.923075`**, control **`0.952868`**, deviations **`−3.0188e-02 / −3.1267e-02`**. The box moved the estimator **+2.956 pp** and the control **+3.045 pp** — i.e. it moved *both routes together* — leaving their difference at 3.02–3.13% against the pre-decided 2.5% threshold, a move of **−0.104 pp**: the wrong direction, and 5× smaller than the 0.5 pp (i) required. That the box itself behaved normally is corroborated independently: this fixture's σ = 0 control reads `0.952868` at padding 0.10 against 3b-xi's *ungapped* reaction route at `0.949744`, and `0.922423` against `0.919676` at 0.08 — a stable `+0.27/+0.31` pp gapped/ungapped offset under enlargement. Negative control, recomputed against this box's own reference: the uncorrected wedge-only estimator gives ratio `0.5181`, deviation `−0.4819`, 15× the threshold. **`REACTION_CONSISTENCY_TOLERANCE` stays at 0.03** — the review's authorized re-size to 0.05 was conditional on convergence, which did not happen, so it is not taken, and nothing under `src/` changed. What this leaves open, for a review: the corrected terminal-to-terminal estimator and the σ = 0 reaction control differ by ≈3% *for a reason that is not truncation and not the wedge limits* — the two problems differ in that the production loop is gapped and σ = 800 S/m while the control's is closed and lossless, and that difference is now the named suspect rather than an assumed nuisance. `MUTUAL_TOLERANCE` unmoved at 0.10. These two tests stay red and unchanged. |
| **Scoped 2026-08-07 (03:00 review)** | Both successors are now full §7 plans and queued (§9 items 1–2). **3b-x**: correct `_gap_arc_quadrature` to terminal-to-terminal limits read off the port facet tags, gate the corrected estimator against the *same-fixture* reaction-route `Im Z₁₂` (≤ 3%; the ωM₁₂ comparison stays printed and tracked here, expected ≈−10.6%), delete the σ-monotonicity test whose premise 3b-ix refuted, and land the branch. **3b-xi**: padding sweep {0.08, 0.10, 0.12} on the ungapped fixture — if the deficit does not shrink monotonically with padding, the PEC-box attribution dies and this entry escalates to the weekly review. `MUTUAL_TOLERANCE` untouched by either. |
| **Progress 2026-08-08** | Step 3b-xiii ran that σ ladder and lands on **disposition (mixed) — but the informative result is that the experiment's premise is disproved, not that the bands were missed.** Parked on `attempt/PORT-1-step3bxiii-20260808T005500Z` (`82bfb40`, carrying the full 3b-ix → 3b-xii lineage); log `20260808T004346Z_PORT-1-step3bxiii-ladder-b-n2.log` (`-n 2`, standard, **344.6 s**, 20 passed + the known consistency gate red). Fixture identity first and it byte-reproduces the branch's record **exactly** — estimator `0.894543 / 0.894022`, control(σ = 0) `0.922423`, deviation `−3.0224e-02` — so nothing geometric moved. The ladder, σ on the wire ∪ gap-box footprints of both loops through the same DG0 material map the production solves use: `control(σ=0) = 0.922423`, `control(σ=200) = 0.496614`, `control(σ=800) = 0.107556` × ωM₁₂. It is **monotone decreasing** (the new ordering gate — the intermediate rung must lie between the endpoints, or the ladder measures noise rather than σ — passes) but lands nowhere near either 0.7 pp band: the σ = 800 rung sits **78.7 pp** from the estimator and **81.5 pp** from control(σ = 0), on a 2.81 pp endpoint spread. **Why:** a *closed* lossy loop is a shorted turn. The induced circulating current reaches `|I_cond/I′| = 0.412` at σ = 200 and **0.865** at σ = 800, and its back-field cancels most of the mutual EMF the reaction integral reads. σ and closed-vs-gapped are therefore **not separable on this control** — the knob the review believed independent is confounded with the very difference it was meant to isolate. The ≈3% deviation is **untouched**: the loss-vs-gap question stays open and now needs a *gapped* control (the production loop driven at σ → 0), not a lossy closed one; that escalation is the weekly review's. **Not tuned:** `REACTION_CONSISTENCY_TOLERANCE` stays 0.03, nothing re-pointed, `MUTUAL_TOLERANCE` unmoved at 0.10. One `src/` change, unrelated to `PORT-1` and parked with the branch: `_validate_material_map_tags` tested rank-local `cell_tags.values`, so a material map over the two 1 mm gap boxes — globally valid — raised on one rank of two while the other entered the solve and hung in the first collective until the ceiling (`20260808T003238Z_PORT-1-step3bxiii-ladder-n2.log`: 16 errors, a 246 s session costing 601 s). The tag set is now reduced with `mesh.comm.allgather`. These two tests stay red and unchanged. |
| **Scoped 2026-08-07 (18:00 review)** | 3b-xii's disposition (ii) successor is now a full §7 plan and queued (§9 item 1): **step 3b-xiii**, the σ ladder {0, 200, 800} S/m on the closed-footprint control — nothing geometric moves, so σ becomes the only difference between the routes and the ≈3% deviation gets an owner: loss or gap. Landing the parked branch (`attempt/PORT-1-step3bxii-20260807T170000Z`) is pre-authorized only under the (loss) disposition with the consistency bound unchanged at 0.03; the (gap) and (mixed) outcomes park and escalate to the weekly review. `attempt/PORT-1-step3bxb-…` was deleted this review as strictly superseded by content (+852/−3, verified — the ancestry test returns false only because the 12:00 run squashed the lineage). `MUTUAL_TOLERANCE` untouched regardless of outcome. |
| **Scoped 2026-08-08 (03:00 review)** | 3b-xiii's successor is now a full §7 plan and queued (§9 item 1): **step 3b-xiv**, the non-degenerate half of the same sweep — hold the **production gapped** fixture fixed and run its σ down {800, 200, 0}; a lossless *gapped* loop has no shorted-turn current, so σ and gap separate in this direction. **Measurement only: every disposition parks and reports** — branch landing and any gate re-pointing stay with the weekly review (2026-08-10), per 3b-xiii's escalation; this step exists so that review adjudicates with the ladder in hand. Branch hygiene: `attempt/PORT-1-step3bxii-20260807T170000Z` deleted this review (verified strict ancestor of the 3b-xiii branch); `attempt/PORT-1-step3bxiii-20260808T005500Z` (`82bfb40`) is the one live lineage. The rank-safety `_validate_material_map_tags` fix that rode along is separately scoped for `main` as `OPS-13`. |

### ✅ RETIRED 2026-09-29 (`PORT-4` ✅, 12:00 slot, §9 item 77) — 6. Rank-dependent: single-port excitation

> **Retired.** `PORT-4` replaced the placeholder test named below with one
> real-solve `@complex_only` test on the 116 085-cell birdcage: unequal
> terminations (P2 200 Ω, P4 25 Ω) ⇒ ‖Y − Yᵀ‖_F/‖Y‖_F = 8.513e-15 ≤ 1e-12
> (`20260929T170418Z_PORT-4-W1.log:1900`, 16 passed, 91 s, `-n 2`). The test
> this entry names no longer exists. Defect (1), the rank-local tagging, survives
> only in `_build_test_problem`, which the three argument-validation tests
> still use. They raise before the tag check (`ports/excitation.py:243–276`),
> so the defect can no longer turn a test red. Defect (2) was fixed 2026-08-13.

| | |
|---|---|
| **Test** | `tests/solver/test_single_port_excitation.py::test_single_port_excitation_returns_finite_estimates` |
| **Symptom** | Passes at 1 rank, fails at `mpiexec -n 8` — **verified 2026-08-08**, and the failure is a `ValueError` raised out of `src/`, not a test assertion: `missing required port tags: [21, 22]` / `[12, 21, 22]` from `ports/definitions.py:99`, on 8/8 ranks. Also red at `-n 4` (`missing required port tags: [22]`, 4/4 ranks) — the entry's `-n 8` is the count that was tried, not the threshold. |
| **Cause** | **Diagnosed 2026-08-08 (09:00 run) by `OPS-14` — two independent defects, and neither one alone explains the red.** See the diagnosis row below. *(Original, superseded:)* Not diagnosed. Two rank-local bugs of this family were already fixed on 2026-07-27 (`tests/solver/test_two_torus.py` asserted on rank-local `cell_tags.values`; `test_helmholtz_v2.py` did a per-rank collision search). Suspect the same pattern — a quantity that is rank-local being treated as global. |
| **Tool** | `tests/mesh/helpers.py::global_cell_tag_set()` exists for the tag case. `post.evaluation.evaluate_vector_field_parallel()` for the point-location case. |
| **Verified pre-existing at** | `ce92e8c` and earlier |
| **Scoped 2026-08-08 (03:00 review)** | `OPS-14` (§7, queued §9 item 4) owns the diagnosis — the last never-diagnosed entry in this file. **Treat the symptom line above as unverified:** the last two never-diagnosed entries both turned out to be misrecorded (entry 2's assertion was wrong; entry 4's green CI signal was a partition artifact), so the first step is a three-way `-n 1/2/8` reproduction that captures the *actual* failure. Prime suspect: the fixture's `tags[cell_indices % 4]` over rank-local indices on a 12-cell cube — at `-n 8` ranks see strict subsets of the tag set. Pre-registered disposition: a defect wholly inside the `PORT-0` placeholder is **not fixed** (entry 3's "resolve in `PORT-1`" logic) — this entry would then be re-pointed at `PORT-1` rather than retired. |
| **Diagnosed 2026-08-08 (09:00 run, `OPS-14`)** | **Two defects, independent, and each is individually sufficient to keep this red at `-n 4` and `-n 8`.** Measured with `scripts/probes/ops14_rank_probe.py` at `-n 1/2/4/8` (logs `20260808T140412Z`/`140424Z`/`140426Z`/`140427Z_OPS-14-table-n{1,2,4,8}.log`; reproduction `20260808T140044Z`/`140055Z`/`140056Z_OPS-14-repro-n{1,2,8}.log`). **(1) The fixture** (`tests/solver/test_single_port_excitation.py:21-26`) builds `tags[cell_indices % 4]` over **rank-local** indices, so the *global* tag set is itself rank-count dependent: `{11,12,21,22}` at `-n 1/2`, **`{11,12,21}` at `-n 4`** (per-tag global cell counts `4/4/4/0`), **`{11,12}` at `-n 8`** (`8/4/0/0`). Tag 22 exists on no rank at `-n 4`, so the `ValueError` is *correct behaviour* — the mesh really lacks the tag. The same defect makes the placeholder's own output rank-dependent while green: `P1.I` = `3.000000e-03+4.263898e-04j` at `-n 1` against `4.000000e-03+5.685197e-04j` at `-n 2` (+33.3%), because `support` counts tagged cells. **(2) `ports/excitation.py:249`** hands rank-local `problem.cell_tags.values` to `validate_required_port_tags_exist`, so the check is not collective — a rank owning no cell of a terminal region raises while others return. Counterfactual B (fixture tags taken over *global* cell numbering — global per-tag counts then exactly `3/3/3/3` at every rank count) still raises on 4/4 ranks at `-n 4` and 8/8 at `-n 8`, which isolates (2) from (1). Counterfactual A (collective validator argument, fixture untouched) still raises at `-n 4/8`, which isolates (1) from (2). **Anchor — the cross-rank identity:** under counterfactual B the estimates are **byte-identical** at `-n 1` and `-n 2` (`P1.I = 3.000000000000e-03+4.263897544510e-04j`, `P2.V = 5.000000000000e-02`, coupling `1.000000000000e-01`), where production diverges by 33.3% across the same two rank counts. |
| **Defect (2) fixed 2026-08-13 (13:30 run, `PORT-1` step 4) — entry stays red** | `run_placeholder_port_coupling_case` now reduces the tag set with `comm.allgather` before `validate_required_port_tags_exist`, exactly as that function's docstring prescribes. Not a courtesy: step 4's negative control runs the retiring heuristic on the two-torus fixture at `-n 2`, where each rank owns one port's cells, and defect (2) made that call raise on one rank while the other returned. Verified by the control's own green run (`20260813T183606Z_PORT-1-step4-packagegate.log`, 7 passed 153.9 s). **Defect (1) — the fixture tagging over rank-local cell indices — is untouched, so this entry stays open and the test stays red at `-n 4`/`-n 8`.** |
| **Disposition 2026-08-08 — not fixed, re-pointed at `PORT-1`** | The pre-registered rule applies as written: both defects are wholly inside the `PORT-0` placeholder that `PORT-1` deletes — (2) is a line of `run_placeholder_port_coupling_case`, and (1) is a fixture that exists only to exercise it. A shared-machinery survey found no other non-collective tag read left in `src/`: `core/time_harmonic.py:162` was fixed by `OPS-13`, `post/sar.py:184` already reduces with `allreduce(..., op=MPI.MAX)`, `io/mesh.py:1711` with `SUM`. The only change landed is a hazard warning in `validate_required_port_tags_exist`'s docstring (behaviour unchanged) so `PORT-1`'s real caller does not repeat (2). **This entry leaves with `PORT-1`, like entry 3.** |

### ✅ RETIRED 2026-09-30 (04:30 slot, §9 item 96: the balance is now dimensionless against an optional `reference_speed_m_per_s`, `None` with no warning when none is given; 1.639e-4 at ω/\|k\| on `TH-6`'s wave, `20260930T094718Z_POST-2-balance-W1.log:134–138`; `mri:1` reads `status: OK` / `n/a` and its guide is corrected, `20260930T094943Z_POST-2-balance-mri1.log:407–425`) — was OPEN 2026-09-30 (`POST-2` W1, §9 item 92; filed by pre-registration) — `compute_field_consistency_diagnostics`'s `mean_balance_rel_diff` compares \|E\| in V/m with \|B\| in T, so it reads ≈ 1 and fires the "strongly imbalanced" warning **on every physical field**, including an exact one

| Field | Value |
|---|---|
| **Test** | `tests/validation/test_field_consistency_metrics.py::test_post2_diagnostics_read_the_lossy_plane_wave_closed_form` (printed record, outside the gate; the test passes) |
| **Literal symptom** | On `TH-6`'s plane wave at 24³: `mean_balance_rel_diff=0.999999962617 vs cf 1-\|k\|/omega=0.999999962610` and `'strongly imbalanced' warning present: True`; 12³ the same (`0.999999962635`, `True`) (`20260930T033548Z_POST-2-W1.log:116–117, 128–129`) — on a field whose `e_to_b_mean_ratio` matches ω/\|k\| to 0.016 %. |
| **Commit verified at** | `10e004da` (tests-only working tree, no `src/` change) |
| **Cause** | `post/consistency.py:54` forms `abs(e_mean − b_mean)/max(e_mean, b_mean)` across two different units; for any physical field \|E\|/\|B\| is of order a phase velocity (≈ 2.7e7 m/s here), so the metric is 1 − \|B\|/\|E\| ≈ 1 − 3.7e-8 and the 0.95 threshold is always exceeded. The `mri:1` guide (`examples/mri/01_coil_phantom_fields.md` ≈ `:222–238`) attributes its `1.000000` reading to "an ungated proxy drive" — that attribution is wrong: an exact solve reads the same. |
| **Scope** | Outside `POST-2`'s ✅ by its Done-when; `e_to_b_mean_ratio` and the span ratios are unaffected. |
| **Resolves with** | The review's call (not licensed by `POST-2`, tests only): drop the metric, or compare \|E\| against a velocity-scaled \|B\| (e.g. c·\|B\| or the medium's ω/\|k\|·\|B\|); and correct `mri:1`'s wording. **Ruled 2026-09-30 03:00 review, queued as §9 item 96:** the metric is replaced — dimensionless against an optional `reference_speed_m_per_s`, reported as `None` with no warning when no speed is given (a coil + phantom field has no single one); the 0.95 warning default does not move; `mri:1`'s guide and docstring are re-recorded. The entry retires with the commit that lands both. |

### ✅ RETIRED 2026-09-30 (04:30 slot, §9 item 95: tolerance floored at 4 ulp of the carrier, last point pinned to `stop`, steps under 16 ulp raise; probe `duplicates=0 missing=0` asserted, 2 000-case `Decimal` family 0 failures at ≤ 1.456 ulp, `20260930T094514Z_PORT-6-floor-W1.log:69–78`) — was OPEN 2026-09-30 (`PORT-6` W1, §9 item 86; non-blocking by pre-registration) — `plan_frequency_sweep` has a **validity floor on the coarse step**: below `step ≈ ulp(stop)/1e-9` (≈ 30 Hz for carriers under 268 MHz, ≈ 60 Hz for 268–537 MHz) an on-grid `stop` can be appended as a **near-duplicate** of the last lattice point

| Field | Value |
|---|---|
| **Test** | `tests/ports/test_frequency_sweep_planner.py::test_port6_regime_probe_sub_60hz_steps_printed_only` (printed record, never asserted; the test passes) — renamed `…::test_port6_regime_probe_sub_60hz_steps` and asserted by the retiring commit |
| **Literal symptom** | `PORT6_REGIME carrier=297200000.71 step=0.1 n=7 stop=decimal len=9 expected=8 duplicates=1 missing=0 min_gap=np.float64(5.960464477539063e-08) last_two=(297200001.40999997, 297200001.41)`; the same at `n=37`; summary `cases=48 duplicates=2 missing=0` (`20260930T020839Z_PORT-6-W1.log:69–71`). |
| **Commit verified at** | `8f161af3` (tests-only working tree, no `src/` change) |
| **Cause** | `_STEP_RELATIVE_TOL * step` (1e-10 Hz at a 0.1 Hz step) is below one ulp of `stop` (5.96e-08 Hz at 297 MHz). When the decimal-typed `stop` rounds one ulp away from the planner's float `start + n*step`, the endpoint test and the merge both see two distinct points 1 ulp apart. Only the decimal-stated `stop` form tripped (2 of 24); `stop = start + n*step` computed in float is clean (0 of 24), and no point is ever **missing**. |
| **Scope** | Steps < 60 Hz are outside `PORT-6`'s ✅ by its Done-when; the integer-Hz family (steps ≥ 100 Hz) is bitwise-exact on 500/500 cases (`:62`). |
| **Resolves with** | A merge/endpoint tolerance floored at a few ulp of the carrier, e.g. `max(1e-9*step, 4*np.spacing(stop))`, re-checked against this probe. Not licensed by `PORT-6` (tests only). **Queued 2026-09-30 03:00 review as §9 item 95**, with two additions the drafter's host emulation found necessary (predictions, not records): the merge must pin the last cluster to `stop` (else `grid[-1] != stop` in ≈ 1 of 2 000 sub-60 Hz cases), and steps under 16 ulp of the carrier must raise. The merge site, which this probe never reaches, is the larger hole (predicted 315 of 2 000 wrong lengths on a wider family). |

### ✅ RETIRED 2026-09-30 (13:30 slot, §9 item 111: the break decision goes through `_ladder_budget_allows`, an `MPI.MIN` allreduce of `remaining − 2·last_cost`; `test_ladder_budget_decision_is_collective` reads `[False, False]` on straddling inputs against the un-reduced control's `[True, False]`, and the gate's rungs are unchanged at 0.846150 / 0.916742 / 0.966977, `20260930T183517Z_GEO-15-budget-collective.log:4340–4347`) — was Latent (has not fired): rank-local ladder-budget break in `test_birdcage_conductor_sizing.py` can desync collectives

**Found:** 2026-08-16, 10:30 review audit of `GEO-15` step 1 (subagent
auditor), at commit `94becb5`. **Not a failure yet** — recorded because the
mode it enables is a wedge, not a wrong number.

**Symptom (potential):** `tests/mesh/test_birdcage_conductor_sizing.py`
lines ≈149–158 break out of the sizing ladder when
`remaining < 2 × previous_rung_time`, using per-rank `time.perf_counter()`
with **no reduction**. If ranks ever straddle the threshold, one rank
`break`s while the other enters a collective mesh call — deadlock until the
container-side `timeout -k 30` kills the job. It cannot corrupt a reported
number (all rungs completed in the closing run, 41 s against a 300 s
budget), but a slower box or a bigger ladder could fire it.

**Cause:** wall-clock is rank-local state used in a collective control-flow
decision — the same class as the `cell_tags.values` rule in CLAUDE.md.

**Resolves with:** any commit that next touches this test — reduce the
decision (`comm.allreduce(remaining, MPI.MIN)` or decide on rank 0 and
`bcast`). `EX-21` must not copy the pattern into the example.

---

## Batch moved 2026-10-03

### ✅ RETIRED 2026-10-03 06:00 slot (§9 item 197 green re-run on the rule (f) band) — was OPEN 2026-10-02 (`TH-20` step 0a, §9 item 176 W1, landed by §9 item 182) — `test_th20_q1_to_q6_quarter_identities` is **red only under `FEM_EM_TH20_STEP0=1`**: Q3 reads **\|4·A₄₀₁,q/A₄₀₁,f − 1\| = 1.498e-4 > 5e-5**; Q0 and Q7 are green, and the default suites skip the test

| | |
| --- | --- |
| **Symptom** | `tests/validation/test_th20_quarter_f_human_mesh.py::test_th20_q1_to_q6_quarter_identities` (`FEM_EM_TH20_STEP0=1`, `-n 1`, real build): `(Q3) facets 401 19836, 501 1825, 502 1811; 401 on planes 0; untagged exterior 0; A401 q 5.195046113e-02 f 2.077707230e-01 \|4q/f-1\| 1.498e-04 (band 5e-05, predicted <= 3e-5)` (`20261002T185456Z_TH-20-step0a-W1.log:8160`), `failures: ['Q3 cavity \|4q/f-1\| 1.498e-04']` (`:8194`). Status 1, 204 s. Same window, green: Q0 (full hole 369 387 cells), Q1, Q2, Q4–Q6; the six Q7 refusals pass in the default suite. |
| **Verified at** | `dea66f77` plus the uncommitted step 0a module (landed by §9 item 182 from `attempt/TH-20-step0a-20261002T191039Z`, `src/fem_em_solver/io/mesh.py` and the test module only), by execution: one window. |
| **Cause** | **Not diagnosed.** The 5e-5 band (`F_HUMAN_WALL_AREA_RTOL`) was registered on a different comparison — hole vs solid on the same full sizing (`test_wf12_f_human_hole.py:83–106`). Candidate (untested at entry): triangulation — the quarter grades over 64 surfaces against the full mesh's 240 and carries 4 × 19 836 = 79 344 wall facets against 78 066 (+1.64 %); the facet model \|δ_f\|·(4N_q/N_f − 1) = 1.33e-4 is 89 % of the reading. §9 item 182 Part 2 (`tests/mesh/probe_th20_quarter_wall_area.py`) measures it. |
| **Scope** | `TH-20` step 0a only; the test is skipped unless the flag is set. No `main` default suite is affected. Q3, its 5e-5 band and `F_HUMAN_WALL_AREA_RTOL` are unedited. |
| **Resolves with** | §9 item 182's probe → the 03:00 review's rule (f) ruling → the 0a re-run (attempt 3). Retire with the 0a re-run's green commit. |
| **Ruling** | **2026-10-03 03:00 daily review (rule (f); the 10-03 weekly ratifies or reverts), `docs/planning/chunks/TH-20.md` "Ruling, 2026-10-03 03:00":** item 182's ladder read TRIANGULATION; Q3 re-registered under a new constant `TH20_QUARTER_WALL_AREA_RTOL` = 1.69e-3 (1.9 × the worst rung) plus the exact ring-port sub-identity Q3-b; `F_HUMAN_WALL_AREA_RTOL` unedited. Item 197 re-ran 0a green: `20261003T111707Z_TH-20-step0a-r3-W1.log` (`-n 1`, Status 0) and `20261003T112030Z_TH-20-step0a-r3-W2.log` (`-n 2`, controls asserted, Status 0), Q3 1.498e-4, Q3-b 1.1e-16 with 772 / 772 facets. |

### ✅ RETIRED 2026-10-02 13:30 slot (§9 item 175 green re-run, κ off per the operator ruling) — was OPEN 2026-10-02 (`WF-12` step 2′ W1, 12:00 implementer slot, §9 item 175) — the **κ-corrected** copper F-human 32×32 reads **σ_max 1.000002129 > `MATRIX_PASSIVITY_CEILING` 1.0**; the uncorrected copper 32×32 on the same fixture (H12b) reads 0.999993339

| | |
| --- | --- |
| **Symptom** | `tests/validation/test_wf12_f_human_copper_c_sweep.py::test_wf12_step2p_copper_c_sweep` (`FEM_EM_WF12_STEP2P=1`, `-n 8`, complex; module parked on `attempt/WF-12-step2p-20261002T172400Z`): `(M) 32x32 reciprocity 7.848739e-15 ...; sigma_max 1.000002129408 (ceiling 1.0, ASSERTED), 1-sigma_max -2.129e-06` (`20261002T171702Z_WF-12-step2p-W1.log:6359`), `AssertionError: 1.0000021294084984` (`:6367`, `:6432`). Status 1, 304 s. Same window, green: (A0) 369 387 cells rel 0; (iv) none; (R) P_src rel 2.111e-11, P_surf rel 1.505e-11 (band 1e-6); (P) 2.412e-13 (band 1e-6); NEG P_surf-blind 123.6× (bar 10×); κ_Cu 4.534948e-03 (σ = 800: 4.590991e-03); 32×32 1 factorising + 31 held; reciprocity 7.85e-15. |
| **Verified at** | `409f93bc` plus the uncommitted step 2′ test (now on the attempt branch), by execution: one window. |
| **Cause** | **Not diagnosed.** Measured facts only: the same copper fixture's uncorrected 32×32 is passive with margin 1 − σ_max = +6.66e-6 (`20261002T094751Z_WF-12-step2-W2.log`); applying the pooled κ = 4.53e-3 width correction moves σ_max by ≈ +8.8e-6, which crosses 1 only because copper leaves so little margin. On σ = 800 the κ-corrected 32×32 read σ_max 0.999923038 (`20261001T110957Z_WF-10-step1c-W2.log`), margin 7.7e-5, so the same shift would not have shown. Candidate (untested): the κ-corrected terminal current is not an energy-consistent port quantity, so the corrected S carries an O(κ · coupling) non-passive excess — the same correction `PORT-21` holds as a hypothesis. ‖S_Cu − S_800‖_F/‖S_800‖_F = 8.69e-02 (κ-corrected, H8a order, printed). |
| **Scope** | `WF-12` step 2′ only. No `main` test is affected; H12b (uncorrected) is green. The per-C identity residuals were not measured (the assert fires before the sweep). |
| **Ruling** | **Operator 2026-10-02 (interactive):** κ is not applied on exact-rectangle sheets (F-human's ring sheets are exact OCC rectangles); the passivity ceiling stands. The measurement the entry asked for is already on record: uncorrected σ_max 0.999993339 (passive) vs κ-corrected 1.000002129 on the same fixture. Retire this entry with item 175's green re-run (κ off). |
| **Resolves with** | A review or operator ruling, with a measurement: e.g. σ_max of the uncorrected vs κ-corrected 32×32 on both conductors from the same window, or the κ-corrected matrix's power balance per drive (P_src vs ΣP_sheet + P_vol + P_surf on corrected sheets). The passivity ceiling is not to be moved in a slot. |
| **Retired 2026-10-02 13:30 slot (§9 item 175)** | Re-run with κ off (`KAPPA_APPLIED = None`, κ printed as a diagnostic only), passivity ceiling unchanged: `(M) 32x32 reciprocity 6.912287e-15 ...; sigma_max 0.999993339099 (ceiling 1.0, ASSERTED), 1-sigma_max 6.661e-06` (`20261002T183158Z_WF-12-step2p-W1.log:6359`, Status 0, 411 s) — equal to H12b's uncorrected value; W2 green as well (`20261002T183901Z_WF-12-step2p-W2.log:6356`, Status 0). The κ-corrected non-passive excess stays a record (the 12:00 log), not a defect on `main`. |

### ✅ RETIRED 2026-10-01 21:00 slot (§9 item 158, on the 18:00 review's rule (f) ruling) — was OPEN 2026-10-01 (`WF-12` step 1 W1, 06:00 implementer slot, §9 item 129) — on the F-human rung as a hole mesh, the meshed cavity wall misses the solid route's conductor interface by **1.155e-05 relative against `INTERFACE_AREA_RTOL` 1e-5**; every other identity on the same window is green

| | |
| --- | --- |
| **Symptom** | `tests/validation/test_wf12_f_human_hole.py::test_wf12_step1_hole_mesh_identities` (`FEM_EM_WF12_STEP1=1`, `-n 1`, real build): `(iii) cavity 2.077707230e-01 vs solid conductor interface 2.077683226e-01 ... rel 1.155e-05 (band 1e-05)` (`20261001T112035Z_WF-12-step1-W1.log:16561`), `AssertionError ... assert 1.1553101716854641e-05 < 1e-05` (`:16571`). Same window, green: the hole meshes (369 387 cells, build 110.8 s); the solid control 507 266 cells (ratio 1.005200 to 504 642, band 0.01) with tag 1 present, hole < solid; (i) the 32 longitudinal ring sheets = chord · w to 2.2e-15 (band 1e-9); (ii) hole CAD + coil = solid CAD to 0.0 (band 1e-9), hole/(box − coil) − 1 = −3.984e-08, solid/box − 1 = −3.960e-08 (band 1e-7); (iv) census has no tag 1, has air, phantom and all 64 ring-port halves, sums to the global count (`:16559–16562`). |
| **Verified at** | `a53e1002` plus the uncommitted `WF-12` step 1 test, by execution: one window, 246 s, Status 1. |
| **Cause** | **Not diagnosed.** The band was set on the two-torus (`TH-15` step 2a read 3.7e-7, band ≈ 30× that) and gated on the 4-leg F-small birdcage; at F-human scale (82 tags touching the solid conductor, a 0.208 m² wall of tori and cylinders) the hole and solid routes are two independent triangulations of the same curved CAD surface, so a 1e-5-level area difference between inscribed meshes is a candidate explanation — not measured. The cavity is *larger* than the interface by 2.4e-06 m². The comparand includes the conductor's faces against the 16 leg-port boxes (tags 101–116), which F-small did not have. **Measured 2026-10-01 16:30 slot (§9 item 147, `mesh-probe`, `tests/mesh/probe_wf12_cavity_wall_area.py`) — pre-registered reading "mixed", not classified:** at `conductor_resolution` 1.6e-3 / 2.4e-3 / 3.2e-3 m (control reproduces W1 to 10 digits, \|Σ_class − A\| ≤ 2.1e-15) D = A_hole − A_solid = **+2.40e-06 / +4.84e-06 / −5.51e-06 m²**, D/(A_CAD − A_hole) +1.41e-3 / +1.83e-3 / −1.35e-3 — the sign flip fails both the triangulation and the face-mismatch tests. 100 % of D sits in the **air-facing wall (tag 2)**; no wall facet borders the leg-port boxes 101–116 on either route, and the ring-port wall is identical on both routes in area and facet count at every rung. D's sign follows the air-class facet-count difference (+32 / +68 / −50); \|D\| grows with the CAD shortfall (A/A_CAD − 1 = −8.14e-3 / −1.26e-2 / −1.95e-2, ≈ h¹ not the predicted h²). Untested interpretation: two independent triangulations of the same air-facing surface, error of either sign, not a missing or extra face. `20261001T214637Z_WF-12-probe-wall-h1p6-hole.log` … `20261001T215509Z_WF-12-probe-wall-h3p2-solid.log`. |
| **Scope** | `WF-12` step 1 only (anchor (iii)). W2 (the copper drive and its power identity) was not run; nothing on `main` is affected. |
| **Resolves with** | The item's negative-result rule ((iii) fails ⇒ generator finding, BLOCKED; next step a `mesh-probe` meshability probe). A band change, if any, is the review's or the operator's with a measurement — e.g. both routes' wall area against each other at two conductor resolutions, and the hole wall against the CAD surface area. **Ruled 2026-10-01 18:00 review (§9 rule (f); the 10-03 weekly ratifies or reverts; full text `docs/planning/chunks/WF-12.md` after item 147):** the sign flip excludes a fixed face mismatch, so the reading is triangulation scatter between two independent meshings of the air-facing wall. (iii) is re-registered on this fixture as `F_HUMAN_WALL_AREA_RTOL` = 5e-5 (measured \|D\|/A 1.155e-5 / 2.341e-5 / 2.684e-5), plus an exact sub-identity, ring-port wall equal ≤ 1e-8 (measured 5.4e-10). `INTERFACE_AREA_RTOL` 1e-5 is unchanged for the two-torus and F-small. This entry retires with §9 item 158's commit if (iii), (iii-b) and (P) are green there. |
| **Retired 2026-10-01 21:00 slot (§9 item 158)** | All three green on the committed test: (iii) 1.155e-05 ≤ `F_HUMAN_WALL_AREA_RTOL` 5e-5 and (iii-b) ring-port wall 5.148e-10 ≤ 1e-8 with 3076 / 3076 facets (`20261002T020254Z_WF-12-step1-W1b.log:16561–16562`, `-n 1`, Status 0, 245 s); (P) 2.412e-13 ≤ 1e-6 on the one Cu drive (`20261002T020713Z_WF-12-step1-W2.log`, `-n 8`, Status 0, 199 s). The band is the fixture's, pending the 10-03 weekly's ratification; if it reverts, this entry reopens. |

### ✅ RETIRED 2026-10-03 (22:30 CDT slot, §9 item 186; `GEO-1`) — was: 🔴 OPEN 2026-08-27 (`OPS-26` step 2 leg (a), second slot) — `test_birdcage_volumes_partition_the_box` aborts in gmsh with **the same "Invalid boundary mesh (overlapping facets)"** — this is the **third** geometry to carry that string, and the first on `birdcage_port_domain`'s own production path

> **RETIRED 2026-10-03 (§9 item 186), under the 08-26 ruling (b) as extended by
> the 2026-10-02 18:00 review (the weekly ratifies).** The test now builds on the
> graded control `GRADED_CONTROL_RESOLUTION = 4.8e-3`: 33 185 cells, V_mesh / V_box
> and Σ tags / V_mesh 1.000000000000, four port boxes 8.000000e-07, conductor 0.8104,
> phantom 0.9705, identical at `-n 2` / `-n 1` (`20261003T035654Z_GEO-1-W1.log`,
> `20261003T035842Z_GEO-1-W2.log`, Status 0). The ungraded build is now asserted
> to raise on every rank (`test_ungraded_build_still_hits_the_documented_floor`;
> rank 0 "overlapping facets", now on surfaces 59 / 79). The limitation itself
> stays OPEN in the 2026-08-25 entry.

> **MEASURED 2026-08-28 (`GEO-23` step 1, 09:00 slot) — geometry-deterministic,
> and this site is the family's positive control.** Red at **both** widths with
> the same string: `-n 1` **Status 1, 4 s** (`1 failed, 2 passed in 2.73s`,
> `20260828T140313Z_GEO-23-step1b-birdcagepart-n1.log`) and `-n 2` **Status 1,
> 5 s** (`20260828T140326Z_…-n2-ports.log`). Note the `-n 2` status: unlike the
> three sibling sites, this one **does not deadlock** — `birdcage_port_domain`
> re-raises the rank-0 gmsh throw as `RuntimeError: birdcage_port_domain
> geometry generation failed on rank 0` on *every* rank, so the command footers
> in 5 s where the unwrapped siblings burn 120 s each. That makes this entry
> the evidence that the family's deadlock is a **raise-path** property, not a
> geometry one, and wrapping the throw is a `GEO-23` step-2 lever that touches
> no mesh, band or record. The two adjacent tests in the module stayed green at
> both widths. Not laddered here — `GEO-21` step 2 already laddered this
> generator and `GEO-23` must not re-record it. Entry stays OPEN.
>
> **OWNER ASSIGNED 2026-08-27, 03:00 review: `GEO-23`** step 1 (b)/(c) —
> taken together with the two sibling entries, as this entry asks.
>
> **Where this fired.**
> `tests/mesh/test_birdcage_port_tags.py::test_birdcage_volumes_partition_the_box`,
> **real** build, `mpiexec -n 2`. Found by the `OPS-26` step 2 execution
> census (2026-08-26 21:00 implementer slot), first in the batch run
> (`20260827T022114Z_OPS-26-step2a-real4-mesh.log`, FAILED at 28% on both
> ranks) and then isolated for its traceback
> (`20260827T022935Z_OPS-26-step2a-mesh-red-tb.log`, **Status 1, 4 s**,
> `1 failed in 2.54s`).
>
> **Symptom, verbatim.**
>
> ```
> Exception: Invalid boundary mesh (overlapping facets) on surface 59 surface 79
> ```
>
> raised from `/usr/local/lib/gmsh.py:2189` inside
> `MeshGenerator.birdcage_port_domain`
> (`src/fem_em_solver/io/mesh.py:3245` → re-raised at `:3275`, wrapped at
> `:3276` as `RuntimeError: birdcage_port_domain geometry generation failed
> on rank 0`). The fragment line printed immediately before the abort is
>
> ```
> [birdcage-mesh] fragment volumes=26 conductor=1.030097e-04(20p) air=1.006440e-02(1p)
> phantom=2.261947e-04(1p) port_P1=8.000000e-07(1p) port_P2=8.000000e-07(1p)
> port_P3=8.000000e-07(1p) port_P4=8.000000e-07(1p)
> ```
>
> so the OCC fragment succeeds and the failure is in 2-D meshing, exactly as
> in the two entries below.
>
> **Why this matters — the string now spans three generators.** `GEO-21`
> filed it on the **open birdcage** conductor sizing, the 19:30 slot filed it
> on the **coil+phantom** generator, and this is `birdcage_port_domain` with
> ports and a phantom. Three different call paths, one symptom. The
> single-generator readings ("a coarse-resolution floor", "a fixture-specific
> sizing") no longer cover the observations; a shared 0.11 gmsh cause is now
> the more economical hypothesis. **Stated as a hypothesis from three shared
> error strings, not as a measurement** — nothing here bisects a resolution
> or attributes a cause.
>
> **Not a kill artifact.** The batch run that found it followed
> `20260827T022014Z_OPS-26-step2a-real3-cheap.log`, which exited **Status 0**,
> so the "do not trust a failure that follows a killed run" rule below does
> not apply; the isolated re-run then reproduced it in 2.54 s with an
> assertion-free gmsh exception, not a `dolfinx/jit.py` `RuntimeError`.
>
> **Adjacent, and green:** the other two tests in the same module pass
> (`test_birdcage_port_layout_diagnostics_match_the_closed_forms`,
> `test_birdcage_port_layout_rejects_too_small_or_overlapping_port_regions`),
> as do `tests/mesh/test_coil_phantom_conforming.py` (2/2) and
> `tests/materials/test_phantom_material_model.py` (3 green, 1 skipped) —
> two of the four coil+phantom consumers the 19:30 slot named as candidates
> to red the same way did **not** red. That narrows the blast radius but does
> not diagnose it.
>
> **Filed, not fixed, not re-recorded**, per the census item's own rule.
> Disposition belongs to a `mesh`-owning chunk, which should take all three
> entries together rather than one at a time.

## RETIRED 2026-09-20 (`OPS-59` (b), 07:30 slot) — `tests/unit/test_setup_figure_title.py::` the call-site count pin (`EXPECTED_EXAMPLE_CALL_SITES = 18`) was stale on `main`

**Retired** by replacing the literal with the identity it stood for. The red
was executed first, as the entry asked: `assert 21 == 18`
(`20260920T130047Z_OPS-59-prechange.log:763–764`) — **21**, not the 20 this
entry predicted, because `75082a0` (`ANS-6`, the same slot's item 14) added a
21st call site and a 21st figure two hours earlier. The test now computes
three counts and asserts them equal —
`call_sites=21 committed_setup_pngs=21 census_ok=21`
(`20260920T130232Z_OPS-59-postchange.log:576`), green at `-n 2`, 8.2 s — so an
`EX-57` slot no longer touches it. The title-length scan is unchanged.
Original entry follows.

**Found** by the 2026-09-20 03:00 review's `OPS-53` audit (the closure itself
PASSes: at `da9db81` the tree had exactly 18 `write_setup_figure(` call sites
under `examples/`). The same day's later slots added `mesh:6` (`6f75424`) and
`mesh:7` (`32f28d8`), so `git grep -l 'write_setup_figure(' -- 'examples/**/*.py'`
on `main` returns **20** and the equality at `:157` cannot hold. The review
tried to execute the module and could not — the docker socket is denied in
the review's sandbox (`20260920T080639Z_OPS-53-callsite-pin.log`, Status 1,
0 s, no test ran) — so this entry is a reading of the source, and the first
window of the fix records the red (or refutes this entry). **Cause:** the
2026-09-19 review's own item text asked for "18 call sites expected"; an
equality on a count the standing figure task increments is the `OPS-44`
artifact pin's disease (entry below) a second time. The title-length scan
the test exists for is unaffected. **Fix:** `OPS-59` (b) (§7; §9 item 16) —
replace the literal by the identity it stood for (call sites == committed
`*_setup.png` == census `ok`). Until then every `EX-57` slot makes it one
worse; do not bump the literal in passing.

## ✅ RETIRED 2026-09-21 (`OPS-59` (a), 06:00 implementer slot) — ~~2026-09-19 — `tests/unit/test_doc_reference_exit_codes.py::test_the_in_tree_exemption_cannot_silently_widen` is red on `main`: the `OPS-44` artifact pin was never widened for the `EX-57` setup figures~~

**Retired:** `ans6_…/metrics.json` declared in `LISTED_EXAMPLE_ARTIFACTS` on
the 2026-09-21 03:00 review's ruling; both modules green at `-n 2` (27 passed,
8.2 s), `pinned=28 (listed=6 by_rule=22) checker=28 git_ls_files=28`, extras
none (`20260921T111239Z_OPS-59.log:354,383`).

**Seen** in `20260919T234019Z_OPS-ans-folder-rename.log` (1 failed / 18
passed, `-n 2`, 7 s) while verifying the benchmark-folder rename; **not
caused by it** — every "extra item in the left set" is an
`examples/*/figures/*_setup.png`, and the test file was last edited
2026-09-09 (`1f78649`) while the first setup PNG landed 2026-09-13
(`40ec08b`). The test's own docstring says a chunk that commits a new
artifact must add it to `COMMITTED_EXAMPLE_ARTIFACTS`; the `EX-57` items
never did, so the pin has drifted by one entry per figure since. **Fix
(chunk opened 2026-09-20 03:00 review: `OPS-59` (a), §9 item 16):**
a small item that re-pins the set from `git ls-files` (the test prints the
list) and, better, makes the `EX-57` per-item template say so — the figure
task will keep committing PNGs for ≈30 more slots. Not a physics gate; the
docrefs checker itself passes (`dead=0`).

**NARROWED, STILL OPEN — 2026-09-20 (`OPS-59` (a), 07:30 slot).** The setup
figures are no longer the cause: the pin is now `LISTED_EXAMPLE_ARTIFACTS`
(5, declared by hand) **plus** `_setup_figures_admitted_by_rule()` (21,
derived from `check_example_setup_figures.py`'s own script↔figure pairing),
so `pinned=26 (listed=5 by_rule=21)` against `checker=27 git_ls_files=27`
(`20260920T130232Z_OPS-59-postchange.log:358`) and the extras fell **22 → 1**
(`…prechange.log:736–758` → `…postchange.log:386`). The one survivor is
**`examples/ansys_benchmarks/ans6_copper_birdcage_four_port_10_64_128MHz/metrics.json`**,
a **non-setup** artifact committed undeclared by `75082a0` (`ANS-6` runnable
half, the same slot's item 14). `OPS-59` pre-registered exactly that case as a
negative result — "name it, leave it red, stop — that one is a real unpinned
artifact" — so the implementer did not declare it and the test is still red
(1 failed / 26 passed, `-n 2`, 8.2 s). The rule half is controlled: a
non-setup PNG and an orphan `*_setup.png` are both refused by the rule and
both turn the identity red if listed (`…postchange.log:420`). **Resolves
with** a review declaring that one path in `LISTED_EXAMPLE_ARTIFACTS` under
the `ANS-1` "an `ans:` case commits its own `metrics.json`" rule — `OPS-44`'s
late declaration of the `ANS-2`/`ANS-4` files is the precedent — plus one
≈ 8 s window, which closes `OPS-59` and retires this entry.

## ✅ RETIRED 2026-09-20 (`PORT-20` steps 1–3, 04:30 + 09:00 implementer slots) — ~~🟡 OPEN 2026-09-19 (`ANS-3` adjudication, weekly review) — on the gap-voltage / impressed-current route `run_n_port_sparameter_sweep` reports an `S` that is **not the 50 Ω S-matrix**; its `z_matrix` is right~~

| | |
| --- | --- |
| **Symptom** | On the two-torus fixture (`EX-20`, `ans:3`) the tabulated `S₂₁` and the `S₂₁` implied by the *same run's* `Z` through `z_to_s(Z, 50 Ω)` differ by **≈ 0.016 in magnitude** (tabulated \|S₂₁\| = 0.0216; Z-implied ≈ 0.0377 by hand from the tabulated `Z`, to be printed exactly by `PORT-20` — both our numbers, `examples/ansys_benchmarks/ans3_two_torus_gap_ports_10MHz/COMPARISON.md`). No test is red: reciprocity and passivity pass on either matrix. |
| **Found by** | The `ANS-3` AED comparison, adjudicated 2026-09-19 (PROJECT_PLAN §10): our `Z` agrees with the independent code on the primary row; our tabulated `S` does not, our Z-implied `S` does. Qualitative verdict only — numbers in the gitignored `docs/private/ans3-ans4-adjudication-2026-09-19.md`. |
| **Cause** | Read, not yet measured: on this route each drive impresses a current across one gap and leaves the other port **open**, so the undriven port's incident wave `a_i ≠ 0` and `S_ij = b_i/a_j` (`ports/sparameters.py::_assemble_sparameter_matrix`) is not a column of the S-matrix. The assembly is exact only when every undriven port is terminated in the reference — the lumped-sheet route at `Z_p = z0`, which is what `PORT-9`/`PORT-11`/`ANS-4` use and what the docstring says. On the current route the "terminated transimpedance" *is* the open-circuit `Z`, and `z_to_s(Z)` is the S-matrix. |
| **Not affected** | Every lumped-sheet-route 4×4 / 32×32 (`PORT-9`, `PORT-11`, `PORT-13`, `ANS-4`, `TH-14`, `TH-15`), `PORT-1`'s gate of record (the mutual ratio is read from `Im Z₂₁`), `MAT-6`, `ANS-1`, `ANS-2`. |
| **Affected** | Any S entry, `‖S‖₂`, or Touchstone export produced on the gap-voltage route: `EX-20` (`ports:2`), `ans:3`'s S table and its `‖S‖₂ = 0.864809` record, `PORT-1` step 4's passivity figure. They are self-consistent records of a quantity that is not S. |
| **Fix** | `PORT-20` (§7). Do not "fix in passing": the reproduction records above move when it lands and are re-recorded by that chunk only. |
| **Retirement** | `PORT-20` steps 1–2 landed 2026-09-20 (the route converts its open-circuit `Z`: `S = z_to_s(Z)`, agreeing with `ports/circuit.py`'s independent implementation at 9.4e-17 relative, closed form at 0.000e+00 on synthetic 2-/3-ports); **step 3, this slot** (`20260920T140828Z_PORT-20.log`, `ports:2`, 181 s, Status 0; `20260920T141141Z_PORT-20.log`, the `PORT-1` step-4 gate module, 18 passed / 190 s) re-recorded the four S-derived constants and removed step 2's strict `xfail`. Measured on the corrected route: `‖S − Sᵀ‖/‖S‖ = 3.112130445718013e-05`, `‖S‖₂ = 0.8613568944845725`; the superseded power-wave digits miss that same run by **1.647e-05 absolute** in the symmetry residual and **3.992e-03 relative** in `‖S‖₂` (gate module: 1.646e-05 and 3.453e-03) — the predicted ≈ 0.016-in-`\|S₂₁\|` size, asserted as each record's negative control. Nothing about `Z` moved: the mutual-ratio records reproduce at 2.10e-10 / 2.06e-10 relative against an unchanged 1e-6 band. **One residue, entry below:** `ans:3`'s *generated* `metrics.json` / `COMPARISON.md` still carry the pre-correction S table. |

## ✅ RETIRED 2026-09-27 (03:00 daily review; fixed by `4719e82`, 2026-09-26 09:00 slot) — was: 2026-09-20 — 🟡 OPEN (`PORT-20` step 3, 09:00 implementer slot) — `ans:3`'s tracked `metrics.json` / `COMPARISON.md` still carry the **pre-`PORT-20`** S table; the example cannot be re-run in a slot without rewriting the gitignored private comparison

| | |
| --- | --- |
| **Symptom** | `examples/ansys_benchmarks/ans3_two_torus_gap_ports_10MHz/metrics.json` and `COMPARISON.md` were generated on 2026-09-19 by the pre-`PORT-20` power-wave assembly, so their `S` block (and the `‖S‖₂` / symmetry entries echoed there) is a record of the quantity `PORT-20` corrected. No test is red: the case's *records* are imported from `examples/ports/02_package_sparameter_sweep.py`, which step 3 re-recorded, so the next run of `ans:3` will gate on the corrected digits — the stale numbers are in the generated artifacts only. Our columns only; the AED columns of the tracked file are blank by construction. |
| **Found by** | `PORT-20` step 3 (§9 item 17), which scoped `ans:3`'s our-half columns and then hit its own pre-registered trap. |
| **Cause** | `03_two_torus_gap_ports_10MHz.py:725-727` writes the gitignored `COMPARISON_private.md` (AED columns filled) **unconditionally whenever `aed_results/` is present**, and it is present on this box. The step's instructions forbid a slot from opening or rewriting that file, so the run that would refresh the tracked artifacts cannot be issued in a scheduled slot. |
| **Not affected** | Everything `PORT-20` step 3 did re-record: `examples/ports/02_package_sparameter_sweep.py`'s four constants and `tests/validation/test_port_package_sparameters.py`'s two, both green this slot. The `ANS-3` adjudication itself (its verdict rests on `Z`, which did not move, and on the S disagreement this correction explains). |
| **Unblock condition** | Either (a) the operator (or an interactive session with the operator present) re-runs `ans:3` and lets the private file be regenerated from the unchanged `aed_results/`, or (b) a small `src`-free change gives the script an opt-out (`--no-private` / `FEM_EM_ANS_SKIP_PRIVATE=1`) so a slot can refresh the tracked half alone — then re-run `ans:3` and commit the regenerated `metrics.json` / `COMPARISON.md`. Until then, read the S rows of that case's tracked artifacts as superseded. |
| **Retired by** | `4719e82` (`EX-57` fallback batch, `ans:3` setup figure, 2026-09-26 09:00 slot) re-ran `ans:3` at `-n 2` and committed the regenerated tracked artifacts: `‖S‖₂` **0.861356894** and `‖S − Sᵀ‖/‖S‖` **3.112129e-05** (`COMPARISON.md`), which are `PORT-20` step 3's re-recorded 0.8613568944845725 / 3.112130445718013e-05 above; `Z₂₂` moved by 2e-7 relative only (`20260926T142337Z_EX-57.log`, Status 0). The executor's "matches to 9–10 figures" described `Z`, not `S`; the S move is this entry's intended correction. **How it was unblocked, not how it should have been:** neither (a) nor (b) — the slot ran the script with `aed_results/` present, so the gitignored `COMPARISON_private.md` was regenerated from the unchanged `aed_results/` (09-26 09:25 file time; no slot read it). The same happened to `ans:1` (07:30 slot) and `ans:2` (09:00 slot). Operator informed (dashboard); the `ans:4` fallback item carries opt-out (b) so it does not recur. |

---

## Moved 2026-10-05 (`POST-10` step 1, in its fixing commit)

### ✅ RETIRED 2026-10-05 (`POST-10` step 1, 15:00 CDT implementer slot, §9 item 254: `parent=(parent_mesh, parent_cells)` on the owner rule, green at -n 1 / 2 / 4) — was: 2026-10-04 — 🟡 OPEN (`ANS-9` step 1a, §9 item 218; filed by the 10:30 daily review) — `evaluate_vector_field_parallel(owner="lowest_original_cell")` raises `RuntimeError: Mixed topology unsupported.` on a `create_submesh` mesh

| | |
| --- | --- |
| **Symptom** | `examples/ansys_benchmarks/ans9_coil_driven_gel_heating_10MHz/09_coil_driven_gel_heating_10MHz.py`'s first W1 crashed at the owner-rule call on the phantom submesh's temperature field (`20261004T142026Z_ANS-9-step1a-W1.log:1949–1978`). The script switched to the default owner (θ is continuous P1, so the owner choice cannot change a value) and the run of record is green (`…142653Z_ANS-9-step1a-W1.log`). |
| **Cause** | Not diagnosed. `_original_cell_index` (`src/fem_em_solver/post/evaluation.py:19–40`) reads `mesh.topology.original_cell_index`, which the DolfinX build refuses on a submesh's topology. |
| **Consequence** | The owner rule — CLAUDE.md's rank-safe point-evaluation rule for fields discontinuous across cells — is unusable on any `create_submesh` mesh. Continuous fields are unaffected. **Watch:** ~~`MAT-19` … on `MAT-12`'s full-mesh sphere it is unaffected~~ — corrected 2026-10-05 (`MAT-19` step 1, §9 item 243): `MAT-12`'s thermal field is on a `create_submesh` mesh too (`_sphere_submesh`). `MAT-19` routes around the rule by integration (`thermal/probe.py`: an owned-cell ball average, one allreduce, no point owner), so it does not exercise or fix this entry. Any other point evaluation on a submesh is still affected. |
| **Resolves with** | A `post/evaluation.py` fix mapping submesh cells to parent original indices (**opened 2026-10-05 10:30 review as `POST-10`, §9 item 254**), or a documented restriction. Remove this entry in the commit that fixes it. |

---

## Batch moved 2026-10-05

### ✅ RETIRED 2026-10-04 as moot (03:00 daily review, on the 2026-10-03 weekly's "retire it as moot or fix the test, your call") — was 🟡 OPEN 2026-10-03 (`xxl` window, `PORT-21` step 3b, `xl-pending.md` entry 23; read by the 10:30 daily review) — `test_the_step3_knobs_reach_the_solve` is **red only under the step-3 knobs at `FEM_EM_PORT21_RESOLUTION=0.005`, degree 2**: the degree-1 reference JSON it compares against was never written, so the degree-knob assertion never runs

| | |
| --- | --- |
| **Symptom** | `tests/validation/test_port21_feed_sensitivity.py::test_the_step3_knobs_reach_the_solve` (`-n 16`, complex, `FEM_EM_PORT21_DEGREE=2`, `FEM_EM_PORT21_RESOLUTION=0.005`, frequencies 10 / 64 MHz, variants `control_prime a`): `Failed: no degree-1 reading at logs/port21/f_10.000MHz_p1_h0.005.json to compare control' against: run W1 (knobs unset) at this frequency first` (`20261003T070008Z_PORT-21-step3b.log:6198`, `:6262`). 1 failed, 15 passed, 4 skipped on all 16 ranks; Status 1, 9361 s. Same window, green: the imported identity gates at both frequencies (reciprocity ≤ 1.47e-13, σ_max ≤ 0.999995519, spreads ≤ 0.032 %, `:6071–6090`); the resolution-knob inequality (592 550 ≠ 116 118 cells, `:6145`); the (a) "knob reached the solve" floor (`:6109–6110`). |
| **Verified at** | `341fe583`, by execution (the `xxl` window itself); `log-pathologist` ruling, 2026-10-03 10:30 review. |
| **Cause** | **Diagnosed:** a missing prerequisite, not a gate red. The test reads `_readings_path(f, 1, resolution)` (`test_port21_feed_sensitivity.py:902–907`), i.e. a degree-1 run *at the same resolution knob*. No run has ever solved degree 1 at h = 0.005 (only `p1_h0.0075` and `p1_hdefault` exist under the gitignored `logs/port21/`), and entry 23's command is degree 2 only — the red was predictable at commissioning. The failure message's remedy is also wrong: "knobs unset" writes `p1_hdefault`, which this lookup cannot use. |
| **Scope** | Only step-3-knob windows at a resolution with no degree-1 reading. Entry 22 (h = 0.0075, queued 10-05) has its `p1_h0.0075` references and is not affected. The default suites skip the test (`_knobs_set()` false). No band moved; the κ readings the window persisted (`:3195`, `:5993`) are unaffected. |
| **Resolves with** | A degree-1 window with `FEM_EM_PORT21_RESOLUTION=0.005` at 10 and 64 MHz (writes the references), followed by any degree-2 window at that resolution, plus the message text corrected to name the resolution knob. Whether that is worth a window after the operator's 2026-10-02 κ ruling (κ off on exact rectangles; this family is ragged-sheet, κ-corrected) is the weekly's (chain F, F7 / F9). **Retired as moot 2026-10-04 (03:00 daily review):** the 2026-10-03 weekly closed the ragged-sheet κ-corrected family — entry 22 superseded, no reference window commissioned, the licence "does not extend any ragged κ-corrected family" — so no window will ever run the step-3 knobs at h = 0.005 and the test cannot be reached by any commissioned command; the default suites skip it. The test and its wrong remedy text are unchanged; if the family is ever reopened, this entry returns with it (text in git, `85612f71`). |

## RETIRED 2026-10-06 (`ANS-7` step 7′, §9 item 262) — was 2026-10-04 — 🟡 OPEN (`ANS-7` step 7, §9 item 217, chain T21; filed by the 10:30 daily review) — `test_ans7_exact_strip.py::test_t21_t2_inmodel_zero` is red by pre-registration: the in-model zero of Im `Y_in` with the k = 1-block `C_mode` sits 4.82e-6 above 64.000 MHz, against `BANDS["f0_hz"]` (1e-6, i.e. 64 Hz)

**Retired 2026-10-06** (`ANS-7` step 7′, §9 item 262, chain T21′, 06:00 implementer slot) under the weekly ruling "T21 / T22 — \"tuned\" means the in-model zero": driven at C_in (item 217's root), the degree-2 in-model zero sits at f₀,A/64 MHz − 1 = −2.536e-6, inside the pre-registered grid resolution `T21P_ZERO_RESOLUTION` 1.287e-5 (`test_ans7_exact_strip.py::test_t21p_inmodel_zero_cin`, `20261006T110217Z_ANS-7-item262-W1-cin.log`). The C_mode readings (4.821e-6 / 1.390e-6) stay records; `BANDS` is unchanged and `test_t21_t2_inmodel_zero` skips under the `cin` knob.


| | |
| --- | --- |
| **Symptom** | Only under the opt-in `FEM_EM_ANS7_EXACT_TUNE=inmodel` (unset ⇒ the T21 tests skip; the default suites never reach it): `assert abs(r["f0_hz"] / TUNE_TARGET_HZ - 1.0) <= BANDS["f0_hz"]` fails with f₀,A = 64 000 308.55 Hz, \|f₀,A/64 MHz − 1\| = 4.821e-6 (`20261004T140738Z_ANS-7-item217-W2-inmodel.log`; 14 passed, 16 skipped, 1 failed). Verified against `84da27b7` (landed by path on `main` by this review). |
| **Pre-registered** | The brief carried `PREDICTED RED: (T2) range [1e-6, 5e-5]`; the reading is inside it and every other anchor is green — (C), (T1), (AB) ≤ 1.29e-5 at `REDUCTION_BAND`, (N1), (W) (`…140347Z_ANS-7-item217-W1-cmode.log`, `…140738Z…`). A record, not a defect report on the solver. |
| **Cause (not diagnosed past the slot's reading)** | `C_mode` is built from the k = 1 block of the 4×4 alone; the full terminated network's own tuning capacitor `C_in` differs by (C_mode − C_in)/C_in = −1.495e-5 (C_mode,Cu 63.589141 pF, C_in 63.590091 pF). The slot's decomposition: +7.35e-6 from the construction minus 2.53e-6 of Route A − Route B gap. The k-block off-diagonal reads 4.56e-6 against a predicted ≤ 1e-6. |
| **Measured at 127.74 MHz (T22, `ANS-7` step 8, §9 item 228; landed as a record by the 2026-10-04 18:00 daily review)** | The same assertion under `FEM_EM_ANS7_EXACT_TUNE_MHZ=127.74` is red inside its pre-registered `[1e-6, 5e-5]`: f₀,A = 127 740 177.57 Hz, \|f₀,A/127.74 MHz − 1\| = **1.390e-6** (`20261004T204206Z_ANS-7-item228-W4-inmodel.log`; 14 passed, 1 failed). Every other anchor green (`…202814Z…-W1-nodes.log`, `…203339Z…-W2-nodes.log`, `…203817Z…-W3-cmode.log`): (T1) Route B at C_mode +3.653e-6, (AB) ≤ 1.62e-5 at `REDUCTION_BAND`, (N1) ≤ 2.2e-12, (W) 1.07e-13. Split: construction c = +3.653e-6, route gap r = f₀,A/f₀,B − 1 = **−2.263e-6** (−2.53e-6 at 64 MHz) — **r is flat in f and alone exceeds the 1e-6 band**, so neither `C_mode` nor `C_in` (estimated δf/f −2.263e-6) puts the in-model zero inside `BANDS`. Code landed by path from `attempt/ANS-7-step8-20261004T2049Z` (`22de12f4`); the tests skip unless the TUNE knob is set. |
| **Resolves with** | The weekly's ruling on T21 / T22's target (whether "tuned" is read on `C_mode` with its construction residual, or on `C_in`, given that the Route A − Route B gap alone is ≈ 2.3–2.5e-6 at both frequencies). **The band does not move** (CLAUDE.md hard rules). |

---

## Batch moved 2026-10-07

### ✅ RETIRED 2026-10-07 (`OPS-68`) — was `MAG-18` anchor (ii) is unreachable as pre-registered: the magnetostatic solve's own cross-width floor is ≈1e-7, not 1e-10 (2026-08-22)

> **Retired 2026-10-07 by `OPS-68` (§9 item 307).** The entry's own condition ("Re-registering (ii) at the measured floor … is a review decision; `MAG-18` stays 🟡 until it is made") is met: the body's 2026-08-23 03:00 ruling re-registered (ii) at ≤ 1e-6, the logged `-n 2` / `-n 4` pair satisfies it (`20260823T003327Z_MAG-18-record-probe.log.gz`, `20260823T003406Z_MAG-18-record-n4.log.gz`, both `Status: 0`; landed with `d494d814`), and `MAG-18` is ✅ (§7 :1667). The ≈1e-7 cross-width floor figure stays recorded in the body below. Body kept verbatim below.

> **Not a failing test** — no assertion is red. This records a
> pre-registered done-when clause that measurement showed could not be met
> by any statistic, so the next reader does not spend a slot re-measuring
> it.
>
> **The clause.** `MAG-18` anchor (ii), PROJECT_PLAN §7 and §9 item 1: "the
> h = 0.0025 value at `-n 2` and `-n 4` agree to **1e-10 relative** — a
> reduced integral has no sampler; this is the control that the new
> statistic lacks the defect the old one had."
>
> **Measured** (`20260823T003327Z_MAG-18-record-probe.log`, `-n 2`, 31 s;
> `20260823T003406Z_MAG-18-record-n4.log`, `-n 4`, 26 s; same 145 884-cell
> mesh in both, `Status: 0`):
>
> | statistic | `-n 2` | `-n 4` | relative |
> |---|---|---|---|
> | `E_Ω` (new) | 1.0728835983e-01 | 1.0728836764e-01 | **7.28e-08** |
> | 10-point, n_points = 8 (retired) | 15.802788% | 15.802785% | 1.9e-07 |
>
> **Cause: the linear solve, not the statistic.** `MagnetostaticSolver`
> runs `ksp_type=preonly, pc_type=lu` — a *direct* factorization, whose
> pivot/elimination order follows the mesh partition, so the computed `A`
> itself differs at roundoff amplified by the gauge-penalty system's
> conditioning. The retired sampled statistic moves the same way on the
> same two runs, which is the discriminator: a defect of the *sampler*
> would not appear in an assembled integral, and this appears in both.
> ≈1e-7 is therefore the solve's cross-width reproducibility floor and no
> functional of this field can read below it.
>
> **What it does and does not mean.** It does **not** revive the sampler
> objection: anchor (ii) existed to show `E_Ω` has no sample-count
> dependence, and it has none — the 34% swing the old statistic showed
> under `n_points` has no counterpart here. It **does** mean the 1e-10
> number was written without knowing the solver's floor.
>
> **Nothing was loosened in-slot.** `E_OMEGA_RECORD_BAND` is the
> separately pre-registered **record** band 1e-4, not a relaxed 1e-10.
> Re-registering (ii) at the measured floor (1e-6 would be 100× the
> observation and still 1e5× tighter than the record band) is a review
> decision; `MAG-18` stays 🟡 until it is made. A cheap alternative the
> review may prefer: assert cross-width agreement *directly* by running
> the record test at both widths in one harness invocation rather than
> against a hard-coded constant.
>
> **✅ RESOLVED 2026-08-23 (03:00 review) — (ii) re-registered at ≤ 1e-6
> relative**, 14× the measured 7.28e-08 and five orders below the 34%
> sampler swing it was commissioned to exclude; the logged `-n 2`/`-n 4`
> pair satisfies it and `MAG-18` is ✅ (PROJECT_PLAN §7). The
> pre-registration error was in the clause, not the solver. Entry kept
> for the floor figure: ≈1e-7 cross-width on a direct-LU magnetostatic
> solve is now a *known* floor for any future rank-independence clause,
> and the 2026-08-23 run-to-run entry below puts the single-width floor
> at ≈1e-9–1e-10.
>
> **The floor confirmed on 0.11 (2026-08-23, 19:30 slot).** The re-gate
> re-measured the same pair on the image `main` now boots: `-n 2`
> 1.0617170177e-01 vs `-n 4` 1.0617175341e-01, **4.86e-07 relative** —
> 6.7× the 0.7.2 observation, still inside the re-registered 1e-6, and
> the two same-width runs agree to 1.86e-08. So the floor is a property
> of the direct-LU-under-partition solve and not of one image, and the
> 1e-6 clause has ≈2× headroom on 0.11, not 14×. A future clause tighter
> than 1e-6 on this fixture would be pre-registering against the solver
> again.

### ✅ RETIRED 2026-10-07 (`OPS-67`) — was 🟡 OPEN 2026-09-03 (`PORT-13` step 1, 22:30 implementer slot) — `GEO-20`'s ring-gap port sheets are **transverse** sections, so the lumped-sheet port model has no `h` on any ring port: the 32-ring-port layout cannot be terminated or driven as `PORT-9` terminates a leg port

> **Retired 2026-10-07 by `OPS-67` (§9 item 303).** The Resolves-with condition is met: `GEO-26` added `ring_sheet_orientation="longitudinal"` (step 1, `41cd5c84`), step 2 recorded the 16-leg rung (`1ad8ba3c`, `20260903T170351Z_GEO-26.log`, the terminal-class STOP that step 3 resolved), and step 3 registered the 16-leg rung's terminal-area band green (`2445b9ed`, `20260904T034441Z_GEO-26.log`, 4 passed, Status 0). `GEO-26` is ✅ and `PORT-13` is ✅ (steps 1–3, the 32×32 on the 32-ring-port layout reciprocal / passive / C16), so the transverse-sheet blocker no longer stands. Body kept verbatim below.

| | |
|---|---|
| **Symptom** | `PORT-13` step 1 (one single-port solve on `mesh:9`) could not be started. On the 4-leg ring-gapped rung (110 786 cells, `GEO-20` step 1's record at ratio 1.000000) every one of the 8 ring sheets spans **7.69e-18 … 1.43e-17 m** along its own `φ̂` — the drive direction — and exactly `1.000000000000 w` along `û` and `ẑ` (`w = 1.000000000e-02 m`), area `1.000000000000 w²` (`20260903T033437Z_PORT-13.log:6959–6966`, Status 0 / 29 s). No solve ran; the parked measurement is `tests/mesh/test_birdcage_ring_sheet_orientation.py` on `attempt/PORT-13-20260903T033437Z` (`30756cf`). |
| **Cause** | `io/mesh.py:3713–3735` builds the ring sheet's four corners at one `φ_c` (varying only `u` and `z`) and records its normal as `φ̂` — the gap's mid-section, "the w x w rectangle" (`:3057–3058`). The `GEO-18` leg pattern was carried over by the wrong analogy: the leg sheet contains the *drive direction* (`ẑ`) and the radial direction, normal `φ̂`; the ring's drive direction *is* `φ̂`, so its analogue must contain `φ̂` and lie in a plane whose normal is transverse. The port model (`ports/lumped.py:110,253,283–289`) needs a sheet spanning the gap along the drive: `R_s = Z_p·w/h`, `E_src = V/h`, `I = (1/R_s)∫E·ĥ dS / h`. On the transverse sheet `ĥ = φ̂` is the facet normal, so `E·ĥ` is the *normal* trace of an H(curl) field on an interior facet — not single-valued for Nédélec elements (`_restrict`'s `'+'` choice is justified only for the tangential trace, `lumped.py:195–206`) — and `w = A/h` is undefined at `h = 0`. **The model would not raise:** `gap_height_m` is caller-supplied (`lumped.py:148,322`), so a solve with an arbitrary `h` would run and produce numbers with no meaning (`log-pathologist` ruling, 2026-09-03 03:00 review). |
| **Consequence** | `PORT-13` is 🚫 (§7 row, §9). Every `GEO-20` identity is unaffected — the transverse section is a valid mesh object and its records (area, C4/C16 spreads, cell counts 110 786 / 265 621) stand; `EX-35` is unchanged. No leg-port result moves (`PORT-9` / `PORT-11` / `WF-6` use leg sheets, whose `h_bbox/dz − 1 < 1e-9` is gated in `test_birdcage_port_sheets.py:274`). |
| **Resolves with** | **`GEO-26`** (opened 2026-09-03 03:00 review, §9 item 1): a `ring_sheet_orientation="longitudinal"` mode emitting the planar rectangle in the plane `u = R` per ring gap — normal `û(φ_c)`, chord `2R·tan α` along `φ̂`, `w` along `ẑ` — with the default unchanged so nothing above moves; the parked module becomes its negative control. This entry retires with `GEO-26` step 2's 16-leg record (the control `PORT-13` step 1 re-opens on), which must point its retired-by row here. |

### ✅ RETIRED 2026-10-07 (`OPS-68`) — was `check_example_doc_references.py` exit codes: 0 clean, 1 real defect, 2 staleness only (`OPS-19`, 2026-08-16 — contract, not an issue)

> **Retired 2026-10-07 by `OPS-68` (§9 item 307).** Resolved 2026-08-16 by `OPS-19` step 1 (`c786f5fb`, `20260816T213312Z_OPS-19-step1-rerun.log.gz`: 8 passed, `Status: 0`); `OPS-19` is ✅ (§7 :1328). The contract itself (`exit != 1`, staleness as `exit 2`) remains stated in §9's standing rules (PROJECT_PLAN.md :3684), so retiring this entry loses no instruction. Body kept verbatim below.

**Resolved 2026-08-16 (`OPS-19` step 1).** Staleness no longer owns the exit
code, so a chunk touching examples can read its companion docrefs log by its
status alone:

| code | meaning | who caused it |
|---|---|---|
| 0 | every reference resolves, every artifact fresh, guide pass green | — |
| 1 | **hard** violation: dead reference, missing guide, missing heading | usually the chunk in the slot |
| 2 | staleness only: references resolve, guides green, some artifact older than `--max-age-s` | the backlog — nobody has re-run those examples |

The last line of the checker's output states the split for a caller that does
not want to parse the body: `RESULT: dead=<n> guide=<n> stale=<n>
stale_severity=<fail|report> exit=<code>`. `--stale-severity fail` restores
the pre-`OPS-19` all-or-nothing reading (staleness exits 1); the default is
`report`. `--max-age-s` (`OPS-15`'s 48 h) did **not** move.

**On `main` at the time of the split**: `dead=0 guide=0 stale=24
stale_severity=report exit=2` — the 24 stale `paraview_output/` files (aged
105–141 h, magnetostatics/MRI examples) are still there, and regenerating
them is still compute rather than a documentation fix. What changed is that
they no longer show up as a failure. Guide pass green at 21/21 examples,
0 pending. Verified in `20260816T213312Z_OPS-19-step1-rerun.log` (8 passed,
1.91 s); the contract is pinned by `tests/unit/test_doc_reference_exit_codes.py`,
whose negative control is a guide naming an artifact no run ever wrote — that
must still exit 1 after the split, or the checker was switched off rather than
sharpened.

*(History: the pre-split behaviour — exit 1 on every invocation, "gate on the
guide-pass violation count, not on exit 0" — cost `EX-20` and `ANS-3` a
red-but-benign companion log each on 2026-08-16, which is what commissioned
`OPS-19`. The `EX-18` heading violations once filed here, three missing
headings in `examples/ports/01_two_torus_port_pair.md`, were fixed 2026-08-16.)*

### ✅ RETIRED 2026-10-07 (`OPS-67`; resolved 2026-09-01 by `TH-13` step 3a‴) — `test_coil_loading_degree2.py` no longer returns inside its 570 s ceiling (exit 124 at 571 s, **twice** on 2026-08-31): the full module ran in 543 s on 2026-08-18, a 27 s margin something ate (`TH-13` steps 3a / 3a″; the original "mesh no longer builds" diagnosis corrected by the 10:30 review and then refuted by measurement at 17:00Z — **the cost is the degree-2 pair alone, ≥ 524 s of the 571 s**; the mesh is 4.3 s and the degree-1 phase is green in 46.8 s)

> **Retired 2026-10-07 by `OPS-67` (§9 item 303).** This is the original open text of the entry the archive already carries as "✅ RETIRED 2026-09-01" (`known-issues-archive.md` ≈ :3278); it is **not** byte-identical to that copy (`20261007T185802Z_OPS-67-diff.log`: 83 live lines against 40 archived, the archive holding the retirement summary), so the body is kept here verbatim. Resolved by `TH-13` step 3a‴ (`f494dc75`): the case was shrunk, not the ceiling raised, and all three windows returned footers — `20260901T033335Z_TH-13-step3a3-original-probe.log` (exit 0, 49 s), `20260901T033434Z_TH-13-step3a3-degree2-loaded.log` (374 s) and `20260901T034059Z_TH-13-step3a3-degree2-free.log` (405 s), each with exactly its pre-registered degree-2 identity red.

**Verified at the step-3a working tree, 04:38 implementer slot,
`20260831T093807Z_TH-13-step3a-coil-degree2-regression.log`** (exit **124**,
571 s, `-n 8`, complex build, `TH12_STEP2_MODE=full`).

| | |
|---|---|
| **Tests** | `tests/validation/test_coil_loading_degree2.py` — the whole module; it never got past `test_the_mesh_is_the_mat6_step3_baseline`. |
| **Symptom** | The run is killed by the container-side `timeout -k 30 570` with the log showing collection, the instant `SKIPPED` of `test_the_memory_exponent_measured_on_the_fine_rung`, and then `test_the_mesh_is_the_mat6_step3_baseline` still in progress on all eight ranks. No assertion was reached and no number moved; `-k 30` worked — the container was healthy afterwards, zero stray `python3`, `memory.max` 137438953472. **Corrected 2026-08-31 10:30 review:** pytest ran **without `-s`**, so the log *cannot* show progress inside the test — and that test consumes the module-scoped `degree_rows` fixture, which is the mesh **plus** the degree-1 control solve **plus** the 61.94 GiB degree-2 factorization. "Still inside the mesh build" was never established by this log. |
| **Cause** | **Not fully diagnosed, but the mesh is excluded (10:30 review, by reading).** The identical mesh — same generator, same `near 0.005` parameters, same 138 490 cells — built in **4.5 s** at 11:02Z the *same morning* (`20260831T110240Z_ANS-5-step1b-ans1.log`) and in **4.8 s** on 2026-08-27 (`20260827T183121Z_OPS-26-step2f-richardson.log`), both post-`OPS-18`, so the gmsh-4.15.2 mesh-time hypothesis is refuted for this mesh and ≈565 s of the killed window belong to something after or alongside it. What remains: on 2026-08-18 (`20260818T200059Z_TH-12-step2-full.log`) the *entire* module completed in **543 s** at the same `-n 8` on the same ceiling — a **27 s (5%) margin** that ordinary machine load or factorization scatter can eat, and the killed run's captured stdout cannot apportion. |
| **Not** | **Not** to be disposed of by raising the ceiling: CLAUDE.md §5.1 says overrun ⇒ kill and shrink, and this module has no cheap variant — `TH12_STEP2_MODE=probe` skips the degree-2 solve but still builds this same mesh first, so it does not fit either. Not evidence about the degree-2 identity bound, which was not reached. |
| **Consequence right now** | `TH-13` step 3a owed two regression re-runs; `test_dodd_deeds_projected_drive.py` ran green (15 passed / 79 s) and **this one did not run at all**, so the claim "the two degree-2 coil identity tests still fail at 1e-9 exactly as before" is *unverified on this commit*. The direct evidence that the default path is unmoved is control (b)'s bit-identity (0.000e+00 relative), the dodd-deeds module, and the fact that the diff's default branch is the old code verbatim — strong, but not the same thing as the measurement. |
| **Fix** | The 10:30 review discharged the entry's option (a) by reading — the two logs above price the generator at ≈5 s — and scoped the live fix as **`TH-13` step 3a″ (§9 item 1)**: one instrumented re-run with `-s` at the unchanged 570 s ceiling. Green (exit 1 with exactly the two degree-2 identity reds) closes this entry and re-verifies the owed claim; a second exit 124 lands with a phase timeline, and the disposition (cached mesh / coarser rung / module split — each moves a record) becomes a review decision made with data. |

**Step 3a″ reading, 12:00 implementer slot 2026-08-31 — a SECOND exit 124 at
the same 571 s, and the phase timeline now apportions the window: the cost is
entirely the degree-2 pair, and the pre-degree-2 phase is green and cheap.**
Two runs, both `-n 8`, complex build, ceiling unchanged at `timeout -k 30 570`:

- `20260831T170038Z_TH-13-step3a2-coil-degree2-rerun.log` — the item's command
  verbatim plus `-s`, `TH12_STEP2_MODE=full`. **Exit 124 at 571 s**, the same
  wall as 04:38. The `-s` prints land: **mesh 4.3 s, 138 490 cells** (record
  138 490; `near 0.005`, skin depth 15.92 mm, 3.18 cells/δ), then the full cost
  probe (162 558 → 881 476 DOFs, 5.42×; degree-1 summed peak RSS 7.06 GiB;
  projection 47.51 GiB against the 102.40 GiB threshold; **VERDICT: under
  cap**), then `SIGTERM` inside the degree-2 pair. The mesh-time hypothesis is
  now **refuted by direct measurement on the failing run itself**, not only by
  reading other logs.
- `20260831T171059Z_TH-13-step3a2-coil-degree2-probe-phase.log` — the same
  module at `TH12_STEP2_MODE=probe`, which runs mesh + the degree-1 pair + the
  probe and stops before degree 2. **8 passed, 6 skipped, exit 0, 49 s**
  (pytest's own 46.8 s), standard tier.

**The partition.** Pre-degree-2 phase = **46.8 s** (mesh 4.1 s + degree-1
solves **20.6 s + 20.5 s** + probe). So in the killed run the **degree-2 pair
alone consumed ≥ 524 s of the 571 s and did not finish**, whereas on 2026-08-18
the whole module finished in 543 s — i.e. the degree-2 pair then took ≈ 496 s.
The regression is **≥ 5.6% on the degree-2 factorization only**; the mesh
(4.1–4.3 s vs the 4.5 / 4.8 s records) and the degree-1 pair are unchanged.
Nothing "ate" a mesh; the module simply has no margin at degree 2 and this box
is now on the wrong side of it.

**This corrects the `Not` row above.** `TH12_STEP2_MODE=probe` was asserted
there to "not fit either" because it still builds the same mesh. Measured, it
fits with 3.7× to spare (**49 s** against a 180 s standard ceiling) — the mesh
was never the expensive part. That makes a probe-mode row a real standing
regression guard for everything except the two degree-2 identities.

**What this re-verifies, and what it does not.** Green on this commit at
degree 1, on the default (`project_source=True`) path, after `TH-13` step 3a:
the mesh test asserts **138 490 cells**; the degree-1 control reproduces its
record at **+1.5838% vs +1.5834% ΔR deviation → +0.00039 pp** against the
0.01 pp floor; the degree-1 complex-power identity residuals read **8.4704e-15
(loaded) / 3.7068e-15 (free)** against the unloosened 1e-9 bound; `P_loss`
**+1.3876226e-01 W** loaded vs **+0.0000000e+00 W** free; the cost probe
prices degree 2 under cap. So step 3a moved **no degree-1 coil number** — that
part of the owed claim is now measured, not inferred. **Still unverified:** the
two `[loaded-2]` / `[free-2]` degree-2 identity reds at 1e-9. They stay red on
`main` by assumption, and no run since 2026-08-18 has observed them.

**This entry stays open.** Per §9 item 1 the ceiling was **not** raised and the
`full` command was **not** retried a third time. The disposition is a review
decision, and it now has data: (a) split the module so the degree-2 pair is its
own heavy-tier run with a real margin; (b) cache the mesh (saves 4 s — worthless
on this evidence); (c) a coarser degree-2 rung, which moves the +1.5834%
record; or (d) accept the degree-2 identities as unobservable in a scheduled
slot and gate the module in probe mode. The measured partition says (a) or (d);
(b) is dead.

**Ruled 2026-08-31 18:00 review: option (a), as `TH-13` step 3a‴ (§9 item 3).**
The degree-2 pair moves to a new `test_coil_loading_degree2_pair.py` and runs
**one σ-half per window** (`TH12_DEGREE2_HALF=loaded|free`; mesh + degree-1
row + cost probe + one degree-2 solve, ≈ 310–330 s at `-n 8`) under
`timeout -k 30 600` — the binding limit for a scheduled slot is the 660 s
Bash-tool window a foregrounded harness run must return inside, which is why
"its own heavy-tier run" alone would not do. The original module defaults to
probe mode (49 s) as the standing guard. **The 570 s ceiling is not raised; the
case is shrunk.** This entry closes with the commit in which both halves return
footers inside their ceilings — expected `exit 1` with exactly one degree-2
identity red each, which is the owed observation — and the degree-2 identity
entry is re-pointed at the new module in the same commit. A third exit 124
(one factorization alone > 550 s) is a stop, and the disposition — then (d) —
returns to the review.

### ✅ RETIRED 2026-10-07 (`OPS-68`) — was Four defects surfaced by replacing finiteness-only tests with real anchors (`OPS-17` step 2, 2026-08-17)

> **Retired 2026-10-07 by `OPS-68` (§9 item 307).** All four defects carry their resolution in the body and each fix's log is green: (1) `GEO-17` step 1 (`81138ba4`, `20260820T110549Z_GEO-17-step1-final.log.gz`, 3 passed, `Status: 0`; `GEO-17` ✅ §7 :1842); (2) `MAG-17` step 1 (`fe1cfb43`, `20260820T123823Z_MAG-17-step1-final2.log.gz`, 6 passed, `Status: 0`; `MAG-17` ✅ §7 :1666); (3) `POST-5` step 4 (`6da9897a`, `20260819T201005Z_POST-5-step4-smoke-final.log.gz`, 12 passed, `Status: 0`); (4) `POST-5` step 1 (`6044a61f`, `20260818T215117Z_POST-5-step1-negcontrol.log.gz`, 8 passed, `Status: 0`); `POST-5` ✅ (§7 :2585). Body kept verbatim below.

**Found:** 2026-08-17, `OPS-17` step 2, at commit `197142f`. All four were
invisible until the tests that exercise these paths were given quantitative
anchors — which is the whole premise of the chunk. **None is being fixed here:**
`OPS-17` is test hygiene, and each of these belongs to the subsystem it lands
in. The first three are carried in the tree as `pytest.mark.xfail(strict=True)`
with the measurement in the docstring, so a fix reports as XPASS rather than
silently passing; the fourth is worked around in one test.

**1. ~~`coil_phantom_domain` region-resolution policy shrinks the meshed coil
volumes by ≈22% while asking for a *finer* mesh.~~ ✅ RESOLVED 2026-08-20
(`GEO-17` step 1, 06:00 implementer slot) — the per-region sizes were never
applied at all: `coil_phantom_domain` walked volume → surfaces → curves →
points through `gmsh.model.getBoundary` with its default `combined=True`,
whose result for a volume's closed shell of surfaces is empty, so every region
collected zero CAD points and `mesh.setSize` was never called
(`20260820T110127Z_GEO-17-step1-diag.log`: `air: 0 pts -> NO SIZE SET` for all
four regions at both sizings). Only the global `CharacteristicLengthMin/Max`
clamps survived — uniform `[0.015, 0.015]`, policy `[0.010, 0.020]` — so the
policy run meshed the coil at the *air's* 0.020 ceiling. **The hypothesis
below is refuted:** no field won the interface, because no field existed. The
sizes are now carried by a `Min` over per-volume `Constant` size fields, and
the volumes move the way refinement requires (final log
`20260820T110549Z_GEO-17-step1-final.log`, `-n 2`, 13 s): coil_1
1.191750413e-04 → **1.319468693e-04 m³** (+10.72%, meshed/CAD 0.754685 →
**0.835563**), coil_2 → 1.316573175e-04 (+10.79%, 0.752565 → 0.833730),
phantom → 4.990112950e-04 (+0.94%, 0.983531 → 0.992751), and the air — the one
region the policy *coarsens* — is the one region that loses volume (−0.26%).
The uniform column reproduces the table below to every printed digit, gated in
the test as the negative control. The test is a plain gate again, renamed
`test_region_resolution_policy_refines_the_tagged_volumes_toward_cad`; its 5%
band was replaced with the sign-of-refinement identity and the meshed/CAD
recovery bounds, because the old band's premise (region sizing "must not move
the geometry") is false for a curved region — see that test's docstring.**
`tests/mesh/test_mesh_tag_integrity.py::test_region_resolution_policy_does_not_move_the_tagged_volumes`
(xfail). Measured at `-n 2`, `20260817T111054Z_OPS-17-step2-mesh-n2.log`,
uniform `h = 0.015` against `coil_resolution=0.012, phantom_resolution=0.010,
air_resolution=0.020`:

| tag | uniform [m³] | policy [m³] | Δ |
| --- | --- | --- | --- |
| 1 `coil_1` | 1.191750413e-04 | 9.333354960e-05 | −21.68% |
| 2 `coil_2` | 1.188402981e-04 | 9.195675344e-05 | −22.62% |
| 3 `phantom` | 4.943767949e-04 | 4.880940997e-04 | −1.27% |
| 4 `air` | 1.143560787e-02 | 1.149461560e-02 | +0.52% |

CAD torus volume is `2π²Rr² = 1.579137e-04 m³`, so uniform recovers 75.5% of
each coil and the policy mesh 59.1%. **The sign is the defect:** a linear-tet
mesh inscribes a curved surface, so refining can only move meshed volume *up*
toward CAD. Every region given a finer size lost volume, and the air took up
exactly what the curved regions lost. **Cause: not diagnosed.** Hypothesis —
the region size fields replace rather than refine the surface sizing on shared
curved interfaces, so the coarse air field (0.020) wins on the coil and phantom
boundaries. **Resolves with:** a `GEO` chunk on `coil_phantom_domain`'s sizing
path. Both meshes still partition their own volume to 1e-9, so this is fidelity,
not conformity.

**2. ~~The Coulomb-gauge Lagrange multiplier does not vanish for a
divergence-free source.~~ ✅ RESOLVED 2026-08-20 (`MAG-17` step 1, 07:30
implementer slot) — candidate (a) confirmed, candidate (b) excluded: **the
anchor was wrong, not the constraint block.** The h-ladder ran
(`20260820T123307Z_MAG-17-step1-ladder.log`, `-n 2`, 95 s), reproducing the
record at its own mesh and refining twice:

| h | cells | multiplier spread |
| --- | --- | --- |
| 0.0050 | 29 190 | 7.836781e+00 |
| 0.0035 | 82 819 | 3.052022e+00 |
| 0.0025 | 208 049 | 1.438617e+00 |

Fitted log-log rate **2.4476** (pairwise 2.645 / 2.234) against the
pre-registered bands rate ≥ 0.7 ⇒ DISCRETE-SOURCE, |rate| < 0.3 ⇒
ASSEMBLY-DEFECT. The verdict is DISCRETE-SOURCE with four sigfigs of margin,
and superlinearly so: `p` is absorbing the interpolated `J`'s discrete
divergence, which is a mesh residual that converges away, not a defect in how
the constraint is assembled. The consequence is that the `OPS-17` anchor
("spread → 0 to solver tolerance") **cannot hold on any single mesh** — the
right anchor is a convergence rate. The strict xfail is retired and the claim
moved to
`tests/solver/test_gauge_multiplier_convergence.py::test_multiplier_spread_converges_for_a_divergence_free_source`,
a plain gate on the ladder at the unmoved pre-registered 0.7 (band deliberately
not tightened to the measurement — it stays the discriminator it was designed
as). The negative control holds in the same run: the incompatible straight
wire stays at its recorded 2.083064e+02 scale, > 10× the loop's base-h spread
(recorded separation 26.6×), so the multiplier has not stopped discriminating
compatible from incompatible sources. Final run
`20260820T123823Z_MAG-17-step1-final2.log`, 6 passed, 97 s.**
`tests/solver/test_gauge_lagrange.py::test_gauge_multiplier_vanishes_for_a_divergence_free_source`
(xfail). Measured at `-n 2`, `20260817T111217Z_OPS-17-step2-solver-n2.log`:
spread **7.836781e+00** on a closed current loop (azimuthal J: `div J = 0`
inside, `J·n = 0` on both the torus surface and the outer sphere) against the
solver-tolerance residual the theory requires. The multiplier is not dead — it
reads **2.083064e+02** on the deliberately incompatible straight wire, 26.6×
larger — but "vanishes for a compatible source" is false as written.
**Cause: not diagnosed.** Two candidates, separable with one h-ladder:
(a) benign — the source enters through an interpolated `J` whose discrete
divergence is O(h), and the fixture is coarse (h = 0.005 against a 0.003 wire
radius), in which case the spread falls with refinement; (b) a real defect in
how the constraint is assembled, in which case it is h-independent.
**Resolves with:** a `MAG`/`OPS` chunk on the gauge formulation. Note `max|A|`
is not a usable normaliser here — `test_lagrange_removes_the_null_space`
requires the LAGRANGE `max|A|` to sit six orders below the penalty solve's.

**3. ~~Real Poynting power does not balance on the time-harmonic smoke fixture,
and the boundary flux has the wrong sign.~~ ✅ RESOLVED 2026-08-19 (`POST-5`
step 4, 15:00 implementer slot) — see the closing block at the end of this
entry; the test is a plain gate again and passes.**
`tests/solver/test_time_harmonic_smoke.py::test_time_harmonic_smoke_solve_conserves_real_power`
(xfail). Measured at `-n 2`,
`20260817T112448Z_OPS-17-step2-th-smoke2-n2.log`: dissipated
`½∫σ|E|²dV = +1.199162e-06 W` against net inward flux
`−∮½Re(E×H̄)·n̂dS = −2.008179e-07 W`, relative imbalance **116.7465%** against
a pre-stated 25%. Power leaves through the boundary while the medium
dissipates, which the identity forbids for any solution of Maxwell's equations
regardless of boundary condition. **Cause: not diagnosed.** Candidates:
(a) resolution — 0.16 m domain at h = 0.03 m is ≈9 cells per in-medium
wavelength (λ = c/(f√78) = 0.266 m), and the boundary leg is a curl trace, the
least accurate quantity a degree-1 N1curl solution carries;
`tests/validation/test_poynting_balance.py` needs a refined mesh to reach 5%
and gates the *convergence* of this imbalance for that reason; (b) the source —
an axial current terminating on the end caps, so `J·n ≠ 0` there, the same
incompatibility as defect 2. **Resolves with:** a `TH`/`POST` chunk; an
h-ladder on this fixture distinguishes the two in one command.

> **Sharpened 2026-08-18, `POST-5` step 1 — candidate (a) is excluded and so
> is a third candidate (c).** The h-ladder ran
> (`20260818T215101Z_POST-5-step1-ladder2.log`, `-n 2`, 5 s):
>
> | h | cells | dissipated [W] | net inward [W] | imbalance |
> |---|---|---|---|---|
> | 0.030 | 1 405 | 1.199162e-06 | −2.008179e-07 | 116.7465% |
> | 0.020 | 2 590 | 1.154337e-06 | −1.778362e-07 | 115.4059% |
> | 0.015 | 4 661 | 1.479920e-06 | −2.134447e-07 | 114.4227% |
>
> Fitted rate in h **0.0290** against the pre-registered ≥ 0.7, and the net
> inward flux is negative on every rung including the finest. **It is not
> resolution.** Candidate (c), a flipped outward measure, is excluded
> exactly: `∮x·n̂dS / 3|Ω| = 1.000000000000` on this fixture with the same
> `dx`/`ds` pair the balance uses
> (`test_smoke_fixture_boundary_measure_is_outward_oriented`). Candidate (b)
> — the drive's `J·n ≠ 0` on the end caps — is what remains, together with a
> defect in the boundary leg's assembly. The xfail is **not** rescoped to
> convergence; its reason string now carries these numbers. `POST-5` step 2
> in PROJECT_PLAN §7 scopes the next discriminator (a closed azimuthal
> source on the same fixture).

> **Sharpened again 2026-08-19, `POST-5` step 2 — candidate (b) is excluded
> too; the defect is in the boundary leg's assembly.** The closed-drive
> discriminator ran on the coarse rung
> (`20260819T051150Z_POST-5-step2-closed-drive2.log`, `-n 2`, **4 s**):
>
> | drive | dissipated [W] | net inward [W] | sign | imbalance | blind diss [W] |
> |---|---|---|---|---|---|
> | axial (record) | 1.199162e-06 | −2.008179e-07 | − | 116.7465% | 0.000000e+00 |
> | closed azimuthal | 4.778876e-09 | −2.849722e-10 | − | 105.9632% | 0.000000e+00 |
>
> The azimuthal drive `J = (−y, x, 0)/a` is a *closed* source on this fixture
> — `div J = 0` pointwise, `J·n = 0` on both end caps and on the rod's own
> lateral surface, so the tag restriction adds no surface divergence either —
> and it is P1-exact, so no interpolation error is folded in
> (`_azimuthal_current` in the smoke module). Both halves of the
> pre-registered SOURCE band fail: the imbalance moves only
> 116.7465% → **105.9632%** (a 10.8 pp move, against a band that required
> < 25%) and the flux **sign does not turn positive**. The axial drive re-run
> in the same session reproduces the step-1 record to `rtol=1e-6` on all
> three numbers (asserted, not eyeballed), so the drive is the only thing
> that changed between the two rows. **Verdict: ASSEMBLY.** The source's
> `J·n ≠ 0` incompatibility is real (defect 2 measures it) but it is not what
> makes this identity fail. **Resolves with:** the boundary-leg probe scoped
> as **`POST-5` step 3** (named in the step-2 §7 entry; queued 2026-08-19
> 03:00 review) — scoring the curl trace `−∮½Re(E×H̄)·n̂dS` against the
> `TH-6` lossy plane wave, where both legs have closed forms, which separates
> a wrong `H = ∇×E/(−jωμ)` reconstruction from a wrong facet assembly. The
> xfail keeps its 25% band and `strict=True`.

> **Resolved as a defect in the *identity*, not the code — 2026-08-19,
> `POST-5` step 3.** Step 2's ASSEMBLY verdict is **overturned by direct
> measurement**: the boundary leg is sound, and what is wrong is that
> `poynting_power_balance` scores the **source-free** identity on a
> **driven** fixture.
>
> *Leg 1 — the boundary leg against its own closed form*
> (`20260819T123438Z_POST-5-step3.log`, `-n 2`;
> `test_each_leg_scored_against_its_own_closed_form`). On the `TH-6` lossy
> plane wave both legs have closed forms —
> `P_flux = ½βL²(1−e^{−2αL})/(ωμ₀μᵣ)`, `P_diss = ½σL²(1−e^{−2αL})/(2α)`,
> equal identically because `k² = k₀²ε_c` gives `2αβ = ωμ₀σ`
> (asserted separately at `rtol=1e-12`: both sides 7.060162290693e+02;
> analytic value 1.241101e-04 W):
>
> | rung | cells | flux leg [W] | flux err | volume leg [W] | volume err |
> |---|---|---|---|---|---|
> | 12³ | 10 368 | 1.140318e-04 | 8.1205% | 1.241984e-04 | 0.0711% |
> | 24³ | 82 944 | 1.190042e-04 | **4.1141%** | 1.241317e-04 | 0.0174% |
>
> The boundary leg is inside the pre-registered 10% band on the fine rung and
> falls at rate `log₂(8.1205/4.1141) = 0.981` — clean O(h) for a degree-1
> N1curl curl trace. `H = ∇×E/(−jωμᵣμ₀)` and the facet assembly are
> **correct**; there is no factor and no conjugation error to find.
>
> *Leg 2 — the full three-term balance on this fixture*
> (`20260819T124405Z_POST-5-step3-source.log`, `-n 2`, 4 s;
> `test_the_missing_impressed_source_term_accounts_for_the_smoke_imbalance`).
> With an impressed `J`, Poynting's theorem reads
> `−∮½Re(E×H̄)·n̂dS = ½∫σ|E|²dV + ½Re∫E·J̄dV`, and the helper omits the
> second term:
>
> | drive | dissipated [W] | net inward [W] | source ½Re∫E·J̄ [W] | two-term | three-term |
> |---|---|---|---|---|---|
> | axial | 1.199162e-06 | −2.008179e-07 | −1.199162e-06 | 116.7465% | **16.7465%** |
> | azimuthal | 4.778876e-09 | −2.849722e-10 | −4.778876e-09 | 105.9632% | **5.9632%** |
>
> Both inside the pre-registered 25% band (the xfail's own), so the omitted
> term **is** the imbalance. **Read the third column carefully**: the source
> term equals `−dissipated` to all seven printed digits on both drives. That
> is not a coincidence and not evidence about the flux — under the *natural*
> boundary condition this fixture uses (`TimeHarmonicBoundaryCondition.NATURAL`,
> the default; the `TH-6` fixture is PEC-with-Dirichlet-data and source-free,
> which is why its 5% gate is honest) the weak form tested with `v = Ē` has no
> boundary term, so `½∫σ|E|² + ½Re∫E·J̄ = 0` holds *algebraically* in the
> discrete solution. The three-term residual is therefore exactly the boundary
> flux over the scale — i.e. 16.7% / 6.0% is the discretisation error of the
> curl trace on a ≈9-cells-per-wavelength gmsh mesh, entirely consistent with
> leg 1's 8.1% at 10 368 cells on a *structured* box.
>
> **Consequence.** The wrong sign is not forbidden: `−∮½Re(E×H̄)·n̂dS` alone
> has no sign law when a source is present inside — only the three-term sum
> does. The original entry's "which the identity forbids for any solution of
> Maxwell's equations regardless of boundary condition" is **wrong as
> written**; it is true only for a source-free domain. **Resolves with:**
> `POST-5` step 4 (scoped, not executed) — teach `poynting_power_balance` the
> impressed-source term and re-gate. Nothing was changed in this step: the
> xfail keeps its 25% band and `strict=True` and still XFAILs.

> **✅ CLOSED 2026-08-19, `POST-5` step 4 (15:00 implementer slot).** The fix is
> in the helper, not in the solve: `poynting_power_balance` now accepts the
> impressed `current_density` (and the `source_measure` it was assembled on),
> assembles `source_power_w = ½Re∫E·J̄dV`, and scores `relative_imbalance` on
> the three-term statement when a drive is given. Omitting the drive changes
> nothing — the source-free two-term identity is still what a source-free
> domain is scored against, and it stays reachable in both cases as
> `two_term_relative_imbalance` / `two_term_power_scale_w`, so the step-1
> h-ladder journal above keeps reconciling.
>
> `test_time_harmonic_smoke_solve_conserves_real_power` is a **plain gate**
> again — the `xfail(strict=True)` is gone and it PASSES against the *unmoved*
> 25% band (`20260819T201005Z_POST-5-step4-smoke-final.log`, `-n 2`, complex
> build, 12 passed / exit 0 / 8 s):
>
> | reading | value |
> |---|---|
> | dissipated ½∫σ\|E\|²dV | 1.199162e-06 W |
> | net inward −∮½Re(E×H̄)·n̂dS | −2.008179e-07 W |
> | source ½Re∫E·J̄dV | −1.199162e-06 W |
> | three-term residual (gated, band 25%) | **16.7465%** |
> | two-term reading (step-1 record 116.7465%) | 116.7465% |
> | σ-blind three-term control | 83.2535% |
>
> Every one of these reproduces the step-3 record; the test asserts them at
> `rtol=1e-6` (powers) and `atol=1e-6` (the two imbalances, which the record
> carries only to 4 decimals as a percentage) rather than printing them.
> The σ-blind control is now scored on the three-term residual too, where its
> arithmetic ceiling is `1/0.167465 = 5.97×`; the pre-registered replacement
> for the old (never-met) 10× factor is **3.0×**, and it must also be *rejected*
> by the very band the honest solve passes — 83.2535% is 4.97× the honest
> reading and well outside 25%.
>
> **Negative control, on the fixture where J = 0**
> (`20260819T200651Z_POST-5-step4-negcontrol.log`, `-n 2`, 15 passed / exit 0 /
> 152 s): `test_zero_impressed_current_leaves_the_source_free_balance_untouched`
> solves the `TH-6` plane wave at 12³ and scores it twice — with no drive and
> with `J = fem.Constant(msh, [0,0,0])` (a `Constant` rather than a literal, so
> the integral is genuinely *assembled* rather than folded away). The source
> term is **exactly `0.0` W**, asserted `== 0.0`, and all seven other returned
> quantities are asserted **bit-identical** between the two calls: 8.185716%
> both ways, which is the step-3 12³ rung's 8.1857% unmoved. All the `POST-3`
> gates in that file (5% MVP, piecewise σ, μᵣ-field, and the three blind
> controls) are green in the same log.

> **JIT trap, 2026-08-19 (`POST-5` step 2).** The step's first window
> (`20260819T050314Z`, exit 124 at 400 s) died with rank 1 parked in
> `MPI_Bcast` — the dolfinx cold-JIT signature — and left **one 0-byte
> `.c` in `/root/.cache/fenics`** (created 7 s into the run, i.e. long
> before the form that was waiting on it was reached). Removing that single
> entry (`rm /root/.cache/fenics/*<hash>*`) and re-running the identical
> command took **2.94 s of pytest**. So a 0-byte cache entry is not only the
> *symptom* of a killed compile, it is a live lock the rest of that same run
> can block on: when a run stalls in `MPI_Bcast`, look for `find
> /root/.cache/fenics -size 0` **before** concluding the case is too big for
> the window or that a form's quadrature degree is unpinned.

**4. ~~`poynting_power_balance` raises on a scalar `sigma=0.0`~~ — FIXED
2026-08-18 (`POST-5` step 1).** The scalar branch is now wrapped in
`fem.Constant(msh, dolfinx.default_scalar_type(σ))`, so the integral keeps
its domain and `sigma=0.0` assembles to **exactly 0.000000e+00 W** (verified
at three mesh sizes, `20260818T215101Z_POST-5-step1-ladder2.log`; asserted
`== 0.0`). The `SIGMA_BLIND = 1e-12 * SIGMA` workaround is deleted and the
control is now a real zero. `tests/validation/test_poynting_balance.py` is
unmoved by the wrap — 8 passed, scalar-vs-DG0-field paths still equal at
`rtol=1e-12` (`20260818T215117Z_POST-5-step1-negcontrol.log`). Original
report follows.

**`poynting_power_balance` raises on a scalar `sigma=0.0`, the σ-blind
negative control its own docstring advertises.**
`src/fem_em_solver/post/power_balance.py:137`. `0.5 * 0.0 * ufl.inner(E, E)`
folds to a domain-less UFL zero and `* ufl.dx` then raises
`ValueError: This integral is missing an integration domain`. First hit at
`20260817T112414Z_OPS-17-step2-th-smoke-n2.log`. **Cause:** the scalar branch
passes a bare Python float rather than wrapping it in `fem.Constant(msh, ...)`,
so UFL constant-folds the whole integrand away. `sigma=0.0` is an intended
input — the module docstring calls it "what makes the σ-blind negative control
possible". Existing callers pass a non-zero scalar or a `sigma_field`, so
nothing was previously red. **Worked around**, not fixed: the smoke test above
uses `SIGMA_BLIND = 1e-12 * SIGMA`. **Resolves with:** wrapping the scalar in
`fem.Constant` — a one-line `POST` fix plus a re-run of
`tests/validation/test_poynting_balance.py`.

### ✅ RETIRED 2026-10-07 (`OPS-67`; repaired 2026-08-28 by `OPS-28`) — `PORT-1` step 4's `allgather` reduction broke `_DummyComm`, and one orientation test regressed silently on 2026-08-13 (`OPS-17` step 3, 2026-08-17)

> **Retired 2026-10-07 by `OPS-67` (§9 item 303).** Its own Resolves-with — "any earlier commit that repairs the double", with entry 3 re-symptomed in the same commit — was met by `OPS-28` (`fd1479ca`): `_DummyComm` gained `allgather`, and entry 3 was re-dated with a measured correction in that commit. Bracketed by `20260828T033037Z_OPS-28-red-baseline.log` (the `AttributeError` on both names) and `20260828T033055Z_OPS-28-gate.log` (no `AttributeError`; `…induced_voltage_sign` green; the remaining two reds are entry 3's, Status 1 by those alone). The archive's `✅ RETIRED 2026-08-28` entry (`known-issues-archive.md` ≈ :588) is the **same defect** found a second time by the `OPS-26` census, under a different heading and body — so this live copy was a stale duplicate, not a distinct open issue. The orientation module has since been rewritten without the double (`PORT-8`, `5418f66b`). Body kept verbatim below.

**Found:** 2026-08-17, `OPS-17` step 3, at commit `e211356`, in the first
*completed* real-mode suite leg since 2026-08-13
(`20260817T201248Z_OPS-17-step3-real-nonvalidation-n2.log`, `-n 2`, 218 s,
3 failed / 134 passed / 32 skipped / 2 xfailed, both ranks identical).

| | |
|---|---|
| **Tests** | `tests/ports/test_port_orientation_sensitivity.py::test_port_orientation_flip_changes_induced_voltage_sign` — **not previously known-red**<br>`tests/ports/test_port_orientation_sensitivity.py::test_port_orientation_flip_changes_off_diagonal_sparameter_sign` — known-red under entry 3, but **not for entry 3's recorded reason** |
| **Symptom** | `AttributeError: '_DummyComm' object has no attribute 'allgather'` raised out of `src/fem_em_solver/ports/excitation.py:258`. The failure is inside `src/`, *before* any test assertion executes. |
| **Cause** | Diagnosed, one line. `PORT-1` step 4 (2026-08-13) added `problem.mesh.comm.allgather(...)` at `excitation.py:258` to reduce the tag set before `validate_required_port_tags_exist` — the documented fix for entry 6 defect (2). The file's test double (`tests/ports/test_port_orientation_sensitivity.py:16`) implements only `rank` and a static `allreduce`, so the new collective call has nothing to dispatch to. |
| **Why it went unseen for four days** | `PORT-1` step 4's own gate ran the *two-torus* negative control, not this file, and every slot since has used targeted per-file runs. No completed suite leg was paid for until this one. The `--ignore=tests/validation` leg above is the cheapest thing that catches this class (218 s). |
| **Two separate consequences** | (i) `test_port_orientation_flip_changes_induced_voltage_sign` is absent from entry 3's list of two tests, so it was green before 2026-08-13 — this is a **silent regression**, not a pre-existing red. (ii) Entry 3's symptom line (`assert np.all(np.abs(diagonal) > 0.0)` on a zero diagonal) is **stale for the orientation test**: that assertion is now unreachable. It still describes `tests/ports/test_sparameter_assembly.py::test_n_port_sweep_assembles_finite_matrix_with_expected_shape` correctly (`test_sparameter_assembly.py:104: AssertionError: assert False`), which is the third failure in the same leg. |
| **Fix** | **Deliberately not fixed by `OPS-17`** — that chunk is test hygiene, and entry 3's standing disposition is that these tests live and die with `PORT-1`'s retirement of the `PORT-0` placeholder. The mechanical repair is two lines in the test double (`allgather = staticmethod(lambda v: [v])`, or use a real `MPI.COMM_SELF`); whoever retires `PORT-1` should decide whether the double survives at all. Do not read the `AttributeError` as evidence about the placeholder's arithmetic — it never runs. |
| **Resolves with** | `PORT-1`'s retirement commit, or any earlier commit that repairs the double. **Entry 3 must be re-symptomed in the same commit**, since its recorded cause now applies to only one of its two tests. |

### ✅ RETIRED 2026-10-07 (`OPS-68`) — was `EX-49`'s negative-control ceiling-first floor (≥2× the band) was a plan-time estimate that the measured fixture falls 1.2% short of — the primary claim is untouched (2026-09-05, implementer slot)

> **Retired 2026-10-07 by `OPS-68` (§9 item 307).** The entry's own Resolves-with reads "Nothing further needed; this is a closed observation, not an open chunk", and the 2026-09-05 18:00 ruling in its body upheld the demotion. `EX-49` is ✅ and audited PASS (§7 :2989; `c3cd9171`, `20260905T201151Z_EX-49.log.gz`, `Status: 0`). The printed 1.9754× margin stays recorded in the body below. Body kept verbatim below.

| | |
| --- | --- |
| **Chunk** | `EX-49` (`examples/ports/13_birdcage_asymmetric_drive.py`), the negative control for the linear (`w_lin = (1,1,0,0)/√2`) drive's C4 covariance. |
| **Verified at** | `20260905T201151Z_EX-49.log`, Status 0, 73 s at `-n 2`, complex build, 2026-09-05. |
| **What the §7 item pre-registered** | "assert ≥ 2× the band, none larger claimed" — reasoning from `POST-6`'s gate module cw-sense/mirror control on the same fixture, which reads 95.1975%. That control is a *different comparison* (sense-swap and mirror, not a same-drive rotation) and was never measured for a two-port linear drive's own rotation-covariance before this item ran. |
| **What is measured** | The linear drive's own C4 covariance (`|B1+|_lin` at points rotated 90° vs unrotated, on the same 51 phantom centroids `ports:7`/`POST-6` use): **9.8768% = 1.9754× the imported `C4_COVARIANCE_BAND` (5%)**. The primary claim — `> C4_COVARIANCE_BAND`, i.e. the linear drive is *not* C4-invariant — holds comfortably. The pre-registered ≥2× (10.0%) floor misses by 1.2% relative. |
| **Why this is not loosened** | `C4_COVARIANCE_BAND` (an imported, gated band) is untouched. The 2× multiple was never an imported band or a measured record from any gate module — it was this item's own plan-time arithmetic, carried over from an unrelated comparison. The example now measures and prints the margin rather than asserting a number it does not reach; asserting it anyway would have been the failing-test-loosening this repo forbids in the other direction (silently passing on a number the run does not support). |
| **Consequence** | None for `EX-49`'s closure — the actual negative-control claim in the §7 item (`> C4_COVARIANCE_BAND`) is asserted and holds. Noted here only so a future run's margin (currently 1.9754×) is checked against this one rather than silently drifting, and so the "≥2×" language in the §7/§9 text is read as measured-not-guaranteed. |
| **Resolves with** | Nothing further needed; this is a closed observation, not an open chunk. Remove this entry only if a future review wants the record gone — it costs nothing to leave standing as provenance for the printed margin. **Review ruling 2026-09-05 18:00 — the demotion stands.** The ≥2× floor was a plan-time prediction carried from a different comparison (`POST-6`'s cw-sense/mirror control), never a measured record on this comparison and never an imported band; the gated claim `> C4_COVARIANCE_BAND` is asserted and holds at 1.9754×. The never-loosen rule protects bands and records that were once green, and this was neither. What the episode did expose is a wording defect in the item, not in the run: §9 standing rule (e) (added this review) now requires every negative-control factor to be labelled **asserted** (a prior measurement of the *same* comparison on the *same* fixture backs it) or **predicted** (printed beside the measured factor, never asserted), and an executor that meets a failing pre-registered assertion whose label is missing stops and reports the negative result rather than deciding in-slot. |

## 2026-09-19 — ✅ RETIRED 2026-10-05 (weekly review, T6 decided on F18a / F18b) — was 🟡 OPEN (`ANS-4` step 3 / 3c private readouts, weekly review) — the order-matched `ANS-4` comparison **disagrees on the self class at 64 MHz and 10 MHz**; the 128 MHz AGREE stands

| | |
| --- | --- |
| **Symptom** | None in CI — every imported `PORT-9`/`PORT-11` gate is green on every degree-2 rung (`xl-ledger.md` rows 2026-09-16/17/18). The disagreement is against the external code only. |
| **Ruling** | 2026-09-13's pre-registered private rule for the 64 MHz order-matched rung **fired on the self class** (couplings inside it); 3c's public rule also fired — the 10 MHz degree 1 → 2 move (4.41 / 1.35 / 1.08 %) is in the same class as the Larmor moves, so the 09-06 "10 MHz AGREE", read at degree 1, was partly coincidental. 3b's rule did not fire: on the C4-congruent cut the 128 MHz degree-2 classes sit within step 2d's own spreads of step 2d's (≤ 0.083 % against ≤ 0.142 %), and the 128 MHz AGREE stands as written. Numbers: gitignored `docs/private/ans3-ans4-adjudication-2026-09-19.md`. |
| **What is and is not known** | The order-matched disagreement is concentrated in `S₁₁` and **falls** with frequency, which a fixed-`h` discretisation error does not do; `ANS-2`'s 09-18 ruling independently named the lumped-port feed as its residual at 10 MHz. But at 10 and 64 MHz only one degree-2 rung exists, so neither figure is yet an `h`-converged value. |
| **Consequence** | No absolute `S₁₁` (hence no `Z_in`, no tuned-`S₁₁` absolute, no match claim) is licensed at 10 or 64 MHz at either order. Coupling classes: AGREE at all three frequencies. The production element order is not changed on this evidence. |
| **Next** | `xl` entries 7 (3e, 64 MHz degree-2 `h`-ladder, queued 09-22) and 11 (3g, the same at 10 MHz) decide whether the rungs are converged; then `PORT-21` (§7) and `ANS-6`'s PEC-coil column (which removes the resolved conductor interior from both codes). |
| **Public sensitivity table (2026-09-21, `PORT-21` step 1b)** | On `main`: `tests/validation/test_port21_feed_sensitivity.py` (gated on `FEM_EM_PORT21_FREQUENCIES_MHZ`); the table as printed is in `docs/testing/logs/20260921T094925Z_PORT-21.log:5795–5810` (10 / 64 / 128 MHz, `-n 2`, degree 1). Self-class `\|ΔS\|/\|S\|`: (a) width-told-as-`w/(1+κ)` 3.411e-02 / 1.552e-02 / 6.327e-03 (falls with f); (c′) halved conductor grading on the C4-congruent cut, against its own congruent control′, 1.885e-02 / 1.538e-02 / 1.679e-02 (non-monotonic). Our own convention sensitivities only — no absolute `S₁₁` is licensed by it. |
| **Re-read 2026-09-26 (weekly; `xl` 3e 09-22, 3g 09-24; `ANS-6` adjudicated)** | **Still OPEN, and re-read as unconverged on our side rather than as a model difference.** The degree-2 ladders at 10 and 64 MHz move monotonically **toward** the order-matched AED column on all three classes; the finest-rung gaps are 5.7 / 1.8 / 1.5 % (10 MHz) and 2.5 / 1.5 / 1.1 % (64 MHz) and each is only ≈ 1.6× its ladder's last step. The ladder refines the global `resolution` with `conductor_resolution` pinned, so the port neighbourhood is never refined and the "ratio ≳ 1.5" reading is withdrawn (§10 epitaph 2026-09-26). `ANS-6`'s copper and PEC columns carry the same gap, which clears the resolved conductor interior. Variant (a) of the table above has the gap's frequency shape and points toward AED (a hypothesis until measured at degree 2). **Consequence unchanged:** no absolute `S₁₁` / `Z_in` at 10 or 64 MHz at either order. **Next:** §10 chain F — F6 / F7 (variant (a) at degree 2 on the finest congruent rung), F8 (the conductor-side ladder), `OPS-61` (Palace, a public second reference). Numbers: gitignored `docs/private/ans6-ans4-ans2-adjudication-2026-09-26.md`. |
| **Re-read 2026-09-30 (weekly; `xl` 3h 09-27, 3i 09-28, 3j 09-30, `ANS-6` 2c 09-29; `PORT-21` step 3's degree-2 window)** | **Still OPEN; a leading hypothesis now has one measured point.** The three ladders on the congruent cut reproduce their cut-off twins: the cut is not a term in the h-statement, and no ladder of this family converges against its gap (last change on `S₁₁` 3.56 / 1.49 / 0.63 % against gaps of 5.7 / 2.5 / 1.1 %) — the family is closed to further windows. The hole-mesh degree-2 rung at 64 MHz sits 3.1 / 1.8 / 1.5 % out, the same on Cu and PEC. **The κ reading:** at degree 2 on the ×1 congruent cut at 10 MHz, variant (a) — the sheet width told as `w/(1+κ)` with the in-run κ = 5.21 %, no AED input — takes the self-class miss against the order-matched AED column from 15.0 % to 0.8 %, collinear with the gap (`20260927T110931Z_PORT-21-step3-W3.log:2396,2914,2942`, read privately). One frequency, one mesh, a 4-digit print: a hypothesis, tested by §10 F11 (64 / 128 MHz on the ×1 cut), F12 (`xl` entry 22, h = 0.0075) and F7 (`xxl` entry 23, the finest rung). **Consequence unchanged:** no absolute `S₁₁` / `Z_in` at 10 or 64 MHz; the production order is not flipped (§10 T6 rule). Numbers: gitignored `docs/private/ans7-ans6-adjudication-2026-09-30.md` (the κ arithmetic is in the session's unlogged read; the gap levels above are its whole public content). |
| **Re-read 2026-10-03 (weekly, F9; transcribed by the 2026-10-04 03:00 daily)** | **Still OPEN — retires at T6, not before.** On exact SPEC port rectangles with κ off (operator ruling 2026-10-02; item 192), our degree 2 on the ×1 mesh is within **0.73 / 0.44 / 0.34 %** (self class, 10 / 64 / 128 MHz) of AED First Order and within 0.55 % on every coupling class; Palace reads the same mesh to 0.0000 % of ours at matched order (F17). The self-class gap was the ragged sheet. **One mesh:** `xl` 10-02 (`ANS-4` 4b) moved the ragged degree-2 self class +1.92 % under conductor halving, so T6 (i′) / (i″) — F18a / F18b on exact rectangles — decide. **Consequence unchanged** until then. Numbers: gitignored `docs/private/f9-t13-ans8-adjudication-2026-10-03.md`. |
| **Retired 2026-10-05 (weekly, T6 decided on F18a / F18b)** | On exact SPEC rectangles, κ off, our degree 2 is within **1 %** of AED First Order on every class at 10 / 64 / 128 MHz on three meshes — ×1 (0.73 % worst class), h = 0.0075 (`xl` entry 27, 0.97 %) and conductor ×0.5 (items 212 / 225, 0.77 %) — and equals Palace at matched order. The 09-19 disagreement was the ragged port sheet at degree 2 and element order at degree 1; `PORT-21` ✅, §2.2 rewritten. Not retired as an h-converged value: two h-rungs, and at 10 MHz the self class moves 0.24 % *away* from the column between them — the three-rung statement is `xl` entry 30 (§10 F25). AED figures: `docs/private/f18-adjudication-2026-10-05.md`; §10 rulings 2026-10-05. |

## RETIRED 2026-10-08 (`WF-17` API pages, §9 item 331, 12:00 slot) — was 2026-10-08 — 🟡 OPEN (found by `PORT-24` step 1, §9 item 308; not fixed in passing) — the committed API reference is stale for `load_study` / `run_study_file`: three `tests/package/test_wf17_api_reference.py` tests fail on clean `main`

- **Tests:** `test_wf17_pages_are_what_the_generator_emits` (`AssertionError: first differing file: load_study.md`), `test_wf17_page_envelope_equals_registry` (`assert ['load_study', 'run_study_file'] == []`), `test_wf17_api_planted_controls` (`assert ['Study.md', 'load_study.md', 'run_study_file.md'] == ['Study.md']`). Printed: `WF-17 (G): gen1 vs gen2 diff=[]; gen1 vs committed diff=['load_study.md', 'run_study_file.md']`.
- **Verified at:** `b2174873` (clean tree), `docs/testing/logs/20261008T094800Z_PORT-24-step1-wf17-baseline.log` — 3 failed, 3 passed.
- **Cause (inferred from history, not diffed page by page):** the pages were last generated at `069dc26c` (WF-17 step 2, 2026-10-06), and `dc9d651b` (WF-1 step 3, YAML `thermal:` section) changed `src/fem_em_solver/study_yaml.py` afterwards. The generator is deterministic (gen1 = gen2), so the committed pages are behind the source; the generator is fine.
- **Resolves it:** whichever item next owns the WF-1 / WF-17 lane re-runs `scripts/docs/generate_api_reference.py` and commits the regenerated pages. `PORT-24` does not touch `study_yaml.py` (its `ring_modes` module has no API page).

**Retired 2026-10-08** (§9 item 331, 12:00 implementer slot): the same three tests red before regeneration (`20261008T171048Z_WF-17-item331-Wpre.log`, 3 failed, 9 passed); `generate_api_reference.py` changed only `load_study.md` (+1 line) and `run_study_file.md` (+5 lines) (`20261008T171101Z_WF-17-item331-gen.log`); then 12 passed (`20261008T171107Z_WF-17-item331-W1.log`) and the guide link pass reads `dead=0` (`20261008T171115Z_WF-17-item331-W2-docrefs-guide.log`). No `src/` or test edit.

---

## Batch moved 2026-10-10

### ✅ RETIRED 2026-10-09 (`OPS-70` fix, §9 item 347) — was 🔴 OPEN 2026-10-08 (`OPS-70` step 1 W7, 21:00 implementer slot, §9 item 335): `OPS-70: tests/package/test_wf16_study_sweep.py::test_wf16_sweep_api_and_refusals red on H`

> Retired by item 347 (2026-10-09 21:00 slot) under the 2026-10-09 18:00 review ruling "pin, not API": the (K′) pin now encodes `Study.__init__`'s current signature (`src/fem_em_solver/study.py:100–110`, `conductor=None` keyword-only after `density_kg_m3`) and additionally asserts `conductor`'s `None` default and `KEYWORD_ONLY` kind; no `src/` edit, no band moved. `20261010T020048Z_OPS-70-item347-W7.log:5641` — `16 passed in 130.54s`, 132 s at `-n 2`, `FEM_EM_WF16_DEGREE=1`; sweep (R′1) 6.601e-11 / 9.557e-11 / 9.449e-11 (band 1e-9), (G′) reciprocity 3.190e-15 / 4.509e-16 / 5.601e-16 at 10 / 64 / 128 MHz, every four-port node green. Negative control: the unamended pin's red on the same signature is `20261009T022354Z_OPS-70-W7-wf16-deg1-fp-sweep.log:5645`.

- **Test id:** `tests/package/test_wf16_study_sweep.py::test_wf16_sweep_api_and_refusals` ((K′) signature pin, `:160–169`).
- **Literal symptom:** `AssertionError: assert ('self', 'geometry', 'material', 'excitation', 'frequency_hz', 'degree', 'density_kg_m3', 'conductor') == ('self', 'geometry', 'material', 'excitation', 'frequency_hz', 'degree', 'density_kg_m3')` (`20261009T022354Z_OPS-70-W7-wf16-deg1-fp-sweep.log:5645`; 1 failed, 15 passed, 124 s, `-n 2`, `FEM_EM_WF16_DEGREE=1`). Every numeric node in the window is green: four-port (R) / (T2) / `PORT-11` gates, sweep (R′1) 6.601e-11, (G′) reciprocity 3.190e-15; controls (c1) 0.0491 / (c2) 0.2492 at their priors.
- **Commit verified at:** `H` = `9b50dab5`. Last green log of this node: `20261006T033609Z_WF-16-s2-W1-deg1.log` (commit `5404e634`). `git log --oneline 5404e634..9b50dab5 -- src/fem_em_solver/study.py src/fem_em_solver/cli.py`: `f7315c43`, `384b0a88`, `b9ad9e94`, `ee87e99b`, `df25001f`, `24728153`, `e763fa93`, `1f303fd6`, `e4b2f445`, `19def312`.
- **Cause (diagnosed, not fixed):** `df25001f` (`WF-16` step 7, §9 item 324) added the keyword `conductor` to `Study.__init__`; the sweep module's (K′) test pins the constructor's parameter tuple and was not updated or re-run in that slot (step 7 re-ran only the step-1 four-port window). The defect is the stale pin or the API addition — a ruling for the fix chunk, not this item. No in-slot fix (`OPS-70` ruling (6)). `OPS-70` step 1 stopped at W7; W8–W20 not run.

## ✅ RETIRED 2026-10-08 (OPS-69, §9 item 336) — was 🟡 OPEN 2026-10-08 (found by `EX-89`, §9 item 332, 12:00 slot; not fixed in passing) — the examples docrefs census reads `dead=2`, `exit=1` on clean `main`

> Retired by `OPS-69` (21:00 implementer slot): SPEC (i) reworded to name no file; `20261009T021429Z_OPS-69-W1.log:96` reads `RESULT: dead=0 guide=0 stale=56 stale_severity=report exit=2` (W0 control `20261009T021416Z_OPS-69-W0.log:98`: `dead=2 … exit=1`); `20261009T021433Z_OPS-69-W2.log:538` — `20 passed in 11.49s`, tree anchor and both planted-dead cases green.

- **Census:** `check_example_doc_references.py` (default root) — `RESULT: dead=2 guide=0 stale=56 … exit=1`, identical before and after item 332 (`20261008T172117Z_EX-89-item332-pre-docrefs.log:36–38`, `20261008T173603Z_EX-89-item332-post-docrefs.log`). Last `dead=0`: `20261008T004224Z_EX-102-post-docrefs.log`.
- **Dead references:** `ans9_aed_results.json` (`examples/ansys_benchmarks/ans9_coil_driven_gel_heating_10MHz/SPEC.md:36`) and `ans9_corrected_rows.json` (`…/SPEC.md:38`) — "no run produced it … and it is not committed under examples/".
- **Cause (from history):** `4d4d93c5` (operator AED run recorded for `ANS-9`, 2026-10-08) added those filenames to the tracked SPEC. They name gitignored private result files, so the census can never find them.
- **Consequence:** every `EX-*` item's "docrefs `exit != 1`" done-when clause is unmeetable until fixed; item 332 committed with the census unchanged by it.
- **Resolves it:** the review rewords the SPEC lines so they name no file (e.g. "the gitignored `aed_results/`") or adds the two names to the census allowlist with the private-file reason; then one docrefs window reads `dead=0`.

---

## Batch moved 2026-10-17

## 2026-10-10 — ✅ RETIRED 2026-10-14 (weekly review; `ANS-9` re-adjudicated on the degree-2 exact-rectangle column, every row AGREE, probe D within 2.5 %) — was 🟡 OPEN (`ANS-9` step 1 adjudication, weekly review; §10 chain P21) — the coil-driven gel phantom's **far-side probe D is 11 % colder in our thermal solve than in AED's** (band 10 %), while the integral rows and the near-feed probes agree within 2 %

- **Test id:** none red — a benchmark finding. Our column: `ans:9`
  (`examples/ansys_benchmarks/ans9_coil_driven_gel_heating_10MHz/`,
  `20261004T142653Z_ANS-9-step1a-W1.log`), degree-1 EM on `ANS-2`'s 10 MHz
  fixture (the pre-2026-10-03 port sheets), P1 thermal, CN Δt = 31.35 s.
- **Symptom:** against the operator's AEDT 2026 R1 run of 2026-10-06 (HFSS
  First Order + IcepakFEA transient; corrected rows, see the SPEC's status
  block) rows 1 / 2 / 3 / 6 agree within 0.8 % / 0.5 % / 1.4 % (bands 5 %),
  probes A–C within 2 % and the B > D > C ordering holds, but **probe D
  (−0.020, 0, 0 — facing the leg opposite the driven port) reads 11 % below
  AED's**, outside the pre-registered 10 %. AED values stay private
  (`docs/private/t13-t29-ans9-adjudication-2026-10-10.md`).
- **Reading (not a diagnosis):** diffusion is the same physics in both codes
  and every integral quantity agrees, so the gap is in the EM source on the
  far side of the coil — the field set by current induced in the three
  50 Ω-terminated legs. Candidates, in order: the EM order mismatch (AED
  First Order against our degree 1 — the self-class signature `ANS-4` showed
  at degree 1); the ragged pre-10-03 port sheets of `ANS-2`'s fixture (T6's
  convention is exact rectangles at degree 2); the termination model.
- **Two AED-side departures stated with the verdict:** one load step per
  report time with a fixed 13 automatic substeps (≈ 125 s steps on the
  1254 s run; AED's own energy balance then over-counts the surface heat —
  the Simpson figure on the field saves is the adjudicated row 6), and the
  operator's offline volume / area correction of the harvested means.
- **Next:** `ANS-9` step 1b (§10 P21) prints our local SAR at AED's nine
  monitor points on today's column and on the degree-2 exact-rectangle
  fixture, re-runs the thermal half on the latter and re-reads rows 2–4
  under a pre-registered rule; the 10-14 weekly re-adjudicates D. Nothing
  about the thermal hand-off is in question.
- **Measured (step 1b (a), §9 item 371a, 2026-10-11, our values only):**
  local SAR at A–E per 1 W in the phantom, on `ans:9`'s own degree-1 field
  (d1), item 192's exact-rectangle mesh at degree 1 (d1x) and at degree 2
  (d2x). (D/C) = 205.6 / 10.2 / 435.6, so **ρ = (D/C)_d2x / (D/C)_d1 =
  2.119 ⇒ "no branch"** under the weekly's pre-registered rule. C
  (0, 0.020, 0) sits near a field null (0.0129 / 0.303 / 0.0070 W/kg), so
  the ratio is mostly C's resolution. SAR_D itself is 2.651 / 3.094 /
  3.064 W/kg: +15.6 % from d1 to d2x, of which +16.7 % is the sheets and
  phantom mesh step (d1 → d1x) and −1.0 % the order step (d1x → d2x).
  (`20261011T184301Z_ANS-9-step1b-W1-d1.log`,
  `20261011T184526Z_ANS-9-step1b-W2-d2x.log`.) The entry stays open for
  the 10-14 weekly.
- **Measured (step 1b (b), §9 item 371b, 2026-10-11, our values only):**
  `ans:9`'s thermal half on item 192's coarse (537-cell) phantom with the
  d2x and d1x sources, all balances green ((R1′) ≤ 4.4e-16, (E) ≤ 3.9e-12,
  (R3′) 4.9e-6 / 5.2e-6). θ₄₀ at D is 0.7233 K (d2x) vs 0.6202 K (d1x):
  **+16.6 % from the order step alone on the same thermal mesh**, while
  ⟨θ⟩₄₀ moves +0.28 % and A / B / C / E ≤ 1.9 %. A coarse-phantom reading
  (≈ 15 mm cells), not a re-record of `ANS-9`'s column
  (`20261011T185145Z_ANS-9-step1b-W2b-d2x-thermal.log`). The entry stays
  open for the 10-14 weekly.
- **Measured (step 1c, §9 item 396, 2026-10-12, our values only):** route B
  (`GEO-27` step 2's builder) reproduces item 192 at x1 (S to 7.3e-11 / 6.6e-11,
  point SAR to ≤ 2e-9, θ₄₀ to ≤ 6.5e-7). Per rung x1 / h005 / h0025 (tag-3
  537 / 8 498 / 58 960): **θ_D order effect (d2x / d1x) 1.1663 / 1.0964 /
  1.0794**; point-SAR_D order effect 0.9903 / 1.1111 / 1.1929; ⟨θ⟩₄₀ order
  effect 1.0028 / 1.0005 / 1.0001. h_p moves against P21 (x1), d2x / d1x:
  SAR_D −2.62 % / −13.2 % (h005), −2.75 % / −19.3 % (h0025); θ_D +4.20 % /
  +10.8 %, +4.64 % / +13.1 %; D/C +19.2 % / ×270, +26.8 % / ×8.9; D/A +1.50 % /
  +3.00 %, +2.34 % / +20.4 %. κ_s rose with refinement (d2x 1.34 / 1.54 / 1.61),
  against the predicted fall. Readings for the weekly; no branch read
  (`20261012T100019Z_ANS-9-step1c-W1-x1.log`, `20261012T100223Z_ANS-9-step1c-W2-h005.log`,
  `20261012T100436Z_ANS-9-step1c-W3-h0025.log`). The entry stays open.
- **Ruled 2026-10-12 (weekly; numbers in
  `docs/private/f27-f26-h15c-ans9-adjudication-2026-10-12.md`):** the gap is
  discretisation on our side, not the far-leg termination model — at degree 1
  θ_D rises 13 % between the 537-cell phantom and h_p = 0.0025 m, and the order
  adds a further 8 % on the refined phantom; **the degree-2, h_p = 0.0025 m
  candidate column sits within 2.5 % of AED at probe D (band 10 %)**, while
  the degree-1 refined column is still outside. P21's D/C instrument was
  dominated by the near-null probe C and decided nothing (§10 epitaph).
  **Next:** §10 P24 re-records `ans:9`'s column at degree 2 on exact
  rectangles with the h_p = 0.0025 m phantom (route B, every SPEC row); the
  2026-10-14 weekly re-adjudicates all rows against the operator's existing
  AED column (no new AED run). The entry retires on that adjudication.
- **Measured (step 1d, §9 item 410, 2026-10-13, our values only):** `ans:9`'s
  primary column is now degree 2 on exact rectangles, h_p = 0.0025 m (route B,
  201 072 / 58 960 cells), every identity green. Rows 2–4, degree 2 against
  the carried degree-1 ragged column: ⟨θ⟩ at 313.5 / 627 / 940.5 / 1254 s
  0.318123 / 0.614370 / 0.892184 / 1.153415 K vs 0.318096 / 0.614315 /
  0.892104 / 1.153314 (×1.00009); row 3 slope 1.058563e-3 K/s both (×1.000000);
  θ₄₀ A 0.889137 / B 3.514465 / C 0.683153 / **D 0.756873** / E 0.959953 K vs
  0.905798 / 3.564854 / 0.685305 / **0.697444** / 0.958639 (×0.9816 / 0.9859 /
  0.9969 / **1.0852** / 1.0014); B > D > C holds
  (`20261013T093826Z_ANS-9-step1d-W2-h0025.log`). The entry stays open for the
  10-14 weekly's adjudication.
- **Retired 2026-10-14 (weekly; numbers in the gitignored
  `docs/private/ans9-readjudication-2026-10-14.md`):** on step 1d's column
  every row AGREES with the operator's existing AED column — probe D within
  2.5 % (band 10 %), probes A–C within 2.3 %, B > D > C held, ⟨θ⟩ within
  0.8 %, slope within 0.5 %, surface heat within 1.3 % (bands 5 %). The
  finding was ours: 13 % phantom resolution plus 8 % element order on a
  column recorded below the production convention; no AED-side cause, no
  port-model cause. **Lesson (§10 epitaph 2026-10-12, now binding):** a
  benchmark column is adjudicated at the production convention only — the
  degree-1 near-feed probes A / B "agreed" within 1.4 % and moved 1.9 /
  1.4 pp under the order move. `ANS-9` ✅.

## 2026-10-17 — ✅ RETIRED 2026-10-17 (opened and retired in the same commit, `WF-13` step 9, §9 item 507, 04:30 CDT implementer slot) — `femem run --estimate` on a thermal-only ball YAML built the ball mesh, then raised an uncaught `AttributeError` on every rank (status 1)

- **Symptom:** `mpiexec -n 2 python3 -m fem_em_solver run --estimate tests/package/wf1/thermal_ball_duty_a.yaml` (complex build, `f3ca6d0a`, unmodified `src/`) printed, on both ranks, `AttributeError: 'ThermalBallGeometry' object has no attribute 'sphere_tag'` at `study.py:1687` (`phantom_tag=self.geometry.sphere_tag`) in `_study_estimate`. Status 1, 11 s (`20261017T093035Z_WF-13-s9-W0-observe.log:35–69`).
- **Cause:** `_study_estimate` had no refusal for a thermal-only Study. Its non-four-port branch built the ball mesh and then read `sphere_tag`, which `ThermalBallGeometry` does not have. `cli._estimate` catches `(ValueError, TypeError)` only, so the `AttributeError` escaped.
- **Fix (`study.py` only, +8 lines; `cli.py` untouched):** a pre-mesh `_refuse_thermal_only(study, method)` beside `_refuse_f_human` raises `TypeError("estimate is not defined on a thermal-only Study (ThermalBallGeometry has no EM solve; use temperature())")`. It is called directly after `_refuse_f_human(self, "estimate")`. Item 480's handler then gives rc 2.
- **Verification:** W1 `test_wf13_cli_estimate.py` 8 passed, 31 s (`20261017T093117Z_WF-13-s9-W1.log`). The new node `test_wf13_cli_estimate_thermal_ball_refusal` gave rc `[2, 2]`, rank 0's one `femem: error:` line equal to the API message, no output dir, and `ThermalBallGeometry.build` calls `[0, 0]`. Control (helper patched to a no-op): `AttributeError` matching `sphere_tag` with build calls `[1, 1]`. RG1 F-human refusal 5 passed (`20261017T093204Z_WF-13-s9-RG1.log`). The real-build thermal-ball re-gates are at item 500's digits: RG2 3 passed / 3 skipped (`20261017T093213Z_WF-13-s9-RG2-thermal-duty.log`), RG3 2 passed / 3 skipped (`20261017T093406Z_WF-13-s9-RG3-thermal-wall.log`).

---

## Batch moved 2026-10-19

## 2026-09-30 — ✅ RETIRED 2026-10-17 (weekly review, T41: the `f₀` clause was degree-1 discretisation — T13, 2026-10-10, within 0.11 % / 0.03 % at degree 2 on exact rectangles — and the Q / `R_in` clause moves to the 2026-10-17 entry below at its h-limit) — was 🟡 OPEN (`ANS-7` v2 adjudication, weekly review; §10 chain T10b) — the tuned F-small resonance sits **1.2–1.4 % above the order-matched AED columns and 1.9 % above AED's First Order** in `f₀` (capacitor equivalent 2.5–3.9 %), with Q 6–10 % and `R_in` 7–11 % low; mode identity and the copper contribution AGREE

| | |
| --- | --- |
| **Symptom** | None in CI — `ANS-7`'s Route A = Route B band, the single − → + zero and the k = ±1 share (0.999994 / 0.999998) are green (`20260929T100040Z_ANS-7-item68-W3-A.log`, `20260929T094009Z_ANS-7-item69-W2-B2.log:2819–2830`). The gap is against the external code only. |
| **What agrees** | Both codes find one k = ±1 parallel resonance in the 62–67 MHz window; each code's two routes agree with each other; the conductor's part of 1/Q agrees to 5 % (ours higher), `Q_Cu/Q_PEC` to 1.1–2.4 %, the Cu − PEC shift of `f₀` to 5 % of AED First Order. |
| **What does not** | `f₀` (ours high, the same on the Cu and PEC columns — so not the conductor model), 14–19 linewidths at Q ≈ 1000. Q and `R_in` (ours low); the split puts most of it in the phantom's part of 1/Q, 8–11 % high. AED's own Zero → First Order move on `f₀` is 0.5 %, away from ours. |
| **What is and is not known** | Our side is one rung: degree 1 on the 80 181-cell hole mesh, port width uncorrected. No order or h move of `f₀` / Q exists on our side, and our degree 1 → 2 move on this mesh family's S is 4–7 % at 64 MHz. The `f₀` gap is the size and sign of the open `S₁₁` self-class gap (entry 2026-09-19 above: our self term low); whether it is the same term is unmeasured. |
| **Consequence** | **No tuned-capacitor, `f₀`, Q or `Z_in` figure of ours is quoted as agreeing with AED.** "Tuned" stays an internal statement (circuit = field, eigen = circuit) on one fixture. `ANS-7` stays 🟡. Verdict: AGREE on the mode and the copper contribution; absolute rows INCONCLUSIVE at degree 1 — not DISAGREE (one rung; §10 epitaph 2026-09-19), not AGREE. |
| **Next** | §10 chain T11 (the sparse-grid Route B, then its degree-2 order move), T11b (the κ-corrected port width), T12 (`xl` entry 24, the h = 0.005 hole rung), T16 (the phantom's 1/Q under `phantom_resolution`); the weekly re-reads on those (T13). Numbers: gitignored `docs/private/ans7-ans6-adjudication-2026-09-30.md`. |
| **Re-read 2026-10-03 (weekly, T13; transcribed by the 2026-10-04 03:00 daily)** | **Still OPEN — retires at T13 on the h-rung, not before.** At degree 2 on exact rectangles, κ off (item 196), `f₀` sits **0.11 % below** AED First Order on both columns (1.1 / 1.9 linewidths, Cu / PEC); Q within 1.9 % / 0.3 %, `R_in` within 1.6 % / 0.1 %; our order move on `f₀` is −1.985 %, the exact-vs-ragged sheet move 0.005 %. **The `f₀` gap was degree-1 discretisation**, not the sheet; the phantom's part of 1/Q is within 0.3 % at degree 2. The conductor's part stays 4.4 % high at both orders (inside the 5 % already quoted). One mesh: the h-rung is T12′ (`xl` entry 28, or entry 24) with T19's intermediate rung. Numbers: gitignored `docs/private/f9-t13-ans8-adjudication-2026-10-03.md`. |

## 2026-10-13 — ✅ RETIRED 2026-10-17 (weekly review: `IMP-6`'s step-2 family is frozen on the row after I29 — the rim truncation below ρ_r/3 stands as a documented builder limit, worked around by `oblate_rim="edge"`; no gate reads the truncating path) — was 🟡 OPEN (`IMP-6` step 2, §9 item 408, 22:30 implementer slot; ruled by the 03:00 daily review) — the private oblate-spheroid builder **truncates the rim when the tip Balls are refined below ρ_r/3**: air mask-cell midpoints sit inside the true spheroid

- **Test id:** the parked `tests/validation/test_imp6_oblate_rim_ladder.py` (`attempt/IMP-6-item408-20261013T0339Z`). The rim-ball Φ reads −13.78 % on rung ρ_r/4 (`20261013T033339Z_IMP-6-step2-W1b.log`) and −20.07 % on ρ_r/5 (`…033502Z_…-W1c.log`), against `MAT4_BAND` 10 %. β and the cavity check stay green on both rungs. Rung ρ_r/3 reproduces step 1 to every digit (`…033249Z_…-W1a.log`).
- **What is measured** (`20261013T033815Z_IMP-6-step2-W0diag.log`): the number of air midpoints inside the body is 0 / 31 / 93 on the three rungs. The deepest sits at a − 0.132 / a − 0.163 mm, past the focal ring at a − 0.121 mm. Step 1's max-edge h_rim/ρ_r jumps from 0.665 to 2.720 to 2.987 while the median falls. The FEM therefore sees a blunted edge.
- **Not affected:** step 1's and `IMP-2`'s default paths (0 midpoints inside the body), and `EX-111`.
- **Resolves with:** §9 item 426, a rule (c) mesher fix on the private builder: the rim constrained explicitly first, the Netgen pass dropped only if that is not enough. If no licensed variant clears ρ_r/4, the measured negative is banked on the `IMP-6` row instead.
- **Item 426 (2026-10-13, 15:00 CDT slot): the variant that avoids it is `spheroid_in_box_domain(oblate=True, oblate_rim="edge")`** — the body built by revolving a half-ellipse, so the equator is a model edge shared by the two cavity faces (2 faces, printed). Inside midpoints 0 / 0 / 0 on rungs ρ_r/3, /4, /5, min SICN 0.3605 / 0.3994 / 0.3466, cells 33 334 / 53 959 / 89 530 (`20261013T200400Z_IMP-6-step2-a2-W0-edge.log:168, :286, :417`); air OCC mass rel 0.000e+00 (exact CAD). The **default path is unchanged and still truncates** — (C-fid) re-reads 31 / 93 on ρ_r/4 / ρ_r/5 (`20261013T200330Z_IMP-6-step2-a2-W0-ctl.log:172, :312`), and step 1, step 3 and `EX-111` keep their meshes (re-gate `20261013T200442Z_IMP-6-step2-a2-W0r.log`: 39 142 / 32 676 cells, SICN 0.2667 / 0.3707). The entry therefore stays open as a documented property of the default builder: **any new oblate fixture finer than `resolution_tip` = ρ_r/3 uses `oblate_rim="edge"`**. Variants (1b) `tube` and (2) `balls-nonetgen` are coded but were not run (the brief's order stops at the first passing variant).

### 🟡 OPEN 2026-10-19 (`POST-8` step 3b-ii-e1, 09:00 implementer slot, §9 item 547): `th:4` (`04_evanescent_waveguide_decay.py`) is **red at `-n 4` on `main` as it stands** — its transverse-profile assertion reads 2.25 % RMS against `PROFILE_RMS_MAX`, where `-n 2` prints 0.124 % *(heading corrected 2026-10-19 10:30 review: the slot's 0.200 % is the August `-n 2` record)*

| | |
| --- | --- |
| **Symptom** | `mpiexec -n 4` of the bare, unedited example: `AssertionError: the exported transverse profile is 2.25% RMS from sin(pi x/a) after peak normalisation`, Status 1, 5 s (`20261019T140617Z_POST-8-step3b-ii-e1-th4-baseline-n4.log`). The edited example under the spy reads the same 2.25 % (`20261019T140524Z_POST-8-step3b-ii-e1-th-n4.log`); the other four th examples and `mag:1` / `mag:5` ran green there. |
| **Verified at** | `main` at `f9b8dce4` plus the item's keyword edits, by execution; the baseline had the `th:4` edit stashed. `-n 2` is green with (O) owner indices exact. |
| **Cause** | Not diagnosed. The profile samples the CG1 export `e_re` / `e_imag` at z = L/2, a lattice plane of the box mesh. At `-n 2` every one of the 25 probes is multi-incidence with owner = legacy and delta 0 (so the rule does not move it); the width dependence is the partition dependence of the CG1 interpolation of the DG split (step 1's mechanism, which the owner rule cannot reach), not the owner rule. **Corrected by the 2026-10-19 10:30 review (log-pathologist ruling):** there is no DG split here. The export is `e_cg.interpolate(fields.e_complex)` (`:250`), a direct N1curl → CG1 vertex interpolation. Owner and legacy agree to ≤ 1.016e-15 at the 25 profile probes, not bitwise at all of them (`…-th4-n2.log:91` on the branch), and the conclusion stands. The N1curl readings are width-free: γ 37.650399, L2 4.406648e-02 at both widths. `CG1_VS_NEDELEC_MAX` moves too, 0.027 % at `-n 2` → 0.490 % at `-n 4` against its 0.5 % bound. The mechanism (last-writer vertex values depend on the partition) is plausible but not established; a CG1 vertex dump at two widths would establish it. Today's `-n 2` reads (0.124 % / 0.027 %) also differ from August's (0.200 % / 0.117 %, `20260810T020325Z_EX-7-run1.log.gz:53`), cause unidentified. |
| **Scope** | `th:4` cannot be run as a `-n 4` example; item 547's `-n 4` / `-n 1` ladder for it fails on the example's own band. No band touched, nothing re-recorded. |
| **Resolves with** | The review re-poses 547: either `th:4` at record width `-n 2` only (as item 548's examples run), or the profile band is ruled width-dependent. The executor sets neither. **10:30 review:** item 547 is re-posed as item 549, with `th:4` at `-n 2` only. That holds the item but does not resolve this entry. The entry resolves when the weekly rules on the content change: sample the profile and the cross-check from the N1curl field under the owner rule, or state the CG1 bands' width. **Ruled 2026-10-19 (weekly, §10): the content change is licensed** — `POST-8` step 3b-ii-e3 re-samples `th:4`'s transverse profile and the CG1 cross-check from the N1curl field under the owner rule, re-recorded under the `(1*)` licence with `PROFILE_RMS_MAX` and `CG1_VS_NEDELEC_MAX` unmoved; the CG1 figures (0.124 % / 0.027 % at `-n 2`, August's 0.200 % / 0.117 %) stay in the guide as a labelled `-n 2` record. **This entry retires when `th:4` is green at `-n 1 / 2 / 4` on the re-sampled profile.** The CG1 export's width dependence itself is step 1's mechanism (entry 2026-09-27) and is not re-opened here. |

---

## Batch moved 2026-10-24

### ✅ RETIRED 2026-10-22 (`POST-8` 3b-iii, §9 item 587; retirement row below) — was 🟡 OPEN 2026-09-27 (03:00 daily review, from the 2026-09-26 09:00 slot's disclosure; `log-pathologist` ruling) — `ans:2`'s **pointwise** SAR / `|E|` sample (Row 1 of its comparison) **depends on the MPI width**: −19 % … +112 % in SAR between `-n 2` and `-n 4` on the same code and mesh, while every reduced quantity agrees to ≈ 12 figures

| | |
|---|---|
| **Symptom** | The four named points (±x, ±y at z = 0) × four drives. `-n 4` gives P1 row `5.878191741e-07 …` in every `-n 4` run 09-09 … 09-25 (`20260909T140558Z_ANS-2-step1.log.gz:1902`, `20260925T110348Z_OPS-63.log:1909–1912`; tracked `metrics.json:49–94`, `mpi_ranks: 4`). `-n 2` gives `4.858355042e-07 …` in all three 09-26 runs, including one before the setup-figure edit (`20260926T140325Z_EX-57.log:1902–1905`, `20260926T140900Z_EX-57.log:1909–1912`, `20260926T141326Z_EX-57.log:1902–1905`). The SAR change is −18.8 % … +4.4 % at the driven points and −42 % … +112 % elsewhere. Inside one run the four C4-symmetric driven points already scatter 15.73 % at `-n 4` (`OPS-63.log:1944`) and 13.3 % at `-n 2`. **Unchanged across widths:** whole-phantom power (5e-12 relative, `metrics.json:100` vs `140900Z:1914`), accepted power, the 10 g / 1 g C4 pairs, the `PORT-9` σ_max and class spreads, and the cell count (199 920). |
| **Verified at** | **`0b00091` (`POST-8` step 1, 2026-09-27): same-commit pair measured.** `-n 2` reproduces its logged P1 row to 8.94e-11 relative and `-n 4` its own to 4.95e-10 (≤ 1e-9); P1 phantom power within 4.80e-12 / 6.62e-14 of `metrics.json:100` at `-n 2` / `-n 4` (`20260927T112115Z_POST-8-step1-n2.log:1902–1913`, `20260927T112501Z_POST-8-step1-n4.log:1904–1915`). The sample is reproducible at a fixed width and moves between widths on one commit. Earlier reading: `500a575` (`-n 2`) against the `-n 4` runs through `d6289ac`, no same-commit pair. A diff review excludes a code change: `d6289ac` touches no pointwise-path line, and `evaluation.py` last changed on 09-21, before the 09-25 `-n 4` run that still reproduces the tracked row. **`POST-8` step 2 (2026-10-02, §9 item 177):** the owner rule (`owner="lowest_original_cell"`, opt-in) returns `ans:2`'s 16 raw N1curl SARs within 1.37e-10 / 6.55e-11 / 3.03e-10 of step 1's per-cell values at `-n 2 / 4 / 1` (cross-width pairwise max 2.98e-10); the default CG1 route still reproduces its own width's logged row (8.89e-11 / 5.89e-10) and still moves ≥ 12.9 % across widths (`20261002T200331Z_POST-8-step2-W1.log`, `20261002T200654Z_POST-8-step2-W2.log`, `20261002T201013Z_POST-8-step2-W3.log`). Still OPEN: no caller adopts the rule until step 3. **`POST-8` step 3a (2026-10-03, §9 item 188):** the same 16 entries as raw owner-rule, owner-cell average and 1 g mass average agree with the `-n 2` npz at `-n 4` / `-n 1` to ≤ 2.21e-10 / 2.15e-10 / 1.07e-10 relative (1e-9 asserted); the CG1 route still separates by ≥ 12.9 % (`20261003T033225Z_POST-8-step3a-W1.log`, `…033703Z_…W2.log`, `…034100Z_…W3.log`). Still OPEN until 3b adopts the rule and re-issues Row 1. **`POST-8` step 3b-i (2026-10-18, §9 item 526):** `post.sar.point_sar` now defaults to the owner rule (`owner=None` = legacy) and `phantom_fields.py:163` is asserted owner-compliant; on step 2's cube the DG1-split N1curl SAR is width-invariant to 4.087e-16 at `-n 1 / 2 / 4`, while the legacy path moves the vertex SAR by 11.9 % / 12.2 % at `-n 4` / `-n 1` (`20261018T125312Z_POST-8-step3b-W1.log`, `…125322Z_…W2.log`, `…125326Z_…W3.log`). `ans:2`'s CG1 route is untouched; still OPEN until 3b-iii re-issues Row 1. **`POST-8` step 3b-ii-t1 (2026-10-18, §9 item 536):** the three time-harmonic test modules sample under the owner rule; on `TH-12`'s coarse sphere the 25 probes are single-incidence at degree 1 and 2 (rule bitwise = legacy there), owner indices are exact and the sampled E width-invariant to 1.108e-11 at `-n 4 / 1` (1e-9 asserted) (`20261018T213254Z_POST-8-step3b-ii-t1-W1.log`, `…213419Z_…W2.log`, `…213438Z_…W3.log`). Still OPEN until 3b-iii. **`POST-8` step 3b-ii-t2 (2026-10-19, §9 item 537):** the coil-phantom DG0 centerline / symmetry samples and the straight-wire radial sampler use the owner rule; all 91 probes are single-incidence (rule bitwise = legacy), coil-phantom owner indices are exact at `-n 1 / 2 / 4`, and both modules re-run green at their bands (`20261019T021325Z_POST-8-step3b-ii-t2-W1.log`, `…021402Z_…W2.log`, `…021435Z_…W3.log`, `…021728Z_…W4-wire.log`). **`POST-8` step 3b-ii-e1 (2026-10-19, §9 item 549):** seven example scripts sample under the owner rule (`mag:1` on the CG1 `B_lag`); each runs as `__main__` under a spy at `-n 2` with its own assertions green, 651 of 706 probes single-incidence (rule bytewise = legacy there), and owner indices exact at `-n 4` / `-n 1` for all but `th:4`, which runs at `-n 2` only (the entry of 2026-10-19 above) (`20261019T170143Z_POST-8-step3b-ii-e1b-mag1-n2.log` … `20261019T170341Z_POST-8-step3b-ii-e1b-th9-n2.log`, `…170352Z_…mag-n4.log`, `…170410Z_…th-n4.log`, `…170448Z_…mag-n1.log`, `…170530Z_…th-n1.log`). **`POST-8` step 3b-ii-e2a (2026-10-20, §9 item 551):** `mag:2` and `mag:6` (CG1, keyword only) sample under the owner rule on ghost-free meshes (`cells_ghost` 0 on both ranks); each runs as `__main__` under the spy at `-n 2` with its own assertions green, 25/25 and 10/10 probes single-incidence, 0 byte mismatches against the legacy owner (`20261020T003156Z_POST-8-step3b-ii-e2a-mag2-n2.log`, `20261020T003421Z_POST-8-step3b-ii-e2a-mag6-n2.log`). **`POST-8` step 3b-ii-e2d (2026-10-20, §9 item 568 Part A):** `mag:2`, `mag:4`, `mag:6`, `ports:12` (keyword edits; `mag:4` and `ports:12` on ghosted meshes) sample under the owner rule at record width, each example's assertions green, (S-sb) 0 at the single-incidence rows on every call, cross-batch (S) printed (`mag:4` 10 mismatches ≤ 2 ULP, `ports:12` 6 ≤ 32 ULP); e1 re-gates green on the `_ncomp` fix (`20261020T213252Z_POST-8-step3b-ii-e2d-mag2-n2.log`, `…213520Z_…mag6-n2.log`, `…213751Z_…mag4-n2.log`, `…213924Z_…ports12-n4.log`, `…214139Z_…e1-real-n2.log`, `…214153Z_…e1-complex-n2.log`). Still OPEN until 3b-iii. **`POST-8` step 3b-ii-e3b (2026-10-21, §9 item 580):** `th:4` (DG split + N1curl) and `mag:1` (DG `B`) sample under the owner rule, `mri:1` keyword; (S-sb) 0 at every width, `th:4` green at `-n 1 / 2 / 4`, (O) exact (`20261021T110150Z_POST-8-step3b-ii-e3b-th4-n2.log`, `…110229Z_…th4-n4.log`, `…110303Z_…th4-n1.log`, `…110135Z_…mag1-n2.log`, `…110214Z_…mag1-n4.log`, `…110242Z_…mag1-n1.log`, `…110204Z_…mri1-n2.log`). Still OPEN until 3b-iii. **`POST-8` step 3b-ii-e4 (2026-10-22, §9 item 586):** `ports:1` samples under the owner rule with its four records re-based `(1*)` at `-n 2`, bands unmoved; green at `-n 1 / 2 / 4`, (S-sb) 0, (O) exact (`20261022T110209Z_POST-8-step3b-ii-e4-ports1-n2.log`, `…110501Z_…n4.log`, `…110738Z_…n1.log`). Still OPEN until 3b-iii. |
| **Cause** | **Localized by `POST-8` step 1 (measured, 2026-09-27); the weekly rules on it.** Each named point lies inside one owned cell at both widths. The raw N1curl `E` evaluated in that cell is width-invariant: at all 16 (point, drive) entries its SAR agrees between `-n 2` and `-n 4` to ≤ 1.0e-9 relative, the resolution of the 10-figure print (e.g. `-n 2` / `-n 4` point +x, P1: 5.514140460e-07 W/kg at both, same `E` digits, global cell index 85735 vs 68816 after renumbering). The CG1 sample does not equal it at either width (`-n 2` 4.858355042e-07, `-n 4` 5.878191741e-07 W/kg). In that cell the CG1 nodal `|E|²` differs from the cell's own N1curl value at the vertex by ratio 0.83–0.98 (`-n 2` log:1925–1929, P1). So the width dependence enters through the N1curl → CG1 nodal interpolation, not the solve or a reduction. The "bracket" readout the item pre-registered does not apply: each point has one incident cell, so there is no per-cell range to bracket with. **Original hypothesis (code reading, not log):** the example samples a CG1 interpolant of the degree-1 N1curl field (`02_birdcage_coil_driven_sar_10MHz.py:327–333`). The normal component of N1curl E jumps across every facet, so a CG1 vertex takes the value of whichever incident cell writes last, and then the owning rank's value after `scatter_forward`. Both depend on the partition. `evaluate_vector_field_parallel` itself is collective and rank-safe (`post/evaluation.py:109–135`), so this is **not a missing reduction**. |
| **Scope** | `ans:2`'s Row 1 pointwise figures, and the **public relative level of the 09-26 weekly's `ANS-2` step-4 pointwise statement** (the "within 6.3 %" was read at one width). Integrals, 10 g / 1 g averages, S / Z and every gate are unaffected. Other examples interpolate N1curl → CG1 too (`th:1`–`th:9`, `ports:1`–`5`), mostly for ParaView export only; whether any of them *samples points* from it is `POST-8` step 1's census. The CG1 ParaView bundles probably carry the same partition dependence (an inference, not checked). |
| **Retired** | **2026-10-22 (`POST-8` 3b-iii, §9 item 587).** `ans:2`'s tracked Row 1 is re-issued as the owner-cell average (3a's route (ii); `metrics.json["pointwise_rule"]` = `"owner_cell"`), asserted against 3a's `-n 2` record at 3a's `WIDTH_RTOL` = 1e-9 at every width: max relative miss 8.329e-11 / 1.838e-10 / 1.890e-10 at `-n 2 / 1 / 4`, owner indices exact `[150135, 162464, 142550, 169932]`, (Im) 0, the legacy CG1 row kept as a labelled record and asserted to separate from the owner-cell row by more than 1e-8 (1.190e-02 / 1.413e-02 / 2.022e-03). Windows `20261022T123323Z_POST-8-step3b-iii-ans2-n2.log` (251 s), `…123751Z_…ans2-n1.log` (386 s), `…124425Z_…ans2-n4.log` (195 s); consumer `…124753Z_…ans9-consumer.log` green (7 passed, 2 skipped, 130 s). The CG1 route itself is unchanged and still width-dependent; it is a labelled record, never published. |
| **Resolves with** | `POST-8` (opened 2026-09-27 03:00 review). Step 1: a same-code control at both widths, plus the per-incident-cell N1curl values at each named point. Then a partition-independent sampling rule (the weekly's call). This entry retires with the step that lands that rule, and `ans:2`'s Row 1 is re-issued under it. |

---

## Batch moved 2026-10-26

### ✅ RETIRED 2026-10-24 as a record (weekly review: `TH-15` step 2 was killed 2026-10-07 and `TH-15` closed ✅ 2026-10-19, so "stays open until step 2 closes" can never fire; the module never reached `main`; the point-sampled `_path_voltage` reading it diagnosed is `POST-8`'s owner-rule subject, ✅ 2026-10-22, `ports:1` re-recorded; nothing on `main` is red from it — the text below stays verbatim as the two-torus record) — was 🟡 OPEN 2026-09-07 (`TH-15` step 2 proper, 16:30 implementer slot) — the PEC-hole two-torus `Z` is **exactly loss-free** (`Re Z ≡ 0`) but **2.07% asymmetric**, which lands `‖SᴴS − I‖_F` at **7.538037e-03** against 1e-9 and the mutual at **0.828497 (−17.15%)** against the imported 10%

| | |
| --- | --- |
| **Test** | `tests/validation/test_two_torus_pec_hole_ports.py::test_pec_hole_network_is_reciprocal_and_unitary` and `::test_pec_hole_mutual_matches_the_closed_form`. **Not on `main`** — the module lives only on `attempt/TH-15-step2proper-20260907T213739Z` (`144feff`). Listed here because the readings are a finding a review must dispose of before `TH-15` step 2 can close. |
| **Verified at** | `aa5ffa0` + the parked worktree. `20260907T213308Z_TH-15.log` — `-n 4`, complex build, `tests/environment` first, `-s -v --tb=short`, `timeout -k 30 500`: **2 failed / 14 passed in 215.30 s**, `Status: 1`, elapsed **217 s** (`:1904, :1907–1908`). Collect-only smoke first, 4 s (`20260907T213256Z_TH-15.log`, 5 tests). |
| **Symptom** | `Z` on the `as_hole=True` mesh with `pec_facet_tags=(1, 301)` and σ = 0 everywhere is `[[6.51916156j, 1.05007456j], [1.02878956j, 6.42566429j]]` Ω (`:624–626`). `max_ij |Re Z_ij|/|Z_ij|` = **0.000000e+00** — the lossless identity is exact, asserted green at `LOSSLESS_BAND` 1e-9 (`:627`). But `Z₁₂ ≠ Z₂₁`: 1.05007456 vs 1.02878956, **2.07%**. `‖SᴴS − I‖_F` = **7.538037e-03**, `σ_max` = **1.002865051123** (`:633–634`); `Im Z₂₁/(ωM₁₂)` = **0.828497**, −17.15% against the imported `MUTUAL_TOLERANCE` 10% (`:636`), where the solid gap route reads 0.909618 and the conduction route 0.939822. |
| **What is green in the same window** | Cells **161 461 / 161 461 = 1.000000** at the imported `CELL_COUNT_BAND`, 24.22 s to mesh (`:590, :616`); `‖S − Sᵀ‖/‖S‖` = **4.286714e-04** inside the imported `S_SYMMETRY_BAND` 1e-3 (`:632`); the σ = 800 S/m solid control through the identical sweep, `Re Z₁₁` = **+3.771673e+00 Ω > 0** asserted with `|Re Z₁₁|/|Z₁₁|` = **0.469192** against the *predicted* order 0.5 (`:1572`); the printed open-circuit reading `\|I_disp,undriven\|/I_drive` = **1.463859e-06** on both drives (`:606–615`). |
| **Cause — the asymmetry, one mechanism for both misses; not dissipation** | `Re Z` is identically zero, so no loss term enters the assembled system and the gauge-penalty / ABC hypothesis the item pre-registered is **excluded by measurement**. For purely imaginary `Z = jX` with real `Z₀`, `S = (jX − Z₀)(jX + Z₀)⁻¹` is unitary exactly iff `X = Xᵀ`; 2.07% of asymmetry lands 7.5e-3 of non-unitarity, which is the reading to the order. `‖S − Sᵀ‖/‖S‖` passes on the *same* matrix only because the S off-diagonals (≈0.02) are diluted by the diagonal (≈0.97) in `‖S‖`, so S-reciprocity is ≈50× weaker than Z-reciprocity on this network — worth knowing wherever `S_SYMMETRY_BAND` is used as a reciprocity gate on a near-totally-reflecting port pair. The asymmetry is present upstream of the reduction: `V_undriven` = 1.0224 V under the P1 drive vs 1.0437 V under the P2 drive, 2.08% apart on a nominally mirror-symmetric fixture (`:593, :595`), while the driven currents agree to 1.5e-4 (0.9937978 / 0.9939479). **Not diagnosed further**: whether the asymmetry is the distributed mesh's or the displacement route's two gap tags was not measured in this slot. |
| **Consequence** | `TH-15` step 2 does not close; the row stays 🟡. Nothing on `main` is red from this — the module is not on `main`. |
| **Resolves with** | A review ruling on (a) the source of the 2.07% — re-run at `-n 2` for a width control, and print the *solid*'s own `Z₁₂/Z₂₁` on this same route (step 2d printed S-reciprocity 5.2613e-04 but never Z-symmetry), which separates "the fixture/route is 2% asymmetric everywhere" from "the hole is"; and (b) whether the −17.15% mutual is the PEC cavity wall genuinely excluding flux the *filament* `ωM₁₂` comparand counts — in which case the comparand, not the band, is what must change (rule (f)). **Neither `LOSSLESS_BAND` nor `MUTUAL_TOLERANCE` may be widened**, and (i)/(ii) are already green. |
| **Ruled 2026-09-07 18:00 review** | (a) The asymmetry is read as a per-port *current-reading calibration* `c₁/c₂` — for a lossless reciprocal network `Z = Zᵀ` exactly, so `Z₁₂/Z₂₁` can only carry the ratio of the fractions of each port's true terminal current its facet integral captures; the displacement route reads one facet set per gap and would see a 2% CAD/mesh difference in `A_gap/g` one-to-one, the conduction route would not. The `-n 2` re-run is *not* commissioned (a MUMPS solve moves at ≈1e-12 with width); the discriminator is the solid on both routes plus the two gap tags' area and length. (b) The mutual miss is the comparand, rule (f) class: a PEC tube links the flux through its *inner-edge* disc (`a − r_w = 35 mm`), and the receiver term alone at dipole order is −13% against the measured hole/solid 0.911; the symmetrised hole mutual is still −16.3%, so the asymmetry is not the mutual's miss. Re-registered comparands `_mutual_inductance(a, a − r_w, d)` and `(a − r_w, a − r_w, d)` are printed by **`TH-15` step 2e (§9 item 2)** together with the symmetrised-`S` unitarity identity asserted at 1e-9; a later review ratifies the bracket. Bands unmoved. Entry stays open until step 2 closes. |
| **Step 2e executed 2026-09-07 21:00 slot — mechanism (a) refuted, the bracket lands** | `20260908T020506Z_TH-15.log` (on the branch, `4275308`), `-n 4`, **2 failed / 18 passed**, `Status: 1`, elapsed **306 s** (`:1973–1974`); the two reds are this entry's, byte-identical. Symmetrised-`S`: `‖S_symᴴ S_sym − I‖_F` = **9.362447e-18**, `σ_max` = 1.000000000000 (`:1582–1583`) — but `result.s_matrix` is wave-assembled, and `z_to_s(Z_raw)`'s own non-unitarity is only **1.183730e-03** (`:1584`), so the 7.538037e-03 is 6.4× what the asymmetry explains. `Z₁₂/Z₂₁` asymmetry: hole 2.026999e-02, solid displacement **2.236184e-02**, solid **conduction 2.229415e-02** (`:1608–1610`) — the conduction route was predicted ≲ 0.5% and reads the same 2.2%; the two ports' `A_gap/g` agree to twelve digits (`:1617–1623`). The hole's symmetrised mutual against `M(a, a − r_w, d)`: **0.999879** (`:1630–1637`); `Z₂₁` alone 0.989641. |
| **Ruled 2026-09-08 03:00 review** | (1) `M(a, a − r_w, d)` is **ratified under rule (f)** as the hole mutual's comparand: a PEC tube carries no flux through its own cross-section, so the EMF is the same on every loop of its surface and the linked flux is the inner-edge disc's; the source-side filament at `a` is the disclosed omission. (2) The **wave-vs-`Z` distinction is real** — the wave-assembled `S` carries ≈ 6.4e-3 of non-unitarity beyond what the `Z` asymmetry explains, unattributed; step 2 cannot close on the symmetrised identity. (3) **No solid mutual reading stands**: 0.909618 is step 2d's *unboxed* fixture, 0.865226 the hole module's *boxed* solid; the same-mesh hole/solid ratio is 0.9576. (4) The asymmetry survives both current routes and both fixtures with identical gap geometry, so it is the **voltage** reading or the field: **`TH-15` step 2f (§9 item 2)** re-points the mutual assert at the ratified comparand (10% band unmoved), prints `_path_voltage` beside a gap-averaged `(g/V_gap)∫E·d̂ dV` on both ports and both fixtures with `Z` rebuilt on the average, and prints `‖S_wave − z_to_s(Z_raw)‖_F`. Unitarity and lossless bands unmoved; entry stays open until step 2 closes. |
| **Step 2f executed 2026-09-08 06:00 slot — the asymmetry is the `_path_voltage` *reading*, and the mutual half of this entry is green** | `20260908T111059Z_TH-15.log` (on the branch, `10b3af1`), `-n 4`, complex build, **1 failed / 21 passed in 348.49 s**, `Status: 1`, elapsed **350 s** (`:1806`, footer `:1957–1958`). The single red is this entry's raw-`S` unitarity gate, `7.538037e-03 > 1e-09`, byte-identical (`:1665, 1730`). **(a) Mutual, now green:** `Im Z₂₁ = 1.028789564e+00` Ω against the ratified `M(a, a − r_w, d) = 1.654508076658e-08` H (`ωM = 1.039558` Ω) reads **0.989641 (−1.04%)** inside the imported, unmoved `MUTUAL_TOLERANCE` 10% (`:1613`); the superseded `M(a, a, d)` ratio 0.828497 is printed beside (`:1614`). The header's "−17.15% against the imported 10%" is therefore **superseded by ruling (1)**, not loosened — the band never moved, the comparand was re-registered under rule (f). **(b) The `Z` asymmetry is the voltage reading** (`:1629–1646`): rebuilding `Z` on the gap-averaged `V̄` collapses `\|Z₁₂ − Z₂₁\|/\|Z₁₂\|` from **2.026999e-02 to 1.510620e-04** on the hole (134×) and **2.236184e-02 to 1.925413e-04** on the boxed solid (116×), with the probe's re-solve equal to the sweep's to round-off (5.07e-08 / 8.07e-08, `:1637, 1645`) and the sweep's own `Z` reproducing the `V_path` asymmetry to nine digits (`:1636, 1644`) — so nothing but the reading changed. The two **undriven** `V̄` values agree across drives to **ten digits** (6.789932166e-01 / 6.789932160e-01, `:1631–1632`) where the `V_path` arc samples differ by 2%. Predicted `\|V_path − V̄\|/\|V̄\| ~ 1–2%` is **refuted** — measured 49–54% undriven and 100.2% driven, because `V̄` is uncalibrated in magnitude (the driven gap volume contains the impressed source and the gap box is longer than `g`); only its reciprocity is meaningful here. **(c)** `‖S_wave − z_to_s(Z_raw)‖_F` = **2.915842e-02**, 4.6× the predicted 6.4e-3 and **essentially all off-diagonal** (8.43e-04 / 2.08e-02 / 2.04e-02 / 8.43e-04, `:1651–1657`) — ruling (2)'s 6.4e-3 was a difference of non-unitarity residuals, not this norm. **Consequence:** the entry stays open on the unitarity gate alone; the `_path_voltage` fix chunk is now licensed by a measurement, and its open question is calibration (restrict the average to the `g`-long gap slab, exclude the impressed source on the driven port) and whether that one change also moves the 2.9e-2 off-diagonal gap and the 7.5e-3 non-unitarity. No band moved. |
| **Ruled 2026-09-08 10:30 review — licensed twice, not yet specifiable; step 2g measures the calibration** | `OPS-41` (`7693a24`, same day) independently attributes the package fixture's 1e-4 width sensitivity to the point-sampled `V` (3.249e-04) against ≤ 3.6e-09 on `I`, so the `src/` replacement of `_path_voltage` is licensed by two fixtures and two routes. But 2f's `V̄` is uncalibrated by 49–100% for a geometric reason: the gap box is 10.4 mm square around a 10 mm wire (1.38× the section) and 13.95 mm long against an 11.95 mm chord (`_gap_half_extents`), and `V̄` was normalised by the box length over the tag's volume. **`TH-15` step 2g (§9 item 4)** prints, beside 2f's whole-tag reading A, the wire-footprint-restricted (B) and footprint-plus-chord-restricted (C) averages from a DG0 indicator on cell midpoints, with `Z` and the mutual rebuilt on each; reading A's reciprocity ≤ 1e-3 is asserted (backed by 1.510620e-04 / 1.925413e-04), C's agreement with `V_path` on the undriven ports (predicted ≤ 5%), C's reciprocity (≤ 1e-3) and C's mutual (inside 10%) are printed as predictions. The `src/` chunk is scoped from that table. The 03:00 "6.4e-3" is restated as a residual difference; the 2.915842e-02 off-diagonal wave-vs-`Z` gap stays unattributed. Unitarity gate red and unmoved; entry stays open until step 2 closes. |
| **Step 2g executed 2026-09-08 16:30 slot — the calibrated reading is B, not C; C's negative-result clause fires** | `20260908T213405Z_TH-15.log` (on the branch, `6f68956`), `-n 4`, complex build, `tests/environment` first, `-s`, `timeout -k 30 560`: **1 failed / 23 passed in 383.64 s**, `Status: 1`, elapsed **386 s**. The single red is this entry's raw-`S` unitarity gate, `7.538037e-03 > 1e-09`, byte-identical to 2f. **Anchor (iii) asserted and green:** reading A's rebuilt `Z` reciprocity **1.510600e-04** (hole) / **1.925424e-04** (solid) against the new `READING_A_RECIPROCITY_BAND` = 1e-3 (`:1666–1667`), reproducing 2f's 1.510620e-04 / 1.925413e-04. **Predictions (1), (2), (3) all fail for C and all but hold for B**, which the item did not separately predict. Undriven `\|V_path − V̄_X\|/\|V̄_X\|` (`:1674–1675, 1686–1687`): `miss_A` 50.577 / 53.716 / 49.341 / 52.725%, **`miss_B` 2.736 / 7.414 / 2.535 / 6.267%**, `miss_C` 114.867 / 307.792 / 108.537 / 265.311%. Rebuilt `Z` (`:1677–1679, 1689–1691`): reciprocity A 1.5106e-04 / 1.9254e-04, **B 2.434274e-02 / 1.412485e-02**, C 8.594116e-01 / 7.135082e-01; mutual against the ratified `M(a, a − r_w, d)` A −34.28% / −30.80%, **B −3.67% / +0.80% (inside the unmoved 10%)**, C −53.94% / −50.45%. Indicated volumes (`:1671–1672, 1683–1684`): `V_A` = 754.689 mm³, `V_B` ≈ 553.4–556.3, `V_C` ≈ 475.3–477.3 mm³, i.e. `V_C`/CAD 938.947 mm³ = **0.5069 (−49.3%)** — the item's CAD comparand looks 2× too large for this tag: `V_A` is exactly half the naive box (108.16 mm² × 13.955 mm = 1509.4 mm³) and the fixture's own `A_gap = V_gap/g` = 5.408e-05 m² is exactly half `(2(r_w + GAP_OVERHANG))²` (`:1603–1604`), so the gap **cell tag appears to be a half-domain**; against the half comparand 469.5 mm³, `V_C` is **+1.39%**, inside the predicted 5%. That reading is the executor's arithmetic on printed numbers, not a measurement, and no review has ruled it. Driven-port miss (no prediction): A 100.23%, B ≈ 102.2%, C ≈ 100.5%. `‖S_wave − z_to_s(Z_X)‖_F` = 2.81 / 2.83 / 2.82 (hole), 2.61 / 2.63 / 2.62 (solid) — all O(1) against 2f's 2.915842e-02, which was `z_to_s` of the *sweep's* `Z_raw`, not of a `V̄` route, so prediction (5) is uninformative as posed. **Mechanism suggested, untested:** the chord slab drops ≈ 14% of the volume (`V_C/V_B` ≈ 0.859) but ≈ 52% of the reading, so the field is concentrated in the `GAP_BURIAL` overhang beyond the chord and C amputates it asymmetrically. **Consequence:** C's reciprocity 8.6e-01 / 7.1e-01 is worse than the item's 1e-2 negative-result threshold, so §9 item 4 is marked 🚫 per its own clause; but the item's deliverable — the specification of the `src/` `_path_voltage` replacement — is *served*, and it points at the footprint restriction **B**. No band moved, no `src/` change, nothing loosened; entry stays open until step 2 closes. |
| **Ruled 2026-09-09 03:00 review — reading B is the basis, and the half-domain factor 2 is settled by geometry before any `src/` chunk is written** | 2g's own negative-result clause fired on reading C (reciprocity 8.594116e-01 / 7.135082e-01 against the item's 1e-2 threshold), so §9 item 4 was correctly marked 🚫 — but the deliverable was served by a reading the item did not separately predict: **B**, the wire-footprint restriction, reproduces the ratified `M(a, a − r_w, d)` at **−3.67% / +0.80%** inside the unmoved 10%, with undriven misses 2.736 / 7.414 / 2.535 / 6.267%. That is the calibration a `_path_voltage` replacement needs. What blocks specifying it is **arithmetic, not physics**: the executor's own reading — `V_C`/CAD = 0.5069, `V_A` exactly half the naive box, `A_gap` exactly half `(2(r_w + GAP_OVERHANG))²` — is that the gap **cell tag is a half-domain**, and that review explicitly recorded it as arithmetic on printed numbers, unmeasured and unruled. The factor 2 divides into **every** `V̄` normalisation on this fixture, so a `src/` chunk written now would bake it in. **`TH-15` step 2h (§9 item 2, `mesh-probe`, no solve)** measures the gap tag's volume, extents and cell count directly on the built mesh against both CAD candidates — the gap box and the wire-footprint-times-chord 938.6 mm³ — which differ by a factor 2, so a 5% assertion cannot straddle them and decides the question by itself; the two ports' tags must agree to ≤ 1e-3 as the symmetry control. No band moved, no `src/`, `TH-15` stays 🟡, both branches kept. The `src/` specification is written from step 2h's three numbers by the next review. Entry stays open until step 2 closes. |
| **Step 2h executed 2026-09-09 06:00 slot (`mesh-probe`, 🧪) — CONFIRMED: the gap cell tag is exactly half the gap box, and the split is the generator's own, documented in `src/`** | Two geometry-only windows, no solve, real build, `-n 2`, `timeout -k 30 300`, `-s`, **character-identical in every measured digit**: `20260909T110257Z_TH-15-step2h.log` (`Status: 1`, elapsed **58 s**) and `20260909T110426Z_TH-15-step2h.log` (`Status: 1`, **56 s**); `Status: 1` is the anchor's own negative-result exit, not a crash. Identical on both fixtures (hole, boxed solid) and both ports (`:513–519, 1438–1444`): `V_tag` = 7.546891363338e-07 m³ = **754.689136 mm³**, **`V_tag`/box = 0.500000** (CAD box 1509.378273 mm³), `V_tag`/(π r_w²·chord) = **0.803761** (938.947478 mm³), and `A_gap = V_tag/g` = 5.408000000000e-05 m² = 54.080000 mm² = **0.500000** of the box cross-section 108.160000 mm². Extents from owned cell midpoints (hole P1): `x` 1.017787e-02 and `y` 1.388758e-02 span the full box, **`z` spans 5.067458e-03 — half of 10.4 mm** — over `z ∈ [−2.5099e-02, −2.0032e-02]`, the half below the torus-1 centre plane at `z = −0.02`. Owned cells 12 585 / 12 632 (hole), 13 661 / 13 648 (solid, reproducing `OPS-39`'s `-n 1` census). **Negative control green:** `|V_P1 − V_P2|/|V_P1|` = **2.525310e-15** / **8.698290e-15** against 1e-3 (`:521, 1446`), so this is the mesh and not the probe. **The anchor's negative-result clause fired — `V_tag` matches neither named CAD candidate inside 5% (`:1450–1453`) — but the third candidate is not a guess.** With `emit_port_sheet=True` the generator splits each gap box at its mid-plane into two cell groups, `101`/`111` (gap 1 below/above) and `102`/`112`, and its docstring already states that a caller selecting the gap volume by tag "must take **both halves**" (`src/fem_em_solver/io/mesh.py:1176–1181, 1425–1426, 1455, 1542, 1582–1583`); the generator's fragment census printed above the mesh reads `gap_1 = gap_2 = gap_1_upper = gap_2_upper = 7.546891e-07` against `gap_box_analytic = 1.509378e-06` (`:40`) — four equal halves, two per port. So the whole gap box is `101 ∪ 111` (`102 ∪ 112`) and `GAP_TAGS = (101, 102)` selects half. **Consequence, for the review to rule:** `_tag_volume(GAP_TAGS[k])` reads half the gap box, so every `V̄` normalised by box length over `V_tag` on this fixture carries a factor 2; 2g's `V_C`/CAD = 0.5069 is `V_C` against the *full*-box-scale comparand, and against the half comparand 469.5 mm³ it is **+1.39%**. Nothing solved, no `src/` change, no test edited, no band moved, no record written; deliverable is `scripts/probes/th15_gap_tag_geometry.py`. `TH-15` stays 🟡, both branches kept. Entry stays open until step 2 closes. |

## ✅ RETIRED 2026-10-24 (weekly review, on `xl` entry 43's first branch; was 🟡 OPEN 2026-10-14 — `IMP-3` step 7 / I15 part 3, `xl` entry 33, 2026-10-14 02:00; ruled by the 2026-10-14 weekly review) — the rod's **pointwise tip ratio η_rod(T1) is not resolved by the cap-ball ladder**: r/4 → r/(4√2) → r/8 reads 3.482328 → 3.556262 → 3.530461, the third rung moving back

- **Retired 2026-10-24 (weekly):** `xl` entry 43 (I30, `20261023T070008Z_IMP-3-step10-I30.log:2654`) read η_rod(T1) = 3.539570 at the 4 mm ball, r/4 grain, degree 2, with every I14 identity green — inside the "Next" clause's [3.48, 3.56], so the spread is read as the ball-edge effect and this entry retires by its own clause. **Caveats recorded with the retirement (§10 rulings and epitaph 2026-10-24):** the interval was the spread of the three rungs it brackets, so "inside" separates little; at T1 the b4 mesh is finer than r/8's (edge 3.383e-4 against 5.082e-4 m); the 1 g ball mass and the degree-2 dof match were printed, not asserted, on that path. **Reading of record:** the pointwise tip ratio at degree 2 is 3.48–3.56 across four meshes (±1.1 %), a record on `IMP-3`'s row, not a converged value; no finer rung is bought (10-14 epitaph). The 1 g ball ratio Q (4.4617 / 4.4610) is the quantity (7c) will compare.

- **Test id:** none red — a pre-registered decision rule's "otherwise"
  branch. Module `tests/validation/test_imp3_tip_ladder_1g.py`
  (`FEM_EM_IMP3_STEP7=r8`, `-n 16`; `20261014T070008Z_IMP-3-step7-I15-r8.log`,
  13 passed, 511 s, 80.55 GiB, 342 072 cells). I14's identities (G1)–(G3),
  (i), (ii), (PP), (M), (E), (V), (D) are green on every rung at degree 2,
  and the 1 g ball ratio Q moves 0.28 % — the finding is on the pointwise
  tip ratio only.
- **Symptom:** with Δ₁ = η(r/(4√2)) − η(r/4) = +0.073934 and Δ₂ = η(r/8) −
  η(r/(4√2)) = −0.025801 the signs are opposite, so entry 33's rule ("same
  sign and \|Δ₂\| ≤ 0.7 \|Δ₁\| ⇒ a three-rung bracket; otherwise ⇒ the tip is
  not resolved by the cap ball alone") read its second branch. The spread
  of the three readings is 2.1 % of the ratio; the near-T1 edge ratios per
  rung were unequal (≈ 1.087, ≈ 1.264), so no order is claimed.
- **Reading (not a diagnosis):** the near-T1 edge — the mesh edge the
  pointwise sample sits on — lies outside the 2 mm cap ball on every rung,
  so the ladder refined cells around the point without refining the cell
  the point is read from; the ratio then follows that cell's grain, which
  the ladder did not control. The 1 g ball ratio, an integral, is
  unaffected.
- **Consequence:** `IMP-3`'s row carries the three values as a record, no
  bracket; no finer cap-ball rung is bought (§10 epitaph 2026-10-14).
- **Next:** §10 I17 — a cap-ball radius knob on `implant_rod` (`tip_ball_radius`,
  default 2 mm) with the near-T1 edge asserted inside a 4 mm ball; then I30 —
  I14's arm on the 4 mm ball at the r/4 grain, the ratio printed beside the
  three-rung record. Inside [3.48, 3.56] reads the spread as the ball-edge
  effect and retires this entry; outside, the pointwise tip ratio is
  recorded as unresolved at this grain on the row and the entry stays as
  the record. No band moves either way.
- **Item 464 (2026-10-15, I17 done):** with `tip_ball_radius` = 4 mm the
  130 near-T1 cells lie wholly inside the ball (largest vertex distance
  3.151970e-03 m from the cap centre; at 2 mm none of the 11 do), near-T1
  mean edge 3.383011e-04 m (6.983224e-04 m at 2 mm, ratio 2.0642); the
  build is 351 478 cells, over `ROD_CELL_GUARD_W2` = 320 000, priced
  564.7 s at `-n 8` (e = 1.37) — **I30 is not slot-feasible; it goes to the
  weekly as an `xl` candidate.** Entry not retired
  (`20261015T123700Z_IMP-3-I17-W2-n2.log`).

---

## Batch moved 2026-10-31

### ✅ RETIRED 2026-10-27 (`OPS-72`, §9 item 665, 07:30 CDT implementer slot — full `tests/environment` 44.90 s at `-n 1`, 44.73 s at `-n 2`; the environment-first ruling is lifted by the next daily review, not here) — was 🟡 OPEN 2026-10-26 (`MAT-3` step 1b, 21:00 CDT implementer slot, §9 item 646): `tests/environment` alone now costs ≈ 145 s at `-n 1`, so a standard-tier (180 s) window with the full `tests/environment` first cannot fit a solve

| | |
| --- | --- |
| **Symptom** | Exit 124 at the ceiling with no assertion failing. In the 2-rank stream, rank 1 reports the current node PASSED while rank 0 is still inside it at SIGTERM. That looks like a deadlock, but it is lag. `MAT-3` step 1b W1 (full `tests/environment` + two plane-wave modules, `timeout -k 30 180`) was killed at 181 s, and step 1's own Debye module behind the full `tests/environment` was killed at 180 s. That module passed at 100 s on 2026-09-29 (`20260929T183252Z_MAT-3-step1.log.gz`). |
| **Verified at** | `b1b3fd28`, by execution, harness logs on `attempt/MAT-3-step1b-20261026T0217Z`. Bisect at `-n 2`, each window `test_complex_mode.py` + one env module + the Debye module (33 s): `+ test_scikit_rf_version.py` green in 35 s (`20261026T021046Z_MAT-3-step1b-bisectA.log`); `+ test_dolfinx_version.py` green in 39 s (`…021129Z_…-bisectB.log`); `+ test_dolfinx_api_migration.py` killed at 120 s with rank 0 still inside `test_negative_control_detects_reverted_migrations` (`…021215Z_…-bisectC.log`). Priced alone at `-n 1`: **3 passed in 144.42 s**, of which `test_negative_control_detects_reverted_migrations` 78.76 s, `test_no_unmigrated_call_sites` 39.01 s, `test_filed_survivors_outside_the_gated_roots_are_unchanged` 26.63 s (`20261026T021500Z_MAT-3-step1b-envprice.log`). |
| **Cause** | Cost growth, not a defect in any gate. `scripts/testing/check_dolfinx_api_migration.py` parses and resolves every `.py` under `src` + `tests` twice: once for the baseline and once over a temp copy for the negative control. It also sweeps `examples` + `scripts`. Both ranks run all of it redundantly. The trees have grown since `OPS-26` priced the module as smoke (`tests/` is 1 600+ files, 38 MB). Not measured: how much the 2-rank contention adds over the `-n 1` price. |
| **Scope** | Any brief that puts the full `tests/environment` first in a standard-tier window, or in a heavy window sized with ≈ 0 s for it. That includes `TH-5` step 3c (§9 item 644, killed at 300 s with the same rank-0-lag signature). `test_complex_mode.py` alone costs < 1 s, the 10:30 smoke rule's choice. No test is red. |
| **Resolves with** | A review ruling: either briefs put only `tests/environment/test_complex_mode.py` first (the regression the CLAUDE.md "environment first" rule is for) and run the migration sweep in its own window, or an `OPS` chunk makes the sweep rank-0-only or cheaper. Do not edit the checker in passing. |
| **Ruled 2026-10-26 (03:00 daily review)** | **The workaround binds; the cost stands, so the entry stays open.** A solve window puts `tests/environment/test_complex_mode.py` first. The full `tests/environment` runs once per item that adds or edits a `.py`, in its own window `W-env` (`-n 1`, `timeout -k 30 300`, priced at 146 s), and never shares a window with a solve. Briefs written before the ruling are read this way (§9 digest; journal `docs/planning/reviews/2026-10-26-daily.md` §4). Whether to make the sweep cheaper (an `OPS` chunk) is returned to the 10-26 weekly. **Correction to the Symptom row** (`log-pathologist`, 2026-10-26 03:00, on `TH-5` 3c's W1): "rank 1 reports the current node PASSED" is an interleave misreading. The stream holds 13 `PASSED` tokens per rank, so neither rank finished the solve node. The trailing `PASSED` on rank 1's node line is rank 0's last migration verdict. The real cost is both ranks running the sweep, not rank-0 lag. Count verdict tokens per rank, never by line adjacency. |
| **Owner (2026-10-26 weekly)** | The 03:00 ruling is ratified and binds until the cost is removed. The cost gets a chunk, **`OPS-72`** (PROJECT_PLAN §7, lane 1): the sweep runs once (rank 0, result broadcast), the planted negative control still red, every `tests/environment` verdict unchanged, the full `tests/environment` at `-n 1` and `-n 2` ≤ 60 s measured against this entry's 144.42 s. This entry retires on that measurement, and the daily then lifts the environment-first clause from the briefs. Never edited in passing. |
| **Retired 2026-10-27 (`OPS-72`, §9 item 665)** | W0's profile (`20261027T124634Z_OPS-72-W0-profile.log`) put 52 of 108 profiled s of one `src` + `tests` sweep in `enclosing_bindings` (a full subtree re-walk per call site). The checker now precomputes each scope's contained node ids once and walks each file's tree once (internal; `RESULT:` line identical before/after, `files=582 calls=1471 apis=39 methods=17 violations=0 uncheckable=0 shadowed=0`, `…124558Z_OPS-72-W-rb.log` / `…124919Z_OPS-72-W-ra.log`). The gate module computes the `src` + `tests` sweep once, on rank 0, and broadcasts it; the negative control takes it as `baseline=` and still sweeps the whole temp copy. No assertion, floor or reversion edited. W1 `-n 1`: 12 passed in 44.90 s, same 12 IDs as `20261027T003839Z_WF-15-s8-Wenv.log` (`…124941Z_OPS-72-W1-n1.log`); W2 `-n 2`: 12 passed per rank (24 `PASSED` tokens, 0 failed) in 44.73 s (`…125038Z_OPS-72-W2-n2.log`). |

## ✅ RETIRED 2026-10-26 (weekly review, on the 10-24 condition's first clause — H33's width law at both orders; re-filed in PROJECT_PLAN §2.2 (4) as a record; H34 withdrawn, §10 epitaph) — was 🟡 OPEN 2026-10-11 (`ANS-8` step 2 / H21, §9 item 369; re-titled 2026-10-22 by §10 H32, §9 item 584, on §10 H30's reading; weekly ruling 2026-10-12 (3) names chain F, ownership of the re-titled finding is the weekly's per the 2026-10-22 03:00 review) — a finding on the lumped-sheet model: **the field's dissipation on a lumped feed sheet exceeds the port's lumped Norton term by the sheet-field variance X_s = CS_s, and the excess does not shrink with order** — ΣX_s / P_acc 8.316822e-03 at degree 1 (118.2437 pF), 2.149715e-02 at degree 2 (115.5209 pF), on the tuned copper F-human

- **Test id:** the parked stage V3 of `ans:8` (`attempt/ANS-8-v3-20261011T1420Z`), assertion (P_acc-2): **8.317e-3 against `DISCRETE_IDENTITY_RTOL` 1e-6** (`20261011T141650Z_ANS-8-v3-V3.log:6197`). Every other anchor in the window is green: (P_w) 5.055e-13, (P_acc-1) 6.808e-16, and the reproductions of items 358 / 359 (`:6177–6226`). The same gap was printed by `MAT-15` step 5 (item 361, 8.75e-3 of the field net, `20261011T003520Z_MAT-15-step5-W1.log:6297`).
- **What is measured:** the field's own energy balance closes, (P_w) at 5e-13. The lumped line cut equals the wave form, (P_acc-1). The field source power equals the lumped Norton source power, terminal + Σ ½ Re(1/Z_th)\|V\|², to every printed digit, so the terms linear in E agree. The whole gap is the two feed sheets' field dissipation ½ Re(Y_s)∫\|n̂×E\|² dA in excess of ½ Re(1/Z_th)\|V\|²: 1.329e-4 W on P17 and 1.214e-4 W on P21 (`20261011T142050Z_ANS-8-v3-diag.log:37–39`, the `log-pathologist`'s arithmetic).
- **Reading (not a diagnosis):** this is the Cauchy–Schwarz excess ½ Re(Y_s)(∫\|n̂×E\|² − \|∫E·ĥ\|²/A) that a field varying across the sheet produces (`tests/validation/test_birdcage_power_identity.py:272–293`), i.e. a property of the discrete sheet, not of the ESR term (which cancels) or of the solve. Unmeasured.
- **Consequence:** no assertion equates the field's net feed power with the terminal power. Per-W figures state their normaliser: `ans:8`'s v3 rows use the terminal line cut (the SPEC's definition) and print the field-partition one beside it, and P20 (§9 item 380) does the same (2026-10-11 03:00 ruling 2). Row 5's fractions of the line cut sum to 0.9917, and the residual is this gap.
- **Mechanism, measured 2026-10-11 (§9 item 386, `20261011T170627Z_ANS-8-v3b-V3.log:6204–6208`, printed, never asserted):** on each feed sheet the Cauchy–Schwarz excess ½ Re(Y_s)(∫\|n̂×E\|² dA − \|∫E·ĥ dA\|²/A_s) equals the field-minus-lumped excess X_s: CS_s / X_s = 1.000000000 on P17 (1.329324349e-4 W) and P21 (1.213651853e-4 W); ΣX_s = the gap, 2.542976202e-4 W. Split per sheet: in-plane component across ĥ 9.402e-5 / 8.562e-5 W, drive-component variance 3.892e-5 / 3.574e-5 W (≈ 70 / 30). P_src = Norton lumped power to 3.3e-16; A_s = w·h to 10 digits (exact rectangles). The gap is the field's non-uniformity over the sheet, which the lumped terminal reduction cannot see; the transverse part dominates.
- **Next:** the split above landed inside its pre-registered 1 % window. Whether the gap is a documented sheet property or an issue, and which normalisation the AED column reads, is the 10-12 weekly's ruling. Degree 2 (`xl` entry 35) is the order move.
- **Ruled 2026-10-12 (weekly):** (1) the gap is a **documented property of the discrete feed sheet** — the Cauchy–Schwarz excess of a field that varies across the sheet, equal to the measured excess to 1e-9 on both sheets, with the field's own energy balance closed at 5e-13 — not a defect in the solve, the ESR term or the port model; no assertion is added or changed. (2) **The normaliser every per-watt row of `ans:8` step 2 reads is the lumped line cut, P_acc,SPEC = ½ Re Σ_feeds V_port · conj(I_line) = incident − reflected at the two 50 Ω lumped ports** — the quantity HFSS reports as accepted power — with the field-partition figure printed beside it (SPEC v3 note added tonight). (3) The entry **retires when `xl` entry 35 prints the excess at degree 2**: a sheet-discretisation property shrinks with the order; if it does not shrink, the entry is re-titled as a finding on the lumped-sheet model and chain F owns it.
- **Checked 2026-10-18 (03:00 daily review, `log-pathologist` ruling; the 10-17 weekly's "retires when the daily confirms" check): NOT PRINTED — the entry stays 🟡 OPEN.** `20261016T070008Z_WF-10-step4-H20.log` (entry 35) prints only the field-side partitions (P_RC) and (P_w) at degree 2 (`:6827`, `:6830–6832`, `:13385–13390`). The lumped terminal side ½ Re(1/Z_th)\|V\|² and the per-sheet CS_s / X_s never appear; a search for the v3b diagnostic's terms over the log and both modules returns nothing. Entry 35's own Readout never asked for the excess, so ruling (3)'s trigger cannot fire from that window. It also ran at the degree-1 node, where the degree-2 coil is detuned (Im Y_in +1.346e-02 S, `:6836`). **For the weekly:** re-pose the trigger onto a degree-2 window that carries the `ans:8` v3b diagnostic (`20261011T170627Z_ANS-8-v3b-V3.log:6204–6208`) at the degree-2 root, e.g. §10 H30 after entry 41. Under rule (6) the daily may not change a retirement condition.
- **Re-posed 2026-10-19 (weekly):** ruling (3)'s trigger moves from entry 35 to **§10 H30** — the `ans:8` v3 rows 3–5 at the degree-2 root from entry 41's record, which carries the v3b diagnostic (CS_s / X_s per feed sheet, ΣX_s / P_acc) at degree 2. The entry **retires if the excess shrinks with order (ΣX_s / P_acc < 8.3e-3 at degree 2)**; if it does not shrink, the entry is re-titled a finding on the lumped-sheet model. Nothing else in ruling (3) moves.
- **Re-titled 2026-10-22 (§10 H32, §9 item 584; the re-pose's "does not shrink" branch, executed):** §10 H30 (§9 item 582) printed, at the degree-2 root from `xl` entry 41's record (`logs/wf10/h9p-node-115.5209-d2.json`, written at `20261022T070008Z_WF-10-step4-H28.log:13525`), X_s 3.262132656e-04 / 2.994809404e-04 W on P17 / P21 and **ΣX_s / P_acc 2.149715e-02** (`20261022T093054Z_ANS-8-H30-W1.log:87–89`, printed, never asserted) against 8.3e-3 — 2.58× the degree-1 8.316822e-03 (`20261021T004407Z_ANS-8-H30a-W1.log:85`). Two points (degree 1 and degree 2, one mesh), not a trend. CS_s is not in entry 41's record; X_s = CS_s is the sheet model's own identity (CS_s / X_s = 1.000000000 at degree 1, `20261011T170627Z_ANS-8-v3b-V3.log:6204–6208`). Ruling (1)'s reading (a documented property of the discrete feed sheet, not of the solve, the ESR term or the port model) and ruling (2)'s normaliser stand; no assertion is added or changed. Ruling (3) names chain F; which chain owns the re-titled finding is the weekly's (2026-10-22 03:00 review ruling on item 584), and the entry stays under `ANS-8` / known-issues. This re-title sets no new retirement condition — that is the weekly's.
- **Owner and retirement condition set 2026-10-24 (weekly, §10 rulings):** **chain H owns the finding** (chain F closed with T6 / F27; ruling (3)'s name is superseded). The excess grows with order because a finite sheet's edge field is better resolved — what a property of the sheet *model* looks like — but two points on one mesh are not a trend, so the question is posed on the axis it lives on. **§10 H33:** X_s / P_acc on the F-small exact-rectangle four-port at three sheet widths (w₀ / 2, w₀, 2 w₀) at degree 1 and 2, CS_s / X_s = 1 at `DISCRETE_IDENTITY_RTOL` asserted, the width exponent per order printed. **§10 H34:** Palace's lumped port on `OPS-61`'s shared mesh at degree 2 — its field-side dissipation against its reported port power, the same arithmetic. **Retires when H33 shows X_s / P_acc monotone in the width at both orders (then re-filed in §2.2 as a record with the width law), or when H34 shows Palace's excess equals ours within 5 % on the same sheet;** if neither, the entry stays a finding and `ANS-8` step 2's adjudication reads AED's accepted-power definition against ours (SPEC v4 names the line cut). Phase 6 exited 2026-10-24 with this finding travelling beside it (§6).
- **Retired 2026-10-26 (weekly, §10 rulings) — the first clause fired.** H33a / H33b (§9 items 624 / 625; the daily's re-pose of the widths to w₀/4, w₀/2, w₀ ratified — 2w₀ cannot be built) measured ΣX_s / P_acc on the F-small exact-rectangle four-port at 64 MHz: **4.274e-04 / 9.958e-05 / 4.135e-05 at degree 1 (LSQ exponent −1.685)** and **8.826e-04 / 1.683e-04 / 4.838e-05 at degree 2 (≈ −2.09)**, with CS_s / X_s = 1 at ≤ 2.6e-10, (P_w) ≤ 7.5e-15, (P_acc-1) ≤ 5.1e-13 on every sheet and the record width reproducing F18a's S at (R0) 9.557e-11 (`20261025T095010Z_ANS-4-H33a-W1.log:5662–5665`; `20261025T100927Z` / `100355Z` / `095803Z_ANS-4-H33b-W{3,2,1}.log` TABLE lines; `ports:23`). Monotone, falling, in the width at both orders and larger at degree 2 at every width — a documented sheet-resolution property with a width law, which is what the 10-24 condition named. Rulings (1) and (2) stand: no assertion equates field and terminal feed power; every per-watt figure states its normaliser (the lumped line cut) with the field partition beside it. H34 is withdrawn (no Palace field export exists; `OPS-61`'s mesh is not H33's fixture). The F-human values (8.3e-3 / 2.1e-2) are the same law on sheets narrower relative to their field; not re-measured.

## ✅ RETIRED 2026-10-26 (weekly review, by the entry's own first branch — C13 read the class's absolute spread at 0.888× the adjacent class's: a relative artefact of a small entry, re-filed as a record) — was 🟡 OPEN 2026-10-24 (`GEO-38` C3, `xxl` entry 44, 2026-10-24 02:00; ruled by the 2026-10-24 weekly review, decision-rule branch 2) — F-clinical-1.5T's degree-1 untuned 64×64 **breaks its C32 × mirror symmetry in one class of 34**: the "other ring, one step" class spreads **0.9262 %** against the 0.5 % `ADJACENT_SPREAD_BAND` the readout named, while the other 33 classes are ≤ 0.3217 % and reciprocity reads 1.887e-14

- **Test id:** none red — the module (`tests/validation/test_geo38_f_clinical_matrix.py`, `FEM_EM_GEO38_C3=1`, `-n 16`) asserts the classes at `OPPOSITE_SPREAD_BAND` (5 %) and passed on every rank (`20261024T070008Z_GEO-38-C3.log:14784`, worst 9.2621e-03); the readout entry 44 was bought against pre-registered `ADJACENT_SPREAD_BAND` (0.5 %). Item 519 disclosed the substitution in the module docstring (`:24–28`) before the window; no review ruled on it. **Ruled 2026-10-24 (weekly): a disclosed substitution is a question, not a ruling (§10 epitaph) — the readout's band is the band of record, so the window reads branch 2 of its decision rule: "a class-spread red with reciprocity green ⇒ the fixture's symmetry is the finding, no band moves".** The module's 5 % assertion stands as the module's gate; it is not the gate of record for any class a readout names at 0.5 %.
- **What is green in the same window:** reciprocity 1.887168e-14 against `RECIPROCITY_BAND` 1e-3 (`:14748`); σ_max 0.999970196163 ≤ 1 (`:14749`); `TH-14`'s copper identity worst rel 2.177e-12 against 1e-6, NEG ≥ 1.6173e-04 (`:14747`); 1 solve + 63 held drives, 1 147 366 cells, U₁ 1 525 193, 1 233 s, 46.85 GiB. The circuit layer's view (C10, `20261024T093056Z_GEO-38-C10-W1.log`): anti(k = 1) 1.000000 and off-block 1.079e-03 — the k = 1 block is C32-clean to ≈ 1e-3 while one small-entry class is not.
- **Reading (not a diagnosis):** F-human's degree-1 32×32 held 18 C16 × mirror classes at ≤ 0.45 % (`PORT-13`) and the degree-2 32×32 is C16-symmetric (entry 36); a relative spread of 0.93 % on a coupling between a port on one ring and the port one step over on the *other* ring — a small-magnitude entry — may be the mesh's asymmetry seen through a small denominator, or the fixture's (the leg-box corner, the shield can, the box placement). Unmeasured which.
- **Consequence:** no gate on F-clinical quotes a class spread at 0.5 % until this is read; C4 (entry 46, the in-model root on the two-port feed) does not read the class spreads and proceeds; the 3 T twin's 64×64 (entry 47) prints its worst class against 0.5 % under the same ruling.
- **Resolves with:** §10 **C13** (lane 1, no solve) — the class's anatomy from the saved `logs/geo38/c3-cu-s64-d1.npz`: per class the member \|S\| mean, absolute and relative spread, the "other ring, one step" members listed with their (slot, ring) pairs, and the same statistics on F-human's degree-1 32×32 beside them. The 10-26 weekly writes the retirement condition from it: a relative artefact of a small entry (absolute spread below the adjacent class's) ⇒ re-filed as a record; an absolute spread above the adjacent class's ⇒ the fixture's symmetry, read at degree 2 by C5 after C5a prices it. **Neither band moves; no finer degree-1 64×64 is bought before C13 reads** (§10 declined, 10-03 epitaph).
- **Retired 2026-10-26 (weekly, §10 rulings) — the first branch.** C13 (§9 item 622, `20261025T093429Z_GEO-38-C13-W1.log:67–70`, `:105`, `:234–237`, `:257–261`): (RR) / (CT) green on 34 / 34 lines; the "other ring, one step" class's **absolute** spread (max − min of \|S\|) is 5.318e-4, **0.8883×** the "same ring, one step" class's 5.987e-4, on a class mean of 5.74e-2 — the 0.9262 % is a small denominator, not a larger asymmetry; F-human's degree-1 32×32 reads the same way (worst class "other ring, 2 steps" 3.2748 % relative at absolute ratio 0.6181). No band moves; the module's 5 % assertion stays the module's gate. **Standing reading:** a readout that names a relative class band on a small-entry class names C13's statistic beside it — entry 47's rule (the 3 T twin) is annotated so, and C15 prints the anatomy on the 3 T matrix; a known-issues row is opened only if a class's absolute spread exceeds the adjacent class's.

---

## Batch moved 2026-11-07

### ✅ RETIRED 2026-11-05 (`OPS-81`, §9 item 828, 12:00 implementer slot — the scope docstrings and the guard landed; the positive control raised and the callers' gates held) — was 🟡 OPEN 2026-11-03 (`GEO-36` I57, §9 item 801, 13:30 implementer slot; filed by the 18:00 daily review): the birdcage power-identity helper sums volume loss over tags 3 and 1 only, so on a mesh with other lossy cells the (P) partition omits them

| | |
| --- | --- |
| **Symptom** | On I41's kept-rod F-small build (182 371 cells, tag 403 = the rod as σ_wire cells), the (P) identity read P1 2.882e-06 against `DISCRETE_IDENTITY_RTOL` 1e-6. The residual is exactly tag 403's loss (1.505772888e-08 W of P_src 5.225076965e-03 W); with P_vol over all cells it closes to ≤ 9.6e-15 on every drive (`20261103T195044Z_GEO-36-I57-diag-W1.log:1919–1942`, on `attempt/GEO-36-I57-diag-20261103T1952Z`). The same omission at σ_gel gave the gel twin's P1 2.853e-07, green but under-summed (`20261103T120306Z_GEO-36-I57-W1.log:1974–1977`). |
| **Verified at** | `main` after item 801, by execution: `-n 2`, 124 s; quadrature flat in degree, algebraic residual 6.4e-13. |
| **Cause** | `_loss_power_w` (`tests/validation/test_port_drive_superposition.py:205–221`), used by `_exact_shares`, integrates σ\|E\|² over `PHANTOM_CELL_TAG` (3) and `CONDUCTOR_CELL_TAG` (1) only — written for fixtures whose only lossy volumes are those two. Entry 48's identity uses the all-cells integral (`test_geo36_rod_f_clinical_drive.py:381`), which is why it was green. |
| **Scope** | No test on `main` is red. Any caller of `_exact_shares` / `_loss_power_w` on a mesh with lossy cells outside tags 3 and 1 (a kept rod or disc on tag 403) under-sums its partition. Hole builds (tag 402 is a facet) are unaffected. **Callers audited 2026-11-05 (`OPS-80`, §9 item 817, code reading, `20261105T013028Z_OPS-80-grep.log`: 46 lines in 16 files, 20 calls, table in `docs/planning/chunks/OPS-80.md`): no other (K) two-tag gate.** Every direct call runs on a hole build or a build without an implant (U). Only `_rb._arm` (`test_imp3_rod_birdcage.py:329`) reaches a kept build, and only through three callers. `test_imp3_tip_q_entry48_def.py` gates (P) over all cells and keeps the two-tag form as an asserted control. `test_imp3_rod_heating_matched.py:171` (I41 arm M, 403 → gel) and `examples/implants/15_rod_heating_matched_phantom.py:164` compute two-tag shares but read no (P). |
| **Resolves with** | I57's re-run with the all-cells form (§9 item 804, `docs/planning/chunks/GEO-36.md` "Review ruling 2026-11-03 18:00"), plus an audit of the other callers on kept builds; whether the shared helper is widened is the weekly's call. Do not edit the helper in passing — other modules import it. **First corrected caller (2026-11-04, §9 item 810):** `tests/validation/test_imp3_tip_q_entry48_def.py` gates (P) over all cells, ≤ 8.3e-15 on every arm (`20261104T104820Z_GEO-36-I57-W1d.log:1968–1994`); the helper is unedited and the other callers are still unaudited, so the entry stays open. **The callers' audit is queued** as `OPS-80` (§9 item 817, 2026-11-04 18:00 review; code reading only, `docs/planning/chunks/OPS-80.md`). **Weekly ruling 2026-11-04: the helper is not widened.** `OPS-80` found no other (K) two-tag gate; `_exact_shares` returns the phantom / conductor shares that 20 callers read as a partition, and widening its sum would change that partition's meaning everywhere while fixing nothing red. The rule: every kept-cell fixture gates (P) on the all-cells integral (entry 48's form, I57's). **This entry retires with `OPS-81`** (lane 1): the helper's docstring states its two-tag scope and names the all-cells form, and — if the case object carries its material map — `_loss_power_w` raises on a lossy tag outside (1, 3); the docstring alone otherwise. |
| **Retired by** | `OPS-81` (the commit landing §9 item 828): the two-tag scope stated, the all-cells form named; `_assert_two_tag_scope` raises on a lossy tag outside (1, 3) unless the case carries `TWO_TAG_ACK_KEY` (positive control `20261105T182800Z_OPS-81-W0.log:63`, 12 passed per rank); the acknowledging wrappers and the (U) caller stay green (`20261105T182816Z_OPS-81-W1-entry48.log:2028` 14 passed / 3 skipped per rank, control 2.8818e-06 at `:2019`; `20261105T183330Z_OPS-81-W2-i41.log:3791` 16 passed per rank; `20261105T183847Z_OPS-81-W3-port16.log:2007` 18 passed per rank, 0 failed, no refusal); the helper is not widened (weekly 2026-11-04). |

---

## Batch moved 2026-11-09

## ✅ RETIRED 2026-11-09 (§9 item 918, 04:30 implementer slot — the premise tests read empty tmp log roots; three W0s green, the control red as registered) — was 🟡 OPEN (2026-11-09 03:00 daily review; read from the code, not observed in a window) — the `xl` readers' "(Pend)" premise tests go red the night their window's log lands

- **Test ids:**
  - `tests/validation/test_geo38_c33_field_strength_reader.py::test_c44_pend_and_knob` and `::test_c44_n_both` — **red on `main` since `3b5b7e1f`**, the commit that tracked `20261109T080008Z_GEO-38-C6-3T.log` (`xl` entry 58);
  - `tests/validation/test_geo38_c39_thermal_reader.py::test_c46_pend_and_knob` — from entry 59's log, Tuesday 11-10 02:00;
  - `tests/validation/test_geo36_c40_rod_frequency_reader.py::test_c45_pend_and_knob` / `::test_c45_n_both` — from entry 61's log, 11-11 or 11-12.
- **Literal symptom (predicted from the code; no window has run these nodes since 02:13):**
  - `assert _on_root(find_log, ENTRY58_GLOB) is None` (`test_geo38_c33_field_strength_reader.py:450`) now finds the tracked log.
  - In `test_c44_n_both`, `_c33_body("both")` no longer raises (`:443`). It writes `c33-1p5t-3t.json` under `/workspace/logs/geo38` (gitignored), and then `assert not … .exists()` fails (`:445`).
- **Cause:** items 905–907 (C44–C46) pinned "entry N's log is absent" against the real tree as a premise. That premise is true only until the window runs.
- **Consequence:** no gate or band is affected. A slot that runs one of these modules unknobbed sees a red that is not its own. Item 917 (C33's `both` window) cannot be green until this is fixed.
- **Resolves with §9 item 918:** the premise and n_both tests monkeypatch `TRACKED_LOG_DIR` / `HELD_LOG_DIR` to empty tmp roots, so each control still fires on an absent entry by construction. C44's module also asserts the real tree reads entry 58 as "tracked". No reader, gate or control is dropped. This entry retires with 918's three green W0s.
- **Verified at:** `3b5b7e1f` (code read; `find_log`'s default roots are read at call time, `test_geo38_c33_field_strength_reader.py:90–96`).
- **Retired by:** §9 item 918. A helper `_empty_log_roots` in each of the three modules monkeypatches `TRACKED_LOG_DIR` / `HELD_LOG_DIR` to two empty `tmp_path` dirs (made on rank 0, bcast); the five premise / n_both tests call it, and C44 gains `test_c44_real58` on the real tree ("tracked", `20261109T080008Z_GEO-38-C6-3T.log`). W0s green per rank at `-n 2`: C44 + ans11 page 30 passed / 1 skipped (`20261109T103136Z_GEO-38-item918-W0-C44.log`), C45 15 / 1 (`20261109T103147Z_GEO-38-item918-W0-C45.log`), C46 13 / 1 (`20261109T103155Z_GEO-38-item918-W0-C46.log`). Control: the C44 copy with the pend test's patch removed is red on that node alone, naming entry 58's base name (`20261109T103211Z_GEO-38-item918-Wctl.log:106`, 1 failed / 7 passed at `-n 1`). (Dry): `logs/geo38` and `logs/geo36` listings identical before and after (`20261109T103130Z_GEO-38-item918-ls-pre.log`, `20261109T103220Z_GEO-38-item918-ls-post.log`); no `c33-1p5t-3t.json` written. **Second premise, same slot (item 917):** the n_both / pend tests also asserted the reader's host JSON *absent*, which is false once that reader's own `both` window has run (`20261109T103446Z_GEO-38-C33-W1-n1.log:143–147`, red, kept). The three tests now compare rank 0's `(mtime_ns, bytes)` before and after; W0-r2s and C33's W1-r2 / W1-n1-r2 are green (`20261109T103539Z_GEO-38-item918-W0-C44-r2.log`, `20261109T103547Z_GEO-38-item918-W0-C45-r2.log`, `20261109T103550Z_GEO-38-item918-W0-C46-r2.log`, `20261109T103601Z_GEO-38-C33-W1-n1-r2.log`).

---

## Batch moved 2026-11-11

### RETIRED 2026-11-10 (18:00 daily review's ruling, `WF-17.md` § "Review ruling (2026-11-10 18:00 daily": the print equals the node's own named record `WF15_DEG2_RECORD` at every digit, and the degree-2 sphere SAR's `.10e` print scatters 4–10 units in its last place between same-commit, same-route `-n 2` runs — `log-pathologist`; a solve output's print is gated against its named constant by its band, not against the last run's print; step 39 closes on these windows; provisional, the 11-11 weekly ratifies or reverts) — was 🟡 OPEN 2026-11-10 (`WF-17` step 39 / U148, §9 item 943, 13:30 implementer slot; observed, not fixed): the degree-2 first-study node's (R) SAR `.10e` print reads 8.0116577245e-08 against its prior 8.0116577252e-08 — a (b) record seven units off in its last place, every band green

| | |
| --- | --- |
| **Symptom** | `test_wf17_first_study_tutorial.py::test_wf17_first_study_block_reproduces_wf15_record` (`-n 2`, complex) passes, but its (R) line reads `SAR = 8.011657724500328e-08 (8.0116577245e-08); record 8.0116577245e-08; rel 4.097e-14` (`20261110T193150Z_WF-17-s39-W-first.log:138`). The prior was `8.011657725185266e-08 (8.0116577252e-08); … rel 8.553e-11` (`20261010T124220Z_OPS-70-S2-W7-wf17-first-study-deg2.log.gz:162`). The value now sits closer to `WF15_DEG2_RECORD` (8.0116577245e-08, unchanged at `tests/package/test_wf17_first_study_tutorial.py:72`). The other records are at their priors: (K) `[]`; (Deg) 2 / 2 / 2; cells 74020 and N1curl dofs 482158; (G) 0.0413 % (`:137`, `:139`). Every band is green. |
| **Verified at** | `main` at `2977d6b4` by execution: 294 s, `[capture] rc=0`, 5 passed per rank. W-rep (`20261110T193043Z_WF-17-s39-W-four-rep.log`, 6 passed per rank, 56 s) ran on the same tree just before it. Its (c4) degree-1 SAR printed 8.29155812850687e-08, against 8.291558128523365e-08 in item 938's W-four (`20261110T121910Z_WF-17-s39-W-four.log:1919`) on the same code. Its distance to `MAT4_DEG1_RECORD` printed 6.113e-11, against 6.312e-11 there and 6.299e-11 at the 10-06 prior. The `.10e` print, 8.2915581285e-08, was unchanged. |
| **Cause** | Not established. The evidence points to run-to-run nondeterminism, not a merge. W-rep's degree-1 digits differ from item 938's on identical code at ≈ 2e-13 relative. Across every tracked log, the degree-2 sphere SAR's `.10e` print takes **11 distinct values, from 8.0116577233e-08 to 8.0116577256e-08** (≈ 2.9e-11 relative). The commonest is 8.0116577245e-08, seen 20 times (`zgrep -ho "8\.0116577[0-9]*e-08" docs/testing/logs/*`). Whether those logs share mesh, width and route was not checked. So the degree-2 `.10e` print does not hold steady at its eleventh digit run to run. |
| **Scope** | Records only. No band is red and no test on `main` is red. `WF-17` is held at 🟡 per the brief's rule and the 10:30 ruling ("two or more units is BLOCKED"). |
| **Resolves with** | The review's ruling on whether the degree-2 (R) `.10e` print is a record under the 10:30 ruling's (b), or a reading because it is not reproducible at its last digit. If it is a reading, step 39 closes on this window. Fix nothing in passing. |

### RETIRED 2026-11-10 (10:30 daily review's ruling, `WF-17.md` § "Review ruling (2026-11-10 10:30 daily": the moved value is a band-gated residual, not a record — `log-pathologist`, no asserted or fixed-precision record moved; the step continues as item 943 / 944 / 945) — was 🟡 OPEN 2026-11-10 (`WF-17` step 42 / U151, §9 item 940, 09:00 implementer slot; observed, not fixed): the S12 copper node's (g) reciprocity residual prints 6.478e-16 against its prior 6.315e-16, a record digit off its prior with every band green

| | |
| --- | --- |
| **Symptom** | `test_wf17_copper_yaml_tutorial.py::test_wf17_s12_copper_block_reproduces_ans6` (`-n 2`, complex, `FEM_EM_WF17_S12=1`) passes, but its (g) line reads `reciprocity 6.478e-16 (ASSERTED <= 1e-03)` (`20261110T150154Z_WF-17-s42-Wcu.log:1785`) against the prior 6.315e-16 (`20261025T003309Z_WF-17-s12-Wcu.log.gz:1785`). The new value equals the `RGCu:2718` record the node cites on the same line. Every other record is unchanged: cells 80798, cu vs record 4.093e-11 / 8.673e-11 / 1.376e-10, σ max 0.999851545099, (c1) 6.422e-04 / 2.812e-04 / 4.863e-04 (`:1783–1788`). Band green. |
| **Verified at** | `main` at `225ed227`, by execution: 62 s, 5 passed per rank. The six W0 windows of the same item are green, 8 passed per rank each (`20261110T150037Z` … `20261110T150135Z_WF-17-s42-W0-*.log`). Wesr, Wyaml, Wtune, Wsar, Wb1, Wsweep and Wprobe were **not run**, because the brief's negative result says stop. |
| **Cause** | Not established. It is the same class as the entries for items 939 and 938 below: round-off of order 1e-16 moving on a green route. The (g) spreads also moved, at ≈ 1e-15 relative; they are printed, not records. |
| **Scope** | Records only. No band is red and no test on `main` is red. `WF-17` stays at 🟡 (held by item 938). |
| **Resolves with** | The same ruling as item 938's entry. Fix nothing in passing. |

### RETIRED 2026-11-10 (10:30 daily review's ruling, `WF-17.md` § "Review ruling (2026-11-10 10:30 daily": the moved value is a band-gated residual, not a record — `log-pathologist`, no asserted or fixed-precision record moved; the step continues as item 943 / 944 / 945) — was 🟡 OPEN 2026-11-10 (`WF-17` step 41 / U150, §9 item 939, 07:30 implementer slot; observed, not fixed): the S8 probes node's (U′) round-off residual prints 1.613e-16 against its prior 2.213e-16 — a record digit off its prior with every band green

| | |
| --- | --- |
| **Symptom** | `test_wf17_four_port_outputs_tutorial.py::test_wf17_s8_probes_block_reproduces_u70` (`-n 2`, complex, `FEM_EM_WF17_S8=1`) passes, but its (U′) line reads `max \|E\|^2 sigma/(2 rho) vs SAR rel 1.613e-16 (ASSERTED <= 1e-10; U70 read 2.213e-16)` (`20261110T133826Z_WF-17-s41-W34.log:3629`) against the prior 2.213e-16 (`20261021T003511Z_WF-17-s8-W34.log.gz:3629`). Everything else the node prints is unchanged: (R) 9.557e-11, the 20 per-port SAR values at `.6e` (`:3625–3628`), (c-disp) 6.3865e-02 / 2.1232e-01 / 8.3482e-02 / 5.6184e-02. Band green (1.613e-16 ≤ 1e-10). |
| **Verified at** | `main` at `808ae4c2` by execution: 77 s, Status 0, 5 passed per rank. Earlier windows of the same item green with every named record at its prior: W0-s7 / s8 / s9 (`20261110T133034Z_WF-17-s41-W0-s7.log` 9, `20261110T133051Z_WF-17-s41-W0-s8.log` 9, `20261110T133059Z_WF-17-s41-W0-s9.log` 8 passed), W29 (`20261110T133121Z_WF-17-s41-W29.log`), W31 (`20261110T133144Z_WF-17-s41-W31.log`), W32 (`20261110T133359Z_WF-17-s41-W32.log`), W33 (`20261110T133621Z_WF-17-s41-W33.log`). W35 (S9) was **not run**: the brief's negative result says stop. |
| **Cause** | Not established; the same class as item 938's entry below (round-off at or below 1e-13 relative moving on a green route). Two further round-off prints, not brief records, also moved: W31's real-build duty probe `values_k` (printed beside only) at ≈ 1e-14 relative (`…W31.log:74`), and W33's builder reciprocity diagnostic `max\|Sij-Sji\|` 7.576e-16 / 1.201e-15 against 3.554e-16 / 3.554e-16 (`…W33.log:3621`). |
| **Scope** | Records only; no band is red, no test on `main` is red. `WF-17` stays at 🟡 (held by item 938). |
| **Resolves with** | The same ruling as item 938's entry: whether a round-off residual print (≈ 1e-16) is a record under U146's rule, or a repeat run on the same tree. Fix nothing in passing. |

### RETIRED 2026-11-10 (10:30 daily review's ruling, `WF-17.md` § "Review ruling (2026-11-10 10:30 daily": the moved value is a band-gated residual, not a record — `log-pathologist`, no asserted or fixed-precision record moved; the step continues as item 943 / 944 / 945) — was 🟡 OPEN 2026-11-10 (`WF-17` step 39 / U148, §9 item 938, 06:00 implementer slot; observed, not fixed): the planted degree-1 control's distance to `MAT4_DEG1_RECORD` prints 6.312e-11 against its prior 6.299e-11 — a record digit off its prior with every band green

| | |
| --- | --- |
| **Symptom** | `test_wf17_first_study_tutorial.py::test_wf17_planted_degree1_block_misses_the_record` (`-n 2`, complex) passes, but its printed `(c4)` line reads `degree-1 SAR = 8.291558128523365e-08 (8.2915581285e-08); … vs MAT4 deg-1 record 6.312e-11` (`20261110T121910Z_WF-17-s39-W-four.log:1919`) against the prior `8.291558128522265e-08 (8.2915581285e-08); … 6.299e-11` (`20261006T004414Z_WF-17-W1-planted-deg1.log.gz:161`). The `.10e` SAR print and the 3.494e-02 miss to the degree-2 record are unchanged; the full-precision SAR moved by ≈ 1.3e-13 relative. `MAT4_DEG1_RECORD` (`tests/package/test_wf15_study_lossy_sphere.py:87`) is unchanged. Every band is green (6.312e-11 ≤ `L0_BAND` 1e-09). The four-port node in the same window reproduces its (R4) prior 9.557e-11 and cells 117206 exactly (`:1842`). |
| **Verified at** | `main` at `4e45beb5`-plus-journal by execution: 53 s, Status 0, 6 passed per rank. W0-a (`20261110T121824Z_WF-17-s39-W0a.log`, 11 passed) and W0-b (`20261110T121836Z_WF-17-s39-W0b.log`, 26 passed / 2 skipped) green. W-first (the degree-2 node) was **not run**: the brief's negative result says stop. |
| **Cause** | Not established. The item's brief reads any record off its prior as "a merge moved a printed digit". A related observation from the same slot: `WF-15` step 19's W-ctl (`20261110T121308Z_WF-15-s19-W-ctl.log`) printed the sphere duty-cycle (Δ) 6.466e-16 against 6.608e-16 and (R-A_t) 3.4708293772e-02 against …773e-02, on a path the slot's edit cannot reach — so round-off at the 1e-13 level on the sphere solve has moved since those priors were logged. Whether this is a merge to the solver path or run-to-run nondeterminism at fixed width is not measured. |
| **Scope** | Records only; no band is red, no test on `main` is red. `WF-17` is held at 🟡 per the brief's rule (U146's rule, item 928). |
| **Resolves with** | The review's ruling: either the record is re-read as the `.10e` print (then step 39 re-runs W-first and closes), or a repeat run of the same node on the same tree shows whether the digit is deterministic. Fix nothing in passing. |

---

## Batch moved 2026-11-14

### RETIRED 2026-11-14 (`ANS-11` C37, §9 item 983, 12:00 CST implementer slot — `test_c19_absent` retired in favour of `test_c37_custody`; the default render reads entry 62's tracked log, never the host JSON) — was 🟡 OPEN 2026-11-14 (`ANS-13` C54, 06:00 implementer slot, §9 item 959; found, not fixed): since `xxl` entry 62 wrote its degree-2 JSONs, the `ans:11` page's default render carries degree-2 cells and `test_c19_absent` is red

| | |
| --- | --- |
| **Symptom** | `tests/unit/test_ans11_record_page.py::test_c19_absent` (`:191–200`) is red. It asserts that `read_c19()` returns `None` for both of entry 62's JSONs, and both are now present (`20261114T121512Z_ANS-13-C54-W1.log:139`). The run reads 1 failed / 84 passed / 1 skipped per rank, and every other node is green. `./scripts/run_examples.sh -e ans:11 -n 1` still prints every (K) `[PASS]`. It then writes a `COMPARISON.md` with 16 changed lines, the degree-2 cells filled, while `metrics.json` stays byte-identical (`20261114T121606Z_ANS-13-C54-ans11-page.log:117`, `read_c19: …/geo38/c6-node-185.7629-d2.json present; …/geo36/i36-node-185.7629-d2.json present`). The slot restored the file with `git checkout`, so `main`'s copy is unchanged. |
| **Verified at** | The parked branch `attempt/ANS-13-C54-20261114T1217Z` (`db34a506`), which is `main` `04387052` plus item 959's thermal case. The default EM path is untouched by that diff; the change in behaviour comes from the gitignored `/workspace/logs/` files, not the tree. |
| **Cause** | `xxl` entry 62 (`3185acfd`, `20261114T080009Z_GEO-38-C19.log`, Status 0) wrote `logs/geo38/c6-node-185.7629-d2.json` and `logs/geo36/i36-node-185.7629-d2.json`. `test_c19_absent` is the pre-window premise node, and `read_c19` is the page's reader for C37 (entry 62's degree-2 rows). Both were built to change behaviour once the files exist. No code defect. |
| **Scope** | Any window that runs `test_ans11_record_page.py`, or that runs `-e ans:11` and then expects `git diff` on `ans:11`'s `COMPARISON.md` to be empty. This blocks item 959's (K) gate as written. It may reach item 979 (C33 loads the `ans:11` page) and `mri:41` (which loads C33); neither was run this slot. A slot that runs `-e ans:11` must restore `COMPARISON.md` unless it is landing C37. |
| **Resolves with** | C37 (entry 62's degree-2 rows on the `ans:11` page), which retires or re-poses `test_c19_absent` and lands the degree-2 `COMPARISON.md`. After that, item 959 re-poses with (K) against the C37-committed files. Ordering is the review's call. **Ruled 2026-11-14 10:30 review:** C37 is §9 item 983. It moves the default render onto entry 62's tracked log, so a clone regenerates the same file, and keeps the host JSONs as custody (J62) only. It retires `test_c19_absent` in favour of `test_c37_custody`. C54 re-poses as item 984, in a later slot than 983 (`ANS-11.md` / `ANS-13.md` § "§9 item 983" / "§9 item 984"). |
| **Retirement (2026-11-14, §9 item 983)** | `render_comparison(m, None)` now fills the degree-2 cells from `m["degree2_c19"]`, the block `degree2_block()` parses from `20261114T080009Z_GEO-38-C19.log`; `main()` no longer passes `read_c19()` to the render, and runs (J62) on it instead. W1 `20261114T181121Z_ANS-11-C37-W1.log` (`-n 2`, ans:11 / ans:12 / ans:13 page tests + docrefs): 76 passed per rank, 0 failed, `[capture] rc=0`. Closing page `20261114T181245Z_ANS-11-C37-page-r.log`: `git diff` on `ans:11`'s two files empty against window 2's (`20261114T181052Z_ANS-11-C37-page.log`). ANS-12's `COMPARISON.md` unchanged after `20261114T181222Z_ANS-11-C37-ans12-page.log`. Item 984 is unblocked. |

### RETIRED 2026-11-14 (`GEO-36` I85c, §9 item 978, 07:30 CST implementer slot — the literal re-recorded at 20 with `FEM_EM_GEO38_C21_DEGREE` asserted as a member; W1 / W2 green) — was 🟡 OPEN 2026-11-13 (`GEO-36` I87′, 04:30 implementer slot, §9 item 968; found, not fixed): I85's knob-parsing test pins the foreign-guard count at 19, and C52 made it 20

| | |
| --- | --- |
| **Symptom** | `tests/validation/test_geo36_rod_f_clinical_3t_heating.py::test_i85_k1_12_knob_parsing` is red at its `:206`: `assert len(foreign) == len(set(foreign)) == 19` reads `assert 20 == 19` on both ranks (`20261113T103903Z_GEO-36-I87p-Wbase.log:58`, `:74`). Every other node in the module is green in item 968's Wreg, including K12-R ("(T-rep) I32b's literals missing none (22 lines; ASSERTED)", `20261113T103546Z_GEO-36-I87p-Wreg.log:478`). |
| **Verified at** | Clean `main` at `cebdb8fc`, by execution: `-n 2` complex, that node alone, 3 s, Status 1 (Wbase). The same red shows in Wreg with item 968's edit applied (1 failed / 22 passed / 10 skipped per rank, `:526`, `:542`). |
| **Cause** | I85's `_foreign_envs()` (`:106–108`) is C32's `_foreign_envs()` plus `HEAT_ENV` and `S11_DRIVE_ENV`. C52 (`7ffc523d`, §9 item 929, 2026-11-12) added `FEM_EM_GEO38_C21_DEGREE` to C32's list in `test_geo36_rod_f_clinical_3t_drive.py:62`, after I85 (`f7e6d193`, 11-10) recorded 19. Item 929's windows did not run I85's module. The count literal is stale, and no guard has lost a name. |
| **Scope** | That node only. I85's gated test (entry 67's knob) and K12-R are unaffected. Any window that runs I85's module reads 1 failed per rank until this is resolved. |
| **Resolves with** | The review's ruling: re-record the literal at 20 in its own step, or derive it from C32's list. It is not fixed in passing, because item 968's scope is one other module. **Ruled 2026-11-13 10:30 daily** (`GEO-36.md` § "Review note (2026-11-13 10:30 daily)"): re-recorded at 20 with `FEM_EM_GEO38_C21_DEGREE` named in an assertion, not derived — §9 item 978, takeable once `xxl` entry 62 has run (I85's module imports entry 62's `test_geo36_rod_f_clinical_drive.py`). This entry leaves with that item's commit. |
| **Retirement (2026-11-14, §9 item 978)** | `:206` reads `== 20`, and a new line asserts `"FEM_EM_GEO38_C21_DEGREE" in foreign`. Nothing is derived from C32's list, and the uniqueness check and the four named members are unchanged. W1 (the node alone, `-n 2`): 1 passed per rank, Status 0, 5 s (`20261114T133050Z_GEO-36-I85c-W1.log`). The K1-12 print lists 20 names. Set against the 19 in `20261110T105744Z_GEO-36-I85-W0b.log:51`, the difference by string is exactly `{FEM_EM_GEO38_C21_DEGREE}`. W2 (`test_complex_mode.py` plus the whole module, `-n 2`): 9 passed / 1 skipped per rank (the knob-gated body), 0 failed, Status 0, 28 s (`20261114T133105Z_GEO-36-I85c-W2.log:203`, `:206`). The control is the logged `assert 20 == 19` at `cebdb8fc` (`20261113T103903Z_GEO-36-I87p-Wbase.log:58`, `:74`). |

### RETIRED 2026-11-14 (`GEO-38` C33c, §9 item 979, 07:30 CST implementer slot — `:265` formats dz `+.6f`; (J58) 0 mismatches, P89′ 0 string-only) — was 🟡 OPEN 2026-11-12 (`MAT-21` step 49 / P89, 06:00 implementer slot, §9 item 963): `c33.custody_json` reports a false mismatch on any c\* with positive dz — its `c_star.dz` format drops the sign the C6a log prints

| | |
| --- | --- |
| **Symptom** | `tests/validation/test_mat21_lattice_ball_superposition_closed_form.py::test_k_and_ks` (on `attempt/MAT-21-P89-20261112T1210Z`) is red: `c33.custody_json(p58, HOST_JSON_DIR / "c6-node-44.0525.json")` returns 76 fields, 1 mismatch, `('c_star.dz', '0.050000', '+0.050000')` (`20261112T120951Z_MAT-21-P89-W1.log:54`, `:99`). The 18 lattice fields and every other field match. |
| **Verified at** | `main` at `ecfb0953` plus the parked module, by execution, `-n 1` complex, 7 s. |
| **Cause** | `tests/validation/test_geo38_c33_field_strength_reader.py` `custody_json` builds `f"{j['c_star'][2]:.6f}"`, while the C6a log prints dz signed (`20261109T080008Z_GEO-38-C6-3T.log:12256`: "dz +0.050000 m"), and the page parser keeps the token as printed. Entry 52's c\* has dz −0.05, whose `.6f` carries its sign, so C33 — which ran custody on entry 52 only (`_c33_body`) — never met a positive dz. A string-format defect, not a value mismatch. |
| **Scope** | Any caller of `c33.custody_json` on a side whose c\* has dz > 0 (entry 58, the 3 T node JSON). C33's own gated test is unaffected. |
| **Resolves with** | The review's ruling on §9 item 963: a fix of `custody_json` in its own step (format `{:+.6f}` or a float comparison of dz), or P89's custody compared without the sign. Nothing to widen. **Ruled 2026-11-12 10:30 daily** (`MAT-21.md` § "Review ruling 2026-11-12 10:30"): P89′ (§9 item 972) splits the helper's mismatches into string-only and value mismatches in its own module and asserts that there are no value mismatches. That clears P89's red but **does not fix the helper**, so this entry stays open until a step fixes `custody_json` itself. Seven modules import it, including the `ans:12` / `ans:13` pages and `mri:41`, so that step re-runs C33's own window and is not queued yet. **Queued 2026-11-13 18:00 daily as §9 item 979 (`GEO-38` C33c, `GEO-38.md` § "§9 item 979"):** by reading, only C33's own module and P89′'s call `custody_json` (the others import the module for `_load_page`), so the fix is `:265` as `{:+.6f}` plus (J58) / (N-sign), run in C33's window at knob `1` and P89′'s window, after `xxl` 62. |
| **Retirement (2026-11-14, §9 item 979)** | `custody_json` now builds `f"{j['c_star'][2]:+.6f}"`. **W1** `20261114T133421Z_GEO-38-C33c-W1.log` (C33's module at knob `1`, `-n 2`): 15 passed per rank, Status 0. (J58) reads `76 fields, mismatches []` (`:176`); (N-sign), a dz-negated copy, reads exactly `[('c_star.dz', '-0.050000', '+0.050000')]` (`:182`); (J52) is unchanged at 0 (`:196`). **W2** `20261114T133435Z_MAT-21-P89-W1c.log` (P89′'s W1b command, `-n 1`): 9 passed, and both (K) custody lines read `0 string-only + 0 value` (`:52`, `:54`). |

### RETIRED 2026-11-14 (`GEO-34` step 1d / U12′, §9 item 985, 12:00 CST implementer slot — the "Resolves with" clause's additive read option landed as `"highest_dim_only": false`; key-absent `:164` call unchanged; W1 / W2a / W2b green) — was 🟡 OPEN 2026-10-08 (`GEO-34` step 1d probe U12a, 16:30 implementer slot, §9 item 333): the importer's default `importShapes` read silently drops a free face that sits beside a solid in the same STEP file

| | |
| --- | --- |
| **Symptom** | `tests/mesh/probe_geo34_step1d_export_forms.py`, form F-comp: a unit box plus a free 0.8 × 0.6 rectangle written to one STEP file. Read with the default `importShapes(path)` (signature `importShapes(fileName, highestDimOnly=True, format='')`), it returns only `[(3, 1)]`: 6 faces and no 0.48 m² face. Read with `highestDimOnly=False`, it returns 7 faces including the 0.48 face. Two separate files merged in one session keep the face under both reads, because `highestDimOnly` applies per file (`20261008T214725Z_GEO-34-step1d-probe.log`, with a bit-identical repeat). |
| **Verified at** | `main` after item 330, by execution: one rank, 5 s, Status 0, twice. |
| **Cause** | Gmsh's documented default: `highestDimOnly=True` keeps only each file's highest-dimension entities. |
| **Scope** | No test on `main` is red. It affects `src/fem_em_solver/io/step_import.py:164` (`_build` calls `importShapes` with the defaults) and `scripts/maintenance/install_cadquery_image.sh:111–120` (step 4b's `len(faces) >= 7` read uses the default). Any single STEP that carries a sheet as a free face beside a solid would lose the sheet silently on import. A face embedded in the solid's topology is not affected by this read. |
| **Resolves with** | U12′'s importer change: an additive read option, or a STEP whose sheet is embedded in the solid. The weekly decides (ruling Q4, 2026-10-08 10:30 review). Do not change `:164` in passing. |
| **2026-11-14 (interactive)** | The operator's CadQuery image install hit exactly this at step 4b (6 faces). The script's own read now passes `highestDimOnly=False` (its separate `-noenv` bug fixed too). On CadQuery's `cq_sample` assembly STEP that read gives 1 volume, 7 faces, and **both product names survive** (`…/AirBox/AirBox`, `…/PortSheet/PortSheet`): the named-sheet route U12′ needs exists. 4c (versions, 180 packages, both DolfinX builds) and 4d (`20261114T093413Z_GEO-34-cadquery-image-4d.log`, 5 passed, Status 0, 104 s) were run by hand on the new image. `step_import.py:164` is unchanged. |
| **Retired 2026-11-14 (§9 item 985)** | W0b (`20261114T181729Z_GEO-34-U12p-W0b-read.log`) re-read this on the committed CadQuery fixture: default read six unnamed faces, no `PortSheet` (`:60–67`); `highestDimOnly=False` keeps face 7 `Shapes/geo34_step1d/PortSheet/PortSheet` (`:136`). The importer's new additive key `"highest_dim_only": false` passes the option; with the key absent the default call is unchanged, and the (n-hDO) control asserts the default read still raises naming `PortSheet` (`20261114T181927Z_GEO-34-U12p-W1.log`). `install_cadquery_image.sh` was already fixed interactively. |

## ✅ RETIRED 2026-11-15 (§9 item 997, 19:30 CST implementer slot) — was 🟡 OPEN 2026-11-15 (`WF-17` U163, §9 item 992, 13:30 CST implementer slot) — the 1.5 T estimate-page custody nodes raise `FileNotFoundError`: the weekly sweep gzipped two ledger-cited `xl` logs they read as plain text

- **Test ids:** `tests/package/test_wf17_estimate_measured_15t.py::test_wf17_s52_ledger_strings_and_custody` and `::test_wf17_s52_c52_custody_fails_on_entry52`, red on both ranks at `-n 2`.
- **Symptom:** `FileNotFoundError: [Errno 2] No such file or directory: '/workspace/docs/testing/logs/20261102T080008Z_GEO-36-step9-I36.log'` (`20261115T193357Z_WF-17-s53-W0.log:304–305`). The rest of the module is green, including (K15) / (Sec15) / (Order) (`:133`).
- **Cause:** the 2026-11-15 weekly sweep (`f7f8f460`, 01:45 CST) gzipped `20261102T080008Z_GEO-36-step9-I36.log` (entry 48) and `20261103T080009Z_GEO-38-C6.log` (entry 52) under the 7-day retention rule. `_custody` in `tests/package/test_wf17_estimate_measured_3t.py:58–76` opened `LOGS / p["log"]` as plain text with no `.gz` fallback. `housekeeping.py:138–143` gzips every log at 7 days, cited or not; citation exempts a log from deletion only.
- **Wider scope (18:00 daily review, code sweep at `7db96c4c`):** `_parse_c52` (`test_mat21_eddy_source_probe_points_closed_form.py:182`) read the same gzipped C6 log as plain text (four p80 / p81 nodes inferred red); eight more nodes (3 T `_custody` users, six `test_s50_*` closed-form page nodes, the lattice `_c39`) were due red at the 2026-11-22 sweep.
- **Resolution (§9 item 997):** `_custody` (new `logs=LOGS` default), `_parse_c52` (iterating `io.StringIO`), the closed-form page's `read_one` / `is_file` / two `splitlines` reads and the lattice `_c39` now read through `tests/log_records.read_log_text(errors="strict")` / `resolve_log_path`. The pre-edit W0 (`20261116T013050Z_WF-17-i997-W0.log`) measured 6 failed / 41 passed per rank: the two 15t nodes **plus the four p80 / p81 nodes the review inferred, found red and fixed here**, all `FileNotFoundError`. W1 (`20261116T013304Z_WF-17-i997-W1.log`) 55 passed per rank (110 `PASSED`, 0 `FAILED`), 16 s, with the new `tests/package/test_wf17_log_gz_readers.py`: opposite-form copies of all ten table logs read text-identical to the tracked form and to `Path.read_text` semantics, custody `[]` in both forms (1.5 T and 3 T), the other readers parse identically, locale `utf-8`, and the absent / truncated (`EOFError`) / CRC-flipped (`BadGzipFile`) controls raise. Exempting cited logs from gzip remains the weekly's call.
