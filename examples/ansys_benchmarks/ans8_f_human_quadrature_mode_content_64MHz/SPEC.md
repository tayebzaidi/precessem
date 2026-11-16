# ANS-8 — The tuned F-human birdcage, two-port 90° quadrature: is the non-k = 1 content seen by AED too?

**Commissioned 2026-10-01 by operator directive (interactive session), the
SPEC written in that session acting for the weekly at the operator's request**
— PROJECT_PLAN §7 `ANS-8` and §10 chain H (H8b ruling 2026-10-01). Seventh
Ansys benchmark case and the first at human scale. Step 1 (this SPEC) is the
mode content of the tuned coil under the clinical drive; step 2 (tuned
`|B₁⁺|` and the 10 g SAR hotspot) is reserved until H9 lands.

**Why this case.** Our F-human birdcage, tuned by the k = 1 block with a
capacitor on every gap and driven at two gaps 90° apart (the drive the
operator's clinical coils use), puts only **≈ 0.79** of the feed-ring gap
current in the k = ±1 (imaging) mode, with a **≈ 0.17 k = 0** leak; the other
ring reads ≈ 0.95 (`WF-10` step 1d, item 133,
`20261001T140319Z_WF-10-step1d-W1.log`). Our circuit layer predicts ≈ 0.99 for
an opposite-pair or four-port drive, so the leak is attributed to the
two-port geometry. Some non-k = 1 content is physical for this drive. The
question for AED is whether **an independent solver sees the same mode
content on the same setup** — if it does, the ≈ 0.79 / 0.17 is the drive's
physics; if it does not, the defect is in our fixture or port model.

## Status

- [x] SPEC v1 (2026-10-01).
- [x] **SPEC v2 (2026-10-01, operator directive, interactive session)** — the
  own-`C_mode` row (2b), the symmetric-drive controls, row 4 redefined, the
  First Order and copper columns commissioned. See "What v2 changes".
- [x] Runnable half (`ans:8`, our side — see "Our side" below; landed 2026-10-01, §9 item 141, against v1's rows).
  **v2 rows 2b / 2c / 4a landed 2026-10-02** (§9 item 161, `ring_mode_rows_v2`, no new solve).
- [x] Runnable half, v3 rows (2026-10-11, §9 item 386; item 369 was attempt 1).
- [x] **SPEC v4 (2026-10-18, §9 item 513, §10 H29)** — step 2's capacitor at the degree-2 root, C = 115.5209 pF; the degree-1 values are a record.
- [x] Operator replication in AED, **v1 column: Zero Order, σ = 800**
  (2026-10-01; private results in `aed_results/`).
- [ ] Operator replication in AED, **v2 columns**: copper Zero Order, σ = 800
  First Order, copper First Order (started 2026-10-01).
- [ ] Operator replication in AED, **v4 column** (C = 115.5209 pF) (step 2, one session — see "Operator session brief")
- [ ] Adjudication (the weekly after both halves land).

## What v2 changes, and why

The Zero Order σ = 800 AED half showed two defects of v1's *design*, both
visible without quoting an AED value:

1. **Row 2 at the fixed C confounds mode content with detuning.** v1 put C
   at *our* `C_mode` on both codes. Wherever the codes' `C_mode` differ, the
   other code is driven off its own k = 1 resonance. The two-port share and
   leak move steeply with that detuning: on the AED sweep the fixed-C share
   changes by several percent per MHz near 64 MHz. So a fixed-C disagreement
   cannot separate "different mode content" from "different `C_mode`". The
   readout v1 pre-registered compares exactly that quantity. **v2 adds row
   2b: every drive evaluated with each code's matrix terminated at its *own*
   `C_mode`**, so both codes sit on their own k = 1 resonance. Row 2 is kept
   as a record of the detuning sensitivity.
2. **v1's row 4 has no root to find.** The − → + zero of `Im Y_in(P17)` with
   C across P17 is not the k = 1 resonance: the input admittance is
   dominated by P17's own capacitor (`ωC ≈ 45 mS` at 64 MHz), and on the AED
   sweep `Im Y_in` stays positive across 60–68 MHz. **v2 row 4 is the k = 1
   block resonance**, the f where `C_mode(f) = C`.

**Blindness, stated honestly.** Row 2b and the controls were written *after*
the Zero Order σ = 800 AED column was seen. For that column they are post
hoc. They are pre-registered — blind — for the three v2 columns (copper
Zero Order, σ = 800 First Order, copper First Order) and for our side's
runnable half. The weekly weighs the v1 column with that caveat.

## Geometry (SI units, metres)

The F-human fixture is `birdcage_port_domain` with `_params(0.15,
scale_sizing=False)` (`tests/validation/test_birdcage_f_human_rung.py`), the
four-leg F-small layout scaled by **s = 0.15 / 0.07 = 2.142857142857**, the
conductor cross-sections *not* scaled. All values below are that call's,
evaluated.

| Part | Definition |
|---|---|
| Air box | `x, y ∈ [−0.2374285714, +0.2374285714]`, `z ∈ [−0.2142857143, +0.2142857143]` (radial extent = R + r_leg + 0.0171428571 + 0.0642857143; z extent = coil_length / 2 + 0.0642857143). |
| Legs | **16** solid cylinders, radius **r_leg = 0.006**, axis ∥ ẑ at radius **R = 0.15**, azimuths **θ_n = n · 22.5°**, n = 0…15, spanning **z ∈ [−0.15, +0.15]** (unbroken — this is a high-pass coil). |
| End rings | Two tori, major radius **0.15**, minor radius **0.004**, centred on the axis at **z = ±0.1178571429** (leg spacing 0.2357142857). The legs pass through the rings and extend 0.028 beyond them; the conductor is the **union** of legs and rings. |
| Ring gaps | On **both** rings, at the 16 mid-azimuths **φ_j = 11.25° + j · 22.5°**, j = 0…15, remove the arc `\|φ − φ_j\| ≤ α`, **α = g / (2R) = 0.0571428571 rad** (g = 0.0171428571), with **radial** cut planes at `φ_j ± α`: **32 gaps**, each leaving two planar circular terminal faces of radius 0.004. |
| Gap sheets | In each gap one planar rectangle in the **tangent plane at radius R** (normal r̂(φ_j)), centred at (R cos φ_j, R sin φ_j, z_ring): **2R·tan α = 0.0171615404** along φ̂(φ_j) × **0.010** along ẑ. Its two φ̂-edges lie in the radial cut planes and cross both terminal faces. Current / integration line along **+φ̂**, from the cut plane at φ_j − α to the one at φ_j + α, at radius R, z = z_ring. |
| Phantom | Solid cylinder, radius **0.0642857143**, **z ∈ [−0.0857142857, +0.0857142857]**, on the axis. |

**Port numbering (ours, used in every row):** gap sheet (ring, j) is port
**P(17 + j)** on the **bottom** ring (z = −0.1178571429) and **P(33 + j)** on
the top ring (z = +0.1178571429), j = 0…15. **Feeds: P17 (φ = 11.25°) and P21
(φ = 101.25°)**, both on the bottom ring, 4 slots = 90° apart.

*Not modelled, by design:* the 16 floating leg-port boxes our generator still
emits in the isocenter plane (`GEO-35` — inert vacuum, no port or material);
the AED model omits them.

## Materials

| Region | σ [S/m] | εᵣ | μᵣ |
|---|---|---|---|
| Air box (everything not listed) | 0 | 1 | 1 |
| Coil (legs ∪ rings) | **800**, **solve inside** (volume conductor — exactly our fixture, as `ANS-4` did) | 1 | 1 |
| Phantom | **0.5** | **78** | 1 |

All six box faces **Perfect E**. **Copper column (commissioned v2):** the coil
as **copper**, Finite Conductivity σ = 5.8e7, Solve Inside off (`ANS-6`'s
treatment). Our copper side follows `WF-12` (step 1 BLOCKED 2026-10-01), so
this column is reference data until then. It is still useful now: the share,
leak and `C_mode` move from σ = 800 to copper *within AED* say how much of the
σ = 800 fixture's mode content is the volume conductor.

## Excitation — Route B (primary): 32 ports + NumPy

All 32 gap sheets are **50 Ω lumped ports** (no de-embedding, renormalised to
50 Ω), the S / Z matrix exported at full precision. Everything below is
circuit algebra on that matrix, run identically on both codes' 32×32:

- `Z = 50 (I + S)(I − S)⁻¹`, `Y_net = Z⁻¹`.
- **Capacitor on every gap: C = 113.11 pF exactly (1.1311e-10 F)**, ideal —
  our H8b `C_mode` (1.131080769863e-10 F, item 133) rounded; fixed, not tuned
  on the AED side.
- **Terminations:** every non-feed port j: admittance `jωC`. Each feed f
  (P17, P21): a Norton source `I_s,f = V_f / 50` in parallel with
  `1/50 + jωC` (the drive across its gap's capacitor through 50 Ω — the
  parallel drive of `ANS-7` v2).
- Solve `(Y_net + Y_load) V = I_s`; the **gap currents** are the network port
  currents `I = Y_net V` (the current crossing each gap through whatever
  terminates it).
- **Drives:** (i) **P17 alone** (`V_17 = 1`, `V_21 = 0`); (ii) **ccw
  quadrature** `V_17 = 1`, `V_21 = −j` (the port at +90° azimuth lags —
  our `quadrature_phase_weights(…, "ccw")`, the sense that co-rotates with
  `B₁⁺`); (iii) **cw**, `V_21 = +j`.

**Route A (optional, if time allows):** the same model with **Lumped RLC**
(capacitance only, C, written **with its unit**, `113.11pF`) on the 30
non-feed gaps and the two feeds as 50 Ω ports, each with C added across it
in post exactly as `ANS-7` Route A; its P17 / P21 2×2 checks Route B.

## Frequency and mesh

- **Adaptive at 64 MHz** to **ΔS ≤ 0.005** (32 ports make the 0.002 of
  earlier cases expensive), **≥ 2 consecutive converged passes**, max 10
  passes, ≥ 2 elements across each gap along φ̂; report if the cap is hit
  instead of converging (the `ANS-7` Zero Order Route A designs did).
- **Primary rows at 64 MHz** on the adaptive mesh.
- **Secondary (AED only for now):** a **discrete** sweep 60 → 68 MHz in
  0.5 MHz steps (17 points, every point solved) on the 64 MHz mesh, to locate
  the terminated network's resonance under drive (i). Our side's sweep comes
  later on a licensed sparse grid (as `ANS-7` step 3).
- **Basis:** **Zero Order** (= our degree 1) primary. **First Order** (= our
  degree 2) **commissioned v2**, both columns, per the `ANS-5` order
  correspondence. Mixed Order forbidden. Confirm unknowns / tet from the
  matrix statistics. First Order answers whether the `C_mode` gap (row 1) is
  a degree-1 effect on our side: the `ANS-7` F-small finding (our `f₀`
  1.2–1.9 % high at degree 1, §10 chain T) predicts the AED First Order
  `C_mode` moves *away* from, or stays at, the Zero Order one. It never moves
  toward ours by construction; that is the question.
- **The sweep runs on Zero Order only** (v2): 17 discrete points on the
  adaptive mesh are a large share of a Zero Order run, and First Order
  carries ≈ 6× the unknowns per tet.
  Row 4 for the First Order columns comes from the one-point form 4a below,
  which the Zero Order sweep checks.
- Solve resources: **20 cores** (operator 2026-10-01; was 12, 2026-09-29), **1 TB RAM ceiling**. This is the Windows AED machine, not the shared Linux box, whose 12-core cap is unchanged.

## Quantities to export (all digits AED prints; do not round)

| Quantity | Definition |
|---|---|
| `S`, `Z` (32×32) at 64 MHz | Route B, raw, Touchstone `.s32p` + JSON |
| `S`, `Z` on the 17-point sweep | Route B, if the sweep is run |
| Loss partition at 64 MHz, drive (ii) | `P_coil / (P_coil + P_phantom)`, both integrals, with the drive applied by the AED source settings that reproduce drive (ii) (or computed from per-port field solutions superposed). **v2:** also `P_coil` by power balance, `P_net − P_phantom` with `P_net = Σ (\|a\|² − \|b\|²)` over the 32 ports at the same source scaling (PEC box, lossless air). This is the only coil-loss route for copper, where Solve Inside is off, and on σ = 800 the balance residual against the coil integral is an identity check. |
| Solve metadata | tets, passes, final ΔS, converged yes / no, wall time, order, unknowns / tet |

## Primary adjudication rows (identical algebra, both codes, 64 MHz)

1. **`C_mode` from the k = 1 block.** DFT the 32×32 `Z` over the 16 slots
   into 2×2 blocks `Ẑ_k` (bottom / top ring); `λ₁` = the eigenvalue of `Ẑ₁`
   whose eigenvector is ring-antisymmetric (carries leg current);
   `C_mode = 1 / (ω Im λ₁)`. Print the off-block norm (C16 breaking).
2. **k-spectrum of the gap currents** for drives (i), (ii), (iii), **per
   ring**: `ĝ_k = Σ_j I_j e^{−2πi jk/16}` over that ring's 16 ports;
   **k = ±1 share** `(|ĝ₁|² + |ĝ₁₅|²) / Σ_k |ĝ_k|²`, **k = 0 leak**
   `|ĝ₀|² / Σ_k |ĝ_k|²`, and the full `|ĝ_k|²` fractions.
3. **Input admittance at the feeds** under drive (ii): `V_f / I_f` per feed,
   and the 2×2 feed-port matrix with every other gap in C.
2b. **(v2) Rows 2 and 3 at each code's own `C_mode`:** the same drives, every
   gap — the feeds' included — at **that matrix's own** row-1 `C_mode`
   instead of 113.11 pF. This is the mode-content comparison with both coils
   on their own k = 1 resonance.
2c. **(v2) Controls, at both C:** the **opposite pair** (`V_17 = 1`,
   `V_25 = −1`, P25 the 180° gap on the feed ring) and the **four-port ccw**
   (`V_17, V_21, V_25, V_29 = 1, −j, −1, +j`). Both carry no k = 0 component
   by symmetry, so their k = 0 leak is a C16 / D16-breaking gauge. Their
   k = ±1 share (≈ 0.99 on our circuit, item 133) is the drive-independent
   ceiling that the two-port share is read against.
4. *Secondary, redefined v2:* **the k = 1 block resonance** `f₀` — the f
   where `C_mode(f) = 113.11 pF`, linearly interpolated on the sweep.
   **4a:** the one-point form `f₀ ≈ 64 MHz · √(C_mode / C)` (λ₁ taken
   ∝ ω), which is all a single-frequency matrix allows. It is used for the
   First Order columns and for our side; the Zero Order sweep checks it. The
   k = 2 block's `C` (the `TH-17` / H8b negative-control mode) as before.
   *(v1's − → + zero of `Im Y_in(P17)` is kept in the JSON as a diagnostic
   only; see "What v2 changes".)*

**Readout, pre-registered (v2).** The primary reading is **row 2b**: the
feed-ring k = ±1 share and k = 0 leak under drive (ii), each code at its own
`C_mode`, agree at a stated relative level ⇒ the leak is the two-port drive's
physics and H8b's fixture stands. Disagreement there ⇒ our fixture or port
model. The controls (2c) must agree on their ≈ 0.99 ceiling for the reading
to stand. Row 2 (fixed C) is read only as the detuning record: its gap
between the codes should be what row 1's `C_mode` gap and the share's slope
along the AED sweep predict. **Row 1, `C_mode`, and row 4, `f₀`,** are read
beside `ANS-7`'s (our `f₀` 1.2–1.9 % high at degree 1 on F-small), and the
First Order columns say whether the gap is order-driven. No absolute claim
follows from any reading. *(v1's readout, kept for the record: share / leak
at the fixed C agree at a stated level ⇒ physics; disagreement beyond the
`C_mode` gap's own size ⇒ fixture or port model. v2 supersedes it for the
reason in "What v2 changes" (1).)*

## Our side (`ans:8`, the runnable half)

`examples/ansys_benchmarks/ans8_f_human_quadrature_mode_content_64MHz/` —
imports H3's / `WF-10` step 1d's F-human construction (never re-implements
it), solves our 32×32 at 64 MHz at degree 1 with the κ-corrected sheets
(`PORT-14` step 3; κ ≈ 4.59e-3 on this fixture), and runs **the same
algebra** (rows 1–3) through one shared function that the private writer
also applies to the AED matrix, so the two halves cannot drift. Its anchors
reproduce item 133's records (C_mode, the circuit shares 0.790993 / 0.952414
for drive (ii)) at their imported bands. `COMPARISON.md` carries our rows
with the AED columns blank; `COMPARISON_private.md` (gitignored) fills them.

**v2 additions to our side** (a follow-up to §9 item 141, which was written
against v1 and may land first): rows 2b, 2c and 4a are pure algebra on the
same 32×32, so they extend `ring_mode_rows` with no new solve. Our
counterpart to the First Order column is a degree-2 32×32. That is a new
solve, priced by the weekly, and none exists today. Our copper counterpart
is `WF-12` step 3.

## Reference values (ours)

None are written here: the runnable half regenerates them. For orientation
only: `C_mode` ≈ 113.1 pF; drive (ii) feed-ring k = ±1 share ≈ 0.79, k = 0
leak ≈ 0.17, other ring ≈ 0.95; opposite-pair / four-port drives ≈ 0.99
(circuit, printed in item 133's log).

**Sense convention, checked 2026-10-01.** Item 133's printed "two-port
quadrature (0/4)" drove `(V_17, V_21) = (1, +j)` — the +90° port *leads*,
which is **cw** under `quadrature_phase_weights`; this SPEC's drive (ii)
**ccw** is `(1, −j)`. The two senses differ by ≈ 0.002 in the feed-ring
share on our matrix; both are rows, so compare like with like. The AED
script's rows functions, run on our held 32×32
(`logs/wf10/h8a-held-s32.npz`, item 128), reproduce item 133's `C_mode`, the
k = 2 control's C, the C16 off-block norm and the P17-alone shares to every
printed digit — the algebra and the port / ring / phase conventions are
verified end to end before AED is run.

## Out of scope

No matching network; no capacitor ESR; no tuning on the AED side; no
`|B₁⁺|`, SAR or compliance quantity (step 2, after H9); no 32-leg coil
(`GEO-33` / `WF-11`); no copper claim from this step (copper is the
optional reference column).

## Step 2 (v3, commissioned 2026-10-10 by the weekly review; §10 H21 / chain H) — the tuned two-port drive with capacitor ESR: mode content, `|B₁⁺|` and the 10 g hotspot

**Precondition met 2026-10-10:** `WF-10` step 1e found the in-model Larmor
root of the copper F-human with ESR (C_in,F-human = 118.2437 pF, resolution
0.10 pF) and step 2 read `|B₁⁺|` and the 10 g hotspot under the two-port
90° drive at that node (items 354 / 358 / 359). This step is the human-scale
deliverable of §1 — construct → tune → drive → safety quantities — compared
with AED on the operator's existing F-human project. **One AED session.**

**v4 (2026-10-18, §10 H29): the degree-2 root.** At degree 2 (the
production order) the degree-1 node above is detuned. `WF-10`'s H26 found
the in-model root at **C_in(d2) = 115.5209 pF** (P17 (1, +1j) secant,
bracket [115.50, 116.00] pF, **resolution 0.50 pF**, four-reading spread
0.036 pF; `20261017T070008Z_WF-10-step1e-H26.log:8674`). The step-2 fixture
capacitor is that value, because a benchmark column is adjudicated at the
production convention only (the 10-12 ruling).

**Fixture (both codes):** step 1's geometry, materials and 32 ring-gap
sheets, **copper** conductor (Finite Conductivity on the AED side, Leontovich
on ours), **every one of the 32 gaps carries the same capacitor
C = 115.5209 pF with a 0.2 Ω series ESR** (binding ruling 2026-10-06; in
HFSS as the parallel equivalent at 64 MHz: R_p = 1/(ω²C²r) = **2 317 Ω ‖ 115.52 pF**, Q_C = 1/(ωCr) = 108); the two feed gaps at 0° and 90° are
50 Ω lumped ports driven in quadrature (sense: `ccw` = (1, −1j) at slots
0 / 4, the convention of `quadrature_phase_weights`), the other 30 gaps are
the capacitor sheets. Frequency **64 MHz** only. First Order (operator order
ruling); a Zero Order column optional. No matching network. **State in the
results file whether AED was re-tuned**: it must not be — the capacitor is
fixed at our root, so that row 1 measures the root, not AED's own tune.

**Quantities to export (all printed digits), per √W of accepted power
(P_acc = incident − reflected at the two feeds):**

*Accepted power, defined (weekly ruling 2026-10-12, known-issues 2026-10-11):*
P_acc is the sum over the two feed ports of HFSS's **accepted power at the
lumped port** — incident minus reflected at the 50 Ω port plane, the
default per-port "accepted power" HFSS reports. The feed gap's capacitor and
its 0.2 Ω ESR are load, not feed loss: their dissipation is inside P_acc and
appears in row 5's ESR fraction. Our side normalises by the same line cut
(½ Re Σ_feeds V_port · conj(I_line), printed as P_acc,SPEC) and prints beside
it the field's own feed-sheet partition, which exceeds the line cut by
≈ 0.8 % on this mesh at degree 1 — a documented property of the discrete feed
sheet (the Cauchy–Schwarz excess of a field varying across the sheet), not a
quantity the AED column is asked for. Row 2 on the AED side is step 1's row
2b algebra on AED's own 32-port S at C = 115.5209 pF with the two feeds
terminated; on our side it is the field route's gap currents at the root
(the two agree on our fixture to the (Lin) identity, 6.6e-16).

| Row | Quantity |
|---|---|
| 1 | Im Y_in at each feed at 64 MHz with C = 115.5209 pF (sign and value): the root check |
| 2 | feed-ring k = ±1 share and k = 0 leak of the two-port drive (step 1's row 2b algebra on the 32-port S at this C) |
| 3 | `\|B₁⁺\|` at five probes, metres: O (0, 0, z_c); P₁ (0.045, 0, z_c); P₂ (0, 0.045, z_c); P₃ (−0.045, 0, z_c); P₄ (0, −0.045, z_c) — z_c the coil's mid-plane; point values |
| 4 | the 10 g-averaged SAR maximum in the phantom and its location (C95.3 cube or sphere — state which) |
| 5 | power fraction: phantom loss / P_acc, copper loss / P_acc, capacitor (ESR) loss / P_acc |
| 6 | `\|B₁⁺\|` / `\|B₁⁻\|` at O (the drive-sense check) |

**Readout (pre-registered):**
- Row 1: Im Y_in within one grid resolution of zero at both feeds ⇒ the root
  AGREES; a definite sign ⇒ the roots differ, and the gap in C (from AED's
  own Im Y_in slope) is the finding.
- Row 2: k = ±1 share and k = 0 leak within **1 %** (step 1 agreed within
  0.6 %) is AGREE.
- Row 3: within **10 %** per probe, with the P₁ … P₄ ordering preserved, is
  AGREE; O is the headline figure.
- Row 4: within **15 %** on the value and within one ball radius on the
  location is AGREE (a maximum on two meshes; our phantom is 0.015 m,
  ≈ 1.8 cells per ball, and §10 H18 refines it — the SPEC does not wait).
- Row 5: phantom fraction within **5 %**, copper and ESR fractions within
  **10 %**, is AGREE.
- Row 6: a ratio ≫ 1 in both codes (ours 21.5) is the sense check; < 1 on one
  side is a convention finding, not a disagreement.

**Our side:** our adjudicated column is `ans:8`'s v3 rows at **degree 2 at
C_in(d2)**, from `xl` entry 41's degree-2 node record at 115.5209 pF
(§10 H30; entry 41's output, not a file in this case directory). The degree-1 rows at 118.2437 pF (§10 H21, `COMPARISON.md`'s
step-2 section, which labels itself "degree 1 … C_in = 118.2437 pF") are
**a labelled record, not the adjudicated column**. Lossless step-1 columns
stay as they are, labelled lossless.

### Operator session brief (step 2 — one AED session)

Each bullet executes the step-2 text above; it adds no row, readout or value.

- **Project.** Open the existing F-human project and copy the step-1 design
  to a new design. Never overwrite a step-1 design.
- If a design was already built at the v3 capacitor ("Degree-1 record"
  below), copy the step-1 design again; never re-use it for v4.
- **Conductor.** Copper: Finite Conductivity σ = 5.8e7 S/m, Solve Inside off
  ("Materials", copper column; "Fixture (both codes)" above).
- **Gaps.** All 32 gaps carry **C = 115.5209 pF in parallel with R_p = 2 317 Ω**
  (the parallel equivalent of the 0.2 Ω ESR at 64 MHz, "Fixture (both
  codes)"). Enter C with all its digits and its unit, `115.5209pF` (step 1's
  Route A lesson, "Excitation").
  - The 30 non-feed gaps are Lumped RLC sheets.
  - The two feed gaps (P17 at 0°, P21 at 90°) are 50 Ω lumped ports. **The
    same R_p ‖ C must be present across each feed gap in the solved field.**
    Rows 3–5 need it, because the feed capacitor's dissipation is inside
    P_acc ("Accepted power, defined"). A capacitor added "in post", as step
    1's Route A did, does not change the field.
  - How the feed-gap R_p ‖ C is realised is the operator's choice; state it
    in the results file.
- **Drive.** Equal magnitudes, P21 lagging P17 by 90° (`ccw` = (1, −1j),
  "Fixture (both codes)"). Row 6 is the sense check.
- **Solve.** 64 MHz only, First Order. Adaptive as in "Frequency and mesh".
  **No re-tune**; say so in the results file.
- **Export.** Rows 1–6 with all printed digits, per "Accepted power,
  defined", plus the solve metadata ("Quantities to export"). Results go
  only to the gitignored `aed_results/`. No AED number ever enters a tracked
  file.

### Degree-1 record (v3 / v3.1 values — do not enter)

v3's capacitor was C = 118.2437 pF (resolution 0.10 pF), the degree-1
root, with R_p = 2 212 Ω and Q_C = 105 at 64 MHz and r = 0.2 Ω. Our
degree-1 column at that node (`COMPARISON.md` step 2) is kept as a record;
no AED column is commissioned at it. (v4 enters C in the high-Q convention
v3 used; the exact parallel C_p = C·Q²/(1 + Q²) is 115.5109 pF, 0.0100 pF
below, 0.02 of the 0.50 pF bracket.)

## Revision history

- **v3 (2026-10-10, weekly review):** step 2 as above — the tuned two-port
  drive with 0.2 Ω ESR on every gap at our in-model root; rows 1–6 and their
  readout pre-registered; one AED session on the existing project.
- **v3.1 (2026-10-13, §9 item 411):** operator session brief added; no row, readout, fixture value or definition changed.
- **v4 (2026-10-18, §9 item 513, §10 H29):** step 2's capacitor moves to the degree-2 root C_in(d2) = 115.5209 pF (R_p ‖ C = 2 317 Ω ‖ 115.52 pF, Q_C = 108, high-Q form, 64 MHz); the operator brief re-pointed; the degree-1 values kept as a labelled record; no row, readout or definition changed.

- **v1 (2026-10-01, operator directive, interactive session):** step 1 as
  above.
- **v2 (2026-10-01, operator directive, interactive session, after the Zero
  Order σ = 800 AED column):** row 2b (own `C_mode`) becomes the primary
  readout and row 2 becomes the detuning record; controls 2c added; row 4
  redefined as the k = 1 block resonance with the one-point form 4a; the
  `P_coil` power-balance route added; the copper and First Order columns
  commissioned; the sweep restricted to Zero Order. The AED script carries
  all of it (`aed/`, gitignored). Blindness caveat in "What v2 changes".
