# ANS-13 — The ASTM-F2182-box-loaded F-clinical birdcages' gel heating at 64 and 128 MHz: HFSS SAR as the source of a one-way transient thermal solve (Icepak), per watt deposited in the box

**Commissioned 2026-11-11 by the weekly planning review (PROJECT_PLAN §10
chain P / chain C — the thermal `ANS` case posed by `xl` entries 59 and 60's
pre-registered decision rules, `docs/testing/xl-pending.md` entries 59 / 60).**
Phase 8's second external check and its first at clinical scale: `ANS-9`
(✅ 2026-10-14) checked the EM → thermal hand-off on F-small's small cylinder
at 10 MHz; this case checks the same route on the production fixture, the
standard's gel box and both Larmor frequencies. **It is cheap for the
operator because the EM half is byte-for-byte `ANS-11` (1.5 T) and `ANS-12`
(3 T)** — the projects W10 / W11 already ask for — plus one thermal design
each, as `ANS-9` added one to `ANS-2`.

Our degree-1 columns exist (`xl` entries 59 and 60, records, below); **the
adjudication columns are the degree-2 ones** (`xl-pending.md` entries 69 and
70, PENDING; the 2026-10-12 epitaph: a benchmark column is recorded at the
production convention before it is adjudicated). Every value below is taken
from the code at the cited line or from the two windows' tracked logs.

## Open before the operator starts

1. **The 1.5 T half is ready to build once `ANS-11` is built:** the node is
   `ANS-11` v2's C = 185.7629 pF (R_p = 896.1 Ω ‖ C at 64 MHz). Nothing
   thermal depends on the node.
2. **The 3 T half waits on `ANS-12`'s node (SPEC v2)** — the 3 T degree-2
   root, `xl` entry 64 (Sunday 2026-11-15) then its refinement. Build the
   thermal design on the `ANS-12` project now; do not solve it before v2.
3. **The gel's thermal constants are the same at both field strengths and
   its 128 MHz σ / ε_r are provisional** (dashboard W3; `ANS-12` item 4).
   The 3 T column is re-posed if the F2182 clause moves them.
4. **The per-watt convention is θ per watt deposited in the box.** Our
   windows print θ per *accepted* watt (P_acc, the two feeds' lumped line
   cut); the two differ by the factor P_box / P_acc, which the EM design
   gives and `ANS-11` / `ANS-12` row 5 already adjudicate. Fixing the
   deposited power at 1 W, as `ANS-9` did, keeps the feed systematic out of
   the thermal comparison. `CEM43` is **not** exported per watt: at both
   field strengths our per-watt `cem43` sits at the 20 °C baseline
   (2821.5 × 0.25²³ ≈ 4.01e-11 s; entries 59 / 60's ledger rows), so a
   per-watt dose is the 1 W history's dose and scales with nothing.

## Status

- [x] SPEC v1 (2026-11-11, skeleton; the 1.5 T half final, the 3 T node TBD
      with `ANS-12`).
- [ ] Runnable half (`ans:13`): a record page that imports P61's and C31's
      modules (`tests/validation/test_geo38_f_clinical_box_heating.py`,
      `tests/validation/test_geo38_f_clinical_3t_box_heating.py`) and C46's
      reader, never re-implements them; `COMPARISON.md` with our rows, AED
      columns blank (§10 C53 / C54).
- [x] Our degree-1 columns (`xl` entry 59, 2026-11-10; entry 60, 2026-11-11)
      — **records** (below).
- [ ] Our degree-2 columns (entries 69 / 70, `xxl`) — **the adjudication
      columns**; their body is §10 P87 (the degree knob on P61's / C31's
      modules).
- [ ] Operator replication in AED (HFSS + IcepakFEA), First Order EM.
- [ ] Adjudication (weekly, on entries 69 / 70).

## EM half — `ANS-11` and `ANS-12`, unchanged

Geometry, materials, the 64 ESR capacitor sheets, the two 50 Ω feeds P33 /
P41, the drive, frequency, order and convergence settings: **exactly
`ANS-11`'s SPEC (64 MHz) and `ANS-12`'s SPEC (128 MHz).** If this file and
either of those ever disagree on the EM half, theirs is authoritative. Use
the **First Order** designs.

**Excitation for the thermal source: drive (ii), ccw = (V_33, V_41) =
(1, −1j)** — the sense that co-rotates with `B₁⁺` and the one both windows
report first — with the source **scaled so that the power absorbed in the
gel box is exactly 1 W** (`∫_box ½ σ |E|² dV = 1 W`; in HFSS: edit sources so
the volume-loss integral over the box reads 1 W, and record the accepted
power P_acc that gives it). The cw sense (1, +1j) is a second column at c\*
only, as in `ANS-11` row 3.

## Thermal half

| Item | Value |
|---|---|
| Thermal domain | **the gel box only** — 0.42 (x) × 0.09 (y) × 0.65 (z) m, centred on the origin (`tests/validation/test_mat21_box_closed_form.py:34`; V = 0.02457 m³, Robin area 0.7386 m²; entry 59 `:12315`) — no air, no coil, no container wall |
| Thermal conductivity k | **0.6 W/(m·K)**, uniform, temperature-independent (`K_GEL`, `tests/validation/test_mat15_coil_driven_heating_phantom.py:82`) |
| Volumetric heat capacity ρc | **4.18e6 J/(m³·K)** (ρ = 1000 kg/m³, c = 4180 J/(kg·K); `RHO_C`, `test_mat21_box_closed_form.py`) |
| Perfusion / flow | none (gel) |
| Heat source | the EM loss density `½ σ |E|²` of the ccw drive above, mapped onto the thermal mesh (one-way; σ does not depend on T) |
| Boundary, all six box faces | convection (Robin): **h = 10 W/(m²·K)** to ambient **T∞ = 20 °C** (`H_GEL`, `test_mat15_coil_driven_heating_phantom.py:83`); no radiation |
| Initial condition | uniform **20 °C** |
| Duration t\* | **2821.5 s** (= 0.2 · C²ρc/k with C = 0.045 m, `T_STAR`, `tests/validation/test_mat21_box_fem.py:72`), source on throughout |
| Time step | ≤ 70.5375 s (40 steps or finer; ours is Crank–Nicolson, N = 40); state the scheme and step used |

If the thermal solver cannot take a convection coefficient on a bare solid
face without a fluid region, state what was used instead; do not substitute
a fixed-temperature wall silently. (`ANS-9`'s replication used AEDT 2026 R1
IcepakFEA through PyAEDT `Mechanical`, Transient Thermal, with one setup per
report time and automatic substepping — state the same details here.)

## Quantities to export (rise θ = T − 20 °C, in K, per watt deposited in the box; both field strengths)

| Row | Quantity |
|---|---|
| 1 | Total power deposited in the box by the mapped source, W (should read 1.000; the mapping's conservation) |
| 2 | **Point rise θ(t\*) at seven probes** (metres, fixture frame: x across the bore, y vertical, z along the bore axis): **c\*** (1.5 T: (−0.105, 0, −0.050); 3 T: (+0.105, 0, +0.050) — each field strength's own SAR10g lattice argmax under ccw), **P0** centre (0, 0, 0), **P1** (+0.21, 0, 0) and **P2** (−0.21, 0, 0) on the long-side walls, **P3** (0, 0, +0.325) and **P4** (0, 0, −0.325) on the end walls, **P5** (0, +0.045, 0) on the top face (`test_geo36_rod_f_clinical_heating.py:155–157`, `_p37_fixture`, P37's points). cw at c\* only. |
| 3 | Initial slope of the volume-mean rise, K/s, from the first time step (calorimetry: equals 1 W / (ρc·V) = 9.737e-06 K/s) |
| 4 | Heat leaving through the six Robin faces **at t\***, W, as a fraction of the deposited watt (ours prints the instantaneous balance "stored + surface = 1 W" at t\*) |
| 5 | Volume-mean rise ⟨θ⟩(t\*), K (a record until our degree-2 column prints it — §10 P87 adds the print) |
| 6 | **P_box / P_acc** under the same ccw drive from the EM design — the factor that turns rows 2–5 into per-accepted-watt figures (it is `ANS-11` / `ANS-12` row 5's gel fraction) |
| 7 | Solve metadata: thermal solver and version, element order of both solves, thermal mesh size, scheme and step, wall time |

Row 2 may be exported as the field value at the point or as the mean over a
ball of radius ≤ 5 mm centred on it (`ANS-9`'s probe-volume clause) — state
which; the wall probes P1–P5 sit on the Robin faces.

## Readout (pre-registered; adjudicated on entries 69 / 70 only)

- Rows 3 and 4 are integral quantities: agreement within **5 %** is AGREE
  (`ANS-9`'s level).
- Row 2 is pointwise on a source hottest near the long-side walls: within
  **10 %** at c\* and P0–P5, with the ordering **{P1, P2} > {P3, P4} > P5 >
  P0** preserved in both codes (ours at both field strengths; P1 / P2 and
  P3 / P4 are near-ties and are not ordered against each other), is AGREE.
  **Tie clause at c\*** (as `ANS-11` v2.1): if AED's SAR10g argmax is the
  other member of c\*'s mirror pair, row 2 is read at AED's own c\* and at
  ours, both exported.
- Row 6 against ours is `ANS-11` / `ANS-12` row 5's question (gel fraction
  within 5 %); it is recorded here, adjudicated there.
- Row 5 is a record until our column prints it; row 7 is a record.
- A disagreement on rows 3 / 4 with row 1 at 1 W on both sides is a finding
  about one code's thermal solve or source mapping: a known-issues entry and
  a diagnosis chunk follow. A row-2 disagreement at one probe with rows 3 / 4
  AGREE is read first against our own degree-1 → degree-2 move at that
  probe (`ANS-9` probe D's lesson, 2026-10-12 epitaph).
- No compliance, tissue, homogeneity or C95.3 claim follows at any outcome;
  the gel is 20 °C with no container wall.

## Our columns (degree 1, records; θ per **accepted** watt as logged, with the factor P_box / P_acc beside it)

Both windows, `-n 8`, every executed assertion green — EM gates ((A0), (A0m),
(U1), (R1), (P) / NEG, (P_RC), coverage (i)) and thermal ((B), (E) ≤ 1e-9,
(E′), the global (S2)). The readout's hotspot-ball slope and "(S0)" were
never asserted under those names (the 11-11 weekly's ruling); the DG0 1 g
tip read is known-issues 2026-11-06's scope.

| Quantity, ccw (1, −1j) | 1.5 T — entry 59, `20261110T080008Z_GEO-38-P33.log` (852 s, 41.70 GiB, node 193.5996 pF) | 3 T — entry 60, `20261111T080008Z_GEO-38-P33-3T.log` (732 s, 32.75 GiB, node 44.0525 pF) |
|---|---|---|
| P_box / P_acc (the factor; log "per-accepted-watt factor") | 5.285093370e-02 (`:12325`) | 1.952048974e-01 (`:12275`) |
| θ(t\*) at c\* per P_acc watt [K/W] | 1.611905088e-03 at (−0.105, 0, −0.050) (`:12326`) | 6.516151132e-03 at (+0.105, 0, +0.050) (`:12276`) |
| θ(t\*) at P0 … P5 per P_acc watt [K/W] | 1.591340145e-04, 4.445002899e-03, 4.110246040e-03, 7.476681123e-04, 8.068176753e-04, 2.497948048e-04 (`:12327–12332`) | 1.317955001e-03, 1.272564088e-02, 1.213484346e-02, 2.939765335e-03, 3.087193253e-03, 2.365388045e-03 (`:12277–12282`) |
| cw (1, +1j): θ(t\*) at c\* per P_acc watt; factor | 1.266521350e-03; 4.165210008e-02 (`:12344`, `:12343`) | 5.749650374e-03; 1.785856659e-01 (`:12294`, `:12293`) |
| 1 g ball mean θ(t\*) at c\* at 1 W deposited [K] | 3.061023711e-02 (`:12317`) | 3.340195106e-02 (`:12267`) |
| Surface heat at t\* per deposited watt (row 4) | 1.686884596e-01 W (`:12324`, (E′)) | 1.678897963e-01 W (`:12274`, (E′)) |
| 360 s running-mean maximum at c\* per P_acc watt | 1.507840749e-03 K (still rising at t\*) | 6.095554712e-03 K (still rising at t\*) |

Per deposited watt, row 2 is the logged value divided by the factor (review
arithmetic: c\* 3.050e-02 K at 1.5 T, 3.338e-02 K at 3 T — the point value
beside the 1 g ball mean above); the `ans:13` page (C53) prints the exact
quotients. **The 3 T / 1.5 T ratio of θ(c\*, t\*) per accepted watt is
4.0425 ccw / 4.5397 cw** (C39, `20261111T103024Z_GEO-38-C39-W1.log:373`,
`:384`), against the SAR10g lever 4.1014 (C33) and the closed-form pictures'
3.6935 (P83 — the power-share part alone).

## Our side

The adjudicated columns are **entries 69 (1.5 T) and 70 (3 T)** — P61's and
C31's modules at degree 2 at the degree-2 roots (`xxl`; PENDING on §10 P87
and on the roots), read through the same rows; the degree-1 table above is a
labelled record, as `ANS-9`'s degree-1 column is. `COMPARISON.md` carries our
rows with the AED columns blank (C54).

**AED / Ansys values are private (CLAUDE.md hard rule): they go in the
gitignored `aed_results/` and `COMPARISON_private.md` only; the qualitative
verdict and a relative agreement level are all that reaches a tracked file.**

**Out of scope:** no implant (`GEO-36` steps 10 / 12, entries 63 / 67, are
the implant-heating windows and have no `ANS` case until W3's reference is
in the repo); no container wall (`MAT-21`'s P48 family is closed form + FEM
only); no duty cycle, no two-way σ(T), no perfusion; no dose row; no
compliance statement; no claim at degree 1.

## Revision history

- **v1 (2026-11-11, weekly planning review; skeleton):** EM half by
  reference to `ANS-11` / `ANS-12`; thermal half, seven probes, rows and
  readout as above; the degree-1 record columns from entries 59 / 60. The
  3 T node is `ANS-12` v2's.
- **v1.0.1 (2026-11-12, 03:00 daily review; errata, for the next weekly to
  ratify or revert):** three transcription corrections against our own
  tracked logs, no AED figure, no band, nothing the operator builds changes.
  The 3 T cw factor read `1.735230713e-01`, a token in no log; the log
  prints `1.785856659e-01` at `20261111T080008Z_GEO-38-P33-3T.log:12293`
  (C46 read it green, `20261111T103024Z_GEO-38-C39-W1.log:300`), and the
  citation `:12304` is corrected to `:12293`. V and Robin area cite entry
  59 `:12315`, not `:12317`. P37's points are `_p37_fixture` at
  `test_geo36_rod_f_clinical_heating.py:155–157`, not `:443–445`. Found by
  the C53 drafter (§9 item 958).
