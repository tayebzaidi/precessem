# ANS-8 step 1 — comparison (our side; AED column blank)

Rows from `fem_em_solver.ports.ring_modes.ring_mode_rows` on our 32×32 (F-human, 64 MHz, degree 1, κ = 4.590991e-03, 507266 cells), C = 1.1311e-10 F on every gap, feeds P17 / P21 (Norton, 50 Ω). AED values are never written here (PROJECT_PLAN §5.4).

| Row | ours | AED |
|---|---|---|
| (1) C_mode from k = 1 [pF] | 113.1081 | |
| (1) C from k = 2 [pF] | 75.4544 | |
| (1) off-block norm | 4.182887e-03 | |
| (1) k = 1 eigenpair ring-antisymmetric fraction (v2) | 1.000000 | |
| (2) rows below at C = 113.11 pF: detuning record (SPEC v2) | | |
| (2) p17_alone feed-ring (bottom) k = ±1 share | 0.789842 | |
| (2) p17_alone feed-ring (bottom) k = 0 leak | 0.173562 | |
| (2) p17_alone other ring (top) k = ±1 share | 0.951117 | |
| (2) ccw feed-ring (bottom) k = ±1 share | 0.789066 | |
| (2) ccw feed-ring (bottom) k = 0 leak | 0.173623 | |
| (2) ccw other ring (top) k = ±1 share | 0.950224 | |
| (2) cw feed-ring (bottom) k = ±1 share | 0.790978 | |
| (2) cw feed-ring (bottom) k = 0 leak | 0.173454 | |
| (2) cw other ring (top) k = ±1 share | 0.952417 | |
| (3) ccw V/I at P17 [Ω] | +1.891406e+01 +2.364137e+01j | |
| (3) ccw V/I at P21 [Ω] | +3.978820e+01 +2.568597e+01j | |
| (3) feed 2×2 Y [S], others in C | +2.140297e-02-1.911774e-02j, +7.453775e-03+1.262029e-03j; +7.453775e-03+1.262029e-03j, +2.144217e-02-1.911256e-02j | |
| (2b) own C_mode [pF] (from row 1, every gap incl. feeds) | 113.108077 | |
| (2b) p17_alone feed-ring (bottom) k = ±1 share | 0.789857 | |
| (2b) p17_alone feed-ring (bottom) k = 0 leak | 0.173544 | |
| (2b) p17_alone other ring (top) k = ±1 share | 0.951115 | |
| (2b) ccw feed-ring (bottom) k = ±1 share | 0.789081 | |
| (2b) ccw feed-ring (bottom) k = 0 leak | 0.173606 | |
| (2b) ccw other ring (top) k = ±1 share | 0.950222 | |
| (2b) cw feed-ring (bottom) k = ±1 share | 0.790993 | |
| (2b) cw feed-ring (bottom) k = 0 leak | 0.173436 | |
| (2b) cw other ring (top) k = ±1 share | 0.952414 | |
| (3b) ccw V/I at P17, own C [Ω] | +1.891461e+01 +2.363992e+01j | |
| (3b) ccw V/I at P21, own C [Ω] | +3.978927e+01 +2.567924e+01j | |
| (3b) feed 2×2 Y [S], others in own C | +2.140529e-02-1.911715e-02j, +7.454496e-03+1.261502e-03j; +7.454496e-03+1.261502e-03j, +2.144449e-02-1.911197e-02j | |
| (2c) opposite, C = 113.11 pF: feed ring share / leak; other ring share / leak | 0.988773 / 2.177e-08; 0.987501 / 9.456e-07 | |
| (2c) four_ccw, C = 113.11 pF: feed ring share / leak; other ring share / leak | 0.987648 / 1.300e-08; 0.986318 / 9.309e-07 | |
| (2c) opposite, own C_mode: feed ring share / leak; other ring share / leak | 0.988772 / 2.177e-08; 0.987500 / 9.455e-07 | |
| (2c) four_ccw, own C_mode: feed ring share / leak; other ring share / leak | 0.987648 / 1.300e-08; 0.986317 / 9.307e-07 | |
| (4a) f₀ = f·√(C_mode / C) [Hz] (one-point, λ₁ ∝ ω) | 63999455.96 | |

Loss partition: not computed on our side in step 1 (left blank).

### Per-ring k-spectra (ours)

| Drive | ring | k = ±1 share | k = 0 leak | \|k\| folds 0…8 |
|---|---|---|---|---|
| ccw | bottom | 0.789066 | 0.173623 | 0.1736 0.7891 0.0227 0.0068 0.0029 0.0020 0.0013 0.0011 0.0005 |
| ccw | top | 0.950224 | 0.000396 | 0.0004 0.9502 0.0298 0.0090 0.0039 0.0027 0.0017 0.0015 0.0007 |
| cw | bottom | 0.790978 | 0.173454 | 0.1735 0.7910 0.0228 0.0054 0.0029 0.0016 0.0013 0.0011 0.0006 |
| cw | top | 0.952417 | 0.000418 | 0.0004 0.9524 0.0299 0.0073 0.0039 0.0021 0.0017 0.0015 0.0008 |
| p17_alone | bottom | 0.789842 | 0.173562 | 0.1736 0.7898 0.0227 0.0062 0.0029 0.0018 0.0013 0.0011 0.0005 |
| p17_alone | top | 0.951117 | 0.000414 | 0.0004 0.9511 0.0298 0.0083 0.0039 0.0024 0.0017 0.0015 0.0007 |

### Per-ring k-spectra at own C_mode (2b, ours, v2)

| Drive | ring | k = ±1 share | k = 0 leak | \|k\| folds 0…8 |
|---|---|---|---|---|
| ccw | bottom | 0.789081 | 0.173606 | 0.1736 0.7891 0.0227 0.0068 0.0029 0.0020 0.0013 0.0011 0.0005 |
| ccw | top | 0.950222 | 0.000396 | 0.0004 0.9502 0.0298 0.0090 0.0039 0.0027 0.0017 0.0015 0.0007 |
| cw | bottom | 0.790993 | 0.173436 | 0.1734 0.7910 0.0228 0.0054 0.0029 0.0016 0.0013 0.0011 0.0006 |
| cw | top | 0.952414 | 0.000418 | 0.0004 0.9524 0.0299 0.0073 0.0039 0.0021 0.0017 0.0015 0.0008 |
| p17_alone | bottom | 0.789857 | 0.173544 | 0.1735 0.7899 0.0227 0.0062 0.0029 0.0018 0.0013 0.0011 0.0005 |
| p17_alone | top | 0.951115 | 0.000414 | 0.0004 0.9511 0.0298 0.0083 0.0039 0.0024 0.0017 0.0015 0.0007 |

Not claimed: any AED value, verdict or agreement level; sweep row 4; Route A; copper. The weekly adjudicates share and leak at a stated relative level.

Not claimed (v2): the row-4 sweep interpolation (4a is the one-point form, λ₁ taken ∝ ω, its first-order error (1 − p)·δ/2 unmeasured); the degree-2 32×32; copper; P_coil; Route A; adjudication.

## Step 2 (SPEC v3) — tuned two-port drive with ESR, degree 1 (ours; AED column blank)

The tables above are **lossless (step 1)** (σ = 800 S/m coil, C = 113.11 pF, no ESR). This section: H9′'s copper F-human (hole mesh, Leontovich, κ off), degree 1, the 0.015 m phantom, C_in = 118.2437 pF and r = 0.2 Ω on all 32 gaps, ccw (1, −j) at P17, P21, 64 MHz, 369387 cells. AED values are never written here (PROJECT_PLAN §5.4).

**Accepted power (SPEC v3, weekly ruling 2026-10-12).** Per-W / per-√W figures are per W of P_acc, the sum over the two feed ports of the accepted power at the 50 Ω lumped port: incident − reflected at the port plane, HFSS's default per-port accepted power. Ours is the line cut ½ Re Σ_feeds V_port·conj(I_line) = 3.057630e-02 W at unit Thevenin drive (relative difference from the wave form (‖v_s‖² − ‖S v_s‖²) / (8 z0): 3.4e-16). The feed gaps' capacitor and 0.2 Ω ESR are load inside P_acc: their dissipation is part of row 5's ESR fraction. The field's feed-sheet partition exceeds the line cut by 0.83 %: the discrete feed sheet's Cauchy–Schwarz excess (CS / X = 1.000000000 at P17, 1.000000000 at P21), a documented sheet property (weekly ruling 2026-10-12; known-issues 2026-10-11), not a defect and not a quantity the AED column is asked for. The per-partition values are printed beside each row.

| Row | ours | AED |
|---|---|---|
| (1) Y_in ccw at P17 [S] | +1.010006e-02 +1.017991e-05j | |
| (1) Y_in ccw at P21 [S] | +1.184993e-02 -2.164135e-04j | |
| (2) feed ring (bottom) k = ±1 share (in-model sheet currents) | 0.993845 | |
| (2) feed ring (bottom) k = 0 leak | 4.739399e-03 | |
| (2) other ring (top) k = ±1 share (in-model sheet currents) | 0.998448 | |
| (2) other ring (top) k = 0 leak | 1.022353e-05 | |
| (3) \|B₁⁺\| at O per √W [T/√W] | 1.668096e-06 (per W of the field partition: 1.675076e-06) | |
| (3) \|B₁⁺\| at P1 per √W [T/√W] | 1.652840e-06 (per W of the field partition: 1.659756e-06) | |
| (3) \|B₁⁺\| at P2 per √W [T/√W] | 1.623936e-06 (per W of the field partition: 1.630731e-06) | |
| (3) \|B₁⁺\| at P3 per √W [T/√W] | 1.610663e-06 (per W of the field partition: 1.617403e-06) | |
| (3) \|B₁⁺\| at P4 per √W [T/√W] | 1.648538e-06 (per W of the field partition: 1.655437e-06) | |
| (4) peak 10 g SAR per W (sphere, lattice argmax) [W/kg/W] | 1.587143e-01 (per W of the field partition: 1.600454e-01) | |
| (4) its location (x, y, z) [m] | (0.044135, 0.008779, 0.000000) | |
| (5) phantom fraction of P_acc | 0.268112 (of the field partition: 0.270361) | |
| (5) copper (Leontovich) fraction of P_acc | 0.019926 (of the field partition: 0.020094) | |
| (5) ESR fraction of P_acc (30 gaps + 2 feed gaps) | 0.703644 (of the field partition: 0.709546) | |
| (5) sum of the fractions of P_acc; residual = feed-sheet Cauchy–Schwarz excess (documented sheet property, weekly 2026-10-12; known-issues 2026-10-11) | 0.991683; residual 0.008317 | |
| (6) \|B₁⁺\| / \|B₁⁻\| at O | 20.7092 | |

Row 2 is `ring_modes.k_spectrum` on the field's sheet currents; the circuit reduction of the stored 32×32 is not our row (H12c identity red at 1.58). c* is a lattice argmax (18 points), not a global maximum.

Not claimed (v3): any AED value, verdict or agreement level; the degree-2 column; f₀, Q, match or tuning accuracy; c* as a global maximum; homogeneity; compliance.

## Step 2 at the degree-2 root (SPEC v4: rows 4–5 from the node record, degree 2; ours; AED column blank)

This block is the degree-2 column at the root, labelled per SPEC v4: degree 2, C_in = 115.5209 pF on all 32 gaps, ccw (1, −j) at P17, P21, 64 MHz, rows from the solved node record (knob d2root). P_acc,SPEC is the v3 section's line cut. AED values are never written here (PROJECT_PLAN §5.4).

| Row | ours | AED |
|---|---|---|
| P_acc,SPEC [W] | 2.910591023475e-02 (degree 1 record: 3.057629563251e-02) | |
| (4) peak 10 g SAR per W of P_acc (sphere, lattice argmax) [W/kg/W] | 0.1562286 (degree 1 record: 0.1587143) | |
| (5) phantom fraction of P_acc | 0.260954 (degree 1 record: 0.268112) | |
| (5) copper (Leontovich) fraction of P_acc | 0.020083 (degree 1 record: 0.019926) | |
| (5) ESR (30 gaps + 2 feed gaps) fraction of P_acc | 0.697465 (degree 1 record: 0.703644) | |
| (5) sum of the fractions of P_acc | 0.978503 (degree 1 record: 0.991683) | |
| X_s at feed sheet P17 [W] | 3.262132656e-04 (degree 1 record: 1.329324349e-04) | |
| X_s at feed sheet P21 [W] | 2.994809404e-04 (degree 1 record: 1.213651853e-04) | |
| ΣX_s / P_acc | 2.149715e-02 (degree 1 record: 8.316822e-03) | |
| row 2 (degree 2, **circuit route — Norton feeds on entry 36's saved 32×32 at C_mode(d2, Cu) = 114.5221 pF; in-model row not in record**) | ccw feed-ring share 0.968088 / k = 0 leak 0.025825 (lossless); 0.948363 / 0.041953 (r = 0.2 Ω, the production convention); source `docs/testing/logs/20261018T095559Z_ANS-8-H25-W1.log:142–143,165–166` | |

Rows 1–3, row 6 and CS_s are not in the node record at degree 2 (§10 H30 as narrowed by the 2026-10-21 weekly).

Not claimed (degree 2): any AED value, verdict or agreement level; rows 1–3 / 6 / CS_s; an order-convergence statement; f₀, Q, match or tuning accuracy; compliance.

Row 2 is now carried at degree 2 from the circuit route (§10 H35), labelled as that route at H25's own C_mode(d2, Cu) = 114.5221 pF, not at the 115.5209 pF root; rows 1, 3, 6 and CS_s stay not in record.
