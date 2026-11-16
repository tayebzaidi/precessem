# ANS-12 — comparison (our degree-1 record; degree-2 and AED columns blank)

F-clinical-3T loaded with the ASTM F2182 gel box, two-port ccw drive (V_33, V_41) = (1, −1j) at 128 MHz, C = 44.0525 pF + 0.2 Ω on all 64 gaps (entry 57's degree-1 root; not the node to enter, SPEC v2), raw values at V_th = 1 V per feed. **Degree 1 is a record, never adjudicated; the adjudication column is degree 2** (entry 66, pending; its node is SPEC v2's). Every value is the printed token of `20261109T080008Z_GEO-38-C6-3T.log` at the cited line (parsed by `11_f_clinical_1p5t_box_loaded_64MHz.py --case ans12`; no solve). AED values are never written here (PROJECT_PLAN §5.4). Rows follow [SPEC.md](SPEC.md).

| Row | Quantity | ours, degree 1 (entry 58, record) | ours, degree 2 (entry 66, pending) | AED First Order |
|---|---|---|---|---|
| 1 | S₁₁ | 3.6878613057e-01-6.7206972791e-02j (`:12228`) | | |
| 1 | S₁₂ = S₂₁ | -1.5727100565e-02+4.2880439168e-02j | | |
| 1 | S₂₂ | 2.0011568377e-01+6.0905743336e-02j | | |
| 1 | Y_in ccw P33 [S] | +8.297670657e-03 +1.031943733e-03j (Im +; `:12232`) | | |
| 1 | Y_in ccw P41 [S] | +1.451253979e-02 -1.347374887e-03j (Im −) | | |
| 2 | 18-point lattice mean `\|B₁⁺\|`_ccw [T] (the comparison row) | 1.17493e-08 (arithmetic on `:12238–12255`) | | |
| 2 | the same per √W of P_acc [T/√W] | 8.5816e-08 | | |
| 3 | c\* (largest of the 18 SAR10g) | x+105mm dz+50mm (r 0.105000 m, dz +0.050000 m, az 0.000°; `:12256`) | | |
| 3 | SAR10g(c\*) [W/kg]; per W of P_acc | 1.748477085e-04; 9.327634e-03 (`:12259`) | | |
| 3 | cw (1, +1j) ball at c\* [W/kg] (record) | 1.427632775e-04 | | |
| 3 | box-wide 10 g maximum and location (record) | not read in entry 58 (the C6b reader was off) | | |
| 4 | 400-centroid gel mean `\|B₁⁺\|`_ccw [T] (second record) | 1.562098621e-08 over 400 / 400 tag-3 centroids (`:12237`) | | |
| 4 | the same per √W of P_acc [T/√W] | 1.140944e-07 | | |
| 5 | P_box / P_src / P_acc [W] | 3.659141247e-03 / 5.347172304e-02 / 1.874513035e-02 | | |
| 5 | gel fraction P_box / P_acc | 1.952049e-01 | | |
| 5 | copper / capacitor-ESR fractions under ccw | not printed by entry 58 (per-column partitions only, `:12229–12231`) | | |
| 6 | mean `\|B₁⁺\|` / mean `\|B₁⁻\|` over the gel (sense) | 1.674561 (`:12237`) | | |

## Row 2 / 3 lattice (ours, degree 1; raw at V_th = 1 V, and per P_acc)

| Point | log line | `\|B₁⁺\|`_ccw [T] | per √W P_acc [T/√W] | SAR10g [W/kg] | per W P_acc [W/kg/W] |
|---|---|---|---|---|---|
| x-105mm dz-250mm | `:12238` | 9.020981381e-09 | 6.58885e-08 | 1.014256845e-04 | 5.41078e-03 |
| x-105mm dz-150mm | `:12239` | 1.474242393e-08 | 1.07677e-07 | 1.247755343e-04 | 6.65642e-03 |
| x-105mm dz-50mm | `:12240` | 1.712663981e-08 | 1.25092e-07 | 1.714031706e-04 | 9.14388e-03 |
| x-105mm dz+50mm | `:12241` | 1.747062147e-08 | 1.27604e-07 | 1.674857587e-04 | 8.93489e-03 |
| x-105mm dz+150mm | `:12242` | 1.460517029e-08 | 1.06675e-07 | 1.249718628e-04 | 6.66690e-03 |
| x-105mm dz+250mm | `:12243` | 8.644831167e-09 | 6.31411e-08 | 9.714595086e-05 | 5.18246e-03 |
| x+0mm dz-250mm | `:12244` | 7.447523171e-09 | 5.43961e-08 | 1.313788917e-04 | 7.00869e-03 |
| x+0mm dz-150mm | `:12245` | 7.651445279e-09 | 5.58855e-08 | 3.791187182e-05 | 2.02249e-03 |
| x+0mm dz-50mm | `:12246` | 7.420152519e-09 | 5.41962e-08 | 6.935343981e-06 | 3.69981e-04 |
| x+0mm dz+50mm | `:12247` | 7.351443042e-09 | 5.36943e-08 | 7.045243663e-06 | 3.75844e-04 |
| x+0mm dz+150mm | `:12248` | 7.356691699e-09 | 5.37327e-08 | 3.746219347e-05 | 1.99850e-03 |
| x+0mm dz+250mm | `:12249` | 7.355732300e-09 | 5.37257e-08 | 1.257707735e-04 | 6.70952e-03 |
| x+105mm dz-250mm | `:12250` | 8.849241010e-09 | 6.46341e-08 | 1.009852525e-04 | 5.38728e-03 |
| x+105mm dz-150mm | `:12251` | 1.534803656e-08 | 1.12101e-07 | 1.270816313e-04 | 6.77945e-03 |
| x+105mm dz-50mm | `:12252` | 1.819489693e-08 | 1.32894e-07 | 1.734192670e-04 | 9.25143e-03 |
| x+105mm dz+50mm | `:12253` | 1.821288890e-08 | 1.33025e-07 | 1.748477085e-04 | 9.32763e-03 |
| x+105mm dz+150mm | `:12254` | 1.552658197e-08 | 1.13405e-07 | 1.263617979e-04 | 6.74105e-03 |
| x+105mm dz+250mm | `:12255` | 9.162569629e-09 | 6.69226e-08 | 9.877796131e-05 | 5.26953e-03 |

## Row 3 tie set

Lattice balls with SAR10g within 1 % of the lattice maximum (the 11-09 weekly's tie clause, `ANS-11` v2.1; SPEC `:150`), and the nearest outsider:

| Point | SAR10g [W/kg] | vs the maximum [%] | role |
|---|---|---|---|
| x+105mm dz+50mm | 1.748477085e-04 | 0.00 | maximum |
| x+105mm dz-50mm | 1.734192670e-04 | -0.82 | member |
| x-105mm dz-50mm | 1.714031706e-04 | -1.97 | nearest outsider |

## Provenance

- Source: `docs/testing/logs/20261109T080008Z_GEO-38-C6-3T.log` (`xl` entry 58, `-n 8`; price line `:12261`: build 448.3 s, window 774.2 s).
- P_acc re-derived from the printed S and V_s through C6's `_p_acc_spec`: 1.874513e-02 W (ccw); the cw weights give 1.8962e-02 W (N-sense, printed only).
- Not claimed: the degree-2 rows, any adjudication or AED number, the node to enter (SPEC v2).
