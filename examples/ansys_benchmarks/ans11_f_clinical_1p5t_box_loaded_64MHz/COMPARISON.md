# ANS-11 — comparison (our degree-1 and degree-2 records; AED column blank)

F-clinical-1.5T loaded with the ASTM F2182 gel box, two-port ccw drive (V_33, V_41) = (1, −1j) at 64 MHz, C = 193.5996 pF + 0.2 Ω on all 64 gaps (entry 52's node), raw values at V_th = 1 V per feed. **Degree 1 is a record, never adjudicated. The degree-2 column is C19's record** (`xxl` entry 62, `-n 16`, degree 2, at SPEC v2's node 185.7629 pF + 0.2 Ω, parsed from `20261114T080009Z_GEO-38-C19.log` at the cited lines), **still un-adjudicated** (Done-when (2) is the weekly's). Every value is the printed token of `20261103T080009Z_GEO-38-C6.log` at the cited line (parsed by `11_f_clinical_1p5t_box_loaded_64MHz.py`; no solve). AED values are never written here (PROJECT_PLAN §5.4). Rows follow [SPEC.md](SPEC.md).

| Row | Quantity | ours, degree 1 (entry 52, record) | ours, degree 2 (entry 62, record) | AED First Order |
|---|---|---|---|---|
| 1 | S₁₁ | -4.8167140712e-01+7.9527084190e-04j (`:12281`) | -4.2727123652e-01-1.2826446214e-03j (`:12379`) | |
| 1 | S₁₂ = S₂₁ | 4.6496927594e-03+1.0179009021e-01j | 4.2062703482e-03+9.9390915599e-02j | |
| 1 | S₂₂ | -5.0211295667e-01-2.1994222327e-02j | -4.5534678749e-01-2.0906895611e-02j | |
| 1 | Y_in ccw P33 [S] | +4.450128950e-02 +4.009154872e-04j (Im +; `:12285`) | +3.950924174e-02 +4.859866152e-04j (Im +; `:12383`) | |
| 1 | Y_in ccw P41 [S] | +8.079211370e-02 +4.413545175e-03j (Im +) | +6.970848681e-02 +3.364730927e-03j (Im +) | |
| 2 | 18-point lattice mean `\|B₁⁺\|`_ccw [T] (the comparison row) | 1.49792e-08 (arithmetic on `:12289–12306`) | 1.47886e-08 (arithmetic on `:12388–12405`) | |
| 2 | the same per √W of P_acc [T/√W] | 6.0833e-08 | 6.0552e-08 | |
| 3 | c\* (largest of the 18 SAR10g) | x-105mm dz-50mm (r 0.105000 m, dz -0.050000 m, az 180.000°; `:12307`) | x-105mm dz-50mm (r 0.105000 m, dz -0.050000 m, az 180.000°; `:12406`) | |
| 3 | SAR10g(c\*) [W/kg]; per W of P_acc | 1.378896513e-04; 2.274234e-03 (`:12310`) | 1.346502340e-04; 2.257378e-03 (`:12410`) | |
| 3 | cw (1, +1j) ball at c\* [W/kg] (record) | 1.059925022e-04 | 1.005190493e-04 | |
| 3 | box-wide 10 g maximum and location (record) | not read in entry 52 (the C6b reader was off) | x+157.5mm y+0.0mm dz+0.0mm: 3.151138314e-04 (`:12407`) | |
| 4 | 400-centroid gel mean `\|B₁⁺\|`_ccw [T] (second record) | 2.090267723e-08 over 400 / 400 tag-3 centroids (`:12288`) | 2.063058069e-08 over 400 / 400 tag-3 centroids (`:12387`) | |
| 4 | the same per √W of P_acc [T/√W] | 8.488944e-08 | 8.447147e-08 | |
| 5 | P_box / P_src / P_acc [W] | 3.204417752e-03 / 8.266643116e-02 / 6.063124201e-02 | 3.112805571e-03 / 8.413089485e-02 / 5.964895108e-02 | |
| 5 | gel fraction P_box / P_acc | 5.285093e-02 | 5.218542e-02 | |
| 5 | copper / capacitor-ESR fractions under ccw | not printed by entry 52 (per-column partitions only, `:12282–12284`) | | |
| 6 | mean `\|B₁⁺\|` / mean `\|B₁⁻\|` over the gel (sense) | 2.963374 (`:12288`) | 2.963373 (`:12387`) | |

## Row 2 / 3 lattice (ours, degree 1; raw at V_th = 1 V, and per P_acc)

| Point | log line | `\|B₁⁺\|`_ccw [T] | per √W P_acc [T/√W] | SAR10g [W/kg] | per W P_acc [W/kg/W] |
|---|---|---|---|---|---|
| x-105mm dz-250mm | `:12289` | 5.319736310e-09 | 2.16044e-08 | 7.896193072e-05 | 1.30233e-03 |
| x-105mm dz-150mm | `:12290` | 1.664794126e-08 | 6.76102e-08 | 1.153727714e-04 | 1.90286e-03 |
| x-105mm dz-50mm | `:12291` | 2.230736624e-08 | 9.05941e-08 | 1.378896513e-04 | 2.27423e-03 |
| x-105mm dz+50mm | `:12292` | 2.255456765e-08 | 9.15981e-08 | 1.368277819e-04 | 2.25672e-03 |
| x-105mm dz+150mm | `:12293` | 1.691890098e-08 | 6.87106e-08 | 1.118693083e-04 | 1.84508e-03 |
| x-105mm dz+250mm | `:12294` | 4.981574558e-09 | 2.02310e-08 | 7.397515403e-05 | 1.22008e-03 |
| x+0mm dz-250mm | `:12295` | 5.345956882e-09 | 2.17109e-08 | 1.106314017e-04 | 1.82466e-03 |
| x+0mm dz-150mm | `:12296` | 1.526969839e-08 | 6.20129e-08 | 7.066740319e-05 | 1.16553e-03 |
| x+0mm dz-50mm | `:12297` | 2.084487448e-08 | 8.46547e-08 | 1.150657127e-05 | 1.89780e-04 |
| x+0mm dz+50mm | `:12298` | 2.082999222e-08 | 8.45943e-08 | 9.807225805e-06 | 1.61752e-04 |
| x+0mm dz+150mm | `:12299` | 1.528202406e-08 | 6.20630e-08 | 6.505140565e-05 | 1.07290e-03 |
| x+0mm dz+250mm | `:12300` | 5.246062790e-09 | 2.13052e-08 | 1.038150798e-04 | 1.71224e-03 |
| x+105mm dz-250mm | `:12301` | 5.140520390e-09 | 2.08766e-08 | 8.277082975e-05 | 1.36515e-03 |
| x+105mm dz-150mm | `:12302` | 1.895998609e-08 | 7.69998e-08 | 1.160653596e-04 | 1.91428e-03 |
| x+105mm dz-50mm | `:12303` | 2.519918694e-08 | 1.02338e-07 | 1.378032615e-04 | 2.27281e-03 |
| x+105mm dz+50mm | `:12304` | 2.486936097e-08 | 1.00999e-07 | 1.349823646e-04 | 2.22628e-03 |
| x+105mm dz+150mm | `:12305` | 1.842139994e-08 | 7.48125e-08 | 1.106521593e-04 | 1.82500e-03 |
| x+105mm dz+250mm | `:12306` | 5.485954992e-09 | 2.22794e-08 | 7.651316756e-05 | 1.26194e-03 |

## Row 2 / 3 lattice (ours, degree 2; raw at V_th = 1 V, and per P_acc)

| Point | degree | log line | `\|B₁⁺\|`_ccw [T] | per √W P_acc [T/√W] | SAR10g [W/kg] | per W P_acc [W/kg/W] |
|---|---|---|---|---|---|---|
| x-105mm dz-250mm | 2 | `:12388` | 5.268208021e-09 | 2.15706e-08 | 7.679531417e-05 | 1.28745e-03 |
| x-105mm dz-150mm | 2 | `:12389` | 1.642481869e-08 | 6.72511e-08 | 1.121139963e-04 | 1.87956e-03 |
| x-105mm dz-50mm | 2 | `:12390` | 2.206869697e-08 | 9.03598e-08 | 1.346502340e-04 | 2.25738e-03 |
| x-105mm dz+50mm | 2 | `:12391` | 2.227693412e-08 | 9.12124e-08 | 1.332203910e-04 | 2.23341e-03 |
| x-105mm dz+150mm | 2 | `:12392` | 1.672608459e-08 | 6.84846e-08 | 1.082819001e-04 | 1.81532e-03 |
| x-105mm dz+250mm | 2 | `:12393` | 4.926130961e-09 | 2.01699e-08 | 7.265538569e-05 | 1.21805e-03 |
| x+0mm dz-250mm | 2 | `:12394` | 5.265810204e-09 | 2.15607e-08 | 1.076638402e-04 | 1.80496e-03 |
| x+0mm dz-150mm | 2 | `:12395` | 1.508172096e-08 | 6.17518e-08 | 6.831118241e-05 | 1.14522e-03 |
| x+0mm dz-50mm | 2 | `:12396` | 2.058504417e-08 | 8.42850e-08 | 1.104530526e-05 | 1.85172e-04 |
| x+0mm dz+50mm | 2 | `:12397` | 2.060088920e-08 | 8.43499e-08 | 9.301285785e-06 | 1.55934e-04 |
| x+0mm dz+150mm | 2 | `:12398` | 1.507667565e-08 | 6.17311e-08 | 6.346880841e-05 | 1.06404e-03 |
| x+0mm dz+250mm | 2 | `:12399` | 5.216555944e-09 | 2.13591e-08 | 1.013868587e-04 | 1.69973e-03 |
| x+105mm dz-250mm | 2 | `:12400` | 5.061057760e-09 | 2.07224e-08 | 7.874806630e-05 | 1.32019e-03 |
| x+105mm dz-150mm | 2 | `:12401` | 1.869455812e-08 | 7.65445e-08 | 1.127575632e-04 | 1.89035e-03 |
| x+105mm dz-50mm | 2 | `:12402` | 2.484351396e-08 | 1.01721e-07 | 1.340265160e-04 | 2.24692e-03 |
| x+105mm dz+50mm | 2 | `:12403` | 2.452617941e-08 | 1.00422e-07 | 1.323408381e-04 | 2.21866e-03 |
| x+105mm dz+150mm | 2 | `:12404` | 1.814979097e-08 | 7.43139e-08 | 1.084982031e-04 | 1.81895e-03 |
| x+105mm dz+250mm | 2 | `:12405` | 5.402883121e-09 | 2.21220e-08 | 7.429744646e-05 | 1.24558e-03 |

## Degree 1 → degree 2 (records)

Move = degree-2 value / degree-1 value − 1, from the two tracked logs' tokens (a ratio of our own numbers; no AED figure). Per-watt rows are the order rule's.

| Quantity | per-watt row | move |
|---|---|---|
| |S11| | no | -1.129e-01 |
| |S22| | no | -9.305e-02 |
| |Y_in ccw P33| | no | -1.121e-01 |
| |Y_in ccw P41| | no | -1.375e-01 |
| 18-point lattice mean |B1+|_ccw | no | -1.272e-02 |
| 18-point lattice mean |B1+|_ccw per sqrt(P_acc) | yes | -4.624e-03 |
| SAR10g(c*) | no | -2.349e-02 |
| SAR10g(c*) per P_acc | yes | -7.412e-03 |
| cw (1,+1j) ball at c* | no | -5.164e-02 |
| 400-centroid mean |B1+|_ccw | no | -1.302e-02 |
| 400-centroid mean |B1+|_ccw per sqrt(P_acc) | yes | -4.924e-03 |
| P_box | no | -2.859e-02 |
| P_src | no | +1.772e-02 |
| P_acc | no | -1.620e-02 |
| gel fraction P_box / P_acc | yes | -1.259e-02 |
| sense mean |B1+| / mean |B1-| | no | -3.375e-07 |

- Order rule (the 10-24 exit's, |move| > 2 % on a per-watt row; a printed verdict, not an assertion): rows over: none. **The 2 % order rule does not fire.**
- c\*: x-105mm dz-50mm → x-105mm dz-50mm (degree 1 → degree 2). The box-wide g\* is degree-2 only (entry 52 ran no grid): x+157.5mm y+0.0mm dz+0.0mm (0.157500, 0.000000, 0.000000) m, `:12407`; no claim about a maximum outside the grid.
- Degree-2 S₁₂ and S₂₁ strings (`:12379`): 4.2062703482e-03+9.9390915599e-02j / 4.2062703483e-03+9.9390915598e-02j — they differ in the last printed digit; reciprocity is the window's asserted (R1) on that line, not a string identity.

## Provenance

- Source: `docs/testing/logs/20261103T080009Z_GEO-38-C6.log` (`xl` entry 52, `-n 8`; price line `:12312`: build 460.8 s, window 862.8 s).
- Degree-2 source: `docs/testing/logs/20261114T080009Z_GEO-38-C19.log` (`xxl` entry 62, `-n 16`; price line `:12412`: build 480.2 s, window 4301.9 s). Entry 62's host JSON is a custody check only (the page's (J62)); this file never reads it.
- P_acc re-derived from the printed S and V_s through C6's `_p_acc_spec`: 6.063124e-02 W (ccw); the cw weights give 6.0953e-02 W (N-sense, printed only).
- Not claimed: any adjudication or AED number (Done-when (2)), the rod rows (entry 62's call 2), a setup figure, and any reading of the order rule's verdict for other clinical figures (the weekly's). The moves compare entry 52 (degree 1, C = 193.5996 pF) with entry 62 (degree 2, C = 185.7629 pF): each move carries the node change as well as the order change.
