# ANS-13 — comparison (our degree-1 records; degree-2 and AED columns blank)

F-clinical box heating at 1.5 T (64 MHz) and 3 T (128 MHz): the gel box driven by the coil's ccw (1, −1j) SAR, every value **per deposited watt** (SPEC item 4) unless the row says otherwise. **Degree 1 is a record, never adjudicated; the adjudication columns are degree 2** (entries 69 / 70, PENDING on P87). Every value is the printed token of entry 59's or entry 60's log at the cited line (parsed by `11_f_clinical_1p5t_box_loaded_64MHz.py --case ans13` through C46's reader; no solve). AED values are never written here (PROJECT_PLAN §5.4). Rows follow [SPEC.md](SPEC.md).

| Row | Quantity | ours, degree 1, 1.5 T (entry 59, record) | ours, degree 1, 3 T (entry 60, record) | ours, degree 2 (entries 69 / 70, pending) | AED First Order + IcepakFEA |
|---|---|---|---|---|---|
| 1 | Total power deposited in the box by the mapped source [W] (ccw (B) int_sub q, per deposited watt) | 9.999999999999984e-01 (`:12316`) | 9.999999999999993e-01 (`:12266`) | | |
| 2 | θ(t\*) at P0 centre per deposited watt [K] (ccw) | +3.010997221e-03 (`:12318`) | +6.751649256e-03 (`:12268`) | | |
| 2 | θ(t\*) at P1 +y wall per deposited watt [K] (ccw) | +8.410452925e-02 (`:12319`) | +6.519119676e-02 (`:12269`) | | |
| 2 | θ(t\*) at P2 -y wall per deposited watt [K] (ccw) | +7.777054731e-02 (`:12320`) | +6.216464656e-02 (`:12270`) | | |
| 2 | θ(t\*) at P3 +x wall per deposited watt [K] (ccw) | +1.414673422e-02 (`:12321`) | +1.505989539e-02 (`:12271`) | | |
| 2 | θ(t\*) at P4 -x wall per deposited watt [K] (ccw) | +1.526591148e-02 (`:12322`) | +1.581514241e-02 (`:12272`) | | |
| 2 | θ(t\*) at P5 top face per deposited watt [K] (ccw) | +4.726402872e-03 (`:12323`) | +1.211746260e-02 (`:12273`) | | |
| 2 | θ(t\*) at c\* per deposited watt [K] (ccw; per P_acc / factor; 1 g ball beside) | 3.049908e-02 (`:12326` / `:12325`); 1 g ball 3.061023711e-02 (`:12317`) | 3.338108e-02 (`:12276` / `:12275`); 1 g ball 3.340195106e-02 (`:12267`) | | |
| 2 | cw (1, +1j): θ(t\*) at c\* per deposited watt [K] (per P_acc / factor); cw factor | 3.040714e-02 (`:12344` / `:12343`); factor 4.165210008e-02 | 3.219548e-02 (`:12294` / `:12293`); factor 1.785856659e-01 | | |
| 3 | Initial volume-mean slope [K/s]: closed form 1 W / (ρc·V_sub) (the logs print (S2) = abs(slope / closed form − 1), not the slope) | 9.7369e-06 (V_sub 2.457000000000e-02 m³, `:12315`); (S2) 1.1895e-05 (`:12324`) | 9.7369e-06 (V_sub 2.457000000000e-02 m³, `:12265`); (S2) 1.1729e-05 (`:12274`) | | |
| 4 | Surface heat at t\* per deposited watt [W] (ccw (E′)) | 1.686884595994e-01 (`:12324`) | 1.678897962921e-01 (`:12274`) | | |
| 5 | Volume-mean rise ⟨θ⟩(t\*) [K] | P87 prints it | P87 prints it | | |
| 6 | P_box / P_acc (the factor), ccw | 5.285093370e-02 (`:12325`) | 1.952048974e-01 (`:12275`) | | |

## Provenance

- Entry 59 (1.5 T): `docs/testing/logs/20261110T080008Z_GEO-38-P33.log` (node 193.5996 pF, `-n 8`, 852 s; CN N = 40; P61 block `:139–12356`).
- Entry 60 (3 T): `docs/testing/logs/20261111T080008Z_GEO-38-P33-3T.log` (node 44.0525 pF, `-n 8`, 732 s; CN N = 40; P61 block `:141–12307`).
- Not claimed: any degree-2 value or reader (entries 69 / 70, P87); row 3's slope value (only (S2) is printed); rows 5 and 7 as rows; an adjudication or any AED number; a SPEC edit.
