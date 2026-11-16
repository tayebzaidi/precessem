# Clinical-scale body-coil fixtures (F-clinical 1.5 T / 3 T)

This page collects the geometry of the two clinical-scale birdcage fixtures
of `GEO-38` (PROJECT_PLAN §10 chain C), the sources each number comes from,
and the rulings that fixed the modelling choices. Every number in the
geometry table below is checked against the code by
`tests/package/test_geo38_clinical_fixture_page.py`, and every number in the
heating section by `tests/package/test_wf17_clinical_heating_page.py`, and
every number in the rod-in-box section by
`tests/package/test_wf17_clinical_rod_records_page.py`, and every number in
the closed-form gel-box section by
`tests/package/test_wf17_clinical_closed_form_records_page.py`.

## What the two fixtures are

F-clinical-1.5T (64 MHz) and F-clinical-3T (128 MHz) are **nominal clinical
body-coil geometries, not scanner models**: 32-leg high-pass birdcages inside
a closed PEC can, with the ASTM F2182 gel box as the load. F-clinical is the
production fixture for new human-scale work; the Phase 6 exit stays on
F-human (PROJECT_PLAN §10, the clinical-scale coil fixtures paragraph).

Both are built by the same generator. The 3 T table is derived from the
1.5 T table by overriding five values (rung diameter, rung span, leg spacing,
shield radius, air padding); everything else is inherited.

## Geometry

Lengths in metres unless marked. "Record" rows are version-tagged measured
values ((1*) licence), not bands.

| Quantity | F-clinical-1.5T | F-clinical-3T |
|---|---|---|
| frequency | 64 MHz | 128 MHz |
| legs | 32 | 32 |
| ring-gap sheets (port census) | 64 | 64 |
| rung (conductor) diameter | 0.710 m | 0.704 m |
| RF shield diameter | 0.745 m | 0.743 m |
| shield (can) length | 0.954 m | 0.700 m |
| rung length (ring-centreline separation) | 0.450 m | 0.516 m |
| leg span | 0.5143 m | 0.5800 m |
| ring gap (arc) | 17.14 mm | 17.14 mm |
| leg width / ring tube radius | 12 mm / 4 mm | 12 mm / 4 mm |
| leg-box section | 4 × 4 × 21.43 mm | 4 × 4 × 21.43 mm |
| leg-box corner clearance to the can | 5.685 mm | 7.685 mm |
| load | ASTM F2182 box 650 × 420 × 90 mm, 650 mm on the axis, centred (C2); C1 stand-in cylinder 0.10969 m × 0.65 m, volume-matched | ASTM F2182 box as at 1.5 T (C8 box window); C8 r1 / r2 on the same stand-in cylinder |
| sizing rungs (resolution, conductor_resolution r1 / r2) | 0.015 m, 1.6 mm / 2.4 mm | 0.015 m, 1.6 mm / 2.4 mm |
| cells r1 / r2 (records) | 1 147 366 / 912 804 | 1 057 895 / 808 126 |
| U₁ / U₂ at r1 (records) | 1 525 193 / 7 885 674 | 1 425 765 / 7 335 988 |
| gmsh at r1 (measured) | 295.1 s, `-n 2`, 2026-10-15 12:07 CDT (`20261015T170721Z_GEO-38-s1a-W1-r1-a2.log`) | 286.5 s, `-n 2`, 2026-10-18 12:04 CDT (`20261018T170415Z_GEO-38-C8-W1-r1.log`) |

Notes on the rows:

- **Shield (can) length** is not a typed parameter. The generator places the
  can at ±z_extent with z_extent = max(coil_length/2, leg_spacing/2 +
  ring_minor_radius, phantom_height/2) + air_padding; the 3 T table computes
  its `air_padding` so that the can is the sourced 0.70 m.
- **Rung length at 3 T** (review ruling (b), 2026-10-18 03:00): the sourced
  0.58 m is the **leg span** (`coil_length`), asserted at 1e-9; the
  ring-centreline separation keeps the 1.5 T absolute leg overhang
  (0.0321428571 m past each ring), so it reads 0.516 m. The 1.5 T convention
  (sourced 0.45 m = ring-centreline separation, leg span 0.5143 m) is not
  carried over.
- **Box cells.** With the ASTM box in place: 1.5 T 1 149 150 cells at
  `-n 1 / 2 / 4` (C2, `20261015T185056Z_GEO-38-s2-W3-n4.log`); 3 T
  1 061 517 cells at `-n 1 / 2 / 4` (C8b `-n 1`,
  `20261019T172435Z_GEO-38-C8b-W1-box-n1.log`; C8b `-n 2`,
  `20261019T173044Z_GEO-38-C8b-W2-box-n2.log`; C8 `-n 4`,
  `20261018T171406Z_GEO-38-C8-W3-box-n4.log`).
- **Can-cap clearance at 3 T** (review ruling (c)): the box ends sit 25 mm
  from the can caps; the closed can stays, and the clearance is named for the
  weekly.

## Sources

The three sources found 2026-10-14 (`docs/planning/chunks/GEO-38.md`,
"Nominal geometry"):

- https://mri-q.com/uploads/3/4/5/7/34572113/body_coil_distriubtion_mrm.20177.pdf
  — "Efficient high-frequency body coil for high-field MRI" (MRM): the 3 T
  high-pass quadrature body coil, 58 cm long, 70.4 cm in diameter, inside a
  shield 74.3 cm in diameter and 70 cm long. Supplies the 3 T rung diameter,
  shield diameter, shield length and rung length.
- https://pmc.ncbi.nlm.nih.gov/articles/PMC12137786/ — a 70 cm-bore 3 T
  scanner modelled with a 32-rung high-pass birdcage. Supplies the leg count
  (32) and the high-pass topology.
- https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/11442125
  — US patent 11 442 125: 1.5 T rungs 45 cm long on a carrier of 65 cm inner
  diameter, the RF shield within ≈ 2 cm of the end-ring diameter. Supplies
  the 1.5 T rung length and the shield-to-ring spacing.

**Modelling choices, not sourced:**

- the 1.5 T rung diameter 0.710 m, set just outside the 0.70 m patient space
  rather than at the patent's 65 cm carrier;
- the 1.5 T shield length, which follows from the air-padding rule
  (padding = (0.03/0.07) · ring_radius) rather than from a source;
- the leg overhang past each ring (0.0321428571 m);
- the leg-box section (4 × 4 × 21.43 mm);
- the closed PEC can itself (a modelling stand-in for the RF shield);
- the 3 T rung-length reading (0.58 m as the leg span, ruling (b)).

The weekly may refine a value only with a cited source (PROJECT_PLAN §10).

## The rulings

- **The human operator's port-census ruling (2026-10-15):** 32 legs ⇒ every
  ring gap of both rings is a lumped sheet, 64 sheets and a 64×64 untuned
  S-matrix; tuning and SAR by two ports in quadrature with the 64 gaps as
  series-RC sheets at the 0.2 Ω ESR; 16 legs considered and rejected because
  it contradicts the cited 32-rung source.
- **The 2026-10-15 03:00 review's provisional choices:** rung length = the
  ring-centreline separation; the shield as a closed PEC can (new
  keyword-only `shield_radius`); F-human's absolute conductor sizing; C1's
  load a volume-matched gel cylinder; the leg-box dy cut to 4 mm.
- **The 2026-10-15 10:30 corner ruling:** leg-box dx 4 mm; the design corner
  clearance ≥ 5 mm; the `shield_radius` refusal computes the leg-box corner
  radius over the actual box azimuths; `shield_radius` is refused together
  with `leg_gap_length`.
- **All of the above ratified** by the 2026-10-17 weekly.
- **The 2026-10-15 18:00 ruling** concerns the build cache for the C1 (b)
  column (`-n 4` build, column alone at `-n 8`), not geometry.
- **Item 518's rulings (2026-10-18 03:00 review):** (b) the 3 T rung length
  is the weekly's 0.58 m, asserted at 1e-9 on the rung (leg) span, the chosen
  parameter (`coil_length`) printed; (c) the closed PEC can stays and the
  25 mm cap-to-box clearance is printed and named for the weekly.

## How it is built

The 1.5 T table is the single source of truth; the 3 T table is a derived
dict. Both are passed to `MeshGenerator.birdcage_port_domain`, whose
keyword-only `shield_radius` adds the closed PEC can and whose `phantom_box`
places the ASTM box as tag 3.

```text
F_CLINICAL_15T_KWARGS   tests/validation/test_geo38_f_clinical_build.py
F_CLINICAL_3T_KWARGS    tests/validation/test_geo38_f_clinical_3t.py
MeshGenerator.birdcage_port_domain(..., shield_radius=..., phantom_box=...)
gate modules            tests/validation/test_geo38_f_clinical_build.py   (C1 (a))
                        tests/validation/test_geo38_f_clinical_box.py     (C2)
                        tests/validation/test_geo38_f_clinical_3t.py      (C8)
```

The meshing example that builds the 1.5 T fixture with the box is
[the `mesh:15` guide](../../examples/meshing/15_f_clinical_coil_astm_box.md).

## Clinical heating records (entries 59 / 60)

The two `xl` windows `20261110T080008Z_GEO-38-P33.log` (1.5 T, 64 MHz) and
`20261111T080008Z_GEO-38-P33-3T.log` (3 T, 128 MHz) ran the P61 body on the
unwalled ASTM box: two-port quadrature drive at the degree-1 tuning node
(C = 193.5996 pF at 1.5 T, 44.0525 pF at 3 T), both senses (ccw and cw), and
the temperature rise θ(t\*) at t\* = 2821.500 s, **per accepted watt** (P_acc
at the lumped line cut, the 2026-10-12 ruling). c\* is the hotspot point in
metres (fixture frame); P1 and P2 are the two box-wall probes. Temperatures
are in K per accepted watt.

| Field, sense | c\* | θ(t\*) at c\* | P1 +y wall | P2 -y wall |
|---|---|---|---|---|
| 1.5 T, ccw | (-0.105, 0.0, -0.05) | 1.611905088e-03 | +4.445002899e-03 | +4.110246040e-03 |
| 1.5 T, cw | (-0.105, 0.0, -0.05) | 1.266521350e-03 | +3.070816668e-03 | +3.262978681e-03 |
| 3 T, ccw | (0.105, 0.0, 0.05) | 6.516151132e-03 | +1.272564088e-02 | +1.213484346e-02 |
| 3 T, cw | (0.105, 0.0, 0.05) | 5.749650374e-03 | +1.042249186e-02 | +1.076902943e-02 |
| 3 T / 1.5 T, ccw | each its own | 4.042515e+00 | not read | not read |
| 3 T / 1.5 T, cw | each its own | 4.539718e+00 | not read | not read |

The two ratio rows are the 3 T / 1.5 T quotient of θ(t\*) at each field's
own c\*, as `GEO-38` C39 printed them in
`20261111T103024Z_GEO-38-C39-W1.log`.

- **degree 1.** Both records are degree-1 fields; no degree-2 or h-converged
  heating figure exists at clinical scale.
- **1.5 T gel constants.** The 3 T window uses the 1.5 T values of σ, ε_r and
  the thermal constants (provisional); the ratio is not a measured
  frequency dependence of the gel.
- **20 °C baseline.** The window also prints a `cem43` per accepted watt; it
  is the unheated 20 °C baseline, not a dose per watt (a thermal dose is
  nonlinear in power), so no `cem43` figure is quoted here.
- **No container wall.** The box is unwalled; the container-wall thermal
  model is not in these records.
- **No compliance claim.** These are internal records with no external
  check; nothing here is an ASTM F2182 or safety statement.

## Closed-form gel-box records (MAT-21 steps 48 / 49)

The numbers below are printed by `20261112T120219Z_MAT-21-P88-W1.log`
(`MAT-21` step 48), `20261112T182020Z_MAT-21-P89-W1b.log` (`MAT-21` step 49)
and `20261110T103343Z_MAT-21-P83-W1.log` (the power share), beside the three
ratios `GEO-38` C39 printed in `20261111T103024Z_GEO-38-C39-W1.log`.

Step 48 heats the gel ball by a deposited power P_dep and compares the
thermal-dose excess over the unheated baseline with its linear term:

| P_dep (W) | excess / D0 | excess / linear term − 1 |
|---|---|---|
| 0.1 | 2.38132e-04 | 1.6194e-04 |
| 1 | 2.38479e-03 | 1.6212e-03 |
| 10 | 2.41996e-02 | 1.6392e-02 |
| 100 | 2.81864e-01 | 1.8384e-01 |

The excess per watt at 100 W is 1.1836 times that at 0.1 W, so no single
number per watt carries it.

Step 49 places the 18-ball lattice at C6's points, weighted by entries
52 / 58's SAR10g, and predicts the 3 T / 1.5 T lever as the power share times
SP, the ratio of θ at each field's own c\* per deposited watt:

| Quantity | Value |
|---|---|
| SP | 1.1041889 |
| power share (3 T / 1.5 T) | 3.6935 |
| predicted lever = share × SP | 4.078321 |
| C39 ccw point | 4.042515e+00 |
| C39 ccw ball | 4.030354e+00 |
| C39 cw point (cw heated at ccw's c\*) | 4.539718e+00 |

Step 49's printed verdict reads "lattice resolves the excess"; the planning
review reads it. The question it answers is set out in
[the `mat:42` example](../../examples/materials/42_f_clinical_box_heating_field_strength_lever.md).

- **Closed form.** No FEM field enters either record.
- **Not the clinical SAR.** The 18-ball lattice at C6's points, weighted by
  entries 52 / 58's SAR10g, is a stand-in, not a fit.
- **Per deposited watt.** Both records are per deposited watt of the ball or
  lattice source, never per accepted watt; P_dep, P_box and P_acc are not
  mixed.
- **1.5 T gel constants.** Both records use the 1.5 T gel values.
- **No dose per watt.** Nothing here is a compliance figure, a dose band, a
  damage threshold or a tissue model; 100 W is not a realistic drive.
- **A record.** Provisional: a planning review reads the verdict, and no band
  is set here.

## Rod-in-box records (entries 48 / 61)

The two `xl` windows `20261102T080008Z_GEO-36-step9-I36.log` (1.5 T, 64 MHz)
and `20261112T080008Z_GEO-36-step11-I36-3T.log` (3 T, 128 MHz) solved the
F-clinical coil around the ASTM box with the 30 mm rod (arm M, σ_wire
800 S/m) and with its gel twin (arm R), under the two-port quadrature drive
(1, +1j) at each field's own in-model degree-1 node (C = 193.5996 pF at
1.5 T, 44.0525 pF at 3 T). Q = SAR_1g(M) / SAR_1g(R) on the four tip balls;
η_rod(T1) = |E_M| / |E_R| at the tip point T1; the rod share is arm M's rod
power over its volume power.

| Quantity | 1.5 T (entry 48) | 3 T (entry 61) | 3 T / 1.5 T |
|---|---|---|---|
| Q top | 4.460568 | 4.854396 | 1.088291e+00 |
| Q bottom | 4.428040 | 4.884059 | 1.102984e+00 |
| Q top_apex | 4.413534 | 4.748393 | 1.075871e+00 |
| Q bottom_apex | 4.341688 | 4.775105 | 1.099827e+00 |
| η_rod(T1) | 3.642324 | 3.407173 | 9.354393e-01 |
| rod share (arm M) | 3.740924e-05 | 7.434370e-05 | 1.987308e+00 |

The ratio column is the 3 T / 1.5 T quotient as `GEO-36` C40 printed it in
`20261112T104117Z_GEO-36-C40-W1-n1.log`. The gel twin's rod share is a
numerical floor and is never quoted as a ratio; Im Y_in is not quoted (its
sign differs between the two nodes).

- **degree 1.** Both records are degree-1 fields; no degree-2 or h-converged
  rod-in-box figure exists at clinical scale.
- **1.5 T gel constants.** The 3 T window uses the 1.5 T values of σ and
  ε_r for the gel (provisional); the ratio column is not a measured frequency
  dependence of the gel.
- **Box containment is not asserted.** Neither window gates the field's
  containment in the box.
- **No published comparand.** No external or published value exists for
  these quantities on this fixture.
- **A record.** Entry 61's decision rule is read by a planning review, not on
  this page; no band is set here.

## Not claimed

Beyond the heating records and the rod-in-box records (internal, degree 1)
and the closed-form records (internal, closed form) above: no
port gate, S-matrix, tuning, resonance, B₁⁺ or SAR figure is quoted on this
page; not a scanner model; no ASTM compliance claim; the closed can is not claimed to
be the physical shield.

## Erratum

Resolved. The `mesh:15` guide
(`examples/meshing/15_f_clinical_coil_astm_box.md`) gave the can as "length
0.745 m", which was its diameter. §9 item 543 (an `EX-117` follow-on)
corrected it to **0.954285714 m** (2·z_extent,
`20261015T171259Z_GEO-38-s1a-W2-r2-a2.log`), together with the stale
"pending the weekly's ratification" text for choices the 2026-10-17 weekly
had ratified. The table above always carried the right value.
