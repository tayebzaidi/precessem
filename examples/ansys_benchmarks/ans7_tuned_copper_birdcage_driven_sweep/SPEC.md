# ANS-7 — The tuned copper four-leg birdcage: a capacitor on every gap, the drive across one of them, `Z_in(f)` across the k = 1 resonance

**Commissioned 2026-09-26 (weekly planning review)** — PROJECT_PLAN §10
chain T10. **Revised 2026-09-28 (v2, operator directive, interactive
session) before any AED number landed** — see "Revision history" at the end.
Sixth Ansys benchmark case, on the physics gated by `PORT-22` ✅ (the tuned
resonance located by a driven sweep, circuit layer = field solve within 1e-3),
`PORT-14` / `PORT-15` ✅ (lumped RLC sheets and the circuit layer), `TH-14` ✅
(copper by the Leontovich surface impedance; externally checked by `ANS-6`)
and `TH-17` step 3 (the four-capacitor coil tuned by the k = 1 block, the
eigenpair at 64.33 MHz).

**Why this case.** Every earlier AED case terminates the undriven ports in
50 Ω. A birdcage is built with a tuning capacitor across **every** gap and
driven across **one of those capacitors** (through a matching network the
field comparison does not need); its resonance is the coil's own k = 1 mode,
and the drive sees it as a **parallel** resonance — a peak in `R_in`. This
case fixes that capacitor and asks both codes where the k = 1 resonance is,
how sharp it is and what the drive sees — agreement on the capacitor a user
would have to fit, in the configuration a user builds. **The capacitor is
fixed, not tuned, on the AED side**: no optimisation loop.

## Status

- [ ] Runnable half (`ANS-7`, §7) — re-registered for v2: our side measures
      the same configuration (all four gaps at `C`, the drive across C₁)
      through `PORT-15`'s circuit layer on the four-port and through the
      in-model route, and confirms the k = 1 resonance sits inside the fine
      window below.
- [x] Operator replication in AED — v2, **landed 2026-09-29** (results in the
      gitignored `aed_results/`). The first v2 upload (09-28) was discarded:
      the PyAEDT call stored the Lumped RLC capacitance without a unit and
      AEDT treated the three Route A capacitors as open gaps (Route A's bare
      port matched Route B's four-port with gaps 2–4 open). The script now
      writes the unit explicitly and guards against that failure mode; the
      09-29 rerun, from scratch, passes both guards in every Route A design.
      Both routes, both columns, Zero and First Order Route A, as SPEC'd.
      Interactive first-pass check (2026-09-29): SpecGrid / FineWindow grids
      exact; rows (a)–(c) recompute from the exported `Z_in`; Route B's Z/S
      and reciprocity identities hold to round-off; the k = 1 resonance is
      present in every design. **For the adjudication:** the two Zero Order
      Route A designs did not meet the ΔS target within the pass cap (the
      Route B and First Order designs did), so AED's own Zero Order route
      A − B difference carries unconverged mesh error; AED's adaptive meshes
      are several times finer than ours at Zero Order.
- [x] Adjudication — **2026-09-30 weekly review** (PROJECT_PLAN §10 rulings;
      known-issues 2026-09-30; numbers in the gitignored
      `docs/private/ans7-ans6-adjudication-2026-09-30.md`). **AGREE** on the
      mode (one k = ±1 parallel resonance in the window in both codes, share
      ≥ 0.9999), on Route A = Route B, and on the copper contribution (the
      conductor's part of 1/Q within 5 %, `Q_Cu/Q_PEC` within 1.1–2.4 %, the
      Cu − PEC shift of `f₀` within 5 % of First Order). **INCONCLUSIVE** on
      the absolute rows: our `f₀` is 1.2–1.4 % above the order-matched
      columns and 1.9 % above First Order (row (d): 2.5–3.9 % in C), Q 6–10 %
      and `R_in` 7–11 % low — our side is one degree-1 rung. Re-read when
      chain T11 / T12 supply the order- and h-matched rung. No tuning or
      match claim is licensed by this case yet (re-read at degree 2 below).
- [x] Re-adjudication at degree 2 on exact rectangles, κ off — **`ANS-7` ✅ (2026-10-17 weekly, T41).** f₀ AGREE with its h-statement (T13, 2026-10-10: within 0.11 % of AED First Order at ×1, 0.03 % at h = 0.005, first order in h). Q and `R_in` converge in h (p = 1.6451 / 1.6786 on h = 0.015 / 0.01 / 0.005 m) to a level ≈ 3 % below AED First Order: INCONCLUSIVE as an absolute Q, pending AED's own pass history (T48). Rungs, rates and limits: `COMPARISON.md` § Verdict.
- [x] The runnable half's box above: done 2026-09-29 (§7 `ANS-7`).

## Geometry (SI units)

**Identical to `ANS-6`** —
`../ans6_copper_birdcage_four_port_10_64_128MHz/SPEC.md` § Geometry (itself
`ANS-4`'s): box `x, y ∈ [−0.120, +0.120]`, `z ∈ [−0.100, +0.100]`; two
end-ring tori, four gapped legs, phantom cylinder radius 0.030,
`z ∈ [−0.040, +0.040]`.

## Materials and the two coil columns

Identical to `ANS-6`: phantom σ = **0.5** S/m, εᵣ = **78**; air elsewhere;
coil **column Cu** — Finite Conductivity boundary, σ = **5.8e7** S/m, Solve
Inside off; coil **column PEC** — Perfect E, Solve Inside off. All six box
faces Perfect E.

## Capacitors and the drive

The four gap sheets are `ANS-4`'s, unchanged (interior-half sheet, width
0.007 × height 0.008, current / integration line along **+ẑ**).

**Every gap carries the same capacitor, C = 66.29 pF exactly
(6.629e-11 F)**, ideal (no ESR — see Out of scope). The value is `TH-17`
step 3's `C_mode` (the k = 1 block tuned to 64 MHz on the F-small hole mesh,
`6.628769508e-11 F`), rounded; the eigenpair there sits at ≈ 64.33 MHz.

**The drive is across C₁** (the leg at azimuth 0°), in parallel with it:
`Z_in` is the impedance seen at the gap-1 terminals **with C₁ in place**.

Two routes, each run by both codes; their agreement is a check each code
owes itself:

| Route | Model | `Z_in` |
|---|---|---|
| **A — in-model** (primary) | Gap 1: a 50 Ω **lumped port** (bare). Gaps 2–4: **Lumped RLC** boundaries, capacitance only, C. | `Z_in = (1/Z₁₁ + jωC)⁻¹` — C₁ added across the port terminals exactly. The port and C₁ share one gap, so the field is the same as with C₁ in the model (an ideal element across the terminals changes only how the source current divides) and the adaptive mesh is refined on the resonant structure. |
| **B — four-port + circuit** (secondary) | All four gaps **50 Ω lumped ports** (`ANS-6`'s four-port). | Ports 2–4 terminated in C and C₁ placed across port 1, both in the circuit domain: `Z_bare = Z₁₁ − Z₁ₜ(Z_tt + Z_C I)⁻¹Zₜ₁`, `Z_in = (1/Z_bare + jωC)⁻¹`, `Z_C = 1/(jωC)`. The termination may be done in an AED Circuit design or, as the operator's script does, in NumPy on the exported 4×4 (say which). |

## Frequencies

Discrete sweeps (every point solved, no interpolating or fast sweep), both on
the mesh adapted **at 64 MHz** to ΔS ≤ 0.002, seeded with ≥ 2 elements across
each gap along ẑ:

- **SpecGrid:** 44 → 84 MHz in 2 MHz steps (21 points) — the overview. It
  contains the k = 1 resonance (≈ 64.3 MHz) and, on our side, the k = 2
  resonance near 78 MHz; say what you see.
- **FineWindow:** 62.0 → 67.0 MHz in 0.1 MHz steps (51 points) — the
  adjudication rows are read here.

## Basis order and mesh

Per the `ANS-5` ruling, Route A runs each coil column at **Zero Order**
(= our degree 1) and **First Order** (= degree 2). Route B at Zero Order
(First Order optional). Mixed Order forbidden. Confirm unknowns-per-tet from
the matrix statistics.

## Quantities to export (all digits AED prints; do not round)

| Quantity | Definition | Units |
|---|---|---|
| `Z_in(f)` both sweeps, per route, column, order | complex, C₁ across the drive | Ω |
| `Z₁₁(f)` of the bare port (Route A) and the 4×4 `Z`, `S` (Route B) | raw, before the C₁ / termination algebra | Ω, — |
| Loss partition `P_coil / (P_coil + P_phantom)` (Cu) and the two integrals | Route A, gap-1 port driven, at 64 MHz and at the SpecGrid point nearest `f₀` | —, W |
| Gap-current k = ±1 share | Route B, at the fine-window point nearest `f₀`: the four gap currents with gap 1 driven, DFT over the legs, (\|ĝ₁\|² + \|ĝ₃\|²) / Σ\|ĝ_k\|² | — |
| Solve metadata | tets, passes, final ΔS, wall time, order, unknowns/tet, mesh re-use | — |

## Primary adjudication rows

Computed identically for both codes on the **FineWindow**, from
`Y_in = 1/Z_in` (a parallel resonance is a zero of `Im Y_in` crossing
− → +):

(a) **`f₀`** — the lowest − → + zero of `Im Y_in`, by linear interpolation
between the bracketing points (the k = 1 mode; confirmed by the Route B
k-share); (b) **`R_in(f₀)` = 1 / Re Y_in(f₀)**, `Re Y_in` interpolated the same
way; (c) **loaded Q** = `f₀ · (d Im Y_in / df) / (2 Re Y_in)` at `f₀`, bracket
secant; (d) the **capacitor equivalent of the frequency gap**
`ΔC/C = −2 Δf₀/f₀`. Secondary: rows (a)–(c) on the SpecGrid; Route A vs
Route B per code; the Cu − PEC difference of `f₀` and Q; the k = 2 resonance
if present.

## Reference values (ours)

None are written here: they are regenerated by the runnable half and appear
in its `COMPARISON.md`, never transcribed. For orientation only: `TH-17`
step 3's eigenpair at `C_mode` sits at ≈ 64.33 MHz on the PEC hole mesh, the
next mode near 78.4 MHz.

## Out of scope

No matching network (the unmatched `Z_in` is compared), no capacitor ESR or
series resistor (a real coil has both; the comparison needs the ideal
element — an ESR column is a possible v3), no quadrature drive, no
mode-spectrum claim beyond rows (a)–(c) and the k = 2 observation, no B₁⁺ or
SAR, no 16-leg or F-human layout (`ANS-8`), no optimisation on the AED side.

## Revision history

- **v1 (2026-09-26, weekly):** C = 15.57 pF on gaps 2–4 only, a **bare**
  50 Ω port on gap 1, rows on the zero of `Im Z_in` — a **series** resonance
  of a three-capacitor network (≈ 66.3 MHz on our side, `PORT-22` step 2 /
  §9 item 50). Well-posed as a code comparison, but not how a birdcage is
  built or driven: gap 1 carries no capacitor, the network is not
  C4-symmetric and carries no degenerate k = 1 pair (the 09-26 weekly's
  epitaph already noted that 15.57 pF on all four gaps puts the lowest mode
  near 129.6 MHz). **Withdrawn 2026-09-28** by the operator before any AED
  result landed; item 50's measurement stays in §7 as a record of v1.
- **v2 (2026-09-28, operator):** a capacitor on every gap at `TH-17`'s
  `C_mode`, the drive across C₁ (parallel resonance), Route A in-model +
  Route B four-port / circuit, fine window 62–67 MHz, rows on `Y_in`, the
  Route B k-share as mode identification.
