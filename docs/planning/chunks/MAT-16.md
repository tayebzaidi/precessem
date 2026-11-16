# MAT-16 — PROPOSED, not opened (for the weekly)

**Status: a proposal, with no §7 row and nothing queued.** The 2026-10-01 18:00 daily review drafted this lane-3 chunk to fill the restock shortfall and then declined to open it, for this reason: §9 rule (5) lets a daily open a Phase 7/8 chunk only "until the weekly has written the §10 subgoals and §7 chunks for them". The 2026-09-30 weekly has since written Phase 8 subgoals (8a)–(8e), and thermal dose is not among them. Daily-review step 5 also limits a daily to chunks *within* the current phase's subgoals. Opening, re-homing it (Phase 8, or Phase 9 (9a)'s "temperature via `thermal/`" results) or declining is the weekly's call. Nothing in §10 excludes dose: `MAT-13` / `MAT-14` / `EX-73` exclude it from their own scope only.

## Proposed §7 row (draft)

**CEM43 thermal dose as a post-processor on a temperature history — Phase 8.** `fem_em_solver/thermal/dose.py: cem43(times_s, temperatures_c, *, rule="trapezoid", r_above=0.5, r_below=0.25)` returns minutes:

CEM43 = (1/60)∫ R(T)^(43 − T) dt, with R = 0.5 for T ≥ 43 °C and 0.25 below (Sapareto & Dewey 1984, *Int. J. Radiat. Oncol. Biol. Phys.* 10:787).

`rule="interval_mean"` is the paper's discrete sum Σ Δtᵢ R^(43 − T̄ᵢ). The function is numpy only. The closed forms go into `utils/analytical.py` (additive).

**Done-when, step 1** (numpy, real build, `-n 1`). Write κ = −ln R and G(T) = expm1(κ(T)(T − 43))/κ(T). G is continuous at 43 °C, and G′ = R^(43 − T).
- **(C1) Constant T.** Five exact binary cases, and an unevenly sampled constant history, to ≤ 1e-13.
- **(C2) Ramp.** The ramp equals [G(T_end) − G(T₀)]/(60a), against the explicit two-piece split ≤ 1e-13 and against Gauss–Legendre ≤ 1e-12 (37 → 47 °C at 0.01 K/s: 37.26933 min). The cooling ramp equals the heating ramp.
- **(C3) Exponential approach.** Via Ei, against Gauss–Legendre ≤ 1e-11 on five cases. The late-time slope is asserted ≤ 1e-9.
- **(C4) Convergence.** The trapezoid rule on the 43 °C-crossing ramp has an LSQ rate in [1.9, 2.1], and so does `interval_mean`.
- **(K) The 43 °C kink.** The ramp from 43 − ε to 43 + ε reads 1 − (ln 2/4)ε within 1e-6 at ε = 1e-3.

Controls, each asserted:
- R swapped: misses by ×4 on the constant cases and ×8.2894 on the ramp.
- Kelvin input: misses by a log₂ factor of 273.15 / 275.15 / 273.194.
- The antiderivative without its 43 °C anchor: misses the ramp by exactly 1.2023 min.

**Out of scope:** tissue-damage thresholds and any biological or regulatory claim; Arrhenius damage; spatial dose maps; ASTM F2182.

**Step 2** (not drafted for queueing): CEM43 of `MAT-12`'s FEM centre history beside the closed-form series. Dose amplifies θ errors exponentially: `MAT-12`'s ≈ 3.5 % SAR offset on an 8 K rise becomes ≈ 21 % in CEM43 above 43 °C. Setting that band is the weekly's.

## Drafting notes kept for the weekly

- **Measured class.** Lane-3 closed form, median 8 min (n = 5). Predicted ≈ 35 min: a ≈ 1 min window, 15 fixed, ≈ 20 for code.
- **Traps the drafter found.**
  - Compute G with `np.expm1`. The naive `exp − 1` loses ≈ 1e-7 relative at ε = 1e-9, and (K) catches that.
  - Use one 43-anchored G on every piece. A per-piece `f/κ` is exactly the unanchored control.
  - `e^{−t/τ}` underflows past t/τ ≈ 745, so refuse a duration > 200τ.
  - Do not export from `thermal/__init__.py`. It imports dolfinx, and an export would force rule (c) windows on `MAT-9` … `MAT-14`.
  - Inputs are absolute °C, not θ.
- **Pre-derived numbers.** Trapezoid errors on N = 12·2^j are ≈ +2.88e-2 / 6.98e-3 / 1.80e-3 / 4.36e-4 / 1.12e-4. The kink-cell term makes the pairwise rates alternate (2.04 / 1.96), so only the LSQ slope is asserted. The exponential case (37 → T∞ 45, τ 300 s, 1800 s) is ≈ 72.45 min by hand.

**Closing note (2026-10-25).** The (K) kink and `expm1` hazards drafted here were read against `MAT-20`'s `cem43` (§9 item 638, `docs/planning/chunks/MAT-20.md#§9-item-638`): both held at `DOSE_BAR`; this draft is not queued.
