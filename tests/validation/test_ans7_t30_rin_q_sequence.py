"""ANS-7 T30 (§9 item 364): the three-rung Cu R_in / Q / f0 h-sequence at C_F
from logged data (lane-1 revalidation, no solve).

Re-derives each rung's five-node Route B Cu rows (degree 2, κ off, exact
rectangles, F-small, 64 MHz band, at ``C_F`` -- untuned, not the C_in rows)
from the stored node JSONs through the production ``sparse_rows`` path,
asserts they reproduce the logged digits within ``NODE_S_RTOL`` (imported),
and prints the successive changes Δ₁ = v(0.01) − v(×1), Δ₂ = v(0.005) − v(0.01)
and their ratio ρ = Δ₂/Δ₁ against the §10 T30 rule (|Δ₂| ≤ ½|Δ₁| → bracket
[x(0.005) − |Δ₂|, x(0.005) + |Δ₂|], entry 34's form, 364-Q1 ruling; else
"unconverged at h = 0.005 -- F31 is the next rung").  No band, tolerance or
target is set here.  Brief: docs/planning/chunks/ANS-7.md "### §9 item 364".
"""

import pytest

from tests.complex_mode import complex_only
from tests.validation.test_ans7_exact_strip import (
    ANS7_EXACT_STRIP_ROUTE_B_ROWS_P2,
    DEFAULT_NODES_MHZ,
    ROWS,
    STEP6_RUNG_ROWS,
    _assert_tie,
    _load_nodes,
    _say,
    _tie_from_json,
    _uj,
    _x1_cu5_rows,
)
from tests.validation.test_ans7_sparse_route_b import (
    ASSEMBLY5_NODES_MHZ,
    NODE_S_RTOL,
    sparse_rows,
)

RUNGS = ("x1", "h0.01", "h0.005")
RUNG_RES = {"x1": None, "h0.01": 0.01, "h0.005": 0.005}

# Records, not bands: the logged Cu five-node Route B rows at C_F (full repr).
#   x1:     20261004T112945Z_ANS-7-item214-W3-rung.log:2259-2261 ((S5) 5-node;
#           also 20261008T070008Z_ANS-7-step4p.log:4171-4173)
#   h0.01:  20261004T112945Z_ANS-7-item214-W3-rung.log:2218
#   h0.005: 20261008T070008Z_ANS-7-step4p.log:4128 (xl entry, T12')
T30_LOGGED_ROWS = {
    "x1": {
        "f0_hz": 62701081.939501435,
        "r_in_ohm": 18441.75458382276,
        "q_loaded": 977.2613692833095,
    },
    "h0.01": {
        "f0_hz": 62726846.66023456,
        "r_in_ohm": 18313.602153804146,
        "q_loaded": 970.4389973340836,
    },
    "h0.005": {
        "f0_hz": 62751784.313079365,
        "r_in_ohm": 18223.232609192608,
        "q_loaded": 965.5456118046693,
    },
}
# Cell counts: x1 item 196 (80 798); h0.01 item 214 W3 log; h0.005 step4p.log
# and docs/testing/xl-ledger.md row "ANS-7-step4p".
T30_CELLS = {"x1": 80798, "h0.01": 140270, "h0.005": 523880}
# Resolution labels; x1 = 0.015 is test_ans7_exact_strip.py's own label.  Used
# only for the printed order estimate.
T30_H_M = (0.015, 0.01, 0.005)


def _rung_inputs():
    """{rung: (recs, nodes_mhz)} for the JSONs each rung's rows came from."""
    out, missing = {}, []
    rows_x1, recs_x1 = _x1_cu5_rows()
    if rows_x1 is None:
        missing += list(recs_x1)
    else:
        out["x1"] = (recs_x1, DEFAULT_NODES_MHZ)
    for rung in RUNGS[1:]:
        recs, miss = _load_nodes(ASSEMBLY5_NODES_MHZ, 2, RUNG_RES[rung], ("cu",))
        if recs is None:
            missing += miss
        else:
            out[rung] = (recs, ASSEMBLY5_NODES_MHZ)
    return out, missing, rows_x1


@pytest.fixture(scope="module")
def t30():
    inputs, missing, rows_x1 = _rung_inputs()
    if missing:
        return {"missing": missing}
    rows = {"x1": rows_x1["cu"]["rows"]}
    for rung in RUNGS[1:]:
        recs = inputs[rung][0]
        rows[rung] = sparse_rows(
            ASSEMBLY5_NODES_MHZ, {"cu": [_uj(r["z"]["cu"]) for r in recs]}, columns=("cu",)
        )["cu"]["rows"]
    return {"missing": [], "inputs": inputs, "rows": rows}


def _need_inputs(t):
    if t["missing"]:
        pytest.fail(f"T30 inputs absent (BLOCKED; never re-solve): {t['missing']}")


@complex_only
def test_t30_inputs(t30):
    """(F): degree 2, κ None, the rung's resolution, columns, one cell count equal
    to T30_CELLS, the tie per file, f_hz within 1 Hz of its node."""
    _need_inputs(t30)
    for rung in RUNGS:
        recs, nodes = t30["inputs"][rung]
        res = RUNG_RES[rung]
        assert len(recs) == len(nodes), (rung, len(recs))
        assert all(r["degree"] == 2 for r in recs), (rung, [r["degree"] for r in recs])
        assert all(r["kappa"] is None for r in recs), rung
        assert {r.get("resolution") for r in recs} == {res}, (
            rung,
            [r.get("resolution") for r in recs],
        )
        if res is not None:
            assert all(list(r["columns"]) == ["cu"] for r in recs), (
                rung,
                [r["columns"] for r in recs],
            )
        cells = {r["cells"] for r in recs}
        assert cells == {T30_CELLS[rung]}, (rung, cells)
        for m, r in zip(nodes, recs):
            _assert_tie(r["cells"], _tie_from_json(r["tie"]), res)
            assert abs(r["f_hz"] - m * 1e6) <= 1.0, (rung, m, r["f_hz"])
        _say(
            f"[t30] (F) {rung}: {len(recs)} JSONs, degree 2, kappa None, resolution {res}, "
            f"cells {T30_CELLS[rung]}, nodes {tuple(nodes)} MHz -- ok"
        )


@complex_only
def test_t30_reproduces_logged_rows(t30):
    """(R-rep): each of the nine (rung, quantity) rows reproduces the logged
    value within NODE_S_RTOL (predicted exactly 0)."""
    _need_inputs(t30)
    fails = []
    for rung in RUNGS:
        for k in ROWS:
            got = t30["rows"][rung][k]
            want = T30_LOGGED_ROWS[rung][k]
            d = abs(got / want - 1.0)
            _say(
                f"[t30] (R-rep) {rung} {k}: {got!r} vs logged {want!r}: rel {d:.3e} "
                f"(ASSERTED <= NODE_S_RTOL {NODE_S_RTOL:.0e}; PREDICTED 0: "
                f"{'met' if d == 0.0 else 'NOT met'})"
            )
            if d > NODE_S_RTOL:
                fails.append((rung, k, got, want, d))
    for k in ROWS:
        rec7 = ANS7_EXACT_STRIP_ROUTE_B_ROWS_P2["cu"][k]
        v5 = t30["rows"]["x1"][k]
        _say(
            f"[t30] printed only: x1 7-node record {k} {rec7!r} vs 5-node {v5!r}: rel "
            f"{abs(rec7 / v5 - 1.0):.3e} (PREDICTED 9.2e-10 R_in / 6.6e-10 Q, as (S5))"
        )
        s6 = STEP6_RUNG_ROWS[k]
        _say(f"[t30] printed only: STEP6_RUNG_ROWS {k} {s6!r} vs h0.01 {t30['rows']['h0.01'][k]!r}")
    assert not fails, fails


@complex_only
def test_t30_rung_swap_control(t30):
    """(N-swap), asserted: each rung's recomputed R_in and Q against the other
    two rungs' logged values separates by > 100 x NODE_S_RTOL."""
    _need_inputs(t30)
    floor = 100.0 * NODE_S_RTOL
    fails = []
    for rung in RUNGS:
        for other in RUNGS:
            if other == rung:
                continue
            for k in ("r_in_ohm", "q_loaded"):
                d = abs(t30["rows"][rung][k] / T30_LOGGED_ROWS[other][k] - 1.0)
                _say(
                    f"[t30] (N-swap) recomputed {rung} {k} vs logged {other}: rel {d:.4e} "
                    f"(ASSERTED > {floor:.0e})"
                )
                if not d > floor:
                    fails.append((rung, other, k, d))
    assert not fails, fails


def _order(rho, h=T30_H_M, lo=0.1, hi=6.0):
    """p with (h1^p - h2^p)/(h0^p - h1^p) = rho, bisection on [lo, hi]; None if
    no sign change."""

    def g(p):
        return (h[1] ** p - h[2] ** p) / (h[0] ** p - h[1] ** p) - rho

    a, b = lo, hi
    ga, gb = g(a), g(b)
    if ga * gb > 0:
        return None
    for _ in range(200):
        c = 0.5 * (a + b)
        gc = g(c)
        if ga * gc <= 0:
            b, gb = c, gc
        else:
            a, ga = c, gc
    return 0.5 * (a + b)


@complex_only
def test_t30_sequence_record(t30):
    """(A) asserted: Δ₁ + Δ₂ = v(0.005) − v(×1) within NODE_S_RTOL relative
    (bookkeeping).  Printed record: Δ₁, Δ₂, ρ, the §10 branch, the order."""
    _need_inputs(t30)
    preds = {
        "r_in_ohm": (-128.152430, -90.369545, 0.70517, 1.68),
        "q_loaded": (-6.822372, -4.893386, 0.71726, 1.64),
        "f0_hz": (25764.72, 24937.65, 0.96790, 1.06),
    }
    fails = []
    for k in ("r_in_ohm", "q_loaded", "f0_hz"):
        v0, v1, v2 = (t30["rows"][r][k] for r in RUNGS)
        d1, d2 = v1 - v0, v2 - v1
        tot = v2 - v0
        ident = abs((d1 + d2) - tot) / abs(tot)
        rho = d2 / d1
        p = _order(rho)
        pd1, pd2, prho, pp = preds[k]
        _say(
            f"[t30] (A) {k}: Delta1 + Delta2 - (v(0.005) - v(x1)) rel {ident:.3e} "
            f"(ASSERTED <= NODE_S_RTOL {NODE_S_RTOL:.0e})"
        )
        _say(
            f"[t30] RECORD {k}: v(x1) {v0!r}, v(0.01) {v1!r}, v(0.005) {v2!r}; "
            f"Delta1 {d1!r}, Delta2 {d2!r}, rho = Delta2/Delta1 {rho!r} "
            f"(PREDICTED Delta1 {pd1}, Delta2 {pd2}, rho {prho}); relative steps "
            f"{d1 / v0:+.4%}, {d2 / v1:+.4%}"
        )
        _say(
            f"[t30] order p (printed; not the §10 rule) on h {T30_H_M}: "
            f"{'n/a (no root in [0.1, 6])' if p is None else f'{p:.4f}'} (PREDICTED ~{pp})"
        )
        if k == "f0_hz":
            _say("[t30] f0 rho printed for comparison with the T13 ruling (first order in h)")
        elif abs(d2) <= 0.5 * abs(d1):
            _say(
                f"[t30] BRANCH {k}: |Delta2| <= 0.5 |Delta1| -> bracket (entry 34's form, "
                f"364-Q1) [{v2 - abs(d2)!r}, {v2 + abs(d2)!r}] (PREDICTED branch: unconverged "
                f"-> NOT met)"
            )
        else:
            _say(
                f"[t30] BRANCH {k}: |Delta2|/|Delta1| = {abs(d2) / abs(d1):.5f} > 0.5 -> "
                f"unconverged at h = 0.005 -- F31 is the next rung (PREDICTED branch: met)"
            )
        if ident > NODE_S_RTOL:
            fails.append((k, ident))
    assert not fails, fails
