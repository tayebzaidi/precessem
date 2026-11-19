"""`ANS-7` lane 1 (§9 item 514, §10 T47): the public verdict section of ``ans:7``'s COMPARISON.md.

Render-only, on scratch copies of the committed ``metrics.json`` / ``COMPARISON.md``:

(X)  every verdict figure is traced to its source line (T30's log, PROJECT_PLAN.md, SPEC.md) and
     occurs in the committed verdict section;
(L)  the printed orders p are recomputed from T30's logged values by bisection (to 5e-5), and the
     Richardson limits / percentages formatted in the section equal the recomputation; the coarse
     pair gives the same limit to 1e-3 relative;
(R)  ``_write_comparison`` on the committed metrics reproduces the committed file byte for byte,
     the verdict heading once, after the header paragraph and before ``## Degree 2``;
(I)  stand-alone ``_write_comparison_esr128`` / ``_write_comparison_palace`` leave it unchanged;
(B)  no "Palace", no private path, no "unconverged", no numeric token with >= 8 significant digits;
(R0) the repo's files are untouched.
Negative controls (asserted): (N-p1/N-p2) Richardson at p = 1 / 2 changes the printed Q limit;
(N-cite) p = 1.6451 cited one line early fails (X); (N-digit) "963.25" -> "963.3" differs from the
render; (N-old) the pinned 248416c9 writer renders no verdict heading and differs.
"""

from __future__ import annotations

import gzip
import importlib.util
import json
import re
import shutil
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
CASE = REPO / "examples" / "ansys_benchmarks" / "ans7_tuned_copper_birdcage_driven_sweep"
SCRIPT = CASE / "07_tuned_copper_birdcage_driven_sweep.py"
PIN_SHA = "248416c9"
# Mirrored depth: the module body resolves ``parents[3]``.
OLD_SCRIPT = REPO / "logs" / "t47_pin" / SCRIPT.relative_to(REPO)
LOGS = REPO / "docs" / "testing" / "logs"
T30 = "20261011T094910Z_ANS-7-item364-T30.log"
PLAN = "PROJECT_PLAN.md"
SPEC = "SPEC.md"
H = (0.015, 0.01, 0.005)

# (cell as it occurs in the section, source, 1-based line, check, log needle or None)
# Checks: "v" the needle (default: the cell) occurs verbatim on the line.
TRACE = (
    ("80 798 cells", T30, 51, "v", "cells 80798"),
    ("140 270 cells", T30, 52, "v", "cells 140270"),
    ("523 880 cells", T30, 53, "v", "cells 523880"),
    ("h = 0.015 / 0.01 / 0.005 m", T30, 86, "v", "on h (0.015, 0.01, 0.005)"),
    ("p = 1.6451 (Q)", T30, 90, "v", "1.6451"),
    ("1.6786 (R_in)", T30, 86, "v", "1.6786"),
    ("p = 1.0625", T30, 94, "v", "1.0625"),
    ("BRANCH lines", T30, 87, "v", "BRANCH r_in_ohm"),
    ("BRANCH lines", T30, 91, "v", "BRANCH q_loaded"),
    (
        "within 0.11 % of AED First Order at ×1",
        PLAN,
        316,
        "v",
        "within 0.11 % of AED First Order at",
    ),
    ("0.03 % at h = 0.005", PLAN, 317, "v", "within 0.03 % at h = 0.005"),
    ("Q sits 3.1 %", PLAN, 4682, "v", "Q −3.1 %"),
    ("R_in 2.7 %", PLAN, 4682, "v", "`R_in` −2.7 %"),
    ("≈ 3 % below", PLAN, 4683, "v", "≈ 3 %"),
    ("INCONCLUSIVE", PLAN, 4686, "v", "INCONCLUSIVE"),
    ("T48", PLAN, 4687, "v", "T48"),
    ("share ≥ 0.9999", SPEC, 50, "v", "≥ 0.9999"),
    ("within 5 %", SPEC, 51, "v", "within 5 %"),
)
RECORD = re.compile(
    r"RECORD (\w+): v\(x1\) ([-+0-9.e]+), v\(0\.01\) ([-+0-9.e]+), "
    r"v\(0\.005\) ([-+0-9.e]+); Delta1 ([-+0-9.e]+), Delta2 ([-+0-9.e]+), "
    r"rho = Delta2/Delta1 ([-+0-9.e]+)"
)


# PROJECT_PLAN.md is edited by every review, so its rows are anchored: the recorded line (as at
# 2026-10-18, 316 / 4682) is re-located by a unique anchor string, the offset kept.
PLAN_ANCHORS = {
    316: ("exact rectangles at degree 2, `f₀` sits within 0.11 % of AED First Order at", 0),
    317: ("exact rectangles at degree 2, `f₀` sits within 0.11 % of AED First Order at", 1),
    4682: ("**Against AED:** T13 read Q −3.1 %", 0),
    4683: ("**Against AED:** T13 read Q −3.1 %", 1),
    4686: ("**Against AED:** T13 read Q −3.1 %", 4),
    4687: ("**Against AED:** T13 read Q −3.1 %", 5),
}


def _plan_line(line: int) -> int:
    anchor, off = PLAN_ANCHORS[line]
    ls = _lines(PLAN)
    hits = [i for i, s in enumerate(ls) if anchor in s]
    assert len(hits) == 1, (anchor, hits)
    return hits[0] + 1 + off


def _lines(name: str) -> list[str]:
    if name == PLAN:
        p = REPO / PLAN
    elif name == SPEC:
        p = CASE / SPEC
    else:
        p = LOGS / name
        if not p.exists():
            return gzip.open(str(p) + ".gz", "rt").read().split("\n")
    return p.read_text().split("\n")


def _holds(cell, src, line, mode, needle) -> bool:
    assert mode == "v"
    if src == PLAN:
        line = _plan_line(line)
    return (needle or cell) in _lines(src)[line - 1]


def _record(line: int):
    mt = RECORD.search(_lines(T30)[line - 1])
    assert mt, line
    return (
        mt.group(1),
        [float(mt.group(k)) for k in (2, 3, 4)],
        float(mt.group(5)),
        float(mt.group(6)),
        float(mt.group(7)),
    )


def _printed_p(line: int) -> float:
    mt = re.search(r"\): ([0-9.]+) \(PREDICTED", _lines(T30)[line - 1])
    assert mt, line
    return float(mt.group(1))


def _solve_p(rho: float) -> float:
    def g(p):
        return (H[1] ** p - H[2] ** p) / (H[0] ** p - H[1] ** p) - rho

    lo, hi = 0.1, 5.0
    assert g(lo) * g(hi) < 0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if g(lo) * g(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def _limit(v, p, pair=(1, 2)):
    a, b = pair
    ha, hb = H[a] ** p, H[b] ** p
    return v[b] - (v[a] - v[b]) * hb / (ha - hb)


def _sp(x: float, fmt: str) -> str:
    """Format and typeset a thousands space in the integer part (18441.75 -> 18 441.75)."""
    s = format(x, fmt)
    sign = s[0] if s[0] in "+-" else ""
    body = s[len(sign) :]
    ip, _, fp = body.partition(".")
    if len(ip) > 3:
        ip = ip[:-3] + " " + ip[-3:]
    return sign + ip + ("." + fp if fp else "")


def _pct(x: float) -> str:
    s = f"{x:+.2f}"
    return s.replace("-", "−")


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def ans7():
    return _load("ans7_writer_verdict_ut", SCRIPT)


@pytest.fixture(scope="module")
def ans7_old():
    if not OLD_SCRIPT.exists():
        pytest.fail(f"{OLD_SCRIPT} missing: run `git show {PIN_SHA}:<script> > {OLD_SCRIPT}` first")
    mod = _load(f"ans7_writer_{PIN_SHA}_verdict_ut", OLD_SCRIPT)
    return mod


def _render(mod, tmp: Path, metrics: dict) -> bytes:
    shutil.copy(CASE / "COMPARISON.md", tmp / "COMPARISON.md")
    shutil.copy(CASE / "metrics.json", tmp / "metrics.json")
    old = mod.CASE_DIR
    mod.CASE_DIR = tmp
    try:
        mod._write_comparison(metrics)
    finally:
        mod.CASE_DIR = old
    return (tmp / "COMPARISON.md").read_bytes()


def _section(text: str, mod) -> str:
    i = text.index(mod.VERDICT_HEAD)
    j = text.index("\n## Degree 2", i)
    return text[i:j]


def _expected_cells():
    """Every computed cell the section must carry, from T30's logged values."""
    cells = []
    _, r, _, _, rho_r = _record(85)
    _, q, _, _, rho_q = _record(89)
    _, f, d1f, d2f, rho_f = _record(93)
    cells.append(" / ".join(_sp(x, ".2f") for x in r) + " Ω")
    cells.append(" / ".join(f"{x:.2f}" for x in q))
    cells.append(" / ".join(f"{x * 1e-6:.5f}" for x in f) + " MHz")
    cells.append(f"{d1f * 1e-3:+.2f} / {d2f * 1e-3:+.2f} kHz")
    cells.append(f"ratios {rho_q:.3f} (Q) / {rho_r:.3f} (R_in)")
    cells.append(f"ratio {rho_f:.3f}")
    # first / second order on the named steps
    first = (H[1] - H[2]) / (H[0] - H[1])
    second = (H[1] ** 2 - H[2] ** 2) / (H[0] ** 2 - H[1] ** 2)
    cells.append(f"first-order quantity gives {first:.2f} and a second-order one {second:.2f}")
    pq, pr, pf = _printed_p(90), _printed_p(86), _printed_p(94)
    qi, ri, fi = _limit(q, pq), _limit(r, pr), _limit(f, pf)
    cells.append(
        f"Q∞ ≈ {qi:.2f} ({_pct(100 * (qi / q[2] - 1))} % from h = 0.005, "
        f"{_pct(100 * (qi / q[0] - 1))} % from ×1)"
    )
    cells.append(f"R_in∞ ≈ {_sp(ri, '.0f')} Ω ({_pct(100 * (ri / r[2] - 1))} %)")
    cells.append(
        f"Richardson limit ≈ {fi * 1e-6:.3f} MHz ({_pct(100 * (fi / f[2] - 1))} % from "
        f"h = 0.005)"
    )
    return cells


def test_verdict_traced(ans7):
    committed = (CASE / "COMPARISON.md").read_text()
    sec = _section(committed, ans7)
    assert sec == "\n".join(ans7.VERDICT_SECTION)
    bad = [t for t in TRACE if not _holds(*t)]
    for t in TRACE:
        print(f"(X) {t[0]!r:46} {t[1]}:{t[2]} {'ok' if t not in bad else 'FAIL'}")
    assert not bad, bad
    for cell, *_ in TRACE:
        assert cell in sec, cell
    for cell in _expected_cells():
        print(f"(X) computed cell {cell!r} in section: {cell in sec}")
        assert cell in sec, cell
    print(f"(X) {len(TRACE)} traced + {len(_expected_cells())} computed cells hold")
    # (N-cite): p = 1.6451 cited at T30:89 (one line early) fails (X).
    assert not _holds("p = 1.6451 (Q)", T30, 89, "v", "): 1.6451")
    print("(N-cite) 1.6451 at T30:89 fails (X): as asserted")


def test_richardson_recomputed(ans7):
    sec = _section((CASE / "COMPARISON.md").read_text(), ans7)
    for line, pline, rel in ((85, 86, 1e-3), (89, 90, 1e-3), (93, 94, 1e-3)):
        name, v, d1, d2, rho = _record(line)
        assert abs(d1 + d2 - (v[2] - v[0])) <= 1e-9 * abs(v[2] - v[0])
        assert abs(d2 / d1 - rho) <= 1e-12
        p = _solve_p(rho)
        pp = _printed_p(pline)
        fine = _limit(v, pp)
        coarse = _limit(v, pp, (0, 1))
        print(
            f"(L) {name}: p solved {p:.6f} printed {pp} |dp| {abs(p - pp):.2e}; "
            f"limit {fine:.6f} (fine pair) {coarse:.6f} (coarse pair) rel {abs(coarse / fine - 1):.2e}"
        )
        assert abs(p - pp) <= 5e-5
        assert abs(coarse / fine - 1) <= rel
    _, q, *_ = _record(89)
    # (N-p1) / (N-p2): the Q limit is rate-sensitive at the printed digit.
    for pc in (1.0, 2.0):
        s = f"Q∞ ≈ {_limit(q, pc):.2f}"
        print(f"(N-p{pc:.0f}) Richardson at p = {pc:.0f}: {s!r} in section: {s in sec}")
        assert s not in sec


def test_verdict_rendered(ans7, ans7_old, tmp_path):
    committed = (CASE / "COMPARISON.md").read_bytes()
    metrics_text = (CASE / "metrics.json").read_text()
    text = committed.decode()

    # (R)
    new = _render(ans7, tmp_path, json.loads(metrics_text))
    n = new.decode().count(ans7.VERDICT_HEAD)
    print(f"\n(R) rendered {len(new)} B, committed {len(committed)} B, verdict heading count {n}")
    assert n == 1
    assert new == committed
    hdr = text.index("No AED value enters this file")
    assert hdr < text.index(ans7.VERDICT_HEAD) < text.index("## Degree 2")
    assert text.endswith("\n" + ans7.CLOSING_LINE + "\n")

    # (I)
    for fn, arg in (
        (ans7._write_comparison_esr128, (json.loads(metrics_text),)),
        (ans7._write_comparison_palace, ()),
    ):
        d = tmp_path / fn.__name__
        d.mkdir()
        shutil.copy(CASE / "COMPARISON.md", d / "COMPARISON.md")
        old = ans7.CASE_DIR
        ans7.CASE_DIR = d
        try:
            fn(*arg)
        finally:
            ans7.CASE_DIR = old
        same = (d / "COMPARISON.md").read_bytes() == committed
        print(f"(I) stand-alone {fn.__name__}: byte-identical {same}")
        assert same

    # (B)
    sec = _section(text, ans7)
    for bad in ("Palace", "docs/private", "aed_results", "COMPARISON_private", "unconverged"):
        assert bad not in sec, bad
    scrub = re.sub(r"\d{8}T\d{6}Z_\S*?\.log", "", sec)
    tokens = re.findall(r"\d[\d ]*(?:\.\d+)?", scrub)
    sig = max(len(t.replace(" ", "").replace(".", "").lstrip("0")) for t in tokens)
    print(f"(B) no Palace / private path / 'unconverged'; max significant digits {sig}")
    assert sig < 8

    # (N-digit)
    assert "Q∞ ≈ 963.25" in sec
    dd = tmp_path / "ndigit"
    dd.mkdir()
    (dd / "COMPARISON.md").write_text(text.replace("Q∞ ≈ 963.25", "Q∞ ≈ 963.3"))
    assert (dd / "COMPARISON.md").read_bytes() != new
    print("(N-digit) 963.25 -> 963.3 copy differs from the render: as asserted")

    # (N-old)
    od = tmp_path / "old"
    od.mkdir()
    oldr = _render(ans7_old, od, json.loads(metrics_text))
    n_old = oldr.decode().count(ans7.VERDICT_HEAD)
    print(f"(N-old) {PIN_SHA} writer: {len(oldr)} B, verdict heading count {n_old}")
    assert n_old == 0
    assert oldr != committed

    # (R0)
    assert (CASE / "COMPARISON.md").read_bytes() == committed
    assert (CASE / "metrics.json").read_text() == metrics_text
    print("(R0) committed COMPARISON.md / metrics.json unchanged")
