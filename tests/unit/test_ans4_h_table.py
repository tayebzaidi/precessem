"""`ANS-4` step 7 — the public h-table (§9 item 417, §10 F35). Render-only, no solve.

(X)  `parse_h_table_logs` on the three tracked logs equals the committed
     ``h_table_s.json`` string for string; the ×1 rows equal
     ``ANS4_EXACT_STRIP_S_P2_X1`` exactly; cells 117 206 / 282 198 / 593 478, ×1 equal
     to ``ANS4_EXACT_STRIP_CELL_RECORD``.
(L)  the moves Δ(×1 → 0.0075) and Δ(×1 → 0.005), computed with the module's convention
     (|m_fine − m_coarse| / |m_coarse|), reproduce the logs' own printed moves
     (``…step6b.log:3726–3728``, ``…step6c.log:3840–3842``) — parsed, never retyped —
     to the printed digits.
(W)  the weekly's public figures (``docs/testing/xl-pending.md:1601``,
     PROJECT_PLAN.md:4525): Δ₂ per class, the self ratios, the p = 2 prediction, and
     entry 30's first branch (``xl-pending.md:1609``).
(R)  the tracked writer, given the committed ``metrics.json`` and ``h_table_s.json``,
     reproduces the committed ``COMPARISON.md`` byte for byte, section exactly once.
(R0) the committed files are unchanged by the test.
Negative controls (rule (e)): (N-digit-doc), (N-digit-data), (N-old).
"""

from __future__ import annotations

import json
import re
import sys
import types
from pathlib import Path

import numpy as np
import pytest

from tests.log_records import read_log_text, resolve_log_path

REPO = Path(__file__).resolve().parents[2]
CASE = REPO / "examples" / "ansys_benchmarks" / "ans4_birdcage_four_port_10_64_128MHz"
SCRIPT = CASE / "04_birdcage_four_port_10_64_128MHz.py"
PIN = REPO / "logs" / "f35_pin" / "04_pre417.py"  # host `git show c7b1b682:<script>`
LOGS = REPO / "docs" / "testing" / "logs"

# (W) comparands — the weekly's public figures, the only typed numbers
# (xl-pending.md:1601; PROJECT_PLAN.md:4525). Δ₂ in %, per frequency 10 / 64 / 128.
WEEKLY_D2_PCT = {
    "self": (0.0144, 0.0569, 0.0662),
    "adjacent": (0.0088, 0.0235, 0.0472),
    "opposite": (0.0161, 0.0611, 0.1110),
}
WEEKLY_SELF_RATIO = (0.057, 0.243, 0.301)
WEEKLY_P2_PREDICTION = 0.185
#: `xl` entry 30's decision rule, first branch (xl-pending.md:1609): max self Δ₂ ≤ 0.10 %.
#: The rule's own threshold, used to reproduce the weekly's reading — not a new band.
F27_SELF_CHANGE_CEILING = 1e-3
PCT_TOL = 5e-5   # %-points: the printed 4 dp
RATIO_TOL = 5e-4  # the printed 3 dp
FREQS = (10, 64, 128)
CLASSES = ("self", "adjacent", "opposite")


def _exec_as(name, source: Path, as_path: Path):
    mod = types.ModuleType(name)
    mod.__file__ = str(as_path)
    sys.modules[name] = mod
    exec(compile(source.read_text(), str(source), "exec"), mod.__dict__)
    return mod


@pytest.fixture(scope="module")
def mod():
    return _exec_as("ans4_h_table_live_ut", SCRIPT, SCRIPT)


def _logs(mod, scratch: Path):
    """The script's logs. One retention has gzipped is read as a scratch copy under
    its listed name, which the numbers file records (`OPS-73`)."""
    out = []
    for r, p in mod.H_TABLE_LOGS:
        found = resolve_log_path(REPO / p)
        if found.name != Path(p).name:
            copy = scratch / Path(p).name
            copy.write_text(read_log_text(found))
            found = copy
        out.append((r, found))
    return out


def _printed_moves(log: Path, tag: str):
    """{f: {class: pct}} from the log's own `S class move` print lines."""
    text = read_log_text(log)
    pat = re.compile(r"\[ANS-4 step6a\] knobbed vs item 192 x1 \(ANS4_EXACT_STRIP_S_P2_X1\) "
                     r"(\d+) MHz S class move self ([\d.]+)%\s+adjacent ([\d.]+)%\s+opposite ([\d.]+)%")
    out = {}
    for m in pat.finditer(text):
        out.setdefault(int(m.group(1)), {"self": float(m.group(2)),
                                         "adjacent": float(m.group(3)),
                                         "opposite": float(m.group(4))})
    print(f"({tag}) printed moves parsed from {log.name}: {out}")
    return out


def test_h_table(mod, tmp_path):
    from tests.validation.test_ans4_exact_strip_records import ANS4_EXACT_STRIP_S_P2_X1
    from tests.validation.test_ops61_exact_strip_ports import ANS4_EXACT_STRIP_CELL_RECORD
    from tests.validation.test_port21_feed_sensitivity import _class_means, _relative_move

    json_path = CASE / "h_table_s.json"
    committed_json = json_path.read_bytes()
    committed_md = (CASE / "COMPARISON.md").read_bytes()
    metrics_before = (CASE / "metrics.json").read_bytes()
    numbers = json.loads(committed_json)

    # (X)
    parsed = mod.parse_h_table_logs(_logs(mod, tmp_path))
    regen = (json.dumps(parsed, indent=1, sort_keys=True) + "\n").encode()
    print(f"\n(X) re-extracted numbers file identical to committed: {regen == committed_json}")
    assert regen == committed_json
    cells = [r["cells"] for r in parsed["rungs"]]
    print(f"(X) cells {cells}; record {ANS4_EXACT_STRIP_CELL_RECORD}")
    assert cells == [117206, 282198, 593478]
    assert cells[0] == ANS4_EXACT_STRIP_CELL_RECORD
    x1 = parsed["rungs"][0]
    for f in FREQS:
        got = mod._h_matrix(x1["freqs"][str(f)])
        ref = np.asarray(ANS4_EXACT_STRIP_S_P2_X1[str(f)], dtype=np.complex128)
        print(f"(X) x1 {f} MHz rows equal ANS4_EXACT_STRIP_S_P2_X1 exactly: {np.array_equal(got, ref)}")
        assert np.array_equal(got, ref)

    # (L)
    means = [{f: _class_means(mod._h_matrix(r["freqs"][str(f)])) for f in FREQS}
             for r in numbers["rungs"]]
    for fine, log in ((1, LOGS / "20261005T070009Z_ANS-4-step6b.log"),
                      (2, LOGS / "20261011T070008Z_ANS-4-step6c.log")):
        printed = _printed_moves(log, "L")
        for f in FREQS:
            mv = _relative_move(means[fine][f], means[0][f])
            for c in CLASSES:
                d = abs(100 * mv[c] - printed[f][c])
                print(f"(L) x1->rung{fine} {f} MHz {c}: {100 * mv[c]:.6f} % vs printed "
                      f"{printed[f][c]:.4f} % (|d| {d:.2e})")
                assert d <= PCT_TOL

    # (W)
    vals = mod.h_table_values(numbers)
    for i, f in enumerate(FREQS):
        for c in CLASSES:
            d = abs(100 * vals[f][c]["d2"] - WEEKLY_D2_PCT[c][i])
            print(f"(W) {f} MHz {c} Δ₂ {100 * vals[f][c]['d2']:.6f} % vs weekly "
                  f"{WEEKLY_D2_PCT[c][i]} (|d| {d:.2e}); ratio {vals[f][c]['ratio']:.5f}")
            assert d <= PCT_TOL
        dr = abs(vals[f]["self"]["ratio"] - WEEKLY_SELF_RATIO[i])
        assert dr <= RATIO_TOL
    p2 = (0.0075**2 - 0.005**2) / (0.015**2 - 0.0075**2)
    print(f"(W) p = 2 prediction {p2:.6f} vs {WEEKLY_P2_PREDICTION}")
    assert abs(p2 - WEEKLY_P2_PREDICTION) <= RATIO_TOL
    max_self = max(vals[f]["self"]["d2"] for f in FREQS)
    print(f"(W) entry 30 branch: max self Δ₂ {max_self:.6e} <= {F27_SELF_CHANGE_CEILING}")
    assert max_self <= F27_SELF_CHANGE_CEILING

    # (R)
    (tmp_path / "metrics.json").write_bytes(metrics_before)
    old = mod.CASE_DIR
    mod.CASE_DIR = tmp_path
    try:
        mod._write_comparison(json.loads(metrics_before), h_table=json_path)
    finally:
        mod.CASE_DIR = old
    render = (tmp_path / "COMPARISON.md").read_bytes()
    n_head = render.decode().count(mod.H_TABLE_HEADING)
    print(f"(R) render {len(render)} B vs committed {len(committed_md)} B: "
          f"identical={render == committed_md}; heading count {n_head}")
    assert n_head == 1
    assert render == committed_md

    # (R0)
    assert json_path.read_bytes() == committed_json
    assert (CASE / "COMPARISON.md").read_bytes() == committed_md
    assert (CASE / "metrics.json").read_bytes() == metrics_before

    # (N-digit-doc)
    planted = committed_md.decode().replace("| 0.0569 |", "| 0.0579 |", 1).encode()
    assert planted != committed_md
    print(f"(N-digit-doc) 0.0569 -> 0.0579: planted differs from render = {planted != render}")
    assert planted != render

    # (N-digit-data)
    bad_text = committed_json.decode().replace("+7.347852893e-02", "+7.357852893e-02", 1)
    assert bad_text != committed_json.decode()
    bad = json.loads(bad_text)
    print(f"(N-digit-data) planted numbers file equals re-extraction = "
          f"{bad_text.encode() == regen}")
    assert bad_text.encode() != regen
    bad_d2 = 100 * mod.h_table_values(bad)[64]["self"]["d2"]
    print(f"(N-digit-data) 64 MHz self Δ₂ {bad_d2:.6f} % (predicted ≈ 0.0528 %, printed); "
          f"|d| from 0.0569 = {abs(bad_d2 - 0.0569):.2e} > {PCT_TOL}")
    assert abs(bad_d2 - WEEKLY_D2_PCT["self"][1]) > PCT_TOL

    # (N-old) the pinned pre-change writer renders no h-table section
    if not PIN.exists():
        pytest.fail(f"{PIN} missing: pin it by host `git show c7b1b682:<script>` first")
    old_mod = _exec_as("ans4_h_table_pin_ut", PIN, SCRIPT)
    (tmp_path / "old").mkdir()
    (tmp_path / "old" / "metrics.json").write_bytes(metrics_before)
    old_mod.CASE_DIR = tmp_path / "old"
    old_mod._write_comparison(json.loads(metrics_before))
    old_render = (tmp_path / "old" / "COMPARISON.md").read_bytes()
    n_old = old_render.decode().count(mod.H_TABLE_HEADING)
    print(f"(N-old) pinned render heading count {n_old}; differs from committed = "
          f"{old_render != committed_md}")
    assert n_old == 0
    assert old_render != committed_md
