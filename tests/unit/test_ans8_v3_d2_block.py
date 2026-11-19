"""`ANS-8` H31 Part B (§9 item 597, §10 H31): H30's degree-2 rows at the root as a
new ``COMPARISON.md`` block after the v3 section, insertions only, rendered by the
tracked writer from the committed ``metrics.json``.

Render-only, on the committed files (nothing written to the repo):

(R2)       ``comparison_text(committed metrics.json)`` equals the committed
           ``COMPARISON.md`` byte for byte;
(Ins)      the committed ``COMPARISON.md`` starts with ``e5d287d4``'s bytes and the
           remainder begins ``"\\n## Step 2 at the degree-2 root"``; the committed
           metrics with ``rows_v3_d2`` popped deep-equals ``8f50751e``'s;
(P_acc-1)  ``check_v3_power(rows_v3_d2, DISCRETE_IDENTITY_RTOL)``;
(Deg)      degree 2, ``c_in_pf == C_IN_F_HUMAN_D2_PF``, ``node_knob == "d2root"``;
(Dig)      every number in the block's "ours" cells (degree-2 values) occurs verbatim
           in a ``[H30a H30]`` line of the tracked H30 log; every degree-1 value beside
           it is ``rows_v3``'s own float in the same format;
(Row2)     §9 item 740: the row 2 row's four values and ``c_mode_pf`` equal
           ``row2_circuit_route_h25``'s strings; that row is H25's, not in (Dig).
Negative controls (asserted): (N-absent) the render without ``rows_v3_d2`` equals
``e5d287d4``'s bytes; (N-digit) a planted digit breaks (R2); (N-d1) the block rendered
with ``rows_v3``'s degree-1 row 4 in place of H30's fails (Dig); (N-row) the key popped
renders neither the row nor the closing sentence; (N-digit′) a planted digit in the key
renders a page unequal to the committed one.
"""

from __future__ import annotations

import copy
import importlib.util
import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

from tests.log_records import read_log_text
from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL
from tests.validation.test_wf10_f_human_tuned_b1_sar import C_IN_F_HUMAN_D2_PF

REPO = Path(__file__).resolve().parents[2]
CASE = REPO / "examples" / "ansys_benchmarks" / "ans8_f_human_quadrature_mode_content_64MHz"
SCRIPT = CASE / "08_f_human_quadrature_mode_content.py"
CASE_REL = CASE.relative_to(REPO).as_posix()
COMPARISON_PRE_SHA = "e5d287d4"  # last commit of COMPARISON.md before this item
METRICS_PRE_SHA = "8f50751e"  # last commit of metrics.json before this item
H30_LOG = REPO / "docs" / "testing" / "logs" / "20261022T093054Z_ANS-8-H30-W1.log"
D2_HEAD = "## Step 2 at the degree-2 root"
NUM = re.compile(r"[-+]?\d+\.\d+(?:e[-+]\d+)?")


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def _git_show(sha: str, rel: str) -> bytes:
    return subprocess.run(
        ["git", "-c", f"safe.directory={REPO}", "-C", str(REPO), "show", f"{sha}:{rel}"],
        check=True,
        capture_output=True,
    ).stdout


@pytest.fixture(scope="module")
def ans8():
    return _load("ans8_v3d2_ut", SCRIPT)


def _block_cells(text: str):
    """(degree-2 cell, degree-1 parenthetical or None) for each block table row."""
    block = text[text.index(D2_HEAD) :]
    out = []
    for ln in block.splitlines():
        if not ln.startswith("| ") or ln.startswith("| Row") or ln.startswith("|---"):
            continue
        if ln.startswith("| row 2"):  # H25's digits, not [H30a H30] numbers: gated by (Row2)
            continue
        ours = ln.split(" | ")[1]
        m = re.match(r"^(.*?)(?: \(degree 1 record: (.*)\))?$", ours)
        out.append((m.group(1), m.group(2)))
    return out


def _dig(text: str, h30_lines: list[str]):
    """Degree-2 numbers absent from every ``[H30a H30]`` line."""
    hay = "\n".join(h30_lines)
    missing = []
    for d2_cell, _ in _block_cells(text):
        for tok in NUM.findall(d2_cell):
            if tok not in hay:
                missing.append(tok)
    return missing


def test_ans8_v3_d2_block(ans8):
    committed = (CASE / "COMPARISON.md").read_bytes()
    metrics = json.loads((CASE / "metrics.json").read_text())
    d2 = metrics["rows_v3_d2"]
    new = ans8.comparison_text(metrics)

    # (R2)
    print(
        f"(R2) render == committed COMPARISON.md: {new.encode() == committed} "
        f"({len(committed)} bytes)"
    )
    assert new.encode() == committed

    # (Row2) §9 item 740: the rendered row's four values and c_mode_pf are the key's strings
    key = d2["row2_circuit_route_h25"]
    row = next(ln for ln in new.splitlines() if ln.startswith("| row 2"))
    want = [
        key["c_mode_pf"],
        key["lossless"]["ccw_feed_share"],
        key["lossless"]["ccw_feed_leak_k0"],
        key["esr_0p2_ohm"]["ccw_feed_share"],
        key["esr_0p2_ohm"]["ccw_feed_leak_k0"],
    ]
    got = re.findall(r"\b\d+\.\d+\b", row.replace("r = 0.2 Ω", ""))
    print(f"(Row2) rendered {got}, key {want}")
    assert got == want

    # (Ins)
    pre_cmp = _git_show(COMPARISON_PRE_SHA, f"{CASE_REL}/COMPARISON.md")
    rest = committed[len(pre_cmp) :].decode()
    pre_metrics = json.loads(_git_show(METRICS_PRE_SHA, f"{CASE_REL}/metrics.json"))
    popped = copy.deepcopy(metrics)
    popped.pop("rows_v3_d2")
    print(
        f"(Ins) starts with {COMPARISON_PRE_SHA} ({len(pre_cmp)} bytes): "
        f"{committed.startswith(pre_cmp)}; remainder begins {rest[:32]!r}; metrics less "
        f"rows_v3_d2 == {METRICS_PRE_SHA}: {popped == pre_metrics}"
    )
    assert committed.startswith(pre_cmp)
    assert rest.startswith("\n" + D2_HEAD)
    assert popped == pre_metrics

    # (P_acc-1)
    r1, r2, ok = ans8.check_v3_power(d2, DISCRETE_IDENTITY_RTOL)
    print(f"(P_acc-1) r1 {r1:.3e} (band {DISCRETE_IDENTITY_RTOL:g} ASSERTED); r2 {r2:.4e} printed")
    assert ok

    # (Deg)
    print(
        f"(Deg) degree {d2['degree']}, c_in_pf {d2['c_in_pf']} (expected {C_IN_F_HUMAN_D2_PF}), "
        f"node_knob {d2['node_knob']!r}"
    )
    assert d2["degree"] == 2
    assert d2["c_in_pf"] == C_IN_F_HUMAN_D2_PF
    assert d2["node_knob"] == "d2root"

    # (Dig)
    h30 = [ln for ln in read_log_text(H30_LOG).splitlines() if ln.startswith("[H30a H30]")]
    cells = _block_cells(new)
    missing = _dig(new, h30)
    n_tok = sum(len(NUM.findall(c)) for c, _ in cells)
    print(
        f"(Dig) {len(cells)} block rows, {n_tok} degree-2 numbers, {len(h30)} [H30a H30] "
        f"lines; missing {missing}"
    )
    assert len(h30) >= 5 and len(cells) == 9 and n_tok == 9
    assert not missing
    v3 = metrics["rows_v3"]
    p1 = v3["power_w"]["p_acc_spec"]
    m1 = v3["gap"]["mechanism"]
    want_d1 = [
        f"{p1:.12e}",
        f"{v3['row4']['sar10g_per_w']:.7f}",
        *(f"{v3['row5'][k]:.6f}" for k in ("phantom", "copper", "esr", "sum")),
        *(f"{m1[n]['x_w']:.9e}" for n in v3["feeds"]),
        f"{sum(m1[n]['x_w'] for n in v3['feeds']) / p1:.6e}",
    ]
    got_d1 = [c for _, c in cells]
    print(f"(Dig) degree-1 cells {got_d1}")
    assert got_d1 == want_d1

    # (D)-compatibility: the block carries none of the v3 section's single-count strings
    block = new[new.index(D2_HEAD) :]
    assert "0.83 %" not in block
    assert "accepted power at the 50 Ω lumped port" not in block
    assert all(ln.endswith("| |") for ln in block.splitlines() if ln.startswith("| ("))
    assert not block.startswith("## Step 2 (SPEC v3)")


def test_ans8_v3_d2_negative_controls(ans8):
    metrics = json.loads((CASE / "metrics.json").read_text())
    committed = (CASE / "COMPARISON.md").read_bytes()

    # (N-absent)
    absent = copy.deepcopy(metrics)
    absent.pop("rows_v3_d2")
    pre_cmp = _git_show(COMPARISON_PRE_SHA, f"{CASE_REL}/COMPARISON.md")
    same = ans8.comparison_text(absent).encode() == pre_cmp
    print(f"(N-absent) render without rows_v3_d2 == {COMPARISON_PRE_SHA}: {same} (ASSERTED True)")
    assert same

    # (N-digit)
    text = committed.decode()
    i = text.index("2.910591023475e-02", text.index(D2_HEAD))
    planted = text[:i] + "3" + text[i + 1 :]
    print(
        f"(N-digit) planted digit equals render: {planted.encode() == committed} (ASSERTED False)"
    )
    assert planted.encode() != committed

    # (N-d1)
    h30 = [ln for ln in read_log_text(H30_LOG).splitlines() if ln.startswith("[H30a H30]")]
    bad = copy.deepcopy(metrics)
    bad["rows_v3_d2"]["row4"]["sar10g_per_w"] = metrics["rows_v3"]["row4"]["sar10g_per_w"]
    missing = _dig(ans8.comparison_text(bad), h30)
    print(f"(N-d1) degree-1 row 4 in the degree-2 cell: missing {missing} (ASSERTED ['0.1587143'])")
    assert missing == ["0.1587143"]

    # (N-row) §9 item 740: the key popped renders neither the row nor the closing sentence
    no_row = copy.deepcopy(metrics)
    no_row["rows_v3_d2"].pop("row2_circuit_route_h25")
    text_nr = ans8.comparison_text(no_row)
    print(
        f"(N-row) row {'| row 2' in text_nr}, sentence {'Row 2 is now carried' in text_nr}, "
        f"equals committed {text_nr.encode() == committed} (ASSERTED False ×3)"
    )
    assert "| row 2" not in text_nr and "Row 2 is now carried" not in text_nr
    assert text_nr.encode() != committed

    # (N-digit′) a planted digit in the key renders a page unequal to the committed one
    planted_key = copy.deepcopy(metrics)
    planted_key["rows_v3_d2"]["row2_circuit_route_h25"]["lossless"]["ccw_feed_share"] = "0.968089"
    same_key = ans8.comparison_text(planted_key).encode() == committed
    print(f"(N-digit′) planted key digit equals committed: {same_key} (ASSERTED False)")
    assert not same_key
