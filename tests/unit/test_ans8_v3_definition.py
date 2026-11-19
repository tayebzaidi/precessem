"""`ANS-8` step 2 to the operator (§9 item 411, §10 H22): SPEC v3's accepted-power
definition in ``COMPARISON.md`` and the arithmetic of SPEC.md's operator brief.

Render-only, on scratch copies of the committed ``metrics.json`` / ``COMPARISON.md``:

(R)     ``comparison_text`` on the committed metrics equals the committed ``COMPARISON.md``
        byte for byte;
(P)     everything before the ``## Step 2 (SPEC v3)`` heading equals the pinned pre-change
        writer's render (``8f50751e``, copied under ``logs/ans8_writer_8f50751e/`` at the
        script's own relative path by a host ``git show`` before the window);
(D)     the v3 section carries the definition paragraph once (accepted power at the 50 Ω
        lumped port, Cauchy–Schwarz, the rendered ``rel_line_cut`` "0.83 %"), no longer the
        old "feed lines" clause, and every v3 table row ends ``| |`` (AED cells empty);
(Rp)    R_p = 1/(ω²C²r) and Q_C = 1/(ωCr) from the imported C_in and ESR reproduce the
        SPEC's 2 212 Ω / 105, and the operator subsection prints 118.2437 and 2 212; the
        exact series→parallel C factor 1/(1 + Q_C⁻²) − 1 is printed only;
(R0)    the repo's ``COMPARISON.md``, ``metrics.json`` and ``SPEC.md`` are untouched.
Negative controls (asserted): (N-digit) one planted digit in a tmp copy differs from the
render; (N-old) the pinned writer's render differs and lacks "Cauchy–Schwarz" in its v3
section; (N-Rp) r = 0.02 Ω does not give 2 212 Ω.
"""

from __future__ import annotations

import importlib.util
import json
import math
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from tests.validation.test_wf10_f_human_copper_esr_root import ESR_OHM
from tests.validation.test_wf10_f_human_tuned_b1_sar import C_IN_F_HUMAN_D2_PF, C_IN_F_HUMAN_PF

REPO = Path(__file__).resolve().parents[2]
CASE = REPO / "examples" / "ansys_benchmarks" / "ans8_f_human_quadrature_mode_content_64MHz"
SCRIPT = CASE / "08_f_human_quadrature_mode_content.py"
# The pinned copy mirrors the script's depth (item 400's precedent).
OLD_SCRIPT = REPO / "logs" / "ans8_writer_8f50751e" / SCRIPT.relative_to(REPO)
V3_HEAD = "## Step 2 (SPEC v3)"
BRIEF_HEAD = "### Operator session brief (step 2 — one AED session)"
F_HZ = 64.0e6


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def ans8():
    return _load("ans8_v3def_ut", SCRIPT)


@pytest.fixture(scope="module")
def ans8_old():
    if not OLD_SCRIPT.exists():
        pytest.fail(f"{OLD_SCRIPT} missing: run `git show 8f50751e:<script> > {OLD_SCRIPT}` first")
    return _load("ans8_writer_8f50751e_ut", OLD_SCRIPT)


def _v3(text: str) -> str:
    i = text.index(V3_HEAD)
    return text[i:]


def test_comparison_definition(ans8, ans8_old, tmp_path):
    committed = (CASE / "COMPARISON.md").read_bytes()
    metrics_text = (CASE / "metrics.json").read_text()
    spec_bytes = (CASE / "SPEC.md").read_bytes()
    shutil.copy(CASE / "metrics.json", tmp_path / "metrics.json")
    shutil.copy(CASE / "COMPARISON.md", tmp_path / "COMPARISON.md")
    metrics = json.loads((tmp_path / "metrics.json").read_text())

    # (R)
    new = ans8.comparison_text(metrics)
    print(
        f"\n(R) rendered {len(new.encode())} B, committed {len(committed)} B, "
        f"equal {new.encode() == committed}"
    )
    assert new.encode() == committed

    # (P) and (N-old)
    old = ans8_old.comparison_text(json.loads(metrics_text))
    pre_new, pre_old = new[: new.index(V3_HEAD)], old[: old.index(V3_HEAD)]
    print(
        f"(P) prefix before the v3 heading: new {len(pre_new)} chars, pinned {len(pre_old)} "
        f"chars, equal {pre_new == pre_old}"
    )
    assert pre_new == pre_old
    print(
        f"(N-old) pinned 8f50751e render differs: {old.encode() != committed}; its v3 section "
        f"has Cauchy–Schwarz: {'Cauchy–Schwarz' in _v3(old)} (ASSERTED False)"
    )
    assert old.encode() != committed
    assert "Cauchy–Schwarz" not in _v3(old)

    # (D)
    v3 = _v3(new)
    rel = metrics["rows_v3"]["gap"]["rel_line_cut"]
    pct = f"{100.0 * rel:.2f} %"
    n_def = v3.count("accepted power at the 50 Ω lumped port")
    print(
        f"(D) definition count {n_def}; rel_line_cut {rel!r} rendered '{pct}' count "
        f"{v3.count(pct)}; Cauchy–Schwarz count {v3.count('Cauchy–Schwarz')}; old clause "
        f"present {'incident − reflected at the two 50 Ω feed lines' in new}"
    )
    assert pct == "0.83 %"
    assert n_def == 1
    assert v3.count(pct) == 1
    assert "Cauchy–Schwarz" in v3
    assert "incident − reflected at the two 50 Ω feed lines" not in new
    rows = [ln for ln in v3.splitlines() if ln.startswith("| (")]
    print(
        f"(D) {len(rows)} v3 table rows, all AED cells empty: "
        f"{all(ln.endswith('| |') for ln in rows)}"
    )
    assert rows and all(ln.endswith("| |") for ln in rows)

    # (N-digit)
    planted = tmp_path / "COMPARISON.md"
    text = planted.read_text()
    i = text.index(V3_HEAD)
    head, tail = text[:i], text[i:]
    assert tail.count("0.268112") >= 1
    planted.write_text(head + tail.replace("0.268112", "0.268113", 1))
    p_bytes = planted.read_bytes()
    first = next(a for a, b in zip(p_bytes.decode().splitlines(), new.splitlines()) if a != b)
    print(
        f"(N-digit) planted copy differs from the render: {p_bytes != new.encode()}; "
        f"first differing line: {first}"
    )
    assert p_bytes != new.encode()

    # (R0)
    assert (CASE / "COMPARISON.md").read_bytes() == committed
    assert (CASE / "metrics.json").read_text() == metrics_text
    assert (CASE / "SPEC.md").read_bytes() == spec_bytes
    print("(R0) repo COMPARISON.md / metrics.json / SPEC.md unchanged")


def _rp(c_f: float, r: float) -> tuple[float, float]:
    w = 2.0 * math.pi * F_HZ
    return 1.0 / (w * w * c_f * c_f * r), 1.0 / (w * c_f * r)


def test_operator_brief_arithmetic():
    c_f = C_IN_F_HUMAN_PF * 1e-12
    r_p, q_c = _rp(c_f, ESR_OHM)
    factor = 1.0 / (1.0 + q_c**-2) - 1.0
    print(
        f"\n(Rp) C {C_IN_F_HUMAN_PF} pF, r {ESR_OHM} Ω: R_p {r_p:.6f} Ω (round "
        f"{round(r_p)}), Q_C {q_c:.6f} (round {round(q_c)}); exact parallel-C factor "
        f"1/(1 + Q_C⁻²) − 1 = {factor:.3e} (printed only; ΔC {factor * C_IN_F_HUMAN_PF:.3e} pF "
        f"vs the 0.10 pF root grid)"
    )
    assert round(r_p) == 2212
    assert round(q_c) == 105
    spec = (CASE / "SPEC.md").read_text()
    sub = spec[spec.index(BRIEF_HEAD) : spec.index("## Revision history")]
    print(
        f"(Rp) brief subsection {len(sub)} chars; '118.2437' {'118.2437' in sub}; "
        f"'2 212' {'2 212' in sub}"
    )
    assert spec.count(BRIEF_HEAD) == 1
    assert "118.2437" in sub and "2 212" in sub

    # (N-Rp)
    r_bad, _ = _rp(c_f, 0.02)
    print(f"(N-Rp) r = 0.02 Ω gives R_p {r_bad:.3f} Ω (round {round(r_bad)} ≠ 2212 ASSERTED)")
    assert round(r_bad) != 2212


# --- SPEC v4 (§9 item 513, §10 H29): step 2's capacitor at the degree-2 root -----------
PRE_V4_SHA = "e5d287d4"  # last commit of SPEC.md before v4 (git log -1, read 2026-10-18)
SPEC_REL = "examples/ansys_benchmarks/ans8_f_human_quadrature_mode_content_64MHz/SPEC.md"
RECORD_HEAD = "### Degree-1 record (v3 / v3.1 values — do not enter)"


def _brief(text: str) -> str:
    """The operator brief alone: to the degree-1 record if present, else to the history."""
    end = RECORD_HEAD if RECORD_HEAD in text else "## Revision history"
    return text[text.index(BRIEF_HEAD) : text.index(end)]


def _brief_v4_offences(sub: str) -> list[str]:
    bad = [s for s in ("118.2437", "2 212") if s in sub]
    bad += [f"missing {s!r}" for s in ("115.5209", "`115.5209pF`", "2 317") if s not in sub]
    return bad


def test_operator_brief_arithmetic_v4():
    """SPEC v4: C = C_in(d2), R_p ‖ C and Q_C recomputed, the degree-1 values a record.

    Disclosed (review ruling 2026-10-18 03:00, Q1): the degree-1 record subsection sits
    between the brief and ``## Revision history``, i.e. inside the span item 411's
    ``test_operator_brief_arithmetic`` reads, so that test stays unedited and green. The
    brief proper (``_brief``) is asserted free of the degree-1 values here (B4).
    """
    spec_bytes = (REPO / SPEC_REL).read_bytes()
    spec = spec_bytes.decode()

    # (Rp4) closed forms on the imported constants
    w = 2.0 * math.pi * F_HZ
    c_d2, c_d1 = C_IN_F_HUMAN_D2_PF * 1e-12, C_IN_F_HUMAN_PF * 1e-12
    r_p, q_c = _rp(c_d2, ESR_OHM)
    r_p_exact = ESR_OHM * (1.0 + q_c * q_c)
    cp_factor = q_c * q_c / (1.0 + q_c * q_c) - 1.0
    r_p_d1, _ = _rp(c_d1, ESR_OHM)
    ratio, ratio_cf = r_p / r_p_d1, (C_IN_F_HUMAN_PF / C_IN_F_HUMAN_D2_PF) ** 2
    print(
        f"\n(Rp4) C {C_IN_F_HUMAN_D2_PF} pF, r {ESR_OHM} Ω, ω {w:.6e}: R_p high-Q "
        f"{r_p:.6f} Ω (round {round(r_p)}), exact r(1+Q²) {r_p_exact:.6f} Ω (round "
        f"{round(r_p_exact)}), difference {r_p_exact - r_p:.12f} Ω vs r {ESR_OHM}; Q_C "
        f"{q_c:.6f} (round {round(q_c)}); C_p/C − 1 = {cp_factor:.4e} (printed only, "
        f"predicted −8.631e-5; ΔC {cp_factor * C_IN_F_HUMAN_D2_PF:.4f} pF vs the 0.50 pF "
        f"bracket); R_p(d2)/R_p(d1) {ratio:.9f} vs (C_d1/C_d2)² {ratio_cf:.9f}"
    )
    assert C_IN_F_HUMAN_D2_PF == 115.5209
    assert round(r_p) == 2317
    assert round(q_c) == 108
    assert round(r_p_exact) == 2317
    assert abs((r_p_exact - r_p) - ESR_OHM) <= 1e-9 * ESR_OHM
    assert abs(ratio / ratio_cf - 1.0) <= 1e-12

    # (N-Rp4)
    r_bad, _ = _rp(c_d2, 0.02)
    r_blind, _ = _rp(c_d1, ESR_OHM)
    print(
        f"(N-Rp4) r = 0.02 Ω: R_p {r_bad:.3f} Ω (round {round(r_bad)} ≠ 2317 ASSERTED); "
        f"v3-blind C {C_IN_F_HUMAN_PF} pF: R_p {r_blind:.3f} Ω (round {round(r_blind)} == "
        f"2212 ASSERTED); separation 2317/2212 = {2317 / 2212:.4f}"
    )
    assert round(r_bad) != 2317
    assert round(r_blind) == 2212 and round(r_blind) != 2317

    # (B4) the brief proper
    sub = _brief(spec)
    bad = _brief_v4_offences(sub)
    print(f"(B4) brief {len(sub)} chars; offences {bad}")
    assert bad == []

    # (F4) fixture paragraph and row 1
    fix = spec[spec.index("**Fixture (both codes):**") : spec.index("**Quantities to export")]
    row1 = [ln for ln in spec.splitlines() if ln.startswith("| 1 | Im Y_in")]
    print(
        f"(F4) fixture {len(fix)} chars: 'C = 115.5209 pF' {'C = 115.5209 pF' in fix}, "
        f"'2 317 Ω ‖ 115.52 pF' {'2 317 Ω ‖ 115.52 pF' in fix}, '= 108' {'= 108' in fix}, "
        f"'C = 118.2437 pF' {'C = 118.2437 pF' in fix}; row 1 lines {len(row1)}"
    )
    assert "C = 115.5209 pF" in fix
    assert "2 317 Ω ‖ 115.52 pF" in fix
    assert "= 108" in fix
    assert "C = 118.2437 pF" not in fix
    assert len(row1) == 1 and "C = 115.5209 pF" in row1[0]

    # (Rec) the labelled degree-1 record and the history line
    rec = spec[spec.index(RECORD_HEAD) : spec.index("## Revision history")]
    v4_lines = [ln for ln in spec.splitlines() if ln.startswith("- **v4 ")]
    print(
        f"(Rec) record heads {spec.count(RECORD_HEAD)}; record has '118.2437' "
        f"{'118.2437' in rec}, '2 212' {'2 212' in rec}; v4 history lines {len(v4_lines)}"
    )
    assert spec.count(RECORD_HEAD) == 1
    assert "118.2437" in rec and "2 212" in rec
    assert len(v4_lines) == 1

    # (N-pin) the pre-v4 SPEC fails (B4)
    old = subprocess.run(
        [
            "git",
            "-c",
            f"safe.directory={REPO}",
            "-C",
            str(REPO),
            "show",
            f"{PRE_V4_SHA}:{SPEC_REL}",
        ],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    old_sub = _brief(old)
    old_bad = _brief_v4_offences(old_sub)
    print(
        f"(N-pin) {PRE_V4_SHA} SPEC: record head present {RECORD_HEAD in old}; brief "
        f"offences {old_bad} (first {old_bad[0] if old_bad else None!r})"
    )
    assert RECORD_HEAD not in old
    assert "118.2437" in old_sub and "115.5209" not in old_sub
    assert old_bad

    # (R0') the test leaves SPEC.md untouched
    assert (REPO / SPEC_REL).read_bytes() == spec_bytes
    print("(R0') repo SPEC.md unchanged")
