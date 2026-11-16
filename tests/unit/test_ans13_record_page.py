"""`ANS-13` — the `ans:13` record page (§9 item 958, §10 C53;
docs/planning/chunks/ANS-13.md "§9 item 958").

Our degree-1 columns of the F-clinical box-heating SPEC, read from ``xl`` entries
59 / 60's tracked logs through C46's reader; no solve, no mesh.

(X) C46's ``check_block`` on both blocks; custody; SPEC row 2's coordinates; every
    record-table token matched once with its cite; the thermal table.
(D) (D-dep) 24 checks; (D-c*); (D-ratio); (D-lever); (D-cal).
(R) the record table rendered from the logs == ``SPEC.md:140–148`` bytes.
(R-cmp) absent until item 959 (C54) lands.
Controls (asserted): (N-digit-log), (N-swap), (N-ratio), (N-doc).
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
CASE = REPO / "examples" / "ansys_benchmarks" / "ans13_f_clinical_box_heating_64_128MHz"
SCRIPT = CASE / "13_f_clinical_box_heating_64_128MHz.py"


def _load():
    name = "ans13_record_page"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def page():
    return _load()


@pytest.fixture(scope="module")
def ctx(page):
    return page.read_all()


def _fails(results):
    return [(n, d) for n, ok, d in results if not ok]


def _show(tag, res):
    for n, ok, d in res:
        print(f"[ANS-13 {tag}] {'PASS' if ok else 'FAIL'} {n}: {d}")


def test_x_extraction(page, ctx):
    res = page.check_x(ctx)
    _show("X", res)
    assert not _fails(res), _fails(res)
    b, x = ctx["blocks"], ctx["extras"]
    # The brief's lines (docs/planning/chunks/ANS-13.md, §9 item 958, "Reads").
    assert (b["59"]["head"]["line"], b["59"]["end"]) == (139, 12356)
    assert (b["60"]["head"]["line"], b["60"]["end"]) == (141, 12307)
    assert [x["59"]["senses"]["ccw"]["vsub"]["line"], x["60"]["senses"]["ccw"]["vsub"]["line"]] \
        == [12315, 12265]
    assert {k: (x[k]["elapsed"]["line"], x[k]["peak"]["line"], x[k]["cmd"]["line"])
            for k in ("59", "60")} == {"59": (12457, 12452, 12), "60": (12375, 12370, 12)}
    assert [x["59"]["norm"]["line"], x["60"]["norm"]["line"]] == [12312, 12262]
    # Entry 60's cw probes sit after the cut at :12298–12303 (item 955's note).
    assert [p["line"] for p in b["60"]["senses"]["cw"]["p61"]["probes"].values()] \
        == list(range(12298, 12304))
    # C46's checks are all in, both blocks, both senses.
    names = [n for n, _, _ in res]
    for want in ("(Node)", "(D-f)", "(P_gel)", "(D-scale)", "(Cem)", "(RA)", "(T*)", "(Probes)"):
        assert sum(want in n for n in names) >= 2, want
    cov = [d for n, _, d in res if "every SPEC record-table token" in n]
    print(f"[ANS-13 X] coverage: {cov}")
    assert cov and cov[0].startswith("52 of 52"), cov


def test_x_raw_multiplicity(page, ctx):
    """Reading 1: matched once to its row's field; the raw multiplicity is printed."""
    sp = ctx["spec"]
    for k in ("59", "60"):
        text = page._read_log(ctx["paths"][k])
        for name, tok in (("factor", sp["factor"][f"f{k}"]), ("c*", sp["cstar"][f"t{k}"]),
                          ("ball", sp["ball"][f"b{k}"])):
            print(f"[ANS-13 X-mult] entry {k} {name} {tok}: {text.count(tok)} occurrences")
            assert text.count(tok) >= 1
    print(f"[ANS-13 X-mult] C39 table copies: {[t['head'] for t in ctx['c39']['tables']]}")
    assert [t["head"] for t in ctx["c39"]["tables"]] == [371, 407]


def test_d_derived(page, ctx):
    res = page.check_d(ctx)
    _show("D", res)
    names = [n for n, _, _ in res]
    assert sum(n.startswith("(D-dep)") for n in names) == 24
    for want, count in (("(D-c*)", 2), ("(D-ratio)", 4), ("(D-lever)", 2), ("(D-cal)", 2)):
        assert sum(n.startswith(want) for n in names) == count, want
    assert not _fails(res), _fails(res)


def test_r_regeneration(page, ctx):
    res = page.check_r(ctx)
    _show("R", res)
    assert not _fails(res), _fails(res)
    assert ctx["spec"]["table"]["line"] == 140


def test_r_cmp_conditional(page):
    res, msg = page.check_r_cmp()
    print(f"[ANS-13 R-cmp] {msg}")
    if page.COMPARISON_MD.exists():
        assert msg == "live" and not _fails(res), (msg, res)
    else:
        assert msg == "absent — item 959 (C54)" and res == []


def test_n_digit_log(page, tmp_path):
    """(N-digit-log): the 4th significant digit of entry 60's ccw P3 per-P_acc token
    (`:12280`) changed in a tmp copy: (D-dep) red on that probe alone; (X) and (R) red."""
    src = page.entry_log("60")
    lines = page._read_log(src).splitlines(keepends=True)
    i = 12280 - 1
    old, new = "per P_acc watt +2.939765335e-03 K", "per P_acc watt +2.930765335e-03 K"
    assert "probe P3 +x wall, sense (1,-1j)" in lines[i] and old in lines[i]
    lines[i] = lines[i].replace(old, new)
    bad = tmp_path / src.name
    bad.write_text("".join(lines))
    c = page.read_all(log60=bad)
    d_f = [n for n, _ in _fails(page.check_d(c))]
    x_f = [n for n, _ in _fails(page.check_x(c))]
    r_f = [n for n, _ in _fails(page.check_r(c))]
    print(f"[ANS-13 N-digit-log] (D) {d_f}; (X) {x_f}; (R) {r_f}")
    assert [n for n in d_f if n.startswith("(D-dep)")] == [
        "(D-dep) 3 T ccw P3 +x wall: θ per P_acc / factor == I32 per-deposited-watt"]
    assert "(X) probes.p60[3]" in x_f
    assert any("(D-scale) probe P3 +x wall" in n for n in x_f)
    assert r_f == ["(R) record table rendered from the logs == SPEC bytes"]


def test_n_swap(page, ctx):
    """(N-swap): entry 59's block read as 3 T: (Node) red; (D-c*) red at ≈ 9.4 % against
    a ≈ 1.6e-4 bound."""
    c = page.read_all(swap=True)
    x_f = [n for n, _ in _fails(page.check_x(c))]
    d = page.check_d(c)
    dc = [(n, det) for n, ok, det in d if n.startswith("(D-c*)") and not ok]
    print(f"[ANS-13 N-swap] (X) node fails {[n for n in x_f if '(Node)' in n]}; (D-c*) {dc}")
    assert any("3 T (entry 60) (Node)" in n for n in x_f)
    assert len(dc) == 2
    rels = [float(det.split("rel ")[1].split(" ")[0]) for _, det in dc]
    assert max(rels) > 0.09 and min(rels) > 0.08, rels


def test_n_ratio(page, ctx, tmp_path):
    """(N-ratio): one digit of one C39 ratio changed in a tmp copy: (D-ratio) red."""
    lines = page._read_log(page.C39_LOG).splitlines(keepends=True)
    i = 373 - 1
    assert lines[i].rstrip("\n").endswith("| 4.042515e+00")
    lines[i] = lines[i].replace("| 4.042515e+00", "| 4.042516e+00")
    bad = tmp_path / page.C39_LOG.name
    bad.write_text("".join(lines))
    c = dict(ctx, c39=page.parse_c39(bad))
    d_f = [n for n, _ in _fails(page.check_d(c))]
    print(f"[ANS-13 N-ratio] (D) fails {d_f}")
    assert "(D-ratio) C46 ratios(side(59), side(60)) .6e == C39's strings" in d_f


def test_n_doc(page, ctx, tmp_path):
    """(N-doc): one digit in a tmp SPEC copy's `:144` makes (R) red."""
    text = page.SPEC_PATH.read_text()
    lines = text.split("\n")
    old = "7.476681123e-04"
    assert old in lines[144 - 1] and text.count(old) == 1
    spec_bad = tmp_path / "SPEC.md"
    spec_bad.write_text(text.replace(old, "7.476691123e-04"))
    c = dict(ctx, spec=page.parse_spec13(spec_bad))
    r_f = [n for n, _ in _fails(page.check_r(c))]
    print(f"[ANS-13 N-doc] (R) fails {r_f}")
    assert r_f == ["(R) record table rendered from the logs == SPEC bytes"]
