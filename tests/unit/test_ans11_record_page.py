"""`ANS-11` — the `ans:11` record page (§9 item 833; docs/planning/chunks/ANS-11.md).

Our degree-1 column of the F-clinical-1.5T box-loaded SPEC, read from entry 52's
tracked C6 log; no solve, no mesh.

(X) extraction: lattice labels == C6's `_lattice()`; SPEC-quoted values (parsed from
    SPEC.md, never retyped) == their log tokens; the SPEC's 36 six-digit lattice
    entries == ``f"{float(tok):.6e}"``; the parsed lines are the brief's.
(D) derived from the parsed tokens, at half a unit in the last printed digit.
(R) the writer on the committed ``metrics.json`` reproduces the committed
    ``COMPARISON.md`` byte for byte; the parse reproduces ``metrics.json``.
Controls (asserted): (N-digit-log), (N-drop), (N-doc).  (N-sense) is printed only
(rule (e): predicted from hand arithmetic, no prior measurement).

C36 (§9 item 871): the entry-62 reader ``read_c19`` (absent on this tree; stubs on
``tmp_path`` only, never under ``/workspace/logs/geo38`` or ``/geo36``), the degree-2
mapping, (K) SPEC-vs-COMPARISON strings, (D-node) / (D-v1) / (D-v2) node figures, and
controls (N-K), (N-D), (N-stub-key), (N-stub-node), (N-stub-deg).  The node-C figures
are printed only (rule (e)).

C37 (§9 item 983): the degree-2 column parsed from entry 62's tracked log
(``degree2_block``), (D2) log-internal identities, (X2) degree-2 cells carry their
tokens, (J62) custody against the host JSON when present, the moves and the order
rule (printed verdict); controls (N-J62), (N-X2), (N-grid), (N-node).
``test_c19_absent`` is retired (premise expired) in favour of ``test_c37_custody``.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
CASE = REPO / "examples" / "ansys_benchmarks" / "ans11_f_clinical_1p5t_box_loaded_64MHz"
SCRIPT = CASE / "11_f_clinical_1p5t_box_loaded_64MHz.py"


def _load():
    name = "ans11_record_page"
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
def parsed(page):
    return page.parse_c6_log(page.SOURCE_LOG)


@pytest.fixture(scope="module")
def spec(page):
    return page.parse_spec(page.SPEC_PATH)


def _fails(results):
    return [(n, d) for n, ok, d in results if not ok]


def test_x_extraction(page, parsed, spec):
    from tests.validation.test_geo38_f_clinical_box_loaded import BRIEF_LATTICE_N

    res = page.check_extraction(parsed, spec)
    for n, ok, d in res:
        print(f"[ANS-11 X] {'PASS' if ok else 'FAIL'} {n}: {d}")
    assert not _fails(res), _fails(res)
    assert len(parsed["lattice"]) == BRIEF_LATTICE_N
    # The brief's cited lines (docs/planning/chunks/ANS-11.md, §9 item 833 "Code").
    lines = {
        k: parsed[k]["line"]
        for k in (
            "node",
            "s2",
            "y_in_ccw",
            "y_in_cw",
            "means",
            "c_star",
            "normalisers",
            "wrote",
            "price",
        )
    }
    assert lines == {
        "node": 12063,
        "s2": 12281,
        "y_in_ccw": 12285,
        "y_in_cw": 12286,
        "means": 12288,
        "c_star": 12307,
        "normalisers": 12310,
        "wrote": 12311,
        "price": 12312,
    }, lines
    assert [r["line"] for r in parsed["lattice"]] == list(range(12289, 12307))


def test_d_derived(page, parsed, spec):
    d = page.derived(parsed)
    res = page.check_derived(parsed, spec, d)
    for n, ok, det in res:
        print(f"[ANS-11 D] {'PASS' if ok else 'FAIL'} {n}: {det}")
    # (N-sense), printed only (rule (e)).
    print(
        f"[ANS-11 N-sense] P_acc cw {d['p_acc_cw']:.4e} W vs ccw {d['p_acc_ccw']:.4e} W "
        f"(printed, not asserted)"
    )
    assert not _fails(res), _fails(res)


def test_r_regeneration(page, parsed):
    committed = json.loads(page.METRICS_JSON.read_text())
    fresh = json.loads(
        json.dumps(page.build_metrics(parsed, page.derived(parsed)), indent=1, sort_keys=True)
    )
    assert fresh == committed
    render = page.render_comparison(committed)
    text = page.COMPARISON_MD.read_text()
    print(f"[ANS-11 R] render {len(render.encode())} B vs committed {len(text.encode())} B")
    assert render == text


def test_n_digit_log(page, spec, tmp_path):
    src = page._read_log(page.SOURCE_LOG).splitlines(keepends=True)
    target = "SAR10g 1.378896513e-04 W/kg, mass"
    idx = [i for i, ln in enumerate(src) if "[GEO-38 C6a] lattice x-105mm dz-50mm:" in ln]
    assert len(idx) == 1 and target in src[idx[0]]
    src[idx[0]] = src[idx[0]].replace(target, "SAR10g 1.379896513e-04 W/kg, mass")
    bad = tmp_path / "c6_digit.log"
    bad.write_text("".join(src))
    p = page.parse_c6_log(bad)
    dres = page.check_derived(p, spec, page.derived(p))
    xres = page.check_extraction(p, spec)
    d_names = [n for n, _ in _fails(dres)]
    x_names = [n for n, _ in _fails(xres)]
    print(f"[ANS-11 N-digit-log] (D) fails {d_names}; (X) fails {x_names}")
    assert "(D) c* SAR string == its lattice line's string" in d_names
    assert x_names


def test_n_drop(page, tmp_path):
    src = page._read_log(page.SOURCE_LOG).splitlines(keepends=True)
    idx = [i for i, ln in enumerate(src) if "[GEO-38 C6a] lattice x+0mm dz+50mm:" in ln]
    assert len(idx) == 1
    del src[idx[0]]
    bad = tmp_path / "c6_drop.log"
    bad.write_text("".join(src))
    with pytest.raises(ValueError, match=r"17 != BRIEF_LATTICE_N 18"):
        page.parse_c6_log(bad)


def test_n_doc(page, tmp_path):
    text = page.COMPARISON_MD.read_text()
    target = "2.963374"
    assert text.count(target) == 1
    bad = tmp_path / "COMPARISON.md"
    bad.write_text(text.replace(target, "2.963384"))
    render = page.render_comparison(json.loads(page.METRICS_JSON.read_text()))
    assert bad.read_text() != render
    assert text == render


# ---------------------------------------------------------------------------
# C36 (§9 item 871).
# ---------------------------------------------------------------------------
def _bcast_dir(tmp_path, writer):
    """Rank 0 writes the stubs under its ``tmp_path``; every rank reads that dir."""
    try:
        from mpi4py import MPI

        comm = MPI.COMM_WORLD
    except ImportError:  # pragma: no cover
        comm = None
    rank = comm.rank if comm is not None else 0
    d = str(tmp_path)
    if rank == 0:
        writer(tmp_path)
    if comm is not None:
        d = comm.bcast(d, root=0)
        comm.Barrier()
    return d


def _c6_stub(page, c_pf=185.7629):
    from tests.validation.test_geo38_f_clinical_box_loaded import _lattice

    labels = [lab for lab, _ in _lattice()]
    rec = {
        "c_pf": c_pf,
        "feeds": ["P33", "P41"],
        "az_feeds": [0.1, 0.2],
        "s2": [[[-0.41, 0.012], [0.031, 0.12]], [[0.031, 0.12], [-0.43, -0.021]]],
        "y_in": {"ccw": [[0.041, 0.0013], [0.077, -0.0021]], "cw": [[0.05, 0.001], [0.06, 0.002]]},
        "mean_b_plus": {"ccw": 2.2e-08, "cw": 7.1e-09},
        "mean_b_minus": {"ccw": 7.3e-09, "cw": 2.1e-08},
        "lattice": [
            {
                "label": lb,
                "c": [0.0, 0.0, 0.0],
                "b_plus_ccw": 1.0e-08 + 1e-10 * i,
                "sar10g": 1.0e-04 + 1e-06 * i,
            }
            for i, lb in enumerate(labels)
        ],
        "c_star": [-0.105, 0.0, -0.05],
        "star_label": labels[2],
        "sar_star": 1.4e-04,
        "sar_star_cw": 1.1e-04,
        "p_box": 3.1e-03,
        "p_src": 8.2e-02,
        "p_acc": 6.1e-02,
        "parts": [],
        "grid": {"star_label": "g017", "sar_star": 1.6e-04},
    }
    assert set(rec) == set(page.K8_JSON_KEYS) | {"grid"}
    return rec


def _write(dirpath, name, rec):
    (Path(dirpath) / name).write_text(json.dumps(rec))


# C37 (§9 item 983) retired ``test_c19_absent``: its premise ("entry 62 has not run")
# expired when entry 62 wrote its JSONs (known-issues 2026-11-14).  ``test_c37_custody``
# replaces it: the tracked-log route asserted always, (J62) when the host JSON exists.
@pytest.fixture(scope="module")
def block(page):
    return page.degree2_block()


def test_c37_custody(page, block):
    c19 = page.read_c19()
    assert c19["paths"] == [
        "/workspace/logs/geo38/c6-node-185.7629-d2.json",
        "/workspace/logs/geo36/i36-node-185.7629-d2.json",
    ]
    m = json.loads(page.METRICS_JSON.read_text())
    text = page.COMPARISON_MD.read_text()
    # The tracked-log route, always: the committed render is the c19=None render.
    assert m["degree2_c19"]["source_log"] == page.C19_LOG_BASENAME
    assert page.render_comparison(m) == text
    assert json.loads(json.dumps(block)) == m["degree2_c19"]
    cust = page.custody_c19(block, c19)
    if cust is None:
        print(f"[ANS-11 C37 J62] {c19['paths'][0]} absent; custody not run")
        return
    n, mism, pform = cust
    print(f"[ANS-11 C37 J62] n_fields {n}, mismatches {mism}, print_form_only {pform}")
    assert mism == []
    assert c19["i36"] is not None  # validated by read_c19 (degree 2, node); no row


def test_c37_lines(page, block):
    p2 = block["parsed"]
    lines = {
        k: p2[k]["line"]
        for k in (
            "node",
            "s2",
            "y_in_ccw",
            "y_in_cw",
            "means",
            "c_star",
            "normalisers",
            "wrote",
            "price",
        )
    }
    print(f"[ANS-11 C37 lines] {lines}; grid :{block['grid']['line']}")
    assert lines == {
        "node": 12160,
        "s2": 12379,
        "y_in_ccw": 12383,
        "y_in_cw": 12384,
        "means": 12387,
        "c_star": 12406,
        "normalisers": 12410,
        "wrote": 12411,
        "price": 12412,
    }, lines
    assert [r["line"] for r in p2["lattice"]] == list(range(12388, 12406))
    assert block["grid"]["line"] == 12407 and block["grid"]["n"] == "90"
    assert p2["node"]["c_pf"] == page.C19_NODE_PF


def test_c37_d2(page, block):
    res = page.check_d2(block)
    for n, ok, det in res:
        print(f"[ANS-11 C37 D2] {'PASS' if ok else 'FAIL'} {n}: {det}")
    s = block["parsed"]["s2"]
    print(f"[ANS-11 C37 D2] (printed) S12 {s['s12']} vs S21 {s['s21']}")
    assert not _fails(res), _fails(res)
    assert len(res) == 7


def test_c37_x2(page, block):
    text = page.COMPARISON_MD.read_text()
    miss = page.x2_misses(text, block)
    print(f"[ANS-11 C37 X2] misses {miss}")
    assert miss == []


def test_c37_moves_rule(page, block):
    for mv in block["moves"]:
        print(
            f"[ANS-11 C37 move] {mv['row']}: {mv['move']}{' (per-watt)' if mv['per_watt'] else ''}"
        )
    rule = block["order_rule"]
    print(f"[ANS-11 C37 order rule] over {rule['rows_over']} => {rule['sentence']}")
    assert len(block["moves"]) == 16
    assert len(rule["per_watt_rows"]) == 4
    # The verdict is printed, not asserted; only its consistency with the moves is.
    over = [
        mv["row"]
        for mv in block["moves"]
        if mv["per_watt"] and abs(float(mv["move"])) > page.ORDER_RULE_REL
    ]
    assert rule["rows_over"] == over
    assert rule["sentence"] == (page.ORDER_RULE_FIRES if over else page.ORDER_RULE_QUIET)
    assert rule["sentence"] in page.COMPARISON_MD.read_text().lower()


def _c19_copy(page, tmp_path, name, edit):
    def w(d):
        src = page._read_log(page.C19_LOG).splitlines(keepends=True)
        (Path(d) / name).write_text("".join(edit(src)))

    return Path(_bcast_dir(tmp_path, w)) / name


def test_c37_n_j62(page, block, tmp_path):
    c19 = page.read_c19()
    if c19["c6"] is None:
        pytest.skip("host C6 JSON absent (a clone); (N-J62) needs it")
    rec = json.loads(json.dumps(c19["c6"]))
    lab = rec["lattice"][4]["label"]
    rec["lattice"][4]["b_plus_ccw"] *= 1.0 + 1e-6
    d = _bcast_dir(tmp_path, lambda p: _write(p, "c6-node-185.7629-d2.json", rec))
    bad = page.read_c19(c6_dir=d, i36_dir=d)
    n, mism, _ = page.custody_c19(block, bad)
    print(f"[ANS-11 C37 N-J62] {lab} x (1 + 1e-6): mismatches {[m[0] for m in mism]}")
    assert [m[0] for m in mism] == [f"lattice.{lab}.b1p_ccw"]


def test_c37_n_x2(page, block, tmp_path):
    # One digit, the third after the point (rel ≈ 1.7e-4): a last-digit change (rel ≈ 1.7e-9)
    # would sit inside (D2)'s 1e-8 identity tolerance.
    tok = block["parsed"]["normalisers"]["p_acc"]
    new = tok[:4] + str((int(tok[4]) + 1) % 10) + tok[5:]
    assert new != tok and len(new) == len(tok)

    def edit(src):
        i = block["parsed"]["normalisers"]["line"] - 1
        assert src[i].count(f"P_acc (lumped line cut, ruling 2026-10-12) {tok} W") == 1
        src[i] = src[i].replace(
            f"P_acc (lumped line cut, ruling 2026-10-12) {tok} W",
            f"P_acc (lumped line cut, ruling 2026-10-12) {new} W",
        )
        return src

    bad = page.degree2_block(_c19_copy(page, tmp_path, "c19_pacc.log", edit))
    d_names = [n for n, _ in _fails(page.check_d2(bad))]
    x = page.x2_misses(page.COMPARISON_MD.read_text(), bad)
    print(f"[ANS-11 C37 N-X2] P_acc {tok} -> {new}: (D2) fails {d_names}; (X2) misses {x}")
    assert any("_p_acc_spec" in n for n in d_names)
    assert (page.Q_P, new) in x


def test_c37_n_grid(page, tmp_path):
    def edit(src):
        i = 12407 - 1
        assert "(C6b) grid n 90 (72 new balls" in src[i]
        return src[: i + 1] + [src[i]] + src[i + 1 :]

    with pytest.raises(ValueError, match=r"grid result lines: 2 matches"):
        page.parse_c19_grid(_c19_copy(page, tmp_path, "c19_grid2.log", edit))


def test_c37_n_node(page, tmp_path):
    def edit(src):
        i = 12160 - 1
        assert "(Z) node 185.7629 pF" in src[i]
        src[i] = src[i].replace("(Z) node 185.7629 pF", "(Z) node 193.5996 pF")
        return src

    with pytest.raises(ValueError, match=r"node 193\.5996 pF != C19_NODE_PF 185\.7629"):
        page.degree2_block(_c19_copy(page, tmp_path, "c19_node.log", edit))


def test_c19_stub(page, spec, tmp_path):
    def w(d):
        _write(d, "c6-node-185.7629-d2.json", _c6_stub(page))
        _write(d, "i36-node-185.7629-d2.json", {"c_pf": 185.7629, "degree": 2})

    d = _bcast_dir(tmp_path, w)
    c19 = page.read_c19(c6_dir=d, i36_dir=d)
    assert c19["c6"] is not None and c19["i36"] == {"c_pf": 185.7629, "degree": 2}
    cells = page.degree2_cells(c19)
    empty = [q for q, c in cells.items() if not c]
    print(f"[ANS-11 C19-stub] {len(cells) - len(empty)} cells filled; empty {empty}")
    assert empty == [page.Q_CU]
    m = json.loads(page.METRICS_JSON.read_text())
    text = page.COMPARISON_MD.read_text()
    stub_render = page.render_comparison(m, c19)
    assert stub_render != text
    got, ref = page.comparison_cells(stub_render), page.comparison_cells(text)
    assert {k: v[0] for k, v in got["main"].items()} == {k: v[0] for k, v in ref["main"].items()}
    assert got["lattice"] == ref["lattice"]
    assert {k: v[1] for k, v in got["main"].items() if v[1]} == {
        q: c for q, c in cells.items() if c
    }
    kres = page.check_k(spec, stub_render)
    assert not _fails(kres), _fails(kres)


def test_k(page, spec):
    res = page.check_k(spec, page.COMPARISON_MD.read_text())
    for n, ok, det in res:
        print(f"[ANS-11 K] {'PASS' if ok else 'FAIL'} {n}: {det}")
    assert not _fails(res), _fails(res)
    assert len(res) == 18
    assert len(page.K_NOT_CARRIED) == 6


def test_d_node(page, parsed):
    sn = page.parse_spec_node(page.SPEC_PATH)
    res, printed = page.check_d(sn, parsed)
    for n, ok, det in res:
        print(f"[ANS-11 D-node] {'PASS' if ok else 'FAIL'} {n}: {det}")
    for line in printed:
        print(f"[ANS-11 node-C, printed not asserted] {line}")
    assert not _fails(res), _fails(res)
    names = [n for n, _, _ in res]
    assert sum(n.startswith("(D-v1)") for n in names) == 9
    assert sum(n.startswith("(D-v2)") for n in names) == 6
    assert sn["node_14"][0] == sn["node_75"][0] == "185.7629"


def test_n_k(page, parsed, tmp_path):
    text = page.SPEC_PATH.read_text()
    target = "−4.8167140712e-01 + 7.9527084190e-04j"
    assert text.count(target) == 1
    bad = tmp_path / "SPEC.md"
    bad.write_text(text.replace(target, "−4.8167140713e-01 + 7.9527084190e-04j"))
    sp = page.parse_spec(bad)
    k_names = [n for n, _ in _fails(page.check_k(sp, page.COMPARISON_MD.read_text()))]
    x_names = [n for n, _ in _fails(page.check_extraction(parsed, sp))]
    print(f"[ANS-11 N-K] (K) fails {k_names}; (X) fails {x_names}")
    assert k_names == ["(K) s2.s11"]
    assert x_names == ["(X) S s11"]


def test_n_d(page):
    sn = page.parse_spec_node(page.SPEC_PATH)
    d = page.node_figures(193.5996e-12)
    print(f"[ANS-11 N-D] v2 R_p string {sn['v2_rp'][0]} vs R_p at 193.5996 pF {d['r_p']:.1f}")
    assert f"{d['r_p']:.1f}" == "825.0"
    assert not page._at(d["r_p"], sn["v2_rp"][0])


def test_n_stub_key(page, tmp_path):
    rec = _c6_stub(page)
    del rec["sar_star_cw"]
    d = _bcast_dir(tmp_path, lambda p: _write(p, "c6-node-185.7629-d2.json", rec))
    with pytest.raises(ValueError, match=r"c6-node-185\.7629-d2\.json.*missing \['sar_star_cw'\]"):
        page.read_c19(c6_dir=d, i36_dir=d)


def test_n_stub_node(page, tmp_path):
    rec = _c6_stub(page, c_pf=193.5996)
    d = _bcast_dir(tmp_path, lambda p: _write(p, "c6-node-185.7629-d2.json", rec))
    with pytest.raises(ValueError, match=r"c6-node-185\.7629-d2\.json: key c_pf 193\.5996"):
        page.read_c19(c6_dir=d, i36_dir=d)


def test_n_stub_deg(page, tmp_path):
    d = _bcast_dir(
        tmp_path, lambda p: _write(p, "i36-node-185.7629-d2.json", {"c_pf": 185.7629, "degree": 1})
    )
    with pytest.raises(ValueError, match=r"i36-node-185\.7629-d2\.json: key degree 1"):
        page.read_c19(c6_dir=d, i36_dir=d)


# ---------------------------------------------------------------------------
# C43 (§9 item 880): SPEC v2's node section and the guide's section 4.
# ---------------------------------------------------------------------------
GUIDE = CASE / "11_f_clinical_1p5t_box_loaded_64MHz.md"
C43_HEADING = "## 4. SPEC v2's node (the capacitor the operator enters)"
C43_NUM = r"[-+−]?\d+\.\d+(?:e[-+]\d+)?"


def _guide_section4(text: str) -> str:
    import re

    assert text.count(C43_HEADING) == 1
    body = text.split(C43_HEADING, 1)[1]
    m = re.search(r"^## ", body, flags=re.M)
    return body[: m.start()] if m else body


def _guide_misses(text: str, lines: list) -> list:
    import re

    out = "\n".join(lines)
    return [t for t in re.findall(C43_NUM, _guide_section4(text)) if t not in out]


def test_c43_node_section(page):
    lines = page.node_section_lines(page.parse_spec_node(page.SPEC_PATH))
    for ln in lines:
        print(f"[ANS-11 C43] {ln}")
    verdicts = [ln for ln in lines if ln.startswith(("[PASS]", "[FAIL]"))]
    assert len(verdicts) == 7, verdicts  # node line + R_p, Q_C, C_p, -Im Z_RC, Re/-Im Z_th
    assert all(ln.startswith("[PASS]") for ln in verdicts), verdicts
    out = "\n".join(lines)
    for s in (
        "896.1 Ω",
        "66.9 vs",
        "185.7215 pF",
        "0.2 − 13.386935347j",
        "3.505604063 − 12.39875443j",
        "185.7629 pF",
    ):
        assert s in out, s
    text = GUIDE.read_text()
    assert text.count(C43_HEADING) == 1
    toks = __import__("re").findall(C43_NUM, _guide_section4(text))
    misses = _guide_misses(text, lines)
    print(f"[ANS-11 C43] guide section 4: {len(toks)} numeric tokens, misses {misses}")
    assert toks and not misses, misses


def test_c43_n_guide(page, tmp_path):
    text = GUIDE.read_text()
    sec = _guide_section4(text)
    assert sec.count("185.7215") == 1
    bad = text.replace(sec, sec.replace("185.7215", "185.7216"))
    lines = page.node_section_lines(page.parse_spec_node(page.SPEC_PATH))
    misses = _guide_misses(bad, lines)
    print(f"[ANS-11 N-guide] misses {misses}")
    assert misses == ["185.7216"]


def test_c43_n_node(page, tmp_path):
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    text = page.SPEC_PATH.read_text()
    bad = Path(comm.bcast(str(tmp_path) if comm.rank == 0 else None, root=0)) / "SPEC.md"
    if comm.rank == 0:
        bad.write_text(text.replace("C = 185.7629 pF", "C = 193.5996 pF"))
    comm.Barrier()
    lines = page.node_section_lines(page.parse_spec_node(bad))
    rp = [ln for ln in lines if "(C43 v2) R_p" in ln]
    print(f"[ANS-11 N-node] {rp}")
    assert len(rp) == 1 and rp[0].startswith("[FAIL]") and "SPEC v2 896.1" in rp[0]
    # a FAIL line shows the derived value at <= 7 digits; read it at the SPEC's one decimal
    shown = __import__("re").search(r"at the root = (\d+\.\d+)", rp[0]).group(1)
    assert f"{float(shown):.1f}" == "825.0", shown


# ---------------------------------------------------------------------------
# C49 (§9 item 931): the "ans12" case — ANS-12's COMPARISON.md from entry 58.
# ---------------------------------------------------------------------------
import re  # noqa: E402

SPEC_3T = (
    REPO / "examples" / "ansys_benchmarks" / "ans12_f_clinical_3t_box_loaded_128MHz" / "SPEC.md"
)


@pytest.fixture(scope="module")
def p58(page):
    return page.parse_c6_log(page.SOURCE_LOG_3T)


def _render_3t(page, parsed, case="ans12"):
    return page.render_comparison(
        page.build_metrics(parsed, page.derived(parsed), case=case), None, case=case
    )


def _x3t_misses(page, text, parsed) -> list:
    """(X3T): every degree-1 cell carries entry 58's parsed token string(s) for its row."""
    c = page.comparison_cells(text)
    m, d = c["main"], page.derived(parsed)
    s, y, mn, cs, nm = (parsed[k] for k in ("s2", "y_in_ccw", "means", "c_star", "normalisers"))
    want = {
        page.Q_S11: [s["s11"], f":{s['line']}"],
        page.Q_S12: [s["s12"], s["s21"]],
        page.Q_S22: [s["s22"]],
        page.Q_Y33: [y["p33_re"], y["p33_im"], f":{y['line']}"],
        page.Q_Y41: [y["p41_re"], y["p41_im"]],
        page.Q_LAT: [f"{d['lattice_mean_b1p']:.5e}"],
        page.Q_LAT_P: [f"{d['lattice_mean_b1p_per_sqrt_p_acc']:.4e}"],
        page.Q_CSTAR: [cs["label"], cs["r_m"], cs["dz_m"], cs["az_deg"], f":{cs['line']}"],
        page.Q_SAR: [cs["sar10g"], nm["sar_per_p_acc"], f":{nm['line']}"],
        page.Q_CW: [cs["sar10g_cw"]],
        page.Q_GRID: ["entry 58"],
        page.Q_MEAN: [mn["b1p_ccw"], f"{mn['n_used']} / {mn['n_total']}", f":{mn['line']}"],
        page.Q_MEAN_P: [nm["b1_per_sqrt_p_acc"]],
        page.Q_P: [nm["p_box"], nm["p_src"], nm["p_acc"]],
        page.Q_GEL: [nm["p_box_over_p_acc"]],
        page.Q_CU: ["entry 58", ":12229–12231"],
        page.Q_SENSE: [mn["sense"]],
    }
    assert len(want) == 17
    miss = [(q, t) for q, toks in want.items() for t in toks if t not in m.get(q, ("",))[0]]
    lines = {
        mm.group(1): mm.group(2)
        for mm in re.finditer(r"^\| (x[+-]\d+mm dz[+-]\d+mm) \| `:(\d+)` \|", text, re.M)
    }
    for r in parsed["lattice"]:
        got = c["lattice"].get(r["label"])
        if got != (r["b1p_ccw"], r["sar10g"]) or lines.get(r["label"]) != str(r["line"]):
            miss.append((r["label"], got, lines.get(r["label"])))
    if len(c["lattice"]) != 18:
        miss.append(("lattice rows", len(c["lattice"])))
    return miss


def test_c49_k_ans11_default_byte_identical(page, parsed):
    """(K): the default render of entry 52 == the committed COMPARISON.md; metrics.json too."""
    m = page.build_metrics(parsed, page.derived(parsed))
    js = json.dumps(m, indent=1, sort_keys=True) + "\n"
    text = page.COMPARISON_MD.read_text()
    print(
        f"[ANS-11 C49 K] metrics {len(js.encode())} B vs {len(page.METRICS_JSON.read_bytes())} B; "
        f"render {len(page.render_comparison(m).encode())} B vs {len(text.encode())} B"
    )
    assert js == page.METRICS_JSON.read_text()
    assert page.render_comparison(json.loads(js)) == text
    assert page.render_comparison(m, None, case="ans11") == text


def test_c49_x3t(page, p58):
    text = page.COMPARISON_MD_3T.read_text()
    miss = _x3t_misses(page, text, p58)
    print(f"[ANS-12 C49 X3T] misses {miss}")
    assert p58["log"] == page.SOURCE_LOG_3T_BASENAME and p58["node"]["c_pf"] == "44.0525"
    assert miss == []


def test_c49_tie(page, p58):
    t = page.tie_set(p58)
    spec = SPEC_3T.read_text()
    sm = re.search(
        r"Row 3's tie set on this column: x\+105 mm dz\+50 mm \((\S+)\), x\+105 mm\s+"
        r"dz−50 mm \((\S+), (−?\d+\.\d+) %\), x−105 mm dz−50 mm \((\S+), (−?\d+\.\d+) %,",
        spec,
    )
    assert sm, "SPEC tie sentence"
    pcts = [g.replace("−", "-") for g in (sm.group(3), sm.group(5))]
    print(f"[ANS-12 C49 Tie] members {t['members']}; outsider {t['outsider']}; SPEC {pcts}")
    assert [(a, c) for a, _, c in t["members"]] == [
        ("x+105mm dz+50mm", "0.00"),
        ("x+105mm dz-50mm", pcts[0]),
    ]
    assert (t["outsider"][0], t["outsider"][2]) == ("x-105mm dz-50mm", pcts[1])
    assert pcts == ["-0.82", "-1.97"]
    text = page.COMPARISON_MD_3T.read_text()
    sec = text.split("## Row 3 tie set", 1)[1].split("## Provenance", 1)[0]
    rows = re.findall(
        r"^\| (x[+-]\d+mm dz[+-]\d+mm) \| (\S+) \| (-?\d+\.\d+) \| ([a-z ]+) \|$", sec, re.M
    )
    assert rows == [
        (a, b, c, r)
        for (a, b, c), r in zip(
            t["members"] + [t["outsider"]], ["maximum", "member", "nearest outsider"]
        )
    ]
    assert "## Row 3 tie set" not in page.COMPARISON_MD.read_text()


def test_c49_blank(page):
    rows = [
        ln for ln in page.COMPARISON_MD_3T.read_text().splitlines() if re.match(r"^\| \d+ \| ", ln)
    ]
    cells = [[c.strip() for c in ln[1:-1].split(" |")] for ln in rows]
    assert len(cells) == 17 and all(len(c) == 5 for c in cells)
    assert all(c[3] == "" and c[4] == "" for c in cells), [c for c in cells if c[3] or c[4]]


def test_c49_r3t(page, p58):
    render = _render_3t(page, p58)
    text = page.COMPARISON_MD_3T.read_text()
    print(f"[ANS-12 C49 R3T] render {len(render.encode())} B vs committed {len(text.encode())} B")
    assert render == text
    assert "[SPEC.md](SPEC.md)" in text


def test_c49_n_case(page, parsed, p58):
    wrong = _render_3t(page, parsed, case="ans12")
    miss = _x3t_misses(page, wrong, p58)
    print(f"[ANS-12 C49 N-case] entry 52 as ans12: {len(miss)} X3T misses")
    assert miss
    assert _render_3t(page, p58, case="ans11") != page.COMPARISON_MD_3T.read_text()


def test_c49_n_tie(page, tmp_path):
    src = page._read_log(page.SOURCE_LOG_3T).splitlines(keepends=True)
    idx = [i for i, ln in enumerate(src) if "[GEO-38 C6a] lattice x+105mm dz-50mm:" in ln]
    assert len(idx) == 1
    tok = re.search(r"SAR10g (\S+) W/kg, mass", src[idx[0]]).group(1)
    src[idx[0]] = src[idx[0]].replace(f"SAR10g {tok} W/kg", "SAR10g 1.714192670e-04 W/kg")
    bad = tmp_path / "c6_3t_tie.log"
    bad.write_text("".join(src))
    t = page.tie_set(page.parse_c6_log(bad))
    print(f"[ANS-12 C49 N-tie] {tok} -> 1.714192670e-04: members {t['members']}")
    assert len(t["members"]) == 1


def test_c49_n_doc(page, p58, tmp_path):
    text = page.COMPARISON_MD_3T.read_text()
    tok = p58["normalisers"]["p_acc"]
    assert text.count(tok) == 1
    last = tok[tok.index("e") - 1]
    bad = tmp_path / "COMPARISON.md"
    bad.write_text(text.replace(tok, tok.replace(f"{last}e", f"{(int(last) + 1) % 10}e")))
    assert bad.read_text() != _render_3t(page, p58)
    assert text == _render_3t(page, p58)


# ---------------------------------------------------------------------------
# C54 (§9 item 959): ANS-13's thermal case — (X13), (Probes), (Blank), (R13), (K)
# and the controls (N-case), (N-digit), (N-doc).  X13 reads the logs' text directly
# (its own needles), independent of the page's regexes.
# ---------------------------------------------------------------------------
import ast  # noqa: E402

N_ROWS_13 = 13  # row 1, six probes + c* + cw (row 2), rows 3, 4, 5, 6


@pytest.fixture(scope="module")
def logs13():
    from tests.validation.test_geo38_c39_thermal_reader import ENTRY59_GLOB, ENTRY60_GLOB, find_log

    out = {}
    for k, g in (("59", ENTRY59_GLOB), ("60", ENTRY60_GLOB)):
        hit = find_log(g)
        assert hit is not None and hit[1] == "tracked", (k, hit)
        out[k] = hit[0]
    return out


def _on(lines, cite, *needles):
    return all(n in lines[cite - 1] for n in needles)


def _x13_misses(page, text, log59, log60) -> list:
    """(X13): every degree-1 cell's token sits on its cited line of entry 59's / 60's
    log (the right sense and quantity); derived cells recomputed from those tokens."""
    from tests.validation.test_mat21_box_closed_form import RHO_C

    logs = {"59": page._read_log(log59).splitlines(), "60": page._read_log(log60).splitlines()}
    miss = []
    for (row, q), (c59, c60, _d2, _aed) in page.comparison_cells_thermal(text).items():
        for k, c in (("59", c59), ("60", c60)):
            L = logs[k]
            cites = [int(n) for n in re.findall(r"`:(\d+)`", c)]
            toks = re.findall(r"[+-]?\d\.\d+e[+-]\d+", c)
            ok = False
            try:
                pm = re.match(r"θ\(t\\\*\) at (P\d .+?) per deposited watt \[K\] \(ccw\)$", q)
                if row == 1:
                    ok = _on(L, cites[0], f"[gel twin ccw] (B) int_sub q {toks[0]} vs")
                elif row == 2 and pm:
                    ok = _on(
                        L,
                        cites[0],
                        f"[gel twin ccw] probe {pm.group(1)} (",
                        f"sense (1,-1j): theta(t*) {toks[0]} K",
                    )
                elif row == 2 and q.startswith(("θ(t\\*) at c\\*", "cw ")):
                    word, tup = ("cw", "+1j") if q.startswith("cw ") else ("ccw", "-1j")
                    hot, src = cites[0], cites[1]
                    th = re.search(
                        r"\[GEO-38 P61\] hotspot c\* .*sense \(1,"
                        + re.escape(tup)
                        + r"\): theta\(t\*\) per P_acc watt (\S+) K",
                        L[hot - 1],
                    ).group(1)
                    fa = re.search(
                        r"\("
                        + word
                        + r", heated\): .*factor P_gel/P_acc \("
                        + word
                        + r"\) (\S+)\s*$",
                        L[src - 1],
                    ).group(1)
                    ok = toks[0] == f"{float(th) / float(fa):.6e}"
                    if word == "ccw":
                        ok = ok and _on(
                            L,
                            cites[2],
                            "[gel twin ccw] 1 g ball",
                            f"sense (1,-1j): theta(t*) {toks[1]} K",
                        )
                    else:
                        ok = ok and toks[1] == fa
                elif row == 3:
                    ok = (
                        _on(L, cites[0], "[gel twin ccw] sense (1,-1j):", f"V_sub {toks[1]} m^3")
                        and toks[0] == f"{1.0 / (RHO_C * float(toks[1])):.4e}"
                        and _on(L, cites[1], "[gel twin ccw] CN N=", f"(S2) {toks[2]} (band")
                    )
                elif row == 4:
                    ok = _on(L, cites[0], "[gel twin ccw] CN N=", f"+ surface {toks[0]} W:")
                elif row == 5:
                    ok = c == "P87 prints it"
                elif row == 6:
                    ok = _on(
                        L, cites[0], "[GEO-38 P61] sources, sense (1,-1j) (ccw, heated)"
                    ) and L[cites[0] - 1].rstrip().endswith(f"(ccw) {toks[0]}")
            except (IndexError, AttributeError, ValueError):
                ok = False
            if not ok:
                miss.append((row, q, k))
    return miss


def test_c54_x13(page, logs13):
    text = page.COMPARISON_MD_13.read_text()
    cells = page.comparison_cells_thermal(text)
    miss = _x13_misses(page, text, logs13["59"], logs13["60"])
    print(f"[ANS-13 C54 X13] {len(cells)} rows; misses {miss}")
    assert len(cells) == N_ROWS_13
    assert sorted({r for r, _ in cells}) == [1, 2, 3, 4, 5, 6]
    assert miss == []


def test_c54_probes(page):
    import tests.validation.test_geo38_f_clinical_box_heating as p61

    cells = page.comparison_cells_thermal(page.COMPARISON_MD_13.read_text())
    labs = [
        m.group(1)
        for r, q in cells
        if r == 2 and (m := re.match(r"θ\(t\\\*\) at (P\d .+?) per deposited watt", q))
    ]
    print(f"[ANS-13 C54 Probes] {labs} vs {list(p61._probes())}")
    assert labs == list(p61._probes())


def test_c54_blank(page):
    cells = page.comparison_cells_thermal(page.COMPARISON_MD_13.read_text())
    assert all(v[2] == "" and v[3] == "" for v in cells.values())
    assert [v[:2] for (r, _), v in cells.items() if r == 5] == [("P87 prints it", "P87 prints it")]


def test_c54_r13(page):
    render = page.render_ans13()
    text = page.COMPARISON_MD_13.read_text()
    print(f"[ANS-13 C54 R13] render {len(render.encode())} B vs committed {len(text.encode())} B")
    assert render == text
    assert "[SPEC.md](SPEC.md)" in text


def test_c54_k_cases_untouched(page, parsed, p58):
    """(K): the thermal case sits beside CASES; the EM renders equal the committed files."""
    assert sorted(page.CASES) == ["ans11", "ans12"] and sorted(page.THERMAL_CASES) == ["ans13"]
    m = page.build_metrics(parsed, page.derived(parsed))
    assert json.dumps(m, indent=1, sort_keys=True) + "\n" == page.METRICS_JSON.read_text()
    assert page.render_comparison(m) == page.COMPARISON_MD.read_text()
    assert _render_3t(page, p58) == page.COMPARISON_MD_3T.read_text()


def test_c54_lazy_imports(page):
    """No top-level import of C46, P61, MAT-21, C33 or the ans:13 page (the closure C33
    and mri:41 load is unchanged; no cycle)."""
    tree = ast.parse(Path(page.__file__).read_text())
    mods = []
    for node in tree.body:
        if isinstance(node, ast.ImportFrom):
            mods.append(node.module or "")
        elif isinstance(node, ast.Import):
            mods += [a.name for a in node.names]
    bad = [
        m
        for m in mods
        if any(
            s in m for s in ("c39_thermal_reader", "box_heating", "mat21_box", "c33_field", "ans13")
        )
    ]
    assert bad == []


def test_c54_n_case(page, logs13):
    swapped = page.render_ans13(log59=logs13["60"], log60=logs13["59"])
    miss = _x13_misses(page, swapped, logs13["59"], logs13["60"])
    print(f"[ANS-13 C54 N-case] 59 / 60 swapped: {len(miss)} X13 misses")
    assert miss
    assert swapped != page.COMPARISON_MD_13.read_text()


def test_c54_n_digit(page, logs13, tmp_path):
    src = page._read_log(logs13["60"]).splitlines(keepends=True)
    ln = src[12271 - 1]
    assert "[gel twin ccw] probe P3 " in ln
    tok = re.search(r"theta\(t\*\) (\S+) K", ln).group(1)
    pos = [i for i, ch in enumerate(tok) if ch.isdigit()][3]
    new = tok[:pos] + str((int(tok[pos]) + 1) % 10) + tok[pos + 1 :]
    src[12271 - 1] = ln.replace(f"theta(t*) {tok} K", f"theta(t*) {new} K")
    bad = tmp_path / "p33_3t_digit.log"
    bad.write_text("".join(src))
    miss = _x13_misses(page, page.COMPARISON_MD_13.read_text(), logs13["59"], bad)
    print(f"[ANS-13 C54 N-digit] :12271 {tok} -> {new}: misses {miss}")
    assert [(r, k) for r, _, k in miss] == [(2, "60")]
    assert miss[0][1].startswith("θ(t\\*) at P3 ")


def test_c54_n_doc(page, tmp_path):
    text = page.COMPARISON_MD_13.read_text()
    i = text.index("| 4 | ")
    j = next(n for n in range(i + 30, len(text)) if text[n].isdigit())
    bad = tmp_path / "COMPARISON.md"
    bad.write_text(text[:j] + str((int(text[j]) + 1) % 10) + text[j + 1 :])
    assert bad.read_text() != page.render_ans13()
    assert text == page.render_ans13()
