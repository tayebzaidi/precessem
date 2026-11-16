"""`ANS-12` — the `ans:12` record page (§9 item 930, §10 C48;
docs/planning/chunks/ANS-12.md "§9 item 930").

Our degree-1 column of the F-clinical-3T box-loaded SPEC, read from entry 58's
tracked C6a log through C33's reader and the `ans:11` page's parser; no solve, no mesh.

(X) extraction: lattice labels == `_lattice()`; node == ``ENTRY58_NODE_PF``; every
    SPEC token (parsed from SPEC.md, never retyped) == its parsed field; the 36
    six-digit lattice entries; the header gate tokens on their C6a lines; every SPEC
    token matched exactly once.
(D) C33's ``check_side`` on entry 58, (D-lat), (D-latP).
(R) the lattice table rendered from the log == SPEC.md's bytes; the row-3 tie set.
(R-cmp) absent until item 931 (C49) lands.
Controls (asserted): (N-digit-log), (N-drop), (N-doc), (N-tie).  (N-sense) is printed
only (rule (e)).
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
CASE = REPO / "examples" / "ansys_benchmarks" / "ans12_f_clinical_3t_box_loaded_128MHz"
SCRIPT = CASE / "12_f_clinical_3t_box_loaded_128MHz.py"


def _load():
    name = "ans12_record_page"
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
        print(f"[ANS-12 {tag}] {'PASS' if ok else 'FAIL'} {n}: {d}")


def test_x_extraction(page, ctx):
    res = page.check_x(ctx)
    _show("X", res)
    assert not _fails(res), _fails(res)
    p = ctx["p58"]
    # The brief's lines (docs/planning/chunks/ANS-12.md, §9 item 930, reading 2).
    lines = {k: p[k]["line"] for k in ("node", "s2", "y_in_ccw", "y_in_cw", "means", "c_star",
                                       "normalisers", "wrote", "price")}
    assert lines == {"node": 12010, "s2": 12228, "y_in_ccw": 12232, "y_in_cw": 12233,
                     "means": 12237, "c_star": 12256, "normalisers": 12259, "wrote": 12260,
                     "price": 12261}, lines
    assert ctx["c21"]["line"] == 12235
    assert [r["line"] for r in p["lattice"]] == list(range(12238, 12256))
    assert {n: line for n, line, _ in ctx["gates"]}["elapsed_s"] == 12350


def test_x_raw_multiplicity(page, ctx):
    """Reading 1: tokens matched once to their field; the raw multiplicity is printed."""
    text = ctx["page"]._read_log(ctx["log_path"])
    sp = ctx["spec"]
    for name, tok in (("centroid mean", sp["means"]["b1p_ccw"]),
                      ("c* SAR", sp["c_star"]["sar10g"]),
                      ("S12", sp["s2"]["s12"])):
        print(f"[ANS-12 X-mult] {name} {tok}: {text.count(tok)} occurrences in the log")
    assert text.count(sp["means"]["b1p_ccw"]) >= 2


def test_d_derived(page, ctx):
    res = page.check_d(ctx)
    _show("D", res)
    names = [n for n, _, _ in res]
    for want in ("(D-Pacc)", "(D-box)", "(D-sar)", "(D-cen)", "(D-argmax)", "(D-cstar)",
                 "(D-S12)", "(Node)", "(D-lat)", "(D-latP)"):
        assert any(want in n for n in names), want
    d = ctx["d"]
    print(f"[ANS-12 N-sense] P_acc cw {d['p_acc_cw']:.4e} W vs ccw {d['p_acc_ccw']:.4e} W, "
          f"ratio {d['p_acc_cw'] / d['p_acc_ccw']:.4f} (printed, not asserted)")
    assert not _fails(res), _fails(res)


def test_r_regeneration(page, ctx):
    res = page.check_r(ctx)
    _show("R", res)
    assert not _fails(res), _fails(res)
    t = page.tie_set(ctx["p58"])
    assert [a for a, _, _ in t["members"]] == ["x+105mm dz+50mm", "x+105mm dz-50mm"]
    assert t["outsider"][0] == "x-105mm dz-50mm"


def test_r_cmp_conditional(page, ctx):
    res, msg = page.check_r_cmp(ctx)
    print(f"[ANS-12 R-cmp] {msg}")
    if page.COMPARISON_MD.exists():
        assert msg == "live" and not _fails(res), (msg, res)
    else:
        assert msg == "absent — item 931 (C49)" and res == []


def test_n_digit_log(page, tmp_path):
    lines = page._page()._read_log(page.entry58_log()).splitlines(keepends=True)
    key = "[GEO-38 C6a] lattice x+105mm dz+50mm:"
    idx = [i for i, ln in enumerate(lines) if key in ln]
    assert len(idx) == 1 and idx[0] + 1 == 12253
    assert "SAR10g 1.748477085e-04 W/kg" in lines[idx[0]]
    lines[idx[0]] = lines[idx[0]].replace("SAR10g 1.748477085e-04", "SAR10g 1.749477085e-04")
    bad = tmp_path / "c6_3t_digit.log"
    bad.write_text("".join(lines))
    c = page.read_all(log_path=bad)
    d_f = [n for n, _ in _fails(page.check_d(c))]
    x_f = [n for n, _ in _fails(page.check_x(c))]
    r_f = [n for n, _ in _fails(page.check_r(c))]
    print(f"[ANS-12 N-digit-log] (D) {d_f}; (X) {x_f}; (R) {r_f}")
    assert any("(D-cstar)" in n for n in d_f)
    assert "(X) SPEC 36 six-digit lattice entries == f'{float(tok):.6e}'" in x_f
    assert "(R) lattice table rendered from the log == SPEC bytes" in r_f


def test_n_drop(page, tmp_path):
    lines = page._page()._read_log(page.entry58_log()).splitlines(keepends=True)
    idx = [i for i, ln in enumerate(lines) if "[GEO-38 C6a] lattice x+0mm dz-50mm:" in ln]
    assert len(idx) == 1
    del lines[idx[0]]
    bad = tmp_path / "c6_3t_drop.log"
    bad.write_text("".join(lines))
    with pytest.raises(ValueError, match="lattice lines: 17 != BRIEF_LATTICE_N 18"):
        page._page().parse_c6_log(bad)


def test_n_doc(page, ctx, tmp_path):
    text = page.SPEC_PATH.read_text()
    old = "| 0 | 7.447523e-09 / 1.313789e-04 |"
    assert text.count(old) == 1
    spec_bad = tmp_path / "SPEC.md"
    spec_bad.write_text(text.replace(old, "| 0 | 7.447523e-09 / 1.313779e-04 |"))
    c = dict(ctx, spec=page.parse_spec_3t(spec_bad))
    r_f = [n for n, _ in _fails(page.check_r(c))]
    print(f"[ANS-12 N-doc] (R) fails {r_f}")
    assert "(R) lattice table rendered from the log == SPEC bytes" in r_f


# ---------------------------------------------------------------------------
# C58 (§9 item 989): the entry-66 degree-2 reader; (K) on entry 62, (Pend) on empty
# tmp roots (never on the real tree), (Stub66) and the controls.  Rank 0 writes every
# tmp file and broadcasts the path (C45's ``_on_root``, own copy).
# ---------------------------------------------------------------------------
def _comm():
    from mpi4py import MPI

    return MPI.COMM_WORLD


def _on_root(fn, *args, **kwargs):
    """Run ``fn`` on rank 0 only, bcast its result or its exception to every rank."""
    comm = _comm()
    out = None
    if comm.rank == 0:
        try:
            out = ("ok", fn(*args, **kwargs))
        except BaseException as exc:  # re-raised on every rank below
            out = ("err", exc)
    out = comm.bcast(out, root=0)
    if out[0] == "err":
        raise out[1]
    return out[1]


def _root_dir(tmp_path):
    return Path(_on_root(lambda: str(tmp_path)))


def _ledger_header(page):
    return [ln for ln in page.XXL_LEDGER.read_text().splitlines()
            if ln.startswith("| Date (UTC) |") or ln.startswith("|---|")]


def _write(path, text):
    def w():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        return str(path)
    return Path(_on_root(w))


def _e62_lines(page):
    return page._page()._read_log(page.TRACKED_LOG_DIR / page.ENTRY62_LOG_BASENAME).splitlines(
        keepends=True)


STUB66 = "20990101T000000Z_GEO-38-C19-3T.log"


def _stub_row(name=STUB66):
    return f"| 2099-01-01 | GEO-38-C19-3T | `{name}` | 16 | stub | — | — | stub |"


def test_c58_k(page):
    res = page.degree2_fields(page.TRACKED_LOG_DIR / page.ENTRY62_LOG_BASENAME)
    mf, blk = page.fields_from_metrics()
    out = page.check_k_d2(res, mf, blk)
    _show("C58 K", out)
    print(f"[ANS-12 C58 K] first {res['lines'][0]}, last {res['lines'][-1]}; grid_lines "
          f"{res['grid_lines']}")
    assert not _fails(out), _fails(out)
    assert len(res["fields"]) == page.D2_FIELDS_N == 53
    assert res["node"] == ("185.7629", 12160) and res["grid_lines"] == 1


def test_c58_pend(page, tmp_path):
    root = _root_dir(tmp_path)
    logs = root / "logs"
    _on_root(lambda: logs.mkdir())
    hdr = _ledger_header(page)
    assert len(hdr) == 2, hdr
    led = _write(root / "ledger_hdr.md", "\n".join(hdr) + "\n")
    assert page.entry66_log(led, logs) is None
    row27 = [ln for ln in page.XXL_LEDGER.read_text().splitlines()
             if ln.startswith("| 2026-11-14 | GEO-38-C19 |")]
    assert len(row27) == 1
    led2 = _write(root / "ledger_row27.md", "\n".join(hdr + row27) + "\n")
    assert page.entry66_scan(led2) == (0, [])
    assert page.entry66_log(led2, logs) is None


def test_c58_stub66(page, tmp_path):
    root = _root_dir(tmp_path)
    led = _write(root / "ledger.md", "\n".join(_ledger_header(page) + [_stub_row()]) + "\n")
    _write(root / "logs" / STUB66, "".join(_e62_lines(page)))
    path = page.entry66_log(led, root / "logs")
    assert path == root / "logs" / STUB66
    res = page.degree2_fields(path)
    ref = page.degree2_fields(page.TRACKED_LOG_DIR / page.ENTRY62_LOG_BASENAME)
    print(f"[ANS-12 C58 Stub66] {len(res['fields'])} fields, node {res['node']}")
    assert len(res["fields"]) == 53 and res["fields"] == ref["fields"]
    assert res["node"][0] == "185.7629" and res["log"] == STUB66


def test_c58_n_twin(page, tmp_path):
    root = _root_dir(tmp_path)
    logs = root / "logs"
    _on_root(lambda: logs.mkdir())
    other = "20990108T000000Z_GEO-38-C19-3T.log"
    led = _write(root / "twin.md", "\n".join(_ledger_header(page) + [_stub_row(), _stub_row(other)]))
    with pytest.raises(ValueError, match="name 2 distinct logs"):
        page.entry66_log(led, logs)
    _write(logs / STUB66, "".join(_e62_lines(page)))
    led2 = _write(root / "repeat.md", "\n".join(_ledger_header(page) + [_stub_row(), _stub_row()]))
    assert page.entry66_scan(led2) == (2, [STUB66])
    assert page.entry66_log(led2, logs) == logs / STUB66


def test_c58_n_absent(page, tmp_path):
    root = _root_dir(tmp_path)
    logs = root / "logs"
    _on_root(lambda: logs.mkdir())
    led = _write(root / "ledger.md", "\n".join(_ledger_header(page) + [_stub_row()]))
    with pytest.raises(ValueError, match=STUB66.replace(".", r"\.")):
        page.entry66_log(led, logs)


def test_c58_n_digit(page, tmp_path):
    lines = _e62_lines(page)
    key = "[GEO-38 C6a] lattice x-105mm dz-250mm:"
    idx = [i for i, ln in enumerate(lines) if key in ln]
    assert len(idx) == 1 and idx[0] + 1 == 12388
    assert lines[idx[0]].count("7.679531417e-05") == 1
    lines[idx[0]] = lines[idx[0]].replace("7.679531417e-05", "7.679531418e-05")
    bad = _write(_root_dir(tmp_path) / "c19_digit.log", "".join(lines))
    mf, blk = page.fields_from_metrics()
    res = page.degree2_fields(bad)
    want = dict(mf)
    mism = [n for n, t in res["fields"] if want[n] != t]
    print(f"[ANS-12 C58 N-digit] mismatches {mism}")
    assert mism == ["lattice.x-105mm dz-250mm.sar10g"]
    k = {n: ok for n, ok, _ in page.check_k_d2(res, mf, blk)}
    assert k["(K) every field token == metrics.json degree2_c19 string"] is False


def test_c58_n_grid2(page, tmp_path):
    lines = _e62_lines(page)
    assert page._page()._GRID_RE.search(lines[12406])
    lines.insert(12407, lines[12406])
    bad = _write(_root_dir(tmp_path) / "c19_grid2.log", "".join(lines))
    with pytest.raises(ValueError, match="grid result lines: 2"):
        page.degree2_fields(bad)


def test_c58_nogrid(page, tmp_path):
    lines = _e62_lines(page)
    assert page._page()._GRID_RE.search(lines[12406])
    del lines[12406]
    bad = _write(_root_dir(tmp_path) / "c19_nogrid.log", "".join(lines))
    res = page.degree2_fields(bad)
    out = page.render_d2(res)
    print("\n".join(out[:3]))
    assert len(res["fields"]) == 51 and res["grid_lines"] == 0
    assert not [n for n, _ in res["fields"] if n.startswith("grid.")]
    assert "[ANS-12 page] degree 2: grid: not run" in out


def test_n_tie(page, ctx, tmp_path):
    src = page._page()._read_log(page.entry58_log())
    assert src.count("1.734192670e-04") == 1
    bad = tmp_path / "c6_3t_tie.log"
    bad.write_text(src.replace("1.734192670e-04", "1.714192670e-04"))
    p = page._page().parse_c6_log(bad)
    t = page.tie_set(p)
    print(f"[ANS-12 N-tie] members {t['members']}; outsider {t['outsider']}")
    assert len(t["members"]) == 1 and t["members"][0][0] == "x+105mm dz+50mm"
    r_f = [n for n, _ in _fails(page.check_r(dict(ctx, p58=p)))]
    assert "(R) tie set (SAR10g within 1 % of the lattice max) == SPEC's members" in r_f
