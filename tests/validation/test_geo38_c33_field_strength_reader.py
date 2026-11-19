"""`GEO-38` C44 — C33's 1.5 T / 3 T field-strength reader (§9 item 905;
docs/planning/chunks/GEO-38.md "§9 item 905").  No solve, no mesh.

C33 reads entries 52 and 58: the 3 T / 1.5 T ratios of the lattice mean and the
centroid mean per sqrt(P_acc), of SAR10g(c*) per P_acc and of P_box / P_acc, the two
hotspot locations, Im Y_in at the two roots.  The host JSONs C33 names live under the
gitignored ``/logs/`` and are not reproducible from a clone, so this reader reads the
**tracked** C6a logs through the ``ans:11`` page's parser (loaded under the same
``sys.modules`` name as ``tests/unit/test_ans11_record_page.py``, never a re-implemented
regex); the host JSON is a custody cross-check only.

Anchors: (X52) parse of entry 52's log, node 193.5996; (D52) the spec-free derived
checks on entry 52; (J52) the host JSON's fields format to their log strings; (Self)
a renamed 3 T copy of entry 52's log reads with ratios 1; (Pend) entry 58 absent.
Controls: (N-digit), (N-node), (N-two), (N-json), (N-both).

Knob ``FEM_EM_GEO38_C33``: unset skip; ``1`` entry 52 required, entry 58 optional
(nothing written); ``both`` both required, rank 0 writes ``c33-1p5t-3t.json`` after
every check holds; anything else ``ValueError``.  No log token is a literal here: every
number checked is parsed at run time.

C33c (§9 item 979): ``custody_json`` formats the JSON's c* dz signed (``+.6f``), as
C6a prints it, so a positive dz (entry 58's) no longer reads a false mismatch;
(J58) and (N-sign) are its nodes, gated by the same knob.
"""

from __future__ import annotations

import importlib.util
import json
import math
import os
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
_PAGE_SCRIPT = (
    REPO
    / "examples"
    / "ansys_benchmarks"
    / "ans11_f_clinical_1p5t_box_loaded_64MHz"
    / "11_f_clinical_1p5t_box_loaded_64MHz.py"
)

C33_ENV = "FEM_EM_GEO38_C33"
#: Entry 58's XL_CHUNK (xl-pending.md entry 58 command).
ENTRY58_GLOB = "*_GEO-38-C6-3T.log"
#: Entry 57's crossing at four decimals, as entry 58's queue file carries it.
ENTRY58_NODE_PF = "44.0525"
ENTRY52_NODE_PF = "193.5996"
TRACKED_LOG_DIR = REPO / "docs" / "testing" / "logs"
HELD_LOG_DIR = REPO / "logs" / "leak-held"
HOST_JSON_DIR = Path("/workspace/logs/geo38")
C33_JSON_NAME = "c33-1p5t-3t.json"


# ---------------------------------------------------------------------------
# Loader (the unit test's, verbatim in effect: one shared module object).
# ---------------------------------------------------------------------------
def _load_page():
    name = "ans11_record_page"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _PAGE_SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


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


# ---------------------------------------------------------------------------
# Pure readers.
# ---------------------------------------------------------------------------
def find_log(glob, dirs=None):
    """``(path, "tracked" | "held")`` of the one log matching ``glob`` (either suffix)
    under ``docs/testing/logs/`` then ``logs/leak-held/`` (read only), or ``None``.
    Two matches raise ``ValueError``."""
    if dirs is None:
        dirs = ((TRACKED_LOG_DIR, "tracked"), (HELD_LOG_DIR, "held"))
    hits = []
    for d, kind in dirs:
        d = Path(d)
        if not d.is_dir():
            continue
        for p in sorted(set(d.glob(glob)) | set(d.glob(glob + ".gz"))):
            hits.append((p, kind))
    if len(hits) > 1:
        raise ValueError(f"{glob}: {len(hits)} matches {[str(p) for p, _ in hits]}")
    if not hits:
        return None
    p, kind = hits[0]
    # Cite the base name; the page's `_read_log` falls back to the `.gz` itself.
    return (Path(str(p).removesuffix(".gz")), kind)


def side(parsed):
    """The C33 fields of one side: strings as printed, plus floats from ``derived``."""
    page = _load_page()
    d = page.derived(parsed)
    nm, cs = parsed["normalisers"], parsed["c_star"]
    return {
        "log": parsed["log"],
        "node_pf": parsed["node"]["c_pf"],
        "lattice_mean_per_sqrt_p_acc": d["lattice_mean_b1p_per_sqrt_p_acc"],
        "lattice_mean_b1p": d["lattice_mean_b1p"],
        "centroid_mean_per_sqrt_p_acc": nm["b1_per_sqrt_p_acc"],
        "centroid_mean_b1p": parsed["means"]["b1p_ccw"],
        "sar_per_p_acc": nm["sar_per_p_acc"],
        "c_star": {k: cs[k] for k in ("label", "r_m", "dz_m", "az_deg", "sar10g", "sar10g_cw")},
        "y_in": {
            sense: {
                f"{p}_{part}": parsed[f"y_in_{sense}"][f"{p}_{part}"]
                for p in ("p33", "p41")
                for part in ("re", "im")
            }
            for sense in ("cw", "ccw")
        },
        "p_box_over_p_acc": nm["p_box_over_p_acc"],
        "p_acc": nm["p_acc"],
    }


def _half_ulp_rel(tok):
    """Half a unit in the printed token's last digit, relative to its value."""
    t = tok.strip()
    mant, _, exp = t.lower().partition("e")
    dec = mant.split(".")[1] if "." in mant else ""
    abs_ulp = 0.5 * 10.0 ** ((int(exp) if exp else 0) - len(dec))
    return abs_ulp / abs(float(t))


def check_side(parsed, node_pf, name):
    """Spec-free derived checks on one side: ``[(check, ok, detail)]``, every check
    named with ``name``."""
    page = _load_page()
    import numpy as np

    res = []
    nm, cs, means = parsed["normalisers"], parsed["c_star"], parsed["means"]

    def rel_check(tag, val, tok, inputs):
        bound = _half_ulp_rel(tok) + sum(_half_ulp_rel(t) for t in inputs)
        err = abs(val - float(tok)) / abs(float(tok))
        res.append(
            (
                f"{name} {tag}",
                err <= bound,
                f"{val:.6e} vs printed {tok}: rel {err:.3e} (<= {bound:.3e})",
            )
        )

    s = parsed["s2"]
    s2 = np.array([[complex(s["s11"]), complex(s["s12"])], [complex(s["s21"]), complex(s["s22"])]])
    p_acc = page._p_acc_spec(s2, page.CCW, complex(parsed["node"]["v_s"]), page.Z0_OHM)
    err = abs(p_acc - float(nm["p_acc"])) / abs(float(nm["p_acc"]))
    res.append(
        (
            f"{name} (D-Pacc)",
            err <= 1e-8,
            f"_p_acc_spec {p_acc:.6e} vs printed {nm['p_acc']}: rel {err:.3e} (<= 1e-08)",
        )
    )
    rel_check(
        "(D-box)",
        float(nm["p_box"]) / float(nm["p_acc"]),
        nm["p_box_over_p_acc"],
        (nm["p_box"], nm["p_acc"]),
    )
    rel_check(
        "(D-sar)",
        float(cs["sar10g"]) / float(nm["p_acc"]),
        nm["sar_per_p_acc"],
        (cs["sar10g"], nm["p_acc"]),
    )
    rel_check(
        "(D-cen)",
        float(means["b1p_ccw"]) / math.sqrt(float(nm["p_acc"])),
        nm["b1_per_sqrt_p_acc"],
        (means["b1p_ccw"], nm["p_acc"]),
    )
    d = page.derived(parsed)
    res.append(
        (
            f"{name} (D-argmax) argmax(lattice SAR10g) == c* label",
            d["argmax_label"] == cs["label"],
            f"argmax {d['argmax_label']} vs c* {cs['label']}",
        )
    )
    lat = {r["label"]: r for r in parsed["lattice"]}.get(cs["label"])
    res.append(
        (
            f"{name} (D-cstar) c* SAR string == its lattice line's string",
            lat is not None and lat["sar10g"] == cs["sar10g"],
            f"c* {cs['sar10g']} vs lattice {lat and lat['sar10g']}",
        )
    )
    res.append(
        (
            f"{name} (D-S12) S12 string == S21 string",
            s["s12"] == s["s21"],
            f"{s['s12']} vs {s['s21']}",
        )
    )
    res.append(
        (
            f"{name} (Node) node token == {node_pf}",
            parsed["node"]["c_pf"] == node_pf,
            f"log node {parsed['node']['c_pf']} vs {node_pf}",
        )
    )
    return res


def require_side(parsed, node_pf, name):
    res = check_side(parsed, node_pf, name)
    for n, ok, d in res:
        print(f"[GEO-38 C44] {'PASS' if ok else 'FAIL'} {n}: {d}")
    bad = [(n, d) for n, ok, d in res if not ok]
    if bad:
        raise AssertionError(f"{name}: {bad}")
    return res


def ratios(s15, s3):
    """3 T / 1.5 T ratios as floats (print at ``.6e``)."""
    return {
        k: float(s3[k]) / float(s15[k])
        for k in (
            "lattice_mean_per_sqrt_p_acc",
            "centroid_mean_per_sqrt_p_acc",
            "sar_per_p_acc",
            "p_box_over_p_acc",
        )
    }


def _cstar_str(sd):
    c = sd["c_star"]
    return f"{c['label']} (r {c['r_m']} m, dz {c['dz_m']} m, az {c['az_deg']} deg)"


def print_table(s15, s3=None):
    pend = "pending (xl entry 58)"
    r = ratios(s15, s3) if s3 is not None else None
    print(
        f"[GEO-38 C44] C33 table: 1.5 T ({s15['log']}, node {s15['node_pf']} pF) | "
        f"3 T ({s3['log'] + ', node ' + s3['node_pf'] + ' pF' if s3 else pend}) | ratio 3T/1.5T"
    )
    rows = (
        (
            "lattice mean |B1+|_ccw / sqrt(P_acc) [T/sqrt(W)]",
            "lattice_mean_per_sqrt_p_acc",
            lambda v: f"{v:.6e}",
        ),
        ("centroid mean |B1+|_ccw / sqrt(P_acc) [T/sqrt(W)]", "centroid_mean_per_sqrt_p_acc", str),
        ("SAR10g(c*) per P_acc [W/kg/W]", "sar_per_p_acc", str),
        ("P_box / P_acc", "p_box_over_p_acc", str),
    )
    for title, key, fmt in rows:
        print(
            f"[GEO-38 C44]   {title}: {fmt(s15[key])} | "
            f"{fmt(s3[key]) if s3 else pend} | {f'{r[key]:.6e}' if r else pend}"
        )
    print(f"[GEO-38 C44]   c* (hotspot): {_cstar_str(s15)} | {_cstar_str(s3) if s3 else pend}")
    for sense in ("cw", "ccw"):  # cw P33 first: C4's root series
        for p in ("p33", "p41"):
            k = f"{p}_im"
            print(
                f"[GEO-38 C44]   Im Y_in {sense} {p.upper()} [S]: {s15['y_in'][sense][k]} | "
                f"{s3['y_in'][sense][k] if s3 else pend} (side by side, never a ratio)"
            )
    return r


def custody_json(parsed, path):
    """Rank 0 reads the host JSON when present; each C33 field formatted as its log
    token is printed must equal that token's string.  ``(present, n_fields, mismatches)``."""

    def _read():
        p = Path(path)
        return json.loads(p.read_text()) if p.exists() else None

    j = _on_root(_read)
    if j is None:
        print(f"[GEO-38 C44] (J) host JSON absent: {path}")
        return (False, 0, [])
    f9, fs = "{:.9e}".format, "{:+.9e}".format
    pairs = []
    s = parsed["s2"]
    for (i, k), key in zip(((0, 0), (0, 1), (1, 0), (1, 1)), ("s11", "s12", "s21", "s22")):
        re_, im_ = j["s2"][i][k]
        pairs.append((f"s2.{key}", f"{re_:.10e}{im_:+.10e}j", s[key]))
    for sense in ("ccw", "cw"):
        for ip, p in enumerate(("p33", "p41")):
            for ipart, part in enumerate(("re", "im")):
                pairs.append(
                    (
                        f"y_in.{sense}.{p}.{part}",
                        fs(j["y_in"][sense][ip][ipart]),
                        parsed[f"y_in_{sense}"][f"{p}_{part}"],
                    )
                )
    pairs.append(("mean_b_plus.ccw", f9(j["mean_b_plus"]["ccw"]), parsed["means"]["b1p_ccw"]))
    pairs.append(("mean_b_minus.ccw", f9(j["mean_b_minus"]["ccw"]), parsed["means"]["b1m_ccw"]))
    if len(j["lattice"]) != len(parsed["lattice"]):
        pairs.append(("lattice.count", str(len(j["lattice"])), str(len(parsed["lattice"]))))
    for jl, pl in zip(j["lattice"], parsed["lattice"]):
        pairs.append((f"lattice[{pl['label']}].label", jl["label"], pl["label"]))
        pairs.append((f"lattice[{pl['label']}].b_plus_ccw", f9(jl["b_plus_ccw"]), pl["b1p_ccw"]))
        pairs.append((f"lattice[{pl['label']}].sar10g", f9(jl["sar10g"]), pl["sar10g"]))
    cs = parsed["c_star"]
    pairs.append(("star_label", j["star_label"], cs["label"]))
    pairs.append(("sar_star", f9(j["sar_star"]), cs["sar10g"]))
    pairs.append(("sar_star_cw", f9(j["sar_star_cw"]), cs["sar10g_cw"]))
    pairs.append(("c_star.dz", f"{j['c_star'][2]:+.6f}", cs["dz_m"]))
    pairs.append(("c_star.r", f"{math.hypot(j['c_star'][0], j['c_star'][1]):.6f}", cs["r_m"]))
    for k in ("p_box", "p_src", "p_acc"):
        pairs.append((k, f9(j[k]), parsed["normalisers"][k]))
    bad = [(n, a, b) for n, a, b in pairs if a != b]
    print(
        f"[GEO-38 C44] (J) host JSON {Path(path).name}: {len(pairs)} fields, "
        f"{len(bad)} mismatches{'; first ' + repr(bad[0]) if bad else ''}"
    )
    return (True, len(pairs), bad)


def _c33_mode():
    raw = os.environ.get(C33_ENV, "").strip()
    if raw == "":
        return None
    if raw in ("1", "both"):
        return raw
    raise ValueError(f"{C33_ENV} must be unset, '1' or 'both', got {raw!r}")


def _parse(path):
    return _on_root(_load_page().parse_c6_log, path)


def _read_entry52():
    page = _load_page()
    p52 = _parse(page.SOURCE_LOG)
    print(
        f"[GEO-38 C44] (X52) parsed {p52['log']}: node {p52['node']['c_pf']} pF "
        f"(line :{p52['node']['line']}), {len(p52['lattice'])} lattice lines"
    )
    assert p52["node"]["c_pf"] == ENTRY52_NODE_PF, p52["node"]
    require_side(p52, ENTRY52_NODE_PF, "1.5 T (entry 52)")
    return p52


def _c33_body(mode):
    """The gated body; returns the summary dict.  Writes only at ``both``."""
    p52 = _read_entry52()
    present, n, bad = custody_json(p52, HOST_JSON_DIR / f"c6-node-{ENTRY52_NODE_PF}.json")
    assert present, "(J52) host JSON absent on this host"
    assert not bad, f"(J52) {len(bad)} of {n} mismatches, first {bad[0]}"
    hit = _on_root(find_log, ENTRY58_GLOB)
    s15, s3, p58 = side(p52), None, None
    if hit is None:
        print("[GEO-38 C44] entry 58: pending (no *_GEO-38-C6-3T.log tracked or held)")
        if mode == "both":
            raise AssertionError("entry 58 (xl, *_GEO-38-C6-3T.log) absent: C33 needs both")
    else:
        path58, kind = hit
        print(f"[GEO-38 C44] entry 58: {path58.name} ({kind})")
        p58 = _parse(path58)
        require_side(p58, ENTRY58_NODE_PF, "3 T (entry 58)")
        s3 = side(p58)
    r = print_table(s15, s3)
    out = {
        "entry52": s15,
        "entry58": s3,
        "ratios": {k: f"{v:.6e}" for k, v in r.items()} if r else None,
    }
    if mode == "both":

        def _write():
            HOST_JSON_DIR.mkdir(parents=True, exist_ok=True)
            p = HOST_JSON_DIR / C33_JSON_NAME
            p.write_text(json.dumps(out, indent=1))
            return str(p)

        print(f"[GEO-38 C44] wrote {_on_root(_write)}")
    return out


# ---------------------------------------------------------------------------
# Coding tests (no knob).
# ---------------------------------------------------------------------------
def _renamed_copy(tmp_dir, node_to=ENTRY58_NODE_PF, cstar_digit=False):
    """Rank 0: entry 52's log copied as a 3 T name, node token (`:node`) and the
    ``wrote`` path rewritten; optionally one c* SAR digit changed on the c* line."""
    page = _load_page()
    p52 = page.parse_c6_log(page.SOURCE_LOG)
    lines = page._read_log(page.SOURCE_LOG).split("\n")
    i_node, i_wrote, i_cs = (p52[k]["line"] - 1 for k in ("node", "wrote", "c_star"))
    old = p52["node"]["c_pf"]
    assert f"node {old} pF" in lines[i_node] and f"c6-node-{old}.json" in lines[i_wrote]
    lines[i_node] = lines[i_node].replace(f"node {old} pF", f"node {node_to} pF", 1)
    lines[i_wrote] = lines[i_wrote].replace(f"c6-node-{old}.json", f"c6-node-{node_to}.json", 1)
    if cstar_digit:
        tok = p52["c_star"]["sar10g"]
        mant = tok[0] + "." + str((int(tok[2]) + 1) % 10) + tok[3:]  # one digit
        assert f"SAR10g {tok} W/kg (ccw)" in lines[i_cs]
        lines[i_cs] = lines[i_cs].replace(
            f"SAR10g {tok} W/kg (ccw)", f"SAR10g {mant} W/kg (ccw)", 1
        )
    d = Path(tmp_dir)
    d.mkdir(parents=True, exist_ok=True)
    out = d / "20990101T000000Z_GEO-38-C6-3T.log"
    out.write_text("\n".join(lines))
    return str(out)


def _tmp_dir(tmp_path):
    return _on_root(lambda: str(tmp_path))


def test_c44_self_renamed_copy(tmp_path):
    """(Self): the 3 T path on the line forms entry 58 prints."""
    page = _load_page()
    p52 = _parse(page.SOURCE_LOG)
    d = _tmp_dir(tmp_path)
    path = _on_root(_renamed_copy, Path(d) / "self")
    hit = _on_root(find_log, ENTRY58_GLOB, ((Path(d) / "self", "tracked"),))
    assert hit is not None and hit[1] == "tracked" and hit[0].name == Path(path).name, hit
    p3 = _parse(hit[0])
    assert p3["node"]["c_pf"] == ENTRY58_NODE_PF
    require_side(p3, ENTRY58_NODE_PF, "3 T (self copy)")
    s15, s3 = side(p52), side(p3)
    r = print_table(s15, s3)
    assert all(f"{v:.6e}" == "1.000000e+00" for v in r.values()), r
    assert s3["c_star"] == s15["c_star"]
    assert s3["y_in"] == s15["y_in"]


def test_c44_n_digit(tmp_path):
    """(N-digit): one c* SAR digit changed fails (D-sar) and the c*-vs-lattice string."""
    d = _tmp_dir(tmp_path)
    path = _on_root(_renamed_copy, Path(d) / "ndigit", ENTRY58_NODE_PF, True)
    p3 = _parse(path)
    res = check_side(p3, ENTRY58_NODE_PF, "3 T (N-digit)")
    failed = {n for n, ok, _ in res if not ok}
    print(f"[GEO-38 C44] (N-digit) failed: {sorted(failed)}")
    assert "3 T (N-digit) (D-sar)" in failed
    assert "3 T (N-digit) (D-cstar) c* SAR string == its lattice line's string" in failed
    with pytest.raises(AssertionError, match=r"3 T \(N-digit\)"):
        require_side(p3, ENTRY58_NODE_PF, "3 T (N-digit)")


def test_c44_n_node(tmp_path):
    """(N-node): the 3 T name with the node left at entry 52's makes (Node) raise."""
    d = _tmp_dir(tmp_path)
    page = _load_page()
    path = _on_root(
        _renamed_copy, Path(d) / "nnode", page.parse_c6_log(page.SOURCE_LOG)["node"]["c_pf"]
    )
    p3 = _parse(path)
    with pytest.raises(AssertionError, match=r"\(Node\)"):
        require_side(p3, ENTRY58_NODE_PF, "3 T (N-node)")


def test_c44_n_two(tmp_path):
    """(N-two): two files matching the glob raise ``ValueError``."""
    d = Path(_tmp_dir(tmp_path)) / "ntwo"

    def _make():
        d.mkdir(parents=True, exist_ok=True)
        (d / "20990101T000000Z_GEO-38-C6-3T.log").write_text("a\n")
        (d / "20990102T000000Z_GEO-38-C6-3T.log").write_text("b\n")

    _on_root(_make)
    with pytest.raises(ValueError, match="2 matches"):
        _on_root(find_log, ENTRY58_GLOB, ((d, "tracked"),))


def test_c44_n_json(tmp_path):
    """(N-json): ``sar_star`` changed in its 8th significant digit fails (J52)."""
    page = _load_page()
    p52 = _parse(page.SOURCE_LOG)
    src = HOST_JSON_DIR / f"c6-node-{ENTRY52_NODE_PF}.json"
    d = Path(_tmp_dir(tmp_path)) / "njson"

    def _make():
        j = json.loads(src.read_text())
        v = j["sar_star"]
        e = math.floor(math.log10(abs(v)))
        j["sar_star"] = v + 10.0 ** (e - 7)  # one unit in the 8th significant digit
        d.mkdir(parents=True, exist_ok=True)
        (d / src.name).write_text(json.dumps(j))
        return str(d / src.name)

    path = _on_root(_make)
    present, n, bad = custody_json(p52, path)
    assert present and [b[0] for b in bad] == ["sar_star"], bad


def _empty_log_roots(monkeypatch, tmp_path):
    """Point ``find_log``'s default roots at two empty tmp dirs (rank 0 makes them,
    bcast), so a premise test reads entry 58 absent by construction (§9 item 918)."""

    def _make():
        t, h = Path(tmp_path) / "tracked-empty", Path(tmp_path) / "held-empty"
        t.mkdir(parents=True, exist_ok=True)
        h.mkdir(parents=True, exist_ok=True)
        return str(t), str(h)

    t, h = _on_root(_make)
    monkeypatch.setitem(globals(), "TRACKED_LOG_DIR", Path(t))
    monkeypatch.setitem(globals(), "HELD_LOG_DIR", Path(h))


def _json_stamp(p):
    """Rank 0's ``(mtime_ns, bytes)`` of ``p``, or ``None`` if absent: "nothing written"
    holds whether or not C33's own window has already written the file (§9 item 917)."""
    return _on_root(lambda: (p.stat().st_mtime_ns, p.read_bytes()) if p.exists() else None)


def test_c44_n_both(monkeypatch, tmp_path):
    """(N-both): at ``both`` on empty log roots the gated body fails, naming entry 58."""
    _empty_log_roots(monkeypatch, tmp_path)
    before = _json_stamp(HOST_JSON_DIR / C33_JSON_NAME)
    monkeypatch.setenv(C33_ENV, "both")
    mode = _c33_mode()
    assert mode == "both"
    with pytest.raises(AssertionError, match="entry 58"):
        _c33_body(mode)
    assert _json_stamp(HOST_JSON_DIR / C33_JSON_NAME) == before


def test_c44_pend_and_knob(monkeypatch, tmp_path):
    """(Pend) entry 58 absent on empty log roots; the knob's shape."""
    _empty_log_roots(monkeypatch, tmp_path)  # item918-ctl
    assert _on_root(find_log, ENTRY58_GLOB) is None
    monkeypatch.setenv(C33_ENV, "bogus")
    with pytest.raises(ValueError):
        _c33_mode()


def test_c44_real58():
    """(Real58): on the real tree entry 58's log is tracked (landed at ``3b5b7e1f``)."""
    hit = _on_root(find_log, ENTRY58_GLOB)
    print(f"[GEO-38 C44] (Real58) {hit}")
    assert hit is not None and hit[1] == "tracked", hit
    assert hit[0].name == "20261109T080008Z_GEO-38-C6-3T.log", hit


def _parse_entry58():
    hit = _on_root(find_log, ENTRY58_GLOB)
    assert hit is not None and hit[0].name == "20261109T080008Z_GEO-38-C6-3T.log", hit
    p58 = _parse(hit[0])
    assert p58["node"]["c_pf"] == ENTRY58_NODE_PF, p58["node"]
    return p58


def test_c33c_j58():
    """(J58), C33c: entry 58's host JSON formats to its log strings, 76 fields, none off."""
    if _c33_mode() is None:
        pytest.skip(f"{C33_ENV} unset")
    p58 = _parse_entry58()
    present, n, bad = custody_json(p58, HOST_JSON_DIR / f"c6-node-{ENTRY58_NODE_PF}.json")
    print(f"[GEO-38 C33c] (J58) present {present}, {n} fields, mismatches {bad} (ASSERTED)")
    assert present, "(J58) host JSON absent on this host"
    assert n == 76 and bad == [], (n, bad)


def test_c33c_n_sign(tmp_path):
    """(N-sign), C33c: entry 58's JSON with c*'s dz negated fails on c_star.dz alone."""
    if _c33_mode() is None:
        pytest.skip(f"{C33_ENV} unset")
    p58 = _parse_entry58()
    src = HOST_JSON_DIR / f"c6-node-{ENTRY58_NODE_PF}.json"
    d = Path(_tmp_dir(tmp_path)) / "nsign"

    def _make():
        j = json.loads(src.read_text())
        j["c_star"][2] = -j["c_star"][2]
        d.mkdir(parents=True, exist_ok=True)
        (d / src.name).write_text(json.dumps(j))
        return str(d / src.name)

    path = _on_root(_make)
    present, n, bad = custody_json(p58, path)
    print(f"[GEO-38 C33c] (N-sign) present {present}, {n} fields, mismatches {bad} (ASSERTED)")
    assert present and bad == [("c_star.dz", "-0.050000", "+0.050000")], bad


# ---------------------------------------------------------------------------
# The gated test.
# ---------------------------------------------------------------------------
def test_geo38_c33_reader():
    mode = _c33_mode()
    if mode is None:
        pytest.skip(f"{C33_ENV} unset")
    before = _on_root(
        lambda: sorted(p.name for p in HOST_JSON_DIR.glob("*")) if HOST_JSON_DIR.is_dir() else []
    )
    out = _c33_body(mode)
    if mode == "1":
        after = _on_root(
            lambda: (
                sorted(p.name for p in HOST_JSON_DIR.glob("*")) if HOST_JSON_DIR.is_dir() else []
            )
        )
        assert after == before, "wrote under logs/geo38 at knob 1"
        assert out["entry58"] is None or out["ratios"] is not None
        print("[GEO-38 C44] knob 1: nothing written under logs/geo38")
