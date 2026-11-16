"""Example (`EX-195`): the finite-radius two-tip twin at the probe -- I86's closed form beside the
zero-radius (point) twin, the line and I83's cylinder, on the 16 / 17 / 18 / 20 mm ladder.

A closed-form record page: numpy / scipy only, **no mesh and no solve**. `IMP-8` I86 (item 934)
answered, in closed form and infinite gel, whether the rod's radius or the tip weighting dominates
the two-tip twin's probe read. This page evaluates I86's twin Phi_2,a (two uniform cylinders of
radius a = ``ROD_RADIUS`` and length h = 2a centred on z = +-L/2, P/2 each; the length is ratified
by the 11-11 weekly) on I80's four rungs, beside I80's zero-radius two-tip twin ``_phi_tips`` (the
plan's "point twin" and "line twin" both mean it; this page calls it the zero-radius (point) twin
and calls ``_phi`` "the line") and I83's cylinder ``_phi_cyl``.

The page's expressions are I86's own, written out operand for operand (disclosed), so the floats
are the same. The comparand strings are **parsed from I86's tracked log** (``I86_LOG``), never
typed: a local ``re.search`` per line on lines carrying ``[imp8-i86]``. Everything else is imported
by name from the gate modules (the `ANS-1` rule); no ``test_*`` function is imported or called.

Asserted (collected first): (X) custody -- the log exists, ``Status: 0``, parse counts 1 / 1 / 1 /
4 / 4 / 4 / 4, rung keys and z; (K) the weekly's four strings at 16 mm; (K-rec) every other logged
string of I86's record block, by string; (Lim) and (Lim-a) at ``B1_BAND``; (C-twin) 1 - cylinder
twin / point twin at 16 mm, continuous 900 s, > ``INTERIOR_BAND``; (C-a10) I86's own negative;
(N-len) the h = a string resolves the ratified length.

Run it::

    ./scripts/run_examples.sh -e imp:29 -n 1 -t 120

Output: none beyond stdout (no mesh, no field). The setup figure is a byte-identical copy of
`imp:28`'s render; this page has no box.
"""

from __future__ import annotations

import re
import sys
import time
from functools import partial
from pathlib import Path

import numpy as np
from mpi4py import MPI

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from tests.mesh.test_birdcage_implant_rod import ROD_RADIUS  # noqa: E402
from tests.validation.test_imp4_pec_cylinder_closed_form import B1_BAND  # noqa: E402
from tests.validation.test_imp8_entry48_lumped import _tips_steady  # noqa: E402
from tests.validation.test_imp8_probe_cylinder_source_closed_form import _phi_cyl  # noqa: E402
from tests.validation.test_imp8_probe_f2182_closed_form import DUTY, _p  # noqa: E402
from tests.validation.test_imp8_probe_ladder_closed_form import LADDER_MM, RUNGS, _lab  # noqa: E402
from tests.validation.test_imp8_probe_twin_cylinder_closed_form import (  # noqa: E402
    A_LIM,
    _twin_a,
    _twin_a_steady,
)
from tests.validation.test_imp8_rod_line_source_closed_form import (  # noqa: E402
    L,
    MM,
    T_PRINT,
    _phi,
    _phi_tips,
)
from tests.validation.test_imp8_rod_line_source_duty_closed_form import SCALE, T_STAR  # noqa: E402
from tests.validation.test_imp8_rod_tips_duty_closed_form import _duhamel  # noqa: E402
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND  # noqa: E402

#: I86's closing log (item 934, real build, -n 2, Status 0), resolved from the repo root.
I86_LOG = "docs/testing/logs/20261110T120247Z_IMP-8-item934-W1.log"

# Local, disclosed line patterns (one re.search per line; never ^-anchored, never line adjacency).
_RE_CTRL = re.compile(r"cylinder twin (\S+) point twin (\S+) K/W")
_RE_CTRL_SS = re.compile(r"\(ctrl\) steady (\S+) vs (\S+):")
_RE_CTRL_D = re.compile(r"\(ctrl\) D=0\.5 (\S+) vs (\S+):")
_RE_RUNG = re.compile(
    r"rung (\d) \(z (\d+) mm.*SCALE\*Φ2a 900 s ([^,]+), D=0\.5 ([^,]+), steady ([^;]+);")
_RE_RF = re.compile(r"radius factor Φ2a/Φ2: 900 s ([^,]+), D=0\.5 ([^,]+), steady (\S+) ")
_RE_TW = re.compile(
    r"tip weighting twin_a/line_a: 900 s ([^,]+), steady (\S+) "
    r"\(point twin/line here ([^;]+); I80 logged ([^)]+)\)")
_RE_HA = re.compile(r"h=a rung (\d) .*Φ2a/Φ2 900 s ([^,]+),")


def _fmt(x) -> str:
    return f"{x:.7g}"


def parse_i86(path: Path) -> dict:
    """Parse I86's tracked log into its comparand strings (UTF-8: Φ, θ)."""
    text = path.read_text(encoding="utf-8")
    out = dict(ctrl=[], ctrl_ss=[], ctrl_d=[], rung=[], rf=[], tw=[], ha=[], status=None)
    for line in text.splitlines():
        if "[imp8-i86]" in line:
            for key, rx in (("ctrl", _RE_CTRL), ("ctrl_ss", _RE_CTRL_SS), ("ctrl_d", _RE_CTRL_D),
                            ("rung", _RE_RUNG), ("rf", _RE_RF), ("tw", _RE_TW), ("ha", _RE_HA)):
                m = rx.search(line)
                if m:
                    out[key].append(m.groups())
    m = re.search(r"## Exit\s*\n- Status: (\d+)", text)
    out["status"] = int(m.group(1)) if m else None
    return out


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()
    _p("=" * 78)
    _p("EX-195 -- the finite-radius two-tip twin at the probe (IMP-8 I86; I80, I83)")
    _p("=" * 78)
    _p(f"[gate] imported: B1_BAND {B1_BAND:.0e}, INTERIOR_BAND {INTERIOR_BAND:.0e}, "
       f"T_PRINT {T_PRINT:g} s, T_STAR {T_STAR:g} s, DUTY {DUTY:g}, "
       f"ROD_RADIUS {ROD_RADIUS / MM:g} mm, L {L / MM:g} mm, A_LIM {A_LIM:.3e} m, "
       f"SCALE {SCALE:.7g} K/W per unit Phi, LADDER_MM {LADDER_MM}")
    _p("[note] closed form only: infinite gel, no box, no FEM, no rod conduction; every digit is "
       "I86's, reproduced by string. 'point twin' = I80's zero-radius two-tip twin _phi_tips.")

    # ---- custody: parse I86's log -------------------------------------------------------
    log_path = _REPO_ROOT / I86_LOG
    log_exists = log_path.is_file()
    rec = parse_i86(log_path) if log_exists else None
    if rec is None:
        _p(f"[X] {I86_LOG} NOT FOUND")
        gates = {"X_log_exists": False}
        _p(f"[gates] {gates}")
        raise AssertionError(gates)
    counts = (len(rec["ctrl"]), len(rec["ctrl_ss"]), len(rec["ctrl_d"]), len(rec["rung"]),
              len(rec["rf"]), len(rec["tw"]), len(rec["ha"]))
    keys = [(int(k), int(z)) for k, z, *_ in rec["rung"]]
    want_keys = [(k + 1, int(round(z / MM))) for k, z in enumerate(RUNGS)]
    ha_keys = [int(k) for k, _ in rec["ha"]]
    _p(f"[X] {I86_LOG}: Status {rec['status']}; parse counts {counts} (want 1/1/1/4/4/4/4); "
       f"rung keys {keys}")
    gate_x = (
        rec["status"] == 0
        and counts == (1, 1, 1, 4, 4, 4, 4)
        and keys == want_keys == [(1, 16), (2, 17), (3, 18), (4, 20)]
        and ha_keys == [1, 2, 3, 4]
    )
    if not gate_x:
        _p(f"[gates] {{'X_custody': False}}")
        raise AssertionError("(X) custody red; refusing to compare against a malformed parse")

    # ---- the closed forms (I86's expressions, operand for operand) ----------------------
    k_ok_all, k_bad = True, []
    n_cmp = 0

    def cmp(tag, got, want):
        nonlocal k_ok_all, n_cmp
        n_cmp += 1
        ok = got == want
        if not ok:
            k_ok_all = False
            k_bad.append((tag, got, want))
        return ok

    rows = []
    for k, z in enumerate(RUNGS):
        cont = SCALE * _twin_a(0.0, z, T_PRINT)
        duty = _duhamel(_twin_a, 0.0, z, T_STAR, DUTY)
        ss = SCALE * _twin_a_steady(z)
        rf = _twin_a(0.0, z, T_PRINT) / _phi_tips(0.0, z, T_PRINT)
        rf_d = duty / _duhamel(_phi_tips, 0.0, z, T_STAR, DUTY)
        rf_ss = _twin_a_steady(z) / _tips_steady(0.0, z)
        tw = _twin_a(0.0, z, T_PRINT) / _phi_cyl(0.0, z, T_PRINT)
        tw_ss = _twin_a(0.0, z, np.inf) / _phi_cyl(0.0, z, np.inf)
        pt_line = _phi_tips(0.0, z, T_PRINT) / _phi(0.0, z, T_PRINT)
        pt_cont = SCALE * _phi_tips(0.0, z, T_PRINT)
        pt_ss = SCALE * _tips_steady(0.0, z)
        pt_d = _duhamel(_phi_tips, 0.0, z, T_STAR, DUTY)
        line = SCALE * _phi(0.0, z, T_PRINT)
        cyl = SCALE * _phi_cyl(0.0, z, T_PRINT)
        rows.append(dict(cont=cont, duty=duty, ss=ss, rf=rf, rf_d=rf_d, rf_ss=rf_ss, tw=tw,
                         tw_ss=tw_ss, pt_line=pt_line, pt_cont=pt_cont, pt_ss=pt_ss, pt_d=pt_d,
                         line=line, cyl=cyl))
        r_log = rec["rung"][k]
        f_log, t_log = rec["rf"][k], rec["tw"][k]
        cmp(f"rung{k + 1}.cont", _fmt(cont), r_log[2].strip())
        cmp(f"rung{k + 1}.duty", _fmt(duty), r_log[3].strip())
        cmp(f"rung{k + 1}.ss", _fmt(ss), r_log[4].strip())
        cmp(f"rung{k + 1}.rf", _fmt(rf), f_log[0].strip())
        cmp(f"rung{k + 1}.rf_d", _fmt(rf_d), f_log[1].strip())
        cmp(f"rung{k + 1}.rf_ss", _fmt(rf_ss), f_log[2].strip())
        cmp(f"rung{k + 1}.tw", _fmt(tw), t_log[0].strip())
        cmp(f"rung{k + 1}.tw_ss", _fmt(tw_ss), t_log[1].strip())
        cmp(f"rung{k + 1}.pt_line", _fmt(pt_line), t_log[2].strip())
        cmp(f"rung{k + 1}.pt_line_I80", _fmt(pt_line), t_log[3].strip())

    # control strings (rung 1)
    c_cyl_log, c_pt_log = rec["ctrl"][0]
    ss_cyl_log, ss_pt_log = rec["ctrl_ss"][0]
    d_cyl_log, d_pt_log = rec["ctrl_d"][0]
    r0 = rows[0]
    cmp("ctrl.cyl", _fmt(r0["cont"]), c_cyl_log)
    cmp("ctrl.pt", _fmt(r0["pt_cont"]), c_pt_log)
    cmp("ctrl.ss_cyl", _fmt(r0["ss"]), ss_cyl_log)
    cmp("ctrl.ss_pt", _fmt(r0["pt_ss"]), ss_pt_log)
    cmp("ctrl.d_cyl", _fmt(r0["duty"]), d_cyl_log)
    cmp("ctrl.d_pt", _fmt(r0["pt_d"]), d_pt_log)

    # (K), the weekly's four strings at 16 mm, each against the log
    k_four = {
        "cont 60.34935-form": _fmt(r0["cont"]) == c_cyl_log == rec["rung"][0][2].strip(),
        "point twin": _fmt(r0["pt_cont"]) == c_pt_log,
        "rf": _fmt(r0["rf"]) == rec["rf"][0][0].strip(),
        "tw": _fmt(r0["tw"]) == rec["tw"][0][0].strip(),
    }
    _p(f"[K] the weekly's four strings at 16 mm, page vs I86 log: cylinder twin "
       f"{_fmt(r0['cont'])} vs {c_cyl_log}; point twin {_fmt(r0['pt_cont'])} vs {c_pt_log}; "
       f"rf {_fmt(r0['rf'])} vs {rec['rf'][0][0].strip()}; tw {_fmt(r0['tw'])} vs "
       f"{rec['tw'][0][0].strip()}; all equal: {all(k_four.values())}")
    _p(f"[K-rec] {n_cmp} string comparisons of I86's record block and control lines; "
       f"all equal: {k_ok_all}; mismatches {k_bad}")

    # ---- (Lim), (Lim-a) -------------------------------------------------------------------
    lim_worst, lima_worst = 0.0, 0.0
    small = partial(_twin_a, a=A_LIM)
    for k, z0 in enumerate(RUNGS):
        for sg in (1.0, -1.0):
            z = sg * z0
            lim_worst = max(lim_worst, abs(_twin_a(0.0, z, np.inf) / _twin_a_steady(z) - 1.0))
            for got, ref in (
                (small(0.0, z, T_PRINT), _phi_tips(0.0, z, T_PRINT)),
                (small(0.0, z, np.inf), _phi_tips(0.0, z, np.inf)),
                (_duhamel(small, 0.0, z, T_STAR, DUTY), _duhamel(_phi_tips, 0.0, z, T_STAR, DUTY)),
            ):
                lima_worst = max(lima_worst, abs(got / ref - 1.0))
    _p(f"[Lim] |Phi2a(0,z,inf) / g-form - 1|, every rung and sign: worst {lim_worst:.3e} "
       f"(I86 logged 4.441e-16; B1_BAND {B1_BAND:.0e})")
    _p(f"[Lim-a] a = A_LIM, 900 s / inf / D=0.5 vs the point twin, every rung and sign: worst "
       f"{lima_worst:.3e} (I86 logged 4.441e-15; B1_BAND {B1_BAND:.0e})")

    # ---- controls ------------------------------------------------------------------------
    c_twin = 1.0 - r0["cont"] / r0["pt_cont"]
    c_twin_ss = 1.0 - r0["ss"] / r0["pt_ss"]
    c_twin_d = 1.0 - r0["duty"] / r0["pt_d"]
    _p(f"[C-twin] 16 mm, continuous theta(900 s): cylinder twin {_fmt(r0['cont'])} vs point twin "
       f"{_fmt(r0['pt_cont'])} K/W; 1 - ratio {c_twin:.4e} (asserted > INTERIOR_BAND "
       f"{INTERIOR_BAND:g}, < 1; I86 logged 4.4198e-02)")
    _p(f"[C-twin] steady {_fmt(r0['ss'])} vs {_fmt(r0['pt_ss'])}: 1 - ratio {c_twin_ss:.4e} "
       f"(printed; I86 logged 4.0796e-02)")
    _p(f"[C-twin] D=0.5 {_fmt(r0['duty'])} vs {_fmt(r0['pt_d'])}: 1 - ratio {c_twin_d:.4e} "
       f"(below the band; not the weekly's control; I86 logged 3.9031e-03)")
    a10 = abs(_twin_a(0.0, RUNGS[0], T_PRINT, a=ROD_RADIUS / 10.0)
              / _phi_tips(0.0, RUNGS[0], T_PRINT) - 1.0)
    _p(f"[C-a10] I86's own negative: |Phi2a(a/10) / Phi2 - 1| at 16 mm, 900 s {a10:.4e} "
       f"(asserted > B1_BAND {B1_BAND:.0e}; I86 logged 8.571e-04)")

    # (N-len): h = a is a sensitivity; the h = 2a string resolves the ratified length
    ha_page = []
    for z in RUNGS:
        ha_page.append(_fmt(_twin_a(0.0, z, T_PRINT, length=ROD_RADIUS) / _phi_tips(0.0, z, T_PRINT)))
    ha_log = [h[1].strip() for h in rec["ha"]]
    cont_ha = _fmt(SCALE * _twin_a(0.0, RUNGS[0], T_PRINT, length=ROD_RADIUS))
    n_len_ok = ha_page[0] == ha_log[0] and cont_ha != c_cyl_log
    _p(f"[N-len] h = a at 16 mm, 900 s: ratio to the point twin {ha_page[0]} vs I86 log "
       f"{ha_log[0]}; SCALE*twin(h=a) {cont_ha} != {c_cyl_log} (the h = 2a string): {n_len_ok}")

    # ---- the record (printed only) --------------------------------------------------------
    _p("[record] placement and time quoted from memory; W3 confirms or corrects. K per watt per "
       "rod, infinite gel, tips h = 2a centred on +-L/2.")
    for k, r in enumerate(rows):
        _p(f"[record] {_lab(k)} (+{LADDER_MM[k]:g} mm): twin_a 900 s {_fmt(r['cont'])}, D=0.5 "
           f"{_fmt(r['duty'])}, steady {_fmt(r['ss'])}; point twin 900 s {_fmt(r['pt_cont'])}, "
           f"D=0.5 {_fmt(r['pt_d'])}, steady {_fmt(r['pt_ss'])}; line SCALE*Phi {_fmt(r['line'])}; "
           f"I83 cylinder SCALE*Phi_a {_fmt(r['cyl'])}")
        _p(f"[record]   radius factor Phi2a/Phi2: 900 s {_fmt(r['rf'])}, D=0.5 {_fmt(r['rf_d'])}, "
           f"steady {_fmt(r['rf_ss'])}; tip weighting twin_a/line_a: 900 s {_fmt(r['tw'])}, steady "
           f"{_fmt(r['tw_ss'])}; point twin/line {_fmt(r['pt_line'])}")
    rf, tw = r0["rf"], r0["tw"]
    dom = "the tip weighting" if abs(np.log(tw)) > abs(np.log(rf)) else "the radius"
    _p(f"[answer] 16 mm, continuous 900 s: radius factor Phi2a/Phi2 {_fmt(rf)} (moves the read "
       f"{100.0 * (1.0 - rf):.3g} %) against tip weighting twin_a/line_a x {_fmt(tw)}: {dom} "
       f"dominates")
    _p("[beyond] rungs beyond the face: from 17 mm on, Phi2a/Phi2 > 1 (the h = 2a cylinder's "
       "quadrupole excess), beside the line's Phi_a/Phi < 1 (I83); a record, no claim:")
    for k, r in enumerate(rows):
        _p(f"[beyond] rung {k + 1} (z {RUNGS[k] / MM:g} mm): twin radius factor 900 s {_fmt(r['rf'])} "
           f"(> 1: {r['rf'] > 1.0})")
    _p("[h=a] h = 2a ratified by the 11-11 weekly; h = a is a sensitivity. Phi2a/Phi2 at 900 s, "
       "page vs I86 log:")
    for k in range(4):
        _p(f"[h=a] rung {k + 1} ({RUNGS[k] / MM:g} mm): {ha_page[k]} vs {ha_log[k]} "
           f"(equal: {ha_page[k] == ha_log[k]})")

    # ---- assertions, collected first ------------------------------------------------------
    gates = {
        "X_custody": gate_x,
        "K_four_strings": all(k_four.values()),
        "K_rec_all_strings": k_ok_all,
        "Lim": lim_worst <= B1_BAND,
        "Lim_a": lima_worst <= B1_BAND,
        "C_twin": INTERIOR_BAND < c_twin < 1.0,
        "C_a10": a10 > B1_BAND,
        "N_len": n_len_ok,
    }
    _p(f"[gates] {gates}")
    _p(f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    assert all(gates.values()), gates
    _p("\nAll imported gates hold.")


if __name__ == "__main__":
    try:
        main()
    except BaseException:  # one rank raising must not leave the other waiting in a collective
        import traceback

        traceback.print_exc()
        sys.stdout.flush()
        sys.stderr.flush()
        MPI.COMM_WORLD.Abort(1)
