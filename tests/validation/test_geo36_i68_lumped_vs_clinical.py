"""`GEO-36` §10 I68′ (§9 item 872): the lumped reading against the clinical FEM
on stated normalisations — a reader over a tracked log, no solve, no mesh.

Brief: docs/planning/chunks/GEO-36.md, "§9 item 872 — 2026-11-08".  The FEM
side is entry 55's tracked log (`ENTRY55_LOG`); there is no saved JSON for it
(`_run_heat` writes none), and `logs/geo36/*.json` is host-only, so it is never
read.  The lumped side is I58's module (`test_imp8_entry48_lumped`), imported by
name; no `test_*` name is imported here.

Asserted: (K) every parsed token string equals its transcribed literal, each
pattern matching exactly once; (K-d) the log's own ratios re-derived at the
print's last digit, and MAT-21's t* reproduced as the log's string; (N) I58's
records reproduced from its own module, with one ball radius on both sides.
Controls: (N-c) I58's line kernel read at entry 55's t* instead of 900 s misses
the 900 s record by more than (N)'s bound; (K-c) a one-digit change to a copy of
line `:12512` fails (K).

Printed (a record, no band): θ_FEM(t*) per watt in the gel, per watt in the rod
(θ_FEM / (G·P_rod)) and per watt of P_vol (θ_FEM / (G·P_vol)); "per accepted
watt: not in a tracked record"; I58 per rod watt at 900 s and at t*; and FEM /
lumped per rod watt, labelled with what differs.

Not claimed: a like-for-like comparison (the FEM has no rod source, the lumped
model no box and no gel SAR), a FEM read at 900 s, P_acc, any tissue or
published comparand.
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path

from mpi4py import MPI

from tests.validation.test_imp8_entry48_lumped import (
    N_BQ,
    P_ROD_LOG,
    PER_WATT,
    _a_1g,
    _ball_average_axis,
    _line,
    _tips,
)
from tests.validation.test_imp8_rod_line_source_closed_form import L, T_PRINT
from tests.validation.test_mat21_box_fem import T_STAR as T_STAR_BOX

ENTRY55_LOG = "20261106T080008Z_GEO-36-step10-I32.log"
LOG_DIR = Path(__file__).resolve().parents[2] / "docs/testing/logs"

# --- The literals, transcribed from ENTRY55_LOG (line numbers of the .log). ---
# :12497
ROD_SHARE = "3.740924e-05"
P_VOL = "2.538288262e-03"
P_GEL_SPLIT = "2.538193307e-03"
P_ROD = "9.495544579e-08"
# :12499
TWIN_P_GEL = "2.538800824240e-03"
G_TWIN = "3.938867478110e+02"
GP_TWIN = "1.000000000000e+00"
# :12501–12502
TWIN_TOP = "3.082473767e-02"
TWIN_BOTTOM = "3.090736898e-02"
BALL_MM = "6.2035"
# :12509
T_STAR_LOG = "2821.500"
DT_LOG = "70.537500"
# :12510
GP_ROD = "9.997607068690e-01"
# :12512–12513
ROD_TOP = "3.356701564e-02"
ROD_BOTTOM = "3.368366738e-02"
# :12521
RATIO_TOP = "1.088964"
RATIO_BOTTOM = "1.089826"

# I58's records (`20261103T211525Z_IMP-8-item786-W1.log:63–66`), per rod watt at
# T_PRINT = 900 s: the 1 g ball average and the rod-surface point.
I58_BALL_LINE = "7.307753"
I58_BALL_TIPS = "12.93694"
I58_SURF_LINE = "12.48967"
I58_SURF_TIPS = "63.16158"
N_BOUND = 1.0e-6  # one unit of the 7th significant digit of 7.307753

_H = r"^\[GEO-36 I32\] "
_S = r"sense \(1,\+1j\)"
_BALL = _H + r"\[{arm}\] 1 g ball \(a (\S+) mm\) at {tip}, " + _S + r": theta\(t\*\) (\S+) K,"
PATTERNS = {
    "share": (_H + r"rod share of P_vol, " + _S
              + r": (\S+) \(P_vol (\S+) W, gel (\S+), rod (\S+); printed\)$",
              ("ROD_SHARE", "P_VOL", "P_GEL_SPLIT", "P_ROD")),
    "twin_g": (_H + r"\[gel twin\] " + _S + r": P_gel (\S+) W; G (\S+) 1/W; G\*P (\S+) W;",
               ("TWIN_P_GEL", "G_TWIN", "GP_TWIN")),
    "twin_top": (_BALL.format(arm="gel twin", tip="top"), ("BALL_MM", "TWIN_TOP")),
    "twin_bottom": (_BALL.format(arm="gel twin", tip="bottom"), ("BALL_MM", "TWIN_BOTTOM")),
    "twin_cn": (_H + r"\[gel twin\] CN N=40 to t\* (\S+) s, dt (\S+) s ", ("T_STAR_LOG", "DT_LOG")),
    "rod_g": (_H + r"\[rod\] " + _S + r": P_gel \S+ W; G \S+ 1/W; G\*P (\S+) W;", ("GP_ROD",)),
    "rod_top": (_BALL.format(arm="rod", tip="top"), ("BALL_MM", "ROD_TOP")),
    "rod_bottom": (_BALL.format(arm="rod", tip="bottom"), ("BALL_MM", "ROD_BOTTOM")),
    "ratio": (_H + r"tip theta\(t\*\) ratio rod / gel twin, " + _S + r": top (\S+), bottom (\S+) ",
              ("RATIO_TOP", "RATIO_BOTTOM")),
}


def _p(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _log_text():
    path = LOG_DIR / ENTRY55_LOG
    if path.exists():
        return path.read_text(encoding="utf-8", errors="replace")
    with gzip.open(path.with_suffix(".log.gz"), "rt", encoding="utf-8", errors="replace") as fh:
        return fh.read()


def _parse(text):
    """Token strings by literal name, and a list of (K) failures."""
    tokens, failures = {}, []
    for key, (pat, names) in PATTERNS.items():
        hits = re.findall(pat, text, flags=re.MULTILINE)
        if len(hits) != 1:
            failures.append(f"{key}: {len(hits)} matches")
            continue
        vals = hits[0] if isinstance(hits[0], tuple) else (hits[0],)
        for name, val in zip(names, vals, strict=True):
            tokens.setdefault(name, []).append(val)
    for name, vals in tokens.items():
        for v in vals:
            if v != globals()[name]:
                failures.append(f"{name}: log {v!r} != literal {globals()[name]!r}")
    return {n: v[0] for n, v in tokens.items()}, failures


def _entry55_record():
    return _parse(_log_text())


def _lumped_ball(f, t, sign):
    return _ball_average_axis(lambda r, z: f(r, z, t), sign * 0.5 * L, _a_1g(), 2 * N_BQ) * PER_WATT


def test_k_record_and_controls():
    rec, failures = _entry55_record()
    _p(f"[GEO-36 I68'] (K) {ENTRY55_LOG}: {len(PATTERNS)} patterns, {len(rec)} tokens; "
       f"failures {failures or 'none'} (ASSERTED, string equality)")
    assert not failures, failures
    # (K-c): a copy of :12512 with one digit changed, in memory.
    text = _log_text()
    bad = text.replace(f"at top, sense (1,+1j): theta(t*) {ROD_TOP} K", "at top, sense (1,+1j): "
                       f"theta(t*) {ROD_TOP[:-5]}7{ROD_TOP[-4:]} K", 1)
    assert bad != text
    _, bad_failures = _parse(bad)
    _p(f"[GEO-36 I68'] (K-c) one digit changed on :12512: failures {bad_failures} (ASSERTED non-empty)")
    assert bad_failures


def test_kd_rederived():
    rec, failures = _entry55_record()
    assert not failures, failures
    for tip, rod, twin, ratio in (("top", ROD_TOP, TWIN_TOP, RATIO_TOP),
                                  ("bottom", ROD_BOTTOM, TWIN_BOTTOM, RATIO_BOTTOM)):
        d = abs(float(rod) / float(twin) - float(ratio))
        _p(f"[GEO-36 I68'] (K-d) {tip}: theta_rod/theta_twin vs {ratio}: |d| {d:.3e} (band 1e-06, ASSERTED)")
        assert d <= 1.0e-6
    d = abs(float(P_ROD) / float(P_VOL) - float(ROD_SHARE))
    _p(f"[GEO-36 I68'] (K-d) rod/P_vol vs {ROD_SHARE}: |d| {d:.3e} (band 5e-12, ASSERTED)")
    assert d <= 5.0e-12
    t_s = f"{T_STAR_BOX:.3f}"
    _p(f"[GEO-36 I68'] (K-d) MAT-21 T_STAR {t_s} s vs log {T_STAR_LOG} s (ASSERTED string equality)")
    assert t_s == T_STAR_LOG


def test_n_lumped_reproduced_and_record():
    a_s = f"{_a_1g() * 1e3:.4f}"
    _p(f"[GEO-36 I68'] (N) I58 ball radius {a_s} mm vs entry 55 {BALL_MM} mm (ASSERTED string equality)")
    assert a_s == BALL_MM
    lumped = {}
    for tip, sign in (("top", 1.0), ("bottom", -1.0)):
        ln = _lumped_ball(_line, T_PRINT, sign)
        tp = _lumped_ball(_tips, T_PRINT, sign)
        _p(f"[GEO-36 I68'] (N) {tip} ball at {T_PRINT:g} s: line {ln:.7g} vs {I58_BALL_LINE}, "
           f"tips {tp:.7g} vs {I58_BALL_TIPS} K/W (ASSERTED string equality)")
        assert f"{ln:.7g}" == I58_BALL_LINE
        assert f"{tp:.7g}" == I58_BALL_TIPS
        ln_t = _lumped_ball(_line, T_STAR_BOX, sign)
        miss = abs(ln_t - float(I58_BALL_LINE))
        _p(f"[GEO-36 I68'] (N-c) {tip} ball line read at t* {T_STAR_BOX:.1f} s: {ln_t:.7g} K/W; "
           f"|miss| vs {I58_BALL_LINE} {miss:.3e} (> {N_BOUND:.0e}, ASSERTED)")
        assert miss > N_BOUND
        lumped[tip] = ln_t

    # The record.  G and the powers are the log's; conversions are new strings.
    g = float(G_TWIN)
    gp_rod, gp_vol = g * float(P_ROD), g * float(P_VOL)
    a = f"ball a {BALL_MM} mm"
    _p(f"[GEO-36 I68'] conversions ({a}): G {G_TWIN} 1/W (:12499); P_rod {P_ROD_LOG} W, P_vol {P_VOL} W "
       f"(:12497); G*P_rod {gp_rod:.7g} W in the rod per W in the gel; G*P_vol {gp_vol:.7g} W")
    _p(f"[GEO-36 I68'] per accepted watt: not in a tracked record ({a})")
    for arm, gp_gel, vals in (("gel twin", GP_TWIN, (("top", TWIN_TOP), ("bottom", TWIN_BOTTOM))),
                              ("rod", GP_ROD, (("top", ROD_TOP), ("bottom", ROD_BOTTOM)))):
        for tip, th_s in vals:
            th = float(th_s)
            per_rod, per_vol = th / gp_rod, th / gp_vol
            _p(f"[GEO-36 I68'] [{arm}] {tip} ({a}): theta_FEM(t*) {th_s} K at G*P {gp_gel} W in the gel "
               f"(the print); per rod watt {per_rod:.7g} K/W; per P_vol watt {per_vol:.7g} K/W "
               f"(nearest tracked normalisation)")
            _p(f"[GEO-36 I68'] [{arm}] {tip} ({a}): lumped per rod watt — 900 s ball line {I58_BALL_LINE}, "
               f"ball tips {I58_BALL_TIPS}, surface line {I58_SURF_LINE}, surface tips {I58_SURF_TIPS}; "
               f"t* ball line {lumped[tip]:.7g} K/W")
            _p(f"[GEO-36 I68'] [{arm}] {tip} ({a}): FEM/lumped per rod watt — vs 900 s ball line "
               f"{per_rod / float(I58_BALL_LINE):.7g}, ball tips {per_rod / float(I58_BALL_TIPS):.7g}, "
               f"surface line {per_rod / float(I58_SURF_LINE):.7g}, surface tips "
               f"{per_rod / float(I58_SURF_TIPS):.7g}; vs t* ball line {per_rod / lumped[tip]:.7g} "
               f"(differs: FEM SAR in the gel, rod unheated (tag 403 outside the thermal domain); lumped "
               f"P_rod in the rod's line and tips; {T_STAR_LOG} s against {T_PRINT:g} s except the t* "
               f"column; walled box against infinite medium)")
