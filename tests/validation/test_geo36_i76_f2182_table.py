"""`GEO-36` §10 I79 (§9 item 893): I76's F2182 table reader, with entry 63's
column blank, stub-proved — a reader over tracked logs, no solve, no mesh.

Brief: docs/planning/chunks/GEO-36.md, "§9 item 893 — 2026-11-08".  The clinical
column reads the print form of `test_geo36_rod_f_clinical_heating.py:336–339`;
it is proved here on I32b's stub log (`STUB_LOG`), and entry 63's own log
(`ENTRY63_LOG`) is `None` until I76 sets it.  No glob: a reader whose output
changes when a file appears is not a record.  `logs/geo36/*.json` is host-only
and never read.  No `test_*` name is imported.

Asserted: (K-s) the F2182 pattern and I68′'s `share` / `twin_g` / `rod_g`
patterns each match exactly once on the stub, every token equal to its literal;
(K-b) I75's three logs, every token equal to its literal at its stated count;
(K-d) I75's means and the stub's rod share re-derived at the print's last digit;
(N) I78's twelve records reproduced from its own module; (P) the placement
identity — the clinical probe sits at I78's `Z_C` (17 mm), not its `Z_P`.
Controls: (K-c) a one-digit change to a copy of stub line `:119` fails (K-s);
(N-c) I78's line read at `Z_P` misses the `Z_C` record by more than (N)'s bound.

Printed (a record, no band): the I76 table — I75's bracket rows at 16 mm, I78 at
16 and 17 mm, the clinical column blank, and the stub dry run of that column
labelled "stub, form proof, not a record".

Not claimed: any clinical number (entry 63's), I76 itself, a like-for-like
comparison (the bracket heats the rod body with 1 W, the clinical FEM heats the
gel by SAR with the rod unheated), a band, the F2182 placement or time.
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path

from tests.mesh.test_birdcage_implant_rod import ROD_LENGTH, ROD_RADIUS
from tests.validation.test_geo36_i68_lumped_vs_clinical import PATTERNS, _H, _S, _p
from tests.validation.test_imp8_probe_f2182_closed_form import SHAPES, Z_C, Z_P, _reads
from tests.validation.test_imp8_rod_line_source_closed_form import MM
from tests.validation.test_imp8_rod_probe_f2182_fem import ARMS
from tests.validation.test_imp8_rod_probe_f2182_fem import Z_P as Z_P_BRACKET_MM

LOG_DIR = Path(__file__).resolve().parents[2] / "docs/testing/logs"

#: Entry 63's log (`xl` Friday 11-13).  I76 sets it; never globbed.
#: Set by I76′ (§9 item 987) from `docs/testing/xl-ledger.md:83`'s basename.
ENTRY63_LOG = "20261113T080008Z_GEO-36-step10-I32-prime.log"
ENTRY63_TIPS = ("top", "bottom")
BLANK = "blank: entry 63 not run (xl Friday 11-13; I76 sets ENTRY63_LOG)"

# --- I76′: entry 63's tokens, copied from ENTRY63_LOG (line numbers of the .log). ---
# :12497
ENTRY63_SHARE = ("3.740924e-05", "2.538288262e-03", "2.538193307e-03", "9.495544579e-08")
# :12500
ENTRY63_TWIN_G = ("2.538800824240e-03", "3.938867478111e+02", "1.000000000000e+00")
# :12515
ENTRY63_ROD_GP = "9.997607068692e-01"
# :12502–12503 (gel twin), :12517–12518 (rod): (point, offset mm, t s, theta K, valid)
ENTRY63_PTS = {"top": "0.105, 0.0, 0.017", "bottom": "0.105, 0.0, -0.017"}
ENTRY63_F2182 = {
    ("gel twin", "top"): (ENTRY63_PTS["top"], "1.0", "900", "+9.364225563e-03", "True"),  # :12502
    ("gel twin", "bottom"): (ENTRY63_PTS["bottom"], "1.0", "900", "+9.382406725e-03", "True"),  # :12503
    ("rod", "top"): (ENTRY63_PTS["top"], "1.0", "900", "+1.488186427e-02", "True"),  # :12517
    ("rod", "bottom"): (ENTRY63_PTS["bottom"], "1.0", "900", "+1.501090199e-02", "True"),  # :12518
}
ENTRY63_ROD_LINE = 12517  # (K-c63)'s line
ENTRY63_STATUS = (
    "entry 63 RED as a window (ledger 2026-11-13; :12532): EM gates green; (E) at dt1/10 and "
    "rod (S2-tip) above band as readings (I32c, item 986; the asserted set is re-bought by xl "
    "entry 71); theta / probe records equal to entry 55's at every printed digit"
)
P63_BOUND = 1.0e-6  # m: two half-units of the print's 6-decimal rounding (heating.py:336)

# --- The stub, I32b's dry run of the print form (line numbers of the .log). ---
STUB_LOG = "20261108T181211Z_GEO-36-I32b-L-W0b.log"
STUB_TIPS = ("a", "b")
# :103
STUB_SHARE = ("1.041667e-02", "4.500000000e-03", "4.453125000e-03", "4.687500000e-05")
# :106
STUB_TWIN_G = ("1.125000000000e-03", "8.888888888889e+02", "1.000000000000e+00")
# :117
STUB_ROD_GP = "3.958333333333e+00"
# :108–109 (gel twin), :119–120 (rod): (point, offset mm, t s, theta K, valid)
STUB_PTS = {"a": "0.03, -0.03, 0.002", "b": "-0.03, 0.03, 0.002"}
STUB_F2182 = {
    ("gel twin", t): (STUB_PTS[t], "1.0", "900", "+2.121079410e-01", "True") for t in STUB_TIPS
} | {
    ("rod", t): (STUB_PTS[t], "1.0", "900", "+8.484317694e-01", "True") for t in STUB_TIPS
}
STUB_ROD_LINE = 119  # (K-c)'s line

# --- I75 (`:120–122` of each log): +16, -16 (theta, D=0.5), then the means. ---
I75_LOGS = {
    "ti": "20261108T194323Z_IMP-8-item881-W1.log",
    "gel": "20261108T195853Z_IMP-8-item881-W2b.log",  # W2b; the crashed W2 is never read
    "hole": "20261108T200636Z_IMP-8-item881-W3.log",
}
I75 = {
    "ti": ((("+16", "1.064117e+01", "2.109553e+00"), ("-16", "1.049960e+01", "2.099602e+00")),
           ("1.057039e+01", "2.104578e+00")),
    "gel": ((("+16", "9.220827e+00", "1.979610e+00"), ("-16", "9.102961e+00", "1.970586e+00")),
            ("9.161894e+00", "1.975098e+00")),
    "hole": ((("+16", "8.815610e+00", "1.895292e+00"), ("-16", "8.750725e+00", "1.888582e+00")),
             ("8.783168e+00", "1.891937e+00")),
}

# --- I78 (`20261108T133157Z_IMP-8-item874-W1.log:58–59`, `:61–62`), .7g: cont, D=0.5, ss. ---
I78 = {
    ("line", "Z_P"): ("9.648662", "1.978384", "15.18156"),
    ("line", "Z_C"): ("6.810428", "1.822237", "12.25753"),
    ("twin", "Z_P"): ("63.14", "2.498341", "68.45374"),
    ("twin", "Z_C"): ("29.97039", "2.451249", "35.22961"),
}
N_BOUND = 1.0e-6  # one unit of the 7th significant digit of 6.810428
P_BOUND = 1.0e-15  # m

F2182 = (_H + r"\[{arm}\] F2182-shaped probe {tip} \(([^)]*)\) \((\S+) mm beyond the apex\), "
         + _S + r": theta\((\S+) s\) (\S+) K \(valid (\S+);")
I75_PROBE = (r"^\[imp8-i75\] probe \(0, ([+-]16)\) mm arm {arm}: theta\(900 s\) (\S+) K/W; "
             r".*D=0\.5 superposed theta\(900 s\) (\S+) K/W$")
I75_MEAN = r"^\[imp8-i75\] probe arm {arm}: mean (\S+) K/W, .*D=0\.5 mean (\S+)$"


def _read_log(basename):
    """A tracked log, or its `.log.gz` when housekeeping has rotated it (never gunzipped)."""
    path = LOG_DIR / basename
    if path.exists():
        return path.read_text(encoding="utf-8", errors="replace")
    with gzip.open(path.with_suffix(".log.gz"), "rt", encoding="utf-8", errors="replace") as fh:
        return fh.read()


def _once(pat, text, key, failures):
    hits = re.findall(pat, text, flags=re.MULTILINE)
    if len(hits) != 1:
        failures.append(f"{key}: {len(hits)} matches")
        return None
    return hits[0] if isinstance(hits[0], tuple) else (hits[0],)


def _clinical_rows(text, tips):
    """Per (arm, tip): (point, offset, t, theta, valid); and share, twin G, rod G*P tokens."""
    failures, rows = [], {}
    for arm in ("gel twin", "rod"):
        for tip in tips:
            hit = _once(F2182.format(arm=re.escape(arm), tip=re.escape(tip)), text, f"{arm}/{tip}", failures)
            if hit is not None:
                rows[(arm, tip)] = hit
    share = _once(PATTERNS["share"][0], text, "share", failures)
    twin_g = _once(PATTERNS["twin_g"][0], text, "twin_g", failures)
    rod_g = _once(PATTERNS["rod_g"][0], text, "rod_g", failures)
    return rows, share, twin_g, rod_g, failures


def _stub_failures(text):
    rows, share, twin_g, rod_g, failures = _clinical_rows(text, STUB_TIPS)
    for key, want in STUB_F2182.items():
        if key in rows and rows[key] != want:
            failures.append(f"{key}: log {rows[key]!r} != literal {want!r}")
    for name, got, want in (("share", share, STUB_SHARE), ("twin_g", twin_g, STUB_TWIN_G),
                            ("rod_g", rod_g, (STUB_ROD_GP,))):
        if got is not None and tuple(got) != want:
            failures.append(f"{name}: log {got!r} != literal {want!r}")
    return failures


def _i75_failures(arm, text):
    failures = []
    probes = re.findall(I75_PROBE.format(arm=arm), text, flags=re.MULTILINE)
    if len(probes) != 2:
        failures.append(f"{arm} probe: {len(probes)} matches")
    elif tuple(probes) != I75[arm][0]:
        failures.append(f"{arm} probe: log {probes!r} != literal {I75[arm][0]!r}")
    mean = _once(I75_MEAN.format(arm=arm), text, f"{arm} mean", failures)
    if mean is not None and tuple(mean) != I75[arm][1]:
        failures.append(f"{arm} mean: log {mean!r} != literal {I75[arm][1]!r}")
    return failures


def test_ks_stub_reads_the_entry63_print_form():
    failures = _stub_failures(_read_log(STUB_LOG))
    _p(f"\n[GEO-36 I79] (K-s) stub {STUB_LOG}: 4 F2182 rows + share/twin_g/rod_g, "
       f"failures {failures or 'none'} (ASSERTED exact)")
    assert not failures


def test_kc_one_digit_changed_fails():
    lines = _read_log(STUB_LOG).splitlines()
    orig = lines[STUB_ROD_LINE - 1]
    assert orig.startswith("[GEO-36 I32] [rod] F2182-shaped probe a ")
    lines[STUB_ROD_LINE - 1] = orig.replace("+8.484317694e-01", "+8.484317695e-01")
    assert lines[STUB_ROD_LINE - 1] != orig
    failures = _stub_failures("\n".join(lines))
    _p(f"[GEO-36 I79] (K-c) :{STUB_ROD_LINE} theta last digit 4->5: failures {failures} (ASSERTED non-empty)")
    assert failures


def test_kb_i75_logs_read_as_recorded():
    allf = []
    for arm in ARMS:
        f = _i75_failures(arm, _read_log(I75_LOGS[arm]))
        _p(f"[GEO-36 I79] (K-b) {arm} {I75_LOGS[arm]}: failures {f or 'none'} (ASSERTED exact)")
        allf += f
    assert not allf


def test_kd_rederived_at_last_digit():
    worst = []
    for arm in ARMS:
        (p, m), means = I75[arm]
        for k, label in ((1, "theta"), (2, "D=0.5")):
            avg = 0.5 * (float(p[k]) + float(m[k]))
            s = means[k - 1]
            unit = 10.0 ** (int(s.split("e")[1]) - 6)
            d = abs(avg - float(s))
            _p(f"[GEO-36 I79] (K-d) {arm} {label} mean {s} vs (a+b)/2: |diff| {d:.2e} (bound {unit:.0e}, ASSERTED)")
            worst.append(d <= unit)
    share, p_vol, _g, p_rod = STUB_SHARE
    d = abs(float(p_rod) / float(p_vol) - float(share))
    _p(f"[GEO-36 I79] (K-d) stub P_rod/P_vol vs {share}: |diff| {d:.2e} (bound 5e-09, ASSERTED)")
    assert all(worst)
    assert d <= 5.0e-9


def test_n_i78_reproduced():
    bad = []
    for name, (kern, _kd, steady) in SHAPES.items():
        for lab, z in (("Z_P", Z_P), ("Z_C", Z_C)):
            got = tuple(f"{v:.7g}" for v in _reads(kern, steady, z))
            _p(f"[GEO-36 I79] (N) {name} {lab}: {got} vs {I78[(name, lab)]} (ASSERTED string)")
            if got != I78[(name, lab)]:
                bad.append((name, lab, got))
    assert not bad


def test_nc_line_at_zp_in_clinical_row_misses():
    kern, _kd, steady = SHAPES["line"]
    cont = _reads(kern, steady, Z_P)[0]
    miss = abs(cont - float(I78[("line", "Z_C")][0]))
    _p(f"[GEO-36 I79] (N-c) line read at Z_P {cont:.7g} in the 17 mm row vs Z_C record "
       f"{I78[('line', 'Z_C')][0]}: |miss| {miss:.4g} K/W (> {N_BOUND:.0e}, ASSERTED)")
    assert miss > N_BOUND


def test_p_placement_identity():
    rows, *_rest, failures = _clinical_rows(_read_log(STUB_LOG), STUB_TIPS)
    assert not failures
    offsets = {r[1] for r in rows.values()}
    assert len(offsets) == 1
    off = float(offsets.pop())
    z_clin = 0.5 * ROD_LENGTH + ROD_RADIUS + off * MM
    d_c = abs(z_clin - Z_C)
    d_p = abs(Z_P_BRACKET_MM * MM - Z_P)
    _p(f"[GEO-36 I79] (P) clinical probe L/2 + a + {off:g} mm = {z_clin / MM:.6g} mm vs I78 Z_C "
       f"{Z_C / MM:.6g} mm: |diff| {d_c:.1e} m; I75 Z_P {Z_P_BRACKET_MM:g} mm vs I78 Z_P "
       f"{Z_P / MM:.6g} mm: |diff| {d_p:.1e} m (bound {P_BOUND:.0e}, ASSERTED)")
    assert d_c <= P_BOUND
    assert d_p <= P_BOUND


def _clinical_cells(text, tips, label):
    rows, share, twin_g, rod_g, failures = _clinical_rows(text, tips)
    assert not failures, failures
    g = float(twin_g[1])
    p_rod = share[3]
    gp_rod = g * float(p_rod)
    _p(f"[GEO-36 I79] {label}: G {twin_g[1]} 1/W, P_rod {p_rod} W, G*P_rod {gp_rod:.7g} W in the rod "
       f"per W in the gel (twin G, as I68' converts)")
    for arm in ("rod", "gel twin"):
        for tip in tips:
            pt, off, t, th, valid = rows[(arm, tip)]
            _p(f"[GEO-36 I79] {label} [{arm}] {tip} ({pt}), {off} mm beyond the capsule apex = I78 Z_C: "
               f"theta({t} s) {th} K (valid {valid}); per rod watt {float(th) / gp_rod:.7g} K/W "
               f"(gel SAR, rod unheated, ASTM box; not a matched source)")


def test_printed_table():
    _p("[GEO-36 I79] the I76 table — theta(900 s), K per W deposited in the rod; a record, no band")
    for arm in ARMS:
        (p, m), means = I75[arm]
        ln = I78[("line", "Z_P")]
        _p(f"[GEO-36 I79] I75 {arm} (16 mm, 1 mm beyond a flat end; 1 W in the rod, self-heated; "
           f"{I75_LOGS[arm]}:120–122): theta mean {means[0]}, D=0.5 mean {means[1]}; "
           f"FEM/line {float(means[0]) / float(ln[0]):.7g}, D=0.5 {float(means[1]) / float(ln[1]):.7g}")
    for name in ("line", "twin"):
        for lab, mm in (("Z_P", "16"), ("Z_C", "17")):
            c, d, s = I78[(name, lab)]
            _p(f"[GEO-36 I79] I78 {name} at {lab} = {mm} mm (infinite gel, 1 W per rod): theta {c}, "
               f"D=0.5 {d}, ss {s}")
    if ENTRY63_LOG is None:
        _p(f"[GEO-36 I79] clinical rod / gel twin at 17 mm (= I78 Z_C), per rod watt: {BLANK}")
    else:
        _clinical_cells(_read_log(ENTRY63_LOG), ENTRY63_TIPS, f"entry 63 {ENTRY63_LOG}")
    _clinical_cells(_read_log(STUB_LOG), STUB_TIPS, "stub, form proof, not a record")


# --- I76′ (§9 item 987): entry 63's column, read through the reader above. ---


def _entry63_failures(text):
    rows, share, twin_g, rod_g, failures = _clinical_rows(text, ENTRY63_TIPS)
    for key, want in ENTRY63_F2182.items():
        if key in rows and rows[key] != want:
            failures.append(f"{key}: log {rows[key]!r} != literal {want!r}")
    for name, got, want in (("share", share, ENTRY63_SHARE), ("twin_g", twin_g, ENTRY63_TWIN_G),
                            ("rod_g", rod_g, (ENTRY63_ROD_GP,))):
        if got is not None and tuple(got) != want:
            failures.append(f"{name}: log {got!r} != literal {want!r}")
    return failures


def test_i76_k63_entry63_reads_as_recorded():
    failures = _entry63_failures(_read_log(ENTRY63_LOG))
    _p(f"\n[GEO-36 I76'] (K-63) {ENTRY63_LOG}: 4 F2182 rows + share/twin_g/rod_g, "
       f"failures {failures or 'none'} (ASSERTED exact)")
    assert not failures


def test_i76_kc63_one_digit_changed_fails():
    lines = _read_log(ENTRY63_LOG).splitlines()
    orig = lines[ENTRY63_ROD_LINE - 1]
    assert orig.startswith("[GEO-36 I32] [rod] F2182-shaped probe top ")
    lines[ENTRY63_ROD_LINE - 1] = orig.replace("+1.488186427e-02", "+1.488186428e-02")
    assert lines[ENTRY63_ROD_LINE - 1] != orig
    failures = _entry63_failures("\n".join(lines))
    _p(f"[GEO-36 I76'] (K-c63) :{ENTRY63_ROD_LINE} theta last digit 7->8: failures {failures} "
       f"(ASSERTED non-empty)")
    assert failures


def test_i76_p63_placement():
    import tests.validation.test_imp8_probe_cylinder_source_closed_form as i83

    rows, *_rest, failures = _clinical_rows(_read_log(ENTRY63_LOG), ENTRY63_TIPS)
    assert not failures, failures
    bad = []
    for arm in ("gel twin", "rod"):
        top = rows[(arm, "top")][0].split(", ")
        bot = rows[(arm, "bottom")][0].split(", ")
        half = 0.5 * (float(top[2]) - float(bot[2]))
        d = abs(half - Z_C)
        offs = sorted({rows[(arm, t)][1] for t in ENTRY63_TIPS})
        _p(f"[GEO-36 I76'] (P63) {arm}: top ({rows[(arm, 'top')][0]}), bottom ({rows[(arm, 'bottom')][0]}); "
           f"(z_top - z_bottom)/2 {half / MM:.6g} mm vs I78 Z_C {Z_C / MM:.6g} mm: |diff| {d:.1e} m "
           f"(bound {P63_BOUND:.0e}); offsets {offs} (ASSERTED)")
        if top[:2] != bot[:2] or d > P63_BOUND or offs != ["1.0"]:
            bad.append(arm)
    d_l = abs(i83.RUNGS[1] - Z_C)
    _p(f"[GEO-36 I76'] (L17) I81 RUNGS[1] {i83.RUNGS[1] / MM:.6g} mm vs Z_C {Z_C / MM:.6g} mm: "
       f"|diff| {d_l:.1e} m (bound {P_BOUND:.0e}, ASSERTED)")
    assert not bad, bad
    assert d_l <= P_BOUND


def test_i76_readout_beside_ladder_and_phi_a():
    import tests.validation.test_imp8_probe_cylinder_source_closed_form as i83

    _p(f"[GEO-36 I76'] {ENTRY63_STATUS}")
    rows, share, twin_g, _rod_g, failures = _clinical_rows(_read_log(ENTRY63_LOG), ENTRY63_TIPS)
    assert not failures, failures
    gp_rod = float(twin_g[1]) * float(share[3])
    per_w = {key: float(r[3]) / gp_rod for key, r in rows.items()}
    for (arm, tip), v in per_w.items():
        _p(f"[GEO-36 I76'] clinical [{arm}] {tip} at 17 mm (= I78 Z_C): theta(900 s) {rows[(arm, tip)][3]} K; "
           f"per rod watt {v:.7g} K/W (G*P_rod {gp_rod:.7g} W; gel SAR, rod unheated, ASTM box; "
           f"not a matched source)")
    for tip in ENTRY63_TIPS:
        tip_ratio = float(rows[("rod", tip)][3]) / float(rows[("gel twin", tip)][3])
        _p(f"[GEO-36 I76'] clinical rod / gel twin F2182 theta ratio {tip}: {tip_ratio:.7g} (printed)")
    for arm in ARMS:
        _p(f"[GEO-36 I76'] I75 {arm} 16 mm (1 W in the rod body, flat end): theta mean {I75[arm][1][0]}, "
           f"D=0.5 mean {I75[arm][1][1]}")
    for arm, rec in i83.FEM.items():
        for k, z in enumerate(i83.RUNGS):
            _p(f"[GEO-36 I76'] I81 {arm} |z| {z / MM:.6g} mm (1 W in the rod, I55's cube; {rec['log']}): "
               f"+/- mean {rec['mean'][k]}, D=0.5 mean {rec['dmean'][k]} K/W")
    closed = [i83._closed(z) for z in i83.RUNGS]
    for z, (phi, dphi) in zip(i83.RUNGS, closed):
        _p(f"[GEO-36 I76'] I83 |z| {z / MM:.6g} mm (uniform cylinder, infinite gel, 1 W): SCALE*Phi_a {phi:.7g}, "
           f"D=0.5 {dphi:.7g} K/W")
    c, d, s = I78[("line", "Z_C")]
    _p(f"[GEO-36 I76'] I78 line at Z_C = 17 mm (infinite gel, 1 W): theta {c}, D=0.5 {d}, ss {s}")
    rod_mean = 0.5 * (per_w[("rod", "top")] + per_w[("rod", "bottom")])
    for arm, rec in i83.FEM.items():
        m17 = float(rec["mean"][1])
        _p(f"[GEO-36 I76'] 17 mm ratio: clinical rod per rod watt (tip mean) {rod_mean:.7g} / I81 {arm} "
           f"{rec['mean'][1]} = {rod_mean / m17:.7g} (printed, a record, no band)")
    _p(f"[GEO-36 I76'] 17 mm ratio: clinical rod per rod watt (tip mean) {rod_mean:.7g} / SCALE*Phi_a "
       f"{closed[1][0]:.7g} = {rod_mean / closed[1][0]:.7g} (printed, a record, no band)")
    # (N-row): I81's gel 16 mm mean placed in the 17 mm cell misses the 17 mm mean.
    m16s, m17s = i83.FEM["gel"]["mean"][0], i83.FEM["gel"]["mean"][1]
    unit = 10.0 ** (int(m17s.split("e")[1]) - 6)
    miss = abs(float(m16s) - float(m17s))
    _p(f"[GEO-36 I76'] (N-row) I81 gel 16 mm mean {m16s} in the 17 mm cell vs {m17s}: |miss| {miss:.7g} K/W "
       f"(> {unit:.0e}, ASSERTED); separation 1 - m17/m16 {1.0 - float(m17s) / float(m16s):.6f}")
    assert miss > unit
