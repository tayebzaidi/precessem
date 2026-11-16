"""Example (`EX-194`): the field-strength lever on the heating -- P83's closed-form prediction
beside C39's FEM reading of `xl` entries 59 / 60, per accepted watt.

A **record page, no solve and no mesh** (the ``mri:41`` precedent).  A user asks how much more
the gel heats at 3 T than at 1.5 T per accepted watt.  ``MAT-21`` step 46 (P83) wrote the
closed-form pictures before entry 60 ran: every picture is mirror-even, so the prediction is the
accepted-power share alone, s3 / s1.5.  ``GEO-38`` C39 then read the two FEM heating windows
(entry 59: 1.5 T, 64 MHz, node 193.5996 pF; entry 60: 3 T, 128 MHz, node 44.0525 pF, degree 1,
each at its own root and its own c*).  This page sets the two on one page and prints their
quotient, the heating's shape part, beside the SAR's.

Everything is read from tracked logs: entries 59 / 60 through the C46 module's ``find_log`` /
``parse_heat_log`` / ``check_block`` / ``side`` / ``ratios`` (module ``c46``); C39's closing
``-n 1`` log and P83's log by local ``re.search`` (disclosed).  No ``test_*`` function, no gated
body, no host JSON, no gate-module edit, no ``src/`` edit.

Asserted (all collected first, then asserted):

* (X)       each ``xl`` log holds one block, node strings equal the node constants, the probe
            names equal ``c46._probes()``; C39's real (non-STUB) table is found once with 22 rows;
            P83's rows are found with counts 1 / 4 / 4 / 3 (+ the C33 ratio tuple once);
* (D)       ``check_block`` green on both sides under both senses at C46's own half-unit bounds;
* (R39)     the anchor: each of the 22 page ratios at ``.6e`` equals C39's closing-log string
            for the same sense and readout (string equality);
* (Share)   the ccw accepted-power factor ratio f3 / f1.5 from the P61 tokens equals P83's
            full-token power share within the sum of the three tokens' half-units;
* (P83-Q)   each of P83's four THETA rows has (3 T token / 1.5 T token) equal to the share within
            its three tokens' half-units;
* (N-sense) control: the cw factor ratio in place of the ccw turns (Share) red;
* (N-plant) control: one digit planted in entry 60's ccw ball token misses (R39) on the ccw ball
            row alone while (D) stays green;
* (N-swap)  control: the sides swapped miss (R39) on all 22 rows and turn (Node) red.

Printed only (<= 7 significant digits): C39's 22-row table, the ccw measured / predicted
quotients, the heating and SAR shape parts, the ccw hotspot ratio over the SAR lever, the cw
rows without a P83 comparand, the ``cem43`` ratio with its caveat.

Run it::

    ./scripts/run_examples.sh -e mat:42 -n 1 -t 30

Not claimed: a 3 T gel claim (the 1.5 T gel's constants at both field strengths), which
closed-form picture the heating follows, a thermal or box-wide maximum, a dose band, ``cem43`` as
a lever, a decision rule or ``ANS`` thermal case, degree 2, tissue, compliance, an absolute or
external figure.
"""

from __future__ import annotations

import contextlib
import copy
import io
import re
import sys
from pathlib import Path

from mpi4py import MPI

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from tests.validation import test_geo38_c39_thermal_reader as c46  # noqa: E402
from tests.validation.test_geo38_c39_thermal_reader import (  # noqa: E402
    ENTRY59_GLOB,
    ENTRY59_NODE_PF,
    ENTRY60_GLOB,
    ENTRY60_NODE_PF,
    SENSES,
    _half_ulp_rel,
    _read_log,
    check_block,
    find_log,
    parse_heat_log,
    ratios,
    side,
)

LOG_DIR = _REPO_ROOT / "docs" / "testing" / "logs"
#: C39's closing ``-n 1`` log (named in item 953's Outcome; the ``-n 2`` log prints the table twice).
C39_LOG = LOG_DIR / "20261111T103037Z_GEO-38-C39-W1-n1.log"
#: P83's log (``MAT-21`` step 46).
P83_LOG = LOG_DIR / "20261110T103343Z_MAT-21-P83-W1.log"
BASE59 = "20261110T080008Z_GEO-38-P33.log"
BASE60 = "20261111T080008Z_GEO-38-P33-3T.log"
TAG = "[EX-194]"

_RX_ROW = re.compile(r"\[GEO-38 C46\]\s+(ccw|cw) (.+?): (\S+) \| (\S+) \| (\S+)\s*$")
_RX_SHARE = re.compile(r"\(K\) full-token power share s3/s1\.5 (\S+) vs C33")
_RX_THETA = re.compile(r"MAT21-46 THETA (.+?) \| (\S+) \| (\S+) \| (\S+) \|")
_RX_DOSE = re.compile(r"MAT21-46 DOSE ratio 3T/1\.5T (.+?): excess (\S+) \(")
_RX_SHAPE = re.compile(r"MAT21-46 BALL shape part 3T/1\.5T (.+?): (\S+) \(")
_RX_C33 = re.compile(r"c33\.ratios at \.6e \('([^']+)', '([^']+)', '([^']+)', '([^']+)'\) vs "
                     r"read_c33_rows")
_RX_BOUND = re.compile(r"rel (\S+) \(<= (\S+)\)")


def plant_digit(tok, place=4):
    """``tok`` with its ``place``-th significant digit raised by one (mod 10)."""
    mant, e, exp = tok.partition("e")
    digits = [i for i, c in enumerate(mant) if c.isdigit()]
    chars = list(mant)
    i = digits[place - 1]
    chars[i] = str((int(chars[i]) + 1) % 10)
    return "".join(chars) + e + exp


def read_c39_table(text):
    """``{(sense, key): (s15, s3, ratio)}`` strings of C39's real table, and its header line.

    The log also holds STUB tables (the stub-dry test); the block read is the one non-STUB header
    naming both real ``xl`` logs, and exactly one such header must exist."""
    lines = text.splitlines()
    heads = [i for i, ln in enumerate(lines)
             if "C39 table:" in ln and "STUB" not in ln and BASE59 in ln and BASE60 in ln]
    if len(heads) != 1:
        raise ValueError(f"C39 real table headers naming {BASE59} and {BASE60}: {len(heads)}")
    rows = {}
    for ln in lines[heads[0] + 1:heads[0] + 23]:
        m = _RX_ROW.search(ln)
        if not m:
            raise ValueError(f"C39 table row does not parse: {ln[:80]!r}")
        key = (m.group(1), m.group(2))
        if key in rows:
            raise ValueError(f"C39 table row repeated: {key}")
        rows[key] = (m.group(3), m.group(4), m.group(5))
    return rows, lines[heads[0]]


def read_p83(text):
    """P83's share token, THETA / DOSE / shape rows and the C33 ratio tuple (strings)."""
    share = _RX_SHARE.findall(text)
    theta = [m for m in _RX_THETA.findall(text) if m[0] != "picture"]
    dose = _RX_DOSE.findall(text)
    shape = _RX_SHAPE.findall(text)
    c33 = _RX_C33.findall(text)
    return {"share": share, "theta": theta, "dose": dose, "shape": shape, "c33": c33}


def rel_bound(toks):
    return sum(_half_ulp_rel(t) for t in toks)


def factor_ratio(b15, b3, sense="ccw"):
    t15 = b15["senses"][sense]["p61"]["src"]["factor"]
    t3 = b3["senses"][sense]["p61"]["src"]["factor"]
    return float(t3) / float(t15), (t15, t3)


def share_gate(b15, b3, share_tok, sense="ccw"):
    """(Share): ``(rel, bound)`` of f3 / f1.5 against P83's full-token share."""
    r, toks = factor_ratio(b15, b3, sense)
    return abs(r / float(share_tok) - 1.0), rel_bound((*toks, share_tok)), r


def quiet_check(block, node, names, name):
    """``check_block`` with its progress prints (10-digit tokens) swallowed."""
    with contextlib.redirect_stdout(io.StringIO()):
        return check_block(block, node, names, name)


def _on_root(fn):
    """Run ``fn`` on rank 0, bcast its result or its exception text to every rank."""
    comm = MPI.COMM_WORLD
    out = None
    if comm.rank == 0:
        try:
            out = ("ok", fn())
        except BaseException:  # re-raised on every rank below, as text
            import traceback

            out = ("err", traceback.format_exc())
    out = comm.bcast(out, root=0)
    if out[0] == "err":
        raise RuntimeError("rank 0 failed:\n" + out[1])
    return out[1]


def build():
    res, out = [], []

    def gate(name, ok, detail):
        res.append((name, bool(ok), detail))

    def say(s):
        out.append(f"{TAG} {s}")

    f6 = lambda v: f"{v:.6e}"  # noqa: E731
    names = list(c46._probes())

    hit59, hit60 = find_log(ENTRY59_GLOB), find_log(ENTRY60_GLOB)
    if hit59 is None or hit60 is None:
        raise ValueError(f"entry 59 / 60 log absent: {hit59} / {hit60}")
    (path59, kind59), (path60, kind60) = hit59, hit60
    b59, b60 = parse_heat_log(path59), parse_heat_log(path60)
    rows39, head39 = read_c39_table(_read_log(C39_LOG))
    p83 = read_p83(_read_log(P83_LOG))
    say(f"entry 59 {b59['log']} ({kind59}), entry 60 {b60['log']} ({kind60})")
    say(f"C39 closing log {C39_LOG.name}; P83 log {P83_LOG.name}")

    # (X)
    gate("(X) logs are the two named base names",
         (b59["log"], b60["log"]) == (BASE59, BASE60), f"{b59['log']} | {b60['log']}")
    gate("(X) node strings", b59["head"]["node"] == ENTRY59_NODE_PF
         and b60["head"]["node"] == ENTRY60_NODE_PF,
         f"{b59['head']['node']} / {b60['head']['node']} pF")
    gate("(X) probe names == c46._probes() in order (both sides, both senses)",
         all(list(b["senses"][s]["p61"]["probes"]) == names for b in (b59, b60) for s in SENSES),
         f"{names}")
    gate("(X) C39 real table found once, 22 rows", len(rows39) == 22, f"{len(rows39)} rows")
    counts = (len(p83["share"]), len(p83["theta"]), len(p83["dose"]), len(p83["shape"]),
              len(p83["c33"]))
    gate("(X) P83 rows found with counts 1 / 4 / 4 / 3 (+ C33 tuple 1)",
         counts == (1, 4, 4, 3, 1), f"{counts}")
    if counts != (1, 4, 4, 3, 1):
        raise ValueError(f"P83 row counts {counts}")
    share_tok = p83["share"][0]

    # (D)
    worst = {"(D-f)": (0.0, 1.0), "(D-scale)": (0.0, 1.0)}
    for b, node, nm in ((b59, ENTRY59_NODE_PF, "1.5 T (entry 59)"),
                        (b60, ENTRY60_NODE_PF, "3 T (entry 60)")):
        for n, ok, det in quiet_check(b, node, names, nm):
            gate(f"(D) {n}", ok, "ok" if ok else "FAILED")
            m = _RX_BOUND.search(det)
            for tag in worst:
                if tag in n and m and float(m.group(1)) / float(m.group(2)) >= \
                        worst[tag][0] / worst[tag][1]:
                    worst[tag] = (float(m.group(1)), float(m.group(2)))

    s15, s60 = side(b59), side(b60)
    r = ratios(s15, s60)

    # (R39)
    miss = []
    for sn in SENSES:
        for k, v in r[sn].items():
            row = rows39.get((sn, k))
            if row is None or f"{v:.6e}" != row[2]:
                miss.append((sn, k))
    gate("(R39) 22 page ratios == C39 closing-log strings", not miss and len(rows39) == 22,
         f"{22 - len(miss)} of 22 equal" + (f"; MISS {miss}" if miss else ""))

    # (Share), (P83-Q)
    rel_s, bnd_s, f_ccw = share_gate(b59, b60, share_tok)
    gate("(Share) ccw f3/f1.5 == P83 power share", rel_s <= bnd_s,
         f"rel {rel_s:.3e} (<= {bnd_s:.3e})")
    for pic, t15, t3, tr in p83["theta"]:
        rel = abs(float(t3) / float(t15) / float(share_tok) - 1.0)
        bnd = rel_bound((t15, t3, share_tok))
        gate(f"(P83-Q) THETA {pic}: 3T/1.5T == share", rel <= bnd, f"rel {rel:.3e} (<= {bnd:.3e})")

    # Controls
    rel_c, bnd_c, f_cw = share_gate(b59, b60, share_tok, sense="cw")
    gate("(N-sense) cw f3/f1.5 in place of ccw: (Share) red", rel_c > bnd_s,
         f"residual {rel_c:.3e} vs bound {bnd_s:.3e} ({rel_c / bnd_s:.2e}x)")
    b60x = copy.deepcopy(b60)
    ball = b60x["senses"]["ccw"]["i32"]["ball"]
    ball["theta"] = plant_digit(ball["theta"], 4)
    rx = ratios(s15, side(b60x))
    missx = [(sn, k) for sn in SENSES for k, v in rx[sn].items()
             if f"{v:.6e}" != rows39[(sn, k)][2]]
    key_ball = next(k for k in s15["ccw"] if k.startswith("ball theta"))
    gate("(N-plant) one digit planted in entry 60's ccw ball token: (R39) misses the ccw ball row "
         "alone", missx == [("ccw", key_ball)],
         f"{len(missx)} row(s) miss: {[(s, k[:20]) for s, k in missx]}; ccw ball rel "
         f"{abs(rx['ccw'][key_ball] / r['ccw'][key_ball] - 1.0):.3e}")
    dx = [ok for _n, ok, _d in quiet_check(b60x, ENTRY60_NODE_PF, names, "planted")]
    gate("(N-plant) (D) stays green on the planted copy (check_block does not read the ball)",
         all(dx), f"{sum(dx)} of {len(dx)} checks green")
    r_sw = ratios(s60, s15)
    miss_sw = [(sn, k) for sn in SENSES for k, v in r_sw[sn].items()
               if f"{v:.6e}" != rows39[(sn, k)][2]]
    gate("(N-swap) swapped sides miss (R39) on all 22 rows", len(miss_sw) == 22,
         f"{len(miss_sw)} of 22 rows miss")
    node_sw = [(n, ok) for n, ok, _ in quiet_check(b60, ENTRY59_NODE_PF, names, "swap")
               if "(Node)" in n]
    gate("(N-swap) check_block(b60, ENTRY59_NODE_PF) (Node) red",
         len(node_sw) == 1 and not node_sw[0][1], f"{node_sw}")

    # ---- Printed only (<= 7 significant digits) ----
    say(f"D bounds, worst error / bound: (D-f) {worst['(D-f)'][0]:.3e} / {worst['(D-f)'][1]:.3e}; "
        f"(D-scale) {worst['(D-scale)'][0] / worst['(D-scale)'][1]:.3f} of its bound")
    say(f"C39 table, 3T/1.5T (1.5 T node {ENTRY59_NODE_PF} pF | 3 T node {ENTRY60_NODE_PF} pF), "
        "per accepted watt (K per W; cem43 in s):")
    for sn in SENSES:
        for k in s15[sn]:
            say(f"  {sn} {k}: {f6(s15[sn][k])} | {f6(s60[sn][k])} | ratio {f6(r[sn][k])}")

    def key(prefix):
        ks = [k for k in s15["ccw"] if k.startswith(prefix)]
        if len(ks) != 1:
            raise ValueError(f"readout key {prefix!r}: {ks}")
        return ks[0]

    k_hot, k_ball = key("hotspot point"), key("ball theta")
    k_win, k_exc, k_cem = key("running mean"), key("cem43 excess"), key("cem43 per")
    share = float(share_tok)
    shape_p83 = {pic: tok for pic, tok in p83["shape"]}
    sar_lever = float(p83["c33"][0][2])
    shape_sar = float(next(iter(shape_p83.values())))
    say(f"P83 power share s3/s1.5 (full token): {share:.7f}; P83 prediction for every "
        "mirror-even picture is the share alone (shape factor 1)")
    for lab, k in (("hotspot point", k_hot), ("ball", k_ball), ("360 s window", k_win),
                   ("cem43 excess", k_exc)):
        v = r["ccw"][k]
        say(f"  ccw {lab}: measured {v:.6f} / predicted {share:.6f} = {v / share:.6f}")
    shape_heat = r["ccw"][k_hot] / share
    say(f"  heating shape part (ccw hotspot ratio / power share) {shape_heat:.6f}; SAR shape "
        f"part (P83's three BALL rows, all equal) {shape_sar:.6f}")
    say(f"  ccw hotspot ratio / (share x SAR shape part) {r['ccw'][k_hot] / (share * shape_sar):.6f}")
    say(f"  SAR lever (C33) {sar_lever:.6f}; ccw hotspot ratio / SAR lever "
        f"{r['ccw'][k_hot] / sar_lever:.6f}")
    say("  cem43 excess ratios in P83: " + ", ".join(f"{p} {float(t):.5g}" for p, t in p83["dose"])
        + f"; measured ccw {r['ccw'][k_exc]:.5g}")
    say("  cw rows have no P83 comparand: P83's share is the ccw partition, the cw arm is heated "
        f"at ccw's c*, and cw's f3/f1.5 is {f_cw:.5g} (ccw {f_ccw:.7g})")
    for lab, k in (("hotspot point", k_hot), ("ball", k_ball), ("360 s window", k_win)):
        say(f"  cw {lab} ratio {r['cw'][k]:.6f} (no comparand)")
    say(f"  cem43 per P_acc ratio: ccw {r['ccw'][k_cem]:.6f}, cw {r['cw'][k_cem]:.6f} -- printed "
        "only: a cem43 per accepted watt is the 1 W history's dose, about D0 plus an excess, not "
        "a scalable per-watt dose and never a lever")
    return res, out


def main() -> int:
    comm = MPI.COMM_WORLD
    res, out = _on_root(build)
    failed = [n for n, ok, _ in res if not ok]
    if comm.rank == 0:
        for line in out:
            print(line)
        for n, ok, det in res:
            print(f"{TAG} [{'PASS' if ok else 'FAIL'}] {n}: {det}")
        print(f"{TAG} {len(res) - len(failed)} of {len(res)} gates green"
              + (f"; FAILED {failed}" if failed else ""))
    return 1 if failed else 0


if __name__ == "__main__":
    try:
        code = main()
    except BaseException:
        import traceback

        traceback.print_exc()
        sys.stderr.flush()
        MPI.COMM_WORLD.Abort(1)
        raise
    sys.exit(code)
