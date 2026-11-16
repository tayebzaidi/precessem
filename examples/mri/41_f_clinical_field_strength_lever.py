"""Example (`EX-193`): the field-strength lever on one page -- F-clinical at 1.5 T and 3 T.

A **record page, no solve and no mesh** (the ``ans:11`` precedent).  It reads two
tracked ``GEO-38`` C6a logs -- ``xl`` entry 52 (1.5 T, 64 MHz, node 193.5996 pF) and
entry 58 (3 T, 128 MHz, node 44.0525 pF), both degree 1, each at its own unloaded
in-model root -- and C33's closing log, and sets the two field strengths side by side.
The new angle is the **factoring of the lever per accepted watt**: the SAR10g(c*) ratio
is (how much more of the accepted power the gel box takes, P_box / P_acc) x (how much
more SAR each box watt makes at c*), and the |B1+| ratio factors the same way with a
square root.  The hotspot "move" is also measured as a margin over the runner-up ball.

Every parser, helper and constant is imported by name from the gate modules (the ANS-1
rule): ``tests/validation/test_geo38_c33_field_strength_reader.py`` (``_load_page``,
``find_log``, ``side``, ``check_side``, ``ratios``, ``_half_ulp_rel``, the node
constants) and, through its ``_load_page()``, the ``ans:11`` page (``parse_c6_log``,
``derived``, ``_read_log``, ``SOURCE_LOG``).  No ``test_*`` function and no gated body
is called; no log regex of the C6a lines is re-implemented.  There is no gate-module
edit and no ``src/`` edit.

Asserted (all collected first, then asserted):

* (X)      both logs parse; node strings equal the node constants; 18 lattice lines per
           side; the labels agree in order across the sides;
* (D)      ``check_side`` green on both sides at C33's own half-unit bounds;
* (R33)    the anchor: each of C33's four ratios at ``.6e`` equals the string C33's
           closing log printed, and both c* labels equal its c* row (string equality);
* (Dec-S)  R(SAR10g(c*) per P_acc) = R(per box watt) x R(P_box / P_acc), at the sum of
           the six tokens' half-units;
* (Dec-B)  R(mean |B1+| per sqrt(P_acc)) = R(per sqrt(box W)) x sqrt(R(P_box / P_acc));
* (Lat)    the 18 parsed |B1+| tokens / sqrt(P_acc) reproduce ``side()``'s lattice mean;
* (Arg)    the runner-up ball is not c* on either side and the margin is > 0;
* (N-norm) control: R(per P_src) in place of R(per box watt) turns (Dec-S) red;
* (N-plant) control: one digit planted (4th significant place) in the 3 T per-box-watt
           string turns (Dec-S) red while (R33) stays green;
* (N-swap) control: the sides swapped miss (R33) on all four rows and turn (Node) red.

Printed only: C33's table with its two factors, the 18-row per-point table, the c*
labels / runner-ups / margins, the SAR ratio at a fixed ball, SAR10g per |B1+|^2, the
P_acc / P_src share, and Im Y_in under both senses side by side (never a ratio).

Run it through the example runner (the ``mri:`` group sources the complex build)::

    ./scripts/run_examples.sh -e mri:41 -n 2 -t 120

Not claimed: a tuning, absolute, compliance, homogeneity or external number; degree 2;
a box-wide maximum; a quasi-static law; any band on a ratio.
"""

from __future__ import annotations

import copy
import math
import re
import sys
from pathlib import Path

from mpi4py import MPI

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from tests.validation.test_geo38_c33_field_strength_reader import (  # noqa: E402
    ENTRY52_NODE_PF,
    ENTRY58_GLOB,
    ENTRY58_NODE_PF,
    _half_ulp_rel,
    _load_page,
    check_side,
    find_log,
    ratios,
    side,
)

#: C33's closing log (tracked; gzip-tolerant through the page's ``_read_log``).
C33_LOG = _REPO_ROOT / "docs" / "testing" / "logs" / "20261109T103601Z_GEO-38-C33-W1-n1-r2.log"
#: Local copies of C33's ``print_table`` row titles, keyed by ``ratios()`` key (disclosed).
C33_TITLES = {
    "lattice_mean_per_sqrt_p_acc": "lattice mean |B1+|_ccw / sqrt(P_acc) [T/sqrt(W)]",
    "centroid_mean_per_sqrt_p_acc": "centroid mean |B1+|_ccw / sqrt(P_acc) [T/sqrt(W)]",
    "sar_per_p_acc": "SAR10g(c*) per P_acc [W/kg/W]",
    "p_box_over_p_acc": "P_box / P_acc",
}
_LABEL = r"(x[+-]\d+mm dz[+-]\d+mm)"
TAG = "[EX-193]"


# ---------------------------------------------------------------------------
# Pure helpers.
# ---------------------------------------------------------------------------
def read_c33_rows(text, log52, log58):
    """``{key: (s15, s3, ratio)}`` strings of C33's four table rows and the two c* labels.

    C33's closing log also holds its self-test table (a renamed copy of entry 52), so the
    block read is the one whose header names both real logs; exactly one such header."""
    lines = text.splitlines()
    heads = [i for i, ln in enumerate(lines) if "C33 table:" in ln and log52 in ln and log58 in ln]
    if len(heads) != 1:
        raise ValueError(f"C33 table headers naming {log52} and {log58}: {len(heads)}")
    text = "\n".join(lines[heads[0]:heads[0] + 10])
    rows = {}
    for key, title in C33_TITLES.items():
        rx = re.compile(r"\[GEO-38 C44\]\s+" + re.escape(title) + r": (\S+) \| (\S+) \| (\S+)\s*$",
                        re.M)
        hits = rx.findall(text)
        if len(hits) != 1:
            raise ValueError(f"C33 row {title!r}: {len(hits)} matches")
        rows[key] = hits[0]
    rx = re.compile(r"\[GEO-38 C44\]\s+c\* \(hotspot\): " + _LABEL + r" \([^)]*\) \| " + _LABEL)
    hits = rx.findall(text)
    if len(hits) != 1:
        raise ValueError(f"C33 c* row: {len(hits)} matches")
    return rows, hits[0]


def rel_bound(tokens):
    return sum(_half_ulp_rel(t) for t in tokens)


def ratio_of(p52, p58, key):
    return float(p58["normalisers"][key]) / float(p52["normalisers"][key])


def dec_s(p52, p58, factor_key="sar_per_box_w"):
    """(Dec-S): ``(rel, bound)`` of R(sar_per_p_acc) against R(factor) x R(P_box/P_acc)."""
    lhs = ratio_of(p52, p58, "sar_per_p_acc")
    rhs = ratio_of(p52, p58, factor_key) * ratio_of(p52, p58, "p_box_over_p_acc")
    toks = [p[k] for p in (p52["normalisers"], p58["normalisers"])
            for k in ("sar_per_p_acc", factor_key, "p_box_over_p_acc")]
    return abs(lhs - rhs) / abs(lhs), rel_bound(toks), lhs, rhs


def dec_b(p52, p58):
    lhs = ratio_of(p52, p58, "b1_per_sqrt_p_acc")
    rhs = ratio_of(p52, p58, "b1_per_sqrt_box_w") * math.sqrt(ratio_of(p52, p58, "p_box_over_p_acc"))
    toks = [p[k] for p in (p52["normalisers"], p58["normalisers"])
            for k in ("b1_per_sqrt_p_acc", "b1_per_sqrt_box_w", "p_box_over_p_acc")]
    return abs(lhs - rhs) / abs(lhs), rel_bound(toks), lhs, rhs


def plant_digit(tok, place=4):
    """``tok`` with its ``place``-th significant digit raised by one (mod 10)."""
    mant, e, exp = tok.partition("e")
    digits = [i for i, c in enumerate(mant) if c.isdigit()]
    i = digits[place - 1]
    chars = list(mant)
    chars[i] = str((int(chars[i]) + 1) % 10)
    return "".join(chars) + e + exp


def runner_up(parsed, d):
    """``(label, sar_per_p_acc, c*_label)`` of the second-largest lattice SAR10g."""
    order = sorted(range(len(d["sar_per_p_acc"])), key=lambda i: float(parsed["lattice"][i]["sar10g"]))
    top, second = order[-1], order[-2]
    return parsed["lattice"][second]["label"], float(parsed["lattice"][second]["sar10g"]), \
        parsed["lattice"][top]["label"], float(parsed["lattice"][top]["sar10g"])


def _on_root(fn):
    """Run ``fn`` on rank 0, bcast its result or its exception to every rank."""
    comm = MPI.COMM_WORLD
    out = None
    if comm.rank == 0:
        try:
            out = ("ok", fn())
        except BaseException:  # re-raised on every rank below, as text (always picklable)
            import traceback

            out = ("err", traceback.format_exc())
    out = comm.bcast(out, root=0)
    if out[0] == "err":
        raise RuntimeError("rank 0 failed:\n" + out[1])
    return out[1]


# ---------------------------------------------------------------------------
# The record (rank 0).
# ---------------------------------------------------------------------------
def build():
    page = _load_page()
    res, out = [], []

    def gate(name, ok, detail):
        res.append((name, bool(ok), detail))

    def say(s):
        out.append(f"{TAG} {s}")

    p52 = page.parse_c6_log(page.SOURCE_LOG)
    hit = find_log(ENTRY58_GLOB)
    path58, kind = hit if hit is not None else (None, "absent")
    if hit is None:
        raise ValueError(f"entry 58 log {ENTRY58_GLOB} absent")
    p58 = page.parse_c6_log(path58)
    c33_text = page._read_log(C33_LOG)
    rows33, (cs52_c33, cs58_c33) = read_c33_rows(c33_text, p52['log'], p58['log'])
    say(f"entry 52 {p52['log']}, entry 58 {p58['log']} ({kind}), C33 {C33_LOG.name}")

    # (X)
    labels52 = [r["label"] for r in p52["lattice"]]
    labels58 = [r["label"] for r in p58["lattice"]]
    gate("(X) node strings", p52["node"]["c_pf"] == ENTRY52_NODE_PF
         and p58["node"]["c_pf"] == ENTRY58_NODE_PF,
         f"{p52['node']['c_pf']} / {p58['node']['c_pf']} pF")
    gate("(X) 18 lattice lines per side, labels equal in order",
         len(labels52) == 18 and len(labels58) == 18 and labels52 == labels58,
         f"{len(labels52)} / {len(labels58)} lines")
    gate("(X) entry 58 log is tracked", kind == "tracked", kind)

    # (D)
    for p, node, nm in ((p52, ENTRY52_NODE_PF, "1.5 T (entry 52)"), (p58, ENTRY58_NODE_PF, "3 T (entry 58)")):
        for n, ok, det in check_side(p, node, nm):
            if "(D-cstar)" in n or "(D-S12)" in n:  # their details echo 10-digit log tokens
                det = "strings equal" if ok else "strings differ"
            elif "(D-Pacc)" in n:  # keep the residual, drop the echoed 10-digit token
                det = "rel " + det.split(": rel ")[-1]
            gate(f"(D) {n}", ok, det)

    s15, s3 = side(p52), side(p58)
    d52, d58 = page.derived(p52), page.derived(p58)
    r = ratios(s15, s3)

    # (R33)
    for key in C33_TITLES:
        mine = f"{r[key]:.6e}"
        gate(f"(R33) {key}", mine == rows33[key][2], f"{mine} vs C33 printed {rows33[key][2]}")
    gate("(R33) c* labels == C33 c* row",
         (s15["c_star"]["label"], s3["c_star"]["label"]) == (cs52_c33, cs58_c33),
         f"{s15['c_star']['label']} | {s3['c_star']['label']} vs {cs52_c33} | {cs58_c33}")

    # (Dec-S), (Dec-B)
    rel_s, bnd_s, lhs_s, rhs_s = dec_s(p52, p58)
    gate("(Dec-S) R(SAR/P_acc) = R(SAR/P_box) x R(P_box/P_acc)", rel_s <= bnd_s,
         f"{lhs_s:.6e} vs {rhs_s:.6e}: rel {rel_s:.3e} (<= {bnd_s:.3e})")
    rel_b, bnd_b, lhs_b, rhs_b = dec_b(p52, p58)
    gate("(Dec-B) R(B/sqrt(P_acc)) = R(B/sqrt(P_box)) x sqrt(R(P_box/P_acc))", rel_b <= bnd_b,
         f"{lhs_b:.6e} vs {rhs_b:.6e}: rel {rel_b:.3e} (<= {bnd_b:.3e})")

    # (Lat)
    import numpy as np
    for nm, p, s in (("1.5 T", p52, s15), ("3 T", p58, s3)):
        mean = float(np.mean([float(x["b1p_ccw"]) for x in p["lattice"]])) \
            / math.sqrt(float(p["normalisers"]["p_acc"]))
        err = abs(mean - s["lattice_mean_per_sqrt_p_acc"]) / s["lattice_mean_per_sqrt_p_acc"]
        gate(f"(Lat) {nm} lattice mean from the 18 tokens", err <= 1e-12,
             f"{mean:.6e} vs {s['lattice_mean_per_sqrt_p_acc']:.6e}: rel {err:.3e} (<= 1e-12)")

    # (Arg)
    ru, marg = {}, {}
    for nm, p, d in (("1.5 T", p52, d52), ("3 T", p58, d58)):
        lab2, sar2, lab1, sar1 = runner_up(p, d)
        ru[nm] = lab2
        marg[nm] = (sar1 - sar2) / sar2
        gate(f"(Arg) {nm} runner-up ball is not c*, margin > 0",
             lab2 != p["c_star"]["label"] and lab1 == p["c_star"]["label"] and marg[nm] > 0,
             f"runner-up {lab2}, margin {marg[nm]:.6e}")

    # Controls
    rel_n, bnd_n, _, _ = dec_s(p52, p58, factor_key="sar_per_p_src")
    gate("(N-norm) R(per P_src) in place of R(per box watt): (Dec-S) red", rel_n > bnd_s,
         f"residual {rel_n:.6e} vs bound {bnd_s:.3e}")
    p58x = copy.deepcopy(p58)
    tok = p58["normalisers"]["sar_per_box_w"]
    p58x["normalisers"]["sar_per_box_w"] = plant_digit(tok, 4)
    rel_p, bnd_p, _, _ = dec_s(p52, p58x)
    r_x = ratios(s15, side(p58x))
    r33_ok = all(f"{r_x[k]:.6e}" == rows33[k][2] for k in C33_TITLES)
    gate("(N-plant) one digit planted in 3 T per-box-watt: (Dec-S) red", rel_p > bnd_p,
         f"residual {rel_p:.6e} vs bound {bnd_p:.3e}")
    gate("(N-plant) (R33) stays green on the planted copy (independent gates)", r33_ok,
         "all four ratio strings still equal C33's")
    r_sw = ratios(s3, s15)
    miss = [k for k in C33_TITLES if f"{r_sw[k]:.6e}" != rows33[k][2]]
    gate("(N-swap) swapped sides miss (R33) on all four rows", len(miss) == 4,
         f"{len(miss)} of 4 rows miss")
    node_sw = [(n, ok) for n, ok, _ in check_side(p58, ENTRY52_NODE_PF, "swap") if "(Node)" in n]
    gate("(N-swap) check_side(p58, ENTRY52_NODE_PF) (Node) red",
         len(node_sw) == 1 and not node_sw[0][1], f"{node_sw}")

    # ---- Printed only (<= 7 significant digits) ----
    f6 = lambda v: f"{v:.6e}"  # noqa: E731
    say(f"C33 table, 3T/1.5T (1.5 T node {s15['node_pf']} pF | 3 T node {s3['node_pf']} pF):")
    for key, title in C33_TITLES.items():
        say(f"  {title}: {f6(float(s15[key]))} | {f6(float(s3[key]))} | ratio {f6(r[key])}")
    r_sbw = ratio_of(p52, p58, "sar_per_box_w")
    r_bbw = ratio_of(p52, p58, "b1_per_sqrt_box_w")
    r_pb = r["p_box_over_p_acc"]
    say(f"  SAR lever factored: R(per P_acc) {f6(r['sar_per_p_acc'])} = R(per box watt) {f6(r_sbw)} "
        f"x R(P_box/P_acc) {f6(r_pb)}")
    say(f"  |B1+| lever factored: R(centroid per sqrt P_acc) {f6(r['centroid_mean_per_sqrt_p_acc'])} "
        f"= R(per sqrt box W) {f6(r_bbw)} x sqrt(R(P_box/P_acc)) {f6(math.sqrt(r_pb))}")
    say(f"  P_acc / P_src: 1.5 T {f6(float(p52['normalisers']['p_acc']) / float(p52['normalisers']['p_src']))}"
        f" | 3 T {f6(float(p58['normalisers']['p_acc']) / float(p58['normalisers']['p_src']))}")
    say(f"  c* (hotspot): {s15['c_star']['label']} | {s3['c_star']['label']}")
    for nm in ("1.5 T", "3 T"):
        say(f"  {nm} runner-up ball {ru[nm]}, c* margin (c* - runner-up)/runner-up {f6(marg[nm])}")
    lab15 = s15["c_star"]["label"]
    i15 = labels52.index(lab15)
    fixed = d58["sar_per_p_acc"][i15] / d52["sar_per_p_acc"][i15]
    say(f"  SAR10g per P_acc ratio at the fixed ball {lab15} (1.5 T's c*): {f6(fixed)}; "
        f"at each side's own c*: {f6(r['sar_per_p_acc'])}")
    say(f"  SAR10g(c*) per P_acc ratio / (|B1+| ratio)^2: lattice "
        f"{f6(r['sar_per_p_acc'] / r['lattice_mean_per_sqrt_p_acc'] ** 2)}, centroid "
        f"{f6(r['sar_per_p_acc'] / r['centroid_mean_per_sqrt_p_acc'] ** 2)}")
    say("  per-point table: label | |B1+|ccw/sqrt(P_acc) 1.5 T | 3 T | ratio | SAR10g/P_acc 1.5 T | 3 T | ratio")
    for i, lab in enumerate(labels52):
        b15, b3 = d52["b1p_per_sqrt_p_acc"][i], d58["b1p_per_sqrt_p_acc"][i]
        q15, q3 = d52["sar_per_p_acc"][i], d58["sar_per_p_acc"][i]
        say(f"  {lab} | {f6(b15)} | {f6(b3)} | {f6(b3 / b15)} | {f6(q15)} | {f6(q3)} | {f6(q3 / q15)}")
    for sense in ("cw", "ccw"):  # cw P33 first: C4's root series, C33's order
        for pp in ("p33", "p41"):
            k = f"{pp}_im"
            say(f"  Im Y_in {sense} {pp.upper()} [S]: {float(s15['y_in'][sense][k]):+.6e} | "
                f"{float(s3['y_in'][sense][k]):+.6e} (side by side, never a ratio)")
    return res, out


def main() -> int:
    comm = MPI.COMM_WORLD
    _load_page()  # on every rank: the page's imports may be collective
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
