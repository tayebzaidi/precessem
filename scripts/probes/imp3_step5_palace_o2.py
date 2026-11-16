"""`IMP-3` step 5 (§9 item 252, chain I13) — Palace order 2 against ours at degree 2.

Two subcommands, both serial in the complex ``fem-em-solver`` container:

* ``configs`` (W1) — (X0) re-runs item 242's census on ``imp3x_off.msh`` /
  ``imp3x_rod.msh``, loads ``imp3p_off.json`` / ``imp3p_rod.json`` read-only,
  checks Order 1 / PEC / Postprocessing / Output, and writes
  ``imp3q_<arm>.json`` with ``Solver.Order = 2`` and Output
  ``postpro_imp3q_<arm>`` (mode ``"x"``).  (X0q) prints every key path where
  ``imp3q_off.json`` differs from ``logs/ops61/f05x_o2.json`` (Mesh / Output /
  Postprocessing excluded); predicted empty, printed not asserted.
* ``read`` (W4) — (H) header bijection on each ``postpro_imp3q_*``, (P) via
  ``_p_check``, (RP2) ``imp3q_off`` S against ``logs/ops61/postpro_f05x_o2/
  port-S.csv`` within ``REUSE_REPRODUCTION_RTOL``, Palace's own order 1 -> 2
  shifts against ``palace_imp3p_fields.json``; parsed arrays saved to
  ``logs/imp3/palace_imp3q_fields.json`` (mode ``"x"``).

Loads item 242's probe by importlib; edits no earlier probe, never writes to
``logs/ops61/`` or to any ``imp3x_*`` / ``imp3p_*`` artifact.
"""

from __future__ import annotations

import copy
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

ARMS = ("off", "rod")
ARM_PEC = {"off": [1], "rod": [1, 402]}
FIELDS_JSON = "palace_imp3q_fields.json"
F05X_O2 = REPO / "logs" / "ops61" / "f05x_o2.json"
F05X_O2_S = REPO / "logs" / "ops61" / "postpro_f05x_o2" / "port-S.csv"


def _step4():
    path = REPO / "scripts" / "probes" / "imp3_step4_palace_fields.py"
    spec = importlib.util.spec_from_file_location("imp3_step4_palace_fields", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _diff_paths(a, b, prefix=""):
    """Every key path where nested dict/list ``a`` differs from ``b``."""
    out = []
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            p = f"{prefix}.{k}" if prefix else str(k)
            if k not in a or k not in b:
                out.append(f"{p} (only in {'imp3q' if k in a else 'f05x_o2'})")
            else:
                out += _diff_paths(a[k], b[k], p)
    elif isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
        for i, (x, y) in enumerate(zip(a, b)):
            out += _diff_paths(x, y, f"{prefix}[{i}]")
    elif a != b:
        out.append(f"{prefix}: {a!r} != {b!r}")
    return out


def configs():
    from mpi4py import MPI

    assert MPI.COMM_WORLD.size == 1, "configs is serial: run it at mpiexec -n 1"
    s4 = _step4()
    out = s4._step3().OUT_DIR
    tag = "[IMP-3 5 W1]"
    meta = json.loads((out / "imp3x_meta.json").read_text())
    fails = []
    for arm in ARMS:
        ntet, attrs = s4._msh_census(out / f"imp3x_{arm}.msh")
        want = {1, 21, 22, 23, 24} | ({402} if arm == "rod" else set())
        ok = ntet == meta[arm]["cells"] and attrs == want
        print(f"{tag} (X0) imp3x_{arm}.msh: {ntet} tets (meta {meta[arm]['cells']}), "
              f"triangle attributes {sorted(attrs)} (want {sorted(want)}): "
              f"{'PASS' if ok else 'FAIL'}", flush=True)
        if not ok:
            fails.append(f"X0 {arm}")
    bases = {}
    for arm in ARMS:
        cfg = json.loads((out / f"imp3p_{arm}.json").read_text())
        order = cfg["Solver"]["Order"]
        pec = cfg["Boundaries"]["PEC"]["Attributes"]
        post = sorted(cfg["Domains"]["Postprocessing"])
        outp = cfg["Problem"]["Output"]
        ok = (order == 1 and pec == ARM_PEC[arm] and post == ["Energy", "Probe"]
              and outp == f"postpro_imp3p_{arm}")
        print(f"{tag} (X0) imp3p_{arm}.json Order {order}, PEC {pec}, Postprocessing {post}, "
              f"Output {outp}: {'PASS' if ok else 'FAIL'}", flush=True)
        if not ok:
            fails.append(f"X0 cfg {arm}")
        bases[arm] = cfg
    assert not fails, f"(X0) red: {fails}"
    targets = [out / f"imp3q_{a}.json" for a in ARMS]
    for t in targets:
        assert not t.exists(), f"{t} exists: never overwrite"
    for arm in ARMS:
        cfg = copy.deepcopy(bases[arm])
        cfg["Solver"]["Order"] = 2
        cfg["Problem"]["Output"] = f"postpro_imp3q_{arm}"
        path = out / f"imp3q_{arm}.json"
        with open(path, "x") as fh:
            fh.write(json.dumps(cfg, indent=2) + "\n")
        print(f"{tag} wrote logs/imp3/{path.name}: Order {cfg['Solver']['Order']}, Output "
              f"{cfg['Problem']['Output']}, PEC {cfg['Boundaries']['PEC']['Attributes']}, "
              f"mesh {cfg['Model']['Mesh']}", flush=True)
    q = json.loads((out / "imp3q_off.json").read_text())
    ref = json.loads(F05X_O2.read_text())
    for c in (q, ref):
        c.get("Model", {}).pop("Mesh", None)
        c.get("Problem", {}).pop("Output", None)
        c.get("Domains", {}).pop("Postprocessing", None)
    diffs = _diff_paths(q, ref)
    print(f"{tag} (X0q) imp3q_off.json vs logs/ops61/f05x_o2.json differing key paths "
          f"(Mesh/Output/Postprocessing excluded; predicted empty): {len(diffs)}", flush=True)
    for d in diffs:
        print(f"{tag} (X0q)   {d}", flush=True)
    print(f"{tag} ANCHORS PASS: (X0)", flush=True)


def read():
    from tests.validation.test_port_birdcage_ring_matrix import REUSE_REPRODUCTION_RTOL

    s4 = _step4()
    s3 = s4._step3()
    o61 = s3._ops61()
    out = s3.OUT_DIR
    tag = "[IMP-3 5 W4]"
    target = out / FIELDS_JSON
    assert not target.exists(), f"{target} exists: never overwrite"
    names = s4.PROBE_NAMES
    fails = []
    e, w, s = {}, {}, {}
    for arm in ARMS:
        d = out / f"postpro_imp3q_{arm}"
        ph, hdr = s4.parse_probe(d / "probe-E.csv")
        print(f"{tag} (H) {arm} probe-E.csv header ({len(hdr)} cols): {hdr}", flush=True)
        ew, hdr = s4.parse_energy(d / "domain-E.csv")
        print(f"{tag} (H) {arm} domain-E.csv header ({len(hdr)} cols): {hdr}", flush=True)
        e[arm], w[arm] = ph, ew
        name = f"imp3q_{arm}"
        s[arm], _ident, f = o61._p_check(name, d / "port-S.csv", tag)
        fails += [f"{name}: {x}" for x in f]
        for k, nm in enumerate(names):
            print(f"{tag} {arm} {nm} P1 E = " + ", ".join(o61._fmt(z) for z in ph[k, :, 0])
                  + f"  |E| {np.linalg.norm(ph[k, :, 0]):.9e}", flush=True)
        print(f"{tag} {arm} W_elec(phantom) per exc " + ", ".join(f"{v:.12e}" for v in ew),
              flush=True)
    s_ref = o61._read_port_s(F05X_O2_S)[0]
    rp = float(np.max(np.abs(s["off"] - s_ref)) / np.max(np.abs(s_ref)))
    ok = rp <= REUSE_REPRODUCTION_RTOL
    print(f"{tag} (RP2) imp3q_off vs logs/ops61/postpro_f05x_o2 entrywise max|ΔS|/max|S| "
          f"{rp:.3e} <= {REUSE_REPRODUCTION_RTOL}: {'PASS' if ok else 'FAIL'}", flush=True)
    if not ok:
        fails.append(f"RP2 {rp}")
    # Palace's own order 1 -> 2 shifts (item 242's order-1 arrays, read-only).
    o1 = json.loads((out / s4.FIELDS_JSON).read_text())
    e1 = {a: np.array(o1["E"][a])[..., 0] + 1j * np.array(o1["E"][a])[..., 1] for a in ARMS}
    w1 = {a: np.array(o1["W_elec"][a]) for a in ARMS}
    mag2 = {a: np.linalg.norm(e[a], axis=1) for a in ARMS}
    mag1 = {a: np.linalg.norm(e1[a], axis=1) for a in ARMS}
    for arm in ARMS:
        dl = mag2[arm][:, 0] / mag1[arm][:, 0] - 1
        print(f"{tag} Palace δ_p (o2/o1 - 1, P1) {arm}: "
              + "  ".join(f"{nm} {dl[k]:+.6e}" for k, nm in enumerate(names)), flush=True)
    for j in range(s4.N_EXC):
        r2 = mag2["rod"][:, j] / mag2["off"][:, j]
        r1 = mag1["rod"][:, j] / mag1["off"][:, j]
        print(f"{tag} Palace η_rod^o2 P{j + 1}: "
              + "  ".join(f"{nm} {r2[k]:.6f}" for k, nm in enumerate(names)), flush=True)
        print(f"{tag} Palace η_rod^o2/η_rod^o1 - 1 P{j + 1}: "
              + "  ".join(f"{nm} {r2[k] / r1[k] - 1:+.6e}" for k, nm in enumerate(names)),
              flush=True)
        print(f"{tag} Palace ΔW(rod/off) P{j + 1}: o2 {w['rod'][j] / w['off'][j] - 1:+.6e}"
              f"  o1 {w1['rod'][j] / w1['off'][j] - 1:+.6e}", flush=True)
    with open(target, "x") as fh:
        json.dump({
            "order": 2,
            "probe_names": list(names),
            "points": [list(p) for p in s4.PROBE_POINTS],
            "E": {a: [[[[z.real, z.imag] for z in e[a][k, c, :]] for c in range(3)]
                      for k in range(6)] for a in e},
            "W_elec": {a: w[a].tolist() for a in w},
        }, fh)
    print(f"{tag} wrote logs/imp3/{FIELDS_JSON}", flush=True)
    assert not fails, f"(P)/(RP2) red: {fails}"
    print(f"{tag} (H) PASS, (P) PASS on imp3q_off and imp3q_rod, (RP2) PASS", flush=True)


if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in ("configs", "read"):
        raise SystemExit("usage: imp3_step5_palace_o2.py configs|read")
    {"configs": configs, "read": read}[sys.argv[1]]()
