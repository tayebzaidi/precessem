"""`IMP-3` step 4 (§9 item 242, chain I12) — the near-tip field, ours against Palace.

Two subcommands, both serial in the complex ``fem-em-solver`` container:

* ``configs`` (W1) — (X0) re-reads item 235's exported ``imp3x_off.msh`` /
  ``imp3x_rod.msh`` (tet counts against ``imp3x_meta.json``, triangle attribute
  sets), checks the three ``imp3x_*.json`` PEC lists, and writes
  ``imp3p_<arm>.json`` (item 235's configs + ``Domains.Postprocessing`` Energy
  on the phantom and six Probe points) and the energy-only fallbacks
  ``imp3e_<arm>.json``.  Every target opened with mode ``"x"``.
* ``read`` (W5) — (H) header bijection on ``probe-E.csv`` / ``domain-E.csv``,
  (P) via ``_p_check`` on ``imp3p_off`` / ``imp3p_rod``, (RP) reproduction of
  item 235's S, Palace's own η printed; parsed arrays saved to
  ``logs/imp3/palace_imp3p_fields.json`` (mode ``"x"``).

Imports item 235's probe by importlib; never edits it, never writes to
``logs/ops61/`` or to any ``imp3x_*`` artifact.
"""

from __future__ import annotations

import copy
import csv
import importlib.util
import json
import re
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

ARM_RUN = {"off": "imp3x_off", "rod": "imp3x_rod", "rodnat": "imp3x_rodnat"}
ARM_MESH = {"off": "off", "rod": "rod", "rodnat": "rod"}
ARM_PEC = {"off": [1], "rod": [1, 402], "rodnat": [1]}

# Probe points (metres), the brief's table; identical in the test module.
PROBE_NAMES = ("T1", "T2", "T3", "S1", "F1", "F2")
PROBE_POINTS = (
    (0.015023, 0.000017, 0.017531),
    (0.015023, 0.000017, 0.019031),
    (0.015023, 0.000017, 0.022031),
    (0.017523, 0.000017, 0.000031),
    (-0.014977, 0.000017, 0.000031),
    (0.000023, -0.017983, 0.010031),
)
N_EXC = 4
FIELDS_JSON = "palace_imp3p_fields.json"


def _step3():
    path = REPO / "scripts" / "probes" / "imp3_step3_palace_rod.py"
    spec = importlib.util.spec_from_file_location("imp3_step3_palace_rod", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _msh_census(path):
    """(tet count, triangle attribute set) from a Gmsh 2.2 ASCII file."""
    ntet, tri_attrs = 0, set()
    in_el = False
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("$Elements"):
                in_el = True
                next(fh)
                continue
            if ln.startswith("$EndElements"):
                break
            if not in_el:
                continue
            parts = ln.split()
            etype = int(parts[1])
            if etype == 4:
                ntet += 1
            elif etype == 2:
                tri_attrs.add(int(parts[3]))
    return ntet, tri_attrs


def configs():
    from mpi4py import MPI

    from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG

    assert MPI.COMM_WORLD.size == 1, "configs is serial: run it at mpiexec -n 1"
    s3 = _step3()
    out = s3.OUT_DIR
    tag = "[IMP-3 4 W1]"
    meta = json.loads((out / "imp3x_meta.json").read_text())
    fails = []
    for arm in ("off", "rod"):
        ntet, attrs = _msh_census(out / f"imp3x_{arm}.msh")
        want = {1, 21, 22, 23, 24} | ({402} if arm == "rod" else set())
        ok = ntet == meta[arm]["cells"] and attrs == want
        print(f"{tag} (X0) imp3x_{arm}.msh: {ntet} tets (meta {meta[arm]['cells']}), "
              f"triangle attributes {sorted(attrs)} (want {sorted(want)}): "
              f"{'PASS' if ok else 'FAIL'}", flush=True)
        if not ok:
            fails.append(f"X0 {arm}")
    for arm, run in ARM_RUN.items():
        cfg = json.loads((out / f"{run}.json").read_text())
        pec = cfg["Boundaries"]["PEC"]["Attributes"]
        ok = pec == ARM_PEC[arm] and cfg["Model"]["Mesh"] == f"imp3x_{ARM_MESH[arm]}.msh"
        print(f"{tag} (X0) {run}.json PEC {pec} (want {ARM_PEC[arm]}), mesh "
              f"{cfg['Model']['Mesh']}: {'PASS' if ok else 'FAIL'}", flush=True)
        if not ok:
            fails.append(f"X0 cfg {arm}")
    assert not fails, f"(X0) red: {fails}"
    targets = [out / f"imp3{k}_{a}.json" for k in "pe" for a in ARM_RUN]
    for t in targets:
        assert not t.exists(), f"{t} exists: never overwrite"
    for arm, run in ARM_RUN.items():
        base = json.loads((out / f"{run}.json").read_text())
        energy = [{"Index": 1, "Attributes": [PHANTOM_CELL_TAG]}]
        probes = [{"Index": k + 1, "Center": list(p)} for k, p in enumerate(PROBE_POINTS)]
        for kind in "pe":
            cfg = copy.deepcopy(base)
            post = {"Energy": energy}
            if kind == "p":
                post["Probe"] = probes
            cfg["Domains"]["Postprocessing"] = post
            cfg["Problem"]["Output"] = f"postpro_imp3{kind}_{arm}"
            path = out / f"imp3{kind}_{arm}.json"
            with open(path, "x") as fh:
                fh.write(json.dumps(cfg, indent=2) + "\n")
            print(f"{tag} wrote logs/imp3/{path.name}: Output {cfg['Problem']['Output']}, "
                  f"PEC {cfg['Boundaries']['PEC']['Attributes']}, Postprocessing "
                  f"{sorted(post)}", flush=True)
    print(json.dumps(json.loads((out / "imp3p_rod.json").read_text())["Domains"], indent=2),
          flush=True)
    print(f"{tag} ANCHORS PASS: (X0)", flush=True)


# --------------------------------------------------------------------------- read

def _rows(path):
    with open(path) as fh:
        rows = list(csv.reader(fh))
    header = [h.strip() for h in rows[0]]
    data = [float(v) for v in rows[1]]
    return header, data


def parse_probe(path):
    """(6, 3, 4) complex E from probe-E.csv, columns mapped by header.

    Returns (array, header).  Raises on any column that does not map, or a map
    that is not a bijection onto 6 probes x 3 components x Re/Im x 4 exc.
    """
    header, data = _rows(path)
    # Palace v0.18.1 form (W2 header): ``Re{E_x[<probe>][<excitation>]} (V/m)``.
    re_c = re.compile(r"^(Re|Im)\{E_([xyz])\[(\d+)\]\[(\d+)\]\} \(V/m\)$")
    seen = {}
    for h, v in zip(header, data):
        if h.startswith("f (GHz)"):
            continue
        m = re_c.match(h)
        if not m:
            raise ValueError(f"(H) unparsed probe column {h!r}; header {header!r}")
        part, comp, k, exc = m.group(1), m.group(2), int(m.group(3)), int(m.group(4))
        key = (k, "xyz".index(comp), part, exc)
        if key in seen:
            raise ValueError(f"(H) duplicate probe column {h!r}; header {header!r}")
        seen[key] = v
    excs = sorted({k[3] for k in seen})
    want = {(k, c, part, e) for k in range(1, 7) for c in range(3)
            for part in ("Re", "Im") for e in excs}
    if set(seen) != want or len(excs) != N_EXC:
        raise ValueError(f"(H) probe map not a bijection onto 6x3x2x{N_EXC} "
                         f"(excitations {excs}); header {header!r}")
    e = np.zeros((6, 3, N_EXC), dtype=complex)
    for (k, c, part, ex), v in seen.items():
        j = excs.index(ex)
        e[k - 1, c, j] += v if part == "Re" else 1j * v
    return e, header


def parse_energy(path):
    """(4,) phantom electric energy from domain-E.csv (domain 1, per excitation)."""
    header, data = _rows(path)
    re_c = re.compile(r"^E_elec\[(\d+)\]\[(\d+)\] \(J\)$")
    seen = {}
    for h, v in zip(header, data):
        m = re_c.match(h)
        if m:
            key = (int(m.group(1)), int(m.group(2)))
            if key in seen:
                raise ValueError(f"(H) duplicate energy column {h!r}")
            seen[key] = v
    # Which index is the domain and which the excitation is decided from the
    # header (one Energy domain => that index set is {1}).
    a = {k[0] for k in seen}
    b = {k[1] for k in seen}
    if a == {1} and b == set(range(1, N_EXC + 1)) and len(seen) == N_EXC:
        w = np.array([seen[(1, e)] for e in range(1, N_EXC + 1)])
    elif b == {1} and a == set(range(1, N_EXC + 1)) and len(seen) == N_EXC:
        w = np.array([seen[(e, 1)] for e in range(1, N_EXC + 1)])
    else:
        raise ValueError(f"(H) energy map not a bijection onto 1 domain x {N_EXC}; "
                         f"matched {sorted(seen)}; header {header!r}")
    return w, header


def read():
    from tests.validation.test_port_birdcage_ring_matrix import REUSE_REPRODUCTION_RTOL

    s3 = _step3()
    o61 = s3._ops61()
    out = s3.OUT_DIR
    tag = "[IMP-3 4 W5]"
    target = out / FIELDS_JSON
    assert not target.exists(), f"{target} exists: never overwrite"
    fails = []
    e, w, s = {}, {}, {}
    for arm in ARM_RUN:
        d = out / f"postpro_imp3p_{arm}"
        ph, hdr = parse_probe(d / "probe-E.csv")
        print(f"{tag} (H) {arm} probe-E.csv header ({len(hdr)} cols): {hdr}", flush=True)
        ew, hdr = parse_energy(d / "domain-E.csv")
        print(f"{tag} (H) {arm} domain-E.csv header ({len(hdr)} cols): {hdr}", flush=True)
        e[arm], w[arm] = ph, ew
        name = f"imp3p_{arm}"
        s[arm], _ident, f = o61._p_check(name, d / "port-S.csv", tag)
        if arm != "rodnat":
            fails += [f"{name}: {x}" for x in f]
        else:
            print(f"{tag} imp3p_rodnat (P) printed only (control)", flush=True)
        s_old = o61._read_port_s(s3.csv_of(ARM_RUN[arm]))[0]
        rp = float(np.max(np.abs(s[arm] - s_old)) / np.max(np.abs(s_old)))
        ok = rp <= REUSE_REPRODUCTION_RTOL
        print(f"{tag} (RP) imp3p_{arm} vs item 235 imp3x_{arm} entrywise max|ΔS|/max|S| "
              f"{rp:.3e} <= {REUSE_REPRODUCTION_RTOL}: {'PASS' if ok else 'FAIL'}", flush=True)
        if not ok:
            fails.append(f"RP {arm} {rp}")
        for k, nm in enumerate(PROBE_NAMES):
            print(f"{tag} {arm} {nm} P1 E = " + ", ".join(o61._fmt(z) for z in ph[k, :, 0])
                  + f"  |E| {np.linalg.norm(ph[k, :, 0]):.9e}", flush=True)
        print(f"{tag} {arm} W_elec(phantom) per exc " + ", ".join(f"{v:.12e}" for v in ew),
              flush=True)
    mag = {a: np.linalg.norm(e[a], axis=1) for a in e}  # (6, 4)
    for lab, num, den in (("η_wall", "rod", "rodnat"), ("η_rod", "rod", "off")):
        r = mag[num] / mag[den]
        for j in range(N_EXC):
            print(f"{tag} Palace {lab} P{j + 1}: "
                  + "  ".join(f"{nm} {r[k, j]:.6f}" for k, nm in enumerate(PROBE_NAMES)),
                  flush=True)
    for j in range(N_EXC):
        print(f"{tag} Palace ΔW P{j + 1}: rod/rodnat - 1 {w['rod'][j] / w['rodnat'][j] - 1:+.6e}"
              f"  rod/off - 1 {w['rod'][j] / w['off'][j] - 1:+.6e}", flush=True)
    with open(target, "x") as fh:
        json.dump({
            "probe_names": list(PROBE_NAMES),
            "points": [list(p) for p in PROBE_POINTS],
            "E": {a: [[[[z.real, z.imag] for z in e[a][k, c, :]] for c in range(3)]
                      for k in range(6)] for a in e},
            "W_elec": {a: w[a].tolist() for a in w},
        }, fh)
    print(f"{tag} wrote logs/imp3/{FIELDS_JSON}", flush=True)
    assert not fails, f"(P)/(RP) red: {fails}"
    print(f"{tag} (H) PASS, (P) PASS on imp3p_off and imp3p_rod, (RP) PASS on all three",
          flush=True)


if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in ("configs", "read"):
        raise SystemExit("usage: imp3_step4_palace_fields.py configs|read")
    {"configs": configs, "read": read}[sys.argv[1]]()
