"""`GEO-34` step 1e (§9 item 996) — Netgen optimiser ladder on the 1d fixture.

Measurement only: not collected, asserts nothing, imported by nothing.  One
invocation = one rung ``(arm, h)`` from argv, so a SEGV inside
``optimize("Netgen")`` kills exactly one rung's process and the harness
records that rung's status.

Replicates ``import_step``'s gmsh read directly (never calls ``import_step``,
whose named-surface branch skips Netgen by design):
``importShapes(highestDimOnly=False)``, ``removeAllDuplicates``, synchronise.

Arms:
  E  the STEP face named ``PortSheet`` fragmented into the box (1d's
     embedded face, the configuration that crashed in
     ``20261114T181903Z_GEO-34-U12p-W1.log``).
  N  the same read, the named face removed before meshing (box only).
  S  the named face removed and replaced by an OCC rectangle spanning the box
     (x 0..1, y 0..1 at z 0.47), fragmented in: two volumes sharing the face.

Sizing: ``Mesh.MeshSizeMin = Mesh.MeshSizeMax = h`` and, as the importer does,
``mesh.setSize`` on every point to h.  Then ``generate(3)``, a pre-optimise
census (tets, minSICN min / p1, count minSICN <= 0, gamma min, getEmbedded),
flush, ``optimize("Netgen")``, the post census, and a final ROW block.
Quality metric: gmsh ``getElementQualities(..., "minSICN")`` and ``"gamma"``.

Usage:  python3 tests/mesh/probe_geo34_step1e_netgen_ladder.py --arm E --h 0.08
"""

import argparse
import hashlib
import time

import gmsh
import numpy as np

STEP = "tests/mesh/data/geo34_step1d_port_sheet.step"


def p(msg):
    print(f"[1e] {msg}", flush=True)


def named_faces():
    out = []
    for _, t in gmsh.model.getEntities(2):
        n = gmsh.model.getEntityName(2, t)
        if n and (n == "PortSheet" or n.split("/")[-1] == "PortSheet"):
            out.append(t)
    return out


def census(tag):
    types, etags, _ = gmsh.model.mesh.getElements(3)
    tets = []
    for ty, et in zip(types, etags):
        if ty == 4:  # 4-node tetrahedron
            tets.extend(int(e) for e in et)
    other = [int(ty) for ty in types if ty != 4]
    n = len(tets)
    if n == 0:
        p(f"{tag} census: tets=0 other_types={other}")
        return {"tets": 0, "sicn_min": float("nan"), "sicn_p1": float("nan"),
                "n_le0": 0, "gamma_min": float("nan")}
    sicn = np.asarray(gmsh.model.mesh.getElementQualities(tets, "minSICN"))
    gamma = np.asarray(gmsh.model.mesh.getElementQualities(tets, "gamma"))
    res = {"tets": n, "sicn_min": float(sicn.min()),
           "sicn_p1": float(np.percentile(sicn, 1)),
           "n_le0": int(np.count_nonzero(sicn <= 0.0)),
           "gamma_min": float(gamma.min())}
    p(f"{tag} census: tets={n} minSICN min={res['sicn_min']:.6e} "
      f"p1={res['sicn_p1']:.6e} n(minSICN<=0)={res['n_le0']} "
      f"gamma min={res['gamma_min']:.6e} other_elem_types={other}")
    for _, v in gmsh.model.getEntities(3):
        p(f"{tag} getEmbedded(3,{v})={gmsh.model.mesh.getEmbedded(3, v)}")
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", required=True, choices=["E", "N", "S"])
    ap.add_argument("--h", required=True, type=float)
    a = ap.parse_args()

    with open(STEP, "rb") as fh:
        p(f"{STEP} sha256={hashlib.sha256(fh.read()).hexdigest()}")
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    p(f"gmsh.__version__={getattr(gmsh, '__version__', '?')}")
    for opt in ("General.Version", "General.BuildInfo", "General.BuildOptions"):
        try:
            p(f"{opt}={gmsh.option.getString(opt)!r}")
        except Exception as exc:  # printed, not fatal
            p(f"{opt} unavailable: {exc!r}")
    p(f"rung arm={a.arm} h={a.h}")
    gmsh.logger.start()
    t0 = time.perf_counter()
    try:
        gmsh.model.add("geo34_1e")
        gmsh.model.occ.importShapes(STEP, highestDimOnly=False)
        gmsh.model.occ.removeAllDuplicates()
        gmsh.model.occ.synchronize()
        faces = named_faces()
        p(f"post-heal vols={gmsh.model.getEntities(3)} named PortSheet faces={faces}")
        vols = gmsh.model.getEntities(3)
        if a.arm == "E":
            out, out_map = gmsh.model.occ.fragment(vols, [(2, s) for s in faces])
            gmsh.model.occ.synchronize()
            p(f"E fragment out={out} out_map={out_map}")
        else:
            gmsh.model.occ.remove([(2, s) for s in faces], recursive=True)
            gmsh.model.occ.synchronize()
            p(f"{a.arm}: removed named faces {faces}; remaining faces="
              f"{len(gmsh.model.getEntities(2))} vols={gmsh.model.getEntities(3)}")
            if a.arm == "S":
                r = gmsh.model.occ.addRectangle(0.0, 0.0, 0.47, 1.0, 1.0)
                out, out_map = gmsh.model.occ.fragment(gmsh.model.getEntities(3), [(2, r)])
                gmsh.model.occ.synchronize()
                p(f"S fragment out={out} out_map={out_map} "
                  f"vol masses={[gmsh.model.occ.getMass(3, v) for _, v in gmsh.model.getEntities(3)]}")
        p(f"pre-mesh entities: vols={len(gmsh.model.getEntities(3))} "
          f"faces={len(gmsh.model.getEntities(2))}")

        gmsh.option.setNumber("Mesh.MeshSizeMin", a.h)
        gmsh.option.setNumber("Mesh.MeshSizeMax", a.h)
        gmsh.model.mesh.setSize(gmsh.model.getEntities(0), a.h)
        tg = time.perf_counter()
        gmsh.model.mesh.generate(3)
        t_gen = time.perf_counter() - tg
        p(f"generate(3) wall={t_gen:.3f} s")
        pre = census("PRE")
        p("calling optimize('Netgen') now")
        to = time.perf_counter()
        gmsh.model.mesh.optimize("Netgen")
        t_opt = time.perf_counter() - to
        p(f"optimize('Netgen') returned wall={t_opt:.3f} s")
        post = census("POST")
        illegal = [m for m in gmsh.logger.get() if "illegal" in m.lower()]
        p(f"logger lines mentioning 'illegal': {illegal}")
        emb = any(gmsh.model.mesh.getEmbedded(3, v) for _, v in gmsh.model.getEntities(3))
        p("ROW arm | h | tets pre | minSICN min pre | n<=0 pre | embedded | "
          "illegal lines | crashed | tets post | minSICN min post | elapsed s")
        p(f"ROW {a.arm} | {a.h} | {pre['tets']} | {pre['sicn_min']:.6e} | "
          f"{pre['n_le0']} | {'Y' if emb else 'N'} | {illegal or '-'} | N | "
          f"{post['tets']} | {post['sicn_min']:.6e} | {time.perf_counter() - t0:.3f}")
    finally:
        gmsh.logger.stop()
        gmsh.finalize()


if __name__ == "__main__":
    main()
