"""`TH-16` step 3a probe: build the quadrant birdcage once, print the outcome.

Measurement only — no assertion.  Run under the harness at ``-n 1``.
"""

from tests.mesh.test_birdcage_port_sheets import _build

if __name__ == "__main__":
    try:
        mesh, cells, facets, diag, elapsed = _build(True, symmetry_quadrant=True)
        print("[probe] cells", mesh.topology.index_map(3).size_global, "elapsed", elapsed)
    except Exception as exc:  # noqa: BLE001
        print("[probe] raised:", exc)
