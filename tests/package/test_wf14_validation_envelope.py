"""WF-14 step 1 (§9 item 221): the validation envelope.

A committed registry (``validation_envelope.toml``) maps every public name of
``fem_em_solver.__all__`` to the gate test nodes that cover it. Asserted:

(S) the surface equals the committed snapshot (``public_surface.txt``);
(R) the registry keys equal the surface, no entry has empty ``nodes``;
(C) every registered node is among the IDs ``pytest --collect-only`` finds,
    with collection ``rc == 0``;
(N) the registry restates no number (node IDs only).

Negative controls (planted, asserted by construction):
(c1) a planted unregistered public name is reported exactly, both as not in
     the snapshot and as having no registry entry (pure checker and the real
     ``__all__`` path via monkeypatch);
(c2) a planted uncollectable node yields exactly that violation against the
     same collected IDs.

Not claimed: that a registered node exercises its name or is green today
(collectability only), the sub-package surfaces, docs rendering.

Collection runs in ONE subprocess (never a nested ``pytest.main`` in this
session; never parsing verbosity-dependent text), with cwd /workspace so the
pyproject ini (``--import-mode=importlib``) and rootdir apply, and with MPI
launcher variables stripped from the child environment.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tomllib
from pathlib import Path

import pytest

import fem_em_solver

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
REGISTRY_PATH = HERE / "validation_envelope.toml"
SNAPSHOT_PATH = HERE / "public_surface.txt"
PLANTED_NAME = "_wf14_planted_name"
PLANTED_NODE = "tests/validation/test_lossy_plane_wave.py::test_wf14_no_such_node"
_STRIP_PREFIXES = ("PMI_", "HYDRA_", "MPIR_", "OMPI_", "PMIX_")

_COLLECT_SCRIPT = r"""
import json, sys
import pytest

class C:
    def __init__(self):
        self.ids = []
    def pytest_collection_finish(self, session):
        self.ids = [i.nodeid for i in session.items]

c = C()
rc = pytest.main(["--collect-only", "-p", "no:cacheprovider", *sys.argv[2:]], plugins=[c])
with open(sys.argv[1], "w") as fh:
    json.dump({"rc": int(rc), "ids": c.ids}, fh)
"""


# ----------------------------------------------------------------- loaders
def _load_registry() -> dict[str, dict]:
    with REGISTRY_PATH.open("rb") as fh:
        data = tomllib.load(fh)
    return dict(data.get("names", {}))


def _load_snapshot() -> list[str]:
    lines = SNAPSHOT_PATH.read_text().splitlines()
    return [ln.strip() for ln in lines if ln.strip()]


def _surface() -> list[str]:
    names = sorted(set(fem_em_solver.__all__))
    for n in names:
        getattr(fem_em_solver, n)  # a dangling __all__ entry raises here
    return names


# ------------------------------------------------------------ pure checkers
def _surface_violations(surface, snapshot, registry) -> list[str]:
    surface_set, snapshot_set = set(surface), set(snapshot)
    out: list[str] = []
    for n in sorted(surface_set - snapshot_set):
        out.append(f"not in snapshot: {n}")
    for n in sorted(snapshot_set - surface_set):
        out.append(f"in snapshot but gone: {n}")
    for n in sorted(surface_set - set(registry)):
        out.append(f"no registry entry: {n}")
    for n in sorted(set(registry) - surface_set):
        out.append(f"registry entry for no surface name: {n}")
    for n in sorted(registry):
        if not list(registry[n].get("nodes", [])):
            out.append(f"empty nodes: {n}")
    return out


def _node_matches(node: str, collected: set[str], prefixes: tuple[str, ...]) -> bool:
    if node in collected:
        return True
    if "[" not in node:
        return any(c.startswith(node + "[") for c in prefixes)
    return False


def _node_violations(registry, collected_ids) -> list[str]:
    collected = set(collected_ids)
    parametrised = tuple(c for c in collected if "[" in c)
    out: list[str] = []
    for name in sorted(registry):
        for node in registry[name].get("nodes", []):
            if not _node_matches(node, collected, parametrised):
                out.append(f"uncollectable node for {name}: {node}")
    return out


def _registry_files(registry) -> list[str]:
    return sorted({node.split("::", 1)[0] for e in registry.values() for node in e.get("nodes", [])})


# ---------------------------------------------------------------- fixtures
@pytest.fixture(scope="module")
def registry():
    return _load_registry()


@pytest.fixture(scope="module")
def snapshot():
    return _load_snapshot()


@pytest.fixture(scope="module")
def collection(registry, tmp_path_factory):
    out = tmp_path_factory.mktemp("wf14") / "collected.json"
    files = _registry_files(registry)
    env = {k: v for k, v in os.environ.items() if not k.startswith(_STRIP_PREFIXES)}
    stripped = sorted(k for k in os.environ if k.startswith(_STRIP_PREFIXES))
    proc = subprocess.run(
        [sys.executable, "-c", _COLLECT_SCRIPT, str(out), *files],
        cwd=str(ROOT),
        env=env,
        capture_output=True,
        text=True,
        timeout=150,
    )
    result = json.loads(out.read_text()) if out.exists() else {"rc": None, "ids": []}
    result["files"] = files
    result["stdout_tail"] = proc.stdout[-4000:]
    result["stderr_tail"] = proc.stderr[-4000:]
    print(
        f"\nWF-14 collection: rc={result['rc']} files={len(files)} "
        f"collected_ids={len(result['ids'])} env_stripped={stripped}"
    )
    return result


# ------------------------------------------------------------------ anchors
def test_surface_matches_snapshot(snapshot):
    """(S)"""
    surface = _surface()
    print(f"\nWF-14 (S): surface={len(surface)} snapshot={len(snapshot)}")
    assert surface == snapshot


def test_every_public_name_has_an_entry(registry, snapshot):
    """(R)"""
    surface = _surface()
    v = _surface_violations(surface, snapshot, registry)
    print(f"\nWF-14 (R): registry keys={len(registry)} violations={v}")
    assert v == []
    assert sorted(registry) == surface


def test_every_entry_names_a_collectable_node(registry, collection):
    """(C)"""
    if collection["rc"] != 0:
        print(collection["stdout_tail"])
        print(collection["stderr_tail"])
    assert collection["rc"] == 0
    v = _node_violations(registry, collection["ids"])
    n_nodes = sum(len(e["nodes"]) for e in registry.values())
    print(f"\nWF-14 (C): registered node refs={n_nodes} violations={v}")
    assert v == []

    # printed, never asserted: does any named node's file mention the name?
    mentioned = 0
    for name in sorted(registry):
        files = {node.split("::", 1)[0] for node in registry[name]["nodes"]}
        hit = any(name in (ROOT / f).read_text() for f in files)
        mentioned += hit
        print(f"WF-14 textual mention: {name}: {'yes' if hit else 'no'}")
    print(f"WF-14 textual mention tally: {mentioned} / {len(registry)}")


def test_registry_restates_no_number(registry):
    """(N)"""
    for name, entry in registry.items():
        assert set(entry) == {"nodes"}, (name, sorted(entry))
        assert isinstance(entry["nodes"], list), name
        for node in entry["nodes"]:
            assert isinstance(node, str), (name, node)
            try:
                float(node)
            except ValueError:
                pass
            else:
                pytest.fail(f"{name}: node {node!r} parses as a number")


# ---------------------------------------------------------------- controls
def test_c1_planted_unregistered_name_is_reported(registry, snapshot, monkeypatch):
    """(c1) pure checker, then the real __all__ path."""
    expected = [f"not in snapshot: {PLANTED_NAME}", f"no registry entry: {PLANTED_NAME}"]
    pure = _surface_violations(_surface() + [PLANTED_NAME], snapshot, registry)
    print(f"\nWF-14 (c1) pure: {pure}")
    assert pure == expected

    monkeypatch.setattr(fem_em_solver, "__all__", [*fem_em_solver.__all__, PLANTED_NAME])
    monkeypatch.setattr(fem_em_solver, PLANTED_NAME, object(), raising=False)
    surface = _surface()
    assert surface != snapshot  # the snapshot check fails on the planted name
    real = _surface_violations(surface, snapshot, registry)
    print(f"WF-14 (c1) real __all__ path: {real}")
    assert real == expected


def test_c2_planted_uncollectable_node_is_reported(registry, collection):
    """(c2) same collected IDs, one planted node."""
    assert collection["rc"] == 0
    planted = {k: {"nodes": list(v["nodes"])} for k, v in registry.items()}
    planted["TimeHarmonicSolver"]["nodes"].append(PLANTED_NODE)
    v = _node_violations(planted, collection["ids"])
    print(f"\nWF-14 (c2): {v}")
    assert v == [f"uncollectable node for TimeHarmonicSolver: {PLANTED_NODE}"]
