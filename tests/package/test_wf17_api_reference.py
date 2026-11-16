"""WF-17 step 2 (§9 item 274, chain U15): the API reference with the envelope.

``scripts/docs/generate_api_reference.py`` writes one page per name in
``fem_em_solver.__all__`` under ``docs/guide/api/`` plus ``index.md``; each
page renders its WF-14 envelope entry between ``<!-- envelope:begin -->`` /
``<!-- envelope:end -->``. Asserted (all exact):

(N) one page per public name: page stems == surface == snapshot; the index
    links each name exactly once; no case-insensitive collision;
(G) the committed pages are byte-identical to ``generate(tmp)``, and a second
    generation is identical to the first (determinism);
(E) the node list parsed from each committed page equals the registry's, in
    registry order;
(C) every node parsed from the pages collects (collect-only subprocess, WF-14
    helpers imported by name), rc == 0, no violations;
(D) ``docrefs --docs-root docs/guide --check-links`` reads dead == 0 with the
    pages included; the ``Link pass:`` line reads >= 36 pages, >= 102 links.

Negative controls (asserted, by construction; OPS-19 precedent):
(c1) a removed ``Study.md`` -> exactly ["no page: Study"];
(c2) a synthetic tmp tree with one dead link and one fenced link-shaped
     string reads dead == 1 exactly;
(c3) a planted uncollectable node in Study's parsed list -> exactly that one
     violation against the same collected IDs;
(c4) a hand-edited ``Study.md`` -> the tree comparator reports ["Study.md"].

Not claimed: a built/hosted site, docstring quality, sub-package surfaces,
that any node is green today, anchor-fragment resolution.
"""

from __future__ import annotations

import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from tests.package.test_wf14_validation_envelope import (
    _COLLECT_SCRIPT,
    _STRIP_PREFIXES,
    _load_registry,
    _load_snapshot,
    _node_violations,
    _surface,
)
from tests.unit.test_doc_reference_exit_codes import run_checker

ROOT = Path(__file__).resolve().parent.parent.parent
API_DIR = ROOT / "docs" / "guide" / "api"
GENERATOR = ROOT / "scripts" / "docs" / "generate_api_reference.py"
PLANTED_NODE = "tests/package/test_no_such_module.py::test_nothing"
LINK_PASS_RE = re.compile(r"Link pass: (\d+) link\(s\) checked over (\d+) page\(s\), (\d+) dead\.")
NODE_LINE_RE = re.compile(r"^- `([^`]+)`")


def _load_generator():
    spec = importlib.util.spec_from_file_location("wf17_generate_api_reference", GENERATOR)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


gen = _load_generator()
render = gen.render
generate = gen.generate


# ------------------------------------------------------------ pure helpers
def parse_page_nodes(text: str) -> list[str]:
    """Node IDs between the envelope markers of one committed page."""
    begin = text.index(gen.ENVELOPE_BEGIN)
    end = text.index(gen.ENVELOPE_END, begin)
    block = text[begin + len(gen.ENVELOPE_BEGIN) : end]
    out = []
    for line in block.splitlines():
        m = NODE_LINE_RE.match(line.strip())
        if m:
            out.append(m.group(1))
    return out


def page_stems(api_dir: Path) -> set[str]:
    return {p.stem for p in Path(api_dir).glob("*.md")} - {"index"}


def page_violations(stems: set[str], surface: list[str]) -> list[str]:
    out = [f"no page: {n}" for n in sorted(set(surface) - stems)]
    out += [f"page for no public name: {s}" for s in sorted(stems - set(surface))]
    return out


def page_registry(api_dir: Path, names) -> dict[str, dict]:
    return {n: {"nodes": parse_page_nodes((Path(api_dir) / f"{n}.md").read_text())} for n in names}


def _sha_tree(api_dir: Path) -> str:
    import hashlib

    h = hashlib.sha256()
    for p in sorted(Path(api_dir).glob("*.md")):
        h.update(p.name.encode())
        h.update(p.read_bytes())
    return h.hexdigest()


@pytest.fixture(scope="module")
def page_reg():
    return page_registry(API_DIR, _surface())


@pytest.fixture(scope="module")
def page_collection(page_reg, tmp_path_factory):
    out = tmp_path_factory.mktemp("wf17api") / "collected.json"
    files = sorted({n.split("::", 1)[0] for e in page_reg.values() for n in e["nodes"]})
    env = {k: v for k, v in os.environ.items() if not k.startswith(_STRIP_PREFIXES)}
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
        f"\nWF-17 api collection: rc={result['rc']} files={len(files)} "
        f"collected_ids={len(result['ids'])}"
    )
    return result


# ------------------------------------------------------------------ anchors
def test_wf17_every_public_name_has_a_page():
    """(N)"""
    surface = _surface()
    stems = page_stems(API_DIR)
    print(f"\nWF-17 (N): pages={len(stems)} (+ index) surface={len(surface)} tree sha256={_sha_tree(API_DIR)}")
    assert page_violations(stems, surface) == []
    assert stems == set(surface) == set(_load_snapshot())
    lowered = [s.lower() for s in stems]
    assert len(lowered) == len(set(lowered)), "case-insensitive page collision"
    index = (API_DIR / "index.md").read_text()
    for n in surface:
        assert index.count(f"]({n}.md)") == 1, n
    no_doc = [n for n in surface if "*(no docstring)*" in (API_DIR / f"{n}.md").read_text()]
    print(f"WF-17 (N) printed: names with no docstring = {len(no_doc)} {no_doc}")


def test_wf17_pages_are_what_the_generator_emits(tmp_path):
    """(G)"""
    t1, t2 = tmp_path / "g1", tmp_path / "g2"
    generate(t1)
    generate(t2)
    d12 = gen.diff_trees(t1, t2)
    dc = gen.diff_trees(t1, API_DIR)
    print(f"\nWF-17 (G): gen1 vs gen2 diff={d12}; gen1 vs committed diff={dc}; files={len(list(t1.glob('*.md')))}")
    assert d12 == []
    assert dc == [], f"first differing file: {dc[0] if dc else None}"


def test_wf17_page_envelope_equals_registry(page_reg):
    """(E)"""
    registry = _load_registry()
    mismatches = [n for n in sorted(page_reg) if page_reg[n]["nodes"] != list(registry[n]["nodes"])]
    total = sum(len(e["nodes"]) for e in page_reg.values())
    print(f"\nWF-17 (E): page node refs={total} mismatches={mismatches}")
    assert sorted(page_reg) == sorted(registry)
    assert mismatches == []


def test_wf17_every_page_node_collects(page_reg, page_collection):
    """(C)"""
    if page_collection["rc"] != 0:
        print(page_collection["stdout_tail"])
        print(page_collection["stderr_tail"])
    print(f"\nWF-17 (C): files={page_collection['files']}")
    print(f"WF-17 (C): collected ids={len(page_collection['ids'])}")
    assert page_collection["rc"] == 0
    v = _node_violations(page_reg, page_collection["ids"])
    print(f"WF-17 (C): violations={v}")
    assert v == []


def test_wf17_docrefs_guide_zero_dead_with_api_pages(capsys):
    """(D)"""
    rc, counts = run_checker("--docs-root", "docs/guide", "--check-links")
    out = capsys.readouterr().out
    m = LINK_PASS_RE.search(out)
    with capsys.disabled():
        print(out)
        print(f"WF-17 (D): rc={rc} counts={counts} link_pass={m.group(0) if m else None}")
    assert m is not None
    n_links, n_pages, n_dead = (int(g) for g in m.groups())
    assert counts["dead"] == 0
    assert n_dead == 0
    n_names = len(_surface())
    assert n_pages >= n_names + 2
    assert n_links >= n_names * 3


# ---------------------------------------------------------------- controls
def test_wf17_api_planted_controls(tmp_path, page_reg, page_collection, capsys):
    """(c1)-(c4), asserted by construction."""
    surface = _surface()

    # (c1) a removed page
    c1 = tmp_path / "c1"
    shutil.copytree(API_DIR, c1)
    (c1 / "Study.md").unlink()
    v1 = page_violations(page_stems(c1), surface)
    print(f"\nWF-17 (c1): {v1}")
    assert v1 == ["no page: Study"]

    # (c2) synthetic tree: one live, one dead, one fenced link-shaped string
    c2 = tmp_path / "c2" / "api"
    c2.mkdir(parents=True)
    (c2 / "index.md").write_text("# index\n")
    (c2 / "Page.md").write_text(
        "# Page\n\n[back](index.md)\n\n[x](NoSuchName.md)\n\n```text\n[y](also_missing.md)\n```\n"
    )
    rc2, counts2 = run_checker("--docs-root", str(c2), "--check-links")
    out2 = capsys.readouterr().out
    with capsys.disabled():
        print(out2)
        print(f"WF-17 (c2): rc={rc2} dead={counts2['dead']}")
    assert counts2["dead"] == 1

    # (c3) an uncollectable node planted in one page's parsed list
    assert page_collection["rc"] == 0
    planted = {k: {"nodes": list(v["nodes"])} for k, v in page_reg.items()}
    planted["Study"]["nodes"].append(PLANTED_NODE)
    v3 = _node_violations(planted, page_collection["ids"])
    print(f"WF-17 (c3): {v3}")
    assert v3 == [f"uncollectable node for Study: {PLANTED_NODE}"]

    # (c4) a hand-edited page
    c4 = tmp_path / "c4"
    shutil.copytree(API_DIR, c4)
    text = (c4 / "Study.md").read_text()
    head, sep, tail = text.partition("## Docstring\n")
    assert sep
    fence_at = tail.index("\n```", tail.index("```") + 3)
    (c4 / "Study.md").write_text(head + sep + tail[:fence_at] + "X" + tail[fence_at:])
    g = tmp_path / "c4gen"
    generate(g)
    v4 = gen.diff_trees(c4, g)
    print(f"WF-17 (c4): {v4}")
    assert v4 == ["Study.md"]
