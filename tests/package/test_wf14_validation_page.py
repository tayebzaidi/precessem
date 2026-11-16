"""WF-14 step 2 (§9 item 375, chain U41): the generated validation page.

``scripts/docs/generate_validation_page.py`` renders ``docs/guide/validation.md``
from the `WF-14` registry: per public name, its API page link, the gating chunk
IDs its nodes' modules name, and the module-level bands each node function
references. Asserted (all exact):

(N) the page's ``###`` names == ``_surface()`` == registry keys, each once,
    each linking ``api/<Name>.md`` exactly once;
(E) per name, the page's node list == the registry's, in registry order;
(G) the committed page is byte-identical to a fresh generation, and two
    generations are identical;
(I) every chunk ID on the page is a ``PROJECT_PLAN.md`` §7 row;
(B) every rendered band name resolves to a top-level assignment in its linked
    file (counts printed);
(D) ``docrefs --docs-root docs/guide --check-links`` reads dead == 0, and the
    ``Link pass:`` page count is the API-only tree's plus exactly 1.

Negative controls (asserted, by construction):
(c1) a planted unregistered entry point -> the generator refuses, naming
     exactly that name, and writes nothing;
(c2) a node in a scratch module with no ID -> (F2) naming that node;
(c3) an override to ``NOPE-999`` -> (F3) naming it;
(c4) a one-byte edit of a tmp copy of the page -> ``--check`` exits 1;
(c5) a tmp copy of the page with one planted dead link -> dead == 1.

Not claimed: that a node is green today, bands asserted in fixtures or
helpers, sub-package surfaces, a hosted site.
"""

from __future__ import annotations

import ast
import importlib.util
import re
import shutil
from pathlib import Path

import pytest

import fem_em_solver
from tests.package.test_wf14_validation_envelope import (
    PLANTED_NAME,
    _load_registry,
    _surface,
)
from tests.unit.test_doc_reference_exit_codes import run_checker

ROOT = Path(__file__).resolve().parent.parent.parent
GUIDE = ROOT / "docs" / "guide"
PAGE = GUIDE / "validation.md"
GENERATOR = ROOT / "scripts" / "docs" / "generate_validation_page.py"
LINK_PASS_RE = re.compile(r"Link pass: (\d+) link\(s\) checked over (\d+) page\(s\), (\d+) dead\.")
SECTION_RE = re.compile(r"^### (\S+)$", re.MULTILINE)
ROW_RE = re.compile(r"^\| \[`([^`]+)`\]\(")
BAND_RE = re.compile(r"`([A-Z][A-Z0-9_]*)(?: = [^`]*| \(computed\))` \(\[[^\]]+\]\(\.\./\.\./([^)]+)\)\)")


def _load_generator():
    spec = importlib.util.spec_from_file_location("wf14_generate_validation_page", GENERATOR)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


gen = _load_generator()


# ------------------------------------------------------------ pure helpers
def parse_sections(text: str) -> dict[str, str]:
    heads = list(SECTION_RE.finditer(text))
    out: dict[str, str] = {}
    for i, m in enumerate(heads):
        end = heads[i + 1].start() if i + 1 < len(heads) else len(text)
        assert m.group(1) not in out, f"duplicate section {m.group(1)}"
        out[m.group(1)] = text[m.end() : end]
    return out


def parse_page_nodes(section: str) -> list[str]:
    return [m.group(1) for line in section.splitlines() if (m := ROW_RE.match(line))]


def _top_level_names(repo_rel: str) -> set[str]:
    tree = ast.parse((ROOT / repo_rel).read_text(encoding="utf-8"))
    out = set()
    for stmt in tree.body:
        if isinstance(stmt, ast.Assign):
            out |= {t.id for t in stmt.targets if isinstance(t, ast.Name)}
        elif isinstance(stmt, ast.AnnAssign) and isinstance(stmt.target, ast.Name):
            out.add(stmt.target.id)
    return out


@pytest.fixture(scope="module")
def page_text():
    return PAGE.read_text(encoding="utf-8")


# ------------------------------------------------------------------ anchors
def test_wf14_page_names_equal_surface(page_text):
    """(N)"""
    surface = _surface()
    secs = parse_sections(page_text)
    print(f"\nWF-14 s2 (N): sections={len(secs)} surface={len(surface)} registry={len(_load_registry())}")
    assert list(secs) == surface
    assert set(secs) == set(_load_registry())
    for n in surface:
        assert page_text.count(f"](api/{n}.md)") == 1, n


def test_wf14_page_nodes_equal_registry(page_text):
    """(E)"""
    registry = _load_registry()
    secs = parse_sections(page_text)
    mismatches = [n for n in secs if parse_page_nodes(secs[n]) != list(registry[n]["nodes"])]
    total = sum(len(parse_page_nodes(s)) for s in secs.values())
    print(f"\nWF-14 s2 (E): page node rows={total} mismatches={mismatches}")
    assert mismatches == []


def test_wf14_page_is_what_the_generator_emits(tmp_path):
    """(G)"""
    g1, g2 = tmp_path / "g1.md", tmp_path / "g2.md"
    stats = gen.generate(g1)
    gen.generate(g2)
    same12 = g1.read_bytes() == g2.read_bytes()
    same_c = g1.read_bytes() == PAGE.read_bytes()
    print(f"\nWF-14 s2 (G): gen1==gen2 {same12}; gen1==committed {same_c}; bytes={len(g1.read_bytes())}; stats={stats}")
    assert same12
    assert same_c
    assert gen.main(["--out", str(PAGE), "--check"]) == 0


def test_wf14_page_chunk_ids_are_plan_rows(page_text):
    """(I)"""
    ids = set(gen.ID_RE.findall(page_text))
    rows = gen.plan_rows()
    missing = sorted(ids - rows)
    print(f"\nWF-14 s2 (I): distinct chunk IDs on page={len(ids)} {sorted(ids)}; not a row={missing}")
    assert ids
    assert missing == []


def test_wf14_page_bands_resolve(page_text):
    """(B)"""
    secs = parse_sections(page_text)
    with_band = without = computed = 0
    unresolved = []
    n_bands = 0
    for sec in secs.values():
        for line in sec.splitlines():
            if not ROW_RE.match(line):
                continue
            bands = BAND_RE.findall(line)
            if bands:
                with_band += 1
            else:
                assert gen.NO_BAND in line, line
                without += 1
            computed += line.count("(computed)`")
            for name, rel in bands:
                n_bands += 1
                if name not in _top_level_names(rel):
                    unresolved.append((name, rel))
    print(
        f"\nWF-14 s2 (B): node rows with >=1 band={with_band}, with none={without}, "
        f"band refs={n_bands}, computed={computed}, unresolved={unresolved}"
    )
    assert n_bands > 0
    assert unresolved == []


def test_wf14_docrefs_guide_zero_dead_with_page(capsys):
    """(D)"""
    rc, counts = run_checker("--docs-root", "docs/guide", "--check-links")
    out = capsys.readouterr().out
    m = LINK_PASS_RE.search(out)
    api_only = sum(1 for p in GUIDE.rglob("*.md") if p != PAGE and not p.name.endswith("_private.md"))
    with capsys.disabled():
        print(out)
        print(f"WF-14 s2 (D): rc={rc} counts={counts} link_pass={m.group(0) if m else None} api-only pages={api_only}")
    assert m is not None
    n_links, n_pages, n_dead = (int(g) for g in m.groups())
    assert counts["dead"] == 0
    assert n_dead == 0
    assert PAGE.is_file()
    assert n_pages == api_only + 1


# ---------------------------------------------------------------- controls
def test_wf14_page_c1_planted_unregistered_name(tmp_path, monkeypatch):
    """(c1) the weekly's control."""
    monkeypatch.setattr(fem_em_solver, "__all__", [*fem_em_solver.__all__, PLANTED_NAME])
    monkeypatch.setattr(fem_em_solver, PLANTED_NAME, object(), raising=False)
    with pytest.raises(gen.GenerationError) as exc:
        gen.build()
    msg = str(exc.value)
    print(f"\nWF-14 s2 (c1): {msg}")
    assert msg == f"(F1) public name(s) with no registry entry: {[PLANTED_NAME]}"
    out = tmp_path / "v.md"
    assert gen.main(["--out", str(out)]) == 1
    assert not out.exists()


def _scratch_module(tmp_path: Path) -> str:
    p = tmp_path / "test_scratch_no_id.py"
    p.write_text('"""A scratch module whose docstring names no chunk."""\n\n\ndef test_x():\n    pass\n')
    return f"{p.as_posix()}::test_x"


def test_wf14_page_c2_c3_attribution_refusals(tmp_path):
    """(c2) (F2) and (c3) (F3)."""
    registry = {k: {"nodes": list(v["nodes"])} for k, v in _load_registry().items()}
    names = sorted(registry)
    node = _scratch_module(tmp_path)
    registry[names[0]]["nodes"].insert(0, node)

    with pytest.raises(gen.GenerationError) as e2:
        gen.build(names=names, registry=registry)
    print(f"\nWF-14 s2 (c2): {e2.value}")
    assert str(e2.value) == f"(F2) node whose module yields no chunk ID: {node}"

    overrides = {**gen.CHUNK_OVERRIDES, node.split("::", 1)[0]: ("NOPE-999",)}
    with pytest.raises(gen.GenerationError) as e3:
        gen.build(names=names, registry=registry, overrides=overrides)
    print(f"WF-14 s2 (c3): {e3.value}")
    assert str(e3.value) == f"(F3) chunk ID(s) with no §7 row: ['NOPE-999'] (from {node})"


def test_wf14_page_c4_hand_edit_is_seen(tmp_path):
    """(c4)"""
    copy = tmp_path / "validation.md"
    shutil.copyfile(PAGE, copy)
    assert gen.main(["--out", str(copy), "--check"]) == 0
    data = bytearray(copy.read_bytes())
    i = data.index(b"# Validation envelope") + 2
    data[i] = ord("v")
    copy.write_bytes(bytes(data))
    rc = gen.main(["--out", str(copy), "--check"])
    print(f"\nWF-14 s2 (c4): --check rc on one-byte edit = {rc}")
    assert rc == 1


def test_wf14_page_c5_planted_dead_link(tmp_path, capsys, page_text):
    """(c5) a tmp copy of the page, its link targets mirrored, one dead link planted."""
    guide = tmp_path / "docs" / "guide"
    (guide / "api").mkdir(parents=True)
    (tmp_path / "tests").symlink_to(ROOT / "tests", target_is_directory=True)
    for n in _surface():
        (guide / "api" / f"{n}.md").write_text(f"# {n}\n")
    (guide / "validation.md").write_text(page_text)
    rc0, counts0 = run_checker("--docs-root", str(guide), "--check-links")
    (guide / "validation.md").write_text(page_text + "\n[planted](api/NoSuchName.md)\n")
    rc1, counts1 = run_checker("--docs-root", str(guide), "--check-links")
    capsys.readouterr()
    with capsys.disabled():
        print(f"\nWF-14 s2 (c5): rc {rc0}/{rc1}; unplanted copy dead={counts0['dead']}; planted copy dead={counts1['dead']}")
    assert counts0["dead"] == 0
    assert counts1["dead"] == 1
