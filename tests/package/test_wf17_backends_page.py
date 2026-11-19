"""`ARCH-1` step 1c (§9 item 529, §10 U62): the "backends" guide page, tested without a solve.

The page (`docs/guide/backends.md`) documents the opt-in Palace backend of
`ARCH-1` step 1 (U16): its three stages, its refusals, the binding
matched-order cross-check condition and the U16 nodes that hold it.  Nothing
here meshes or solves (the solves are `tests/package/test_arch1_palace_backend.py`):
one block per tag and one ``solve`` call with ``backend="palace"`` (X); the
Palace command block equals ``PalaceRunRequired(...).argv`` built from
``backends.palace``'s default rank count, and the tolerance line equals
``LEG_D0_REPRODUCTION_BAND`` formatted at ``.0e`` (E); every cited node is
registered under ``Study`` / ``run_study_file`` and the validation / API links
resolve (L); the guide tree reads 0 dead links (D).  Three negative controls
(c1)-(c3) show each comparator can fail.  Brief: docs/planning/chunks/ARCH-1.md.
"""

from __future__ import annotations

import ast
import inspect
import re
import shutil
import tomllib
from pathlib import Path

from fem_em_solver.backends.palace import PalaceRunRequired, solve_palace
from tests.package.test_wf17_first_study_tutorial import extract_tagged_block
from tests.unit.test_doc_reference_exit_codes import run_checker
from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND

ROOT = Path(__file__).resolve().parents[2]
GUIDE = ROOT / "docs" / "guide"
PAGE = GUIDE / "backends.md"
REGISTRY = ROOT / "tests" / "package" / "validation_envelope.toml"
REGISTRY_NAMES = ("Study", "run_study_file")
TAGS = (
    ("python", "palace-backend"),
    ("bash", "palace-backend-cli"),
    ("text", "palace-command"),
    ("text", "palace-tolerance"),
)
NODE_RE = re.compile(r"`(tests/[A-Za-z0-9_/.]+\.py::[A-Za-z0-9_]+)`")
LINK_RE = re.compile(r"\]\(([^)\s]+)\)")
TOL_RE = re.compile(r"^Frobenius-relative S tolerance: (\S+)$")
LINK_PASS_RE = re.compile(r"Link pass: (\d+) link\(s\) checked over (\d+) page\(s\), (\d+) dead\.")


def _page() -> str:
    return PAGE.read_text()


def _differing_lines(a: str, b: str) -> list[int]:
    la, lb = a.splitlines(), b.splitlines()
    n = max(len(la), len(lb))
    return [
        i for i in range(n) if (la[i] if i < len(la) else None) != (lb[i] if i < len(lb) else None)
    ]


def _expected_command() -> str:
    ranks = inspect.signature(solve_palace).parameters["palace_ranks"].default
    req = PalaceRunRequired(["palace", "-np", str(int(ranks)), "case.json"], ROOT)
    return " ".join(req.argv) + "\n"


def _registered_nodes() -> set[str]:
    reg = tomllib.loads(REGISTRY.read_text())["names"]
    return {node for name in REGISTRY_NAMES for node in reg[name]["nodes"]}


def _unregistered(nodes, registered) -> list[str]:
    return [n for n in nodes if n not in registered]


def test_arch1_u62_x_blocks_and_python_call():
    """(X) one block per tag; one .solve(backend="palace", case_dir=...) call; CLI run twice."""
    text = _page()
    blocks = {tag: extract_tagged_block(text, tag, lang) for lang, tag in TAGS}
    calls = [
        n
        for n in ast.walk(ast.parse(blocks["palace-backend"]))
        if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute) and n.func.attr == "solve"
    ]
    print(f"\n[ARCH-1 U62] (X) tags={len(blocks)}, solve calls={len(calls)}")
    assert len(calls) == 1
    kw = {k.arg: k.value for k in calls[0].keywords}
    assert isinstance(kw.get("backend"), ast.Constant) and kw["backend"].value == "palace"
    assert "case_dir" in kw
    cli = blocks["palace-backend-cli"]
    runs = [line for line in cli.splitlines() if line.startswith("femem run ")]
    print(
        f"[ARCH-1 U62] (X) femem run lines={len(runs)}; with --backend palace={sum('--backend palace' in r for r in runs)}"
    )
    assert len(runs) == 2
    assert all("--backend palace" in r for r in runs)
    assert cli.count("--backend palace") == 2


def test_arch1_u62_e_command_and_tolerance_rendered_and_c1():
    """(E) palace-command == PalaceRunRequired(...).argv; tolerance == LEG_D0 at .0e; (c1) -np 4 differs on one line."""
    text = _page()
    block = extract_tagged_block(text, "palace-command", "text")
    expected = _expected_command()
    print(f"\n[ARCH-1 U62] (E) page command={block.strip()!r} expected={expected.strip()!r}")
    assert block == expected

    tol_block = extract_tagged_block(text, "palace-tolerance", "text")
    lines = [line for line in tol_block.splitlines() if line.strip()]
    assert len(lines) == 1
    m = TOL_RE.match(lines[0])
    assert m is not None
    tol = f"{LEG_D0_REPRODUCTION_BAND:.0e}"
    print(f"[ARCH-1 U62] (E) page tolerance={m.group(1)} LEG_D0_REPRODUCTION_BAND -> {tol}")
    assert m.group(1) == tol

    mutated = block.replace("-np 8", "-np 4")
    diff = _differing_lines(mutated, expected)
    print(f"[ARCH-1 U62] (c1) -np 4 copy differs on lines {diff}")
    assert mutated != expected
    assert len(diff) == 1


def test_arch1_u62_l_nodes_registered_links_resolve_and_c2():
    """(L) cited nodes registered under Study / run_study_file; links resolve; (c2) a planted node is named."""
    text = _page()
    nodes = NODE_RE.findall(text)
    registered = _registered_nodes()
    missing = _unregistered(nodes, registered)
    print(
        f"\n[ARCH-1 U62] (L) nodes cited={len(nodes)} registered under {REGISTRY_NAMES}; unregistered={missing}"
    )
    assert len(nodes) == 6
    assert len(set(nodes)) == len(nodes)
    assert missing == []

    planted = "tests/package/test_arch1_palace_backend.py::test_u16_no_such_node"
    caught = _unregistered(nodes + [planted], registered)
    print(f"[ARCH-1 U62] (c2) planted node -> unregistered={caught}")
    assert caught == [planted]

    links = LINK_RE.findall(text)
    val_heads = {
        line.lstrip("#").strip().lower()
        for line in (GUIDE / "validation.md").read_text().splitlines()
        if line.startswith("### ")
    }
    required = {"validation.md", "api/Study.md", "api/run_study_file.md"}
    targets = {lk.split("#")[0] for lk in links}
    print(f"[ARCH-1 U62] (L) links={len(links)}; required present={required <= targets}")
    assert required <= targets
    for lk in links:
        if lk.startswith("http"):
            continue
        path, _, anchor = lk.partition("#")
        assert (PAGE.parent / path).resolve().exists(), lk
        if path == "validation.md" and anchor:
            assert anchor in val_heads, lk
    assert "validation.md#study" in links and "validation.md#run_study_file" in links
    assert (
        "../../examples/ansys_benchmarks/ans4_birdcage_four_port_10_64_128MHz/COMPARISON.md"
        in links
    )


def test_arch1_u62_d_guide_zero_dead(capsys):
    """(D) run_checker --docs-root docs/guide --check-links reads dead == 0."""
    rc, counts = run_checker("--docs-root", "docs/guide", "--check-links")
    out = capsys.readouterr().out
    m = LINK_PASS_RE.search(out)
    with capsys.disabled():
        print(
            f"\n[ARCH-1 U62] (D) rc={rc} {m.group(0) if m else 'no link-pass line'} counts={counts}"
        )
    assert m is not None
    assert int(m.group(3)) == 0
    assert counts["dead"] == 0


def test_arch1_u62_c3_planted_dead_link(tmp_path, capsys):
    """(c3) a tmp copy of the guide tree with one planted dead link reads dead == 1."""
    guide = tmp_path / "docs" / "guide"
    shutil.copytree(GUIDE, guide)
    for top in ("tests", "examples", "src", "docs/planning", "scripts"):
        src = ROOT / top
        dst = tmp_path / top
        if src.exists() and not dst.exists():
            dst.parent.mkdir(parents=True, exist_ok=True)
            dst.symlink_to(src, target_is_directory=True)
    rc0, counts0 = run_checker("--docs-root", str(guide), "--check-links")
    page = guide / PAGE.name
    page.write_text(page.read_text() + "\n[planted](api/NoSuchName.md)\n")
    rc1, counts1 = run_checker("--docs-root", str(guide), "--check-links")
    capsys.readouterr()
    with capsys.disabled():
        print(
            f"\n[ARCH-1 U62] (c3) rc {rc0}/{rc1}; copy dead={counts0['dead']}; planted copy dead={counts1['dead']}"
        )
    assert counts0["dead"] == 0
    assert counts1["dead"] == 1
