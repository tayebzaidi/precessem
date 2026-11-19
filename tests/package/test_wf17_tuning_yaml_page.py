"""`WF-17` step 6: the "Tuning from a YAML file" guide page, tested without a solve.

The page (`docs/guide/tuning_from_yaml.md`) documents the `WF-1` `tune:`
section and the `WF-16` ``Study.tune`` read it drives.  Nothing here meshes or
solves (the solve is ``test_wf1_yaml_tune`` and ``mri:21``'s run): the page's
YAML block must equal the tracked file byte for byte (Y=), its key list must
equal ``_tune_summary``'s keys on a stub (S), its printed ``C_mode`` must be
the TH-17 record times item 429's measured (a), formatted at ``.4e`` from
imported constants (E), every test node it names must be registered in
``validation_envelope.toml`` under the Study / YAML names (L), and the guide
tree must read 0 dead links (D).  Four negative controls (c1)-(c4) show each
comparator can fail.
"""

from __future__ import annotations

import ast
import hashlib
import re
import shutil
import tomllib
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from mpi4py import MPI

from fem_em_solver.study import TuneResult
from fem_em_solver.study_yaml import _tune_summary
from tests.package.test_wf1_yaml_tune import ITEM429_A_PRINT, TUNE_FILE
from tests.package.test_wf16_study_tune import TH17_STEP3_C_MODE_RECORD_F
from tests.package.test_wf17_first_study_tutorial import extract_tagged_block
from tests.unit.test_doc_reference_exit_codes import run_checker
from tests.validation.test_lossy_sphere_fullwave import FREQUENCY_64_HZ

ROOT = Path(__file__).resolve().parents[2]
GUIDE = ROOT / "docs" / "guide"
PAGE = GUIDE / "tuning_from_yaml.md"
REGISTRY = ROOT / "tests" / "package" / "validation_envelope.toml"
REGISTRY_NAMES = ("Study", "load_study", "run_study_file", "StudyConfigError")
TAGS = (
    ("yaml", "f-small-hole-pec-tune"),
    ("python", "tune-yaml"),
    ("bash", "tune-yaml-cli"),
    ("text", "tune-yaml-output"),
    ("text", "tune-summary-keys"),
)
# U49's own YAML reading (20261014T110326Z_WF-1-s11-W1.log:1846), printed beside, never asserted.
U49_YAML_C_MODE_PRINT = 6.626796421e-11
# Item 429's k = 2 offset (+48.95 %, docs/planning/chunks/WF-16.md:984), for (c4) only.
ITEM429_K2_FACTOR = 1.4895
NODE_RE = re.compile(r"`(tests/[A-Za-z0-9_/.]+\.py::[A-Za-z0-9_]+)`")
LINK_RE = re.compile(r"\]\(([^)\s]+)\)")
OUTPUT_RE = re.compile(r"^C_mode = (\S+) F$")
LINK_PASS_RE = re.compile(r"Link pass: (\d+) link\(s\) checked over (\d+) page\(s\), (\d+) dead\.")


def _page() -> str:
    return PAGE.read_text()


def _differing_lines(a: str, b: str) -> list[int]:
    la, lb = a.splitlines(), b.splitlines()
    n = max(len(la), len(lb))
    return [
        i for i in range(n) if (la[i] if i < len(la) else None) != (lb[i] if i < len(lb) else None)
    ]


def _expected_output() -> str:
    return f"{TH17_STEP3_C_MODE_RECORD_F * (1 + ITEM429_A_PRINT):.4e}"


def _registered_nodes() -> set[str]:
    reg = tomllib.loads(REGISTRY.read_text())["names"]
    return {node for name in REGISTRY_NAMES for node in reg[name]["nodes"]}


def _unregistered(nodes, registered) -> list[str]:
    return [n for n in nodes if n not in registered]


def test_wf17_s6_x_blocks_and_python_call():
    """(X) one block per tag; one run_study_file call on the YAML's name; no Study(."""
    text = _page()
    blocks = {tag: extract_tagged_block(text, tag, lang) for lang, tag in TAGS}
    src = blocks["tune-yaml"]
    calls = [
        n
        for n in ast.walk(ast.parse(src))
        if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)
    ]
    rsf = [c for c in calls if c.func.id == "run_study_file"]
    studies = [c for c in calls if c.func.id == "Study"]
    print(
        f"\n[WF-17 s6] (X) tags={len(blocks)}, run_study_file calls={len(rsf)}, Study calls={len(studies)}"
    )
    assert len(rsf) == 1
    assert isinstance(rsf[0].args[0], ast.Constant) and rsf[0].args[0].value == TUNE_FILE.name
    assert studies == [] and "Study(" not in src
    assert f"femem run {TUNE_FILE.name}" in blocks["tune-yaml-cli"]


def test_wf17_s6_y_yaml_block_byte_equal_and_c1():
    """(Y=) and (c1): k: 1 -> k: 2 differs on exactly one line."""
    block = extract_tagged_block(_page(), "f-small-hole-pec-tune", "yaml")
    tracked = TUNE_FILE.read_bytes()
    sha_b = hashlib.sha256(block.encode()).hexdigest()
    sha_t = hashlib.sha256(tracked).hexdigest()
    print(f"\n[WF-17 s6] (Y=) block sha256 = {sha_b}\n[WF-17 s6] (Y=) file  sha256 = {sha_t}")
    assert block.encode() == tracked
    mutated = block.replace("  k: 1\n", "  k: 2\n")
    diff = _differing_lines(mutated, tracked.decode())
    print(f"[WF-17 s6] (c1) k: 2 copy differs on lines {diff}")
    assert mutated.encode() != tracked
    assert len(diff) == 1


def test_wf17_s6_s_keys_equal_tune_summary():
    """(S) the page's key list == set(_tune_summary(stub, {"k": 1}, comm))."""
    stub = SimpleNamespace(
        tune=lambda **kw: TuneResult(
            k=1,
            frequency_hz=FREQUENCY_64_HZ,
            z_k=np.zeros(4, complex),
            c_mode_farad=1.0,
            off_diagonal_residual=0.0,
        )
    )
    keys = set(_tune_summary(stub, {"k": 1}, MPI.COMM_WORLD))
    block = extract_tagged_block(_page(), "tune-summary-keys", "text")
    listed = [line.split()[0] for line in block.splitlines() if line.strip()]
    print(f"\n[WF-17 s6] (S) page keys={listed}\n[WF-17 s6] (S) _tune_summary keys={sorted(keys)}")
    assert len(listed) == len(set(listed))
    assert set(listed) == keys


def test_wf17_s6_e_output_from_imported_constants_and_c4():
    """(E) the output line == f"{TH17 record * (1 + item 429 a):.4e}"; (c4) the k = 2 value differs."""
    block = extract_tagged_block(_page(), "tune-yaml-output", "text")
    lines = [line for line in block.splitlines() if line.strip()]
    assert len(lines) == 1
    m = OUTPUT_RE.match(lines[0])
    assert m is not None
    expected = _expected_output()
    k2 = f"{TH17_STEP3_C_MODE_RECORD_F * ITEM429_K2_FACTOR:.4e}"
    print(
        f"\n[WF-17 s6] (E) page={m.group(1)} expected={expected} "
        f"(U49 YAML reading {U49_YAML_C_MODE_PRINT:.9e} -> {U49_YAML_C_MODE_PRINT:.4e}, printed only)"
        f"\n[WF-17 s6] (c4) k = 2 value {k2} != page: {k2 != m.group(1)}"
    )
    assert m.group(1) == expected
    assert k2 != m.group(1)


def test_wf17_s6_l_nodes_registered_links_resolve_and_c2():
    """(L) cited nodes registered; validation/api links resolve (anchors on headings); (c2) a planted node is named."""
    text = _page()
    nodes = NODE_RE.findall(text)
    registered = _registered_nodes()
    missing = _unregistered(nodes, registered)
    print(
        f"\n[WF-17 s6] (L) nodes cited={len(nodes)} registered under {REGISTRY_NAMES}; unregistered={missing}"
    )
    assert len(nodes) == 3
    assert missing == []

    planted = "tests/package/test_wf1_yaml_tune.py::test_wf1_s11_no_such_node"
    caught = _unregistered(nodes + [planted], registered)
    print(f"[WF-17 s6] (c2) planted node -> unregistered={caught}")
    assert caught == [planted]

    links = LINK_RE.findall(text)
    val_heads = {
        line.lstrip("#").strip().lower()
        for line in (GUIDE / "validation.md").read_text().splitlines()
        if line.startswith("### ")
    }
    required = {"validation.md", "api/Study.md", "api/run_study_file.md"}
    targets = {lk.split("#")[0] for lk in links}
    print(f"[WF-17 s6] (L) links={len(links)}; required present={required <= targets}")
    assert required <= targets
    for lk in links:
        path, _, anchor = lk.partition("#")
        assert (PAGE.parent / path).resolve().exists(), lk
        if path == "validation.md" and anchor:
            assert anchor in val_heads, lk
    for ex in ("19_study_tune_mode_capacitance.md", "21_femem_run_tune_section.md"):
        assert f"../../examples/mri/{ex}" in links


def test_wf17_s6_d_guide_zero_dead(capsys):
    """(D) run_checker --docs-root docs/guide --check-links reads dead == 0."""
    rc, counts = run_checker("--docs-root", "docs/guide", "--check-links")
    out = capsys.readouterr().out
    m = LINK_PASS_RE.search(out)
    with capsys.disabled():
        print(
            f"\n[WF-17 s6] (D) rc={rc} {m.group(0) if m else 'no link-pass line'} counts={counts}"
        )
    assert m is not None
    assert int(m.group(3)) == 0
    assert counts["dead"] == 0


def test_wf17_s6_c3_planted_dead_link(tmp_path, capsys):
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
            f"\n[WF-17 s6] (c3) rc {rc0}/{rc1}; copy dead={counts0['dead']}; planted copy dead={counts1['dead']}"
        )
    assert counts0["dead"] == 0
    assert counts1["dead"] == 1
