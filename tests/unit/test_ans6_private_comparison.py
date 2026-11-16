"""`ANS-6` step 3 (§9 item 60, chain F6b): the private-mode comparison writer.

Everything here runs on the tracked ``metrics.json`` and a **synthetic** AED
JSON written into ``tmp_path`` — never the real ``aed_results/`` (the module's
``CASE_DIR`` / ``AED_RESULTS_PATH`` are repointed before any call, and
``ans:6`` itself is never run).

Anchors (asserted):

1. One sentinel ``1.<iii>7351928`` (11 significant digits, unique; the
   imaginary part of a complex value is ``2.<iii>7351928``, so re != im) per
   AED key. Every key the private table renders appears exactly once, in its
   (section, frequency, column, entry, order) cell; every key it does not
   render appears nowhere. Counts printed.
2. The tracked render equals the committed ``COMPARISON.md`` byte for byte and
   holds no sentinel.
3. A scratch git repo holding the case, ``aed_results/`` and the private file
   ignored as in ``.gitignore``: ``check_private_leak.py --root <scratch>
   --audit`` exits 0 with ``clean:`` and > 0 reference figures (``SKIPPED``
   fails).
4. ``FEM_EM_ANS_SKIP_PRIVATE=1``: no private file, and the note names the knob.

Negative control (asserted): the private render committed as the tracked
``COMPARISON.md`` in the scratch repo makes the checker exit 1 with the
location ``COMPARISON.md``; the hit count is printed beside the rendered
sentinel-token count it is predicted to approximate.
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

REPO = Path(__file__).resolve().parents[2]
CASE_REL = Path("examples/ansys_benchmarks/ans6_copper_birdcage_four_port_10_64_128MHz")
CASE = REPO / CASE_REL
SCRIPT = CASE / "06_copper_birdcage_four_port_10_64_128MHz.py"
CHECKER = REPO / "scripts" / "testing" / "check_private_leak.py"

LABELS = ("10 MHz", "64 MHz", "128 MHz")
COLUMNS = ("cu", "pec")
ORDERS = ("zero_order", "first_order")
CLASS_NAMES = ("self", "adjacent", "opposite")
S_KEYS = tuple(f"S{i}{j}" for i in range(1, 5) for j in range(1, 5))
META_KEYS = ("elements", "basis_order", "passes", "delta_s", "solve_time",
             "port_model", "coil_treatment")
#: `.gitignore:134–135`, verbatim.
GITIGNORE = ("examples/ansys_benchmarks/*/aed_results/\n"
             "examples/ansys_benchmarks/*/COMPARISON_private.md\n")
#: `scripts/testing/test_private_leak.py:63–67`.
GIT_ENV = {
    "GIT_CONFIG_GLOBAL": os.devnull, "GIT_CONFIG_NOSYSTEM": "1",
    "GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@example.invalid",
    "GIT_COMMITTER_NAME": "t", "GIT_COMMITTER_EMAIL": "t@example.invalid",
}


def _env() -> dict:
    e = dict(os.environ)
    for k in ("PRIVATE_LEAK_ROOT", "PRIVATE_LEAK_SECRET_GLOBS", "GIT_DIR",
              "GIT_INDEX_FILE", "GIT_WORK_TREE"):
        e.pop(k, None)
    e.update(GIT_ENV)
    return e


def _git(repo: Path, *args: str) -> str:
    return subprocess.run(["git", "-C", str(repo), *args], env=_env(), check=True,
                          capture_output=True, text=True).stdout


def _checker(repo: Path) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(CHECKER), "--root", str(repo), "--audit"],
                          cwd=repo, env=_env(), capture_output=True, text=True, timeout=60)


class Sentinels:
    """``1.<iii>7351928`` per key, in insertion order; complex im = ``2.<iii>7351928``."""

    def __init__(self):
        self.by_key: dict[tuple, str] = {}
        self.complex_keys: set[tuple] = set()

    def _next(self, where: tuple) -> int:
        i = len(self.by_key)
        assert i < 1000
        self.by_key[where] = f"1.{i:03d}7351928"
        return i

    def cplx(self, where):
        i = self._next(where)
        self.complex_keys.add(where)
        return {"re": float(f"1.{i:03d}7351928"), "im": float(f"2.{i:03d}7351928")}

    def real(self, where):
        return float(f"1.{self._next(where):03d}7351928")

    def text(self, where):
        return f"synthetic {where[-1].replace('|', '/')} 1.{self._next(where):03d}7351928 (not AED)"


def _synthetic_aed() -> tuple[dict, Sentinels]:
    """The PyAEDT script's key layout (`aed/ans6_hfss_pyaedt.py:36–52`), synthetic values."""
    sen = Sentinels()
    out: dict = {"_provenance": "synthetic, ANS-6 step 3 test", "_meta": {}, "_reading": ""}
    for column in COLUMNS:
        for order in ORDERS:
            blk = f"{column}_{order}"
            block = {}
            for label in LABELS:
                for s in S_KEYS:
                    block[f"{label}|{s}"] = sen.cplx((blk, f"{label}|{s}"))
                for c in CLASS_NAMES:
                    block[f"{label}|class_{c}"] = sen.cplx((blk, f"{label}|class_{c}"))
                block[f"{label}|Z11"] = sen.cplx((blk, f"{label}|Z11"))
                block[f"{label}|identities"] = sen.text((blk, f"{label}|identities"))
                block[f"{label}|power"] = sen.text((blk, f"{label}|power"))
                for k in ("accepted_w", "p_phantom_w", "p_coil_w", "loss_partition"):
                    block[f"{label}|{k}"] = sen.real((blk, f"{label}|{k}"))
            for k in META_KEYS:
                block[k] = sen.text((blk, k))
            out[blk] = block
    diff = {}
    for order in ORDERS:
        for label in LABELS:
            for c in CLASS_NAMES:
                diff[f"{order}|{label}|{c}"] = sen.cplx(("_column_difference", f"{order}|{label}|{c}"))
            diff[f"{order}|{label}|_rel_norm"] = sen.real(("_column_difference", f"{order}|{label}|_rel_norm"))
    out["_column_difference"] = diff
    return out, sen


def _expected_cell(where: tuple):
    """(section substring, row prefix cells, column index) the key must land in, or None."""
    blk, key = where
    if blk == "_column_difference":
        order, label, cls = key.split("|")
        if cls == "_rel_norm":
            return None
        return ("## Column difference", (label, cls), 3 if order == "zero_order" else 4)
    column, order = blk.split("_", 1)
    zo = order == "zero_order"
    if key in ("elements", "basis_order", "port_model"):
        row = {"elements": "Elements", "basis_order": "Basis order", "port_model": "Port model"}[key]
        return ("## Solve metadata", (row,), 2 if zo else 3)
    if "|" not in key:
        return None
    label, entry = key.split("|")
    if entry.startswith("class_"):
        section = f"## S-matrix, C4 classes, column {'Cu' if column == 'cu' else 'PEC'}"
        return (section, (label, column, entry[len("class_"):]), 4 if zo else 5)
    if entry in S_KEYS and column == "cu":
        return (f"### {label}", (entry,), 2 if zo else 3)
    if entry == "identities":
        section = "## Identities (computable" if column == "cu" else "## Identities — column PEC"
        return (section, (label,), 4 if zo else 5)
    if entry == "loss_partition" and column == "cu":
        return ("## Loss partition", (label,), 3 if zo else 4)
    return None


def _cells(line: str) -> list[str]:
    return [c.strip() for c in line.strip().strip("|").split("|")]


def _locate(text: str, sentinel: str) -> list[tuple[str, list[str], int]]:
    """Every (section heading, row cells, column index) holding ``sentinel``."""
    hits, section = [], ""
    for line in text.split("\n"):
        if line.startswith("#"):
            section = line
        if sentinel in line:
            cells = _cells(line)
            for idx, cell in enumerate(cells):
                if sentinel in cell:
                    hits.append((section, cells, idx))
    return hits


@pytest.fixture(scope="module")
def ans6():
    spec = importlib.util.spec_from_file_location("ans6_private_writer_under_test", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture
def scratch(tmp_path, ans6, monkeypatch):
    """A scratch repo holding the case; the module writes into it and nowhere else."""
    repo = tmp_path / "scratch"
    case = repo / CASE_REL
    (case / "aed_results").mkdir(parents=True)
    _git(repo, "init", "-q", "-b", "main")
    (repo / ".gitignore").write_text(GITIGNORE)
    shutil.copyfile(CASE / "metrics.json", case / "metrics.json")
    aed, sen = _synthetic_aed()
    aed_path = case / "aed_results" / "ans6_aed_results.json"
    aed_path.write_text(json.dumps(aed, indent=2))
    monkeypatch.setattr(ans6, "CASE_DIR", case)
    monkeypatch.setattr(ans6, "AED_RESULTS_PATH", aed_path)
    monkeypatch.delenv(ans6.SKIP_PRIVATE_ENV, raising=False)
    m = json.loads((case / "metrics.json").read_text())
    return repo, case, m, sen


def test_private_writer_sentinels_tracked_bytes_and_leak_checker(ans6, scratch):
    repo, case, m, sen = scratch
    tracked = ans6._write_comparison(m)
    note = ans6._write_private_if_available(m)
    private = case / "COMPARISON_private.md"
    assert tracked == case / "COMPARISON.md"
    assert private.exists(), note
    print(f"[ANS-6 step 3] private note: {note}")

    # -- anchor 2: tracked render byte-identical, no sentinel ------------------
    committed = (CASE / "COMPARISON.md").read_bytes()
    rendered = tracked.read_bytes()
    assert rendered == committed, "tracked render differs from the committed COMPARISON.md"
    tracked_text = rendered.decode()
    assert "7351928" not in tracked_text
    print(f"[ANS-6 step 3] tracked render == committed COMPARISON.md: "
          f"{len(rendered)} bytes, byte-identical, 0 sentinels")

    # -- anchor 1: every sentinel once, in its cell; unrendered keys absent ---
    text = private.read_text()
    rendered_keys, absent_keys, tokens = 0, 0, 0
    for where, s in sen.by_key.items():
        expect = _expected_cell(where)
        count = text.count(s)
        if expect is None:
            assert count == 0, f"{where} not rendered but found {count}x"
            absent_keys += 1
            continue
        assert count == 1, f"{where}: sentinel appears {count}x"
        (section, cells, idx), = _locate(text, s)
        want_section, want_prefix, want_idx = expect
        assert want_section in section, (where, section)
        assert idx == want_idx, (where, idx, want_idx)
        assert all(cells[k] == p or cells[k].startswith(p + " ")
                   for k, p in enumerate(want_prefix)), (where, cells)
        rendered_keys += 1
        if where in sen.complex_keys:
            # the imaginary part sits right beside it, in the same cell
            assert text.count(f"2.{s[2:]}") == 1, where
            tokens += 2
        else:
            tokens += 1
    assert rendered_keys > 0
    print(f"[ANS-6 step 3] sentinels: {len(sen.by_key)} keys; {rendered_keys} rendered "
          f"exactly once in their cells; {absent_keys} unrendered keys absent; "
          f"{tokens} sentinel digit tokens in the private render")

    # -- anchor 3: scratch repo, ignore rules, --audit clean -------------------
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "scratch: ANS-6 case, tracked render")
    files = _git(repo, "ls-files").split()
    assert str(CASE_REL / "COMPARISON.md") in files
    assert not any("aed_results" in f or "COMPARISON_private" in f for f in files), files
    clean = _checker(repo)
    out = clean.stdout + clean.stderr
    print(f"[ANS-6 step 3] checker on the tracked render: exit {clean.returncode}; "
          f"{out.strip()}")
    assert "SKIPPED" not in out
    assert clean.returncode == 0 and out.startswith("clean:"), out
    n_ref = int(re.search(r"carry none of the (\d+) reference figures", out).group(1))
    assert n_ref > 0

    # -- negative control: private render committed as the tracked file -------
    (case / "COMPARISON.md").write_text(text)
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "scratch: leaking writer variant")
    leak = _checker(repo)
    out = leak.stdout + leak.stderr
    hits = re.search(r"(\d+) match\(es\)", out)
    print(f"[ANS-6 step 3] NEGATIVE CONTROL checker on the leaking variant: exit "
          f"{leak.returncode}; hits {hits.group(1) if hits else None} "
          f"(predicted ~ {tokens} sentinel tokens); reference figures {n_ref}")
    assert leak.returncode == 1, out
    assert f"{CASE_REL / 'COMPARISON.md'}:" in out
    assert "7351928" not in out  # locations only, never the figure


def test_skip_private_knob(ans6, scratch, monkeypatch, capsys):
    repo, case, m, _ = scratch
    monkeypatch.setenv(ans6.SKIP_PRIVATE_ENV, "1")
    assert ans6.AED_RESULTS_PATH.exists()
    note = ans6._write_private_if_available(m)
    print(f"[ANS-6 step 3] with {ans6.SKIP_PRIVATE_ENV}=1: {note}")
    assert not (case / "COMPARISON_private.md").exists()
    assert "FEM_EM_ANS_SKIP_PRIVATE" in note
