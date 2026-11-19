"""`ANS-7` lane 1 (§9 item 400): ``_write_comparison`` keeps the 127.74 MHz ESR section.

Render-only, on scratch copies of the committed ``metrics.json`` / ``COMPARISON.md``:

(R)     the writer on the committed metrics reproduces the committed ``COMPARISON.md``
        byte for byte, the 127.74 MHz ESR section included (count 1);
(R0)    the repo's ``COMPARISON.md`` is untouched by the test;
(N-old) the pinned pre-change writer (``b01dc790``, copied under ``logs/ans7_writer_b01dc790/``
        at the script's own relative path by a host ``git show`` before the window) drops the
        section and so differs;
(N-key) with ``tuned_larmor_esr["127.74"]`` deleted the new writer adds no section and
        equals (N-old)'s output byte for byte.
"""

from __future__ import annotations

import importlib.util
import json
import shutil
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
CASE = REPO / "examples" / "ansys_benchmarks" / "ans7_tuned_copper_birdcage_driven_sweep"
SCRIPT = CASE / "07_tuned_copper_birdcage_driven_sweep.py"
# The pinned copy mirrors the script's depth: its module body resolves ``parents[3]``.
OLD_SCRIPT = REPO / "logs" / "ans7_writer_b01dc790" / SCRIPT.relative_to(REPO)
HEAD128 = "## Tuned with capacitor ESR at 127.74 MHz (§10 T35 / T36)"


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def ans7():
    return _load("ans7_writer_ut", SCRIPT)


@pytest.fixture(scope="module")
def ans7_old():
    if not OLD_SCRIPT.exists():
        pytest.fail(f"{OLD_SCRIPT} missing: run `git show b01dc790:<script> > {OLD_SCRIPT}` first")
    return _load("ans7_writer_b01dc790_ut", OLD_SCRIPT)


def _render(mod, tmp: Path, metrics: dict) -> bytes:
    shutil.copy(CASE / "COMPARISON.md", tmp / "COMPARISON.md")
    shutil.copy(CASE / "metrics.json", tmp / "metrics.json")
    old = mod.CASE_DIR
    mod.CASE_DIR = tmp
    try:
        mod._write_comparison(metrics)
    finally:
        mod.CASE_DIR = old
    return (tmp / "COMPARISON.md").read_bytes()


def test_writer_reproduces_committed_comparison(ans7, ans7_old, tmp_path):
    committed = (CASE / "COMPARISON.md").read_bytes()
    metrics_text = (CASE / "metrics.json").read_text()

    # (R)
    new = _render(ans7, tmp_path, json.loads(metrics_text))
    n_new = new.decode().count(HEAD128)
    print(
        f"\n(R) rendered {len(new)} B, committed {len(committed)} B, "
        f"127.74 MHz heading count {n_new}"
    )
    assert n_new == 1
    assert new == committed

    # (N-old)
    old_dir = tmp_path / "old"
    old_dir.mkdir()
    old = _render(ans7_old, old_dir, json.loads(metrics_text))
    n_old = old.decode().count(HEAD128)
    print(f"(N-old) b01dc790 writer: {len(old)} B, 127.74 MHz heading count {n_old}")
    assert n_old == 0
    assert old != committed

    # (N-key)
    key_dir = tmp_path / "nokey"
    key_dir.mkdir()
    m = json.loads(metrics_text)
    del m["degree2_exact"]["tuned_larmor_esr"]["127.74"]
    nokey = _render(ans7, key_dir, m)
    n_key = nokey.decode().count(HEAD128)
    print(
        f"(N-key) new writer without the 127.74 key: {len(nokey)} B, heading count {n_key}, "
        f"equals (N-old): {nokey == old}"
    )
    assert n_key == 0
    # §9 item 460 (licensed by the 2026-10-15 03:00 review ruling, item-417 ruling-(2) pattern):
    # the new writer adds the unconditional Palace section, which b01dc790 predates; cut it
    # (present exactly once) before the comparison, every other byte still compared.
    nokey_text = nokey.decode()
    assert nokey_text.count(ans7.PALACE_HEAD) == 1
    i = nokey_text.index("\n" + ans7.PALACE_HEAD)
    j = i + len("\n".join(ans7.PALACE_SECTION)) + 1
    nokey_text = nokey_text[:i] + nokey_text[j:]
    # §9 item 514 (licensed by the 2026-10-18 03:00 review ruling R2): the new writer adds the
    # verdict section (present exactly once) and re-points the closing line (ruling R1 (A)),
    # both of which b01dc790 predates; cut the section and restore the old closing line (present
    # exactly once) before the comparison, every other byte still compared.
    assert nokey_text.count(ans7.VERDICT_HEAD) == 1
    vs = "\n".join(ans7.VERDICT_SECTION) + "\n"
    assert nokey_text.count(vs) == 1
    nokey_text = nokey_text.replace(vs, "")
    old_close = (
        "Not claimed: any absolute `Z_in`, a match; no h-statement (T12′); no absolute "
        "f₀ / Q / R_in."
    )
    assert nokey_text.count("\n" + ans7.CLOSING_LINE + "\n") == 1
    nokey_text = nokey_text.replace("\n" + ans7.CLOSING_LINE + "\n", "\n" + old_close + "\n")
    nokey = nokey_text.encode()
    assert nokey == old

    # (R0)
    assert (CASE / "COMPARISON.md").read_bytes() == committed
    assert (CASE / "metrics.json").read_text() == metrics_text
    print("(R0) committed COMPARISON.md / metrics.json unchanged")
