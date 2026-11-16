"""`ANS-7` lane 1 (§9 item 460, §10 T45): the public "Palace at matched order" section.

Render-only, on scratch copies of the committed ``metrics.json`` / ``COMPARISON.md``:

(X)  every figure of the section traces to a tracked read log line (verbatim, rounded, or the
     max of three verbatim per-point values), and occurs in the committed section;
(R)  ``_write_comparison`` on the committed metrics reproduces the committed ``COMPARISON.md``
     byte for byte, the Palace heading once, after the 127.74 MHz ESR section and before the
     closing "Not claimed" line;
(I)  a stand-alone ``_write_comparison_esr128`` (the XT128R path) leaves the committed file
     byte-identical — the Palace section survives its cut-to-tail;
(B)  no bridge to AED (§5.4 (iv)): "AED" only in the two licensed phrases, no private path;
(R0) the repo's ``COMPARISON.md`` / ``metrics.json`` are untouched.
Negative controls (asserted): (N-cite) a citation moved by one line fails (X); (N-digit) a
one-digit edit of the committed file differs from the render; (N-old) the pinned pre-change
writer ``3c359dd0`` (copied under ``logs/t45_pin/`` at the script's relative path by a host
``git show`` before the window) renders no Palace section and differs.
"""

from __future__ import annotations

import gzip
import importlib.util
import json
import re
import shutil
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
CASE = REPO / "examples" / "ansys_benchmarks" / "ans7_tuned_copper_birdcage_driven_sweep"
SCRIPT = CASE / "07_tuned_copper_birdcage_driven_sweep.py"
PIN_SHA = "3c359dd0"
# Mirrored depth: the module body resolves ``parents[3]``.
OLD_SCRIPT = REPO / "logs" / "t45_pin" / SCRIPT.relative_to(REPO)
LOGS = REPO / "docs" / "testing" / "logs"
HEAD128 = "## Tuned with capacitor ESR at 127.74 MHz (§10 T35 / T36)"
# §9 item 514 (ruling R1 (A), R2 licence, 2026-10-18 03:00 review): the re-pointed closing line.
TAIL = ("\nNot claimed: any absolute `Z_in`, a match; an absolute Q or R_in (INCONCLUSIVE, T48); "
        "the tables above are one mesh each — the h-statements are the verdict section's.\n")

L421 = "20261013T141248Z_OPS-61-item421-W4-read.log"
L420 = "20261013T140547Z_OPS-61-item420-W6-read.log"
L436 = "20261013T215115Z_OPS-61-item436-W6-read.log"
L444 = "20261014T094041Z_OPS-61-item444-W7-read.log"
L448 = "20261014T171034Z_OPS-61-item448-W6-read.log"
LP64 = "20261013T140240Z_OPS-61-item420-W3-palace-t39-esr.log"
LPH = "20261014T093349Z_OPS-61-item444-W3-palace-t43-f64.0.log"
CMP = "COMPARISON.md"  # the committed file itself (ours, metrics-rendered)
# §9 item 514 (R2 licence): CMP line 54 -> 63, the verdict section adds 9 lines above it.

# (cell, log, line, mode). Modes: "v" verbatim substring; "q" / "r" rounded .4f / .3f of the
# "<X> vs <Y>" pair on the Q / R_in line, cell "P / O"; "max" cell is max of the three per-point
# ties on lines line..line+2, each verbatim.
TRACE = (
    ("+6.530e-11", L421, 48, "v"), ("-4.871e-09", L421, 49, "v"), ("-7.456e-10", L421, 50, "v"),
    ("2.580e-10", L421, 38, "max"),
    ("2.565e-10", L420, 46, "v"), ("1.994e-09", L420, 48, "v"), ("parallel", L420, 49, "v"),
    ("164.5434 / 164.5434", L421, 49, "q"), ("-2.583046e-06", L421, 48, "v"),
    ("-2.583111e-06", CMP, 63, "v"), ("3175.979 / 3175.979", L421, 50, "r"),
    ("+6.073e-11", L436, 51, "v"), ("-2.301744e-06", L436, 51, "v"), ("-2.301805e-06", L436, 51, "v"),
    ("+5.444e-09", L436, 52, "v"), ("254.4693 / 254.4693", L436, 52, "q"),
    ("+2.906e-09", L436, 53, "v"), ("9790.704 / 9790.704", L436, 53, "r"),
    ("2.373e-10", L436, 39, "max"), ("2.308e-10", L436, 42, "v"), ("5.648e-09", L436, 43, "v"),
    ("+6.437e-11", L444, 51, "v"), ("-2.588508e-06", L444, 51, "v"), ("-2.588573e-06", L444, 51, "v"),
    ("+1.829e-09", L444, 52, "v"), ("164.2158 / 164.2158", L444, 52, "q"),
    ("+6.507e-10", L444, 53, "v"), ("3166.748 / 3166.748", L444, 53, "r"),
    ("2.591e-10", L444, 39, "max"), ("2.499e-10", L444, 42, "v"), ("1.717e-09", L444, 43, "v"),
    ("-5.462608e-09 / -5.461672e-09", L444, 54, "v"), ("-1.991402e-03 / -1.991409e-03", L444, 54, "v"),
    ("-2.906654e-03 / -2.906656e-03", L444, 54, "v"),
    ("+5.912e-11", L448, 51, "v"), ("-2.307258e-06", L448, 51, "v"), ("-2.307317e-06", L448, 51, "v"),
    ("+1.027e-08", L448, 52, "v"), ("252.3012 / 252.3012", L448, 52, "q"),
    ("+4.872e-09", L448, 53, "v"), ("9698.113 / 9698.113", L448, 53, "r"),
    ("2.408e-10", L448, 39, "max"), ("2.093e-10", L448, 42, "v"), ("4.356e-09", L448, 43, "v"),
    ("-5.514397e-09 / -5.512787e-09", L448, 54, "v"), ("-8.519899e-03 / -8.519904e-03", L448, 54, "v"),
    ("-9.456957e-03 / -9.456959e-03", L448, 54, "v"),
    ("v0.18.1", LP64, 45, "v"), ("v0.18.1", LPH, 45, "v"),
    ("80798", LP64, 68, "v"), ("140270", LPH, 68, "v"),
    ("ND (p = 2): 570332", LP64, 104, "v"), ("ND (p = 2): 968016", LPH, 104, "v"),
)
# Cells whose section spelling differs from the log's (thousands spaces, "v" prefix kept).
SECTION_SPELLING = {"80798": "80 798", "140270": "140 270", "ND (p = 2): 570332": "570 332",
                    "ND (p = 2): 968016": "968 016", "parallel": "read as parallel"}
PAIR = re.compile(r"(?:Q_Pal|R_in,Pal) ([-+0-9.e]+)(?: ohm)? vs (?:Q_A|R_in,A) ([-+0-9.e]+)")
TIE = re.compile(r"\|S11_Pal - S11_A\|/\|S11_A\| ([0-9.e+-]+)")


def _lines(name: str) -> list[str]:
    if name == CMP:
        return (CASE / CMP).read_text().split("\n")
    p = LOGS / name
    if p.exists():
        return p.read_text().split("\n")
    return gzip.open(str(p) + ".gz", "rt").read().split("\n")


def _holds(cell: str, log: str, line: int, mode: str) -> bool:
    ls = _lines(log)
    s = ls[line - 1]
    if mode == "v":
        return cell in s
    if mode in ("q", "r"):
        mt = PAIR.search(s)
        if not mt:
            return False
        fmt = ".4f" if mode == "q" else ".3f"
        return f"{float(mt.group(1)):{fmt}} / {float(mt.group(2)):{fmt}}" == cell
    if mode == "max":
        ties = []
        for k in range(3):
            mt = TIE.search(ls[line - 1 + k])
            if not mt:
                return False
            ties.append(mt.group(1))
        return max(ties, key=float) == cell
    raise ValueError(mode)


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def ans7():
    return _load("ans7_writer_palace_ut", SCRIPT)


@pytest.fixture(scope="module")
def ans7_old():
    if not OLD_SCRIPT.exists():
        pytest.fail(f"{OLD_SCRIPT} missing: run `git show {PIN_SHA}:<script> > {OLD_SCRIPT}` first")
    return _load(f"ans7_writer_{PIN_SHA}_ut", OLD_SCRIPT)


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


def _section(text: str, head: str) -> str:
    i = text.index("\n" + head)
    j = text.index(TAIL, i)
    return text[i:j]


def test_palace_section_traced(ans7):
    committed = (CASE / "COMPARISON.md").read_text()
    section = _section(committed, ans7.PALACE_HEAD)
    # The inserter writes the section plus one newline before the tail's own leading newline.
    assert section == "\n".join(ans7.PALACE_SECTION) + "\n"
    bad = [t for t in TRACE if not _holds(*t)]
    for t in TRACE:
        print(f"(X) {t[0]!r:34} {t[1]}:{t[2]} [{t[3]}] {'ok' if t not in bad else 'FAIL'}")
    assert not bad, bad
    for cell, *_ in TRACE:
        assert SECTION_SPELLING.get(cell, cell) in section, cell
    print(f"(X) {len(TRACE)} traced cells hold and occur in the committed section")
    # (N-cite): Q +1.027e-08 cited one line late (448 W6:53) fails.
    assert not _holds("+1.027e-08", L448, 53, "v")
    print("(N-cite) +1.027e-08 at 448 W6:53 fails (X): as asserted")


def test_palace_section_rendered(ans7, ans7_old, tmp_path):
    committed = (CASE / "COMPARISON.md").read_bytes()
    metrics_text = (CASE / "metrics.json").read_text()
    text = committed.decode()

    # (R)
    new = _render(ans7, tmp_path, json.loads(metrics_text))
    n_pal = new.decode().count(ans7.PALACE_HEAD)
    print(f"\n(R) rendered {len(new)} B, committed {len(committed)} B, Palace heading count {n_pal}")
    assert n_pal == 1
    assert new == committed
    assert text.index(HEAD128) < text.index(ans7.PALACE_HEAD) < text.rindex(TAIL)
    assert text.endswith(TAIL)

    # (I)
    idir = tmp_path / "esr128"
    idir.mkdir()
    shutil.copy(CASE / "COMPARISON.md", idir / "COMPARISON.md")
    old_dir = ans7.CASE_DIR
    ans7.CASE_DIR = idir
    try:
        ans7._write_comparison_esr128(json.loads(metrics_text))
    finally:
        ans7.CASE_DIR = old_dir
    after = (idir / "COMPARISON.md").read_bytes()
    print(f"(I) stand-alone _write_comparison_esr128: byte-identical {after == committed}")
    assert after == committed

    # (B)
    sec = _section(text, ans7.PALACE_HEAD)
    stripped = (sec.replace("No AED column: no AED ESR case exists", "")
                .replace("any AED value", ""))
    print(f"(B) 'AED' outside the licensed phrases: {stripped.count('AED')}")
    assert "AED" not in stripped
    for bad in ("docs/private", "aed_results", "COMPARISON_private"):
        assert bad not in sec

    # (N-digit)
    edited = text.replace("+6.530e-11", "+6.540e-11", 1).encode()
    assert edited != committed and edited != new
    print("(N-digit) +6.530e-11 -> +6.540e-11 differs from the render: as asserted")

    # (N-old)
    odir = tmp_path / "old"
    odir.mkdir()
    old = _render(ans7_old, odir, json.loads(metrics_text))
    n_old = old.decode().count(ans7.PALACE_HEAD)
    print(f"(N-old) {PIN_SHA} writer: {len(old)} B, Palace heading count {n_old}, "
          f"differs {old != committed}")
    assert n_old == 0
    assert old != committed

    # (R0)
    assert (CASE / "COMPARISON.md").read_bytes() == committed
    assert (CASE / "metrics.json").read_text() == metrics_text
    print("(R0) committed COMPARISON.md / metrics.json unchanged")
