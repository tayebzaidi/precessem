"""§9 item 997 — the retention-safe log readers read the same text in either form.

`scripts/maintenance/housekeeping.py` gzips every tracked log at 7 days, cited
or not, so a reader that opened ``docs/testing/logs/<name>.log`` as plain text
went red on the 11-15 sweep (known-issues 2026-11-15). Item 997 routes
``_custody`` (3 T module, imported by the 1.5 T one), ``_parse_c52`` (P80,
imported by P81), the closed-form records page's ``read_one`` and the lattice
module's ``_c39`` through ``tests/log_records.read_log_text(errors="strict")``.

Gates (all asserted):

  (L)       ``locale.getpreferredencoding(False)`` is UTF-8 (the plain branch of
            ``read_log_text`` passes no ``encoding=``; the old readers passed
            ``encoding="utf-8"``).
  (Id-b)    for every log in the item's table, a copy in the *opposite* form in
            ``tmp_path`` reads byte-identical text to the tracked one, and both
            equal ``TextIOWrapper(BytesIO(raw), "utf-8", newline=None).read()``
            — what ``Path.read_text`` does.
  (Id-c)    ``_custody(prices, logs=tmp) == _custody(prices, logs=LOGS)`` and
            both are all-empty, for the 1.5 T (48 / 52 / 59) and 3 T
            (58 / 60 / 61) ledger entries.
  (Id-p)    each other edited reader, with its path constant monkeypatched to
            the opposite-form copy, parses what it parses on the real path.
  (N-absent) a log in neither form raises ``FileNotFoundError`` naming both.
  (N-trunc) entry 48's ``.gz`` truncated to half raises in ``_custody``.
  (N-crc)   the same with one CRC32-trailer byte flipped raises
            ``gzip.BadGzipFile``.

No copy is ever written under ``docs/testing/logs/``; every opposite-form copy
lives in the rank's own ``tmp_path``.
"""

from __future__ import annotations

import gzip
import io
import locale
import zlib
from pathlib import Path

import pytest
from mpi4py import MPI

from tests.log_records import read_log_text, resolve_log_path
from tests.package import test_wf17_clinical_closed_form_records_page as cf
from tests.package.test_wf17_estimate_measured_3t import LOGS, _custody
from tests.package.test_wf17_estimate_table_page import gen
from tests.validation import test_mat21_eddy_source_probe_points_closed_form as p80
from tests.validation import test_mat21_lattice_ball_superposition_closed_form as lat

TAG = "[WF-17 i997]"

# The item's table (WF-17.md § "§9 item 997"), every basename it names.
TABLE = (
    "20261102T080008Z_GEO-36-step9-I36.log",  # entry 48
    "20261103T080009Z_GEO-38-C6.log",  # entry 52 / P80's C52_LOG
    "20261110T080008Z_GEO-38-P33.log",  # entry 59
    "20261109T080008Z_GEO-38-C6-3T.log",  # entry 58
    "20261111T080008Z_GEO-38-P33-3T.log",  # entry 60
    "20261112T080008Z_GEO-36-step11-I36-3T.log",  # entry 61
    "20261112T120219Z_MAT-21-P88-W1.log",
    "20261112T182020Z_MAT-21-P89-W1b.log",
    "20261110T103343Z_MAT-21-P83-W1.log",
    "20261111T103024Z_GEO-38-C39-W1.log",
)


def _p0(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _raw_plain_bytes(name: str) -> bytes:
    plain, gz = LOGS / name, LOGS / f"{name}.gz"
    if plain.exists():
        return plain.read_bytes()
    if gz.exists():
        return gzip.decompress(gz.read_bytes())
    pytest.fail(f"{name}: neither form on disk")


def _opposite(name: str, tmp: Path) -> Path:
    """Write ``name`` into ``tmp`` in the form opposite to the tracked one; return
    the plain-form path the readers are given."""
    plain, gz = LOGS / name, LOGS / f"{name}.gz"
    tmp.mkdir(parents=True, exist_ok=True)
    if plain.exists():
        (tmp / f"{name}.gz").write_bytes(gzip.compress(plain.read_bytes()))
        assert not (tmp / name).exists()
    elif gz.exists():
        (tmp / name).write_bytes(gzip.decompress(gz.read_bytes()))
        assert not (tmp / f"{name}.gz").exists()
    else:
        pytest.fail(f"{name}: neither form on disk")
    return tmp / name


def _form(path: Path) -> str:
    return resolve_log_path(path).name.removeprefix(path.name) or ".log"


def test_i997_l_locale_is_utf8():
    enc = locale.getpreferredencoding(False)
    _p0(f"\n{TAG} (L) locale.getpreferredencoding(False) = {enc!r}")
    assert enc.replace("-", "").upper() == "UTF8", enc


def test_i997_id_b_byte_identity(tmp_path):
    tmp = tmp_path / "opp"
    for name in TABLE:
        real = LOGS / name
        opp = _opposite(name, tmp)
        raw = _raw_plain_bytes(name)
        want = io.TextIOWrapper(io.BytesIO(raw), encoding="utf-8", newline=None).read()
        a = read_log_text(opp, errors="strict")
        b = read_log_text(real, errors="strict")
        if a != b or b != want:
            la, lb, lw = a.splitlines(True), b.splitlines(True), want.splitlines(True)
            i = next(
                (
                    k
                    for k in range(max(len(la), len(lb), len(lw)))
                    if not (k < len(la) and k < len(lb) and k < len(lw) and la[k] == lb[k] == lw[k])
                ),
                None,
            )
            pytest.fail(
                f"{name}: first differing line {i}: "
                f"{la[i] if i is not None and i < len(la) else None!r} / "
                f"{lb[i] if i is not None and i < len(lb) else None!r} / "
                f"{lw[i] if i is not None and i < len(lw) else None!r}"
            )
        _p0(
            f"{TAG} (Id-b) {name}: tracked {_form(real)}, tmp {_form(opp)}; "
            f"{len(b)} chars, {len(b.splitlines())} lines; equal"
        )


@pytest.mark.parametrize(
    "label,entries,keys",
    [
        ("1.5 T", gen.MEASURED_15T, ["48", "52", "59"]),
        ("3 T", gen.MEASURED_3T, ["58", "60", "61"]),
    ],
)
def test_i997_id_c_custody_either_form(tmp_path, label, entries, keys):
    prices = gen.ledger_prices(entries=entries)
    tmp = tmp_path / "opp"
    for p in prices:
        _opposite(p["log"], tmp)
    real = _custody(prices, logs=LOGS)
    opp = _custody(prices, logs=tmp)
    forms = {p["entry"]: (_form(LOGS / p["log"]), _form(tmp / p["log"])) for p in prices}
    _p0(f"{TAG} (Id-c) {label}: forms (tracked, tmp) {forms}; custody tracked {real}, tmp {opp}")
    assert opp == real
    assert real == {k: [] for k in keys}, real


def test_i997_id_p_other_readers(tmp_path, monkeypatch):
    tmp = tmp_path / "opp"
    # P80 / P81: _parse_c52 on C52_LOG (entry 52, gz on disk now).
    want_c52 = p80._parse_c52()
    monkeypatch.setattr(p80, "C52_LOG", _opposite(p80.C52_LOG.name, tmp))
    got_c52 = p80._parse_c52()
    _p0(
        f"{TAG} (Id-p) _parse_c52: lat {len(want_c52[0])} rows, star {len(want_c52[1])}; "
        f"opposite form equal {got_c52 == want_c52}"
    )
    assert got_c52 == want_c52
    assert want_c52[0] and want_c52[1]
    # Closed-form page: read_one(…, STATUS_RE) on all four, and expected().
    names = ("P88_LOG", "P89_LOG", "P83_LOG", "C39_LOG")
    want_st = [cf.read_one(getattr(cf, n), cf.STATUS_RE) for n in names]
    want_exp = cf.expected()
    for n in names:
        monkeypatch.setattr(cf, n, _opposite(getattr(cf, n).name, tmp))
    got_st = [cf.read_one(getattr(cf, n), cf.STATUS_RE) for n in names]
    got_exp = cf.expected()
    _p0(
        f"{TAG} (Id-p) closed-form page: status {[g for g, _ in got_st]}; "
        f"expected() equal {got_exp == want_exp}"
    )
    assert got_st == want_st
    assert got_exp == want_exp
    # Lattice module: _c39().
    want_c39 = lat._c39()
    monkeypatch.setattr(lat, "C39_LOG", _opposite(lat.C39_LOG.name, tmp))
    got_c39 = lat._c39()
    _p0(f"{TAG} (Id-p) _c39: keys {sorted(want_c39)}; opposite form equal {got_c39 == want_c39}")
    assert got_c39 == want_c39


def test_i997_n_absent_raises():
    row = [{"entry": "x", "log": "20990101T000000Z_ABSENT.log", "seconds": "1", "gib": "1.00"}]
    with pytest.raises(FileNotFoundError) as ei:
        _custody(row)
    msg = str(ei.value)
    _p0(f"{TAG} (N-absent) {type(ei.value).__name__}: {msg}")
    assert "20990101T000000Z_ABSENT.log " in msg + " "
    assert "20990101T000000Z_ABSENT.log.gz" in msg


def _entry48(tmp: Path):
    prices = [p for p in gen.ledger_prices(entries=gen.MEASURED_15T) if p["entry"] == "48"]
    assert len(prices) == 1, prices
    name = prices[0]["log"]
    gz = LOGS / f"{name}.gz"
    assert gz.exists() and not (LOGS / name).exists(), "entry 48 is expected gz-only on disk"
    tmp.mkdir(parents=True, exist_ok=True)
    return prices, name, gz.read_bytes()


def test_i997_n_trunc_raises(tmp_path):
    tmp = tmp_path / "trunc"
    prices, name, data = _entry48(tmp)
    (tmp / f"{name}.gz").write_bytes(data[: len(data) // 2])
    with pytest.raises((EOFError, OSError, zlib.error)) as ei:
        _custody(prices, logs=tmp)
    _p0(
        f"{TAG} (N-trunc) {len(data) // 2} of {len(data)} bytes: "
        f"{type(ei.value).__name__}: {ei.value}"
    )


def test_i997_n_crc_raises(tmp_path):
    tmp = tmp_path / "crc"
    prices, name, data = _entry48(tmp)
    bad = bytearray(data)
    bad[-6] ^= 0xFF
    (tmp / f"{name}.gz").write_bytes(bytes(bad))
    with pytest.raises(gzip.BadGzipFile) as ei:
        _custody(prices, logs=tmp)
    _p0(f"{TAG} (N-crc) byte -6 flipped: {type(ei.value).__name__}: {ei.value}")
