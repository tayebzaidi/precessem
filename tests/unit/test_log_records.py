"""`OPS-73` (§9 item 724): the retention-safe log reader.

Anchors (asserted):
  (H)  `read_log_text` returns the identical text for a plain file and for a
       ``.log.gz`` of the same bytes.
  (A)  with neither present it raises `FileNotFoundError` naming both paths.
Negative control (asserted): a bare `Path.read_text` on a ``.gz``-only path
raises `FileNotFoundError`, so the test really is on the state that was red.
"""

from __future__ import annotations

import gzip

import pytest

from tests.log_records import glob_log_paths, read_log_text, resolve_log_path

TEXT = "[tag] first line\n- Status: 0 é\n\nlast line without a newline"


def _plain(tmp_path):
    p = tmp_path / "plain" / "20261101T000000Z_X.log"
    p.parent.mkdir()
    p.write_bytes(TEXT.encode())
    return p


def _gz_only(tmp_path):
    p = tmp_path / "gz" / "20261101T000000Z_X.log"
    p.parent.mkdir()
    with gzip.open(f"{p}.gz", "wb") as fh:
        fh.write(TEXT.encode())
    return p


def test_h_plain_and_gz_read_identically(tmp_path):
    plain, gz = _plain(tmp_path), _gz_only(tmp_path)
    assert not gz.exists()
    a, b = read_log_text(plain), read_log_text(gz)
    print(f"\n[OPS-73] (H) plain {len(a)} chars, .gz {len(b)} chars, identical: {a == b}")
    assert a == b == TEXT
    assert a.splitlines() == b.splitlines() == TEXT.splitlines()
    assert read_log_text(f"{gz}.gz") == TEXT
    assert read_log_text(str(gz)) == TEXT


def test_h_plain_wins_when_both_exist(tmp_path):
    p = _plain(tmp_path)
    with gzip.open(f"{p}.gz", "wb") as fh:
        fh.write(b"not this one")
    assert resolve_log_path(p) == p
    assert read_log_text(p) == TEXT


def test_a_both_absent_refused_naming_both(tmp_path):
    p = tmp_path / "absent.log"
    with pytest.raises(FileNotFoundError) as exc:
        read_log_text(p)
    print(f"\n[OPS-73] (A) {exc.value}")
    assert str(p) in str(exc.value) and f"{p}.gz" in str(exc.value)


def test_nc_bare_read_fails_on_gz_only(tmp_path):
    p = _gz_only(tmp_path)
    with pytest.raises(FileNotFoundError):
        p.read_text()


def test_glob_finds_either_form_once(tmp_path):
    p = _plain(tmp_path)
    with gzip.open(tmp_path / "plain" / "20261102T000000Z_X.log.gz", "wb") as fh:
        fh.write(TEXT.encode())
    with gzip.open(f"{p}.gz", "wb") as fh:
        fh.write(TEXT.encode())
    found = glob_log_paths(tmp_path / "plain" / "2026*_X.log")
    assert found == [str(p), str(tmp_path / "plain" / "20261102T000000Z_X.log")]
    assert glob_log_paths(tmp_path / "plain" / "1999*_X.log") == []
