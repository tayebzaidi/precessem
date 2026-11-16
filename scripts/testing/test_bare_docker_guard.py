#!/usr/bin/env python3
"""Asserting cases for `bash_guard.py` rule 3 (bare docker / harness calls).

The DENY cases are the command shapes that actually lost windows on
2026-09-29/30 (slot sessions and spawned executors: `permission denied …
docker.sock`, 0-1 s); the ALLOW cases are the bare calls that ran green, and
mentions that must never trip the rule. Run as a file — the cases live here so
this command line does not itself trip the guard:

    python3 scripts/testing/test_bare_docker_guard.py
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
GUARD = ROOT / "scripts" / "automation" / "hooks" / "bash_guard.py"

H = 'scripts/testing/run_and_log.sh X "docker compose exec -T fem-em-solver bash -lc \'cd /workspace && timeout -k 30 60 python3 a.py\'"'

CASES = [
    # --- must DENY: the shapes that lost windows ---
    ("DENY", "cd prefix (EX-66 executor, 13:30)", "cd /home/u/fem-em-solver; " + H),
    ("DENY", "cd && prefix", "cd /home/u/fem-em-solver && " + H),
    ("DENY", "; date suffix (EX-64 executor, 06:00)", H + "; date"),
    ("DENY", "; tail suffix", H + "; tail -6 docs/testing/logs/x.log | cut -c1-200"),
    ("DENY", "date -u; prefix (TH-5 13:30)", "date -u; " + H),
    ("DENY", "harness piped to tail", H + " | tail -20"),
    ("DENY", "docker ps piped to head", "docker compose -f docker/docker-compose.yml ps 2>&1 | head -5"),
    ("DENY", "docker with redirect", "docker build -f docker/palace.Dockerfile -t t docker > logs/b.log 2>&1"),
    ("DENY", "git status && docker ps (preflight)", "git status --short && docker compose -f docker/docker-compose.yml ps"),
    ("DENY", "absolute harness path", "/home/u/fem-em-solver/" + H),
    ("DENY", "./ harness path", "./" + H),
    ("DENY", "env assignment prefix", "FOO=1 " + H),
    # --- must ALLOW: bare calls and mere mentions ---
    ("ALLOW", "bare harness call", H),
    ("ALLOW", "bare docker compose ps", "docker compose -f docker/docker-compose.yml ps"),
    ("ALLOW", "bare docker build", "docker build --progress=plain -f docker/palace.Dockerfile -t t docker"),
    ("ALLOW", "grep mentioning docker", 'grep -n "docker" docs/testing/attempts.md | head'),
    ("ALLOW", "grep for run_and_log.sh", "grep -rn run_and_log.sh scripts | head -3"),
    ("ALLOW", "commit message mentioning docker", 'git commit -m "docker: note; harness | fine"'),
    ("ALLOW", "plain compound without docker", "cd /tmp && ls | wc -l; date"),
    ("ALLOW", "unbalanced quote fails open", "echo 'oops"),
]


def verdict(cmd: str) -> tuple[str, str]:
    r = subprocess.run([sys.executable, str(GUARD)],
                       input=json.dumps({"tool_input": {"command": cmd}}),
                       capture_output=True, text=True)
    if not r.stdout.strip():
        return "ALLOW", ""
    return "DENY", json.loads(r.stdout)["hookSpecificOutput"]["permissionDecisionReason"]


def main() -> int:
    bad = 0
    for want, label, cmd in CASES:
        got, reason = verdict(cmd)
        ok = got == want
        bad += not ok
        print(f"{'ok  ' if ok else 'FAIL'} {want:5} {label:40} -> {got}"
              + (f"  [{reason[reason.find('because'):][:70]}]" if got == "DENY" else ""))
    print(f"\n{len(CASES) - bad}/{len(CASES)} cases as expected")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
