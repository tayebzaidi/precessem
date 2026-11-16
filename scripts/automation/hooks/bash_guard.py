#!/usr/bin/env python3
"""PreToolUse guard for Bash commands (wired in .claude/settings.json).

Enforces the CLAUDE.md compute hard rules mechanically, for every session
including the headless cron runs:

  1. Rank ceiling: any `mpiexec/mpirun -n/-np N` with N > 12 is denied —
     except the XL tier (operator directive 2026-09-05): N <= 16 is allowed
     when the command is a `docker compose … exec` against the
     `fem-em-solver-xl` service (a mention of the name alone is not a run), goes through
     `scripts/testing/run_and_log.sh`, carries a `timeout -k 30 <s>` of at
     most 7200 s, and `docs/testing/xl-ledger.md` shows no XL run in the last
     7 days. The same four conditions gate *any* command against the XL
     service, whatever its rank count.
  2. Harness routing: pytest invocations must go through
     scripts/testing/run_and_log.sh (collection-only queries are exempt).
  3. Bare docker / harness calls (2026-09-30, interactive session): a command
     that runs `docker` or `run_and_log.sh` must be exactly that one call —
     first word `docker` or `scripts/testing/run_and_log.sh`, no `cd`
     prefix, no `;` `&&` `||` `|` `&`, no redirect. `.claude/settings.json`
     exempts only whole commands matching `docker *` /
     `scripts/testing/run_and_log.sh *` from the sandbox; anything else runs
     sandboxed, where /var/run/docker.sock is unreachable, and the window
     dies in 0-1 s with `permission denied … docker.sock` (every such
     failure on 09-29/09-30, spawned executors included, had a prefix, a
     `; date` / `; tail` suffix or a pipe). Denying up front turns a wasted,
     logged window into an immediate, explained retry.

Fail-open on malformed input: a broken guard must not brick every Bash call.
"""
import datetime as dt
import json
import os
import re
import shlex
import sys
from pathlib import Path

RANK_CEILING = 12
# Two big-compute tiers, gated the same way and budgeted separately
# (operator directive 2026-09-10). `xl` moved from one run per 7 days to a
# budget of three, because the constraint turned out to be wall-clock and not
# memory; `xxl` is the Saturday window for what a 2 h slot cannot hold.
_DOCS = Path(__file__).resolve().parents[3] / "docs" / "testing"
TIERS = {
    "xl": {
        "service": "fem-em-solver-xl",
        "ranks": 16,
        "timeout_s": 14400,         # 4 h (operator directive 2026-09-13; was 2 h — ANS-4 step 2d ran 25 s past 7200 and finished only because the deadline was lifted by hand)
        "per_week": 6,              # trailing 7 days (operator directive 2026-09-13: nightly Sun-Fri; was 3 from 09-10, 1 from 09-05)
        "ledger": _DOCS / "xl-ledger.md",
    },
    "xxl": {
        "service": "fem-em-solver-xxl",
        "ranks": 16,
        "timeout_s": 28800,         # 8 h
        "per_week": 1,
        "ledger": _DOCS / "xxl-ledger.md",
    },
}
# Longest service name first: "fem-em-solver-xl" is not a substring of
# "fem-em-solver-xxl", but matching order is still made explicit so a future
# rename cannot silently send an xxl command down the xl path.
TIER_ORDER = ("xxl", "xl")

XL_SERVICE = TIERS["xl"]["service"]          # kept: referenced in messages
XL_RANK_CEILING = TIERS["xl"]["ranks"]
XL_TIMEOUT_CEILING_S = TIERS["xl"]["timeout_s"]
XL_INTERVAL_DAYS = 7
XL_LEDGER = TIERS["xl"]["ledger"]
XL_OVERRIDE = Path(__file__).resolve().parents[1] / "xl-override.env"

# Shared-box courtesy for the XL tier (operator directive 2026-09-09). An XL
# window holds up to 16 of this 36-core box's cores for up to two hours, so
# starting one onto a box somebody else is already using is the single rudest
# thing this project can do. Refuse when the cores are not there and say when
# to come back; the run is postponed, never silently degraded.
XL_LOAD_HEADROOM_CORES = 4.0
XL_LOAD_OVERRIDE_ENV = "FEM_EM_XL_IGNORE_LOAD"
XL_BOX_OK_ENV = "FEM_EM_XL_BOX_OK"


def deny(reason: str) -> None:
    print(json.dumps({
        "hookSpecificOutput": {
            "hookEventName": "PreToolUse",
            "permissionDecision": "deny",
            "permissionDecisionReason": reason,
        }
    }))
    sys.exit(0)


def ledger_dates(ledger: Path) -> "list[dt.date]":
    """Every run date in a tier's ledger (`| YYYY-MM-DD | ...`)."""
    try:
        text = ledger.read_text(encoding="utf-8")
    except OSError:
        return []
    return [dt.date.fromisoformat(m)
            for m in re.findall(r"^\|\s*(\d{4}-\d{2}-\d{2})\s*\|", text, flags=re.M)]


def runs_in_trailing_week(ledger: Path) -> int:
    """Rows inside the trailing 7 days that actually consumed the box.

    The budget is on **box time**, not on attempts. `run_and_log.sh` writes the
    row when a window starts and fills the elapsed column when it ends, so a
    window that never started carries 0 — and on 2026-09-10 two of five rows
    were exactly that (docker unreachable, 0 s, no compute). Charging those
    against a six-per-week budget would price an infrastructure failure like
    a two-hour solve.

    This cannot be gamed the way "re-run until it works" could: you cannot
    produce a 0-second row for a run that did any work, because the harness
    measures the elapsed time itself.
    """
    cutoff = dt.date.today() - dt.timedelta(days=XL_INTERVAL_DAYS - 1)
    try:
        text = ledger.read_text(encoding="utf-8")
    except OSError:
        return 0
    n = 0
    for row in re.findall(r"^\|\s*(\d{4}-\d{2}-\d{2})\s*\|(.*)$", text, flags=re.M):
        date_s, rest = row
        try:
            when = dt.date.fromisoformat(date_s)
        except ValueError:
            continue
        if when < cutoff:
            continue
        cells = [c.strip() for c in rest.split("|")]
        # date | chunk | log | ranks | cells | peak | elapsed | readout
        elapsed = cells[5] if len(cells) > 5 else ""
        if elapsed in ("", "-", "\u2014") or elapsed == "0":
            continue          # never started: no box time spent
        n += 1
    return n


def last_xl_run() -> "dt.date | None":
    """Date of the most recent `xl` run. Kept for callers and tests."""
    dates = ledger_dates(XL_LEDGER)
    return max(dates) if dates else None


def xl_interval_override() -> "str | None":
    """The reason string of an in-date operator override, else None.

    The 7-day interval is a budget the weekly review spends, and normally
    nothing may buy a second run inside it.  This is the one documented door:
    an **operator-authorised, dated, self-expiring** file that says who
    authorised the exception and until when.  It expires by date here, so an
    override cannot outlive its reason, and it never touches the ledger — the
    record of what ran stays complete and a spent slot stays spent.

    Malformed or absent file ⇒ no override, and the ordinary interval rule
    applies.  Fail-closed on purpose: the safe direction for a rule that
    protects a shared box is to deny.
    """
    try:
        text = XL_OVERRIDE.read_text(encoding="utf-8")
    except OSError:
        return None
    fields = {}
    for line in text.splitlines():
        line = line.split("#", 1)[0].strip()
        if "=" in line:
            k, v = line.split("=", 1)
            fields[k.strip()] = v.strip().strip('"').strip("'")
    until = fields.get("XL_OVERRIDE_UNTIL", "")
    reason = fields.get("XL_OVERRIDE_REASON", "") or "no reason recorded"
    try:
        expiry = dt.date.fromisoformat(until)
    except ValueError:
        return None
    if dt.date.today() > expiry:
        return None
    return f"{reason} (operator override, expires {expiry.isoformat()})"


def xl_load_check(cmd: str) -> "str | None":
    """Reason to postpone an XL window on box load, or None to allow.

    Uses the 1-minute load average against the core count, reserving the ranks
    this command asks for plus a headroom margin. Fails **open** on anything
    unreadable — a guard that cannot measure the box must not become the reason
    no work can run.

    Honest limitation, worth knowing before trusting a refusal: load average is
    a decaying mean and does not distinguish our own load from anyone else's,
    so for a minute or two after one of our own runs ends it will read high and
    postpone the next one. That is the safe direction, and
    ``FEM_EM_XL_IGNORE_LOAD=1`` is the operator's escape hatch when the load is
    known to be our own and decaying.
    """
    raw = os.environ.get(XL_LOAD_OVERRIDE_ENV, "").strip().lower()
    if raw and raw not in ("0", "false", "no", "off"):
        return None
    try:
        with open("/proc/loadavg", encoding="utf-8") as fh:
            load1 = float(fh.read().split()[0])
        ncpu = float(os.cpu_count() or 0)
        if ncpu <= 0:
            return None
    except Exception:
        return None
    m = re.search(r"mpi(?:exec|run)\s+(?:-n|-np)\s+(\d+)", cmd)
    ranks = float(m.group(1)) if m else float(XL_RANK_CEILING)
    free = ncpu - load1
    needed = ranks + XL_LOAD_HEADROOM_CORES
    if free < needed:
        return (
            f"the box is busy: 1-min load {load1:.1f} of {ncpu:.0f} cores leaves "
            f"{free:.1f} free, and this window wants {ranks:.0f} ranks plus "
            f"{XL_LOAD_HEADROOM_CORES:.0f} cores of headroom ({needed:.0f}). An XL "
            "window holds its cores for up to two hours on a shared machine, so it "
            "waits rather than oversubscribing. Re-check `cat /proc/loadavg` and "
            f"start when the load is under {ncpu - needed:.1f}. If the load is our "
            f"own and still decaying, {XL_LOAD_OVERRIDE_ENV}=1 overrides — operator "
            "judgement only, and never to push past somebody else's work."
        )
    return None


def check_tier(cmd: str, tier: str) -> None:
    """The tier's four conditions; deny on the first that fails."""
    spec = TIERS[tier]
    name, ceiling = spec["service"], spec["timeout_s"]

    if "run_and_log.sh" not in cmd:
        deny(f"Commands against {name} must go through scripts/testing/run_and_log.sh "
             f"({tier.upper()} tier, PROJECT_PLAN §5.1).")

    timeouts = [int(v) for v in re.findall(r"timeout\s+-k\s+\d+\s+(\d+)", cmd)]
    if not timeouts:
        deny(f"Commands against {name} must carry `timeout -k <n> <s>` with s <= "
             f"{ceiling} ({tier.upper()} tier ceiling {ceiling // 3600} h).")
    if max(timeouts) > ceiling:
        deny(f"timeout {max(timeouts)} s exceeds the {tier.upper()} tier ceiling of {ceiling} s.")

    busy = xl_load_check(cmd)
    if busy is not None:
        deny(f"{tier.upper()} window postponed — {busy}")

    # NOTE (2026-09-09): a mandatory per-run "box confirmed free" token used to
    # live here. Removed the same day — see the scheduling note in §5.1.

    used = runs_in_trailing_week(spec["ledger"])
    budget = spec["per_week"]
    if used >= budget and xl_interval_override() is None:
        recent = sorted(set(ledger_dates(spec["ledger"])))[-budget:]
        deny(f"The {tier.upper()} budget for the trailing {XL_INTERVAL_DAYS} days is spent: "
             f"{used} of {budget} runs already recorded ({', '.join(d.isoformat() for d in recent)}, "
             f"{spec['ledger'].name}). Operator directive 2026-09-10. Wait for the oldest to age "
             f"out, shrink the case to a tier that fits, or have an operator authorise an "
             f"exception in {XL_OVERRIDE.name} (dated, self-expiring); a scheduled session may "
             "not write one.")


def check_xl(cmd: str) -> None:
    """Back-compatible entry point: the `xl` tier's conditions."""
    check_tier(cmd, "xl")


SANDBOX_EXEMPT_FIRST_WORDS = ("docker", "scripts/testing/run_and_log.sh")
_ASSIGNMENT = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*=")


def unsandboxed_call_problem(cmd: str) -> "str | None":
    """Rule 3: why a docker / harness call would run sandboxed, or None if it is bare.

    Only command-position words count (the first word of the line and of each
    part after a shell operator, past any VAR=value assignments), so a quoted
    mention — a grep pattern, a commit message — never trips it.
    """
    try:
        lex = shlex.shlex(cmd, posix=True, punctuation_chars=True)
        lex.whitespace_split = True
        tokens = list(lex)
    except ValueError:
        return None                                   # unbalanced quotes: fail open
    if not tokens:
        return None
    operators = [t for t in tokens if t and set(t) <= set(";&|<>()")]
    command_words, at_start = [], True
    for t in tokens:
        if t in operators:
            at_start = True
            continue
        if at_start and _ASSIGNMENT.match(t):
            continue
        if at_start:
            command_words.append(t)
            at_start = False
    runs_docker = any(w == "docker" or w.endswith("run_and_log.sh") for w in command_words)
    if not runs_docker:
        return None
    first = command_words[0]
    if operators:
        return f"it chains shell operators {sorted(set(operators))}"
    if first not in SANDBOX_EXEMPT_FIRST_WORDS:
        return f"its first word is `{first}`, not `docker` or `scripts/testing/run_and_log.sh`"
    if _ASSIGNMENT.match(tokens[0]):
        return "it starts with a VAR=value assignment"
    return None


def main() -> None:
    try:
        payload = json.load(sys.stdin)
        cmd = payload.get("tool_input", {}).get("command", "") or ""
    except Exception:
        sys.exit(0)  # fail open

    # Rule 3: docker and the harness only work as bare, whole commands.
    problem = unsandboxed_call_problem(cmd)
    if problem:
        deny(
            "This docker / run_and_log.sh call would run inside the sandbox, where "
            "/var/run/docker.sock is unreachable, and fail with `permission denied` — "
            f"because {problem}. Only a whole command starting with `docker` or "
            "`scripts/testing/run_and_log.sh` (repo-relative, from the repo root) is "
            "sandbox-exempt. Re-issue it bare: no `cd … &&` prefix, no `; date` / "
            "`; tail` suffix, no pipe or redirect. Read the log in a separate call."
        )

    # "Targets the XL service" means an actual exec against it (a compute
    # command), not any mention of its name — heredocs, greps and commit
    # messages that name the service are not XL runs. `--profile xl up/stop`
    # is service lifecycle, not compute, and is not gated either.
    targets_tier = None
    for _tier in TIER_ORDER:
        if re.search(r"docker\s+compose[^;&|]*\bexec\b[^;&|]*" + re.escape(TIERS[_tier]["service"]), cmd):
            targets_tier = _tier
            break
    targets_xl = targets_tier is not None

    # Rule 1: rank ceiling, absolute — applies even inside harness commands.
    for m in re.finditer(r"mpi(?:exec|run)\s+(?:-n|-np)\s+(\d+)", cmd):
        n = int(m.group(1))
        if n > RANK_CEILING and not targets_xl:
            deny(
                f"mpiexec -n {n} exceeds the {RANK_CEILING}-rank ceiling for this "
                "project (CLAUDE.md hard rules; shared 36-core box). Use the "
                f"smallest rank count that fits the tier. (Up to {XL_RANK_CEILING} ranks "
                f"exist only in the weekly XL slot against {XL_SERVICE} — PROJECT_PLAN §5.1.)"
            )
        tier_ranks = TIERS[targets_tier]["ranks"] if targets_tier else XL_RANK_CEILING
        if n > tier_ranks:
            deny(f"mpiexec -n {n} exceeds the {(targets_tier or 'xl').upper()} tier's "
                 f"{tier_ranks}-rank ceiling.")

    if targets_tier is not None:
        check_tier(cmd, targets_tier)

    # Rule 2: pytest runs go through the logging harness. Matching is scoped
    # to actual invocations (container exec, mpiexec, or a host-side pytest at
    # a command boundary) so grep/git/commit-message mentions never trip it.
    if "run_and_log.sh" not in cmd:
        collect_only = "--collect-only" in cmd or "--co" in cmd or "--version" in cmd
        in_container = "docker compose" in cmd and "exec" in cmd and re.search(r"\bpytest\b", cmd)
        host_side = re.search(r"(?:^|[;&|]\s*)(?:timeout\s+\S+\s+)?(?:mpi(?:exec|run)[^;&|]*\s)?(?:python3?\s+-m\s+)?pytest\b", cmd)
        if (in_container or host_side) and not collect_only:
            deny(
                "pytest must run through the logging harness: "
                "scripts/testing/run_and_log.sh <CHUNK-ID> \"docker compose exec -T "
                "fem-em-solver bash -lc '...'\" (CLAUDE.md hard rules). "
                "Collection-only checks (--collect-only) are exempt."
            )

    sys.exit(0)


if __name__ == "__main__":
    main()
