#!/usr/bin/env bash
# Scheduled XL window (operator directive 2026-09-09). Cron, 02:00 local.
#
# The box is shared with work this sandbox cannot see — on 2026-09-09 the
# Windows host was at 100% CPU with two solves while this VM read load 0.84 and
# every container 0.00%. Rather than try to detect that (which cannot be done
# from here), XL windows run at 02:00 when the machine is quiet.
#
# **No Claude session.** An XL window is one deterministic command, and running
# it from an agent is what made it hard: a foreground harness call is capped at
# 660 s and an implementer slot is killed at 65 min, so no scheduled *session*
# can hold a 2 h (since 2026-09-13, 4 h) window. A plain cron script has neither limit. That is the
# whole reason this file exists.
#
# What runs is the lexically first `*.env` file in docs/testing/<tier>-queue.d/
# (operator directive 2026-09-15: a multi-entry FIFO, so the daily review can
# queue every READY entry of docs/testing/xl-pending.md at once and the
# Tuesday / Thursday windows — nights with no review to refill a one-slot
# file — are no longer lost). The legacy one-slot docs/testing/<tier>-queue.env
# is still honoured, first, when non-empty. Empty queue means nothing runs and
# this exits quietly — so the entry can sit in cron every night and cost
# nothing. The launcher checks the tier's ledger budget BEFORE taking an
# entry, so a spent budget leaves the queue intact for the night it reopens.
set -uo pipefail

# Which tier this invocation runs. `xl` = 2 h / 512 GiB / 3 per week (Sun-Fri
# 02:00); `xxl` = 8 h / 754 GiB / 1 per week (Saturday 02:00). Operator
# directive 2026-09-10. Each tier has its own queue file, its own ledger and
# its own compose service; the guard enforces the ceilings and the budget.
TIER="${1:-xl}"
case "$TIER" in
  xl)  SERVICE="fem-em-solver-xl"  ;;
  xxl) SERVICE="fem-em-solver-xxl" ;;
  *)   echo "usage: xl-run.sh [xl|xxl]" >&2; exit 2 ;;
esac

# Overridable so scripts/testing/test_xl_window_isolation.sh can run this script
# end to end in a scratch repo; cron sets nothing.
REPO="${FEM_EM_REPO:-/home/user/fem-em-solver}"
# Overridable because $HOME is read-only inside the agent sandbox; cron has a
# writable one. `housekeeping.sh` carries the same escape hatch.
LOCK="${FEM_EM_XL_LOCK:-$HOME/.fem-em-$TIER.lock}"   # its own lock, NOT the automation flock:
                                      # a 2 h XL window must not starve the
                                      # 02:15 weekly or the 03:00 review, both
                                      # of which are documentation-only.
LOGDIR="$REPO/logs/automation"
QUEUE="$REPO/docs/testing/$TIER-queue.env"   # legacy one-slot file; moved from scripts/automation/ 2026-09-13: Edit(scripts/automation/**) is on the ask list, so a headless daily review (the XL clerk, daily-review.md step 6b) could not write it there
QUEUE_DIR="$REPO/docs/testing/$TIER-queue.d"  # FIFO: NN-<chunk>.env files, same XL_CHUNK / XL_COMMAND format (2026-09-15)
GUARD_DIR="$REPO/scripts/automation/hooks"
export PATH="$HOME/.local/bin:/usr/local/bin:/usr/bin:/bin"

mkdir -p "$LOGDIR"
LOG="$LOGDIR/$(date -u +%Y%m%dT%H%M%SZ)_${TIER}-run.log"
exec >>"$LOG" 2>&1
echo "$(date -u) xl-run starting: tier=$TIER service=$SERVICE"

# Distinguish "cannot create the lock" from "another window holds it". They had
# the same message and the same exit 0, so an unwritable lock path looked
# exactly like a healthy skip — which is how a failed launch would hide.
if ! exec 9>"$LOCK" 2>/dev/null; then
  echo "$(date -u) FAILED: cannot open lock file $LOCK (set FEM_EM_XL_LOCK to a writable path)"
  exit 1
fi
if ! flock -n 9; then
  echo "$(date -u) another $TIER window holds the lock; skipping"
  exit 0
fi

# Pick the entry: the legacy one-slot file if it is non-empty, else the first
# file in the FIFO directory. ENTRY is what gets consumed afterwards.
ENTRY=""
XL_CHUNK=""; XL_COMMAND=""; XL_OUTPUTS=""
if [[ -f "$QUEUE" ]]; then
  # shellcheck source=/dev/null
  source "$QUEUE" || true
  [[ -n "${XL_CHUNK:-}" && -n "${XL_COMMAND:-}" ]] && ENTRY="$QUEUE"
fi
if [[ -z "$ENTRY" && -d "$QUEUE_DIR" ]]; then
  for f in "$QUEUE_DIR"/*.env; do
    [[ -f "$f" ]] || continue
    XL_CHUNK=""; XL_COMMAND=""; XL_OUTPUTS=""
    # shellcheck source=/dev/null
    source "$f" || true
    if [[ -n "${XL_CHUNK:-}" && -n "${XL_COMMAND:-}" ]]; then ENTRY="$f"; break; fi
    echo "$(date -u) skipping malformed queue entry $(basename "$f") (XL_CHUNK/XL_COMMAND unset)"
  done
fi
CHUNK="${XL_CHUNK:-}"
CMD="${XL_COMMAND:-}"
# Optional (OPS-64, 2026-09-28): space-separated repo-relative paths the
# window writes into the tree and that belong in its commit (a regenerated
# metrics file, say). Undeclared tree outputs are left alone and named in this
# log — and then make the next window refuse on a dirty tree, which is why a
# queue entry that writes into the tree must declare what it writes.
OUTPUTS="${XL_OUTPUTS:-}"
if [[ -z "$ENTRY" ]]; then
  echo "$(date -u) queue empty; nothing to run"
  exit 0
fi
QUEUED_AHEAD=$(( $(ls "$QUEUE_DIR"/*.env 2>/dev/null | wc -l) ))
echo "$(date -u) queue entry $(basename "$ENTRY") -> chunk $CHUNK ($QUEUED_AHEAD file(s) in $(basename "$QUEUE_DIR"))"

# Budget preflight (same arithmetic as bash_guard.py: charged rows in the
# trailing 7 days of the tier's ledger). A spent budget must NOT consume the
# entry — the guard inside run_and_log.sh would refuse the command anyway,
# and clearing the queue on a refusal is how a window gets silently lost.
USED="$(cd "$GUARD_DIR" && python3 -c "import bash_guard as g; t=g.TIERS['$TIER']; print(g.runs_in_trailing_week(t['ledger']), t['per_week'])" 2>/dev/null)"
if [[ -n "$USED" ]]; then
  set -- $USED
  if [[ "$1" -ge "$2" ]]; then
    echo "$(date -u) $TIER budget spent ($1 of $2 charged rows in the trailing 7 days); entry left in the queue for the night it reopens"
    exit 0
  fi
  echo "$(date -u) $TIER budget: $1 of $2 used in the trailing 7 days"
else
  echo "$(date -u) WARNING: could not read the $TIER budget from bash_guard.py; relying on the harness guard"
fi

cd "$REPO" || exit 1
if [[ -n "$(git status --porcelain)" ]]; then
  echo "$(date -u) tree dirty; refusing to run an XL window on it"
  git status --short
  exit 0
fi
if [[ "$(git rev-parse --abbrev-ref HEAD)" != "main" ]]; then
  echo "$(date -u) not on main; skipping"
  exit 0
fi

# Declared-output preflight (OPS-64), BEFORE the ledger or the queue is
# touched: a declared output this script could not commit is refused now,
# not discovered after a 4 h window. Refused: an absolute path, a path that
# climbs out of the repo, and a gitignored path (`git add` would skip it and
# the window's product would silently stay out of its commit).
for p in $OUTPUTS; do
  WHY=""
  if [[ "$p" == /* ]]; then WHY="an absolute path"
  elif [[ "/$p/" == *"/../"* ]]; then WHY="outside the repo"
  elif git check-ignore -q --no-index -- "$p"; then WHY="gitignored, so it could not be committed"
  fi
  if [[ -n "$WHY" ]]; then
    echo "$(date -u) REFUSED before the window: declared output '$p' is $WHY; entry $(basename "$ENTRY") left in the queue, ledger untouched"
    exit 1
  fi
done
[[ -n "$OUTPUTS" ]] && echo "$(date -u) declared outputs: $OUTPUTS"

# Environment preflight, BEFORE anything that consumes the slot. The harness
# appends the ledger row the moment the window starts, so a failure after that
# point spends the week for nothing. Measured 2026-09-10: run by hand from an
# agent session, docker was permission-denied from inside this script (it works
# from a direct agent command, not from a script), the window died in 0 s, and
# it had already taken a ledger row and cleared the queue. Cron is outside that
# sandbox and unaffected — but the check is cheap and the failure was silent.
if ! docker compose -f docker/docker-compose.yml ps >/dev/null 2>&1; then
  echo "$(date -u) FAILED: docker unreachable from this context; queue and ledger untouched"
  exit 1
fi

# §5.1: restart the service immediately before the window so the container's
# `memory.peak` belongs to this run — it is per container lifetime and cannot
# be reset on this kernel.
docker compose -f docker/docker-compose.yml --profile "$TIER" up -d --force-recreate "$SERVICE"
sleep 5

# Orphan count in the container's namespace (OPS-64), taken by the launcher
# itself so it does not depend on the queued command carrying it (CLAUDE.md:
# a dead wrapper does not stop the compute — measured 2026-09-09). Prints an
# integer, or nothing when the count cannot be read.
count_python() {
  docker compose -f docker/docker-compose.yml --profile "$TIER" exec -T "$SERVICE" \
    bash -c 'pgrep -c python3; true' 2>/dev/null | tr -dc '0-9\n' | grep -m1 -E '^[0-9]+$'
}
ORPH_BEFORE="$(count_python)"
echo "$(date -u) container python3 processes before the window: ${ORPH_BEFORE:-unreadable}"

echo "$(date -u) running $TIER chunk $CHUNK"
# The box is quiet by schedule, so the guard's interactive box-confirmation is
# satisfied here by the clock rather than by a human at 02:00.
FEM_EM_XL_BOX_OK=1 scripts/testing/run_and_log.sh "$CHUNK" "$CMD"
STATUS=$?
echo "$(date -u) run_and_log exit=$STATUS"

ORPH_AFTER="$(count_python)"
docker compose -f docker/docker-compose.yml --profile "$TIER" stop "$SERVICE"
STILL_UP="$(docker compose -f docker/docker-compose.yml --profile "$TIER" ps --status running -q "$SERVICE" 2>/dev/null)"
if [[ -z "$ORPH_BEFORE" || -z "$ORPH_AFTER" ]]; then
  ORPH_VERDICT="not established (count unreadable: before '${ORPH_BEFORE}', after '${ORPH_AFTER}')"
elif [[ "$ORPH_AFTER" -eq 0 ]]; then
  ORPH_VERDICT="none"
else
  ORPH_VERDICT="ORPHANS: $ORPH_AFTER python3 process(es) outlived the window's command; the service stop below ends them"
fi
echo "$(date -u) orphans: $ORPH_VERDICT"

# Append the launcher's own record to the harness log this window commits
# (the newest untracked log for this chunk), after the harness's ## Exit.
SAFE_CHUNK="${CHUNK//[^a-zA-Z0-9._-]/_}"
WINDOW_LOG="$(git ls-files --others --exclude-standard -- "docs/testing/logs/*_${SAFE_CHUNK}.log" | sort | tail -1)"
if [[ -n "$WINDOW_LOG" ]]; then
  {
    echo ""
    echo "## Launcher (xl-run.sh, OPS-64)"
    echo "- Orphans before (container python3 processes): ${ORPH_BEFORE:-unreadable}"
    echo "- Orphans after (container python3 processes): ${ORPH_AFTER:-unreadable}"
    echo "- Orphans: $ORPH_VERDICT"
    if [[ -z "$STILL_UP" ]]; then echo "- Service stopped: yes"; else echo "- Service stopped: NO (still running: $STILL_UP)"; fi
    echo "- Declared outputs: ${OUTPUTS:-none}"
  } >> "$WINDOW_LOG"
else
  echo "$(date -u) WARNING: no harness log for $CHUNK found; launcher record kept in this log only"
fi

# Consume the queue entry either way: a failed XL window has still spent its
# slot (§5.1), and re-running it unattended tomorrow is exactly what the
# once-per-interval rule exists to prevent. Re-queue deliberately or not at all.
if [[ "$ENTRY" == "$QUEUE" ]]; then
  : > "$QUEUE"
  echo "$(date -u) legacy queue file cleared"
else
  git rm -q --cached "$ENTRY" 2>/dev/null || true
  rm -f "$ENTRY"
  echo "$(date -u) queue entry $(basename "$ENTRY") consumed; $(ls "$QUEUE_DIR"/*.env 2>/dev/null | wc -l) left"
fi

# Stage ONLY what this window owns (2026-09-19; was `git add -A`). The window
# holds its own lock, not the automation flock, so the 03:00 review and the
# implementer slots run while it does, and `git add -A` at the end of a 4 h /
# 8 h window would have committed whatever they had in flight under this
# script's name. Every window since 09-16 has committed exactly these paths:
# its harness log(s), the tier ledger, the test-results row, and the queue
# entry's removal (already staged above) — plus, since OPS-64, the paths the
# queue entry declares in XL_OUTPUTS. The commit is path-limited too, so
# another session's *staged* files stay out of it. test-results.md is staged
# whole: a concurrent session's uncommitted row rides along, which is an index
# line, not work. Anything else dirty is left alone and named in this log.
OWN=("docs/testing/$TIER-ledger.md" docs/testing/test-results.md "${ENTRY#"$REPO"/}")
while IFS= read -r f; do
  [[ -n "$f" ]] && OWN+=("$f")
done < <(git ls-files --others --exclude-standard -- "docs/testing/logs/*_${SAFE_CHUNK}.log")
for p in $OUTPUTS; do
  if [[ -e "$p" ]]; then OWN+=("$p"); else echo "$(date -u) NOTE: declared output $p was not written"; fi
done
for f in "${OWN[@]}"; do
  [[ -e "$f" ]] && git add -- "$f"
done
LEFT="$(git status --porcelain | grep -v -F -f <(printf '%s\n' "${OWN[@]}") || true)"
if [[ -n "$LEFT" ]]; then
  echo "$(date -u) NOTE: paths not owned by this window left out of its commit (another session's work in flight?):"
  echo "$LEFT"
fi

# The tree was on main at 02:00; a session may have moved it since (implementer
# runs park work on attempt/*). Never commit a window onto somebody's branch.
NOW_ON="$(git rev-parse --abbrev-ref HEAD)"
if [[ "$NOW_ON" != "main" ]]; then
  echo "$(date -u) FAILED: HEAD is on $NOW_ON, not main; the window's record is left staged, uncommitted"
elif git diff --cached --quiet -- "${OWN[@]}"; then
  echo "$(date -u) nothing to commit"
else
  # A session committing at the same moment holds .git/index.lock; retry
  # rather than leave the window's record uncommitted on main.
  COMMITTED=0
  for attempt in 1 2 3 4 5; do
    if git -c user.name="fem-em xl-run" -c user.email="xl-run@localhost" commit -q -m "chore($TIER): scheduled $TIER window — $CHUNK $(date -u +%Y-%m-%d)

Ran by scripts/automation/xl-run.sh $TIER at the 02:00 slot. Harness log, ledger row
and test-results row are in this commit; the readout still needs a human or a
review to interpret it, and the ledger's last four columns are filled by hand." -- "${OWN[@]}"; then
      COMMITTED=1
      break
    fi
    echo "$(date -u) commit attempt $attempt failed; retrying in 20 s"
    sleep 20
  done
  if [[ "$COMMITTED" == 1 ]]; then
    echo "$(date -u) committed $(git rev-parse --short HEAD)"
  else
    echo "$(date -u) FAILED: could not commit the window's record; it is staged on main"
  fi
fi
exit "$STATUS"
