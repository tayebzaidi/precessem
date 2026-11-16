# Implementer run protocol (Opus, scheduled)

Run by `scripts/automation/implementer-run.sh` via cron, **twelve times a
day, every day** — four after each daily review: 04:30 / 06:00 / 07:30 /
09:00 (after 03:00), 12:00 / 13:30 / 15:00 / 16:30 (after 10:30), 19:30 /
21:00 / 22:30 / 00:00 (after 18:00); no 21:00 run on Monday, Wednesday and Saturday,
the weekly review's window (full schedule restored 2026-09-28, operator
directive; canonical crontab `scripts/automation/crontab`). One
run = one attempt at the **top "On deck" item** in PROJECT_PLAN.md §9, then
the next if it commits before minute 30 (step 2). Sessions never overlap:
the grid spacing exceeds the 65-minute hard kill, and all runs share one
`flock`.
Since 2026-09-12 (operator directive, enacted by the 18:00 daily review) a
run that finishes its item early takes the next one — see step 2.

## Timebox

60 minutes wall clock, externally enforced (the wrapper kills the session at
65). Note the start time immediately; **start no new implementation work
after minute 45** — the final 15 minutes are for documentation and leaving a
clean tree. **A second (or later) item may be started only before minute
30**, and only from a clean, committed tree (step 2); the minute-45 rule and
every per-command ceiling apply to it unchanged. The per-command compute budget is unchanged and non-negotiable:
§5.1 tiers, 20-minute hard ceiling per compute command, `mpiexec -n 12` max,
`timeout` at the tier ceiling, shared machine. The one exception is an on-deck
item the weekly review marked **`xl`** (§5.1 XL tier: ≤ 512 GiB / 16 ranks /
14400 s against `fem-em-solver-xl`, six per week since 2026-09-13, ledger-gated) — run exactly
the command it specifies, bring the service up before and stop it after, and
never invent an XL run yourself.

## Steps

1. Preflight: `git status` must be clean. If it is not, do **no** chunk work:
   append an `anomaly` entry to `docs/testing/attempts.md` describing what
   you found, commit only that, and stop. Confirm the container is Up per
   CLAUDE.md.

   **Exception — already-journaled documentation drift.** If a *prior* run's
   attempts.md anomaly entry already describes this exact dirty tree, you may
   land it and proceed with chunk work, provided ALL of the following hold
   (verify each; do not assume):
   - the current diff is **byte-identical** to the one the prior entry
     journaled (same files, same `git diff` content);
   - it touches **documentation only** — nothing under `src/`, `tests/`,
     `scripts/`, and no §7 status or done-when change in PROJECT_PLAN.md;
   - the edits are internally consistent and complete (they read as a
     finished change, not half of one).

   Commit the diff **by itself** first (`docs: land dirty tree journaled
   <UTC timestamp of the anomaly entry>`), note in your own attempts.md
   entry that you landed it and why the conditions held, then continue at
   step 2. If any condition fails — new or different dirtiness, anything
   non-documentation, no prior journal entry — the original rule applies
   unchanged: anomaly entry, commit only that, stop. Never discard or stash
   a human's edits; landing an already-journaled doc diff is the only
   permitted action.

   **Second encounter — park it and proceed.** If a prior run's attempts.md
   anomaly entry already journaled a dirty tree and the tree is *still* dirty
   now, the tree is stuck and stopping again just burns the rest of the day's
   slots. Whatever the diff contains:
   - commit it **as-is** to a new branch `recovered/<UTC-timestamp>` (create
     the branch, commit there, return to a clean `main`);
   - journal in your attempts.md entry what you parked, the branch name, and
     the timestamp of the prior entry that made this the second encounter;
   - then continue at step 2 and do chunk work normally.

   **A launcher's refused record — leak-held on first encounter (weekly
   ruling 2026-11-07).** If the dirty tree is exactly a scheduled `xl` /
   `xxl` window's record left staged by its launcher — the new harness log
   under `docs/testing/logs/`, its `test-results.md` row, one ledger row
   and a consumed queue file, nothing else — and the launcher's log
   (`logs/automation/*xl-run.log`) shows its commit refused by the
   `pre-commit` leak hook, do not stop and do not park: apply §9's
   leak-held rule yourself, on the first encounter. Move the harness log
   unchanged to the gitignored `logs/leak-held/`, commit the rest through
   the hook as `chore(<tier>): scheduled <tier> window — <chunk> <date>
   (log held by the leak guard)`, journal which line number the hook named
   (never the token or its digits) and that the daily review reads the
   held copy, then continue at step 2. Never `--no-verify`, never edit the
   checker; a tree that does not match this shape exactly follows the
   rules above. (Four slots stopped on one such record on 2026-11-07.)

   Parking preserves the content — nothing is discarded or stashed, and the
   branch is recoverable with one `git checkout`. The first encounter still
   stops, so a human editing interactively is never interrupted mid-change;
   only a tree that survived a full slot unattended gets moved. Never delete a
   `recovered/*` branch; the daily review disposes of it.
2. Take the FIRST item under "On deck" in §9 that is not marked done or
   blocked. Do not choose a different item for any reason. **The item's
   brief** (operator directive 2026-10-01) is the `### §9 item N` section of
   `docs/planning/chunks/<CHUNK>.md` that its two-line §9 entry links — read
   it whole before acting; an older item without a link carries its brief in
   §9 itself. Mark the item done / blocked in the **§9 entry**, and append
   the outcome line under the brief's heading. If every item is
   done or blocked, fall back to the chunk named in §9's "obvious next entry"
   sentence, scoped to one run; note in attempts.md that you used the
   fallback. If that sentence names nothing, append an attempts.md entry
   saying so and stop.

   **Then, if the item is done early, take the next one.** After an item's
   outcome commit (step 4) has landed and `git status --porcelain` is empty,
   if the clock reads **before minute 30** of the slot, return to this step
   and take the next open item — the same rule, the same order, and it may
   repeat while the clock allows. Measured 2026-09-11/12: slots were
   stopping at minutes 5–20 of 60 because one run meant one item, while the
   compute-safety rationale ("one executor at a time") is satisfied by
   sequential items just as well. Three guards, all load-bearing: a second
   item **never starts on a dirty tree** (the first item's commit is what
   licenses it — an incomplete first item parks per step 4 and the slot then
   stops, it does not move on — **except a pre-registered negative**,
   operator directive 2026-10-01: if the item's brief carries
   `PREDICTED RED: <assertion> range [a, b]` lines, **every** red assertion
   is one of them and lands inside its range, and every other assertion is
   green, the item is parked and journaled exactly as an incomplete item
   (step 4, outcome `blocked`, "predicted red, inside range" in the entry),
   and once `main` is clean the slot **continues** under the minute-30 rule;
   an unpredicted red, a red outside its range, or any other failure stops
   the slot as before); an item that says "depends on item N
   landing" is taken only if N's commit is on `main`; and each item gets its
   own outcome commit and its own attempts.md entry (step 5), never a merged
   one. Nothing here changes the ban on concurrent executors or the
   per-command tier ceilings.
3. Execute the chunk following `.claude/agents/implementer.md` (read it
   first) and the chunk's §7 entry, which carries the implementation plan.
   Specialist executors exist for three chunk classes — an EX-* example
   chunk goes to `example-runner`, a measurement-only mesh/resolution probe
   to `mesh-probe`, a version-bump record sweep to `record-reconciler`
   (.claude/agents/); everything else stays with the implementer agent.
   `plan-navigator` and `log-pathologist` are available read-only for
   lookups and disputed prior logs. Never spawn `auditor` on your own
   closure — auditing is the review's job.

   Three rules on spawning, all load-bearing:
   - **Only this repo's agents exist.** `implementer-run.sh` denies every
     built-in agent type by name, so `general-purpose`/`Explore`/`Plan` are
     unavailable by construction. Never route chunk execution around the
     specialist definitions — the executor's prompt *is* the accumulated
     trap list, and a generic agent would rediscover them at your expense.
   - **One executor at a time, never concurrent.** `example-runner` and
     `mesh-probe` issue compute; two of them in parallel would double-book
     the 12-core budget the guard hook enforces per-command but cannot
     enforce across simultaneous sessions. One chunk per run makes this
     natural — keep it that way.
   - **A delegated chunk is still yours.** The executor's report is
     evidence; you own the commit, the attempts.md entry, and the §4 claim,
     exactly as if you had run it. If its report and the logs disagree, the
     logs win.
   - **Running an example: the closing window is the script through the
     harness** (11-02 weekly ruling, PROJECT_PLAN §10; `OPS-74`; `EX-170`
     … `EX-172` closed on it). `./scripts/run_examples.sh` is host-side and
     is *not* in `sandbox.excludedCommands`, so a scheduled slot invoking it
     bare gets `permission denied ... docker API` from its inner
     `docker compose` call (2026-08-29, 08-30, 09-01). Wrapping it in
     `docker compose exec` is the other half of the same trap and exits
     127 — never do that. The ruling's three clauses:

     1. This, as one bare, whole command, **is the closing window**:

        ```
        scripts/testing/run_and_log.sh <ID> "./scripts/run_examples.sh -e <group>:<n> -n <n> -t <s>"
        ```

        The harness call is sandbox-exempt, so the script's inner
        `docker compose` reaches the socket.
     2. A setup figure needs its own earlier window (with
        `FEM_EM_SETUP_FIGURES=1`) that closes nothing. The closing window
        is the last.
     3. A guide edited after its closing window admits only digits from
        the closing log, or arithmetic on them, and says so ("by arithmetic
        from the log"). A phrase with no log behind it is removed or marked
        "a reading, not measured".

     Always write the selector with its group (`-e mesh:3`, never `-e 3`) —
     hand-built selectors have cost slots. `--dry-run` touches no socket and
     runs fine sandboxed; use it to see the inner command, not as a required
     step. This rule binds a spawned `example-runner` too; say so in its
     spawn prompt.
   - **Executors run in the foreground, and their harness windows carry
     the 660 000 ms timeout.** Spawn with `run_in_background: false`, say
     in the spawn prompt that the "Never run the harness with
     `run_in_background`" rule below binds the executor's own Bash calls
     (foreground, `timeout` 660000 ms, container-side `timeout -k 30`
     sized to return a footer), and never end your turn while an executor
     is in flight: the headless CLI waits at most 600 s for background
     tasks and then terminates the whole session (`Background tasks still
     running after 600s; terminating`), SIGKILLing the executor and any
     harness window it holds. The 2026-09-01 00:00 slot — the first able
     to delegate — lost `EX-36` leg (mesh) this way: `example-runner`
     returned with its third `./run_examples.sh` window still running
     (killed host-side 267 s in, mid-`mesh:8`, footerless; the container
     process ran on and wrote artifacts nobody witnessed), the slot
     "resumed" it and ended its turn, and the session was terminated at
     minute 19 with no journal, no commit and three untracked logs. A
     launcher backstop (`export CLAUDE_CODE_PRINT_BG_WAIT_CEILING_MS=0`
     in implementer-run.sh, so a slot that backgrounds anyway is bounded
     by the 65-minute kill rather than the 600 s ceiling) is proposed to
     the operator, not landed — reviews cannot write `scripts/automation/`.
     Until it is, this rule is the only protection.
4. Outcome:
   - **Complete** (§4-done: verification executed, quantitative assertion,
     harness log + elapsed time recorded): commit code + tests + logs +
     §7 status flip together on `main`, marking the On-deck item done in the
     same commit. A chunk whose §7 row points to a history file
     (`docs/planning/chunks/<ID>.md`, `OPS-46`) gets its step narrative
     appended to that file; the status flip and the ≤ 2-sentence state line
     stay in the row.
   - **Incomplete** (out of time, or blocked): park ALL code changes on a
     branch `attempt/<CHUNK-ID>-<UTC-timestamp>` (commit there, return to a
     clean `main`). On `main`, commit only the attempts.md entry plus any §7
     annotation (🟡/🚫 with the blocker named). **Never leave `main` red or
     dirty** — a half-applied change on main costs the next run its
     preflight.
5. Always append an entry to `docs/testing/attempts.md` (append-only): UTC
   timestamp, chunk ID, outcome (`complete|incomplete|blocked|anomaly`),
   what was tried, measured numbers, harness log filenames, branch name if
   parked, and a one-line hypothesis for the next attempt. The daily review
   is the reader — write for it.

## Non-negotiables

Everything in `.claude/agents/implementer.md` applies unchanged: all
verification through `run_and_log.sh`, tier ceilings with kill-and-shrink on
overrun, never loosen assertions, rank-safety reductions, no work on ⚠️
subsystems, known-issues.md discipline for unrelated failures.

In tracked Markdown (plan annotations, attempts.md, known-issues, example
guides), write "approximately" as `≈`, never `~`: GFM renders a single
`~text~` as strikethrough, so two approximate tildes in one paragraph cross
out everything between them (swept 2026-09-28). Literal tildes go in
backticks or as `\~`.

## Working inside the permission allowlist

Scheduled sessions run against `.claude/settings.json` with no human to
answer prompts — a denied tool call is simply denied. Consequences:

- All compute goes through `scripts/testing/run_and_log.sh ...` or
  `docker compose exec ...` at the **top level** of the Bash command. Host
  `python3`/`pytest` are not allowed and the host lacks dolfinx anyway.
- **Never run the harness with `run_in_background`, and never end your turn
  while a harness command is running.** You are a headless `claude -p`
  session: ending the turn exits the CLI (exit 0), which SIGKILLs the
  backgrounded harness — footerless log, untracked artifacts, no journal
  entry, and the tree left dirty for the next slot. Three consecutive slots
  paid this on 2026-08-10/11 (19:30, 22:30, 00:00 — the whole night) on the
  same item, and a spawned `example-runner` paid it again on 2026-09-01
  00:00 (step 3, fourth rule). Run harness commands in the **foreground** with the Bash tool's
  `timeout` parameter at its 660000 ms maximum, and size the
  **container-side** `timeout` so the command returns a footer inside that
  window (≤ ≈590 s of container time for a command with ≈1 min of setup).
  Always write the container-side timeout as **`timeout -k 30 <s>`** — a
  plain TERM does not reliably stop an `mpiexec` job (MAT-6 step 10,
  2026-08-12: `timeout 590` fired and the ranks ran on for ≈1 700 s,
  wedging the container). If a container wedges anyway (`exec` hangs,
  `restart`/`kill` report "did not receive an exit event"), the recovery
  is `docker compose -f docker/docker-compose.yml up -d --force-recreate`,
  then verify `memory.max` and zero stray `python3` before continuing. A
  compute step that cannot finish inside one foreground window at ≤ 12 ranks
  is too big for a scheduled slot: shrink the case or journal it for the
  review — do not background it.
- Read files with the Read/Grep/Glob tools, not `cat`/`sed`/`head` — the
  generic readers are not allowlisted, deliberately.
- Multi-line commit messages: use a **literal multi-line `-m` string**
  (a quoted argument containing newlines). The older `git commit -F <file>`
  route does not work in scheduled sessions — the 2026-09-01 19:30 slot found
  `Write` denied for both `$TMPDIR` and `.git/`, so there is no writable
  target for the message file (corrected by the 2026-09-02 03:00 review;
  `CLAUDE.md`'s `-F` advice is for interactive sessions, where Write is
  available). Command substitution like `git commit -m "$(cat ...)"` is
  treated as injection by the permission layer and will be denied.
- If a genuinely needed command is denied, do not fight it: document the
  denial in your attempts.md entry so the daily review can propose an
  allowlist change to the human.

- AED benchmark numbers are private (CLAUDE.md hard rules): read the gitignored `aed_results/` / `COMPARISON_private.md` / `docs/private/` files if you need context, but never copy a number from them into a tracked file, a journal entry or a commit message.


## Standing execution rules (from §9, moved verbatim 2026-10-24)

*Operator ruling 2026-10-24 (interactive session): these blocks used to sit in PROJECT_PLAN §9 above the queue. They bind every slot, every executor a slot spawns, and every review that writes a brief. They keep their letters — "rule (d)", "rule (i)", "rule (j)" cite this section.*

**⚠️ Standing constraint on the compose allow — read before editing that
file.** `docker-compose.yml` line 9 is `- ..:/workspace`, so write access
to it is write access to *what the container mounts from the host*. The
operator granted this knowingly and narrowly (2026-08-22). **Edit only
`environment:` keys. Do not touch `volumes:`, do not add a mount, do not
widen a path, and do not change the memory limit (**128 G**, raised from
64 G by operator directive 2026-08-24) — in this or any future chunk.**
*(The `fem-em-solver-xl` service added 2026-09-05 is the one operator-
authorized exception — a second service with its own 512 G limit under the
`xl` profile; the constraint applies to it verbatim, and chunks do not
edit it either. The `fem-em-palace` service added 2026-09-28 — `OPS-61`'s
Palace image, `palace` profile, the same `..:/workspace` mount and 128 G —
is the second operator-authorized exception, on the same terms.)* A
chunk that believes it needs a mount change is a **blocked finding for the
operator**. The `Edit(docker/.claude/**)` caution stands for the same
reason: a nested `.claude/` is a settings-override surface. One surviving
mechanic: `git checkout` cannot swap `docker/Dockerfile` /
`docker-compose.yml` in this sandbox — bind-mounted, "Device or resource
busy", a *silent* wrong-content switch — so any chunk that must move them
uses the Edit tool and verifies `git status --porcelain`.

**The docker-socket "runner trap" is diagnosed (`fcffa701`, operator's
interactive session, 2026-09-30):** every `permission denied …
/var/run/docker.sock` since 08-29 — slots and spawned executors alike — was
a **non-bare** call: a `cd …` / `date -u;` prefix, a `; date` / `| tail`
suffix, a pipe or a redirect. Only whole commands matching `docker *` /
`scripts/testing/run_and_log.sh *` (repo-relative; an absolute path is
denied) are sandbox-exempt, and `bash_guard.py` rule 3 now denies the other
shapes before they run. So: every harness call is the whole Bash command;
`./run_examples.sh … --dry-run` only prints the per-example `docker compose
exec` lines — run each through `run_and_log.sh` as its own call, never
chained with `&&`. *(Old paragraph archived: `plan-archive.md` § "§9
runner-trap paragraph".)*

**Standing rules for every item below (from the 10:30 review, kept):**
(a) an `EX-*` item may add an *additive* return key or keyword to a gate
fixture's helper it imports, provided the gate module is re-run green
through the harness in the same slot and the diff is disclosed; (b) guide
artifact references carry the full filename; (c) a `src/` change beyond an
item's letter follows ruling (a) above — disclosed, its pre-existing gate
re-run green in the same slot; **(d) (added 2026-09-05 10:30) a slot that
parks its item on `attempt/*`, or otherwise leaves it neither done nor
runnable, marks the §9 item BLOCKED with the unblock condition in the
*same* commit as the record** — the next slot reads §9 first, and an
unmarked item costs it the ten minutes the 09:00 slot spent on
`e49cb67`; **(e) (added 2026-09-05 18:00) every negative-control factor
in an item is labelled *asserted* or *predicted*: asserted only when a
prior measurement of the *same* comparison on the *same* fixture backs
it (cite the log line), predicted otherwise — and a predicted factor is
printed beside the measured one, never asserted. An executor that meets
a failing pre-registered assertion whose label is missing stops and
reports the negative result; it does not decide in-slot** (`EX-49`,
2026-09-05 15:00, ruled by the 18:00 review; `MAT-6` step 11, 2026-09-06
07:30, ruled by the 10:30 review); **(f) (added 2026-09-06 18:00) a band
an item states against an *analytic* comparand that the executed
fixture's own CAD or mesher cannot meet is re-registered by measurement
on the same quantity where both sides are the fixture's numbers, said
so in the docstring and the journal, and the review ratifies or reverts
it** (`TH-15` step 3a, 16:30 slot — ratified above, ruling (2)); it is
never widened silently and never on a quantity that was already green.
**(g) (added 2026-09-08 10:30)** every pytest window below runs with `-s`
— two windows (≈ 9 min) were lost this interval to swallowed prints
(`TH-15` step 2f, `OPS-41`), each re-run only because the executor
noticed; a log without the readings is a window not spent. **(h) (added
2026-09-13 18:00)** an executor that re-registers a red anchor's
*comparand* in-slot (`TH-15` step 3, 15:00) does so as rule (f) allows only
when the new comparand measures the *same physics* as the old one and
says so; a substitute that is true by construction on the fixture (an
identity with one side removed by the mesh) is a **record**, and the item
stays at its pre-registered result — red — for the review to rule on.
**(i) (added 2026-09-13 18:00)** the window a step is closed on runs the
module **as committed**: if the module is edited after its last green
window (even a comment), the window is re-run before the commit, or the
row says it was not. **(j) (added 2026-09-14 03:00)** a review that writes
"**Status it can move:** `X` → ✅" first reads `X`'s §7 Done-when and either
matches it or amends the row in the same review commit; an executor that
meets a mismatch closes on the *row* (holds 🟡) and flags it, as the 21:00
slot did with `TH-14` — the §9 letter never outranks the row.


**Every window below that uses §5.1 durable capture copies the idiom
verbatim, trailing `; exit $rc` included.** Since `OPS-45` the footer honours
a `[capture] rc=` line only when it is the *last* output line. **Every pytest
window runs with `-s`** (rule (g)). **Every window runs inside the container
through `run_and_log.sh`** — host `python3` is denied.

**Every item:**
- The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. **A slot that spawns an executor puts that sentence in the spawn prompt verbatim.**
- **An executor that loses the docker socket mid-item stops and reports.**
- **No edit to a script, test or page after its last window** (rule (i)).
- Every Bash call that runs a window passes `timeout: 660000`. Run pytest with `-s`. Never delete a log.
- An item that parks marks itself BLOCKED in the same commit (rule (d)).
- **A window killed at its tier ceiling is not a result.** After a kill, run `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` first.
- **Two heavy windows never run at once.**
- **Read the chunk file's latest "Review ruling(s)" paragraph (the newest dated one; no item is open after the 2026-10-17 18:00 review) before the brief.**
- **Any item that adds lines to `WF-14`'s registry regenerates the API pages and `docs/guide/validation.md` and runs the page gates in the same slot** (`WF-1.md`, 03:00 ruling Q3; item 394's Wgen / Wval / Wapi / Wdoc windows).
- No slot opens `docs/private/` or `aed_results/`.
- **A smoke-tier (30 s) window puts only `tests/environment/test_complex_mode.py` on the path, never the full `tests/environment`** (10:30 review ruling, item 364).
