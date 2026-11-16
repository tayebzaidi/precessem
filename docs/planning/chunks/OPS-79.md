# OPS-79 — the two tracked closing-window texts amended to the 11-02 ruling

*History file for the §7 `OPS-79` row (`OPS-46`). Opened by the 2026-11-03 03:00 daily review from the 11-02 weekly's closing-window ruling (PROJECT_PLAN §10, "Rulings, 2026-11-02", the closing-window bullet), which names this chunk. Lane 1, docs.*

### §9 item 792 — 2026-11-03 (lane 1, ops / docs; no compute beyond a dry-run)

**What.** Two tracked texts still teach emit-then-harness:
- `scripts/run_examples.sh:39–48`: "Until that exclusion is added, the supported path for a scheduled slot is emit-then-harness …"
- `docs/automation/implementer-run.md:133–151`: "Running an example: emit, then harness — do not improvise."

Replace each with the ruling's three clauses, citing it ("11-02 weekly ruling, PROJECT_PLAN §10; `OPS-74`; `EX-170` … `EX-172` closed on it"):
1. `scripts/testing/run_and_log.sh <ID> "./scripts/run_examples.sh -e <group>:<n> -n <n> -t <s>"`, as one bare, whole command, **is the closing window.** The harness call is sandbox-exempt, so the script's inner `docker compose` reaches the socket.
2. A setup figure needs its own earlier window (with `FEM_EM_SETUP_FIGURES=1`) that closes nothing. The closing window is the last.
3. A guide edited after its closing window admits only digits from the closing log, or arithmetic on them, and says so ("by arithmetic from the log"). A phrase with no log behind it is removed or marked "a reading, not measured".

**Keep** the two traps that are still true: never wrap the script in `docker compose exec` (Status 127), and always write `-e group:N`, never `-e N`. `--dry-run` stays documented as a way to see the inner command, not as the required path. **Do not touch** `.claude/agents/example-runner.md`. Its `:57–72` carry the same text, but that file waits for the operator's worktree session (W4); say so in the commit.

**Anchor: an identity, not a band.**
- **(a) The script's diff is comment-only:** every changed line in `git diff -U0 scripts/run_examples.sh` begins with `#`.
- **(b) Its behaviour is unchanged:** `./scripts/run_examples.sh -e mat:2 --dry-run` prints byte-identical output from `HEAD`'s copy (`git show HEAD:scripts/run_examples.sh` into a repo-root scratch file) and from the edited file.
- **(c)** `bash -n scripts/run_examples.sh` exits 0.

**Negative control.** Run the same comparison as (b) on a scratch copy with one character changed in a non-comment line. The outputs must differ.

**Windows.** Smoke, no container, `timeout -k 30 60`. Use one bare harness window per check, e.g. `scripts/testing/run_and_log.sh OPS-79-W1 "timeout -k 30 60 bash -n scripts/run_examples.sh"`. A check that needs a pipe or a comparison goes into a repo-root scratch script run by `bash <scratch>.sh` inside the harness. Remove the scratch files before the commit.

**Class:** ops / docs. **Cost:** ≈ 8.8 measured / 20 predicted min.

**Traps already paid for.** The bash guard denies harness calls that are not bare (no `cd … &&`, no pipe outside the quoted inner command).

**Scope.** Two texts. **Not claimed:** the agent file.

**Negative result.** If (b) differs, a non-comment line moved: revert it and re-run.

**Status it can move.** `OPS-79` ⬜ → ✅. **Depends on:** nothing.

**Outcome (22:30 slot, 2026-11-03): DONE, `OPS-79` ✅.** Both texts replaced with the ruling's three clauses, each citing "11-02 weekly ruling, PROJECT_PLAN §10; `OPS-74`; `EX-170` … `EX-172` closed on it"; the Status-127 and `-e group:N` traps kept; `--dry-run` documented as a way to see the inner command only. (a) `20261103T223106Z_OPS-79-W2.log`: 30 changed lines, 0 not beginning with `#`. (b) `20261103T223107Z_OPS-79-W3.log`: `-e mat:2 --dry-run` from HEAD's copy and the edited file both exit 0, sha256 `8f6c149583654b83770ff4e829ae28ac7d4f044de307178dcb0cd0e8822790d5` on both, `cmp` identical. (c) `20261103T223105Z_OPS-79-W1.log`: `bash -n` Status 0. Control `20261103T223108Z_OPS-79-W4.log`: `NPROC=2` → `NPROC=3` in a scratch copy changes 2 output lines (`mpiexec -n 3`), PASS. Each window 1 s. One deviation from the brief's letter: the HEAD and control copies were written into `scripts/` (not the repo root) because the script derives `ROOT_DIR` from `dirname "$0"/..` — a repo-root copy would resolve a different root and differ for that reason alone; the comparison driver was a repo-root scratch script. All scratch files removed. `.claude/agents/example-runner.md` untouched (its `:57–72` wait for the operator's worktree session, W4).
