# OPS-54 — XL / XXL windows run against an immutable checkout

### §9 item 487 — 2026-10-15 (lane 1)

**Step 1: the scratch-repo clause of the Done-when, behind an opt-in key, with the default path byte-identical.** Implementer; smoke; no MPI, no container, no compute.

**Why it is queueable now.**
- The §7 row says "do not queue until a measured window has actually overlapped a slot, or the weekly decides the risk warrants it". The weekly has decided twice:
  - "`OPS-54` — queue it; it is the right use of a drained slot" (`plan-archive.md:34768`);
  - "`OPS-54` step 1 if assertable in a slot", in the current lane-1 list (`PROJECT_PLAN.md` §10, "Lane 1 and lane 2, targets for the daily").
- The earlier "the operator's (a compose mount)" (`plan-archive.md:35131`) is answered by keeping the mount change a default-preserving interpolation. The review rules on that reading.
- **Rule (j):** the Done-when has two clauses: (i) the scratch-repo test, (ii) one real XL window through the new path. This step meets (i) only.

**Construction (additive).**
- `xl-run.sh` reads a new queue-entry key, `XL_IMMUTABLE_CHECKOUT=1`, beside `XL_OUTPUTS`. Unset leaves today's path untouched.
- When the key is set, `xl-run.sh`:
  - records the window's start `HEAD` sha;
  - runs `git worktree add --detach <gitignored path, e.g. .xl-worktree/> <sha>`;
  - exports `FEM_EM_XL_WORKSPACE` to that path before its existing `up -d --force-recreate`;
  - runs the queued command against that view;
  - after the window, copies the declared `XL_OUTPUTS` and any in-tree raw-capture files back into the main tree before its existing staging step (`xl-run.sh:199–245`);
  - removes the worktree (`git worktree remove --force`, then `git worktree prune`).
- In `docker/docker-compose.yml` the xl and xxl services' `..:/workspace` becomes `${FEM_EM_XL_WORKSPACE:-..}:/workspace`. Edit it with the Edit tool (see traps).
- The gitignore entry goes in `.gitignore`. The scratch test copies that file (`test_xl_window_isolation.sh:53`).

1. **Anchor (byte identity).** In a new window 7 of `scripts/testing/test_xl_window_isolation.sh`, the queued probe:
   - reads `src.txt` in its own view (`$PWD` under the stub, i.e. the worktree);
   - commits a change to `src.txt` on `main` in `$REPO` mid-window (a real `git commit`, not an edit);
   - re-reads its view.

   Asserted:
   - the view's sha256 after the concurrent commit equals `git show <start-sha>:src.txt`'s sha256;
   - `main`'s `src.txt` does carry the change;
   - after the window, the worktree path is gone and `git worktree list` has one entry;
   - a declared output written inside the view is present in the xl-run commit.
2. **Negative control.** Window 7′: the same probe and the same concurrent commit with the key **unset**. Assert that the view's sha256 **differs** from the start blob, i.e. the live tree moved under the window, which is the defect. This is a planted boolean control, asserted as OPS-64's windows 5–6 are. Rule (e) governs factors, and no factor is claimed. **Regression:** windows 1–6 stay green unmodified. This is the evidence that the default path is unchanged.
3. **Tier, ranks, wall clock.** Smoke, no ranks. Priced from `20260920T004915Z_OPS-XL-ISOLATION.log` (23 s) and `OPS-64`'s three windows of 24 / 24 / 47 s (`docs/testing/test-results.md:2430–2432`). Predicted ≈ 45 slot-minutes, mostly script work. The window is one bare whole command: `scripts/testing/run_and_log.sh OPS-54-step1 "timeout -k 30 120 bash scripts/testing/test_xl_window_isolation.sh"`.
4. **Traps.**
   - `bash_guard.py` denies any Bash text that names the XL service (weekly-review mechanics). Edit compose and read it with Read/Edit/Grep tools, never with a shell command naming the service. The test file's existing strings are fine, because the harness command does not contain them.
   - A git worktree's `.git` is a *file* holding an absolute host path. In-container `git` cannot resolve it. The harness runs host-side, so step 1 is unaffected. Name this as a step-2 check, not a step-1 claim.
   - `mktemp` must honour `$TMPDIR` (`:37` already does). `/tmp/claude` is denied in headless sessions.
   - The real `xl-run.sh` is live at 02:00 tonight. Never run it against the real repo or docker. Never write to `docs/testing/xl-queue.d/` (implementers never commission). The key must default off, so tonight's window is byte-identical in behaviour.
   - Durable-capture raw files and `--capture-orphan` read repo-relative paths. Under the key they land in the worktree, so copy them back before the worktree is removed, or the window loses its record.
   - Never pipe the test's output inside the harness command.
5. **Scope boundary.**
   - Clause (ii) is not attempted: no real XL window, and no flip of the default to the key. Both belong to the weekly, which commissions the step-2 window with `XL_IMMUTABLE_CHECKOUT=1` in its queue entry.
   - No change to `bash_guard.py`'s budget arithmetic, to the XXL cadence, or to any compute.
6. **Negative result.**
   - Window 7 red (the view moved, or the worktree left behind), or any of windows 1–6 red: report and stop. Park the work on `attempt/OPS-54-step1-…`, mark the item BLOCKED in the same commit, and leave `main`'s `xl-run.sh` and compose untouched.
   - If 7′ does **not** differ, the probe is not reading the window's view. The control is broken and nothing is claimed.
7. **Status it can move.** `OPS-54` ⬜ → 🟡 on clause (i). ✅ waits on clause (ii), the weekly's real window.

### Review ruling 2026-10-15 18:00 (item 487 — NOT queued)

**Declined; the draft above is banked for the operator.** It edits `docker/docker-compose.yml`'s `volumes:` (`..:/workspace` → `${FEM_EM_XL_WORKSPACE:-..}:/workspace`). §9's standing constraint on the compose allow says: "Edit only `environment:` keys. Do not touch `volumes:`, do not add a mount, do not widen a path". It also says "A chunk that believes it needs a mount change is a **blocked finding for the operator**". A default-preserving interpolation still makes the mounted host path settable from the environment, which is a widening. Only the operator can authorize it. The dashboard carries it as the operator's (carried item "`OPS-54` is yours (a compose mount)"). Number 487 is spent and is not reused.
