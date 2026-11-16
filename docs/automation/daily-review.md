# Daily review protocol (scheduled)

*Model: `scripts/automation/review-model.env` (dated override, self-expiring —
see CLAUDE.md § Scheduled automation). High effort (operator directive
2026-09-15; was medium).*

Run by `scripts/automation/daily-review.sh` via cron **three times a day,
every day** (03:00 / 10:30 / 18:00 local; full schedule restored 2026-09-28,
operator directive — the 2026-09-13 wind-down is over; canonical crontab
`scripts/automation/crontab`), each followed by four implementer runs (04:30 /
06:00 / 07:30 / 09:00, 12:00 / 13:30 / 15:00 / 16:30, 19:30 / 21:00 / 22:30 /
00:00; no 21:00 run on Monday, Wednesday and Saturday, when the weekly holds that
window). "Since the last review" therefore spans about eight hours, and the
restock floor is one four-slot interval. One
session, documentation work only — **no solves, no meshing**; reading harness
logs is fine. You are maintaining the plan, not executing it.

**Subagents are available to you, and web tools are not.** Because this session
never solves or meshes, a subagent costs tokens rather than cores and does not
touch the 12-core compute budget. Use them for the mechanical, parallel, read-
heavy work — one `auditor` (.claude/agents/auditor.md) per newly-✅ chunk in
step 3, a `log-pathologist` ruling on any disputed or status-changing log
reading in step 4, `plan-navigator` for plan/known-issues lookups instead of
paging those files yourself, an `Explore` sweep when
step 5 asks whether the backlog still reaches §10 — and keep your own context
for the judgement calls: what a negative result means, and what to queue next.
A subagent's report is evidence, not a verdict; if one says a chunk passes §4,
the log it cites is what you cite back in the commit. Do not delegate steps 2,
6, 7, or 8 — disposition, queue order, the dashboard refresh, and the commit
are yours.

## Steps

1. Establish what happened since the last review:
   - `git log` since the previous `docs(plan): daily review` commit (or 24 h
     if none exists)
   - new rows in `docs/testing/test-results.md` and logs in
     `docs/testing/logs/`
   - new entries in `docs/testing/attempts.md`
   - `git branch --list 'attempt/*'` for parked incomplete work
   - `git status --porcelain -uno` — dirty tracked files at review time mean
     every implementer run since they appeared has been tripping preflight

2. **Clear any stalled tree, and dispose of `recovered/*` branches.** A dirty
   tracked tree older than one implementer cycle (90 min) is an outage, not a
   curiosity, and this review is the scheduled actor responsible for ending
   it. Read the diff and the attempts.md anomaly entries about it, then
   resolve it now: commit the changes (accurate message, its own commit) if
   they describe reality, or revert them if they do not, and record which
   you did and why in the review commit. Documentation-only diffs that a
   prior run journaled as an anomaly should normally have been landed by the
   next implementer run (implementer-run.md step 1 exception); if one is
   still sitting here, also note why that didn't happen. Never leave the
   tree dirty at the end of the review.

   `git branch --list 'recovered/*'` lists trees a *second* implementer
   encounter parked rather than stopped for (implementer-run.md step 1). Each
   one is a change nobody has adjudicated: read it, then land it, fold it into
   a §7 entry, or delete the branch — and say which in the review commit. An
   accumulating `recovered/*` list means the tree keeps going dirty from a
   source nobody has found; name that in the commit rather than clearing it
   silently.

3. **Audit every chunk whose status changed to ✅ since the last review**
   against PROJECT_PLAN.md §4: does a harness log exist, was the verification
   executed by the agent itself, is at least one assertion quantitative
   (closed form / convergence rate / conservation, reciprocity, or symmetry
   identity), is elapsed time recorded? Demote anything non-compliant to 🧪
   with a dated note. Do not re-run anything. A chunk's evidence is its §7
   row **plus** its history file `docs/planning/chunks/<ID>.md` when the row
   points to one (`OPS-46`); the auditor reads both. Delegate each audit to the
   `auditor` agent; treat its DEMOTE/PASS as evidence and re-cite its
   log:line evidence yourself in the review commit. Through 2026-09-03:
   re-verify one cited claim per agent report (a random digit trace) before
   acting on it; this clause expires on its own after that date.

4. For each incomplete attempt (attempts.md entries + `attempt/*` branches):
   diagnose from the logs and the parked diff; rescope the chunk's §7 entry —
   smaller case, sharper implementation plan, or split into two chunks — and
   record the diagnosis in the entry. Delete an attempt branch only when its
   useful content is fully captured in the plan. Before banking a diagnosis
   whose log reading is surprising or changes a status, get a
   `log-pathologist` ruling; UNCOUNTABLE means diagnose from the parked diff
   only.

5. Assess against §10 success criteria: does the existing backlog still lead
   to the mission? If a gap exists, add new chunk entries (stable IDs,
   §4-compliant done-whens, implementation plans meeting the rubric below).
   If no gap exists, do not invent work. **Scope boundary:** the §6 phase map
   and §10's long-horizon roadmap (phases, subgoals, dated assessments)
   belong to the weekly planning review (docs/automation/weekly-review.md) —
   add chunks *within* the current phase's subgoals; do not restructure
   phases or edit the roadmap here.

   **Example chunks (§5.4 ramp):** for each chunk that newly closed a
   quantitative gate since the last review (the step-3 list, post-audit),
   check whether an existing example already demonstrates that capability;
   if not, add a standalone example chunk to §7 — sized for one implementer
   run, executed via `./run_examples.sh`, producing combined-XDMF that opens
   in ParaView, and demonstrating the capability from an angle no existing
   example covers (geometry, materials, drive, or output quantity). Example
   chunks are never riders on physics chunks and never target ungated
   capability.

6. Refresh **"On deck"** in §9. **The restock floor counts *measured*
   slot-minutes (operator directive 2026-10-01; was predicted minutes from
   2026-09-13):** the queue must carry **≥ 160 measured slot-minutes** of
   open, unblocked work — four slots × ≈ 40 min, the work a slot can
   actually consume, since it starts items only before minute 30 — **and
   ≥ 5 items** (the spare). Predicted minutes ran 5–10× high (09-30: six
   items predicted at 245 min landed in 33; 17 of 36 slots short of work),
   so they no longer count toward the floor; each item still states its
   predicted cost (rubric element 3), which sets its tier and timebox.
   **An item's measured slot-minutes** are its class's median over the last
   ≥ 5 landed items of that class, measured from the slots' own commits: an
   item's duration is its outcome-commit time minus the slot's start (first
   item) or minus the previous item's outcome commit (take-next). Classes:
   chain step (solve), chain step (no solve), lane-1 revalidation, lane-2
   feature step, lane-3 closed form, example (`example-runner`), ops / docs.
   A class with < 5 landed items uses its predicted minutes × the last
   interval's overall measured / predicted ratio until it calibrates. Write
   the class table (class, n, median minutes) in the review journal every
   run and the running measured total beside the list. Since 2026-09-12 a slot that commits its item with
   a clean tree before minute 30 takes the next one (implementer-run.md
   step 2), which is why five items is not a floor any more: on 09-12 a
   five-item queue of cheap serial moves lasted two slots and three further
   slots drained with every review alive. Each item is still sized for one
   run (≤ 1 h wall clock, ≤ 20 min per compute command, ≤ 12 ranks) and
   ordered; independent items first. If the ready items do not reach the
   floor, list what exists, **state the shortfall in minutes**, and say so —
   step 5 still forbids inventing work, and a visible shortfall is the
   signal the weekly needs. Each listed item meets the rubric below; an item
   that cannot yet state its anchor is not ready to queue, and writing that
   anchor is itself the better queue item.

   **Where an item lives (operator directive 2026-10-01).** The full brief —
   every rubric element, the spawn prompt, commands, `PREDICTED RED` lines —
   goes in the chunk's history file `docs/planning/chunks/<CHUNK>.md` under
   a heading `### §9 item N — <date> (<chain step / lane>)` (create the file
   if the chunk has none; the §7 row already points there or gains the
   pointer). **§9 carries two lines per item:** line 1 —
   `N. <status> — <chunk> <step> (<chain step / lane>; <executor>; <tier>,
   -n <ranks>; ≈ <measured> measured / <predicted> predicted min;
   depends on <item | nothing>)`; line 2 — the anchor in one clause and the
   link `→ docs/planning/chunks/<CHUNK>.md#§9-item-N`. When an item lands,
   its §9 entry collapses to the single DONE / BLOCKED line the digest
   already uses; the brief stays in the chunk file as history. §9 itself
   keeps only the rules, the rulings that bind the next slots, the entries
   and the running totals — **target ≤ 25 KB**. *Migration:* the first
   review after this directive moves every open item's body from §9 into its
   chunk file verbatim (the zero-loss check: each moved body byte-identical
   at its destination) and leaves the two-line entries.

   **§7 rows follow the same split (operator directive 2026-10-01): each
   row ≤ 2 KB** — ID, short title, glyph, a ≤ 2-sentence current-state line,
   the Done-when, the tier. When you amend a row (a ruling, a re-scope, an
   audit result), write the narrative into `docs/planning/chunks/<ID>.md`
   and change only the glyph / state line / Done-when in the row. The
   weekly runs the census and moves any row that has grown past 2 KB
   (weekly-review.md step 6).

   **The enumerated-chain licence (operator directive 2026-09-18; §9 rule
   (4)).** "Forbids inventing work" does not forbid *executing the weekly's
   plan*. When the ready items do not reach the floor, you may write and
   queue — at the full rubric grain, anchor and control and status-it-can-
   move included — any **numbered** step the weekly has enumerated in a §10
   chain (today the "step count to a tuned birdcage" list, §10 Phase 6)
   whose serial prerequisites have landed on `main`. The weekly supplied the
   question and the validation target; you supply the execution detail
   (fixture, files, tier by measurement or cost probe, imported bands). You
   may not add a step to a chain, reorder one, or change a step's target; a
   step whose target you cannot make assertable is queued as "write the
   anchor" (a documentation item) and named in the dashboard for the
   weekly. Say in the item which chain step it executes. This is what
   keeps four slots a day from draining into figures while a written plan
   sits idle between weeklies.

   **The backlog-lane licence (operator directive 2026-09-28; §9 rule
   (5)).** When the ready chain steps do not reach the floor, also write
   and queue, at the same full rubric grain, items from the three standing
   lanes §9 rule (5) defines — **lane 1** revalidation of `🧪` / `⚠️`
   chunks, **lane 2** feature-ladder Tier B (`TH-16`, then `TH-18`, then
   `TH-5`), **lane 3** Phases 7 and 8 and their materials (`MAT-3`,
   `MAT-5`, Pennes gel regime, implant geometry). Fill the shortfall
   round-robin across the lanes; each item is independent of the chain
   critical path and never displaces a ready chain step above it. Here
   "forbids inventing work" means you do not invent *targets*: a lane item's
   anchor is the chunk's own done-when, a closed form (Jin or another cited
   text, chapter/equation named), or a cited published result. A lane-3
   chunk you open because the weekly has not yet written one is named in
   the dashboard for the weekly to ratify. **A drained queue is a supply
   problem, never a reason to recommend fewer slots:** if chains and all
   three lanes are genuinely empty, say which lane is empty and why in a
   🔴 dashboard item.

   **The chain-extension licence (operator directive 2026-10-02; §9 rule
   (6)).** When the ready chain steps and the lanes still do not reach the
   floor — or a chain's last enumerated step has landed or been parked — you
   may **write the next numbered steps of a live §10 chain** yourself, at
   full rubric grain, toward that chain's **stated subgoal and validation
   target** (the weekly's; you never write or change one). Tag each step in
   §10 *"(daily, YYYY-MM-DD — provisional)"* and queue it like any chain
   step. **At most six provisional steps per chain** may stand between
   weeklies (raised from three, operator ruling 2026-10-24 — W8); the weekly ratifies or reverts them. You may also **re-pose a
   blocked step** — a different fixture, split windows, a measurement that
   separates two readings — when the re-pose keeps the step's target and
   its family's attempt count (the four-attempt cap still applies; a
   re-pose is the next attempt). Never under this licence: set or move a
   band, tolerance, done-when or validation target; open a new chain,
   subgoal or phase; ratify your own provisional ruling. Bands you rule
   provisionally stay provisional until the weekly ratifies (as `WF-12`'s
   and `MAT-15`'s did, 2026-10-01). Do the supply work in the review that
   finds the shortfall — usually 18:00 and 03:00 — so the 10:30 run stays
   light when the morning queue holds.

   **The setup-figure filler (operator directive 2026-09-13, `EX-57`).**
   After the real items, **always append one `EX-57` item** — the recurring
   background task the operator asked for — naming the scripts that
   `python3 scripts/testing/check_example_setup_figures.py --next 3` prints,
   cut at the first one outside the leading script's group directory
   (run it; do not guess the order). **Batched since 2026-09-25** (operator
   directive): one `example-runner` run draws all of them, so its start-up
   reads are paid once, not once per figure. Cost it at the examples'
   recorded windows + 15 min + 10 min per figure after the first (≈ 20
   slot-min for one, ≈ 40 for three), tier the heaviest example's,
   `example-runner`, independent, "**Status it can move:** `EX-57` census
   `missing` − k" (k = scripts named). Its
   done-when is the §7 `EX-57` entry's per-item list; copy nothing else.
   When the census prints nothing (`missing=0`), flip `EX-57` to ✅ in §7
   and stop appending. This item counts toward the slot-minutes floor but
   never displaces a physics item above it, and it is also the
   drained-queue fallback in §9, so a slot that empties the queue draws the
   next figure instead of stopping.

   **The family cap (operator directive 2026-09-12, four attempts).** A
   *step family* is a numbered step plus all of its lettered sub-steps
   (`4`, `4a`–`4k` are one family). Count the attempts since the last
   sub-step in the family that landed a gate — an asserted, pre-registered
   anchor that moved a §7 status, retired a known-issues entry or registered
   a record; a green run of *imported* gates does not count. **At four, the
   family is frozen**: this review must either re-scope the question at a
   higher altitude as a *new numbered step* or bank the measured negative and
   close the question in the §7 row and known-issues — it may not queue a
   fifth letter. Relettering does not reset the count; a landed gate does.
   The rule replaced §9's "fails twice" rule, which never fired because every
   reletter was formally a new item (`WF-6` step 4 ran twelve).

   **Prefer independent items over a dependency chain.** Four runs will take
   items 1–4 in order without waiting for each other's results, so an item
   that only works if the item above it landed will fail for reasons that are
   not about that item. Where the critical path is genuinely serial, say so
   explicitly in the item ("depends on item 1 landing; if it did not, skip to
   item 3") rather than leaving the run to discover it.

6b. **XL clerk and XL backlog floor (operator directives 2026-09-13 and
   2026-09-15).** The XL / XXL windows are the cheapest compute the project
   has — cron scripts, no tokens — and they are meant to run **every
   night**. This step keeps them fed. Read
   `docs/testing/xl-pending.md` and both ledgers, then do all of the
   following, in order.

   1. **Mark runs.** When a ledger has gained a row since the last review,
      mark the matching pending entry `RUN <log>` and remove nothing else.
   2. **Queue every READY entry, not just the top one.** For each tier,
      copy each `READY` entry whose prerequisite has landed on `main` into
      its own file in `docs/testing/<tier>-queue.d/` — named
      `<intended run date>-<chunk>.env`, `XL_CHUNK` / `XL_COMMAND` (and
      `XL_OUTPUTS` when the entry declares one, `OPS-64`) copied
      **verbatim** (the Write tool, not a shell redirect — the guard trips
      on XL command text in a shell command) — up to the tier's budget:
      `xl` six charged rows per trailing 7 days of `xl-ledger.md`, one
      window per night Sun–Fri; `xxl` one per trailing 7 days, Saturday.
      Run dates are consecutive tier nights from the next 02:00. Mark each
      entry `QUEUED <today> for <run date>`. The launcher runs the files
      in lexical order, re-checks the budget itself and deletes each file
      after its window, so queueing ahead over several nights is
      exactly the point. Never alter a command, never queue on a dirty
      tree. The legacy one-slot `<tier>-queue.env` stays empty.
   3. **Check the floor.** After queueing, count the `xl` entries that are
      `READY` or `QUEUED` for a future night, and the `xxl` ones. **Floor:
      ≥ 4 `xl` and ≥ 1 `xxl` ahead** — enough to cover every window until
      the next weekly. Write the count and any shortfall in the review
      commit and in the dashboard's Automation-health section as *XL:
      N of 6 windows used in the trailing 7 days, M entries ahead (floor
      4)*, and when the floor is not met add a 🔴 Waiting-on-you item
      that says so.
   4. **Fill the shortfall under the daily licence.** When the floor is not
      met you may **write new pending entries yourself**, in two classes
      only, and must still write them in full (tier, chunk and step, exact
      command, price, readout, decision rule — the same rubric as the
      weekly's):
      - **Priced-family variants:** the same script and mesh as an entry
        or ledger row that already carries a *measured* price at this tier,
        varying only frequency, polynomial degree, drive or port set, rank
        count, or one named env knob. The measured price of the family
        member is the price you cite, with the scaling stated (a degree or
        rung change is a variant only when a measured point at the target
        order or rung exists; otherwise it is a cost probe).
      - **Cost probes:** the first window on an unpriced case, whose
        readout *is* the price (§5.1: "the first XL window on it is a cost
        probe whose readout says so"). Fixed brackets in the readout, no
        gate, no band.
      Prefer variants that answer a §10 question the weekly has already
      posed (the ANS-4 degree-2 rung at the other Larmor frequency; the
      F-human degree-1 solve at the full port set; resolution-ladder rungs
      the heavy tier cannot hold; the `TH-16` symmetry-plane cost). New
      physics questions, new fixtures and anything that moves a band stay
      the weekly's to commission; say in the entry which licence class it
      is under. An entry you write is `READY` or `PENDING PREREQUISITE`
      like any other and is queued by rule 2 in the same review.
   5. **Prerequisites.** A `PENDING PREREQUISITE` entry's prerequisite is an
      ordinary §9 item — queue it in step 6 like any other, and write into
      that item that the implementer landing it marks the pending entry
      `READY` in the same commit (so a landing does not wait for the next
      review).

6c. **Keep §9 current-state only** (added 2026-09-19, operator directive —
   every agent reads this section before it starts, and on that date it was
   71 KB). §9 answers one question: *what does the next slot do?* Three
   rules, applied every review, in the same commit as the restock:
   - **A DONE item collapses to one line** the review after it lands:
     `N. DONE <date> (<slot>) — <chunk/step>: <result in ≤ 25 words>;
     <commit>; <log>`. Its execution text has served its purpose — the
     closure narrative is already in the §7 row / `docs/planning/chunks/<ID>.md`
     and the journal entry. Move nothing; delete the execution text (git
     keeps it) unless it holds a ruling recorded nowhere else, in which case
     move that paragraph verbatim to the chunk's history file first.
   - **This review's own journal** — interval recap, audit reports, rulings,
     §10 assessment, restock and XL-clerk notes — is written to
     `docs/planning/reviews/YYYY-MM-DD-daily.md` (create the directory on
     first use), **not** into §9. §9 carries a digest of **≤ 15 lines**:
     the rulings that bind the next slots, one line each, and the path of
     the full journal. The dashboard (step 7) digests the same file.
   - **Standing rules stay, dated one-offs go.** A standing rule that a
     slot must obey stays in §9 verbatim. A paragraph that records *why* a
     rule exists, or an incident, moves verbatim to
     `docs/planning/plan-archive.md` with a one-line pointer.
   Target, not a gate: §9 under ≈25 KB with a full queue. If it is over,
   say so in the journal and say which rule was not enough.

7. **Refresh the status dashboard.** Rewrite `docs/status/dashboard.md` from
   what steps 1–6 established — Waiting-on-you first, then the §2 digest
   (only when §2 changed), recent activity, automation health, on-deck
   summary. Keep it a digest: no content that exists only there.

   **Then regenerate the cold-start index (added 2026-09-19):**

   ```
   scripts/testing/run_and_log.sh PLAN-INDEX "timeout -k 30 30 python3 scripts/maintenance/plan_index.py"
   ```

   It rewrites the generated open-chunk block of `docs/status/START_HERE.md`
   from the §7 rows you have just edited (through the harness: a bare
   `python3` is not on the headless allowlist). Never edit that block by
   hand, and put no status, number or ruling into the rest of that file — it
   points, it does not record. Commit it with the dashboard.

   The artifact republish is **interactive-only**: the Artifact tool is not
   available in headless scheduled sessions, so do not attempt it here —
   the file update is this step's whole deliverable. The published copy at

   `https://claude.ai/code/artifact/d5040a1e-ae6c-42dd-8e11-2330e0b9bbc8`

   is refreshed by the next interactive session (pass that URL as `url` so
   the link stays stable), and lags `docs/status/dashboard.md` until one
   runs — the operator confirmed this arrangement 2026-08-11 after the
   artifact was found six days stale. Anything blocked on the human
   operator goes at the
   top of Waiting-on-you — the dashboard is the only alerting channel; do
   not send push notifications. Dashboard staleness alone does not justify
   a commit — fold the refresh into a commit the other steps already
   earned, or skip it this interval.

   **Weekly usage-limit outages are accepted (operator ruling 2026-10-09).**
   When launcher logs read "You've hit your weekly limit", report the lost
   sessions in Automation health as information — count and window — not as
   a Waiting-on-you decision. The allowance resets Friday 12:00 local, so the
   outage falls on Thursday night / Friday morning, where no weekly runs, and
   the two-slot dirty-tree containment covers a session cut off mid-run.
   **Re-raise it to Waiting-on-you only if** a weekly review is lost to it,
   the first failure moves earlier than Thursday 18:00, a mid-run cut-off
   leaves a `recovered/*` branch, or the outage exceeds one 12-hour window.

8. Commit everything as `docs(plan): daily review YYYY-MM-DD`. If nothing
   needs changing, **commit nothing** — §5.2 explicitly prohibits audit-note
   commits, and that rule exists because of a 35-commit pile of them.

## Rubric: what a queueable item states

This is the standard the strongest items have already met — `PORT-1` step 1
named its closed form, its meshed-vs-nominal current correction, and the two
traps that had each cost a run, which is why it returned a decisive negative
inside one slot instead of a confused half-result. Below that bar, the
implementer spends its hour rescoping instead of measuring. Every item added in
step 5 or listed in step 6 states all six. An item whose chunk class has a
specialist executor (`example-runner` for EX-* example chunks, `mesh-probe`
for measurement-only mesh/resolution probes, `record-reconciler` for
version-bump record sweeps) says so in its first line:

1. **The anchor** — the specific closed form, conservation/reciprocity
   identity, or convergence rate the item will assert against, named with the
   symbol or the function that computes it (`utils/analytical.py`, a paper's
   kernel). "Check that it works" is not an anchor, and an item without one
   cannot close a chunk under §4.
2. **The negative control** — what a blind or broken solver returns on the same
   fixture, and the separation to assert. Compute the *ceiling* before naming a
   factor: `POST-3` step 2's blind imbalance saturates just under 100%, so
   1/0.1185 = 8.4× is arithmetically the most that fixture can show and a 10×
   bar would have been unreachable, not merely unmet.
3. **Tier, ranks, and expected wall clock** — smoke 30 s / standard 180 s /
   heavy 1200 s (never `xl`: that tier is the weekly review's to commission,
   §5.1), the rank count, and a cost estimate taken from a prior
   measurement where one exists (a probe's solve time, a comparable fixture).
   An item nobody has costed is an item that overruns.
4. **The traps already paid for** — name the failures this project has already
   bought, so the run does not buy them twice: `ufl.max_value` and any
   ordering comparison (`<=`) on complex-typed operands do not compile in the
   complex build — the error surfaces either as a UFL `ComplexComparisonError`
   or as a swallowed FFCx "root node" failure, and three fixture
   `current_density` callables carried it (`OPS-22`, 2026-08-19; regularise
   inside the `sqrt` instead) — and **real operands do not save a comparison
   that contains a `sqrt`**: UFL types `Sqrt` as complex whatever its
   argument, so `proj > sqrt(x² + y²)·c` on a real `SpatialCoordinate` raised
   the same error and cost a 400 s window (`WF-6` step 4c, 2026-09-08; compare
   squares — `proj > 0 ∧ proj² > (x² + y²)·c²`); a killed compile leaves an **orphan `.c`** in
   `/root/.cache/fenics` — a stem with its `.c` but no `.cpython-*.so`,
   either a 0-byte stub (2026-08-18/19, three windows) or a partial file
   (44 KB, `MAT-15` item 160, 2026-10-01, missed by the old `-size 0`
   sweep) — that is a *live lock*: later runs (even in the same session)
   stall in `MPI_Bcast` or fail with "JIT compilation timed out, probably
   due to a failed previous compile". Every complete stem has four files
   (`.c`, `.c.cached`, `.o`, `.cpython-312-…so`). **Sweep only when no
   compute is running** (`pgrep -c python3` in the container reads 0 — a
   live compile legitimately has `.c` without `.so`), as one bare harness
   call:
   `scripts/testing/run_and_log.sh <ID>-cache-sweep "docker compose exec -T fem-em-solver bash -lc 'cd /root/.cache/fenics && for c in libffcx_*.c; do s=\${c%.c}; ls \$s.*.so >/dev/null 2>&1 || { echo orphan \$s; rm -f \$s.c \$s.c.cached \$s.o; }; done; echo sweep done'"`
   — delete orphan stems only, never clear the cache wholesale (the
   targeted delete is also the diagnostic: a real defect re-creates its
   orphan, a cache artifact does not); `cell_tags.values` and
   `assemble_scalar` are rank-local; pytest captures prints without `-s`;
   `-k a or b` splits into stray argv inside an already-quoted container
   command; a headless session that backgrounds a harness run and ends its
   turn exits the CLI and SIGKILLs the harness (footerless log, no journal —
   three slots on 2026-08-10/11): harness runs go foreground, Bash-tool
   timeout 660000 ms, container-side `timeout` sized to return a footer
   inside that window — and the same rule binds a **spawned executor's**
   windows and the slot that spawned it: on 2026-09-01 00:00
   `example-runner` returned with a `./run_examples.sh` window still
   running (killed host-side 267 s in, footerless, the container process
   orphaned), and the slot that ended its turn waiting for it was
   terminated by the CLI's 600 s background ceiling with no journal —
   spawn executors foreground, state the rule in the spawn prompt; the
   container-side `timeout` needs `-k 30` — a plain
   TERM does not reliably stop an `mpiexec` job, and an overrun can wedge
   the container (MAT-6 step 10, 2026-08-12; recovery is
   `docker compose up -d --force-recreate`); piping pytest through
   `grep -v` (or anything) inside the harness command makes the log
   footer record the pipe's exit status, not pytest's — two OPS-17
   step-2 footers show exit 0 over a failing and a killed run
   (2026-08-17); filter after the fact, never in the pipeline; a
   `SpatialCoordinate`-bearing facet integral on a gmsh mesh without
   `metadata={"quadrature_degree": …}` can send FFCx into a compile that
   does not finish in nine minutes, and each killed window poisons that
   form's cache entry (`rm /root/.cache/fenics/*<hash>*` recovers; pin
   the degree — `POST-5` step 1, 2026-08-18, two windows); landing a
   parked branch by **path checkout** of an append-only record
   (`git checkout attempt/… -- docs/testing/test-results.md`) replaces
   the file and silently deletes every row `main` appended since the
   branch — `8d4cf58` (`TH-15` step 3b, 2026-09-09) dropped nine rows
   across four chunks, restored by the 2026-09-10 03:00 review; take
   code by path, append record rows by hand, and check
   `git diff --stat` shows only insertions on record files; a negative
   control that copies `HEAD:` of the file its change edits goes stale the
   moment that change lands — `test_orphan_guard.sh`'s control was red on
   every HEAD after `d10a940` and cost `OPS-43` (a)'s regression a window
   (2026-09-10 16:30); pin the pre-change commit (`<sha>^`) instead; an
   orphan-rank check spelled `pgrep -f "python3 -m pytest"` is **denied by
   `bash_guard.py`** (the string matches its pytest-routing rule) —
   `pgrep -c python3` inside the container is the check that runs
   (`TH-19` steps 1–2, 2026-09-10 22:30); a durable-capture command that
   drops the idiom's trailing **`; exit $rc`** exits with `cat`'s status, so
   the harness footer and test-results row read **0 over a failed pytest** —
   `WF-6` step 4h's red window (`[capture] rc=1`, `Status: 0`, 2026-09-11
   09:00); copy §5.1's idiom verbatim. Since `OPS-45` (2026-09-11 12:00) the
   footer honours a **final** `[capture] rc=` line, but an rc line followed
   by any further output is still ignored by design, so the trailing
   `; exit $rc` stays mandatory. A DOLFINx `.eval` / point sample's
   **bytes depend on the batch it was evaluated in**: re-sampling four
   captured batches as one 46-point batch differed in the last ULP at some
   probes and turned a W4 red (`POST-8` 3b-ii-t2, item 537, 2026-10-18
   21:00, `20261019T020649Z_POST-8-step3b-ii-t2-W4-wire.log`). A bitwise
   comparison compares call by call on the same batch. Add to this list as
   runs discover more.
5. **The scope boundary** — what the item does *not* close, stated so the
   implementer holds the chunk at 🟡 rather than over-claiming. `POST-3` step 1
   correctly stayed 🟡 because a scalar-σ identity does not gate the coil+
   phantom case, which is where it would earn its keep.
6. **What a negative result means** — the disposition if the measurement comes
   back wrong or zero. The answer is always *report the measurement and stop*,
   never fabricate a gate around it or loosen a bound to swallow it; the item
   should say which artifact captures it (a §7 annotation, a known-issues
   entry, an `attempt/*` branch). **A pre-registered negative** (operator
   directive 2026-10-01): when the item is queued *expecting* an assertion to
   fail — the confirming measurement of a ruling, as item 133 was — write it
   in this fixed form on its own line of the brief:
   `PREDICTED RED: <assertion id> range [<a>, <b>]` (one line per predicted
   red; a range, never a point). Only a red that lands inside its stated
   range, with every other assertion green, lets the slot continue
   (implementer-run.md step 2); without that line a red stops the slot as
   before.
7. **The status its result can move** *(added 2026-09-12, operator
   directive)* — the §7 chunk row, §2 sentence or known-issues entry the
   measurement can change, and what result would change it (a glyph flip, a
   retired entry, a registered record, a re-pointed §2 clause). An item that
   cannot name one is a curiosity, not work, and the review declines to queue
   it — `WF-6` steps 4h–4k all diagnosed a rung the review had already
   removed from the ladder, four slots that could not move a status by
   construction. This is a filter applied while queueing, not a regret
   afterwards.

Prefer items where a negative result is still informative — those convert a
failed hour into a finding. `PORT-1` step 1 measured exactly-zero mutual
coupling and thereby found the unfragmented mesh that had been quietly
corrupting three other fixtures.

## Constraints

- AED benchmark numbers are private (CLAUDE.md hard rules): read the gitignored `aed_results/` / `COMPARISON_private.md` / `docs/private/` files if you need context, but never copy a number from them into a tracked file, a journal entry or a commit message.
- Never loosen a test bound or a done-when to make history look better.
- **Write "approximately" as `≈`, never `~`,** in any tracked Markdown
  (plan, journals, dashboard, known-issues, chunk files). GitHub-flavoured
  Markdown renders a single `~text~` as strikethrough, so two approximate
  tildes in one paragraph cross out everything between them (swept
  2026-09-28, 1 800 lines). Literal tildes go in backticks or as `\~`;
  `~~…~~` stays reserved for deliberate retractions.
- known-issues.md discipline applies: failures observed but not fixed get an
  entry; entries leave only with the commit that fixes them.
- Your session transcript is not durable. Anything worth keeping goes into
  the repo in this session's commit.


## Standing §9 rules (moved verbatim from PROJECT_PLAN §9, 2026-10-24)

*Operator ruling 2026-10-24 (interactive session): §9 holds the queue; the standing blocks that govern how the reviews fill it live here, verbatim, and keep their numbers — "§9 rule (N)", "§9 item 5" and "ladder A1–C4" cite this section. Amend them here, not in §9.*

### §9 head items (history pointers, the feature ladder)

*(Items 1–3 — the birdcage-port lineage to `PORT-9` ✅, the 64 MHz h → 0
bracket, the `PORT-4`…`8` ordering — are history, archived verbatim
2026-09-28 in `docs/planning/plan-archive.md` § "§9 head items 1–3".)*

4. *(Operator directive 2026-09-04, the conductor lineage `TH-15` → `TH-14` → `ANS-6` — its ordering is spent: `TH-14` ✅, `ANS-6` AGREE on the loss partition, `TH-15` 🟡 on its own row, all now driven by the §10 chains; archived verbatim 2026-10-01 in `plan-archive.md` § "§9 head items 4 and 6".)*
5. **Operator directive 2026-09-04 (interactive session) — the Ansys
   feature ladder.** The boundary conditions and excitations HFSS users
   rely on for the §1 workflow that this repo lacked, ranked once by the
   operator so that no review has to re-derive the priority. **Reviews
   queue the ladder in this order; an item is skipped only when blocked
   (serial dependency not ✅, or a measured cost the box cannot pay), and
   every skip is recorded in this table with its reason, in the row's
   Note column.** Tier A is cheap linear algebra on solves that already
   pass their gates and comes *before* the conductor lineage of item 4;
   Tier B is real FEM work and comes after it; Tier C is phase-gated and
   is not queued until its phase opens.

   | Rank | Chunk | HFSS feature | Serial on | Note |
   |---|---|---|---|---|
   | **A1** | `POST-6` | Edit Sources (multi-port amplitude/phase) | — | first item; package-level version of `WF-6` step 2 |
   | **A2** | `PORT-14` | Lumped RLC boundary | — | capacitors in the model |
   | **A3** | `PORT-15` | HFSS + Circuit link | `PORT-13` ✅ | tuning at zero field-solve cost. Queued in ladder order 2026-09-04 18:00 review as §9 item 3 — **step 1 is the part that is serial on nothing**: the ladder-network closed form and the termination-reduction algebra, gated by pure-numpy identities (no FEM); gate (i) against `PORT-14`'s in-model solve waits on A2 |
   | — | `TH-15` → `TH-14` → `ANS-6` | Perfect E on bodies → Finite Conductivity → the copper AED case | — | item 4's lineage, between the tiers. `TH-15` step 0 (`mesh-probe`, 🧪 by the §3 rule until step 1 gates it) queued 2026-09-04 18:00 review as §9 item 4 |
   | **B1** | `TH-5` | Radiation boundary | — | existing ⬜ chunk, now dated by this ladder. **Skip recorded 2026-09-23 03:00 review:** not queued — Tier B follows item 4's conductor lineage, which has not closed (`TH-15` 🟡, `ANS-6` 🟡 awaiting the 09-26 adjudication); B3 (`TH-17`) went ahead only as the weekly's chain step T5. Whether B1 opens before the lineage closes is the weekly's call. **Released onto lane 2 by rule (5), 2026-09-28; 🟡 2026-09-30 (step 1a, §9 item 94); step 1b re-posed on a sphere and green (items 107, 114: Γ_sph(π) to 1.83e-3, band 6e-3; the weekly ratifies the re-pose); step 2's bar is above its predicted ceiling — the weekly's (10:30 review)** |
   | **B2** | `TH-16` | Perfect E / H symmetry planes | `GEO-25` report | the F-human memory lever. The `GEO-25` report landed 2026-09-04 16:30 (cost exponent ≈ 0.84 at fixed sizing, not r³; 504 642 cells at 0.15 m) — the weekly re-dates the 62 GiB figure from it. **✅ 2026-09-30 (§9 item 93, audited PASS 10:30 review) on the re-anchored gate: quarter = its mirror-built full 4×4 to 1.5e-12; quarter sweeps 0.36 GiB against the full's 1.24 on F-small at degree 1. Its F-human use is a chain step the weekly writes** |
   | **B3** | `TH-17` | Eigenmode solution | `TH-15`, `PORT-14` | Phase 6's first target |
   | **B4** | `TH-18` | Layered Impedance | `TH-14` | last in Tier B. **✅ 2026-09-30 (§9 item 87) on the σ = 10 S/m amended gate — the weekly ratifies** |
   | **C1** | `WF-9` | (implant TF excitation, ISO 10974) | Phase 7 open | scope `WF-7` as its consumer |
   | **C2** | `PORT-17` | Wave Port | a benchmark that needs it | low *(renumbered from `PORT-16` 2026-09-09 03:00 review — ID collision with the closed accounting-gap chunk)* |
   | **C3** | `MAT-3` | frequency-dependent materials | Phase 7 open | existing ⬜ chunk |
   | **C4** | Phase 9 AMR | adaptive passes | Phase 6 closed | not a BC; the largest workflow gap |

   What is deliberately *not* on the ladder, because it already exists:
   PEC outer box, lumped ports (`PORT-9`), incident plane wave (`TH-6`),
   projected current drives, the 4-port fixed quadrature (`WF-6` step 2).
6. *(Operator directive 2026-09-05, the XL tier — now §5.1; archived verbatim, same section.)*

### §9 rules (1)–(6) and the binding rulings

**Three rules enacted 2026-09-12 18:00 review (operator directives of
2026-09-12, `247290d`):** (1) **a slot that commits its item with a clean
tree before minute 30 takes the next open item** (implementer-run.md step 2 —
never on a dirty tree, never an item whose stated dependency has not landed,
one outcome commit and one journal entry per item); (2) **a step family is
capped at four attempts** since its last sub-step that landed a gate — a
numbered step plus all its lettered sub-steps is one family, relettering
does not reset the count, a landed gate does — and at four the family is
frozen: a review re-scopes the question as a new numbered step or banks the
measured negative, and never queues a fifth letter (this replaces the former
"items that fail twice get rescoped" sentence, which never fired); (3)
**every item below names the status its result can move** (daily-review.md
rubric element 7) — an item that cannot is not queued; (4) **enumerated-chain
licence (operator directive 2026-09-18):** when the queue is below the floor,
the daily review may write and queue, at full rubric grain, any *numbered*
step that the weekly has enumerated in a §10 chain (today: the "step count to
a tuned birdcage" list) whose serial prerequisites have landed — the weekly
supplied the question and the validation target, the review supplies the
execution detail; it may not add, reorder or re-target chain steps, and a
step whose target it cannot make assertable is queued as "write the anchor"
and returned to the weekly. Families frozen at
enactment: `WF-6` step 4 (4h–4k since 4g), `ANS-4` step 2 (2a‴, 2e, 2f, 2g
since 2a″), `PORT-14` step 2 (2, 2b–2e; no gate landed) — each carries its
ruling in §7.

**Binding ruling — sheets and κ (operator, 2026-10-02).** **Port and lumped-element sheets are exact rectangles; κ is not applied to them** (operator ruling 2026-10-02, interactive session). Evidence: on exact-rectangle ports Palace and ours agree to 0.0000 % at matched order and κ falls to ≈ 1.6e-5 / 1.8e-5 (degree 1 / 2, items 159 / 166) against ≈ 1.1e-2 on the old ragged sheets — κ measured the ragged sheet geometry, not the solver; and the κ-corrected copper F-human 32×32 is non-passive (σ_max 1.000002129) where the uncorrected one is passive (0.999993339, item 175). So: (i) **every new fixture or step builds its port / RLC / impedance sheets as exact CAD rectangles** (OCC faces — `leg_sheet_strip_fraction`, the longitudinal ring sheets), never a sub-sheet triangulated out of a box; (ii) **the κ width correction is off on exact-rectangle sheets** (`width_correction_kappa=None`), with every gate at its existing band — the passivity ceiling stands; (iii) records measured on the old ragged sheets stay valid as records; moving existing fixtures (the `ANS-4` / `ANS-6` / `ANS-7` F-small half-sheets) to exact rectangles is the weekly's to plan, re-recorded under the (1*) licence, no band moved; κ remains available as a *diagnostic* (a non-zero κ on an exact sheet is a reading, not a correction).

**Binding ruling — capacitor ESR (operator, 2026-10-06, interactive session).** **Every lumped
capacitor carries a 0.2 Ω series resistance (ESR)**: tuning capacitors on ring gaps, and any
RLC sheet that models a capacitor. It is frequency-independent and the same on every capacitor.
Why: real RF capacitors are lossy. A lossless C overstates the Q of a copper coil, whose own
conductor loss is small enough that the capacitor ESR dominates it.
(i) **Circuit layer:** a capacitor termination is `Z = r + 1/(jωC)` with r = 0.2 Ω. A feed is
driven across its gap's (r + C) branch. Uniform r adds r to every mode's eigen-impedance, so
`C_mode` (from Im λ₁) is unchanged and Re λ shifts by r.
(ii) **In-model:** a series-RC sheet impedance carrying the same r and C.
(iii) **AED:** HFSS's Lumped RLC boundary is a **parallel** R ‖ L ‖ C. A series (r, C) maps
exactly at one frequency to R_p = r(1 + Q²), C_p = C·Q²/(1 + Q²), with Q = 1/(ωCr). At
113 pF and 64 MHz, Q ≈ 110, R_p ≈ 2.4 kΩ and C_p is C to 8e-5. Entering 0.2 Ω as the HFSS R
shorts the gap. On a sweep, R_p follows 1/f² away from the conversion frequency, so a SPEC
names that frequency.
**Scope:** new steps and fixtures. Lossless-C results stay valid as records, labelled lossless.
**Not a gate change:** r never re-reads or rescues a red. H12c's identity question was posed on
lossless C and stands as posed. Whether H12c / H9 are re-posed with r (both sides at once) is
the 10-07 weekly's decision, together with which in-flight chains and SPECs adopt it (H9,
`ANS-8` step 2, the ANS SPEC template, an `esr_ohm` argument on `ring_mode_rows` and the
terminations).

**Rule (6), the chain-extension licence (operator directive 2026-10-02,
interactive session).** Slots consume a weekly's enumerated steps in about two
days, then sit idle until the next weekly (09-30 weekly: 17 of 36 slots short
of work). The daily review may therefore write the **next numbered steps of a
live §10 chain** toward that chain's stated subgoal and validation target,
tagged *(daily, date — provisional)*, at most **six per chain** between
weeklies (raised from three by operator ruling 2026-10-24, interactive
session — dashboard W8: twelve idle periods traced to chains at the cap; the
dailies' provisional steps have been ratified without reversion), and **re-pose blocked steps** within the four-attempt family cap.
It never sets or moves a band, tolerance, done-when or validation target,
never opens a chain, subgoal or phase, and never ratifies its own provisional
ruling; the weekly ratifies or reverts (daily-review.md step 6,
weekly-review.md step 3).

**Rule (5), the backlog lanes (operator directive 2026-09-28, interactive
session).** The queue ran dry on 09-27 and 09-28 while 60 chunks were open:
rule (4) is the only licence, and every enumerated chain sits on one serial
critical path gated by `xl` windows, weekly rulings and operator steps. From
2026-09-28, when the ready chain steps do not reach the floor, the daily
review **also** writes and queues, at full rubric grain (anchor, control,
tier by measurement, the status it can move), items from three standing
lanes. Each lane item must be independent of the chain critical path, and a
review fills the shortfall round-robin across the lanes, not from one:

- **Lane 1 — revalidation.** A `🧪` chunk's never-executed test run
  through the harness with a quantitative assertion (→ ✅, or a recorded
  negative plus a known-issues entry), or a `⚠️` chunk's placeholder test
  replaced by an assertion against the real solve (§3: never by loosening
  it; the ⚠️ no-new-features rule still holds). Today's list (pruned 2026-09-30 10:30
  review; closed since 09-28: `PORT-7`, `TH-4`, `PORT-8`, `GEO-3`, `PORT-4`,
  `MAT-1`, `PORT-6`, `TH-3`, `POST-2`, `PORT-5`, `WF-3`, `GEO-5`):
  `GEO-1`/`2`/`4`/`6`/`7`, `GEO-27`…`32`, `PORT-2`/`3`/`18`,
  `OPS-7`/`8`, `WF-1`/`4` (`🧪`); `TH-2` (`⚠️`; `WF-2` ✅ 2026-10-09, §9 item 346, provisional). **Every listed
  member has been examined and rejected, with reasons**
  (`docs/planning/reviews/2026-09-29-daily-1800.md` §6,
  `…2026-09-30-daily.md` §6): the lane's 🧪 / ⚠️ list is exhausted unless
  the weekly opens new members. What the 10:30 review queued in its place
  is the same kind of work on other rows: a closing window on the current
  image for a 🟡 chunk whose steps closed on the old one (`POST-3`), an
  audit's unmeasured control (`MAT-5` (D′)), and open known-issues entries
  whose own "Resolves with" a slot can execute. A probe-only chunk (§3:
  measurement steps are `🧪` until a gate lands) is revalidated by writing
  its gate.
- **Lane 2 — feature ladder Tier B, released now.** `TH-16` (symmetry
  planes), then `TH-18` (layered impedance), then `TH-5` (radiation
  boundary) — item 5's table, no longer held behind the conductor lineage's
  formal closure. Each is gated on its closed form (Jin, cited) before any
  birdcage use.
- **Lane 3 — Phases 7 and 8 in parallel, and their materials.** Phases 7
  (implants) and 8 (Pennes bioheat) are **open** from 2026-09-28 and do not
  wait for Phase 6; they consume `WF-7`'s SAR and the F-small phantom, not
  the tuned coil. `MAT-3` (Debye / Cole–Cole) and `MAT-5`
  (temperature-dependent σ) join this lane (ladder C3 released). Phase 8
  starts at the gel regime — conduction plus a prescribed source, gated on
  analytic solutions (§10 Phase 8); Phase 7 at parametric implant geometry
  in the phantom, gated on a closed form or a published measurement. Until
  the weekly has written the §10 subgoals and §7 chunks for them, the daily
  may open a Phase 7/8 chunk itself **only** where the validation target is
  a closed form or a cited published result, and names it in the dashboard
  for the weekly to ratify.

A lane item never displaces a ready chain step above it. Slots do not
write items, so the lanes work through the reviews: with three dailies a
day, each restocks its four-slot interval from chains and lanes, and a slot
that still finds the queue drained stops and journals (the `EX-57` fallback
is exhausted, `missing=0`). A drain is reported as a 🔴 dashboard item
naming the empty lanes — a supply failure, never a reason to cut slots.
