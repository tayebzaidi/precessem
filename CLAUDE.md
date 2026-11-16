# FEM-EM Solver — agent onboarding

FEniCSX/DolfinX FEM toolkit targeting the MRI-RF-safety slice of Ansys
Electronics Desktop (HFSS + circuit solver first; Pennes bioheat long-term):
construct birdcage coil + gelled saline phantom (± implant) simulations, tune
the coil at 64/128 MHz, extract B1+/SAR/S-parameters. Scope rescoped
2026-08-04 — see PROJECT_PLAN.md §1; parity claims are per-workflow, never
"HFSS parity".
Magnetostatics is validated against closed forms. The time-harmonic path is a
real complex curl-curl solve validated against the analytic lossy plane wave
(`TH-1` closed 2026-07-31; decay/phase constants to < 0.06%, `TH-6`) and,
at the Larmor frequencies themselves, against the lossy-sphere series on an
imposed field (`TH-10` closed 2026-08-13: 3.643% / 1.826% at 64 / 128 MHz,
power to 3.629%; the 128 MHz figure is 1.769% on the 0.11 image `main`
boots — re-recorded with its mesh by `OPS-18`, 2026-08-22). Coil loading is gated only in the eddy-current regime —
`MAT-6` closed 2026-07-31 with ΔR matching Dodd–Deeds to 1.58% at 10 MHz,
σ = 100 S/m; the **coil-at-Larmor** case is measured flat in f across
10–64 MHz on one XL window (`TH-11` step 5d, 2026-09-09) — a record, not a
gate (2026-09-13 weekly). SAR is gated against the lossy-sphere closed form to 3.5%
(`MAT-4` step 1, 2026-08-03) on an **imposed uniform field**, never on a
coil. S-parameters are field-derived through the package (`PORT-1` closed
2026-08-15; feed systematics named on the two-torus fixture), and the
gapped 4-leg birdcage carries four lumped-sheet ports whose 4×4 passes
reciprocity / passivity / C4-symmetry gates at 10, 64 and 128 MHz
(`PORT-9` ✅ 2026-08-25, `PORT-11` ✅ 2026-08-26) — **self-consistency
identities on one fixture, not an absolute-accuracy, resonance or tuning
claim**. **On the F-small birdcage the degree-2 S-matrix and the tuned
resonance are h-statements** (AED values private; relative levels only, per
the 2026-09-24 ruling). On **exact-rectangle ports with κ off** (operator
ruling 2026-10-02) the degree-2 `ANS-4` 4×4 is within **0.73 %** of AED First
Order on every class at 10 / 64 / 128 MHz and h-converged to ≈ 0.1 % on three
mesh rungs (T6 decided 2026-10-05, F27 2026-10-12 weekly); Palace — the public
external reference (`OPS-61`) — matches ours to 0.0000 % at matched order.
The tuned `f₀` (`ANS-7`, T13, 2026-10-10 weekly) is within **0.11 %** of AED
First Order on the base mesh and 0.03 % at h = 0.005, converging first order
in h. `R_in` and Q also converge in h (T41, 2026-10-17 weekly: successive
ratios ≈ 0.71 on three named steps), to limits ≈ 3 % below AED First Order
that Palace shares to six digits — so `ANS-7` is ✅ with the absolute Q
**INCONCLUSIVE**, not DISAGREE. **Licensed:** the degree-2 exact-rectangle
S / `S₁₁`, `f₀`, and Q / `R_in` as h-converged levels on F-small; **the
F-small side of the Phase 6 exit is decided.** The
old `S₁₁` / `f₀` gaps were the ragged port sheet (which κ had been
correcting) plus degree-1 discretisation. The two-torus current-drive route's
S defect is fixed (`PORT-20` ✅ 2026-09-20). Copper via Leontovich is
externally checked: `ANS-6` ✅ (2026-10-12 weekly) within 0.33 % of AED First
Order on every class at three frequencies at degree 2.
Human scale: a 32-port F-human birdcage has symmetry-gated `|B₁⁺|` and a
located 10 g SAR hotspot at 64 MHz (`WF-7` ✅ 2026-09-25 — degree 1, σ =
800 S/m conductor). Its copper twin (`WF-12`, Leontovich hole mesh) is
passive and loss-gated; under the clinical two-port 90° drive its mode
content **agrees with AED within 0.6 %** and the degree-2 `C_mode` / `f₀`
within 0.7 % (`ANS-8` ✅ on step 1, 2026-10-12 weekly; step 2, with the
capacitor ESR in the solved field, is the operator's). `WF-10` ✅
(2026-10-11) puts the copper F-human on its in-model Larmor root with
r = 0.2 Ω ESR and feed-bisector mirror gates at degree 1; at degree 2 the
root moves to C_in = 115.5209 pF (H26, 2026-10-17) and `ANS-8` step 2's AED
session is set there (SPEC v4). **Phase 6 exited 2026-10-24** (entry 41 read
green): the human-scale tuned figures are quoted at degree 2 only, at that
root, as internal tuning statements; **no absolute, external or compliance
claim at human scale** — `ANS-8` step 2 (W2) is the external reading. New human-scale
work moves to **F-clinical** (0.70 m bore, 32 legs / 64 ring-gap sheets; 1.5 T
first, then the 3 T twin; operator directive 2026-10-14, §10 chain C). Beside Phase 6 run Phases 7–9:
Phase 8 heating of the gel ball in closed form and FEM — conduction,
SAR-driven, Pennes, duty-cycled, container-wall, CEM43 dose (`MAT-9` …
`MAT-20`), and coil-driven SAR heating the gel phantom (`MAT-15` ✅), with
**the first external thermal check**: `ANS-9` (HFSS + IcepakFEA, gel,
one-way) is ✅ (2026-10-14 weekly): on the refined degree-2 column it agrees
on every row — the hand-off integrals within 1.3 %, probes A–D within 2.5 %
(probe D's earlier gap was our degree-1 discretisation). Phase 7 PEC implants against exact series (sphere,
spheroid tip, cylinder, coated lead — `IMP-1` … `IMP-5`) and a PEC rod inside
the coil-driven phantom (`IMP-3`, cross-checked against Palace at matched
order) — **no tissue model and no external implant check yet**. Phase 9 (a
usable package: `Study` API, `femem run` YAML, tutorials as tests) opened
2026-10-03; the rename to **PrecessEM** (`WF-18`) is placed as one
interactive atomic commit.
Read PROJECT_PLAN.md §2 before trusting any S-parameter or
coil-loading/SAR figure. Anything
that solves in the frequency domain needs the complex DolfinX build
(`source /usr/local/bin/dolfinx-complex-mode`); real mode raises.

## Read these before working

- **docs/status/START_HERE.md** — read this FIRST (≈10 KB): what to read for
  which job, where everything lives and how big it is, and a generated
  one-line index of every open chunk with its plan line number. The files
  below total > 2 MB; none of them is meant to be read whole — grep a chunk
  ID and read the span.
- **PROJECT_PLAN.md** — the single source of truth. §2 honest current state,
  §4 definition of done, §5 compute budget + Docker + logging harness,
  §7 chunk backlog (stable IDs; several entries carry ready-to-execute
  implementation plans).
- **docs/testing/known-issues.md** — check BEFORE debugging any failing test.
  A set of failures on main predates you; that file says which.
- **docs/references/** — Jin, *The FEM in Electromagnetics* 3rd ed. is
  available as searchable per-chapter markdown under
  `docs/references/jin-fem-3e/` (gitignored, operator-provided; start at
  `INDEX.md`). Consult it before re-deriving formulation, feed/port-model,
  boundary-condition, or solver theory; cite chapter/section/equation
  numbers in plan annotations. Absent on a fresh clone — see
  `docs/references/README.md`.

## Hard rules

- **Shared 36-core box, 12 cores max for this project.** Declare a tier
  (smoke 30 s / standard 180 s / heavy 1200 s), wrap runs in `timeout` at the
  ceiling, never exceed 20 minutes for a single compute command.
  Overrun ⇒ kill and shrink the case; never just raise the timeout.
  **Two exceptions, mechanically gated** (PROJECT_PLAN §5.1): the **`xl`** tier —
  six runs per trailing 7 days (nightly Sun–Fri, 2026-09-13) at ≤ 512 GiB / 16 ranks / 4 h against
  `fem-em-solver-xl`, 02:00 Sun–Fri — and the **`xxl`** tier — one run per
  trailing 7 days at ≤ 754 GiB / 16 ranks / 8 h against `fem-em-solver-xxl`,
  02:00 Saturday. Both are commissioned by the weekly planning review —
  pre-registered in `docs/testing/xl-pending.md` with the exact command,
  measured price and readout — queued into `docs/testing/<tier>-queue.d/`
  (one file per window, drained in order by cron) by the daily review acting
  as clerk when the budget allows (2026-09-13), run by cron with no Claude
  session, and recorded in `docs/testing/<tier>-ledger.md`. Since
  2026-09-15 the daily review also keeps a backlog floor (≥ 4 `xl`, ≥ 1
  `xxl` ahead) and may fill it with priced-family variants and cost probes
  (daily-review.md step 6b). Implementers never commission one.
- **All verification runs in Docker through the logging harness** (service must
  be Up — `docker compose -f docker/docker-compose.yml ps`):

  ```
  scripts/testing/run_and_log.sh <CHUNK-ID> "docker compose exec -T fem-em-solver \
    bash -lc 'cd /workspace && PYTHONPATH=/workspace/src timeout -k 30 180 \
    mpiexec -n 2 python3 -m pytest <paths> -v --tb=short'"
  ```

  The `-k 30` is mandatory: a plain `timeout` TERM does not reliably stop an
  `mpiexec` job (MAT-6 step 10, 2026-08-12 — a `timeout 590` run burned cores
  for ≈1 700 s and wedged the container; recovery is
  `docker compose -f docker/docker-compose.yml up -d --force-recreate`, see
  known-issues).

- `mpiexec -n 12` is the hard ceiling (16 only inside the XL slot above); use the smallest rank count that fits the
  tier, and keep `-n 2` for anything a rank-local bug could hide in (that is the
  only width where a missing reduction is visible in CI).
- **A window over ≈10 minutes must survive its own wrapper** (§5.1, measured
  2026-09-09): redirect container-side output to a file under the gitignored
  `/logs/` first and echo it back, never a bare pipe or `tee`. If a window is
  killed, **check for orphaned ranks before anything else** — the wrapper
  dying does not stop the compute, and eight ranks once kept running on
  260 GiB with nothing consuming their output.
- A chunk is ✅ only per §4: verification executed by the agent itself, at
  least one quantitative assertion (closed form, convergence rate, or a
  conservation/reciprocity identity), elapsed time recorded. Finiteness-only
  checks close nothing.
- Never loosen a failing assertion to make a test pass; a failing analytic
  comparison is evidence about the test as much as the code (PROJECT_PLAN §7,
  MAG table, defect 5).
- **Port / lumped-element sheets are exact CAD rectangles, and the κ width
  correction is not applied to them** (operator ruling 2026-10-02; PROJECT_PLAN
  §9 binding rulings): new fixtures never triangulate a sub-sheet out of a box;
  κ stays a diagnostic. Ragged-sheet records stay valid as records.
- Rank-safety: `cell_tags.values`, `assemble_scalar`, local max/min are
  rank-local — reduce before asserting. Point evaluation goes through
  `post.evaluation.evaluate_vector_field_parallel`, never
  `f.eval(points, np.arange(n))`.
- **Ansys (AED) benchmark numbers are private and never enter a tracked
  file** — licence terms restrict disclosure. They live only in the gitignored
  `examples/ansys_benchmarks/*/aed_results/`, `*/COMPARISON_private.md` and
  `docs/private/`; quote none of them (R/X/ΔZ, tet counts, passes, energy
  errors, timings) in plan annotations, journals or commit messages — the
  qualitative verdict (agree / disagree / inconclusive) goes public, and
  (operator ruling 2026-09-24) so may a **relative agreement level** — "agrees
  to within x %", our-vs-AED gap in % — though never the AED value itself.
- **Palace is the public external reference** (operator directive
  2026-10-02, PROJECT_PLAN §5.4): its numbers are published at full detail
  (Palace version, config and mesh cited; matched element order on a shared
  mesh), but no tracked file ever sets a Palace number against an AED one —
  the three-way reading stays in the gitignored private files.
- In tracked Markdown write "approximately" as `≈`, never `~` — GFM renders
  a single `~text~` as strikethrough, so two approximate tildes in one
  paragraph cross out everything between them (swept 2026-09-28). `~~…~~`
  is for deliberate retractions only.
- Commit code, tests, harness logs, and PROJECT_PLAN/known-issues updates
  together.
- The project permission allowlist (.claude/settings.json) funnels all compute
  through `run_and_log.sh` / `docker compose` and blocks host package managers,
  network commands, and reads of files outside this repo. Prefer Read/Grep/Glob
  tools over shell readers. For multi-line commit messages use
  `git commit -F <file>` — `$(...)` substitution in commands is denied.

## Chunk execution

Implementation work is delegated per chunk to the `implementer` agent
(.claude/agents/implementer.md, pinned to Opus) — one chunk ID per invocation:
"Use the implementer agent to execute MAG-13." Planning, chunk design, and
review of returned diffs stay in the main (Fable) session.

## Scheduled automation

System cron runs headless sessions via `scripts/automation/` — **full
schedule restored 2026-09-28** (operator directive; the 2026-09-13 wind-down
is over). Every day: three daily reviews, each followed by a four-slot
implementer interval (twelve slots a day, eleven on weekly evenings), plus
three weekly reviews (Mon / Wed / Sat). Canonical crontab: `scripts/automation/crontab`.

- **Plan review** (Opus 5.5, high effort, 03:00 / 10:30 / 18:00 local, every day) — audits
  results against §4, rescopes failed attempts, disposes of `recovered/*`
  branches, tops the §9 "On deck" queue to the slot-minutes floor (160 **measured**
  minutes since 2026-10-01; item briefs live in `docs/planning/chunks/`, §9
  holds two-line entries; a slot continues past a `PREDICTED RED` that lands
  in range; since 2026-10-02 it may also write up to six (three until the operator's
  2026-10-24 ruling) provisional next steps per live §10 chain and re-pose blocked steps, never a target — §9
  rule (6), ratified by the weekly) — since
  2026-09-18 by executing the weekly's enumerated §10 chain steps at full
  rubric grain (§9 rule (4)) and, since 2026-09-28, three backlog lanes —
  🧪/⚠️ revalidation, feature-ladder Tier B, Phases 7/8 (open in parallel
  with Phase 6) — (§9 rule (5)), never by inventing validation targets.
  The weekly keeps ≥ 4–5 lanes live and enough steps to reach the next
  weekly; a drained queue is never a reason to cut slots. Protocol:
  docs/automation/daily-review.md.
- **Weekly planning review** (Fable 5.1, high effort, **Monday, Wednesday
  and Saturday 21:00 local**) —
  owns the long horizon with brutal realism: §6 phase map and §10 roadmap
  (phases → subgoals, dated pace-based assessments), examples/ health, and
  commissioning/adjudicating `examples/ansys_benchmarks/` cases the human
  operator replicates in Ansys Electronics Desktop (PROJECT_PLAN §5.4).
  Never edits §9. Three times weekly since 2026-10-02 (gaps 48 / 48 / 72 h;
  twice weekly from 2026-09-28; once weekly
  2026-09-13 … 09-26 — what the plan calls "the 09-16 weekly" is the 09-19
  one; at 21:00, not the old Sunday/Wednesday 02:15, because that session
  twice exhausted the account's 5-hour allowance and took the next block
  with it); "this week" in that
  protocol means since the last weekly-review commit, so pace arithmetic
  uses the measured elapsed time. Protocol: docs/automation/weekly-review.md.
- **Implementer runs** (Opus 5.5, 00:00 / 04:30 / 06:00 / 07:30 / 09:00 /
  12:00 / 13:30 / 15:00 / 16:30 / 19:30 / 21:00 / 22:30 local, every day;
  no 21:00 run on Monday, Wednesday and Saturday, the weekly's window) — each takes the top On-deck item inside a 1-hour timebox, and the
  next one if it commits before minute 30 (take-next, 2026-09-12); incomplete
  work is parked on `attempt/*` branches and journaled in
  docs/testing/attempts.md, never left on main. A dirty tree stops the first
  run that meets it and is parked on `recovered/*` by the second, so an
  outage costs two slots rather than the rest of the day. Protocol:
  docs/automation/implementer-run.md.
- **XL / XXL windows** (02:00, cron scripts, no Claude session, zero tokens)
  are unchanged — the cheapest compute the project has.

The daily review also maintains a status dashboard for the human operator
(`docs/status/dashboard.md`, republished as a Claude artifact — URL in
daily-review.md step 7); its Waiting-on-you section is the only alerting
channel — scheduled sessions never send push notifications. A PreToolUse
hook (`scripts/automation/hooks/bash_guard.py`, wired in
.claude/settings.json) mechanically denies `mpiexec` rank counts above 12,
pytest runs that bypass the logging harness, and (since 2026-09-30) any
`docker` / `run_and_log.sh` call that is not a **bare, whole command** — no
`cd … &&` prefix, no `; date` / `| tail` suffix, no redirect — because only
whole commands matching `docker *` / `scripts/testing/run_and_log.sh *` are
sandbox-exempt; anything else is denied the docker socket.

**Review model.** Both reviews take their model from
`scripts/automation/review-model.env` — a default, standing per-review
choices (since 2026-09-28, operator: **daily on `claude-opus-5-5`, weekly on
`claude-fable-5-1`**), plus an optional **dated** override the launcher
expires by itself, so a borrowed model cannot outlive its reason and no
revert has to be remembered. No override is set now; the launcher logs which
model it chose on every run. Implementer runs do not read
that file and always run on Opus (`claude-opus-5-5` since 2026-09-22, pinned
in `implementer-run.sh`).

If you are one of these scheduled sessions, your protocol document is
authoritative; read it before acting. If you are an interactive session,
expect these runs to exist: check `git log` and attempts.md before assuming
the tree is as you left it.
