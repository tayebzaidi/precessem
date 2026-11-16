# XL / XXL queue directories

`xl-queue.d/` and `xxl-queue.d/` are the FIFO queues `scripts/automation/xl-run.sh`
reads at 02:00 (operator directive 2026-09-15). Each `*.env` file holds one
window in the `XL_CHUNK` / `XL_COMMAND` format of the legacy `<tier>-queue.env`
(still honoured first when non-empty). Files run in lexical order, so name them
`<intended run date>-<chunk>.env`; the launcher checks the ledger budget before
taking a file, runs the first valid one, and deletes it afterwards whether the
window succeeded or failed. Entries are copied **verbatim** from
`xl-pending.md` by the daily review (daily-review.md step 6b), the weekly, or
the operator — never authored here.

**Optional `XL_OUTPUTS`** (`OPS-64`, 2026-09-28): space-separated repo-relative
paths the window writes *into the tree* that belong in its commit. The launcher
refuses the entry before the window (queue and ledger untouched) if a declared
path is absolute, outside the repo, or gitignored, and stages the declared
paths that exist after it. An undeclared tree output is left uncommitted and
makes the next window and the next slot stop on a dirty tree, so a
pre-registered entry that writes into the tree must declare it. Output under
the gitignored `/logs/` (the usual raw-log route) needs no declaration. The
launcher also records the container's python3 count before and after every
window in a `## Launcher` section of the committed harness log, so an entry
no longer needs its own `[orphans-before]` / `[orphans-after]` lines.
