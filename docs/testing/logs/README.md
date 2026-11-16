# Harness logs: a representative subset

Every verification run in this project goes through
`scripts/testing/run_and_log.sh`, which writes one log recording the commit,
environment, exact command, full output, exit status and elapsed time. The
private development repository holds about 4,700 of them (about 74 MB). This
public snapshot keeps 67:

- **every log cited in `PROJECT_PLAN.md` §2.1** ("Validated, with the number
  that licenses it"), so each headline validation claim can be traced to the
  run that produced it;
- **the oversight logs:** the private-data leak audits (`OPS-46`, `OPS-47`,
  `OPS-58`, `LEAK-CHECK-*`) and the `MAT-6` step 10 runs, where a plain
  `timeout` failed to stop an MPI job and led to the mandatory `timeout -k 30`
  rule (see `AI_WORKFLOW.md`).

Logs are gzipped; read them with `zcat` or `zless`. Elsewhere in the plan,
the chunk files, the attempt journal and `test-results.md`, log names are
cited that are not included here. Those citations are records of runs in the
private repository. Tooling that checks every cited log exists (the
doc-reference census) will report them as missing in this snapshot.
