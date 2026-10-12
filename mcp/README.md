# Rheome MCP server — Rheome used as a product

Other projects use Rheome the way an outside user would: through this MCP server over its public API, never by
reading or importing its source (lab architecture, Part C5). Every call records the Rheome version that answered
it and its inputs, so any result made through it can say which release made it.

## Tools

| tool | what it does | MATLAB |
|:--|:--|:--|
| `version` | the Rheome install the server runs: release (`Contents.m`), `git describe`, commit, dirty flag; the server's version | no |
| `analyses` | the analyses `rheome.scale.run` knows, and the "ported" default set (cached per Rheome commit) | once per version |
| `run` | `rheome.scale.run`: one subject, chosen analyses → `out_dir/<subject>/` (`metrics.csv`, `timing.csv`, `provenance.json`, each analysis's tables) | a job |
| `reduce` | `rheome.scale.reduce`: per-subject folders → group `distribution.csv`, `status.csv`, `anchor.csv`, `cohorts.csv`, `twins.csv` | a job |
| `modes_discover` | `rheome.modes.discover`: participants' `modes_fields.mat` → the frozen, hashed centroids of the resting-state pattern atlas | a job |
| `modes_reduce` | `rheome.modes.reduce`: the keep rule over a cohort's `modesfit.csv` → `modes_catalogue.csv`, `modes_tests.csv` | a job |
| `job`, `jobs` | a job's status (queued / running / done / error), Rheome version, MATLAB release, the function's resolved path, outputs, error, log tail | no |
| `readout` | a run's tables: without `table` the files in a folder (and its `provenance.json`); with `table` its columns and the rows matching `where` (column → value or list) | no |
| `calls` | the provenance log: every tool call, its inputs, server and Rheome versions | no |

The atlas's per-participant stages (`modesfields`, `modesfit`, `jointfields`, …) are analyses of `run`
(`analyses=["modesfit"]`, `env={"RHEOME_MODES_CENTROIDS": …}`); `options` passes any other name-value of the
function and `env` any `RHEOME_*` variable an analysis needs (`RHEOME_DATA`, `RHEOME_G1`, …). The function's own
help in the release documentation says what each option and output means.

A typical session: `version` → `analyses` → `run` (a job id) → `job` until `done` → `readout` with
`table="metrics"` and `where={"analysis": "resolution"}` → cite `rheome.describe` in the result.

## Provenance

- `STATE/jobs/<job>/call.json` — function, positional args, options, `RHEOME_*` env, sha256 and size of every
  input file, the Rheome version (release, describe, commit, dirty), server version, host, time.
- `STATE/jobs/<job>/provenance.json` — the same plus MATLAB release, `which()` of the function (the install that
  really ran), seconds, ok / error message; copied into the job's output folder as `rheome_mcp_<job>.json`.
- `STATE/calls.jsonl` — one line per tool call, reads included.

## One MATLAB run at a time

Jobs run as `matlab -batch` one at a time (the hub's rule, 2026-10-11): a single worker starts the oldest queued
job when none runs. The queue is on disk, so it survives a server restart. On the hub `mem-guard.sh` also pauses
this server's run while another agent's batch run is going.

## Configuration

| variable | default | |
|:--|:--|:--|
| `RHEOME_ROOT` | the folder above `mcp/` | the Rheome install: a production clone (`~/workspace/library/software/rheome`) or an unpacked release. A development checkout (`research/code`) is refused |
| `RHEOME_MCP_STATE` | `~/workspace/.rheome-mcp` | jobs and the call log |
| `RHEOME_MCP_ROOTS` | `~/workspace` | `:`-separated folders every `out_dir` and `readout` path must be inside |
| `RHEOME_MATLAB` | `matlab` on PATH, else MATLAB R2023b | the MATLAB launcher |
| `RHEOME_MCP_ALLOW_DEV` | unset | `1` allows a development checkout (tests only) |

Register it with a client (Claude Code `.mcp.json`; [uv](https://docs.astral.sh/uv/) installs the `mcp` package):

```json
{ "mcpServers": { "rheome": {
    "command": "uv",
    "args": ["run", "--quiet", "--script", "/Users/you/workspace/library/software/rheome/mcp/rheome_mcp.py"] } } }
```

## Tests

```sh
uv run --with "mcp>=2.2,<3" --with pytest pytest mcp/tests                 # fake MATLAB, seconds
RHEOME_MCP_MATLAB_TEST=1 uv run --with "mcp>=2.2,<3" --with pytest \
    pytest mcp/tests -k "real or analyses"                                  # the real release through MATLAB
```

`mcp/` is not part of the `.mltbx` (the toolbox build packages `+rheome`, `doc` and the root files only).
