#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["mcp>=2.2,<3"]
# ///
"""
Rheome MCP server — Rheome used as a product, through its public API, the way an outside user would.

A project that needs Rheome (its scale-up analyses, the resting-state pattern atlas, their readouts) calls
these tools instead of reading or importing Rheome's source. Every call is recorded with the Rheome version
that answered it and its inputs (provenance), so a result can always say which release made it.

Tools
  version      the Rheome install the server runs (version, git describe, commit, dirty), and the server's version
  analyses     the analyses rheome.scale.run knows, and the "ported" default set
  run          rheome.scale.run: one subject, the chosen analyses -> per-subject tables       (a MATLAB job)
  reduce       rheome.scale.reduce: per-subject folders -> group distributions               (a MATLAB job)
  modes_discover  rheome.modes.discover: group discovery of the pattern atlas, frozen + hashed (a MATLAB job)
  modes_reduce    rheome.modes.reduce: the keep rule over a cohort's back-fits                 (a MATLAB job)
  job / jobs   a job's status, log tail, outputs and provenance; the job list
  readout      the tables a run wrote (metrics.csv, timing.csv, modes_catalogue.csv, ...), filtered by column
  calls        the provenance log: every tool call with its inputs and the Rheome version

MATLAB jobs run ONE AT A TIME (the hub's rule, Diellor 2026-10-11): a single worker starts the oldest queued
job when none is running. The queue lives on disk (STATE/jobs), so it survives a server restart.

Environment
  RHEOME_ROOT            the Rheome install to run: a production clone (library/software/rheome, main or a tag)
                         or an unpacked release; default: the folder above this file
  RHEOME_MCP_ALLOW_DEV   "1" allows a RHEOME_ROOT under research/code (tests only); otherwise refused
  RHEOME_MCP_STATE       jobs and the call log            (default: ~/workspace/.rheome-mcp)
  RHEOME_MCP_ROOTS       ":"-separated folders that out_dir / readout paths must be inside (default: ~/workspace)
  RHEOME_MATLAB          the matlab launcher               (default: matlab on PATH, else MATLAB_R2023b)
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import secrets
import shutil
import socket
import subprocess
import threading
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

SERVER_VERSION = "0.1.0"
HERE = Path(__file__).resolve().parent
ROOT = Path(os.path.expanduser(os.environ.get("RHEOME_ROOT", str(HERE.parent)))).resolve()
STATE = Path(os.path.expanduser(os.environ.get("RHEOME_MCP_STATE", "~/workspace/.rheome-mcp")))
JOBS = STATE / "jobs"
CALLS = STATE / "calls.jsonl"
ROOTS = [Path(os.path.expanduser(p)).resolve() for p in os.environ.get("RHEOME_MCP_ROOTS", "~/workspace").split(":") if p]
MATLAB = os.environ.get("RHEOME_MATLAB") or shutil.which("matlab") or "/Applications/MATLAB_R2023b.app/bin/matlab"

OPTION = re.compile(r"^[A-Za-z][A-Za-z0-9]*$")
ENVVAR = re.compile(r"^RHEOME_[A-Z0-9_]+$")
NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
JOBID = re.compile(r"^[0-9]{8}T[0-9]{12}Z-[0-9a-f]{6}$")
_lock = threading.Lock()


# ---------------------------------------------------------------- the install and its version

def check_root() -> None:
    """Refuse a development checkout: products are used from production copies or releases (C5)."""
    if not (ROOT / "+rheome").is_dir():
        raise RuntimeError(f"RHEOME_ROOT {ROOT} holds no +rheome package")
    if "/research/code/" in f"{ROOT}/" and os.environ.get("RHEOME_MCP_ALLOW_DEV") != "1":
        raise RuntimeError(f"RHEOME_ROOT {ROOT} is a development checkout; point it at library/software/rheome "
                           "or an unpacked release (RHEOME_MCP_ALLOW_DEV=1 is for tests)")


def _git(*a: str) -> str:
    r = subprocess.run(["git", "-C", str(ROOT), *a], capture_output=True, text=True)
    return r.stdout.strip() if r.returncode == 0 else ""


def rheome_version() -> dict:
    """The release (Contents.m) and, in a git clone, exactly which commit -- the provenance of every call."""
    m = re.search(r"% Version (\d+\.\d+\.\d+)", (ROOT / "Contents.m").read_text(errors="replace")) \
        if (ROOT / "Contents.m").exists() else None
    v = {"release": m.group(1) if m else None, "root": str(ROOT)}
    commit = _git("rev-parse", "HEAD")
    if commit:
        v.update(commit=commit, describe=_git("describe", "--tags", "--always"),
                 dirty=bool(_git("status", "--porcelain", "--untracked-files=no")))
    return v


# ---------------------------------------------------------------- provenance

def now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def log_call(tool: str, inputs: dict, result: Any) -> None:
    STATE.mkdir(parents=True, exist_ok=True)
    rec = {"at": now(), "tool": tool, "server": SERVER_VERSION, "rheome": rheome_version(), "inputs": inputs,
           "job": result.get("job") if isinstance(result, dict) else None,
           "error": result.get("error") if isinstance(result, dict) else None}
    with _lock, CALLS.open("a") as f:
        f.write(json.dumps(rec) + "\n")


def inside(path: str, must_exist: bool = False) -> Path:
    """Validate a caller's path: absolute, and inside RHEOME_MCP_ROOTS."""
    p = Path(os.path.expanduser(path))
    if not p.is_absolute():
        raise ValueError(f"{path}: give an absolute path")
    p = p.resolve()
    if not any(p == r or r in p.parents for r in ROOTS):
        raise ValueError(f"{p} is outside the allowed roots {[str(r) for r in ROOTS]}")
    if must_exist and not p.exists():
        raise ValueError(f"{p} does not exist")
    return p


def fingerprint(value: Any) -> Any:
    """sha256 + size of every input that is an existing file (folders: listed, not hashed)."""
    out = {}
    for v in value if isinstance(value, list) else [value]:
        if isinstance(v, str) and v.startswith("/") and Path(v).is_file():
            h = hashlib.sha256()
            with open(v, "rb") as f:
                for b in iter(lambda: f.read(1 << 20), b""):
                    h.update(b)
            out[v] = {"sha256": h.hexdigest(), "bytes": Path(v).stat().st_size}
    return out


# ---------------------------------------------------------------- the job queue (one MATLAB run at a time)

def submit(function: str, args: list, options: dict, env: dict, out_dir: Optional[Path], returns: bool = False) -> dict:
    check_root()
    for k in options:
        if not OPTION.match(k):
            raise ValueError(f"option {k!r} is not a MATLAB name-value name")
    for k, v in env.items():
        if not ENVVAR.match(k) or not isinstance(v, str):
            raise ValueError(f"env {k!r}: only RHEOME_* variables with string values")
    # microseconds: ids sort in submission order, which is the queue's order
    jid = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ") + "-" + secrets.token_hex(3)
    d = JOBS / jid
    d.mkdir(parents=True)
    files = {}
    for a in [*args, *options.values(), *env.values()]:
        files.update(fingerprint(a))
    call = {"job": jid, "function": function, "args": args, "options": options or None, "env": env,
            "returns": returns, "out_dir": str(out_dir) if out_dir else None, "root": str(ROOT),
            "rheome": rheome_version(), "server": SERVER_VERSION, "host": socket.gethostname(),
            "submitted": now(), "input_files": files, "matlab_launcher": MATLAB}
    (d / "call.json").write_text(json.dumps(call, indent=2))
    (d / "status").write_text("queued")
    pump()
    return {"job": jid, "status": (d / "status").read_text(), "function": "rheome." + function,
            "rheome": call["rheome"]}


def _alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except OSError:
        return False
    try:  # a zombie child of this server is not alive
        return os.waitpid(pid, os.WNOHANG) == (0, 0)
    except ChildProcessError:
        return True


def _finish(d: Path) -> None:
    """A run whose MATLAB process is gone: done / error from result.json; copy the provenance to its out_dir."""
    r = d / "result.json"
    res = json.loads(r.read_text()) if r.exists() else {"ok": False, "message": "MATLAB exited without a result (see log)"}
    (d / "status").write_text("done" if res.get("ok") else "error")
    call = json.loads((d / "call.json").read_text())
    call.update(finished=now(), matlab=res.get("matlab"), resolved=res.get("resolved"), seconds=res.get("seconds"),
                ok=bool(res.get("ok")), message=res.get("message") or "")
    (d / "provenance.json").write_text(json.dumps(call, indent=2))
    if call.get("out_dir"):
        Path(call["out_dir"]).mkdir(parents=True, exist_ok=True)
        (Path(call["out_dir"]) / f"rheome_mcp_{call['job']}.json").write_text(json.dumps(call, indent=2))


def pump() -> None:
    """Finish the running job if its MATLAB is gone; start the oldest queued one if nothing runs."""
    if not JOBS.is_dir():
        return
    with _lock:
        running = False
        for d in sorted(JOBS.iterdir()):
            if (d / "status").exists() and (d / "status").read_text() == "running":
                if _alive(int((d / "pid").read_text())):
                    running = True
                else:
                    _finish(d)
        if running:
            return
        queued = [d for d in sorted(JOBS.iterdir()) if (d / "status").exists() and (d / "status").read_text() == "queued"]
        if not queued:
            return
        d = queued[0]
        call = json.loads((d / "call.json").read_text())
        env = {**os.environ, **call["env"]}
        q = lambda s: str(s).replace("'", "''")
        cmd = [MATLAB, "-batch", f"addpath('{q(HERE)}'); rheome_mcp_call('{q(d)}')"]
        with open(d / "log.txt", "w") as log:
            p = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, env=env,
                                 cwd=str(d), start_new_session=True)
        (d / "pid").write_text(str(p.pid))
        (d / "started").write_text(now())
        (d / "status").write_text("running")


def _worker() -> None:
    while True:
        try:
            pump()
        except Exception:  # noqa: BLE001 -- the worker must outlive one bad job folder
            pass
        time.sleep(5)


def job_info(jid: str, log_lines: int = 20) -> dict:
    if not JOBID.match(jid) or not (JOBS / jid).is_dir():
        raise ValueError(f"no job {jid!r}")
    pump()
    d = JOBS / jid
    call = json.loads((d / "call.json").read_text())
    out = {"job": jid, "status": (d / "status").read_text(), "function": "rheome." + call["function"],
           "submitted": call["submitted"], "started": (d / "started").read_text() if (d / "started").exists() else None,
           "out_dir": call["out_dir"], "rheome": call["rheome"]}
    if (d / "provenance.json").exists():
        p = json.loads((d / "provenance.json").read_text())
        out.update({k: p.get(k) for k in ("finished", "seconds", "matlab", "resolved", "message")})
        if (d / "result.json").exists():
            out["value"] = json.loads((d / "result.json").read_text()).get("value")
        if call["out_dir"] and Path(call["out_dir"]).is_dir():
            out["outputs"] = sorted(str(f.relative_to(call["out_dir"])) for f in Path(call["out_dir"]).rglob("*") if f.is_file())[:200]
    if (d / "log.txt").exists() and log_lines:
        out["log_tail"] = (d / "log.txt").read_text(errors="replace").splitlines()[-log_lines:]
    return out


def wait(jid: str, timeout_s: float) -> dict:
    t = time.time() + timeout_s
    while time.time() < t:
        info = job_info(jid, 0)
        if info["status"] in ("done", "error"):
            return job_info(jid)
        time.sleep(1)
    return job_info(jid)


# ---------------------------------------------------------------- tools (plain functions; registered below)

def tool_version() -> dict:
    """The Rheome install this server runs and the server's own version. Cite `rheome.describe` (or `release`)
    in any result made through these tools."""
    check_root()
    return {"rheome": rheome_version(), "server": SERVER_VERSION, "matlab_launcher": MATLAB, "state": str(STATE),
            "roots": [str(r) for r in ROOTS]}


def tool_analyses(timeout_s: float = 300) -> dict:
    """The analyses rheome.scale.run knows (in run order) and the "ported" set it runs by default. Answered by
    MATLAB once per Rheome version, then cached."""
    v = rheome_version()
    cache = STATE / f"analyses-{v.get('commit') or v.get('release')}.json"
    if cache.exists():
        return json.loads(cache.read_text())
    a = wait(submit("scale.analyses", [], {}, {}, None, returns=True)["job"], timeout_s)
    p = wait(submit("scale.analyses", ["ported"], {}, {}, None, returns=True)["job"], timeout_s)
    if a["status"] != "done" or p["status"] != "done":
        return {"error": "MATLAB did not answer in time or failed", "jobs": [a, p]}
    res = {"all": a["value"], "ported": p["value"], "rheome": v,
           "help": "per-analysis outputs and inputs: rheome.scale.run / rheome.scale.analyses help in the release docs"}
    cache.write_text(json.dumps(res))
    return res


def tool_run(subject: str, out_dir: str, analyses: Optional[list[str]] = None, proto_dir: Optional[str] = None,
             options: Optional[dict] = None, env: Optional[dict] = None) -> dict:
    """rheome.scale.run for one subject: writes out_dir/<subject>/ (metrics.csv, timing.csv, provenance.json and
    each analysis's tables). analyses: names from `analyses` (default: the ported set). proto_dir: the Brainstorm
    protocol to import from (omit when the subject is already cached). options: any other rheome.scale.run
    name-value (Cohort, Dataset, DurationS, NoiseStudy, FlowTiles, Modality, Edf, Events, Import, Sub).
    env: RHEOME_* variables some analyses need (RHEOME_DATA, RHEOME_G1, RHEOME_MODES_CENTROIDS, ...).
    Returns a job id at once; follow it with `job`."""
    if not NAME.match(subject):
        raise ValueError(f"subject {subject!r}: letters, digits, . _ - only")
    od = inside(out_dir)
    od.mkdir(parents=True, exist_ok=True)
    opts = {**(options or {}), "OutDir": str(od)}
    if analyses:
        opts["Analyses"] = list(analyses)
    if proto_dir:
        opts["ProtoDir"] = str(Path(os.path.expanduser(proto_dir)).resolve())
    return submit("scale.run", [subject], opts, env or {}, od / subject)


def tool_reduce(in_dir: str, out_dir: str, options: Optional[dict] = None) -> dict:
    """rheome.scale.reduce: in_dir holds one folder per subject (run's out_dir); writes metrics_long.csv,
    status.csv, distribution.csv, anchor.csv, cohorts.csv, twins.csv to out_dir. options: Anchor, AnchorDir,
    Reference, TwinSide. Returns a job id."""
    i, o = inside(in_dir, must_exist=True), inside(out_dir)
    o.mkdir(parents=True, exist_ok=True)
    return submit("scale.reduce", [str(i), str(o)], options or {}, {}, o)


def tool_modes_discover(files: list[str], out_dir: str, options: Optional[dict] = None) -> dict:
    """rheome.modes.discover: the participants' modes_fields.mat (from run with analyses ["modesfields"]) ->
    the group's frozen, hashed centroids in out_dir. options: Variant, Bands, NPer, Ks, KMin, NSplits, RThresh,
    NullSources, Seed, Workers, ... (see the release's help). Returns a job id."""
    fs = [str(inside(f, must_exist=True)) for f in files]
    o = inside(out_dir)
    o.mkdir(parents=True, exist_ok=True)
    return submit("modes.discover", [fs, str(o)], options or {}, {}, o)


def tool_modes_reduce(fit_files: list[str], centroids: str, out_dir: str, options: Optional[dict] = None) -> dict:
    """rheome.modes.reduce: the keep rule over a cohort's modesfit.csv files against the frozen centroids ->
    modes_catalogue.csv, modes_tests.csv, reduce_provenance.json in out_dir. options: Cohort, Alpha, N1R,
    StandingRad, Fdr, NBoot, Seed, N3Margin. No label makes a propagation call (travelling is unresolved at
    >= 100 mm). Returns a job id."""
    fs = [str(inside(f, must_exist=True)) for f in fit_files]
    c, o = inside(centroids, must_exist=True), inside(out_dir)
    o.mkdir(parents=True, exist_ok=True)
    return submit("modes.reduce", [fs, str(c), str(o)], options or {}, {}, o)


def tool_job(job: str, log_lines: int = 20) -> dict:
    """A job: status (queued / running / done / error), its Rheome version, MATLAB release, the function's
    resolved path (proof of which install ran), seconds, outputs, error message and the log tail."""
    return job_info(job, log_lines)


def tool_jobs(status: Optional[str] = None, limit: int = 50) -> dict:
    """The jobs, newest first; status filters (queued, running, done, error)."""
    pump()
    out = []
    for d in sorted(JOBS.iterdir(), reverse=True) if JOBS.is_dir() else []:
        s = (d / "status").read_text() if (d / "status").exists() else "?"
        if status and s != status:
            continue
        c = json.loads((d / "call.json").read_text())
        out.append({"job": d.name, "status": s, "function": "rheome." + c["function"], "out_dir": c["out_dir"]})
        if len(out) >= limit:
            break
    return {"jobs": out}


def tool_readout(path: str, table: Optional[str] = None, where: Optional[dict] = None, limit: int = 200) -> dict:
    """Read a run's outputs. path: an output folder (or a .csv). Without table: the tables (.csv) and files there,
    and provenance.json if present. With table (e.g. "metrics", "timing", "modes_catalogue"): its columns and the
    rows matching `where` (column -> value, or list of values), at most `limit`; `total` counts every match."""
    p = inside(path, must_exist=True)
    if p.is_dir() and not table:
        files = sorted(str(f.relative_to(p)) for f in p.rglob("*") if f.is_file())
        out = {"path": str(p), "tables": [f[:-4] for f in files if f.endswith(".csv")], "files": files[:500]}
        if (p / "provenance.json").exists():
            out["provenance"] = json.loads((p / "provenance.json").read_text())
        return out
    f = p if p.is_file() else p / (table if table.endswith(".csv") else table + ".csv")
    f = inside(str(f), must_exist=True)
    want = {k: {str(x) for x in (v if isinstance(v, list) else [v])} for k, v in (where or {}).items()}
    rows, total = [], 0
    with f.open(newline="") as fh:
        r = csv.DictReader(fh)
        missing = set(want) - set(r.fieldnames or [])
        if missing:
            raise ValueError(f"{f.name} has no column(s) {sorted(missing)}; columns: {r.fieldnames}")
        for row in r:
            if all(row[k] in vs for k, vs in want.items()):
                total += 1
                if len(rows) < limit:
                    rows.append(row)
    return {"table": str(f), "columns": r.fieldnames, "total": total, "rows": rows}


def tool_calls(limit: int = 20, tool: Optional[str] = None) -> dict:
    """The provenance log, newest first: every tool call with its inputs, the server and Rheome versions."""
    lines = CALLS.read_text().splitlines() if CALLS.exists() else []
    recs = [json.loads(x) for x in reversed(lines)]
    return {"calls": [r for r in recs if not tool or r["tool"] == tool][:limit]}


TOOLS = {"version": tool_version, "analyses": tool_analyses, "run": tool_run, "reduce": tool_reduce,
         "modes_discover": tool_modes_discover, "modes_reduce": tool_modes_reduce, "job": tool_job,
         "jobs": tool_jobs, "readout": tool_readout, "calls": tool_calls}
ACTS = {"run", "reduce", "modes_discover", "modes_reduce"}


def call(name: str, **kw) -> dict:
    """One tool call, logged with its inputs and the Rheome version whatever the outcome."""
    try:
        res = TOOLS[name](**kw)
    except (ValueError, RuntimeError) as e:
        res = {"error": str(e)}
    if name != "calls":
        log_call(name, kw, res)
    return res


def build_server():
    import functools
    from mcp.server.mcpserver import MCPServer
    from mcp.types import ToolAnnotations
    srv = MCPServer("rheome", version=SERVER_VERSION, instructions=(
        "Rheome (the multiscale geometry of human cortical dynamics) as a product: run its analyses, the "
        "resting-state pattern atlas (modes) and read their tables through these tools, never from its source. "
        "Start with `version` and `analyses`; run / reduce / modes_* return a job id (MATLAB runs one at a time) "
        "to follow with `job`; read results with `readout`. Every call is recorded with the Rheome version and "
        "inputs (`calls`); cite the version in anything made from them."))
    for name, fn in TOOLS.items():
        @functools.wraps(fn)
        def tool(*a, _n=name, **kw):
            return call(_n, **kw)
        tool.__name__ = name
        read = name not in ACTS
        srv.tool(name=name, annotations=ToolAnnotations(readOnlyHint=read, destructiveHint=False,
                                                        idempotentHint=read, openWorldHint=False))(tool)
    return srv


if __name__ == "__main__":
    check_root()
    threading.Thread(target=_worker, daemon=True).start()
    build_server().run()
