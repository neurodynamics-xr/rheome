"""Tests of the Rheome MCP server. Unit tests use a fake MATLAB (fake_matlab.py); the MATLAB test runs the real
release functions and is opt-in (RHEOME_MCP_MATLAB_TEST=1: one `matlab -batch` per job, one at a time).

    uv run --with "mcp>=2.2,<3" --with pytest pytest mcp/tests
"""
import asyncio
import importlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
sys.path.insert(0, str(HERE.parent))


@pytest.fixture
def srv(tmp_path, monkeypatch):
    monkeypatch.setenv("RHEOME_ROOT", str(REPO))
    monkeypatch.setenv("RHEOME_MCP_ALLOW_DEV", "1")
    monkeypatch.setenv("RHEOME_MCP_STATE", str(tmp_path / "state"))
    monkeypatch.setenv("RHEOME_MCP_ROOTS", str(tmp_path))
    if os.environ.get("RHEOME_MCP_MATLAB_TEST") != "1":
        monkeypatch.setenv("RHEOME_MATLAB", str(HERE / "fake_matlab.py"))
    else:
        monkeypatch.delenv("RHEOME_MATLAB", raising=False)
    import rheome_mcp
    return importlib.reload(rheome_mcp)


def settle(m, jid, timeout=600):
    t = time.time() + timeout
    while time.time() < t:
        j = m.call("job", job=jid)
        if j["status"] in ("done", "error"):
            return j
        time.sleep(0.2)
    raise AssertionError(f"job {jid} did not finish: {j}")


def test_version_records_release_and_commit(srv):
    v = srv.call("version")["rheome"]
    assert v["release"] and v["root"] == str(REPO)
    assert len(v["commit"]) == 40 and "describe" in v


def test_development_checkout_refused(srv, monkeypatch):
    if "/research/code/" not in str(REPO) + "/":
        pytest.skip("checkout not under research/code")
    monkeypatch.delenv("RHEOME_MCP_ALLOW_DEV")
    r = srv.call("run", subject="sub01", out_dir=str(Path(srv.ROOTS[0]) / "o"))
    assert "development checkout" in r["error"]


def test_paths_outside_roots_and_relative_refused(srv):
    assert "outside the allowed roots" in srv.call("reduce", in_dir="/etc", out_dir=str(srv.ROOTS[0] / "o"))["error"]
    assert "absolute path" in srv.call("readout", path="relative/dir")["error"]
    assert "outside" in srv.call("run", subject="s", out_dir=str(srv.ROOTS[0]) + "/../x")["error"]


def test_bad_names_refused(srv, tmp_path):
    assert "subject" in srv.call("run", subject="../x", out_dir=str(tmp_path / "o"))["error"]
    assert "name-value" in srv.call("run", subject="s", out_dir=str(tmp_path / "o"), options={"a;b": 1})["error"]
    assert "RHEOME_" in srv.call("run", subject="s", out_dir=str(tmp_path / "o"), env={"PATH": "/x"})["error"]


def test_run_records_provenance_in_job_and_out_dir(srv, tmp_path):
    proto = tmp_path / "proto"; proto.mkdir()
    r = srv.call("run", subject="sub01", out_dir=str(tmp_path / "out"), analyses=["resolution"],
                 proto_dir=str(proto), options={"Cohort": "norm"}, env={"RHEOME_TESTVAR": "seen"})
    assert r["function"] == "rheome.scale.run" and r["rheome"]["commit"]
    j = settle(srv, r["job"])
    assert j["status"] == "done" and j["resolved"].endswith("scale.run")
    assert any("RHEOME_TESTVAR seen" in x or "env seen" in x for x in j["log_tail"])
    p = json.loads((tmp_path / "out" / "sub01" / f"rheome_mcp_{r['job']}.json").read_text())
    assert p["args"] == ["sub01"] and p["options"]["Analyses"] == ["resolution"]
    assert p["options"]["OutDir"] == str(tmp_path / "out") and p["options"]["ProtoDir"] == str(proto)
    assert p["rheome"]["commit"] and p["server"] == srv.SERVER_VERSION and p["ok"] is True


def test_input_files_are_hashed(srv, tmp_path):
    f = tmp_path / "fit.csv"; f.write_text("a\n1\n")
    c = tmp_path / "centroids.mat"; c.write_bytes(b"x")
    r = srv.call("modes_reduce", fit_files=[str(f)], centroids=str(c), out_dir=str(tmp_path / "red"))
    p = json.loads((srv.JOBS / r["job"] / "call.json").read_text())
    assert set(p["input_files"]) == {str(f), str(c)} and len(p["input_files"][str(f)]["sha256"]) == 64
    settle(srv, r["job"])


def test_one_matlab_at_a_time(srv, tmp_path, monkeypatch):
    if os.environ.get("RHEOME_MCP_MATLAB_TEST") == "1":
        pytest.skip("timing test for the fake MATLAB")
    monkeypatch.setenv("FAKE_MATLAB_SLEEP", "1.5")
    (tmp_path / "in").mkdir()
    a = srv.call("reduce", in_dir=str(tmp_path / "in"), out_dir=str(tmp_path / "g1"))
    b = srv.call("reduce", in_dir=str(tmp_path / "in"), out_dir=str(tmp_path / "g2"))
    assert a["status"] == "running" and b["status"] == "queued"
    assert srv.call("job", job=b["job"])["status"] == "queued"
    assert [j["status"] for j in srv.call("jobs")["jobs"]] == ["queued", "running"]
    assert settle(srv, a["job"])["status"] == "done"
    assert settle(srv, b["job"])["status"] == "done"


def test_error_is_reported(srv, tmp_path):
    jid = srv.submit("fail", [], {}, {}, None)["job"]
    j = settle(srv, jid)
    assert j["status"] == "error" and j["message"] == "boom"


def test_readout_lists_and_filters(srv, tmp_path):
    d = tmp_path / "sub01"; d.mkdir()
    (d / "metrics.csv").write_text("subject,analysis,metric,band,value,unit\n"
                                   "sub01,resolution,fwhm,alpha,12,mm\nsub01,resolution,fwhm,beta,9,mm\n"
                                   "sub01,bandsnr,snr,alpha,3,dB\n")
    (d / "provenance.json").write_text('{"commit": "abc"}')
    ls = srv.call("readout", path=str(d))
    assert ls["tables"] == ["metrics"] and ls["provenance"]["commit"] == "abc"
    r = srv.call("readout", path=str(d), table="metrics", where={"analysis": "resolution"}, limit=1)
    assert r["total"] == 2 and len(r["rows"]) == 1 and r["rows"][0]["band"] == "alpha"
    r = srv.call("readout", path=str(d / "metrics.csv"), where={"band": ["alpha", "beta"], "metric": "fwhm"})
    assert r["total"] == 2
    assert "no column" in srv.call("readout", path=str(d), table="metrics", where={"nope": 1})["error"]


def test_every_call_is_logged_with_version_and_inputs(srv, tmp_path):
    srv.call("version")
    srv.call("readout", path=str(tmp_path))
    log = srv.call("calls", limit=5)["calls"]
    assert [c["tool"] for c in log] == ["readout", "version"]
    assert log[0]["inputs"] == {"path": str(tmp_path)} and log[0]["rheome"]["commit"]


def test_analyses_cached_per_version(srv):
    a = srv.call("analyses")
    assert "resolution" in a["all"] and set(a["ported"]) <= set(a["all"])
    n = len(list(srv.JOBS.iterdir()))
    assert srv.call("analyses")["all"] == a["all"] and len(list(srv.JOBS.iterdir())) == n


def test_mcp_surface(srv):
    s = srv.build_server()
    tools = asyncio.run(s.list_tools())
    names = {t.name for t in tools}
    assert names == set(srv.TOOLS)
    run = next(t for t in tools if t.name == "run")
    assert {"subject", "out_dir", "analyses"} <= set(run.input_schema["properties"])
    assert run.annotations.read_only_hint is False


@pytest.mark.skipif(os.environ.get("RHEOME_MCP_MATLAB_TEST") != "1", reason="real MATLAB: RHEOME_MCP_MATLAB_TEST=1")
def test_real_reduce_through_matlab(srv, tmp_path):
    """rheome.scale.reduce on two synthetic subjects, run by the real release through the MCP path."""
    for s, v in (("sub01", 1.0), ("sub02", 3.0)):
        d = tmp_path / "subjects" / s; d.mkdir(parents=True)
        (d / "metrics.csv").write_text("subject,cohort,dataset,analysis,metric,band,value,unit\n"
                                       f"{s},norm,t,resolution,fwhm,alpha,{v},mm\n")
        (d / "timing.csv").write_text("analysis,seconds,peak_rss_GB,status,message\nresolution,1,0.1,ok,\n")
    r = srv.call("reduce", in_dir=str(tmp_path / "subjects"), out_dir=str(tmp_path / "group"))
    j = settle(srv, r["job"], timeout=900)
    assert j["status"] == "done", j
    assert j["resolved"].startswith(str(REPO)) and "R20" in j["matlab"]
    d = srv.call("readout", path=str(tmp_path / "group"), table="distribution", where={"metric": "fwhm"})
    assert d["total"] == 1 and float(d["rows"][0]["median"]) == 2.0


def test_stdio_handshake(srv):
    """The server as a client sees it: initialize over stdio, list the tools, call one."""
    from mcp import ClientSession, StdioServerParameters
    from mcp.client.stdio import stdio_client

    async def main():
        p = StdioServerParameters(command=sys.executable, args=[str(HERE.parent / "rheome_mcp.py")], env=dict(os.environ))
        async with stdio_client(p) as (r, w), ClientSession(r, w) as s:
            init = await s.initialize()
            tools = await s.list_tools()
            v = await s.call_tool("version", {})
            return init, tools, json.loads(v.content[0].text)

    init, tools, v = asyncio.run(asyncio.wait_for(main(), 60))
    assert init.server_info.name == "rheome" and init.server_info.version == srv.SERVER_VERSION
    assert {t.name for t in tools.tools} == set(srv.TOOLS)
    assert v["rheome"]["release"] and v["server"] == srv.SERVER_VERSION
