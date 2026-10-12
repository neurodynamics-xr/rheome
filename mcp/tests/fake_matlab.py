#!/usr/bin/env python3
"""Stands in for `matlab -batch "...; rheome_mcp_call('<jobDir>')"` in the unit tests: no MATLAB needed.
FAKE_MATLAB_SLEEP seconds of work; a function named "fail" errors; scale.analyses returns a list."""
import json, os, re, sys, time
job = re.search(r"rheome_mcp_call\('(.*)'\)", sys.argv[-1]).group(1).replace("''", "'")
c = json.load(open(os.path.join(job, "call.json")))
time.sleep(float(os.environ.get("FAKE_MATLAB_SLEEP", "0")))
ok = c["function"] != "fail"
print("fake matlab ran", c["function"], "env", os.environ.get("RHEOME_TESTVAR", ""))
r = {"ok": ok, "message": "" if ok else "boom", "matlab": "fake", "resolved": c["root"] + "/+rheome/" + c["function"],
     "seconds": 0.0, "value": ["resolution", "bandsnr"] if c["returns"] else None}
json.dump(r, open(os.path.join(job, "result.json"), "w"))
