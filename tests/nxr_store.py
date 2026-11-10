"""A synthetic schema-50 nxr dataset folder for tNxrReaders: dataset.sqlite from nxr's own model.sql (foreign keys off: only the rows the readers
follow are written), a subject and a template each with a Hodge-Dirac axis, one acquisition, and the subject's stored registration
hodgeDirac_to_<template> (T float32, blosc-compressed when asked), as nxr's store_mode_registration writes it. Nothing here is real data.

    uv run --frozen --project <nxr checkout> python -I nxr_store.py <out folder> [compress]
"""
import json
import sqlite3
import sys
from pathlib import Path

import numpy as np
import zarr
from nxr.io.schema import model_sql

out = Path(sys.argv[1])
compress = len(sys.argv) > 2
out.mkdir(parents=True, exist_ok=True)
cx = sqlite3.connect(out / "dataset.sqlite")
cx.executescript(model_sql())
cx.execute("PRAGMA foreign_keys = OFF")
for sid, name in (("s1", "subj"), ("s2", "@tmpl")):
    cx.execute("INSERT INTO subject (id, dataset_id, name, path, status) VALUES (?, 'd', ?, ?, 'complete')", (sid, name, f"{name.strip('@')}.nxr.zarr"))
    cx.execute("INSERT INTO manifold (id, subject_id, name, type, n_vertices, topology_id, status) VALUES (?, ?, 'hodgeDirac', 'points', 6, 't', 'complete')", ("ax" + sid, sid))
cx.execute("INSERT INTO acquisition (id, subject_id, label, ses, task, run) VALUES ('a1', 's1', 'ses-02_task-rest_run-02', '02', 'rest', '02')")
T = np.arange(12, dtype=np.float32).reshape(3, 4) / 7
rec = {"rows": "template modes", "cols": "subject modes", "sampling": "shepard", "neighbours": 6, "power": 3.0, "flow": "fields",
       "from": {"package": "aa", "K": 4, "mu": [1.0, 2.0, 3.0, 4.0], "mesh": "m1"}, "to": {"package": "bb", "K": 3, "mu": [1.5, 2.5, 3.5], "kind": [0, 2, 0], "mesh": "m2"}}
cx.execute("INSERT INTO operator (id, subject_id, name, path, kind, from_manifold_id, to_manifold_id, layout, n_rows, n_cols, params_json, status, acquisition_id) "
           "VALUES ('o1', 's1', 'hodgeDirac_to_tmpl', 'hodgeDirac_to_tmpl', 'registration', 'axs1', 'axs2', 'dense', 3, 4, ?, 'complete', 'a1')", (json.dumps(rec),))
cx.commit()
cx.close()
z = out / "subj.nxr.zarr" / "hodgeDirac_to_tmpl"
a = zarr.create_array(store=str(z), shape=T.shape, dtype="float32", chunks=(2, 2), zarr_format=3, fill_value=0, compressors="auto" if compress else None)
a[...] = T
