"""nxr-convert's own chain on one subject of an nxr store, for tNxrStoreParity: modes_input (its lookups and reads, unchanged) ->
glue_mesh -> prepared -> solve_hodge. The store is opened READ-ONLY (a sqlite URI, no writer lock); outputs are little-endian f64 files.

    uv run --frozen --project <nxr-convert> python -I nxr_convert_chain.py <store> <subject> <n_modes> <out>
"""
import sqlite3
import sys
from pathlib import Path

import numpy as np
from nxr_convert.eigen.glue import glue_mesh, mesh_key
from nxr_convert.eigen.hodge_dirac import modes_input, solve_hodge
from nxr_convert.eigen.intrinsic import prepared


class _Db:
    def __init__(self, path: Path):
        self.cx = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
        self.cx.row_factory = sqlite3.Row

    def one(self, q, *args):
        r = self.cx.execute(q, args).fetchone()
        return dict(r) if r else None


class _Subject:
    """what modes_input reads of a subject: db.one, id, row, at"""
    def __init__(self, store: Path, name: str):
        self.db = _Db(store / "dataset.sqlite")
        self.row = self.db.one("SELECT * FROM subject WHERE name = ? AND status = 'complete'", name)
        if self.row is None:
            raise KeyError(name)
        self.id, self.root = self.row["id"], store / self.row["path"]

    def at(self, path: str) -> Path:
        return self.root / path


store, name, n, out = Path(sys.argv[1]), sys.argv[2], int(sys.argv[3]), Path(sys.argv[4])
inp = modes_input(_Subject(store, name))
V, G, keep, _ = glue_mesh(inp["P"], inp["F"], inp["wall"])
m = prepared(V, G)
sol = solve_hodge(m, n)
out.mkdir(parents=True, exist_ok=True)
for k, a in {"P": inp["P"], "F": inp["F"], "wall": inp["wall"], "V": V, "G": G, "keep": keep, "iF": m.F,
             "el": sol["exact_lambda"], "cl": sol["coexact_lambda"], "em": sol["exact_modes"], "cm": sol["coexact_modes"]}.items():
    np.ascontiguousarray(a, "<f8").tofile(out / f"{k}.f64")
(out / "key.txt").write_text(mesh_key(G))
