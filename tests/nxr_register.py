"""nxr's sphere registration on a synthetic sphere, for tRegister: the reference rheome.register is held to.

    uv run --frozen --project <nxr checkout> python -I nxr_register.py <folder>

Reads X.f64 [n x 3] (source sphere), F.f64 [nF x 3] (0-based faces), Y.f64 [m x 3] (points), writes for each of the three samplings
R_<sampling>.f64 (dense [m x n], row-major), and for shepard its nxr.register.sphere.shepard_weights as S_idx.f64 / S_w.f64 (0-based).
"""
import sys
from pathlib import Path

import numpy as np
from nxr.register.sphere import registration, shepard_weights

d = Path(sys.argv[1])
rd = lambda n, w: np.fromfile(d / f"{n}.f64", "<f8").reshape(-1, w)
X, F, Y = rd("X", 3), rd("F", 3).astype(int), rd("Y", 3)
hs, hd = np.full(len(X), "L"), np.full(len(Y), "L")
for s in ("barycentric", "nearest", "shepard"):
    R = registration(X, F, hs, Y, hd, sampling=s, neighbours=6, power=3)
    np.ascontiguousarray(R.toarray(), "<f8").tofile(d / f"R_{s}.f64")
i, w = shepard_weights(X / np.linalg.norm(X, axis=1, keepdims=True), Y / np.linalg.norm(Y, axis=1, keepdims=True), 6, 3.0)
np.ascontiguousarray(i, "<f8").tofile(d / "S_idx.f64")
np.ascontiguousarray(w, "<f8").tofile(d / "S_w.f64")
