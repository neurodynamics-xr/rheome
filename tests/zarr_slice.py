"""zarr-python's read of a part of an array, for tNxrStoreParity: the reference rheome.io.read.zarr is held to on a real store.

    uv run --frozen --project <nxr-convert> python -I zarr_slice.py <array folder> <out .f64> [first last per dimension, 1-based inclusive]
"""
import sys

import numpy as np
import zarr

a = zarr.open_array(sys.argv[1], mode="r")
r = [int(v) for v in sys.argv[3:]]
sl = tuple(slice(r[2 * d] - 1, r[2 * d + 1]) if 2 * d + 1 < len(r) else slice(None) for d in range(a.ndim))
np.ascontiguousarray(a[sl], "<f8").tofile(sys.argv[2])
