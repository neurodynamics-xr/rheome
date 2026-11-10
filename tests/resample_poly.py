"""scipy's resample_poly on rows, for tResample: the reference rheome.filters.resample is held to (nxr-convert's call: the ratio of
the two rates in mHz, reduced; default window and padding).

    uv run --frozen --project <nxr-convert> python -I resample_poly.py <rows.f64 [C x T] row-major> <C> <fs> <to> <out.f64>
"""
import sys
from fractions import Fraction

import numpy as np
from scipy.signal import resample_poly

src, C, fs, to, out = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), sys.argv[5]
x = np.fromfile(src, "<f8").reshape(C, -1)
r = Fraction(round(to * 1000), round(fs * 1000))
np.ascontiguousarray(resample_poly(x, r.numerator, r.denominator, axis=1), "<f8").tofile(out)
