"""Golden cases for rheome.geom private/pysort: Python's sorted() with nsp's atlas-ladder comparator
(key[x] - key[y], NaN read as 0) on keys with NaN / +-Inf / ties. Run with CPython 3.12 (the
reference's), then in MATLAB: pysortWriteGolden('pysortGolden.txt') -> pysortGolden.mat."""
import math, random, sys
from functools import cmp_to_key

assert sys.version_info[:2] == (3, 12), sys.version
rng = random.Random(20261014)
out = []
for case in range(300):
    n = rng.randint(64, 600) if case % 3 else rng.randint(2, 63)
    pool = [rng.random() for _ in range(max(1, n // rng.randint(1, 8)))]
    pn, pi = rng.choice([0.05, 0.2, 0.5]), rng.choice([0.0, 0.05, 0.2])
    key = []
    for _ in range(n):
        u = rng.random()
        key.append(math.nan if u < pn else (rng.choice([math.inf, -math.inf]) if u < pn + pi else rng.choice(pool)))
    def cmp(x, y):
        c = key[x] - key[y]
        c = 0.0 if math.isnan(c) else c
        return -1 if c < 0 else (1 if c > 0 else 0)
    perm = sorted(range(n), key=cmp_to_key(cmp))
    out.append(" ".join([str(n)] + [repr(k) for k in key] + [str(p + 1) for p in perm]))
open(sys.argv[1], "w").write("\n".join(out) + "\n")
