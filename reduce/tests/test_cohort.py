"""reduce/cohort.py: CI coverage, the gauge tile map from flat and tar.gz inputs, relative power from sums or rows.
python3 -m unittest discover -s reduce/tests"""
import csv
import io
import os
import sys
import tarfile
import tempfile
import unittest

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import cohort  # noqa: E402

GH = "hemi,depth,node_id,n_vertices,t11,t22,t12,anisotropy,axis_deg,normal_share,trace,colat_median_deg\n"


def _gt(t12, colat):
    a = np.hypot(0.2, 2 * t12) / 1.0
    return GH + "L,3,8,100,0.6,0.4,%g,%g,%g,0.2,1,%g\nL,3,9,100,0.5,0.5,0,0,0,0.3,1,10\n" % (t12, a, cohort._axis(.6, .4, t12), colat)


class Cohort(unittest.TestCase):
    def test_ci_covers_the_median(self):
        rng = np.random.default_rng(0)
        hit = [lo <= 0 <= hi for lo, hi in (cohort.median_ci(rng.standard_cauchy(30))[2:] for _ in range(2000))]
        self.assertGreaterEqual(np.mean(hit), 0.94)
        self.assertTrue(np.isnan(cohort.median_ci([1, 2, 3])[2]))   # n < 6: no interval
        self.assertEqual(cohort.median_ci([1, 2, 3, np.nan])[:2], (3, 2.0))

    def test_gauge_flat_and_tar(self):
        with tempfile.TemporaryDirectory() as d:
            ind = os.path.join(d, "in")
            os.makedirs(ind)
            for i in range(7):            # 4 flat files, 3 tarballs, dataset from metrics
                sub, txt = "sub-%02d" % i, _gt(0.01 * i, 100)
                met = "subject,cohort,dataset,analysis,metric,band,value,unit\nx,unlabelled,omega,a,m,L,1,u\n"
                if i < 4:
                    open(os.path.join(ind, sub + "_gaugetensor.csv"), "w").write(txt)
                    open(os.path.join(ind, sub + "_metrics.csv"), "w").write(met)
                else:
                    with tarfile.open(os.path.join(ind, sub + ".tar.gz"), "w:gz") as tf:
                        for n, t in ((sub + "_gaugetensor.csv", txt), ("./" + sub + "_metrics.csv", met)):
                            b = t.encode(); ti = tarfile.TarInfo(n); ti.size = len(b); tf.addfile(ti, io.BytesIO(b))
            out = cohort.gauge(os.path.join(d, "out"), [ind])
            r8 = [r for r in out if r["node_id"] == 8][0]
            self.assertEqual((r8["dataset"], r8["n"]), ("omega", 7))
            self.assertAlmostEqual(r8["t12_median"], 0.03)
            self.assertAlmostEqual(r8["t12_ci_lo"], 0.0)       # n = 7: k = 1, the extremes
            self.assertAlmostEqual(r8["t12_ci_hi"], 0.06)
            self.assertAlmostEqual(r8["axis_deg"], cohort._axis(.6, .4, .03))
            self.assertFalse(r8["polar"])
            r9 = [r for r in out if r["node_id"] == 9][0]
            self.assertTrue(r9["polar"] and np.isnan(r9["axis_deg"]) and np.isnan(r9["axis_R"]))
            self.assertTrue(os.path.isfile(os.path.join(d, "out", "gaugecohort.csv")))

    def test_relative_from_sums_equals_rows(self):
        hdr = "level,code,hemisphere,field,modality,unit,member,scale_mm,area_m2,sum_ax,sum_ax2,mean,ms\n"
        sums = {"delta": 1.0, "theta": 2.0, "alpha": 3.0, "beta": 4.0}
        a = hdr + "".join("0,1,0,meg_%s,MEG,u,0,0,1,%g,1,%g,1\n" % (b, v, v) for b, v in sums.items())
        b = a + "".join("0,1,0,meg_%s_rel,MEG,fraction,0,0,1,NaN,NaN,%g,NaN\n" % (k, v / 10) for k, v in sums.items())
        for t in (a, b):
            rel, _ = cohort.relative(list(csv.DictReader(io.StringIO(t))))
            self.assertAlmostEqual(rel[(0, 1, 0)]["alpha"], 0.3)
            self.assertAlmostEqual(sum(rel[(0, 1, 0)].values()), 1.0)

    def test_multimodal_pools_members_across_own_scales(self):
        hdr = "level,code,hemisphere,field,modality,unit,member,scale_mm,area_m2,sum_ax,sum_ax2,mean,ms\n"
        with tempfile.TemporaryDirectory() as d:
            ind = os.path.join(d, "in")
            os.makedirs(ind)
            for i in range(7):            # each subject's member 1 has its own scale_mm, as with own trees
                rows = "".join("%d,%d,0,meg_%s,MEG,u,%d,%g,1,%g,1,1,1\n" % (lv, c, b, m, 0 if m == 0 else 8 + i, v + c * (b == "alpha"))
                               for lv, cs in ((0, [1]), (3, [1, 2, 3, 4])) for c in cs for m in (0, 1)
                               for b, v in (("alpha", 1.0 + i), ("beta", 2.0)))
                rows += "".join("3,%d,0,pet_x,PET,u,%d,%g,1,1,1,%g,%g\n" % (c, m, 0 if m == 0 else 8 + i, c, c)
                                for c in (1, 2, 3, 4) for m in (0, 1))
                open(os.path.join(ind, "sub-%02d_multimodal.csv" % i), "w").write(hdr + rows)
            _, _, C = cohort.multimodal(os.path.join(d, "out"), [ind])
            m1 = [r for r in C if r["quantity"] == "relative_power_cortex" and r["member"] == 1]
            self.assertEqual(len(m1), 2)                       # alpha, beta: one row each, not one per subject
            self.assertEqual([(r["quantity"], r["member"], r["n"]) for r in C if r["n"] != 7], [])
            self.assertAlmostEqual(m1[0]["scale_mm"], 11.0)    # median of 8..14


if __name__ == "__main__":
    unittest.main()
