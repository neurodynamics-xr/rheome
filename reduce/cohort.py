"""Cohort reducers of Rheome's per-subject tables: median and a distribution-free 95 % CI per tile or quantity.

    python3 reduce/cohort.py gauge      OUT_DIR IN_DIR [IN_DIR ...]   -> OUT_DIR/gaugecohort.csv
    python3 reduce/cohort.py multimodal OUT_DIR IN_DIR [IN_DIR ...]   -> OUT_DIR/multimodal_{relative,subjects,cohort}.csv

IN_DIR holds rheome.scale.run's per-subject folders (<sub>/gaugetensor.csv, <sub>/multimodal.csv), nsp's flat
<sub>_<table>.csv, or nsp's promoted <sub>.tar.gz (read in place, never unpacked); <sub>_metrics.csv / metrics.csv
gives the dataset, which is the group (cf-atlas labels every cohort "unlabelled").

gauge (cohort tile map, MS1 4.4): per (dataset, hemi, depth, node_id) of rheome.scale.measure_atlas "gaugetensor"
-- template tiles reached through the registered sphere, so tile k is the same in everyone -- the median over
subjects and its CI of anisotropy, normal_share and the trace-normalised tensor t11 t22 t12; axis_deg and
anisotropy_of_median from the median tensor; axis_R = |mean exp(2i axis)| over subjects. A tile whose median
colatitude is within POLAR_DEG of a pole has axis_deg and axis_R NaN (the frame turns inside it; see
rheome.scale.reducegauge, whose mean-tensor table this complements).

multimodal (MS1 5.5/5.6): per subject of rheome.scale.measure_multimodal, (a) relative band power, band / sum of
the bands, per tile, level and member, from the stored sums (exact; taken from the meg_<band>_rel rows when the table
has them); (b) the whole-cortex relative band power (level 0); (c) Spearman rho across tiles, levels >= 3, between
each relative band and each non-MEG field (PET, fibres, later BOLD) and between non-MEG fields (mean at member 0
and for MEG, ms -- the coefficient energy -- for the signed members j > 0). Then per quantity the cohort median and
CI. ⚠ Each person's tiles are their OWN rheome.geom.tree, so tile codes are not the same place across people: the
cohort step reduces per-subject quantities, never tiles.

CI: the order-statistic interval for the median, [x_(k), x_(n-k+1)], k = binom.ppf(0.025, n, 0.5): coverage >= 95 %
for any continuous distribution; NaN when n < 6. Descriptive only: no test across tiles, cohorts or covariates.

Author: Diellor Basha, 2026
"""
import csv
import io
import os
import sys
import tarfile
from collections import defaultdict

import numpy as np
from scipy.stats import binom, spearmanr

POLAR_DEG = 35.0
GAUGE_METRICS = ["anisotropy", "normal_share", "t11", "t22", "t12"]


def median_ci(x, alpha=0.05):
    """(n, median, lo, hi) of the finite values of x; the order-statistic CI of the median."""
    x = np.sort(np.asarray([v for v in x if np.isfinite(v)], float))
    n = x.size
    if n == 0:
        return 0, np.nan, np.nan, np.nan
    k = int(binom.ppf(alpha / 2, n, 0.5))
    lo, hi = (x[k - 1], x[n - k]) if k >= 1 else (np.nan, np.nan)
    return n, float(np.median(x)), float(lo), float(hi)


def _rows(text):
    return list(csv.DictReader(io.StringIO(text)))


def subjects(dirs, table):
    """Yield (subject, dataset, rows of <table>.csv) for every subject under dirs (folders, flat files, tar.gz)."""
    for root in dirs:
        for dp, _, files in os.walk(root, followlinks=True):
            for f in sorted(files):
                p = os.path.join(dp, f)
                if f.startswith("._"):
                    continue
                if f.endswith(".tar.gz"):
                    with tarfile.open(p) as tf:
                        got = {}
                        for m in tf.getmembers():
                            b = os.path.basename(m.name)
                            if m.isfile() and (b.endswith("_%s.csv" % table) or b.endswith("_metrics.csv")):
                                got[b] = tf.extractfile(m).read().decode()
                    for b, txt in got.items():
                        if b.endswith("_%s.csv" % table):
                            sub = b[: -len("_%s.csv" % table)]
                            yield sub, _dataset(got.get(sub + "_metrics.csv")), _rows(txt)
                elif f == table + ".csv" or f.endswith("_%s.csv" % table):
                    sub, pre = (os.path.basename(dp), "") if f == table + ".csv" else (f[: -len("_%s.csv" % table)], f[: -len("%s.csv" % table)])
                    m = os.path.join(dp, pre + "metrics.csv")
                    yield sub, _dataset(open(m).read() if os.path.isfile(m) else None), _rows(open(p).read())


def _dataset(metrics_text):
    if not metrics_text:
        return ""
    r = _rows(metrics_text)
    return r[0].get("dataset", "") if r else ""


def _write(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def _axis(t11, t22, t12):
    return 0.5 * np.degrees(np.arctan2(2 * t12, t11 - t22))


# ---------------------------------------------------------------- gauge: the cohort tile map
def gauge(out_dir, dirs, polar_deg=POLAR_DEG):
    T = defaultdict(lambda: defaultdict(list))  # (dataset, hemi, depth, node) -> column -> values over subjects
    subs = set()
    for sub, ds, rows in subjects(dirs, "gaugetensor"):
        subs.add((ds, sub))
        for r in rows:
            c = T[(ds, r["hemi"], int(r["depth"]), int(r["node_id"]))]
            for k in GAUGE_METRICS + ["axis_deg", "colat_median_deg"]:
                c[k].append(float(r[k]) if r[k] not in ("", "NaN", "nan") else np.nan)
    if not T:
        raise SystemExit("no gaugetensor.csv under %s" % " ".join(dirs))
    out = []
    for key in sorted(T):
        c = T[key]
        row = dict(zip(["dataset", "hemi", "depth", "node_id"], key))
        row["n"] = len(c["t11"])
        med = {}
        for k in GAUGE_METRICS:
            _, med[k], lo, hi = median_ci(c[k])
            row.update({k + "_median": med[k], k + "_ci_lo": lo, k + "_ci_hi": hi})
        d = np.hypot(med["t11"] - med["t22"], 2 * med["t12"])
        row["anisotropy_of_median"] = d / (med["t11"] + med["t22"])
        a = np.radians([v for v in c["axis_deg"] if np.isfinite(v)])
        row["axis_deg"] = _axis(med["t11"], med["t22"], med["t12"])
        row["axis_R"] = float(abs(np.mean(np.exp(2j * a)))) if a.size else np.nan
        row["colat_median_deg"] = float(np.nanmedian(c["colat_median_deg"]))
        row["polar"] = bool(row["colat_median_deg"] < polar_deg or row["colat_median_deg"] > 180 - polar_deg)
        if row["polar"]:
            row["axis_deg"] = row["axis_R"] = np.nan
        out.append(row)
    os.makedirs(out_dir, exist_ok=True)
    _write(os.path.join(out_dir, "gaugecohort.csv"), out)
    per = defaultdict(int)
    for ds, _ in subs:
        per[ds] += 1
    print("[cohort gauge] %s; %d tile rows (%d polar)" % (", ".join("%s %d" % kv for kv in sorted(per.items())), len(out), sum(r["polar"] for r in out)))
    return out


# ---------------------------------------------------------------- multimodal
def relative(rows):
    """{(level, code, member): {band: rel}} and scale_mm per member, from meg_<band>_rel rows or the stored sums."""
    rel, scale = defaultdict(dict), {}
    have = any(r["field"].endswith("_rel") for r in rows)
    raw = defaultdict(dict)
    for r in rows:
        key = (int(r["level"]), int(r["code"]), int(r["member"]))
        scale[key[2]] = float(r["scale_mm"])
        f = r["field"]
        if r["modality"] != "MEG":
            continue
        if have and f.endswith("_rel"):
            rel[key][f[4:-4]] = float(r["mean"])
        elif not have:
            raw[key][f[4:]] = float(r["sum_ax"])
    for key, b in raw.items():
        tot = sum(b.values())
        rel[key] = {k: v / tot for k, v in b.items()}
    return rel, scale


def multimodal(out_dir, dirs, min_level=3):
    R, S = [], []  # per-subject relative tiles; per-subject quantities
    for sub, ds, rows in subjects(dirs, "multimodal"):
        rel, scale = relative(rows)
        other = defaultdict(dict)  # field -> (level, code, member) -> value
        for r in rows:
            if r["modality"] != "MEG":
                m = int(r["member"])
                other[r["field"]][(int(r["level"]), int(r["code"]), m)] = float(r["mean"] if m == 0 else r["ms"])
        for (lv, code, m), b in sorted(rel.items()):
            for band, v in b.items():
                R.append(dict(dataset=ds, subject=sub, level=lv, code=code, member=m, scale_mm=scale[m], band=band, rel=v))
                if lv == 0:
                    S.append(dict(dataset=ds, subject=sub, quantity="relative_power_cortex", level=0, member=m,
                                  scale_mm=scale[m], field_a="meg_" + band + "_rel", field_b="", value=v, n_tiles=1))
        levels = sorted({k[0] for k in rel if k[0] >= min_level})
        bands = sorted({b for v in rel.values() for b in v})
        names = sorted(other)
        for lv in levels:
            for m in sorted(scale):
                codes = sorted(k[1] for k in rel if k[0] == lv and k[2] == m)
                vec = {"meg_" + b + "_rel": [rel[(lv, c, m)][b] for c in codes] for b in bands}
                vec.update({f: [other[f][(lv, c, m)] for c in codes] for f in names})
                pairs = [(a, f) for a in ("meg_" + b + "_rel" for b in bands) for f in names] + \
                        [(names[i], names[j]) for i in range(len(names)) for j in range(i + 1, len(names))]
                for a, f in pairs:
                    rho = spearmanr(vec[a], vec[f])[0] if len(codes) > 2 else np.nan
                    S.append(dict(dataset=ds, subject=sub, quantity="rho_tiles", level=lv, member=m, scale_mm=scale[m],
                                  field_a=a, field_b=f, value=float(rho), n_tiles=len(codes)))
    if not S:
        raise SystemExit("no multimodal.csv under %s" % " ".join(dirs))
    G = defaultdict(list)
    for s in S:
        G[(s["dataset"], s["quantity"], s["level"], s["member"], s["scale_mm"], s["field_a"], s["field_b"])].append(s["value"])
    C = []
    for key in sorted(G):
        n, med, lo, hi = median_ci(G[key])
        C.append(dict(zip(["dataset", "quantity", "level", "member", "scale_mm", "field_a", "field_b"], key),
                      n=n, median=med, ci_lo=lo, ci_hi=hi))
    os.makedirs(out_dir, exist_ok=True)
    _write(os.path.join(out_dir, "multimodal_relative.csv"), R)
    _write(os.path.join(out_dir, "multimodal_subjects.csv"), S)
    _write(os.path.join(out_dir, "multimodal_cohort.csv"), C)
    print("[cohort multimodal] %d subjects, %d per-subject quantities, %d cohort rows"
          % (len({(s["dataset"], s["subject"]) for s in S}), len(S), len(C)))
    return R, S, C


if __name__ == "__main__":
    if len(sys.argv) < 4 or sys.argv[1] not in ("gauge", "multimodal"):
        raise SystemExit(__doc__.split("\n\n")[1])
    {"gauge": gauge, "multimodal": multimodal}[sys.argv[1]](sys.argv[2], sys.argv[3:])
