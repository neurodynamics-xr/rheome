"""tSurfaceLadder golden: nsp's atlas ladder (nxr_convert/atlas/ladder.py) on ladderWriteMeshes' meshes.
python3.12 -I ladderCaptureGolden.py <ladder.py> <meshes.mat> <out.mat>  (numpy + scipy)"""
import importlib.util, sys
import numpy as np, scipy.io as sio
spec = importlib.util.spec_from_file_location("ladder", sys.argv[1]); L = importlib.util.module_from_spec(spec)
sys.modules["ladder"] = L; spec.loader.exec_module(L)
m = sio.loadmat(sys.argv[2], squeeze_me=True)["meshes"]
labels = np.empty(len(m), dtype=object)
for i, s in enumerate(m):
    lad = L.surface_ladder(s["V"].item(), s["F"].item().astype(np.int64) - 1)
    labels[i] = np.stack(lad.labels, 1).astype(np.double)
    print(i, labels[i].shape)
sio.savemat(sys.argv[3], {"meshes": m, "labels": labels, "python": sys.version.split()[0], "numpy": np.__version__}, do_compression=True)
