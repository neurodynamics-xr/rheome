# FlipOut reference answers (written by the released viewer's nxr-compute WASM; read by Rheome, nxr and the Editor)

`flipout-reference.json`: geometry-central's FlipOut geodesic (Sharp & Crane 2020) as the viewer draws it,
`new Manifold(V, F).tracePath(a, b)` in nxr-compute's WASM (`tracePathSurface`: Dijkstra edge path,
`iterativeShorten()`, `getPathPolyline()` on the input mesh). Generators: `meshes.py` (the meshes and the pairs,
seeded) then `node ref.mjs` beside `nxr_compute.js`/`.wasm` of the released viewer. `tFlipOut` checks
`rheome.geom.draw.flipOut` against it.

Keys: `bumpy642`, `bumpy2562` (icospheres of 3 and 4 subdivisions, radius 0.05 m, folded by
r = 1 + a sin(6 theta) cos(5 phi); `vertices` in metres and `faces` inline, 0-based) and `fsaverage5_lh_pial`
(no mesh inline: `$NXR_FSAVERAGE5/surf/lh.pial`, FreeSurfer binary, vertices / 1000 to metres, faces as stored).
Each has `paths`: `from`, `to` (0-based), `length` (m, FlipEdgeNetwork::length), `exact` (the MMP distance, for
context), `nPoints`, `positions`, `vertices` (0-based), `weights` (each flat n x 3, row-major). A point is a vertex
(`[v v v]`, `[1 0 0]`), an edge crossing (`[a b a]`, `[1-t t 0]`, a -> b the edge's first halfedge) or a face point
(`[a b c]`, barycentric).
