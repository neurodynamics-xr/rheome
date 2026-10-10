% GEOM  Mesh generation and geometry utilities.
%
% Reference geometries for validation (built, not read from disk) and intrinsic-geometry
% primitives ported from geometry-central (via nxr-compute) for measuring/averaging on the sheet.
%
% Functions:
%   rheome.geom.icosphere     - geodesic sphere by icosahedron subdivision (FreeSurfer-style ico)
%   rheome.geom.geodesic      - heat-method geodesic distance on a surface (+ reusable solver)
%                        ⚠ short range only: asymmetric on the reference cortex (87.6 vs 43.7 mm on one
%                        pair) and saturating at ~49 edges on a sphere. Read speeds with edgegraph.
%   rheome.geom.edgegraph     - the mesh as a Dijkstra graph, edges + intrinsic diagonals: +1.4% against the
%                        analytic great circle at any range, symmetric. The ruler for a speed.
%   rheome.geom.karcher_mean  - intrinsic (geodesic Fréchet) mean of points on the surface
%
% THE SPATIAL PYRAMID (2026-09-25):
%   rheome.geom.tree   - recursive spectral bisection of a surface: a NESTED, DISJOINT partition,
%                 the cortical twin of rheome.sensors.tree. Sums over its vertices roll up the tree
%                 exactly as sums over samples roll up the dyadic time grid. ⚠ sparse eigs,
%                 not dense eig: 20484 vertices is not a dense eigenproblem. Measured on the
%                 reference cortex: 127 nodes to depth 6 in 1 s, each level halving area exactly
%                 and dividing the diameter by sqrt(2), leaves partitioning all 20484.
%   rheome.geom.sphereframe - the group gauge: the registered sphere's meridian pushed to the cortex
%   rheome.geom.spherepatches - ico patches on the registered sphere (same place in every subject), and
%                               where the polar gauge turns too much to pool (.excluded -> use chart x)
%                 (poles at sphere.reg's +-z, the same anatomical place in every subject).
%   rheome.geom.tiles  - one depth of the tree as a graph: vertex -> tile lookup, membership, and the
%                 tile adjacency P'*A*P (cut-edge counts). .heap says whether node k's children
%                 are 2k and 2k+1 -- true only if no branch stopped early.
%
%   rheome.geom.tilemean - pool a vertex map over tiles (area-weighted mean, or median). Linear, so it
%                 rolls up the tree. ⚠ a tile mean keeps a band-pass peak only up to the peak's
%                 core radius: sigma-25 mm tracked on 23-33 mm tiles, lost at 47 mm and coarser.
%
% THE JOINT BOOKKEEPING LADDER (2026-09-29):
%   rheome.geom.jointladder - one row per observation window of the eigenmode-frequency plane, written in
%                 space and time: tile size (depth) x octave window x frame rate, with the feature
%                 sizes and speeds (m/s = the rate) each window can hold. A row is a SELECTION; what
%                 is measured in it is chosen per question and stored at its address.
%   rheome.geom.jointcell   - pick a row: coarsest tile holding the feature, coarsest rate following its
%                 speed (envelope) or the carrier (phase); refuses what no row can hold.
%
% See also: rheome.detect.blobscale (the smooth, overlapping half: measurement rather than
% bookkeeping), rheome.operators.laplace_beltrami, rheome.sensors.tree
%
% Author: Diellor Basha, 2026
