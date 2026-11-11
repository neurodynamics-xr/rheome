% GEOM  Mesh generation, surface geometry, and the surface's bookkeeping.
%
% Reference geometries for validation (built, not read from disk), the tile and ladder bookkeeping of a surface, and
% intrinsic-geometry primitives ported from geometry-central (via nxr-compute).
%
% BY JOB (2026-11-10): a function that MEASURES, SELECTS or DRAWS lives in a sub-package named for its job, with the
% same names as nxr and the Cortical Flow Editor; mesh structure and bookkeeping stay here in rheome.geom.
%   rheome.geom.measure  returns numbers               (distances, tile means, the intrinsic mean)
%   rheome.geom.select   returns a set of vertices     (the cortex)
%   rheome.geom.draw     returns a curve on the surface (the FlipOut geodesic)
% The old names (geodesic, heatdistance, tilemean, karcher_mean, cortexmask, flipOut, geodesicpath) forward and warn once
% (rheome:deprecated).
%
% Sub-packages:
%   rheome.geom.measure.heatDistance      - heat-method geodesic distance on a surface (+ reusable solver)
%   rheome.geom.measure.heatDistanceEigen - the heat method in the LBO eigenbasis (no solve)
%   rheome.geom.measure.karcherMean       - intrinsic (geodesic Fréchet) mean of points on the surface
%   rheome.geom.measure.tileMean          - pool a vertex map over tiles (area-weighted mean, or median). Linear,
%                        so it rolls up the tree. ⚠ a tile mean keeps a band-pass peak only up to the peak's core
%                        radius: sigma-25 mm tracked on 23-33 mm tiles, lost at 47 mm and coarser.
%   rheome.geom.select.cortex             - which vertices are cortex (the medial wall excluded)
%   rheome.geom.draw.flipOut              - the FlipOut geodesic between two vertices (Sharp & Crane 2020): a port
%                        of geometry-central's, the line the viewer, nxr and the Editor draw (to 1e-9 m)
%   rheome.geom.draw.line                 - a geodesic between two vertices on the faces (calls draw.flipOut)
%   rheome.geom.draw.flipMesh             - the mesh FlipOut works on, built once (pass it to line / flipOut)
%
% Vectors on a fold:
%   rheome.geom.GeometryKit   - the vector heat method's transport and the connection Laplacian's gauge: lift / embed tangent vectors, carry them to an anchor
%                        (nxr GeometryKit). Vectors from different places are carried before they are added.
%
% Functions:
%   rheome.geom.icosphere     - geodesic sphere by icosahedron subdivision (FreeSurfer-style ico)
%   rheome.geom.edgegraph     - the mesh as a Dijkstra graph, edges + intrinsic diagonals: +1.4% against the
%                        analytic great circle at any range, symmetric. The ruler for a speed.
%   rheome.geom.glue          - the two hemispheres made one closed surface by anatomy (medial walls cut, a bridge
%                        grown): the editor's glueMesh, bit-identical to nxr.geom.glue
%   rheome.geom.meshkey       - the editor's mesh fingerprint '<nF>:<fnv1a>' a stored basis or kernel records
%
% THE SPATIAL PYRAMID (2026-09-25):
%   rheome.geom.tree   - recursive spectral bisection of a surface: a NESTED, DISJOINT partition,
%                 the cortical twin of rheome.sensors.tree. Sums over its vertices roll up the tree
%                 exactly as sums over samples roll up the dyadic time grid. ⚠ sparse eigs,
%                 not dense eig: 20484 vertices is not a dense eigenproblem. Measured on the
%                 reference cortex: 127 nodes to depth 6 in 1 s, each level halving area exactly
%                 and dividing the diameter by sqrt(2), leaves partitioning all 20484.
%                 ⭐ Method="Ladder" / Codes= give the ATLAS's tiles instead (VR1 §2.9: Fiedler bisection
%                 is a different, unequal tree, ARI 0.36-0.40, area CV 25 % at depth 8).
%   rheome.geom.surfaceladder - the atlas's dyadic ladder (equal-mass cuts across the geodesic diameter),
%                 bit-exact to nsp's atlas ladder.py: identical on 10 VR1 hemispheres at all 12 levels.
%   rheome.geom.sphereframe - the group gauge: the registered sphere's meridian pushed to the cortex
%   rheome.geom.spherepatches - ico patches on the registered sphere (same place in every subject), and
%                               where the polar gauge turns too much to pool (.excluded -> use chart x)
%   rheome.geom.spiralorder   - a 1-D order of vertices or tiles along a pole-to-pole spiral in the same gauge:
%                               the same order in every subject (sphere.reg); meridional motion reads as
%                               monotone motion along it
%                 (poles at sphere.reg's +-z, the same anatomical place in every subject).
%   rheome.geom.tiles  - one depth of the tree as a graph: vertex -> tile lookup, membership, and the
%                 tile adjacency P'*A*P (cut-edge counts). .heap says whether node k's children
%                 are 2k and 2k+1 -- true only if no branch stopped early.
%
%   rheome.geom.measure.tileMean - pools a vertex map over these tiles (above).
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
