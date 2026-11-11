% OPERATORS  Discrete differential operators on a triangle surface mesh.
%
% Pure-MATLAB reproductions of the surface operators Brainstorm computes, built from
% just the mesh (vertices + faces).
%
% Functions (per-hemisphere, mesh-only):
%   rheome.operators.mass              - mass matrix (mixedVoronoi default / galerkin / barycentric)
%   rheome.operators.laplace_beltrami  - cotangent Laplace-Beltrami stiffness (+ mass); intrinsic Delaunay mesh + mixed-Voronoi mass by default (mesh=, mass=, mollify=, flipEps=)
%   rheome.operators.connection_laplacian - Levi-Civita connection Laplacian (complex tangent); intrinsic by default, mesh='embedded' option
%   rheome.operators.gauge             - a per-vertex tangent FRAME: a reference direction field, by
%                                 vector diffusion, the smoothest eigenvector (Knoppel 2013) or a
%                                 trivial connection (parallel, singular only at chosen poles).
%                                 ⚠ connection_laplacian's own .e1/.e2 are NOT a usable gauge --
%                                 79 deg between neighbours -- so anything reading a per-component
%                                 value needs this instead
%   rheome.operators.trivial_connection - Crane 2010 trivial connection on the dual complex: indices k_f
%                                 summing to chi, face curvature WRAPPED (the chart parts removed)
%   rheome.operators.face_gradient     - per-face constant-gradient primitives (div/curl)
%   rheome.operators.dirac_frame       - relative Dirac quaternion operator (default: tau blend)
%   rheome.operators.dirac_intrinsic   - Crane's intrinsic quaternionic Dirac itself: vertex quaternions -> faces
%   rheome.operators.dirac_intrinsic_sq - intrinsic (edge-based) quaternionic Dirac squared, D' MF D
%   rheome.operators.dirac_extrinsic   - relative (extrinsic) Dirac energy block E = D' *_F D
%   rheome.operators.spin_dirac        - spin Diracs on quaternion spinors per face (extrinsic, Hoffmann-Ye, Ye)
%   rheome.operators.vertex_to_face    - vertex spinors to faces, A (4nF x 4nV)
%   rheome.operators.intrinsic_delaunay - intrinsic Delaunay triangulation (flips, mollified lengths, mixed-Voronoi mass, signposts) and its cotangent Laplacian
%   rheome.operators.hodge_dirac       - Hodge-Dirac d + delta on 0/1/2-forms (K, M, d0, d1, stars) over the intrinsic Delaunay triangulation
%   rheome.operators.connection_laplacian_idt - connection (vector) Laplacian on the intrinsic Delaunay triangulation
%   rheome.operators.weak_differential - weak-form curl and divergence as [nV x 3nV] operators
%   rheome.operators.harmonicfill      - a field on the medial wall replaced by its harmonic extension from the rim
%
% WHOLE-BRAIN connectome operators (from tractography fibers; the connectome couples the
% two hemispheres, so these are single-block whole-brain, not per-hemisphere):
%   rheome.operators.connectome            - vertex connectome W from fiber endpoints (+smoothing, CC)
%   rheome.operators.connectome_laplacian  - symmetric normalized graph Laplacian  I - D^-1/2 W D^-1/2
%   rheome.operators.lb_connectome         - LBO bridged by the connectome  K + gamma*(D-W)
%
% THE ARROWS (2026-09-25):
%   rheome.operators.registry - every operator as a typed arrow: the field type it takes, the one
%                        it gives, and the shape of the matrix as a function of the domain's
%                        element counts. docs/OPERATOR-REGISTRY.md is canonical for the
%                        names; this table is canonical for the arrows, and the tests build
%                        each mesh operator on an icosphere and check the declared shape.
%   rheome.operators.check    - is this field admissible here, and what comes back
%
% See also: rheome.domain.of, rheome.fieldtype.registry
%
% Author: Diellor Basha, 2026
