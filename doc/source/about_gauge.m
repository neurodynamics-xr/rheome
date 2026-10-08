%% The gauge: a frame with two chosen singularities
% *Explanation.* Why the toolbox fixes a tangent frame on the cortex, why it
% must have singularities, and why it puts them at two chosen poles.
%
%% Why a frame is needed
% Divergence, rotation and the Helmholtz potentials do not need a frame.
% Reading a flow per component, or by its orientation, does, and so does
% planting a field of given winding (|rheome.flow.seedvortex|). A frame is a pair
% of orthonormal tangent vectors $e_1, e_2$ at every vertex.
%
%% No smooth frame exists
% On a closed surface no tangent field is smooth everywhere. By the
% Poincare-Hopf theorem the indices of its singularities sum to the Euler
% characteristic, $\chi = 2$ for a hemisphere mesh. The frames that come
% for free are worse than that: the connection Laplacian's own half-edge
% frame turns by 79.1 degrees between neighbouring vertices on the
% paper's reference cortex, and vector diffusion, the smoothest generic
% choice, scatters its singularities over 20 faces.
%
%% A trivial connection with two poles
% The toolbox fixes the gauge by choice: a *trivial connection* (Crane,
% Desbrun and Schroeder, 2010), the least-norm change of discrete
% Levi-Civita transport that is flat everywhere except at prescribed faces.
% Two singularities of index +1 are prescribed, at two reference poles
% (|rheome.operators.gauge(V, F, Method="trivial")|,
% |rheome.operators.trivial_connection|). The face curvature that enters the solve
% is the principal value of the holonomy, so the integer parts carried by
% the connection's charts are not prescribed as extra singularities. The
% frame then has exactly two singularities, at the poles, and is parallel in
% its own connection everywhere else. On the reference cortex it is also the
% smoothest of the toolbox's frames (20.2 degrees between neighbours); on a
% sphere with the poles at the poles it is the geographic frame.
%
%% What depends on the gauge
% The topological charge of a flow singularity is read with the Levi-Civita
% connection and does not depend on the gauge. Anything defined through the
% frame does, near a pole: a field of winding $m$ built in the gauge reads
%
% $$m + (1 - m)\,k$$
%
% at a frame singularity of index $k$. A +1 field reads +1 whatever $k$
% is; a saddle ($m = -1$) built next to a +1 pole reads +1. Because the
% poles are fixed and known, a singularity at a pole is a property of the
% reference frame, never a finding. The paper places the poles where the
% signal is weakest and excludes their one-ring from every frame-dependent
% readout.
%
%% See also
% <reference/rheome.operators.gauge.html rheome.operators.gauge>,
% <reference/rheome.operators.trivial_connection.html rheome.operators.trivial_connection>,
% <reference/rheome.flow.seedvortex.html rheome.flow.seedvortex>.
%
% _Written for Rheome @COMMIT@._
