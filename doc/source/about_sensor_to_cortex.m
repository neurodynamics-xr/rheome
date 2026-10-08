%% From sensors to cortex: the surface flow operator
% *Explanation.* This page explains what the toolbox computes between a MEG
% recording and a map of cortical flow, and why each step is the one it is.
% The method is specified in the methods paper (in preparation); where the
% paper is provisional, so is this page.
%
%% The object: a current vector at every vertex
% Source imaging returns, at each cortical vertex, a current vector $J$ in
% A·m. MEG never records a single neuron: what a vertex carries is the sum
% of the moments of an activated patch of cortex, and on folded cortex such
% a patch spans columns pointing in different directions. Activation across
% the patch is not synchronous, so the summed moment is in general *tilted*
% off the vertex normal, and the tilt can change within a cycle.
%
% Two consequences shape the toolbox. The normal and tangential parts of
% $J$ are two *readings* of one tilted moment, not two generators: nothing
% in the toolbox attributes the normal part to pyramidal cells and the
% tangential part to horizontal connections. And the source estimate must be
% *unconstrained*, three components per vertex; an orientation-constrained
% estimate keeps only the normal projection, flips its sign across a sulcus,
% and turns a rotation of the moment into an amplitude change. The flow path
% refuses a constrained leadfield.
%
%% Step 1: a plain whitened minimum norm
% With $G$ the unconstrained leadfield ($C$ sensors by $3V$ source rows),
% $W$ the noise whitener and $WG = USV^\top$, the imaging kernel is
%
% $$K = V\,\mathrm{diag}\!\left(\frac{\Lambda s}{\Lambda s^2+1}\right)U^\top W,
%   \qquad J(t) = K\,b(t),$$
%
% with $\Lambda$ set by the assumed signal-to-noise ratio (|SnrFixed|,
% default 3) and the noise covariance regularised with a ridge of 0.1
% (|rheome.inverse.mne|). The amplitude measure keeps $J$ a physical current;
% dSPM and sLORETA are statistics, not currents, and are not
% differentiated.
%
% No source basis sits between the leadfield and the inverse. The
% relative-Dirac eigenbasis inverse (|Method = 'dirac'|) is available as an
% alternative, but it band-limits the current before the inverse is solved,
% so it is not the default, and the paper does not compare it with the
% minimum norm on any downstream number.
%
%% The analytic signal is formed at the sensors
% $K$ is real and does not change in time, so it commutes with the Hilbert
% transform: $\mathcal{H}(Kb) = K\,\mathcal{H}(b)$. Envelope and phase of
% any flow quantity therefore cost one complex product per frame, and the
% transform runs on $C$ channels instead of $3V$ source rows.
%
%% Step 2: divergence and rotation on the mesh
% On each triangle the gradient of a hat function is
% $\nabla\varphi_i = (N_f \times e_i)/2A_f$. Stacking these gives the
% ambient divergence and rotation of the vertex field,
%
% $$DJ = W_A\,(G_xJ_x + G_yJ_y + G_zJ_z), \qquad
%   RJ = W_A\,[(\nabla\times J)\cdot N_f],$$
%
% with $W_A$ the area-weighted average from faces back to vertices
% (|rheome.differential.divergence|, |rheome.differential.curl|). A positive divergence is
% a source; a positive rotation turns counter-clockwise seen from outside.
%
%% Two readings of the divergence
% On a curved surface the ambient divergence is not the surface divergence
% of the tangential part. They differ by the curvature coupling of the normal
% reading,
%
% $$DJ = \mathrm{div}_{\mathcal{M}} J_t + 2H\,(J\cdot n),$$
%
% with $H$ the mean curvature and $n$ the outward normal. The toolbox
% computes both. |rheome.flow.divergence| and |rheome.flow.curl| return the ambient
% operators $D$ and $R$, curvature term included. The weak forms, and with
% them the Helmholtz potentials |rheome.flow.potential| and |rheome.flow.stream| and their
% band maps, pair $J$ with test fields that lie in the surface, so they
% depend on the tangential reading alone.
%
% *Which of the two is the primary reading is not settled in this
% release.* The methods paper states it as a choice that must be declared,
% and its current choice is provisional. The ambient $D$ and $R$ are, up
% to sign, two components of the intrinsic surface Dirac operator applied
% to the whole moment written as a pure-imaginary quaternion (the identity
% holds to $3\times10^{-16}$); a reading of the whole moment as one
% quaternion-valued quantity is future work. Whichever you read, say which
% it is.
%
%% Step 3: one matrix per quantity
% Every stage is linear, so each flow quantity is one fixed matrix applied
% to the sensor vector. In the Laplace-Beltrami eigenbasis $\Phi_{LB}$
% (eigenvalues $\Lambda$, mass $M$):
%
% $$F_{\mathrm{div}} = \Phi_{LB}^\top M D K, \quad
%   F_{\mathrm{rot}} = \Phi_{LB}^\top M R K, \quad
%   F_\Phi = \Lambda^{-1}\odot\Phi_{LB}^\top B_{\mathrm{div}} K, \quad
%   F_\Psi = \Lambda^{-1}\odot\Phi_{LB}^\top B_{\mathrm{rot}} K,$$
%
% each of size modes by channels (|coeffOperator|), and the same products
% without $\Phi_{LB}^\top M$ give vertex maps (|vertexOperator|). Energy,
% enstrophy and helicity are quadratic forms $b^\top Q b$ (|gram|). The
% kernels depend on the anatomy, the head model, the noise covariance and
% the regularisation, not on the recording, so |rheome.flow.context| and
% |rheome.flow.build| compute them once per participant.
%
% Fusion is exact by linearity: in the paper the fused kernels equal
% reconstruct-then-differentiate to within $3\times10^{-13}$. It saves time
% and memory. *It adds no information:* every flow map is a linear image of
% the minimum-norm estimate and inherits its point spread (see
% <about_resolution_floor.html The resolution floor>).
%
%% See also
% <about_helmholtz_hodge.html The Helmholtz-Hodge split>,
% <about_optical_flow.html Surface optical flow>,
% <GettingStarted.html Getting Started>.
%
% _Written for Rheome @COMMIT@._
