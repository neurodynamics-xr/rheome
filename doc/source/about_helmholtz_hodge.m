%% The Helmholtz-Hodge split
% *Explanation.* Why the toolbox reads sources and vortices from two scalar
% potentials, and what that reading can and cannot separate.
%
%% Two potentials for one current
% On a surface a tangent vector field has a divergence, where it springs
% and sinks, and a scalar rotation, where it turns. The Helmholtz-Hodge
% decomposition writes the current as
%
% $$J = \nabla\Phi + N\times\nabla\Psi + (J\cdot n)\,n + h,$$
%
% an irrotational part, the gradient of a scalar potential $\Phi$ (maxima
% are sources, minima sinks); a solenoidal part, the rotated gradient of a
% stream function $\Psi$ (its extrema are vortices); the normal reading;
% and a harmonic residual $h$. The potentials solve two Poisson problems,
% $L\Phi = B_{\mathrm{div}}J$ and $L\Psi = B_{\mathrm{rot}}J$, with $L$ the
% Laplace-Beltrami stiffness and $B_{\mathrm{div}}$, $B_{\mathrm{rot}}$ the
% weak (integrated) divergence and rotation (|rheome.differential.helmholtz|,
% |rheome.flow.potential|, |rheome.flow.stream|).
%
% The weak forms put the derivative on smooth test functions instead of on
% the data, so the potentials are smooth even where the current is not. The
% Poisson solve is diagonal in the Laplace-Beltrami eigenbasis,
% $\hat\Phi_k = \lambda_k^{-1}(\Phi_{LB}^\top B_{\mathrm{div}}J)_k$, with
% the constant mode of each hemisphere dropped, which is why the fused
% potential kernels are one matrix each. The potentials see the tangential
% reading of the current only (see <about_sensor_to_cortex.html From
% sensors to cortex>).
%
%% Why not read per-vertex divergence and rotation
% A per-vertex divergence or rotation map has no declared scale: its finest
% structure is set by the mesh and, after the inverse, by the point spread.
% On folded cortex neighbouring normals differ by tens of degrees, and the
% pointwise operators leak between divergence and rotation. The potentials
% are scalars, smooth, and can be filtered by scale exactly.
%
%% Helmholtz band maps
% |rheome.differential.helmholtzbands| filters $\Phi$ and $\Psi$ with a tight
% Laplace-Beltrami wavelet bank (log-itersine, two voices per octave,
% $\sum_m g_m^2 = 1$). Because the bank is tight and the Helmholtz stiffness
% is the Laplace-Beltrami $L$, the band energies are exact
% $\lambda$-weighted Dirichlet energies and sum to the total. A *Helmholtz
% band map* is a potential filtered to one band; its extremum is the
% read-out location of a source or vortex of that scale.
%
% In the methods paper, planted sources and vortices were recovered in the
% correct band in 90 of 90 cases on a real cortex, from 315 to 23 mm, with
% under 1% of the energy in the wrong part. Through MEG the instrument sets
% the limit, not the operator (see <about_resolution_floor.html The
% resolution floor>).
%
%% The parts are not orthogonal on a mesh
% The four parts are orthogonal in the continuum but not on a discrete
% surface. Band filtering controls the leak between the parts; it does not
% make them orthogonal. In the paper the harmonic residual is 4-9% of the
% energy down to 46 mm and 14-18% at 23 mm. Report it with every
% decomposition (|harmFrac| in the output of |rheome.differential.helmholtzbands|).
%
%% Composition is not a property of the cortex
% How a field divides between irrotational and solenoidal parts after the
% inverse is set largely by the instrument: every minimum-norm estimate is
% a combination of leadfield rows, and the paper finds planted vortices and
% spirals with the same composition as an empty-room recording. Read where
% the band maps put their extrema, not the share of energy in each part.
%
%% See also
% <howto_flow_maps.html How to compute flow maps>,
% <reference/rheome.differential.helmholtzbands.html rheome.differential.helmholtzbands>,
% <reference/rheome.differential.helmholtz.html rheome.differential.helmholtz>.
%
% _Written for Rheome @COMMIT@._
