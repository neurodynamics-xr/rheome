%% The resolution floor
% *Explanation.* What the MEG array and the minimum-norm inverse let the
% flow maps say, and the rules of reading that follow. The numbers on this
% page are those of the methods paper (in preparation).
%
%% Fusion adds no information
% Every flow map is a linear image of the minimum-norm estimate, so it
% inherits the estimate's point spread. The operator is not the limit; the
% instrument is.
%
%% The resolution matrix
% For the amplitude kernel $K$ and the leadfield $G$, the resolution matrix
% is $R = KG$. Its eigenvalues are the Wiener gains of the inverse, so its
% trace is the number of independent source directions the data determine.
% |rheome.inverse.resolution| reduces $R$ to lengths:
%
% * $r_{50}$, the geodesic radius holding half the point-spread power of a
%   point source, with the peak localisation error and the fraction of
%   power on the wrong hemisphere;
% * the modulation-transfer function over normally oriented modes,
%   summarised by its centroid. Its half-maximum is not used: it jumps
%   between modes of a non-monotone curve.
%
%%
% A measured conversion, $r_{50} = 0.196\,\lambda$, turns an $r_{50}$ into
% the wavelength at which a structure is localised. A noise-normalised
% kernel (dSPM, sLORETA) is no longer a resolution matrix: resolution is
% measured for the amplitude estimate only.
%
%% What the paper measures
% Across several hundred resting recordings from two cohorts, the median
% $r_{50}$ is 48.8 and 47.3 mm, with interquartile ranges under 4 mm, and it
% roughly doubles from the shallowest to the deepest cortex. With a
% recording's own background as noise, planted sources and vortices are
% located to 7-22 mm at wavelengths of 130 mm and above, down to -5 dB;
% vortices also at 92 mm above 0 dB; and nothing at 65 mm or finer, even
% without noise. These scales are the coarse units of the dyadic tree
% (depths 0-3), not anatomical parcels, most of which are smaller than the
% separation floor.
%
% The instrument also manufactures some quantities: the irrotational and
% solenoidal composition of an estimate, the size of a blob below about
% 20 dB, and vortex counts near the noise floor.
%
%% Rules of reading
% The methods paper states these rules for reading the operator above the
% floor:
%
% # Read flow from Helmholtz band maps at a declared scale of at least
%   about 130 mm (vortices to 92 mm at positive SNR), never from per-vertex
%   divergence or rotation.
% # Fix a source's band in advance; a vortex's band may be chosen blind.
% # Treat rings of divergence or rotation around a compact activation as
%   the point spread until a test shows otherwise.
% # Do not read irrotational/solenoidal composition, blob size below about
%   20 dB, or vortex count and size near the noise floor as properties of
%   the cortex.
% # Report a region's flux only where its own-region fraction
%   (|rheome.flow.crosstalk|) reaches 25%, and use sensitivity (|rheome.flow.sensitivity|)
%   as a mask, never as a divisor.
% # Quote $r_{50}$ and the MTF centroid, not the MTF half-maximum.
% # Exclude the reference poles' one-ring from every frame-dependent
%   readout; a singularity at a pole is the frame's (see
%   <about_gauge.html The gauge>).
%
%% What is not yet known
% The method states a fixed floor; it cannot flag, reading by reading, that
% a value is unresolved. Planted-recovery numbers come from one reference
% recording and one head model. The planted fields are synthetic fields on
% the analysis mesh, not tilted patches of phase-delayed dipoles, so the
% size and rate of the tilt on human cortex remain unmeasured. Treat the
% numbers above as the paper's, at its commit, not as constants of the
% toolbox; |rheome.scale.run| measures them for your own recordings.
%
%% See also
% <about_scale.html Scale: local to global, fast to slow>,
% <reference/rheome.inverse.resolution.html rheome.inverse.resolution>,
% <reference/rheome.flow.crosstalk.html rheome.flow.crosstalk>.
%
% _Written for Rheome @COMMIT@._
