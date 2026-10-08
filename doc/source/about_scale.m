%% Scale: local to global, fast to slow
% *Explanation.* Every reading the toolbox makes has a scale in space and a
% scale in time. This page explains how each is defined, how the two meet in
% the per-participant measures of |rheome.scale.run|, and why a scale must be
% declared rather than read off a map.
%
%% Space: from local to global
% The Laplace-Beltrami eigenmodes of a hemisphere, $L\phi_k = \lambda_k
% M\phi_k$, are ordered from global to local, and each carries a wavelength
% $\ell_k = 2\pi/\sqrt{\lambda_k}$. Scale enters as wavelets on that
% spectrum,
%
% $$W_m x = \Phi_{LB}\, g_m(\Lambda)\, \Phi_{LB}^\top M x,$$
%
% one diagonal per scale on a shared basis, so a filtered kernel is still
% one matrix (|graphfilterbank|). The bank used for the Helmholtz readout
% has two voices per octave; on the reference cortex of the methods paper
% its centre wavelengths run 315, 265, 183, 130, 92, 65, 46, 32 and 23 mm.
%
% The same octave ladder appears as a dyadic tree of the cortex itself
% (|rheome.geom.tree|): each level halves the area, and the median node diameter
% runs from 378 mm for the whole hemisphere through 267, 192 and 135 mm at
% depths 1-3. A wavelet at a given scale stands for the mixture of cortical
% sources the sensors pick up at that scale; tiles of the tree are the units
% that tracking and pattern statistics are pooled over.
%
% How fine the ladder can go is set by the basis: |rheome.flow.context| prints the
% finest usable member ($\sigma_m \ge 3.08/\sqrt{\lambda_{\max}}$) for the
% number of modes it loaded. That number moves with the basis, which is why
% the number of modes is pinned rather than left to whatever is cached.
%
%% Time: from fast to slow
% In time, the toolbox works in octaves: ten from 0.125 to 128 Hz, and
% constant-Q frames (|timefilterbank|) in which each band has its own rate.
% Each octave has its own signal-to-noise ratio at the sensors, measured
% against an empty-room recording, and the inverse is regularised for that
% band's own SNR ($\mathrm{SnrFixed} = 10^{\mathrm{dB}/20}$). Envelopes and
% phases are taken within one band, at the sensors.
%
%% Where space and time meet
% The instrument's resolution is the same in every band only in its shape.
% In the paper's two cohorts the point-spread radius varies by about 1.35
% across octaves, worsening above about 16 Hz where the sensor SNR
% collapses, while the number of independent locations the data determine
% varies far more. More signal buys more independent locations, not finer
% wavelengths. So a faster band is not a finer one, and a slower band is not
% a coarser one: the two axes are set separately.
%
%% The per-participant measures
% |rheome.scale.run| runs, for one participant, the measures that place a
% recording on these axes:
%
% * |resolution|: the point spread of the inverse, in millimetres;
% * |bandsnr|: the sensor SNR of each octave against the noise run;
% * |bandresolution|: the point spread at each octave's own SNR;
% * |periodicflow|: apparent motion of the alpha envelope, total and
%   periodic, over tiles of rest;
% * |grouptrack|: tile tracking of the alpha envelope over the dyadic tree
%   at several depths and frame rates, against a surrogate.
%
%%
% It writes compact tables and no figures, so that it can run once per
% participant on many participants; |rheome.scale.reduce| then gives the group
% distributions. The group summary is descriptive.
%
%% Declare the scale
% A per-vertex map has no scale of its own: its finest structure is the
% mesh's and, after the inverse, the point spread's. Every reading should
% name its band in space and in time. For sources, choose the spatial band
% before looking at the data; the band of largest energy is about three
% bands too fine in the paper's tests (for vortices it is a valid choice).
%
%% See also
% <about_resolution_floor.html The resolution floor>,
% <howto_scale_run.html How to run the per-participant measures>,
% <reference/rheome.scale.run.html rheome.scale.run>, <reference/rheome.geom.tree.html rheome.geom.tree>.
%
% _Written for Rheome @COMMIT@._
