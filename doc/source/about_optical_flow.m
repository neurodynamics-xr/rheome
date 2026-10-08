%% Surface optical flow: how a pattern moves
% *Explanation.* The toolbox has two flow pipelines, and they answer
% different questions. This page explains the second one and how it relates
% to the first.
%
%% Flow of the current, motion of a pattern
% The fused kernels (|rheome.flow.divergence|, |rheome.flow.curl|, |rheome.flow.potential|,
% |rheome.flow.stream|) differentiate the *current* $J$. They are linear, exact,
% and need one frame: where current springs, sinks and turns at one
% instant.
%
% Optical flow estimates how an *activity pattern* moves between frames.
% |rheome.flow.activation| reduces the current to one non-negative value per vertex
% and frame, and |rheome.flow.apparent| fits a velocity field $v$ to its change by
% Horn-Schunck optical flow on the surface (smoothness weight
% $\alpha = 1$ by default), then applies the same divergence and rotation
% operators to $v$.
%
% The two measure different things. A rotating current can sit inside a
% stationary envelope (rotation of $J$ large, rotation of $v$ zero), and a
% travelling bump can be irrotational everywhere (rotation of $J$ zero,
% divergence of $v$ large). On real resting alpha the paper finds the two
% unrelated (r = -0.11), which is expected and not a discrepancy.
%
%% A comparator and an input
% The method works on any vector field on the cortex. An optical-flow
% velocity field is one such field: it can be passed through the same
% operators, and it is a comparator for questions about propagation. Its
% results are never scored against, or pooled with, the divergence and
% rotation of the current.
%
%% What the estimate is, and is not
% * *Kinematic, not physical.* Cortical activity is generated locally; it
%   does not flow. $v$ is the apparent velocity of the pattern. Where
%   activity appears or vanishes instead of moving, the residual shows up
%   as divergence of $v$.
% * *Aperture-limited.* Only the component of $v$ along the intensity
%   gradient comes from the data; the smoothness prior supplies the rest,
%   so $\alpha$ sets how much of $v$ is data. Vary it and report the range.
% * *Envelope, not instantaneous magnitude.* The magnitude of a
%   band-limited current pulses at twice the carrier frequency and touches
%   zero; optical flow on it tracks the carrier. Use the analytic envelope
%   in a narrow band (the default of |rheome.flow.activation|), formed at the
%   sensors.
% * *Nonlinear and two frames.* Nothing in this pipeline can be fused into
%   a kernel, so it is far slower than the first; work on one hemisphere.
%
%%
% The scalar Horn-Schunck on the surface recovers a translating bump on a
% sphere with the right direction and 103% of its speed; the paper reports
% its vector variant as not validated. In both cohorts of the paper the
% apparent motion of resting alpha is slow, about 22 mm/s.
%
%% See also
% <about_sensor_to_cortex.html From sensors to cortex>,
% <reference/rheome.flow.apparent.html rheome.flow.apparent>,
% <reference/rheome.flow.activation.html rheome.flow.activation>.
%
% _Written for Rheome @COMMIT@._
