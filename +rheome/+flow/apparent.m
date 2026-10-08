function R = apparent(A, S, opts)
% FLOW.APPARENT  Apparent motion of an activation pattern, and its divergence and curl.
%
%   R = rheome.flow.apparent(A, S)
%   R = rheome.flow.apparent(A, S, Rate=300, Band=[8 16], Alpha=1)
%
% The second flow pipeline, and it answers a DIFFERENT question from rheome.flow.curl. Given a
% scalar activation map A [nV x nT] (rheome.flow.activation), estimate the velocity field v of the
% pattern by manifold Horn-Schunck optical flow, then apply the same differential operators
% to v that rheome.flow.curl and rheome.flow.divergence apply to the current.
%
% ⭐⭐ THE TWO PIPELINES MEASURE DIFFERENT THINGS AND MUST NOT BE COMPARED CASUALLY.
%   rheome.flow.curl / rheome.flow.divergence : curl and div of the CURRENT J. An instantaneous property
%       of the field at one time -- where current circulates, where it springs and sinks.
%       Linear in the data, exactly fused from the inverse, no time derivative anywhere.
%   rheome.flow.apparent               : curl and div of the VELOCITY v of the amplitude pattern.
%       A kinematic property of how a pattern MOVES between frames. Nonlinear in the data
%       (optical flow is a variational estimate), and it needs at least two frames.
% A rotating current can sit in a stationary envelope (curl J large, curl v zero) and a
% travelling bump can be irrotational everywhere (curl J zero, div v large). Both are real;
% they are simply different questions, and the vortex-over-the-alpha-phase question is the
% SECOND one.
%
% ⚠ OPTICAL FLOW IS KINEMATIC, NOT PHYSICAL. Cortical activity is generated locally, it does
% not advect: v is the APPARENT velocity of the pattern, in the sense in which a wave of
% stadium spectators moves. Where brightness constancy breaks -- activity created or
% destroyed rather than moved -- the residual appears as div v, so a source-sink motif here
% can be either genuine propagation converging or an amplitude change the model cannot
% express. That ambiguity is intrinsic to the method, not to this implementation.
%
% ⚠ IT IS AN APERTURE-LIMITED ESTIMATE. Horn-Schunck recovers only the velocity component
% along the intensity gradient (the aperture problem); the rest is supplied by the smoothness
% prior, so `Alpha` is not a cosmetic knob -- it sets how much of v is data and how much is
% regularisation. Vary it and report the range (.alpha is recorded for that reason).
%
% ⚠ v IS PER FRAME, NOT PER SECOND. rheome.dynamics.opticalflow_scalar differences adjacent frames,
% so |v| is metres per FRAME; .speed here is multiplied by Rate and is metres per second, and
% .speedPerFrame keeps the raw one. Getting this wrong scales every speed by the sample rate.
%
% ⚠ COST IS THE REASON TO WORK ON ONE HEMISPHERE. Each frame solves a sparse 2nV x 2nV
% system, so the whole cortex at 20484 vertices costs about four times a hemisphere and
% there is nothing to gain: the hemispheres are disconnected, the cached basis is already
% per-hemisphere, and 12% of a point source's power lands on the wrong one anyway
% (rheome.inverse.resolution). Pass B.L.S and the matching rows of the kernel.
%
% INPUTS:
%   A      [nV x nT] scalar activation (rheome.flow.activation)
%   S      surface struct .Vertices .Faces .VertNormals .nV -- the SAME mesh A lives on
%   opts (name-value):
%     Rate   frame rate of A in Hz. ⚠ REQUIRED for .speed to be in m/s; without it .speed is
%            returned empty and only .speedPerFrame is filled.
%     Band   [fLo fHi] recorded for provenance
%     Alpha  Horn-Schunck smoothness weight (1)
%     Frames restrict to these frame indices before estimating (all)
%
% OUTPUT (struct R):
%   .velocity   [3nV x nT-1] ambient tangent velocity, interleaved
%   .divergence [nV x nT-1]  div v: where the pattern converges (-) or emerges (+)
%   .vorticity  [nV x nT-1]  curl v: rotation of the pattern, CCW positive from outside
%   .speed      [nV x nT-1]  |v| in m/s          .speedPerFrame  |v| in m/frame
%   .rate .band .alpha .nV .nT .seconds  provenance of the estimate
%
% See also: rheome.flow.activation, rheome.flow.curl, rheome.flow.divergence, rheome.dynamics.opticalflow_scalar,
%           rheome.dynamics.flow_readout, rheome.differential.curl
%
% Author: Diellor Basha, 2026

    arguments
        A      double
        S      (1,1) struct
        opts.Rate   (1,1) double = NaN
        opts.Band         double = []
        opts.Alpha  (1,1) double {mustBePositive} = 1
        opts.Frames       double = []
    end
    if ~isreal(A) || any(A(:) < 0)
        error('flow:apparent:map', ...
            'A must be a real non-negative activation map; use rheome.flow.activation to build it.');
    end
    if size(A,1) ~= S.nV
        error('flow:apparent:mesh', ...
            'A has %d vertices but the surface has %d. Same mesh, or the gradient is nonsense.', ...
            size(A,1), S.nV);
    end
    if ~isempty(opts.Frames), A = A(:, opts.Frames); end
    if size(A,2) < 2
        error('flow:apparent:frames', 'Optical flow needs at least two frames; got %d.', size(A,2));
    end

    t0 = tic;
    v  = rheome.dynamics.opticalflow_scalar(A, S, struct('alpha', opts.Alpha));
    rd = rheome.dynamics.flow_readout(v, S);

    R = struct();
    R.velocity      = v;
    R.divergence    = rd.divergence;
    R.vorticity     = rd.vorticity;
    R.speedPerFrame = rd.speed;
    if isfinite(opts.Rate), R.speed = rd.speed * opts.Rate; else, R.speed = []; end
    R.rate = opts.Rate;  R.band = opts.Band;  R.alpha = opts.Alpha;
    R.nV = S.nV;  R.nT = size(A,2);  R.seconds = toc(t0);
end
% Author: Diellor Basha, 2026
