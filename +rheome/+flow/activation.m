function [A, info] = activation(J, opts)
% FLOW.ACTIVATION  A current field reduced to a scalar activation map, for optical flow.
%
%   A = rheome.flow.activation(J)                         % analytic envelope (the default)
%   A = rheome.flow.activation(J, Method="norm")           % instantaneous |J|
%   [A, info] = rheome.flow.activation(J, Method="analytic", Rate=300, Band=[8 16])
%
% Takes the unconstrained current J [3nV x nT] -- the raw minimum-norm estimate, interleaved
% [x1;y1;z1;x2;...] -- and returns one non-negative scalar per vertex per frame [nV x nT],
% which is what a brightness-constancy optical flow needs.
%
% ⚠⚠ THE INSTANTANEOUS NORM IS THE WRONG MAP, AND IT FAILS IN A SPECIFIC WAY. For a
% band-limited real signal J ~ a(t)cos(2 pi f t + phi), |J| goes as |cos| -- it PULSES AT 2f
% AND TOUCHES ZERO TWICE PER CYCLE. ⚠ that zero is exact only for a LINEARLY polarised
% current (one phase for all three components, i.e. a dipole of fixed orientation whose
% amplitude oscillates, which is the physiological case). Give the components independent
% phases and the current is elliptical and never passes through zero -- measured on a
% synthetic, 0.178 of peak instead of ~0. Real resting alpha behaves like the linear case. Brightness constancy (dI/dt + v.grad I = 0) then has to
% explain the carrier's own rise and fall, which no velocity field can do, so the solver
% returns whatever minimises the smoothness term at the zero crossings. The apparent motion
% it reports is the carrier, not the pattern. MEASURED on one 2 s resting alpha tile, left
% hemisphere, 300 Hz: the per-vertex power at 2*f0 relative to DC is 210x larger for the norm
% than for the envelope over the 100 strongest vertices (280x over the strongest 1000), and
% the norm falls to a median of 3.4% of each vertex's own peak within the tile against 21%
% for the envelope, with 2.5% of vertices dropping below 1% of peak against none. .carrierRatio
% and .dipRatio report both, per call, so the choice is never a matter of taste.
%
% ⚠ DO NOT MEASURE THE CARRIER ON THE SPATIAL MEAN. The 2f pulsation has a different phase at
% every vertex, so averaging over the cortex cancels it and both methods score zero -- which
% is what the first version of this check did, and it made the norm look fine.
%
% ⭐ USE THE ANALYTIC ENVELOPE (Method="analytic", the default). Hilbert-transform the
% band-limited signal along time, so J becomes complex and |J| is the ENVELOPE a(t) with the
% carrier removed. That is the quantity whose motion across the cortex is the thing being
% asked about -- how an alpha patch travels over the phase of the oscillation -- and it is
% smooth on the timescale optical flow assumes.
%
% ⚠ HILBERT BELONGS AT THE SENSORS, NOT HERE, WHEN YOU CAN PUT IT THERE. The inverse is a
% real linear time-invariant map, so hilbert(K*F) == K*hilbert(F) exactly, and the sensor
% version costs nCh = 270 transforms instead of 3nV = 61452. This function accepts an
% ALREADY COMPLEX J and simply takes its modulus, which is the cheap path; pass a real J and
% it does the transform here, correctly but 200x more expensively. Either way the answer is
% the same to floating point (asserted in tFlowApparent).
%
% ⚠ THE HILBERT TRANSFORM NEEDS A NARROW BAND TO MEAN ANYTHING. The analytic signal of a
% broadband field has no interpretable envelope: band-limit first (the alpha octave, 8-16 Hz,
% is band 6 of the store's ladder) and say so in Band, which is recorded and never applied --
% this function does not filter.
%
% INPUTS:
%   J     [3nV x nT] interleaved current; real (a band-limited field) or complex (already
%         analytic, e.g. a timefilterbank sub-band or hilbert applied at the sensors)
%   opts (name-value):
%     Method  "analytic" (default) | "norm"
%     Rate    the frame rate of J in Hz; recorded, and used for .samplesPerCycle
%     Band    [fLo fHi] the band J was limited to; recorded, and used for .samplesPerCycle
%     Normalise  scale A to unit maximum (false). ⚠ rheome.dynamics.opticalflow_scalar normalises
%                internally anyway, so leaving this false keeps A in A.m for reporting.
%
% OUTPUT:
%   A     [nV x nT] non-negative scalar activation
%   info  .method .rate .band .f0 .samplesPerCycle .carrierRatio .dipRatio .wasComplex
%         .carrierRatio  per-vertex power at 2*f0 over power at DC, median over the 100
%                        strongest vertices: two orders of magnitude larger for a norm
%         .dipRatio      median over vertices of min(A)/max(A) within the window: an envelope
%                        stays up (0.21 measured), a norm dips to the floor (0.034)
%
% See also: rheome.flow.apparent, rheome.dynamics.opticalflow_scalar, rheome.flow.rateplan, rheome.select.ladder
%
% Author: Diellor Basha, 2026

    arguments
        J
        opts.Method    (1,1) string {mustBeMember(opts.Method, ["analytic","norm"])} = "analytic"
        opts.Rate      (1,1) double = NaN     % ⚠ NOT mustBePositive: NaN means "not stated",
                                              % and a validator rejects its own default
        opts.Band            double = []
        opts.Normalise (1,1) logical = false
    end
    if mod(size(J,1), 3) ~= 0
        error('flow:activation:shape', 'J must be [3nV x nT] interleaved; got %d rows.', size(J,1));
    end
    nV = size(J,1)/3;  nT = size(J,2);
    wasComplex = ~isreal(J);

    if opts.Method == "analytic" && ~wasComplex
        if nT < 8
            error('flow:activation:short', ...
                'The analytic signal needs more than %d frames; band-limit and pass more.', nT);
        end
        J = hilbert(J.').';                      % ⚠ hilbert works down COLUMNS: transpose
    end

    % |J| per vertex: the 3-vector modulus, and for a complex J the modulus of the complex
    % 3-vector, which is the envelope of that vector's magnitude.
    A = sqrt(abs(J(1:3:end,:)).^2 + abs(J(2:3:end,:)).^2 + abs(J(3:3:end,:)).^2);

    if opts.Normalise, A = A / max(A(:) + eps); end

    f0 = NaN;
    if ~isempty(opts.Band), f0 = sqrt(opts.Band(1) * opts.Band(2)); end
    [cr, dip] = i_carrier(A, opts.Rate, f0);
    info = struct('method', opts.Method, 'rate', opts.Rate, 'band', opts.Band, 'f0', f0, ...
                  'samplesPerCycle', opts.Rate / f0, 'wasComplex', wasComplex, ...
                  'nV', nV, 'nT', nT, 'carrierRatio', cr, 'dipRatio', dip);
end

function [r, dip] = i_carrier(A, fs, f0)
% PER-VERTEX power at the 2f carrier relative to DC, median over the strongest vertices, and
% how close each vertex's map comes to zero within the window. ⚠ per vertex, NOT on the
% spatial mean: the pulsation's phase varies across the cortex and the mean cancels it.
    r = NaN;  dip = median(min(A, [], 2) ./ max(max(A, [], 2), eps));
    if ~isfinite(fs) || ~isfinite(f0) || size(A,2) < 16, return; end
    p = mean(A.^2, 2);  [~, o] = sort(p, 'descend');
    v = o(1:min(100, numel(o)));
    n = size(A,2);  f = (0:n-1) * fs / n;
    [~, i2] = min(abs(f - 2*f0));
    X = abs(fft(A(v,:), [], 2)).^2 / n^2;
    r = median(X(:, i2) ./ max(X(:,1), eps));
end
% Author: Diellor Basha, 2026
