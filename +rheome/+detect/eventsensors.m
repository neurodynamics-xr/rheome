function E = eventsensors(track, G, Gain, K, X, opts)
% DETECT.EVENTSENSORS  Which sensor samples and channels carry a tracked event.
%
%   E = rheome.detect.eventsensors(track, G, Gain, K, X, FrameSamples=fs)
%   E = rheome.detect.eventsensors(track, G, Gain, K, X, FrameSamples=fs, Fraction=0.5)
%       track   one element of rheome.detect.tiletrack / tilepath .tracks (.frames .tiles)
%       G       rheome.geom.tiles (.P [nV x nTile] membership)
%       Gain    [C x 3nV] (free) or [C x nV] (fixed) leadfield on the selected channels
%       K       [3nV x C] or [nV x C] imaging kernel: J = K * X, the estimate the tracker read
%       X       [C x nS] sensor data (raw, band-passed or analytic)
%       FrameSamples  [nF x 2] first and last sample of each tracker frame (frame f -> row f)
%
% ⭐ THE EVENT'S PART OF THE SENSOR DATA. At each frame the track sits on a tile; its sensor signal is
% the source estimate on that tile, forwarded:  x_ev(t) = Gain(:,tile) * K(tile,:) * x(t). That is
% exactly what the recording contributes to what the tracker followed, in the recording's own units,
% so it can be overlaid on the traces. ⚠ It is not a unique decomposition: neighbouring tiles' estimates
% overlap through the point spread, so x_ev of two tiles do not add to x.
% A channel CARRIES the event at a frame when it is in the smallest set holding Fraction of that frame's
% x_ev energy; E.mask marks those (channel, sample) cells.
%
% OUTPUT (struct E)
%   .Xevent [C x nS] the event's sensor signal (zero outside its frames)
%   .mask   [C x nS] logical, channel carries the event at that sample
%   .channels {nStep x 1} carrying channels per step (strongest first)   .samples [nStep x 2]
%   .power  [C x nStep] mean |x_ev|^2 per channel and step   .footprint [C x 1] summed over steps
%   .tiles  [nStep x 1] tile of each step   .fraction
%
% See also: rheome.show.eventsensors, rheome.detect.tiletrack, rheome.detect.tilepath, rheome.geom.tiles
%
% Author: Diellor Basha, 2026

    arguments
        track (1,1) struct
        G (1,1) struct
        Gain double
        K double
        X double
        opts.FrameSamples (:,2) double
        opts.Fraction (1,1) double {mustBeInRange(opts.Fraction, 0, 1, "exclude-lower")} = 0.5
    end
    [C, nS] = size(X);  nV = size(G.P, 1);
    if size(Gain, 1) ~= C || size(K, 2) ~= C || size(Gain, 2) ~= size(K, 1)
        error('detect:eventsensors:size', 'Gain [C x R], K [R x C] and X [C x nS] must agree (C = %d).', C);
    end
    per = size(Gain, 2) / nV;
    if per ~= 1 && per ~= 3, error('detect:eventsensors:size', 'Gain has %d columns; expects nV or 3nV (nV = %d).', size(Gain, 2), nV); end
    nStep = numel(track.frames);
    E.tiles = track.tiles(:);  E.samples = opts.FrameSamples(track.frames, :);  E.fraction = opts.Fraction;
    E.Xevent = zeros(C, nS, 'like', X);  E.mask = false(C, nS);
    E.power = zeros(C, nStep);  E.channels = cell(nStep, 1);
    for k = 1:nStep
        v = find(G.P(:, E.tiles(k)));
        r = reshape((v(:)' - 1) * per + (1:per)', [], 1);
        s = E.samples(k, 1):E.samples(k, 2);
        xe = Gain(:, r) * (K(r, :) * X(:, s));
        E.Xevent(:, s) = xe;
        p = mean(abs(xe).^2, 2);  E.power(:, k) = p;
        [ps, o] = sort(p, 'descend');
        n = find(cumsum(ps) >= opts.Fraction * sum(ps), 1);
        E.channels{k} = o(1:n);
        E.mask(o(1:n), s) = true;
    end
    E.footprint = sum(E.power, 2);
end

% Author: Diellor Basha, 2026
