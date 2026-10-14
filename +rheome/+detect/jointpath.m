function out = jointpath(Y, G, opts)
% DETECT.JOINTPATH  Track-before-detect across frequency: the best walk on the joint (tile x octave) graph in a window.
%
%   out = rheome.detect.jointpath(Y, G, Baseline=b, LogF=log2fc)
%   out = rheome.detect.jointpath(Y, G, Null=N, LogF=log2fc, SampleRate=18.75, NumPaths=3)   % N = jointpathnull(...)
%
% Y is a band-envelope field over ONE observation window, [nTile x nBand x nF]: tile means (rheome.geom.tiles)
% of each band's envelope, on a common frame grid (the time tile). G is rheome.geom.tiles. The joint graph's
% nodes are (tile, band); its edges are the tile graph inside a band and the band chain at one tile (the
% Cartesian product). A path is a walk on it: each frame it stays, steps to an adjacent tile, OR steps to an
% adjacent band. The best path is found exactly by rheome.detect.tilepath (the same Viterbi, local alignment)
% run on the product graph, so a pattern whose frequency drifts is followed across bands instead of being
% cut into one track per band.
%
% ⭐ BANDS MUST BE COMMENSURATE BEFORE A PATH CAN CROSS THEM. Envelopes fall ~1/f, so raw values would send
% every path to the lowest band. Norm (from a null, rheome.detect.jointpathnull) divides each (tile, band)
% series by its null median: a value is then "times this place-and-band's typical level", and one Baseline
% serves every band.
%
% ⭐ THE FREQUENCY READOUT IS A CENTROID, NOT THE PATH'S BAND. Bands are LogF steps apart (half an octave by
% default), so the band index is quantised; along the path, each frame's log2 frequency is the centroid of
% LogF over the path's band and its two neighbours at the path's tile, weighted by the normalised value above
% the baseline. .sweepOctS is the least-squares slope of that centroid over time (octaves per second).
%
% ⚠ One step per frame on the product graph: a pattern that changes tile AND band in one frame takes two
%   frames. The frame rate must make that rare (as for space: rheome.geom.jointcell's speedTrackMS).
% ⚠ The common frame grid oversamples the low bands' envelopes; the null (a surrogate of the same data on
%   the same grid) carries that correlation, so the threshold and false rate are honest for the grid used.
%
% INPUTS
%   Y           [nTile x nBand x nF] envelope tile values     G  rheome.geom.tiles (.A; .D for trackstats)
%   Null        a rheome.detect.jointpathnull struct: sets Norm and Baseline (overrides both)
%   Norm        [nTile x nBand] divisor (ones)               Baseline   subtracted after Norm (0)
%   LogF        [1 x nBand] log2 centre frequency of each band (1:nBand)
%   StepCost    per tile or band change (0)                  NumPaths   disjoint paths, best first (1)
%   SampleRate  frames per second (1)                        Span       "local" | "full" (rheome.detect.tilepath)
%
% OUTPUT (struct out) -- rheome.detect.tilepath's, so rheome.detect.trackstats reads the tile walk; per track also
%   .bands      [n x 1] band index per frame            .logf   [n x 1] centroid log2 frequency per frame
%   .sweepOctS  slope of .logf over time (oct/s; NaN under 3 frames)   .octSpan  max - min of .logf (octaves)
%   .tiles      the tile per frame (the band removed)   .values  the NORMALISED value per frame
%
% See also: rheome.detect.jointpathnull, rheome.detect.tilepath, rheome.detect.trackstats, rheome.scale.measure_freqtrack
%
% Author: Diellor Basha, 2026

    arguments
        Y double
        G struct
        opts.Null struct = struct([])
        opts.Norm double = []
        opts.Baseline (1,1) double = 0
        opts.LogF double = []
        opts.StepCost (1,1) double {mustBeNonnegative} = 0
        opts.NumPaths (1,1) double {mustBeInteger, mustBePositive} = 1
        opts.SampleRate (1,1) double {mustBePositive} = 1
        opts.Span (1,1) string {mustBeMember(opts.Span, ["local","full"])} = "local"
    end
    [nTile, nB, nF] = size(Y);
    if ~isempty(opts.Null), opts.Norm = opts.Null.norm;  opts.Baseline = opts.Null.baseline; end
    if isempty(opts.Norm), opts.Norm = ones(nTile, nB); end
    if isempty(opts.LogF), opts.LogF = 1:nB; end
    if ~isequal(size(opts.Norm), [nTile nB])
        error('detect:jointpath:norm', 'Norm must be [nTile x nBand] = [%d x %d].', nTile, nB);
    end
    Z = Y ./ opts.Norm;
    Ab = spdiags(ones(nB, 2), [-1 1], nB, nB);                      % the band chain
    J = struct('A', kron(speye(nB), sparse(double(G.A > 0))) + kron(Ab, speye(nTile)));   % node = tile + nTile*(band-1)
    out = rheome.detect.tilepath(reshape(Z, nTile*nB, nF), J, Baseline=opts.Baseline, StepCost=opts.StepCost, ...
                                 NumPaths=opts.NumPaths, SampleRate=opts.SampleRate, Span=opts.Span);
    lf = opts.LogF(:);  [out.tracks.bands] = deal([]);  [out.tracks.logf] = deal([]);
    [out.tracks.sweepOctS] = deal(NaN);  [out.tracks.octSpan] = deal(NaN);
    for q = 1:out.nTracks
        x = out.tracks(q);  node = x.tiles;
        k = mod(node - 1, nTile) + 1;  b = (node - k)/nTile + 1;
        c = zeros(numel(k), 1);
        for i = 1:numel(k)
            nb = max(b(i) - 1, 1):min(b(i) + 1, nB);
            w = max(squeeze(Z(k(i), nb, x.frames(i))) - opts.Baseline, 0);  w = w(:);
            if sum(w) > 0, c(i) = sum(w.*lf(nb))/sum(w); else, c(i) = lf(b(i)); end
        end
        sw = NaN;  if numel(c) >= 3, p = polyfit(x.frames(:)/opts.SampleRate, c, 1);  sw = p(1); end
        out.tracks(q).tiles = k;  out.tracks(q).bands = b;  out.tracks(q).logf = c;
        out.tracks(q).hops = sum(diff(k) ~= 0);  out.tracks(q).sweepOctS = sw;  out.tracks(q).octSpan = max(c) - min(c);
    end
end

% Author: Diellor Basha, 2026
