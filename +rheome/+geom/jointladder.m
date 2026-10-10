function L = jointladder(T, opts)
% GEOM.JOINTLADDER  The joint (space x time) bookkeeping table: one row per observation window.
%
%   L = rheome.geom.jointladder(T)                          % T = rheome.geom.tree(S, L=L, M=M, MaxDepth=8)
%   L = rheome.geom.jointladder(T, Fs=600, FrequencyLimits=[1 64], SamplesPerCycle=30)
%
% ⭐⭐ A ROW IS A SELECTION, NOT A MEASUREMENT. Each row is an observation window of the joint
% eigenmode-frequency plane, written in the spatial and temporal units it is easier to reason in:
% a tile of diameter D (a rheome.geom.tree depth, standing for a band of Laplace-Beltrami eigenvalues), a
% window of length W (a constant-Q octave's support, standing for a band of frequencies), and a
% frame interval dt (a dyadic level of the time grid). The row says what CAN be measured there --
% which feature sizes, which frequencies, which speeds. What IS measured (envelope power, a peak's
% trajectory, a group velocity, a phase velocity) is chosen per question and stored at the row's
% address. Power rolls up exactly to coarser rows; a trajectory or a speed does not, and lives at
% its home row, as the spectral peak does in the time atlas.
%
% ⭐ SPEED IS THE RATE, IN METRES PER SECOND. A window whose tiles are D across and whose frames are
% dt apart follows at most one tile per frame, so its speed is D/dt; cells on the same diagonal
% share it. Three speeds are reported:
%   speedMaxMS    D / dt           one tile per frame: the ceiling
%   speedTrackMS  D / (MinDwell*dt) the fastest a tile tracker actually followed: MinDwell = 2.5
%                                   frames per tile, measured by planting
%   speedMinMS    D / W            one tile in the whole window: slower than this is "stationary"
%
% ⭐ THE FEATURE-SIZE COLUMN IS MEASURED, NOT ASSUMED: a tile of diameter D holds features of
% sigma >= D / TileSigma, and which ratio is right depends on the tracker.
%   TileSigma = sqrt(2)  tile <= the core RADIUS. What a PER-FRAME tile peak needs: a tile mean keeps a
%                        band-pass peak only up to that size (measured by planting: sigma 25 mm
%                        tracked on 23-33 mm tiles at 10 dB, lost on 47 mm and coarser).
%   TileSigma = 1        tile <= sigma, one depth finer. What TRACK-BEFORE-DETECT prefers
%                        (rheome.detect.tilepath): the path spreads over more of the peak's tiles, and its
%                        best cell was one depth finer than the sqrt(2) rule's (measured by planting,
%                        confirmed on an independent seed).
% The octave of eigenvalues a depth stands for is sigma in [sigmaMin, 2*sigmaMin], lambda = 2/sigma^2
% (the mexhat's peak, t = sigma^2/2), reported as lambdaLo/lambdaHi. Two tree depths are one spatial
% octave, because area halves per depth (rheome.geom.ladder).
%
% ⭐ THE WINDOW IS THE CONSTANT-Q SUPPORT. An octave [f, 2f] needs k/f seconds (rheome.ingest.support: k = 15
% at 4 voices), rounded up to a dyadic length anchored at 1 s: alpha's 8-16 Hz octave gets 2 s and
% 22.6 cycles, as in the time atlas. Frame levels run from the record's rate down to one frame per
% window; `isCarrier` marks the level nearest SamplesPerCycle frames per cycle of the band's centre
% (in log), which is where a PHASE measurement has to live. An ENVELOPE measurement can sit at any
% coarser level whose speedTrack still exceeds the motion -- the roll-up.
%
% ⚠ `nodeAdmissible` is the instrument's bound, not the ladder's: a tile smaller than 2*r50 (104 mm on
% a 270-channel CTF array, rheome.inverse.resolution) cannot be told from its neighbour through the inverse. Planted maps on
% the cortex are not bound by it; source estimates are.
%
% INPUTS
%   T                rheome.geom.tree table (depths >= 1 are used; depth 0 is the whole surface)
%   Fs               record rate, Hz (600)
%   FrequencyLimits  [lo hi] Hz, octaves anchored at Anchor ([1 64])
%   Anchor           octave anchor, Hz (1)
%   Voices           the time bank's voices, for the support constant k (4)
%   SamplesPerCycle  the carrier rate target (30)
%   MinDwell         frames per tile the tracker needs (2.5, measured)
%   R50mm            the inverse's localisation radius (52, rheome.inverse.resolution)
%   TileSigma        tile diameter / smallest feature sigma a depth holds (1; sqrt(2) for per-frame
%                    tile peaks -- see above)
%
% OUTPUT (table L, one row per depth x octave x frame level)
%   depth nTiles tileMM sigmaMinMM lambdaLo lambdaHi
%   fLo fHi fCenter windowS level frameHz framesPerCycle framesPerWindow isCarrier
%   speedMaxMS speedTrackMS speedMinMS nodeAdmissible
%
% See also: rheome.geom.jointcell, rheome.geom.ladder, rheome.select.ladder, rheome.ingest.support, rheome.geom.tiles
%
% Author: Diellor Basha, 2026

    arguments
        T table
        opts.Fs (1,1) double {mustBePositive} = 600
        opts.FrequencyLimits (1,2) double {mustBePositive} = [1 64]
        opts.Anchor (1,1) double {mustBePositive} = 1
        opts.Voices (1,1) double {mustBePositive} = 4
        opts.SamplesPerCycle (1,1) double {mustBePositive} = 30
        opts.MinDwell (1,1) double {mustBePositive} = 2.5
        opts.R50mm (1,1) double {mustBePositive} = 52
        opts.TileSigma (1,1) double {mustBePositive} = 1
    end
    k = rheome.ingest.support(1, Voices=opts.Voices).support_times_fc;

    dep = unique(T.depth);  dep = dep(dep >= 1);
    D = arrayfun(@(d) median(T.diameter(T.depth == d)), dep);          % metres
    nTile = arrayfun(@(d) sum(T.depth == d), dep);

    j0 = floor(log2(opts.FrequencyLimits(1) / opts.Anchor));
    j1 = ceil(log2(opts.FrequencyLimits(2) / opts.Anchor)) - 1;
    fLo = opts.Anchor * 2.^(j0:j1);  fHi = 2 * fLo;  fC = sqrt(fLo .* fHi);
    W = 2.^ceil(log2(k ./ fLo));                                       % dyadic, anchored at 1 s

    rows = {};
    for b = 1:numel(fLo)
        lMax = floor(log2(W(b) * opts.Fs));                            % one frame per window
        lv = 0:lMax;  fr = opts.Fs ./ 2.^lv;
        [~, ic] = min(abs(log(fr / (opts.SamplesPerCycle * fC(b)))));
        for i = 1:numel(dep)
            for q = 1:numel(lv)
                sm = D(i) / opts.TileSigma;                        % smallest sigma held
                rows{end+1, 1} = [dep(i) nTile(i) 1e3*D(i) 1e3*sm ...
                    1/(2*sm^2) 2/sm^2 fLo(b) fHi(b) fC(b) W(b) lv(q) fr(q) fr(q)/fC(b) ...
                    W(b)*fr(q) q == ic D(i)*fr(q) D(i)*fr(q)/opts.MinDwell D(i)/W(b) ...
                    1e3*D(i) > 2*opts.R50mm]; %#ok<AGROW>
            end
        end
    end
    L = array2table(vertcat(rows{:}), 'VariableNames', {'depth','nTiles','tileMM','sigmaMinMM', ...
        'lambdaLo','lambdaHi','fLo','fHi','fCenter','windowS','level','frameHz','framesPerCycle', ...
        'framesPerWindow','isCarrier','speedMaxMS','speedTrackMS','speedMinMS','nodeAdmissible'});
    L.isCarrier = logical(L.isCarrier);  L.nodeAdmissible = logical(L.nodeAdmissible);
    L.Properties.UserData = struct('Fs', opts.Fs, 'k', k, 'MinDwell', opts.MinDwell, ...
        'SamplesPerCycle', opts.SamplesPerCycle, 'R50mm', opts.R50mm, 'TileSigma', opts.TileSigma);
end

% Author: Diellor Basha, 2026
