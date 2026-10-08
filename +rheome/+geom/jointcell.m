function [row, why] = jointcell(L, opts)
% GEOM.JOINTCELL  Pick the observation window for a measurement from the joint ladder.
%
%   row = rheome.geom.jointcell(L, SigmaMM=25, Hz=10, SpeedMS=0.25)                   % envelope / group
%   row = rheome.geom.jointcell(L, SigmaMM=25, Hz=10, SpeedMS=0.25, Measure="phase")  % at the carrier
%   [row, why] = rheome.geom.jointcell(...)          % why: a sentence per decision, or the refusal
%   row = rheome.geom.jointcell(L, SigmaMM=25, Hz=10, SpeedMS=0.5, SnrDB=10, PathMM=120)   % + detectability
%
% L is rheome.geom.jointladder. The pick is three decisions, each taken as COARSE as the question allows,
% because a coarser window pools more samples -- more averaging, fewer tiles, fewer frames:
%
%   octave   the one holding Hz
%   depth    the COARSEST depth whose tile still holds the feature: sigmaMinMM <= SigmaMM, i.e.
%            tileMM <= TileSigma * SigmaMM with the ladder's own TileSigma (rheome.geom.jointladder)
%   level    Measure="envelope" (default): the COARSEST frame level whose speedTrackMS still
%            covers SpeedMS -- the roll-up -- AND that still holds MinFrames frames per window.
%            Measure="phase": the carrier level, whatever the speed, because a phase is only
%            defined at ~SamplesPerCycle frames per cycle.
%
% ⚠ MinFrames IS WHAT STOPS A STILL FEATURE ROLLING UP TO NOTHING. Without it SpeedMS = 0 took the
% coarsest level, ONE frame per window (0.59 frames/s for alpha), where no track or path can exist.
%
% ⚠ A REFUSAL IS AN ANSWER. If no depth is fine enough for the feature, or no level fast enough for
% the speed, `row` is empty and `why` says which: a sigma-10 mm peak on a depth-8 tree (finest tile
% 23 mm > 14 mm) is refused, and plant_tiletrack_omega.m found it untrackable at every depth.
%
% ⚠ SpeedMS must also exceed the window's speedMinMS to be told from stationary; below it the
% cell is still returned but `why` says the motion is below one tile per window.
%
% ⭐⭐ DETECTABILITY IS A PROPERTY OF THE EVENT, NOT OF THE WINDOW. Given SnrDB, the row gains
%     eventS   = min(PathMM / SpeedMS, windowS)     how long the event exists in the window
%     nTau     = eventS / TauS                       independent noise samples it spans
%     evidence = 10^(SnrDB/10) * nTau                what track-before-detect can accumulate
%     detectable = evidence >= EvidenceMin
% TauS defaults to the octave's envelope correlation time 1/(fHi - fLo) (0.125 s for alpha).
% EvidenceMin = 45.6 is where rheome.detect.tilepath's hit rate reaches 0.8: fitted on one planted run
% (seed 13) and CHECKED on an independent one (seed 29) -- cells predicted detectable hit 0.99, the
% rest 0.38 (ladder_rules_omega.m). At 10 dB that is 4.6 correlation times (~0.57 s of alpha); at
% 20 dB 0.46. A fast mover over a short path is brief, and no rate or tile rescues a brief event:
% frame rate moved the hit rate by <= 0.2 at 10 dB. ⚠ SnrDB is the planted tests' FIELD SNR (peak
% over noise RMS of the pooled-to-be map), and EvidenceMin holds for their smooth, one-cycle noise;
% recalibrate both on real data.
%
% See also: rheome.geom.jointladder, rheome.geom.tiles, rheome.detect.tilepath, rheome.detect.tiletrack
%
% Author: Diellor Basha, 2026

    arguments
        L table
        opts.SigmaMM (1,1) double {mustBePositive}
        opts.Hz (1,1) double {mustBePositive}
        opts.SpeedMS (1,1) double {mustBeNonnegative} = 0
        opts.Measure (1,1) string {mustBeMember(opts.Measure, ["envelope","phase"])} = "envelope"
        opts.MinFrames (1,1) double {mustBePositive} = 8
        opts.SnrDB = []
        opts.PathMM (1,1) double {mustBePositive} = Inf
        opts.TauS = []
        opts.EvidenceMin (1,1) double {mustBePositive} = 45.6
    end
    row = L([], :);  why = strings(0, 1);
    band = L.fLo <= opts.Hz & opts.Hz < L.fHi;
    if ~any(band)
        why = sprintf("%.3g Hz is outside the ladder's octaves [%g %g) Hz.", opts.Hz, min(L.fLo), max(L.fHi));
        return
    end
    B = L(band, :);
    ts = i_tilesigma(L);
    ok = B.tileMM <= ts * opts.SigmaMM + 1e-9;
    if ~any(ok)
        why = sprintf(['No depth is fine enough: a sigma-%g mm feature needs tiles <= %.1f mm ' ...
            '(TileSigma %.3g); the finest here are %.1f mm.'], opts.SigmaMM, ts*opts.SigmaMM, ts, min(B.tileMM));
        return
    end
    d = max(B.depth(ok & B.tileMM == max(B.tileMM(ok))));
    B = B(B.depth == d, :);
    why(end+1) = sprintf("octave %g-%g Hz, window %g s", B.fLo(1), B.fHi(1), B.windowS(1));
    why(end+1) = sprintf("depth %d: %.1f mm tiles, the coarsest <= %.3g sigma = %.1f mm", d, ...
        B.tileMM(1), ts, ts*opts.SigmaMM);
    if opts.Measure == "phase"
        row = B(B.isCarrier, :);
        why(end+1) = sprintf("phase: carrier level %d, %.4g frames/s, %.1f per cycle", row.level, ...
            row.frameHz, row.framesPerCycle);
    else
        fast = B.speedTrackMS >= opts.SpeedMS & B.framesPerWindow >= opts.MinFrames;
        if ~any(fast)
            row = L([], :);
            why = sprintf(['Too fast: %.3g m/s exceeds the fastest trackable speed at depth %d ' ...
                '(%.3g m/s at %.4g frames/s).'], opts.SpeedMS, d, max(B.speedTrackMS), max(B.frameHz));
            return
        end
        row = B(find(fast, 1, 'last'), :);                 % levels ascend: last = coarsest
        why(end+1) = sprintf("envelope: level %d, %.4g frames/s, tracks up to %.3g m/s", row.level, ...
            row.frameHz, row.speedTrackMS);
    end
    if ~isempty(opts.SnrDB)
        tau = opts.TauS;  if isempty(tau), tau = 1 / (row.fHi - row.fLo); end
        ev = row.windowS;  if opts.SpeedMS > 0, ev = min(opts.PathMM / 1e3 / opts.SpeedMS, row.windowS); end
        row.eventS = ev;  row.nTau = ev / tau;  row.evidence = 10^(opts.SnrDB / 10) * row.nTau;
        row.detectable = row.evidence >= opts.EvidenceMin;
        if row.detectable
            why(end+1) = sprintf("detectable: %.2f s event = %.1f tau, evidence %.3g >= %.3g", ev, ...
                row.nTau, row.evidence, opts.EvidenceMin);
        else
            why(end+1) = sprintf(['⚠ NOT detectable: %.2f s event = %.1f tau at %g dB, evidence %.3g < %.3g ' ...
                '(needs %.2f s at this SNR)'], ev, row.nTau, opts.SnrDB, row.evidence, opts.EvidenceMin, ...
                opts.EvidenceMin * tau / 10^(opts.SnrDB / 10));
        end
    end
    if opts.SpeedMS > 0 && opts.SpeedMS < row.speedMinMS
        why(end+1) = sprintf("⚠ %.3g m/s is below one tile per window (%.3g m/s): indistinguishable from still", ...
            opts.SpeedMS, row.speedMinMS);
    end
    why = strjoin(why, "; ");
end

function ts = i_tilesigma(L)
    ts = sqrt(2);                                        % ladders built before TileSigma existed
    u = L.Properties.UserData;
    if isstruct(u) && isfield(u, 'TileSigma'), ts = u.TileSigma; end
end

% Author: Diellor Basha, 2026
