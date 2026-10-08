function out = chirality(out, vorticity)
% DETECT.CHIRALITY  Annotate tracked ridges with the handedness of the field they sit in.
%
%   out = rheome.detect.chirality(out, vorticity)
%
% Detection runs on |W|, which discards sign, so a ridge knows where and how big a motif is but
% not which way it turns. This samples the SIGNED vorticity at each ridge's own (vertex, frame)
% and attaches the result, so handedness is carried rather than recovered later.
%
% Per ridge it adds:
%   .chirality      [nFrames x 1]  +1 counter-clockwise, -1 clockwise, at each step of the ridge
%   .chiralityMean  mean over the ridge, in [-1 1]
%   .nFlips         sign changes along the ridge
%   .halfCycles     nFlips + 1 -- how many half cycles the ridge spans
%
% ⚠ READ .nFlips BEFORE .chiralityMean. A ridge tracked on the ENVELOPE persists across sign
% reversals by design, so it spans several half cycles and its mean handedness tends to zero --
% that is the envelope working, not an absence of chirality. A ridge with nFlips = 0 is
% single-handed; one with nFlips = 5 has turned both ways three times and its mean says little.
%
% INPUTS:
%   out        from rheome.detect.ridges
%   vorticity  [nVert x nT] SIGNED field on the same grid the ridges were tracked on
%              (e.g. obs.vorticity from rheome.flow.observables)
%
% See also: rheome.detect.ridges, rheome.flow.observables
%
% Author: Diellor Basha, 2026

    if isempty(out.ridges), return; end
    [nV, nT] = size(vorticity);
    for k = 1:numel(out.ridges)
        v = out.ridges(k).v(:);  tt = out.ridges(k).t(:);
        if any(v > nV) || any(tt > nT)
            error('detect:chirality:range', ...
                'ridge %d indexes vertex %d / frame %d but the field is [%d x %d].', ...
                k, max(v), max(tt), nV, nT);
        end
        s = sign(vorticity(sub2ind([nV nT], v, tt)));
        out.ridges(k).chirality     = s;
        out.ridges(k).chiralityMean = mean(s);
        out.ridges(k).nFlips        = sum(diff(s) ~= 0);
        out.ridges(k).halfCycles    = sum(diff(s) ~= 0) + 1;
    end
    out.summary.chiralityMean = [out.ridges.chiralityMean];
    out.summary.nFlips        = [out.ridges.nFlips];
    out.summary.singleHanded  = mean([out.ridges.nFlips] == 0);
end

% Author: Diellor Basha, 2026
