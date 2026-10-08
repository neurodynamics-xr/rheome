function g = grid(nT, fs, cfg)
% INGEST.GRID  The dyadic frame grid of a record: frames per level, their centres and extents.
%
%   g = rheome.ingest.grid(nT, fs, cfg)
%
% Pure arithmetic, no data. Level 0 tiles the record in frames of F = FrameFloor*fs samples
% anchored at sample 1 (the store's Time starts at 0). Level L frames are 2^L level-0
% frames. The last frame at every level may be partial; it keeps its NOMINAL centre and
% extent, and carries its real sample count in the statistics (rheome.ingest.reduce .n). Lmax is
% the first level whose single frame covers the record, so every parent at every level has
% exactly two children, the second possibly with n = 0.
%
% ⚠ Centres are of the nominal frame, not of the samples it holds. A query box is answered
% against the grid; how many samples a tile actually saw is a stored statistic, not a
% geometry.
%
% INPUTS:
%   nT   samples in the record
%   fs   sampling rate (Hz)
%   cfg  rheome.ingest.config
% OUTPUT (struct g):
%   .F        samples per level-0 frame          .K0      frames at level 0
%   .Lmax     top level (one frame)              .K       [1 x Lmax+1] frames per level
%   .tExtent  [1 x Lmax+1] seconds               .tCenter {1 x Lmax+1} seconds
%   .nT .fs .FrameFloor
%
% See also: rheome.ingest.config, rheome.ingest.reduce, rheome.ingest.rollup
%
% Author: Diellor Basha, 2026

    arguments
        nT  (1,1) double {mustBeInteger, mustBePositive}
        fs  (1,1) double {mustBePositive}
        cfg (1,1) struct
    end

    F = cfg.FrameFloor * fs;
    if abs(F - round(F)) > 1e-9 * max(1, F)
        error('ingest:grid:frameFloor', ...
              'FrameFloor*fs = %.6g is not an integer (FrameFloor %.4g s at %.6g Hz).', ...
              F, cfg.FrameFloor, fs);
    end
    F  = round(F);
    K0 = ceil(nT / F);
    Lmax = max(0, ceil(log2(K0)));
    L  = 0:Lmax;

    g = struct();
    g.F          = F;
    g.K0         = K0;
    g.Lmax       = Lmax;
    g.K          = ceil(K0 ./ 2.^L);
    g.tExtent    = cfg.FrameFloor * 2.^L;
    g.tCenter    = arrayfun(@(l) ((0:g.K(l+1)-1) + 0.5) * g.tExtent(l+1), L, 'UniformOutput', false);
    g.nT         = nT;
    g.fs         = fs;
    g.FrameFloor = cfg.FrameFloor;
end
% Author: Diellor Basha, 2026
