function out = extract(obj, pageIndex)
% EXTRACT  Tabulate one page: [window x band x scale x ROI].
%
%   out = extract(fe, pageIndex)
%
% One flowpage per band (the CWT, which dominates), then per window the coefficient
% covariance and one inner product per (ROI, scale).
%
% ⚠ WINDOWS TILE THE CORE ONLY. A window straddling the margin would mix this page's samples
% with the neighbours' AND with cone-contaminated ones, and nothing downstream could tell.
% Whole windows only; a trailing partial window is dropped rather than being reported short.
%
% OUTPUT (struct out):
%   .Energy        [nWin x nBand x nScale x nROI]  area-weighted, summed over the window
%   .Density       the same, divided by WindowSamples and ROIArea -- the comparable form
%   .TotalEnergy   [nWin x nBand]           whole cortex, all scales
%   .ScaleEnergy   [nWin x nBand x nScale]  whole cortex, per scale
%   .ModeSpectrum  [nWin x nBand x Ks]  energy per LBO mode -- the full-resolution scale
%                  axis, for a centroid the 7-member bank is too coarse to resolve
%   .ROIArea       [1 x nROI]   .WindowSamples [nWin x nBand]   .Lambda [Ks x 1]
%   .tStart .tEnd  [nWin x 1] seconds on the recording's own axis
%   .WindowFrames  {nWin x nBand} frame indices used
%   .Bands .ScaleMM .ROILabel .ROIHemi .PageIndex .NumWindows .NumScales
%
% See also: totable, rheome.flowfeatures
%
% Author: Diellor Basha, 2026

    B  = obj.Bundle;
    nB = numel(obj.Bands);
    nS = obj.NumScales;
    nR = obj.Atlas.nScout;
    Ks = obj.Ks_;

    pages = cell(1, nB);
    for b = 1:nB
        pages{b} = rheome.flowpage(B, 'Band', obj.BandLimits{b}, 'PageIndex', pageIndex, ...
            'FreqLimits', B.freqLimits, 'VoicesPerOctave', B.voicesOct, ...
            'SamplesPerCycle', obj.SamplesPerCycle, 'Method', obj.Method);
    end

    % ---- the window grid, defined in TIME so it is identical across bands ----
    % Bands run at different rates, so a grid in samples would not line up between them and
    % window k of alpha would not be window k of beta.
    fp0  = pages{1};
    tCore = fp0.Time(fp0.Core);
    t0 = tCore(1);  t1 = tCore(end);
    % ⚠ THE CORE SPANS ONE SAMPLE MORE THAN (t1 - t0). A 10 s core at 400 Hz has 4000
    % samples covering 9.9975 s between its first and last TIMESTAMP, so flooring that
    % dropped a whole window per page -- 10% of the recording, silently.
    coreDur = (t1 - t0) + 1/fp0.Rate;
    nWin = floor(coreDur / obj.WindowSec + 1e-9);
    if nWin < 1
        error('flowfeatures:window', ...
            'Core is %.3f s but WindowSec is %g -- no whole window fits.', t1 - t0, obj.WindowSec);
    end
    tS = t0 + (0:nWin-1).' * obj.WindowSec;
    tE = tS + obj.WindowSec;

    out = struct();
    out.Energy        = zeros(nWin, nB, nS, nR);
    out.ScaleEnergy   = zeros(nWin, nB, nS);
    out.TotalEnergy   = zeros(nWin, nB);
    out.WindowSamples = zeros(nWin, nB);
    out.WindowFrames  = cell(nWin, nB);
    % Per-mode energy, summed over the window. FREE -- it is the diagonal of the same Gwin
    % the ROI traces use. Needed because a centroid read off the 7-member BANK is pinned
    % near the bank's own centre of mass (measured: 55-58 mm of range against 83-91 mm from
    % the modes), so the bank cannot resolve a scale change that the modes can.
    out.ModeSpectrum  = zeros(nWin, nB, Ks);

    for b = 1:nB
        fp = pages{b};
        H  = double(fp.ScaleGains);                       % [Ks x nS]
        C  = double(fp.Coefficients);
        for iw = 1:nWin
            idx = find(fp.Time >= tS(iw) & fp.Time < tE(iw) & fp.Core);
            out.WindowFrames{iw, b}  = idx;
            out.WindowSamples(iw, b) = numel(idx);
            if isempty(idx), continue; end

            Cw   = C(:, idx);
            Gwin = real(Cw * Cw');                        % [Ks x Ks] Hermitian -- the statistic.
                                                          % G_roi is real symmetric, so the
                                                          % imaginary part cannot contribute
            for g = 1:nS
                Mg = single((H(:,g) * H(:,g).') .* Gwin);
                % ⚠ CAST THE SMALL SIDE. The ROI Gram stack is [nROI x Ks^2] = 348 MB in
                % double; casting it here instead of Mg reallocated that on every one of
                % nWin*nScale*nBand calls and took the page from 2 s to 71 s.
                out.Energy(iw, b, g, :) = double(obj.ROIGram_ * Mg(:));
                % whole-cortex counterpart: Phi is orthonormal in the mass inner product, so
                % the field energy is the coefficient energy -- no vertex pass needed
                out.ScaleEnergy(iw, b, g) = sum(abs(H(:,g) .* Cw).^2, 'all');
            end
            out.TotalEnergy(iw, b)     = sum(abs(Cw).^2, 'all');
            out.ModeSpectrum(iw, b, :) = real(diag(Gwin));
        end
    end

    kc = double(reshape(centerWavenumbers(B.gfb), 1, []));
    out.Lambda      = double(B.Lambda(:));
    out.ROIArea     = full(obj.Atlas.Membership * double(B.wVert(:))).';
    out.Density     = out.Energy ./ max(out.WindowSamples, 1) ./ reshape(out.ROIArea, 1,1,1,[]);
    out.tStart      = tS;   out.tEnd = tE;
    out.Bands       = obj.Bands;
    out.ScaleMM     = 1000 * 2*pi ./ kc;
    out.ROILabel    = obj.Atlas.Label;
    out.ROIHemi     = obj.Atlas.Hemi;
    out.PageIndex   = pageIndex;
    out.NumWindows  = nWin;
    out.NumScales   = nS;
    out.Rate        = cellfun(@(p) p.Rate, pages);
end

% Author: Diellor Basha, 2026
