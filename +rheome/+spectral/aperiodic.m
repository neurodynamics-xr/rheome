function ap = aperiodic(P, f, opts)
% SPECTRAL.APERIODIC  Fit the aperiodic (1/f-like) background of power spectra.
%
%   ap = rheome.spectral.aperiodic(P, f)
%   ap = rheome.spectral.aperiodic(P, f, opts)
%
% Fits the specparam / FOOOF aperiodic model (Donoghue et al. 2020) to each COLUMN of P by
% iteratively down-weighting the points that sit above the current fit -- i.e. the peaks -- and
% refitting on what remains. The peaks themselves are never fitted: only the background is, which
% is all rheome.spectral.gains needs.
%
% MODEL (log10 power vs log10 frequency):
%   no knee   L(f) = b - chi*log10(f)                   a line in log-log
%   knee      L(f) = b - log10(k + f^chi)               flattens below f ~ k^(1/chi)
%
% INPUTS:
%   P     [nF x nS] LINEAR power, one spectrum per column (nS spectra; e.g. one per mode)
%   f     [nF x 1]  frequencies (Hz), strictly positive -- DC must be excluded by the caller
%   opts  .knee     false (default) | true      fit the knee term
%         .nIter    peak-rejection iterations (default 3)
%         .thresh   reject points above fit + thresh*sd of the residual (default 2)
%         .range    [fmin fmax] restrict the fit (default: all of f)
%
% OUTPUT (struct ap):
%   .Pap      [nF x nS] fitted aperiodic power, LINEAR units, evaluated at every f
%   .exponent [1 x nS]  chi
%   .offset   [1 x nS]  b
%   .knee     [1 x nS]  k (zero when opts.knee is false)
%   .r2       [1 x nS]  goodness of fit in log space, on the retained points
%   .nUsed    [1 x nS]  points retained after peak rejection
%
% The exponent returned per column is chi(lambda) when the columns are spatial modes.
%
% See also: rheome.spectral.gains, rheome.spectral.decompose, rheome.spectral.carrier
%
% Author: Diellor Basha, 2026

    if nargin < 3, opts = struct(); end
    if ~isfield(opts,'knee')   || isempty(opts.knee),   opts.knee   = false; end
    if ~isfield(opts,'nIter')  || isempty(opts.nIter),  opts.nIter  = 3;     end
    if ~isfield(opts,'thresh') || isempty(opts.thresh), opts.thresh = 2;     end
    if ~isfield(opts,'range')  || isempty(opts.range),  opts.range  = [];    end

    f = double(f(:));
    if any(f <= 0)
        error('spectral:aperiodic:dc', 'f must be strictly positive (exclude DC before calling).');
    end
    if size(P,1) ~= numel(f)
        error('spectral:aperiodic:size', 'P has %d rows but f has %d entries.', size(P,1), numel(f));
    end
    P  = double(P);
    nF = numel(f);  nS = size(P,2);

    inFit = true(nF,1);
    if ~isempty(opts.range)
        inFit = f >= opts.range(1) & f <= opts.range(2);
    end
    if nnz(inFit) < 4
        error('spectral:aperiodic:range', 'fewer than 4 frequency points in the fit range.');
    end

    lf  = log10(f);
    lP  = log10(max(P, realmin));

    ap.Pap      = zeros(nF, nS);
    ap.exponent = zeros(1, nS);
    ap.offset   = zeros(1, nS);
    ap.knee     = zeros(1, nS);
    ap.r2       = zeros(1, nS);
    ap.nUsed    = zeros(1, nS);

    for s = 1:nS
        y    = lP(:,s);
        keep = inFit;
        b = 0; chi = 0; k = 0;

        for it = 1:opts.nIter
            if opts.knee
                [b, chi, k] = i_fit_knee(f(keep), y(keep));
                yhat_keep   = b - log10(k + f(keep).^chi);
            else
                A         = [ones(nnz(keep),1), -lf(keep)];
                coef      = A \ y(keep);
                b = coef(1);  chi = coef(2);  k = 0;
                yhat_keep = b - chi*lf(keep);
            end
            if it == opts.nIter, break; end
            % reject peaks: points sitting ABOVE the fit by more than thresh*sd
            r  = y(keep) - yhat_keep;
            sd = std(r);
            if ~(sd > 0), break; end
            idx  = find(keep);
            drop = idx(r > opts.thresh*sd);
            keep(drop) = false;
            if nnz(keep) < 4, keep = inFit; break; end   % refuse to over-prune
        end

        if opts.knee
            yhat_all  = b - log10(k + f.^chi);
            yhat_keep = yhat_all(keep);
        else
            yhat_all  = b - chi*lf;
            yhat_keep = yhat_all(keep);
        end
        res = y(keep) - yhat_keep;
        tot = y(keep) - mean(y(keep));

        ap.Pap(:,s)    = 10.^yhat_all;
        ap.exponent(s) = chi;
        ap.offset(s)   = b;
        ap.knee(s)     = k;
        ap.r2(s)       = 1 - sum(res.^2)/max(sum(tot.^2), eps);
        ap.nUsed(s)    = nnz(keep);
    end
end

% Knee model: grid the exponent, solve the offset in closed form, keep the best.
% Cheap and derivative-free -- the knee is weakly identified and a fine search buys nothing.
function [b, chi, k] = i_fit_knee(f, y)
    chiGrid = 0.1:0.1:6;
    kGrid   = [0, 10.^(-2:0.5:3)];
    best = inf;  b = 0;  chi = 1;  k = 0;
    for c = chiGrid
        fc = f.^c;
        for kk = kGrid
            t   = log10(kk + fc);
            bb  = mean(y + t);              % least-squares offset for this (chi,k)
            err = sum((y - (bb - t)).^2);
            if err < best, best = err; b = bb; chi = c; k = kk; end
        end
    end
end

% Author: Diellor Basha, 2026
