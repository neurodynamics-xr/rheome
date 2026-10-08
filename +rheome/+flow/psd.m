function out = psd(F, sfreq, coeffOperator, opts)
% FLOW.PSD  Welch PSD of flow-kernel mode coefficients (for fitting, not for filtering).
%
%   out = rheome.flow.psd(F, sfreq, coeffOperator)
%   out = rheome.flow.psd(F, sfreq, coeffOperator, opts)
%
% Averages periodograms over overlapping Hamming windows, following bst_psd: demean each window,
% taper, FFT with no zero-padding, keep the one-sided spectrum, apply the kernel to the FFT
% coefficients, convert to power, average across windows.
%
% WHY THIS EXISTS SEPARATELY FROM rheome.flow.joint. rheome.flow.joint returns a SINGLE periodogram, because it
% must stay complex and invertible -- it is what the pipeline filters and synthesises from. A
% single periodogram is chi-squared with 2 degrees of freedom per bin (~100%% variability), and it
% does NOT improve with a longer recording: more bins, each just as noisy.
%
% MEASURED, ON subject01, 30 s, 1-45 Hz (800 modes):
%
%                 min chi  max chi  med chi   med r2
%   single FFT       0.64     1.49     0.99    0.267
%   Welch (29 win)   0.60     1.37     0.92    0.841
%
% The EXPONENT is barely affected -- IQR 0.147 vs 0.140, no negative values either way. Two
% reasons: the log-periodogram bias for 2 dof is a CONSTANT (-gamma), so it shifts the fitted
% offset and leaves the slope alone; and least squares over ~1300 bins already averages the
% scatter. The regression does the averaging Welch would have done.
%
% What breaks is R^2: 0.267 on the periodogram is not a poor fit, it is chi-squared scatter
% dominating the residual. Using R^2 as a quality gate on a raw periodogram rejects good fits.
%
% So use this for the FIT (interpretable R^2, 15x fewer bins, faster) and rheome.flow.joint for the
% FILTERING. Pass as opts.fitPower to rheome.spectral.decompose. Exponents already estimated from a
% periodogram do not need revisiting.
%
% INPUTS:
%   F              [C x nT] sensor time series
%   sfreq          sampling rate (Hz)
%   coeffOperator  [Ks x C] flow kernel
%   opts  .winLength  window length in SECONDS (default 2)
%         .overlap    fractional overlap (default 0.5)
%         .band       [f1 f2] restrict the returned bins (default: all positive)
%
% OUTPUT (struct out):
%   .P      [Ks x nF] mean power per mode per bin
%   .f      [1 x nF]  bin frequencies (Hz), DC excluded
%   .nWin   windows averaged      .df  bin spacing (Hz)
%   .dof    2*nWin -- the chi-squared degrees of freedom the estimate carries
%
% See also: rheome.flow.joint, rheome.spectral.decompose, rheome.spectral.aperiodic
%
% Author: Diellor Basha, 2026

    if nargin < 4, opts = struct(); end
    if ~isfield(opts,'winLength') || isempty(opts.winLength), opts.winLength = 2;   end
    if ~isfield(opts,'overlap')   || isempty(opts.overlap),   opts.overlap   = 0.5; end
    if ~isfield(opts,'band')      || isempty(opts.band),      opts.band      = [];  end

    [~, nT] = size(F);
    Lwin = round(opts.winLength * sfreq);  Lwin = Lwin - mod(Lwin,2);
    if Lwin < 8 || Lwin > nT
        error('flow:psd:win', 'window of %d samples does not fit %d samples of data.', Lwin, nT);
    end
    Lov  = round(Lwin * opts.overlap);
    nWin = floor((nT - Lov) / (Lwin - Lov));
    NFFT = Lwin;                                   % no zero-padding, as bst_psd
    fAll = sfreq/2 * linspace(0, 1, NFFT/2+1);
    win  = 0.54 - 0.46*cos(2*pi*(0:Lwin-1)/(Lwin-1));   % hamming
    npg  = sum(win.^2);                                  % noise power gain

    K = size(coeffOperator,1);
    S1 = zeros(K, NFFT/2+1);
    for iw = 1:nWin
        idx  = (1:Lwin) + (iw-1)*(Lwin - Lov);
        Fw   = F(:,idx) - mean(F(:,idx), 2);
        Ff   = fft(Fw .* win, NFFT, 2);
        TF   = Ff(:, 1:NFFT/2+1) * sqrt(2 / (sfreq * npg));
        TF(:, [1 end]) = TF(:, [1 end]) / sqrt(2);       % DC and Nyquist are not doubled
        S1 = S1 + abs(coeffOperator * TF).^2;
    end
    P = S1 / nWin;

    keep = fAll > 0;
    if ~isempty(opts.band)
        keep = keep & fAll >= opts.band(1) & fAll <= opts.band(2);
    end
    out.P    = P(:, keep);
    out.f    = fAll(keep);
    out.nWin = nWin;
    out.df   = sfreq / NFFT;
    out.dof  = 2 * nWin;
end

% Author: Diellor Basha, 2026
