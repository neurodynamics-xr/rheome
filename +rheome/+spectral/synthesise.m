function X = synthesise(nT, fs, chi, opts)
% SPECTRAL.SYNTHESISE  Synthesise real time series with a prescribed 1/f-like background.
%
%   X = rheome.spectral.synthesise(nT, fs, chi)
%   X = rheome.spectral.synthesise(nT, fs, chi, opts)
%
% Draws circular complex white noise, shapes it by the specparam aperiodic model, forces
% Hermitian symmetry and inverse-transforms. The result is the synthesis counterpart of
% rheome.spectral.aperiodic: fitting the output reproduces the chi that went in.
%
% MODEL (matching rheome.spectral.aperiodic's knee form)
%   S(f) = 1 / (k + f^chi)                    k = 0 gives the plain power law
%
% INPUTS:
%   nT    scalar   samples per series
%   fs    scalar   sampling rate (Hz)
%   chi   scalar   aperiodic exponent
%   opts  .knee    k (default 0)
%         .nS      number of series, one per ROW of X (default 1)
%         .power   [nS x 1] in-band power to match, or [] to return unit-variance rows
%                  matched EXACTLY, on the realisation, by Parseval -- not in expectation
%         .band    [fmin fmax] the band .power refers to (default [1 fs/2])
%         .shape   [nF x 1] LINEAR power to reproduce, in place of the chi/knee model
%         .shapeF  [nF x 1] the frequencies .shape is given at (Hz, strictly positive)
%
% OUTPUT:
%   X     [nS x nT] real series, rows independent
%
% ⚠ The Hermitian mirror is load-bearing. Without it ifft returns a complex series, and
%   taking real() of an unmirrored draw halves the variance at every non-DC bin -- the
%   spectrum keeps its SHAPE, so chi still fits, and only the power is silently wrong.
%   DC and, for even nT, Nyquist must be real for the same reason.
%
% ⭐⭐ .shape REPRODUCES A MEASURED CURVE INSTEAD OF A MODEL, which is what to use when the model
%   is pinned at a bound. rheome.spectral.aperiodic searches the knee on `[0, 10.^(-2:0.5:3)]`, so a knee
%   of exactly 1000 means the grid ran out rather than the fit converging -- real MEG sensor data
%   does that -- and a background synthesised from those parameters is then not free to match the
%   low end. Passing the fitted curve itself as .shape sidesteps the parameterisation entirely.
%   Interpolation is log-log and the end values are HELD outside [min(shapeF), max(shapeF)], so
%   the shape must span every band that will be read back.
%
% ⚠ .power is the realised in-band power, so it is matched exactly. The inline version in
%   synth_background_omega.m instead set each mode's TOTAL variance to its in-band power,
%   which under-fills the band; prefer this one.
%
% ⭐ Shape and level are set independently: chi and knee fix the shape, opts.power fixes the
%   level. Synthesising a background to sit under planted activity needs both.
%
% See also: rheome.spectral.aperiodic, rheome.spectral.gains, rheome.spectral.decompose
%
% Author: Diellor Basha, 2026

    if nargin < 4, opts = struct(); end
    if ~isfield(opts,'knee')  || isempty(opts.knee),  opts.knee  = 0;         end
    if ~isfield(opts,'nS')    || isempty(opts.nS),    opts.nS    = 1;         end
    if ~isfield(opts,'power'), opts.power = [];                               end
    if ~isfield(opts,'band')  || isempty(opts.band),  opts.band  = [1 fs/2];  end

    nT = round(nT);
    assert(nT >= 8, 'rheome.spectral.synthesise: nT must be at least 8, got %d.', nT);
    hasShape = isfield(opts,'shape') && ~isempty(opts.shape);
    if ~hasShape
        assert(isscalar(chi) && isfinite(chi), 'rheome.spectral.synthesise: chi must be a finite scalar.');
    end
    if ~isempty(opts.power)
        assert(numel(opts.power) == opts.nS, ...
            'rheome.spectral.synthesise: opts.power has %d entries for nS = %d.', ...
            numel(opts.power), opts.nS);
    end

    % the one-sided shape, with DC excluded rather than divided by
    ff = (0:nT-1) * fs / nT;
    ff(1) = ff(2);                                  % placeholder; DC is zeroed below
    if isfield(opts,'shape') && ~isempty(opts.shape)
        assert(isfield(opts,'shapeF') && numel(opts.shapeF) == numel(opts.shape), ...
            'rheome.spectral.synthesise: .shape needs a .shapeF of the same length.');
        sf = opts.shapeF(:);  sv = opts.shape(:);
        ok = sf > 0 & sv > 0 & isfinite(sf) & isfinite(sv);
        assert(sum(ok) >= 2, 'rheome.spectral.synthesise: .shape needs two positive finite points.');
        shape = 10.^interp1(log10(sf(ok)), log10(sv(ok)), log10(ff(:)), 'linear', 'extrap');
        % ⚠ hold the ends rather than extrapolate the slope, which runs away over decades
        shape(ff < min(sf(ok))) = max(sv(ok(1)), realmin);
        shape(ff > max(sf(ok))) = sv(find(ok, 1, 'last'));
        shape = reshape(shape, 1, []);
    else
        shape = 1 ./ (opts.knee + ff .^ chi);
    end
    shape(1) = 0;

    isEven = (mod(nT,2) == 0);
    nyq    = nT/2 + 1;                              % valid only when isEven
    hi     = 2 : ceil(nT/2);                        % the freely drawn positive bins

    X = zeros(opts.nS, nT);
    for i = 1:opts.nS
        w = (randn(1,nT) + 1i*randn(1,nT)) .* sqrt(shape);
        w(1) = 0;                                   % no DC
        if isEven, w(nyq) = real(w(nyq)) * sqrt(2); end   % Nyquist is its own conjugate
        w(nT - hi + 2) = conj(w(hi));               % ⚠ the mirror
        x = ifft(w);
        X(i,:) = real(x);                           % imaginary part is now ~1e-17
    end

    % level: unit variance, or rescale so the REALISED in-band power hits the target
    sd = std(X, 0, 2);
    sd(sd == 0) = 1;
    X = X ./ sd;
    if ~isempty(opts.power)
        for i = 1:opts.nS
            X(i,:) = X(i,:) * sqrt(opts.power(i) / i_inband(X(i,:), fs, opts.band));
        end
    end
end

function p = i_inband(x, fs, band)
% Exact in-band power of THIS realisation, by Parseval on the full-length DFT.
% ⚠ Scaling by the expected power of the SHAPE does not work: a single 1/f draw's in-band
%   power varies by a factor of 5 around its expectation at chi = 2, because the series is
%   dominated by a handful of low-frequency bins and the effective DOF is tiny.
    nT = numel(x);
    W  = abs(fft(x)).^2 / nT^2;
    k  = 1:floor(nT/2);                             % positive bins, excluding DC
    f  = k * fs / nT;
    in = f >= band(1) & f <= band(2);
    p  = 2 * sum(W(k(in)));                         % the 2 folds in the negative bins
    p  = max(p, realmin);
end

% Author: Diellor Basha, 2026
