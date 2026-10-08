function s = support(obj, varargin)
% SUPPORT  Measured time support and bandwidth constants of the kernels.
%
%   s = support(tfb)                      on a probe grid (2^16 samples at fs)
%   s = support(tfb, 'Energy', 0.999)     fraction of |psi|^2 the support must hold
%
% A band member's time kernel is the inverse FFT of its gains. Its support is the
% shortest centred interval holding the given energy fraction; times the centre
% frequency it is a constant of the (shape, voices) pair, since the kernels are scale
% covariant. Its half-power width over the centre frequency is 1/Q.
%
% OUTPUT: s.timeTimesFc   support * fc (median over the probe's band members)
%         s.loOverFc, s.hiOverFc   half-power borders over fc
%         s.Q             fc / half-power width
%         s.perMember     [3 x nProbe] the raw values
%
% Author: Diellor Basha, 2026

    p = inputParser;  p.addParameter('Energy', 0.999);  p.parse(varargin{:});
    frac = p.Results.Energy;
    fs = obj.SamplingFrequency;  Np = 2^16;
    probe = rheome.timefilterbank(Np, 'SamplingFrequency', fs, 'VoicesPerOctave', obj.VoicesPerOctave, ...
                           'Anchor', obj.Anchor, 'Oversample', obj.Oversample);
    [H, f] = freqz(probe);
    kinds = probe.Kind_;  band = find(strcmp(kinds, 'band'));
    fc = probe.Anchor * 2.^(probe.K_(band) / probe.VoicesPerOctave);
    vals = zeros(3, numel(band));
    for i = 1:numel(band)
        m = band(i);
        h = zeros(Np, 1);  h(1:numel(f)) = H(m, :);          % one-sided -> analytic kernel
        psi = ifft(h);  e = abs(psi).^2;
        e = fftshift(e);  c = cumsum(e) / sum(e);             % centred at Np/2+1
        lo = find(c >= (1 - frac)/2, 1);  hi = find(c >= 1 - (1 - frac)/2, 1);
        vals(1, i) = (hi - lo) / fs * fc(i);
        g = H(m, :);  half = g >= max(g) / sqrt(2);
        vals(2, i) = f(find(half, 1)) / fc(i);
        vals(3, i) = f(find(half, 1, 'last')) / fc(i);
    end
    ok = vals(1, :) < 0.25 * Np / fs * fc & fc > 8 * fs / Np;   % untruncated, representable
    if nnz(ok) < 3, ok = true(1, numel(band)); end
    s = struct('timeTimesFc', median(vals(1, ok)), 'loOverFc', median(vals(2, ok)), ...
               'hiOverFc', median(vals(3, ok)), 'perMember', vals);
    s.Q = 1 / (s.hiOverFc - s.loOverFc);
end
% Author: Diellor Basha, 2026
