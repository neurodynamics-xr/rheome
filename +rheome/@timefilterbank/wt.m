function C = wt(obj, x)
% WT  Sub-band transform: each member's coefficients at the member's own rate.
%
%   C = wt(tfb, x)     x [N x 1] real -> C struct array [1 x M]:
%       .coef      [Nm x 1] complex, the member's analytic coefficients at .t
%       .t         [Nm x 1] seconds, instants k*N/(Nm*fs) (may fall between samples)
%       .rate      Nm*fs/N, Hz
%       .weight    N/Nm: sum(|coef|^2)*weight is the member's energy (Parseval)
%       .kind      'lowpass' | 'band' | 'highpass'
%
% One FFT of the record, then per member an inverse FFT of ITS BINS ALONE at length
% Nm = the smallest 5-smooth number >= Oversample * (#bins). Placing the bins at the start
% of the short array demodulates the member by its first bin, which rotates the phase by
% a linear ramp and leaves magnitudes exact. Cost is the bank's total bandwidth times the
% duration, independent of the sampling rate; memory per member is Nm.
%
% ⚠ Exactness: |coef| equals |ifft(fft(x).*H_m)| AT THE INSTANTS .t. Those instants
% coincide with the record's samples only when N/Nm is an integer. Tile sums over them
% are sums over a different sampling of the same band-limited series.
%
% See also: freqz, rheome.timefilterbank
%
% Author: Diellor Basha, 2026

    N = obj.SignalLength;  fs = obj.SamplingFrequency;
    x = double(x(:));
    if numel(x) ~= N
        error('timefilterbank:length', 'x has %d samples; the bank was built for %d.', numel(x), N);
    end
    X = fft(x);
    f = obj.Freq_;
    M = obj.NumMembers;
    C = repmat(struct('coef', [], 't', [], 'rate', [], 'weight', [], 'kind', ''), 1, M);
    for m = 1:M
        b = obj.Bins_(m, :);  nb = b(2) - b(1) + 1;
        Nm = tfb_smooth(ceil(obj.Oversample * nb));
        Y = zeros(Nm, 1);
        Y(1:nb) = X(b(1):b(2)) .* gains(obj, m, f(b(1):b(2))).';
        C(m).coef   = ifft(Y) * (Nm / N);
        C(m).t      = (0:Nm-1)' * (N / Nm) / fs;
        C(m).rate   = Nm * fs / N;
        C(m).weight = N / Nm;
        C(m).kind   = obj.Kind_{m};
    end
end
% Author: Diellor Basha, 2026
