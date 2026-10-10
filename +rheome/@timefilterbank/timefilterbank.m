classdef timefilterbank
% TIMEFILTERBANK  A designed constant-Q tight frame on the frequency axis, evaluated per member at its own rate.
%
%   tfb = rheome.timefilterbank(SignalLength, 'SamplingFrequency', fs)
%   tfb = rheome.timefilterbank(SignalLength, 'SamplingFrequency', fs, 'VoicesPerOctave', 4, ...
%                        'Anchor', 1, 'FrequencyLimits', [], 'Oversample', 2)
%
% Analytic kernels g_m(f) on the one-sided DFT grid of a record, itersine windows of two
% voices' width on xi = V*log2(f/Anchor), member k centred at xi = k, plus a low-pass and
% a high-pass member completing the frame so that S(f) = sum_m |g_m|^2 = 2 on the whole
% axis (0, fs/2): a TIGHT frame in the analytic gain-2 convention, under which
% sum_m ||w_m||^2 = sum x^2 for a real signal up to the DC and Nyquist bins.
%
% wt(tfb, x) returns each member's coefficients at the member's OWN rate: the inverse FFT
% of its bins alone, at Oversample x the bin count. Magnitudes are exact at those instants.
% This is the non-stationary Gabor transform in its constant-Q form, in this repo's
% vocabulary. Modelled on graphfilterbank; owns the range, the anchor and the shapes; not
% the data.
%
% ⚠ The grid is ANCHORED AT AN ABSOLUTE FREQUENCY (Anchor, 1 Hz): member k sits at
% Anchor*2^(k/V) for every record and every rate, so all records share scales.
% cwtfilterbank anchors at its top scale and shares only across power-of-two rates.
%
% ⚠ No Nyquist cap is needed: every kernel is compactly supported inside [0, fs/2], so it
% is band-limited on the grid by construction -- the property Morse's top scales lacked.
%
% ⚠ Periodic boundary: the record's ends meet. Coefficients within a member's support of
% either end are edge-conditioned; the support cone labels them (rheome.ingest.cone).
%
% ⚠ VOICES COST TIME LOCALISATION HERE. A member is two voices wide and compactly
% supported in frequency, so its time support (99.9 % of |psi|^2) grows with voices:
% support*fc = 4.8, 8.0, 11.5, 15.0, 22.4, 37.1 at 1, 2, 3, 4, 6, 10 voices, against 2.9
% for Morse at ten under the same definition (Q 1.4 ... 14.5; Morse 4.6). The default
% of four gives Q ~ 6 and keeps an 8 Hz band inside 2 s tiles. The transform's cost does
% not depend on voices (total coefficients = the axis's bandwidth x duration).
%
% rates(tfb) gives those per-member rates from the design alone, without transforming
% anything: on the reference bank every band member runs at 0.7*fc, so a tile of an octave
% costs the same number of the band's own samples wherever the octave sits (rheome.select.ladder).
%
% See also: rates, wt, cwtfilterbank, rheome.graphfilterbank, rheome.jointfilterbank, rheome.ingest.bank
%
% Author: Diellor Basha, 2026

    properties (SetAccess = immutable)
        SignalLength
        SamplingFrequency
        VoicesPerOctave
        Anchor
        Oversample
        MinBins = 8
    end

    properties (SetAccess = private)
        FrequencyLimits              % [fLo fHi] of the BAND members (edge members complete the frame)
    end

    properties (SetAccess = private, Hidden)
        K_          % [1 x M] member index on the voice axis (NaN for edge members)
        Kind_       % {1 x M} 'lowpass' | 'band' | 'highpass'
        Bins_       % [M x 2] first and last one-sided bin of each member's support
        Freq_       % [1 x nBins] one-sided grid, Hz
    end

    properties (Dependent, SetAccess = private)
        NumMembers
    end

    methods
        function obj = timefilterbank(N, varargin)
            p = inputParser;
            p.addParameter('SamplingFrequency', [], @(v) isscalar(v) && v > 0);
            p.addParameter('VoicesPerOctave', 4, @(v) isscalar(v) && v >= 1 && v == round(v));
            p.addParameter('Anchor', 1, @(v) isscalar(v) && v > 0);
            p.addParameter('FrequencyLimits', [], @(v) isempty(v) || (numel(v) == 2 && v(1) > 0 && v(2) > v(1)));
            p.addParameter('Oversample', 2, @(v) isscalar(v) && v >= 1);
            p.parse(varargin{:});
            o = p.Results;
            if isempty(o.SamplingFrequency)
                error('timefilterbank:fs', 'SamplingFrequency is required.');
            end
            if ~(isscalar(N) && N == round(N) && N >= 16)
                error('timefilterbank:length', 'SignalLength must be an integer >= 16.');
            end
            obj.SignalLength      = N;
            obj.SamplingFrequency = o.SamplingFrequency;
            obj.VoicesPerOctave   = o.VoicesPerOctave;
            obj.Anchor            = o.Anchor;
            obj.Oversample        = o.Oversample;
            obj = i_design(obj, o.FrequencyLimits);
        end

        function v = get.NumMembers(obj), v = numel(obj.K_); end
    end
end

function obj = i_design(obj, limits)
    N = obj.SignalLength;  fs = obj.SamplingFrequency;  V = obj.VoicesPerOctave;  a = obj.Anchor;
    nb = floor(N/2) + 1;
    f  = (0:nb-1) * fs / N;
    obj.Freq_ = f;
    T  = N / fs;
    % candidate band members: xi = k, support xi in (k-1, k+1) -> f in (a*2^((k-1)/V), a*2^((k+1)/V))
    fLoRep = obj.MinBins / ((2^(1/V) - 2^(-1/V)) * T);      % >= MinBins bins of support
    if isempty(limits)
        fLo = fLoRep;  fHi = fs / 2;
    else
        fLo = max(limits(1), fLoRep);  fHi = min(limits(2), fs / 2);
    end
    kLo = ceil(V * log2(fLo / a));                          % first member centred >= fLo
    kHi = floor(V * log2(fHi / a)) - 1;                     % last member whose support ends <= fHi: a*2^((k+1)/V) <= fHi
    if kHi < kLo
        error('timefilterbank:range', 'No band member fits between %.4g and %.4g Hz at %d voices.', fLo, fHi, V);
    end
    K = kLo:kHi;
    obj.K_    = [NaN, K, NaN];
    obj.Kind_ = [{'lowpass'}, repmat({'band'}, 1, numel(K)), {'highpass'}];
    obj.FrequencyLimits = a * 2.^([K(1) - 1, K(end) + 1] / V);   % where the band members' support starts/ends
    % bins of support per member
    M = numel(obj.K_);
    B = zeros(M, 2);
    for m = 2:M-1
        k = obj.K_(m);
        lo = a * 2^((k-1)/V);  hi = a * 2^((k+1)/V);
        i1 = find(f > lo, 1);  i2 = find(f < hi, 1, 'last');
        B(m, :) = [i1, i2];
    end
    % low-pass: from bin 1 (DC) to the first band member's centre; high-pass: from the last centre to Nyquist
    fc1 = a * 2^(K(1)/V);  fcE = a * 2^(K(end)/V);
    B(1, :)   = [1, find(f < fc1, 1, 'last')];
    B(end, :) = [find(f > fcE, 1), nb];
    obj.Bins_ = B;
end
% Author: Diellor Basha, 2026
