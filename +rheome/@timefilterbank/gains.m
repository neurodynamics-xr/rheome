function G = gains(obj, m, f)
% GAINS  A member's gain on given frequencies (Hz). Peak in-band response = tfb.Gain.
%
%   G = gains(tfb, m, f)     -> same size as f
%
% Band members: Gain * itersine(V*log2(f/Anchor) - k), zero outside |u| <= 1 (peak Gain).
% Edge members: Gain * sqrt(1 - itersine^2) on their side of the first/last centre, so the frame
% is tight on the whole axis: sum g^2 = Gain^2.
%
% Author: Diellor Basha, 2026

    V = obj.VoicesPerOctave;  a = obj.Anchor;
    f = double(f);
    xi = V * log2(max(f, realmin) / a);
    switch obj.Kind_{m}
        case 'band'
            G = obj.Gain * i_itersine(xi - obj.K_(m));
        case 'lowpass'
            K = obj.K_(2);
            S = i_itersine(xi - K).^2;                     % only the first band member reaches below its centre
            G = obj.Gain * sqrt(max(1 - S, 0)) .* (xi < K & f > 0);
            G(f == 0) = 0;
        case 'highpass'
            K = obj.K_(end-1);
            S = i_itersine(xi - K).^2;
            G = obj.Gain * sqrt(max(1 - S, 0)) .* (xi > K);
            G(f >= obj.SamplingFrequency/2) = 0;            % the Nyquist bin is not analytic
    end
end

function w = i_itersine(u)
    w = zeros(size(u));
    in = abs(u) <= 1;
    w(in) = sin(0.5*pi * cos(0.5*pi*u(in)).^2);
end
% Author: Diellor Basha, 2026
