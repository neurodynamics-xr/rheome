function G = gains(obj, m, f)
% GAINS  A member's gain on given frequencies (Hz). Analytic gain-2 convention.
%
%   G = gains(tfb, m, f)     -> same size as f
%
% Band members: sqrt(2) * itersine(V*log2(f/Anchor) - k), zero outside |u| <= 1.
% Edge members: sqrt(2 - S_bands(f)) on their side of the first/last centre, so the frame
% is tight on the whole axis.
%
% Author: Diellor Basha, 2026

    V = obj.VoicesPerOctave;  a = obj.Anchor;
    f = double(f);
    xi = V * log2(max(f, realmin) / a);
    switch obj.Kind_{m}
        case 'band'
            G = sqrt(2) * i_itersine(xi - obj.K_(m));
        case 'lowpass'
            K = obj.K_(2);
            S = 2 * i_itersine(xi - K).^2;                 % only the first band member reaches below its centre
            G = sqrt(max(2 - S, 0)) .* (xi < K & f > 0);
            G(f == 0) = 0;
        case 'highpass'
            K = obj.K_(end-1);
            S = 2 * i_itersine(xi - K).^2;
            G = sqrt(max(2 - S, 0)) .* (xi > K);
            G(f >= obj.SamplingFrequency/2) = 0;            % the Nyquist bin is not analytic
    end
end

function w = i_itersine(u)
    w = zeros(size(u));
    in = abs(u) <= 1;
    w(in) = sin(0.5*pi * cos(0.5*pi*u(in)).^2);
end
% Author: Diellor Basha, 2026
