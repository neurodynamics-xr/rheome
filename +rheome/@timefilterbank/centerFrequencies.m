function fc = centerFrequencies(obj)
% CENTERFREQUENCIES  Centre of every member, Hz; edge members report the centroid of their support.
%
%   fc = centerFrequencies(tfb)    [1 x M], ascending
%
% Author: Diellor Basha, 2026

    fc = obj.Anchor * 2.^(obj.K_ / obj.VoicesPerOctave);
    f = obj.Freq_;
    for m = [1, obj.NumMembers]
        b = obj.Bins_(m, :);  g2 = gains(obj, m, f(b(1):b(2))).^2;
        fc(m) = sum(f(b(1):b(2)) .* g2) / sum(g2);
    end
end
% Author: Diellor Basha, 2026
