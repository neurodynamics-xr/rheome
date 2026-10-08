function [H, f] = freqz(obj)
% FREQZ  One-sided gains of every member on the record's DFT grid.
%
%   [H, f] = freqz(tfb)    H [M x nBins], f [1 x nBins] Hz
%
% Same convention as cwtfilterbank/freqz: one-sided, real, analytic, peak gain sqrt(2)
% per member and S = sum |H|^2 = 2 on (0, fs/2).
%
% Author: Diellor Basha, 2026

    f = obj.Freq_;
    M = obj.NumMembers;
    H = zeros(M, numel(f));
    for m = 1:M
        b = obj.Bins_(m, :);
        H(m, b(1):b(2)) = gains(obj, m, f(b(1):b(2)));
    end
end
% Author: Diellor Basha, 2026
