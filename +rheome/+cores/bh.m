function q = bh(p)
% CORES.BH  Benjamini-Hochberg q values (nxr bh).
%
%   q = rheome.cores.bh(p)
%
% q = p(sorted) * m / rank, reverse cumulative minimum, clipped to 1, scattered back; m counts every p given, so remove NaN first. Tied p get
% the same q. Ports nxr-analysis src/nxr/vortex_cores.py bh (line 96, commit 0db82e4).
%
% Author: Diellor Basha, 2026
    p = p(:);  m = numel(p);
    [ps, o] = sort(p);
    qs = ps * m ./ (1:m).';
    qs = flipud(cummin(flipud(qs)));
    q = zeros(m, 1);  q(o) = min(qs, 1);
end
