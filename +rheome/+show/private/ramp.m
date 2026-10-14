function C = ramp(n)
% RAMP  n colours for an ORDERED series (members fine -> coarse, modes low -> high) from the house
% sequential map (atlas-tile-code), stopping short of its light end so every line reads on paper.
%
% Author: Diellor Basha, 2026
    tokens();
    M = nxr_cmap('atlas-tile-code');
    C = M(round(linspace(1, 0.8*size(M,1), max(n,1))), :);
end
