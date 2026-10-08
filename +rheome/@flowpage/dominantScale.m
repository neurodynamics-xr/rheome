function [iDom, total] = dominantScale(obj, iFrame, Z)
% DOMINANTSCALE  Which spatial scale dominates at each VERTEX, and the total energy there.
%
%   [iDom, total] = dominantScale(fp, iFrame)      -> both [V x 1]
%   [iDom, total] = dominantScale(fp, iFrame, Z)   % reuse a complexmaps result
%
% ⭐ WHAT THE PER-SCALE PANELS CANNOT SHOW. Separate surfaces separate exactly what you want
% to compare: to ask "is the flow coarse here and fine there" you must hold the cortex fixed
% and let scale vary WITHIN it. This returns, per vertex, the argmax over scales -- rendered
% as hue, with `total` as brightness, it is one picture of where the flow is coarse and where
% it is fine.
%
% ⚠ ARGMAX IS UNSTABLE WHERE THE SPECTRUM IS FLAT. Two scales within a few percent of each
% other will alternate frame to frame and flicker. `total` is returned so the caller can
% darken low-energy vertices, where the argmax means least.
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(Z), Z = complexmaps(obj, iFrame); end
    P = abs(double(Z(:, 2:end))).^2;               % per-scale energy per vertex
    [~, iDom] = max(P, [], 2);
    total = sum(P, 2);
end

% Author: Diellor Basha, 2026
