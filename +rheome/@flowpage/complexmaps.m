function Z = complexmaps(obj, iFrame)
% COMPLEXMAPS  The complex scale-filtered synthesis at one frame, all scales, one GEMM.
%
%   Z = complexmaps(fp, iFrame)      -> [V x (nScale+1)] complex, col 1 = all scales
%
% ⭐ THE SINGLE SOURCE FOR EVERY PER-FRAME CORTICAL PRODUCT. Signed vorticity is real(Z),
% the envelope is abs(Z), and the dominant-scale map is argmax over abs(Z(:,2:end)).^2 --
% all three from ONE [V x Ks]*[Ks x nScale+1] product. Computing them separately means
% streaming the 65 MB basis through cache two or three times per displayed frame, which is
% the whole cost of a redraw: it took the browser from 39 ms to 74 ms before this existed.
%
% See also: maps, map, dominantScale
%
% Author: Diellor Basha, 2026

    if ~isscalar(iFrame) || iFrame < 1 || iFrame > obj.NumFrames || mod(iFrame,1) ~= 0
        error('flowpage:frame', 'Frame must be an integer in 1..%d.', obj.NumFrames);
    end
    c = obj.Coefficients(:, iFrame);
    G = [sum(obj.ScaleGains, 2), obj.ScaleGains];        % [Ks x (nScale+1)]
    Z = obj.Bundle.Phi * (cast(G, 'like', c) .* c);
end

% Author: Diellor Basha, 2026
