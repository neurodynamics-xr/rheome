function M = maps(obj, iFrame, mode)
% MAPS  Every spatial scale at one frame, in ONE GEMM.
%
%   M = maps(fp, iFrame)              -> [V x (nScale+1)], column 1 = all scales summed
%   M = maps(fp, iFrame, 'magnitude')
%
% ⭐ WHY THIS EXISTS RATHER THAN LOOPING map(). Calling map() per scale is nScale separate
% [V x Ks]*[Ks x 1] GEMVs, each streaming the whole 65 MB Phi through cache for one column
% of output. Batching them is Phi * (H .* c), a single [V x Ks]*[Ks x nScale] GEMM that
% streams Phi ONCE. Same arithmetic, measurably faster, and it is what the browser calls on
% every frame.
%
% Both readouts come from complexmaps, so asking for signed and magnitude at the same frame
% costs one product, not two.
%
% See also: complexmaps, map, rheome.flowbrowser/drawCortex
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(mode), mode = 'signed'; end
    if ~isscalar(iFrame) || iFrame < 1 || iFrame > obj.NumFrames || mod(iFrame,1) ~= 0
        error('flowpage:frame', 'Frame must be an integer in 1..%d.', obj.NumFrames);
    end

    Z = complexmaps(obj, iFrame);                        % one GEMM

    switch lower(mode)
        case 'signed',    M = real(Z);
        case 'magnitude', M = abs(Z);
        otherwise
            error('flowpage:mode', 'mode must be ''signed'' or ''magnitude'', got ''%s''.', mode);
    end
end

% Author: Diellor Basha, 2026
