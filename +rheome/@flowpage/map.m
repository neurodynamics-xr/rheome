function m = map(obj, iFrame, iScale, mode)
% MAP  The cortical vorticity map at one frame and one spatial scale.
%
%   m = map(fp, iFrame, iScale)              % signed vorticity, [V x 1]
%   m = map(fp, iFrame, iScale, 'magnitude')
%   m = map(fp, iFrame, 0)                   % summed over the whole scale bank
%
% ONE GEMV: Phi * (h_g .* C(:,t)). This is the interactive path -- ~33 MFLOP, so scrubbing
% time or switching scale is instant once the page exists.
%
% 'signed' is real(.) -- positive = counter-clockwise seen from outside, negative = clockwise.
% Both handednesses live in ONE signed map at every instant; there are not two maps.
% 'magnitude' is abs(.), which is steadier to watch because it does not flicker with the
% carrier phase, but it discards the handedness that makes a vortex PAIR visible.
%
% See also: scaleEnergy, centroid, rheome.flowbrowser
%
% Author: Diellor Basha, 2026

    if nargin < 4 || isempty(mode), mode = 'signed'; end
    if ~isscalar(iFrame) || iFrame < 1 || iFrame > obj.NumFrames || mod(iFrame,1) ~= 0
        error('flowpage:frame', 'Frame must be an integer in 1..%d, got %s.', ...
            obj.NumFrames, mat2str(iFrame));
    end
    if ~isscalar(iScale) || iScale < 0 || iScale > obj.NumScales || mod(iScale,1) ~= 0
        error('flowpage:scale', 'Scale must be an integer in 0..%d, got %s.', ...
            obj.NumScales, mat2str(iScale));
    end

    c = obj.Coefficients(:, iFrame);
    if iScale == 0
        g = sum(obj.ScaleGains, 2);          % the bank's total response
    else
        g = obj.ScaleGains(:, iScale);
    end
    z = obj.Bundle.Phi * (cast(g, 'like', c) .* c);

    switch lower(mode)
        case 'signed',    m = real(z);
        case 'magnitude', m = abs(z);
        otherwise
            error('flowpage:mode', 'mode must be ''signed'' or ''magnitude'', got ''%s''.', mode);
    end
end

% Author: Diellor Basha, 2026
