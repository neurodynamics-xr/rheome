function S = modeSpectrum(obj, iFrame)
% MODESPECTRUM  Energy per LBO mode at one frame -- the raw lambda spectrum.
%
%   S = modeSpectrum(fp, iFrame)     -> [Ks x 1], aligned with fp.Lambda
%
% The unbinned version of scaleEnergy: what the eigenspectrum panel plots. scaleEnergy is
% this seen through the graph bank's smooth windows.
%
% Author: Diellor Basha, 2026

    if ~isscalar(iFrame) || iFrame < 1 || iFrame > obj.NumFrames || mod(iFrame,1) ~= 0
        error('flowpage:frame', 'Frame must be an integer in 1..%d.', obj.NumFrames);
    end
    S = abs(double(obj.Coefficients(:, iFrame))).^2;
end

% Author: Diellor Basha, 2026
