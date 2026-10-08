function E = scaleEnergy(obj, iFrame)
% SCALEENERGY  Energy per spatial scale at one frame -- the globality readout.
%
%   E = scaleEnergy(fp, iFrame)      -> [1 x nScale]
%
% ⭐ THIS IS THE GLOBAL/LOCAL AXIS. Energy in LOW lambda is globally distributed structure;
% energy in HIGH lambda is spatially confined structure. A feature is not merely "a vortex"
% but a vortex occurring while a particular part of the spatial spectrum is occupied, and
% that occupancy is what this returns.
%
% Read in the coefficient domain -- sum_k |h_g(lambda_k) c_k|^2 -- so no vertex field is
% formed. Phi is orthonormal in the mass inner product, so this equals the field energy of
% the scale-filtered map.
%
% See also: centroid, globalIndex, modeSpectrum
%
% Author: Diellor Basha, 2026

    if ~isscalar(iFrame) || iFrame < 1 || iFrame > obj.NumFrames || mod(iFrame,1) ~= 0
        error('flowpage:frame', 'Frame must be an integer in 1..%d.', obj.NumFrames);
    end
    c = double(obj.Coefficients(:, iFrame));
    E = sum(abs(obj.ScaleGains .* c).^2, 1);
end

% Author: Diellor Basha, 2026
