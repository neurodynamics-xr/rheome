function readout = flow_readout(velocityField, surface)
% DYNAMICS.FLOW_READOUT  GCF readouts of a propagation velocity field.
%
%   readout = rheome.dynamics.flow_readout(velocityField, surface)
%
% Applies the (instantaneous) flow operators to the PROPAGATION velocity field v (from
% dynamics.opticalflow_*), giving the Global Cortical-Flow descriptors:
%   .divergence  [nV x T]  div(v) -- propagation source-sink motifs (where activity emerges /
%                          converges as it propagates)
%   .vorticity   [nV x T]  curl(v) -- circulation of the propagation
%   .speed       [nV x T]  |v| per vertex (propagation speed, per-frame units; x sfreq for /s)
%
% Reuses rheome.differential.divergence / rheome.differential.curl on the velocity field.
%
% See also: rheome.dynamics.opticalflow_scalar, rheome.dynamics.opticalflow_vector, rheome.differential.divergence
%
% Author: Diellor Basha, 2026

    readout = struct();
    readout.divergence = rheome.differential.divergence(velocityField, surface);
    readout.vorticity  = rheome.differential.curl(velocityField, surface);
    vx = velocityField(1:3:end, :);  vy = velocityField(2:3:end, :);  vz = velocityField(3:3:end, :);
    readout.speed = sqrt(vx.^2 + vy.^2 + vz.^2);
end

% Author: Diellor Basha, 2026
