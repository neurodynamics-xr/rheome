function kf = field(ctx)
% FLOW.FIELD  Root fused kernel: sensors -> per-vertex current 3-vectors.
%   kf = rheome.flow.field(ctx)   kf.vertexOperator [3V x C];  J(:,t) = kf.vertexOperator * ctx.F(:,t)
%
% Whatever source estimate rheome.flow.context built -- a plain whitened minimum norm by default, or
% the Dirac-mode inverse if it was asked for. In AMPLITUDE measure this IS the imaging kernel
% (physical current, A.m), and it is the root every other flow kernel is built from.
% ctx.Method says which route produced it; the provenance carries it forward.
%
% See also: rheome.flow.context, rheome.flow.divergence, rheome.flow.curl
%
% Author: Diellor Basha, 2026

    kf = struct('vertexOperator', ctx.currentKernel, 'domain', 'vertex', 'out', 'current J [A.m]', ...
        'chNames', {ctx.chNames}, ...
        'provenance', struct('inverse', ctx.Method, 'P', ctx.P, 'measure', 'amplitude'));
end

% Author: Diellor Basha, 2026
