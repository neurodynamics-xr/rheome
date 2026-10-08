function kd = divergence(ctx)
% FLOW.DIVERGENCE  Fused flow operator: sensors -> divergence (sources/sinks).
%   kd = rheome.flow.divergence(ctx)
%
% CANONICAL (eigenbasis): the Laplace-Beltrami coefficients of the divergence, straight from
% sensors -- no per-vertex field materialized:
%   g_div(:,t) = kd.coeffOperator * ctx.F(:,t)          coeffOperator [Ks x C]  (= Phi_LB' M div reconstruct)
% FALLBACK (vertex, for maps/validation only):
%   div(:,t)   = kd.vertexOperator * ctx.F(:,t)          vertexOperator [V x C]
%   or reconstruct from coeffs: kd.scalarModes * (kd.coeffOperator * F).
%
% Exact by linearity: div(reconstruct(diracInverse*F)) = rheome.differential.divergence(currentKernel,S,fg) * F, and its
% Laplace-Beltrami coefficients are Phi_LB' * M * (that vertex kernel). Positive = source, negative = sink.
%
% See also: rheome.flow.field, rheome.flow.curl, rheome.flow.potential, rheome.differential.divergence, OPERATOR-REGISTRY.md
%
% Author: Diellor Basha, 2026

    vertexOperator = rheome.differential.divergence(ctx.currentKernel, ctx.S, ctx.fg);        % [V x C]  (fallback)
    coeffOperator  = ctx.lbo.Phi' * (ctx.lbo.Mass * vertexOperator);                   % [Ks x C] Laplace-Beltrami coeffs (canonical)
    kd = struct('coeffOperator', coeffOperator, 'scalarModes', ctx.lbo.Phi, 'vertexOperator', vertexOperator, ...
        'domain', 'eigenbasis', 'out', 'divergence [sources+/sinks-]', ...
        'chNames', {ctx.chNames}, ...
        'provenance', struct('inverse', ctx.Method, 'analysisBasis','laplace-beltrami', ...
            'P',ctx.P, 'Ks',ctx.Ks, 'op','rheome.differential.divergence'));
end

% Author: Diellor Basha, 2026
