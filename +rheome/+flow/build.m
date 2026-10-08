function K = build(ctx, opts)
% FLOW.BUILD  Assemble all fused flow kernels + a coeff-vs-vertex agreement report.
%   K = rheome.flow.build(ctx [,opts])
%
% Returns a struct of every kernel (each from its dedicated rheome.flow.* builder) plus a report:
%   K.field K.divergence K.curl K.potential K.stream K.energy K.enstrophy K.helicity
%   K.report .Ks .C .V .phiCorr .psiCorr   (coeff-vs-vertex correlation of the potentials)
%
% Apply linear kernels as  q(:,t) = K.<name>.vertexOperator * ctx.F(:,t)  (potentials: coeff via
% K.potential.scalarModes * (K.potential.coeffOperator * F), or full-rank K.potential.vertexOperator * F).
% Grams as  scalar(t) = ctx.F(:,t)' * K.<name>.gram * ctx.F(:,t).
%
% See also: rheome.flow.context, rheome.flow.field, rheome.flow.divergence, rheome.flow.curl, rheome.flow.potential, rheome.flow.stream
%
% Author: Diellor Basha, 2026

    if nargin < 2, opts = struct(); end
    K.field      = rheome.flow.field(ctx);
    K.divergence = rheome.flow.divergence(ctx);
    K.curl       = rheome.flow.curl(ctx);
    K.potential  = rheome.flow.potential(ctx, struct('vertex', true));
    K.stream     = rheome.flow.stream(ctx, struct('vertex', true));
    K.energy     = rheome.flow.energy(ctx);
    K.enstrophy  = rheome.flow.enstrophy(ctx);
    K.helicity   = rheome.flow.helicity(ctx);

    phiC = K.potential.scalarModes * K.potential.coeffOperator;
    psiC = K.stream.scalarModes    * K.stream.coeffOperator;
    K.report = struct('Ks', ctx.Ks, 'C', ctx.C, 'V', ctx.V, ...
        'phiCorr', corr(phiC(:), K.potential.vertexOperator(:)), ...
        'psiCorr', corr(psiC(:), K.stream.vertexOperator(:)));
end

% Author: Diellor Basha, 2026
