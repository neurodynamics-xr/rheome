function ks = stream(ctx, opts)
% FLOW.STREAM  Fused kernel: sensors -> Helmholtz stream function Psi (vortices).
%   ks = rheome.flow.stream(ctx [,opts])
%
% Psi solves  K*Psi = vortw(J)  (vortices = extrema of Psi). Like rheome.flow.potential, the
% COEFFICIENT kernel is diagonal (1/lambda) in the scalar-LBO eigenbasis:
%   Psi_coeff = Lambda^-1 .* (scalarModes' * (Brot * currentKernel))     [Ks x C]
% With opts.vertex (default true) the full-rank vertex kernel (rheome.differential.helmholtz Psi)
% is also returned.
%
% OUTPUT (struct ks): .coeffOperator [Ks x C], .scalarModes [V x Ks], .Lambda [Ks x 1],
%   .vertexOperator [V x C] (if opts.vertex), .domain='coeff', .out, .chNames, .provenance.
%
% See also: rheome.flow.weak, rheome.flow.potential, rheome.differential.helmholtz
%
% Author: Diellor Basha, 2026

    if nargin < 2, opts = struct(); end
    if ~isfield(opts,'vertex') || isempty(opts.vertex), opts.vertex = true; end

    wk  = rheome.flow.weak(ctx);
    Phi = ctx.lbo.Phi(:, 1:ctx.Ks);                 % [V x Ks]
    lam = ctx.lbo.Lambda(1:ctx.Ks);                 % [Ks x 1]
    tol = 1e-6 * max(ctx.lbo.Lambda);
    invLam = 1 ./ lam;  invLam(lam < tol) = 0;      % drop the per-hemisphere near-null constants

    coeffOperator = invLam .* (Phi' * (wk.Brot * ctx.currentKernel));   % [Ks x C]

    ks = struct('coeffOperator', coeffOperator, 'scalarModes', Phi, 'Lambda', lam, ...
        'domain', 'coeff', 'out', 'stream function Psi [vortices]', ...
        'chNames', {ctx.chNames}, ...
        'provenance', struct('inverse', ctx.Method, 'analysisBasis','laplace-beltrami', ...
            'Ks',ctx.Ks,'op','rheome.flow.weak.Brot'));

    if opts.vertex
        H = rheome.differential.helmholtz(ctx.currentKernel, ctx.S);
        ks.vertexOperator = H.Psi;                   % [V x C] full-rank, pre-solved
    end
end

% Author: Diellor Basha, 2026
