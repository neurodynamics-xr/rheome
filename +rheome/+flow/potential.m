function kp = potential(ctx, opts)
% FLOW.POTENTIAL  Fused kernel: sensors -> Helmholtz scalar potential Phi (sources/sinks).
%   kp = rheome.flow.potential(ctx [,opts])
%
% Phi solves the Poisson problem  K*Phi = divw(J)  (sources = maxima, sinks = minima). In the
% scalar-LBO eigenbasis K^+ is diagonal (1/lambda), so the COEFFICIENT kernel is a single
% matrix, band-limited to Ks modes:
%   Phi_coeff = Lambda^-1 .* (scalarModes' * (Bdiv * currentKernel))     [Ks x C]   (near-null modes zeroed)
% With opts.vertex (default true) the full-rank VERTEX kernel is also returned, pre-solved for
% all C channel columns via rheome.differential.helmholtz (still one matvec at runtime).
%
% OUTPUT (struct kp): .coeffOperator [Ks x C], .scalarModes [V x Ks], .Lambda [Ks x 1],
%   .vertexOperator [V x C] (if opts.vertex), .domain='coeff', .out, .chNames, .provenance.
%   Vertex field for display:  Phi(:,t) = kp.scalarModes * (kp.coeffOperator * ctx.F(:,t))   (coeff)
%                       or      Phi(:,t) = kp.vertexOperator * ctx.F(:,t)              (full-rank)
%
% See also: rheome.flow.weak, rheome.flow.stream, rheome.differential.helmholtz
%
% Author: Diellor Basha, 2026

    if nargin < 2, opts = struct(); end
    if ~isfield(opts,'vertex') || isempty(opts.vertex), opts.vertex = true; end

    wk  = rheome.flow.weak(ctx);
    Phi = ctx.lbo.Phi(:, 1:ctx.Ks);                 % [V x Ks]
    lam = ctx.lbo.Lambda(1:ctx.Ks);                 % [Ks x 1]
    tol = 1e-6 * max(ctx.lbo.Lambda);
    invLam = 1 ./ lam;  invLam(lam < tol) = 0;      % drop the per-hemisphere near-null constants

    coeffOperator = invLam .* (Phi' * (wk.Bdiv * ctx.currentKernel));    % [Ks x C]

    kp = struct('coeffOperator', coeffOperator, 'scalarModes', Phi, 'Lambda', lam, ...
        'domain', 'coeff', 'out', 'scalar potential Phi [sources/sinks]', ...
        'chNames', {ctx.chNames}, ...
        'provenance', struct('inverse', ctx.Method, 'analysisBasis','laplace-beltrami', ...
            'Ks',ctx.Ks,'op','rheome.flow.weak.Bdiv'));

    if opts.vertex
        H = rheome.differential.helmholtz(ctx.currentKernel, ctx.S);
        kp.vertexOperator = H.Phi;                   % [V x C] full-rank, pre-solved
    end
end

% Author: Diellor Basha, 2026
