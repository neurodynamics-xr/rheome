function kc = curl(ctx, opts)
% FLOW.CURL  Fused flow operator: sensors -> curl (scalar vorticity).
%   kc = rheome.flow.curl(ctx)
%   kc = rheome.flow.curl(ctx, struct('method','weak'))
%
% CANONICAL (eigenbasis): the Laplace-Beltrami coefficients of the vorticity, straight from
% sensors -- no per-vertex field materialized:
%   g_curl(:,t) = kc.coeffOperator * ctx.F(:,t)          coeffOperator [Ks x C]  (= Phi_LB' M curl reconstruct)
% FALLBACK (vertex, for maps/validation only):
%   vort(:,t)   = kc.vertexOperator * ctx.F(:,t)          vertexOperator [V x C]
%   or reconstruct from coeffs: kc.scalarModes * (kc.coeffOperator * F).
%
% Exact by linearity: curl(reconstruct(diracInverse*F)) = rheome.differential.curl(currentKernel,S,fg) * F, and its
% Laplace-Beltrami coefficients are Phi_LB' * M * (that vertex kernel). Positive = CCW (from
% outside), negative = CW.
%
% ---------------------------------------------------------------------------------------------
% TWO CONSTRUCTIONS OF coeffOperator -- opts.method
%
%   'lumped' (DEFAULT, and what the pipeline uses)
%       Phi' * M * (fg.W * curl_face)
%       rheome.differential.curl computes vorticity per FACE, then area-averages it to VERTICES (fg.W),
%       then the projection integrates it. That per-vertex step is what a surface MAP needs.
%
%   'weak'   (option 2, for comparison -- NOT the pipeline default)
%       Phi' * wd.Curl        with wd = rheome.operators.weak_differential
%       Integrates by parts, int psi (curl J) = -int (N x grad psi) . J, so the derivative sits on
%       the smooth basis function and the face->vertex average never happens. No mass matrix: the
%       pairing already integrates.
%
% ⚠ THESE ARE THE SAME OPERATOR up to that averaging step -- the weak form agrees with the strong
% form integrated per face to 7.2e-16 (curl_weakform.m). Integration by parts is EXACT here: P1
% fields, exact per-face quadrature, closed surface, no boundary term. So 'weak' is not a different
% discretisation, and 'lumped' is not wrong; the difference is one redundant smoothing step.
%
% MEASURED against the analytic curl(zhat x p) = 2z/R on an ico5 sphere:
%       'lumped'  rel err 4.53e-4
%       'weak'    rel err 1.76e-4      <- 2.6x closer
% but at 1-30% noise on J the weak form is a flat 5% WORSE, because the lumping smooths noise while
% it blurs signal. It is a trade. The default is unchanged so that no existing number moves; run
% flow_curl_methods.m to see both on real data before switching.
% ---------------------------------------------------------------------------------------------
%
% See also: rheome.flow.field, rheome.flow.divergence, rheome.flow.stream, rheome.differential.curl,
%           rheome.operators.weak_differential, flow_curl_methods, OPERATOR-REGISTRY.md
%
% Author: Diellor Basha, 2026

    if nargin < 2, opts = struct(); end
    if ~isfield(opts,'method') || isempty(opts.method), opts.method = 'lumped'; end

    % The per-vertex map is ALWAYS the lumped route: a surface render needs a value per vertex, and
    % producing one is exactly what the face->vertex average is for. Only the ANALYSIS path
    % (coeffOperator) is affected by the method.
    vertexOperator = rheome.differential.curl(ctx.currentKernel, ctx.S, ctx.fg);             % [V x C]

    switch lower(opts.method)
        case 'lumped'
            coeffOperator = ctx.lbo.Phi' * (ctx.lbo.Mass * vertexOperator);           % [Ks x C]
        case 'weak'
            wd = rheome.operators.weak_differential(ctx.S.Vertices, ctx.S.Faces, ctx.fg);
            coeffOperator = ctx.lbo.Phi' * (wd.Curl * ctx.currentKernel);             % [Ks x C]
        otherwise
            error('flow:curl:method', 'method must be ''lumped'' or ''weak'', got ''%s''.', opts.method);
    end

    kc = struct('coeffOperator', coeffOperator, 'scalarModes', ctx.lbo.Phi, 'vertexOperator', vertexOperator, ...
        'domain', 'eigenbasis', 'out', 'vorticity [CCW+/CW-]', 'method', lower(opts.method), ...
        'chNames', {ctx.chNames}, ...
        'provenance', struct('inverse', ctx.Method, 'analysisBasis', 'laplace-beltrami', ...
            'P',ctx.P, 'Ks',ctx.Ks, 'op','rheome.differential.curl', 'method',lower(opts.method)));
end

% Author: Diellor Basha, 2026
