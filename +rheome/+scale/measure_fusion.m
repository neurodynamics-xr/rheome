function [T, X] = measure_fusion(name, ctx, K, opts)
% SCALE.MEASURE_FUSION  Fused flow kernels against reconstruct-then-differentiate, on one subject's data.
%
%   [T, X] = rheome.scale.measure_fusion(name)
%   [T, X] = rheome.scale.measure_fusion(name, ctx, rheome.flow.build(ctx), NumFrames=20)
%
% fusion_exactness_omega.m (MS1 Section 8.2) as a per-subject spot check (group plan G14): on NumFrames
% random frames of the cached recording (seed 0), at the default minimum-norm context:
%   divVertex, curlVertex   max|K.(div|curl).vertexOperator*F - differential.(divergence|curl)(currentKernel*F)| / max|staged|
%   divCoeff, curlCoeff     the same for coeffOperator against Phi_LB' M (staged vertex map)
%   helmholtzVertex         max over Phi, Psi of the vertex potentials against differential.helmholtz
%   energy, enstrophy, helicity   |F' Q F - sum(map(F))| / |sum(map(F))|
%   phiCorr, psiCorr        corr of the Ks-mode potentials with the full-rank vertex potentials
%                           (truncation, NOT exactness: it measures how much lies inside the modes)
% Rows (analysis "fusion"): <quantity>_median and <quantity>_worst (max; min for the correlations),
% max_exact_error (the worst of every exactness quantity), Ks, build_s (NaN when K is passed in). X is the per-frame table.
%
% ⚠ For scale, a reference subject (Ks 800): div 9.0e-15 (worst 1.9e-14), curl 2.2e-15, potentials
% 1.3e-13 (2.6e-13), Grams 0.7-3.1e-15; Phi/Psi 0.987/0.998. Expected <= 1e-12 everywhere.
% ⚠ helicity is ~1e-20 in SI units, far below eps: never guard its denominator with eps.
%
% See also: rheome.flow.context, rheome.flow.build
%
% Author: Diellor Basha, 2026

    arguments
        name (1,:) char
        ctx = []
        K = []
        opts.NumFrames (1,1) double {mustBeInteger, mustBePositive} = 20
        opts.Seed (1,1) double = 0
    end
    if isempty(ctx), ctx = rheome.flow.context(name); end
    tBuild = NaN;                                                % NaN: the kernels were built by the caller
    if isempty(K), tb = tic;  K = rheome.flow.build(ctx);  tBuild = toc(tb); end
    rs = rng(opts.Seed);  c = onCleanup(@() rng(rs));
    nT = size(ctx.F, 2);  frames = sort(randperm(nT, min(opts.NumFrames, nT)));
    Mm = ctx.lbo.Mass;  Phi = ctx.lbo.Phi;
    relmax = @(a, b) max(abs(a - b)) / max(abs(b));
    rec = zeros(numel(frames), 11);
    for i = 1:numel(frames)
        Ft = double(ctx.F(:, frames(i)));
        J  = ctx.currentKernel * Ft;
        dS = rheome.differential.divergence(J, ctx.S, ctx.fg);
        cS = rheome.differential.curl(J, ctx.S, ctx.fg);
        H  = rheome.differential.helmholtz(J, ctx.S);
        pF = K.potential.scalarModes * (K.potential.coeffOperator * Ft);
        sF = K.stream.scalarModes    * (K.stream.coeffOperator    * Ft);
        pV = K.potential.vertexOperator * Ft;  sV = K.stream.vertexOperator * Ft;
        eF = Ft' * K.energy.gram * Ft;     eS = sum(K.energy.map(Ft));
        oF = Ft' * K.enstrophy.gram * Ft;  oS = sum(K.enstrophy.map(Ft));
        hF = Ft' * K.helicity.gram * Ft;   hS = sum(K.helicity.map(Ft));
        rec(i, :) = [frames(i), relmax(K.divergence.vertexOperator * Ft, dS), relmax(K.curl.vertexOperator * Ft, cS), ...
            relmax(K.divergence.coeffOperator * Ft, Phi' * (Mm * dS)), relmax(K.curl.coeffOperator * Ft, Phi' * (Mm * cS)), ...
            corr(pF, pV), corr(sF, sV), max(relmax(pV, H.Phi), relmax(sV, H.Psi)), ...
            abs(eF - eS)/abs(eS), abs(oF - oS)/abs(oS), abs(hF - hS)/abs(hS)];
    end
    X = array2table(rec, 'VariableNames', {'frame','divVertex','curlVertex','divCoeff','curlCoeff', ...
        'phiCorr','psiCorr','helmholtzVertex','energy','enstrophy','helicity'});
    q = string(X.Properties.VariableNames(2:end));  m = strings(0,1);  v = [];  u = strings(0,1);
    for k = q
        x = X.(k);  isC = endsWith(k, "Corr");
        w = max(x);  un = "relative";
        if isC, w = min(x);  un = "r"; end
        m = [m; k + "_median"; k + "_worst"];  v = [v; median(x); w];  u = [u; un; un]; %#ok<AGROW>
    end
    ex = X{:, setdiff(q, ["phiCorr" "psiCorr"], 'stable')};
    T = rheome.scale.rows("fusion", [m; "max_exact_error"; "Ks"; "build_s"], [v; max(ex(:)); ctx.Ks; tBuild], ...
                          [u; "relative"; "modes"; "s"]);
end

% Author: Diellor Basha, 2026
