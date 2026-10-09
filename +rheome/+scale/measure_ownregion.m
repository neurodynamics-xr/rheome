function [T, X] = measure_ownregion(name, S, K, opts)
% SCALE.MEASURE_OWNREGION  Readout rule 5 for one subject: which Desikan-Killiany regions own their flux.
%
%   [T, X] = rheome.scale.measure_ownregion(name)
%   [T, X] = rheome.scale.measure_ownregion(name, rheome.scale.sensors(name), rheome.flow.build(ctx))
%   [T, X] = rheome.scale.measure_ownregion(name, S, K, Atlas="Desikan-Killiany", Threshold=0.25)
%
% rheome.flow.crosstalk on the fused divergence and curl vertex kernels (rheome.flow.build on the default
% minimum-norm context) against the subject's own leadfield S.G: per region, the fraction of the
% region's measured flux power that originates inside it. MS1 rule 5 reports a region's flux only where
% this reaches 25% (group plan G3, rule 5).
% Rows (analysis "ownregion"), per kernel (band "curl" | "div"):
%   own_median, own_q25, own_q75   over regions                                          fraction
%   frac_regions_pass              fraction of regions with own >= Threshold ⭐ the rule's number
%   n_regions_pass, n_regions      regions passing, regions with any vertex
%   area_frac_pass                 share of the atlas's cortical area in passing regions  fraction
% X: one row per region -- label, n_vertices, area_m2, own_curl, own_div, dominant_curl, dominant_div.
%
% ⚠ For scale, a reference subject (270 CTF channels): median own fraction 6%, 17 of 68 regions >= 25%,
% covering half the cortical area (rheome.flow.crosstalk).
% ⚠ A MASK, NOT A DIVISOR (rheome.flow.crosstalk): never divide a reading by sensitivity.
%
% See also: rheome.flow.crosstalk, rheome.flow.build, rheome.load.atlas
%
% Author: Diellor Basha, 2026

    arguments
        name (1,:) char
        S = []
        K = []
        opts.Atlas (1,:) char = 'Desikan-Killiany'
        opts.Threshold (1,1) double = 0.25
    end
    if isempty(S), S = rheome.scale.sensors(name); end
    if isempty(K), K = rheome.flow.build(rheome.flow.context(name)); end
    A = rheome.load.atlas(name, opts.Atlas);  Pm = A.Membership;
    a = zeros(size(Pm, 2), 1);                                 % lumped vertex area, from the cached bases
    for h = intersect(fieldnames(S.B), {'L','R'})', a(S.B.(h{1}).gv) = full(sum(S.B.(h{1}).lbo.Mass, 2)); end
    [ownC, domC] = rheome.flow.crosstalk(K.curl.vertexOperator, S.G, Pm);
    [ownD, domD] = rheome.flow.crosstalk(K.divergence.vertexOperator, S.G, Pm);
    X = table(string(A.Label(:)), full(sum(Pm, 2)), double(Pm) * a, ownC, ownD, domC, domD, ...
              'VariableNames', {'label','n_vertices','area_m2','own_curl','own_div','dominant_curl','dominant_div'});
    T = rheome.scale.rows("ownregion", strings(0,1), [], "");
    for k = ["curl" "div"]
        o = X.("own_" + k);  ok = ~isnan(o);  pass = ok & o >= opts.Threshold;
        T = [T; rheome.scale.rows("ownregion", ["own_median" "own_q25" "own_q75" "frac_regions_pass" ...
              "n_regions_pass" "n_regions" "area_frac_pass"], ...
              [median(o(ok)) prctile(o(ok), 25) prctile(o(ok), 75) mean(pass(ok)) nnz(pass) nnz(ok) ...
               sum(X.area_m2(pass)) / sum(X.area_m2(ok))], ...
              ["fraction" "fraction" "fraction" "fraction" "regions" "regions" "fraction"], k)]; %#ok<AGROW>
    end
end

% Author: Diellor Basha, 2026
