function [T, X] = measure_helmholtzbands(name, opts)
% SCALE.MEASURE_HELMHOLTZBANDS  Helmholtz-band recovery on the subject's own cortex (geometry only; no MEG).
%
%   [T, X] = rheome.scale.measure_helmholtzbands(name)
%   [T, X] = rheome.scale.measure_helmholtzbands(name, Plants=5, Hemis=["L" "R"], Seed=11)
%
% Part 1 of the Helmholtz-band test (MS1 Fig. 4A-B) per hemisphere, for the group (plan G9). For every
% interior band m0 of the tight Laplace-Beltrami bank (logitersine, 2 voices per octave; 315-23 mm on a
% 1000-mode reference cortex) and Plants random vertices: a source field grad(psi) and a vortex field
% N x grad(psi), psi the bank's own member m0 at the vertex, read by rheome.differential.helmholtzbands.
% Per plant (X, one row per hemisphere x band x plant):
%   bandSource/bandVortex   argmax band of the irrotational / solenoidal energy (recovered = equals m0)
%   leakSource_pct          solenoidal / irrotational energy of the source, %       (wrong part / right part)
%   leakVortex_pct          irrotational / solenoidal energy of the vortex, %
%   harmSource_pct, harmVortex_pct   harmonic residual share of the field energy, %
%   captureSource, captureVortex     projected / unprojected Dirichlet energy (K modes hold it)
% Summary rows per band ("<wl> mm"), pooled over hemispheres and plants: recovered_source,
% recovered_vortex (fractions), leak_source_pct, leak_vortex_pct, harm_source_pct, harm_vortex_pct
% (medians), n_plants. And per hemisphere (band ""): recovered_all_<H> (fraction of its 2 x bands x
% Plants cases recovered).
%
% ⚠ For scale, a reference subject (left hemisphere, 5 plants): 45/45 sources and 45/45 vortices
% recovered, leak 0.4-0.8% (band medians), harmonic 4-9% to 46 mm, 9-11% at 32 mm, 14-18% at 23 mm.
% ⚠ The baseline is closed form on the sphere (tHelmholtzBands); this is the folded-cortex check.
%
% See also: rheome.differential.helmholtzbands, rheome.differential.helmholtz, rheome.graphfilterbank
%
% Author: Diellor Basha, 2026

    arguments
        name (1,:) char
        opts.Plants (1,1) double {mustBeInteger, mustBePositive} = 5
        opts.Hemis string = ["L" "R"]
        opts.Seed (1,1) double = 11
    end
    B = rheome.load.bases(name);  rs = rng(opts.Seed);  c = onCleanup(@() rng(rs));
    X = table();
    for h = opts.Hemis(:)'
        Hh = B.(h);  S = Hh.S;  lbo = Hh.lbo;  lam = lbo.Lambda(:);  P = lbo.Phi;  Mm = lbo.Mass;
        nV = size(S.Vertices, 1);  nrm = double(S.VertNormals);  nrm = nrm ./ vecnorm(nrm, 2, 2);
        fg = rheome.operators.face_gradient(double(S.Vertices), double(S.Faces));
        gfb = rheome.graphfilterbank(lam, 'Wavelet', 'logitersine', 'VoicesPerOctave', 2);  nM = gfb.NumMembers;
        G = cell2mat(arrayfun(@(m) feval(gain(gfb, m), lam), 1:nM, 'uni', 0));  wl = 1e3 * wavelengths(gfb);
        for m0 = 2:nM-1
            for rep = 1:opts.Plants
                v0 = randi(nV);
                psi = P * (G(:, m0) .* (P' * (Mm * full(sparse(v0, 1, 1, nV, 1)))));
                g = [fg.W*(fg.Gx*psi) fg.W*(fg.Gy*psi) fg.W*(fg.Gz*psi)];  g = g - sum(g.*nrm, 2) .* nrm;
                Jf = [reshape(g', [], 1) reshape(cross(nrm, g, 2)', [], 1)];
                Bd = rheome.differential.helmholtzbands(Jf, S, lbo);
                [~, mi] = max(Bd.Eirr(:, 1));  [~, ms] = max(Bd.Esol(:, 2));
                X = [X; {h, m0, wl(m0), rep, v0, mi, ms, 100*sum(Bd.Esol(:,1))/sum(Bd.Eirr(:,1)), ...
                         100*sum(Bd.Eirr(:,2))/sum(Bd.Esol(:,2)), 100*Bd.harmFrac(1), 100*Bd.harmFrac(2), ...
                         Bd.capture(1,1), Bd.capture(2,2)}]; %#ok<AGROW>
            end
        end
    end
    X.Properties.VariableNames = {'hemi','m0','wavelengthMM','plant','vertex','bandSource','bandVortex', ...
        'leakSource_pct','leakVortex_pct','harmSource_pct','harmVortex_pct','captureSource','captureVortex'};
    T = rheome.scale.rows("helmholtzbands", strings(0,1), [], "");
    okS = X.bandSource == X.m0;  okV = X.bandVortex == X.m0;
    for m0 = unique(X.m0)'
        k = X.m0 == m0;  lab = string(sprintf('%.0f mm', median(X.wavelengthMM(k))));
        T = [T; rheome.scale.rows("helmholtzbands", ["recovered_source" "recovered_vortex" "leak_source_pct" ...
              "leak_vortex_pct" "harm_source_pct" "harm_vortex_pct" "n_plants"], ...
              [mean(okS(k)) mean(okV(k)) median(X.leakSource_pct(k)) median(X.leakVortex_pct(k)) ...
               median(X.harmSource_pct(k)) median(X.harmVortex_pct(k)) nnz(k)], ...
              ["fraction" "fraction" "percent" "percent" "percent" "percent" "plants"], lab)]; %#ok<AGROW>
    end
    for h = opts.Hemis(:)'
        k = X.hemi == h;
        T = [T; rheome.scale.rows("helmholtzbands", "recovered_all_" + h, mean([okS(k); okV(k)]), "fraction")]; %#ok<AGROW>
    end
end

% Author: Diellor Basha, 2026
