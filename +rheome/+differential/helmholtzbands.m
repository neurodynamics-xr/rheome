function B = helmholtzbands(J, S, lbo, opts)
% DIFFERENTIAL.HELMHOLTZBANDS  Source and vortex content of a field, band by band, via its potentials.
%
%   B = rheome.differential.helmholtzbands(J, S, lbo)                     % logitersine, 2 voices per octave
%   B = rheome.differential.helmholtzbands(J, S, lbo, Voices=1, Maps=true, H=Hprecomputed)
%       J    [3nV x nT] ambient vectors (rows x1,y1,z1,...)
%       S    surface: .Vertices .Faces .VertNormals (.Hemi)
%       lbo  Laplace-Beltrami basis on the SAME surface: .Phi [nV x K] .Lambda [K x 1] .Mass .L
%
% The route that separates rotation from divergence on a real cortex:
%   1. rheome.differential.helmholtz  J = grad(Phi) + N x grad(Psi) + (J.n) n + harmonic
%      two Poisson solves -- a scalar potential Phi (sources +, sinks -) and a stream function Psi
%      (vortices at its extrema)
%   2. both scalars filtered by a TIGHT scalar Laplace-Beltrami bank (graphfilterbank on lbo.Lambda,
%      sum_m g_m^2 = 1): Phi_m = lbo.Phi * (g_m .* c_Phi), likewise Psi_m
%
% ⭐ THE BAND ENERGIES ARE EXACT. The field energies are the potentials' Dirichlet energies --
% ||grad Phi||^2 = Phi' K Phi and ||N x grad Psi||^2 = Psi' K Psi, with K the cotan stiffness helmholtz
% solves with, which IS lbo.L (1.5e-16 on a reference subject) -- so in the eigenbasis
%       E_irr(m) = sum_k lambda_k g_m(lambda_k)^2 c_Phi,k^2,    E_sol(m) likewise from c_Psi
% and sum_m E(m) = the projected Dirichlet energy, by tightness. .capture says how much of each
% potential's energy the K modes hold (1000 modes reach 17 mm wavelengths on this cortex).
% ⭐ WHY NOT POINTWISE CURL AND DIV. On a reference subject neighbouring normals differ by ~30 deg and a pure
% source reads curl/div energy 0.13-0.30 per vertex (dirac_parts_omega.m); helmholtz separates the same
% atoms at 91% / 0.7%, because its sources are WEAK (integrated against hat gradients) and the potentials
% are smooth. The filtering is then on scalars, where the tight bank's exactness holds.
%
% ⭐ MEASURED (helmholtz_bands_omega.m): planted bands 315-23 mm recovered exactly on the cortex, 0.4-0.8%
% cross-leak, a mixed source + vortex split with both peaks within 3 mm. ⚠ Through MEG + MNE the limit
% is SCALE (helmholtz_noise_omega.m): sources and vortices at >= 130 mm are located to 7-20 mm down to
% -5 dB, vortices hold to 92 mm above 0 dB, nothing below 65 mm is located even noise-free; a source's
% max-energy band is 3 bands too fine, so choose its band in advance.
%
% OUTPUT (struct B)
%   .Eirr .Esol     [M x nT]  band energies of grad(Phi) and N x grad(Psi)
%   .irrTotal .solTotal [1 x nT]  Phi'K Phi and Psi'K Psi (unprojected)
%   .capture        [2 x nT]  projected / unprojected Dirichlet energy, Phi and Psi
%   .cPhi .cPsi     [K x nT]  mode coefficients of the potentials
%   .G [K x M] gains    .wavelengthMM [1 x M] band centres (graphfilterbank wavelengths)
%   .normalFrac .harmFrac [1 x nT]  energy shares of the normal and harmonic parts (lumped area)
%   .PhiBand .PsiBand [nV x M x nT]  band-limited potentials, when Maps=true
%   .H              the helmholtz struct
%
% See also: rheome.differential.helmholtz, rheome.graphfilterbank, rheome.differential.diracparts
%
% Author: Diellor Basha, 2026

    arguments
        J double
        S (1,1) struct
        lbo (1,1) struct
        opts.Voices (1,1) double {mustBePositive} = 2
        opts.Wavelet (1,1) string = "logitersine"
        opts.Maps (1,1) logical = false
        opts.H = []
    end
    if isempty(opts.H), H = rheome.differential.helmholtz(J, S); else, H = opts.H; end
    lam = lbo.Lambda(:);  P = lbo.Phi;  M = lbo.Mass;  K = lbo.L;
    gfb = rheome.graphfilterbank(lam, 'Wavelet', char(opts.Wavelet), 'VoicesPerOctave', opts.Voices);
    G = cell2mat(arrayfun(@(m) feval(gain(gfb, m), lam), 1:gfb.NumMembers, 'uni', 0));
    B.cPhi = P' * (M * H.Phi);  B.cPsi = P' * (M * H.Psi);
    B.Eirr = (G.^2)' * (lam .* B.cPhi.^2);
    B.Esol = (G.^2)' * (lam .* B.cPsi.^2);
    B.irrTotal = sum(H.Phi .* (K * H.Phi), 1);
    B.solTotal = sum(H.Psi .* (K * H.Psi), 1);
    B.capture = [sum(lam .* B.cPhi.^2, 1) ./ max(B.irrTotal, eps); sum(lam .* B.cPsi.^2, 1) ./ max(B.solTotal, eps)];
    B.G = G;  B.wavelengthMM = 1e3 * wavelengths(gfb);
    a = full(sum(M, 2));  Nv = S.VertNormals ./ max(vecnorm(S.VertNormals, 2, 2), eps);
    Jn = J(1:3:end,:).*Nv(:,1) + J(2:3:end,:).*Nv(:,2) + J(3:3:end,:).*Nv(:,3);
    tot = sum(a .* H.Fmag.^2, 1);
    B.normalFrac = sum(a .* Jn.^2, 1) ./ max(tot, eps);
    B.harmFrac = sum(a .* H.Hmag.^2, 1) ./ max(tot, eps);
    if opts.Maps
        nT = size(J, 2);  nM = gfb.NumMembers;
        B.PhiBand = zeros(size(P, 1), nM, nT);  B.PsiBand = B.PhiBand;
        for t = 1:nT
            B.PhiBand(:, :, t) = P * (G .* B.cPhi(:, t));
            B.PsiBand(:, :, t) = P * (G .* B.cPsi(:, t));
        end
    end
    B.H = H;
end

% Author: Diellor Basha, 2026
