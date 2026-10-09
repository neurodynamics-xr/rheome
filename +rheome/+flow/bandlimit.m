function B = bandlimit(J, S, lbo, opts)
% FLOW.BANDLIMIT  A current field, its divergence and its curl, band-limited by the graph-wavelet bank.
%
%   B = rheome.flow.bandlimit(J, S, lbo)                        % keep wavelengths >= 92 mm
%   B = rheome.flow.bandlimit(J, S, lbo, CutoffMM=130, Voices=2)
%       J    [3nV x nT] ambient current (rows x1,y1,z1,...), e.g. the minimum-norm estimate
%       S    surface: .Vertices (m) .Faces .VertNormals (.Hemi)
%       lbo  Laplace-Beltrami basis on the SAME surface: .Phi .Lambda .Mass .L
%
% ⭐ WHY THE PER-VERTEX FIELD LOOKS NOISY, AND WHY THIS IS SMOOTH. The minimum-norm current is
% smooth at the scale MEG resolves (point-spread r50 ~ 48 mm) but it is a free-orientation 3-vector
% per vertex, and on the folded cortex neighbouring normals differ by ~30 deg: its tangential
% reading flips across every sulcal wall, and the pointwise div/curl (a derivative per vertex) are
% dominated by the mesh, not by the source. Nothing in that chain band-limits. Here the field is
% split into its potentials (rheome.differential.helmholtzbands: J = grad Phi + n x grad Psi + ...),
% the potentials are kept only in the members of the TIGHT logitersine bank whose wavelength is
% >= CutoffMM, h(lambda) = sqrt(sum_{kept m} g_m(lambda)^2), and everything is rebuilt from them:
%       Phi_bl = Phi_LB (h .* cPhi)        div_bl  = Phi_LB (-lambda .* h .* cPhi)
%       Psi_bl = Phi_LB (h .* cPsi)        curl_bl = Phi_LB (-lambda .* h .* cPsi)
%       J_bl   = grad Phi_bl + n x grad Psi_bl          (tangential; normal part and harmonic dropped)
% Because the bank is tight, h <= 1 and h = 1 on the kept band: it removes, never reweights inside.
% div_bl and curl_bl are exact Laplacians of the band-limited potentials (weak form, no per-vertex
% difference), so their sign matches rheome.differential.divergence / curl on a smooth field.
%
% OUTPUT (struct B)
%   .J [3nV x nT] .Div .Curl .Phi .Psi [nV x nT]  band-limited
%   .h [K x 1] gain   .cutoffMM   .keptWavelengthMM [1 x m] band centres kept
%   .finestSigmaMM    spatial width sigma = sqrt(2t) of the finest kept member (rheome.graphfilterbank widths)
%   .Hb               the helmholtzbands struct (unfiltered potentials, band energies)
%
% See also: rheome.differential.helmholtzbands, rheome.flow.smoothness, rheome.graphfilterbank
%
% Author: Diellor Basha, 2026

    arguments
        J double
        S (1,1) struct
        lbo (1,1) struct
        opts.CutoffMM (1,1) double {mustBePositive} = 92
        opts.Voices (1,1) double {mustBePositive} = 2
    end
    Hb = rheome.differential.helmholtzbands(J, S, lbo, Voices=opts.Voices);
    lam = lbo.Lambda(:);  P = lbo.Phi;
    gfb = rheome.graphfilterbank(lam, 'Wavelet', 'logitersine', 'VoicesPerOctave', opts.Voices);
    keep = Hb.wavelengthMM >= opts.CutoffMM;
    if ~any(keep)
        error('flow:bandlimit:cutoff', 'no bank member has a wavelength >= %g mm (coarsest %.0f mm)', ...
              opts.CutoffMM, max(Hb.wavelengthMM));
    end
    h = sqrt(sum(Hb.G(:, keep).^2, 2));
    w = 1e3 * widths(gfb);
    B.Phi = P * (h .* Hb.cPhi);   B.Psi = P * (h .* Hb.cPsi);
    B.Div = P * (-lam .* h .* Hb.cPhi);   B.Curl = P * (-lam .* h .* Hb.cPsi);
    fg = rheome.operators.face_gradient(S.Vertices, S.Faces);  Nf = fg.FaceNormal;
    gp = {fg.Gx*B.Phi, fg.Gy*B.Phi, fg.Gz*B.Phi};  gs = {fg.Gx*B.Psi, fg.Gy*B.Psi, fg.Gz*B.Psi};
    ns = {Nf(:,2).*gs{3} - Nf(:,3).*gs{2}, Nf(:,3).*gs{1} - Nf(:,1).*gs{3}, Nf(:,1).*gs{2} - Nf(:,2).*gs{1}};
    B.J = zeros(size(J));
    for c = 1:3, B.J(c:3:end, :) = fg.W * (gp{c} + ns{c}); end
    B.h = h;  B.cutoffMM = opts.CutoffMM;  B.keptWavelengthMM = Hb.wavelengthMM(keep);
    B.finestSigmaMM = min(w(keep));  B.Hb = Hb;
end

% Author: Diellor Basha, 2026
