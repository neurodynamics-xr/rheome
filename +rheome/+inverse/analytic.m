function R = analytic(Gain, NoiseCovMat, basis, opts)
% INVERSE.ANALYTIC  Inverse in the leadfield's OWN measurable basis, truncated at the noise floor.
%
%   R = rheome.inverse.analytic(Gain, NoiseCovMat, basis)
%   R = rheome.inverse.analytic(Gain, NoiseCovMat, basis, struct('K',44))
%
% ⭐⭐ WHY THIS IS NOT ANOTHER SMOOTHNESS PRIOR. rheome.inverse.dirac reduces 3nV unknowns to K Dirac modes,
% which is a real reduction but the wrong subspace: measured, leadfield rows lie only 0.187 inside the
% 400-mode span, just 10 of 270 measurable directions have a principal-angle cosine above 0.5 to it,
% and the mean cos^2 is 0.041 (§43). Because rank(G*Psi) is still 270 the truncated model fits any
% sensor pattern exactly, so the misalignment never shows as a bad fit -- it shows as a source
% estimate dominated by near-null directions. ⭐ This truncates in the span of the leadfield ITSELF,
% where every retained direction is one the array actually measures.
%
% ⭐ AND THE BASIS IS ANALYTIC. `basis` is the analytic reconstruction of the leadfield
% (rheome.forward.analyticbank), which spans the true row space to 0.9988 and is built from a closed-form
% law -- a sphere centre, the sensor position and orientation, and a few envelope coefficients. So
% the source basis is not extracted numerically from a sampled matrix; it is written down.
%
% ⚠⚠ IT DOES NOT RECOVER THE INVISIBLE PART, AND NOTHING DOES. Only 0.276 of a planted vortex and
% 0.066 of a focal dipole patch lies in the measurable subspace at all. This makes the estimate of
% the visible part over-determined -- about 44 coordinates rise above the noise at SNR 3, against
% 270 measurements -- rather than regularised into existence. It says nothing about the rest.
%
% INPUTS
%   Gain         [nCh x 3nV]         NoiseCovMat  struct with .NoiseCov, as rheome.inverse.mne
%   basis        [nCh x 3nV] the analytic reconstruction; its rows span the source basis
%   opts .K      directions to keep ([] = the discrepancy choice from the noise floor)
%        .Snr    assumed SNR for that choice (3)
%
% OUTPUT (struct R)
%   .ImagingKernel [3nV x nCh]   .K  directions kept   .spectrum  the truncated singular values
%
% See also: rheome.inverse.mne, rheome.inverse.dirac, rheome.forward.analyticbank, rheome.forward.leadfield
%
% Author: Diellor Basha, 2026

    if nargin < 4, opts = struct(); end
    if ~isfield(opts,'Snr') || isempty(opts.Snr), opts.Snr = 3; end
    [nCh, n3] = size(Gain);
    assert(isequal(size(basis), [nCh n3]), ...
        'rheome.inverse.analytic: basis must be [%d x %d], got %s.', nCh, n3, mat2str(size(basis)));

    C = NoiseCovMat.NoiseCov;
    [Uc, Sc] = svd((C+C')/2);
    sc = diag(Sc);  sc = max(sc, max(sc)*1e-8);
    iW = Uc*diag(1./sqrt(sc))*Uc';                 % the whitener, as rheome.inverse.mne

    B  = basis';                                    % [3nV x nCh] the source basis
    GB = iW*Gain*B;                                 % [nCh x nCh] whitened forward in that basis
    [U, Sg, Vv] = svd(GB, 'econ');
    s = diag(Sg);

    K = [];
    if isfield(opts,'K') && ~isempty(opts.K), K = round(opts.K); end
    if isempty(K)
        % ⭐ keep the directions the data determines: a singular value counts when it lifts a
        %   unit-amplitude coefficient above the whitened noise, s > s(1)/(Snr*sqrt(nCh)).
        K = max(1, sum(s > s(1)/(opts.Snr*sqrt(nCh))));
    end
    K = min(K, numel(s));

    R.ImagingKernel = B*(Vv(:,1:K)*diag(1./s(1:K))*U(:,1:K)')*iW;
    R.K = K;
    R.spectrum = s;
    R.Whitener = iW;
end

% Author: Diellor Basha, 2026
