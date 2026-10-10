function frame = frame(family, Nf, lrange, varargin)
% FILTERS.FRAME  Design a spectral-graph wavelet FRAME (filterbank) that tiles the spectrum.
%
%   frame = rheome.filters.frame(family, Nf, lrange)   % lrange = [lmin lmax] over the eigenvalue axis
%
% Builds a bank of M spectral filters g_1(lambda)..g_M(lambda) covering [lmin lmax], the
% multi-member generalization of a single filter. Applied through an eigenbasis
% (rheome.filters.frame_analysis / rheome.filters.frame_synthesis) it is a spectral-graph WAVELET TRANSFORM
% (Hammond/Vandergheynst/Gribonval 2011; GSPBox, Perraudin et al. 2014), using the exact
% eigenbasis Phi*diag(g)*Phi' (not the Chebyshev approximation). Port of bst_eigenwavelet('Design').
%
% FAMILIES:
%   'itersine' : half-cosine TIGHT frame (sum_m g_m^2 = const over the interior -> EXACT
%                reconstruction), M = Nf members.
%   'mexhat'   : Nf log-scaled band-pass mexhat wavelets + 1 low-pass scaling function, M = Nf+1.
%   'heat'     : Nf low-pass heat scales (a diffusion scale-space, coarse -> fine), M = Nf.
%
% INPUTS:
%   family  'itersine' | 'mexhat' | 'heat'
%   Nf      number of wavelet scales. [] (recommended) DERIVES it from the scale range at
%           'perOctave' members per octave -- see below, because a fixed Nf is wrong for the same
%           reason a fixed lambda_max is: the right count depends on how wide the range turned out.
%   lrange  [lmin lmax], OR the full eigenvalue vector (then lmin = the smallest NONZERO eigenvalue)
%
% NAME-VALUE:
%   'lpfactor'  GSPBox's constant, lmin = lmax/lpfactor (default 20). Sets the COARSE end.
%   'sigmaMax'  pin the coarse end in METRES instead. ⭐ Use this: it makes the bank INVARIANT to K,
%               which a physical motif size requires. lpfactor still leaves sigma_max proportional
%               to 1/sqrt(lmax), so a richer basis silently narrows the coarse end.
%   'perOctave' members per octave when Nf is derived (default 3)
%   'tmin'      'basis' (default, t_min = 1/lmax -- tiles the WHOLE spectrum, GSP's choice) |
%               'usable' (t_min = 4.744/lmax -- every member responds properly, but the top of the
%               spectrum is then UNCOVERED, so A collapses and reconstruction loses it. Prefer
%               choosing K so lmax lands where you trust it, over refusing members after the fact.)
%
% ⭐ THE SCALE RANGE MUST NOT DEPEND ON Nf -- this follows GSPBox exactly. gsp_design_mexican_hat
% sets lmin = lmax/lpfactor with lpfactor a CONSTANT, then gsp_wlog_scales spaces Nf-1 scales over
% t in [1/lmax, 2/lmin]. Nf therefore controls only the DENSITY of the sampling, never the reach.
% Ours had t_max = 2/max(lmin, lmax/(4*Nf)) = 8*Nf/lmax, which tied the coarse end to the member
% count: raising K from 400 to 1000 at fixed Nf = 6 shrank the coarsest member 46.6 -> 27.9 mm and
% sent 50% of detections onto the low-pass, where size is Inf.
%
% OUTPUT (struct frame):
%   .Family .Nf .g {1 x M} member gain handles g_m(lambda)  .Lrange [lmin lmax]
%   .t       [1 x M]  the scale parameter of each member (NaN where the family has none, and for
%                     the mexhat scaling function, which is a low-pass rather than a scale)
%   .Sigma   [1 x M]  the member's spatial scale in METRES, sqrt(2*t) -- EXACT. ⚠ REPORT THIS,
%                     not .Centers: a mexhat member peaks at lambda = 1/t and a Gaussian of width
%                     sigma is exp(-lambda*sigma^2/2), so t = sigma^2/2 identically.
%   .Centers [1 x M]  per-member characteristic wavenumber (gain-weighted centroid of sqrt(lambda)).
%                     A DISPLAY summary only -- it is contaminated by where the eigenvalue axis was
%                     truncated, so it is not the scale the member is matched to.
%   .Gamma            SIZE CALIBRATION: sigma_m / sigma_true. 1/sqrt(2) for mexhat (derived below),
%                     NaN for families with no matched-filter interpretation.
%   .SigmaTrue [1xM]  .Sigma / .Gamma -- the vortex size a detection on member m corresponds to.
%                     ⭐ REPORT THIS for motif sizes; .Sigma is the member's own scale.
%   .MassLost [1 x M] fraction of each member's spectral mass falling beyond lmax
%   .Usable  [1 x M]  logical, .MassLost <= 5%. Members below this under-respond.
%   .SigmaFloor       the smallest usable sigma, 3.08/sqrt(lmax)
%
% ⭐ WHY Gamma = 1/sqrt(2) AND NOT 1. Matching the PEAK of the member gain to the peak of a
% sigma-vortex's vorticity spectrum lambda*exp(-lambda*sigma^2/2) gives sigma_m = sigma. But a
% detector maximises the RESPONSE, which integrates the overlap against the mode density -- and by
% Weyl's law that density is uniform in lambda on a surface. So
%
%     Resp(t) ~ int_0^inf (t*lam*e^{-t*lam})(lam*e^{-lam*sigma^2/2}) dlam = 2t/(t + sigma^2/2)^3
%
% is maximised at t = sigma^2/4, hence sigma_m = sqrt(2t) = sigma/sqrt(2). Analytic, not fitted;
% measured 0.735 over sigma = 8-40 mm in the sphere validation.
%
% ⚠ THE FINEST MEMBERS FALL OFF THE END OF THE BASIS. A mexhat member has mass 1/t, of which
% (1+u)exp(-u) lies beyond lmax with u = t*lmax. The 'basis' policy sets t_min = 1/lmax, i.e. u = 1,
% so the finest member loses 2/e = 73.6% -- it sees only the RISING half of its own passband,
% under-responds, and pushes the scale-space maximum one member coarser. Members lose under 5% once
% u >= 4.744, i.e. sigma_m >= 3.08/sqrt(lmax). A warning is raised naming the offending members.
%
% See also: rheome.filters.frame_gains, rheome.filters.frame_bounds, rheome.filters.frame_analysis, rheome.filters.frame_synthesis
%
% Author: Diellor Basha, 2026 (port of bst_eigenwavelet('Design'))

    % ─── COMPATIBILITY SHIM ────────────────────────────────────────────────────────────
    % The numerics now live in @graphfilterbank. This function reproduces the output struct
    % FIELD FOR FIELD so every existing caller is unaffected; rheome.filters.frame_gains,
    % frame_bounds, frame_analysis, frame_synthesis and frame_scalogram consume that struct
    % and needed no change at all. New code should use graphfilterbank directly:
    %     gfb = rheome.graphfilterbank(lrange, 'Wavelet','mexhat', 'VoicesPerOctave',3);
    if nargin < 3 || isempty(lrange), error('filters:frame:lrange', 'lrange = [lmin lmax] is required.'); end
    if nargin < 2, Nf = []; end

    pp = inputParser;
    pp.addParameter('tmin', 'basis');
    pp.addParameter('lpfactor', 20);
    pp.addParameter('sigmaMax', []);
    pp.addParameter('perOctave', 3);
    pp.addParameter('warn', true);
    pp.parse(varargin{:});
    opt = pp.Results;
    U_USABLE = 4.744;                     % (1+u)exp(-u) = 0.05

    family = lower(family);
    if ~any(strcmp(family, {'mexhat','itersine','heat'}))
        error('filters:frame:family', 'unknown family ''%s'' (use mexhat|heat|itersine).', family);
    end
    if ~any(strcmpi(opt.tmin, {'basis','usable'}))
        error('filters:frame:tmin', 'tmin must be ''basis'' or ''usable''.');
    end

    % lmin/lmax exactly as this function has always parsed them, so .Lrange is unchanged.
    if numel(lrange) > 2
        lv   = sort(double(lrange(:)));
        lmax = lv(end);
        nz   = lv(lv > lmax*1e-12);
        lmin = nz(1);
    else
        lmin = lrange(1);  lmax = lrange(2);
    end
    if ~(lmax > 0) || (lmax <= lmin)
        error('filters:frame:lrange', 'invalid lrange = [%g %g].', lmin, lmax);
    end

    args = {'Wavelet', family, 'NumFilters', Nf, 'VoicesPerOctave', opt.perOctave, ...
            'LPFactor', opt.lpfactor, 'FineLimit', lower(opt.tmin)};
    if ~isempty(opt.sigmaMax)
        % sigmaMax pins only the COARSE end; the fine end stays on the tmin policy.
        switch lower(opt.tmin)
            case 'basis',  tminPolicy = 1/lmax;
            case 'usable', tminPolicy = U_USABLE/lmax;
        end
        args = [args, {'ScaleLimits', [tminPolicy, opt.sigmaMax^2/2]}];
    end

    % The class raises graphfilterbank:degenerateFrame on a DIFFERENT criterion (coverage,
    % A <= 1e-6*B) than the two warnings this function has always raised. Suppress it here
    % and re-issue filters:frame:truncated / :sparse below with their original identifiers,
    % which callers assert on or suppress.
    wsRestore = warning('off', 'graphfilterbank:degenerateFrame');
    restoreW  = onCleanup(@() warning(wsRestore));
    gfb = rheome.graphfilterbank(lrange, args{:});

    M   = gfb.NumMembers;
    sig = widths(gfb);
    lost = gfb.MassLost;
    if strcmp(family,'mexhat'), gamma = 1/sqrt(2); else, gamma = NaN; end

    % ⚠ SigmaMax is NOT gfb.SizeLimits(2). 'heat' runs one octave coarser internally and
    % 'itersine' has no scale parameter at all, so take it from the member scales
    % themselves: max() ignores NaN, and returns NaN when every entry is NaN.
    tmaxFam = max(scales(gfb));

    % Nf as this function has always reported it: the WAVELET count, so mexhat's M = Nf+1.
    NfOut = M - strcmp(family,'mexhat');

    frame = struct('Family', family, 'Nf', NfOut, 'g', {i_handles(gfb)}, ...
                   'Lrange', [lmin lmax], 'Centers', centerWavenumbers(gfb), ...
                   't', scales(gfb), 'Sigma', sig, ...
                   'Gamma', gamma, 'SigmaTrue', sig / gamma, ...
                   'MassLost', lost, 'Usable', gfb.Usable, ...
                   'SigmaFloor', gfb.FinestScale, ...
                   'SigmaMax', sqrt(2*tmaxFam), 'LminEff', gfb.SpectralLimits(1), ...
                   'PerOctave', gfb.AchievedVoicesPerOctave);

    % ⚠ NAME THE COMPROMISED MEMBERS. A member losing most of its mass past lmax still returns a
    % number, and that number looks like a size -- which is exactly why this has to be said out loud
    % at design time rather than discovered from a size distribution that piles up at the fine end.
    if opt.warn && any(~frame.Usable)
        bad = find(~frame.Usable);
        warning('filters:frame:truncated', ...
            ['%d of %d members lose >5%% of their mass beyond lambda_max = %.4g (worst %.0f%%, ' ...
             'sigma_m = %.1f mm). They under-respond, so sizes near the fine end read too coarse. ' ...
             'Usable floor is sigma_m >= %.1f mm; pass ''tmin'',''usable'' to start there.'], ...
            numel(bad), M, lmax, 100*max(lost(bad)), 1000*min(sig(bad)), 1000*frame.SigmaFloor);
    end

    % ⚠ IS THE SCALE AXIS SAMPLED DENSELY ENOUGH? With the range now fixed by the spectrum, Nf is
    % purely a density choice, and too few members means a motif whose size falls between two of
    % them is reported at whichever is nearer -- a quantisation error that looks like a measurement.
    % Below ~1.5 members/octave adjacent members barely overlap and the frame bound A collapses.
    if opt.warn && strcmp(family,'mexhat') && frame.PerOctave < 1.5 && M >= 3
        tmn = gfb.ScaleLimits(1);  tmx = gfb.ScaleLimits(2);
        warning('filters:frame:sparse', ...
            ['%d members over %.1f octaves is %.1f per octave -- the scale axis is under-sampled ' ...
             '(adjacent members are %.2fx apart). Pass Nf = [] to derive it (%d at %g/octave), or ' ...
             'raise Nf.'], NfOut, log2(frame.SigmaMax/sig(2)), frame.PerOctave, sig(3)/sig(2), ...
            max(2, round(opt.perOctave*log2(sqrt(tmx/tmn)))+1), opt.perOctave);
    end
end

function h = i_handles(gfb)
    h = cell(1, gfb.NumMembers);
    for m = 1:gfb.NumMembers, h{m} = gain(gfb, m); end
end

% Author: Diellor Basha, 2026
