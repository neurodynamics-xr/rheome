function out = joint(F, sfreq, coeffOperator, band, opts)
% FLOW.JOINT  Apply a flow kernel in the JOINT (spatial-mode x temporal-frequency) domain.
%
%   out = rheome.flow.joint(F, sfreq, coeffOperator, band)
%   out = rheome.flow.joint(F, sfreq, coeffOperator, band, opts)
%
% Transforms the sensor series once, keeps only the POSITIVE in-band bins, and applies the flow
% kernel there:
%
%   C(lambda,omega) = coeffOperator * Btilde(:,omega)          one GEMM, [Ks x nOmega]
%
% WHY IN FREQUENCY. The MEG forward model is quasi-static, so the kernel carries no frequency
% dependence and acts on the SENSOR index alone -- every bin is mapped identically. Discarding
% out-of-band bins BEFORE the multiply is therefore free, and a narrow band is sparse in
% frequency: at 600 s the resolution is ~1/600 Hz, so an 8-13 Hz band is ~3e3 of ~1.8e5 bins.
% The same GEMM in the decimated time domain would run over ~7e4 columns.
%
% Retaining positive frequencies only makes the result the ANALYTIC signal by construction (the
% sensor series is real, so its spectrum is conjugate-symmetric); no Hilbert transform is needed
% anywhere downstream. Use rheome.flow.synth to return to time at any output rate.
%
% ⚠ TAPER THE BAND EDGES -- opts.edgeTaper, ON BY DEFAULT. Keeping a set of bins and discarding the
% rest IS a bandpass, with a RECTANGULAR transfer function: a sinc in time decaying as 1/t. Because
% the FFT treats the record as periodic, that long tail wraps around and contaminates the retained
% band. A raised cosine at each band edge shortens the tail and the problem goes with it.
%
% MEASURED (a reference subject, 8-13 Hz, each route against ITSELF on a longer record, so only edge effects
% differ). Relative error in the interior, beyond a 2.5 s guard:
%
%   truncate, no taper .................. 0.0194
%   edge taper 1 Hz ..................... 0.0004     <- 50x better, and this is the whole fix
%   mirror-pad 3 s only ................. 0.0200     <- padding alone does nothing
%   pad + taper ......................... 0.0003     <- padding adds nothing on top of the taper
%   FIR prefilter (bst-style) + truncate  0.0166     <- 40x WORSE than the taper
%
% So the fix belongs in the joint domain, not in a separate time-domain stage: this is a weighting
% of coefficients like every other filter here. opts.prefilter and opts.pad remain available and
% are documented, but neither is needed.
%
% ⚠ NONE of these fixes the EDGES themselves (~0.25-0.36 within 1 s for every route). The first and
% last samples simply have less context; guard them.
%
% INPUTS:
%   F              [C x nT] sensor time series (channels already selected)
%   sfreq          sampling rate (Hz)
%   coeffOperator  [Ks x C] flow kernel, e.g. kc.coeffOperator from rheome.flow.curl
%   band           [f1 f2] retained band (Hz)
%   opts  .prefilter  false (default) | true | a spec from rheome.filters.firbandpass
%                     apply a declared FIR bandpass before the transform (RECOMMENDED)
%         .fir        options forwarded to rheome.filters.firbandpass (.tranBand .ripple .atten .mode)
%         .pad        SECONDS of mirror padding before the transform (recommended: >= 2/bandwidth).
%                     Makes the circular convolution linear over the original support, which is
%                     what the periodicity artefact actually needs -- no filter shape fixes it.
%         .edgeTaper  Hz of raised cosine at each band edge (the shape fix, in the joint domain)
%         .lambda     eigenvalues of the scalar basis, so the returned .axes is COMPLETE.
%                     Without it the axes carry the frequency side only and rheome.jtv.compatible cannot
%                     check that two representations share a basis.
%         .taper      'none' (default) | 'hann'   taper before the transform
%         .nfft       transform length (default: nT)
%
% OUTPUT (struct out):
%   .C     [Ks x nOmega] complex joint spectrum
%   .f     [1 x nOmega]  bin frequencies (Hz), strictly positive
%   .df    bin spacing (Hz)         .nfft   transform length used
%   .nT    original sample count    .sfreq  original rate
%   .kept  nOmega / (nfft/2) -- the fraction of the positive axis retained
%   .psiTime  [1 x nOmega] the band member as applied (rheome.filters.psitime). It reaches EXACTLY zero
%          at both band edges, which is what makes the truncation to .f lossless rather than
%          approximately so -- the discarded bins lie beyond a gain that is already zero.
%   .edgeTaper  the edge width actually used (Hz)
%   .axes  rheome.jtv.axes object: the frequency axis, sampling, retained bands, half convention and the
%          implied time boundary ('ring' for a plain FFT, 'path' if mirror-padded). Attach the
%          eigenvalue axis with rheome.jtv.axes(...,'lambda',Lambda) and check pairs with rheome.jtv.compatible.
%   .t0    time origin offset: pass as opts.t0 to rheome.flow.synth so the output lands on the original
%          samples. Nonzero whenever .pad was used; ignoring it misaligns the whole record.
%   .fir   the FIR spec if prefiltered, else []
%   .edgeSamples  samples at each end the filter touches (0 if not prefiltered)
%
% See also: rheome.flow.synth, rheome.flow.curl, rheome.spectral.decompose
%
% Author: Diellor Basha, 2026

    if nargin < 5, opts = struct(); end
    if ~isfield(opts,'taper')     || isempty(opts.taper),     opts.taper = 'none'; end
    if ~isfield(opts,'prefilter') || isempty(opts.prefilter), opts.prefilter = false; end
    if ~isfield(opts,'fir')       || isempty(opts.fir),       opts.fir = struct(); end
    if ~isfield(opts,'lambda')    || isempty(opts.lambda),    opts.lambda = []; end
    if ~isfield(opts,'pad')       || isempty(opts.pad),       opts.pad = 0;   end   % seconds of mirror pad
    % DEFAULT ON, BUT CAPPED. A raised-cosine edge is the whole fix for a narrow analysis band, but
    % a fraction-of-bandwidth rule is catastrophic on a WIDE one: 0.2*bandwidth over [1 45] Hz is
    % 8.8 Hz of taper eating the low-frequency end. Cap it.
    % ⚠ USE edgeTaper = 0 FOR ANY PASS WHOSE SPECTRUM WILL BE FITTED. Measured on a reference subject, fitting
    % the aperiodic background over 1-45 Hz: taper 0 Hz -> chi in [0.57 1.44]; taper 0.5 Hz -> chi in
    % [-5.08 -3.76]; taper 8.8 Hz -> [-5.92 -4.59]. NEGATIVE exponents, i.e. power rising with
    % frequency. The lowest bins carry most of the leverage on a 1/f slope and tapering removes them.
    % Better still: use rheome.flow.psd (Welch) for fitting -- it tapers in TIME, which a PSD wants, and has
    % no frequency-domain taper to corrupt.
    if ~isfield(opts,'edgeTaper') || isempty(opts.edgeTaper)
        opts.edgeTaper = min(0.1 * (band(2) - band(1)), 1.0);
    end
    [C, nT] = size(F);
    if size(coeffOperator,2) ~= C
        error('flow:joint:size', 'kernel expects %d channels but F has %d.', size(coeffOperator,2), C);
    end
    if ~isfield(opts,'nfft') || isempty(opts.nfft), opts.nfft = nT; end
    nfft = opts.nfft;

    firSpec = [];  edgeSamples = 0;
    if ~isequal(opts.prefilter, false)
        [F, firSpec] = rheome.filters.firbandpass(F, sfreq, band, opts.fir);
        edgeSamples  = firSpec.edgeSamples;
    end

    switch lower(opts.taper)
        case 'none', w = ones(1, nT);
        case 'hann', w = 0.5 - 0.5*cos(2*pi*(0:nT-1)/(nT-1));
        otherwise,   error('flow:joint:taper', 'unknown taper ''%s''.', opts.taper);
    end

    % ---- mirror-pad so the convolution is LINEAR, not circular ----
    % A whole-record FFT multiply is CIRCULAR convolution: the transform treats the record as
    % periodic, so x(T) ~= x(1) is a discontinuity whose energy spreads across every bin. No
    % filter shape removes it, because it is already mixed into the coefficients. Padding makes
    % circular convolution equal linear convolution over the original support, after which a
    % smooth joint filter is entirely legitimate -- and stays in the joint framework rather than
    % importing a separate time-domain stage.
    nPad = round(opts.pad * sfreq);
    if nPad > 0
        nPad = min(nPad, nT-1);
        F    = [F(:, nPad+1:-1:2), F, F(:, end-1:-1:end-nPad)];
        nfft = size(F,2);
        w    = ones(1, nfft);
        if strcmpi(opts.taper,'hann'), w = 0.5 - 0.5*cos(2*pi*(0:nfft-1)/(nfft-1)); end
    end

    Btil = fft(F .* w, nfft, 2);                       % [C x nfft]
    df   = sfreq / nfft;
    fAll = (0:nfft-1) * df;
    keep = find(fAll > 0 & fAll <= sfreq/2 & fAll >= band(1) & fAll <= band(2));
    if isempty(keep)
        error('flow:joint:band', 'no bins in [%g %g] Hz at df = %g Hz.', band(1), band(2), df);
    end

    % SCALING. rheome.flow.synth evaluates a direct sum with no 1/N, so the 1/nfft of the inverse
    % transform is carried here; the factor 2 completes the ANALYTIC signal, since discarding the
    % negative frequencies of a real signal halves the amplitude. abs(c) is then the envelope in
    % the units of F.
    Bk = Btil(:, keep);
    % ---- the band member psi_time(omega) ----
    % The band edge is a FILTER MEMBER, not a property of the transform, so it lives in
    % rheome.filters.psitime and is applied here like any other gain. Keeping it there rather than inline
    % means the same member can be handed to rheome.jtv.bank, and that a test comparing this route against
    % a time-domain one compares the same FILTER rather than two different ones (which is exactly
    % how flow_joint_validate came to report a 3.4e-2 "equivalence failure" that was a taper the
    % reference did not have).
    psiTime = rheome.filters.psitime(fAll(keep), band, opts.edgeTaper);
    Bk      = Bk .* psiTime;
    out.C     = (2/nfft) * (coeffOperator * Bk);              % [Ks x nOmega]  <-- whole recording
    out.f     = fAll(keep);
    out.psiTime   = psiTime;                                  % the member as applied, [1 x nOmega]
    out.edgeTaper = opts.edgeTaper;
    % ⚠ UNITS. |C|^2 is amplitude-squared PER BIN, not a power spectral DENSITY: the (2/nfft) scaling
    % above makes |C| an envelope in the units of F, which is what synthesis and tracking want. A PSD
    % (rheome.flow.psd, units^2 per Hz) is a DIFFERENT quantity, and the two differ by nfft/(2*fs) -- a factor
    % that GROWS WITH RECORD LENGTH: 4x at 8 s, 300x at 600 s. Comparing them absolutely is a unit
    % error that a short test cannot expose and a long recording turns catastrophic (measured: the
    % periodic component came back identically zero on 600 s). To fit a background on a PSD and apply
    % it here, divide:  fitPower = psd.P / out.psdScale.
    out.psdScale  = nfft / (2*sfreq);                         % |C|^2 * psdScale = PSD (to ~5%)
    out.df    = df;
    out.nfft  = nfft;
    out.nT    = nT;
    out.nPad  = nPad;
    % ⚠ TIME ORIGIN. The padded record begins nPad samples BEFORE the original, so the synthesis
    % must be evaluated from t0 = nPad/sfreq to land on the original samples. Ignoring this
    % silently misaligns everything by the pad length -- measured as a 26%% "interior error" that
    % was pure time shift, not filtering.
    out.t0    = nPad / sfreq;

    % ---- the axes, as a first-class object ----
    % Carry them rather than passing Lambda and f loose: rheome.jtv.compatible can then check that two
    % representations belong together before they are combined. Two per-band runs truncate to
    % DIFFERENT bins and do not share a frequency axis, even though both are "the joint spectrum".
    out.axes = rheome.jtv.axes('f', out.f, 'fs', sfreq, 'nT', nT, 'nfft', nfft, ...
        'lambda', opts.lambda, 'bands', band(:).', 'half', 'positive', 't0', out.t0, ...
        'boundary', ternary(nPad > 0, 'path', 'ring'));
    out.sfreq = sfreq;
    out.kept  = numel(keep) / floor(nfft/2);
    out.fir   = firSpec;
    out.edgeSamples = edgeSamples;
end

function y = ternary(c,a,b), if c, y=a; else, y=b; end, end

% Author: Diellor Basha, 2026
