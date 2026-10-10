function out = decompose(C, f, opts)
% SPECTRAL.DECOMPOSE  Split a coefficient array into periodic and aperiodic parts.
%
%   out = rheome.spectral.decompose(C, f)
%   out = rheome.spectral.decompose(C, f, opts)
%
% The usual entry point: fits the aperiodic background of every ROW of C (one row per spatial
% mode) and applies the complementary gains, so that
%
%   C_per + ... and C_ap  satisfy   |C_per|^2 + |C_ap|^2 = |C|^2   pointwise.
%
% INPUTS:
%   C     [K x nF] complex coefficients, rows = spatial modes, columns = frequency bins
%   f     [nF x 1] frequencies (Hz), strictly positive
%   opts  .fitPower [K x nFb] BROADBAND power to fit the background on (STRONGLY RECOMMENDED)
%         .fitF     [nFb x 1]  its frequencies
%         others passed to rheome.spectral.aperiodic (.knee .nIter .thresh .range)
%
% ⚠ FIT BROADBAND, APPLY NARROW. Without .fitPower the background is fitted on the same bins the
% split is applied to, which fails whenever those bins are a narrow band -- there is not enough
% frequency range to identify a 1/f slope.
%
% ⚠⚠ FIT ON A WELCH PSD (rheome.flow.psd), NOT A SINGLE PERIODOGRAM. The log of a 2-dof periodogram is
% biased LOW by gamma (~0.25 in log10), so the fitted background sits ~1.8x below the truth. That
% bias is CONSTANT in frequency, so it leaves the exponent alone -- but the split depends on the
% ABSOLUTE LEVEL of Pap, not its slope, so it inflates the periodic share. Measured on a reference subject,
% 8-13 Hz: single periodogram -> 77.5%% periodic; Welch (2 s, 50%% overlap) -> 26.4%%. The Welch
% figure is the correct one, and the difference is not a detail: it is the difference between
% "alpha-band flow is mostly rhythm" and "alpha-band flow is mostly background".
%
% ⚠ AND NEVER TAPER A FIT PASS. rheome.flow.joint's edgeTaper removes the lowest bins, which carry most of
% the leverage on a 1/f slope: a 0.5 Hz taper over [1 45] Hz drove chi to -5.08..-3.76 (negative,
% i.e. power rising with frequency). rheome.flow.psd has no frequency taper and is the safe route.
%
% OUTPUT (struct out):
%   .Cper .Cap  [K x nF] complex, the two components
%   .hper .hap  [K x nF] real gains
%   .ap         struct from rheome.spectral.aperiodic, with .exponent [1 x K] = chi(lambda)
%   .partition  max |  |Cper|^2 + |Cap|^2 - |C|^2  | / max|C|^2   -- should be ~1e-15
%
% NOTE ON ORIENTATION. C is [modes x frequency] but rheome.spectral.aperiodic fits COLUMNS, so the power
% is transposed on the way in and the gains transposed back. Handled here so callers need not.
%
% See also: rheome.spectral.aperiodic, rheome.spectral.gains, rheome.spectral.carrier
%
% Author: Diellor Basha, 2026

    if nargin < 3, opts = struct(); end
    f = double(f(:));
    if size(C,2) ~= numel(f)
        error('spectral:decompose:size', 'C has %d columns but f has %d entries.', size(C,2), numel(f));
    end

    Pow = abs(C).^2;                       % [K x nF]

    if isfield(opts,'fitPower') && ~isempty(opts.fitPower)
        % ---- BROADBAND FIT, NARROW-BAND APPLICATION ----
        % The background must be fitted over a wide frequency range and only then evaluated on the
        % retained band. Fitting inside the band is not merely imprecise, it is meaningless: an
        % 8-13 Hz window gave exponents from -2.5 to +7.4 on real data, negative values included.
        if ~isfield(opts,'fitF') || isempty(opts.fitF)
            error('spectral:decompose:fitF', 'opts.fitPower requires opts.fitF.');
        end
        if size(opts.fitPower,1) ~= size(C,1)
            error('spectral:decompose:fitRows', ...
                'fitPower has %d rows but C has %d.', size(opts.fitPower,1), size(C,1));
        end
        ap  = rheome.spectral.aperiodic(opts.fitPower.', opts.fitF, opts);
        Pap = rheome.spectral.evaluate(ap, f);               % [nF x K] on the RETAINED band

        % ⚠ THE FIT AND THE TARGET MUST BE IN THE SAME UNITS, and nothing else checks it. The split is
        % h_ap = sqrt(Pap/P), an ABSOLUTE comparison, so a background fitted on a power spectral
        % DENSITY (rheome.flow.psd, units^2/Hz) cannot be applied to an amplitude-per-bin spectrum
        % (rheome.flow.joint, |C|^2) without conversion. They differ by nfft/(2*fs) -- 4x at 8 s but 300x at
        % 600 s -- so a short test passes while a long recording returns h_per = 0 at every bin.
        % Divide the PSD by joint.psdScale before passing it here.
        above = mean(Pap(:) > reshape(Pow.', [], 1));
        if above > 0.9
            error('spectral:decompose:scale', ...
                ['the fitted background exceeds the observed power at %.0f%% of bins, so the ' ...
                 'periodic residual is empty. This is a UNIT MISMATCH, not physiology: fitPower has ' ...
                 'median %.3e against |C|^2 median %.3e, a factor of %.0f. If fitPower came from ' ...
                 'rheome.flow.psd, pass  fitPower = psd.P / joint.psdScale.'], ...
                100*above, median(opts.fitPower(:)), median(Pow(:)), ...
                median(opts.fitPower(:))/max(median(Pow(:)), realmin));
        end
    else
        ap  = rheome.spectral.aperiodic(Pow.', f, opts);     % fits columns -> pass [nF x K]
        Pap = ap.Pap;
    end
    [hp, ha] = rheome.spectral.gains(Pow.', Pap);            % [nF x K]

    out.hper = hp.';                       % [K x nF]
    out.hap  = ha.';
    out.Cper = out.hper .* C;
    out.Cap  = out.hap  .* C;
    out.ap   = ap;

    lhs = abs(out.Cper).^2 + abs(out.Cap).^2;
    out.partition = max(abs(lhs(:) - Pow(:))) / max(max(Pow(:)), eps);
end

% Author: Diellor Basha, 2026
