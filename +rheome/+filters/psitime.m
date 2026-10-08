function h = psitime(f, band, width)
% FILTERS.PSITIME  The temporal band member psi_time(omega) -- a raised-cosine band-pass.
%
%   h = rheome.filters.psitime(f, band)          % width = min(0.1*bandwidth, 1.0) Hz
%   h = rheome.filters.psitime(f, band, width)   % explicit edge width, Hz;  width = 0 -> brick wall
%
% THE BAND EDGE IS A FILTER MEMBER, NOT A PROPERTY OF THE TRANSFORM. Keeping a set of frequency
% bins and discarding the rest IS a band-pass, and if nothing shapes its edges it is a RECTANGULAR
% one -- a sinc in time decaying as 1/t. Since the FFT treats the record as periodic, that tail
% wraps back into the retained band. A raised cosine at each edge removes it.
%
%   h(f) = 0                                                    f <= band(1)  or  f >= band(2)
%          (1 - cos(pi*(f-band(1))/width))/2                    rising edge
%          1                                                    interior
%          (1 + cos(pi*(f-(band(2)-width))/width))/2            falling edge
%
% h reaches EXACTLY ZERO at both band edges, which is why the subsequent truncation to the
% retained bins is LOSSLESS rather than approximately so: the discarded bins sit beyond a gain
% that is already zero.
%
% MEASURED (subject01, 8-13 Hz, each route against ITSELF on a longer record so that only edge
% effects differ; relative interior error beyond a 2.5 s guard):
%
%   truncate, no shaping ................ 0.0194
%   RAISED-COSINE EDGE (this) ........... 0.0004    <- 50x better than a hard cut
%   mirror-pad 3 s only ................. 0.0200    <- padding alone does nothing
%   pad + shaped edge ................... 0.0003    <- padding adds nothing on top
%   FIR prefilter (Brainstorm-style) .... 0.0166    <- 40x WORSE than the shaped edge
%
% None of these fixes the RECORD edges themselves (~0.25-0.36 within 1 s for every route): the
% first and last samples have less context, so a guard is still required.
%
% ⚠ NEVER APPLY THIS TO A PASS WHOSE SPECTRUM WILL BE FITTED. It removes the lowest bins, which
% carry most of the leverage on a 1/f slope. Measured on subject01 over [1 45] Hz: width 0 Hz
% gives chi in [0.57 1.44]; width 0.5 Hz gives chi in [-5.08 -3.76] -- negative exponents, i.e.
% power RISING with frequency. Estimation passes use rheome.flow.psd (Welch), which tapers in TIME.
%
% INPUTS:
%   f      [1 x nOmega] or [nOmega x 1] frequencies (Hz)
%   band   [f1 f2] pass band (Hz)
%   width  Hz of raised cosine at each edge. Default min(0.1*(f2-f1), 1.0) -- a
%          fraction-of-bandwidth rule alone is catastrophic on a WIDE band (0.2*bandwidth over
%          [1 45] Hz is 8.8 Hz of taper eating the low-frequency end), hence the cap. Clamped to
%          at most half the bandwidth so the two edges cannot overlap.
%
% OUTPUT:
%   h      same shape as f, gains in [0,1]
%
% Use it directly as a psiTime kernel of a joint bank:
%   rheome.jtv.bank(Lambda, f, frame.g, {@(w) rheome.filters.psitime(w/(2*pi), band)}, [], 'type','js')
%
% See also: rheome.flow.joint, rheome.jtv.bank, rheome.filters.firbandpass
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(width), width = min(0.1*(band(2)-band(1)), 1.0); end
    bw = band(2) - band(1);
    if ~(bw > 0), error('filters:psitime:band', 'band must be increasing, got [%g %g].', band(1), band(2)); end
    width = max(0, min(width, bw/2));            % the two edges must not overlap

    ff = double(f);
    h  = double(ff >= band(1) & ff <= band(2));  % brick wall; shaped below

    if width > 0
        lo = ff > band(1) & ff < band(1) + width;
        hi = ff > band(2) - width & ff < band(2);
        h(lo) = 0.5 * (1 - cos(pi * (ff(lo) - band(1)) / width));
        h(hi) = 0.5 * (1 + cos(pi * (ff(hi) - (band(2) - width)) / width));
        h(ff == band(1) | ff == band(2)) = 0;    % exactly zero at the edges
    end
    h = max(h, 0);
end

% Author: Diellor Basha, 2026
