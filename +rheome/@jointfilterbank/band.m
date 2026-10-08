function h = band(f1, f2, width)
% JOINTFILTERBANK.BAND  A raised-cosine temporal band member, as an @(omega) handle.
%
%   h = rheome.jointfilterbank.band(f1, f2)          % width = min(0.1*(f2-f1), 1.0) Hz
%   h = rheome.jointfilterbank.band(f1, f2, width)   % width = 0 -> brick wall
%
% ⭐ THE BAND EDGE IS A FILTER MEMBER, NOT A PROPERTY OF THE TRANSFORM. Keeping a set of
% frequency bins and discarding the rest IS a band-pass, and if nothing shapes its edges
% it is a RECTANGULAR one -- a sinc in time decaying as 1/t. Since the FFT treats the
% record as periodic, that tail wraps back into the retained band. A raised cosine at
% each edge removes it.
%
% h reaches EXACTLY ZERO at both band edges, which is why a subsequent truncation to the
% retained bins is LOSSLESS rather than approximately so: the discarded bins sit beyond a
% gain that is already zero.
%
% ⚠ NEVER APPLY THIS TO A PASS WHOSE SPECTRUM WILL BE FITTED. It removes the lowest bins,
% which carry most of the leverage on a 1/f slope. Pass width = 0 for such a pass.
%
% INPUTS:  f1, f2  pass band in Hz;  width  Hz of raised cosine at EACH edge. The default
%          caps a fraction-of-bandwidth rule at 1 Hz, because 10% of a [1 45] Hz band is
%          4.4 Hz of taper eating the low-frequency end. Clamped to at most half the
%          bandwidth so the two edges cannot overlap.
% OUTPUT:  h  @(omega) -> gains in [0,1], same shape as omega (rad/s in, so it composes
%          directly as a psi_t factor).
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(width), width = min(0.1*(f2-f1), 1.0); end
    bw = f2 - f1;
    if ~(bw > 0)
        error('jointfilterbank:band', 'band must be increasing, got [%g %g].', f1, f2);
    end
    width = max(0, min(width, bw/2));
    h = @(w) i_band(w, f1, f2, width);
end

function g = i_band(omega, f1, f2, width)
    ff = double(omega) / (2*pi);
    g  = double(ff > f1 & ff < f2);              % brick wall; shaped below
    if width > 0
        lo = ff > f1 & ff < f1 + width;
        hi = ff > f2 - width & ff < f2;
        g(lo) = 0.5 * (1 - cos(pi * (ff(lo) - f1) / width));
        g(hi) = 0.5 * (1 + cos(pi * (ff(hi) - (f2 - width)) / width));
    end
    g(ff == f1 | ff == f2) = 0;                  % exactly zero at the edges
    g = max(g, 0);
end

% Author: Diellor Basha, 2026
