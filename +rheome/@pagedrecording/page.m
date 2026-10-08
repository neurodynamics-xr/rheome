function p = page(obj, k)
% PAGE  Read page k: its data, its time axis, and which samples it owns.
%
%   p = page(pr, k)
%
% OUTPUT (struct p):
%   .F        [NumChannels x nSpan] the data actually read, in pr.Precision
%   .t        [1 x nSpan]  time of each column (seconds, on the full axis)
%   .samples  [1 x nSpan]  sample index of each column, on the full axis
%   .core     [1 x nSpan] logical -- the columns this page OWNS
%   .clipped  [pre post]   margin the record could not supply (0 0 for interior pages)
%   .index    k            .numPages  total
%
% ⚠ REDUCE OVER .core, NOT over the whole span. The margin belongs to the neighbouring
% pages too; summing spans counts it twice.
%
% See also: pagerange, read, rheome.pagedrecording
%
% Author: Diellor Basha, 2026

    rng = pagerange(obj, k);
    s1  = rng.span(1);  s2 = rng.span(2);

    p          = struct();
    p.index    = k;
    p.numPages = obj.NumPages;
    p.samples  = s1:s2;
    p.core     = p.samples >= rng.core(1) & p.samples <= rng.core(2);
    p.clipped  = rng.clipped;
    p.t        = obj.Time(p.samples);
    p.F        = read(obj, s1, s2);
end

% Author: Diellor Basha, 2026
