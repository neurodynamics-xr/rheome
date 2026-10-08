function rng = pagerange(obj, k)
% PAGERANGE  The sample bookkeeping for page k, without reading any data.
%
%   rng = pagerange(pr, k)
%
% OUTPUT (struct rng):
%   .core     [c1 c2]  first/last sample this page OWNS (cores partition the record)
%   .span     [s1 s2]  first/last sample actually READ (core widened by Overlap, clipped)
%   .clipped  [pre post]  samples of margin the record could not supply at each end
%
% Separated from page() so the caller can size buffers, plan a loop, or check the
% clipping of the edge pages without paying for a read.
%
% See also: page, rheome.pagedrecording
%
% Author: Diellor Basha, 2026

    if ~isscalar(k) || ~isnumeric(k) || mod(k,1) ~= 0 || k < 1 || k > obj.NumPages
        error('pagedrecording:pageIndex', ...
            'Page index must be an integer in 1..%d, got %s.', obj.NumPages, mat2str(k));
    end

    c1 = (k-1) * obj.PageLength + 1;
    c2 = min(k * obj.PageLength, obj.NumSamples);       % the last page is short

    s1 = max(1,              c1 - obj.Overlap);
    s2 = min(obj.NumSamples, c2 + obj.Overlap);

    % What the record could NOT supply -- not the same as (requested - got) clipped to zero
    % only at the ends: an Overlap wider than the whole record clips at BOTH.
    pre  = obj.Overlap - (c1 - 1);
    post = obj.Overlap - (obj.NumSamples - c2);

    rng = struct('core',    [c1 c2], ...
                 'span',    [s1 s2], ...
                 'clipped', [max(0, pre), max(0, post)]);
end

% Author: Diellor Basha, 2026
