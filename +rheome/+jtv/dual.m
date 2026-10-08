function [Gd, b] = dual(G, opts)
% JTV.DUAL  Canonical dual of a joint time-vertex filterbank.
%
%   [Gd, b] = rheome.jtv.dual(G)
%
% Returns gd_i = g_i / S with S = sum_j |g_j|^2, so that sum_i gd_i * g_i = 1 wherever S > 0 and
% synthesis(analysis(C)) = C exactly. This is the one thing that makes a filterbank a TRANSFORM
% rather than an analysis: without a dual, coefficients can be computed but nothing can be
% reconstructed from them.
%
% ⚠ THE DUAL ONLY EXISTS WHERE THE BANK COVERS. Dividing by S is exactly where a gap in coverage
% becomes fatal, so the bounds are returned alongside and a warning is issued if any of the joint
% plane is uncovered. Fix the bank (more members, wider range) rather than the division.
%
% INPUTS:
%   G     [K x nOmega x Nf] filter gains
%   opts  .tol  treat S <= tol*max(S) as uncovered and set the dual to 0 there (default 1e-12)
%
% OUTPUTS:
%   Gd    [K x nOmega x Nf] canonical dual gains
%   b     the frame bounds from rheome.jtv.bounds
%
% See also: rheome.jtv.bounds, rheome.jtv.analysis, rheome.jtv.synthesis
%
% Author: Diellor Basha, 2026 (after F. Grassi, GSPBox)

    if nargin < 2, opts = struct(); end
    if ~isfield(opts,'tol') || isempty(opts.tol), opts.tol = 1e-12; end
    if ndims(G) < 3, G = reshape(G, size(G,1), size(G,2), []); end

    b = rheome.jtv.bounds(G);
    S = b.S;
    bad = S <= opts.tol * max(b.B, realmin);
    if any(bad(:))
        warning('jtv:dual:uncovered', ...
            ['%.2f%% of the joint plane has no coverage (A = %.3g). The dual is set to zero there ' ...
             'and reconstruction will LOSE that content. Widen the bank rather than the tolerance.'], ...
            100*mean(bad(:)), b.A);
    end
    Sinv = zeros(size(S));
    Sinv(~bad) = 1 ./ S(~bad);
    Gd = G .* Sinv;
end

% Author: Diellor Basha, 2026
