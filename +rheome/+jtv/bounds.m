function b = bounds(G)
% JTV.BOUNDS  Frame operator and bounds of a joint time-vertex filterbank.
%
%   b = rheome.jtv.bounds(G)
%
% For a bank G(:,:,i) of joint filters on the (lambda,omega) grid, the frame operator is
%
%   S(lambda,omega) = sum_i |g_i(lambda,omega)|^2
%
% and the bank is a frame iff A = min S > 0. A = B is TIGHT (reconstruction needs no dual beyond a
% constant). A -> 0 means part of the joint plane is uncovered and the canonical dual does not
% exist there.
%
% INPUT:   G  [K x nOmega x Nf] filter gains
% OUTPUT (struct b):
%   .S          [K x nOmega] the frame operator
%   .A .B       min and max of S
%   .Tightness  B/A (1 = tight)
%   .uncovered  fraction of the grid with S below eps*max(S) -- where the dual fails
%   .worst      [k, j] index of the weakest point
%
% See also: rheome.jtv.dual, rheome.jtv.bank
%
% Author: Diellor Basha, 2026 (after F. Grassi, GSPBox)

    if ndims(G) < 3, G = reshape(G, size(G,1), size(G,2), []); end
    S = sum(abs(G).^2, 3);
    b.S = S;
    b.A = min(S(:));
    b.B = max(S(:));
    b.Tightness = b.B / max(b.A, eps);
    tol = eps * max(b.B, realmin);
    b.uncovered = mean(S(:) <= tol);
    [~, i] = min(S(:));
    [k, j] = ind2sub(size(S), i);
    b.worst = [k, j];
end

% Author: Diellor Basha, 2026
