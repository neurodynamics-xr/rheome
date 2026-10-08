function b = framebounds(obj)
% FRAMEBOUNDS  Frame operator, bounds and coverage of the joint bank.
%
%   b = framebounds(jfb)
%
% For S(lambda,omega) = sum_i |W_i|^2 the bank is a frame iff A = min S > 0. A = B is
% TIGHT. A -> 0 means part of the joint plane is uncovered and the canonical dual does
% not exist there, so iwt will LOSE that content.
%
% ⚠ COVERAGE IS A PRECONDITION for iwt, not a nicety. A band-limited bank is uncovered
% outside its band BY CONSTRUCTION -- that is not a bug, it just means iwt reconstructs
% only what the bank passes.
%
% OUTPUT: b.S [K x nOmega]  b.A  b.B  b.Tightness (B/A)  b.Uncovered  b.Worst [k j]
%
% See also: isframetight, iwt
%
% Author: Diellor Basha, 2026

    S = jfb_frameop(obj);
    b.S = S;
    b.A = min(S(:));
    b.B = max(S(:));
    b.Tightness = b.B / max(b.A, eps);
    b.Uncovered = mean(S(:) <= eps * max([b.B; realmin]));
    [~, i] = min(S(:));
    [k, j] = ind2sub(size(S), i);
    b.Worst = [k, j];
end

% Author: Diellor Basha, 2026
