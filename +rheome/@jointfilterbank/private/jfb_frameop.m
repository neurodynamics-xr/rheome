function S = jfb_frameop(obj)
% JFB_FRAMEOP  S(lambda,omega) = sum_i |W_i|^2  -> [K x nOmega].
%
% Loops ONE MEMBER AT A TIME, so the [K x nOmega x Nf] product is never formed. S itself
% is the one unavoidable 2-D array; it is computed on demand and NOT cached, so a bank
% held in memory stays small.
%
% Author: Diellor Basha, 2026

    lam = obj.Lambda;  om = obj.Omega;
    S = zeros(numel(lam), numel(om));
    for m = 1:obj.NumMembers
        S = S + abs(jfb_member(obj, m, lam, om)).^2;
    end
end

% Author: Diellor Basha, 2026
