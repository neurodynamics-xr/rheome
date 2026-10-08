function Y = wt(obj, C)
% WT  Apply the joint filterbank to a joint spectrum.
%
%   Y = wt(jfb, C)      C [K x nOmega] -> Y [K x nOmega x Nf]
%
% Every member is an elementwise multiply -- the bank costs Nf multiplies and NO
% transforms, because C and the members already live on the same (lambda, omega) grid.
%
% ⚠ The OUTPUT is [K x nOmega x Nf] and is therefore byte-guarded like jointfilters. If
% it refuses, you want a reduction (scalogram) rather than the coefficients.
%
% See also: iwt, scalogram, framebounds
%
% Author: Diellor Basha, 2026

    lam = obj.Lambda;  om = obj.Omega;
    K = numel(lam);  nO = numel(om);  Nf = obj.NumMembers;
    if ~isequal(size(C), [K nO])
        error('jointfilterbank:size', ...
            'C is %s but the bank grid is %s.', mat2str(size(C)), mat2str([K nO]));
    end
    jfb_guard(obj, [K nO Nf], 'the coefficients');

    Y = zeros(K, nO, Nf);
    if ~isreal(C), Y = complex(Y); end
    for m = 1:Nf
        Y(:,:,m) = C .* jfb_member(obj, m, lam, om);
    end
end

% Author: Diellor Basha, 2026
