function C = synthesis(Y, Gd)
% JTV.SYNTHESIS  Reconstruct a joint spectrum from filterbank coefficients.
%
%   C = rheome.jtv.synthesis(Y, Gd)
%
% With Gd the canonical dual of the analysing bank, C = sum_i gd_i .* Y_i recovers the original
% exactly wherever the bank covers.
%
% INPUTS:  Y [K x nOmega x Nf] coefficients;  Gd [K x nOmega x Nf] dual gains
% OUTPUT:  C [K x nOmega]
%
% See also: rheome.jtv.analysis, rheome.jtv.dual
%
% Author: Diellor Basha, 2026 (after F. Grassi, GSPBox)

    if ~isequal(size(Y), size(Gd))
        error('jtv:synthesis:size', 'Y is %s but Gd is %s.', mat2str(size(Y)), mat2str(size(Gd)));
    end
    C = sum(conj(Gd) .* Y, 3);
end

% Author: Diellor Basha, 2026
