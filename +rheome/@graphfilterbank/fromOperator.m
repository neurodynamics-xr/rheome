function gfb = fromOperator(L, varargin)
% GRAPHFILTERBANK.FROMOPERATOR  Build a bank straight from an operator.
%
%   gfb = rheome.graphfilterbank.fromOperator(L)
%   gfb = rheome.graphfilterbank.fromOperator(L, 'Mass', M, 'Order', 40, ...)
%
% Estimates lambda_max by Lanczos (one eigenvalue, no decomposition), builds the bank
% on that range and attaches a Chebyshev transform. GSPBox's ergonomics without its
% coupling -- and the safe path, since it guarantees the polynomial's lmax bounds the
% operator's real spectrum.
%
% Any other name-value pair is forwarded to the graphfilterbank constructor.
%
% See also: rheome.graphtransform.chebyshev, rheome.graphtransform.eigen
%
% Author: Diellor Basha, 2026

    p = inputParser;  p.KeepUnmatched = true;
    p.addParameter('Mass',  []);
    p.addParameter('Order', 40);
    p.parse(varargin{:});
    o    = p.Results;
    rest = namedargs2cell(p.Unmatched);

    opts = struct('tol', 1e-4, 'maxit', 500, 'disp', 0);
    if isempty(o.Mass), lmax = eigs(L, 1, 'largestabs', opts);
    else,               lmax = eigs(L, o.Mass, 1, 'largestabs', opts);
    end
    lmax = abs(double(lmax));

    T   = rheome.graphtransform.chebyshev(L, lmax, 'Order', o.Order, 'Mass', o.Mass);
    gfb = rheome.graphfilterbank(lmax, rest{:}, 'Transform', T);
end

% Author: Diellor Basha, 2026
