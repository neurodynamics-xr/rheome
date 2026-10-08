function T = chebyshev(L, lmax, varargin)
% GRAPHTRANSFORM.CHEBYSHEV  Apply g(L) by polynomial recursion -- no eigenvectors.
%
%   T = rheome.graphtransform.chebyshev(L)                 % lmax estimated by Lanczos
%   T = rheome.graphtransform.chebyshev(L, lmax)
%   T = rheome.graphtransform.chebyshev(L, lmax, 'Order', 40, 'Mass', M)
%
% This is why graphfilterbank requires only a spectral RANGE and not an eigenvalue
% vector: on a large graph you never compute the spectrum. g(lambda) is expanded in
% Chebyshev polynomials on [0, lmax] and evaluated by Clenshaw recursion in L.
%
% For a generalized pencil L phi = lambda M phi, pass 'Mass' and the recursion runs on
% M\L, whose spectrum is lambda.
%
% ⚠ lmax MUST BOUND THE OPERATOR'S WHOLE SPECTRUM, and it is NOT the lambda_max of a
% truncated eigenbasis. The shifted operator L/a - I maps [0, lmax] onto [-1, 1]; any
% true eigenvalue beyond lmax lands outside, where T_k grows like (2x)^k and the
% recursion diverges silently to 1e31 and then to NaN. Measured on an ico4 sphere: the
% top of a 100-mode basis is 9.3e3 while the operator's true lambda_max is 4.7e5, a
% factor of 50 -- so reusing a basis lambda_max here is catastrophic, not merely
% inaccurate. Omit lmax (it is then estimated by Lanczos) or use
% rheome.graphfilterbank.fromOperator, which does that for you. A supplied lmax is
% cross-checked against a cheap power iteration and a warning is raised if it is low.
%
% Supplies a fifth handle, .filter(g, X), which graphfilterbank prefers over the
% forward/inverse sandwich: a polynomial route has no coefficient space at all.
%
% ⚠ THE LOW-PASS MEMBER IS THIS METHOD'S WEAK SPOT. A mexhat bank's scaling function
% is a QUARTIC roll-off near lambda ~ 0.4*lambda_min, and when the operator's true
% lambda_max is orders of magnitude larger, that feature occupies a vanishing fraction
% of the fitted interval. Measured on an ico4 sphere (roll-off at 186, fit over
% [0, 4.7e5], i.e. 0.04% of the range): 35% error at order 160, 5.6% at order 320,
% while every WAVELET member is already exact to three decimals at order 160. It
% converges, but slowly -- use rheome.graphtransform.eigen when the low-pass band matters.
%
% See also: rheome.graphtransform.eigen, rheome.graphfilterbank.fromOperator
%
% Author: Diellor Basha, 2026

    if nargin < 2, lmax = []; end
    p = inputParser;
    p.addParameter('Order', 40);
    p.addParameter('Mass',  []);
    p.parse(varargin{:});
    o = p.Results;

    n  = size(L, 1);
    Mm = o.Mass;
    if isempty(Mm)
        Aop = @(X) L*X;
    else
        dM = decomposition(Mm, 'chol');       % factor once, reuse every recursion step
        Aop = @(X) dM \ (L*X);
    end

    % lmax: estimate it, or verify what was supplied. A power iteration converges to
    % lambda_max FROM BELOW, so lest is a lower bound -- lmax below it is provably wrong.
    lest = i_lmaxest(Aop, n);
    if isempty(lmax)
        lmax = 1.01 * lest;
    elseif lmax < 0.99 * lest
        warning('graphtransform:chebyshev:lmax', ...
            ['lmax = %.4g is below the operator''s spectrum (at least %.4g). The ' ...
             'Chebyshev recursion will DIVERGE, not merely lose accuracy. This is the ' ...
             'signature of reusing a truncated eigenbasis''s lambda_max. Omit lmax to ' ...
             'have it estimated, or use rheome.graphfilterbank.fromOperator.'], lmax, lest);
    end

    T = struct('forward', @(X) X, ...
               'inverse', @(X) X, ...
               'norm',    @(X) abs(X), ...
               'filter',  @(gfun, X) i_cheby(gfun, X, Aop, lmax, o.Order), ...
               'lambda',  [], ...
               'nV',      n, ...
               'rows',    n, ...
               'order',   o.Order, ...
               'lmax',    lmax);
end

function lest = i_lmaxest(Aop, n, nIter)
    % Power iteration: converges to lambda_max from below, so it is a LOWER bound.
    if nargin < 3, nIter = 25; end
    rs = RandStream('twister', 'Seed', 0);      % deterministic, and leaves rng alone
    x  = randn(rs, n, 1);  x = x / norm(x);
    lest = 0;
    for i = 1:nIter
        y = Aop(x);
        nv = norm(y);
        if nv == 0, break; end
        x = y / nv;
        lest = abs(x' * Aop(x));
    end
end

function Y = i_cheby(gfun, X, Aop, lmax, K)
    % Chebyshev coefficients of g on [0, lmax] at the Chebyshev nodes.
    a  = lmax/2;
    j  = (0:K)';
    xj = cos(pi*(j+0.5)/(K+1));                 % nodes on [-1,1]
    gj = gfun(a*(xj+1));  gj = gj(:);
    c  = zeros(K+1,1);
    for k = 0:K
        c(k+1) = (2/(K+1)) * sum(gj .* cos(pi*k*(j+0.5)/(K+1)));
    end
    c(1) = c(1)/2;

    % Clenshaw recursion on the shifted operator  Ltil = L/a - I,  which maps
    % lambda in [0, lmax] onto [-1, 1].
    Ltil = @(Y) Aop(Y)/a - Y;
    b2 = zeros(size(X));  b1 = zeros(size(X));
    for k = K:-1:1
        bk = 2*Ltil(b1) - b2 + c(k+1)*X;
        b2 = b1;  b1 = bk;
    end
    Y = Ltil(b1) - b2 + c(1)*X;
end

% Author: Diellor Basha, 2026
