function G = kleingordon(lambda, t, c, m, varargin)
% FILTERS.KLEINGORDON  Massive (Klein-Gordon) wave kernel (eigenmode x time lag, 'ts').
%
%   G = rheome.filters.kleingordon(lambda, t, c, m)                       % continuum (default)
%   G = rheome.filters.kleingordon(lambda, t, c, m, 'Variant','discrete', 'LambdaMax', lmax)
%
% CONTINUUM (the definition). u_tt = -c^2 L u - m^2 u gives each mode
%       omega_k^2 = c^2 * lambda_k + m^2,      g(lambda,t) = cos( t * sqrt(c^2*lambda + m^2) ).
% m is the REST FREQUENCY (rad/s): even the lambda = 0 mode oscillates at m, and the dispersion
% has a frequency floor. PHYSICAL UNITS (lambda 1/length^2, c length/s, t s), not normalised by lambda_max;
% m enters ONLY through omega^2 (never inside an inverse cosine).
%
% DISCRETE (named variant, never the default): the earlier Chebyshev-angle form
% cos( c * t * acos(1 - 2*(lambda + m)/lmax) ), where the mass shifts the Laplacian's eigenvalue
% before the change of variable (a shifted Laplacian, not a mass term). Port of
% bst_eigfilter_design_kleingordon; there m is the shift "mu" in eigenvalue units.
%
% A (lambda,t) kernel acts on a time LAG; it cannot multiply a joint (lambda,f) spectrum.
%
% INPUTS:
%   lambda [K x 1] (1/length^2)   t [1 x N] time-lag axis (s)
%   c      speed (length/s, default 1)   m rest frequency (rad/s, default 0.1)
% NAME-VALUE:
%   'Variant'    'continuum' (default) | 'discrete'
%   'LambdaMax'  normaliser of the discrete variant (default max(lambda))
% OUTPUT:
%   G      [K x N]
%
% See also: rheome.filters.wave, rheome.filters.dampedwave
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(c), c = 1;   end
    if nargin < 4 || isempty(m), m = 0.1; end
    if ~isempty(varargin) && isnumeric(varargin{1})
        error('filters:kleingordon:arguments', ...
            ['The fifth argument is no longer lmax: the continuum form needs no normaliser. ' ...
             'Use ''Variant'',''discrete'',''LambdaMax'',lmax for the Chebyshev-angle form.']);
    end
    p = inputParser;
    p.addParameter('Variant', 'continuum', @(v) any(strcmpi(v, {'continuum','discrete'})));
    p.addParameter('LambdaMax', []);
    p.parse(varargin{:});
    lam = double(lambda(:));  tt = double(t(:)');
    switch lower(p.Results.Variant)
        case 'continuum'
            G = cos( sqrt(c^2 * max(lam, 0) + m^2) * tt );
        case 'discrete'
            lmax = p.Results.LambdaMax;  if isempty(lmax), lmax = max(lam); end
            sc = 2 / max(lmax, eps);
            th = acos( max(min(1 - sc*(lam + m), 1), -1) );
            G  = cos( c * (th * tt) );
    end
end

% Author: Diellor Basha, 2026
