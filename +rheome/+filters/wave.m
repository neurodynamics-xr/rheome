function G = wave(lambda, t, c, varargin)
% FILTERS.WAVE  Wave propagation kernel (eigenmode x time lag, 'ts' domain).
%
%   G = rheome.filters.wave(lambda, t, c)                        % continuum (default)
%   G = rheome.filters.wave(lambda, t, c, 'Variant','discrete', 'LambdaMax', lmax)
%
% CONTINUUM (the definition). The d'Alembert wave equation u_tt = -c^2 L u on a surface
% with Laplacian eigenvalue lambda gives each mode the frequency  omega_k = c*sqrt(lambda_k),
% so the propagator is
%       g(lambda,t) = cos( t * c * sqrt(lambda) ).
% In PHYSICAL UNITS: lambda in 1/length^2, c in length/s, t in seconds, omega in rad/s. Nothing
% is normalised by lambda_max. On the unit sphere omega_l = c*sqrt(l(l+1))/R, the exact spectrum
% (see rheome.filters.sphere_kernel).
%
% DISCRETE (a named variant, never the default). The Chebyshev change of variable
% theta = acos(1 - 2*lambda/lmax) turns a polynomial filter on [0,lmax] into a cosine
% series, and  g = cos( c * t * theta )  is then a discrete-time propagator whose "frequency"
% c*theta is bounded (theta in [0,pi]) and is NOT c*sqrt(lambda): lambda is normalised by lmax
% and the time unit is relative to it. Port of bst_eigfilter_design_wave. Pass 'LambdaMax'
% (default max(lambda)).
%
% A (lambda,t) kernel acts on a time LAG: it shapes the response to an impulse as a function of
% elapsed time. It cannot multiply a joint (lambda,f) spectrum; transform it along the time
% axis first (rheome.filters.jtv) before using it as a 'js' kernel.
%
% INPUTS:
%   lambda [K x 1] eigenvalues (1/length^2)   t [1 x N] time-lag axis (s)
%   c      wave speed (length/s, default 1)
% NAME-VALUE:
%   'Variant'    'continuum' (default) | 'discrete'
%   'LambdaMax'  spectral normaliser of the discrete variant (default max(lambda)); ignored by 'continuum'
% OUTPUT:
%   G      [K x N]
%
% See also: rheome.filters.dampedwave, rheome.filters.kleingordon, rheome.filters.sphere_kernel, rheome.filters.jtv
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(c), c = 1; end
    [variant, lmax] = i_options(lambda, 'wave', varargin);
    lam = double(lambda(:));  tt = double(t(:)');
    switch variant
        case 'continuum', G = cos( (c * sqrt(max(lam, 0))) * tt );
        case 'discrete',  G = cos( c * (i_theta(lam, lmax, 0) * tt) );
    end
end

function [variant, lmax] = i_options(lambda, who, args)
    if ~isempty(args) && isnumeric(args{1})
        error(['filters:' who ':arguments'], ...
            ['The fourth argument is no longer lmax: the continuum form needs no normaliser. ' ...
             'Use ''Variant'',''discrete'',''LambdaMax'',lmax for the Chebyshev-angle form.']);
    end
    p = inputParser;
    p.addParameter('Variant', 'continuum', @(v) any(strcmpi(v, {'continuum','discrete'})));
    p.addParameter('LambdaMax', []);
    p.parse(args{:});
    variant = lower(p.Results.Variant);  lmax = p.Results.LambdaMax;
    if isempty(lmax), lmax = max(lambda(:)); end
end

function th = i_theta(lam, lmax, shift)
    sc = 2 / max(lmax, eps);
    th = acos( max(min(1 - sc*(lam + shift), 1), -1) );
end

% Author: Diellor Basha, 2026
