function G = gabor(lambda, f, k0, f0, sf, half)
% FILTERS.GABOR  Spatiotemporal Gabor packet (eigenmode x frequency, 'js' domain).
%
%   G = rheome.filters.gabor(lambda, f, k0, f0, sf)
%   G = rheome.filters.gabor(lambda, f, k0, f0, sf, half)
%
% Separable band-pass in BOTH axes: a Gaussian around spatial frequency k0 (= sqrt(lambda))
% times a Gaussian around temporal frequency +/- f0:
%   g(lambda,f) = exp(-(sqrt(lambda)-k0)^2 / 2*sk^2) *
%                 [ exp(-(f-f0)^2/2*sf^2) + exp(-(f+f0)^2/2*sf^2) ],  sk = k0/3.
% UNITS: f, f0 and sf are in HERTZ (f signed); k0 in 1/length (lambda in 1/length^2).
% Selects a spatial-scale x temporal-frequency packet -- a "note" in the joint spectrum.
% Port of bst_eigfilter_design_gabor (Hermitian -> real impulse response).
%
% INPUTS:
%   lambda [K x 1]   f [1 x N] frequency axis (Hz); a full FFT grid unless half = 'positive'
%   k0     centre spatial frequency sqrt(lambda) (default 209)   f0 centre temporal freq (Hz, default 10)
%   sf     temporal bandwidth (Hz, default 2)
%   half   'full' (default) | 'positive': see rheome.filters.signedfrequencies
% OUTPUT:
%   G      [K x N]
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(k0), k0 = 209; end
    if nargin < 4 || isempty(f0), f0 = 10;  end
    if nargin < 5 || isempty(sf), sf = 2;   end
    sk = max(k0/3, eps);
    if nargin < 6 || isempty(half), half = 'full'; end
    ws = rheome.filters.signedfrequencies(f, half);
    sl = exp( -((sqrt(max(double(lambda(:)),0)) - k0).^2) / (2*sk^2) );          % [K x 1]
    gt = exp( -((ws - f0).^2) / (2*sf^2) ) + exp( -((ws + f0).^2) / (2*sf^2) );  % [1 x N]
    G  = sl * gt;
end


% Author: Diellor Basha, 2026
