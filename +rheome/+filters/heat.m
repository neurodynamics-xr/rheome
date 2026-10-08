function g = heat(lambda, t, lmax)
% FILTERS.HEAT  Heat / diffusion low-pass kernel  g(lambda) = exp(-t * lambda).
%
%   g = rheome.filters.heat(lambda, t)          % t scalar -> [K x 1];  t vector -> [K x numel(t)]
%   g = rheome.filters.heat(lambda, t, lmax)    % normalized spectrum: exp(-t * lambda / lmax)
%
% The heat kernel is the Green's function of diffusion on the surface: applied to a
% delta it spreads it into a smooth, strictly positive bump whose width grows with t.
% Low-pass (attenuates high lambda), prior-admissible (g(0)=1). Mirrors Brainstorm's
% bst_eigfilter_design_heat.
%
% INPUTS:
%   lambda [K x 1] eigenvalues (basis.Lambda)
%   t      scalar or vector of diffusion times (larger t = wider, smoother)
%   lmax   optional spectral normalization constant (e.g. max(lambda))
%
% OUTPUT:
%   g      [K x numel(t)] spectral gains, one column per scale
%
% See also: rheome.filters.mexhat, rheome.filters.localize
%
% Author: Diellor Basha, 2026

    lambda = double(lambda(:));
    t      = double(t(:)');
    if nargin >= 3 && ~isempty(lmax) && lmax > 0
        lambda = lambda / lmax;
    end
    if any(t < 0), error('filters:heat:t', 't must be >= 0.'); end
    g = exp(-(lambda * t));   % [K x numel(t)]
end

% Author: Diellor Basha, 2026
