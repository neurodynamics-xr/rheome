function g = mexhat(lambda, t)
% FILTERS.MEXHAT  Mexican-hat band-pass wavelet  g(lambda) = (t*lambda) .* exp(-t*lambda).
%
%   g = rheome.filters.mexhat(lambda, t)   % t scalar -> [K x 1];  t vector -> [K x numel(t)] bank
%
% The Mexican-hat wavelet is band-pass and zero-mean (g(0)=0), so applied to a delta
% it produces a SIGNED atom: a positive core ringed by a negative surround, at a
% spatial scale set by t. Sweeping t gives a wavelet filterbank at multiple scales.
% Mirrors Brainstorm's bst_eigfilter_design_mexhat.
%
% INPUTS:
%   lambda [K x 1] eigenvalues (basis.Lambda)
%   t      scalar or vector of scales (larger t selects LOWER frequencies -> wider atom)
%
% OUTPUT:
%   g      [K x numel(t)] spectral gains, one column per scale
%
% See also: rheome.filters.heat, rheome.filters.localize
%
% Author: Diellor Basha, 2026

    lambda = double(lambda(:));
    t      = double(t(:)');
    if any(t < 0), error('filters:mexhat:t', 't must be >= 0.'); end
    tl = lambda * t;          % [K x numel(t)]
    g  = tl .* exp(-tl);
end

% Author: Diellor Basha, 2026
