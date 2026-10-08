function g = diffgauss(lambda, t1, t2)
% FILTERS.DIFFGAUSS  Difference-of-Gaussians band-pass  g(lambda) = exp(-t1*lambda) - exp(-t2*lambda).
%
%   g = rheome.filters.diffgauss(lambda, t1, t2)   % t1,t2 scalar -> [K x 1]; equal-length vectors -> [K x n] bank
%
% A difference of two heat kernels (t1 < t2): band-pass and zero at DC (g(0)=0), so applied to a
% delta it produces a SIGNED atom -- a positive core with a negative surround -- like the mexhat,
% but built from two Gaussians (the classic DoG, a Ricker/Laplacian-of-Gaussian approximation).
% The pass band sits between the two scales; larger t = wider, lower-frequency atom. Mirrors
% Brainstorm's bst_eigfilter_design_diffgauss.
%
% INPUTS:
%   lambda  [K x 1] eigenvalues (basis.Lambda)
%   t1, t2  the two Gaussian scales with t1 < t2; scalars, or equal-length vectors for a bank
%
% OUTPUT:
%   g       [K x n] spectral gains (n = numel(t1)), one column per (t1,t2) pair
%
% See also: rheome.filters.mexhat, rheome.filters.heat, rheome.filters.localize
%
% Author: Diellor Basha, 2026 (port of bst_eigfilter_design_diffgauss)

    lambda = double(lambda(:));
    t1 = double(t1(:)');  t2 = double(t2(:)');
    if numel(t1) ~= numel(t2)
        error('filters:diffgauss:size', 't1 and t2 must have the same number of elements.');
    end
    if any(t1 < 0), error('filters:diffgauss:t', 't1 must be >= 0.'); end
    if any(t1 >= t2), error('filters:diffgauss:order', 'require t1 < t2 (band-pass, non-negative gain).'); end
    g = exp(-(lambda * t1)) - exp(-(lambda * t2));   % [K x n]
end

% Author: Diellor Basha, 2026
