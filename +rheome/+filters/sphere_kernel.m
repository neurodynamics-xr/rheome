function K = sphere_kernel(gamma, gl)
% FILTERS.SPHERE_KERNEL  Analytic response of a spectral filter on the unit sphere.
%
%   K = rheome.filters.sphere_kernel(gamma, gl)
%
% Closed-form (Legendre-series) response of a rotationally-invariant filter to a point
% source on the UNIT sphere, evaluated at geodesic angle gamma from the source. By the
% spherical-harmonic addition theorem, a filter with gain g_l on degree l (eigenvalue
% lambda_l = l(l+1)) gives
%       K(gamma) = sum_{l=0}^{L} g_l * (2l+1)/(4*pi) * P_l(cos gamma).
% This is the ANALYTIC ground truth our numerical eigenfiltering is compared against:
%   heat    g_l = exp(-l(l+1) t)         (diffusion)
%   wave    g_l = cos(sqrt(l(l+1)) c t)  (propagation at speed c)
%   Poisson g_l = 1/(l(l+1)),  g_0 = 0   (inverse Laplacian / Green's function)
%
% INPUTS:
%   gamma  [N x 1] geodesic angles from the source (radians), = acos(p . q)
%   gl     [L+1 x 1] filter gains per spherical-harmonic degree l = 0,1,...,L
%
% OUTPUT:
%   K      [N x 1] analytic field value at each angle
%
% See also: rheome.filters.heat, rheome.filters.wave, demos.sphere_pde
%
% Author: Diellor Basha, 2026

    gamma = gamma(:);
    gl    = gl(:);
    L     = numel(gl) - 1;
    x     = cos(gamma);

    % Legendre P_l(x) by the standard three-term recurrence, accumulating the series.
    Plm1 = ones(size(x));                       % P_0
    K    = gl(1) * (1/(4*pi)) * Plm1;
    if L >= 1
        Pl = x;                                 % P_1
        K  = K + gl(2) * (3/(4*pi)) * Pl;
        for l = 2:L
            Pnew = ((2*l-1) * x .* Pl - (l-1) * Plm1) / l;
            K    = K + gl(l+1) * ((2*l+1)/(4*pi)) * Pnew;
            Plm1 = Pl;  Pl = Pnew;
        end
    end
end

% Author: Diellor Basha, 2026
