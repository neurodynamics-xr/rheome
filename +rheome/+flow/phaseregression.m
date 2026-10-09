function out = phaseregression(z, S, opts)
% FLOW.PHASEREGRESSION  Source-space phase-gradient estimator: circular-linear regression of phase on geodesic distance in a patch.
%
%   out = rheome.flow.phaseregression(z, S, Rate=fs)
%   out = rheome.flow.phaseregression(z, S, Rate=fs, Centres=v, RadiusMM=20, Neighbours=out0.Neighbours)
%
% The comparator a reader would otherwise use (MS1 G16, estimator iii; the circular-linear regression of
% phase on distance of Zhang & Jacobs 2015 / Kempter 2012, moved onto the cortex). Per centre vertex c:
%   patch   the vertices within RadiusMM of c by GEODESIC distance (rheome.geom.edgegraph);
%   chart   each patch vertex at its geodesic distance d from c, in the direction of its displacement
%           projected onto c's tangent plane (x, y), so |(x, y)| = d -- the regression is on geodesic distance;
%   fit     the wavevector q maximising the amplitude-weighted resultant over the window,
%               R(q) = sum_t | sum_i z_i(t) exp(-i q.(x_i, y_i)) | / sum_t sum_i |z_i(t)|,
%           on a grid of LambdaMM x NumDirections (plus q = 0), refined by a parabola in log wavelength
%           and in angle. R in [0, 1] is the fit's coherence (1 = one plane wave across the patch).
%   speed   omega / |k|, omega the patch's mean phase advance per sample x Rate, k = -q (for
%           z = A exp(i(omega t - k.x)) the resultant peaks at q = -k), velocity = omega k / |k|^2.
% A best fit at q = 0 (no spatial phase gradient: a standing or a stationary pattern) returns speed Inf
% and velocity 0.
%
% INPUTS
%   z   [nV x nT] complex analytic field (band-limited), S surface (.Vertices .Faces .VertNormals)
%   Rate        Hz (required)
%   Centres     vertex indices to fit at (default every vertex)
%   RadiusMM    patch radius (default 20)
%   Neighbours  sparse [nV x nV], (i, c) = geodesic distance + 1e-12 (m) for i within RadiusMM of c; reuse
%               out.Neighbours across calls on the same mesh (it costs one Dijkstra per new centre)
%   LambdaMM    wavelength grid (default 40 log-spaced, 8-1000 mm)   NumDirections (default 24)
% OUTPUT out: .centres [nC x 1]  .k .velocity [nC x 3] (1/m, m/s)  .speed .wavelength (mm) .omega (rad/s)
%             .R [nC x 1]  .Neighbours
%
% ⚠ The grid floors the wavelength at LambdaMM(1); below ~4 mesh edges per wavelength the patch aliases,
%   as rheome.flow.phasegradient does (its .aliased).
%
% See also: rheome.flow.phasegradient, rheome.flow.bstopticalflow, rheome.scale.measure_catalogue
%
% Author: Diellor Basha, 2026

    arguments
        z {mustBeNumeric}
        S struct
        opts.Rate (1,1) double {mustBePositive}
        opts.Centres double = []
        opts.RadiusMM (1,1) double {mustBePositive} = 20
        opts.Neighbours = []
        opts.LambdaMM double = logspace(log10(8), log10(1000), 40)
        opts.NumDirections (1,1) double {mustBeInteger, mustBePositive} = 24
    end
    V = double(S.Vertices);  nV = size(V, 1);
    assert(size(z, 1) == nV, 'flow:phaseregression:size', 'z has %d rows, the mesh %d vertices', size(z, 1), nV);
    N = S.VertNormals ./ max(vecnorm(S.VertNormals, 2, 2), eps);
    c = opts.Centres(:);  if isempty(c), c = (1:nV)'; end
    D = opts.Neighbours;  if isempty(D), D = sparse(nV, nV); end
    todo = c(~any(D(:, c), 1)');
    if ~isempty(todo), D = D + i_neighbours(S, todo, opts.RadiusMM*1e-3, nV); end

    lam = opts.LambdaMM(:) * 1e-3;  nL = numel(lam);  nD = opts.NumDirections;  ang = (0:nD-1) * 2*pi/nD;
    km = 2*pi ./ lam;  kx = km * cos(ang);  ky = km * sin(ang);              % [nL x nD]
    nC = numel(c);  K = nan(nC, 3);  Rr = nan(nC, 1);  om = nan(nC, 1);
    for j = 1:nC
        [P, ~, d] = find(D(:, c(j)));  d = d - 1e-12;
        n = N(c(j), :);  e1 = cross(n, [1 0 0]);  if norm(e1) < 0.1, e1 = cross(n, [0 1 0]); end
        e1 = e1 / norm(e1);  e2 = cross(n, e1);
        u = V(P, :) - V(c(j), :);  x = u * e1';  y = u * e2';  rho = hypot(x, y);
        s = d ./ max(rho, eps);  s(rho < eps) = 0;  x = x .* s;  y = y .* s;
        Z = z(P, :);  a = sum(abs(Z(:)));
        if a == 0, continue, end
        sc = reshape(sum(abs(exp(-1i * (x * kx(:)' + y * ky(:)')).' * Z), 2), nL, nD) / a;
        s0 = sum(abs(sum(Z, 1))) / a;
        [Rb, ib] = max(sc(:));  [iL, iD] = ind2sub([nL nD], ib);
        om(j) = angle(sum(Z(:, 2:end) .* conj(Z(:, 1:end-1)), 'all')) * opts.Rate;
        if s0 >= Rb, K(j, :) = 0;  Rr(j) = s0;  continue, end
        ll = log(lam(iL));
        if iL > 1 && iL < nL, ll = ll + i_vertex(sc(iL-1:iL+1, iD)) * (log(lam(iL+1)) - log(lam(iL))); end
        dd = mod(iD + [-2 0], nD) + 1;                                        % circular neighbours
        th = ang(iD) + i_vertex(sc(iL, [dd(1) iD dd(2)])) * 2*pi/nD;
        q = 2*pi / exp(ll) * (cos(th) * e1 + sin(th) * e2);
        K(j, :) = -q;  Rr(j) = Rb;
    end
    k2 = sum(K.^2, 2);
    vel = om .* K ./ k2;  vel(k2 == 0, :) = 0;
    sp = abs(om) ./ sqrt(k2);
    out = struct('centres', c, 'k', K, 'velocity', vel, 'speed', sp, 'wavelength', 2*pi ./ sqrt(k2) * 1e3, ...
                 'omega', om, 'R', Rr, 'Neighbours', D);
end

function d = i_vertex(y)
% offset of a parabola's vertex through three equally spaced samples, in samples, clamped to +-0.5
    y = y(:);  den = y(1) - 2*y(2) + y(3);  d = 0;
    if den < 0, d = max(min(0.5 * (y(1) - y(3)) / den, 0.5), -0.5); end
end

function D = i_neighbours(S, c, R, nV)
    g = rheome.geom.edgegraph(S);  I = cell(numel(c), 1);  J = I;  W = I;
    for j = 1:numel(c)
        dj = distances(g, c(j))';  k = find(dj <= R);
        I{j} = k;  J{j} = repmat(c(j), numel(k), 1);  W{j} = dj(k) + 1e-12;
    end
    D = sparse(vertcat(I{:}), vertcat(J{:}), vertcat(W{:}), nV, nV);
end

% Author: Diellor Basha, 2026
