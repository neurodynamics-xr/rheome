function M = modes(G, varargin)
% SENSORS.MODES  The sensor graph's spectrum, and the two ways to filter with it.
%
%   M = rheome.sensors.modes(G)                          % exact eigenbasis
%   M = rheome.sensors.modes(G, 'Route', 'chebyshev', 'Order', 40)
%
% ⭐ ON A SENSOR ARRAY THE EXACT ROUTE IS FREE. A few hundred nodes means a dense eig costs
% milliseconds, so the truncation that a large mesh is forced into does not arise here
% and the eigenbasis is the default. Chebyshev is carried anyway because it is the route
% that SCALES, and having both on the same operator is what lets a demo show they agree --
% which is the claim that transfers to a graph of tens of thousands of nodes, where only
% one of them is affordable.
%
% INPUTS:
%   G  a rheome.sensors.graph struct
%   'Route'  'eigen' (default) | 'chebyshev'
%   'Order'  Chebyshev order (default 40)
% OUTPUT:
%   M  .Lambda [N x 1] ascending   .Phi [N x N] orthonormal ([] for chebyshev)
%      .T  a graphtransform struct  .Lmax  .Route
%
% See also: rheome.sensors.graph, rheome.sensors.calibrate, rheome.graphtransform.eigen, rheome.graphtransform.chebyshev
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('Route', 'eigen');
    p.addParameter('Order', 40);
    p.parse(varargin{:});
    o = p.Results;

    route = lower(char(o.Route));
    if ~any(strcmp(route, {'eigen', 'chebyshev'}))
        error('sensors:modes:route', ...
            'Route must be ''eigen'' or ''chebyshev'', got ''%s''.', route);
    end

    n = G.nV;
    L = full(G.L);
    L = (L + L.') / 2;

    switch route
        case 'eigen'
            [Phi, Lam] = eig(L, 'vector');
            [Lam, ix]  = sort(real(Lam), 'ascend');
            Phi        = Phi(:, ix);
            Lam(1)     = 0;                          % the constant mode, exactly
            M = struct('Lambda', Lam, 'Phi', Phi, 'Lmax', max(Lam), 'Route', route, ...
                       'T', rheome.graphtransform.eigen(Phi, speye(n), Lam));

        case 'chebyshev'
            Lam  = sort(real(eig(L)), 'ascend');
            lmax = 1.01 * max(Lam);                  % must BOUND the spectrum, not equal it
            M = struct('Lambda', Lam, 'Phi', [], 'Lmax', lmax, 'Route', route, ...
                       'T', rheome.graphtransform.chebyshev(sparse(L), lmax, 'Order', o.Order));
    end
end

% Author: Diellor Basha, 2026
