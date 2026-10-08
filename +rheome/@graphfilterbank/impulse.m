function u = impulse(obj, m, vertex)
% IMPULSE  Vertex-domain impulse response of one member at a vertex. The atom.
%
%   u = impulse(gfb, m, vertex)      -> [nV x 1]
%
% Seeds a delta at the vertex and pushes it through member m: g_m(L) applied to
% delta_vertex. An impulse is a delta on whatever index set the operator lives on, so
% this is the SAME operation jointfilterbank/impulse performs on the vertex-time product
% graph, and the same one cwtfilterbank/wavelets performs on a path or ring. Only the
% graph differs.
%
% ⚠ A GRAPH HAS NO TRANSLATION INVARIANCE, so there is no single wavelet per scale -- the
% response genuinely differs at every vertex and a location is REQUIRED. There is no
% canonical centring either.
%
% Equivalent to graphfilters(gfb,'Type','vertex','Vertex',vertex) taken at member m. Reach
% for that form to seed several vertices at once; it returns [nV x P x M].
%
% See also: graphfilters, gain, widths, rheome.jointfilterbank/impulse
%
% Author: Diellor Basha, 2026

    if ~isscalar(m) || m < 1 || m > obj.NumMembers || mod(m,1) ~= 0
        error('graphfilterbank:member', ...
            'member must be an integer in 1..%d, got %s.', obj.NumMembers, mat2str(m));
    end

    T  = gfb_transform(obj);
    nV = T.rows;
    if ~isscalar(vertex) || vertex < 1 || vertex > nV || mod(vertex,1) ~= 0
        error('graphfilterbank:vertex', ...
            'vertex must be an integer in 1..%d, got %s.', nV, mat2str(vertex));
    end

    x = zeros(nV, 1);  x(vertex) = 1;

    % ⚠ TWO ROUTES, AND ONLY ONE OF THEM HAS A SPECTRAL DOMAIN. rheome.graphtransform.chebyshev
    % has IDENTITY forward/inverse and does its work through .filter, by polynomial
    % recursion in L. Reaching for forward/inverse there hands back the delta unchanged and
    % multiplies it by a gain vector -- silently wrong rather than an error, and the atom
    % comes out all but zero.
    if isfield(T, 'filter') && ~isempty(T.filter)
        u = T.filter(obj.G_{m}, x);
        return;
    end

    % One member, one inverse transform. The [K x M] gain table is scalar evaluations and
    % is cheap; materialising all M atoms to return one would not be.
    H = graphfilters(obj, 'Lambda', gfb_lambda(obj, T));
    u = T.inverse(H(:, m) .* T.forward(x));
end

% Author: Diellor Basha, 2026
