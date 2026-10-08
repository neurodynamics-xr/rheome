function G = graph(arr, varargin)
% SENSORS.GRAPH  A Gaussian-weighted kNN graph on sensor coordinates -- the GSPBox operator.
%
%   G = rheome.sensors.graph(arr)
%   G = rheome.sensors.graph(arr, 'K', 8, 'Sigma', 4e-4, 'Laplacian', 'combinatorial')
%
%   w_ij = exp(-d_ij^2 / (2*sigma^2))  over the K nearest neighbours,  W = max(W, W'),
%   L = D - W.
%
% ⭐ IT RETURNS THE SURFACE SHAPE. .Vertices, .Faces and .nV are populated exactly as the
% mesh code expects, so rheome.flow.phasegradient, rheome.detect.phasesingularity and rheome.flow.phaselock
% accept a sensor array with no adapter and no special case.
%
% ⚠ SYMMETRISE BY max, NOT BY AVERAGE. kNN is not a symmetric relation -- a boundary sensor
% is often a neighbour of an interior one without the reverse holding. max(W,W') keeps such
% an edge at full weight; averaging halves it, and since alpha (rheome.sensors.calibrate) is linear
% in the weights that biases the calibration low exactly where the array is most irregular.
%
% ⚠ lambda IS DIMENSIONLESS HERE. This operator carries no length units of its own; they
% come from rheome.sensors.calibrate, and no speed may be reported without it.
%
% ⚠ THE NORMALIZED LAPLACIAN BREAKS THE CLOSED FORM. D^-1/2 W D^-1/2 rescales each vertex by
% its own degree, so sum_j w_ij |d_ij|^2 is no longer the quantity that appears in the
% small-k expansion. It is available and its spectrum is bounded by 2, but rheome.sensors.calibrate
% refuses it.
%
% INPUTS:
%   arr  a rheome.sensors.* array struct
%   'K'          neighbours (default 4 for Dim 1, 8 for Dim 2)
%   'Sigma'      kernel width, metres (default arr.Pitch)
%   'Laplacian'  'combinatorial' (default) | 'normalized'
%   'Faces'      build the readout triangulation (default true)
% OUTPUT:
%   G  .Vertices .Faces .nV .W .L .D .Sigma .K .Laplacian .Array
%
% See also: rheome.sensors.modes, rheome.sensors.calibrate, rheome.graphtransform.eigen
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('K',         []);
    p.addParameter('Sigma',     []);
    p.addParameter('Laplacian', 'combinatorial');
    p.addParameter('Faces',     true);
    p.parse(varargin{:});
    o = p.Results;

    lap = lower(char(o.Laplacian));
    if ~any(strcmp(lap, {'combinatorial', 'normalized'}))
        error('sensors:graph:laplacian', ...
            'Laplacian must be ''combinatorial'' or ''normalized'', got ''%s''.', lap);
    end

    n = arr.nCh;
    K = o.K;      if isempty(K),     K     = 4 * arr.Dim;              end
    sg = o.Sigma; if isempty(sg),    sg    = arr.Pitch;                end
    K  = min(double(K), n - 1);
    if ~(sg > 0)
        error('sensors:graph:sigma', 'Sigma must be positive, got %g.', sg);
    end

    [~, ~, D] = sen_metrics(arr.Pos);

    % kNN mask: the K smallest strictly-positive distances per row.
    Dnn = D;  Dnn(1:n+1:end) = Inf;
    [~, ord] = sort(Dnn, 2, 'ascend');
    rows = repmat((1:n).', 1, K);
    cols = ord(:, 1:K);
    val  = exp(-D(sub2ind([n n], rows(:), cols(:))).^2 / (2*sg^2));
    W    = sparse(rows(:), cols(:), val, n, n);
    W    = max(W, W.');                       % NOT (W+W')/2 -- see the header
    W(1:n+1:end) = 0;

    deg = full(sum(W, 2));
    switch lap
        case 'combinatorial'
            L = spdiags(deg, 0, n, n) - W;
        case 'normalized'
            di = 1 ./ sqrt(max(deg, eps));
            L  = speye(n) - spdiags(di, 0, n, n) * W * spdiags(di, 0, n, n);
    end
    L = (L + L.') / 2;                        % kill the last bit of round-off asymmetry

    G = struct();
    G.Vertices  = arr.Pos;
    G.Faces     = [];
    G.nV        = n;
    G.W         = W;
    G.L         = L;
    G.D         = D;
    G.Sigma     = sg;
    G.K         = K;
    G.Laplacian = lap;
    G.Array     = arr;

    if o.Faces
        G.Faces = sen_faces(arr);
    end
end

% Author: Diellor Basha, 2026
