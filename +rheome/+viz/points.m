function P = points(cp, S, opts)
% VIZ.POINTS  Critical points (vortices, sources, sinks, saddles) as viewer point sets, one set a type.
%
%   P = rheome.viz.points(rheome.flow.critical(dbasis, a, S, frames), S)       % at their frames
%   P = rheome.viz.points(rheome.detect.criticalPoints(J, S), S)               % one frame: shown at every frame
%   rheome.view(S, struct('J', J), Points=P)
%
% Each point goes to its nearest vertex. One set a type, so the marks sidecar keys each type (cortex.marks 'critical
% points': "each type"); sources and sinks in flow-source, vortices in flow-vortex (the figure tokens), saddles in the
% viewer's default; the sign is the charge.
%
% INPUTS:
%   cp  a rheome.detect.criticalPoints struct (.pos .type .charge), or rheome.flow.critical's output (.cp .frames)
%   S   the surface the points lie on (.Vertices)
% OPTIONS:
%   Name    the sets' prefix (default "critical"): sets "<Name>: <type>"
%   Tokens  the figure tokens' JSON (default $NXR_FIGURE_TOKENS, else ~/workspace/library/styles/tokens/figure-tokens.json;
%           missing: no colour)
%
% OUTPUT:
%   P : struct array .name .vertices (1-based) .frames (1-based; 0: every frame) .signs .color, for
%       rheome.io.write.bundle 'Points'
%
% NOT YET (the viewer): the tokens' glyphs (+ source, - sink, open circle vortex); the viewer draws plus/minus badges.
%
% Example:
%   S = struct('Vertices', eye(3), 'Faces', [1 2 3]);
%   cp = struct('pos', [0.9 0 0; 0 0 1], 'type', {{'vortex'; 'source'}}, 'charge', [1; 1]);
%   P = rheome.viz.points(cp, S);  assert(isequal([P.vertices], [1 3]))
%
% See also: rheome.flow.critical, rheome.detect.criticalPoints, rheome.io.write.bundle
%
% Author: Diellor Basha, 2026
    arguments
        cp struct
        S struct
        opts.Name = "critical"
        opts.Tokens = ""
    end
    if isfield(cp, 'cp'), C = cp.cp; fr = cp.frames; else, C = cp; fr = zeros(1, numel(cp)); end
    pos = zeros(0, 3); ty = {}; q = zeros(0, 1); f = zeros(0, 1);
    for k = 1:numel(C)
        n = size(C(k).pos, 1);
        pos = [pos; C(k).pos]; ty = [ty; cellstr(C(k).type(:))]; q = [q; C(k).charge(:)]; f = [f; repmat(fr(k), n, 1)]; %#ok<AGROW>
    end
    v = zeros(size(pos, 1), 1);
    for i = 1:numel(v), [~, v(i)] = min(sum((S.Vertices - pos(i, :)).^2, 2)); end
    col = colours(opts.Tokens);
    P = struct('name', {}, 'vertices', {}, 'frames', {}, 'signs', {}, 'color', {});
    for t = unique(lower(ty), 'stable')'
        k = strcmpi(ty, t{1});
        c = '';  if isfield(col, t{1}), c = col.(t{1}); end
        P(end+1) = struct('name', sprintf('%s: %s', opts.Name, t{1}), 'vertices', v(k), 'frames', f(k), ...
            'signs', q(k), 'color', c); %#ok<AGROW>
    end
end

function c = colours(file)
% the figure tokens' concept colours for each type (sources and sinks: flow-source; vortices: flow-vortex)
    c = struct();
    file = char(file);
    if isempty(file), file = getenv('NXR_FIGURE_TOKENS'); end
    if isempty(file), file = fullfile(getenv('HOME'), 'workspace', 'library', 'styles', 'tokens', 'figure-tokens.json'); end
    if ~isfile(file), return; end
    s = jsondecode(fileread(file)).colour.semantic;
    hex = @(k) s.(matlab.lang.makeValidName(k)).x_value;
    c = struct('source', hex('flow-source'), 'sink', hex('flow-source'), 'vortex', hex('flow-vortex'));
end
