classdef tSnapshot < matlab.unittest.TestCase
% TSNAPSHOT  rheome.snapshot: the shot list is checked against the figure tokens before anything is drawn (a stand-in
% viewer folder, no browser); with the real viewer, Chrome and Node (RHEOME_VIEWER or the nxr-cortical-flow clone,
% skipped otherwise), a real cortex is rendered and laid out, the files are print-sized, and the scale bar is drawn
% where the viewer draws its vertices (the camera read from the render is checked end to end).
%
% Author: Diellor Basha, 2026
    properties
        S = struct('Vertices', [0 0 0; 1 0 0; 0 1 0; 0 0 1] * 0.05, 'Faces', [1 2 3; 1 3 4; 1 4 2; 2 4 3]);
        Fake
    end
    methods (TestMethodSetup)
        function fakeViewer(tc)
        % a viewer folder with the files read before rendering: viewer.html and the default colormap library
            d = tempname;  mkdir(fullfile(d, 'colormaps', 'default'));
            fclose(fopen(fullfile(d, 'viewer.html'), 'w'));
            names = {'viridis', 'cividis', 'PuOr_r', 'magma'};
            j = struct('version', 2, 'library', 'default', 'width', 256, 'rows', 4, ...
                'colormaps', struct('name', names, 'index', {0, 1, 2, 3}));
            fid = fopen(fullfile(d, 'colormaps', 'default', 'colormaps.json'), 'w');  fwrite(fid, jsonencode(j));  fclose(fid);
            imwrite(uint8(repmat(0:255, 4, 1, 3)), fullfile(d, 'colormaps', 'default', 'colormaps.png'), 'Alpha', uint8(255 * ones(4, 256)));
            tc.Fake = d;
            tc.addTeardown(@() rmdir(d, 's'));
        end
    end
    methods (Test)
        function refusesBadShots(tc)
            r = struct('a', [-1; 0; 1; 2] * ones(1, 4), 'v', zeros(12, 4));
            run = @(shot) rheome.snapshot(tc.S, r, Shots=shot, Out=tempname, Viewer=tc.Fake, Node='/nonexistent', Chrome=fullfile(tc.Fake, 'viewer.html'));
            tc.verifyError(@() run(struct('Field', 'a')), 'rheome:snapshot:unit');                 % a colourbar needs a unit
            tc.verifyError(@() run(struct('Field', 'b', 'Unit', 'mm')), 'rheome:snapshot:option');
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'Wrong', 1)), 'rheome:snapshot:shot');
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'View', 'sideways')), 'rheome:snapshot:option');
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'Frames', [1 2])), 'rheome:snapshot:shot');      % a strip: one view
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'Frames', [1 2 4], 'View', 'L lateral')), 'rheome:snapshot:shot');   % unequal
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'Frames', 5)), 'rheome:snapshot:shot');
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'CLim', [-1 2])), 'rheome:snapshot:clim');      % signed: symmetric
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'Quantity', 'coupling')), 'rheome:snapshot:colormap');   % no BrBG_r
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'Vector', 'a')), 'rheome:snapshot:option');
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'Name', 'a/b')), 'rheome:snapshot:shot');
            tc.verifyError(@() rheome.snapshot(tc.S, r, Viewer=tc.Fake), 'rheome:snapshot:input');                 % Out
            tc.verifyError(@() rheome.snapshot(tc.S, r, Out=tempname, Viewer=tempname), 'rheome:snapshot:noviewer');
            tc.verifyError(@() rheome.snapshot(tc.S, r, Out=tempname, Viewer=tc.Fake, Tokens=[tempname '.json']), 'rheome:snapshot:tokens');
            in = struct('Vertices', [1 2 3], 'Label', 'a tile');
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'Inset', rmfield(in, 'Label'))), 'rheome:snapshot:inset');   % say what it shows
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'Inset', setfield(in, 'Vertices', [1 9]))), 'rheome:snapshot:inset');
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'Inset', setfield(in, 'Zoom', 1.5))), 'rheome:snapshot:inset');
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'Inset', setfield(in, 'View', 'top'))), 'rheome:snapshot:inset');
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'Frames', 1:3, 'View', 'L lateral', 'Inset', in)), 'rheome:snapshot:inset');
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'Inset', setfield(in, 'Field', 'zz'))), 'rheome:snapshot:option');
            sel = struct('Vertices', [1 2 3], 'Label', 'tile 1');
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'Selection', rmfield(sel, 'Label'))), 'rheome:snapshot:selection');
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'Selection', struct('Centre', 1, 'Label', 'd'))), 'rheome:snapshot:selection');   % no radius
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'Selection', setfield(sel, 'Axis', 'depth'))), 'rheome:snapshot:option');
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'Paths', struct('Points', tc.S.Vertices(1:2, :), 'Label', 'p'))), 'rheome:snapshot:path');   % How
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'Paths', struct('From', 1, 'To', 9, 'Label', 'p'))), 'rheome:snapshot:path');
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'Inset', setfield(in, 'Shape', 'star'))), 'rheome:snapshot:option');
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'Paths', struct('Faces', [1 2], 'Weights', [1 0 0], 'How', 'h', 'Label', 'p'))), 'rheome:snapshot:path');
            tc.verifyError(@() run(struct('Layers', {{{'a', 'x'}}})), 'rheome:snapshot:option');            % a real field has form f only
            tc.verifyError(@() run(struct('Layers', {{{'v', 'x', 'colour'}}})), 'rheome:snapshot:layers');  % a vector form as arrows
            tc.verifyError(@() run(struct('Layers', {{{'a'}, {'v', 'div'}}}, 'Unit', 'mm')), 'rheome:snapshot:layers');   % two colour layers
            tc.verifyError(@() run(struct('Layers', {{{'v', 'div'}}}, 'Unit', '1/s')), 'rheome:snapshot:clim');  % the viewer's form: CLim
            tc.verifyError(@() run(struct('Layers', {{{'a'}}})), 'rheome:snapshot:unit');                   % a colour layer: a colourbar
            tc.verifyError(@() run(struct('Layers', {{{'v'}}}, 'Vector', 'v')), 'rheome:snapshot:layers');
            % every check passed: it gets as far as looking for Node
            tc.verifyError(@() run(struct('Layers', {{{'v', 'x', 'arrows'}}})), 'rheome:snapshot:node');    % arrows only: no unit
            tc.verifyError(@() run(struct('Layers', {{}}, 'Paths', struct('Faces', [1 2], 'Weights', [1 0 0; 0 0.5 0.5], 'How', 'h', 'Label', 'p'))), 'rheome:snapshot:node');
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'Selection', struct('Centre', 1, 'RadiusMM', 80, 'Label', 'disk'), ...
                'Paths', struct('From', 1, 'To', 4, 'Label', 'p'), 'Inset', struct('Label', 'the disk'))), 'rheome:snapshot:node');
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'View', 'R medial', 'Vector', 'v')), 'rheome:snapshot:node');
        end
    end
    methods (Test, TestTags = {'Viewer'})
        function rendersACortex(tc)
            [vdir, S] = needs(tc);
            V = S.Vertices;  c = mean(V);
            x = (V(:, 1) - c(1)) * 1e3;
            r = struct('x', [x x x], 'm', abs([x x x]));
            out = rheome.snapshot(S, r, Out=tempname, Viewer=vdir, Fs=100, Shots={ ...
                struct('Field', 'x', 'Label', 'anterior position', 'Unit', 'mm', 'Name', 'four'), ...
                struct('Field', 'm', 'Label', 'distance', 'Unit', 'mm', 'Reference', 'centre', 'Frames', 1:3, 'View', 'R lateral', 'Name', 'strip')});
            tc.addTeardown(@() rmdir(fileparts(out(1).svg), 's'));
            a = out(1);  b = out(2);
            tc.verifyEqual(a.views, {'L lateral', 'L medial', 'R medial', 'R lateral'});
            tc.verifyEqual(a.quantity, 'signed');  tc.verifyEqual(a.colormap, 'PuOr_r');
            tc.verifyEqual(a.clim, prctile(abs(x), 95) * [-1 1], 'AbsTol', 1e-9);   % the token: +-p95 of |x|, symmetric
            tc.verifyEqual(b.quantity, 'magnitude');  tc.verifyEqual(b.colormap, 'cividis');  tc.verifyEqual(b.clim(1), 0);
            tc.verifyEqual(b.colorbar, 'distance (mm re centre)');
            tc.verifySubstring(b.legend, [char(916) 't = 10 ms']);
            tc.verifyGreaterThanOrEqual([a.dpi b.dpi], 500);
            for o = [a b]
                tc.verifyTrue(isfile(o.svg) && isfile(o.pdf) && isfile(o.png));
                tc.verifyEqual(numel(o.renders), max(numel(o.views), numel(o.times)));
                tc.verifyEqual(o.scalebar.length_mm, 20, 'AbsTol', 1e-9);   % the curve's length along the folded surface
            end
            tc.verifySubstring(fileread(a.svg), 'width="180mm"');          % the theme's double column
            tc.verifySubstring(fileread(b.svg), 'width="180mm"');          % a strip: a row, double
            pdf = fileread(a.pdf);
            tc.verifyTrue(startsWith(pdf, '%PDF'));
            tc.verifyEqual(size(imread(a.png), 2), round(180 / 25.4 * 600), 'AbsTol', 2);   % 600 dpi at 180 mm (the Pergammon theme)
        end
        function insetIsTheRegionAtHigherResolution(tc)
        % the inset is a re-render (more pixels per mm of cortex than the view), and its outline on the view runs where
        % the viewer paints the region
            [vdir, S] = needs(tc);
            V = S.Vertices;  lat = find(V(:, 2) == max(V(:, 2)), 1);
            reg = vecnorm(V - V(lat, :), 2, 2) < 0.012;                         % a 12 mm patch on the lateral surface
            x = (V(:, 1) - mean(V(:, 1))) * 1e3;  fine = x + 3 * sin(V(:, 3) * 900);
            out = rheome.snapshot(S, struct('x', x, 'fine', fine), Out=tempname, Viewer=vdir, Format=[], ScaleBar=0, ...
                Shots=struct('Field', 'x', 'Unit', 'mm', 'Inset', ...
                struct('Vertices', reg, 'Field', 'fine', 'Label', 'patch, finer field', 'Zoom', 4, 'View', 'L lateral')));
            tc.addTeardown(@() rmdir(fileparts(out.svg), 's'));
            tc.verifyGreaterThan(out.inset.magnification, 1.5);                % shown larger than on its view
            tc.verifyGreaterThanOrEqual(out.inset.dpi, 600 - 1e-6);              % never enlarged past its pixels (I1 is checked inside)
            tc.verifyError(@() rheome.snapshot(S, struct('x', x), Out=tempname, Viewer=vdir, Format=[], ScaleBar=0, Shots=struct('Field', 'x', ...
                'Unit', 'mm', 'View', 'L lateral', 'Inset', struct('Vertices', reg, 'Label', 'p'))), 'rheome:snapshot:inset');   % one view at 600 dpi: < 2x
            tc.verifySubstring(out.legend, 'inset: patch, finer field');
            tc.verifySubstring(fileread(out.svg), 'patch, finer field');
            paint = rheome.snapshot(S, struct('m', double(reg)), Out=tempname, Viewer=vdir, Format=[], ScaleBar=0, ...
                Shots=struct('Field', 'm', 'Unit', 'a.u.', 'Quantity', 'magnitude', 'CLim', [0 1], 'Colormap', 'viridis', ...
                'Settings', struct('isolines', 0)));
            tc.addTeardown(@() rmdir(fileparts(paint.svg), 's'));
            [im, ~, al] = imread(paint.renders{1});                             % L lateral, as the inset's view
            w0 = size(imread(out.renders{1}), 2);                               % the renders differ in pixel ratio only
            px = min(max(round((out.inset.px - 0.5) * size(al, 2) / w0 + 0.5), 1), flip(size(al)));
            g = double(im);                                                       % viridis >= ~0.47: green over blue
            hot = g(:, :, 2) >= g(:, :, 3) & al > 0;                              % (lighting scales both alike)
            [r, c] = find(hot);
            near = r < 0.9 * size(hot, 1) | c > 0.2 * size(hot, 2);                  % not the viewer's own colourbar (lower left)
            r = r(near);  c = c(near);
            tc.assertNotEmpty(r, 'the painted region is not visible');
            dOut = sqrt(min((px(:, 1) - c').^2 + (px(:, 2) - r').^2, [], 2));          % drawn -> nearest painted
            dIn = sqrt(min((c - px(:, 1)').^2 + (r - px(:, 2)').^2, [], 2));          % painted -> nearest drawn
            % the drawn region is the region's vertices' nearest-pixel cells; the painted one the interpolated indicator's
            % 1/2 iso-line: both lie half way between the region's vertices and the others, to within a vertex spacing
            tol = sqrt(numel(r) / nnz(reg));
            tc.verifyLessThan(norm(mean(px) - mean([c r])), tol, 'the region drawn is off-centre');
            tc.verifyGreaterThan(mean(dOut <= tol), 0.95, 'the region drawn spills off where the viewer paints it');
            tc.verifyGreaterThan(mean(dIn <= tol), 0.95, 'the region drawn misses part of where the viewer paints it');
        end
        function selectionPathAndCameraAreWhereTheViewerDraws(tc)
        % a geodesic disk selected (filled), a tile as an extent (dashed), a geodesic path, a square inset; the camera
        % returned places a vertex where the render's own camera does, and the path is a geodesic on the folded surface
            [vdir, S] = needs(tc);
            V = S.Vertices;  lat = find(V(:, 2) == max(V(:, 2)), 1);
            near = find(vecnorm(V - V(lat, :), 2, 2) < 0.03 & V(:, 2) > max(V(:, 2)) - 0.015);
            [~, j] = max(vecnorm(V(near, :) - V(lat, :), 2, 2));  far = near(j);
            x = (V(:, 1) - mean(V(:, 1))) * 1e3;
            tile = vecnorm(V - V(far, :), 2, 2) < 0.008;
            out = rheome.snapshot(S, struct('x', x), Out=tempname, Viewer=vdir, Format=[], ScaleBar=0, Shots=struct( ...
                'Field', 'x', 'Unit', 'mm', 'View', 'standard', ...
                'Selection', [struct('Centre', lat, 'RadiusMM', 12, 'Label', 'disk, 12 mm', 'Vertices', [], 'Dashed', false); ...
                              struct('Centre', [], 'RadiusMM', [], 'Label', 'support', 'Vertices', tile, 'Dashed', true)], ...
                'Paths', struct('From', lat, 'To', far, 'Label', 'path of the event'), ...
                'Inset', struct('Label', 'the disk', 'Shape', 'square', 'View', 'L lateral')));
            tc.addTeardown(@() rmdir(fileparts(out.svg), 's'));
            svg = fileread(out.svg);
            tc.verifySubstring(svg, 'stroke-dasharray');                                % the extent
            tc.verifyEqual({out.selections.kind}, {'disk', 'tile'});
            tc.verifyEqual(out.inset.shape, 'square');
            p = out.paths;
            tc.verifyGreaterThan(p.length_mm, 0.99 * p.geodesic_mm - 0.5);             % a path on the surface
            tc.verifyLessThan(p.length_mm, 1.06 * p.geodesic_mm + 0.5);                 % and a geodesic, not a mesh walk
            tc.verifySubstring(out.legend, 'path of the event');
            c = out.cameras(1);                                                          % L lateral
            tc.verifyEqual(c.view, 'L lateral');
            px = c.s * V(lat, :) * c.R' + c.off;                                        % the vertex on the render
            [~, ~, al] = imread(out.renders{1});
            tc.verifyGreaterThan(al(round(px(2)), round(px(1))), 0, 'the camera puts the vertex on the cortex');
            tc.verifyEqual(out.settings{1}.paths, {'shot1_path1'}, 'the viewer draws the shot''s path');
        end
        function arrowsOnlyWithAPathTheViewerDraws(tc)
        % an arrows-only shot (no colour layer: no colourbar) with a path given as faces + barycentric weights: the
        % viewer accepts the layers and draws the path in the render, where the render's camera puts its points
            [vdir, S] = needs(tc);
            V = S.Vertices;  nV = size(V, 1);  lat = find(V(:, 2) == max(V(:, 2)), 1);
            near = find(vecnorm(V - V(lat, :), 2, 2) < 0.03 & V(:, 2) > max(V(:, 2)) - 0.015);
            [~, j] = max(vecnorm(V(near, :) - V(lat, :), 2, 2));  far = near(j);
            [X, info] = rheome.geom.geodesicpath(S, lat, far);
            fc = info.faces([1:end end]);
            A = V(S.Faces(fc, 1), :);  e1 = V(S.Faces(fc, 2), :) - A;  e2 = V(S.Faces(fc, 3), :) - A;
            w = zeros(numel(fc), 3);
            for i = 1:numel(fc), uv = [e1(i, :); e2(i, :)]' \ (X(i, :) - A(i, :))';  w(i, :) = [1 - sum(uv) uv']; end
            w = max(w, 0);  w = w ./ sum(w, 2);
            fg = rheome.operators.face_gradient(V, S.Faces);  N = fg.W * fg.FaceNormal;  N = N ./ vecnorm(N, 2, 2);
            v = repmat([0 0 1], nV, 1);  v = v - sum(v .* N, 2) .* N;   % tangent, upward
            r = struct('v', reshape(v', [], 1));
            out = rheome.snapshot(S, r, Out=tempname, Viewer=vdir, Format="pdf", ScaleBar=0, Shots=struct('View', 'L lateral', ...
                'Layers', {{struct('field', 'v', 'form', 'x', 'as', 'arrows', 'density', 0.05)}}, 'Paths', struct('Faces', fc, 'Weights', w, 'How', 'test: geodesicpath', ...
                'Label', 'trajectory', 'Color', '#ff00ff', 'Width', 2)));
            tc.addTeardown(@() rmdir(fileparts(out.svg), 's'));
            tc.verifyEmpty(out.colormap);  tc.verifyEmpty(out.colorbar);
            tc.verifyTrue(isfile(out.pdf));
            L = out.settings{1}.layers;  if iscell(L), L = L{1}; end
            tc.verifyEqual({L.field, L.as}, {'v', 'arrows'});
            [im, ~, al] = imread(out.renders{1});
            g = double(im);  mag = g(:, :, 1) > 180 & g(:, :, 3) > 180 & g(:, :, 2) < 120 & al > 0;
            [rr, cc] = find(mag);
            tc.assertNotEmpty(rr, 'the path is not drawn in the render');
            c = out.cameras(1);  P = out.paths.points;
            P = interp1(1:size(P, 1), P, linspace(1, size(P, 1), 50 * size(P, 1)));   % along its segments
            px = c.s * P * c.R' + c.off;
            d = min(hypot(cc - px(:, 1)', rr - px(:, 2)'), [], 2);       % each drawn pixel to the path
            head = hypot(cc - px(end, 1), rr - px(end, 2)) < 80;          % the arrowhead: a cone 5 x 4 widths at the end
            tc.verifyGreaterThan(nnz(head), 100, 'no arrowhead at the path''s end');
            % the tube, 2 pt (~11 px at this pixel ratio): within 25 px of the path's line
            tc.verifyGreaterThan(mean(d(~head) < 25), 0.95, 'the drawn path is off where its points project');
            tc.verifyEqual([min(cc) max(cc) min(rr) max(rr)], [min(px(:, 1)) max(px(:, 1)) min(px(:, 2)) max(px(:, 2))], 'AbsTol', 30);
        end
        function barLiesWhereTheViewerDrawsIt(tc)
        % paint the bar's vertices (and their neighbours) in the viewer; the curve drawn from the camera read off the
        % render must run through the painted pixels
            [vdir, S] = needs(tc);
            nV = size(S.Vertices, 1);
            first = rheome.snapshot(S, struct('z', S.Vertices(:, 3) * 1e3), Out=tempname, Viewer=vdir, Format=[], ...
                Shots=struct('Field', 'z', 'Unit', 'mm', 'View', 'L lateral'));
            tc.addTeardown(@() rmdir(fileparts(first.svg), 's'));
            tc.verifySubstring(fileread(first.svg), 'width="88mm"');       % one view: single column
            tc.verifyGreaterThanOrEqual(first.dpi, 500);                    % the render scale raised to print it at 500 dpi
            bar = first.scalebar;
            A = sparse(double(S.Faces), double(S.Faces(:, [2 3 1])), 1, nV, nV);  A = A + A' + speye(nV);
            mark = full(A * sparse(bar.vertices, 1, 1, nV, 1)) > 0;
            second = rheome.snapshot(S, struct('mark', double(mark)), Out=tempname, Viewer=vdir, Format=[], ScaleBar=0, ...
                Shots=struct('Field', 'mark', 'Unit', 'a.u.', 'Quantity', 'magnitude', 'CLim', [0 1], 'Colormap', 'viridis', ...
                'View', 'L lateral', 'Settings', struct('isolines', 0)));
            tc.addTeardown(@() rmdir(fileparts(second.svg), 's'));
            [im, ~, al] = imread(second.renders{1});
            hot = double(im(:, :, 1)) > 150 & double(im(:, :, 2)) > 150 & double(im(:, :, 3)) < 120 & al > 0;   % viridis' yellow
            [r, c] = find(hot);
            tc.assertNotEmpty(r, 'the painted vertices are not visible');
            d = min(hypot(bar.px(:, 1) - c', bar.px(:, 2) - r'), [], 2);
            tc.verifyLessThan(max(d), 8, 'the drawn bar is off the painted vertices (px)');
        end
    end
end

function [vdir, S] = needs(tc)
    vdir = getenv('RHEOME_VIEWER');
    if isempty(vdir), vdir = fullfile(char(java.lang.System.getProperty('user.home')), 'workspace', 'research', 'code', 'nxr-cortical-flow', 'editor', 'dist-viewer'); end
    tc.assumeTrue(isfile(fullfile(vdir, 'viewer.html')), 'needs the viewer (RHEOME_VIEWER)');
    tc.assumeTrue(isfile('/Applications/Google Chrome.app/Contents/MacOS/Google Chrome') || ~isempty(getenv('RHEOME_CHROME')), 'needs Chrome');
    tc.assumeTrue(rheome.load.has('omega_sub0002'), 'needs the omega_sub0002 surface in +data');
    S = rheome.load.surface('omega_sub0002');
end
