classdef tSnapshot < matlab.unittest.TestCase
% TSNAPSHOT  rheome.snapshot: the shot list is checked against the figure tokens before anything is drawn (a stand-in
% viewer folder, no browser); with the real viewer, Chrome and Node (RHEOME_VIEWER or the nxr-cortical-flow clone,
% skipped otherwise), a real cortex is rendered and composed by the viewer's print figure, the files are print-sized,
% and the pixels the viewer reports (the scale bar, a path, a selection's outline) lie where it paints their vertices.
%
% Author: Diellor Basha, 2026
    properties
        S = struct('Vertices', [0 0 0; 1 0 0; 0 1 0; 0 0 1] * 0.05, 'Faces', [1 2 3; 1 3 4; 1 4 2; 2 4 3]);
        Fake
    end
    methods (TestMethodSetup)
        function fakeViewer(tc)
        % a viewer folder with the files read before rendering: viewer.html, the driver and the default colormap library
            d = tempname;  mkdir(fullfile(d, 'colormaps', 'default'));  mkdir(fullfile(d, 'driver'));
            fclose(fopen(fullfile(d, 'viewer.html'), 'w'));
            fclose(fopen(fullfile(d, 'driver', 'snapshot-driver.mjs'), 'w'));
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
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'Frames', [1 2])), 'rheome:snapshot:node');      % a strip of several views: allowed (phase 3)
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'Frames', [1 2 4], 'View', 'L lateral')), 'rheome:snapshot:shot');   % unequal
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'Frames', 5)), 'rheome:snapshot:shot');
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'CLim', [-1 2])), 'rheome:snapshot:clim');      % signed: symmetric
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'Quantity', 'coupling')), 'rheome:snapshot:colormap');   % no BrBG_r
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'Vector', 'a')), 'rheome:snapshot:option');
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'Name', 'a/b')), 'rheome:snapshot:shot');
            tc.verifyError(@() rheome.snapshot(tc.S, r, Viewer=tc.Fake), 'rheome:snapshot:input');                 % Out
            tc.verifyError(@() rheome.snapshot(tc.S, r, Out=tempname, Viewer=tempname), 'rheome:snapshot:noviewer');
            tc.verifyError(@() rheome.snapshot(tc.S, r, Out=tempname, Viewer=tc.Fake, Tokens=[tempname '.json']), 'rheome:viz:tokens');   % (CHANGELOG: was rheome:snapshot:tokens)
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
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'Paths', struct('Points', tc.S.Vertices(1:2, :), 'How', 'by hand', 'Label', 'p'))), ...
                'rheome:snapshot:path');                                                     % only From/To geodesics are drawn now
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'Inset', struct('Label', 'no region'))), 'rheome:snapshot:inset');   % no Vertices, no Selection
            nodrv = tempname;  copyfile(tc.Fake, nodrv);  rmdir(fullfile(nodrv, 'driver'), 's');  tc.addTeardown(@() rmdir(nodrv, 's'));
            tc.verifyError(@() rheome.snapshot(tc.S, r, Out=tempname, Viewer=nodrv), 'rheome:snapshot:driver');   % a viewer without the driver
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'Inset', setfield(in, 'Shape', 'star'))), 'rheome:snapshot:option');
            % every check passed: it gets as far as looking for Node
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'Selection', struct('Centre', 1, 'RadiusMM', 80, 'Label', 'disk'), ...
                'Paths', struct('From', 1, 'To', 4, 'Label', 'p'), 'Inset', struct('Label', 'the disk'))), 'rheome:snapshot:node');
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'View', 'R medial', 'Vector', 'v')), 'rheome:snapshot:node');
        end
        function refusesBadNewOptions(tc)
        % the phase-3 fields (no field, Cortex, Labelling, Outlines, raw pictures, the bar's kind, 'auto' and tile), each refused by name
            r = struct('a', [-1; 0; 1; 2] * ones(1, 4), 'v', zeros(12, 4));
            L = struct('name', {'t', 'h'}, 'codes', {[0; 0; 1; 1], [0; 0; 0; 1]});
            run = @(shot, varargin) rheome.snapshot(tc.S, r, 'Shots', shot, 'Out', tempname, 'Viewer', tc.Fake, 'Node', '/nonexistent', ...
                'Chrome', fullfile(tc.Fake, 'viewer.html'), 'Labellings', L, varargin{:});
            tc.verifyError(@() run(struct('Field', '', 'Unit', 'mm')), 'rheome:snapshot:shot');                   % no field: no unit
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'Labelling', 't')), 'rheome:snapshot:shot');  % a field or a labelling
            tc.verifyError(@() run(struct('Labelling', 'zz')), 'rheome:snapshot:option');
            tc.verifyError(@() rheome.snapshot(tc.S, r, Shots=struct('Labelling', 't'), Out=tempname, Viewer=tc.Fake), 'rheome:snapshot:option');   % no Labellings
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'Cortex', 2)), 'rheome:snapshot:shot');
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'View', {{'L lateral', 'l lateral'}})), 'rheome:snapshot:shot');   % each view once
            tc.verifyError(@() run(struct('Field', '', 'Outlines', struct('Coarser', 1))), 'rheome:snapshot:outlines');         % its labelling
            tc.verifyError(@() run(struct('Field', '', 'Outlines', struct('Labelling', 't', 'Wrong', 1))), 'rheome:snapshot:outlines');
            tc.verifyError(@() run(struct('Field', '', 'Outlines', struct('Labelling', 't', 'Width', 50))), 'rheome:snapshot:outlines');
            tc.verifyError(@() run(struct('Field', '', 'Outlines', struct('Labelling', 't', 'Coarser', 0.5))), 'rheome:snapshot:outlines');
            tc.verifyError(@() run(struct('Field', '', 'Outlines', struct('Labelling', 't', 'Select', 1.5))), 'rheome:snapshot:outlines');
            tc.verifyError(@() run(struct('Field', '', 'Outlines', struct('Labelling', 't', 'SelectStyle', 'dots'))), 'rheome:snapshot:option');
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'StagePx', [10 10])), 'rheome:snapshot:raw');
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'Fit', true)), 'rheome:snapshot:raw');         % fit: raw pictures only
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'Dpi', 300)), 'rheome:snapshot:raw');
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'StagePx', [64 64], 'Fit', struct('Margin', 0.7))), 'rheome:snapshot:raw');
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'StagePx', [64 64], 'Dpi', 30)), 'rheome:snapshot:raw');
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'StagePx', [64 64], 'Inset', struct('Vertices', [1 2 3], 'Label', 'x'))), 'rheome:snapshot:raw');
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm'), 'ScaleBar', 'big'), 'rheome:snapshot:input');
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm'), 'ScaleBarKind', 'curved'), 'rheome:snapshot:option');
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm'), 'ScaleBarTile', 0), 'rheome:snapshot:input');
            tc.verifyError(@() run(struct('Field', 'a', 'Unit', 'mm', 'View', 'L lateral'), 'ScaleBarTile', 2), 'rheome:snapshot:shot');   % 1 render
            % every check passed: as far as Node
            tc.verifyError(@() run(struct('Labelling', 't', 'Outlines', struct('Labelling', 't', 'Coarser', 'h', 'Select', [0 1], 'SelectStyle', 'extent'), ...
                'Cortex', 0.5, 'StagePx', [64 48], 'Fit', struct('Margin', 0.1, 'Shared', true), 'Dpi', 300), 'ScaleBar', 'auto', 'ScaleBarKind', 'plane', ...
                'ScaleBarTile', 4), 'rheome:snapshot:node');
            tc.verifyError(@() run(struct('Field', '', 'Vector', 'v', 'Cortex', false, 'View', {{'L lateral', 'left'}}, 'Frames', 1:3)), 'rheome:snapshot:node');
        end
    end
    methods (Test, TestTags = {'Viewer'})
        function rendersACortex(tc)
            [vdir, S] = needs(tc);
            V = S.Vertices;  c = mean(V);
            x = (V(:, 1) - c(1)) * 1e3;
            r = struct('x', [x x x], 'm', abs([x x x]));
            Out = tempname;
            out = rheome.snapshot(S, r, Out=Out, Viewer=vdir, Fs=100, Shots={ ...
                struct('Field', 'x', 'Label', 'anterior position', 'Unit', 'mm', 'Name', 'four'), ...
                struct('Field', 'm', 'Label', 'distance', 'Unit', 'mm', 'Reference', 'centre', 'Frames', 1:3, 'View', 'R lateral', 'Name', 'strip')});
            tc.addTeardown(@() rmdir(Out, 's'));
            a = out(1);  b = out(2);
            tc.verifyEqual(a.views, {'L lateral', 'L medial', 'R medial', 'R lateral'});
            tc.verifyEqual(a.quantity, 'signed');  tc.verifyEqual(a.colormap, 'PuOr_r');
            tc.verifyEqual(a.clim, prctile(abs(double(single(x))), 95) * [-1 1], 'RelTol', 1e-6);   % the token: +-p95 of |x|, symmetric (the bundle's float32)
            tc.verifyEqual(b.quantity, 'magnitude');  tc.verifyEqual(b.colormap, 'cividis');  tc.verifyEqual(b.clim(1), 0);
            tc.verifyEqual(b.colorbar, 'distance (mm re centre)');
            tc.verifyEqual(b.times, [0 0.01 0.02], 'AbsTol', 1e-12);
            tc.verifySubstring(b.legend, [char(916) 't = 10 ms']);
            tc.verifyGreaterThanOrEqual([a.dpi b.dpi], 500);
            for o = [a b]
                tc.verifyTrue(isfile(o.svg) && isfile(o.pdf) && isfile(o.png));
                tc.verifyTrue(all(cellfun(@isfile, o.renders)));
                tc.verifyEqual(numel(o.renders), max(numel(o.views), numel(o.times)));
                tc.verifyEqual(o.scalebar.length_mm, 20, 'AbsTol', 1e-6);   % the curve's length along the folded surface
                tc.verifyTrue(isfile(fullfile(Out, [o.name '.marks.json'])) && isfile(fullfile(Out, [o.name '.report.json'])));
            end
            tc.verifyEqual(b.renders, fullfile(Out, 'strip', {'r-lateral_f0.png', 'r-lateral_f1.png', 'r-lateral_f2.png'}));   % frames 0-based
            tc.verifyTrue(isfile(fullfile(Out, 'job.json')));
            tc.verifySubstring(fileread(a.svg), 'width="180.000mm"');      % the theme's double column
            tc.verifySubstring(fileread(b.svg), 'width="180.000mm"');      % a strip: a row, double
            pdf = fileread(a.pdf);
            tc.verifyTrue(startsWith(pdf, '%PDF'));
            tc.verifyEqual(size(imread(a.png), 2), round(180 / 25.4 * 600), 'AbsTol', 2);   % 600 dpi at 180 mm (the Pergammon theme)
        end
        function insetIsTheRegionAtHigherResolution(tc)
        % the inset is a re-render (more pixels per mm of cortex than the view), and the region's outline the viewer
        % reports (a dashed extent of the same region, its pixels in the marks) runs where the viewer paints the region
            [vdir, S] = needs(tc);
            V = S.Vertices;  lat = find(V(:, 2) == max(V(:, 2)), 1);
            reg = vecnorm(V - V(lat, :), 2, 2) < 0.012;                         % a 12 mm patch on the lateral surface
            x = (V(:, 1) - mean(V(:, 1))) * 1e3;  fine = x + 3 * sin(V(:, 3) * 900);
            Out = tempname;
            out = rheome.snapshot(S, struct('x', x, 'fine', fine), Out=Out, Viewer=vdir, Format=[], ScaleBar=0, ...
                Shots=struct('Field', 'x', 'Unit', 'mm', 'View', 'L lateral', 'Name', 'inset', ...
                'Selection', struct('Vertices', reg, 'Label', 'patch', 'Dashed', true), ...
                'Inset', struct('Field', 'fine', 'Label', 'patch, finer field', 'Zoom', 4)));
            tc.addTeardown(@() rmdir(Out, 's'));
            tc.verifyGreaterThanOrEqual(out.inset.magnification, 2);            % rule I1: >= 2x the view's pixels per mm
            tc.verifyGreaterThanOrEqual(out.inset.dpi, 600 - 1e-6);              % never enlarged past its pixels
            tc.verifyEqual(out.inset.vertices, find(reg));                       % the region: the selection's (1-based)
            tc.verifyEqual(out.selections.vertices, find(reg));
            tc.verifyTrue(isfile(fullfile(Out, 'inset', 'l-lateral_inset.png')));
            tc.verifySubstring(out.legend, 'inset: patch, finer field');
            tc.verifySubstring(fileread(out.svg), 'patch, finer field');
            mk = jsondecode(fileread(fullfile(Out, 'inset.marks.json')));
            ovs = mk.overlays;  if isstruct(ovs), ovs = num2cell(ovs); end
            k = find(cellfun(@(o) strcmp(o.overlay, 'selection: patch'), ovs), 1);
            tc.assertNotEmpty(k, 'the extent is not in the marks');
            px = pixelsOf(ovs{k}.pixels);
            paint = rheome.snapshot(S, struct('m', double(reg)), Out=tempname, Viewer=vdir, Format=[], ScaleBar=0, ...
                Shots=struct('Field', 'm', 'Unit', 'a.u.', 'Quantity', 'magnitude', 'CLim', [0 1], 'Colormap', 'viridis', ...
                'View', 'L lateral', 'Settings', struct('isolines', 0)));
            tc.addTeardown(@() rmdir(fileparts(paint.svg), 's'));
            [im, ~, al] = imread(paint.renders{1});
            px = px * size(al, 2) / size(imread(out.renders{1}), 2);             % the renders differ in pixel ratio only
            g = double(im);                                                       % viridis >= ~0.47: green over blue
            hot = g(:, :, 2) >= g(:, :, 3) & al > 0;                              % (lighting scales both alike)
            [r, c] = find(hot);
            tc.assertNotEmpty(r, 'the painted region is not visible');
            c = c - 0.5;  r = r - 0.5;                                            % pixel centres, the report's corner origin
            dOut = sqrt(min((px(:, 1) - c').^2 + (px(:, 2) - r').^2, [], 2));          % drawn -> nearest painted
            % the drawn outline is the region's boundary on the mesh; the painted region the interpolated indicator's
            % 1/2 iso-line: both lie half way between the region's vertices and the others, to within a vertex spacing
            tol = sqrt(numel(r) / nnz(reg));
            tc.verifyLessThan(norm(mean(px) - mean([c r])), 2 * tol, 'the outline drawn is off-centre');
            tc.verifyGreaterThan(mean(dOut <= tol), 0.95, 'the outline drawn runs off where the viewer paints the region');
        end
        function selectionPathAndPixelsAreWhereTheViewerDraws(tc)
        % a geodesic disk selected (filled), a tile as an extent (dashed), a geodesic path, a square inset; the path's
        % first pixel is where the viewer paints its first vertex, and the path is a geodesic on the folded surface
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
            tc.verifyEqual(out.selections(2).vertices, find(tile));
            tc.verifyTrue(ismember(lat, out.selections(1).vertices));
            tc.verifyEqual(out.inset.shape, 'square');
            p = out.paths;
            G = rheome.geom.edgegraph(S, Diagonals=false);                              % the exact geodesic (the viewer's) is
            tc.verifyLessThanOrEqual(p.length_mm, 1e3 * distances(G, lat, far) + 1e-6); % no longer than a mesh walk,
            tc.verifyGreaterThanOrEqual(p.length_mm, 1e3 * norm(V(lat, :) - V(far, :)));  % no shorter than the chord,
            tc.verifyEqual(p.length_mm, p.geodesic_mm, 'RelTol', 0.1);                  % and near the heat distance between its ends
            tc.verifyEqual(p.points(1, :), V(lat, :), 'AbsTol', 1e-6);                  % from its first vertex (m)
            tc.verifySubstring(out.legend, 'path of the event');
            tc.verifyEqual(numel(p.px), 4);  tc.verifyEqual(numel(p.visible), 4);
            tc.verifyTrue(p.visible{1}(1), 'the path starts seen on the L lateral view');
            tc.verifyFalse(any(p.visible{3}) || any(p.visible{4}), 'a left-hemisphere path is not drawn on the right''s views');
            % the path's first pixel where the viewer paints its first vertex (and its one ring) on the L lateral view
            nV = size(V, 1);
            A = sparse(double(S.Faces), double(S.Faces(:, [2 3 1])), 1, nV, nV);  A = A + A' + speye(nV);
            mark = full(A * sparse(lat, 1, 1, nV, 1)) > 0;
            paint = rheome.snapshot(S, struct('mark', double(mark)), Out=tempname, Viewer=vdir, Format=[], ScaleBar=0, ...
                Shots=struct('Field', 'mark', 'Unit', 'a.u.', 'Quantity', 'magnitude', 'CLim', [0 1], 'Colormap', 'viridis', ...
                'View', 'L lateral', 'Settings', struct('isolines', 0)));
            tc.addTeardown(@() rmdir(fileparts(paint.svg), 's'));
            [im, ~, al] = imread(paint.renders{1});
            q = p.px{1}(1, :) * size(al, 2) / size(imread(out.renders{1}), 2);       % the renders differ in pixel ratio only
            hot = double(im(:, :, 1)) > 150 & double(im(:, :, 2)) > 150 & double(im(:, :, 3)) < 120 & al > 0;   % viridis' yellow
            [r, c] = find(hot);
            tc.assertNotEmpty(r, 'the painted vertex is not visible');
            tc.verifyLessThan(min(hypot(q(1) - (c - 0.5), q(2) - (r - 0.5))), 4, 'the path''s first pixel is off its painted vertex (px)');
        end
        function barLiesWhereTheViewerDrawsIt(tc)
        % paint the bar's vertices (and their neighbours) in the viewer; the curve the viewer reports must run through
        % the painted pixels
            [vdir, S] = needs(tc);
            nV = size(S.Vertices, 1);
            first = rheome.snapshot(S, struct('z', S.Vertices(:, 3) * 1e3), Out=tempname, Viewer=vdir, Format=[], ...
                Shots=struct('Field', 'z', 'Unit', 'mm', 'View', 'L lateral'));
            tc.addTeardown(@() rmdir(fileparts(first.svg), 's'));
            tc.verifySubstring(fileread(first.svg), 'width="88.000mm"');   % one view: single column
            tc.verifyGreaterThanOrEqual(first.dpi, 500);                    % the render scale raised to print it at 500 dpi
            bar = first.scalebar;
            tc.verifyEqual(bar.view, 'L lateral');
            A = sparse(double(S.Faces), double(S.Faces(:, [2 3 1])), 1, nV, nV);  A = A + A' + speye(nV);
            mark = full(A * sparse(bar.vertices, 1, 1, nV, 1)) > 0;
            second = rheome.snapshot(S, struct('mark', double(mark)), Out=tempname, Viewer=vdir, Format=[], ScaleBar=0, ...
                Shots=struct('Field', 'mark', 'Unit', 'a.u.', 'Quantity', 'magnitude', 'CLim', [0 1], 'Colormap', 'viridis', ...
                'View', 'L lateral', 'Settings', struct('isolines', 0)));
            tc.addTeardown(@() rmdir(fileparts(second.svg), 's'));
            [im, ~, al] = imread(second.renders{1});
            tc.verifyEqual(size(al), size(imread(first.renders{1}), [1 2]));   % the same view at the same pixel ratio
            hot = double(im(:, :, 1)) > 150 & double(im(:, :, 2)) > 150 & double(im(:, :, 3)) < 120 & al > 0;   % viridis' yellow
            [r, c] = find(hot);
            tc.assertNotEmpty(r, 'the painted vertices are not visible');
            d = min(hypot(bar.px(:, 1) - (c' - 0.5), bar.px(:, 2) - (r' - 0.5)), [], 2);
            tc.verifyLessThan(max(d), 8, 'the drawn bar is off the painted vertices (px)');
        end
    end
end

function px = pixelsOf(p)
% a marks overlay's pixels ([x y seen] per point, as JSON lists of runs) as [x y] rows
    if iscell(p), p = cellfun(@(q) reshape(q, [], 3), p(:), 'UniformOutput', false); p = vertcat(p{:});
    else, p = reshape(p, [], 3);
    end
    px = p(:, 1:2);
end

function [vdir, S] = needs(tc)
    vdir = getenv('RHEOME_VIEWER');
    if isempty(vdir), vdir = fullfile(char(java.lang.System.getProperty('user.home')), 'workspace', 'research', 'code', 'nxr-cortical-flow', 'editor', 'dist-viewer'); end
    tc.assumeTrue(isfile(fullfile(vdir, 'viewer.html')), 'needs the viewer (RHEOME_VIEWER)');
    tc.assumeTrue(isfile('/Applications/Google Chrome.app/Contents/MacOS/Google Chrome') || ~isempty(getenv('RHEOME_CHROME')), 'needs Chrome');
    tc.assumeTrue(rheome.load.has('omega_sub0002'), 'needs the omega_sub0002 surface in +data');
    S = rheome.load.surface('omega_sub0002');
end
