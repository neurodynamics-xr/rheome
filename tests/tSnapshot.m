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
            % every check passed: it gets as far as looking for Node
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
