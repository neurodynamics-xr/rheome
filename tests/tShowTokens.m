classdef tShowTokens < matlab.unittest.TestCase
% Every rheome.show figure takes its sizes and colours from the figure tokens
% (library/styles/tokens/figure-tokens.json) and its export passes the lab's figure check
% (library/styles/schematics/figure_check.py, on the PDF: column width, page height, fonts embedded and
% in the house face, no text below 6 pt). Request d12bfeff, 2026-10-14.
% Revision 1 (request e2ea030f): every figure is drawn in each theme (Pergammon, NeuroImage for MS1) and
% checked against it; colourbars carry a token quantity (one colormap per quantity) and sit in a key band;
% the colormaps' descriptive names are on the figure for the legend; no Unicode maths glyphs in +show.
%
% Author: Diellor Basha, 2026

    properties (ClassSetupParameter)
        Theme = {'pergammon', 'neuroimage'};
    end

    properties
        S; out; styles; check; theme
    end

    methods (TestClassSetup)
        function build(t, Theme)
            t.theme = Theme;
            [x, y, z] = sphere(30);  fv = surf2patch(0.07*x, 0.07*y, 0.07*z, 'triangles');
            t.S = struct('Vertices', fv.vertices, 'Faces', fv.faces, 'nV', size(fv.vertices, 1));
            t.out = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'results', 'figure-check', Theme);
            if ~isfolder(t.out), mkdir(t.out); end
            t.styles = getenv('NXR_STYLES');
            if isempty(t.styles), t.styles = fullfile(getenv('HOME'), 'workspace', 'library', 'styles'); end
            t.check = fullfile(t.styles, 'schematics', 'figure_check.py');
            rng(7);
        end
    end

    methods (TestMethodTeardown)
        function shut(~), close all force; end
    end

    methods (Test)
        function eventsensors(t)
            C = 12;  nS = 300;  X = randn(C, nS);  m = false(C, nS);  m(1:4, 100:160) = true;
            E = struct('footprint', rand(C, 1), 'channels', {{(1:4)'}}, 'mask', m);
            h = rheome.show.eventsensors(X, (0:nS-1)/300, E, randn(C, 2), 'NumChannels', 6, 'Visible', 'off', ...
                'Unit', 'fT', 'Theme', t.theme, 'Export', t.stem('eventsensors'));
            t.passes(h, 'eventsensors');
        end

        function eventsensorsSpan(t)
            C = 12;  nS = 300;  X = randn(C, nS);  m = false(C, nS);  m(1:4, 100:160) = true;
            E = struct('footprint', rand(C, 1), 'channels', {{(1:4)'}}, 'mask', m);
            h = rheome.show.eventsensors(X, (0:nS-1)/300, E, randn(C, 2), 'NumChannels', 6, 'Visible', 'off', ...
                'Unit', 'fT', 'EventMode', 'span', 'EventLabel', 'spindle', 'Theme', t.theme, 'Export', t.stem('eventsensors-span'));
            t.verifyNotEmpty(findall(h, 'Tag', 'nxr-event-span'));
            t.passes(h, 'eventsensors-span');
        end

        function activity(t)
            K = 20;  nT = 80;  Psi = orth(randn(t.S.nV, K));  c = randn(K, nT) + 1i*randn(K, nT);
            h = rheome.show.activity(t.S, c, Psi, (0:nT-1)/100, 'Visible', 'off', 'Theme', t.theme, 'Export', t.stem('activity'));
            t.passes(h, 'activity');
        end

        function critical(t)
            types = {'vortex', 'source', 'sink', 'saddle'};
            cp = struct('pos', t.S.Vertices(1:40:160, :), 'type', {types}, 'strength', [1; 2; 3; 4]);
            o = struct('cp', cp, 'frames', [1 2 3], 'counts', [1 1 1 1; 2 0 1 1; 1 2 0 1], 'kinds', {types}, ...
                       'chi', [2 2 2], 'chiExpected', 2);
            h = rheome.show.critical(t.S, o, 'Visible', 'off', 'Theme', t.theme, 'Export', t.stem('critical'));
            t.passes(h, 'critical');
        end

        function filmstrip(t)
            U = sin(40*t.S.Vertices(:, 1) + (1:5)) .* exp(10*t.S.Vertices(:, 3));
            h = rheome.show.filmstrip(t.S, U, 1:5, 'Visible', 'off', 'RowLabel', 'Impulse response', ...
                'Titles', compose('t = %d ms', (0:4)*20), 'Unit', 'nA·m', 'Theme', t.theme, 'Export', t.stem('filmstrip'));
            t.passes(h, 'filmstrip');
        end

        function vectorstrip(t)
            J = randn(3*t.S.nV, 4);
            h = rheome.show.vectorstrip(t.S, J, 1:4, 'Visible', 'off', 'Zoom', 1, 'Theme', t.theme, 'Export', t.stem('vectorstrip'));
            t.passes(h, 'vectorstrip');
        end

        function frame(t)
            lam = linspace(0, 4e5, 200)';
            fr = rheome.filters.frame('itersine', 5, [min(lam(lam > 0)) max(lam)]);
            h = rheome.show.frame(fr, lam, 'Visible', 'off', 'Theme', t.theme, 'Export', t.stem('frame'));
            t.passes(h, 'frame');
        end

        function jointscalogram(t)
            sc = struct('sigma', [NaN 0.01 0.02 0.04], 'f', 1:40, 'energy', rand(4, 40), 'kbar', 100:10:200, ...
                        'fRidge', linspace(5, 20, 11), 'wRidge', rand(1, 11), 'slope', 0.5, 'slopeR2', 0.9, ...
                        'fSupport', [4 22]);
            h = rheome.show.jointscalogram(sc, 'Visible', 'off', 'Theme', t.theme, 'Export', t.stem('jointscalogram'));
            t.passes(h, 'jointscalogram');
        end

        function jointspectrum(t)
            lam = [linspace(0, 4e5, 30) linspace(0, 4e5, 30)]';
            h = rheome.show.jointspectrum(randn(60, 40) + 1i, 1:40, lam, 'Speed', [0.5 2], 'Visible', 'off', ...
                'Theme', t.theme, 'Export', t.stem('jointspectrum'));
            t.passes(h, 'jointspectrum');
        end

        function resolution(t)
            marks = table(["fine"; "mid"; "coarse"], [20; 100; 300], ["instrument"; "floor"; "pyramid"], ...
                          'VariableNames', {'name', 'mm', 'kind'});
            h = rheome.show.resolution(marks, 'Surface', t.S, 'Visible', 'off', 'Patches', {1:50, 1:20}, ...
                'PatchTitles', {'a', 'b'}, 'Fields', rand(t.S.nV, 2), 'FieldTitles', {'c', 'd'}, ...
                'Bands', {10, 80, 'low'; 80, 400, 'high'}, 'Theme', t.theme, 'Export', t.stem('resolution'));
            t.passes(h, 'resolution');
        end

        function ridges(t)
            o = t.ridgeOut();
            h = rheome.show.ridges(t.S, o, 'Visible', 'off', 'Truth', struct('sigma', 0.02, 'speed', 0.5), ...
                'Theme', t.theme, 'Export', t.stem('ridges'));
            t.passes(h, 'ridges');
        end

        function ridgestrip(t)
            U = rand(t.S.nV, 6);
            h = rheome.show.ridgestrip(t.S, U, t.ridgeOut(), 1:5, 'Visible', 'off', 'Zoom', 1, ...
                'Theme', t.theme, 'Export', t.stem('ridgestrip'));
            t.passes(h, 'ridgestrip');
        end

        function specparam(t)
            f = (1:60)';  Pap = 1 ./ f.^[1 1.5 2];  P = Pap .* (1 + exp(-(f - 10).^2/4));
            ap = struct('Pap', Pap, 'exponent', [1 1.5 2]);
            h = rheome.show.specparam(P, f, ap, 'Visible', 'off', 'Theme', t.theme, 'Export', t.stem('specparam'));
            t.passes(h, 'specparam');
        end

        function timetracks(t)
            tt = linspace(0, 1, 400);
            tr = struct('Label', {'heading', 'scale', 'chirality'}, 'Unit', {'rad', 'mm', 'sign'}, ...
                'Spec', {struct('Start', 0, 'To', 4*pi, 'Profile', 'ramp'), struct('Start', 10, 'To', 30, 'Profile', 'ease'), ...
                         struct('Keys', struct('x', {0, 0.5}, 'v', {1, -1}, 'seg', {'hold', 'hold'}))}, ...
                'Values', {4*pi*tt(1:20:end) + 0.2*randn(1, 20), [], []}, 'Time', {tt(1:20:end), [], []});
            [h, V] = rheome.show.timetracks(tt, tr, 'Visible', 'off', 'Playhead', 0.3, 'Theme', t.theme, 'Export', t.stem('timetracks'));
            t.verifyEqual(V(1).driver, ['ramp 0' char(8594) '12.6']);  t.verifyEqual(V(3).keyT, [0 0.5]);
            t.verifyEqual({V.kind}, {'angle', 'value', 'sign'});
            t.verifyEqual(V(3).plant([1 end]), [1 -1]);
            t.passes(h, 'timetracks');
        end

        function phasor(t)
            F = nxr_figure('single', 70, {'A'}, 'Kinds', struct('A', 'plot'), 'Letters', false, 'Theme', t.theme);
            axis(F.ax.A, 'equal', 'off');  F.fig.Visible = 'off';
            g = rheome.show.phasor(F.ax.A, pi/3 + 0.3*randn(80, 1), 'Labels', true, 'Theme', t.theme);
            t.verifyEqual(g.mu, pi/3, 'AbsTol', 0.15);  t.verifyGreaterThan(g.R, 0.8);
            g2 = rheome.show.phasor(F.ax.A, [0; pi/2], 'At', [3 0; 6 0], 'Radius', 1, 'Amplitude', [1; 0.5]);
            t.verifyEqual(g2.R, [1; 1], 'AbsTol', 1e-12);  t.verifyEqual(g2.mu, [0; pi/2], 'AbsTol', 1e-12);
            g3 = rheome.show.phasor(F.ax.A, [0; pi]);
            t.verifyEqual(g3.R, 0, 'AbsTol', 1e-12);                       % opposite phases cancel
            nxr_save(F, t.stem('phasor'));
            t.passes(F.fig, 'phasor');
        end

        function levelstack(t)
            [h, G] = rheome.show.levelstack(4, 0:4, 'Slice', 2, 'Tile', [4 11], 'SpaceMM', 122 ./ 2.^(0:4), ...
                'Visible', 'off', 'Theme', t.theme, 'Export', t.stem('levelstack'));
            t.verifyEqual(arrayfun(@(g) numel(g.t0), G), [1 2 4 8 16]);  % level l: 2^l tiles, each two of level l+1
            t.verifyEqual(G(3).t1 - G(3).t0, ones(1, 4));
            t.verifyError(@() rheome.show.levelstack(4, 0:4, 'Tile', [2 5], 'Visible', 'off'), 'rheome:show:levelstack');
            t.passes(h, 'levelstack');
            fld = arrayfun(@(g) rand(1, numel(g.t0)), G, 'UniformOutput', false);
            h2 = rheome.show.levelstack(4, 0:4, 'Field', fld, 'Label', 'power', 'Unit', 'dB', 'Visible', 'off', 'Theme', t.theme, ...
                'Export', t.stem('levelstack-field'));
            t.passes(h2, 'levelstack-field');
        end

        function noColourLiteralsLeft(t)
            % the figure rules are tokens, not literals: no RGB triplets or hand-picked maps in +show
            d = fullfile(fileparts(fileparts(mfilename('fullpath'))), '+rheome', '+show');
            for f = dir(fullfile(d, '*.m'))'
                code = regexprep(fileread(fullfile(d, f.name)), '%[^\n]*', '');
                t.verifyEmpty(regexp(code, '\[\s*0?\.\d+\s+0?\.\d+\s+0?\.\d+', 'match', 'once'), f.name);
                t.verifyEmpty(regexp(code, '\<(parula|hot|jet|lines)\(', 'match', 'once'), f.name);
                t.verifyEmpty(regexp(code, 'flow-signed-field|''signed''', 'match', 'once'), [f.name ': retired kind']);
                % symbols and unit exponents are LaTeX maths (tokens type.math), never Unicode glyphs
                t.verifyEmpty(regexp(code, '[\x{0370}-\x{03FF}\x{00B2}\x{00B3}\x{00B5}\x{2070}-\x{209F}\x{221A}\x{2211}]', 'match', 'once'), ...
                    [f.name ': a Unicode maths glyph outside a comment']);
            end
        end
    end

    methods
        function s = stem(t, name), s = fullfile(t.out, name); end

        function passes(t, h, name)
            t.verifyTrue(isgraphics(h, 'figure'));
            T = nxr_figure_tokens(t.theme);
            for cb = findall(h, 'Type', 'ColorBar')'      % one colormap per quantity, from the key band
                p = strsplit(cb.Tag, '|');
                t.verifyTrue(numel(p) >= 5 && isfield(T.quantity, p{2}), sprintf('%s: colourbar %s has no token quantity', name, cb.Tag));
            end
            nm = getappdata(h, 'colormapNames');
            t.verifyEqual(isempty(nm), isempty(findall(h, 'Type', 'ColorBar')), [name ': colormap names for the legend']);
            t.verifyEmpty(regexp(strjoin([nm {''}]), '_r\>|PuOr|viridis|cividis', 'once'), name);
            uv = fullfile(getenv('HOME'), '.local', 'bin', 'uv');
            t.assumeTrue(isfile(t.check) && isfile(uv), 'figure_check.py or uv is not on this machine');
            % MATLAB prints a variable font as Courier: the Pergammon face (Inter) needs its static files installed
            st = [dir(fullfile(getenv('HOME'), 'Library', 'Fonts', 'Inter*Regular.ttf')); dir('/Library/Fonts/Inter*Regular.ttf')];
            t.assumeTrue(~strcmp(t.theme, 'pergammon') || ~isempty(st), ...
                'static Inter is not installed: MATLAB cannot embed the Pergammon face (the figure was drawn, not checked)');
            cmd = sprintf(['export PATH=/opt/homebrew/bin:$PATH; cd "%s" && "%s" run -q --with numpy ' ...
                '--with matplotlib --with pillow --with pdfplumber python figure_check.py "%s.pdf" --theme %s'], ...
                fileparts(t.check), uv, t.stem(name), t.theme);
            [st, txt] = system(cmd);
            t.verifyEqual(st, 0, sprintf('%s.pdf fails the figure check:\n%s', name, txt));
        end

        function o = ridgeOut(t)
            R = struct('lifetime', {}, 'speedNet', {}, 'sigmaMode', {}, 'nFrames', {}, 'pos', {}, 'v', {}, 't', {});
            for k = 1:4
                v = (1:5)' * 7 + k*90;
                R(k) = struct('lifetime', 0.01*k, 'speedNet', 0.2*k, 'sigmaMode', 0.01*k, 'nFrames', 5, ...
                              'pos', t.S.Vertices(v, :), 'v', v, 't', (1:5)');
            end
            o = struct('ridges', R, 'summary', struct('sigmas', [0.01 0.02 0.03 0.04], 'scaleHist', [3 5 2 1]));
        end
    end
end
