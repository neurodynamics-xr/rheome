classdef tShowResolution < matlab.unittest.TestCase
% The scale ruler renders headless and refuses a malformed marks table.
%
% Author: Diellor Basha, 2026

    properties
        S; marks
    end

    methods (TestClassSetup)
        function build(t)
            [V, F] = i_ico(2);
            t.S = struct('Vertices', 0.07*V, 'Faces', F, 'nV', size(V,1));
            t.marks = table(["fine"; "mid"; "coarse"], [20; 100; 300], ...
                            ["instrument"; "floor"; "pyramid"], ...
                            'VariableNames', {'name','mm','kind'});
        end
    end

    methods (TestMethodTeardown)
        function shut(~), close all force; end
    end

    methods (Test)
        function theRulerAloneRenders(t)
            h = rheome.show.resolution(t.marks, 'Visible', 'off');
            t.verifyTrue(isgraphics(h, 'figure'));
            ax = findobj(h, 'Type', 'axes');
            t.verifyEqual(numel(ax), 1);
            t.verifyEqual(get(ax, 'XScale'), 'log');     % a ruler of lengths is logarithmic
        end

        function cortexRowsAreAddedForPatchesAndFields(t)
            nV = t.S.nV;
            h = rheome.show.resolution(t.marks, 'Surface', t.S, 'Visible', 'off', ...
                    'Patches', {1:50, 1:20}, 'PatchTitles', {'a','b'}, ...
                    'Fields', [rand(nV,1) rand(nV,1)], 'FieldTitles', {'c','d'});
            t.verifyGreaterThanOrEqual(numel(findobj(h, 'Type', 'axes')), 5);
        end

        function aPerPanelColormapIsHonoured(t)
            % ⚠ the regression: a binary scout under 'hot' is black and the panel reads empty.
            nV = t.S.nV;  f = zeros(nV,1);  f(1:30) = 1;
            blue = [linspace(0.8,0.1,64)', linspace(0.8,0.3,64)', linspace(0.9,0.7,64)'];
            h = rheome.show.resolution(t.marks, 'Surface', t.S, 'Visible', 'off', ...
                    'Fields', f, 'FieldTitles', {'scout'}, 'FieldColormaps', {blue});
            ax = findobj(h, 'Type', 'axes');
            cm = arrayfun(@(a) isequal(size(colormap(a)), [64 3]), ax);
            t.verifyTrue(any(cm));
        end

        function bandsAndFloorsAreDrawn(t)
            h = rheome.show.resolution(t.marks, 'Visible', 'off', ...
                    'Bands', {10, 80, 'low'; 80, 400, 'high'});
            t.verifyEqual(numel(findobj(h, 'Type', 'patch')), 2);
            % the one mark of kind "floor" gets a dashed guide line
            t.verifyGreaterThanOrEqual(numel(findobj(h, 'Type', 'line', 'LineStyle', '--')), 1);
        end

        function aMalformedMarksTableIsRefused(t)
            t.verifyError(@() rheome.show.resolution(table([1;2], 'VariableNames', {'mm'}), ...
                'Visible', 'off'), 'show:resolution:marks');
            t.verifyError(@() rheome.show.resolution(struct('mm', 1), 'Visible', 'off'), ...
                'show:resolution:marks');
        end
    end
end

function [V, F] = i_ico(n)
    p = (1+sqrt(5))/2;
    V = [-1 p 0; 1 p 0; -1 -p 0; 1 -p 0; 0 -1 p; 0 1 p; 0 -1 -p; 0 1 -p; p 0 -1; p 0 1; -p 0 -1; -p 0 1];
    V = V ./ vecnorm(V,2,2);
    F = [1 12 6; 1 6 2; 1 2 8; 1 8 11; 1 11 12; 2 6 10; 6 12 5; 12 11 3; 11 8 7; 8 2 9; ...
         4 10 5; 4 5 3; 4 3 7; 4 7 9; 4 9 10; 5 10 6; 3 5 12; 7 3 11; 9 7 8; 10 9 2];
    for k = 1:n
        nF = [];  M = containers.Map('KeyType','char','ValueType','double');
        for i = 1:size(F,1)
            a = F(i,1); b = F(i,2); c = F(i,3);  m = zeros(1,3);  pr = {[a b],[b c],[c a]};
            for z = 1:3
                key = sprintf('%d_%d', min(pr{z}), max(pr{z}));
                if M.isKey(key), m(z) = M(key); else
                    v = (V(pr{z}(1),:) + V(pr{z}(2),:))/2;  v = v/norm(v);
                    V = [V; v];  m(z) = size(V,1);  M(key) = m(z);   %#ok<AGROW>
                end
            end
            nF = [nF; a m(1) m(3); b m(2) m(1); c m(3) m(2); m(1) m(2) m(3)];   %#ok<AGROW>
        end
        F = nF;
    end
end
% Author: Diellor Basha, 2026
