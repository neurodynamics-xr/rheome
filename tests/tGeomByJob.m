classdef tGeomByJob < matlab.unittest.TestCase
% rheome.geom's surface functions grouped by job (2026-11-10): measure, select, draw. Each old name forwards to the
% new one with the same outputs and warns once per session (rheome:deprecated).
%
% Author: Diellor Basha, 2026

    properties
        S       % icosphere(3), unit radius
    end

    methods (TestClassSetup)
        function sphere(tc)
            [V, F] = rheome.geom.icosphere(3);
            tc.S = struct('Vertices', V, 'Faces', F, 'nV', size(V, 1), 'Cortex', V(:, 3) < 0.9);
        end
    end

    methods (Test)
        function geodesicIsMeasureHeatDistance(tc)
            [dOld, sOld] = tc.firstCallWarnsOnce(@() rheome.geom.geodesic(tc.S, 3), 2);
            [dNew, sNew] = rheome.geom.measure.heatDistance(tc.S, 3);
            tc.verifyEqual(dOld, dNew);  tc.verifyEqual(sOld, sNew);
        end
        function heatdistanceIsMeasureHeatDistanceEigen(tc)
            V = tc.S.Vertices;  F = tc.S.Faces;
            [L, M] = rheome.operators.laplace_beltrami(V, F);
            [Phi, Lam] = eigs(-L, M, 40, 'smallestabs');  lam = abs(diag(Lam));  [lam, o] = sort(lam);
            basis = struct('Phi', Phi(:, o), 'Lambda', lam, 'Mass', M);  fg = rheome.operators.face_gradient(V, F);
            o = struct('t', 50 / max(lam));                                   % fg carries no .Vertices: give t
            pOld = tc.firstCallWarnsOnce(@() rheome.geom.heatdistance(basis, fg, [1 7], o), 1);
            tc.verifyEqual(pOld, rheome.geom.measure.heatDistanceEigen(basis, fg, [1 7], o));
        end
        function karcher_meanIsMeasureKarcherMean(tc)
            [vOld, pOld, dOld] = tc.firstCallWarnsOnce(@() rheome.geom.karcher_mean(tc.S, [1; 5; 9]), 3);
            [vNew, pNew, dNew] = rheome.geom.measure.karcherMean(tc.S, [1; 5; 9]);
            tc.verifyEqual(vOld, vNew);  tc.verifyEqual(pOld, pNew);  tc.verifyEqual(dOld, dNew);
        end
        function tilemeanIsMeasureTileMean(tc)
            P = sparse((1:tc.S.nV)', 1 + (tc.S.Vertices(:, 1) > 0), 1);   % two tiles: x <= 0, x > 0
            G = struct('P', P, 'tileOf', full(P * [1; 2]));  X = tc.S.Vertices;  a = ones(tc.S.nV, 1);
            yOld = tc.firstCallWarnsOnce(@() rheome.geom.tilemean(G, X, a), 1);
            tc.verifyEqual(yOld, rheome.geom.measure.tileMean(G, X, a));
            tc.verifyEqual(rheome.geom.tilemean(G, X, a, 'median'), rheome.geom.measure.tileMean(G, X, a, 'median'));
        end
        function cortexmaskIsSelectCortex(tc)
            cOld = tc.firstCallWarnsOnce(@() rheome.geom.cortexmask(tc.S), 1);
            tc.verifyEqual(cOld, rheome.geom.select.cortex(tc.S));  tc.verifyEqual(cOld, tc.S.Cortex);
        end
        function flipOutIsDrawFlipOut(tc)
            [pOld, iOld] = tc.firstCallWarnsOnce(@() rheome.geom.flipOut(tc.S, 1, 200), 2);
            [pNew, iNew] = rheome.geom.draw.flipOut(tc.S, 1, 200);
            tc.verifyEqual(pOld, pNew);  tc.verifyEqual(iOld.length, iNew.length);
        end
        function geodesicpathIsDrawLine(tc)
            [pOld, iOld] = tc.firstCallWarnsOnce(@() rheome.geom.geodesicpath(tc.S, 1, 200), 2);
            [pNew, iNew] = rheome.geom.draw.line(tc.S, 1, 200);
            tc.verifyEqual(pOld, pNew);  tc.verifyEqual(iOld.length, iNew.length);
            tc.verifyEqual(iOld.faces, iNew.faces);  tc.verifyEqual(iOld.vertices, iNew.vertices);
        end
    end

    methods
        function varargout = firstCallWarnsOnce(tc, f, n)
            % The first call after a clear warns rheome:deprecated; the second does not. Returns the outputs.
            clear functions %#ok<CLFUNC>   the wrappers' persistent "warned"
            varargout = cell(1, n);
            [varargout{:}] = tc.verifyWarning(f, 'rheome:deprecated');
            tc.verifyWarningFree(f);
        end
    end
end
