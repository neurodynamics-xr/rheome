classdef tOpticalflowConnection < matlab.unittest.TestCase
% rheome.dynamics.opticalflow_connection (the smoothness of opticalflow_scalar / flow.apparent) is the Cortical Flow
% Editor's optical flow: parity with the Editor's own output (tests/goldens/opticalflow-connection, made by nxr-neuro's generator
% from the Editor's code) and the Editor test's assertions (sphere bump direction and speed, gauge freedom, stationary -> 0).
%
% Author: Diellor Basha, 2026

    methods (Test)
        function matchesTheEditorsPencilAndVelocity(tc)
            G = i_golden(fullfile(fileparts(mfilename('fullpath')), 'goldens', 'opticalflow-connection', 'golden.txt'));
            tc.assertGreaterThanOrEqual(numel(G), 3);
            for m = 1:numel(G)
                g = G(m);  S = struct('Vertices', g.positions, 'Faces', g.faces + 1, 'nV', g.n);
                % the pencil: C (real 2n x 2n) and M
                [L, ~, Sd] = rheome.operators.connection_laplacian_idt(S.Vertices, S.Faces);
                Cg = sparse(g.C(:,1)+1, g.C(:,2)+1, g.C(:,3), 2*g.n, 2*g.n);
                tc.verifyLessThan(norm(L - Cg, 'fro') / norm(Cg, 'fro'), 1e-10, g.name + ": C");
                tc.verifyEqual(Sd.mass, g.M, 'RelTol', 1e-10, g.name + ": M");
                % the velocity of frame 0: the Editor's dense float64 solve (zv) to 1e-10, its float32 output to float32
                [v, info] = i_flow(g.I, S, g.alpha);
                tc.verifyEqual(info.peak, g.peak, 'RelTol', 1e-15);
                v0 = reshape(v(:,1), 3, [])';
                tc.verifyLessThan(max(abs(v0 - g.zv), [], 'all') / max(abs(g.zv), [], 'all'), 1e-10, g.name + ": frame 0 vs the dense solve");
                vg = reshape(g.velocity', 3*g.n, []);                  % [3n x nT-1]
                tc.verifyLessThan(max(abs(v - vg), [], 'all') / max(abs(vg), [], 'all'), 2e-6, g.name + ": every frame vs the Editor's float32 velocity");
            end
        end

        function sphereBumpHasTheRightDirectionAndSpeed(tc)
            [V, F] = rheome.geom.icosphere(3);  R = 0.1;  V = R*V;  n = size(V, 1);
            S = struct('Vertices', V, 'Faces', F, 'nV', n);
            [I, cen] = i_bump(V, R, 3);
            v = rheome.dynamics.opticalflow_scalar(I, S);                       % the default: connection, alpha 0.01
            [cs, sp] = i_read(V, v(:,1), cen(:,1), R);
            tc.verifyGreaterThan(cs, 0.99);
            tc.verifyGreaterThan(sp, 0.85);  tc.verifyLessThan(sp, 1.1);        % the prior pulls the speed low (Editor: 94.8 % at sub 4)
        end

        function isGaugeFree(tc)
            % the faces listed in another order move every vertex's reference edge; v does not move
            [V, F] = rheome.geom.icosphere(2);  V = 0.1*V;  n = size(V, 1);
            rng(9);  V = V .* (1 + 0.04*(rand(n, 1) - 0.5));
            [I, ~] = i_bump(V ./ vecnorm(V, 2, 2) * 0.1, 0.1, 4);
            a = rheome.dynamics.opticalflow_connection(I, struct('Vertices', V, 'Faces', F), struct('alpha', 0.1));
            b = rheome.dynamics.opticalflow_connection(I, struct('Vertices', V, 'Faces', flipud(F)), struct('alpha', 0.1));
            tc.verifyLessThan(max(abs(a - b), [], 'all') / max(abs(a), [], 'all'), 1e-6);
        end

        function aStationaryPatternGivesNoFlow(tc)
            [V, F] = rheome.geom.icosphere(2);  V = 0.1*V;
            S = struct('Vertices', V, 'Faces', F, 'nV', size(V, 1));
            [I, ~] = i_bump(V, 0.1, 1);
            v = rheome.dynamics.opticalflow_scalar(repmat(I(:,1), 1, 4), S);
            tc.verifyEqual(max(abs(v), [], 'all'), 0);
        end

        function theDefaultsAreTheEditorsAndTheFrameMethodIsGone(tc)
            [V, F] = rheome.geom.icosphere(2);  V = 0.1*V;  n = size(V, 1);
            S = struct('Vertices', V, 'Faces', F, 'nV', n);
            [I, ~] = i_bump(V, 0.1, 3);
            a = rheome.dynamics.opticalflow_scalar(I, S);                      % a call of opticalflow_connection
            tc.verifyEqual(a, rheome.dynamics.opticalflow_connection(I, S));
            R = rheome.flow.apparent(I, S, Rate=100);                           % flow.apparent: alpha 0.01
            tc.verifyEqual(R.alpha, 0.01);  tc.verifyEqual(R.velocity, a);
            tc.verifyFalse(isfield(R, 'regulariser'));
            tc.verifyError(@() rheome.flow.apparent(I, S, Rate=100, Regulariser="frame"), 'MATLAB:TooManyInputs');   % the option is gone
        end
    end
end

function [v, info] = i_flow(I, S, alpha)
    [v, info] = rheome.dynamics.opticalflow_connection(I', S, struct('alpha', alpha));   % golden I is [nT x n]
end

function [I, cen] = i_bump(V, R, nT)
% the Editor test's travelling bump: a Gaussian of geodesic width 0.2 R whose centre turns about y at 3 rad/s, dt = 0.004 s
    n = size(V, 1);  sigma = 0.2*R;  I = zeros(n, nT);  cen = zeros(3, nT);
    for k = 1:nT
        th = 3*(k-1)*0.004;  c = [R*cos(th); 0; -R*sin(th)];  cen(:,k) = c;
        d = R * acos(min(1, max(-1, (V * c) ./ (vecnorm(V, 2, 2) * R))));
        I(:,k) = exp(-d.^2 / (2*sigma^2));
    end
end

function [cs, sp] = i_read(V, v, c, R)
% median direction cosine and speed ratio against the rigid rotation's tangential part, on the vertices within 2.5 sigma
    v = reshape(v, 3, [])';  sigma = 0.2*R;
    d = R * acos(min(1, max(-1, (V * c) ./ (vecnorm(V, 2, 2) * R))));
    on = d < 2.5*sigma;
    tg = 3*0.004 * [V(:,3), zeros(size(V,1),1), -V(:,1)];
    cs = median(sum(v(on,:) .* tg(on,:), 2) ./ max(vecnorm(v(on,:), 2, 2) .* vecnorm(tg(on,:), 2, 2), 1e-30));
    sp = median(vecnorm(v(on,:), 2, 2)) / median(vecnorm(tg(on,:), 2, 2));
end

function G = i_golden(file)
% the generator's text format: "mesh name n nT alpha peak", then "name rows cols" + one line of numbers, per array
    L = splitlines(string(fileread(file)));  G = struct('name', {}, 'n', {}, 'alpha', {}, 'peak', {});  k = 1;
    while k <= numel(L)
        t = split(L(k));
        if t(1) == "mesh"
            m = numel(G) + 1;  G(m).name = t(2);  G(m).n = double(t(3));  G(m).alpha = double(t(5));  G(m).peak = double(t(6));
            k = k + 1;
            while k <= numel(L) && ~startsWith(L(k), "mesh ") && strlength(L(k)) > 0
                h = split(L(k));  r = double(h(2));  c = double(h(3));
                x = sscanf(char(L(k+1)), '%f');  G(m).(char(h(1))) = reshape(x, c, r)';  k = k + 2;
            end
        else, k = k + 1; end
    end
end
