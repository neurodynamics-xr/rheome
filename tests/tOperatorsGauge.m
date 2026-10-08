classdef tOperatorsGauge < matlab.unittest.TestCase
% Tests for rheome.operators.gauge: the frame is a frame, the parallel one is parallel, the smooth one
% is smooth, the trivial one is singular exactly where asked, and the index budget (chi) is
% enforced rather than silently absorbed.
%
% Author: Diellor Basha, 2026
    properties
        V; F; C
    end
    methods (TestClassSetup)
        function build(t)
            [t.V, t.F] = i_sphere(3);              % closed, chi = 2
            t.C = rheome.operators.connection_laplacian(t.V, t.F);
        end
    end
    methods (Test)
        function everyMethodReturnsAnOrthonormalFrame(t)
            for m = ["diffusion" "smoothest" "trivial"]
                g = rheome.operators.gauge(t.V, t.F, Method=m, Connection=t.C);
                t.verifyLessThan(max(abs(vecnorm(g.e1,2,2)-1)), 1e-10, m);
                t.verifyLessThan(max(abs(sum(g.e1.*g.e2,2))), 1e-10, m);
                t.verifyLessThan(max(abs(sum(g.e1.*g.normal,2))), 1e-8, m);
            end
        end
        function trivialConnectionIsExactlyParallel(t)
            g = rheome.operators.gauge(t.V, t.F, Method="trivial", Connection=t.C);
            t.verifyLessThan(g.residual, 1e-9, ...
                'a trivial connection must integrate to a frame with zero transport residual');
        end
        function diffusionIsSmootherThanTheHalfEdgeFrame(t)
            g = rheome.operators.gauge(t.V, t.F, Connection=t.C);
            [ii,jj] = find(triu(t.C.A~=0,1));
            raw = median(acosd(max(-1,min(1,sum(t.C.e1(ii,:).*t.C.e1(jj,:),2)))));
            t.verifyLessThan(g.neighbourAngle, raw/2, ...
                'the whole point is that connection .e1 is not a smooth gauge');
        end
        function diffusionConvergesToTheSmoothestField(t)
            % the claim worth testing is not that t does not matter -- on a coarse sphere it
            % moves the frame by 3.9 deg -- but that large t reaches the connection Laplacian's
            % own lowest eigenvector, which is what makes the cheap solve the right one.
            b = rheome.operators.gauge(t.V, t.F, Time=1e3, Connection=t.C);
            e = rheome.operators.gauge(t.V, t.F, Method="smoothest", Connection=t.C);
            t.verifyLessThan(abs(b.neighbourAngle - e.neighbourAngle), 1, ...
                'diffusion at large t must land on the smoothest field');
        end
        function theTwoResidualsAreNotTheSameQuestion(t)
            % ⚠ .residual is against the gauge's own connection, .residualLC against Levi-Civita.
            % For a trivial connection the first is 0 and the second is not.
            g = rheome.operators.gauge(t.V, t.F, Method="trivial", Connection=t.C);
            t.verifyLessThan(g.residual, 1e-9);
            t.verifyGreaterThan(g.residualLC, 1e-3);
        end
        function poincareHopfHolds(t)
            % ⭐ on a closed surface the winding of ANY tangent field sums to chi = 2
            for m = ["diffusion" "smoothest" "trivial"]
                g = rheome.operators.gauge(t.V, t.F, Method=m, Connection=t.C);
                t.verifyEqual(sum(g.charge), 2, m);
            end
        end
        function singularitiesGoWhereAskedAndNowhereElse(t)
            % ⭐⭐ REGRESSION FOR THE UNWRAPPED SOLVE. .singular is MEASURED from the frame, and it
            % must be exactly the prescription. Solving against d1*argE without the wrap gave index
            % k - n_f and wound on 2560 cortex faces while this list echoed the prescription.
            g = rheome.operators.gauge(t.V, t.F, Method="trivial", Singular=[7 900], ...
                                SingularCharge=[1 1], Connection=t.C);
            t.verifyEqual(sort(g.singular(:)'), [7 900]);
            t.verifyEqual(g.charge(:)', [1 1]);
        end
        function aPairNeedsTheBudgetElsewhere(t)
            g = rheome.operators.gauge(t.V, t.F, Method="trivial", Singular=[7 900 300 1200], ...
                                SingularCharge=[1 -1 1 1], Connection=t.C);
            [s, i] = sort(g.singular(:)');
            t.verifyEqual(s, [7 300 900 1200]);
            c = g.charge(:)';  t.verifyEqual(c(i), [1 1 -1 1]);
        end
        function aBudgetOtherThanChiIsRejected(t)
            % ⚠ the trap the error exists for: it used to solve "successfully" and put the remainder
            % on the pinned face. [+1 -1] was once the REQUIRED form; it is now the wrong one.
            t.verifyError(@() rheome.operators.gauge(t.V, t.F, Method="trivial", Singular=[7 900], ...
                SingularCharge=[1 -1], Connection=t.C), 'operators:gauge:budget');
        end
        function trivialIsNoRougherThanDiffusion(t)
            % ⚠ the old docstring called the parallel frame "ambiently roughest" -- that was the bug
            p = i_poles(t.V, t.F);
            g = rheome.operators.gauge(t.V, t.F, Method="trivial", Singular=p, Connection=t.C);
            d = rheome.operators.gauge(t.V, t.F, Connection=t.C);
            t.verifyLessThan(g.neighbourAngle, d.neighbourAngle + 0.5);
            t.verifyLessThan(g.residualLC, d.residualLC + 0.01);
        end
        function polesAtThePolesGiveTheGeographicFrame(t)
            % ⭐ on the round sphere the parallel frame with +1 at both poles is the east field
            p = i_poles(t.V, t.F);
            g = rheome.operators.gauge(t.V, t.F, Method="trivial", Singular=p, Connection=t.C);
            t.verifyLessThan(i_eastdev(t.V, g.e1, g.normal), 2, 'median deg from east, colat 15-165');
        end
        function matchesNxrComputeTrivialFaceField(t)
            % ⭐ the same Crane 2010 equation solved by nxr-compute (geometry-central) on the PRIMAL
            % complex: a face field with vertex singularities. Poles at a vertex of each polar face;
            % the two fields agree up to a constant rotation away from the poles.
            t.assumeEqual(exist('nxr_compute', 'file'), 3, 'nxr_compute MEX not on the path');
            V = t.V; F = t.F;  p = i_poles(V, F);  pv = F(p, 1);
            g = rheome.operators.gauge(V, F, Method="trivial", Singular=p, Connection=t.C);
            h = nxr_compute('create', V, F);  c = onCleanup(@() nxr_compute('destroy', h));
            G = nxr_compute('geometry', h);
            T = nxr_compute('gauge', h, 'trivial', struct('singVerts', uint32(pv(:)), 'singValues', [1;1]));
            ef = real(T.face.rotation .* G.face.grid);
            N = cross(V(F(:,2),:)-V(F(:,1),:), V(F(:,3),:)-V(F(:,1),:), 2);  N = N./vecnorm(N,2,2);
            ev = g.e1(F(:,1),:) + g.e1(F(:,2),:) + g.e1(F(:,3),:);  ev = ev - sum(ev.*N,2).*N;
            ef = ef - sum(ef.*N,2).*N;
            a = atan2(sum(cross(ev,ef,2).*N,2), sum(ev.*ef,2));
            dv = abs(mod(a - angle(mean(exp(1i*a))) + pi, 2*pi) - pi);
            fc = (V(F(:,1),:)+V(F(:,2),:)+V(F(:,3),:))/3;  far = abs(fc(:,3)./vecnorm(fc,2,2)) < cosd(15);
            t.verifyLessThan(rad2deg(median(dv(far))), 2, 'median deg between the two fields');
        end
    end
end

function [V,F] = i_sphere(n)
    [V,F] = i_ico();
    for k = 1:n
        nE = containers.Map('KeyType','char','ValueType','double');
        Fn = zeros(0,3);
        for f = 1:size(F,1)
            v = F(f,:);  m = zeros(1,3);
            for e = 1:3
                a = v(e);  b = v(mod(e,3)+1);  key = sprintf('%d_%d', min(a,b), max(a,b));
                if ~isKey(nE, key)
                    V(end+1,:) = (V(a,:)+V(b,:))/2;  nE(key) = size(V,1); %#ok<AGROW>
                end
                m(e) = nE(key);
            end
            Fn = [Fn; v(1) m(1) m(3); m(1) v(2) m(2); m(3) m(2) v(3); m(1) m(2) m(3)]; %#ok<AGROW>
        end
        F = Fn;
    end
    V = V ./ vecnorm(V,2,2);
end

function p = i_poles(V, F)
    fc = (V(F(:,1),:) + V(F(:,2),:) + V(F(:,3),:))/3;
    [~, a] = max(fc(:,3));  [~, b] = min(fc(:,3));  p = [a b];
end

function d = i_eastdev(V, e1, N)
    u = V./vecnorm(V,2,2);
    eE = [-u(:,2) u(:,1) 0*u(:,1)];  eE = eE./max(vecnorm(eE,2,2), eps);  eQ = cross(N, eE, 2);
    a = atan2(sum(e1.*eQ,2), sum(e1.*eE,2));  cl = acosd(u(:,3));  w = cl > 15 & cl < 165;
    d = rad2deg(median(abs(mod(a(w) - angle(mean(exp(1i*a(w)))) + pi, 2*pi) - pi)));
end

function [V,F] = i_ico()
    t = (1+sqrt(5))/2;
    V = [-1 t 0; 1 t 0; -1 -t 0; 1 -t 0; 0 -1 t; 0 1 t; 0 -1 -t; 0 1 -t; ...
          t 0 -1; t 0 1; -t 0 -1; -t 0 1];
    V = V ./ vecnorm(V,2,2);
    F = [1 12 6; 1 6 2; 1 2 8; 1 8 11; 1 11 12; 2 6 10; 6 12 5; 12 11 3; 11 8 7; 8 2 9; ...
         4 10 5; 4 5 3; 4 3 7; 4 7 9; 4 9 10; 5 10 6; 3 5 12; 7 3 11; 9 7 8; 10 9 2];
end
