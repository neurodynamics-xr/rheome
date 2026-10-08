classdef tRotorForward < matlab.unittest.TestCase
% Seeded rotors through the forward model, and the basis-truncation fact that came out of it.
%
% Author: Diellor Basha, 2026
    properties
        B; ctx; db
    end
    methods (TestClassSetup)
        function build(t)
            try
                t.B = rheome.load.bases(rheomeTestSubject());
                d = rheome.load.dirac(rheomeTestSubject());
                t.ctx = struct('H', t.B.L, 'd', d, 'st', rheome.load.study(rheomeTestSubject()), ...
                               'er', rheome.load.study(rheomeTestSubject('noise')));
                gv = double(t.B.L.gv(:))';  nVh = size(t.B.L.S.Vertices,1);
                cols = find(d.Hemisphere == 1);
                rows4 = reshape((gv' - 1)*4 + (1:4), [], 1);
                t.db = struct('Phi', d.Phi(rows4, cols), 'Mass', d.Mass(rows4, rows4), ...
                              'nVert', nVh, 'nModes', numel(cols));
            catch
                t.assumeFail('cached test subject / noise study data not present');
            end
        end
    end
    methods (Test)
        function thePatternTimeCourseIsSignedNotRectified(t)
            % ⚠⚠ THE REGRESSION FOR THE WORST BUG IN THIS FILE'S HISTORY. The Pattern route collapsed
            % C to a scalar with vecnorm(C,2,1), which is NON-NEGATIVE: for the oscillator family that
            % turns cos(2*pi*f0*t) into |cos|, whose fundamental is at 2*f0 = 22.6 Hz, OUTSIDE the
            % 8-16 Hz band. The band-pass removed it and every Pattern threshold came back inflated by
            % ~1500x -- which was reported as "a rotor is 2000-5000x harder to see" and was wrong.
            r = rheome.flow.seedrotor(rheomeTestSubject(), WavelengthMM=140, Bases=t.B, Check=false);
            o = rheome.forward.simulate(rheomeTestSubject(), Family="oscillator", Band=[8 16], ...
                    Pattern=r.J, NoiseTrials=10, MomentNAm=1, Context=t.ctx, Verbose=false);
            t.verifyLessThan(min(o.C), 0, 'the time course must take negative values');
        end
        function aRotorIsHarderThanAnOscillatorButNotByOrdersOfMagnitude(t)
            % ⭐ corrected: the ratio is single digits, and both are far below physiological strengths
            r = rheome.flow.seedrotor(rheomeTestSubject(), WavelengthMM=140, Bases=t.B, Check=false);
            a = rheome.forward.simulate(rheomeTestSubject(), Family="oscillator", Band=[8 16], ...
                    Pattern=r.J, NoiseTrials=10, MomentNAm=1, Context=t.ctx, Verbose=false);
            b = rheome.forward.simulate(rheomeTestSubject(), Family="oscillator", Band=[8 16], ...
                    WavelengthMM=140, NoiseTrials=10, MomentNAm=1, Context=t.ctx, Verbose=false);
            rat = a.threshMomentNAm / b.threshMomentNAm;
            t.verifyGreaterThan(rat, 1);
            t.verifyLessThan(rat, 50, 'if this exceeds 50 the rectification bug is back');
            t.verifyLessThan(a.threshMomentNAm, 10, 'a 140 mm rotor is detectable at physiological nAm');
            t.verifyEqual(a.route, "pattern");
            t.verifyEqual(b.route, "modes");
        end
        function aQuadraturePairMakesTheFieldSpin(t)
            % ⭐ a spin about the normal is multiplication of the complex tangent field by exp(i*theta),
            % so z*exp(i*w*t) = cos(w*t)*z + sin(w*t)*(i*z): the vortex and its pi/2-spun partner in
            % quadrature. The forward stays RANK TWO and the chirality reverses every half turn.
            B = t.ctx.H;
            g = rheome.operators.gauge(B.S.Vertices, double(B.S.Faces), Method="diffusion");
            o = rheome.flow.seedvortex(rheomeTestSubject(), WavelengthMM=140, Bases=struct('L',B), ...
                                Gauge=g, Check=false);
            JA = [o.J(1:3:end) o.J(2:3:end) o.J(3:3:end)];
            JB = rheome.filters.spin(JA, g.normal, pi/2);
            P  = [reshape(JA',[],1) reshape(JB',[],1)];
            q = rheome.forward.simulate(rheomeTestSubject(), Band=[8 16], Pattern=P, Duration=2, ...
                    NoiseTrials=10, MomentNAm=1, Context=t.ctx, Verbose=false);
            t.verifyTrue(q.spinning);
            t.verifyEqual(q.spinRateHz, sqrt(8*16), 'RelTol', 1e-9, 'defaults to the band centre');
            t.verifyEqual(rank(q.B, 1e-9*norm(q.B)), 2, 'a quadrature pair gives a rank-2 B');
        end
        function spinningAndStandingAreNearlyEquallyDetectable(t)
            % ⚠ 0.78 against 0.84 nAm -- 8%. So AMPLITUDE alone cannot tell a rotating vortex from a
            % standing one; the discriminator is the quadrature phase between the two topographies.
            B = t.ctx.H;
            g = rheome.operators.gauge(B.S.Vertices, double(B.S.Faces), Method="diffusion");
            o = rheome.flow.seedvortex(rheomeTestSubject(), WavelengthMM=140, Bases=struct('L',B), ...
                                Gauge=g, Check=false);
            JA = [o.J(1:3:end) o.J(2:3:end) o.J(3:3:end)];
            JB = rheome.filters.spin(JA, g.normal, pi/2);
            a = rheome.forward.simulate(rheomeTestSubject(), Band=[8 16], Pattern=reshape(JA',[],1), ...
                    Family="oscillator", Duration=2, NoiseTrials=10, MomentNAm=1, ...
                    Context=t.ctx, Verbose=false);
            b = rheome.forward.simulate(rheomeTestSubject(), Band=[8 16], ...
                    Pattern=[reshape(JA',[],1) reshape(JB',[],1)], Duration=2, NoiseTrials=10, ...
                    MomentNAm=1, Context=t.ctx, Verbose=false);
            t.verifyLessThan(abs(log2(a.threshMomentNAm/b.threshMomentNAm)), 0.5);
        end
        function thePatternRouteIsTheVertexLeadfieldAndIsRankOne(t)
            r = rheome.flow.seedrotor(rheomeTestSubject(), WavelengthMM=140, Bases=t.B, Check=false);
            o = rheome.forward.simulate(rheomeTestSubject(), Pattern=r.J, NoiseTrials=5, MomentNAm=1, ...
                                 Context=t.ctx, Verbose=false);
            t.verifyEqual(o.routeRelErr, 0);            % identity, not a measurement
            t.verifyEqual(rank(o.B, 1e-9*norm(o.B)), 1, 'a separable pattern gives a rank-1 B');
        end
        function theDiracBasisRetainsAboutItsDimensionRatioOfAnyVertexField(t)
            % ⚠⚠ NOT a rotor-specific fact, which is how I first read it. A rotor retains 0.0139 of
            % its energy, a purely NORMAL field of the same scale 0.0201, and the tangential
            % projection changes nothing -- all about 400/(3*10242) = 0.013, the dimension ratio.
            % ⭐ So a field CONSTRUCTED in the vertex domain cannot be tested through the Dirac
            % projection; use the vertex leadfield. The sensor-driven pipeline is unaffected because
            % it synthesises fields that are in the span by construction.
            r = rheome.flow.seedrotor(rheomeTestSubject(), WavelengthMM=140, Bases=t.B, Check=false);
            Nv = t.B.L.S.VertNormals;  Nv = Nv ./ vecnorm(Nv, 2, 2);
            Jn = reshape((r.Psi .* Nv)', [], 1);                    % a purely normal field
            er = norm(rheome.forward.project(r.J, t.db))^2 / norm(r.J)^2;
            en = norm(rheome.forward.project(Jn,  t.db))^2 / norm(Jn)^2;
            t.verifyLessThan(er, 0.05);
            t.verifyLessThan(en, 0.05);
            t.verifyLessThan(abs(log2(er/en)), 1.5, ...
                'a rotor and a normal field must be represented about equally badly');
        end
        function anInBasisFieldRoundTripsInShapeButNotInEnergy(t)
            % ⚠ the ceiling for ANY pure 3-vector: the quaternion REAL part is dropped, so even a
            % field built from the basis retains only ~0.57 of its energy while its shape is exact.
            rng(1);  W = t.db.Phi * randn(t.db.nModes, 1);
            J = reshape([W(2:4:end) W(3:4:end) W(4:4:end)]', [], 1);
            Jp = rheome.forward.project(J, t.db);
            t.verifyGreaterThan(abs(corr(Jp, J)), 0.99, 'the shape must survive');
            t.verifyLessThan(norm(Jp)^2/norm(J)^2, 0.75, 'the energy must not');
        end
    end
end
