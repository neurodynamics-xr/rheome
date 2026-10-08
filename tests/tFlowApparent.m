classdef tFlowApparent < matlab.unittest.TestCase
% The second flow pipeline: activation map, apparent motion, and the claims in the docstrings.
%
% Author: Diellor Basha, 2026

    properties
        S; nV
    end

    methods (TestClassSetup)
        function build(t)
            [V, F] = i_ico(3);                       % 642 vertices
            V = 0.07 * V;
            N = V ./ vecnorm(V, 2, 2);               % a sphere's normals are its positions
            t.S = struct('Vertices', V, 'Faces', F, 'VertNormals', N, 'nV', size(V,1));
            t.nV = size(V,1);
        end
    end

    methods (Test)
        function theNormCarriesTheOscillationAndTheEnvelopeDoesNot(t)
            % ⭐ the claim the default rests on: |J| of a real band-limited field pulses at 2f.
            fs = 300;  f0 = 11;  nT = 600;
            [J, ~] = i_wave(t.nV, nT, fs, f0);
            [~, iN] = rheome.flow.activation(real(J), Method="norm",     Rate=fs, Band=[8 16]);
            [~, iE] = rheome.flow.activation(J,       Method="analytic", Rate=fs, Band=[8 16]);
            t.verifyGreaterThan(iN.carrierRatio, 50 * iE.carrierRatio);
            t.verifyLessThan(iN.dipRatio, 0.2);            % the norm dives toward zero
            % ⭐ and the envelope recovers the PLANTED envelope's own depth, 0.5/1.5 = 1/3,
            % which is the round trip rather than an arbitrary threshold
            t.verifyEqual(iE.dipRatio, 1/3, 'RelTol', 0.05);
            t.verifyGreaterThan(iE.dipRatio, 5 * iN.dipRatio);
        end

        function hilbertAtTheSensorsEqualsHilbertAtTheSource(t)
            % ⚠ the 200x saving the script relies on; it must be exact, not approximate.
            rng(4);  nCh = 12;  nT = 512;
            K = randn(3*t.nV, nCh);
            Fb = filtfilt(ones(1,5)/5, 1, randn(nT, nCh));      % something band-ish
            Aat  = rheome.flow.activation(K * hilbert(Fb).', Rate=300, Band=[8 16]);   % sensors
            Asrc = rheome.flow.activation(K * Fb.', Method="analytic", Rate=300, Band=[8 16]); % source
            t.verifyEqual(Aat, Asrc, 'RelTol', 1e-10);
        end

        function anAlreadyComplexFieldIsNotTransformedTwice(t)
            rng(5);  J = complex(randn(3*t.nV, 64), randn(3*t.nV, 64));
            [A, i1] = rheome.flow.activation(J, Rate=300, Band=[8 16]);
            t.verifyTrue(i1.wasComplex);
            t.verifyEqual(A, abs(hilbert(0) + 0) * 0 + ...
                sqrt(abs(J(1:3:end,:)).^2 + abs(J(2:3:end,:)).^2 + abs(J(3:3:end,:)).^2), ...
                'RelTol', 1e-12);
        end

        function aTranslatingBumpGivesTheSpeedItWasGiven(t)
            % ⭐ the end-to-end check: plant a rigid translation, recover its speed.
            fs = 100;  nT = 12;
            [A, planted] = i_translate(t.S, nT, fs);
            R = rheome.flow.apparent(A, t.S, Rate=fs, Alpha=0.05);
            w = mean(A(:,1:end-1), 2);  w = w / sum(w);          % weight by where the bump is
            got = sum(w .* median(R.speed, 2));
            t.verifyEqual(got, planted, 'RelTol', 0.45);          % Horn-Schunck, aperture-limited
            t.verifyGreaterThan(got, 0);
        end

        function speedPerFrameTimesRateIsSpeed(t)
            % ⚠ the units trap: getting this wrong scales every speed by the sample rate.
            fs = 100;  A = i_translate(t.S, 6, fs);
            R = rheome.flow.apparent(A, t.S, Rate=fs);
            t.verifyEqual(R.speed, R.speedPerFrame * fs, 'RelTol', 1e-12);
            R2 = rheome.flow.apparent(A, t.S);
            t.verifyEmpty(R2.speed);                              % no rate, no m/s
            t.verifyNotEmpty(R2.speedPerFrame);
        end

        function theReadoutsHaveTheRightShapes(t)
            A = i_translate(t.S, 8, 100);
            R = rheome.flow.apparent(A, t.S, Rate=100);
            t.verifySize(R.velocity,   [3*t.nV, 7]);
            t.verifySize(R.divergence, [t.nV, 7]);
            t.verifySize(R.vorticity,  [t.nV, 7]);
            t.verifySize(R.speed,      [t.nV, 7]);
        end

        function theSmoothnessWeightChangesTheAnswerAndIsRecorded(t)
            % ⚠ Alpha is not cosmetic: a speed quoted without it is not a measurement.
            A = i_translate(t.S, 10, 100);
            lo = rheome.flow.apparent(A, t.S, Rate=100, Alpha=0.01);
            hi = rheome.flow.apparent(A, t.S, Rate=100, Alpha=100);
            t.verifyEqual(lo.alpha, 0.01);  t.verifyEqual(hi.alpha, 100);
            t.verifyGreaterThan(abs(median(lo.speed(:)) - median(hi.speed(:))), ...
                                0.05 * median(lo.speed(:)));
        end

        function badInputsAreRefusedByName(t)
            A = i_translate(t.S, 4, 100);
            t.verifyError(@() rheome.flow.apparent(complex(A, A), t.S), 'flow:apparent:map');
            t.verifyError(@() rheome.flow.apparent(-A, t.S), 'flow:apparent:map');
            t.verifyError(@() rheome.flow.apparent(A(1:end-1,:), t.S), 'flow:apparent:mesh');
            t.verifyError(@() rheome.flow.apparent(A(:,1), t.S), 'flow:apparent:frames');
            t.verifyError(@() rheome.flow.activation(randn(3*t.nV+1, 10), Rate=300), 'flow:activation:shape');
            t.verifyError(@() rheome.flow.activation(randn(3*t.nV, 4), Rate=300), 'flow:activation:short');
            % ⚠ Rate is optional and defaults to NaN; that must not trip a validator
            t.verifyWarningFree(@() rheome.flow.activation(randn(3*t.nV, 32), Method="norm"));
        end

        function activationIsNonNegativeAndPerVertex(t)
            rng(6);  J = randn(3*t.nV, 64);
            A = rheome.flow.activation(J, Method="norm");
            t.verifySize(A, [t.nV, 64]);
            t.verifyGreaterThanOrEqual(min(A(:)), 0);
            t.verifyEqual(A(1,1), norm(J(1:3,1)), 'RelTol', 1e-12);
        end
    end
end

function [J, env] = i_wave(nV, nT, fs, f0)
% A band-limited oscillation with a slowly varying envelope, as an analytic 3-vector field.
% ⚠ ONE PHASE PER VERTEX, NOT PER COMPONENT. A dipole of fixed orientation whose amplitude
% oscillates is LINEARLY polarised, and that is the physiological case; giving the three
% components independent phases makes it elliptical, and then |J| never reaches zero (it
% measured 0.178 of peak instead of ~0). The distinction is real: the zero-crossing argument
% for preferring the envelope is exact for linear polarisation and weakens as it becomes
% elliptical. Real resting alpha behaves like the linear case (dipRatio 0.034).
    tt = (0:nT-1) / fs;
    env = 1 + 0.5*sin(2*pi*0.7*tt);
    rng(1);  amp = rand(3*nV, 1) + 0.5;
    ph = repelem(2*pi*rand(nV, 1), 3, 1);
    J = (amp .* exp(1i*ph)) * (env .* exp(1i*2*pi*f0*tt));
end

function [A, speed] = i_translate(S, nT, fs)
% A Gaussian bump swept along the sphere at a known geodesic speed.
    V = S.Vertices;  R = mean(vecnorm(V,2,2));
    dtheta = 0.05;                                   % radians per frame
    speed = dtheta * R * fs;                         % m/s
    A = zeros(S.nV, nT);
    for k = 1:nT
        th = (k-1)*dtheta;
        c = R * [cos(th), sin(th), 0];
        d = vecnorm(V - c, 2, 2);
        A(:,k) = exp(-(d/(0.35*R)).^2);
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
            a=F(i,1); b=F(i,2); c=F(i,3);  m=zeros(1,3);  pr={[a b],[b c],[c a]};
            for z = 1:3
                key = sprintf('%d_%d', min(pr{z}), max(pr{z}));
                if M.isKey(key), m(z)=M(key); else
                    v=(V(pr{z}(1),:)+V(pr{z}(2),:))/2; v=v/norm(v); V=[V;v]; m(z)=size(V,1); M(key)=m(z);
                end
            end
            nF=[nF; a m(1) m(3); b m(2) m(1); c m(3) m(2); m(1) m(2) m(3)];  %#ok<AGROW>
        end
        F = nF;
    end
end
% Author: Diellor Basha, 2026
