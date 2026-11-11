classdef tFilterKernels < matlab.unittest.TestCase
% The (lambda, f) and (lambda, t) kernels of rheome.filters and the joint bank, after the unit and
% definition fixes: frequency is SIGNED HERTZ everywhere; the signed-frequency grid is right on a
% one-sided and on a full grid; wave / dampedwave / kleingordon / diffusion / schrodinger are the
% CONTINUUM forms (physical units) with the Chebyshev-angle form kept as the named 'discrete' variant;
% travwave has a sense; (lambda,t) kernels are refused by the joint bank with the reason.
%
% Author: Diellor Basha, 2026

    methods (Test)

        function jointBankPassesHertzToEveryKernel(tc)
            lam = (0:3).';  f = [-20 -5 0 5 12.5];
            [G, info] = rheome.jtv.bank(lam, f, {@(l) ones(size(l))}, {@(x) x});
            tc.verifyEqual(squeeze(G(1,:,1)), f * 1, 'AbsTol', 0);                  % psiTime got f (Hz)
            [G2, ~] = rheome.jtv.bank(lam, f, [], [], {@(l, x) x + 0*l});
            tc.verifyEqual(G2(2,:,1), f, 'AbsTol', 0);                              % K got f (Hz)
            fx = jfbFixture();
            j = rheome.jointfilterbank(fx.gfb, 'SignalLength', 32, 'SamplingFrequency', 64, ...
                    'TimeFilters', {@(x) x}, 'JointKernels', {@(l, x) 1 + 0*l + 0*x});
            W = jointfilters(j, 1);
            g1 = gain(fx.gfb, 1);
            tc.verifyEqual(W(2, :), j.Frequencies * g1(j.Lambda(2)), 'RelTol', 1e-12);
            tc.verifyError(@() rheome.jtv.bank(lam, f, [], [], [], 'type', 'ts'), 'jtv:bank:tsUnsupported');
            try, rheome.jtv.bank(lam, f, [], [], [], 'type', 'ts'); catch e, msg = e.message; end
            tc.verifySubstring(msg, 'time LAG');
        end

        function speedKernelKeepsASpeedInMetresPerSecond(tc)
            c = 3;  d = 0.5;  K = rheome.jointfilterbank.speedkernel(c, d);
            lam = [4; 25; 100];  f = linspace(0, 40, 4001);
            [~, ipk] = max(K(lam, f), [], 2);
            tc.verifyEqual(f(ipk).', c*sqrt(lam) / (2*pi), 'AbsTol', 0.011);        % f = c sqrt(lambda) / 2pi, Hz
        end

        function signedFrequenciesAreRightOnBothGrids(tc)
            fs = 100;  N = 100;  df = fs/N;
            full = (0:N-1) * df;
            s = rheome.filters.signedfrequencies(full, 'full');
            tc.verifyEqual(s, [0:49, -50:-1] * df);                                  % +fs/2 maps to -fs/2, as fftshift's order
            tc.verifyEqual(s, rheome.filters.signedfrequencies(full));              % 'full' is the default
            one = (0:N/2) * df;                                                      % 0 ... fs/2, 51 bins
            p = rheome.filters.signedfrequencies(one, 'positive');
            tc.verifyEqual(p, one);                                                  % every bin non-negative
            tc.verifyGreaterThanOrEqual(min(p), 0);
            w = rheome.filters.signedfrequencies(one, 'full');                       % what the old rule did to it
            tc.verifyLessThan(min(w), 0);                                            % (garbage: the 'full' reading is wrong here)
            tc.verifyError(@() rheome.filters.signedfrequencies(one, 'both'), 'filters:signedfrequencies:half');
            % the kernels see the same: a 20 Hz bin of a one-sided grid is +20 Hz, not -0.83 Hz
            f1 = (0:N/2) * df;  G = rheome.filters.bandpass(f1, 20, 1, 'positive');
            [~, ik] = max(G);  tc.verifyEqual(f1(ik), 20);
            Gt = rheome.filters.travwave(4, f1, 2*pi*20/2, 1, 1, 'positive');       % ridge at +20 Hz
            [~, ik] = max(Gt);  tc.verifyEqual(f1(ik), 20);
        end

        function waveIsTheContinuumDispersionOnTheSphere(tc)
            R = 0.1;  c = 12;  l = (0:12).';
            lam = l .* (l + 1) / R^2;  t = linspace(0, 0.05, 200);
            G = rheome.filters.wave(lam, t, c);
            omega_l = c * sqrt(l .* (l + 1)) / R;                                   % the exact spectrum
            tc.verifyEqual(G, cos(omega_l * t), 'AbsTol', 1e-12);
            tc.verifyEqual(G(1, :), ones(1, numel(t)));                              % lambda = 0 does not move
        end

        function discreteVariantIsTheOldWaveExactly(tc)
            lam = linspace(0, 9, 30).';  t = 0:0.1:2;  lmax = 9;  a = 1.7;  b = 0.3;  mu = 0.4;
            th = @(x) acos(max(min(1 - 2*x/lmax, 1), -1));
            tc.verifyEqual(rheome.filters.wave(lam, t, a, 'Variant', 'discrete', 'LambdaMax', lmax), cos(a * th(lam) * t), 'AbsTol', 1e-13);
            tc.verifyEqual(rheome.filters.wave(lam, t, a, 'Variant', 'discrete'), cos(a * th(lam) * t), 'AbsTol', 1e-13);   % lmax = max(lambda)
            tc.verifyEqual(rheome.filters.dampedwave(lam, t, a, b, 'Variant', 'discrete'), exp(-b*abs(t)) .* cos(a * th(lam) * t), 'AbsTol', 1e-13);
            tc.verifyEqual(rheome.filters.kleingordon(lam, t, a, mu, 'Variant', 'discrete'), cos(a * acos(max(min(1 - 2*(lam + mu)/lmax, 1), -1)) * t), 'AbsTol', 1e-13);
            tc.verifyError(@() rheome.filters.wave(lam, t, a, lmax), 'filters:wave:arguments');   % the old positional lmax
            tc.verifyError(@() rheome.filters.wave(lam, t, a, 'Variant', 'x'), ?MException);
        end

        function dampedWaveAndKleinGordonAreTheirContinuumForms(tc)
            lam = linspace(0, 400, 25).';  t = linspace(0, 1, 50);  c = 2;  beta = 1.3;  m = 7;
            tc.verifyEqual(rheome.filters.dampedwave(lam, t, c, beta), exp(-beta*t) .* cos(c * sqrt(lam) * t), 'AbsTol', 1e-14);
            KG = rheome.filters.kleingordon(lam, t, c, m);
            tc.verifyEqual(KG, cos(sqrt(c^2 * lam + m^2) * t), 'AbsTol', 1e-14);       % omega^2 = c^2 lambda + m^2
            tc.verifyEqual(KG(1, :), cos(m * t), 'AbsTol', 1e-14);                     % the DC mode rings at the rest frequency
            tc.verifyEqual(rheome.filters.kleingordon(lam, t, c, 0), rheome.filters.wave(lam, t, c), 'AbsTol', 1e-14);
            % m is not a shift of lambda: doubling m changes omega^2 by 3 m^2, whatever lambda is
            KG2 = rheome.filters.kleingordon(lam, t, c, 2*m);
            tc.verifyEqual(KG2, cos(sqrt(c^2 * lam + 4*m^2) * t), 'AbsTol', 1e-14);
        end

        function diffusionAndSchrodingerArePhysicalNotNormalised(tc)
            lam = [0; 2; 50];  t = [0 0.1 0.5];  D = 0.7;
            tc.verifyEqual(rheome.filters.diffusion(lam, t, D), exp(-D * lam * t), 'AbsTol', 1e-15);
            % a different lambda_max (more modes) must not change a mode's decay: no normalisation by the largest eigenvalue
            tc.verifyEqual(rheome.filters.diffusion([lam; 1e4], t, D), [exp(-D * lam * t); exp(-D*1e4*t)], 'AbsTol', 1e-15);
            tc.verifyEqual(rheome.filters.schrodinger(lam, t, 2), exp(-1i * lam * t / 2), 'AbsTol', 1e-15);
            tc.verifyEqual(rheome.filters.schrodinger([lam; 1e4], t, 2), [exp(-1i * lam * t / 2); exp(-1i*1e4*t/2)], 'AbsTol', 1e-15);
        end

        function travwaveHasASenseOnASignedAxis(tc)
            % a complex (rotary) coefficient of mode lambda = k^2: a wave one way is e^{+i w t}, the other e^{-i w t}
            c = 4;  k = 5;  lam = k^2;  fs = 100;  N = 200;  t = (0:N-1) / fs;
            f0 = c * k / (2*pi);   f0 = round(f0 * N / fs) * fs / N;  c = 2*pi*f0 / k;   % put the ridge on a bin
            right = exp(1i * 2*pi*f0 * t);  left = exp(-1i * 2*pi*f0 * t);
            [Cr, f] = rheome.filters.jspectrum(right, fs);  Cl = rheome.filters.jspectrum(left, fs);
            E = @(C, s) sum(abs(C .* rheome.filters.travwave(lam, f, c, 1, s)).^2) / sum(abs(C).^2);
            tc.verifyGreaterThan(E(Cr, +1), 0.99);   tc.verifyLessThan(E(Cr, -1), 1e-12);
            tc.verifyGreaterThan(E(Cl, -1), 0.99);   tc.verifyLessThan(E(Cl, +1), 1e-12);
            tc.verifyGreaterThan(E(Cr, 0), 0.99);    tc.verifyGreaterThan(E(Cl, 0), 0.99);        % sense 0: |f|, blind
            tc.verifyEqual(rheome.filters.travwave(lam, f, c, 1), rheome.filters.travwave(lam, f, c, 1, 0));
            tc.verifyError(@() rheome.filters.travwave(lam, f, c, 1, 2), 'filters:travwave:sense');
            % sense +1 plus sense -1 is sense 0 away from f = 0 (the two ridges overlap only through their tails there)
            Gp = rheome.filters.travwave(lam, f, c, 1, 1);  Gm = rheome.filters.travwave(lam, f, c, 1, -1);
            far = abs(rheome.filters.signedfrequencies(f)) > f0/2;
            D = Gp + Gm - rheome.filters.travwave(lam, f, c, 1, 0);
            tc.verifyLessThan(max(abs(D(far))), 1e-3);   % Gaussian tails of the other ridge
        end

    end
end

% Author: Diellor Basha, 2026
