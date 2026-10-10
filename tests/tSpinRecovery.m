classdef tSpinRecovery < matlab.unittest.TestCase
% Does a planted rotation survive the inverse, and does the inverse invent rotation from a standing
% source? The second question is the one that licenses or forbids every rotation claim from a source
% estimate.
%
% Author: Diellor Basha, 2026
    properties
        K; P; Ppinv; bA; bB; rss; fs; nT; f0
    end
    methods (TestClassSetup)
        function build(t)
            rheomeTestSubject(t, 'subject');   % skips, with the reason, when the cache or the subject is absent
            try
                B  = rheome.load.bases(rheomeTestSubject());
                st = rheome.load.study(rheomeTestSubject());
            catch
                t.assumeFail('cached test subject data not present');
            end
            S = B.L.S;  gv = double(B.L.gv(:))';
            g = rheome.operators.gauge(S.Vertices, double(S.Faces), Method="diffusion");
            ok = strcmpi(st.chan.Type,'MEG') & (st.rec.ChannelFlag(:)'==1);
            G  = double(st.hm.Gain(ok,:));  NC = double(st.ncov.NoiseCov(ok,ok));
            nVall = size(G,2)/3;
            R = rheome.inverse.mne(G, struct('NoiseCov',NC), struct('ChannelTypes',{st.chan.Type(ok)}, ...
                    'InverseMeasure','amplitude','nVert',nVall,'SnrFixed',3));
            t.K = R.ImagingKernel;
            o = rheome.flow.seedvortex(rheomeTestSubject(), WavelengthMM=140, Bases=B, Gauge=g, Check=false);
            A3 = [o.J(1:3:end) o.J(2:3:end) o.J(3:3:end)];
            B3 = rheome.filters.spin(A3, g.normal, pi/2);
            eA = zeros(nVall,3);  eA(gv,:) = A3;   JA = reshape(eA',[],1);
            eB = zeros(nVall,3);  eB(gv,:) = B3;   JB = reshape(eB',[],1);
            t.P = [JA JB];  t.Ppinv = pinv(t.P);
            t.bA = G*JA;  t.bB = G*JB;  t.rss = max(vecnorm(t.P,2,1));
            t.fs = 300;  t.nT = 600;  t.f0 = sqrt(8*16);
        end
    end
    methods
        function ic = icoh(t, cA, cB)
            tt = (0:t.nT-1)/t.fs;  env = hann(t.nT)';
            Bs = (t.bA*(cA(tt).*env) + t.bB*(cB(tt).*env)) * (1e-9/t.rss);
            ts = t.Ppinv * (t.K * Bs);
            h  = hilbert(ts.').';
            c  = mean(h(1,:).*conj(h(2,:))) / ...
                 sqrt(mean(abs(h(1,:)).^2)*mean(abs(h(2,:)).^2));
            ic = imag(c);
        end
    end
    methods (Test)
        function aPlantedRotationSurvivesTheInverse(t)
            ic = t.icoh(@(x) cos(2*pi*t.f0*x), @(x) sin(2*pi*t.f0*x));
            t.verifyGreaterThan(ic, 0.8, 'a planted rotation must come back as one');
        end
        function theSenseOfRotationIsPreserved(t)
            a = t.icoh(@(x) cos(2*pi*t.f0*x), @(x)  sin(2*pi*t.f0*x));
            b = t.icoh(@(x) cos(2*pi*t.f0*x), @(x) -sin(2*pi*t.f0*x));
            t.verifyGreaterThan(a, 0.8);
            t.verifyLessThan(b, -0.8, 'reversing the chirality must reverse the recovered sign');
        end
        function theInverseDoesNotInventRotationFromAStandingSource(t)
            % ⭐⭐ THE CONTROL THAT LICENSES THE WHOLE MEASUREMENT. MNE leakage is a real-valued linear
            % operator, hence INSTANTANEOUS, and imaginary coherency is blind to instantaneous mixing
            % -- so leakage cannot manufacture a quadrature phase. Measured 0.000 noiseless.
            ic = t.icoh(@(x) cos(2*pi*t.f0*x), @(x) zeros(size(x)));
            t.verifyLessThan(abs(ic), 0.05, ...
                'if this ever grows, rotation claims from source estimates are unsupported');
        end
    end
end
