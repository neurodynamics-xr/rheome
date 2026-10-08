classdef tFlowPage < matlab.unittest.TestCase
% @flowpage -- the data model behind the browser. Given a band and a page it holds the LBO
% coefficients of the vorticity and serves per-frame products.
%
% Everything heavy (context, flow kernel, graph bank, pager) is INJECTED as a bundle, so the
% model can be tested against a small synthetic anatomy in milliseconds instead of the ~40 s
% a real rheome.flow.context takes. rheome.flowpage.prepare builds the real bundle.
%
% Author: Diellor Basha, 2026

    properties
        B                       % the injected bundle
        C = 8; V = 40; Ks = 12; nT = 400; fs = 400;
    end

    methods (TestClassSetup)
        function makeBundle(tc)
            d = tempname; mkdir(d);
            tc.addTeardown(@() rmdir(d, 's'));
            f = fullfile(d, 'recording.mat');

            rng(11);
            F = randn(tc.C, tc.nT);                                   %#ok<NASGU>
            Time = (0:tc.nT-1)/tc.fs;  sfreq = tc.fs;                 %#ok<NASGU>
            ChannelFlag = ones(tc.C,1);                               %#ok<NASGU>
            ChannelName = compose('C%d', 1:tc.C);                     %#ok<NASGU>
            ChannelType = repmat({'MEG'}, 1, tc.C);                   %#ok<NASGU>
            nCh = tc.C; nT = tc.nT; Comment = 'synth'; source = '';   %#ok<NASGU>
            save(f, 'F','Time','sfreq','ChannelFlag','ChannelName','ChannelType', ...
                 'nCh','nT','Comment','source', '-v7.3', '-nocompression');

            % a small synthetic anatomy: orthonormal Phi, increasing Lambda
            [Phi, ~] = qr(randn(tc.V, tc.Ks), 0);
            Lambda = linspace(1, 400, tc.Ks)';
            tc.B = struct();
            tc.B.pager  = rheome.pagedrecording(f, 'PageLength', 200, 'Overlap', 50, 'Precision','single');
            tc.B.Kc     = randn(tc.Ks, tc.C);
            tc.B.Phi    = Phi;
            tc.B.Lambda = Lambda;
            tc.B.fs     = tc.fs;
            tc.B.gfb    = rheome.graphfilterbank(max(Lambda), 'NumFilters', 3);
            tc.B.wVert  = ones(tc.V,1);
        end
    end

    methods (Test)

        function bandSelectionAndRateClearSamplesPerCycle(tc)
            % The whole band shares ONE time axis so the browser can scrub it, and the rate
            % must clear samplesPerCycle for the HIGHEST in-band filter -- the tightest one.
            fp = rheome.flowpage(tc.B, 'Band', [10 20], 'PageIndex', 1, 'SamplesPerCycle', 10);
            tc.verifyGreaterThan(fp.NumFilters, 0);
            tc.verifyGreaterThanOrEqual(fp.Rate / max(fp.CenterFrequencies), 10 - 1e-9);
            tc.verifyEqual(numel(fp.Time), size(fp.Coefficients, 2));
        end

        function theAcquisitionRateIsACeilingAndTheBandFallsBackToIt(tc)
            % 30 samples/cycle needs spc*f <= fs. Above that ceiling the recording cannot
            % supply the target and no processing recovers it: the band must fall back to
            % the FULL available rate (decim 1) rather than silently decimating below what
            % it already has. At fs = 400 a 40-80 Hz band at 10/cycle wants 600+ Hz.
            fp = rheome.flowpage(tc.B, 'Band', [40 80], 'PageIndex', 1, 'SamplesPerCycle', 10);
            tc.verifyEqual(fp.Decim, 1);
            tc.verifyEqual(fp.Rate, tc.fs, 'AbsTol', 1e-9);
            tc.verifyLessThan(fp.Rate / max(fp.CenterFrequencies), 10);   % genuinely short
        end

        function coefficientsAreTheKernelAppliedToTheSummedBand(tc)
            % Kc is linear, so summing the in-band CWT coefficients in SENSOR space and
            % applying the kernel once must equal applying it per filter and summing. That
            % identity is what makes this one GEMM instead of nF.
            fp = rheome.flowpage(tc.B, 'Band', [40 80], 'PageIndex', 1, 'SamplesPerCycle', 10);
            tc.verifySize(fp.Coefficients, [tc.Ks, numel(fp.Time)]);
            tc.verifyEqual(double(fp.Coefficients), ...
                tc.B.Kc * double(fp.SensorCoefficients), 'RelTol', 1e-4);
        end

        function mapIsTheScaleFilteredSynthesis(tc)
            fp = rheome.flowpage(tc.B, 'Band', [40 80], 'PageIndex', 1, 'SamplesPerCycle', 10);
            g = 2;  iT = 7;
            expected = real(tc.B.Phi * (fp.ScaleGains(:,g) .* double(fp.Coefficients(:,iT))));
            tc.verifyEqual(double(map(fp, iT, g)), expected, 'RelTol', 1e-4);
            tc.verifySize(map(fp, iT, g), [tc.V 1]);
        end

        function magnitudeIsTheEnvelopeNotTheAbsoluteOfTheSignedMap(tc)
            % ⚠ THESE ARE DIFFERENT QUANTITIES AND THE DIFFERENCE IS THE POINT.
            % 'signed' is Re(z), the INSTANTANEOUS vorticity -- it flickers with the carrier
            % phase. 'magnitude' is |z|, the ENVELOPE -- phase-invariant, which is why it is
            % steady to watch. |z| >= |Re z| always, with equality only where the imaginary
            % part vanishes, so asserting magnitude == abs(signed) would be asserting that
            % the analytic signal is real.
            fp = rheome.flowpage(tc.B, 'Band', [40 80], 'PageIndex', 1, 'SamplesPerCycle', 10);
            sg = map(fp, 7, 2, 'signed');
            mg = map(fp, 7, 2, 'magnitude');
            tc.verifyTrue(any(sg < 0), 'signed vorticity must carry both handednesses');
            tc.verifyTrue(all(mg >= 0));
            tc.verifyGreaterThanOrEqual(double(mg), abs(double(sg)) - 1e-6);
            z = tc.B.Phi * (double(fp.ScaleGains(:,2)) .* double(fp.Coefficients(:,7)));
            tc.verifyEqual(double(mg), abs(z), 'RelTol', 1e-4);
            tc.verifyGreaterThan(max(double(mg) - abs(double(sg))), 1e-3);  % genuinely distinct
        end

        function batchedMapsEqualTheLoopedOnes(tc)
            % maps() is one GEMM where map() is nScale GEMVs -- same arithmetic, and the
            % browser calls it every frame, so a divergence here would be invisible on
            % screen and wrong in every figure.
            fp = rheome.flowpage(tc.B, 'Band', [40 80], 'PageIndex', 1, 'SamplesPerCycle', 10);
            % ⚠ TOLERANCE IS SINGLE-PRECISION SUMMATION ORDER, not slack. The GEMM and the
            % GEMVs accumulate 800-term dot products in different orders, which in single
            % gives ~3e-5 relative on terms with cancellation. 1e-3 still catches anything
            % structural -- a wrong scale, transposed gains, an off-by-one column.
            M = maps(fp, 9);
            tc.verifySize(M, [tc.V, fp.NumScales + 1]);
            for g = 0:fp.NumScales
                tc.verifyEqual(double(M(:, g+1)), double(map(fp, 9, g)), 'RelTol', 1e-3);
            end
            tc.verifyEqual(double(maps(fp, 9, 'magnitude')), abs(double(complexZ(fp, 9))), ...
                'RelTol', 1e-4);
        end

        function bothMethodsAgreeOnRelativeMeasuresButNotAbsolute(tc)
            % 'cwt' sums the in-band wavelets; 'bandpass' decimates and applies one declared
            % zero-phase filter. They must describe the SAME band -- same derived rate, same
            % normalised scale profile -- while differing in ABSOLUTE scale, because summing
            % overlapping wavelets is not unit gain.
            %
            % MEASURED on real data: power time courses agree at Spearman 0.99+, scale
            % centroids at 0.94-0.99, the normalised scale profile to 0.02 percentage
            % points, and the frame gain runs ~19-24x. Anything that reads a RATIO is
            % interchangeable between them; anything that reads an absolute energy is not.
            args = {'Band', [40 80], 'PageIndex', 1, 'SamplesPerCycle', 10};
            fA = rheome.flowpage(tc.B, args{:}, 'Method', 'cwt');
            fB = rheome.flowpage(tc.B, args{:}, 'Method', 'bandpass');
            tc.verifyEqual(fA.Method, 'cwt');
            tc.verifyEqual(fB.Method, 'bandpass');
            tc.verifyEqual(fB.Rate, fA.Rate, 'AbsTol', 1e-9);      % the rate is DERIVED, shared

            EA = mean(scalogram(fA), 2);  EA = EA / sum(EA);
            EB = mean(scalogram(fB), 2);  EB = EB / sum(EB);
            tc.verifyEqual(EB, EA, 'AbsTol', 0.06);                % same profile
            tc.verifyGreaterThan(corr(EA, EB), 0.95);

            gain = mean(bandpower_(fA)) / mean(bandpower_(fB));
            tc.verifyGreaterThan(gain, 1.5);                       % the frame gain is real
        end

        function anUnknownMethodErrors(tc)
            tc.verifyError(@() rheome.flowpage(tc.B, 'Band', [40 80], 'PageIndex', 1, ...
                'SamplesPerCycle', 10, 'Method', 'wavelet'), 'flowpage:method');
        end

        function signedAndEnvelopeComeFromOneSynthesis(tc)
            % complexmaps is the single source: real(Z) is the signed map, abs(Z) the
            % envelope, and dominantScale reads the same Z. If these drifted apart the
            % browser would show three views of subtly different data at the same frame.
            fp = rheome.flowpage(tc.B, 'Band', [40 80], 'PageIndex', 1, 'SamplesPerCycle', 10);
            Z = complexmaps(fp, 9);
            tc.verifySize(Z, [tc.V, fp.NumScales + 1]);
            tc.verifyFalse(isreal(Z));
            tc.verifyEqual(double(maps(fp, 9, 'signed')),    double(real(Z)), 'RelTol', 1e-6);
            tc.verifyEqual(double(maps(fp, 9, 'magnitude')), double(abs(Z)),  'RelTol', 1e-6);
        end

        function scalogramAgreesWithPerFrameScaleEnergy(tc)
            % The vectorised whole-page version must equal the per-frame one it replaces.
            fp = rheome.flowpage(tc.B, 'Band', [40 80], 'PageIndex', 1, 'SamplesPerCycle', 10);
            E = scalogram(fp);
            tc.verifySize(E, [fp.NumScales, fp.NumFrames]);
            for i = [1 7 fp.NumFrames]
                tc.verifyEqual(E(:, i).', scaleEnergy(fp, i), 'RelTol', 1e-6);
            end
            tc.verifyTrue(all(E(:) >= 0));
        end

        function dominantScaleIsAnArgmaxOverScalesPerVertex(tc)
            fp = rheome.flowpage(tc.B, 'Band', [40 80], 'PageIndex', 1, 'SamplesPerCycle', 10);
            [iDom, total] = dominantScale(fp, 9);
            tc.verifySize(iDom, [tc.V 1]);
            tc.verifyTrue(all(iDom >= 1 & iDom <= fp.NumScales));
            tc.verifyTrue(all(total >= 0));
            % passing a precomputed Z must not change the answer -- that reuse is what keeps
            % the browser at one synthesis per frame
            Z = complexmaps(fp, 9);
            [iDom2, total2] = dominantScale(fp, 9, Z);
            tc.verifyEqual(iDom2, iDom);
            tc.verifyEqual(total2, total, 'RelTol', 1e-9);
            % and it must actually be the argmax
            P = abs(double(Z(:, 2:end))).^2;
            [~, ref] = max(P, [], 2);
            tc.verifyEqual(iDom, ref);
        end

        function scaleZeroSumsTheBankRatherThanFilteringOnce(tc)
            fp = rheome.flowpage(tc.B, 'Band', [40 80], 'PageIndex', 1, 'SamplesPerCycle', 10);
            total = zeros(tc.V, 1);
            for g = 1:fp.NumScales, total = total + double(map(fp, 5, g)); end
            tc.verifyEqual(double(map(fp, 5, 0)), total, 'RelTol', 1e-4);
        end

        function scaleEnergyIsNonNegativeAndTracksTheSpectrum(tc)
            fp = rheome.flowpage(tc.B, 'Band', [40 80], 'PageIndex', 1, 'SamplesPerCycle', 10);
            E = scaleEnergy(fp, 5);
            tc.verifySize(E, [1 fp.NumScales]);
            tc.verifyTrue(all(E >= 0));
            tc.verifyGreaterThan(sum(E), 0);
        end

        function centroidLiesInsideTheBanksWavenumberRange(tc)
            fp = rheome.flowpage(tc.B, 'Band', [40 80], 'PageIndex', 1, 'SamplesPerCycle', 10);
            k = centroid(fp, 5);
            kc = centerWavenumbers(fp.GraphBank);
            tc.verifyGreaterThanOrEqual(k, min(kc) - 1e-9);
            tc.verifyLessThanOrEqual(k, max(kc) + 1e-9);
        end

        function globalIndexIsAFractionOverTheWholePage(tc)
            % The running global/local readout: fraction of spectral energy below a lambda
            % cut. A fraction, so it must stay in [0,1] at every frame.
            fp = rheome.flowpage(tc.B, 'Band', [40 80], 'PageIndex', 1, 'SamplesPerCycle', 10);
            gi = globalIndex(fp);
            tc.verifySize(gi, [1 numel(fp.Time)]);
            tc.verifyTrue(all(gi >= -1e-9 & gi <= 1 + 1e-9));
        end

        function bandPowerIsPerFrameAndNonNegative(tc)
            fp = rheome.flowpage(tc.B, 'Band', [40 80], 'PageIndex', 1, 'SamplesPerCycle', 10);
            P = bandpower_(fp);
            tc.verifySize(P, [1 numel(fp.Time)]);
            tc.verifyTrue(all(P >= 0));
        end

        function coreMaskMarksTheSamplesThisPageOwns(tc)
            % The browser must not show margin samples as if they were the page's own --
            % they belong to the neighbouring pages and carry cone contamination.
            fp = rheome.flowpage(tc.B, 'Band', [40 80], 'PageIndex', 1, 'SamplesPerCycle', 10);
            tc.verifySize(fp.Core, [1 numel(fp.Time)]);
            tc.verifyTrue(any(fp.Core) && ~all(fp.Core));
        end

        function frameIndexOutOfRangeErrors(tc)
            fp = rheome.flowpage(tc.B, 'Band', [40 80], 'PageIndex', 1, 'SamplesPerCycle', 10);
            tc.verifyError(@() map(fp, 0, 1), 'flowpage:frame');
            tc.verifyError(@() map(fp, numel(fp.Time)+1, 1), 'flowpage:frame');
        end

        function scaleIndexOutOfRangeErrors(tc)
            fp = rheome.flowpage(tc.B, 'Band', [40 80], 'PageIndex', 1, 'SamplesPerCycle', 10);
            tc.verifyError(@() map(fp, 3, fp.NumScales+1), 'flowpage:scale');
        end

    end
end

function Z = complexZ(fp, iT)
    G = [sum(fp.ScaleGains,2), fp.ScaleGains];
    Z = fp.Bundle.Phi * (cast(G,'like',fp.Coefficients) .* fp.Coefficients(:,iT));
end

% Author: Diellor Basha, 2026
