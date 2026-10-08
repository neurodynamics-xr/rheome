classdef flowpage
% FLOWPAGE  One band, one page: the LBO coefficients of the cortical vorticity.
%
%   B  = rheome.flowpage.prepare('sub01')          % ~40 s, once per anatomy
%   fp = rheome.flowpage(B, 'Band', [8 13], 'PageIndex', 3)
%   m  = map(fp, iFrame, iScale)                    % [V x 1] cortical map
%
% The data model behind @flowbrowser. Everything heavy -- the flow context, the fused curl
% kernel, the graph bank, the pager -- is INJECTED as a bundle, so a page costs only its own
% CWT and the model can be tested against a small synthetic anatomy.
%
% ⭐ WHY SCRUBBING IS INTERACTIVE. Once a page's coefficients C [Ks x nT] exist, a cortical
% map at any frame and scale is ONE GEMV: Phi * (h_g .* C(:,t)), ~33 MFLOP. The expensive
% stages -- the CWT and the sensors->modes GEMM -- are paid once per page, not per frame.
%
% ⭐ ONE GEMM, NOT ONE PER FILTER. The kernel is linear, so summing the in-band CWT
% coefficients in SENSOR space and applying the kernel once equals applying it per filter and
% summing: Kc*(sum_m W_m) = sum_m (Kc*W_m). Since C = 270 < Ks = 800 the sum is also done on
% the smaller index.
%
% ⚠ ONE TIME AXIS PER BAND, and it is DERIVED. Multirate gives each filter its own rate, but
% a browser needs a single axis to scrub. The band therefore takes one decimation, set by its
% HIGHEST filter so every member clears SamplesPerCycle:
%       decim = floor(fs / (SamplesPerCycle * max(cf_inband)))
% Lower filters get more than the target, which is harmless. The rate is a property of the
% band, never a free parameter.
%
% ⚠ THE ACQUISITION RATE IS A CEILING. spc*f <= fs must hold; above it the recording cannot
% supply the target and no processing recovers it. The band then falls back to decim = 1 --
% the full available rate -- rather than decimating below what it already has. At 2400 Hz
% the ceiling is 80 Hz for 30 samples/cycle, so the whole 1-60 Hz range clears it.
%
% ⚠ THE MARGIN IS NOT THE PAGE. .Core marks the samples this page OWNS; the rest are the
% wavelet cone's margin and belong to the neighbours too. A browser that lets you scrub into
% the margin without saying so is showing cone-contaminated frames as data.
%
% PROPERTIES
%   Band  PageIndex  Rate  Time  Core  Coefficients [Ks x nT]  SensorCoefficients [C x nT]
%   ScaleGains [Ks x nScale]  CenterFrequencies  NumFilters  NumScales  GraphBank  Lambda
%
% TWO ROUTES TO THE BAND -- 'Method'
%   'cwt'       sum the in-band wavelets. Exact within the frame and the validated route,
%               but it costs nFilters FULL-SPAN transforms per channel. MEASURED at 36 s
%               per page over five bands: 92% of the whole feature pass, and DELTA is the
%               most expensive band (12.1 s) despite producing the FEWEST output samples,
%               because the transform runs at the acquisition rate whatever the output rate.
%   'bandpass'  decimate to the derived rate, then one declared zero-phase bandpass, then
%               Hilbert. MEASURED ~10-30x faster.
%
% ⚠ THE TWO DIFFER IN ABSOLUTE SCALE. Summing overlapping wavelets is NOT unit gain: the
% frame contributes its own factor (measured 16.6x in energy for a 7-member alpha band
% against 2*sum(x^2)). Relative quantities -- scale profiles, fractions, time courses,
% correlations -- are unaffected; absolute energies are not comparable between methods.
% Arguably the bandpass is the better-defined object: a declared flat response against an
% incidental rippled one.
%
% ⚠ DECIMATE BEFORE FILTERING, not after. A narrow band at the acquisition rate needs a
% ~9900-tap FIR (kaiserord at 2400 Hz for a 0.8 Hz transition), and a low-order IIR there
% sits at 0.3-0.5% of Nyquist where it is numerically fragile. At the band's own derived
% rate the same band is a comfortable few percent of Nyquist.
%
% See also: rheome.flowbrowser, rheome.flow.rateplan, rheome.flow.curl, rheome.graphfilterbank, rheome.filters.firbandpass
%
% Author: Diellor Basha, 2026

    properties (SetAccess = immutable)
        Band
        PageIndex
        Rate
        Time
        Core
        Coefficients
        SensorCoefficients
        ScaleGains
        CenterFrequencies
        Lambda
        GraphBank
        Bundle
        Decim
        SamplesPerCycle
        Method
    end

    properties (Dependent, SetAccess = private)
        NumFilters
        NumScales
        NumFrames
        NumVertices
    end

    methods
        function obj = flowpage(B, varargin)
            p = inputParser;
            p.addParameter('Band',            [8 13]);
            p.addParameter('PageIndex',       1);
            p.addParameter('SamplesPerCycle', 30);
            p.addParameter('VoicesPerOctave', 10);
            p.addParameter('FreqLimits',      [1 60]);
            p.addParameter('Precision',       'single');
            p.addParameter('Method',          'cwt');      % 'cwt' | 'bandpass'
            p.parse(varargin{:});
            o = p.Results;

            obj.Bundle = B;
            obj.Band   = o.Band;
            obj.PageIndex = o.PageIndex;
            obj.SamplesPerCycle = o.SamplesPerCycle;
            fs = B.fs;

            % ---- the band's slice of the MASTER grid ----
            % The grid is anchored at the bank's TOP frequency, so a sub-bank's upper limit
            % must land on a master grid point or its centre frequencies differ from the
            % master bank's and bands are not comparable.
            probeN = ceil(8 * fs / o.FreqLimits(1));
            cfM = centerFrequencies(cwtfilterbank('SignalLength', probeN, ...
                'SamplingFrequency', fs, 'FrequencyLimits', o.FreqLimits, ...
                'VoicesPerOctave', o.VoicesPerOctave));
            inb = cfM(:).' >= o.Band(1) & cfM(:).' <= o.Band(2);
            if ~any(inb)
                error('flowpage:band', 'No master-grid filter inside [%g %g] Hz.', o.Band(1), o.Band(2));
            end
            lims = [o.Band(1), max(cfM(inb))];

            % ---- the page ----
            pr = B.pager;
            pg = page(pr, o.PageIndex);
            nSpan = numel(pg.samples);

            fb  = cwtfilterbank('SignalLength', nSpan, 'SamplingFrequency', fs, ...
                                'FrequencyLimits', lims, 'VoicesPerOctave', o.VoicesPerOctave);
            cf  = centerFrequencies(fb);
            obj.CenterFrequencies = cf(:).';

            % ---- one decimation for the whole band, set by its highest filter ----
            obj.Decim = max(1, floor(fs / (o.SamplesPerCycle * max(cf))));
            obj.Rate  = fs / obj.Decim;
            keep = 1:obj.Decim:nSpan;
            obj.Time = pg.t(keep);
            obj.Core = pg.core(keep);

            % ---- the band's analytic sensor signal, by one of two routes ----
            obj.Method = lower(o.Method);
            switch obj.Method
                case 'cwt'
                    % Sum the in-band wavelets in SENSOR space. Exact within the frame, and
                    % the route everything was validated against -- but it costs
                    % nFilters full-span transforms per channel, which MEASURED at 36 s per
                    % page across five bands and is ~92% of the whole feature pass.
                    nCh  = size(pg.F, 1);
                    Wsum = complex(zeros(nCh, numel(keep), o.Precision));
                    for c = 1:nCh
                        Wc = wt(fb, pg.F(c, :));         % [nF x nSpan]
                        Wsum(c, :) = sum(Wc(:, keep), 1);
                    end
                    obj.SensorCoefficients = Wsum;
                case 'bandpass'
                    obj.SensorCoefficients = i_bandpass(pg.F, fs, lims, obj.Decim, o.Precision);
                otherwise
                    error('flowpage:method', ...
                        'Method must be ''cwt'' or ''bandpass'', got ''%s''.', o.Method);
            end

            % ---- one GEMM to the mode domain ----
            obj.Coefficients = cast(B.Kc, o.Precision) * obj.SensorCoefficients;   % [Ks x nT]

            % ---- the spatial-scale axis: a diagonal in lambda ----
            obj.Lambda    = double(B.Lambda(:));
            obj.GraphBank = B.gfb;
            obj.ScaleGains = graphfilters(B.gfb, 'Lambda', obj.Lambda);   % [Ks x nScale]
        end

        function n = get.NumFilters(obj),  n = numel(obj.CenterFrequencies); end
        function n = get.NumScales(obj),   n = size(obj.ScaleGains, 2);      end
        function n = get.NumFrames(obj),   n = numel(obj.Time);              end
        function n = get.NumVertices(obj), n = size(obj.Bundle.Phi, 1);      end

        m  = map(obj, iFrame, iScale, mode)
        M  = maps(obj, iFrame, mode)
        Z  = complexmaps(obj, iFrame)
        E  = scalogram(obj)
        [iDom, total] = dominantScale(obj, iFrame, Z)
        E  = scaleEnergy(obj, iFrame)
        k  = centroid(obj, iFrame, basis)
        gi = globalIndex(obj, lambdaCut)
        P  = bandpower_(obj)
        S  = modeSpectrum(obj, iFrame)
        disp(obj)
    end

    methods (Static)
        B = prepare(name, varargin)
    end
end

% Author: Diellor Basha, 2026

% ---- decimate to the band's own rate, bandpass there, then form the analytic signal ----
function Z = i_bandpass(F, fs, lims, decim, precision)
    X = double(F).';                                     % [nSpan x nCh]; resample works on cols
    if decim > 1
        X  = resample(X, 1, decim);                      % polyphase, anti-aliased
        fsD = fs / decim;
    else
        fsD = fs;
    end

    lo = max(lims(1), 1e-3);  hi = min(lims(2), 0.45 * fsD);
    d  = designfilt('bandpassiir', 'FilterOrder', 8, ...
                    'HalfPowerFrequency1', lo, 'HalfPowerFrequency2', hi, ...
                    'SampleRate', fsD);
    if ~isstable(d)
        % ⚠ NOT A THEORETICAL GUARD. A narrow band at a high rate puts the poles close
        % enough to the unit circle that order 8 goes unstable; fall back rather than
        % return a plausible-looking divergent result.
        d = designfilt('bandpassiir', 'FilterOrder', 4, ...
                       'HalfPowerFrequency1', lo, 'HalfPowerFrequency2', hi, ...
                       'SampleRate', fsD);
        if ~isstable(d)
            error('flowpage:bandpass', ...
                'No stable IIR for [%g %g] Hz at %g Hz; use Method ''cwt''.', lo, hi, fsD);
        end
    end

    Y = filtfilt(d, X);                                  % zero phase
    Y = hilbert(Y);                                      % analytic, as the CWT route is
    Z = cast(Y.', precision);
end

% Author: Diellor Basha, 2026
