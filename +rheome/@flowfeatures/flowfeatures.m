classdef flowfeatures
% FLOWFEATURES  Windowed, band- and scale-stratified tabulation of the cortical flow tensor.
%
%   fe  = rheome.flowfeatures(B)
%   fe  = rheome.flowfeatures(B, 'Bands', {'alpha','beta'}, 'WindowSec', 1)
%   out = extract(fe, pageIndex)
%   T   = totable(fe, out)
%
% The bookkeeping layer: one row per (time window, band, spatial scale, cortical ROI). This
% is the "integrate" step -- features from every cell of the joint space, tabulated and
% indexed, so that state labelling, recurrence and law fitting become analyses OF THE TABLE
% rather than three separate pipelines over the recording.
%
% ⭐ THE COEFFICIENT COVARIANCE IS A SUFFICIENT STATISTIC. Every quantity here is a quadratic
% form, and a quadratic form commutes with summation over a window:
%
%     sum_{t in win} c* Q c  =  trace( Q * sum_{t in win} c c* )  =  trace( Q * Gwin )
%
% So Gwin [Ks x Ks] is formed once per (window, band) and every ROI x scale energy is a
% single inner product against it. MEASURED: ~0.7 s per page against ~7.6 s for the route
% that materialises vertex maps per frame -- about 50x, and exact rather than approximate.
%
% ⚠ THIS IS THE OPPOSITE VERDICT TO THE PER-FRAME CASE, and the difference is the window.
% Evaluating a form frame by frame, an ROI Gram costs Ks^2 PER ROI against V*Ks once, so the
% vertex route wins. Aggregating over a window, Gwin is formed ONCE and each ROI is then one
% inner product, so the Gram route wins by the length of the window. Same algebra, opposite
% answer, decided by where the summation sits.
%
% ⚠ WHAT IS NOT AVAILABLE THIS WAY. Only quadratic quantities. Peak amplitude within a
% window, envelope shape, phase-locking and anything non-quadratic need the per-frame series
% and are a separate, ~50x more expensive pass. Do not add them here.
%
% ⚠ NORMALISERS ARE CARRIED, NOT APPLIED. Raw energy is not comparable across subjects --
% head position, sensor coverage, individual anatomy -- but which normalisation is right
% depends on the question, and a table that has already divided cannot be undivided.
% TotalEnergy, ScaleEnergy, ROIArea and WindowSamples travel with the data.
%
% ⚠ WindowSamples DIFFERS BETWEEN BANDS. Multirate gives each band its own derived rate, so
% a one-second window holds a different number of samples per band. Comparing summed energy
% across bands without dividing by it compares sample counts, not physiology. `Density` does
% the division (per sample and per unit area).
%
% PROPERTIES
%   Bundle  Bands  BandLimits  WindowSec  Atlas  NumScales
%
% See also: rheome.flowpage, rheome.flowbrowser, docs/2026-08-24-flow-descriptors-design.md
%
% Author: Diellor Basha, 2026

    properties (SetAccess = immutable)
        Bundle
        Bands
        BandLimits
        WindowSec
        Atlas
        NumScales
        SamplesPerCycle
        Method
    end

    properties (Access = private)
        ROIGram_                                % [nROI x Ks^2] single -- vec(Phi' diag(w_r) Phi)
        Ks_
    end

    methods
        function obj = flowfeatures(B, varargin)
            p = inputParser;
            p.addParameter('Bands',     {'delta','theta','alpha','beta','gamma'});
            p.addParameter('WindowSec', 1);
            % ⚠ A SPEED KNOB WITH A REAL COST. The display rate (30/cycle) is set by what
            % it takes to SEE a rotation; a window ENERGY needs far fewer samples, and the
            % high bands dominate the cost purely through sample count -- gamma at 2400 Hz
            % carries 21x more samples than delta at 114 Hz. Lowering this for feature
            % extraction is defensible (a 1 s gamma window still holds hundreds of samples)
            % but it is a different rate from the one the browser shows, so the two are no
            % longer sample-for-sample comparable.
            p.addParameter('SamplesPerCycle', 30);
            % 'bandpass' by default: the table wants window ENERGIES, and the summed-CWT
            % route spends ~92% of the pass on transforms whose per-filter structure the
            % table then sums away. 'cwt' remains for cross-checking.
            p.addParameter('Method', 'bandpass');
            p.parse(varargin{:});
            o = p.Results;

            known = struct('delta',[1 4], 'theta',[4 8], 'alpha',[8 13], ...
                           'beta',[13 30], 'gamma',[30 60]);
            names = cellstr(o.Bands);
            lims  = cell(1, numel(names));
            for k = 1:numel(names)
                if ~isfield(known, names{k})
                    error('flowfeatures:band', ...
                        'Unknown band ''%s''. Known: %s.', names{k}, strjoin(fieldnames(known)', ', '));
                end
                lims{k} = known.(names{k});
            end

            obj.Bundle     = B;
            obj.Bands      = names(:).';
            obj.BandLimits = lims;
            obj.WindowSec  = o.WindowSec;
            obj.SamplesPerCycle = o.SamplesPerCycle;
            obj.Method     = o.Method;
            obj.Atlas      = B.atlas;
            obj.NumScales  = B.gfb.NumMembers;
            obj.Ks_        = numel(B.Lambda);
            obj.ROIGram_   = i_roigrams(B);
        end

        out = extract(obj, pageIndex)
        T   = totable(obj, out)
    end
end

% ---- vec(Phi' diag(w_r) Phi) per ROI, stacked so a whole page of traces is one GEMM ----
% Built once per anatomy: Ks^2 * V flops in total, and Ks^2 * nROI * 4 bytes resident
% (174 MB at 800 modes and 68 scouts, the same order as Phi itself).
function G = i_roigrams(B)
    Phi = double(B.Phi);  w = double(B.wVert(:));
    M   = B.atlas.Membership;  nR = size(M, 1);  Ks = size(Phi, 2);
    G   = zeros(nR, Ks*Ks, 'single');
    for r = 1:nR
        v  = find(M(r, :));
        Pr = Phi(v, :);
        Gr = Pr.' * (w(v) .* Pr);                 % [Ks x Ks], symmetric PSD
        G(r, :) = single(Gr(:).');
    end
end

% Author: Diellor Basha, 2026
