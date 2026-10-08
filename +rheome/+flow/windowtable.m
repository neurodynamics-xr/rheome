function [T, X] = windowtable(name, varargin)
% FLOW.WINDOWTABLE  One table indexed by time window: every measurement under a tile on each axis.
%
%   T        = rheome.flow.windowtable('sub01')
%   [T, X]   = rheome.flow.windowtable(name, WindowSec=2, Bands=[4 8; 8 16], MaxDepth=4)
%
% One row per (time window, band, cortical tile). Columns carry the activation statistics under
% that tile, the sensor band power for that window, the eigenmode-octave shares, the gauge split
% and the flow scalars. `X` holds what does not fit a rectangular table: the full per-group Dirac
% power [nGroup x nWindow x nBand] and the axis definitions.
%
% ⭐ ROLL-UP IS DELIBERATELY NOT REQUIRED HERE. An earlier design restricted this table to
% quadratic quantities so every coarser tile would be an exact sum of finer ones. That buys speed
% and exactness and it excludes peak amplitude, crest, vortex counts and anything else that is not
% a quadratic form. This table includes them: it is a LABELLING table, so the cost of recomputing a
% coarser window is accepted. ⚠ The consequence is real and must not be forgotten: a row at depth 2
% is NOT the sum of its depth-3 children for the non-quadratic columns, and `energy` is the only
% column for which summing children is valid.
%
% ⚠⚠ WHAT THIS TABLE CANNOT LABEL, MEASURED BY PLANTING (plant_scale_omega.m). A source of known
% scale planted through the leadfield, noise and inverse comes back at ~100 mm REGARDLESS of its
% true size once the SNR is at or below 10 dB: slope +0.02, R^2 0.002 across a planted 70-267 mm
% range. So `kCentroidMM` and the octave shares describe WHICH FILTER responded, not how big the
% source was, and they are named and documented as such. ⭐ Location survived the same test
% (43-52 mm median error, correct depth-3 node 55-65%), so `node` IS a measurement.
%
% ⚠ ADMISSIBILITY IS A COLUMN, NOT A FILTER. Tiles below 2*r50 = 104 mm are kept with
% `admissible=false` rather than dropped, because they are the control that shows when a result is
% instrumental -- the 70 mm rung has repeatedly behaved like the admissible ones.
%
% ⚠ DEPTH BEATS EVERYTHING on localisation: at 10 dB the planted-source error runs 13/35/41/97 mm
% by distance-to-nearest-coil quartile and the correct-node rate 98/64/28/28%. `depthMM` is carried
% per tile for exactly this reason; a comparison between two tiles at different depths is a depth
% comparison.
%
% ⚠ ONE WINDOW GRID FOR ALL BANDS, which is NOT the constant-Q convention. `rheome.select.ladder` gives
% each band its own natural level so that cycles per tile is constant (22.6 on this bank). A table
% indexed by a single time window cannot do that, so `nCycles` is carried per row and a band whose
% `nCycles` is small is under-resolved in that window. At WindowSec=2 the 4-8 Hz row holds 8-16
% cycles and delta would hold 4 or fewer -- which is why delta is not in the default band list.
%
% ⭐ WINDOW LENGTH: 1-2 s is the measured sweet spot. The activation decorrelates in ~2 s and
% adjacent 2 s windows are already 90% decorrelated, so rows at that length are near-independent.
% ⚠ Samples WITHIN a window are not: Bartlett gives n_eff = 25 for 599 frames at lag-1 0.99, so
% `nEff` and not `nSamples` is the divisor for any error bar. Both are carried.
%
% ⚠ NORMALISERS ARE CARRIED, NOT APPLIED. A table that has divided cannot be undivided.
%
% NAME-VALUE
%   Tiles ([])      [tStart tEnd] in seconds, one row per window -- pass the STORE's tile spans at
%                   a level (rheome.select.derive returns t_lo/t_hi) so every row joins on (level, k)
%   BandIds ([])    the store's own band ids, carried into a band_id column
%   WindowSec (2)   used only when Tiles is empty
%   Bands ([4 8; 8 16; 16 32; 32 64])   MaxDepth (4)   Hemi ("L")
%   Flow (true)     div/curl/vortex scalars on the window-mean current (the costly extras)
%   MaxWindows (Inf)   stop after this many windows per band
%   Verbose (true)
%
% COLUMNS. ⭐ A `win_` prefix means the column is a property of the WINDOW and is identical across
% every tile in it; everything without the prefix varies by tile. That prefix is load-bearing: an
% earlier version computed the gauge split and the vortex count hemisphere-wide and copied them into
% every row, where they were indistinguishable from measured per-tile values.
%   index : window tStart tEnd band fLo fHi nCycles
%   tile  : node depth nodeDiameterMM nodeAreaM2 depthMM admissible
%   per tile : meanAct peakAct p95Act stdAct crest energy density share normalShare tangShare
%              divRMS curlRMS nVortex
%   per window : win_nSamples win_nEff win_bandPowerSensor win_envMean win_oct1..6
%                win_kCentroidMM (⚠ which filter responded, NOT a source size)
%
% See also: rheome.flowfeatures, rheome.geom.tree, rheome.operators.gauge, plant_scale_omega,
%           docs/2026-09-26-feature-table-design.md
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('WindowSec', 2,  @(x) isscalar(x) && x > 0);
    p.addParameter('Bands', [4 8; 8 16; 16 32; 32 64], @(x) size(x,2)==2);
    p.addParameter('MaxDepth', 4, @isscalar);
    p.addParameter('Hemi', "L");
    p.addParameter('Flow', true, @islogical);
    p.addParameter('Tiles', [], @(x) isempty(x) || size(x,2)==2);   % [tStart tEnd] seconds
    p.addParameter('BandIds', [], @isnumeric);        % the store's own band ids, carried through
    p.addParameter('MaxWindows', Inf, @isscalar);      % for tests and quick looks
    p.addParameter('Verbose', true, @islogical);
    p.parse(varargin{:});
    o = p.Results;  vb = o.Verbose;

    %% geometry, gauge and the tile axis -- once
    B = rheome.load.bases(name);
    H = B.(char(o.Hemi));
    S = H.S;  gv = double(H.gv(:))';  lbo = H.lbo;  nV = size(S.Vertices,1);
    av = full(sum(lbo.Mass,2));
    T0 = rheome.geom.tree(S, L=lbo.L, M=lbo.M, MaxDepth=o.MaxDepth);
    keepN = T0.depth >= 1;                                  % the root is the whole hemisphere
    nodes = T0(keepN,:);  nN = height(nodes);
    C = rheome.operators.connection_laplacian(S.Vertices, double(S.Faces));
    g = rheome.operators.gauge(S.Vertices, double(S.Faces), Connection=C);   % smooth frame, one solve
    if vb, fprintf('%d tiles at depths %d-%d, %.0f-%.0f mm\n', nN, min(nodes.depth), ...
            max(nodes.depth), 1e3*max(nodes.diameter), 1e3*min(nodes.diameter)); end

    %% the instrument
    st = rheome.load.study(name);
    isMEG = strcmpi(st.chan.Type,'MEG');
    iSel = find(isMEG(:) & (st.rec.ChannelFlag(:)==1));
    G = double(st.hm.Gain(iSel,:));  NC = double(st.ncov.NoiseCov(iSel,iSel));
    chT = st.chan.Type(iSel);  fs = st.rec.sfreq;  F = double(st.rec.F(iSel,:));
    Loc = cell2mat(arrayfun(@(c) c.Loc(:,1), st.chan.Channel(iSel), 'UniformOutput', false));
    clear st
    depthMM = 1e3*min(pdist2(S.Vertices, Loc'), [], 2);
    d = rheome.load.dirac(name);  Phi = d.Phi;  nM = d.nModes;
    R = rheome.inverse.dirac(rheome.forward.dirac(G,d), d, struct('NoiseCov',NC), ...
            struct('ChannelTypes',{chT},'InverseMeasure','amplitude'));
    K = R.ImagingKernelMode;
    rx = (gv-1)*4+2;  ry = (gv-1)*4+3;  rz = (gv-1)*4+4;
    PhiC = Phi([rx ry rz], :);                              % the current rows, once

    %% the Dirac length axis and the Weyl octaves
    kb = nan(nM,1);  shr = zeros(nM,1);
    for m = 1:nM
        e = sum(abs(lbo.Phi'*(lbo.Mass*[Phi(rx,m) Phi(ry,m) Phi(rz,m)])).^2, 2);
        if sum(e) > 0, shr(m) = sum(e);  kb(m) = sqrt(sum(e.*lbo.Lambda(:))/sum(e)); end
    end
    wlM = 1e3*2*pi./max(kb,eps);
    okM = isfinite(kb) & kb>0 & shr>0.1*max(shr) & wlM<400;  % ⚠ kb>0 alone lets 6e9 mm modes through
    % ⚠⚠ the eigenmode unit is the DEGENERATE GROUP of 4, never the mode: the four coefficients in a
    % group are a gauge. Blocks of 4 in sorted order, with the spread asserted -- NOT a tolerance,
    % which chains distinct eigenvalues into groups of [4 36 40 44 48 444].
    GS = 4;  assert(mod(nM,GS)==0, 'nModes must be a multiple of 4');
    grp = repelem((1:nM/GS)', GS);  nG = nM/GS;
    spread = accumarray(grp, d.Lambda(:), [nG 1], @(v)(max(v)-min(v))/max(mean(abs(v)),realmin));
    if max(spread) > 1e-4, warning('flow:windowtable:degeneracy', ...
            'blocks of 4 do not match the degeneracy (max spread %.1e)', max(spread)); end
    kg  = accumarray(grp, kb.*shr, [nG 1]) ./ max(accumarray(grp, shr, [nG 1]), realmin);
    wlG = 1e3*2*pi./max(kg,eps);
    okG = accumarray(grp, okM, [nG 1]) > 0 & wlG < 400;
    % ⚠⚠ SPATIAL OCTAVES ARE BINNED ON WAVELENGTH, NOT ON THE MODE INDEX. An earlier version used
    % Weyl's x4-modes-per-octave to put the edges at cumulative mode counts [1 4 16 66 249 663].
    % That is the right COUNT but the wrong AXIS: it defines octaves by how many modes have been
    % passed rather than by the wavelength band they occupy, and rheome.geom.ladder bins on wavenumber, so
    % the two disagreed about which octave a group belonged to. One codebase, two definitions of
    % "spatial octave". Binned on wavelength both agree, and Weyl's x4 then comes out as a RESULT
    % (measured 0.95-0.98 of prediction) instead of being assumed.
    A0 = sum(av);
    wlMax = 1e3*2*sqrt(A0/pi);                          % the hemisphere as an equivalent disc, mm
    nOct = 6;
    edgesW = wlMax * 2.^(0:-1:-nOct);                   % descending wavelength edges
    oct = discretize(wlG, flip(edgesW));                % 1 = coarsest
    oct = nOct + 1 - oct;  oct(~isfinite(oct)) = nOct;
    oct = min(nOct, max(1, oct));

    %% flow operators, if asked
    op = [];
    if o.Flow
        % ⚠ rheome.detect.operator -> rheome.utils.hemisphere reads .Hemi{h}, .HemiLabel{h}, .Comment, .nV,
        % .VertNormals and .SurfaceFile. rheome.load.bases's per-hemisphere S carries only some of them, so
        % a single-hemisphere surface has to be re-dressed as a one-hemisphere cortex.
        Sg = S;  Sg.nV = nV;  Sg.Hemi = {(1:nV)'};  Sg.HemiLabel = {sprintf('Cortex %s', o.Hemi)};
        if ~isfield(Sg,'Comment')     || isempty(Sg.Comment),     Sg.Comment = char(name); end
        if ~isfield(Sg,'SurfaceFile'), Sg.SurfaceFile = ''; end
        if ~isfield(Sg,'VertNormals'), Sg.VertNormals = C.normal; end
        op = rheome.detect.operator(Sg);
    end

    %% per band, per window
    nB = size(o.Bands,1);  rows = {};  Xg = [];
    for bi = 1:nB
        lo = o.Bands(bi,1);  hi = o.Bands(bi,2);  f0 = sqrt(lo*hi);
        dec = max(1, floor(fs/max(4*hi, 60)));  fsd = fs/dec;
        bp = designfilt('bandpassiir','FilterOrder',8,'HalfPowerFrequency1',lo, ...
                        'HalfPowerFrequency2',hi,'SampleRate',fs);
        Fa = hilbert(filtfilt(bp, F.')).';
        Fa = Fa(:, 1:dec:end);
        % ⭐ WINDOWS COME FROM THE STORE'S TILE GRID WHEN ONE IS GIVEN, so a sweep at level L uses
        % exactly the tiles the database already indexes and every row joins on (level, k) with no
        % interpolation. Falling back to a uniform WindowSec is for standalone use.
        if isempty(o.Tiles)
            nPW0 = round(o.WindowSec*fsd);  nWall = floor(size(Fa,2)/nPW0);
            edges = [(0:nWall-1)'*o.WindowSec, (1:nWall)'*o.WindowSec];
        else
            edges = o.Tiles;
        end
        nW = min(size(edges,1), o.MaxWindows);
        if vb, fprintf('band %2g-%2g Hz: %g Hz, %d windows, %.2f-%.2f s (%.1f cycles)\n', ...
                lo, hi, fsd, nW, min(diff(edges,1,2)), max(diff(edges,1,2)), ...
                mean(diff(edges,1,2))*f0); end
        Pg = zeros(nG, nW);
        for w = 1:nW
            i0 = max(1, floor(edges(w,1)*fsd) + 1);
            i1 = min(size(Fa,2), max(i0+1, round(edges(w,2)*fsd)));
            ix = i0:i1;  nPW = numel(ix);
            Fw = Fa(:,ix);
            c  = K*Fw;                                       % [nM x nPW] Dirac coefficients
            pm = mean(abs(c).^2, 2);
            pgw = accumarray(grp, pm, [nG 1]);  Pg(:,w) = pgw;
            J  = PhiC*c;                                     % [3nV x nPW] band-limited current
            A  = sqrt(sum(reshape(abs(J), nV, 3, nPW).^2, 2));
            A  = reshape(A, nV, nPW);                        % |J| per vertex per sample
            aMean = mean(A,2);  aPeak = max(A,[],2);
            % gauge split and flow scalars on the window-MEAN complex current.
            % ⚠ PER-VERTEX, so they can be restricted to a tile below. An earlier version summed
            % these over the whole hemisphere once per window and copied the result into every tile's
            % row: the columns then looked per-tile and were identical down the window, which is the
            % kind of thing a reader trusts silently. Anything that CAN be localised, is.
            Jm = reshape(mean(reshape(J, nV, 3, nPW), 3), nV, 3);
            qn = av.*abs(sum(Jm.*g.normal,2)).^2;
            q1 = av.*abs(sum(Jm.*g.e1,2)).^2;
            q2 = av.*abs(sum(Jm.*g.e2,2)).^2;
            dv = []; cv = []; cpNode = zeros(nV,1);
            if o.Flow
                Jr = reshape(real(Jm)', [], 1);
                dv = rheome.differential.divergence(Jr, S, op.fg);
                cv = rheome.differential.curl(Jr, S, op.fg);
                cp = rheome.detect.criticalPoints(Jr, op.Sh{1}, 'all', op);
                if ~isempty(cp.pos)                          % attribute each one to a vertex
                    [~, nvv] = min(pdist2(cp.pos, S.Vertices), [], 2);
                    cpNode = accumarray(nvv, 1, [nV 1]);
                end
            end
            bpow = mean(abs(Fw(:)).^2);  emean = mean(sqrt(mean(abs(Fw).^2,1)));
            nEff = i_neff(sqrt(mean(abs(Fw).^2,1)));
            oshare = accumarray(oct(okG), pgw(okG), [nOct 1]) / max(sum(pgw(okG)), realmin);
            kc = sum(wlG(okG).*pgw(okG)) / max(sum(pgw(okG)), realmin);
            for k = 1:nN
                v = nodes.members{k};  w_a = av(v);  sa = sum(w_a);
                mA = sum(w_a.*aMean(v))/sa;  pA = max(aPeak(v));
                en = sum(w_a.*mean(A(v,:).^2, 2));
                t3 = sum(qn(v)) + sum(q1(v)) + sum(q2(v));   % ⭐ the gauge split UNDER THIS TILE
                rowk = [w, edges(w,1), edges(w,2), bi, lo, hi, (edges(w,2)-edges(w,1))*f0, ...
                        nPW, nEff, nodes.node_id(k), nodes.depth(k), 1e3*nodes.diameter(k), ...
                        nodes.area(k), median(depthMM(v)), 1e3*nodes.diameter(k) > 104, ...
                        mA, pA, prctile(aMean(v),95), std(aMean(v)), pA/max(mA,realmin), ...
                        en, en/max(nEff,1)/max(nodes.area(k),realmin), NaN, ...
                        sum(qn(v))/max(t3,realmin), (sum(q1(v))+sum(q2(v)))/max(t3,realmin), ...
                        i_rms(dv,v,w_a), i_rms(cv,v,w_a), sum(cpNode(v)), ...
                        bpow, emean, oshare(:)', kc];
                rows{end+1} = rowk; %#ok<AGROW>
            end
        end
        Xg = cat(3, Xg, Pg);
    end

    M = vertcat(rows{:});
    % ⭐ the `win_` prefix marks a column that is a property of the WINDOW, identical across every
    % tile in it. Everything without the prefix varies by tile. This is naming as documentation:
    % without it the replicated columns are indistinguishable from measured ones.
    vn = [{'window','tStart','tEnd','bandIdx','fLo','fHi','nCycles','win_nSamples','win_nEff', ...
           'node','depth','nodeDiameterMM','nodeAreaM2','depthMM','admissible', ...
           'meanAct','peakAct','p95Act','stdAct','crest','energy','density','share', ...
           'normalShare','tangShare','divRMS','curlRMS','nVortex', ...
           'win_bandPowerSensor','win_envMean'}, ...
          arrayfun(@(k) sprintf('win_oct%d',k), 1:nOct, 'UniformOutput', false), ...
          {'win_kCentroidMM'}];
    T = array2table(M, 'VariableNames', vn);
    T.admissible = logical(T.admissible);
    T.band = categorical(arrayfun(@(k) sprintf('%g-%g Hz', o.Bands(k,1), o.Bands(k,2)), ...
                                  T.bandIdx, 'UniformOutput', false));
    % `share` is the only column a parent may take as a sum of children; fill it now
    for bi = 1:nB
        for w = unique(T.window(T.bandIdx==bi))'
            m = T.bandIdx==bi & T.window==w;
            tot = sum(T.energy(m & T.depth==max(T.depth(m))));
            T.share(m) = T.energy(m)/max(tot,realmin);
        end
    end
    % ⭐ energy is the one column for which a parent IS the sum of its children -- assert it, since
    % the whole point of admitting non-additive columns is that only this one keeps the property.
    rollErr = 0;
    for bi = 1:nB
        w1 = T.bandIdx==bi & T.window==1;
        if any(w1)
            e1 = sum(T.energy(w1 & T.depth==1));  e4 = sum(T.energy(w1 & T.depth==max(T.depth(w1))));
            rollErr = max(rollErr, abs(e1-e4)/max(e1,realmin));
        end
    end
    if vb, fprintf('energy roll-up (depth 1 total vs deepest total): max relative error %.2e\n', rollErr); end
    if rollErr > 1e-9, warning('flow:windowtable:rollup', ...
            'energy does not roll up (%.1e); the tree partition may not be disjoint', rollErr); end

    if ~isempty(o.BandIds), T.band_id = reshape(o.BandIds(T.bandIdx), [], 1);
    else,                    T.band_id = zeros(height(T), 1); end
    X = struct('GroupPower', Xg, 'wavelengthMM', wlG, 'octave', oct, 'octaveEdgesMM', edgesW, ...
               'okGroup', okG, ...
               'Bands', o.Bands, 'BandIds', o.BandIds, 'WindowSec', o.WindowSec, 'nodes', nodes, ...
               'gaugeMethod', g.method, 'gaugeSingular', g.singular);
    if vb, fprintf('table: %d rows x %d columns\n', height(T), width(T)); end
end

function r = i_rms(x, v, w)
    if isempty(x), r = NaN; else, r = sqrt(sum(w.*x(v).^2)/max(sum(w),realmin)); end
end

function n = i_neff(e)
% Bartlett: the autocorrelation of the envelope, not the sample count. ⚠ 599 frames at lag-1 0.99
% give n_eff = 25, so using nSamples as a divisor overstates precision by ~5x.
    e = e(:) - mean(e);  n = numel(e);
    if n < 4 || all(e==0), n = numel(e); return; end
    L = min(n-1, floor(n/4));
    ac = arrayfun(@(k) (e(1:n-k)'*e(1+k:n))/max(e'*e,realmin), 1:L);
    n = n / max(1 + 2*sum((1-(1:L)/n).*ac), 1);
end

% Author: Diellor Basha, 2026
