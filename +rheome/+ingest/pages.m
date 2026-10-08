function plan = pages(g, bands, cfg)
% INGEST.PAGES  The page plan: one grid level per band, bands grouped into jobs.
%
%   plan = rheome.ingest.pages(g, bands, cfg)
%
% A PAGE IS A TILE AT A CHOSEN LEVEL, and its halo is the band's time support (design
% §3.12). Each band gets the smallest level whose frame is at least
% 2*tSupport/PageOverhead long (so the halo costs at most PageOverhead of the core), but
% not below MinPageLength and not above Lmax. A band at Lmax runs on the whole record with
% no halo -- its edges are the record's edges, which the cone count already labels.
%
% Bands that share a level and are contiguous in scale form one JOB: one sub-bank per
% page, anchored on the job's top master scale, so every coefficient is the master's.
% A job whose per-channel working set (5x scales x longest span x bytes) would exceed MaxBytes is split
% band by band; a single band that still exceeds it errors with the number.
%
% INPUTS:
%   g      rheome.ingest.grid       bands  rheome.ingest.bank       cfg  rheome.ingest.config
% OUTPUT:
%   plan  table, one row per job: job, level, bands {j...}, scales [top bottom] master
%         indices, haloSamples, coreSamples, nPages, spanMax, bufferBytes (working set)
%
% See also: rheome.ingest.reducepaged, rheome.ingest.build
%
% Author: Diellor Basha, 2026

    arguments
        g     (1,1) struct
        bands table
        cfg   (1,1) struct
    end

    fs = g.fs;  nT = g.nT;
    Lmin = max(0, ceil(log2(cfg.MinPageLength / cfg.FrameFloor)));
    Lmin = min(Lmin, g.Lmax);
    nB = height(bands);
    level = zeros(nB, 1);  halo = zeros(nB, 1);
    for b = 1:nB
        need = 2 * bands.tSupport(b) / cfg.PageOverhead;
        L = find(g.tExtent >= need, 1) - 1;               % index -> level
        if isempty(L), L = g.Lmax; end
        L = min(max(L, Lmin), g.Lmax);
        level(b) = L;
        if L < g.Lmax, halo(b) = ceil(bands.tSupport(b) * fs); end
    end

    bytesPer = (8 + 8 * strcmp(cfg.Precision, 'double')) * 5;   % x5: wt's working set, as in rheome.ingest.build
    rows = {};
    b = 1;
    while b <= nB
        L = level(b);
        e = b;
        while e < nB && level(e+1) == L, e = e + 1; end   % contiguous run at this level
        grp = b:e;
        % greedy split of the run so each job fits MaxBytes
        s = grp(1);
        while s <= grp(end)
            t = s;
            while t < grp(end) && i_bytes(bands, s, t+1, L, halo, g, bytesPer) <= cfg.MaxBytes
                t = t + 1;
            end
            bb = i_bytes(bands, s, t, L, halo, g, bytesPer);
            if bb > cfg.MaxBytes
                error('ingest:pages:bytes', ...
                      'Band %d alone needs %.3g bytes per channel at level %d (span %d samples); MaxBytes is %.3g. Lower MinPageLength or raise MaxBytes.', ...
                      s, bb, L, i_span(s, t, L, halo, g), cfg.MaxBytes);
            end
            rows{end+1} = {numel(rows)+1, L, s:t, [bands.scales{s}(1), bands.scales{t}(end)], ...
                           max(halo(s:t)), min(g.F * 2^L, nT), g.K(L+1), i_span(s, t, L, halo, g), bb}; %#ok<AGROW>
            s = t + 1;
        end
        b = e + 1;
    end
    R = vertcat(rows{:});
    plan = table(cell2mat(R(:,1)), cell2mat(R(:,2)), R(:,3), cell2mat(R(:,4)), cell2mat(R(:,5)), ...
                 cell2mat(R(:,6)), cell2mat(R(:,7)), cell2mat(R(:,8)), cell2mat(R(:,9)), ...
                 'VariableNames', {'job','level','bands','scales','haloSamples','coreSamples', ...
                                   'nPages','spanMax','bufferBytes'});
end

function n = i_span(s, t, L, halo, g)
    n = min(g.F * 2^L, g.nT) + 2 * max(halo(s:t));
    n = min(n, g.nT);
end

function bb = i_bytes(bands, s, t, L, halo, g, bytesPer)
    nS = numel(bands.scales{s}(1):bands.scales{t}(end));   % master indices, top first
    bb = nS * i_span(s, t, L, halo, g) * bytesPer;
end
% Author: Diellor Basha, 2026
