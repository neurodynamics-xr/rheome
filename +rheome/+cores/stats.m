function [R, steps] = stats(C, A, vertArea, Vmm, f0, opts)
% CORES.STATS  Core tracks and track measures of a curl movie (nxr core_stats).
%
%   [R, steps] = rheome.cores.stats(C, A, vertArea, Vmm, f0)
%   [R, steps] = rheome.cores.stats(C, A, vertArea, Vmm, f0, Fs=12, Frames=1:size(C,2), Threshold=[], Link=15)
%
%   C [V x T] curl movie (single; rheome.cores.operator), A adjacency (rheome.cores.adjacency), vertArea [V x 1] mm^2,
%   Vmm [V x 3] mm, f0 the band's centre frequency (Hz).
%
% Threshold-on-curl core detector. Per frame, cores = rheome.cores.patches at rheome.cores.threshold (or Threshold=).
% Tracks: greedy per core in patch order; candidate tracks sorted by distance (stable, ties by track insertion order);
% the first of the same sense, within Link mm (15), not yet used. Unmatched cores start tracks; a frame gap closes all
% tracks. Tracks under 3 frames are dropped from life, path and disp but stay in steps; track ids count every track.
%
% R fields: threshold, cores_per_frame, area_mm2, n_tracks, life_ms, life_cycles, speed_mm_s, net_speed_mm_s, disp_mm,
%   straightness, reversal, persist_half, recur_cycle, jaccard_cycle, jaccard_half, frames, life_hist, speed_hist
%   (20 ms and 20 mm/s bins on 0:20:2000), dens_ccw, dens_cw (single, [V x 1]: fraction of frames above +thr / below -thr),
%   plus patches (struct: frame, sense, centre, area; frames are columns of C, 1-based).
%   speed = path/((n-1)/Fs); medians for area (over all patches), life, speed, disp, straightness; means for reversal,
%   persist_half, recur_cycle, jaccard. H = max(round_half_even(Fs/(2 f0)),1), P = max(round_half_even(Fs/f0),1).
%   A core persists if any core H frames earlier lies within Link mm with the same sense; reverses if with the opposite sense.
%   Jaccard: area-weighted over vertices with a nonzero mask, over frames with a frame lag earlier, skipping empty unions.
% steps [S x 9]: sense, from xyz, displacement xyz, track id, frame (1-based).
%
% Ports nxr-analysis src/nxr/vortex.py core_stats (commit 903fdc2). Not nxr's critical-point route; see rheome.detect.
%
% See also: rheome.cores.patches, rheome.cores.threshold, rheome.cores.operator, rheome.detect.track
%
% Author: Diellor Basha, 2026
    arguments
        C
        A
        vertArea
        Vmm
        f0 (1,1) double
        opts.Fs (1,1) double = 12
        opts.Frames double = []
        opts.Threshold double = []
        opts.Link (1,1) double = 15
        opts.Chunk (1,1) double = 3000
    end
    C = single(C); vertArea = double(vertArea(:)); fs = opts.Fs; link = opts.Link;
    if isempty(opts.Frames), fr = 1:size(C, 2); else, fr = opts.Frames(:).'; end
    C = C(:, fr); nf = numel(fr); nk = size(C, 1);
    if isempty(opts.Threshold), thr = rheome.cores.threshold(C, Chunk=opts.Chunk); else, thr = opts.Threshold; end
    H = max(halfeven(fs / (2*f0)), 1); P = max(halfeven(fs / f0), 1);
    pos = C >= single(thr); neg = C <= -single(thr);
    mask = int8(pos) - int8(neg);
    [cf, cs, cc, ca] = rheome.cores.patches(C, thr, A, vertArea, Vmm);
    starts = [0; cumsum(accumarray(cf, 1, [nf 1]))];       % core q of frame j: starts(j)+1 : starts(j+1)
    ncores = diff(starts);
    life = []; path = []; disp_ = []; steps = zeros(0, 9);
    active = struct('id', {}, 'sg', {}, 'c', {}, 'c0', {}, 'n', {}, 'path', {});
    tid = 0; prev = [];
    for j = 1:nf
        k = fr(j); new = active([]);
        if ~(~isempty(prev) && k == prev + 1 && ~isempty(active))
            [life, path, disp_] = closeAll(active, life, path, disp_); active = active([]);
        end
        cen = reshape([active.c], 3, []).'; sgn = [active.sg]; usedIds = [];
        for q = starts(j)+1 : starts(j+1)
            sg = cs(q); ce = cc(q, :); match = 0;
            if ~isempty(active)
                d = vecnorm(cen - ce, 2, 2).'; d(sgn ~= sg) = inf;
                [ds, ord] = sort(d);                         % stable: ties by insertion order
                for u = 1:numel(ord)
                    if ds(u) > link, break, end
                    if ~any(usedIds == active(ord(u)).id), match = ord(u); break, end
                end
            end
            if match == 0
                tid = tid + 1;
                new(end+1) = struct('id', tid, 'sg', sg, 'c', ce, 'c0', ce, 'n', 1, 'path', 0); %#ok<AGROW>
            else
                t = active(match); usedIds(end+1) = t.id; %#ok<AGROW>
                steps(end+1, :) = [sg, t.c, ce - t.c, t.id, k]; %#ok<AGROW>
                t.path = t.path + norm(ce - t.c); t.c = ce; t.n = t.n + 1; new(end+1) = t; %#ok<AGROW>
            end
        end
        gone = ~ismember([active.id], usedIds);
        [life, path, disp_] = closeAll(active(gone), life, path, disp_);
        active = new; prev = k;
    end
    [life, path, disp_] = closeAll(active, life, path, disp_);

    posOf = zeros(1, max(fr) + 1);                           % frame value -> position (0 = absent)
    posOf(fr + 1) = 1:nf; rec = struct();
    lags = [H P]; nm = {'half', 'cycle'};
    for li = 1:2
        lag = lags(li);
        jj = find(fr >= lag); jj = jj(posOf(fr(jj) - lag + 1) > 0); oo = posOf(fr(jj) - lag + 1);
        Mo = mask(:, oo); Mn = mask(:, jj);
        uni = vertArea.' * double((Mo ~= 0) | (Mn ~= 0)); itr = vertArea.' * double((Mo == Mn) & (Mn ~= 0)); keep = uni > 0;
        same = false(0, 1); opp = false(0, 1);
        for m = 1:numel(oo)
            a_ = oo(m); b_ = jj(m);
            if starts(a_+1) == starts(a_), continue, end
            qa = starts(a_)+1 : starts(a_+1); qb = starts(b_)+1 : starts(b_+1);
            if isempty(qb), near = false(numel(qa), 0); sm = near;
            else
                near = pdist2s(cc(qa, :), cc(qb, :)) <= link; sm = cs(qa) == cs(qb).';
            end
            same = [same; any(near & sm, 2)]; opp = [opp; any(near & ~sm, 2)]; %#ok<AGROW>
        end
        rec.(nm{li}) = {same, opp, itr(keep) ./ uni(keep)};
    end
    dur = (life - 1) / fs; sp = path ./ max(dur, 1e-9); nsp = disp_ ./ max(dur, 1e-9);
    md = @(v) iff(isempty(v), NaN, @() median(v)); mn = @(v) iff(isempty(v), NaN, @() mean(v));
    R = struct('threshold', thr, 'cores_per_frame', mean(ncores), 'area_mm2', md(ca), 'n_tracks', numel(life), ...
        'life_ms', md(life) / fs * 1000, 'life_cycles', md(life) / fs * f0, 'speed_mm_s', md(sp), 'net_speed_mm_s', md(nsp), ...
        'disp_mm', md(disp_), 'straightness', md(disp_ ./ max(path, 1e-9)), ...
        'reversal', mn(rec.half{2}), 'persist_half', mn(rec.half{1}), 'recur_cycle', mn(rec.cycle{1}), ...
        'jaccard_cycle', mn(rec.cycle{3}), 'jaccard_half', mn(rec.half{3}), 'frames', nf);
    R.life_hist = histcounts(life / fs * 1000, 0:20:2000).';
    R.speed_hist = histcounts(sp, 0:20:2000).';
    R.dens_ccw = single(sum(pos, 2) / max(nf, 1)); R.dens_cw = single(sum(neg, 2) / max(nf, 1));
    R.patches = struct('frame', cf, 'sense', cs, 'centre', cc, 'area', ca);
    R.life = life; R.speed = sp; R.nk = nk;
end

function [life, path, disp_] = closeAll(tr, life, path, disp_)
    for i = 1:numel(tr)
        t = tr(i);
        if t.n >= 3, life(end+1) = t.n; path(end+1) = t.path; disp_(end+1) = norm(t.c - t.c0); end %#ok<AGROW>
    end
end

function r = halfeven(x)
% Python round(): halfway cases to the even integer
    r = round(x);
    if abs(x - fix(x)) == 0.5, r = 2 * round(x / 2); end
end

function d = pdist2s(a, b)
    d = sqrt((a(:,1) - b(:,1).').^2 + (a(:,2) - b(:,2).').^2 + (a(:,3) - b(:,3).').^2);
end

function v = iff(c, a, b)
    if c, if isa(a, 'function_handle'), v = a(); else, v = a; end
    else, if isa(b, 'function_handle'), v = b(); else, v = b; end
    end
end
