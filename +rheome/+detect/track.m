function out = track(J, S, opts)
% DETECT.TRACK  Track flow critical points across time into trajectories (Lagrangian view).
%
%   out = rheome.detect.track(J, S, opts)
%
% Per frame, rheome.detect.criticalPoints finds the classified critical points (vortex/source/sink/
% saddle) of the ambient current; this links them across frames into trajectories and measures
% each trajectory's velocity, lifetime, straightness and intrinsic home location
% (rheome.geom.karcher_mean). The Lagrangian complement to the Eulerian velocity field: it follows the
% FEATURES of the flow, not the field.
%
% ASSOCIATION — ByteTrack-style two-stage lifecycle. Ported from the sibling repo
%   nxr-tracking (tracking/meshByteTrack.ts, tracking/meshAssociation.ts; N. colleague, 2026),
% which adapts ByteTrack [Zhang et al., "ByteTrack: Multi-Object Tracking by Associating Every
% Detection Box", ECCV 2022] to a cortical mesh. The idea: do not throw weak detections away —
% match STRONG detections first, then RECOVER still-open tracks from the WEAK ones, and keep a
% missed track "lost" for a few frames before ending it. This turns the flickering, short-lived
% detections of spontaneous data into fewer, longer, coherent trajectories, and it flags
% merge/split events (for flow singularities, pair annihilation/bifurcation). Specialised here to
% point features: the mesh-region cost terms of the source (support Jaccard, footprint size) are
% dropped; the cost is the on-surface distance, with same critical-point TYPE as a hard gate.
%
% INPUTS:
%   J    [3nV x nWin] ambient current over a window
%   S    surface struct (.Vertices, .Faces, .Hemi)
%   opts .types {'vortex','source','sink','saddle'}
%        .minStrength (activation; detections below are ignored; default 0.5*median frame-max)
%        .highStrength (HIGH band; default 1.5*minStrength) | .maxStep (max m/frame; default 5*edge)
%        .metric 'euclidean'(default)|'geodesic' | .dt (s; default 1)
%        .lostTrackBuffer (frames a confirmed track survives while missed; default 2)
%        .minimumConsecutiveFrames (tentative->confirmed; default 2)
%        .minLifetimeFrames (drop tracks with fewer samples; default 3)
%
% OUTPUT (struct out):
%   .tracks  struct array (one per trajectory), fields: type frames pos strength velocity speed
%            lifetime pathLength netDisplacement straightness birthFrame deathFrame state
%            missedFrames lowRecoveryCount splitCount mergeCount centerVertex centerPos.
%            Sorted by lifetime*mean(strength) descending.
%   .summary counts per type + per state + lifetime/speed vectors + nSplit/nMerge/nRecovered
%
% See also: rheome.detect.criticalPoints, rheome.geom.geodesic, rheome.geom.karcher_mean
%
% Author: Diellor Basha, 2026

    if nargin < 3, opts = struct; end
    if ~isfield(opts,'types')  || isempty(opts.types),  opts.types = {'vortex','source','sink','saddle'}; end
    if ~isfield(opts,'metric') || isempty(opts.metric), opts.metric = 'euclidean'; end
    if ~isfield(opts,'dt')     || isempty(opts.dt),     opts.dt = 1; end
    if ~isfield(opts,'lostTrackBuffer')          || isempty(opts.lostTrackBuffer),          opts.lostTrackBuffer = 2; end
    if ~isfield(opts,'minimumConsecutiveFrames') || isempty(opts.minimumConsecutiveFrames), opts.minimumConsecutiveFrames = 2; end
    if ~isfield(opts,'minLifetimeFrames')        || isempty(opts.minLifetimeFrames),        opts.minLifetimeFrames = 3; end
    meanEdgeLength = i_mean_edge(S);
    if ~isfield(opts,'maxStep') || isempty(opts.maxStep), opts.maxStep = 5*meanEdgeLength; end

    geodesicSolver = [];

    % ---- 1. per-frame detections. Pre-computed via opts.frameFeatures (for STREAMED detection over
    %         long recordings -- detect per block, discard J, associate here); else detect from J.
    %         ByteTrack matches every detection, so weak ones are kept. ----
    if isfield(opts,'frameFeatures') && ~isempty(opts.frameFeatures)
        frameFeatures = opts.frameFeatures(:);
        nWin = numel(frameFeatures);
    else
        nWin = size(J, 2);
        operator = rheome.detect.operator(S);
        frameFeatures = cell(nWin, 1);
        for t = 1:nWin
            cp = rheome.detect.criticalPoints(J(:,t), S, opts.types, operator);
            if isempty(cp.strength), frameFeatures{t} = struct('pos',[],'type',{{}},'strength',[]); continue; end
            frameFeatures{t} = struct('pos',cp.pos,'type',{cp.type},'strength',cp.strength);
        end
    end
    frameMaxStrength = cellfun(@i_frame_max, frameFeatures);
    if ~isfield(opts,'minStrength') || isempty(opts.minStrength)
        opts.minStrength = 0.5*median(frameMaxStrength(frameMaxStrength>0));      % activation
    end
    if ~isfield(opts,'highStrength') || isempty(opts.highStrength)
        opts.highStrength = 1.5*opts.minStrength;                                 % HIGH band
    end

    % ---- 2. ByteTrack-style two-stage association with a lost buffer ----
    % track fields: type frames pos strength stage state birthFrame lastFrame consecutive
    %               missed lowRecovery splitCount mergeCount
    tracks = i_empty_tracks();
    nSplit = 0;  nMerge = 0;  nRecovered = 0;
    for t = 1:nWin
        f = frameFeatures{t};
        strength = f.strength;
        isHigh = strength >= opts.highStrength;
        isLow  = strength >= opts.minStrength & ~isHigh;          % below activation -> ignored
        highIdx = find(isHigh);  lowIdx = find(isLow);

        openMask = strcmp({tracks.state}, 'tentative') | strcmp({tracks.state}, 'confirmed') | strcmp({tracks.state}, 'lost');
        openTracks = find(openMask);
        matchedTrack = false(1, numel(tracks));
        matchedDet   = false(1, numel(f.strength));

        % --- Stage 1: open tracks x HIGH detections ---
        [pairT, pairD, sC, mC] = i_associate(tracks, openTracks, f, highIdx, S, opts, geodesicSolver);
        nSplit = nSplit + sC;  nMerge = nMerge + mC;
        for p = 1:numel(pairT)
            tracks = i_apply_match(tracks, pairT(p), f, pairD(p), t, 'high', opts);
            matchedTrack(pairT(p)) = true;  matchedDet(pairD(p)) = true;
        end

        % --- Stage 2: still-open tracks x remaining LOW detections (recovery) ---
        openLeft = openTracks(~matchedTrack(openTracks));
        lowLeft  = lowIdx(~matchedDet(lowIdx));
        [pairT2, pairD2, sC2, mC2] = i_associate(tracks, openLeft, f, lowLeft, S, opts, geodesicSolver);
        nSplit = nSplit + sC2;  nMerge = nMerge + mC2;
        for p = 1:numel(pairT2)
            tracks = i_apply_match(tracks, pairT2(p), f, pairD2(p), t, 'low-recovery', opts);
            matchedTrack(pairT2(p)) = true;  matchedDet(pairD2(p)) = true;
            nRecovered = nRecovered + 1;
        end

        % --- Births: unmatched HIGH detections start tentative tracks ---
        for d = highIdx(~matchedDet(highIdx))'
            tracks(end+1) = i_new_track(f, d, t); %#ok<AGROW>
        end

        % --- Misses: open tracks unmatched this frame (tentative->end, else lost/end) ---
        for k = openTracks(~matchedTrack(openTracks))
            tracks(k) = i_register_miss(tracks(k), t, opts);
        end
    end

    % ---- 3. descriptors (keep tracks with >= minLifetimeFrames samples) ----
    keep = arrayfun(@(tr) numel(tr.frames) >= opts.minLifetimeFrames, tracks);
    tracks = tracks(keep);
    if ~isempty(tracks)
        [~, geodesicSolver] = rheome.geom.geodesic(S, 1, geodesicSolver);   % factor ONCE for all karcher means
    end
    descriptor = i_empty_descriptor();
    for k = 1:numel(tracks)
        tr = tracks(k);
        frameGaps = diff(tr.frames(:))';                              % >1 where a gap was bridged
        steps = diff(tr.pos, 1, 1);
        stepLen = sqrt(sum(steps.^2, 2));
        stepTime = frameGaps(:) * opts.dt;
        tr.velocity = steps ./ max(stepTime, eps);
        % PATH speed: mean |step|/dt. This is an UNSIGNED mean, so position-estimation jitter does
        % NOT cancel -- it biases this upward monotonically, and the bias grows LINEARLY with the
        % tracking rate (v ~ d*f_track). Compare against out.summary.vFloor before reading it as
        % motion; prefer speedNet when quoting a physical speed.
        tr.speed = mean(stepLen ./ max(stepTime, eps));
        tr.lifetime = (tr.frames(end) - tr.frames(1)) * opts.dt;
        tr.pathLength = sum(stepLen);
        tr.netDisplacement = norm(tr.pos(end,:) - tr.pos(1,:));
        tr.straightness = tr.netDisplacement / max(tr.pathLength, eps);
        % NET speed: net displacement over lifetime. Zero-mean jitter cancels, so this is the
        % speed to quote for propagation. Identically speed*straightness.
        tr.speedNet = tr.netDisplacement / max(tr.lifetime, eps);
        tr.birthFrame = tr.frames(1);  tr.deathFrame = tr.frames(end);
        tr.missedFrames = sum(frameGaps - 1);
        [tr.centerVertex, tr.centerPos, ~] = rheome.geom.karcher_mean(S, tr.pos, struct('solver',geodesicSolver));
        descriptor(end+1) = i_public_track(tr); %#ok<AGROW>
    end

    % ---- 4. sort + summarise ----
    if ~isempty(descriptor)
        salience = arrayfun(@(tr) tr.lifetime * mean(tr.strength), descriptor);
        [~, order] = sort(salience, 'descend');
        descriptor = descriptor(order);
    end
    out.tracks = descriptor;
    out.summary.nTracks = numel(descriptor);
    for ty = opts.types
        out.summary.(['n_' ty{1}]) = sum(strcmp({descriptor.type}, ty{1}));
    end
    out.summary.nConfirmed = sum(strcmp({descriptor.state}, 'confirmed') | strcmp({descriptor.state}, 'lost') | strcmp({descriptor.state}, 'ended'));
    out.summary.nSplit = nSplit;  out.summary.nMerge = nMerge;  out.summary.nRecovered = nRecovered;
    out.summary.lifetime = [descriptor.lifetime];
    out.summary.speed    = [descriptor.speed];
    out.summary.speedNet = [descriptor.speedNet];
    out.summary.dt       = opts.dt;
    % ---- speed detection floor -----------------------------------------------------------------
    % Tracks whose NET displacement is under one mean edge went nowhere, yet still report a nonzero
    % PATH speed. That value is position-estimation jitter divided by dt -- i.e. the noise floor of
    % the speed estimate at this tracking rate. Speeds at or below it are not interpretable as
    % motion. It scales as d*f_track, so it rises if you track faster: report it beside any speed.
    if isempty(descriptor)
        out.summary.vFloor = NaN;  out.summary.vFloorJitter = NaN;
    else
        stationary = [descriptor.netDisplacement] < meanEdgeLength;
        if any(stationary)
            out.summary.vFloor = median([descriptor(stationary).speed]);
        else
            out.summary.vFloor = median([descriptor.speed]) * median([descriptor.straightness]);
        end
        out.summary.vFloorJitter = out.summary.vFloor * opts.dt;   % implied per-frame jitter (m)
    end
    out.summary.meanEdgeLength = meanEdgeLength;
end


% ===== association: gated same-type distance, greedy one-per-track/one-per-detection =====
% Returns matched (trackIdx, detIdx) pairs, plus split/merge counts from the viable candidate set
% (a track with >=2 in-gate detections = split; a detection with >=2 in-gate tracks = merge).
function [pairT, pairD, splitCount, mergeCount] = i_associate(tracks, trackIdx, f, detIdx, S, opts, solver)
    pairT = [];  pairD = [];  splitCount = 0;  mergeCount = 0;
    if isempty(trackIdx) || isempty(detIdx), return; end
    trackIdx = trackIdx(:)';  detIdx = detIdx(:)';
    % viable candidates: same type, within the max-step gate
    ci = [];  cj = [];  cost = [];
    for a = 1:numel(trackIdx)
        tr = tracks(trackIdx(a));  lastPos = tr.pos(end,:);  lastType = tr.type;
        for b = 1:numel(detIdx)
            if ~strcmp(f.type{detIdx(b)}, lastType), continue; end
            dstep = i_dist(lastPos, f.pos(detIdx(b),:), S, opts.metric, solver);
            if dstep <= opts.maxStep
                ci(end+1) = a;  cj(end+1) = b;  cost(end+1) = dstep; %#ok<AGROW>
            end
        end
    end
    if isempty(cost), return; end
    % split/merge from the viable set (before greedy pruning), mirroring meshByteTrack.ts
    for a = unique(ci), if numel(unique(cj(ci==a))) >= 2, splitCount = splitCount + 1; end, end
    for b = unique(cj), if numel(unique(ci(cj==b))) >= 2, mergeCount = mergeCount + 1; end, end
    % greedy: smallest cost first, one detection per track and one track per detection
    [~, order] = sort(cost, 'ascend');
    usedT = false(1,numel(trackIdx));  usedD = false(1,numel(detIdx));
    for o = order
        if usedT(ci(o)) || usedD(cj(o)), continue; end
        usedT(ci(o)) = true;  usedD(cj(o)) = true;
        pairT(end+1) = trackIdx(ci(o)); %#ok<AGROW>
        pairD(end+1) = detIdx(cj(o));   %#ok<AGROW>
    end
end

function tracks = i_apply_match(tracks, k, f, d, t, stage, opts)
    tr = tracks(k);
    consecutive = (t == tr.lastFrame + 1);
    tr.frames(end+1)     = t;
    tr.pos(end+1,:)      = f.pos(d,:);
    tr.strength(end+1,1) = f.strength(d);
    tr.stage{end+1}      = stage;
    tr.lastFrame = t;
    if consecutive, tr.consecutive = tr.consecutive + 1; else, tr.consecutive = 1; end
    if strcmp(stage,'low-recovery'), tr.lowRecovery = tr.lowRecovery + 1; end
    if strcmp(tr.state,'lost'), tr.state = 'confirmed'; end                  % reactivate
    if strcmp(tr.state,'tentative') && tr.consecutive >= opts.minimumConsecutiveFrames
        tr.state = 'confirmed';                                             % promote
    end
    tracks(k) = tr;
end

function tr = i_register_miss(tr, t, opts)
    if strcmp(tr.state,'tentative')
        tr.state = 'ended';                                                % tentative-missed
    elseif (t - tr.lastFrame) > opts.lostTrackBuffer
        tr.state = 'ended';                                                % lost-buffer-expired
    else
        tr.state = 'lost';
    end
end

function tr = i_new_track(f, d, t)
    tr = i_empty_tracks();
    tr(1).type = f.type{d};  tr.frames = t;  tr.pos = f.pos(d,:);  tr.strength = f.strength(d);
    tr.stage = {'birth'};  tr.state = 'tentative';  tr.birthFrame = t;  tr.lastFrame = t;
    tr.consecutive = 1;  tr.missed = 0;  tr.lowRecovery = 0;  tr.splitCount = 0;  tr.mergeCount = 0;
end

function tr = i_empty_tracks()
    tr = struct('type',{},'frames',{},'pos',{},'strength',{},'stage',{},'state',{}, ...
        'birthFrame',{},'lastFrame',{},'consecutive',{},'missed',{},'lowRecovery',{}, ...
        'splitCount',{},'mergeCount',{});
end

function d = i_empty_descriptor()
    d = struct('type',{},'frames',{},'pos',{},'strength',{},'velocity',{},'speed',{},'speedNet',{}, ...
        'lifetime',{},'pathLength',{},'netDisplacement',{},'straightness',{}, ...
        'birthFrame',{},'deathFrame',{},'state',{},'missedFrames',{},'lowRecoveryCount',{}, ...
        'centerVertex',{},'centerPos',{});
end

function pub = i_public_track(tr)
    pub.type = tr.type;  pub.frames = tr.frames;  pub.pos = tr.pos;  pub.strength = tr.strength;
    pub.velocity = tr.velocity;  pub.speed = tr.speed;  pub.speedNet = tr.speedNet;  pub.lifetime = tr.lifetime;
    pub.pathLength = tr.pathLength;  pub.netDisplacement = tr.netDisplacement;
    pub.straightness = tr.straightness;  pub.birthFrame = tr.birthFrame;  pub.deathFrame = tr.deathFrame;
    pub.state = tr.state;  pub.missedFrames = tr.missedFrames;  pub.lowRecoveryCount = tr.lowRecovery;
    pub.centerVertex = tr.centerVertex;  pub.centerPos = tr.centerPos;
end

function d = i_dist(posA, posB, S, metric, solver)
    if strcmp(metric,'geodesic')
        [~, vA] = min(sum((S.Vertices - posA).^2, 2));
        distField = rheome.geom.geodesic(S, vA, solver);
        [~, vB] = min(sum((S.Vertices - posB).^2, 2));
        d = distField(vB);
    else
        d = norm(posA - posB);
    end
end

function m = i_frame_max(f)
    if isempty(f.strength), m = 0; else, m = max(f.strength); end
end

function meanEdgeLength = i_mean_edge(S)
    F = double(S.Faces);  V = double(S.Vertices);
    edges = [F(:,[1 2]); F(:,[2 3]); F(:,[3 1])];
    meanEdgeLength = mean(sqrt(sum((V(edges(:,1),:)-V(edges(:,2),:)).^2, 2)));
end

% Author: Diellor Basha, 2026
