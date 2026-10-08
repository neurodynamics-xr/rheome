function file = preview(src, opts)
% INGEST.PREVIEW  A min/max pyramid of a recording, beside it, with no transform at all.
%
%   file = rheome.ingest.preview('sub01')
%   file = rheome.ingest.preview(name, Floor=1/60, Channels="MEG", Overwrite=true)
%
% ⭐ WHAT IT IS FOR. Before deciding what to analyse you want to LOOK: browse every channel
% of a recording, find the stretch worth having, and only then load it. That needs one thing
% and nothing else -- the minimum and maximum of each channel over each tile. Those two merge
% by extremum, so the pair at any level bounds every sample under it, and a band drawn from
% them is an exact preview of the trace at any zoom. No bank, no bands, no moments, no
% coefficients: a reshape and two extrema per chunk.
%
% ⭐ AND IT SETS THE ANALYSIS. Choosing a window in the preview chooses a window LENGTH, and
% a window length is a statement about which bands can be resolved in it -- constant Q means
% a band needs about rheome.ingest.support seconds. So the explorer's selection already names the
% band support of whatever comes next (rheome.ingest.support, printed by the browser).
%
% ⚠ THE FLOOR IS THE WHOLE COST. 8 bytes per channel per level-0 tile, doubled by the
% pyramid. Measured on a 300-channel, 600 s, 600 Hz recording (868 MB):
%
%     floor    level 0    pyramid    of the recording    pixel-exact down to
%     1 s       1.4 MB     2.9 MB        0.2 %              1114 s
%     0.25 s    5.8 MB    11.5 MB        0.7 %               278 s
%     1/60 s   86.4 MB   172.8 MB       10   %                19 s
%     1/150 s 216.0 MB   432.0 MB       25   %                 7 s
%
% The last column is the rule: a preview stays pixel-exact down to (axis pixels x floor)
% seconds. Choose the floor so that meets the window where raw reads become cheap -- about
% 20 s / 1100 px = 1/60 s here, which is why that is the default.
%
% ⚠ IT IS A TILE STORE WITH NO BANDS. rheome.select.open reads it, rheome.select.derive gives min, max and
% ptp from it, and @recordingbrowser opens on it in envelope mode. Everything band-shaped is
% absent by construction, and says so rather than returning zeros.
%
% OPTIONS
%   Floor      level-0 tile in seconds (1/60); Floor*fs must be an integer
%   Channels   "MEG" (good channels of that type) or explicit row indices
%   Store      "default" | "native", as pagedrecording
%   Levels     "all" (the pyramid) | "base" (level 0 only; roll up on read, 13 ms a channel)
%   ChunkSeconds  seconds of every channel read at once (60)
%   Overwrite, Verbose
%
% See also: rheome.ingest.build, rheome.ingest.support, rheome.select.open, rheome.recordingbrowser
%
% Author: Diellor Basha, 2026

    arguments
        src
        opts.Floor        (1,1) double {mustBePositive} = 1/60
        opts.Channels     = "MEG"
        opts.Store        (1,1) string {mustBeMember(opts.Store, ["default","native"])} = "default"
        opts.Levels       (1,1) string {mustBeMember(opts.Levels, ["all","base"])} = "all"
        opts.ChunkSeconds (1,1) double {mustBePositive} = 60
        opts.Overwrite    (1,1) logical = false
        opts.Verbose      (1,1) logical = true
    end

    t0 = tic;
    % ⚠ READ IN DOUBLE. pagedrecording defaults to single, and an extremum of the
    % single-rounded signal is not a bound on the original: containment failed on three
    % quarters of the tiles until this was explicit. rheome.ingest.build reads double for the same
    % reason.
    pr0 = rheome.pagedrecording(src, 'Store', char(opts.Store), 'Precision', 'double');
    nT = pr0.NumSamples;  fs = pr0.SamplingFrequency;

    if isnumeric(opts.Channels)
        iSel = double(opts.Channels(:));
    else
        typ = char(string(opts.Channels));
        iSel = find(strcmp(pr0.ChannelType(:), typ) & pr0.ChannelFlag(:) == 1);
        if isempty(iSel)
            error('ingest:preview:channels', 'No good channels of type ''%s'' in %s.', typ, pr0.File);
        end
    end
    C = numel(iSel);

    cfg = rheome.ingest.config(FrameFloor=opts.Floor, Space=false);
    g = rheome.ingest.grid(nT, fs, cfg);                                 % the same dyadic grid as a tile store
    g.bandsAt = repmat({zeros(1, 0)}, 1, g.Lmax + 1);             % no bands, at any level
    g.channelLevel = 0;  g.envelopeLevel = 0;
    g.spaceSlots = table();  g.chanSpace = [];

    outDir = fileparts(pr0.File);
    tag = '';  if strcmp(opts.Store, "native"), tag = '__native'; end
    file = fullfile(outDir, sprintf('preview__F%g%s.mat', opts.Floor, tag));
    if exist(file, 'file') == 2 && ~opts.Overwrite
        m = matfile(file);  old = m.meta;
        if strcmp(char(old.hash), char(cfg.hash)) && old.C == C
            if opts.Verbose, fprintf('rheome.ingest.preview: %s is current (%d ch, %g s tiles)\n', file, C, opts.Floor); end
            return
        end
        error('ingest:preview:hash', ...
              'A different preview is at %s (floor %g, %d channels). Overwrite=true to replace it.', ...
              file, old.cfg.FrameFloor, old.C);
    end

    % ---- level 0: read time chunks of EVERY channel, reshape, take the two extrema ----
    % ⚠ Whole-channel reads would pull 8 bytes x nT x C; a time chunk is bounded and is also
    % how the store is laid out, so this is one sequential pass over the recording.
    F = g.F;  K0 = g.K0;
    mn = zeros(K0, C, 'single');  mx = zeros(K0, C, 'single');
    n  = zeros(K0, 1, 'uint32');
    chunk = max(F, round(opts.ChunkSeconds * fs / F) * F);        % a whole number of tiles
    for a = 1:chunk:nT
        b = min(a + chunk - 1, nT);
        X = read(pr0, a, b);                                      % [nChAll x nSamp], double
        X = X(iSel, :);
        nk = b - a + 1;
        k0 = (a - 1) / F;
        full = floor(nk / F);
        if full > 0
            Y = reshape(X(:, 1:full*F), C, F, full);
            mn(k0 + (1:full), :) = rheome.ingest.bound(squeeze(min(Y, [], 2))', 'down');
            mx(k0 + (1:full), :) = rheome.ingest.bound(squeeze(max(Y, [], 2))', 'up');
            n(k0 + (1:full)) = uint32(F);
        end
        rest = nk - full*F;
        if rest > 0                                               % the record's last, partial tile
            Y = X(:, full*F+1:end);
            mn(k0 + full + 1, :) = rheome.ingest.bound(min(Y, [], 2)', 'down');
            mx(k0 + full + 1, :) = rheome.ingest.bound(max(Y, [], 2)', 'up');
            n(k0 + full + 1) = uint32(rest);
        end
        if opts.Verbose
            fprintf('  %.0f-%.0f s of %.0f  (%.0f s)\n', (a-1)/fs, b/fs, nT/fs, toc(t0));
        end
    end

    % ---- the pyramid: extremum of the two children, level by level ----
    S = struct('meta', i_meta(pr0, src, fs, nT, C, iSel, cfg, opts), 'bands', i_nobands(), ...
               'frame', i_noframe(), 'grid', g);
    S.L00_n = n;  S.L00_min = mn;  S.L00_max = mx;
    top = g.Lmax;  if opts.Levels == "base", top = 0; end
    for L = 1:top
        [n, mn, mx] = i_pair(n, mn, mx, g.K(L+1));
        S.(sprintf('L%02d_n', L)) = n;
        S.(sprintf('L%02d_min', L)) = mn;
        S.(sprintf('L%02d_max', L)) = mx;
    end
    S.meta.levels = top;
    S.meta.buildSeconds = toc(t0);

    tmp = [file '.tmp'];
    save(tmp, '-struct', 'S', '-v7.3', '-nocompression');
    movefile(tmp, file, 'f');
    if opts.Verbose
        d = dir(file);  r = dir(pr0.File);
        fprintf('rheome.ingest.preview: %s (%.0f MB, %.2f %% of the recording) in %.0f s\n', ...
                file, d.bytes/1e6, 100*d.bytes/r.bytes, S.meta.buildSeconds);
    end
end

% A parent is the extremum of its two children and the sum of their counts. Padding the odd
% tail with the neutral elements keeps the pairing exact rather than special-casing it.
function [n2, mn2, mx2] = i_pair(n, mn, mx, K2)
    K = size(mn, 1);
    if mod(K, 2), n(K+1, :) = 0;  mn(K+1, :) = Inf;  mx(K+1, :) = -Inf; end
    n2  = n(1:2:end, :) + n(2:2:end, :);
    mn2 = min(mn(1:2:end, :), mn(2:2:end, :));
    mx2 = max(mx(1:2:end, :), mx(2:2:end, :));
    n2 = n2(1:K2, :);  mn2 = mn2(1:K2, :);  mx2 = mx2(1:K2, :);
end

function meta = i_meta(pr0, src, fs, nT, C, iSel, cfg, opts)
    meta = struct();
    meta.source = pr0.File;
    if rheome.load.has(src)
        meta.name = char(src);
    else
        [d0, stem] = fileparts(pr0.File);  [~, dn] = fileparts(d0);
        meta.name = [dn '_' stem];
    end
    meta.kind = 'preview';                                        % not a band store, and says so
    meta.fs = fs;  meta.nT = nT;  meta.duration = nT / fs;  meta.C = C;
    meta.iChannel = iSel(:)';
    meta.ChannelName = pr0.ChannelName(iSel);
    meta.ChannelType = pr0.ChannelType(iSel);
    meta.units = 'T';
    meta.cfg = cfg;  meta.hash = cfg.hash;
    meta.paged = false;  meta.plan = table();
    meta.release = version;  meta.commit = '';
    meta.created = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
    meta.tree = [];  meta.groupNodes = [];  meta.groupChannelRows = [];
    meta.sbands = table();  meta.sframe = struct();
    meta.masked = false;  meta.nValid = double(nT) * C;
end

function b = i_nobands()
    b = table('Size', [0 12], ...
        'VariableTypes', {'double','cell','double','double','double','double','double','double','double','logical','double','cell'}, ...
        'VariableNames', {'j','scales','fLo','fHi','fCenter','fExtent','fLoMeasured','fHiMeasured','tSupport','partial','naturalLevel','kind'});
end

function f = i_noframe()
    f = struct('fc', [], 'nF', 0, 'fmin', NaN, 'fmax', NaN, 'f', [], 'S', [], ...
               'A', NaN, 'B', NaN, 'interior', [NaN NaN], 'Q', NaN, 'coiConst', NaN, ...
               'wavelet', 'none', 'bank', 'preview', 'support', struct('timeTimesFc', NaN), 'probeLength', 0);
end
% Author: Diellor Basha, 2026
