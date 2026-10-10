% Author: Diellor Basha, 2026
function out = badchannels(rawDirs, dataRoot, outDir, varargin)
% BADCHANNELS  Pre-registered automatic bad-channel detection on staged Brainstorm raw links.
%
%   out = rheome.qc.badchannels(rawDirs, dataRoot, outDir)            % run, write tables + figures
%   out = rheome.qc.badchannels(rawDirs, dataRoot, outDir, 'DryRun', true)   % list runs + thresholds only
%
%   rawDirs  : cellstr of extracted '@raw*' study folders (each holds data_0raw_*.mat and channel_*.mat),
%              COPIES of the staged protocols -- this function never writes into them.
%   dataRoot : where the BIDS root '/data' of the import container lives on this host.
%   outDir   : output folder (tables, per-channel evidence PNGs).
%
%   Criteria and thresholds are fixed in research/methods/omega-badchannels/PROTOCOL.md, written before
%   the first run. They are the defaults of i_params and must not be tuned on the validation set.
%   ⭐ v2 (PROTOCOL-v2.md, Diellor Basha 2026-10-02): an EMPTY-ROOM recording (rheome.qc.isemptyroom: task-noise,
%   task-emptyroom, sub-emptyroom) is never labelled. Its metrics are still computed and saved, and the
%   rules that WOULD have fired are kept as rawFlag/rawWhy for the record, but flag is all false. Subject
%   recordings are untouched: same criteria, same thresholds, bit-identical flags to v1.
%   ⚠ v1 marked 772 empty-room channels across one resting cohort, 652 of them PSD-only: an empty-room spectrum is
%   sensor noise alone, so a run-relative PSD z there measures sensor noise floors, not artefacts.
%   ⭐ v3 (2026-10-08): the 65–115 Hz rule (PSD_HIGH) never labels a channel on its own, it only
%   corroborates another rule (rheome.qc.badchannelrules). A PSD_HIGH-only channel stays in rawFlag/rawWhy
%   and is counted in runs.csv nRulesFiredUnlabelled. Metrics and every other rule are unchanged from v2.
%   Only Brainstorm type 'MEG' is evaluated; MEG REF / system / EEG channels are never touched.
%   ⚠ process_detectbad and process_detectbad_mad flag SEGMENTS, not channels -- they do not fit here.
%   Requires Brainstorm on the path and started (brainstorm server) for in_fread.
%
% See also: rheome.qc.badchannelrules, rheome.qc.isemptyroom, in_fread, process_detectbad_mad

ip = inputParser;
ip.addParameter('DryRun', false);
ip.addParameter('Seed', 20261001);
ip.parse(varargin{:});
P = i_params();
if ~exist(outDir, 'dir'), mkdir(outDir); end
figDir = fullfile(outDir, 'figures'); if ~exist(figDir, 'dir'), mkdir(figDir); end

flagRows = {}; runRows = {}; ctlRows = {};
rng(ip.Results.Seed);
ctlDir = struct();
for k = numel(rawDirs):-1:1                          % reverse so the FIRST match wins
    sb = strrep(regexp(rawDirs{k}, 'sub-[A-Za-z0-9]+', 'match', 'once'), '-', '_');
    if ~isfield(ctlDir, sb) || contains(rawDirs{k}, 'rest') || ~contains(ctlDir.(sb), 'rest')
        if ~isfield(ctlDir, sb) || contains(rawDirs{k}, 'rest'), ctlDir.(sb) = rawDirs{k}; end
    end
end
fprintf('BADCHAN params: %s\n', jsonencode(P));
for k = 1:numel(rawDirs)
    d = rawDirs{k};
    fd = dir(fullfile(d, 'data_0raw_*.mat')); fc = dir(fullfile(d, 'channel_*.mat'));
    if isempty(fd) || isempty(fc), fprintf('SKIP %s (no raw link/channel file)\n', d); continue; end
    R = load(fullfile(d, fd(1).name), 'F'); sFile = R.F;
    sFile.filename = regexprep(sFile.filename, '^/data', dataRoot);
    if isfield(sFile.header, 'meg4_files')            % ⚠ ctf_seek reads these, not sFile.filename
        sFile.header.meg4_files = regexprep(sFile.header.meg4_files, '^/data', dataRoot);
    end
    ChannelMat = load(fullfile(d, fc(1).name));
    types = {ChannelMat.Channel.Type};
    iMeg = find(strcmp(types, 'MEG'));
    [~, runName] = fileparts(fd(1).name); runName = erase(runName, 'data_0raw_');
    sub = regexp(runName, 'sub-[A-Za-z0-9]+', 'match', 'once');
    tsvBad = i_tsvbad(sFile.filename);
    sfreq = sFile.prop.sfreq; dur = diff(sFile.prop.times);
    fprintf('RUN %-50s MEG=%3d EEG=%2d sfreq=%g dur=%.1fs exists=%d tsvBad={%s} importFlagBad=%d\n', ...
        runName, numel(iMeg), sum(strcmp(types,'EEG')), sfreq, dur, exist(sFile.filename,'file')==2, ...
        strjoin(tsvBad, ','), sum(sFile.channelflag == -1));
    if ip.Results.DryRun, continue; end

    M = i_metrics(sFile, ChannelMat, iMeg, P);
    [flag, why, rawFlag, rawWhy] = rheome.qc.badchannelrules(M, P);
    emptyRoom = rheome.qc.isemptyroom(runName);
    if emptyRoom                                      % v2: empty-room recordings get no label from ANY rule
        flag(:) = false; why(:) = {''};
        fprintf('EMPTYROOM %s: not labelled (v1 rules would have flagged %d: %s)\n', runName, sum(rawFlag), ...
            i_ruleslist({ChannelMat.Channel(iMeg).Name}, rawFlag, rawWhy));
    elseif any(rawFlag & ~flag)                       % v3: PSD_HIGH alone, not labelled
        fprintf('PSDHIGHONLY %s: not labelled: %s\n', runName, ...
            i_ruleslist({ChannelMat.Channel(iMeg).Name}, rawFlag & ~flag, rawWhy));
    end
    names = {ChannelMat.Channel(iMeg).Name};
    nFlag = sum(flag); safety = nFlag > P.maxFrac * numel(iMeg);
    runRows(end+1, :) = {sub, runName, numel(iMeg), nFlag, safety, strjoin(names(flag), ' '), strjoin(tsvBad, ' '), ...
        emptyRoom, sum(rawFlag & ~flag)}; %#ok<AGROW>
    snip = i_snippet(sFile, ChannelMat, iMeg, P);
    for c = find(flag(:))'
        flagRows(end+1, :) = {sub, runName, names{c}, why{c}, M.z_sd(c), M.r_sd(c), M.jumpfrac(c), ...
            M.z_lo(c), M.r_lo(c), M.z_hi(c), M.r_hi(c), ismember(names{c}, tsvBad)}; %#ok<AGROW>
        i_figure(M, snip, c, names, why{c}, sprintf('%s  %s', runName, names{c}), ...
            fullfile(figDir, sprintf('FLAG_%s_%s.png', runName, names{c})), flag);
    end
    % controls: 3 unflagged channels on each subject's first rest run (first run if it has none)
    if strcmp(d, ctlDir.(strrep(sub, '-', '_')))
        pool = find(~flag); pick = pool(randperm(numel(pool), min(3, numel(pool))));
        for c = pick(:)'
            ctlRows(end+1, :) = {sub, runName, names{c}, M.z_sd(c), M.r_sd(c), M.jumpfrac(c), M.z_lo(c), M.z_hi(c)}; %#ok<AGROW>
            i_figure(M, snip, c, names, 'CONTROL (unflagged)', sprintf('%s  %s', runName, names{c}), ...
                fullfile(figDir, sprintf('CTRL_%s_%s.png', runName, names{c})), flag);
        end
    end
    % tsv-marked bads that the detector did not flag: show them too
    for nm = tsvBad(:)'
        c = find(strcmp(names, nm{1}));
        if ~isempty(c) && ~flag(c)
            i_figure(M, snip, c, names, 'channels.tsv BAD, NOT flagged', sprintf('%s  %s', runName, nm{1}), ...
                fullfile(figDir, sprintf('TSV_%s_%s.png', runName, nm{1})), flag);
        end
    end
    save(fullfile(outDir, sprintf('metrics_%s.mat', runName)), 'M', 'flag', 'why', 'names', 'P', 'tsvBad', ...
        'emptyRoom', 'rawFlag', 'rawWhy');
end
if ip.Results.DryRun, out = struct('dryrun', true); return; end
out.runs = i_table(runRows, {'subject','run','nMEG','nFlagged','safetyHold','flagged','tsvBad', ...
    'emptyRoom','nRulesFiredUnlabelled'});
out.flags = i_table(flagRows, {'subject','run','channel','rules','z_logSD','r_SD', ...
    'jumpFrac','z_PSD_1_45','r_PSD_1_45','z_PSD_65_115','r_PSD_65_115','tsvBad'});
out.controls = i_table(ctlRows, {'subject','run','channel','z_logSD','r_SD','jumpFrac','z_PSD_1_45','z_PSD_65_115'});
writetable(out.runs, fullfile(outDir, 'runs.csv'));
writetable(out.flags, fullfile(outDir, 'flags.csv'));
writetable(out.controls, fullfile(outDir, 'controls.csv'));
disp(out.runs); disp(out.flags);
end

% ---------------------------------------------------------------------------------------------
function s = i_ruleslist(names, flag, why)
% 'MLT55:JUMPY MRO12:PSD_LOW+PSD_HIGH' for the flagged channels; '' when none
% ⚠ strcat of a 0x0 and a 1x0 cell errors -- the v2 job died on the first empty-room run with no flags
iF = reshape(find(flag), 1, []);
s = strjoin(arrayfun(@(c) sprintf('%s:%s', names{c}, strrep(why{c}, ' ', '+')), iF, 'UniformOutput', false), ' ');
end

function T = i_table(rows, vars)
% ⚠ cell2table on a 0x0 cell with VariableNames errors -- v1 died here when nothing was flagged anywhere
if isempty(rows), T = cell2table(cell(0, numel(vars)), 'VariableNames', vars);
else, T = cell2table(rows, 'VariableNames', vars); end
end

function P = i_params()
P.version = 'v3'; P.emptyRoom = 'never labelled'; P.psdHigh = 'corroborates only';
P.nWin = 60; P.winSec = 2; P.k = 6;
P.zFlat = -5; P.rFlat = 0.2; P.zNoisy = 5; P.rNoisy = 3;
P.zJump = 5; P.jumpFrac = 0.20;
P.zPsd = 5; P.lrPsd = 0.5; P.lo = [1 45]; P.hi = [65 115];
P.maxFrac = 0.10; P.snipSec = 5;
end

function M = i_metrics(sFile, ChannelMat, iMeg, P)
sf = sFile.prop.sfreq; n = round(P.winSec * sf);
s0 = round(sFile.prop.times(1) * sf); s1 = round(sFile.prop.times(2) * sf);
starts = round(linspace(s0, s1 - n + 1, P.nWin)); starts = unique(starts);
opt = db_template('ImportOptions'); opt.DisplayMessages = 0;
nC = numel(iMeg); nW = numel(starts);
sd = zeros(nC, nW); jmp = zeros(nC, nW);
nfft = round(sf); hw = 0.5 - 0.5*cos(2*pi*(0:nfft-1)'/(nfft-1)); f = (0:nfft/2)' * sf / nfft;
pxx = zeros(numel(f), nC);
for w = 1:nW
    F = in_fread(sFile, ChannelMat, 1, [starts(w), starts(w)+n-1], [], opt);
    X = detrend(F(iMeg, :)')';                       % [nC x n]
    sd(:, w) = std(X, 0, 2); jmp(:, w) = max(abs(diff(X, 1, 2)), [], 2);
    for s = 1:nfft/2:(n - nfft + 1)
        Y = fft(X(:, s:s+nfft-1)' .* hw); pxx = pxx + abs(Y(1:numel(f), :)).^2;
    end
end
M.f = f; M.psd = pxx;
loc = zeros(3, nC);
for c = 1:nC, L = ChannelMat.Channel(iMeg(c)).Loc; loc(:, c) = mean(L(:, 1:min(4, size(L,2))), 2); end
D = sqrt(sum((permute(loc, [2 3 1]) - permute(loc, [3 2 1])).^2, 3)); D(1:nC+1:end) = inf;
[~, ord] = sort(D, 2); M.nb = ord(:, 1:P.k); M.loc = loc;
M.sd = median(sd, 2);
M.z_sd = i_rz(log10(max(M.sd, realmin))); M.r_sd = i_nr(M.sd, M.nb);
zj = zeros(nC, nW); for w = 1:nW, zj(:, w) = i_rz(log10(max(jmp(:, w), realmin))); end
M.jumpfrac = mean(zj > P.zJump, 2);
lo = f >= P.lo(1) & f <= P.lo(2); hi = f >= P.hi(1) & f <= P.hi(2);
plo = mean(log10(max(pxx(lo, :), realmin)), 1)'; phi = mean(log10(max(pxx(hi, :), realmin)), 1)';
M.z_lo = i_rz(plo); M.r_lo = i_nr(10.^plo, M.nb); M.z_hi = i_rz(phi); M.r_hi = i_nr(10.^phi, M.nb);
end

function z = i_rz(x)
x = x(:); m = median(x); s = 1.4826 * median(abs(x - m)); z = (x - m) / max(s, eps);
end

function r = i_nr(x, nb)
x = x(:); r = x ./ max(median(x(nb), 2), realmin);
end

function S = i_snippet(sFile, ChannelMat, iMeg, P)
sf = sFile.prop.sfreq; n = round(P.snipSec * sf);
mid = round(mean(sFile.prop.times) * sf);
opt = db_template('ImportOptions'); opt.DisplayMessages = 0;
F = in_fread(sFile, ChannelMat, 1, [mid, mid + n - 1], [], opt);
S.X = detrend(F(iMeg, :)')'; S.t = (0:n-1) / sf;
end

function i_figure(M, S, c, names, why, ttl, png, flag)
fig = figure('Visible', 'off', 'Position', [0 0 1500 420], 'Color', 'w');
nb = M.nb(c, :); nbc = nb(~flag(nb)); if isempty(nbc), nbc = nb; end   % flagged neighbours excluded from the band
lf = log10(max(M.psd, realmin)); keep = M.f >= 1 & M.f <= 200; f = M.f(keep);
a = subplot(1, 3, 1); hold(a, 'on');
q = prctile(lf(keep, :), [5 95], 2);
fill(a, [f; flipud(f)], [q(:,1); flipud(q(:,2))], [.85 .85 .85], 'EdgeColor', 'none');
fill(a, [f; flipud(f)], [min(lf(keep, nbc), [], 2); flipud(max(lf(keep, nbc), [], 2))], [.6 .75 1], ...
    'EdgeColor', 'none', 'FaceAlpha', .7);
plot(a, f, lf(keep, c), 'r', 'LineWidth', 1.2); set(a, 'XScale', 'log'); xlim(a, [1 200]);
xlabel(a, 'Hz'); ylabel(a, 'log_{10} power (run-relative)');
legend(a, {'array 5–95%', sprintf('%d unflagged neighbours min–max', numel(nbc)), names{c}}, 'Location', 'southwest', 'Box', 'off');
title(a, 'PSD vs neighbours');
b = subplot(1, 3, 2); hold(b, 'on');
show = [c nbc(1:min(2, numel(nbc)))]; col = {'r', 'b', [0 .5 0]}; lab = cell(1, numel(show));
for i = 1:numel(show)                                % each trace on its OWN scale (+-4 SD), SD printed
    sdi = std(S.X(show(i), :)); plot(b, S.t, S.X(show(i), :) / (4*sdi) - 2*(i-1), 'Color', col{i});
    lab{i} = sprintf('%s  SD=%.0f fT', names{show(i)}, sdi*1e15);
end
set(b, 'YTick', -2*(numel(show)-1):2:0, 'YTickLabel', fliplr(lab)); xlabel(b, 's'); ylim(b, [-2*numel(show)+1 1]);
title(b, sprintf('5-s snippet (mid-run), each trace +-4 own SD; array median SD = %.0f fT', median(M.sd)*1e15));
d = subplot(1, 3, 3); hold(d, 'on');
L = M.loc - mean(M.loc, 2); r = sqrt(sum(L(1:2,:).^2)); th = atan2(L(2,:), L(1,:));
el = atan2(r, L(3,:)); x = el .* cos(th); y = el .* sin(th);
scatter(d, x, y, 12, [.6 .6 .6], 'filled'); scatter(d, x(nb), y(nb), 30, 'b', 'filled');
scatter(d, x(c), y(c), 70, 'r', 'filled'); axis(d, 'equal', 'off');
title(d, 'sensor position (top view, nose +x)');
sgtitle(fig, {ttl, sprintf('%s | z(logSD)=%.1f  r(SD)=%.2f  jump=%.2f  z(PSD1-45)=%.1f r=%.2f  z(PSD65-115)=%.1f r=%.2f', ...
    why, M.z_sd(c), M.r_sd(c), M.jumpfrac(c), M.z_lo(c), M.r_lo(c), M.z_hi(c), M.r_hi(c))}, ...
    'Interpreter', 'none', 'FontSize', 10);
print(fig, png, '-dpng', '-r100'); close(fig);
end

function bad = i_tsvbad(meg4)
bad = {};
tsv = regexprep(fileparts(meg4), '_meg\.ds$', '_channels.tsv');
if ~exist(tsv, 'file'), return; end
T = readtable(tsv, 'FileType', 'text', 'Delimiter', '\t');
if ~ismember('status', T.Properties.VariableNames), return; end
bad = T.name(strcmpi(strtrim(string(T.status)), 'bad'))';
bad = cellstr(bad);
end
% Author: Diellor Basha, 2026
