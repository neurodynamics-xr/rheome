function fx = selFixture()
% SELFIXTURE  A synthetic recording store and its tile store in a temp dir, for the select tests.
%   fx.dir .rec .store .fs .nT .C .X .cfg .db
% 20 s at 200 Hz, four MEG channels, a 10 Hz rhythm with three bursts and a 3 Hz rhythm.
% Cached per MATLAB session (persistent), since the build takes a second.
%
% Author: Diellor Basha, 2026

    persistent cached
    if ~isempty(cached) && exist(cached.store, 'file') == 2
        fx = cached;  fx.db = rheome.select.open(fx.store);  return
    end
    fs = 200;  nT = 4000;  C = 4;  t = (0:nT-1)' / fs;
    rng(31);
    X = 1e-13 * randn(nT, C);
    burst = zeros(nT, 1);
    for on = [3 9 15]                                          % three 1 s bursts of alpha
        burst = burst + (t >= on & t < on + 1);
    end
    X(:, 1) = X(:, 1) + 4e-13 * sin(2*pi*10*t) .* burst;
    X(:, 2) = X(:, 2) + 3e-13 * sin(2*pi*3*t);
    X(:, 3) = X(:, 3) + 4e-13 * sin(2*pi*10*t) .* burst .* (t < 10);
    d = tempname;  mkdir(d);
    rec = fullfile(d, 'recording.mat');
    F = X.';  Time = t';  sfreq = fs;                                          %#ok<NASGU>
    ChannelFlag = ones(C, 1);  ChannelName = {'A','B','C','D'};  ChannelType = repmat({'MEG'}, 1, C); %#ok<NASGU>
    nCh = C;  Comment = 'select fixture';  Events = [];  nAvg = 1;  source = struct(); %#ok<NASGU>
    save(rec, 'F', 'Time', 'sfreq', 'ChannelFlag', 'ChannelName', 'ChannelType', ...
         'nCh', 'nT', 'Comment', 'Events', 'nAvg', 'source', '-v7.3', '-nocompression');
    cfg = rheome.ingest.config(ChannelMinTile=0);                                  % channels at every level: the tests' reference
    P = [0 0 0; 0.03 0 0; 0 0.03 0; 0.03 0.03 0];                       % a 2 x 2 array, 3 cm pitch
    store = rheome.ingest.build(rec, cfg, Verbose=false, Groups=P);
    fx = struct('dir', d, 'rec', rec, 'store', store, 'fs', fs, 'nT', nT, 'C', C, 'X', X, 'cfg', cfg, 't', t, 'P', P);
    cached = fx;
    fx.db = rheome.select.open(store);
end
% Author: Diellor Basha, 2026
