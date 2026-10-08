function [B, info] = bases(name, tau, K, Klbo)
% LOAD.BASES  Load the per-hemisphere operators & eigenbases (all families, both hemispheres).
%
%   B = rheome.load.bases(name)                    % Dirac slice for (tau,K)=(0.5,400); richest LBO cache
%   B = rheome.load.bases(name, tau, K)            % choose which cached DIRAC eigenbasis to slice
%   B = rheome.load.bases(name, tau, K, Klbo)      % REQUIRE an LBO basis with Klbo modes per hemisphere
%   [B, info] = rheome.load.bases(...)             % info: which file, which K, lambda_max, sigma floor
%
% Returns the struct cached by rheome.import.bases, and -- for a unified per-hemisphere view across ALL
% THREE operator families -- attaches the Dirac eigenbasis of each hemisphere by slicing the
% whole-cortex Dirac cache (rheome.utils.hemisphere), so it is not duplicated on disk:
%   B.(L|R).lbo    struct(.L .M  .Phi .Lambda .Mass)     Laplace–Beltrami
%   B.(L|R).conn   struct(.C  .Psi .Lambda)              complex Connection Laplacian
%   B.(L|R).dirac  relative-Dirac eigenbasis (if the (tau,K) Dirac cache exists)
%   B.(L|R).gv .S  global vertices + slim hemisphere surface
%
% ⚠ K AND Klbo ARE DIFFERENT THINGS. K selects which DIRAC cache file to slice; Klbo is how many
% LBO modes per hemisphere the returned basis must have. The old signature conflated them: a caller
% asking for K = 400 got whatever LBO K happened to be on disk, silently. Pass Klbo whenever the
% lambda axis matters -- which is whenever a scale, a size or a speed will be read off it -- and a
% mismatch becomes an ERROR rather than a different answer.
%
% RESOLUTION ORDER for the LBO cache:
%   1. bases__K<Klbo>__Kc*.mat        if Klbo was requested
%   2. the keyed file with the LARGEST K, if Klbo was not requested
%   3. legacy unkeyed bases.mat
%
% See also: rheome.import.bases, rheome.load.dirac, rheome.utils.hemisphere
%
% Author: Diellor Basha, 2026

    if nargin < 2 || isempty(tau),  tau  = 0.5; end
    if nargin < 3 || isempty(K),    K    = 400; end
    if nargin < 4,                  Klbo = [];  end

    dsdir = fullfile(rheome.load.root(), char(name));
    f = i_resolve(dsdir, Klbo);
    if isempty(f)
        error('load:bases:missing', ...
            'No bases cached for ''%s''%s. Run  rheome.import.bases(''%s''%s)  first.', ...
            name, i_tern(isempty(Klbo), '', sprintf(' with K = %d', Klbo)), ...
            name, i_tern(isempty(Klbo), '', sprintf(', %d', Klbo)));
    end
    B = getfield(builtin('load', f, 'bases'), 'bases');

    % ---- verify we got the basis that was asked for ----
    gotK = NaN;
    if isfield(B,'meta') && isfield(B.meta,'K'), gotK = B.meta.K; end
    if isnan(gotK)
        labs = intersect(fieldnames(B), {'L','R'}, 'stable');
        if ~isempty(labs), gotK = numel(B.(labs{1}).lbo.Lambda); end
    end
    if ~isempty(Klbo) && gotK ~= Klbo
        error('load:bases:K', ...
            ['requested %d LBO modes per hemisphere but %s holds %d. Nothing downstream would ' ...
             'have caught this -- it would simply have run on a different lambda axis. Build it ' ...
             'with rheome.import.bases(''%s'', %d).'], Klbo, f, gotK, name, Klbo);
    end

    info = struct('file', f, 'K', gotK);
    if isfield(B,'meta') && isfield(B.meta,'lambdaMax')
        info.lambdaMax = B.meta.lambdaMax;  info.sigmaFloor = B.meta.sigmaFloor;
    else
        labs = intersect(fieldnames(B), {'L','R'}, 'stable');
        lm = max(cellfun(@(l) max(B.(l).lbo.Lambda), labs));
        info.lambdaMax = lm;  info.sigmaFloor = 3.08/sqrt(lm);
    end

    % attach the per-hemisphere Dirac slice from the whole-cortex Dirac cache (no duplication)
    dfile = fullfile(dsdir, sprintf('dirac__tau%.2f__K%d.mat', tau, K));
    if ~exist(dfile, 'file')
        g = dir(fullfile(dsdir, sprintf('dirac__tau%.2f__K*.mat', tau)));
        if ~isempty(g), dfile = fullfile(g(1).folder, g(1).name); end
    end
    if exist(dfile, 'file')
        dbG = getfield(builtin('load', dfile, 'dbasis'), 'dbasis');
        Sfile = fullfile(dsdir, 'surface.mat');
        S = getfield(builtin('load', Sfile, 'S'), 'S');
        for h = 1:numel(B.hemi)
            lab = upper(regexp(B.hemi{h}, '[LR]$', 'match', 'once'));
            if isempty(lab) || ~isfield(B, lab), continue; end
            [~, dbh] = rheome.utils.hemisphere(S, dbG, h);
            B.(lab).dirac = dbh;
        end
    end
end

% ----- helpers -----
function f = i_resolve(dsdir, Klbo)
    f = '';
    if ~isempty(Klbo)
        g = dir(fullfile(dsdir, sprintf('bases__K%d__Kc*.mat', Klbo)));
        if ~isempty(g), f = fullfile(g(1).folder, g(1).name); return; end
        % Fall back to the legacy unkeyed cache, which predates the naming and may well BE the
        % requested K. Do not substitute a different keyed file -- the caller's meta.K check either
        % accepts this one or reports the mismatch by name.
        legacy = fullfile(dsdir, 'bases.mat');
        if exist(legacy, 'file'), f = legacy; end
        return;
    end
    g = dir(fullfile(dsdir, 'bases__K*__Kc*.mat'));
    if ~isempty(g)
        kk = zeros(1, numel(g));
        for i = 1:numel(g)
            t = regexp(g(i).name, 'bases__K(\d+)__Kc', 'tokens', 'once');
            if ~isempty(t), kk(i) = str2double(t{1}); end
        end
        [~, i] = max(kk);                          % richest basis available
        f = fullfile(g(i).folder, g(i).name);  return;
    end
    legacy = fullfile(dsdir, 'bases.mat');         % pre-keying cache
    if exist(legacy, 'file'), f = legacy; end
end

function y = i_tern(c,a,b), if c, y=a; else, y=b; end, end

% Author: Diellor Basha, 2026
