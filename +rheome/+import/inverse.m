function file = inverse(name, studyDir, pattern)
% IMPORT.INVERSE  Find a Brainstorm kernel-only inverse in a study + cache it to +data.
%
%   file = rheome.import.inverse(name, studyDir)             % auto: the unconstrained MN kernel
%   file = rheome.import.inverse(name, studyDir, pattern)    % e.g. 'results_dSPM*KERNEL*.mat'
%
% Globs the study folder for a Brainstorm results_*_KERNEL_*.mat, PREFERS an unconstrained
% (nComponents == 3) match -- so it compares like-for-like with the module's unconstrained
% Dirac field -- reads it via rheome.io.read.inverse, and writes +data/<name>/inverse_bst.mat.
% Read back through rheome.load.inverse.
%
% See also: rheome.io.read.inverse, rheome.load.inverse, rheome.import.study
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(pattern), pattern = 'results_MN_*KERNEL*.mat'; end
    dsdir = fullfile(rheome.load.root(), char(name));
    if ~exist(dsdir, 'dir'), mkdir(dsdir); end

    L = dir(fullfile(studyDir, pattern));
    L = L(~startsWith({L.name}, '._'));                        % drop macOS ._ sidecars
    if isempty(L), error('import:inverse:none', 'No %s in %s', pattern, studyDir); end

    chosen = '';                                               % prefer an unconstrained kernel
    for i = 1:numel(L)
        m = matfile(fullfile(studyDir, L(i).name));            % read nComponents only (no kernel load)
        if double(m.nComponents) == 3, chosen = fullfile(studyDir, L(i).name); break; end
    end
    if isempty(chosen), [~, k] = max([L.datenum]); chosen = fullfile(studyDir, L(k).name); end

    inv  = rheome.io.read.inverse(chosen);
    file = fullfile(dsdir, 'inverse_bst.mat');
    builtin('save', file, 'inv', '-v7.3');
    fprintf('rheome.import.inverse[%s]: cached %s (%s, nComp=%d) -> %s\n', ...
        name, inv.Comment, inv.Function, inv.nComponents, file);
end

% Author: Diellor Basha, 2026
