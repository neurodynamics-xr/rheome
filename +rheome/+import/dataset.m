function dataset(name, studyDir, cortexFile, dataName, tau, K)
% IMPORT.DATASET  Import a full dataset (surface + study + Dirac eigenbasis) in one call.
%
%   rheome.import.dataset(name, studyDir, cortexFile, dataName)             % tau=0.5, K=400
%   rheome.import.dataset(name, studyDir, cortexFile, dataName, tau, K)
%   rheome.import.dataset(name, '',       cortexFile, '')                   % geometry-only (no study)
%
% Reads every Brainstorm structure the demos use, computes the LBO operators and the
% relative-Dirac eigenbasis, and caches them under +data/<name>/. Run this ONCE;
% thereafter the demos load instantly via rheome.load.*. Pass an empty studyDir/dataName to
% import a geometry-only dataset (surface + Dirac, no leadfield/recording).
%
% See also: rheome.import.surface, rheome.import.study, rheome.import.dirac, rheome.load.dataset, rheome.load.list
%
% Author: Diellor Basha, 2026

    if nargin < 5 || isempty(tau), tau = 0.5; end
    if nargin < 6 || isempty(K),   K   = 400; end

    fprintf('=== rheome.import.dataset[%s] ===\n', name);
    rheome.import.surface(name, cortexFile);
    if nargin >= 2 && ~isempty(studyDir) && nargin >= 4 && ~isempty(dataName)
        rheome.import.study(name, studyDir, dataName);
    else
        fprintf('  (geometry-only: no study imported)\n');
    end
    rheome.import.dirac(name, tau, K);
    rheome.import.bases(name, K);                        % per-hemisphere LBO + connection operators & eigenbases
    fprintf('=== dataset ''%s'' ready (tau=%.2f, K=%d). Load with rheome.load.dataset(''%s''). ===\n', name, tau, K, name);
end

% Author: Diellor Basha, 2026
