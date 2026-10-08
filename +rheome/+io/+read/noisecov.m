function N = noisecov(NoiseCovFile)
% IO.READ.NOISECOV  Read a Brainstorm noise-covariance file.
%
%   N = rheome.io.read.noisecov(NoiseCovFile)
%
% OUTPUT (struct N):
%   .NoiseCov     [nCh x nCh] noise covariance
%   .FourthMoment [nCh x nCh] (for Ledoit-Wolf shrinkage), else []
%   .nSamples     [nCh x nCh] sample counts, else []
%   .nCh          count
%   .Comment
%
% Author: Diellor Basha, 2026

    raw = load(NoiseCovFile);
    if ~isfield(raw, 'NoiseCov') || isempty(raw.NoiseCov)
        error('io:read:noisecov:noCov', 'No NoiseCov matrix in %s', NoiseCovFile);
    end
    N = struct();
    N.NoiseCov     = double(raw.NoiseCov);
    N.nCh          = size(N.NoiseCov, 1);
    N.FourthMoment = i_get(raw, 'FourthMoment', []);
    N.nSamples     = i_get(raw, 'nSamples', []);
    N.Comment      = i_get(raw, 'Comment', '');
end

function v = i_get(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

% Author: Diellor Basha, 2026
