function endPts = fibers(name)
% IO.READ.FIBERS  Fiber endpoints from a cached Brainstorm Fibers .mat.
%
%   endPts = rheome.io.read.fibers(name)     % name = dataset (reads +data/<name>/fibers.mat) OR a .mat path
%
% Reads a Brainstorm Fibers surface .mat (its .Points field, [nFib x nPtsPerFiber x 3]) and
% returns just the two ENDPOINTS of each streamline, [nFib x 2 x 3], in the surface's
% coordinates -- the input to rheome.operators.connectome. The fibers must already be on / registered
% to the target cortex (the template->subject registration is a Brainstorm step, done once when
% the cache is written); this reader does no registration.
%
% See also: rheome.operators.connectome, rheome.load.dataset
%
% Author: Diellor Basha, 2026

    if exist(name, 'file') == 2
        f = name;
    else
        f = fullfile(rheome.load.root(), char(name), 'fibers.mat');
    end
    if ~exist(f, 'file')
        error('io:read:fibers:missing', ...
            'No fibers cached for ''%s'' (expected %s). Register + export a Fibers .mat first.', name, f);
    end
    Fm = load(f, 'Points');
    P  = double(Fm.Points);
    endPts = P(:, [1 end], :);      % first + last point of each streamline -> [nFib x 2 x 3]
end

% Author: Diellor Basha, 2026
