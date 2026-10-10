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
    % a 1e6-streamline Brainstorm file holds 2.4 GB of Points, its ends 48 MB: read only the two end slices
    % from a v7.3 (HDF5) file; an older MAT-file cannot be sliced and is loaded whole, as before
    fid = fopen(f, 'r');  hdr = fread(fid, [1 20], '*char');  fclose(fid);
    if startsWith(hdr, 'MATLAB 7.3')
        m = matfile(f);  n = size(m, 'Points', 2);
        endPts = double(cat(2, m.Points(:, 1, :), m.Points(:, n, :)));   % first + last point -> [nFib x 2 x 3]
    else
        Fm = load(f, 'Points');  endPts = double(Fm.Points(:, [1 end], :));
    end
end

% Author: Diellor Basha, 2026
