function inv = inverse(name)
% LOAD.INVERSE  Load a cached Brainstorm inverse kernel from +data ([] if none).
%
%   inv = rheome.load.inverse(name)
%
% Returns the rheome.io.read.inverse struct written by rheome.import.inverse, or [] when no Brainstorm
% inverse has been cached for this dataset (so callers can treat the comparison as optional).
%
% See also: rheome.import.inverse, rheome.io.read.inverse, rheome.load.dataset
%
% Author: Diellor Basha, 2026

    inv = [];
    f = fullfile(rheome.load.root(), char(name), 'inverse_bst.mat');
    if exist(f, 'file')
        C = builtin('load', f);
        inv = C.inv;
    end
end

% Author: Diellor Basha, 2026
