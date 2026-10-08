function arr = positions(name, P, varargin)
% SENSORS.POSITIONS  An array from bare coordinates. The general case.
%
%   arr = rheome.sensors.positions('MEG helmet', P)          % P [N x 3] or [N x 2], metres
%   arr = rheome.sensors.positions(name, P, 'Dim', 2, 'Kind', 'meg', 'Labels', lab)
%
% rheome.sensors.linear, .grid, .meg and .eeg each KNOW where their coordinates come from. This
% one does not: it takes them as given and measures the same pitch and aperture from them,
% so everything downstream -- rheome.sensors.graph, .modes, .calibrate, .window -- accepts the
% result unchanged.
%
% ⭐ FOR COORDINATES THAT ALREADY EXIST. Reach for this when the positions are a record
% rather than a construction: a projected plane a figure was drawn in, an array read from
% somewhere with no constructor, a geometry recovered from a published result. It is what
% lets an array be analysed on exactly the coordinates it was previously DISPLAYED on.
%
% ⚠ TWO COLUMNS MEANS A PLANE, AND A PLANE MAY BE A PROJECTION. Padding z with zeros is
% correct for a flat array and is an ORTHOGRAPHIC APPROXIMATION for a curved one: on a
% helmet, separations near the rim come back short. That is a real distortion and it is
% the caller's to declare -- but it is the SAME distortion any analysis of those projected
% coordinates already carries, so it is consistent rather than new.
%
% Dim is inferred from the RANK of the centred coordinates, not from axis alignment: a
% probe inserted obliquely has spread on all three axes and is still a chain.
%
% INPUTS:
%   name     char, display name
%   P        [N x 3] or [N x 2] coordinates in METRES
%   'Dim'    1 or 2; default inferred
%   'Kind'   default 'custom'
%   'Labels' cellstr; default S001..SNNN
% OUTPUT:
%   arr  .Name .Kind .Pos .Dim .Pitch .Aperture .Labels .nCh
%
% See also: rheome.sensors.linear, rheome.sensors.grid, rheome.sensors.meg, rheome.sensors.eeg, rheome.sensors.graph
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('Dim',    []);
    p.addParameter('Kind',   'custom');
    p.addParameter('Labels', {});
    p.parse(varargin{:});
    o = p.Results;

    if ~isnumeric(P) || ~ismatrix(P) || ~any(size(P,2) == [2 3])
        error('sensors:positions:pos', ...
            'P must be [N x 3] or [N x 2] in metres, got %s.', mat2str(size(P)));
    end
    P = double(P);
    if size(P,2) == 2, P(:,3) = 0; end

    dim = o.Dim;
    if isempty(dim), dim = i_dim(P); end

    labels = o.Labels;
    if isempty(labels)
        labels = compose('S%03d', (1:size(P,1)).');
        labels = cellstr(labels(:).');
    end

    arr = sen_array(name, o.Kind, P, dim, labels);
end

function d = i_dim(P)
% Rank of the centred coordinates, by singular value. Two directions carrying real spread
% is a surface; one is a chain.
    Pc = P - mean(P, 1);
    s  = svd(Pc, 'econ');
    if numel(s) < 2 || s(1) <= 0, d = 1; return; end
    d = 1 + double(s(2) > 1e-6 * s(1));
end

% Author: Diellor Basha, 2026
