function arr = eeg(varargin)
% SENSORS.EEG  A scalp-cap array: sparse, quasi-uniform sampling on a curved surface.
%
%   arr = rheome.sensors.eeg()                                  % 64 channels, 90 mm, 120 deg cap
%   arr = rheome.sensors.eeg('NumChannels', 128)
%   arr = rheome.sensors.eeg('ChannelFile', f)                  % real positions from a channel file
%
% ⚠ THE BUILT-IN GEOMETRY IS AN IDEALISATION AND NOT A DIGITISED MONTAGE. It is a
% spherical-Fibonacci sampling of a cap -- quasi-uniform by construction, with generic
% labels. Nothing in this suite depends on millimetre-accurate electrode placement; what
% this array contributes is SPARSE SAMPLING ON A CURVED CAP, which the idealisation
% delivers exactly. Pass 'ChannelFile' when real positions matter.
%
% Quasi-uniformity is the point rather than a convenience: pitch is a single number, and it
% only means anything if the sampling has no clusters or holes for it to average over.
%
% INPUTS (name-value):
%   'NumChannels'  (default 64)
%   'Radius'       sphere radius, metres (default 0.09)
%   'CapAngle'     cap half-extent from the vertex, degrees (default 120)
%   'ChannelFile'  a Brainstorm channel file; overrides the idealisation
%   'Channel'      an rheome.io.read.channel struct already in hand; overrides the idealisation
%   'Name'         display name
% OUTPUT:
%   arr  the rheome.sensors.grid array contract, with .Kind 'eeg' and .Dim 2
%
% See also: rheome.sensors.meg, rheome.sensors.grid, rheome.sensors.window
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('NumChannels', 64);
    p.addParameter('Radius',      0.09);
    p.addParameter('CapAngle',    120);
    p.addParameter('ChannelFile', '');
    p.addParameter('Channel',     []);
    p.addParameter('Name',        'EEG cap');
    p.parse(varargin{:});
    o = p.Results;

    if ~isempty(o.Channel) || ~isempty(o.ChannelFile)
        if ~isempty(o.Channel)
            C = o.Channel;
        else
            C = rheome.io.read.channel(o.ChannelFile);
        end
        % ⚠ EXACT MATCH, NOT A PREFIX. strncmpi(...,'EEG',3) also matches 'EEG REF' and any
        % other 'EEG '-prefixed type. The same three-character prefix match in rheome.sensors.meg
        % silently admitted 26 CTF reference gradiometers -- which sit in a rack behind the
        % head, not over it -- and overstated that array's aperture by 40%.
        keep = find(strcmpi(C.Type, 'EEG'));
        if isempty(keep)
            error('sensors:eeg:noChannels', 'No EEG channels in the channel struct.');
        end
        P = zeros(numel(keep), 3);
        for i = 1:numel(keep)
            P(i,:) = C.Channel(keep(i)).Loc(:,1).';
        end
        arr = sen_array(o.Name, 'eeg', P, 2, C.Name(keep));
        return;
    end

    n = double(o.NumChannels);
    if ~isscalar(n) || n < 4 || mod(n, 1) ~= 0
        error('sensors:eeg:numChannels', 'NumChannels must be an integer >= 4, got %s.', mat2str(n));
    end
    R   = double(o.Radius);
    cap = double(o.CapAngle) * pi/180;

    % Spherical Fibonacci over the cap: uniform in z gives uniform AREA on a sphere, and
    % the golden-angle azimuth avoids the rings that a naive lattice would produce.
    i    = (0:n-1).' + 0.5;
    zc   = 1 - (i/n) * (1 - cos(cap));
    rc   = sqrt(max(1 - zc.^2, 0));
    gold = pi * (3 - sqrt(5));
    th   = gold * i;
    P    = R * [rc.*cos(th), rc.*sin(th), zc];

    L = arrayfun(@(k) sprintf('E%03d', k), 1:n, 'UniformOutput', false);
    arr = sen_array(o.Name, 'eeg', P, 2, L);
end

% Author: Diellor Basha, 2026
