function arr = meg(src, varargin)
% SENSORS.MEG  MEG sensor positions -- the one array in the suite with real coordinates.
%
%   arr = rheome.sensors.meg('subject01')            % a dataset cached in +data
%   arr = rheome.sensors.meg('ChannelFile', f)       % a Brainstorm channel file
%   arr = rheome.sensors.meg('Channel', C)           % an rheome.io.read.channel struct already in hand
%
% Irregularly sampled on a curved surface, which is where a graph earns its keep: there is
% no regular lattice here for an image method to run on, and no interpolation step to
% defend.
%
% ⚠ THE PRIMARY COIL IS THE POSITION. A gradiometer's .Loc carries the pickup coil AND the
% reference coil several centimetres above it. Averaging the columns puts the sensor in
% mid-air halfway to the reference, which quietly inflates both pitch and aperture. Column
% one is the sensor.
%
% ⭐ NOTHING HERE IS A SOURCE CLAIM. This is a graph on a helmet. Wavelengths measured on it
% are millimetres of SENSOR SEPARATION and speeds are metres per second ACROSS THE ARRAY.
% No leadfield is involved and none is implied.
%
% INPUTS:
%   src  a dataset name cached in +data (positional), or use the name-value forms
%   'ChannelFile'  path to a Brainstorm channel file
%   'Channel'      an rheome.io.read.channel struct
%   'Name'         display name
% OUTPUT:
%   arr  the rheome.sensors.grid array contract, with .Kind 'meg' and .Dim 2
%
% See also: rheome.sensors.eeg, rheome.sensors.graph, rheome.load.study, rheome.io.read.channel
%
% Author: Diellor Basha, 2026

    if nargin >= 1 && (ischar(src) || isstring(src)) && ...
            any(strcmpi(src, {'ChannelFile', 'Channel', 'Name'}))
        varargin = [{src}, varargin];  src = '';
    elseif nargin < 1
        src = '';
    end

    p = inputParser;
    p.addParameter('ChannelFile', '');
    p.addParameter('Channel',     []);
    p.addParameter('Name',        'MEG helmet');
    p.parse(varargin{:});
    o = p.Results;

    if ~isempty(o.Channel)
        C = o.Channel;
    elseif ~isempty(o.ChannelFile)
        C = rheome.io.read.channel(o.ChannelFile);
    elseif ~isempty(src)
        st = rheome.load.study(char(src));
        C  = st.chan;
    else
        error('sensors:meg:source', ...
            'Give a cached dataset name, a ''ChannelFile'', or a ''Channel'' struct.');
    end

    keep = find(strcmpi(C.Type, 'MEG'));
    if isempty(keep)
        error('sensors:meg:noChannels', 'No MEG channels in the channel struct.');
    end

    P = zeros(numel(keep), 3);
    for i = 1:numel(keep)
        loc = C.Channel(keep(i)).Loc;
        if isempty(loc)
            error('sensors:meg:noLoc', 'Channel %s has no Loc.', C.Name{keep(i)});
        end
        P(i,:) = loc(:,1).';                 % primary coil, NOT the coil mean
    end

    arr = sen_array(o.Name, 'meg', P, 2, C.Name(keep));
end

% Author: Diellor Basha, 2026
