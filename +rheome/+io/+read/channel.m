function C = channel(ChannelFile)
% IO.READ.CHANNEL  Read a Brainstorm channel file (sensor definitions).
%
%   C = rheome.io.read.channel(ChannelFile)
%
% OUTPUT (struct C):
%   .Channel  [1 x nCh] struct array (Name, Type, Loc, ...), as stored
%   .Type     {1 x nCh} cell of channel-type strings (convenience)
%   .Name     {1 x nCh} cell of channel names (convenience)
%   .nCh      count
%   .Comment
%
% Author: Diellor Basha, 2026

    raw = load(ChannelFile);
    if ~isfield(raw, 'Channel') || isempty(raw.Channel)
        error('io:read:channel:noChannel', 'No Channel struct in %s', ChannelFile);
    end
    C = struct();
    C.Channel = raw.Channel;
    C.nCh     = numel(raw.Channel);
    C.Type    = {raw.Channel.Type};
    C.Name    = {raw.Channel.Name};
    C.Comment = '';
    if isfield(raw, 'Comment'), C.Comment = raw.Comment; end
end

% Author: Diellor Basha, 2026
