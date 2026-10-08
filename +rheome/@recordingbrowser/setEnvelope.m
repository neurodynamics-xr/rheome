function setEnvelope(app, on)
% SETENVELOPE  Draw the min/max band on the top axis instead of the statistic.
%
% ⭐ WHY THIS IS THE RIGHT LEVEL-OF-DETAIL VIEW OF A TRACE. min and max merge by extremum,
% so the pair stored on a tile is a bound that contains every sample in that tile, at every
% level of the pyramid. The band drawn from them is therefore exact: it can be too wide, it
% can never be too narrow, and no transient can hide between the drawn columns. A viewer
% that decimates by sampling has no such guarantee and will step over a spike.
%
% ⚠ IT NEEDS THE ENVELOPE EXEMPTION IN THE STORE. Below grid.channelLevel a store keeps the
% per-channel min and max only when it was built with ChannelEnvelope (the default); an
% older store stops at the channel floor and the envelope goes no finer than that. Group
% scope always reaches level 0, because a node carries rows everywhere.
%
% Author: Diellor Basha, 2026

    app.Envelope = logical(on);
    app.refresh();
end
% Author: Diellor Basha, 2026
