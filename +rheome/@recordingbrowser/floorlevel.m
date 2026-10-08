function [L, why] = floorlevel(db, scope, band, stat, envelope)
% RECORDINGBROWSER.FLOORLEVEL  The finest level the diagonal allows, and which constraint set it.
%
%   [L, why] = rheome.recordingbrowser.floorlevel(db, 'channel', bandId, 'rms')
%
% Three constraints, whichever is coarsest:
%   * the strip needs at least one band, so never below min(bands.naturalLevel)
%   * channel rows exist only from grid.channelLevel up (group rows exist everywhere)
%   * a PER-BAND statistic needs its own band, from that band's natural level up
% A time-domain statistic (rms, crest, ...) and a spectral one (entropy, centroid) need no
% particular band -- they read whatever the level carries -- so only a bandPower or share
% column pins the level to the selected band.
%
% Author: Diellor Basha, 2026

    if nargin < 5, envelope = false; end
    L = 0;  why = 'the record';
    minN = 0;
    if ~isempty(db.bands), minN = min(db.bands.naturalLevel); end   % a preview has no bands
    if minN > L && ~envelope
        L = minN;  why = sprintf('no band is carried below level %d', minN);
    end
    % ⭐ THE ENVELOPE HAS ITS OWN, LOWER FLOOR. min and max are kept below the channel floor
    % (rheome.ingest.config ChannelEnvelope), and the envelope needs no band at all, so an envelope
    % view reaches grid.envelopeLevel while a statistic view stops at grid.channelLevel.
    chanFloor = db.channelLevel;
    if envelope, chanFloor = db.envelopeLevel; end
    if ~strcmp(scope, 'group') && chanFloor > L
        L = chanFloor;
        why = sprintf('channel %s start at level %d (%g s tiles)', i_what(envelope), L, db.grid.tExtent(L+1));
    end
    V = rheome.select.derive();
    perBand = ~envelope && ~isempty(db.bands) && ismember(string(stat), V.name(V.kind == "band"));
    if perBand && ~isempty(band)
        j = find(db.bands.j == band, 1);
        if ~isempty(j) && db.bands.naturalLevel(j) > L
            L = db.bands.naturalLevel(j);
            why = sprintf('band %d (%.3g-%.3g Hz) starts at level %d', band, db.bands.fLo(j), db.bands.fHi(j), L);
        end
    end
end

function s = i_what(envelope)
    if envelope, s = 'min/max'; else, s = 'rows'; end
end
% Author: Diellor Basha, 2026
