function tf = ismeg(types)
% SENSORS.ISMEG  True for the Brainstorm channel types that are MEG measurement channels.
%
%   tf = rheome.sensors.ismeg(chan.Type)
%
% 'MEG' (CTF, the OMEGA / PREVENT-AD cohorts) and 'MEG MAG' / 'MEG GRAD' (MEGIN / Elekta VectorView,
% Cam-CAN). Never 'MEG REF' (CTF reference coils): strncmpi(Type, 'MEG', 3) would take it.
% Keep the types as they are downstream: rheome.inverse.mne whitens per type, so magnetometers (T)
% and planar gradiometers (T/m) must not be relabelled to one 'MEG' type.
%
% Author: Diellor Basha, 2026

    tf = ismember(upper(string(types)), ["MEG" "MEG MAG" "MEG GRAD"]);
end
