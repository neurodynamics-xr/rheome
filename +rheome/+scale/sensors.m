function S = sensors(name)
% SCALE.SENSORS  The channel selection, leadfield, noise covariance and kernel every measure shares.
%
%   S = rheome.scale.sensors(name)
%
% The block that opens resolution_scales_omega.m, band_resolution_omega.m and the flow scripts,
% verbatim in effect: good MEG channels only, the unconstrained leadfield on those rows, the
% Brainstorm noise covariance, and the plain minimum norm with SnrFixed = 3, amplitude measure.
%
% Returns S with .G .ncm .chT .Loc .iSel .nV .Res (rheome.inverse.mne output) .B (rheome.load.bases) .Sf.
%
% See also: rheome.inverse.mne, rheome.inverse.resolution, rheome.scale.measure_resolution
%
% Author: Diellor Basha, 2026

    S.B  = rheome.load.bases(name);
    S.Sf = rheome.load.surface(name);
    st = rheome.load.study(name);
    isMEG = strcmpi(st.chan.Type, 'MEG');
    S.iSel = find(isMEG(:) & (st.rec.ChannelFlag(:) == 1));
    S.G   = double(st.hm.Gain(S.iSel,:));
    S.ncm = struct('NoiseCov', double(st.ncov.NoiseCov(S.iSel,S.iSel)));
    S.chT = st.chan.Type(S.iSel);
    S.Loc = cell2mat(arrayfun(@(c) c.Loc(:,1), st.chan.Channel(S.iSel), 'UniformOutput', false));
    S.nV  = S.Sf.nV;
    S.Res = rheome.inverse.mne(S.G, S.ncm, struct('ChannelTypes', {S.chT}, 'InverseMeasure', 'amplitude', ...
                                           'nVert', S.nV, 'SnrFixed', 3));
end

% Author: Diellor Basha, 2026
