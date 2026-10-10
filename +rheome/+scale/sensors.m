function S = sensors(name, opts)
% SCALE.SENSORS  The channel selection, leadfield, noise covariance and kernel every measure shares.
%
%   S = rheome.scale.sensors(name)
%   S = rheome.scale.sensors(name, Modality="EEG")     % sleep EEG (rheome.scale.importeeg): average reference
%
% The block that opens the single-subject resolution, band and flow analyses,
% verbatim in effect: good MEG channels only, the unconstrained leadfield on those rows, the
% Brainstorm noise covariance, and the plain minimum norm with SnrFixed = 3, amplitude measure.
%
% Returns S with .G .ncm .chT .Loc .iSel .nV .Res (rheome.inverse.mne output) .B (rheome.load.bases) .Sf
% .Ref ([nSel x nSel] the reference operator applied to the gain; identity for MEG).
%
% ⚠ EEG IS AVERAGE-REFERENCED: .G is Ref*G with Ref = I - 1/n over the good EEG channels, so the kernel
% expects average-referenced data -- apply S.Ref to the recording's rows S.iSel before the kernel.
%
% See also: rheome.inverse.mne, rheome.inverse.resolution, rheome.scale.measure_resolution
%
% Author: Diellor Basha, 2026

    arguments
        name (1,:) char
        opts.Modality (1,1) string {mustBeMember(opts.Modality, ["MEG" "EEG"])} = "MEG"
    end
    S.B  = rheome.load.bases(name);
    S.Sf = rheome.load.surface(name);
    st = rheome.load.study(name);
    isMod = strcmpi(st.chan.Type, char(opts.Modality));
    S.iSel = find(isMod(:) & (st.rec.ChannelFlag(:) == 1));
    n = numel(S.iSel);  S.Ref = eye(n);  if opts.Modality == "EEG", S.Ref = S.Ref - 1/n; end
    S.G   = S.Ref * double(st.hm.Gain(S.iSel,:));
    S.ncm = struct('NoiseCov', S.Ref * double(st.ncov.NoiseCov(S.iSel,S.iSel)) * S.Ref');
    S.chT = st.chan.Type(S.iSel);
    S.Loc = cell2mat(arrayfun(@(c) c.Loc(:,1), st.chan.Channel(S.iSel), 'UniformOutput', false));
    S.nV  = S.Sf.nV;
    S.Res = rheome.inverse.mne(S.G, S.ncm, struct('ChannelTypes', {S.chT}, 'InverseMeasure', 'amplitude', ...
                                           'nVert', S.nV, 'SnrFixed', 3));
end

% Author: Diellor Basha, 2026
