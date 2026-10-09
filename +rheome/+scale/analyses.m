function a = analyses(which)
% SCALE.ANALYSES  The analyses the scale-up covers, and which of them the driver runs yet.
%
%   a = rheome.scale.analyses()            % all, in run order
%   a = rheome.scale.analyses("ported")    % only those with a figure-free measure in +scale
%
% The seven reports in results/reports plus the injection-recovery controls, by report:
%   resolution       resolution_article        rheome.inverse.resolution, SnrFixed 3, 300 seeds
%   bandsnr          (input to the next)        per-octave SNR vs the same-session noise run
%   bandresolution   bandresolution_article    resolution at each octave's own SNR
%   flowmap          flowmap_article           alpha_apparent_flow_omega: div/curl from sensors
%   periodicflow     periodicflow_article      alpha_periodic_flow_omega: specparam split, speed
%   grouptrack       grouptrack_article        alpha_grouptrack_omega: Viterbi paths vs swap/null
%   inject           grouptrack_article fig 5  alpha_inject_omega: injection recovery
%   vortex           report_cortical_vortex    planted spin through forward + inverse
%   sensorwavelet    report_sensor_wavelet     leadfield rows as dyadic atoms
%   coefficients     (Prognome input)          rheome.scale.coefficients: tile x scale x time envelopes;
%                                              run only when asked (a ~200 MB file), so not "ported"
%
% See also: rheome.scale.run
%
% Author: Diellor Basha, 2026

    all = ["resolution" "bandsnr" "bandresolution" "flowmap" "periodicflow" "grouptrack" ...
           "inject" "vortex" "sensorwavelet" "coefficients"];
    ported = ["resolution" "bandsnr" "bandresolution" "periodicflow" "grouptrack"];
    if nargin && which == "ported", a = ported; else, a = all; end
end

% Author: Diellor Basha, 2026
