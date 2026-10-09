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
%   fieldsmooth      MS1 Fig. 12               per-vertex vs graph-wavelet band-limited J, div, curl
%   eventsensors     MS1 Figs. 8, 10           the tracked event's samples and channels (after grouptrack)
%   bandperiodic     MS1 Table 5               per-band periodic fraction, oscillation SNR, floors, IAF, halves
%   helmholtzbands   MS1 Fig. 4A-B             planted Helmholtz bands recovered on the subject's cortex
%   ownregion        MS1 rule 5                own-region fraction of div/curl per Desikan-Killiany region
%   geometry         MS1 Figs. 2, 16, 20       atom-tile overlap, gauge singularities, roll-up exactness
%   fusion           MS1 Section 8.2           fused kernels against reconstruct-then-differentiate
%   inject           grouptrack_article fig 5  alpha_inject_omega: injection recovery
%   inject           MS1 G7 (Fig. 8)           alpha_inject_omega per participant: injection recovery
%   trackfactorial   MS1 G8 (Fig. 9)           the 0.52 diagnosis: sphere movers, 2^5 factorial
%                    -- both run only when asked (nsp cf-track), not "ported"
%   vortex           report_cortical_vortex    planted spin through forward + inverse
%   sensorwavelet    report_sensor_wavelet     leadfield rows as dyadic atoms
%   plantfloors      MS1 G1, G9 (Table 3)      planted source/vortex bands x SNR through own MEG: floors
%   movingvortex     MS1 G11 (Fig. 11)         moving vortex core error per readout sigma x SNR
%   composition sizeruler vortexscale rotation detection diracangles
%                    MS1 G2 (section 9.4)      what the noise floor and the inverse manufacture
%                    -- the eight MS1 plant analyses run only when asked (nsp cf-plant-floors), not "ported"
%   patterns         MS1 G5 (Fig. 15C)        catalogue detectors per frame vs 200 phase-randomised surrogates
%                                              and the empty room, alpha and IAF +- 2 Hz, + theta-gamma PAC
%   patternnulls     MS1 G13 (section 6.3)     speed-sweep peakedness, Dirac dispersion ratio, r(curl v, curl J)
%                    -- both run only when asked (nsp cf-patterns), not "ported"
%   catalogue        MS1 G6, G3 rule 11, G16   the 19 catalogue plants x 3 centres per hemisphere through own gain,
%                                              whitened MNE and sensor noise (rest, empty room) at Inf and 10 dB;
%                                              framework vs bst_opticalflow (HornSchunck sweep) vs phase regression
%   catalognulls     MS1 G16                   false-propagation nulls: coherent and phase-lagged generator pairs
%                    -- both run only when asked (nsp cf-plants), not "ported"
%   coefficients     (Prognome input)          rheome.scale.coefficients: tile x scale x time envelopes;
%                                              run only when asked (a ~200 MB file), so not "ported"
%
% See also: rheome.scale.run
%
% Author: Diellor Basha, 2026

    all = ["resolution" "bandsnr" "bandresolution" "flowmap" "periodicflow" "grouptrack" ...
           "fieldsmooth" "eventsensors" "inject" "vortex" "sensorwavelet" ...
           "plantfloors" "movingvortex" "composition" "sizeruler" "vortexscale" "rotation" "detection" "diracangles" ...
           "bandperiodic" "helmholtzbands" "ownregion" "geometry" "fusion" "patterns" "patternnulls" "catalogue" "catalognulls" "trackfactorial" "coefficients"];
    ported = ["resolution" "bandsnr" "bandresolution" "periodicflow" "grouptrack" "fieldsmooth" "eventsensors"];
    if nargin && which == "ported", a = ported; else, a = all; end
end

% Author: Diellor Basha, 2026
