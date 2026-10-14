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
%   flowmap          flowmap_article           div/curl from sensors
%   periodicflow     periodicflow_article      specparam split, speed
%   grouptrack       grouptrack_article        Viterbi paths vs swap/null
%   fieldsmooth      MS1 Fig. 12               per-vertex vs graph-wavelet band-limited J, div, curl
%   eventsensors     MS1 Figs. 8, 10           the tracked event's samples and channels (after grouptrack)
%   bandperiodic     MS1 Table 5               per-band periodic fraction, oscillation SNR, floors, IAF, halves
%   helmholtzbands   MS1 Fig. 4A-B             planted Helmholtz bands recovered on the subject's cortex
%   ownregion        MS1 rule 5                own-region fraction of div/curl per Desikan-Killiany region
%   geometry         MS1 Figs. 2, 16, 20       atom-tile overlap, gauge singularities, roll-up exactness
%   fusion           MS1 Section 8.2           fused kernels against reconstruct-then-differentiate
%   inject           grouptrack_article fig 5  injection recovery
%   inject           MS1 G7 (Fig. 8)           per participant: injection recovery
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
%   aperiodic        MS1 G10 (Fig. 14)         a moving 1/f change planted at 30-100 mm through own MEG: on-patch
%                                              false-path threshold and the periodic split's selectivity;
%                                              run only when asked (nsp cf-aperiodic), not "ported"
%   slowosc          MS1 G15, G16              sleep slow oscillations (EEG, rheome.scale.importeeg): direction vs A->P
%                                              and speed per event; framework, bst_of, phasereg, sensor latency;
%                                              time-reversed, phase-randomised and lagged-source nulls -- run only when asked
%                                              (nsp cf-slowosc), not "ported"
%   slowoscplant     SLEEP-1 TP1 gate V        planted A->P / P->A slow-oscillation movers at 1, 2, 5 m/s in the
%                                              EEG head, on the participant's own N2/N3 background at its median
%                                              event p2p and at the detection threshold; same estimators and nulls;
%                                              run only when asked (nsp cf-slowosc-plant), not "ported"
%   correspondence   MS1 G12 (section 12.2)    DK through the sphere vs raw coordinates; frame vs the shared meridian
%   atlasevents      MS1 G12 (Fig. 10)         alpha events on the group tiles, depth <= 3, vs phase-randomised surrogates
%   gaugetensor      MS1 G12 (section 12.3)    alpha orientation tensor per depth-3 group tile in the shared gauge
%   connectome       MS1 G12 (Fig. 3C)         Destrieux -> DK vs dyadic roll-up; connectome-wavelet widths, gamma sweep
%                    -- the four run only when asked (nsp cf-atlas), not "ported"; rheome.scale.measure_atlas
%   atlasone         MS1 3.4 (F8, F12)          one person's alpha atlas, depths 0-6 on one tree, on the paper's kernel
%   xlevel           MS1 2.4 (F23)              a planted mover and a real alpha event tracked across depths 0-6
%                    -- both rheome.scale.measure_atlas; need a template cortex; run only when asked
%                    (nsp rheome-ms1-panels), not "ported"
%   scaleramp        MS1 2.4                    a peak whose size ramps 10 -> 40 mm, read by rheome.detect.ridges on a
%                                              sphere and the left cortex (rheome.scale.measure_scaleramp); run only
%                                              when asked (nsp rheome-ms1-panels), not "ported"
%   multimodal       MS1 multimodal      MEG band power, PET SUVR, fibre degree through one wavelet bank on
%                                              one dyadic ladder, exact roll-up; run only when asked, not "ported"
%   modesfields      resting-state pattern atlas, stage 1: complex band tile series + nulls N1, N4, N5
%   modesfit         stage 3: back-fit of the frozen centroids, features, N2, N3, N5 (RHEOME_MODES_CENTROIDS)
%                    -- both run only when asked (nsp rheome-modes-fields / rheome-modes-fit), not "ported"
%   jointfields      H4 ATLAS-1: joint graph-time wavelet field + null F, N4, N5 + scale transfer T(s) (opt-in)
%   jointcanary      H4 ATLAS-1 J5: the four plants on null F, no real-data series (opt-in)
%   upstream         AGE-1 TP2, WAVE-OR-LEAKAGE TP3   rheome.scale.measure_upstream: share of slow-band travel up
%                                              the principal gradient (upstreamshare v1) -- U, U_surrogate, U_total,
%                                              U_of and the aperiodic covariates; run only when asked (nsp
%                                              rheome-upstream), not "ported"
%   lambdaedge       TP2/TP3 gate S            rheome.scale.measure_lambdaedge: readout v2's short-lambda edge on
%                                              the head alone, noiseless, no recording read; run only when asked
%                                              (nsp rheome-lambdaedge), not "ported"
%   panelexamples    MS1 Figs. 4-7 panel A      one example each (plant + band maps, PSF, injected mover, aperiodic
%                                              patch + spectra) by the figure's own measure; run only when asked
%                                              (nsp rheome-ms1-examples), not "ported"
%   ambient          MS1 1.4 panel             a planted tilted field, sphere + own cortex: ambient vs tangential
%                                              divergence, the 2H(J.n) term against its closed form (no instrument)
%   cptrack          MS1 2.3 panel             a drifting rotor + a waxing/waning source tracked by detect.track:
%                                              core path, curl and div against the truth (no instrument)
%                    -- both run only when asked (nsp rheome-ms1-panels), not "ported"
%   events           MS1 SPEC1 §8 B            the tracker's detections (depth 8, 18.75 fps) as an event table: risefall
%                                              boundaries, peak frequency, cycles; onset/offset-locked TF vs random
%                                              onsets; intervals and autocorrelation vs surrogate and Poisson
%   eventplants      MS1 SPEC1 §8 A            planted episodes of known timing in own rest, boundaries read back
%                    -- both run only when asked (nsp cf-events), not "ported"
%   freqtrack        MS1 2.5                   planted chirps in own rest tracked across frequency (jointpath), the
%                                              surrogate's false and cross-band false rates; run only when asked
%                                              (nsp cf-chirps), not "ported"
%   coefficients     (Prognome input)          rheome.scale.coefficients: tile x scale x time envelopes;
%                                              run only when asked (a ~200 MB file), so not "ported"
%
% See also: rheome.scale.run
%
% Author: Diellor Basha, 2026

    all = ["resolution" "bandsnr" "bandresolution" "flowmap" "periodicflow" "grouptrack" ...
           "fieldsmooth" "eventsensors" "inject" "vortex" "sensorwavelet" ...
           "plantfloors" "movingvortex" "composition" "sizeruler" "vortexscale" "rotation" "detection" "diracangles" ...
           "bandperiodic" "helmholtzbands" "ownregion" "geometry" "fusion" "patterns" "patternnulls" "catalogue" "catalognulls" "trackfactorial" "aperiodic" "slowosc" "slowoscplant" "correspondence" "atlasevents" "gaugetensor" "connectome" "atlasone" "xlevel" "scaleramp" "multimodal" "modesfields" "modesfit" "jointfields" "jointcanary" "upstream" "lambdaedge" "panelexamples" "ambient" "cptrack" "events" "eventplants" "freqtrack" "coefficients"];
    ported = ["resolution" "bandsnr" "bandresolution" "periodicflow" "grouptrack" "fieldsmooth" "eventsensors"];
    if nargin && which == "ported", a = ported; else, a = all; end
end

% Author: Diellor Basha, 2026
