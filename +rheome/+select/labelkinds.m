function K = labelkinds(kind)
% SELECT.LABELKINDS  Every label and feature kind, with what it can CLAIM: the calibration registry.
%
%   K = rheome.select.labelkinds()            % the whole registry, one row per kind
%   r = rheome.select.labelkinds("alpha_traj_speedMS")      % one row (empty if unregistered)
%
% A label is a claim about the brain, and a claim is only as good as its calibration. Each kind here
% carries the evidence for what it measures THROUGH THIS INSTRUMENT, so a table of labels cannot hold a
% confident number that a calibration has already shown to be meaningless.
%
% ⭐ rheome.select.measure REFUSES a kind whose calibration is "unmeasurable" and WARNS on an unregistered one.
% The canonical case: an envelope or blob SIZE. Blobs injected at 63-377 mm read back 147-183 mm
% through forward and inverse with no noise at all, and 130-144 mm in the recording
% (alpha_envscale_omega.m; plant_scale_omega.m found the same for single snapshots). Writing a size
% label would put the instrument's rendering in a column named for the source.
%
% COLUMNS
%   kind          the name written to measure.kind (or a mergeable relation column)
%   relation      "measure" (note: measured where measured, never merges) or "feature_cortex[_space]"
%                 (mergeable: a parent is the sum of its children)
%   scope         "cortex" | "none"
%   unit          physical unit, or "class" / "count"
%   producer      the function that computes it
%   null          what it was judged against
%   calibration   validated | compressed | pending | unmeasurable
%   evidence      where the calibration was measured, with the number
%
% See also: rheome.select.measure, rheome.select.measures, rheome.flow.labelalpha, docs/2026-09-29-alpha-labelling-design.md
%
% Author: Diellor Basha, 2026

    r = {
    % kind                  relation                scope     unit       producer               null                        calibration     evidence
    'alpha_onFrac'          'measure'               'cortex'  'fraction' 'rheome.flow.labelalpha'      'bimodality (noise 0/16 bimodal)' 'validated'  'rheome.detect.onstate; alpha_occupancy_omega.m'
    'alpha_class'           'measure'               'cortex'  'class'    'rheome.flow.labelalpha'      'time block-shuffle; per-tile circular shift' 'pending' 'duration validated (time null, leakage-free); spatial extent awaits a leakage calibration (alpha_jointstate_omega.m)'
    'alpha_episode_tOn'     'measure'               'cortex'  's'        'rheome.detect.risefall'      'planted double sigmoid'    'validated'     'tTilePath/risefallRecoversPlantedTimes: tOn, tOff within 0.03 s'
    'alpha_episode_tOff'    'measure'               'cortex'  's'        'rheome.detect.risefall'      'planted double sigmoid'    'validated'     'as tOn'
    'alpha_episode_tauR'    'measure'               'cortex'  's'        'rheome.detect.risefall'      'planted double sigmoid'    'validated'     'within 0.04 s on planted shapes'
    'alpha_episode_tauF'    'measure'               'cortex'  's'        'rheome.detect.risefall'      'planted double sigmoid'    'validated'     'within 0.04 s on planted shapes'
    'alpha_episode_plateauS' 'measure'              'cortex'  's'        'rheome.detect.risefall'      'planted double sigmoid'    'validated'     'derived from tOn, tOff, tauR, tauF'
    'alpha_episode_r2'      'measure'               'cortex'  'fraction' 'rheome.detect.risefall'      '-'                         'validated'     'goodness of the fit itself'
    'alpha_traj_netMM'      'measure'               'cortex'  'mm'       'rheome.detect.tilepath'      'per-mode circular-shift surrogate, alpha 0.1/s' 'validated' 'alpha_inject_omega.m: still source 0 mm; movers >= 0.05 m/s at 2x typical amplitude'
    'alpha_traj_speedMS'    'measure'               'cortex'  'm/s'      'rheome.detect.tilepath'      'surrogate, alpha 0.1/s'    'compressed'    'ratio 1.00 / 0.93 at 0.05 / 0.1 m/s; 0.02 m/s reads +65%; a 40 mm amplitude swap reads 0.012-0.024 m/s'
    'alpha_traj_straight'   'measure'               'cortex'  'fraction' 'rheome.detect.tilepath'      'surrogate'                 'validated'     'movers 0.2-0.45 vs real 0.05-0.24 (alpha_inject_omega.m)'
    'alpha_traj_durS'       'measure'               'cortex'  's'        'rheome.detect.tilepath'      'surrogate'                 'validated'     'path duration within its 2 s window'
    'alpha_traj_endNode'    'measure'               'cortex'  'node'     'rheome.detect.tilepath'      'surrogate'                 'validated'     'location recovers to 43-52 mm (plant_scale_omega.m)'
    'alpha_traj_score'      'measure'               'cortex'  'score'    'rheome.detect.tilepath'      'surrogate, alpha 0.1/s'    'validated'     'held-out false rate 0.08-0.12/s (alpha_grouptrack_omega.m)'
    'occupancy'             'feature_cortex'        'cortex'  'fraction' 'rheome.flow.cortexfeatures'  'bimodality'                'validated'     'exact roll-up of a mean of a binary state (tOccupancy)'
    'energy'                'feature_cortex'        'cortex'  '(A m)^2 m^2 s' 'rheome.flow.cortexfeatures' '-'                     'validated'     'quadratic, sums exactly over nodes and tiles'
    'energy_cband'          'feature_cortex_space'  'cortex'  '(A m)^2 m^2 s' 'rheome.flow.cortexfeatures' '-'                     'validated'     'a FILTERED VIEW per spatial octave, not a size; partitions energy only globally (feature table s13)'
    'onset_slowness'        'measure'               'cortex'  's/m'      'rheome.detect.onsetspread'   'surrogate; injected fronts' 'compressed'   'injected 1 / 0.25 m/s read 0.20 / 0.65 s/m, ~5x compressed (alpha_burstspread_omega.m)'
    'envelope_scale'        'measure'               'cortex'  'mm'       '-'                    'injected blobs'            'unmeasurable'  'injected 63-377 mm read 147-183 mm noise-free (alpha_envscale_omega.m)'
    'blob_size'             'measure'               'cortex'  'mm'       'rheome.detect.blobscale'     'planted blobs'             'unmeasurable'  '~100 mm recovered whatever was planted below 10 dB (plant_scale_omega.m)'
    };
    K = cell2table(r, 'VariableNames', {'kind','relation','scope','unit','producer','null','calibration','evidence'});
    K = convertvars(K, @iscell, 'string');
    if nargin
        K = K(K.kind == string(kind), :);
    end
end

% Author: Diellor Basha, 2026
