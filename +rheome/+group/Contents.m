% +GROUP  Do the patterns differ between groups or with a covariate, tile by tile, beyond chance and beyond the map's own smoothness?
%
% Tile-wise inference on the shared group tiles, one cohort at a time, with family-wise error control that
% holds at any spatial dependence between tiles. Calibrated on null labels and planted effects in the real
% cohort tile maps before any real outcome is read.
%
%   glm       - tile-wise GLM: one effect's t, Freedman-Lane permutation p and max-|t| FWER; repeated recordings
%               by the random-intercept model (summary statistics)
%   spin      - spin test of an effect map against reference maps, FWER over the references
%   calibrate - null labels and planted group / covariate effects on the cohort tile maps: FWER, out-of-plant
%               false positives, power, spin false rate (the nsp group job)
%
% Author: Diellor Basha, 2026
