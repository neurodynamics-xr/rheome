function w = window(arr)
% SENSORS.WINDOW  Which wavelengths does this array admit? Computed before any signal.
%
%   w = rheome.sensors.window(arr)
%
% Three array-intrinsic facts, and nothing else, decide what an array can measure:
%
%   FACES     a winding number is summed around a face (a 2-cell); an array that admits
%             no faces cannot carry one AT ALL
%   PITCH     wavelengths below 2*pitch alias -- indistinguishable from longer ones
%   APERTURE  above the aperture the array measures a phase GRADIENT, not a wave: one
%             number, with no way to confirm the pattern is even periodic
%
% ⭐ THE WINDOW IS A PROPERTY OF THE INSTRUMENT. It is known before a signal exists, and no
% amount of SNR or cleverness moves either edge. The Utah and ECoG presets are the same
% geometry class two orders of magnitude apart in pitch, and their windows do not overlap.
%
% ⭐ THE PROBE'S TWO LIMITS ARE DIFFERENT IN KIND, which is why .Supports separates them.
% Winding is UNDEFINED on a chain -- not noisy, not hard, absent. Direction is likewise
% unavailable: a 1D array measures only k along its axis, so a wave crossing at angle theta
% reads speed c/cos(theta), biased HIGH and uncorrectable from the probe alone. But
% omega(lambda) and speed ARE available in that projection, so they stay true.
%
% ⚠ EVERY LENGTH IS SENSOR SEPARATION. LambdaMin and LambdaMax are metres between sensors.
% They are not source extents and imply nothing about what produced the signal.
%
% INPUTS:
%   arr  a rheome.sensors.* array struct
% OUTPUT:
%   w  .LambdaMin .LambdaMax .Octaves .Pitch .Aperture .Dim .HasFaces .Supports .Name
%
% See also: rheome.sensors.calibrate, rheome.sensors.graph, rheome.demos.sensor_limits
%
% Author: Diellor Basha, 2026

    hasFaces = arr.Dim >= 2;
    w = struct();
    w.Name       = arr.Name;
    w.Dim        = arr.Dim;
    w.Pitch      = arr.Pitch;
    w.Aperture   = arr.Aperture;
    w.LambdaMin  = 2 * arr.Pitch;
    w.LambdaMax  = arr.Aperture;
    w.Octaves    = log2(w.LambdaMax / max(w.LambdaMin, eps));
    w.HasFaces   = hasFaces;
    w.Supports   = struct('dispersion', true, ...
                          'speed',      true, ...
                          'direction',  hasFaces, ...
                          'winding',    hasFaces, ...
                          'chirality',  hasFaces);
end

% Author: Diellor Basha, 2026
