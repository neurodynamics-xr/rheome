function disp(obj)
% DISP  One-screen summary of the page in hand.
% Author: Diellor Basha, 2026

    fprintf('  rheome.flowpage\n');
    fprintf('    Band        [%g %g] Hz -> %d filters, %.2f..%.2f Hz\n', ...
        obj.Band(1), obj.Band(2), obj.NumFilters, ...
        min(obj.CenterFrequencies), max(obj.CenterFrequencies));
    fprintf('    Page        %d of %d\n', obj.PageIndex, obj.Bundle.pager.NumPages);
    fprintf('    Rate        %.1f Hz (decim %d, %.1f samples/cycle at the top filter)\n', ...
        obj.Rate, obj.Decim, obj.Rate / max(obj.CenterFrequencies));
    fprintf('    Frames      %d (%d core, %.2f..%.2f s)\n', obj.NumFrames, sum(obj.Core), ...
        obj.Time(1), obj.Time(end));
    fprintf('    Space       %d modes -> %d vertices | %d spatial scales\n', ...
        numel(obj.Lambda), obj.NumVertices, obj.NumScales);
end

% Author: Diellor Basha, 2026
