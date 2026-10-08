function disp(obj)
% DISP  One-screen summary: the full axis, the paging, and the margin cost.
%
% Author: Diellor Basha, 2026

    fprintf('  rheome.pagedrecording\n');
    fprintf('    File          %s\n', obj.File);
    fprintf('    Axis          %d samples @ %g Hz = %.1f s\n', ...
        obj.NumSamples, obj.SamplingFrequency, obj.Duration);
    fprintf('    Channels      %d of %d selected\n', obj.NumChannels, numel(obj.ChannelFlag));
    fprintf('    PageLength    %d samples (%.2f s)\n', ...
        obj.PageLength, obj.PageLength / obj.SamplingFrequency);
    fprintf('    Overlap       %d samples (%.3f s each side)\n', ...
        obj.Overlap, obj.Overlap / obj.SamplingFrequency);
    fprintf('    NumPages      %d\n', obj.NumPages);
    fprintf('    Precision     %s\n', obj.Precision);
    span = obj.PageLength + 2*obj.Overlap;
    fprintf('    Margin cost   read %d per %d core = %.0f%% overhead\n', ...
        span, obj.PageLength, 100 * (span - obj.PageLength) / obj.PageLength);
end

% Author: Diellor Basha, 2026
