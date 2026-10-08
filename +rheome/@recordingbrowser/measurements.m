function M = measurements(app, opts)
% MEASUREMENTS  The measurements on this store, optionally joined to their tiles' statistics.
%
%   M = app.measurements()
%   M = app.measurements(Kind="amplitude", InView=true, Context=["rms","share"])
%
% InView keeps only the measurements whose tile overlaps the current window and whose unit
% is on screen. Context joins each row to the moments of its own tile and, with Ancestors,
% of the node above it -- the state the measurement was made in.
%
% See also: rheome.recordingbrowser/measure, rheome.select.measures, rheome.select.context
%
% Author: Diellor Basha, 2026

    arguments
        app
        opts.Kind      string = string.empty
        opts.InView    (1,1) logical = false
        opts.Context   string = string.empty
        opts.Ancestors double = 0
    end
    args = {};
    if ~isempty(opts.Kind), args = [args, {'Kind', opts.Kind}]; end
    if opts.InView
        args = [args, {'Window', app.Window}];
        if strcmp(app.View, 'channels') && ~isempty(app.Stack.units)
            args = [args, {'Units', app.Stack.units}];
        end
    end
    M = rheome.select.measures(app.Db, args{:});
    if ~isempty(opts.Context) && ~isempty(M)
        M = rheome.select.context(app.Db, M, Stats=opts.Context, Ancestors=opts.Ancestors);
    end
end
% Author: Diellor Basha, 2026
