function p = start(type, d, opts)
% PIPELINE.START  Begin a pipeline at a field of a given type on a given domain.
%
%   p = rheome.pipeline.start("sensorScalar", dSensors, Frames=rheome.domain.index(nT, Name="time"))
%
% The pipeline is a VALUE: a struct carrying the state it is currently in (field type,
% spatial domain, frame domain) and the steps taken so far. Nothing is bound and nothing
% runs until rheome.pipeline.run, so a plan can be built, printed, compared and stored on its own.
%
% See also: rheome.pipeline.then, rheome.pipeline.window, rheome.pipeline.run, rheome.pipeline.describe
%
% Author: Diellor Basha, 2026

    arguments
        type
        d (1,1) struct
        opts.Frames = []
        opts.Name (1,1) string = ""
    end
    f = rheome.fieldtype.describe(type);                      % errors on an unknown type, here
    if ~ismember(string(d.kind), string(f.domain_kind))
        error('pipeline:start:domain', '%s cannot live on a %s domain.', f.id, d.kind);
    end
    fr = opts.Frames;
    if isempty(fr), fr = rheome.domain.index(1, Name="frames"); end
    p = struct('name', char(opts.Name), 'type', string(f.id), 'domain', d, 'frames', fr, ...
               'in', struct('type', string(f.id), 'domain', d, 'frames', fr), ...
               'steps', struct('kind', {}, 'id', {}, 'in', {}, 'out', {}, 'domain_in', {}, ...
                               'domain_out', {}, 'frames_out', {}, 'range', {}, 'note', {}));
end
% Author: Diellor Basha, 2026
