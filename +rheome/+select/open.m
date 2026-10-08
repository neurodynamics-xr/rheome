function db = open(file)
% SELECT.OPEN  A handle on a tile store: metadata, a matfile, a level cache and cost counters.
%
%   db = rheome.select.open(file)
%
% Reads go through rheome.select.level so they can be counted: .cost is a containers.Map (a
% handle, so it accumulates through a value struct) with 'bytes' (array bytes loaded
% from disk), 'reads' (matfile reads) and 'tiles' (tile rows examined by executors).
%
% Author: Diellor Basha, 2026

    arguments
        file (1,:) char
    end
    if exist(file, 'file') ~= 2
        error('select:open:missing', 'No store at %s.', file);
    end
    m = matfile(file);
    db = struct();
    db.file  = file;
    db.m     = m;
    db.meta  = m.meta;
    db.bands = m.bands;
    db.grid  = m.grid;
    db.frame = m.frame;
    db.recording_id = db.meta.name;
    if isfield(db.meta, 'groupNodes'), db.groupNodes = db.meta.groupNodes; else, db.groupNodes = []; end
    if isfield(db.meta, 'tree') && ~isempty(db.meta.tree), db.tree = db.meta.tree; else, db.tree = table(); end
    if isfield(db.meta, 'sbands') && ~isempty(db.meta.sbands), db.sbands = db.meta.sbands; else, db.sbands = table(); end
    if isfield(db.grid, 'channelLevel'), db.channelLevel = db.grid.channelLevel; else, db.channelLevel = 0; end
    % the envelope exemption: min and max may reach below the channel floor (rheome.ingest.config
    % ChannelEnvelope). A store built before it exists has no such rows, so the two floors coincide.
    if isfield(db.grid, 'envelopeLevel'), db.envelopeLevel = db.grid.envelopeLevel; else, db.envelopeLevel = db.channelLevel; end
    if isfield(db.grid, 'spaceSlots'), db.spaceSlots = db.grid.spaceSlots; else, db.spaceSlots = table(); end
    % the HOME level of each band: where its non-mergeable spectral peak is stored
    if isfield(db.grid, 'homeAt'), db.homeAt = db.grid.homeAt; else, db.homeAt = {}; end
    if isfield(db.grid, 'chanSpace'), db.chanSpace = db.grid.chanSpace; else, db.chanSpace = 1:height(db.sbands); end
    % A preview store (rheome.ingest.preview) carries min and max and nothing else: no bank, no
    % bands, no sums. Readers check this rather than discovering it as a missing variable.
    if isfield(db.meta, 'kind'), db.kind = char(db.meta.kind); else, db.kind = 'tile'; end
    db.preview = strcmp(db.kind, 'preview');
    db.cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
    db.cost  = containers.Map({'bytes', 'reads', 'tiles'}, {0, 0, 0});
    [d, n] = fileparts(file);
    db.labelFile = fullfile(d, [n '__labels.mat']);
end
% Author: Diellor Basha, 2026
