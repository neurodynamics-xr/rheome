function [db, choice] = resolve(src, kind)
% RECORDINGBROWSER.RESOLVE  A recording -> an open handle on one of its pyramids.
%
%   db = rheome.recordingbrowser.resolve('sub01')              % the best store it has
%   db = rheome.recordingbrowser.resolve('sub01', 'preview')   % the min/max pyramid
%   db = rheome.recordingbrowser.resolve(path)                         % a named store file
%   db = rheome.recordingbrowser.resolve(db)                           % already open: passed through
%
% ⭐ A RECORDING USUALLY HAS MORE THAN ONE PYRAMID. The tile store carries bands and moments
% at 0.25 s and up; a preview store (rheome.ingest.preview) carries min and max alone, but reaches
% far finer, often 1/60 s. The browser is about the RECORDING, so it picks between them
% rather than making the caller name a file:
%
%   'auto'     the tile store if there is one, else the preview (the default)
%   'tile'     the tile store: bands, moments, every statistic
%   'preview'  the min/max pyramid: the finest one present
%   'native'   a store built from the native-rate recording
%
% A dataset with several tile stores gives the default-rate frame store, newest first,
% because that is the one the rest of the pipeline builds on. `choice` says what was opened.
%
% Author: Diellor Basha, 2026

    if nargin < 2 || isempty(kind), kind = 'auto'; end
    kind = char(kind);
    if ~ismember(kind, {'auto','tile','preview','native'})
        error('recordingbrowser:resolve', ...
            'Store must be ''auto'', ''tile'', ''preview'' or ''native'', got ''%s''.', kind);
    end

    if isstruct(src) && isfield(src, 'grid') && isfield(src, 'm')
        db = src;  choice = i_kindof(db);  return
    end
    src = char(src);
    if exist(src, 'file') == 2
        db = rheome.select.open(src);  choice = i_kindof(db);  return
    end
    d = fullfile(rheome.load.root(), src);
    if exist(d, 'dir') ~= 7
        error('recordingbrowser:resolve', ...
            ['No dataset or store at ''%s''. Give a dataset name under %s, a path to a ' ...
             'store, or an open select handle.'], src, rheome.load.root());
    end

    tiles = i_tiles(src, kind);
    prevs = dir(fullfile(d, 'preview__*.mat'));
    switch kind
        case 'preview'
            if isempty(prevs)
                error('recordingbrowser:resolve', ...
                    'No preview store for ''%s''. Run  rheome.ingest.preview(''%s'')  first.', src, src);
            end
            db = rheome.select.open(i_finest(prevs));  choice = 'preview';
        case {'tile','native'}
            if isempty(tiles)
                error('recordingbrowser:resolve', ...
                    'No %s tile store for ''%s''. Run  rheome.ingest.build(''%s'')  first.', kind, src, src);
            end
            db = rheome.select.open(tiles);  choice = kind;
        otherwise                                       % auto: the richer store wins
            if ~isempty(tiles)
                db = rheome.select.open(tiles);  choice = 'tile';
            elseif ~isempty(prevs)
                db = rheome.select.open(i_finest(prevs));  choice = 'preview';
            else
                error('recordingbrowser:resolve', ...
                    ['''%s'' has no pyramid. Run  rheome.ingest.preview(''%s'')  to look at it, or ' ...
                     'rheome.ingest.build(''%s'')  for the full tile store.'], src, src, src);
            end
    end
end

% The finest preview, since a coarser one can always be rolled up from it but not the other
% way round.
function f = i_finest(prevs)
    flo = zeros(numel(prevs), 1);
    for i = 1:numel(prevs)
        t = regexp(prevs(i).name, 'preview__F([0-9.eE+-]+)', 'tokens', 'once');
        if isempty(t), flo(i) = Inf; else, flo(i) = str2double(t{1}); end
    end
    [~, i] = min(flo);
    f = fullfile(prevs(i).folder, prevs(i).name);
end

function f = i_tiles(src, kind)
    f = '';
    C = rheome.select.catalog(rheome.load.root());
    C = C(C.dataset == string(src), :);
    if isempty(C), return; end
    if strcmp(kind, 'native'), keep = C.store == "native";
    else,                      keep = C.store == "default" & C.bank == "frame";
    end
    if ~any(keep) && ~strcmp(kind, 'native'), keep = C.store == "default"; end
    if ~any(keep), return; end
    C = sortrows(C(keep, :), 'created', 'descend');
    f = char(C.file(1));
end

function k = i_kindof(db)
    if isfield(db, 'preview') && db.preview, k = 'preview'; else, k = 'tile'; end
end
% Author: Diellor Basha, 2026
