classdef pagedrecording
% PAGEDRECORDING  A long sensor recording read one time page at a time.
%
%   pr = rheome.pagedrecording('sub01', 'PageLength', 6000, 'Overlap', 1929)
%   pr = rheome.pagedrecording(recordingFile, ...)
%   p  = page(pr, k)
%
% The object declares the FULL time axis without touching the data; page(pr,k) is the
% only thing that reads. Backed by a v7.3 store whose F is a TOP-LEVEL variable, so
% m.F(:, a:b) is a genuine partial read -- measured 2.4 ms for a 4 s page against 250 ms
% for a full load. See rheome.import.recording, which writes that store.
%
% ⚠ WHY THE STORE IS SEPARATE FROM study.mat. MatFile objects only support '()' indexing
% on TOP-LEVEL variables: `m.rec.F(:, a:b)` fails outright ("MatFile objects only support
% '()' indexing"), and `m.rec` loads the whole 977 MB struct. Nesting F inside a struct
% defeats paging entirely, which is why rheome.import.recording surfaces it as its own variable.
%
% ⚠ CORE vs SPAN -- the distinction the caller must respect. A page is read WIDER than it
% is used. `.samples` is the span actually read; `.core` masks the part this page OWNS.
% The cores of successive pages PARTITION the record exactly once, so summing a statistic
% over cores counts every sample once; summing over spans double-counts the overlap and
% silently inflates every total. The margin exists to absorb the wavelet cone of
% influence: set Overlap from waveletsupport(fb), whose lowest filter fixes the longest
% support (a [1 60] Hz bank at 600 Hz reaches +/- 3.215 s = 1929 samples).
%
% ⚠ EDGE PAGES CANNOT HAVE THEIR FULL MARGIN. The first and last pages clip against the
% record, so their cores carry cone-contaminated samples. `.clipped` = [pre post] reports
% how many samples of margin were unavailable; nonzero means treat that end with care.
%
% PROPERTIES
%   File  NumSamples  SamplingFrequency  Time  Duration
%   ChannelFlag  ChannelName  ChannelType  Channels  NumChannels
%   PageLength  Overlap  Precision  NumPages
%
% See also: rheome.import.recording, page, cwtfilterbank, waveletsupport
%
% Author: Diellor Basha, 2026

    properties (SetAccess = immutable)
        File                                   % the v7.3 store
        NumSamples                             % samples on the full time axis
        SamplingFrequency                      % Hz
        Time                                   % [1 x nT] full time axis (seconds)
        ChannelFlag                            % [nCh x 1] +1 good / -1 bad, as stored
        ChannelName                            % {1 x nCh}
        ChannelType                            % {1 x nCh}
        Channels                               % row indices into the stored F
        Comment
        Source
    end

    properties
        PageLength (1,1) double {mustBePositive, mustBeInteger} = 1
        Overlap    (1,1) double {mustBeNonnegative, mustBeInteger} = 0
        Precision  (1,:) char = 'single'
    end

    properties (Dependent, SetAccess = private)
        NumChannels
        NumPages
        Duration
    end

    properties (Access = private)
        MF_                                    % matfile handle
    end

    methods
        function obj = pagedrecording(src, varargin)
            if nargin < 1 || isempty(src)
                error('pagedrecording:source', ...
                    'A dataset name or a recording-store path is required.');
            end
            f = i_resolve(src, o_store(varargin));
            if ~exist(f, 'file')
                error('pagedrecording:missing', ...
                    ['No recording store at %s. Run  rheome.import.recording(name)  to surface ' ...
                     'the cached recording as a v7.3, chunk-readable store.'], f);
            end

            p = inputParser;
            p.addParameter('PageLength', []);
            p.addParameter('Overlap',    0);
            p.addParameter('Channels',   []);
            p.addParameter('Precision',  'single');
            p.addParameter('Store',      'default');   % 'default' | 'native'
            p.parse(varargin{:});
            o = p.Results;

            obj.File = f;
            obj.MF_  = matfile(f);
            hdr      = i_header(f);

            obj.NumSamples        = hdr.nT;
            obj.SamplingFrequency = hdr.sfreq;
            obj.Time              = hdr.Time(:).';
            obj.ChannelFlag       = hdr.ChannelFlag(:);
            obj.ChannelName       = hdr.ChannelName;
            obj.ChannelType       = hdr.ChannelType;
            obj.Comment           = hdr.Comment;
            obj.Source            = hdr.source;

            if isempty(o.Channels)
                obj.Channels = (1:hdr.nCh).';
            else
                ch = double(o.Channels(:));
                if any(ch < 1 | ch > hdr.nCh | mod(ch,1) ~= 0)
                    error('pagedrecording:channels', ...
                        'Channels must be integer row indices in 1..%d.', hdr.nCh);
                end
                obj.Channels = ch;
            end

            if isempty(o.PageLength), o.PageLength = obj.NumSamples; end
            obj.PageLength = o.PageLength;
            obj.Overlap    = o.Overlap;
            switch lower(o.Precision)
                case {'single','double'}, obj.Precision = lower(o.Precision);
                otherwise
                    error('pagedrecording:precision', ...
                        'Precision must be ''single'' or ''double'', got ''%s''.', o.Precision);
            end
        end

        function n = get.NumChannels(obj), n = numel(obj.Channels); end
        function n = get.NumPages(obj),    n = ceil(obj.NumSamples / obj.PageLength); end
        function d = get.Duration(obj),    d = obj.NumSamples / obj.SamplingFrequency; end

        p   = page(obj, k)
        rng = pagerange(obj, k)
        X   = read(obj, a, b)
        disp(obj)
    end
end

% ---- a dataset NAME resolves through the +data cache; anything else is a path ----
% 'default' is the study-cache store (rheome.import.recording); 'native' is the original
% acquisition rate read straight from BST-BIN (rheome.import.rawrecording) -- on the reference CTF data
% 2400 Hz against 600 Hz, which is what the multirate plan needs.
function f = i_resolve(src, store)
    src = char(src);
    if ~rheome.load.has(src), f = src; return; end
    switch lower(store)
        case 'native',  f = fullfile(rheome.load.root(), src, 'recording_native.mat');
        case 'default', f = fullfile(rheome.load.root(), src, 'recording.mat');
        otherwise
            error('pagedrecording:store', ...
                'Store must be ''default'' or ''native'', got ''%s''.', store);
    end
end

% Store must be read BEFORE inputParser runs, since it selects the file the header
% comes from; parse it out of varargin directly.
function s = o_store(args)
    s = 'default';
    for k = 1:2:numel(args)-1
        if ischar(args{k}) && strcmpi(args{k}, 'Store'), s = char(args{k+1}); end
    end
end

% ---- read every variable EXCEPT F: the header is small, F is the thing we are avoiding ----
function h = i_header(f)
    w    = whos('-file', f);
    want = setdiff({w.name}, {'F'});
    h    = builtin('load', f, want{:});
    for fld = {'Comment', 'source', 'ChannelName', 'ChannelType'}
        if ~isfield(h, fld{1}), h.(fld{1}) = []; end
    end
    if ~isfield(h, 'nT') || isempty(h.nT)
        iF = strcmp({w.name}, 'F');
        h.nT = w(iF).size(2);
    end
    if ~isfield(h, 'nCh') || isempty(h.nCh)
        iF = strcmp({w.name}, 'F');
        h.nCh = w(iF).size(1);
    end
    if ~isfield(h, 'ChannelFlag') || isempty(h.ChannelFlag)
        h.ChannelFlag = ones(h.nCh, 1);
    end
end

% Author: Diellor Basha, 2026
