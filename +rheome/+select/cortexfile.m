function f = cortexfile(db)
% SELECT.CORTEXFILE  The cortical feature sidecar of a store: <store>__cortex.mat, beside it.
%
%   f = rheome.select.cortexfile(db)
%
% The mergeable cortical relations (feature_cortex, feature_cortex_space) and the cortex_node rows live
% in their own file next to the ingest store, which is never modified -- the same rule the label and
% measurement sidecar (<store>__labels.mat) follows. rheome.flow.cortexfeatures writes it; rheome.select.rows reads it.
%
% See also: rheome.flow.cortexfeatures, rheome.select.rows, rheome.select.open
%
% Author: Diellor Basha, 2026

    [d, stem] = fileparts(char(db.file));
    f = fullfile(d, [stem '__cortex.mat']);
end

% Author: Diellor Basha, 2026
