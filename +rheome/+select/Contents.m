% SELECT  The query side of the tile store, as a relational database.
%
% The store is a set of relations (rheome.select.schema, rheome.select.ddl); a query is a value
% (rheome.select.query) that two executors interpret and must agree on: a full scan
% (rheome.select.scan, the reference) and a pruned descent of the pyramid (rheome.select.frames), which
% rheome.select.sql renders as the SQL a database runs. Writes are transactions (rheome.select.label).
% Design: docs/2026-09-22-select-design.md
%
% Schema and relations:
%   rheome.select.schema   - relations, columns, types, keys, foreign keys
%   rheome.select.ddl      - Postgres DDL with the indexes the plans assume
%   rheome.select.open     - a handle on a store: metadata, matfile, level cache, cost counters
%   rheome.select.level    - one level array, cached and counted
%   rheome.select.rows     - any relation materialised as a table (the reference view)
%   rheome.select.catalog  - every store under +data, one row per (dataset, store, bank, config)
%   rheome.select.cortexfile   - <store>__cortex.mat: the cortical mergeable sidecar (rheome.flow.cortexfeatures)
%   rheome.select.rollupcortex - finest-cell sums -> every node and level, exactly (a parent is its children)
%   rheome.select.labelkinds   - what each label kind may claim: calibration and evidence
%
% Queries:
%   rheome.select.query    - Stat/Op/Threshold/Level/Band/Channels/MinDuration as a value
%   rheome.select.scan     - reference executor: full scan at the level
%   rheome.select.frames   - pruned executor: semi-join descent from the top; cost.plan per level
%   rheome.select.sql      - the query as SQL (threshold CTE, survivor chain, hits, runs)
%   rheome.select.bursts   - two-phase: prune on envMax, refine on raw samples with the same members
%   rheome.select.cohort   - a query over a catalogue selection, thresholds per recording
%   (stats spaceEnergy / spaceEnvMax address the wavelength axis: space_band, feature_space,
%    feature_group_space)
%   rheome.select.tree     - descend the sensor tree at a time level (Scope="group", Depth=d);
%                     the position-axis twin of frames; the two compose
%
% Derived statistics (nothing below is stored: all of it is a function of the moments):
%   rheome.select.derive   - per-tile mean/rms/std/crest/spectral entropy/centroid/share/... at a
%                     level, for channels or tree nodes, over a time window; with no
%                     arguments, the vocabulary itself (name, kind, moments needed, formula)
%   rheome.select.coupling - phase-amplitude coupling per tile from the stored sums: mean vector
%                     length, preferred phase and the amplitude behind them, at the pair's
%                     own level or any coarser one (the sums merge, so coarser is exact)
%   rheome.select.named    - the physiological bands (delta..gamma) mapped onto the store's octave
%                     bands, with how well each one is covered and at which level it lives
%   rheome.select.ladder   - the constant-Q ladder: per band, its support, natural level, tile
%                     length, cycles per tile, the rate its members arrive at and the
%                     samples that costs, against the record's samples for the same tile
%   (min and max reach BELOW the channel floor -- rheome.ingest.config ChannelEnvelope -- because
%    they merge by extremum and are stored as outward-rounded bounds: that pair is a
%    min/max mipmap, and an exact envelope of one trace at any zoom is one array read)
%   @recordingbrowser - the app over all of this. Opens a RECORDING, picking between its
%                     tile store and its preview store (Store="auto"|"tile"|"preview"), and
%                     browses the channels at whatever level the zoom asks for: the finest
%                     under a tile budget that is the axis width in pixels, or the level a
%                     BAND implies (FollowBand). The strip draws each band at its own tile
%                     length (StripMode), Envelope draws the min/max band, and raw samples
%                     are read only below the diagonal (recordingbrowser_omega.m)
%
% Two bookkeeping matrices over the same nodes (rheome.select.schema marks which is which):
%   'mergeable'  the feature relations -- dense, one row per (unit, tile), a parent exactly
%                the merge of its children, so bounds hold and queries prune by descending
%   'note'       label and measure -- sparse, one row per (node, kind), whatever was
%                measured where it was measured; nothing merges, so they are scanned
%
% Transactions (labels and measurements share one log):
%   rheome.select.label    - append labels in one atomic, idempotent transaction
%   rheome.select.labels   - read labels and transactions, filtered
%   (Scope="cortex" addresses a rheome.geom.tree node on the far side of the inverse -- the
%    cortex_node dimension -- which is where a time-averaged magnitude, divergence or curl
%    under a cortical tile belongs; rheome.flow.sweep writes them)
%   rheome.select.measure  - write measurements onto nodes (Scope channel | group), same contract
%   rheome.select.measures - read measurements, filtered by kind, node, unit, value or window
%   rheome.select.context  - join measurements to the mergeable matrix at their own node and at
%                     any ancestor level: the state each measurement occurred in
%                     (@recordingbrowser/measure writes one onto the clicked tile, and
%                      measurements() reads them back with that context attached)
%   rheome.select.rollback - remove a transaction's rows atomically
%
% Author: Diellor Basha, 2026
