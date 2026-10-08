% INGEST  Backend data preparation: constant-Q tile summaries of a recording.
%
% Reads a recording once and writes queryable, exactly mergeable statistics on
% constant-Q time-frequency tiles (centre + extent on both axes), so that a later query
% can find the channels, times and frequencies where something of interest may be and
% read raw samples only there. Features are the first family of tile labels.
%
% Design: docs/2026-09-22-ingest-design.md    Inventory: docs/atlas/INVENTORY.md
%
% THE BANK (2026-09-22): the default is @timefilterbank, a designed tight constant-Q frame
% on the log-frequency axis anchored at 1 Hz, evaluated per member at its own rate
% (docs/2026-09-22-timefilterbank-design.md): 11-14 s per subject at 600 Hz, 60 s at
% 2400 Hz, bands partition sum x^2 exactly. cwtfilterbank (Bank="morse") is the reference
% path, 20x slower, with its own paging machinery.
%
% M1 -- the pure functions:
%   rheome.ingest.config   - parameters in canonical order, with a hash
%   rheome.ingest.grid     - dyadic frame grid: frames per level, centres, extents (no data)
%   rheome.ingest.bank     - full-range cwtfilterbank; scales grouped into octave bands with
%                     nominal and measured extents; frame bounds A, B
%   rheome.ingest.reduce   - level-0 statistics from samples and coefficients, under a mask
%   rheome.ingest.rollup   - every coarser level from level 0, by pairing rows
%
% The position pyramid (select design §2b):
%   rheome.ingest.groups      - group rows as exact sums/maxima over the sensor tree's members
%                        (rheome.sensors.tree); build writes them as LNN_g_* when positions exist
%
% The wavelength axis (select design §2c):
%   rheome.ingest.spatial     - a tight log-itersine frame on the sensor graph's spectrum applied
%                        to the field per sample; energy and envMax per tile, channel and
%                        spatial band (LNN_s_*, LNN_gs_*); exact partition at the root
%
% Paging as tiling (design §3.12):
%   rheome.ingest.pages       - one grid level per band from its support; bands grouped into jobs
%   rheome.ingest.reducepaged - level 0 page by page, one master-anchored sub-bank per job
%   rheome.ingest.cone        - the analytic cone of influence, shared by both reduce paths
%
% Preview (no transform at all):
%   rheome.ingest.preview  - a min/max pyramid beside the recording, for looking before loading:
%                     8 bytes per channel per tile, doubled by the pyramid, and nothing
%                     else. rheome.select.open reads it, @recordingbrowser opens on it in envelope
%                     mode. Floor 1/60 s costs 10 %% of the recording and stays pixel-exact
%                     down to a ~19 s window; 0.25 s costs 0.7 %% and reaches 278 s.
%   rheome.ingest.support  - which bands a window of a given length can resolve (support = k/f):
%                     choosing a stretch in the explorer chooses the analysis that fits it
%   rheome.ingest.bound    - cast an extremum to single ROUNDING OUTWARD, so it stays a bound
%
% M2 -- the store:
%   rheome.ingest.build    - recording (pagedrecording) -> +data/<name>/ingest__F*__V*__*.mat,
%                     flat per-level top-level variables + meta/bands/frame/grid
%   rheome.ingest.totable  - a store as long-form rows (recording, channel, level, tCenter,
%                     tExtent, band, fCenter, fExtent, stat, value)
%
% Scripts: ingest_frame_omega.m builds the frame stores (600 Hz x2, 2400 Hz native);
% ingest_build_omega.m / ingest_paged_omega.m are the Morse-path builds and comparisons.
% All print the numbers the notes file records (docs/2026-09-22-ingest-notes.md).
%
% Author: Diellor Basha, 2026
