function T = registry()
% SELECTION.REGISTRY  Every way we index or weigh a manifold's elements, as one type.
%
%   T = rheome.selection.registry()
%   T = rheome.selection.registry(); T(T.merges, :)        % the ones a coarser cell can be summed from
%
% ⭐ THE FOURTH PILLAR. +domain says where a field lives, +fieldtype says what it is,
% rheome.operators.registry says which arrows move it -- and this says how we CUT it. A tile, a window, a
% frequency band, a tree node, an ROI and a wavelet member are all the same kind of thing: a weight
% over the elements of one domain. They differed only in which package happened to implement them
% (+select for tiles, rheome.pipeline.window for frames, @graphfilterbank for kernels, rheome.geom.tree for
% patches), and nothing related them, so the same question had to be re-answered per mechanism.
%
% ⭐⭐ THE COLUMN THAT MATTERS IS `merges`. It is the difference between a partition and a frame, and
% getting it wrong is the most expensive mistake in this codebase:
%   PARTITION (merges)   disjoint and exhaustive, so a coarser cell is EXACTLY the sum of finer ones.
%                        rheome.geom.tree, the dyadic time tiles, the constant-Q bands, the Weyl octaves.
%                        Verified: energy rolls up at 4.3e-16 relative error across tree depths.
%   KERNEL (does not)    smooth and overlapping, a frame you can invert but not add up. One member is
%                        a PROJECTION: measured, a single joint wavelet keeps 0.048 = 1/21 of its
%                        (tile x octave) cell, and the graph bank has A/B = 0.274 against the
%                        temporal bank's exactly 1.000.
%   SUBSET (does not)    a plain index range or an ROI list. ⚠ atlas ROIs are NOT disjoint and NOT
%                        exhaustive, so they never merge even though they look like a partition.
%
% ⚠⚠ AND `merges` IS A PROPERTY OF THE SELECTION, NOT OF THE MEASUREMENT. Even under a partition only
% quantities LINEAR in the accumulated statistic roll up. In rheome.flow.windowtable, `energy` is the only
% column a parent may take as a sum of its children; peak, crest, vortex count and the gauge shares
% all have to be recomputed. A partition permits merging; it does not make every column mergeable.
%
% ⭐ THE PAYOFF: rheome.select.ladder and rheome.geom.ladder are the same function. Both pair a partition of one
% axis with a partition of its conjugate and ask how many units of the first fit in a cell of the
% second -- 22.6 cycles per time tile, 2.01 diameters per wavelength. Same predicate, two axes, and
% the 11x gap between those numbers is a fact about the instrument that only shows once they are
% recognised as the same quantity.
%
% COLUMNS
%   id            the selection, as called
%   fcn           the function that builds it -- must resolve (same contract as rheome.operators.registry)
%   domain_kind   which domain kinds it cuts: complex | graph | points | index
%   axes          the axis or axes it cuts, as a cell of names: {'time'}, {'eigenmode','frequency'}
%   rank          the element rank it selects on
%   kind          partition | kernel | subset
%   merges        can a coarser cell be summed from finer ones (true only for partitions)
%   nested        is there a parent/child hierarchy
%   tight         for kernels: 1 if the frame is tight (A == B), 0 if not, NaN where not applicable.
%                 ⚠ a DOUBLE, not a logical: cell2mat refuses a column mixing NaN with true/false,
%                 which is how this table failed to build the first time.
%   notes
%
% ⭐ A PRODUCT SELECTION IS NOW TYPED, NOT LABELLED. `joint_wavelet` has domain_kind 'product' and two
% entries in `axes`, and rheome.domain.product builds the descriptor it cuts. The earlier version of this
% table carried a stringly `axis2` column because no product domain existed; before that it typed
% joint_wavelet as a plain eigenmode selection, which is the error that prompted all of this -- the one
% genuinely product-valued row got one axis. The invariant is now checkable:
% numel(axes) == 2 <=> domain_kind == 'product'.
%
% See also: rheome.domain.kinds, rheome.fieldtype.registry, rheome.operators.registry, rheome.selection.describe,
%           rheome.selection.check, rheome.select.ladder, rheome.geom.ladder
%
% Author: Diellor Basha, 2026

    r = { ...
    % id, fcn, domain_kind, axes (cell), rank, kind, merges, nested, tight, notes
    'time_tile', 'rheome.ingest.grid', 'index', {'time'}, 'element', 'partition', true, true, NaN, 'dyadic: level L tile is 2^L base. A parent is exactly its two children; rheome.select.level reads one level.'
    'time_window', 'rheome.pipeline.window', 'index', {'time'}, 'element', 'subset', false, false, NaN, 'an explicit span in seconds or samples: the bookkeeping step of a plan'
    'burst_window', 'rheome.select.bursts', 'index', {'time'}, 'element', 'subset', false, false, NaN, 'data-driven spans from an envelope threshold. WARNING overlapping bursts share samples, so consecutive ones are not independent.'
    'freq_band', 'rheome.select.ladder', 'index', {'frequency'}, 'element', 'partition', true, true, NaN, 'constant-Q octaves; the bank is designed so bands partition sum x^2 exactly. 22.6 cycles per tile in every octave.'
    'cortex_tile', 'rheome.geom.tree', 'complex', {'cortex'}, 'vertex', 'partition', true, true, NaN, 'recursive spectral bisection by AREA, so diameter falls by sqrt(2) per depth: two depths per spatial octave.'
    'sensor_group', 'rheome.sensors.tree', 'graph', {'sensor'}, 'vertex', 'partition', true, true, NaN, 'the sensor-side twin of cortex_tile; rheome.ingest.groups writes its rows'
    'space_octave', 'rheome.geom.ladder', 'index', {'eigenmode'}, 'element', 'partition', true, true, NaN, 'Weyl octaves on the eigenvalue line, x4 modes per octave. Counted vs predicted 0.95-0.98.'
    'dirac_group', 'rheome.load.dirac', 'index', {'eigenmode'}, 'element', 'partition', true, false, NaN, 'the 4-fold quaternionic degeneracy. WARNING the GROUP is the invariant unit; a single coefficient is a gauge.'
    'roi', 'rheome.load.atlas', 'complex', {'cortex'}, 'vertex', 'subset', false, false, NaN, 'WARNING atlas scouts are neither disjoint nor exhaustive, so they never merge. Sizes run 37-110 mm.'
    'time_wavelet', 'rheome.timefilterbank', 'index', {'time'}, 'element', 'kernel', false, false, 1, 'constant-Q voices; A = B = 1.0000 exactly. One voice keeps 0.245 of its band.'
    'wavelet_tile', 'rheome.selection.wavelettile', 'product', {'time','frequency'}, 'element', 'kernel', false, true, 1, 'A TILING OF WAVELETS: centre and time SUPPORT where a tile has centre and span. NESTED EXACTLY: one support is two of the next level, since supportCycles is scale-invariant and fc halves. Two strides -- strideSelect = support (redundancy 1, the tile lattice) and strideFrame = 1/rate (redundancy 12.5, needed to invert). NOT mergeable: one voice holds 0.245 of its band, and the tile is 1.50 supports long.'
    'graph_wavelet', 'rheome.graphfilterbank', 'index', {'eigenmode'}, 'element', 'kernel', false, false, 0, 'A = 0.195, B = 0.712, A/B = 0.274 on the cached spectrum: invertible, not tight. One member keeps ~0.196 of its octave.'
    'joint_wavelet', 'rheome.jointfilterbank', 'product', {'eigenmode','frequency'}, 'element', 'kernel', false, false, 0, 'THE ONE PRODUCT SELECTION: the separable psi_g(lambda)*psi_t(omega). One member keeps 0.048 of its (tile x octave) cell.'
    'mexhat_scale', 'rheome.filters.mexhat', 'index', {'eigenmode'}, 'element', 'kernel', false, false, 0, 'a single scale-normalised LoG member. WARNING (t*lambda)exp(-t*lambda) IS already normalised; dividing by its norm destroys the scale selection.'
    'heat_scale', 'rheome.filters.heat', 'index', {'eigenmode'}, 'element', 'kernel', false, false, 0, 'low-pass: a partition of UNITY rather than of energy, so it does not merge either'
    };
    % ⚠ `axes` is a CELL column (one or two axis names), so it cannot go through string() with the
    % others; it is assigned after the table is built.
    T = table(string(r(:,1)), string(r(:,2)), string(r(:,3)), string(r(:,5)), ...
              string(r(:,6)), cell2mat(r(:,7)), cell2mat(r(:,8)), cell2mat(r(:,9)), string(r(:,10)), ...
              'VariableNames', {'id','fcn','domain_kind','rank','kind','merges','nested', ...
                                'tight','notes'});
    T.axes = r(:,4);
    T = movevars(T, 'axes', 'After', 'domain_kind');
end

% Author: Diellor Basha, 2026
