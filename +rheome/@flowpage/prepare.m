function B = prepare(name, varargin)
% FLOWPAGE.PREPARE  Build the anatomy bundle a flowpage needs. Slow; do it once.
%
%   B = rheome.flowpage.prepare('sub01')
%   B = rheome.flowpage.prepare(name, 'Store','native', 'PageSec',10, 'GraphFilters',6, ...)
%
% Assembles everything that depends on the ANATOMY and not on the band or the page: the flow
% context, the fused curl kernel, the LBO basis, the graph (spatial-scale) bank, the atlas
% and a pagedrecording opened on the chosen store. ~40 s, dominated by loading the Dirac and
% LBO bases; every page and band afterwards reuses it.
%
% ⚠ THE MARGIN IS SET BY THE LOWEST FREQUENCY IN FreqLimits, not by the band you will browse.
% The pager is opened once with a margin wide enough for the whole master range so that
% switching band never needs a different page grid -- at 1 Hz that is +/-3.215 s.
%
% OUTPUT (struct B): .pager .Kc [Ks x C] .Phi [V x Ks] .Lambda .gfb .atlas .wVert .fs .chSel
%
% See also: rheome.flowpage, rheome.flowbrowser, rheome.flow.context, rheome.flow.curl
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('Store',        'native');
    p.addParameter('PageSec',      10);
    p.addParameter('FreqLimits',   [1 60]);
    p.addParameter('VoicesPerOctave', 10);
    p.addParameter('GraphFilters', 6);
    p.addParameter('DiracK',       400);
    p.addParameter('LboK',         400);
    p.addParameter('Atlas',        'Desikan-Killiany');
    p.parse(varargin{:});
    o = p.Results;

    fprintf('rheome.flowpage.prepare[%s]: building anatomy bundle (~40 s)\n', name);
    ctx = rheome.flow.context(name, o.DiracK, o.LboK);
    kc  = rheome.flow.curl(ctx);

    pr0   = rheome.pagedrecording(name, 'Store', o.Store);
    chSel = find(strcmp(pr0.ChannelType(:).','MEG') & pr0.ChannelFlag(:).' > 0);
    if ~isequal(chSel(:), ctx.iSel(:))
        error('flowpage:channels', ...
            ['Page channels (%d) do not match rheome.flow.context channels (%d). The flow ' ...
             'operators would be applied to the wrong sensor axis.'], numel(chSel), numel(ctx.iSel));
    end
    fs = pr0.SamplingFrequency;

    % the margin is the master range's, so switching band never regrids the pages
    probeN = ceil(8 * fs / o.FreqLimits(1));
    ov = 0;
    for it = 1:10
        fbp = cwtfilterbank('SignalLength', round(o.PageSec*fs) + 2*ov, 'SamplingFrequency', fs, ...
                            'FrequencyLimits', o.FreqLimits, 'VoicesPerOctave', o.VoicesPerOctave);
        need = ceil(max(waveletsupport(fbp).End) * fs);
        if need <= ov, break; end
        ov = need;
    end

    B = struct();
    B.pager  = rheome.pagedrecording(name, 'Store', o.Store, 'PageLength', round(o.PageSec*fs), ...
                              'Overlap', ov, 'Channels', chSel, 'Precision', 'single');
    B.Kc     = kc.coeffOperator;
    % ⚠ SINGLE, DELIBERATELY. Phi is [20484 x 800]; in double it is 131 MB and every display
    % frame streams all of it through cache, which is the actual cost of a redraw -- the
    % synthesis is memory-bandwidth bound, not flop bound. Single halves it and doubles the
    % frame rate. The double basis stays in ctx for anything that needs it.
    B.Phi    = single(ctx.lbo.Phi);
    B.Lambda = ctx.lbo.Lambda;
    B.fs     = fs;
    B.chSel  = chSel;
    B.wVert  = full(sum(ctx.lbo.Mass, 2));
    B.surface = ctx.S;
    B.gfb    = rheome.graphfilterbank(max(ctx.lbo.Lambda), 'NumFilters', o.GraphFilters, ...
                               'Transform', rheome.graphtransform.eigen(ctx.lbo.Phi, ctx.lbo.Mass, ctx.lbo.Lambda));
    B.atlas  = rheome.load.atlas(name, o.Atlas);
    B.name   = name;
    B.probeN = probeN;
    B.freqLimits = o.FreqLimits;
    B.voicesOct  = o.VoicesPerOctave;

    fprintf('  %d ch, %d modes, %d vertices, %d scales | %d pages of %g s (margin %d = %.2f s)\n', ...
        numel(chSel), numel(B.Lambda), size(B.Phi,1), B.gfb.NumMembers, ...
        B.pager.NumPages, o.PageSec, ov, ov/fs);
end

% Author: Diellor Basha, 2026
