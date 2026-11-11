function [G, info] = bank(Lambda, f, psiGraph, psiTime, K, varargin)
% JTV.BANK  Build a joint time-vertex filterbank from separable and non-separable kernels.
%
%   [G, info] = rheome.jtv.bank(Lambda, f, psiGraph, psiTime)
%   [G, info] = rheome.jtv.bank(Lambda, f, psiGraph, psiTime, K)
%   [G, info] = rheome.jtv.bank(..., 'type','js', 'axes',ax, 'labels',{...})
%
% Follows GSPBox's gsp_jtv_design_dgw:
%
%   W(lambda,f) = K(lambda,f) .* psi_graph(lambda) .* psi_time(f)
%
% and forms every combination, so a bank of Ng graph kernels, Nt time kernels and Nk joint kernels
% has Ng*Nt*Nk members. The separable case is K = 1; a speed-selective filter is a non-separable K,
% e.g. rheome.jointfilterbank.speedkernel(c, d), or @(l,f) exp(-(f - c*sqrt(l)/(2*pi)).^2 / (2*df^2)).
%
% UNITS. The frequency axis f is in HERTZ and SIGNED, and every time factor and joint kernel
% receives f in hertz (never omega = 2*pi*f). A kernel that needs an angular frequency converts
% itself (rheome.jointfilterbank.speedkernel does, so a speed stays in length/s).
%
% FILTER TYPE -- 'js' vs 'ts'. GSPBox distinguishes filters defined in (lambda,f) from those
% defined in (lambda,t), and the distinction is not cosmetic: a 'ts' kernel (diffusion, wave,
% dampedwave, kleingordon, schrodinger) is a function of a time LAG, the response to an impulse as
% elapsed time goes by. It cannot multiply a joint (lambda,f) SPECTRUM: it must be transformed along
% the time axis first (rheome.filters.jtv), and applying it as though it were 'js' filters the wrong
% axis while looking entirely plausible. This bank is built and applied in (lambda,f), so 'js' is the
% only type it evaluates; 'ts' is REFUSED rather than quietly mishandled. The field is recorded on
% info so a caller storing banks can tell them apart.
%
% THE AXES ARE FIRST CLASS. Pass 'axes' (a rheome.jtv.axes) and the bank carries the axes it was designed
% on; rheome.jtv.analysis then checks the spectrum it is handed against them. This catches the failure a
% bare size check does NOT: two different bands can retain the same number of bins, so [K x nOmega]
% agreeing proves nothing about the two representations belonging together.
%
% INPUTS:
%   Lambda    [K x 1] eigenvalues
%   f         [1 x nOmega] frequencies (Hz, signed); passed to the kernels as they are
%   psiGraph  function handle or cell array, each @(lambda) -> gain  (default {@(l) ones(size(l))})
%   psiTime   function handle or cell array, each @(f) -> gain, f in Hz   (default {@(f) ones(size(f))})
%   K         function handle or cell array, each @(lambda,f) -> gain, f in Hz (default {@(l,f) 1})
%
% NAME-VALUE:
%   'type'    'js' (default, (lambda,f)) | 'ts' ((lambda,t) -- REFUSED, see above)
%   'axes'    a rheome.jtv.axes struct; recorded on info and checked by rheome.jtv.analysis
%   'labels'  {1 x Nf} cellstr naming the members, for figures and tables
%
% OUTPUTS:
%   G     [K x nOmega x Nf] gains
%   info  .Nf .Ng .Nt .Nk  .index [Nf x 3] giving (ig, it, ik) per member
%         .type  'js'      .axes  the axes the bank was designed on ([] if not supplied)
%         .labels {1 x Nf} .separable  true iff every joint kernel K is identically 1
%
% ⚠ A BANK IS NOT USEFUL UNTIL IT COVERS. Check rheome.jtv.bounds before rheome.jtv.dual: the dual divides by
% the frame operator, so any uncovered region of the (lambda,omega) plane is unreconstructable.
%
% See also: rheome.jtv.bounds, rheome.jtv.dual, rheome.jtv.axes, rheome.jtv.analysis, rheome.filters.frame, rheome.filters.psitime
%
% Author: Diellor Basha, 2026 (after F. Grassi, GSPBox)

    if nargin < 3 || isempty(psiGraph), psiGraph = {@(l) ones(size(l))}; end
    if nargin < 4 || isempty(psiTime),  psiTime  = {@(f) ones(size(f))}; end
    if nargin < 5 || isempty(K),        K        = {@(l,f) ones(size(l))}; end
    if ~iscell(psiGraph), psiGraph = {psiGraph}; end
    if ~iscell(psiTime),  psiTime  = {psiTime};  end
    if ~iscell(K),        K        = {K};        end

    p = inputParser;
    p.addParameter('type',   'js');
    p.addParameter('axes',   []);
    p.addParameter('labels', {});
    p.parse(varargin{:});
    o = p.Results;

    if ~any(strcmpi(o.type, {'js','ts'}))
        error('jtv:bank:type', 'type must be ''js'' (lambda,f) or ''ts'' (lambda,t), got ''%s''.', o.type);
    end
    if strcmpi(o.type, 'ts')
        error('jtv:bank:tsUnsupported', ...
            ['type ''ts'' defines the kernel in (lambda,t): it acts on a time LAG, as a propagator ' ...
             '(diffusion, wave, ...), and cannot multiply a joint (lambda,f) spectrum, which this bank ' ...
             'is evaluated on. Transform the kernel along the time axis first (rheome.filters.jtv) and pass it as ''js''.']);
    end

    lam = double(Lambda(:));
    om  = double(f(:)).';                       % Hz, signed: the kernels' argument
    [LL, WW] = ndgrid(lam, om);

    Ng = numel(psiGraph);  Nt = numel(psiTime);  Nk = numel(K);
    Nf = Ng*Nt*Nk;
    G  = zeros(numel(lam), numel(om), Nf);
    idx = zeros(Nf, 3);
    sep = true;

    m = 0;
    for ik = 1:Nk
        Kv = K{ik}(LL, WW);
        if isscalar(Kv), Kv = Kv * ones(size(LL)); end
        if any(Kv(:) ~= 1), sep = false; end
        for ig = 1:Ng
            gv = psiGraph{ig}(lam);  gv = gv(:);
            for it = 1:Nt
                tv = psiTime{it}(om);  tv = tv(:).';
                m = m + 1;
                G(:,:,m) = Kv .* (gv * tv);
                idx(m,:) = [ig it ik];
            end
        end
    end

    if isempty(o.labels)
        o.labels = arrayfun(@(mm) sprintf('g%d/h%d/K%d', idx(mm,1), idx(mm,2), idx(mm,3)), ...
            1:Nf, 'UniformOutput', false);
    elseif numel(o.labels) ~= Nf
        error('jtv:bank:labels', 'labels has %d entries but the bank has %d members.', ...
            numel(o.labels), Nf);
    end

    info.Nf = Nf;  info.Ng = Ng;  info.Nt = Nt;  info.Nk = Nk;  info.index = idx;
    info.type = lower(o.type);  info.axes = o.axes;  info.labels = o.labels(:).';
    info.separable = sep;
end

% Author: Diellor Basha, 2026
