function ax = axes(varargin)
% JTV.AXES  The axes of a joint time-vertex representation, as a first-class object.
%
%   ax = rheome.jtv.axes('lambda',Lambda, 'f',f, 'fs',fs, ...)
%
% GSPBox carries the graph, the eigenvalue axis, the time axis and the frequency axis together in
% one structure, and derives everything from it. We had been passing Lambda and f loose between
% functions with nothing checking they belonged together -- which is how a taper meant for an
% analysis pass ended up corrupting a fit pass without anything raising an error.
%
% This is the object those axes live in. Carry it; check it with rheome.jtv.compatible before combining
% two representations.
%
% FIELDS (all optional at construction; whatever is supplied is recorded):
%   .lambda [K x 1]   eigenvalues of the scalar operator      .k  sqrt(lambda), wavenumbers rad/m
%   .nK               number of modes                          .lmax  max(lambda)
%   .f      [1 x Om]  RETAINED frequencies (Hz)                .omega  2*pi*f
%   .nOmega           number of retained bins                  .df  bin spacing (Hz)
%   .bands  [nb x 2]  the bands retained, Hz
%   .half   'positive' | 'full'   which half of the frequency axis is kept
%   .fs .nT .dur .nfft            the source time axis
%   .t0               time-origin offset (nonzero if padded)
%   .boundary 'ring' | 'path'     the implied time boundary condition -- 'ring' for a plain FFT
%                                 (periodic, wraps), 'path' if mirror-padded (Neumann, reflects)
%   .nV               vertices, if a surface is attached
%
% ⚠ .half MATTERS FOR REUSE. We keep the positive half only, so the representation is ANALYTIC and
% complex. GSPBox keeps the full axis and its signals are real. A time kernel designed for one is
% wrong for the other, silently.
%
% See also: rheome.jtv.compatible, rheome.flow.joint
%
% Author: Diellor Basha, 2026

    p = inputParser;  p.KeepUnmatched = true;
    p.addParameter('lambda', []);
    p.addParameter('f', []);
    p.addParameter('fs', []);
    p.addParameter('nT', []);
    p.addParameter('nfft', []);
    p.addParameter('bands', []);
    p.addParameter('half', 'positive');
    p.addParameter('boundary', 'ring');
    p.addParameter('t0', 0);
    p.addParameter('nV', []);
    p.parse(varargin{:});
    o = p.Results;

    ax = struct();
    if ~isempty(o.lambda)
        ax.lambda = double(o.lambda(:));
        ax.k      = sqrt(ax.lambda);
        ax.nK     = numel(ax.lambda);
        ax.lmax   = max(ax.lambda);
    end
    if ~isempty(o.f)
        ax.f      = double(o.f(:)).';
        ax.omega  = 2*pi*ax.f;
        ax.nOmega = numel(ax.f);
        if numel(ax.f) > 1, ax.df = median(diff(sort(ax.f))); else, ax.df = NaN; end
    end
    ax.fs = o.fs;  ax.nT = o.nT;  ax.nfft = o.nfft;
    if ~isempty(o.fs) && ~isempty(o.nT), ax.dur = o.nT / o.fs; else, ax.dur = []; end
    ax.bands = o.bands;  ax.half = o.half;  ax.boundary = o.boundary;
    ax.t0 = o.t0;  ax.nV = o.nV;
end

% Author: Diellor Basha, 2026
