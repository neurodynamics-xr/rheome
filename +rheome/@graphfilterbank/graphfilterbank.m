classdef graphfilterbank
% GRAPHFILTERBANK  Spectral filterbank on a self-adjoint operator.
%
%   gfb = rheome.graphfilterbank(SpectralRange)
%   gfb = rheome.graphfilterbank(SpectralRange, 'Name', value, ...)
%
% A set of M spectral multipliers g_m(L) plus the frame structure S = sum |g_m|^2 they
% induce. Modelled on cwtfilterbank. Owns the SPECTRAL RANGE and the filter shapes; it
% does not own the operator, the eigenvectors, the graph, the vertex count, or the data.
%
% SpectralRange -- the SignalLength analogue -- in any of three forms:
%   lmax scalar     lambda_min derived as lambda_max/LPFactor
%   [lmin lmax]     explicit range
%   Lambda vector   takes max, AND retains the exact spectrum for exact frame bounds
%
% See also: cwtfilterbank, rheome.graphtransform.eigen, rheome.graphtransform.chebyshev
%
% Author: Diellor Basha, 2026

    properties
        Wavelet         (1,:) char   = 'mexhat'
        VoicesPerOctave (1,1) double {mustBeInteger, mustBePositive} = 3
        NumFilters                   = []
        LPFactor        (1,1) double {mustBePositive} = 20
        FineLimit       (1,:) char   = 'basis'
        Transform                    = []
    end

    properties (SetAccess = immutable)
        SpectralRange                 % lambda_max
    end

    properties (Dependent)
        SpectralLimits                % [lo hi] on the lambda axis
        ScaleLimits                   % [tmin tmax]
        SizeLimits                    % [smin smax], sigma = sqrt(2t)
    end

    properties (SetAccess = private, Hidden)
        Lmax_                         % lambda_max
        Lmin_                         % lambda_min as supplied ([] if only lmax given)
        Lambda_        = []           % the exact spectrum, [] unless a vector was passed
        Tmin_          = []           % canonical bank scale range
        Tmax_          = []
        LimitsSet_     = false        % true once a user set any limit view
        G_             = {}           % 1 x M gain handles
        T_             = []           % 1 x M scale params
    end

    properties (Dependent, SetAccess = private)
        HasSpectrum                   % true iff constructed from a Lambda vector
        NumMembers                    % M
        AchievedVoicesPerOctave       % what the member count actually delivers
        CenterWavenumbers             % [1 x M] rad/m
        MassLost                      % [1 x M] fraction of mass beyond lambda_max
        Usable                        % [1 x M] logical, MassLost <= 5%
        FinestScale                   % smallest usable sigma, sqrt(2*4.744/lmax)
    end

    methods
        function obj = graphfilterbank(range, varargin)
            if nargin < 1
                error('graphfilterbank:range', 'A spectral range is required.');
            end
            p = inputParser;
            p.addParameter('Wavelet',        'mexhat');
            p.addParameter('VoicesPerOctave', 3);
            p.addParameter('NumFilters',      []);
            p.addParameter('LPFactor',       20);
            p.addParameter('FineLimit',      'basis');
            p.addParameter('SpectralLimits', []);
            p.addParameter('ScaleLimits',    []);
            p.addParameter('SizeLimits',     []);
            p.addParameter('Transform',      []);
            p.parse(varargin{:});
            o = p.Results;

            [obj.Lmax_, obj.Lmin_, obj.Lambda_] = i_parserange(range);
            obj.SpectralRange = obj.Lmax_;
            obj.Wavelet   = o.Wavelet;
            obj.LPFactor  = o.LPFactor;
            obj.FineLimit = o.FineLimit;
            obj.VoicesPerOctave = o.VoicesPerOctave;
            obj.NumFilters      = o.NumFilters;

            nSet = ~isempty(o.SpectralLimits) + ~isempty(o.ScaleLimits) + ~isempty(o.SizeLimits);
            if nSet > 1
                error('graphfilterbank:limits', ...
                    ['SpectralLimits, ScaleLimits and SizeLimits are three views of ONE range. ' ...
                     'Set exactly one.']);
            end
            obj = i_defaultlimits(obj);
            if ~isempty(o.SpectralLimits), obj.SpectralLimits = o.SpectralLimits; end
            if ~isempty(o.ScaleLimits),    obj.ScaleLimits    = o.ScaleLimits;    end
            if ~isempty(o.SizeLimits),     obj.SizeLimits     = o.SizeLimits;     end
            obj = rebuild(obj);
            obj.Transform = o.Transform;
        end

        function obj = set.Transform(obj, T)
            if ~isempty(T)
                [ok, why] = rheome.graphtransform.validate(T);
                if ~ok
                    error('graphfilterbank:transform', ...
                        'Transform is not usable: %s', strjoin(why, '; '));
                end
            end
            obj.Transform = T;
        end

        function obj = set.Wavelet(obj, v)
            v = lower(v);
            if ~any(strcmp(v, {'mexhat','itersine','heat','logitersine'}))
                error('graphfilterbank:wavelet', ...
                    'Wavelet must be ''mexhat'', ''itersine'', ''heat'' or ''logitersine'', got ''%s''.', v);
            end
            obj.Wavelet = v;
            obj = rebuild(obj);
        end

        function obj = set.FineLimit(obj, v)
            v = lower(v);
            if ~any(strcmp(v, {'basis','usable'}))
                error('graphfilterbank:fineLimit', ...
                    'FineLimit must be ''basis'' or ''usable'', got ''%s''.', v);
            end
            obj.FineLimit = v;
            obj = i_defaultlimits(obj);
            obj = rebuild(obj);
        end

        function obj = set.VoicesPerOctave(obj, v)
            if v > 48
                error('graphfilterbank:voices', ...
                    'VoicesPerOctave must be in [1 48], got %g.', v);
            end
            obj.VoicesPerOctave = v;
            obj = rebuild(obj);
        end

        function obj = set.NumFilters(obj, v)
            if ~isempty(v) && (~isscalar(v) || v < 2 || mod(v,1) ~= 0)
                error('graphfilterbank:numFilters', ...
                    'NumFilters must be [] or an integer >= 2, got %s.', mat2str(v));
            end
            obj.NumFilters = v;
            obj = rebuild(obj);
        end

        function v = get.NumMembers(obj), v = numel(obj.G_); end

        function v = get.CenterWavenumbers(obj), v = centerWavenumbers(obj); end

        function v = get.MassLost(obj)
            % Fraction of each member's spectral mass falling beyond lambda_max. The
            % coi analogue: which coefficients do not mean what they look like.
            M = obj.NumMembers;  v = zeros(1, M);
            switch obj.Wavelet
                case 'mexhat'
                    u = obj.T_ * obj.Lmax_;
                    v = (1 + u) .* exp(-u);
                    % The scaling function is a quartic roll-off, not an exponential,
                    % so integrate it rather than reuse a closed form that no longer
                    % applies.
                    lq = linspace(0, 6*obj.Lmax_, 4096)';  gq = abs(obj.G_{1}(lq));
                    v(1) = trapz(lq(lq>obj.Lmax_), gq(lq>obj.Lmax_)) / max(trapz(lq, gq), eps);
                case 'heat'
                    v = exp(-obj.T_ * obj.Lmax_);
                case {'itersine', 'logitersine'}
                    v = zeros(1, M);            % compactly supported on [0 lmax]
            end
        end

        function v = get.Usable(obj),      v = obj.MassLost <= 0.05;    end
        function v = get.FinestScale(obj), v = sqrt(2*4.744/obj.Lmax_); end

        function v = get.AchievedVoicesPerOctave(obj)
            tf = obj.T_(~isnan(obj.T_));
            if numel(tf) < 2, v = NaN; return; end
            v = (numel(tf)-1) / max(log2(sqrt(max(tf)/min(tf))), eps);
        end

        function obj = set.LPFactor(obj, v)
            obj.LPFactor = v;
            obj = i_defaultlimits(obj);
            obj = rebuild(obj);
        end

        function v = get.ScaleLimits(obj),    v = [obj.Tmin_ obj.Tmax_];         end
        function v = get.SizeLimits(obj),     v = sqrt(2*[obj.Tmin_ obj.Tmax_]); end
        function v = get.SpectralLimits(obj), v = [2/obj.Tmax_, 1/obj.Tmin_];    end
        function v = get.HasSpectrum(obj),    v = ~isempty(obj.Lambda_);         end

        function obj = set.ScaleLimits(obj, v)
            i_checkpair(v, 'ScaleLimits');
            obj.Tmin_ = v(1);  obj.Tmax_ = v(2);  obj.LimitsSet_ = true;
            obj = rebuild(obj);
        end

        function obj = set.SizeLimits(obj, v)
            i_checkpair(v, 'SizeLimits');
            obj.ScaleLimits = v.^2 / 2;
        end

        function obj = set.SpectralLimits(obj, v)
            i_checkpair(v, 'SpectralLimits');
            obj.ScaleLimits = [1/v(2), 2/v(1)];
        end
    end

    methods (Static)
        gfb = fromOperator(L, varargin)
    end

    methods (Access = private)
        function obj = rebuild(obj)
            if isempty(obj.Lmax_) || isempty(obj.Tmax_) || isempty(obj.Wavelet)
                return;                        % partially constructed
            end
            lminEff = 2 / obj.Tmax_;
            d = designbank(obj.Wavelet, obj.NumFilters, obj.Tmin_, obj.Tmax_, ...
                           obj.Lmax_, lminEff, obj.VoicesPerOctave);
            obj.G_ = d.g;
            obj.T_ = d.t;

            % The guardrail is COMPUTED, not a hardcoded per-octave threshold, so it is
            % correct for every family and for any family added later. B/A would be the
            % wrong test: heat legitimately runs B/A ~ 56 (every member is ~1 at lambda=0)
            % while a perfectly healthy mexhat sits at 1.7-5.5. The real failure is
            % A -> 0, which is exactly where coverage fails and the dual stops existing.
            b = framebounds(obj);
            if b.A <= 1e-6 * b.B
                warning('graphfilterbank:degenerateFrame', ...
                    ['This bank does not cover its spectrum: A = %.3g, B = %.3g, %.0f%% of ' ...
                     'the axis is uncovered. The canonical dual does not exist there, so ' ...
                     'iwt will LOSE that content. Widen the scale limits or raise ' ...
                     'VoicesPerOctave (currently %d members, %.1f per octave).'], ...
                    b.A, b.B, 100*b.Uncovered, obj.NumMembers, obj.AchievedVoicesPerOctave);
            end
        end
    end
end

function [lmax, lmin, lambda] = i_parserange(range)
    range = double(range(:));
    lambda = [];
    if numel(range) > 2
        lambda = sort(range);
        lmax = lambda(end);
        nz = lambda(lambda > lmax*1e-12);
        if isempty(nz), lmin = []; else, lmin = nz(1); end
    elseif numel(range) == 2
        lmin = range(1);  lmax = range(2);
    else
        lmin = [];  lmax = range;
    end
    if ~isempty(lmin) && lmin < 0
        error('graphfilterbank:signedSpectrum', ...
            ['The spectral range [%g %g] is SIGNED. A first-order Dirac operator has a ' ...
             'two-sided spectrum: mexhat and heat diverge for lambda < 0, and sqrt(lambda) ' ...
             'is imaginary, so every scale reading would be undefined. Square the operator ' ...
             '(D^2 is PSD) or supply a family defined on a two-sided axis.'], lmin, lmax);
    end
    if ~isfinite(lmax) || lmax <= 0
        error('graphfilterbank:range', 'lambda_max must be finite and positive, got %g.', lmax);
    end
    if ~isempty(lmin) && lmax <= lmin
        error('graphfilterbank:range', 'invalid range [%g %g].', lmin, lmax);
    end
end

function obj = i_defaultlimits(obj)
    % The LPFactor default. This ties the coarse end to lambda_max -- set
    % SpectralLimits/ScaleLimits/SizeLimits explicitly to make the bank reproducible.
    if obj.LimitsSet_, return; end
    if isempty(obj.Lmax_), return; end
    % A supplied lambda_min of 0 carries no scale information (it is the DC mode), so
    % fall back to the LPFactor rule rather than producing t_max = Inf.
    if isempty(obj.Lmin_) || obj.Lmin_ <= 0
        lminEff = obj.Lmax_ / obj.LPFactor;
    else
        lminEff = obj.Lmin_;
    end
    switch lower(obj.FineLimit)
        case 'basis',  obj.Tmin_ = 1/obj.Lmax_;
        case 'usable', obj.Tmin_ = 4.744/obj.Lmax_;
    end
    obj.Tmax_ = 2 / lminEff;
end

function i_checkpair(v, name)
    if numel(v) ~= 2 || ~all(isfinite(v)) || ~(v(2) > v(1)) || any(v <= 0)
        error('graphfilterbank:limits', ...
            '%s must be [lo hi], finite, positive and increasing.', name);
    end
end

% Author: Diellor Basha, 2026
