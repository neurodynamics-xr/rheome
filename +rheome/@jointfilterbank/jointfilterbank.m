classdef jointfilterbank
% JOINTFILTERBANK  Filterbank on the joint time-vertex plane (lambda, omega).
%
%   jfb = rheome.jointfilterbank(gfb,  'SignalLength', N, ...)   % gfb a graphfilterbank
%   jfb = rheome.jointfilterbank(lmax, 'SignalLength', N, ...)   % builds the graph bank for you
%
% Every combination of three factor lists -- GSPBox's dgw composition:
%
%   W_i(lambda,omega) = K_k(lambda,omega) * psi_g,g(lambda) * psi_t,t(omega)
%   Nf = Ng * Nt * Nk
%
% ⭐ NOTHING IS MATERIALISED. The class stores HANDLES and an index, never products. A
% separable member is a rank-1 outer product: storing its factors costs [K]+[nOmega]
% against [K x nOmega] materialised, a 780x saving at cortex scale. jointfilters(jfb)
% therefore materialises only on explicit request and refuses above MaxBytes.
%
% REQUIRED: the graph side AND the temporal grid -- neither derives from the other. This
% is graphfilterbank's required argument plus cwtfilterbank's.
%
% See also: rheome.graphfilterbank, cwtfilterbank, rheome.graphtransform.eigen
%
% Author: Diellor Basha, 2026

    properties
        GraphBank                              % a graphfilterbank; supplies psi_g
        SamplingFrequency (1,1) double {mustBePositive} = 1
        Half              (1,:) char = 'positive'
        Boundary          (1,:) char = 'ring'
        Transform                    = []
        MaxBytes          (1,1) double {mustBePositive} = 2e9
        Bands                        = []      % [nb x 2] Hz
        TimeFilters                  = {}      % {1 x Nt} @(omega)
        JointKernels                 = {}      % {1 x Nk} @(lambda, omega)
        Labels                       = {}
    end

    properties (SetAccess = private, Hidden)
        PsiT_  = {}                            % resolved time factors
        Kern_  = {}                            % resolved joint kernels
        Index_ = []                            % [Nf x 3] (ig, it, ik)
        Sep_   = true
    end

    properties (SetAccess = immutable)
        SignalLength                           % [] if Frequencies was given directly
        Frequencies                            % [1 x nOmega] Hz
    end

    properties (Dependent, SetAccess = private)
        Lambda                                 % [K x 1] from the graph bank
        Omega                                  % [1 x nOmega] rad/s
        NumOmega
        Df
        NumMembers
        NumGraph
        NumTime
        NumJoint
        Index
        Separable
        TimeMembers                            % resolved {1 x Nt} @(omega) -- Bands + TimeFilters
        KernelMembers                          % resolved {1 x Nk} @(lambda,omega)
    end

    methods
        function obj = jointfilterbank(graphSide, varargin)
            if nargin < 1
                error('jointfilterbank:graphSide', ...
                    'A rheome.graphfilterbank or a lambda_max is required.');
            end
            p = inputParser;
            p.addParameter('SignalLength',      []);
            p.addParameter('Frequencies',       []);
            p.addParameter('SamplingFrequency', 1);
            p.addParameter('Half',              'positive');
            p.addParameter('Boundary',          'ring');
            p.addParameter('Transform',         []);
            p.addParameter('MaxBytes',          2e9);
            p.addParameter('Bands',             []);
            p.addParameter('TimeFilters',       {});
            p.addParameter('JointKernels',      {});
            p.addParameter('Labels',            {});
            p.parse(varargin{:});
            o = p.Results;

            if isa(graphSide, 'rheome.graphfilterbank')
                obj.GraphBank = graphSide;
            else
                obj.GraphBank = rheome.graphfilterbank(graphSide);
            end

            obj.SamplingFrequency = o.SamplingFrequency;
            obj.Half              = o.Half;
            obj.Boundary          = o.Boundary;
            obj.MaxBytes          = o.MaxBytes;

            hasN = ~isempty(o.SignalLength);
            hasF = ~isempty(o.Frequencies);
            if hasN == hasF
                error('jointfilterbank:temporalGrid', ...
                    ['Exactly one of ''SignalLength'' or ''Frequencies'' is required. ' ...
                     'SignalLength derives the DFT bin grid; Frequencies gives it ' ...
                     'explicitly in Hz (a band-limited analysis retains only a slice).']);
            end
            if hasN
                N = o.SignalLength;
                if ~isscalar(N) || N < 2 || mod(N,1) ~= 0
                    error('jointfilterbank:temporalGrid', ...
                        'SignalLength must be an integer >= 2, got %s.', mat2str(N));
                end
                obj.SignalLength = N;
                switch lower(o.Half)
                    case 'positive', bins = 0:floor(N/2);
                    case 'full',     bins = 0:N-1;
                    otherwise
                        error('jointfilterbank:half', ...
                            'Half must be ''positive'' or ''full'', got ''%s''.', o.Half);
                end
                obj.Frequencies = bins * (o.SamplingFrequency / N);
            else
                obj.SignalLength = [];
                obj.Frequencies  = double(o.Frequencies(:)).';
            end

            % Inherit the graph bank's transform when none is given: under composition
            % the field-domain route is a property of the graph side, and requiring it
            % twice is a papercut with no upside.
            if isempty(o.Transform) && ~isempty(obj.GraphBank.Transform)
                obj.Transform = obj.GraphBank.Transform;
            else
                obj.Transform = o.Transform;
            end

            obj.Bands        = o.Bands;
            obj.TimeFilters  = o.TimeFilters;
            obj.JointKernels = o.JointKernels;
            obj = rebuild(obj);
            if ~isempty(o.Labels)
                if numel(o.Labels) ~= obj.NumMembers
                    error('jointfilterbank:labels', ...
                        'Labels has %d entries but the bank has %d members.', ...
                        numel(o.Labels), obj.NumMembers);
                end
                obj.Labels = o.Labels(:).';
            end
        end

        function obj = set.Half(obj, v)
            v = lower(v);
            if ~any(strcmp(v, {'positive','full'}))
                error('jointfilterbank:half', ...
                    'Half must be ''positive'' or ''full'', got ''%s''.', v);
            end
            obj.Half = v;
        end

        function obj = set.Boundary(obj, v)
            v = lower(v);
            if ~any(strcmp(v, {'ring','path'}))
                error('jointfilterbank:boundary', ...
                    'Boundary must be ''ring'' or ''path'', got ''%s''.', v);
            end
            obj.Boundary = v;
        end

        function obj = set.Transform(obj, T)
            if ~isempty(T)
                [ok, why] = rheome.graphtransform.validate(T);
                if ~ok
                    error('jointfilterbank:transform', ...
                        'Transform is not usable: %s', strjoin(why, '; '));
                end
            end
            obj.Transform = T;
        end

        function v = get.Lambda(obj)
            g = obj.GraphBank;
            if g.HasSpectrum
                v = g.Lambda_;
            elseif ~isempty(obj.Transform) && ~isempty(obj.Transform.lambda)
                v = obj.Transform.lambda;
            else
                v = linspace(0, g.SpectralRange, 512).';
            end
            v = double(v(:));
        end

        function v = get.NumMembers(obj), v = size(obj.Index_, 1);        end
        function v = get.NumGraph(obj),   v = obj.GraphBank.NumMembers;  end
        function v = get.NumTime(obj),    v = numel(obj.PsiT_);          end
        function v = get.NumJoint(obj),   v = numel(obj.Kern_);          end
        function v = get.Index(obj),      v = obj.Index_;                end
        function v = get.Separable(obj),  v = obj.Sep_;                  end
        function v = get.TimeMembers(obj),   v = obj.PsiT_;              end
        function v = get.KernelMembers(obj), v = obj.Kern_;              end

        function v = get.Omega(obj),    v = 2*pi*obj.Frequencies;     end
        function v = get.NumOmega(obj), v = numel(obj.Frequencies);   end

        function v = get.Df(obj)
            if obj.NumOmega > 1, v = median(diff(sort(obj.Frequencies)));
            else,                v = NaN;
            end
        end
    end

    methods (Access = private)
        function obj = rebuild(obj)
            % Resolve the time factors: explicit handles plus one raised-cosine per band.
            psiT = obj.TimeFilters;
            if isa(psiT, 'function_handle'), psiT = {psiT}; end
            bandsT = {};
            if ~isempty(obj.Bands)
                B = obj.Bands;
                if size(B,2) ~= 2
                    error('jointfilterbank:bands', 'Bands must be [nb x 2] in Hz.');
                end
                bandsT = arrayfun(@(i) rheome.jointfilterbank.band(B(i,1), B(i,2)), ...
                                  1:size(B,1), 'UniformOutput', false);
            end
            psiT = [bandsT, psiT];
            if isempty(psiT), psiT = {@(w) ones(size(w))}; end

            kern = obj.JointKernels;
            if isa(kern, 'function_handle'), kern = {kern}; end
            if isempty(kern), kern = {@(l,w) 1}; end

            obj.PsiT_ = psiT;
            obj.Kern_ = kern;

            Ng = obj.GraphBank.NumMembers;  Nt = numel(psiT);  Nk = numel(kern);
            [IG, IT, IK] = ndgrid(1:Ng, 1:Nt, 1:Nk);
            obj.Index_ = [IG(:), IT(:), IK(:)];

            % Separable iff every kernel is identically 1 on a probe of the real axes.
            lp = obj.Lambda;  wp = obj.Omega;
            lp = lp(round(linspace(1, numel(lp), min(8, numel(lp)))));
            wp = wp(round(linspace(1, numel(wp), min(8, numel(wp)))));
            obj.Sep_ = true;
            for ik = 1:Nk
                Kv = kern{ik}(lp(:), wp(:).');
                if ~all(Kv(:) == 1), obj.Sep_ = false; break; end
            end

            obj.Labels = arrayfun( ...
                @(m) sprintf('g%d/h%d/K%d', obj.Index_(m,1), obj.Index_(m,2), obj.Index_(m,3)), ...
                1:size(obj.Index_,1), 'UniformOutput', false);
        end
    end

    methods (Static)
        h = band(f1, f2, width)
        K = speedkernel(c, d)
    end
end

% Author: Diellor Basha, 2026
