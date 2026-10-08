function varargout = fbd_selftest(what, varargin)
% DEMOS.FBD_SELFTEST  Reach the private demo helpers from the test suite.
%
%   rheome.demos.fbd_selftest('check',  name, measured, expected, tol, unit)
%   rheome.demos.fbd_selftest('report', checks)
%   rheome.demos.fbd_selftest('ternary', cond, a, b)
%   rheome.demos.fbd_selftest('fig',    doExport, pos)
%   rheome.demos.fbd_selftest('sphere',   nSub, K)
%   rheome.demos.fbd_selftest('bump',     S, seedIdx, sigma)
%   rheome.demos.fbd_selftest('sectoral', S, l, phase)
%   rheome.demos.fbd_selftest('rotating', S, degrees, Omega, tvec)
%
% +demos/private is invisible outside the package, so this thin forwarder exists to let
% tests exercise the helpers directly rather than only through a whole demo.
%
% Author: Diellor Basha, 2026

    switch lower(what)
        case 'check',   varargout{1} = fbd_check(varargin{:});
        case 'report',  [varargout{1}, varargout{2}] = fbd_report(varargin{:});
        case 'ternary', varargout{1} = fbd_ternary(varargin{:});
        case 'fig',     varargout{1} = fbd_fig(varargin{:});
        case 'sphere',   varargout{1} = fbd_sphere(varargin{:});
        case 'bump',     varargout{1} = fbd_bump(varargin{:});
        case 'sectoral', varargout{1} = fbd_sectoral(varargin{:});
        case 'rotating', [varargout{1}, varargout{2}] = fbd_rotating(varargin{:});
        otherwise, error('demos:fbd_selftest', 'unknown target ''%s''.', what);
    end
end

% Author: Diellor Basha, 2026
