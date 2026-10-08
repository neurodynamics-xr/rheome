function varargout = sensor_dynamics_selftest(what, varargin)
% DEMOS.SENSOR_DYNAMICS_SELFTEST  Reach sensor_dynamics's local helpers from the test suite.
%
%   rheome.demos.sensor_dynamics_selftest('travelindex', z)
%   rheome.demos.sensor_dynamics_selftest('analytic', U)
%
% The travelling/standing discriminator is the demo's central claim, and it deserves tests
% of its own rather than only being exercised through a whole five-figure run.
%
% Author: Diellor Basha, 2026

    switch lower(what)
        case 'travelindex', varargout{1} = i_travelindex_pub(varargin{:});
        case 'analytic',    varargout{1} = i_analytic_pub(varargin{:});
        otherwise, error('demos:sensor_dynamics_selftest', 'unknown target ''%s''.', what);
    end
end

function z = i_analytic_pub(U)
% The FFT route to the analytic signal, mirroring rheome.demos.sensor_phase.
    U = U - mean(U, 2);
    nT = size(U,2);  F = fft(U, [], 2);  h = zeros(1, nT);  h(1) = 1;
    if mod(nT,2)==0, h(nT/2+1) = 1; h(2:nT/2) = 2; else, h(2:(nT+1)/2) = 2; end
    z = ifft(F .* h, [], 2);
end

function s = i_travelindex_pub(z)
    A = [real(z(:)), imag(z(:))];
    A = A - mean(A,1);
    e = sort(eig(cov(A)), 'descend');
    s = e(2) / max(sum(e), eps);
end

% Author: Diellor Basha, 2026
