function obs = observables(C, f, Psi, fsOut, dur, opts)
% FLOW.OBSERVABLES  The named observables of a flow map, from the joint spectrum.
%
%   obs = rheome.flow.observables(C, f, Psi, fsOut, dur)
%   obs = rheome.flow.observables(C, f, Psi, fsOut, dur, opts)
%
% Synthesises the COMPLEX flow map z(v,t) = Psi * IFFT(C) and returns the quantities that are
% immediately readable from it. These exist at the END OF FLOW MAPPING -- before any analysis
% choice -- which is why they are collected in one place rather than recomputed downstream.
%
%   vorticity  Re z         the signed flow map itself
%   chirality  sign(Re z)   +1 counter-clockwise, -1 clockwise (seen from outside)
%   amplitude  |z|          vortical activity, unsigned
%   phase      angle(z)     WHEN the reversal happens
%
% ON CHIRALITY. Both handednesses coexist in ONE signed map at every instant; there are not two
% maps. The sign of the vorticity IS the handedness and needs no demodulation. What the phase
% reference (rheome.flow.demod) adds is COMPARABILITY ACROSS LOCATIONS: with a common reference, equal
% phase means co-rotating and a difference of pi means counter-rotating. Without it every vertex
% carries its own w0*t ramp and the comparison is meaningless.
%
% MEMORY. z is [nVert x nT] COMPLEX -- 20k vertices at 60 Hz over 600 s is ~5.9 GB. Use opts.keep
% to materialise only what is needed, and prefer short windows chosen with rheome.show.activity (whose
% WHEN panel costs nothing, by Parseval).
%
% INPUTS:
%   C      [K x nOmega] joint spectrum (rheome.flow.joint; optionally through rheome.spectral.decompose /
%          rheome.flow.demod -- pass the shifted f if demodulated)
%   f      [1 x nOmega] bin frequencies (Hz)
%   Psi    [N x K] scalar eigenbasis
%   fsOut  output rate (Hz)      dur  duration (s)
%   opts   .keep  cellstr subset of {'z','vorticity','chirality','amplitude','phase'}
%                 (default: all except 'z')
%          .t0    start time (default 0)
%
% OUTPUT (struct obs): the requested fields, each [N x nT], plus
%   .t [1 x nT]   .ccwFraction  fraction of (vertex,time) samples with positive vorticity
%
% ⚠ .ccwFraction should sit near 0.5: conservation gives the closed-surface integral of curl as
% zero, so positive and negative vorticity must balance. A marked departure indicates a mesh or
% operator problem, not a finding.
%
% See also: rheome.flow.joint, rheome.flow.demod, rheome.flow.synth, rheome.show.activity, rheome.detect.chirality
%
% Author: Diellor Basha, 2026

    if nargin < 6, opts = struct(); end
    if ~isfield(opts,'keep') || isempty(opts.keep)
        opts.keep = {'vorticity','chirality','amplitude','phase'};
    end
    if ~isfield(opts,'t0') || isempty(opts.t0), opts.t0 = 0; end
    keep = @(n) any(strcmpi(opts.keep, n));

    [c, t] = rheome.flow.synth(C, f, fsOut, dur, struct('t0', opts.t0));
    z = Psi * c;                                   % [N x nT] complex -- the flow map

    obs.t = t;
    rz = real(z);
    obs.ccwFraction = mean(rz(:) > 0);
    if keep('z'),         obs.z         = z;        end
    if keep('vorticity'), obs.vorticity = rz;       end
    if keep('chirality'), obs.chirality = sign(rz); end
    if keep('amplitude'), obs.amplitude = abs(z);   end
    if keep('phase'),     obs.phase     = angle(z); end
end

% Author: Diellor Basha, 2026
