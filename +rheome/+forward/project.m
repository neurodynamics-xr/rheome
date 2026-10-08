function [Jproj, C] = project(J, dbasis)
% FORWARD.PROJECT  Band-limit a per-vertex 3-vector field onto the Dirac eigenbasis.
%
%   Jproj      = rheome.forward.project(J, dbasis)
%   [Jproj, C] = rheome.forward.project(J, dbasis)
%
% Analysis then synthesis: embeds the current as a pure-imaginary quaternion field,
% takes its Dirac-mode coefficients under the mass (L²) inner product, and reconstructs.
% The result is the B-orthogonal projection of J onto the span of the kept modes -- i.e.
% the SAME field expressed in the Dirac eigenbasis (smoothed / band-limited, per-vertex
% scale preserved). Use it to represent an external reconstruction (e.g. a Brainstorm MNE
% current) in the module's eigenbasis, isolating the basis effect from the inverse operator.
%
% INPUTS:
%   J       [3nV x m] per-vertex 3-vectors, rows [x1,y1,z1, x2,...]  (m columns = frames)
%   dbasis  Dirac eigenbasis from rheome.eigen.dirac_frame (.Phi, .Mass, .nVert, .nModes)
%
% OUTPUTS:
%   Jproj   [3nV x m] the band-limited field (Phi * C, quaternion vector part)
%   C       [nModes x m] the mode coefficients (analysis: C = Phi' * B * embed(J))
%
% See also: rheome.forward.dirac, rheome.forward.reconstruct
%
% Author: Diellor Basha, 2026

    nV = dbasis.nVert;
    if size(J, 1) ~= 3*nV
        error('forward:project:size', 'J has %d rows but the basis expects 3*nV = %d.', size(J,1), 3*nV);
    end
    Jq = zeros(4*nV, size(J, 2));            % embed as pure-imaginary quaternion (w rows = 0)
    Jq(2:4:end, :) = J(1:3:end, :);          % x -> i
    Jq(3:4:end, :) = J(2:3:end, :);          % y -> j
    Jq(4:4:end, :) = J(3:3:end, :);          % z -> k

    C     = dbasis.Phi' * (dbasis.Mass * Jq);   % analysis (mass inner product)  [nModes x m]
    Jproj = rheome.forward.reconstruct(C, dbasis);     % synthesis (band-limited)        [3nV x m]
end

% Author: Diellor Basha, 2026
