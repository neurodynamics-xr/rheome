% FORWARD  Forward-model spectral transforms for Dirac source mapping.
%
% Change-of-basis between the vertex leadfield and the Dirac eigenbasis. Pure-MATLAB
% ports of Brainstorm's bst_dirac (TRANSFORM / RECONSTRUCT).
%
% ⚠⚠ FOUR OPERATORS, TWO DIRECTIONS, TWO CONVENTIONS -- and the pairs that look interchangeable are
% not. Read this table before chaining anything:
%
%   name                 arrow                                 convention   use it to
%   rheome.forward.leadfield    ambientVertexWorld -> sensorScalar    --           forward a current field
%   rheome.forward.diracgain    coeffCurrent       -> sensorScalar    SYNTHESIS    forward mode coefficients
%   rheome.forward.dirac        coeffCurrent       -> sensorScalar    ANALYSIS     build an inverse
%   rheome.forward.project      ambientVertexWorld -> coeffCurrent    analysis     band-limit a field
%   rheome.forward.reconstruct  coeffCurrent       -> ambientVertexWorld           make a field from modes
%
% ⚠ rheome.forward.diracgain and rheome.forward.dirac share their arrow and differ by a MASS WEIGHTING: measured on
% a reference subject the two differ by a factor of ~1.2e5, which is 1/(mean vertex area) = 9.1e4. Swapping
% them type-checks, runs, returns tesla, gives the right spatial pattern and is 100+ dB wrong.
%
% Functions:
%   rheome.forward.leadfield    - the MEG gain for the good channels, optionally one hemisphere
%   rheome.forward.diracgain    - the SYNTHESIS gain Gmode = G*Phi: one sensor pattern per eigenmode
%   rheome.forward.simulate     - seed a Dirac field at a vertex, forward it, add real empty-room
%                          background, and solve for the detection threshold in nAm
%   rheome.forward.dirac        - project an unconstrained leadfield onto the Dirac modes (ANALYSIS)
%   rheome.forward.reconstruct  - map Dirac mode coefficients back to per-vertex 3-vectors
%
% Author: Diellor Basha, 2026
