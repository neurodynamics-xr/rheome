classdef tTypedFields < matlab.unittest.TestCase
% Domains, field types and operator arrows: the registries, and the code checked against them.
%
% The point of a registry is that it is not prose. docs/OPERATOR-REGISTRY.md is canonical
% for the NAMES; rheome.operators.registry is canonical for the ARROWS, and the case that earns
% its keep is theDeclaredShapesAreTheOnesTheCodeBuilds -- every mesh operator is built on a
% small icosphere and its size compared with the shape the registry declares. If an
% operator changes and the table does not, that fails here rather than three functions
% downstream.
%
% Author: Diellor Basha, 2026

    properties
        V; F; d
    end

    methods (TestClassSetup)
        function mesh(tc)
            [tc.V, tc.F] = rheome.geom.icosphere(2);                       % 162 vertices, 320 faces
            tc.d = rheome.domain.of(struct('Vertices', tc.V, 'Faces', tc.F), Name="ico2");
        end
    end

    methods (Test)

        % ---------------- domains ----------------

        function theDomainKindsDeclareTheirRanks(tc)
            K = rheome.domain.kinds();
            % ⭐ 'product' was added when joint manifolds got a type (rheome.domain.product). It mirrors its
            % first factor's counts, so every existing shape lambda is unaffected -- asserted in
            % tSelectionRegistry/theProductKindExistsAndIsAdditive.
            tc.verifyEqual(sort(K.kind), ["complex"; "graph"; "index"; "points"; "product"]);
            tc.verifyEqual(K.ranks{K.kind == "complex"}, {'vertex','edge','face'});
            tc.verifyEqual(K.ranks{K.kind == "graph"}, {'vertex','edge'});
            tc.verifyEqual(K.ranks{K.kind == "index"}, {'element'});
            tc.verifyEqual(K.ranks{K.kind == "product"}, {'vertex','edge','face','element'});
        end

        function aMeshDescriptorCountsEveryRankAndSatisfiesEuler(tc)
            tc.verifyEqual(tc.d.kind, 'complex');
            tc.verifyEqual(tc.d.nV, size(tc.V, 1));
            tc.verifyEqual(tc.d.nF, size(tc.F, 1));
            tc.verifyEqual(tc.d.nV - tc.d.nE + tc.d.nF, 2);         % a closed sphere
            tc.verifyEqual(rheome.domain.elements(tc.d, 'vertex'), tc.d.nV);
            tc.verifyEqual(rheome.domain.elements(tc.d, 'edge'), tc.d.nE);
            tc.verifyEqual(rheome.domain.elements(tc.d, 'face'), tc.d.nF);
            tc.verifyError(@() rheome.domain.elements(rheome.domain.index(10), 'face'), 'domain:elements:rank');
        end

        function theIdIsOfTheTopologyNotTheGeometry(tc)
            % ⭐ Why: an inflated and a white surface share connectivity, so a field defined
            % on one IS defined on the other -- which is what the flow code relies on when it
            % computes on one mesh and draws on another.
            scaled = rheome.domain.of(struct('Vertices', tc.V * 7 + 1, 'Faces', tc.F));
            tc.verifyEqual(scaled.id, tc.d.id);
            [V2, F2] = rheome.geom.icosphere(3);
            other = rheome.domain.of(struct('Vertices', V2, 'Faces', F2));
            tc.verifyNotEqual(other.id, tc.d.id);
            perm = rheome.domain.of(struct('Vertices', tc.V, 'Faces', tc.F(:, [2 3 1])));
            tc.verifyEqual(perm.id, tc.d.id);                       % winding is not topology
            tc.verifyEqual(rheome.domain.of(tc.d).id, tc.d.id);            % idempotent
        end

        function anIndexDomainRemembersWhatItCameFrom(tc)
            m = rheome.domain.index(400, Name="lb_modes", Of=tc.d);
            tc.verifyEqual(m.kind, 'index');
            tc.verifyEqual(m.nV, 400);
            tc.verifyEqual(m.of, tc.d.id);                          % the surface behind the line
            tc.verifyNotEqual(m.id, rheome.domain.index(400, Name="time").id);
        end

        % ---------------- field types ----------------

        function theFieldRegistryIsWellFormed(tc)
            T = rheome.fieldtype.registry();
            tc.verifyEqual(numel(unique(T.id)), height(T));
            kinds = rheome.domain.kinds().kind;
            for i = 1:height(T)
                tc.verifyTrue(all(ismember(string(T.domain_kind{i}), kinds)), T.id(i));
                tc.verifyTrue(ismember(T.status(i), ["built","planned"]), T.id(i));
                tc.verifyTrue(ismember(T.bundle(i), ["scalar","tangent","ambient","immersion"]), T.id(i));
                tc.verifyTrue(ismember(T.rank(i), ["vertex","edge","face","element"]), T.id(i));
            end
            tc.verifyGreaterThan(sum(T.status == "built"), 5);
        end

        function componentsFollowTheBundle(tc)
            tc.verifyEqual(rheome.fieldtype.components("scalarVertex"), 1);
            tc.verifyEqual(rheome.fieldtype.components("tangentVertex"), 1);   % one complex coordinate
            tc.verifyEqual(rheome.fieldtype.components("ambientVertexWorld"), 3);
            tc.verifyEqual(rheome.fieldtype.components("immersionVertex"), 4);
            T = rheome.fieldtype.registry();
            tc.verifyEqual(unique(T.components(T.bundle == "ambient")), 3);
            tc.verifyEqual(unique(T.components(T.bundle == "immersion")), 4);
            tc.verifyError(@() rheome.fieldtype.describe("nosuchfield"), 'fieldtype:describe:unknown');
        end

        function aTypePlusADomainResolvesTheShapeCollisionThatShapeAloneCannot(tc)
            % ⚠ THE CASE THE REPO ALREADY DOCUMENTS. rheome.inverse.mne warns that a [nCh x 3nV]
            % unconstrained gain and a [nCh x nV] constrained one are indistinguishable by
            % shape whenever nV divides by 3, which is why it must be told the vertex count.
            % A type plus a domain has no such ambiguity: the same 486 rows are a scalar
            % field on 486 vertices OR an ambient field on 162, and saying which is the point.
            small = tc.d;                                            % 162 vertices
            big = rheome.domain.of(struct('Vertices', randn(3*small.nV, 3), 'Faces', [1 2 3]));
            X = randn(3 * small.nV, 4);
            tc.verifyTrue(rheome.fieldtype.validate(X, "ambientVertexWorld", small, Throw=false));
            tc.verifyTrue(rheome.fieldtype.validate(X, "scalarVertex", big, Throw=false));
            tc.verifyFalse(rheome.fieldtype.validate(X, "scalarVertex", small, Throw=false));
            tc.verifyFalse(rheome.fieldtype.validate(X, "ambientVertexWorld", big, Throw=false));
            tc.verifyError(@() rheome.fieldtype.validate(X, "scalarVertex", small), 'fieldtype:validate:shape');
        end

        function aTypeIsRefusedOnADomainThatCannotCarryIt(tc)
            m = rheome.domain.index(400, Name="lb_modes", Of=tc.d);
            tc.verifyTrue(rheome.fieldtype.validate(randn(400, 3), "coeffScalar", m, Throw=false));
            [ok, why] = rheome.fieldtype.validate(randn(400, 3), "scalarVertex", m, Throw=false);
            tc.verifyFalse(ok);
            tc.verifySubstring(why, 'index');
            tc.verifyTrue(rheome.fieldtype.matches("scalarVertex", "scalarVertex"));
            tc.verifyFalse(rheome.fieldtype.matches("scalarVertex", "ambientVertexWorld"));
            tc.verifyFalse(rheome.fieldtype.matches("scalarVertex", "scalarFace"));   % same payload, other rank
        end

        function anAmbientFieldValidatesInEitherLayout(tc)
            % [nF x 3] and [3nF x 1] are the same face vector field written two ways, and the
            % repo uses both (rheome.flow.phasegradient returns [nF x 3], the weak operators eat
            % interleaved [3nV x 1]).
            tc.verifyTrue(rheome.fieldtype.validate(randn(tc.d.nF, 3), "ambientFaceWorld", tc.d, Throw=false));
            tc.verifyTrue(rheome.fieldtype.validate(randn(3*tc.d.nF, 1), "ambientFaceWorld", tc.d, Throw=false));
            tc.verifyFalse(rheome.fieldtype.validate(randn(tc.d.nF, 2), "ambientFaceWorld", tc.d, Throw=false));
        end

        % ---------------- operator arrows ----------------

        function theOperatorRegistryIsWellFormedAndTyped(tc)
            T = rheome.operators.registry();
            V = rheome.fieldtype.registry();
            tc.verifyEqual(numel(unique(T.id)), height(T));
            for i = 1:height(T)
                tc.verifyTrue(ismember(T.in(i), V.id), T.id(i));
                tc.verifyTrue(ismember(T.out(i), V.id), T.id(i));
                tc.verifyTrue(ismember(T.status(i), ["built","planned"]), T.id(i));
                tc.verifyTrue(ismember(T.dom(i), ["mesh","cross"]), T.id(i));
                tc.verifyNotEmpty(char(T.fcn(i)));
            end
            % every built operator names a function that exists. ⚠ exist() does not resolve
            % a package-qualified name; which() does.
            for i = find(T.status == "built" & ~contains(T.fcn, '('))'
                tc.verifyNotEmpty(which(char(T.fcn(i))), T.fcn(i));
            end
        end

        function theDeclaredShapesAreTheOnesTheCodeBuilds(tc)
            % ⭐ THE CASE THAT EARNS THE REGISTRY. Build each mesh operator for real and
            % compare its size with the arrow's declared shape.
            d = tc.d;  V = tc.V;  F = tc.F;                          %#ok<PROPLC>
            got = struct();
            got.mass = rheome.operators.mass(V, F);
            [got.laplace_beltrami, ~] = rheome.operators.laplace_beltrami(V, F);
            fg = rheome.operators.face_gradient(V, F);
            got.face_gradient = fg.Gx;                               % one ambient component
            got.face_average = fg.W;
            wd = rheome.operators.weak_differential(V, F, fg);
            got.weak_curl = wd.Curl;  got.weak_div = wd.Div;
            cl = rheome.operators.connection_laplacian(V, F);
            got.connection_laplacian = cl.A;                         % the builder returns a bundle
            got.dirac_extrinsic = rheome.operators.dirac_extrinsic(V, F);
            got.dirac_intrinsic_sq = rheome.operators.dirac_intrinsic_sq(V, F);
            got.dirac_frame = rheome.operators.dirac_frame(V, F, 0.5);
            T = rheome.operators.registry();
            names = fieldnames(got);
            for i = 1:numel(names)
                r = T(T.id == string(names{i}), :);
                tc.verifyEqual(height(r), 1, names{i});
                want = r.shape{1}(d);
                tc.verifyEqual(size(got.(names{i})), want, ...
                    sprintf('%s: built %s, registry says %s', names{i}, ...
                            mat2str(size(got.(names{i}))), mat2str(want)));
            end
            % ⭐ AND NOTHING IS LEFT OUT: every built mesh row that declares a shape must
            % have been built above, so adding an operator to the registry without building
            % it here fails rather than going unchecked.
            shaped = T.id(T.dom == "mesh" & T.status == "built" & ~cellfun(@isempty, T.shape));
            tc.verifyEqual(sort(string(names)), sort(shaped));
        end

        function anArrowAnswersWhatComesBackAndRefusesWhatDoesNot(tc)
            J = randn(3 * tc.d.nV, 6);
            tc.verifyEqual(rheome.operators.check("weak_curl", J, tc.d), "scalarVertex");
            tc.verifyEqual(rheome.operators.check("face_gradient", randn(tc.d.nV, 6), tc.d), "ambientFaceWorld");
            tc.verifyEqual(rheome.operators.check("face_average", randn(tc.d.nF, 6), tc.d), "scalarVertex");
            tc.verifyError(@() rheome.operators.check("weak_curl", randn(tc.d.nV, 6), tc.d), 'fieldtype:validate:shape');
            tc.verifyError(@() rheome.operators.check("nosuchop", [], tc.d), 'operators:check:unknown');
            tc.verifyError(@() rheome.operators.check("d0", [], tc.d), 'operators:check:planned');
        end

        function theChainFromSensorsToEigencoefficientsTypeChecks(tc)
            % ⭐ THE WHOLE POINT, as a composition: sensors -> current -> vorticity ->
            % coefficients, with every step's output the next step's declared input.
            chain = ["inverse_mne", "weak_curl", "lb_forward"];
            t = "sensorScalar";
            for i = 1:numel(chain)
                r = rheome.operators.registry();
                r = r(r.id == chain(i), :);
                tc.verifyEqual(r.in, t, sprintf('%s does not accept %s', chain(i), t));
                t = r.out;
            end
            tc.verifyEqual(t, "coeffScalar");
            % and the lift back, which is the only reason to leave the coefficients
            r = rheome.operators.registry();
            tc.verifyEqual(r.out(r.id == "lb_inverse"), "scalarVertex");
            tc.verifyEqual(r.in(r.id == "lb_inverse"), "coeffScalar");
        end

        function theFusedFlowOperatorsGoStraightToCoefficients(tc)
            % The registry records what the pipeline actually does: no per-vertex field is
            % materialised, sensors land on coefficients in one product.
            T = rheome.operators.registry();
            for id = ["flow_curl","flow_divergence","flow_potential","flow_stream"]
                r = T(T.id == id, :);
                tc.verifyEqual(r.in, "sensorScalar");
                tc.verifyEqual(r.out, "coeffScalar");
                tc.verifyEqual(r.dom, "cross");
            end
            tc.verifyEqual(T.out(T.id == "flow_field"), "ambientVertexWorld");   % the fallback
        end

    end
end

% Author: Diellor Basha, 2026
