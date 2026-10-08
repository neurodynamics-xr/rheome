classdef tSelectionRegistry < matlab.unittest.TestCase
% The fourth pillar and the convention axis: rheome.selection.registry, and rheome.operators.registry's
% direction/convention/adjoint columns.
%
% Author: Diellor Basha, 2026
    properties
        S; O
    end
    methods (TestClassSetup)
        function build(t)
            t.S = rheome.selection.registry();
            t.O = rheome.operators.registry();
        end
    end
    methods (Test)
        % ---------- the selection registry ----------
        function everySelectionNamesAResolvableFunction(t)
            % the same contract rheome.operators.registry enforces, and the reason rheome.forward.leadfield exists
            for i = 1:height(t.S)
                t.verifyNotEmpty(which(t.S.fcn(i)), t.S.fcn(i));
            end
        end
        function onlyPartitionsMerge(t)
            % ⭐ the load-bearing invariant: merges <=> partition
            t.verifyEqual(t.S.merges, t.S.kind == "partition", ...
                'a kernel or a subset must never be marked mergeable');
        end
        function everyDomainKindIsInTheControlledVocabulary(t)
            t.verifyTrue(all(ismember(t.S.domain_kind, rheome.domain.kinds().kind)));
        end
        function kindsAreFromTheClosedSet(t)
            t.verifyTrue(all(ismember(t.S.kind, ["partition" "kernel" "subset"])));
        end
        function tightIsSetForKernelsAndNotForPartitions(t)
            t.verifyTrue(all(isnan(t.S.tight(t.S.kind ~= "kernel"))));
            t.verifyTrue(all(~isnan(t.S.tight(t.S.kind == "kernel"))));
        end
        function bothConjugatePairsArePresent(t)
            % ⭐ the payoff: rheome.select.ladder and rheome.geom.ladder pair a partition with its conjugate, and
            % both axes of both pairs must be registered for that to be one function
            for a = ["time" "frequency" "cortex" "eigenmode"]
                onAxis = cellfun(@(c) any(strcmp(c, a)), t.S.axes);
                t.verifyTrue(any(onAxis & t.S.merges), ...
                    sprintf('no mergeable selection on the %s axis', a));
            end
        end
        function aProductSelectionIsTypedNotLabelled(t)
            % ⭐ numel(axes) == 2 <=> domain_kind == 'product'. joint_wavelet is the separable
            % psi_g(lambda)*psi_t(omega) and is the only product-valued row; it was first typed as a
            % plain eigenmode selection, which is the error this invariant forecloses.
            n = cellfun(@numel, t.S.axes);
            t.verifyEqual(n == 2, t.S.domain_kind == "product");
            j = t.S.id == "joint_wavelet";
            t.verifyEqual(t.S.axes{j}, {'eigenmode','frequency'});
        end
        function theProductKindExistsAndIsAdditive(t)
            % ⭐ the whole reason the kind was cheap: every existing shape lambda must give the SAME
            % answer on a product as on its first factor, so no arrow needs migrating.
            t.verifyTrue(ismember("product", rheome.domain.kinds().kind));
            dA = rheome.domain.of(struct('Vertices', zeros(50,3), 'Faces', [1 2 3]));
            d  = rheome.domain.product(dA, rheome.domain.index(7, Name="time"));
            O  = rheome.operators.registry();  sh = O.shape(~cellfun(@isempty, O.shape));
            for i = 1:numel(sh)
                t.verifyEqual(sh{i}(d), sh{i}(dA), sprintf('shape lambda %d moved', i));
            end
            t.verifyEqual(rheome.domain.elements(d, 'element'), 50*7);
            t.verifyEqual(rheome.domain.elements(d, 'vertex'), 50);
        end
        function domainSameCatchesAConventionMismatch(t)
            % ⚠ the check @jointfilterbank/iscompatible calls "silent and fatal", now in the type
            dM = rheome.domain.index(9, Name="modes");  dF = rheome.domain.index(5, Name="freq");
            a = rheome.domain.product(dM, dF, Conventions=struct('half','positive'));
            b = rheome.domain.product(dM, dF, Conventions=struct('half','whole'));
            t.verifyTrue(rheome.domain.same(a, a));
            [ok, why] = rheome.domain.same(a, b);
            t.verifyFalse(ok);
            t.verifyTrue(contains(why, "half"));
        end
        function aProductOfTwoComplexesIsRefused(t)
            dA = rheome.domain.of(struct('Vertices', zeros(9,3), 'Faces', [1 2 3]));
            t.verifyError(@() rheome.domain.product(dA, dA), 'domain:product:second');
        end
        function theFrameCheckIsOptInAndWorks(t)
            % ⚠ the hole that let a wrong 'js' frame count through to the arithmetic
            dA = rheome.domain.of(struct('Vertices', zeros(6,3), 'Faces', [1 2 3]));
            X = zeros(6, 4);
            t.verifyTrue(rheome.fieldtype.validate(X, "scalarVertex", dA, Throw=false));
            t.verifyTrue(rheome.fieldtype.validate(X, "scalarVertex", dA, Throw=false, Frames=4));
            [ok, why] = rheome.fieldtype.validate(X, "scalarVertex", dA, Throw=false, Frames=5);
            t.verifyFalse(ok);
            t.verifyTrue(contains(why, "frames"));
        end
        function checkRefusesMergingOnAKernel(t)
            [ok, why] = rheome.selection.check("graph_wavelet", struct('kind','index'), Merging=true);
            t.verifyFalse(ok);
            t.verifyTrue(contains(why, "kernel"));
            t.verifyTrue(rheome.selection.check("cortex_tile", struct('kind','complex'), Merging=true));
        end
        function checkRefusesAMismatchedDomain(t)
            [ok, why] = rheome.selection.check("cortex_tile", struct('kind','index'));
            t.verifyFalse(ok);
            t.verifyTrue(contains(why, "complex"));
        end

        % ---------- direction and convention ----------
        function everyOperatorHasADirection(t)
            t.verifyTrue(all(ismember(t.O.direction, ["forward" "inverse" "endo"])));
        end
        function twoOperatorsWithTheSameArrowMustDifferInConvention(t)
            % ⚠⚠ THE REGRESSION THIS WHOLE COLUMN EXISTS FOR. forward_dirac and leadfield_dirac are
            % both coeffCurrent -> sensorScalar and differ by a mass factor of 1.2e5; swapping them
            % type-checks, runs and is 100+ dB wrong. A direction flag alone cannot separate them.
            % ⚠ only among rows that CARRY a convention. Four flow_* rows share
            % sensorScalar -> coeffScalar legitimately: they are different computations, not an
            % analysis/synthesis ambiguity, and requiring them to differ in convention is wrong.
            cross = t.O(t.O.dom == "cross" & t.O.convention ~= "", :);
            [~, ~, g] = unique(cross(:, {'in','out'}), 'rows');
            for k = unique(g)'
                rows = cross(g == k, :);
                if height(rows) > 1
                    t.verifyEqual(numel(unique(rows.convention)), height(rows), ...
                        sprintf('rows %s share an arrow without distinct conventions', ...
                                strjoin(rows.id', ', ')));
                end
            end
        end
        function sameArrowSameDirectionIsSeparatedByConvention(t)
            % ⭐ the pair that motivated all of this: identical arrow, identical direction, not
            % interchangeable, separated only by convention. And NOT adjoints -- an adjoint swaps
            % in and out, and these do not.
            a = t.O(t.O.id == "leadfield_dirac", :);
            b = t.O(t.O.id == "forward_dirac", :);
            t.verifyEqual(a.in, b.in);       t.verifyEqual(a.out, b.out);
            t.verifyEqual(a.direction, b.direction);
            t.verifyNotEqual(a.convention, b.convention);
            t.verifyEqual(a.adjoint, "");    t.verifyEqual(b.adjoint, "");
        end
        function adjointPairingIsSymmetricAndSwapsTheArrow(t)
            for i = 1:height(t.O)
                a = t.O.adjoint(i);
                if a == "", continue; end
                j = find(t.O.id == a);
                t.verifyNumElements(j, 1, sprintf('%s names a missing adjoint %s', t.O.id(i), a));
                t.verifyEqual(t.O.adjoint(j), t.O.id(i), 'the pairing must be symmetric');
                t.verifyEqual(t.O.in(j),  t.O.out(i), 'an adjoint must swap in and out');
                t.verifyEqual(t.O.out(j), t.O.in(i),  'an adjoint must swap in and out');
                t.verifyNotEqual(t.O.convention(j), t.O.convention(i), ...
                    'an adjoint pair is one analysis and one synthesis');
            end
        end
        function aRegularisedInverseHasNoAdjoint(t)
            % ⚠ inverse_mne is a regularised pseudo-inverse, NOT the adjoint of the leadfield
            t.verifyEqual(t.O.adjoint(t.O.id == "inverse_mne"), "");
            t.verifyEqual(t.O.direction(t.O.id == "inverse_mne"), "inverse");
        end
        function theForwardDirectionOfSourceMappingExists(t)
            % it was absent from the table while every inverse was present
            t.verifyTrue(any(t.O.in == "ambientVertexWorld" & t.O.out == "sensorScalar"));
            t.verifyEqual(t.O.direction(t.O.id == "leadfield"), "forward");
        end
    end
end
