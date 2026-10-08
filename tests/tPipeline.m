classdef tPipeline < matlab.unittest.TestCase
% A plan built from the registries, run against an environment, with its record.
%
% What is asserted is the three promises the layer makes: a chain that does not type-check
% is refused WHEN IT IS WRITTEN, a plan that cannot be bound is refused BEFORE IT RUNS, and
% a plan that does run gives exactly what the hand-written chain gives -- so the metadata is
% the algorithm rather than a description of one.
%
% Author: Diellor Basha, 2026

    properties
        V; F; S; e; nCh = 16; nT = 400; fs = 600; basis; K
    end

    methods (TestClassSetup)
        function build(tc)
            rng(11);
            [tc.V, tc.F] = rheome.geom.icosphere(2);
            tc.S = struct('Vertices', tc.V, 'Faces', tc.F, 'nV', size(tc.V,1), 'nF', size(tc.F,1));
            [L, M] = rheome.operators.laplace_beltrami(tc.V, tc.F);
            [Phi, Lam] = eigs(L, M, 12, 'smallestabs');
            tc.basis = struct('Phi', Phi, 'Lambda', diag(Lam), 'Mass', M);
            tc.K = randn(3*tc.S.nV, tc.nCh) * 1e-9;
            tc.e = rheome.pipeline.env(Surface=tc.S, Inverse=tc.K, Modes=tc.basis, ...
                                Sensors=randn(tc.nCh,3), Name="ico2");
        end
    end

    methods
        function p = chain(tc)
            p = rheome.pipeline.start("sensorScalar", tc.e.domains.sensors, ...
                               Frames=rheome.domain.index(tc.nT, Name="time", Units="s"), Name="curl coeffs");
            p = rheome.pipeline.window(p, [0.1 0.3], Rate=tc.fs);
            p = rheome.pipeline.then(p, "inverse_mne");
            p = rheome.pipeline.then(p, "curl");
            p = rheome.pipeline.then(p, "lb_forward");
        end
    end

    methods (Test)

        function theEnvironmentNamesTheDomainsItHolds(tc)
            tc.verifyTrue(all(isfield(tc.e.domains, {'source','sensors','modes'})));
            tc.verifyEqual(tc.e.domains.source.nV, tc.S.nV);
            tc.verifyEqual(tc.e.domains.sensors.kind, 'points');
            tc.verifyEqual(tc.e.domains.modes.kind, 'index');
            tc.verifyEqual(tc.e.domains.modes.nV, size(tc.basis.Phi, 2));
            tc.verifyEqual(tc.e.domains.modes.of, tc.e.domains.source.id);   % the line knows its surface
            tc.verifyNotEmpty(tc.e.fg);  tc.verifyNotEmpty(tc.e.wd);  tc.verifyNotEmpty(tc.e.L);
        end

        function aKernelOfTheWrongFibreIsRefusedWhenItIsBound(tc)
            % ⚠ The collision rheome.inverse.mne documents, caught at the environment instead.
            tc.verifyError(@() rheome.pipeline.env(Surface=tc.S, Inverse=randn(tc.S.nV, tc.nCh)), ...
                           'pipeline:env:inverse');
        end

        function aPlanIsMetadataAndReadsAsATable(tc)
            p = tc.chain();
            T = rheome.pipeline.describe(p);
            tc.verifyEqual(height(T), 4);
            tc.verifyEqual(T.operator', ["window","inverse_mne","curl","lb_forward"]);
            tc.verifyEqual(T.in(2), "sensorScalar");
            tc.verifyEqual(T.out(end), "coeffScalar");
            tc.verifyEqual(T.domain(end), "modes");           % resolved against an env at run time
            tc.verifySubstring(char(T.note(1)), 'samples 61:180');
            tc.verifyEqual(p.type, "coeffScalar");            % the state after the last step
        end

        function achainThatDoesNotTypeCheckIsRefusedWhereItIsWritten(tc)
            % ⭐ The promise: a mismatch surfaces at the line that wrote it, naming both
            % types and what would have been accepted.
            p = rheome.pipeline.start("sensorScalar", tc.e.domains.sensors);
            p = rheome.pipeline.then(p, "inverse_mne");
            tc.verifyError(@() rheome.pipeline.then(p, "lb_forward"), 'pipeline:then:arrow');
            try
                rheome.pipeline.then(p, "lb_forward");
            catch err
                tc.verifySubstring(err.message, 'ambientVertexWorld');
                tc.verifySubstring(err.message, 'curl');       % what does accept it
            end
            tc.verifyError(@() rheome.pipeline.then(p, "nosuchop"), 'pipeline:then:unknown');
            tc.verifyError(@() rheome.pipeline.then(p, "d0"), 'pipeline:then:planned');
            tc.verifyError(@() rheome.pipeline.start("scalarVertex", tc.e.domains.modes), 'pipeline:start:domain');
        end

        function aPlanThatCannotBeBoundIsRefusedBeforeItRuns(tc)
            p = tc.chain();
            bare = rheome.pipeline.env(Surface=tc.S, Name="no kernel");
            [f, why] = rheome.pipeline.compile(p, bare, Throw=false);
            tc.verifyEmpty(f);
            tc.verifyGreaterThanOrEqual(numel(why), 2);
            tc.verifyTrue(any(contains(why, 'no K')));
            tc.verifyTrue(any(contains(why, 'eigenbasis')));
            tc.verifyError(@() rheome.pipeline.compile(p, bare), 'pipeline:compile:unbound');
            tc.verifyNotEmpty(rheome.pipeline.compile(p, tc.e));
        end

        function runningThePlanEqualsTheHandWrittenChain(tc)
            % ⭐ The metadata IS the algorithm: no separate implementation to drift.
            p = tc.chain();
            X = randn(tc.nCh, tc.nT);
            [g, rec] = rheome.pipeline.run(p, tc.e, X);
            a = floor(0.1*tc.fs) + 1;  b = ceil(0.3*tc.fs);
            J = tc.K * X(:, a:b);
            v = rheome.differential.curl(J, tc.S, tc.e.fg);
            want = tc.basis.Phi' * (tc.basis.Mass * v);
            tc.verifyEqual(g, want);
            tc.verifySize(g, [size(tc.basis.Phi, 2), b - a + 1]);
            tc.verifyEqual(rec.out.type, "coeffScalar");
            tc.verifyEqual(rec.out.domain.id, tc.e.domains.modes.id);
            tc.verifyEqual(rec.out.frames.nV, b - a + 1);
        end

        function theRecordSaysHowTheAnswerWasProduced(tc)
            p = tc.chain();
            [~, rec] = rheome.pipeline.run(p, tc.e, randn(tc.nCh, tc.nT));
            R = rec.table;
            tc.verifyEqual(height(R), 4);
            tc.verifyEqual(R.rows', [tc.nCh, 3*tc.S.nV, tc.S.nV, size(tc.basis.Phi,2)]);
            tc.verifyTrue(all(R.cols == 120));
            tc.verifyEqual(R.out', ["sensorScalar","ambientVertexWorld","scalarVertex","coeffScalar"]);
            tc.verifyEqual(R.domain(2), "ico2/source");        % resolved, not 'source'
            tc.verifyEqual(R.domain(4), "ico2/lb_modes");
            tc.verifyTrue(all(R.seconds >= 0));
            tc.verifyEqual(rec.env, 'ico2');
            tc.verifyEqual(rec.in.type, "sensorScalar");
        end

        function everyIntermediateIsValidatedAsItGoes(tc)
            % A wrong applier must fail at the step that broke the contract.
            bad = tc.e;
            bad.K = tc.K(1:end-3, :);                          % one vertex short
            p = rheome.pipeline.start("sensorScalar", tc.e.domains.sensors);
            p = rheome.pipeline.then(p, "inverse_mne");
            tc.verifyError(@() rheome.pipeline.run(p, bad, randn(tc.nCh, 5)), 'fieldtype:validate:shape');
            try
                rheome.pipeline.run(p, bad, randn(tc.nCh, 5));
            catch err
                tc.verifySubstring(err.message, 'step 1');
                tc.verifySubstring(err.message, 'inverse_mne');
            end
        end

        function theInputItselfIsCheckedAgainstItsDeclaredType(tc)
            p = tc.chain();
            tc.verifyError(@() rheome.pipeline.run(p, tc.e, randn(tc.nCh + 5, tc.nT)), 'fieldtype:validate:shape');
        end

        function aWindowIsAStepAndTravelsIntoTheRecord(tc)
            p = rheome.pipeline.start("sensorScalar", tc.e.domains.sensors, Frames=rheome.domain.index(tc.nT, Name="time"));
            p = rheome.pipeline.window(p, [100 199], Units="samples", Note="a burst");
            tc.verifyEqual(p.frames.nV, 100);
            tc.verifyEqual(p.type, "sensorScalar");            % a selection changes no type
            [Y, rec] = rheome.pipeline.run(p, tc.e, randn(tc.nCh, tc.nT));
            tc.verifySize(Y, [tc.nCh, 100]);
            tc.verifySubstring(char(rec.table.note(1)), 'a burst');
            tc.verifyEqual(rec.out.frames.nV, 100);
            tc.verifyError(@() rheome.pipeline.window(p, [0.1 0.2]), 'pipeline:window:rate');
            tc.verifyError(@() rheome.pipeline.window(p, [50 10], Units="samples"), 'pipeline:window:span');
        end

        function thePlanRunsUnchangedOnAnotherEnvironment(tc)
            % ⭐ Why a plan is metadata: the same chain on another subject is the same value.
            p = tc.chain();
            [V2, F2] = rheome.geom.icosphere(2);
            S2 = struct('Vertices', V2 * 1.3, 'Faces', F2, 'nV', size(V2,1), 'nF', size(F2,1));
            [L2, M2] = rheome.operators.laplace_beltrami(S2.Vertices, S2.Faces);
            [Phi2, Lam2] = eigs(L2, M2, 12, 'smallestabs');
            e2 = rheome.pipeline.env(Surface=S2, Inverse=randn(3*S2.nV, tc.nCh)*1e-9, ...
                              Modes=struct('Phi', Phi2, 'Lambda', diag(Lam2), 'Mass', M2), ...
                              Sensors=randn(tc.nCh,3), Name="other");
            [g2, rec2] = rheome.pipeline.run(p, e2, randn(tc.nCh, tc.nT));
            tc.verifySize(g2, [12, 120]);
            tc.verifyEqual(rec2.env, 'other');
            tc.verifyNotEqual(rec2.out.domain.id, tc.e.domains.modes.id);
        end

    end
end

% Author: Diellor Basha, 2026
