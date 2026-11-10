classdef tNxrReaders < matlab.unittest.TestCase
% rheome.io.read.nxrregistration, and the schema-50 catalogue queries of nxrrecording / nxrkernel, on a synthetic nxr dataset
% (tests/nxr_store.py: dataset.sqlite from nxr's own model.sql, a subject and a template, one acquisition, the subject's stored
% registration hodgeDirac_to_<template>); the refusal of a catalogue that is not schema 50, with nxr's wording. The dataset is built by
% nxr: skipped when $NXR (the nxr checkout) or uv is absent. The compressed T needs numcodecs in MATLAB's Python ($RHEOME_PYTHON).
%
% Author: Diellor Basha, 2026

    properties
        plain, packed
    end

    methods (TestClassSetup)
        function build(tc)
            nxr = getenv('NXR');
            tc.assumeTrue(~isempty(nxr) && isfile(fullfile(nxr, 'pyproject.toml')) && system('command -v uv >/dev/null') == 0, 'nxr ($NXR) or uv not found');
            v = getenv('RHEOME_PYTHON');                               % before MATLAB's first Python call of the session
            if ~isempty(v) && pyenv().Status == "NotLoaded", pyenv(Version=v); end
            here = fileparts(mfilename('fullpath'));
            tc.plain = tempname; tc.packed = tempname;
            tc.addTeardown(@() cellfun(@(d) rmdir(d, 's'), {tc.plain, tc.packed}));
            for d = {tc.plain, ''; tc.packed, ' 1'}'
                [st, out] = system(sprintf('uv run --frozen --project "%s" python -I "%s" "%s"%s 2>&1', nxr, fullfile(here, 'nxr_store.py'), d{1}, d{2}));
                tc.assertEqual(st, 0, out);
            end
        end
    end

    methods (Test)
        function storedRegistrationIsRead(tc)
            G = rheome.io.read.nxrregistration(tc.plain, "subj");
            tc.verifyEqual(G.T, double(single(reshape(0:11, 4, 3)' / 7)), 'AbsTol', 1e-7);
            tc.verifyEqual([string(G.Template) string(G.Sampling) G.Flow], ["tmpl" "shepard" "fields"]);
            tc.verifyEqual([G.Neighbours G.Power], [6 3]);
            tc.verifyEqual([G.From.K G.To.K], [4 3]); tc.verifyEqual(G.To.mu, [1.5; 2.5; 3.5]); tc.verifyEqual(G.To.Kind, [0; 2; 0]);
            tc.verifyEqual([G.From.Mesh G.To.Package], ["m1" "bb"]);
            tc.verifyEqual(rheome.io.read.nxrregistration(tc.plain, "subj", "@tmpl").Name, 'hodgeDirac_to_tmpl');
            tc.verifyError(@() rheome.io.read.nxrregistration(tc.plain, "subj", "other"), 'io:read:nxr:notFound');
            tc.verifyError(@() rheome.io.read.nxrregistration(tc.plain, "@tmpl"), 'io:read:nxr:notFound');
        end

        function compressedRegistrationIsRead(tc)
            try, py.importlib.import_module('numcodecs'); ok = true; catch, ok = false; end
            tc.assumeTrue(ok, 'numcodecs is not importable in MATLAB''s Python (set RHEOME_PYTHON)');
            tc.verifyEqual(rheome.io.read.nxrregistration(tc.packed, "subj").T, rheome.io.read.nxrregistration(tc.plain, "subj").T);
        end

        function recordingAndKernelQueriesRunOnSchema50(tc)
            tc.verifyError(@() rheome.io.read.nxrrecording(tc.plain, "subj", "rest"), 'io:read:nxr:notFound');
            tc.verifyError(@() rheome.io.read.nxrkernel(tc.plain, "subj", "k"), 'io:read:nxr:notFound');
        end

        function otherSchemasAreRefusedWithNxrsWords(tc)
            d = tempname; copyfile(tc.plain, d); tc.addTeardown(@() rmdir(d, 's'));
            for v = [49 48]
                [st, out] = system(sprintf('sqlite3 "%s" "PRAGMA user_version = %d"', fullfile(d, 'dataset.sqlite'), v));
                tc.assumeEqual(st, 0, out);
                try
                    rheome.io.read.nxrregistration(d, "subj"); tc.verifyFail('read a schema ' + string(v));
                catch e
                    tc.verifyEqual(e.identifier, 'io:read:nxr:schema');
                    tc.verifySubstring(e.message, sprintf('schema version %d', v));
                    if v == 49, tc.verifySubstring(e.message, 'nxr migrate --dataset <dataset folder>'); else, tc.verifySubstring(e.message, 're-convert it from its source'); end
                end
            end
        end
    end
end
