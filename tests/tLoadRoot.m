classdef tLoadRoot < matlab.unittest.TestCase
% rheome.load.root says why the +data cache is unreachable (drive unmounted, link dangling), and does not
% object to a cache folder that merely does not exist yet.
    methods (Test)
        function anUnmountedDriveIsNamed(tc)
            import matlab.unittest.fixtures.EnvironmentVariableFixture
            tc.applyFixture(EnvironmentVariableFixture('RHEOME_DATA', '/Volumes/rheome_no_such_drive/rheome-data'));
            tc.verifyError(@() rheome.load.root(), 'load:noroot');
            try, rheome.load.root(); catch err, end
            tc.verifySubstring(err.message, 'not mounted');
        end

        function aDanglingLinkIsNamed(tc)
            import matlab.unittest.fixtures.*
            f = tc.applyFixture(TemporaryFolderFixture);
            link = fullfile(f.Folder, 'data');
            [s, msg] = system(sprintf('ln -s "%s" "%s"', fullfile(f.Folder, 'gone'), link));
            tc.assumeEqual(s, 0, msg);                                       % no symlinks on this OS
            tc.applyFixture(EnvironmentVariableFixture('RHEOME_DATA', link));
            tc.verifyError(@() rheome.load.root(), 'load:noroot');
        end

        function aFreshFolderIsFine(tc)
            import matlab.unittest.fixtures.*
            f = tc.applyFixture(TemporaryFolderFixture);
            d = fullfile(f.Folder, 'not-yet');
            tc.applyFixture(EnvironmentVariableFixture('RHEOME_DATA', d));
            tc.verifyEqual(rheome.load.root(), d);
        end
    end
end
