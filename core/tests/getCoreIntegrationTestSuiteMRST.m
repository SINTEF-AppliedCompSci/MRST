function suite = getCoreIntegrationTestSuiteMRST()
%Undocumented utility for collecting core integration tests.

    import matlab.unittest.TestSuite;
    p = fullfile(fileparts(mfilename('fullpath')), 'integrationTests');
    suite = TestSuite.fromFolder(p);
end
