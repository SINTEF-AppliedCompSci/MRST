function suite = getCoreUnitTestSuiteMRST()
%Undocumented utility for collecting core unit tests.

    import matlab.unittest.TestSuite;
    p = fullfile(fileparts(mfilename('fullpath')), 'unitTests');
    suite = TestSuite.fromFolder(p);
end
