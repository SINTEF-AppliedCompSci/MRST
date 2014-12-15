function suite = getIntegrationTestSuiteMRST()
    mrstModule add ad-unittest
    import matlab.unittest.TestSuite;
    adpath = mrstPath('query', 'ad-unittest');
    p = fullfile(adpath, 'test_sim');
    suite = matlab.unittest.TestSuite.fromFolder(p);
end