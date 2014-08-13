function suite = getUnitTestSuiteMRST()
    mrstModule add ad-unittest
    import matlab.unittest.TestSuite;
    adpath = mrstPath('query', 'ad-unittest');
    
    unitfolders = {'test_models', 'test_utils'};
    suite = [];
    for i = 1:numel(unitfolders)
        p = fullfile(adpath, unitfolders{i});
        suite = [suite, matlab.unittest.TestSuite.fromFolder(p)];
    end
end