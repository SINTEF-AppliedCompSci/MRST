function suite = getExampleIntegrationTestSuiteMRST(modules)
    mrstModule add ad-unittest
    import matlab.unittest.TestSuite;
    suite = TestSuite.fromClass(?MRSTExampleTests);
    if nargin > 0
        if ~iscell(modules)
            modules = {modules};
        end
        active = arrayfun(@(x) any(strcmpi(x.Parameterization(2).Value, modules)),  suite);
        suite = suite(active);
    end
end