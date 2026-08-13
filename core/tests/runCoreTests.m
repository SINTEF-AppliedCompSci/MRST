function results = runCoreTests(varargin)
%Run unit and integration tests for the core module.

    opt = struct('runIntegration', true);
    opt = merge_options(opt, varargin{:});

    suite = getCoreUnitTestSuiteMRST();
    if opt.runIntegration
        suite = [suite, getCoreIntegrationTestSuiteMRST()];
    end

    runner = matlab.unittest.TestRunner.withTextOutput;
    results = runner.run(suite);
end
