function results = runTestsAD(varargin)
    opt = struct(...
        'runUnit', true, ...
        'runIntegration', true, ...
        'runExamples', true, ...
        'writeToDisk',  false ...
        );
    opt = merge_options(opt, varargin{:});
    
    import matlab.unittest.TestRunner;

    ix = 1;
    [names, suites] = deal(cell(opt.runUnit + opt.runIntegration, 1));
    if opt.runUnit
        suites{ix} = getUnitTestSuiteMRST();
        names{ix} = 'unit';
        ix = ix+1;
    end
    
    if opt.runIntegration
        suites{ix} = getIntegrationTestSuiteMRST();
        names{ix} = 'integration';
        ix = ix+1;
    end
    
    if opt.runExamples
        [suite, nm] = getExampleIntegrationTestSuiteMRST(mrstPath(), 'seperateModules', true);
        suites = [suites; suite];
        names = [names; nm];
    end
    
    results = [];
    for i = 1:numel(suites)
        fprintf('Running test suite %d of %d: %s\n', i, numel(suites), names{i});
        runner = TestRunner.withTextOutput;
        res = runner.run(suites{i});
        results = [results, res]; %#ok
        if opt.writeToDisk
            outd = fullfile(mrstPath('query', 'ad-unittest'), 'output', 'TAP');
            if exist(outd, 'dir') == 0
                mkdir(outd);
            end
            mrstModule add ad-unittest
            tapFile = fullfile(outd, [names{i}, '.tap']);
            writeTestsTAP_YAMLISH(res, tapFile)
        end
    end
end
