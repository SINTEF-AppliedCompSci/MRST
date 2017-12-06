function results = runTestsAD(varargin)
    opt = struct(...
        'runUnit', true, ...
        'runIntegration', true, ...
        'runExamples', true, ...
        'writeToDisk',  false ...
        );
    opt = merge_options(opt, varargin{:});
    
    import matlab.unittest.TestRunner;

    suites = {};
    names = {};
    if opt.runUnit
        suites{end+1} = getUnitTestSuiteMRST();
        names{end+1} = 'unit';
    end
    
    if opt.runIntegration
        suites{end+1} = getIntegrationTestSuiteMRST();
        names{end+1} = 'integration';
    end
    
    if opt.runExamples
        suites{end+1} = getExampleIntegrationTestSuiteMRST();
        names{end+1} = 'examples';
    end
    
    results = [];
    for i = 1:numel(suites)
        runner = TestRunner.withTextOutput;
        res = runner.run(suites{i});
        results = [results, res]; %#ok
        if opt.writeToDisk
            outd = fullfile(mrstPath('query', 'ad-unittest'), 'output', 'TAP');
            if exist(outd, 'dir') == 0
                mkdir(outd);
            end
            tapFile = fullfile(outd, [names{i}, '.tap']);
            writeTestsTAP_YAMLISH(res, tapFile)
        end
    end
    

end
