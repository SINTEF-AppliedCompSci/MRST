function output = runTestsAD(varargin)
    opt = struct(...
        'runUnit', true, ...
        'runIntegration', true, ...
        'writeToDisk',  false ...
        );
    opt = merge_options(opt, varargin{:});
    
    import matlab.unittest.TestRunner;
    import matlab.unittest.plugins.TAPPlugin;
    import matlab.unittest.plugins.ToFile;
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
    
    output = [];
    for i = 1:numel(suites)
        runner = TestRunner.withTextOutput;
        if opt.writeToDisk
            outd = fullfile(mrstPath('query', 'ad-unittest'), 'output', 'TAP');
            if exist(outd, 'dir') == 0
                mkdir(outd);
            end
            tapFile = fullfile(outd, [names{i}, '.tap']);
            plugin = TAPPlugin.producingOriginalFormat(ToFile(tapFile));
            runner.addPlugin(plugin);
        end
        output = [output; runner.run(suites{i})]; %#ok
    end
end