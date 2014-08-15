run mrst-core/startup;

% Set up module directories
names = {'autodiff', ...
         'model-io'};
     
names = cellfun(@(x) fullfile(ROOTDIR, '..', ['mrst-', x]), names, ...
                    'UniformOutput', false);
mrstPath('addroot', names{:});

names = {'test-datasets'};
names = cellfun(@(x) fullfile(ROOTDIR, '..', x), names, ...
                    'UniformOutput', false);


mrstPath('addroot', names{:});


mrstModule add ad-unittest
runTestsAD('runIntegration', false, 'runUnit', true, 'writeToDisk', true)
% Disable matlab
exit();