cd('../../../')
ls
run('mrst-core/startup.m');

% Set up module directories
names = {'autodiff', ...
         'model-io'};
     
names = cellfun(@(x) fullfile(ROOTDIR, '..', ['mrst-', x]), names, ...
                    'UniformOutput', false);
mrstPath('addroot', names{:});

names = {'test-datasets'};
names = cellfun(@(x) fullfile(ROOTDIR, '..', x), names, ...
                    'UniformOutput', false);

mrstPath('register', 'agmg', 'hnil-agmg');
mrstPath('addroot', names{:});


mrstModule add ad-unittest ad-core ad-blackoil
runTestsAD('runIntegration', true, 'runUnit', true, 'writeToDisk', true)
% Disable matlab
exit();