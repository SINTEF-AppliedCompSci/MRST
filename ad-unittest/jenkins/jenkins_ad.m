cd('../../../')
ls
run('mrst-core/startup.m');

% Set up module directories
names = {'autodiff', ...
         'internal', ...
         'multiscale', ...
         'visualization', ...
         'model-io', ...
         'solvers'};
     
names = cellfun(@(x) fullfile(ROOTDIR, '..', ['mrst-', x]), names, ...
                    'UniformOutput', false);
mrstPath('addroot', names{:});

names = {'test-datasets'};
names = cellfun(@(x) fullfile(ROOTDIR, '..', x), names, ...
                    'UniformOutput', false);
mrstPath('register', 'co2lab', fullfile(ROOTDIR, 'mrst-co2lab', 'co2lab'));
mrstPath('register', 'agmg', 'hnil-agmg');
mrstPath('addroot', names{:});

% Remove previous tests
tap_folder = fullfile(mrstPath('query', 'ad-unittest'), 'output', 'TAP');
delete(tap_folder);

mrstModule add ad-unittest ad-core ad-blackoil
runTestsAD('runIntegration', true, 'runUnit', true, 'runExamples', true, 'writeToDisk', true)
% Disable matlab
exit();