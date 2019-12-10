pth = fullfile(mrstPath('dg'), 'examples');

examples = {'dgExampleDiscretization.m', 'dgExampleBLDisplacement.m', 'dgExampleIFS.m'};

%% Run all

for i = 1:numel(examples)
    run(fullfile(pth, examples{i}));
end

%% Run selected

ix = [2,3];
for i = ix
    run(fullfile(pth, examples{i}));
end