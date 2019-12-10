pth = fullfile(mrstPath('dg'), 'examples');

examples = {'dgExampleDiscretization.m', 'dgExampleBLDisplacement.m', 'dgExampleIFS.m'};
for i = 1:numel(examples)
    run(fullfile(pth, examples{i}));
end