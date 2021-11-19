mrstModule add example-suite
mrstModule add ad-core ad-props ad-blackoil

mrstVerbose on

%%
names = {'buckley_leverett_wo', 'qfs_wo'};
reports = struct();
for name = names
    test   = RegressionTest(name{1});
    reports.(name{1}) = test.runRegressionTest();
end

%%