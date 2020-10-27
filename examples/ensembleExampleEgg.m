mrstModule add ad-core ad-props ad-blackoil example-suite ensemble ...
    mrst-gui
mrstVerbose on

%%
example = MRSTExample('egg_wo');

%%
generatorFn = @(problem, seed) getDeckEGG('realization', seed-1);
samples = DeckSamples('generatorFn', generatorFn, 'num', 101);

%%
problem = example.getPackedSimulationProblem();
data    = samples.getSample(91, problem);
problem = samples.setSample(data, problem);

simulatePackedProblem(problem);

%%
is_prod = vertcat(example.schedule.control(1).W.sign) < 0;
qoi = WellQoI('wellIndices', is_prod, 'fldname', 'qOs', 'combined', true);

%%
qoiVal = qoi.validateQoI(problem);
qOs    = qoiVal.computeQoI(problem);

%%
ensemble = MRSTEnsemble(example, samples, qoi, 'simulationType', 'background');

%%
close all
ensemble.simulateEnsembleMembers(20, 'plotProgress', true);

%%
ensemble.qoi.plotEnsembleQoI(ensemble);