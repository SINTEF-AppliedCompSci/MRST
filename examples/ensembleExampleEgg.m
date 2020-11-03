%% Egg ensamble example
% In this example, we show how to set up an ensemble from a collection of
% deck files, exemplified by the Egg model [1]

%% Add modules
mrstModule add ad-core ad-props ad-blackoil example-suite ensemble ...
    mrst-gui
mrstVerbose on

%% Set up base problem
example = MRSTExample('egg_wo');
example.plot(example.model.rock, 'log10', true); colormap(pink);

%% Set up samples
% We use the function getDeckEGG as a generatorFn to set up a sample. The
% realizations are numbered from 0 to 100, hence seed-1.
generatorFn = @(problem, seed) getDeckEGG('realization', seed-1);
% The class DeckSamples implements what we need for setting up samples from
% a deck file. The class uses `initEclipseProblemAD` to set up the problem,
% and extra input arguments to this function can be provided in the class
% property initArgs. The property gridFromDeck is used to determine if the
% simulation grid should be constructed from the deck, or taken from the
% base problem.
samples = DeckSamples('generatorFn', generatorFn, 'num', 101);
disp(samples)

%% Set up QoI
% For our QoI, we choose the total oil production rate
is_prod = vertcat(example.schedule.control(1).W.sign) < 0;
qoi = WellQoI('wellIndices', is_prod, 'fldname', 'qOs', 'combined', true);

%% Set up ensemble
ensemble = MRSTEnsemble(example, samples, qoi, ...
               'simulationStrategy', 'background'); % Run in the background

%% Simulate the ensemble members
% We simulate 20 samples. Each time the code in this block is called, 20
% new samples will be simulated. To reset the ensemble data, call
% ensemble.reset();
ensemble.simulateEnsembleMembers('range', 20, 'plotProgress', true);

%% Plot the QoI
ensemble.plotQoI();