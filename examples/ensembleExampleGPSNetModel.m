%% Ensemble of GPSNET models 
% GPSNET models constitutes of a network model representing the reservoir
% and well connections. The model consists of a set of 1D reservoir models
% between the wells. In the graph terminology, the wells are represented as
% edges and the 1D reservoirs connecting the wells are edges.
% 
% In this example we set up an ensemble of such models, and we run ensemble
% simulations with uncertain rock properties (each connection has
% homogeneous rock properties) and well production indices.

mrstModule add ad-core ad-blackoil mrst-gui ad-props ...
    example-suite incomp ensemble dd-models diagnostics

mrstVerbose off

%% Set the name of the base problem and ensemble size
% The base problem contains all properties that are common throughout the
% ensemble

baseProblemName = 'ensemble_base_problem_simple_gpsnet_model';
baseProblemOptions = {};

ensembleSize = 20;

%% Create example
baseExample = MRSTExample(baseProblemName, ...
                          'deleteOldResults', false, ...
                          'plotNetwork', false);

                      


%% Run the base problem if we wish
% Simulate and plot it for illustration:
simulateExample = false;
plotSimulation = false;
rerunBaseProblemFromScratch = true;

problem = baseExample.getPackedSimulationProblem();

if simulateExample
    if rerunBaseProblemFromScratch
        clearPackedSimulatorOutput(problem, 'prompt', false);
    end
    simulatePackedProblem(problem);

    [DD_wellSols, DD_states, DD_reports] = getPackedSimulatorOutput(problem);

    % Animate water saturation
    if plotSimulation
        figure;
        title('water saturation in DD model');
        plotGrid(problem.SimulatorSetup.model.G, 'FaceAlpha', 0, 'EdgeAlpha', 0.1);
        plotWell(problem.SimulatorSetup.model.G, problem.SimulatorSetup.schedule.control.W);
        view(30, 50);
        pause(1);
        hs = []; % handle for saturation plot, empty initially
        for i = 1:size(problem.SimulatorSetup.schedule.step.val, 1)
            hs = plotCellData(problem.SimulatorSetup.model.G, ...
                              DD_states{i}.s(:,1), DD_states{i}.s(:,1) > 0.1);
            drawnow, pause(0.5);
        end

        plotWellSols(DD_wellSols)
    end
end
                      
%% Define samples that give different rock properties for each connection

rockData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    tmpPoro = zeros(1,4);
    for j = 1:4
        while tmpPoro(j) > 0.4 || tmpPoro(j) < 0.2
            tmpPoro(j) = 0.3 + randn(1)*0.1;
        end
    end
    rockData{i}.poro = tmpPoro;
    rockData{i}.perm = rockData{i}.poro.^3.*(1e-5)^2./(0.81*72*(1-rockData{i}.poro).^2);
end

%% Create sample object
rockSamples = NetworkRockSamples('data', rockData, ...
                                 'connectionIndices', baseExample.options.connectionIndices)

%% Create QoI
qoi = WellQoI(...
    'wellNames', {'P1', 'P2'}, ...
    'fldname', {'qOs', 'qWs'});


%% Create the ensemble
ensemble = MRSTEnsemble(baseExample, rockSamples, qoi, ...
    'simulationStrategy', 'parallel', ...
    'maxWorkers', 8, ...
    'reset', true, ...
    'verbose', true);

%% Run ensemble
ensemble.simulateEnsembleMembers();

%% Plot results
ensemble.plotQoI('subplots', true);

                      