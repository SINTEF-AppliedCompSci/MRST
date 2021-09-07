%% Calibration of a GPSNet model for Egg using adjoints
% This example demonstrates how to set up a GPSNet type reduced network for
% the Egg model, a conceptual model of a channelized reservoir. For details
% about the model, see Jansen, J. D., et al. "The egg model -- a geological
% ensemble for reservoir simulation." Geoscience Data Journal 1.2 (2014):
% 192-195. The network topology can be determined in different ways. The
% default setup of the example is to use flow diagnostic postprocessing of
% the reference simulation to determine nonzero interwell connections.
%
% For the calibration, we use the Broyden–Fletcher–Goldfarb–Shanno (BFGS)
% algorithm. This is an iterative line-search method that gradually
% improves an approximation to the Hessian matrix of the mismatch function,
% obtained only from adjoint gradients via a generalized secant method.
%
% This example was first introduced in MRST 2021b.
mrstModule add ad-core ad-blackoil deckformat diagnostics mrst-gui ...
    ad-props incomp optimization network-models example-suite

%% Setup 3D reference model
% The Egg data set consists of a base case and an ensemble of one hundred
% different permeability realizations. (If you do not have the data
% available, you have to download it manually from the internet; see
% mrstDatasetGUI for details.) Herein, we will use realization number 0,
% the base case, but you can also run the example with any number between 1
% and 100.
realization = 0; 
[G, ~, ~, deck] = setupEGG('realization', realization);
[state0, modelRef, scheduleRef, nonlinear] = ...
    initEclipseProblemAD(deck, 'G', G, 'TimestepStrategy', 'none');
modelRef.getPhaseNames();

problem = packSimulationProblem(state0, modelRef, scheduleRef, ...
    ['EGG_realization_',num2str(realization)], 'NonLinearSolver', nonlinear);
simulatePackedProblem(problem);

[wellSolsRef, statesRef] = getPackedSimulatorOutput(problem);
scheduleRef = problem.SimulatorSetup.schedule;
Wref        = scheduleRef.control.W;

%% Create the network
% We start by creating a network that connects injectors and producers.
% This network describes the possible 1D flow paths that each will have an
% associated pore volume and transmissibility.

% To create the network, we only need a single perforation from each well,
% which we somewhat arbitrarily pick from the top.
Wnw = Wref;
for i = 1:numel(Wnw)
    Wnw(i).cells = Wnw(i).cells(7);
end

% Select type and create network 
% The module supplies different ways to do this: all to all, injectors to
% producers (and vice versa), connections specified manually by the user,
% connections determined by a flow diagnostics analysis of the 3D
% geological model, or connections based on flow diagnostics analysis of
% the fine-scale reference simulation. Here, we use the second last
% approach but supply code for all other options as well.

networkType = 'fd_preprocessor';
switch networkType
    case 'all_to_all'
        ntwrk =  Network(Wnw, modelRef.G, 'type', networkType);
    case 'injectors_to_producers'
        ntwrk =  Network(Wnw, modelRef.G, 'type', networkType, ...
                         'injectors', 1:8, 'producers', 9:12);
    case 'user_defined_edges'
        edges = [1 9;2 9;2 10;3 9;3 11;
                 4 9; 4 10; 4 11; 4 12;
                 5 10; 5 12; 6 11; 7 11; 7 12; 8 12];
        ntwrk =  Network(Wnw, modelRef.G, 'type', networkType, ...
                         'edges',edges);
    case 'fd_preprocessor'
         ntwrk = Network(Wnw, modelRef.G, 'type', networkType, ...
                         'problem', problem,                   ...
                         'flow_filter',1*stb/day);
    case 'fd_postprocessor'
        ntwrk = Network(Wnw, modelRef.G, 'type', networkType, ...
                         'problem', problem,                  ...
                         'state_number',40,                   ...
                         'flow_filter', 1*stb/day);
    otherwise
        error('\nNetwork of type %s is not implemented\n', networkType);
end
                     
% Plot the network
% If based on flow diagnostics, we plot the network twice to show the
% relative magnitude of the associated transmissibilities and pore volumes.
figure; ntwrk.plotNetwork()

%% Create the data-driven model
% We subgrid each flow path into ten uniform cells and map the resulting
% network onto a rectangular Cartesian grid having the same number of rows
% as the number of flow paths.  The fluid model is copied from the
% reference model.

% Grid and petrophysics
% The grid is set to have an aspect ratio of [5 1 1] and a volum that
% matches the bulk volume of the reservoir. The constant petrophysical
% properties are dummy values primarily used to compute an initial guess
% for matching parameters that have not be set by the network type.
cellsPerPath = 10;
L    = nthroot(sum(modelRef.operators.pv./modelRef.rock.poro)*25,3);
G    = cartGrid([cellsPerPath, 1, numedges(ntwrk.network)], [L, L/5 ,L/5]*meter^3);
G    = computeGeometry(G);
rock = makeRock(G, 200*milli*darcy, 0.1);

% Fluid model
gravity off
fluid = modelRef.fluid;
model = GenericBlackOilModel(G, rock, fluid, 'gas', false);
model.OutputStateFunctions = {};

% Then we map the Network into the MRST model
ntwrkModel = NetworkModel(model, cellsPerPath, ntwrk.network, Wref);
model      = ntwrkModel.model;
model      = model.validateModel();
W          = ntwrkModel.W;
state0     = initState(G, W , 400*barsa,[0.2, 0.8]);

%% Specify the training problem and parameters to be matched
% We need to specify which of the parameters in the model we wish to match
% and any scaling and box limits that will applied to these in the BFGS
% optimization. These are specified through through the ModelParameter
% class from the optimization module.

% Training problem
% We use the first 40 time steps from the fine-scale simulation to train
% the network model. The remaining time steps will be used for prediction.
trainSchedule = simpleSchedule(scheduleRef.step.val(1:40), 'W', W);
trainProb = struct('model', model, 'schedule', trainSchedule, 'state0', state0);

% Parameter lumping
% In the GPSNet type of models, we only match a single pore volume and a
% single transmissibility for each network edge. We thus need a mapping
% between edges and cells/faces in the Cartesian grid to describe the
% corresponding lumping of parameters in this grid.
[cellEdgeNo, cellIx] = reorganizeIndices(ntwrkModel.Graph.Edges.Cell_Indices);
[faceEdgeNo, faceIx] = reorganizeIndices(ntwrkModel.Graph.Edges.Face_Indices);

% Set parameters                             
parameters =  {};                          
parameters{1} = ModelParameter(trainProb, 'name', 'conntrans', 'scaling', ...
                               'linear', 'relativeLimits', [.001 20]);
parameters{2} = ModelParameter(trainProb, 'name', 'porevolume', 'lumping', ...
                               cellEdgeNo,'subset', cellIx, ...
                               'relativeLimits', [.01 5]);
parameters{3} = ModelParameter(trainProb, 'name', 'transmissibility', ...
                               'lumping', faceEdgeNo, 'subset', ...
                               faceIx, 'scaling', 'log', ...
                               'relativeLimits', [.001 100]);

%% Define the mismatch function
% The mismatch function is defined as a function handle to a library
% function from the optimization module that computes the mismatch between
% a given simulation and a reference state. For an oil-water system, the
% match is computed based on three quantities (water/oil rate and bhp) and
% these must be given an associated weight.
weighting =  {'WaterRateWeight',  day/150, ...
              'OilRateWeight',    day/day, ...
              'BHPWeight',        1/(40*barsa)};   
mismatchFn = @(model, states, schedule, states_ref, tt, tstep, state) ...
    matchObservedOW(model, states, schedule, states_ref,...
                    'computePartials', tt, 'tstep', tstep, weighting{:},...
                    'state',state,'from_states',false);

%% Set parameters defining the initial model
% To define the initial model, we extract the network parameters from the
% initial specification, scale them to the unit interval, and organize them
% into a vector that can be passed onto a library function from the
% optimization module that performs the forward simulation and evaluates
% the mismatch using mismatchFn.

% Extract and scale parameters
% Remember to overwrite pore volumes and transmissibilities if initial
% values for these have been computed by flow diagnostics
values    = applyFunction(@(p)p.getParameterValue(trainProb), parameters);
values{1} =  7*values{1};
if any(strcmp(networkType,{'fd_preprocessor','fd_postprocessor'}))
    values{2} = ntwrk.network.Edges.PoreVolume/cellsPerPath;
    values{3} = ntwrk.network.Edges.Transmissibility;
end
for k = numel(values):-1:1    
    u{k} = parameters{k}.scale(values{k});
end
pvec0 = vertcat(u{:});  
 
%{
% Perform a forward simulation and evaluate the mismatch
% This forward simulation is redundant because the same simulation is
% performed inside the model calibration routine. It is only included here
% to illustrate the mismatch.
[misfitVal0,~,wellSols0,states0] = ...
    evaluateMatch(pvec0,mismatchFn,trainProb,parameters,statesRef,'Gradient','none');
 
% Plot well curves for the reference and the initial model
plotWellSols({wellSolsRef,wellSols0}, ...
    {scheduleRef.step.val,trainProb.schedule.step.val}, ...
    'datasetnames', {'reference','initial'}, 'zoom', true, ...
    'field', 'qOs', 'SelectedWells', 9)
%}

%% Model calibration
% Calibrate the model using the BFGS method with parameters
objh = @(p) evaluateMatch(p, mismatchFn, trainProb, parameters, statesRef);

[v, pOpt, history] = unitBoxBFGS(pvec0, objh,'objChangeTol', 1e-8, ...
    'maxIt', 25, 'lbfgsStrategy', 'dynamic', 'lbfgsNum', 5);

%% Evaluate mismatch over the full simulation schedule
prob = trainProb;
prob.schedule = simpleSchedule(scheduleRef.step.val, 'W', W);
[~,~,wellSolsOpt] = ...
    evaluateMatch(pOpt,mismatchFn,prob,parameters,statesRef,'Gradient','none');
[~,~,wellSols0] = ...
    evaluateMatch(pvec0,mismatchFn,prob,parameters,statesRef,'Gradient','none');

%% Plot well curves
plotWellSols({wellSolsRef,wellSols0,wellSolsOpt}, ...
    {scheduleRef.step.val,prob.schedule.step.val,prob.schedule.step.val}, ...
    'datasetnames',{'reference','initial','calibrated'},'zoom', true, ...
    'field', 'qOs', 'SelectedWells', 9)
legend('reference model','initial DD model','calibrated model')

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
