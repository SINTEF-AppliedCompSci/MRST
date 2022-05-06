function [state, schedulereport] = Simulate(model, ...
    twoPhaseOilWaterModel, schedule, verbose)
%
% DESCRIPTION: simulate the fluid model and the schedule
%
% SYNOPSIS:
%   [state, schedulereport] = Simulate(model, ...
%    twoPhaseOilWaterModel, schedule, verbose)
%
% PARAMETERS:
%   model - struct containing following fields:
%   - state: previous or initial state of the model to build on
%   - simulation: time stepping and grid cells information
%   twoPhaseOilWaterModel - model from MRST package
%   schedule: schedule of the simulation
%   verbose: switch for diagnostic messages
%
% RETURNS:
%   state - simulation state for each time step
%   schedulereport - output report
%
% ----------------------------------
% (c) 2020-2022
% Siroos Azizmohammadi
% Omidreza Amrollahinasab
% Montanuniversität Leoben, Austria
% Chair of Reservoir Engineering
% https://dpe.ac.at/
% ----------------------------------
%
%%

state = model.state;
if isfield(model.simulation,'high_percision_mode')
    mode  = model.simulation.high_percision_mode; 
else
    mode = false;
end

% for USS it is better if we use it for better accuracy
stepSel = StateChangeTimeStepSelector(...
'targetProps', {'s'},...
'targetChangeAbs', 0.05);
if mode 
    solver = NonLinearSolver('timeStepSelector', stepSel);
else
    solver = [];
end

if(schedule.counter == 1)
    [~, state, schedulereport] = simulateScheduleAD(...
        state, twoPhaseOilWaterModel, schedule, ...
        'Verbose', verbose, 'NonLinearSolver', solver, ...
        'afterStepFn', []);
else
    [~, state, schedulereport] = simulateScheduleAD(...
        state{end,1}, twoPhaseOilWaterModel, schedule, ...
        'Verbose', verbose, 'NonLinearSolver', solver, ...
        'afterStepFn', []);
end

% different ways to change the solver
% mrstModule add agmg blackoil-sequential
% psolver = AGMGSolverAD();
% tsolver = GMRES_ILUSolverAD();
% % Set up the sequential model
% seqModel = getSequentialModelFromFI(twoPhaseOilWaterModel, 'pressureLinearSolver', psolver,....
%                                            'transportLinearSolver', tsolver);
% solver = NonLinearSolver();
% solver.LinearSolver = CPRSolverAD('ellipticSolver', AGMGSolverAD());
% pressureSolver = AGMGSolverAD('tolerance', 1e-4);
% linsolve = CPRSolverAD('ellipticSolver', pressureSolver, 'relativeTolerance', 1e-3);
