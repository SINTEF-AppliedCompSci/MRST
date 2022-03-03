function [state, schedulereport] = Simulate(model, twoPhaseOilWaterModel, schedule, verbose)
% <keywords>
%
% Purpose : simulated the schedule with the input modeling parameters
%
% Syntax : 
%   [state, schedulereport] = Simulate(model, twoPhaseOilWaterModel, schedule, verbose)
%
% Input Parameters :
%   model: struct containing modeling parameters
%   twoPhaseOilWaterModel: struct containing modeling params from MRST
%   functions
%   schedule: simulation schedule
%   verbose: boolean to show the diagnostics
%
% Return Parameters :
%   state: pressure and saturation solutions from the simulation
%   schedulereport: detailed reports from the simulation
%
% Description :
%
% Author : 
%    Siroos Azizmohammadi
%    Omidreza Amrollahinasab
%
% History :
% \change{1.0}{09-Nov-2021}{Original}
%
% --------------------------------------------------
% (c) 2021, Siroos Azizmohammadi,
% Omidreza Amrollahinasab
% Chair of Reservoir Engineering, University of Leoben, Austria
% email: info@dpe.ac.at
% url: dpe.ac.at
% --------------------------------------------------
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
            'Verbose', verbose, 'NonLinearSolver', solver);
    else
        [~, state, schedulereport] = simulateScheduleAD(...
            state{end,1}, twoPhaseOilWaterModel, schedule, ...
            'Verbose', verbose, 'NonLinearSolver', solver);
    end
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
