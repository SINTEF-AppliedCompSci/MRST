function [wellSols, states, schedulereport] = simulateScheduleAD(initState, model, schedule, varargin)
% Run a schedule for a non-linear physical model using an automatic differention
%
% SYNOPSIS:
%   wellSols = simulateScheduleAD(initState, model, schedule)
%
%   [wellSols, state, report]  = simulateScheduleAD(initState, model)
%
% DESCRIPTION:
%   This function takes in a valid schedule file (see required parameters)
%   and runs a simulation through all timesteps for a given model and
%   initial state.
%
%   simulateScheduleAD is the outer bookkeeping routine used for running
%   simulations with non-trivial control changes. It relies on the model
%   and (non)linear solver classes to do the heavy lifting.
%
% REQUIRED PARAMETERS:
%
%   initState    - Initial reservoir/model state. It should have whatever
%                  fields are associated with the physical model, with
%                  reasonable values. It is the responsibility of the user
%                  to ensure that the state is properly initialized.
%
%   model        - The physical model that determines jacobians/convergence
%                  for the problem. This must be a subclass of the
%                  PhysicalModel base class.
%
%   schedule     - Schedule containing fields step and control, defined as
%                  follows:
%                         - schedule.control is a struct array containing
%                         fields that the model knows how to process.
%                         Typically, this will be the fields such as .W for
%                         wells or .bc for boundary conditions.
%
%                         - schedule.step contains two arrays of equal size
%                         named "val" and "control". Control is a index
%                         into the schedule.control array, indicating which
%                         control is to be used for the timestep.
%                         schedule.step.val is the timestep used for that
%                         control step.
% 
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   
%  'Verbose'        - Indicate if extra output is to be printed such as
%                     detailed convergence reports and so on. 
%                    
%
%  'OutputMinisteps' - The solver may not use timesteps equal to the
%                      control steps depending on problem stiffness and 
%                      timestep selection. Enabling this option will make
%                      the solver output the states and reports for all
%                      steps actually taken and not just at the control
%                      steps. See also 'convertReportToSchedule' which can
%                      be used to construct a new schedule from these
%                      timesteps.
%
% 'NonLinearSolver' - An instance of the NonLinearSolver class. Consider
%                     using this if you for example want a special timestep
%                     selection algorithm. See the NonLinear Solver class
%                     docs for more information.
%
% 'OutputHandler'   - Output handler class, for example for writing states
%                     to disk during the simulation or in-situ
%                     visualization. See the ResultHandler base class.
%
% 'LinearSolver'    - Class instance subclassed from LinearSolverAD. Used
%                     to solve linearized problems in the NonLinearSolver
%                     class. Note that if you are passing a
%                     NonLinearSolver, you might as well put it there.
%
% RETURNS:
%  wellSols         - Well solution at each control step (or timestep if
%                     'OutputMinisteps' is enabled.
%
%  states           - State at each control step (or timestep if
%                     'OutputMinisteps' is enabled.
%
%  schedulereport   - Report for the simulation. Contains detailed info for
%                     the whole schedule run, as well as arrays containing
%                     reports for the whole stack of routines called during
%                     the simulation.
%
% NOTE:
%   For valid models, see ThreePhaseBlackOilModel or TwoPhaseOilWaterModel.
%
% SEE ALSO:
%   computeGradientAdjointAD, PhysicalModel

%{
Copyright 2009-2014 SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    assert(isa(model, 'PhysicalModel'), ...
        'The model must have PhyiscalModel as its base class!')
    
    validateSchedule(schedule);
    
    opt = struct('Verbose',         mrstVerbose(),...
                 'OutputMinisteps', false, ...
                 'NonLinearSolver', [], ...
                 'OutputHandler',   [], ...
                 'LinearSolver',    []);

    opt = merge_options(opt, varargin{:});

    vb = opt.Verbose;
    %--------------------------------------------------------------------------

    dt = schedule.step.val;
    tm = cumsum(dt);
    dispif(vb, '*****************************************************************\n')
    dispif(vb, '********** Starting simulation: %5.0f steps, %5.0f days *********\n', numel(dt), tm(end)/day)
    dispif(vb, '*****************************************************************\n')
    
    solver = opt.NonLinearSolver;
    if isempty(solver)
        solver = NonLinearSolver('linearSolver', opt.LinearSolver);
    elseif ~isempty(opt.LinearSolver)
        % We got a nonlinear solver, but we still want to override the
        % actual linear solver passed to the higher level schedule function
        % we're currently in
        solver.LinearSolver = opt.LinearSolver;
        assert(isa(solver.LinearSolver, 'LinearSolverAD'), ...
        'Passed linear solver is not an instance of LinearSolverAD class!')
    end
    nSteps = numel(dt);

    [wellSols, states, reports] = deal(cell(nSteps, 1));
    wantStates = nargout > 1;
    wantReport = nargout > 2;

    getWell = @(index) schedule.control(schedule.step.control(index)).W;
    state = initState;
    if ~isfield(state, 'wellSol')
        state.wellSol = initWellSolAD(getWell(1), model, state);
    end
    
    failure = false;
    simtime = zeros(nSteps, 1);
    prevControl = nan;
    for i = 1:nSteps
        fprintf('Solving timestep %d of %d at %s\n', i, nSteps, formatTimeRange(tm(i)));
        currControl = schedule.step.control(i);
        if prevControl ~= currControl 
            W = schedule.control(currControl).W;
            forces = model.getDrivingForces(schedule.control(currControl));
            prevControl = currControl;
        end

        timer = tic();
        
        
        state0 = state;
        state0.wellSol = initWellSolAD(W, model, state);
        
        if opt.OutputMinisteps
            [state, report, ministeps] = solver.solveTimestep(state0, dt(i), model, ...
                                            forces{:}, 'controlId', currControl);
        else
            [state, report] = solver.solveTimestep(state0, dt(i), model,...
                                            forces{:}, 'controlId', currControl);
        end
        t = toc(timer);
        
        if ~report.Converged
            warning('Nonlinear solver aborted, returning incomplete results!');
            failure = true;
            break;
        end
        dispif(vb, 'Completed %d iterations in %2.2f seconds (%2.2fs per iteration)\n', ...
                    report.Iterations, t, t/report.Iterations);
        

        W = updateSwitchedControls(state.wellSol, W);
        
        
        % Handle massaging of output to correct expectation
        if opt.OutputMinisteps
            % We have potentially several ministeps desired as output
            nmini = numel(ministeps);
            ise = find(cellfun(@isempty, states), 1, 'first');
            if isempty(ise)
                ise = numel(states) + 1;
            end
            ind = ise:(ise + nmini - 1);
            states_step = ministeps;
        else
            % We just want the control step
            ind = i;
            states_step = {state};
        end
        wellSols_step = cellfun(@(x) x.wellSol, states_step, 'UniformOutput', false);
        
        wellSols(ind) = wellSols_step;
        
        if ~isempty(opt.OutputHandler)
            opt.OutputHandler{ind} = states_step;
        end
        
        if wantStates
            states(ind) = states_step;
        end
        
        if wantReport
            reports{i} = report;
        end
    end
    
    if wantReport
        reports = reports(~cellfun(@isempty, reports));
        
        schedulereport = struct();
        schedulereport.ControlstepReports = reports;
        schedulereport.ReservoirTime = cumsum(schedule.step.val);
        schedulereport.Converged  = cellfun(@(x) x.Converged, reports);
        schedulereport.Iterations = cellfun(@(x) x.Iterations, reports);
        schedulereport.SimulationTime = simtime;
        schedulereport.Failure = failure;
    end
end

function validateSchedule(schedule)
    assert(isfield(schedule, 'control'));
    assert(isfield(schedule, 'step'));
    
    steps = schedule.step;
    
    assert(isfield(steps, 'val'));
    assert(isfield(steps, 'control'));
    
    assert(numel(steps.val) == numel(steps.control));
    assert(numel(schedule.control) >= max(schedule.step.control))
    assert(min(schedule.step.control) > 0);
end
