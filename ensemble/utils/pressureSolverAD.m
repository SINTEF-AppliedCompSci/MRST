function pressureSolverAD(problem, varargin)
%Undocumented Utility Function

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

    opt = struct('dt'          , nan  , ...
                 'ctrlID'      , 1    , ...
                 'pressureOnly', true , ...
                 'incomp'      , false);
    
    opt = merge_options(opt, varargin{:});
    % Get model
    model = problem.SimulatorSetup.model;
    model = getModel(model, opt);
    problem.SimulatorSetup.model = model;
    % Get single-step schedule
    schedule = problem.SimulatorSetup.schedule;
    schedule = getSchedule(schedule, opt);
    problem.SimulatorSetup.schedule = schedule;
    % Get solver
    solver = problem.SimulatorSetup.NonLinearSolver;
    solver = getSolver(solver, opt);
    problem.SimulatorSetup.NonLinearSolver = solver;
    % Simulate
    simulatePackedProblem(problem);
    % Add boolean indicating that this was a single-step pressure solve
    state = problem.OutputHandlers.states{1};
    state.singlePressureSolve = true;
    problem.OutputHandlers.states{1} = state;
end

%-------------------------------------------------------------------------%
function model = getModel(model, opt)
    if opt.pressureOnly
        require sequential
        model = PressureModel(model);
    end
end

%-------------------------------------------------------------------------%
function schedule = getSchedule(schedule, opt)
    time = sum(schedule.step.val);
    if isnan(opt.dt)
        opt.dt = time*0.1;
    end
    schedule.step.val = opt.dt;
    schedule.step.control = schedule.step.control(opt.ctrlID);
end

%-------------------------------------------------------------------------%
function solver = getSolver(solver, opt)
    if opt.pressureOnly
        lsolver = AMGCLSolverAD('tolerance'    , 1e-4         , ...
                                'maxIterations', 100          , ...
                                'solver'       , 'bicgstab'   , ...
                                'direct_coarse', true         , ...
                                'coarsening'   , 'aggregation');
        solver.LinearSolver = lsolver;
    end
end
