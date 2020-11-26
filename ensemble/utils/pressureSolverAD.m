function pressureSolverAD(problem, varargin)
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