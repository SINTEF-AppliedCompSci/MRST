function flowDiagnosticsSolver(problem, varargin)
    opt = struct('dt'             , []   , ...
                 'LinearSolver'   , []   , ...
                 'diagnosticsType', 'tof');
    [opt, extra] = merge_options(opt, varargin{:});
    if isempty(opt.dt)
        % Set timestep equal to 5% of total simulation time if not given
        opt.dt = sum(problem.SimulatorSetup.schedule.step.val)*0.05;
    end
    % Get pressure model
    model  = PressureModel(problem.SimulatorSetup.model);
    pmodel = model.parentModel.validateModel();
    ndof   = pmodel.G.cells.num*pmodel.getNumberOfComponents();
    % Set solver
    if isempty(opt.LinearSolver) && ndof > 1e4
        opt.LinearSolver = AMGCLSolverAD('tolerance', 1e-4);
    end
    % Solve a single pressure step
    W = problem.SimulatorSetup.schedule.control(1).W;
    [state, report] = standaloneSolveAD(problem.SimulatorSetup.state0, model, opt.dt, ...
                                               'W', W, 'LinearSolver', opt.LinearSolver);
    % Adjust pore volumes to account for multiphase flow (does
    % nothing at the moment)
    model = getReservoirModel(model);
    model = model.validateModel();
    model = adjustPoreVolumes(model);
    % Get pore volume
    pv = model.getProp(state, 'PoreVolume');
    % Compute time of flight
    maxTOF = sum(problem.SimulatorSetup.schedule.step.val)*10;
    timer  = tic();
    D      = computeTOFandTracer(state, model.G, model.rock, ...
             'wells', W, 'maxTOF', maxTOF, extra{:});
    % Use TOF or residence time distribution (more accurate)
    switch opt.diagnosticsType
        case 'tof'
            % Compute flow capacity and storage capacity from tof
            [F, Phi] = computeFandPhi(pv, D.tof);
        case 'rtd'
            WP  = computeWellPairs(state, model.G, model.rock, W, D);
            rtd = computeRTD(state, model.G, pv, D, WP, W, ...
                                             'showWaitbar', false);
            [F, Phi] = computeFandPhiFromDist(rtd, 'sum', true);
    end
    % Compute sweep efficiency
    [Ev, tD] = computeSweep(F, Phi);
    state.flowDiagnostics = struct('D', D, 'Ev', Ev, 'tD', tD, 'pv', pv);
    time = toc(timer);
    % Add FD estimation time to total time
    report.SimulationTime = report.SimulationTime + time;
    % Write results to disk
    problem.OutputHandlers.states{1}  = state;
    problem.OutputHandlers.states{1}  = state;
    problem.OutputHandlers.reports{1} = report;
end

%-------------------------------------------------------------------------%
function model = adjustPoreVolumes(model)
    % TODO: implement pore volume adjustment based on BL speed
end