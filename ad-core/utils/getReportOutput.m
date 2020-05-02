function output = getReportOutput(reports, varargin)
    % Get output from report after call to simulateScheduleAD
    opt = struct('solver', []                   , ... % Sepcify solver (e.g. TransportSolver)
                 'get'   , []                   , ... % Getter. Can be user-defined function
                 'type'  , 'nonlinearIterations');    % Used to pick preimplemented getters (se below)
    opt = merge_options(opt, varargin{:});
    if isfield(reports, 'ControlstepReports')
        reports = reports.ControlstepReports;
    end
    nr = numel(reports);
    if isa(reports, 'ResultHandler')
        nr = reports.numelData();
    end
    opt = checkOptions(opt, reports{1});
    % Loop through timesteps and get output
    out = zeros(nr,4);
    for t = 1:nr
        out(t,:) = getControlStepReportData(reports{t}, opt);
    end
    output = struct('total', out(:,1), 'wasted', out(:,2), 'cuts', out(:,4));
    if ~isempty(opt.solver)
        output.steps = out(:,3);
    end
end

%-------------------------------------------------------------------------%
% Helpers

%-------------------------------------------------------------------------%
function opt = checkOptions(opt, sample)
    % Check that requested solver output exists in reports
    if ~isempty(opt.solver) && ~isfield(sample.StepReports{1}.NonlinearReport{1}, opt.solver)
        warning('Did not find %s output in reports, switching to default', opt.solver)
        opt.solver = [];
    end
    if isempty(opt.get)
        % Predefined getters for nonlinear and linear iterations
        switch opt.type
            case 'nonlinearIterations'
                opt.get = @(report) getNonLinearIterations(report);
            case 'nonlinearSolverTime'
                opt.get = @(report) getNonlinearSolverTime(report);
            case 'linearIterations'
                opt.get = @(report) getLinearIterations(report);
            case 'linearSolverTime'
                opt.get = @(report) getLinearSolverTime(report);
            otherwise
                warning('Unknown output type %s, switching to nonlinear iterations', opt.type);
        end
    end
    try % Test by getting controlstep report data from sample
        getControlStepReportData(sample, opt);
    catch
        error('Invalid report outuput requested');
    end
end

%-------------------------------------------------------------------------%
function output = getControlStepReportData(report, opt)
    % Get total output data from a timestep
    % output(1) = total data for all step reports)
    % output(2) = total data for steps that did not result in covergence
    % output(3) = steps for stepreport if ~isempty(opt.solver)
    % output(4) = timestep cuts for control step
    stepreports = report.StepReports;
    output  = zeros(1,4);
    for i = 1:numel(stepreports)
        out = getStepReportData(stepreports{i}, opt);
        output(1:3) = output(1:3) + out(1:3);
        if numel(out) == 4
            output(4) = output(4) + out(4);
        else
            output(4) = output(4) + ~stepreports{i}.Converged;
        end
    end
end

%-------------------------------------------------------------------------%
function out = getStepReportData(stepreport, opt)
    % Get total iterations from a report step
    out = 0;
    if ~isempty(opt.solver)
        % Solver field is nonempty - get report iterations for the solver
        reports    = stepreport.NonlinearReport;
        solver     = opt.solver;
        opt.solver = [];
        for r = 1:numel(reports)
            out = out + getControlStepReportData(reports{r}.(solver), opt);
        end
        out(3) = numel(reports);
    else
        % Solver field is empty - we have reached the bottom!
        reports = stepreport.NonlinearReport;
        for i = 1:numel(reports)
            out = out + opt.get(reports{i});
        end
        out = [out, out.*stepreport.Converged, nan];
    end
end

%-------------------------------------------------------------------------%
function iterations = getNonLinearIterations(report)
    % Count as one iteration if we acutally solved the Newton step, not
    % just checked for convergence
    iterations = report.Solved;
end

%-------------------------------------------------------------------------%
function time = getNonlinearSolverTime(report)
    % Get total nonlinear solver time, which equals assembly time + time
    % spent by the linear solver
    time = 0;
    time = time + report.AssemblyTime;
    time = time + getLinearSolverTime(report);
    
end

%-------------------------------------------------------------------------%
function iterations = getLinearIterations(report)
    % Get total nonlinear iterations if we acutally solved the Newton step,
    % not just checked for convergence
    iterations = 0;
    if report.Solved
        iterations = report.LinearSolver.Iterations;
    end
end

%-------------------------------------------------------------------------%
function time = getLinearSolverTime(report)
    % Get total nonlinear iterations if we acutally solved the Newton step,
    % not just checked for convergence
    time = 0;
    if report.Solved && isfield(report, 'LinearSolver') && isfield(report.LinearSolver, 'SolverTime')
        time = report.LinearSolver.SolverTime;
        time = time + report.LinearSolver.PostProcessTime;
    end
end