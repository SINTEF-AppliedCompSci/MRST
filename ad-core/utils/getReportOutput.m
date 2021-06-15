function output = getReportOutput(reports, varargin)
    % Get output from report after call to simulateScheduleAD

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

    opt = struct('solver', []                   , ... % Sepcify solver (e.g. TransportSolver)
                 'get'   , []                   , ... % Getter. Can be user-defined function
                 'type'  , 'nonlinearIterations', ... % Used to pick preimplemented getters (se below)
                 'ministeps', true);                  % ministeps
    opt = merge_options(opt, varargin{:});
    if isfield(reports, 'ControlstepReports')
        reports = reports.ControlstepReports;
    end
    nr = numel(reports);
    if isa(reports, 'ResultHandler')
        nr = reports.numelData();
    end
    [opt, scalar] = checkOptions(opt, reports{1});
    % Loop through timesteps and get output
    time=[];
    if(opt.ministeps)
        out = {};
        for t = 1:nr
            ns = numel(reports{t}.StepReports);
            for k = 1:ns
                outs = getStepReportData(reports{t}.StepReports{k}, opt);
                out{end+1,1} = [outs,not(reports{t}.StepReports{k}.Converged)];
                time(end+1) = reports{t}.StepReports{k}.Timestep;
            end
        end         
    else
        out = cell(nr,1);
        for t = 1:nr
            out{t} = getControlStepReportData(reports{t}, opt);
            time(end+1) = reports{t}.StepReports{k}.TimeStep;
        end
    end
    [total, wasted, cuts, steps] = unpack(out, scalar);
    output = struct('total', total, 'wasted', wasted, 'cuts', cuts,'time',cumsum(time));
    if ~isempty(opt.solver)
        output.steps = steps;
    end
end

%-------------------------------------------------------------------------%
% Helpers

%-------------------------------------------------------------------------%
function [opt, scalar] = checkOptions(opt, sample)
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
        v = getControlStepReportData(sample, opt);
    catch
        error('Invalid report outuput requested');
    end
    v = getControlStepReportData(sample, opt);
    scalar = size(v,1) == 1;
end

%-------------------------------------------------------------------------%
function [total, wasted, cuts, steps] = unpack(out, scalar)
    if scalar
        out    = cell2mat(out);
        total  = out(:,1);
        wasted = out(:,2);
        cuts   = out(:,4);
        steps  = out(:,3);
    else
        total  = {cellfun(@(v) v(:,1), out, 'UniformOutput', false)};
        wasted = {cellfun(@(v) v(:,2), out, 'UniformOutput', false)};
        cuts   = {cellfun(@(v) v(:,4), out, 'UniformOutput', false)};
        steps  = {cellfun(@(v) v(:,3), out, 'UniformOutput', false)};
    end
end

%-------------------------------------------------------------------------%
function output = getControlStepReportData(report, opt)
    % Get total output data from a timestep
    % output(1) = total data for all step reports)
    % output(2) = total data for steps that did not result in covergence
    % output(3) = steps for stepreport if ~isempty(opt.solver)
    % output(4) = timestep cuts for control step
    output  = zeros(1,4);
    if ~isfield(report, 'StepReports')
        return
    end
    stepreports = report.StepReports;
    for i = 1:numel(stepreports)
        out = getStepReportData(stepreports{i}, opt);
        if size(output,1) < size(out,1)
            output = zeros(size(out,1), 4);
        end
        output(:,1:3) = output(:,1:3) + out(:,1:3);
        if numel(out) == 4
            output(:,4) = output(:,4) + out(:,4);
        else
            output(:,4) = output(:,4) + ~stepreports{i}.Converged;
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
        out = [out, out.*(~stepreport.Converged), nan(numel(out),1)];
    end
end

%-------------------------------------------------------------------------%
function iterations = getNonLinearIterations(report)
    % Count as one iteration if we actually solved the Newton step, not
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
    % Get total nonlinear iterations if we actually solved the Newton step,
    % not just checked for convergence
    iterations = 0;
    if report.Solved
        if(isfield(report.LinearSolver,'Iterations'))
            iterations = report.LinearSolver.Iterations;
        end
    end
end

%-------------------------------------------------------------------------%
function time = getLinearSolverTime(report)
    % Get total nonlinear iterations if we actually solved the Newton step,
    % not just checked for convergence
    time = 0;
    if report.Solved && isfield(report, 'LinearSolver') && isfield(report.LinearSolver, 'SolverTime')
        time = report.LinearSolver.SolverTime;
        time = time + report.LinearSolver.PostProcessTime;
    end
end
