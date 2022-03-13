function output = getReportOutput(reports, varargin)
%Get output from report after call to simulateScheduleAD
% SYNOPSIS:
%   output = getReportOutput(reports, 'pn1', pv1', ...)
%
% REQUIRED PARAMETERS:
%   reports - Simulation reports as returned (or retrieved from disk) from
%             `simulateScheduleAD`.
%
% OPTIONAL PARAMETERS:
%   solver - Get report output for a solver (e.g., 'TransportSolver').
%            Assumes results are stored as
%            reports{i}.StepReports{j}.NonlinearReport{k}.(solver)
%   get    - Function that takes in a nonlinear report and returns the
%            quantity of interest.
%   type   - String to set a predefined get function. Options are
%               * nonlinearIterations (default)
%               * linearIterations
%               * nonlinearSolverTime
%               * linearSolverTime
%
% RETURNS:
%   output - A structure with the following fields
%               * total : total per timestep or ministep
%               * wasted: total that was part of a timestep/ministep that
%                         did not converge
%               * cuts  : number of timestep cuts (not applicable for
%                         ministeps)
%               * time  : reservoir time
%
% EXAMPLES:
%   % Assuming simulation `reports` is a cell array of simulation reports
%   for each timstep of a simulation, run with `simulateScheduleAD`:
%
%   % Get nonlinear iterations per timestep
%   nlit = getReportOutput(reports);
%
%   % Get linear solver time for each ministep
%   ltime = getReportOutput(reports, 'type', 'linearSolverTime', 'ministeps', true);
%
%   % Get linear solver time for each transport timestep (assuming this was
%   a sequential simulation)
%   lit = getReportOutput(reports, 'solver', 'TransportSolver', 'type', 'linearIterations');

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

    opt = struct('solver'   , []                   , ... % Sepcify solver (e.g. TransportSolver)
                 'get'      , []                   , ... % Getter. Can be user-defined function
                 'type'     , 'nonlinearIterations', ... % Used to pick preimplemented getters (se below)
                 'ministeps', false                );    % Ouput data per ministep instead of timestep
    opt = merge_options(opt, varargin{:});
    if isfield(reports, 'ControlstepReports')
        reports = reports.ControlstepReports;
    end
    % Check input
    assert(~(opt.ministeps && ~isempty(opt.solver)), ['Output per ', ...
        'ministep for a specific nonlinear solver is not implemented yet']);
    if isa(reports, 'ResultHandler')
        assert(~opt.ministeps, ['Output per ministep for '  , ...
            'result-handler report data not implemented -- ', ...
            'please load reports first']                    );
        nr = reports.numelData();
    else
        nr = numel(reports);
    end
    [opt, scalar] = checkOptions(opt, reports{1});
    % Loop through timesteps and get output
    [out, time] = deal(cell(nr,1));
    for t = 1:nr
        [out{t}, time{t}] = getControlStepReportData(reports{t}, opt);
    end
    % Unpack data
    if opt.ministeps, out = vertcat(out{:}); end
    time = vertcat(time{:});
    [total, wasted, cuts, steps] = unpack(out, scalar);
    % Make output struct
    output = struct('total' , total , ...
                    'wasted', wasted, ...
                    'cuts'  , cuts  , ...
                    'time'  , time  );
    if ~isempty(opt.solver), output.steps = steps; end
end

%-------------------------------------------------------------------------%
% Helpers

%-------------------------------------------------------------------------%
function [opt, scalar] = checkOptions(opt, sample)
% Check that requested solver output exists in reports
    assert(isfield(sample, 'StepReports'), ['Each report must be a ', ...
        'struct with at least a `StepReport` field']);
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
    if opt.ministeps, v = v{1}; end
    scalar = size(v,1) == 1;
end

%-------------------------------------------------------------------------%
function [total, wasted, cuts, steps] = unpack(out, scalar)
% Unpack data in out
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
function [output, time] = getControlStepReportData(report, opt)
% Get total output data from a timestep
% output(1) = total data for all step reports)
% output(2) = total data for steps that did not result in covergence
% output(3) = steps for stepreport if ~isempty(opt.solver)
% output(4) = timestep cuts for control step
    if ~isfield(report, 'StepReports')
        output = zeros(1,4);
        time   = 0;
        return
    end
    stepreports = report.StepReports;
    ns          = numel(stepreports);
    if opt.ministeps
        output = {};
        time   = [];
        add    = @(oA, oB) [oA; oB];
    else
        output = zeros(1,4);
        time   = 0;
        add    = @(oA, oB) oA + oB;
    end
    
    for i = 1:ns
        [out, t] = getStepReportData(stepreports{i}, opt);
        if ~opt.ministeps && size(output,1) < size(out,1)
            output = zeros(size(out,1), 4);
        end
        if numel(out) == 3
            out(4) = ~stepreports{i}.Converged;
        end
        output = add(output, out);
        time   = add(time, t);
    end
end

%-------------------------------------------------------------------------%
function [out, time] = getStepReportData(stepreport, opt)
% Get total iterations from a report step
    [out, time] = deal(0);
    if ~isempty(opt.solver)
        % Solver field is nonempty - get report iterations for the solver
        reports    = stepreport.NonlinearReport;
        solver     = opt.solver;
        opt.solver = [];
        for r = 1:numel(reports)
            [o, t] = getControlStepReportData(reports{r}.(solver), opt);
            out  = out + o;
            time = time + t;
        end
        out(3) = numel(reports);
    else
        % Solver field is empty - we have reached the bottom!
        reports = stepreport.NonlinearReport;
        for i = 1:numel(reports)
            out = out + opt.get(reports{i});
        end
        out  = [out, out.*(~stepreport.Converged), nan(numel(out),1)];
        time = stepreport.LocalTime;
    end
end

%-------------------------------------------------------------------------%
function iterations = getNonLinearIterations(report)
% Getter for nonlinear iterations
    % Count as one iteration if we actually solved the Newton step, not
    % just checked for convergence
    iterations = report.Solved;
end

%-------------------------------------------------------------------------%
function time = getNonlinearSolverTime(report)
% Getter for nonlinear solver time
    % Get total nonlinear solver time, which equals assembly time + time
    % spent by the linear solver
    time = 0;
    time = time + report.AssemblyTime;
    time = time + getLinearSolverTime(report);
    
end

%-------------------------------------------------------------------------%
function iterations = getLinearIterations(report)
% Getter for linear iterations
    % Get total nonlinear iterations if we actually solved the Newton step,
    % not just checked for convergence
    iterations = 0;
    if report.Solved
        if(isfield(report.LinearSolver,'Iterations'))
            iterations = report.LinearSolver.Iterations(1);
        end
    end
end

%-------------------------------------------------------------------------%
function time = getLinearSolverTime(report)
% Getter for linear solver time
    % Get total nonlinear iterations if we actually solved the Newton step,
    % not just checked for convergence
    time = 0;
    if report.Solved && isfield(report, 'LinearSolver') && isfield(report.LinearSolver, 'SolverTime')
        time = report.LinearSolver.SolverTime;
        time = time + report.LinearSolver.PostProcessTime;
    end
end
