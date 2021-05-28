function varargout = computeEnsembleSimulationObjective(setupFn, varargin)
% This is an undocumented utility function for OptimizationProblem which
% runs a and saves similations/objectives/adjoints for requested
% realizations

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

opt = struct('outputPath',           '', ...
    'memberIx',             [], ...
    'schedule',             [], ...
    'objective',            [], ...
    'computeGradient',    true, ...
    'clearStatesAfterAdjoint', true);

[opt, objectiveOpts] = merge_options(opt, varargin{:});
if isempty(opt.outputPath)
    opt.outputPath = fullfile(mrstOutputDirectory(), 'tmp');
end

chk = numel(opt.memberIx) > 1 && nargout >0;
assert(~chk, 'Output only provided for single model evaluation');

computeGradient = opt.computeGradient || nargout == 2;

% schedule may be located in outputPath
schedule_input = opt.schedule;
if isempty(schedule_input) && isfile(fullfile(opt.outputPath, 'schedule.mat'))
    tmp = load(fullfile(opt.outputPath, 'schedule.mat'));
    schedule_input = tmp.schedule;
end

% % trajectories (if relevant) may be located in outputPath
% traj = opt.trajectories;
% if isempty(traj) && isfile(fullfile(opt.outputPath, 'trajectories.mat'))
%     tmp = load(fullfile(opt.outputPath, 'trajectories.mat'));
%     traj = tmp.trajectories;
% end

% list of objective functions (possibly empty)
obj = @(model, states, schedule, varargin)opt.objective(model, states, schedule, varargin{:}, objectiveOpts{:});


for k = 1:numel(opt.memberIx)
    ix = opt.memberIx(k);
    
    problem = setupFn(ix);
    problem.Directory = opt.outputPath;
    % update schedule
    if ~isempty(schedule_input)
        problem.SimulatorSetup.schedule = mergeSchedules(problem.SimulatorSetup.schedule, schedule_input);
    end
    %     % update trajectory
    %     if ~isempty(traj)
    %         problem.SimulatorSetup.schedule = updateWellTrajectories(problem.SimulatorSetup.schedule, traj);
    %     end
    % update handlers
    fn = fieldnames(problem.OutputHandlers);
    for kh = 1:numel(fn)
        problem.OutputHandlers.(fn{kh}).dataDirectory = opt.outputPath;
    end
    
    % check if ouput-folder exists
    pth = fullfile(problem.OutputHandlers.states.dataDirectory, problem.OutputHandlers.states.dataFolder);
    if ~isfolder(pth)
        mkdir(pth);
    end
    
    % run simulation
    ok = simulatePackedProblem(problem);
    
    % objectives
    if ~isempty(obj)
        model    = problem.SimulatorSetup.model;
        schedule = problem.SimulatorSetup.schedule;
        state0   = problem.SimulatorSetup.state0;
        ws       = problem.OutputHandlers.wellSols;
        states   = problem.OutputHandlers.states;
        
        value = [];
        if ok
            vals = obj(model, states, schedule);
            value  = sum(vertcat(vals{:}));
        end
        objFileNm = [func2str(opt.objective), '.mat'];
        save(fullfile(ws.getDataPath(), objFileNm), 'value');
        if nargout > 0
            varargout{1} = value;
        end
        % adjoints
        if computeGradient
            gradient = [];
            if ok
                objh = @(tstep, model, state)obj(model, states, schedule, 'ComputePartials', true, 'tStep', tstep,'state', state);
                gradient = computeGradientAdjointAD(state0, states, model, schedule, objh, 'OutputPerTimestep', true);%, 'LinearSolver', linSolve);
                gradient = processAdjointGradients(gradient, ws, 'controlIx', schedule.step.control);
            end
            gradFileNm = ['gradient_', func2str(opt.objective), '.mat'];
            save(fullfile(ws.getDataPath(), gradFileNm), 'gradient');
            if nargout > 1
                varargout{2} = gradient;
            end
            if opt.clearStatesAfterAdjoint
                states.resetData();
            end
        end
        
    end
end
end

function schedule = mergeSchedules(schedule, schedule_tmp)
% merge targets/limits of schedule_tmp into schedule
for kc = 1:numel(schedule.control)
    W  = schedule.control(kc).W;
    Wt = schedule_tmp.control(kc).W;
    for kw = 1:numel(W)
        assert(strcmp(W(kw).type, Wt(kw).type), 'Can''t  merge wells of different types')
        W(kw).val  = Wt(kw).val;
        W(kw).lims = Wt(kw).lims;
    end
    schedule.control(kc).W = W;
end
end
