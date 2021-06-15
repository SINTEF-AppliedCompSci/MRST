function [ws, states, reports, model] = getPackedSimulatorOutput(problem, varargin)
%Get output from a packed simulation problem
%
% SYNOPSIS:
%   [ws, states, reports] = getPackedSimulatorOutput(problem)
%
% REQUIRED PARAMETERS:
%   problem - Problem generated using 'packSimulationProblem'.
%
% OPTIONAL PARAMETERS:
%   readFromDisk - Indicating if states and reporst will be read from disk,
%                  or returned as OutputHandler instances. Reading from
%                  disk can take some time and is recommended when further
%                  changes to output is desired.
%   readWellSolsFromDisk - See above. Applies for wellSols only.
%   readReportsFromDisk  - See above. Applies for reports only.
%
% RETURNS:
%   ws      - Well output.
%   states  - States for each simulated step.
%   reports - Report for each of the time-steps.
% 
% EXAMPLE:
%   demoPackedProblems
%
% SEE ALSO:
%   packSimulationProblem, getMultiplePackedSimulatorOutputs

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

    opt = struct('readFromDisk', true, ...
                 'readWellSolsFromDisk', [], ...
                 'readStatesFromDisk',   [], ...
                 'readReportsFromDisk',  []);
    opt = merge_options(opt, varargin{:});
    
    if nargout > 3
        ctrl = problem.SimulatorSetup.schedule.control(1);
        [~, fstruct] = problem.SimulatorSetup.model.getDrivingForces(ctrl);
        model = problem.SimulatorSetup.model.validateModel(fstruct);
    end
    
    if isempty(opt.readWellSolsFromDisk)
        opt.readWellSolsFromDisk = opt.readFromDisk;
    end
    if isempty(opt.readReportsFromDisk)
        opt.readReportsFromDisk = opt.readFromDisk;
    end
    if isempty(opt.readStatesFromDisk)
        opt.readStatesFromDisk = opt.readFromDisk;
    end
    dt = problem.SimulatorSetup.schedule.step.val;
    nstep = numel(dt);
    
    sh = problem.OutputHandlers.states;
    wh = problem.OutputHandlers.wellSols;
    rh = problem.OutputHandlers.reports;
    
    ndata = sh.numelData();
    ws = cell(wh.numelData(), 1);
    states = cell(sh.numelData(),1);
    reports = cell(rh.numelData(),1);
    wantWells = false;
    for i = 1:numel(problem.SimulatorSetup.schedule.control)
        ctrl = problem.SimulatorSetup.schedule.control(i);
        if isfield(ctrl, 'W') && ~isempty(ctrl.W)
            wantWells = true;
            break
        end
    end

    sn = sprintf('%s (%s)', problem.BaseName, problem.Name);
    if nstep == ndata
        fprintf('Found complete data for %s: %d steps present\n', sn, ndata);
    elseif ndata > nstep
        endstate = states{ndata};
        if isfield(endstate, 'time') && endstate.time > sum(dt)
            warning('Found too much data for %s: %d of %d steps present. Case may have been redefined!\n', sn, ndata, nstep);
        else
            fprintf('Found complete data for %s: %d steps present\n', sn, ndata);
        end
    elseif ndata > 0
        fprintf('Found partial data for %s: %d of %d steps present\n', sn, ndata, nstep);
    else
        fprintf('Did not find data for %s\n', sn);
    end
    wellOutputMissing = wantWells && wh.numelData() == 0;
    ns = sh.numelData();
    for i = 1:ns
        if nargout > 1 && opt.readStatesFromDisk
            states{i} = sh{i};
        end
         if wantWells && opt.readWellSolsFromDisk
            if wellOutputMissing
                if isempty(states{i})
                    ws{i} = sh{i}.wellSol;
                else
                    ws{i} = states{i}.wellSol;
                end
            end
         end
    end
    nw = wh.numelData();
    for i = 1:nw
        if wantWells && opt.readWellSolsFromDisk
            if ~wellOutputMissing
                try
                    ws{i} = wh{i};
                catch
                    ws{i} = states{i}.wellSol;
                end
            end
           
        end
    end
    nr = rh.numelData();
    for i = 1:nr
        if nargout > 2 && opt.readReportsFromDisk
            try
                reports{i} = rh{i};
            catch
                reports{i} = [];
            end
        end
    end

    if ~opt.readStatesFromDisk
        % Just return handlers instead
        states = sh;
    end
    if ~opt.readReportsFromDisk
        reports = rh;
    end
    if ~opt.readWellSolsFromDisk
        ws = wh;
    end
end