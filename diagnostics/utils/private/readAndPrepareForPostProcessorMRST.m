function [G, data, Gs, valid_ix] = readAndPrepareForPostProcessorMRST(problem, steps, info, precomp)
% Utility function to create input data for flow diagnostics postprocessor 
% GUI from MRST simulation output.
%
% SYNOPSIS:
%  [G, data, Gs, valid_ix] = ...
%        readAndPrepareForPostProcessorMRST(problem, steps, info, precomp)
%
% DESCRIPTION:
%   Takes an MRST simulation packed problem and outputs the relevant data 
%   structures required as input to PostProcessDiagnostics flow diagnostics
%   viewer. 
%
% REQUIRED PARAMETERS:
%
%   problem - Simulation problem created using packSimulationProblem(). The
%               simulation should already have been run.
%
%   steps   - Timesteps which are to be displayed in the interactive gui
%
%   info    - structure containing the following fields:
%               .date - start date of the format 
%                [day month year] e.g. [1 1 2000] for the 1st January
%                2000.
%               .time - array of cumulative time for all timesteps in days
%   
%   precomp - Structure optionally containing array of precomputed
%               diagnostics data for each timestep to pass to the GUI. 
%               If empty, diagnostics will be calculated when the GUI is 
%               run.
%
% RETURNS:
%
%   G       - Grid structure with additional .cells.PORV field
%   data    - Data structure containing static, dynamic and computed
%                properties. (Computed properties not yet implemented.)
%   Gs      - Same as G for compatibility with other functions. 
%   valid_ix - index of states which are valid for flow diagnostics
%               calculations.
%   
% SEE ALSO:
%   PostProcessDiagnosticsMRST.m

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


model = problem.SimulatorSetup.model;
schedule = problem.SimulatorSetup.schedule;
G = model.G;
G.cells.PORV = poreVolume(G,model.rock);
Gs = G;

init.PORO.values = model.rock.poro;
if size(model.rock.perm, 2) == 3
    init.PERMX.values = model.rock.perm(:,1);
    init.PERMY.values = model.rock.perm(:,2);
    init.PERMZ.values = model.rock.perm(:,3);
else
    [init.PERMX.values, init.PERMY.values, init.PERMZ.values] = deal(model.rock.perm(:,1));
end
if G.griddim == 3
    init.DEPTH.values = G.cells.centroids(:,3);
else
    init.DEPTH.values = zeros(G.cells.num, 1);
end
init.PORV.values = G.cells.PORV;

% more fields from init might be interesting
data = setStatic([], init, {'PORO', 'PERMX', 'PERMY', 'PERMZ', 'DEPTH', 'PORV'});

% time in days (start and end of restart step)
startday = datenum(info.date(1, [3 2 1]));
data.time.prev = startday + info.time( max(steps-1,1) ) - info.time(1);
data.time.cur  = startday + info.time( steps ) - info.time(1);

statesHandler = problem.OutputHandlers.states;
states = statesHandler(steps);

for i = 1:numel(steps)
    data.wells{i} = schedule.control(schedule.step.control(i)).W;
end
data.wells = data.wells';

if isempty(precomp)
    data.states = states;
else
    data.states = cellfun(@(x)x.states{1}, precomp, 'UniformOutput', false);
end

valid_ix = isValidState(data.states);
if ~all(valid_ix)
    ns = nnz(~valid_ix);
    warning('Current version requires at least one open well.\n Skipping %2.0d of the %2.0d selected restart steps.', ns, numel(valid_ix));
    data.time.prev = data.time.prev(valid_ix);
    data.time.cur  = data.time.cur(valid_ix);
    data.states    = data.states(valid_ix);
    data.wells     = data.wells(valid_ix);
end

% include some more fields from restart later on
data = setDynamic(G, data, valid_ix);

% set empty computed-field
data.computed = repmat(struct('name', '', 'values', [], 'limits', []), [0 1]);

end 

function data = setStatic(~, init, propnames)
for k = 1:numel(propnames)
    if isfield(init, propnames{k})
        v = init.(propnames{k}).values;
        data.static(k) = struct('name',    propnames{k}, ...
                                'values' , v, ...
                                'limits' , [min(v), max(v)]);
    end
end
end

function data = setDynamic(G, data, valid_ix)
if ~any(valid_ix)
    data.dynamic = [];
else
    flds = {'PRESSURE', 'SWAT', 'SOIL', 'SGAS'};
    p = cellfun(@(x)x.pressure, data.states, 'UniformOutput', false);
    p = horzcat(p{:})/barsa;
    data.dynamic(1) = struct('name', 'PRESSURE', ...
        'values' , p, ...
        'limits' , [min(min(p)), max(max(p))]);
    
    for k = 1: size(data.states{1}.s, 2)
        nm = flds{k+1};
        vals = cellfun(@(x)x.s(:,k), data.states, 'UniformOutput', false);
        vals = horzcat(vals{:});
        % take min/max over all steps
        [minv, maxv] = deal(min(min(vals)), max(max(vals)));
        data.dynamic(k+1) = struct('name', nm, ...
            'values' , vals, ...
            'limits' , [minv, maxv]);
    end
    

end
end

function flag = isValidState(states)
flag = true(numel(states), 1);
for k = 1:numel(states)
    ws = states{k}.wellSol;
    stat = and([ws.status], abs([ws.val])>0);
    openPrd = any(and(stat, [ws.sign]<0));
    openInj = any(and(stat, [ws.sign]>0));
    flag(k) = openPrd || openInj;
end
end

