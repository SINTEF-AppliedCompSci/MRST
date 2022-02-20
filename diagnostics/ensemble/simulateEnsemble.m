function h = simulateEnsemble(setupFn, varargin)
%Undocumented Utility Function

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

opt = struct('schedule',                 [], ...
             'handlers',                 [], ...
             'outputDir',                '', ...
             'ensembleSize',            inf, ...
             'memberIx',                 [], ...
             'overWrite',             false, ...       
             'writeToDiskMode',       'all', ...
             'deleteStates',          false); 

[opt, extra] = merge_options(opt, varargin{:});         
memberIx = opt.memberIx;
if isempty(memberIx)
    nMem = min(opt.ensembleSize, 1000); % just something big
    memberIx = (1:nMem)';
end

h = opt.handlers;

for k = 1:numel(memberIx)
    ix = memberIx(k);
    % setup
    tmp = setupFn(ix);
    if isempty(tmp.model), break; end
    if ~isfield(tmp, 'state0') || isempty(tmp.state0)
        dispif(mrstVerbose, 'Initial state not provided, using default');
        state0 = initResSol(model.G, opt.p_ref);
    else
        state0 = tmp.state0;
    end
    [model, W] = deal(tmp.model, tmp.W);
    schedule = [];
    if isfield(tmp, 'schedule')
        schedule = tmp.schedule;
    end
    
    % handle schedule, if opt.schedule is non-empty, merge with W
    if ~isempty(opt.schedule)
        schedule = getSchedule(opt.schedule, W);
    end
    
    simOpts = {};
    if isfield(tmp, 'simOpts')
        simOpts = tmp.simOpts;
    end
    
    % handle result-handlers (if not input)
    h = setupHandlers(h, ix, opt.outputDir, opt.writeToDiskMode);
    
    % check if previously (partly) simulated
    nstep = numel(schedule.step.val);
    ndata = h.ws{ix}.numelData();
    if ndata == nstep
        dispif(mrstVerbose, 'Complete output found for realization %d', ix);
        continue
    end
    
    restartStep = h.states{ix}.numelData()+1;
    if restartStep > 1
        state0 = h.states{ix}.getRestart;
    end
    
    try
        simulateScheduleAD(state0, model, schedule, 'restartStep', restartStep,...
                           'OutputHandler', h.states{ix}, ...
                           'WellOutputHandler', h.ws{ix}, ...
                           'ReportHandler', h.reports{ix}, simOpts{:});
        if opt.deleteStates
            h.states{ix}.resetData;
        end
    catch ex
        fprintf('!!! Simulation resulted in fatal error !!!\n Exception thrown: %s\n', ex.message);
    end
end
end

%--------------------------------------------------------------------------

function schedule_out = getSchedule(schedule, W)
% merge schedule and W -> schedule_out

% allow passing schedule as mat-file
if ischar(schedule)
    assert(isfile(schedule));
    tmp = load(schedule);
    schedule = tmp.schedule;
end
schedule_out = schedule;

flds0 = {'type', 'val', 'sign', 'status', 'cstatus', 'lims'};
for cn = 1:numel(schedule.control)
    control = schedule.control(cn);
    flds = flds0(isfield(control.W, flds0));
    assert(numel(W)==numel(control.W), 'Current well-structure is not compatible with schedule');
    schedule_out.control(cn).W = W;
    for wn = 1:numel(W)
        assert(strcmp(W(wn).name, control.W(wn).name), 'Unexpected mismatch, check well ordering');
        for j = 1:numel(flds)
            schedule_out.control(cn).W(wn).(flds{j}) = control.W(wn).(flds{j});
        end
    end
end
end

%--------------------------------------------------------------------------

function h = setupHandlers(h, ix, outputDir, mode)
if ~isfield(h, 'states') || numel(h.states) < ix
    writeToDisk  = strcmp(mode, 'all');
    h.states{ix} = ResultHandler('dataDirectory', outputDir, 'dataFolder', sprintf('realization%4.4d', ix), ...
                        'writeToDisk', writeToDisk, 'dataPrefix', 'states');
end

if ~isfield(h, 'ws') || numel(h.ws) < ix
    h.ws{ix} = ResultHandler('dataDirectory', outputDir, 'dataFolder', sprintf('realization%4.4d', ix), ...
                             'writeToDisk', true, 'dataPrefix', 'ws');
end    

if ~isfield(h, 'reports') || numel(h.reports) < ix
    writeToDisk  = strcmp(mode, 'all');
    h.reports{ix} = ResultHandler('dataDirectory', outputDir, 'dataFolder', sprintf('realization%4.4d', ix), ...
                                  'writeToDisk', writeToDisk, 'dataPrefix', 'reports');
end  
end
