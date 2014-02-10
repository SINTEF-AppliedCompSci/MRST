function varargout = runMrstADI(initState, G, system, schedule, varargin)
% varargout = runMrstADI(initState, G, system, schedule, varargin)
% Given a schedule and a system, solve for all time steps
%
% SYNOPSIS:
% PARAMETERS:
%
%
% RETURNS:
%
% SEE ALSO:
%   solvefiADI runScheduleADI

%{
Copyright 2009-2014 SINTEF ICT, Applied Mathematics.

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


default_outputDir = fullfile(fileparts(mfilename('fullpath')), 'cache');

opt = struct('Verbose', mrstVerbose, 'writeOutput', false, 'outputName', 'state', 'scaling', [], ...
    'startAt', 1, 'outputDir', default_outputDir, 'Wext',[],'fluxes',false,...
    'dt_min', inf,'force_step',true,'targetIts',6,'report_all',true);

opt = merge_options(opt, varargin{:});

meta = struct('converged'   , false, ...
                  'stopped'     , false, ...
                  'wellschanged', false, ...
                  'relax'       , system.nonlinear.relaxation, ...
                  'stagnate'    , false, ...
                  'iteration'   , 0);
vb = opt.Verbose;
outputStates = nargout > 1;
outputIter = nargout > 2;
if(opt.report_all)
    schedulenew=schedule;
    tstep_all=1;
end



%--------------------------------------------------------------------------
dt = schedule.step.val;
tm = cumsum(dt);
dispif(vb, '*****************************************************************\n')
dispif(vb, '********** Starting simulation: %5.0f steps, %5.0f days *********\n', numel(dt), tm(end)/day)
dispif(vb, '*****************************************************************\n')
%--------------------------------------------------------------------------


if opt.writeOutput
    directory = fullfile(fileparts(mfilename('fullpath')), 'cache');
    %delete existing output
    delete(fullfile(directory, [opt.outputName, '*.mat']));
    % output file-names
    outNm  = @(tstep)fullfile(opt.outputDir, [opt.outputName, sprintf('%05.0f', tstep)]);
end


%--------------------------------------------------------------------------
% output initState
if opt.startAt == 1
    state = initState;
    if opt.writeOutput, save(outNm(0), 'state');end
else
    rf = load(outNm(opt.startAt-1));
    state = rf.state;
end

%--------------------------------------------------------------------------
% collect all wellsols in cell wellsols
wellSols = cell(numel(dt), 1);

if outputStates
    states = cell(numel(dt)+1, 1);
    states{1} = state;
end

iter = zeros(numel(dt),1);

%--------------------------------------------------------------------------
% default is to report all steps
if ~isfield(schedule.step, 'repStep')
    schedule.step.repStep = true(numel(dt), 1);
end

T_all=sum(schedule.step.val);
t=0;
prevControl = -1;
timero = tic;
repStep = 0;
W = [];
dt_prev_ok=min(dt(1),opt.dt_min);
for tstep = opt.startAt:numel(dt)
    dispif(vb, 'Global time step %5.0f of %d\n', tstep, numel(dt));
    control = schedule.step.control(tstep);

    if(control==0)
        W=[];
    else
        W=schedule.W{control};
    end
    if control~=prevControl
        if(isfield(state,'wellSol'))
            state = rmfield(state,'wellSol');
        end
        state.wellSol=initWellSolLocal(W, state);
    end
    not_converged = true;
    dispif(vb, sprintf('Global time step length: %g day.\n', convertTo(dt(tstep), day)))

    if(opt.force_step)
        state0 = state;
        [state, its] = solvefiADI(state, dt(tstep), W, G, system);
        t=t+dt(tstep);
    else
        % local timestepping
        t_loc=0;
        dt_loc=dt_prev_ok;
        while t_loc< dt(tstep)
            state0 = state;
            dispif(vb, sprintf('Local time step length: %g day.\n', convertTo(dt_loc, day)))
            assert(dt_loc>0)
            [state, its, converged] = solvefiADI(state, dt_loc, W, G, system);
            % figure(99),clf;
            % plotCellData(G,state.s(:,1)),plotWell(G,W)
            if(~converged.converged)
                if(dt_loc>opt.dt_min)
                    warning('Cutting timestep')
                    dt_loc=dt_loc/2;
                    state=state0;
                    continue;
                else
                    warning('Time step cutting not successfull, continue anyway')
                end
            end
            t=t+dt_loc;
            if(opt.report_all)
                % report all steps and produce schedual accordingly
                iter(tstep_all) = its;
                wellSols{tstep_all} = state.wellSol;
                wellSols{tstep_all} = addWellInfo(wellSols{tstep_all}, W);
                if outputStates
                    states{tstep_all + 1} = state;
                end
                if(opt.fluxes)
                    [eqs, flux] = system.getEquations(state0, state, dt_loc, G, W, system.s, system.fluid,'bc',[],'fluxes',true);
                    states{tstep_all +1}.fluxes=flux;
                end
                if opt.writeOutput && schedule.step.repStep(tstep)
                    repStep = repStep + 1;
                    save(outNm(tstep), 'state', 'W');
                end
                schedulenew.step.control(tstep_all)=control;
                schedulenew.step.val(tstep_all)=dt_loc;
                dispif(~opt.Verbose, 'Step %4g of %3.2f (Used %3g iterations)\n', tstep_all, t/T_all, its);
                tstep_all=tstep_all+1;
            end
            dt_history=[];
            [dt_new, dt_history] = simpleStepSelector(dt_history, dt_loc, its,...
                'targetIts', opt.targetIts, ...
                'dt_min', opt.dt_min, ...
                'dt_max', dt(tstep), ...
                'stepModifier', 1.5);
            %[dt, dt_history] = simpleStepSelector(dt_history, dt_loc, its, varargin);
            t_loc=t_loc+dt_loc;
            if(t_loc<dt(tstep))
                if(its<10)
                    dt_prev_ok=dt_loc;
                end
            else
               if(its<10)
                  if(dt_loc>dt(tstep)/2)
                       dt_prev_ok = dt(tstep);
                  end
               end
            end
            dt_loc=min(dt_new,dt(tstep)-t_loc);
        end
    end
    if(~opt.report_all || opt.force_step)
        iter(tstep) = its;
        wellSols{tstep} = state.wellSol;
        wellSols{tstep} = addWellInfo(wellSols{tstep}, W);
        if outputStates
            states{tstep + 1} = state;
        end
        if(opt.fluxes)
            [eqs,flux] = system.getEquations(state0, state, dt(tstep), G, W, system.s, system.fluid,'bc',[],'fluxes',true,'stepOptions', system.stepOptions);
            states{tstep +1}.fluxes=flux;
        end
        %prevControl = control;
        if opt.writeOutput && schedule.step.repStep(tstep)
            repStep = repStep + 1;
            save(outNm(tstep), 'state', 'W');
        end
        dispif(~opt.Verbose, 'Step %4g of %4g (Used %3g iterations)\n', tstep, numel(dt), its);
    end

end

dispif(vb, '************Simulation done: %7.2f seconds ********************\n', toc(timero))
varargout{1} = wellSols;
if outputStates
    varargout{2} = states;
end

if outputIter
    varargout{3} = iter;
end
if opt.report_all
   varargout{4} = schedulenew;
end
end

function wellSol = addWellInfo(wellSol, W)
%nm = fieldnames(W);
nm = {'name', 'sign'};
for k = 1:numel(nm)
    for wnum = 1:numel(W)
        wellSol(wnum).(nm{k}) = W(wnum).(nm{k});
    end
end
end


