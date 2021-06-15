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


default_outputDir = fullfile(fileparts(mfilename('fullpath')), 'cache');

opt = struct('Verbose', mrstVerbose, 'writeOutput', false, ...
             'outputName', 'state', 'scaling', [], 'startAt', 1, ...
             'outputDir', default_outputDir, 'Wext',[],'fluxes',false,...
             'dt_min', inf,'force_step',true,'targetIts',6, ...
             'report_all',true, 'bc', [], 'tstep_hook_fn', []);

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
end

total_stepcount = 1;


%--------------------------------------------------------------------------
dt = schedule.step.val;
tm = cumsum(dt);
dispif(vb, '*****************************************************************\n')
dispif(vb, '********** Starting simulation: %5.0f steps, %5.0f days *********\n', numel(dt), tm(end)/day)
dispif(vb, '*****************************************************************\n')
%--------------------------------------------------------------------------


if opt.writeOutput | (opt.startAt ~= 1)
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
dt_prev_ok=min(dt(opt.startAt),opt.dt_min);
for tstep = opt.startAt:numel(dt)
    dispif(vb, 'Global time step %5.0f of %d\n', tstep, numel(dt));
    control = schedule.step.control(tstep);
    dt_already_reported = false;

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
    prevControl = control;
    not_converged = true;
    dispif(vb, sprintf('Global time step length: %g day.\n', convertTo(dt(tstep), day)))

    if(opt.force_step)
        state0 = state;
        [state, its] = solvefiADI(state, dt(tstep), W, G, system, 'bc', opt.bc);
        t=t+dt(tstep);
    else
        % local timestepping
        t_loc=0;
        dt_loc=min(dt_prev_ok,dt(tstep));
        while t_loc< dt(tstep)
            state0 = state;
            dispif(vb, sprintf('Local time step length: %g day.\n', convertTo(dt_loc, day)))
            assert(dt_loc>0)
            [state, its, converged] = solvefiADI(state, dt_loc, W, G, system, ...
                                                 'bc', opt.bc);
            % figure(99),clf;
            % plotCellData(G,state.s(:,1)),plotWell(G,W)
            if(~converged.converged)
                if(dt_loc>opt.dt_min)
                    warning('Cutting timestep')
                    dt_loc=dt_loc/2;
                    state=state0;
                    continue;
                else
                    warning('Time step cutting not successful, continue anyway')
                end
            end
            t=t+dt_loc;
            t_loc=t_loc+dt_loc;
            fprintf('t = %g days\n', convertTo(t, day));
            if(opt.report_all)
                % report also for local timesteps and produce schedule accordingly
                [iter, wellSols, states] = ...
                    poststepHousekeeping(its, total_stepcount, state0, state, system, G, W, iter, ...
                                         dt_loc, wellSols, states, outputStates, opt.fluxes);

                schedulenew.step.control(total_stepcount)=control;
                schedulenew.step.val(total_stepcount)=dt_loc;
                dispif(~opt.Verbose, 'Step %4g of %3.2f (Used %3g iterations)\n', ...
                       total_stepcount, t/T_all, its);
                
                total_stepcount = total_stepcount + 1;
                dt_already_reported = true; % no need to report again at the
                                            % end of the global loop
            end
            dt_history=[];
            [dt_new, dt_history] = simpleStepSelector(dt_history, dt_loc, its,...
                'targetIts', opt.targetIts, ...
                'dt_min', opt.dt_min, ...
                'dt_max', dt(tstep), ...
                'stepModifier', 1.5);
            %[dt, dt_history] = simpleStepSelector(dt_history, dt_loc, its, varargin);
            if(t_loc<dt(tstep))
                if(its<system.nonlinear.maxIterations)
                    dt_prev_ok=dt_loc;
                end
            else
               if(converged.converged)
                  if(dt_prev_ok >= dt(tstep)/1.5 || dt_new >= dt(tstep)/1.5 )
                  %if(dt_loc>dt(tstep)/2)
                       dt_prev_ok = dt(tstep);
                  end
               end
            end
            dt_loc=min(dt_new,dt(tstep)-t_loc);
        end
    end

    if(~dt_already_reported)
        [iter, wellSols, states] = poststepHousekeeping(its, total_stepcount, state0, ...
                                                        state, system, G, W, ...
                                                        iter, dt(total_stepcount), ...
                                                        wellSols, states, ...
                                                        outputStates, opt.fluxes);

        dispif(~opt.Verbose, 'Step %4g of %4g (Used %3g iterations)\n', tstep, numel(dt), its);
        total_stepcount = total_stepcount + 1;
    end

    % Stuff to do once per 'tstep' (regardless of whether local stepping is used)
    if opt.writeOutput && schedule.step.repStep(tstep)
        repStep = repStep + 1;
        save(outNm(tstep), 'state', 'W');
    end
    % Running hook function if provided
    if ~isempty(opt.tstep_hook_fn)
        opt.tstep_hook_fn(state, tstep);
        drawnow;
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

% ----------------------------------------------------------------------------
function wellSol = addWellInfo(wellSol, W)
    nm = {'name', 'sign'};
    for k = 1:numel(nm)
        for wnum = 1:numel(W)
            wellSol(wnum).(nm{k}) = W(wnum).(nm{k});
        end
    end
end

% ----------------------------------------------------------------------------
function [iter, wellSols, states] = ...
        poststepHousekeeping(its, step, state0, state, system, G, W, iter, ...
                             dt, wellSols, states, outputStates, outputFluxes)

    iter(step) = its;
    wellSols{step} = addWellInfo(state.wellSol, W);
    if outputStates
        states{step+1} = state;
    end
    if outputFluxes
        [eqs, flux] = system.getEquations(state0, state, dt, G, W, system.s, ...
                                                  system.fluid, 'bc', opt.bc, ...
                                                  'fluxes', true, 'stepOptions', ...
                                                  system.stepOptions);
        states{step+1}.fluxes = flux;
    end
    
end


