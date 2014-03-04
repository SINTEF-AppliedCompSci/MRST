function varargout = ...
   runScheduleADI(initState, G, rock, system, schedule, varargin)
% Given a schedule and a system, solve for all time steps
%
% SYNOPSIS:
%  [wellSols, states, its] = runScheduleADI(initState, G, rock, system, schedule)
%  [wellSols, states, its] = runScheduleADI(initState, G, rock, system, ...
%                                           schedule, 'pn', pv, ...)
% PARAMETERS:
%   initState - The initial state at t = 0;
%
%   G         - A valid grid. See grid_structure.
%
%   rock      - A valid rock structure. Should contain an Nx1 array
%               'poro' containing cell wise porosity values. A permeability
%               field is not *needed* for all the ad-fi solvers as they can
%               work directly with transmissibilities, but it is
%               highly recommended to supply them in either a Nx1 or Nx3
%               array. N is here equal to G.cells.num.
%
%  system     - System configuration as defined by initADISystem.
%
%  schedule   - Schedule (usually found in the deck.SCHEDULE field from
%               the return value of readEclipseDeck from the deckformat
%               module). This fully defines the well configurations for all
%               timesteps.
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%
%   Verbose - If verbose output should be outputted. Defaults to
%             mrstVerbose.
%
%   writeOutput - Save output to the cache folder. This can be practical
%                 when states becomes too big to solve in memory or when
%                 running adjoint simulations.
%
%   outputName  - The string which prefixes .mat files written if
%                 writeOutput is enabled. Defaults to 'state'.
%
%
% RETURNS:
%   wellSols - Well solution struct for each timestep. Cellarray of size
%              Ntx1.
%
%   states (OPTIONAL) - State solution struct for each timestep. Cellarray
%                       of size Ntx1. Note that as this can be come
%                       prohibitively big for long simulations this should
%                       be only outputted if neded.
%
%   its (OPTIONAL) - Nonlinear iteration count for each timestep.
%
%
% SEE ALSO:
%   solvefiADI

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

opt = struct('Verbose'       , mrstVerbose      , ...
             'writeOutput'   , false            , ...
             'outputName'    , 'state'          , ...
             'scaling'       , []               , ...
             'startAt'       , 1                , ...
             'outputDir'     , default_outputDir, ...
             'plotCallback'  , [],  ...
             'outputNameFunc', []);

opt = merge_options(opt, varargin{:});

vb = opt.Verbose;
outputStates =      nargout > 1;
outputIter =        nargout > 2;
outputConvergence = nargout > 3;

%--------------------------------------------------------------------------

dt = schedule.step.val;
tm = cumsum(dt);
dispif(vb, '*****************************************************************\n')
dispif(vb, '********** Starting simulation: %5.0f steps, %5.0f days *********\n', numel(dt), tm(end)/day)
dispif(vb, '*****************************************************************\n')

%--------------------------------------------------------------------------

if opt.writeOutput
    %delete existing output
    % delete(fullfile(directory, [opt.outputName, '*.mat']));
    % output file-names

    if ~isdir(opt.outputDir),
       [success, msg, id] = mkdir(opt.outputDir);
       if ~ success,
          error(id, 'Failed to create output diretory ''%s'': %s', ...
                opt.outputDir, msg);
       end
    end

    if isempty(opt.outputNameFunc)
       outNm  = @(tstep)fullfile(opt.outputDir, [opt.outputName, sprintf('%05.0f', tstep)]);
    else
       outNm  = @(tstep)fullfile(opt.outputDir, opt.outputNameFunc(tstep));
    end
end


%--------------------------------------------------------------------------
% output initState
state = initState;
if opt.writeOutput, save(outNm(0), 'state');end;

%--------------------------------------------------------------------------
% collect all wellsols in cell wellsols
wellSols = cell(numel(dt), 1);

if outputStates
    states = cell(numel(dt)+1, 1);
    initState.wellSol = initWellSolLocal([], state);
    states{1} = initState;
end

iter = zeros(numel(dt),1);

%--------------------------------------------------------------------------
% default is to report all steps
if ~isfield(schedule.step, 'repStep')
    schedule.step.repStep = true(numel(dt), 1);
end

prevControl = -1;
timero = tic;
repStep = 0;
convergence = [];
useMrstSchedule = isfield(schedule.control(1), 'W');
if useMrstSchedule
    nw = arrayfun(@(x)numel(x.W), schedule.control);
    uniformSchedule = all(nw == nw(1));
    if uniformSchedule
        W_temp = schedule.control(1).W;
        for k=1:numel(W_temp), W_temp(k).status = true; end
        wellSol_init = initWellSolLocal(W_temp, state);  % initial guess (udated after each time-step)
        for k=1:numel(W_temp), W_temp(k).status = false; end
        wellSol_zero = initWellSolLocal(W_temp, state);  % default 0-well-sol
    end
end


for tstep = 1:numel(dt)
    dispif(vb, 'Time step %5.0f of %d\n', tstep, numel(dt));
    control = schedule.step.control(tstep);
    if control ~= prevControl
        if control == 0, % when is control == 0 ?
           W = processWellsLocal(G, rock, [], 'createDefaultWell', true);
        else
           if ~useMrstSchedule
               W = processWellsLocal(G, rock, schedule.control(control), ...
                                 'Verbose', opt.Verbose, ...
                                 'DepthReorder', false);
           else
               W = schedule.control(control).W;
           end
           openWells = vertcat(W.status);
           assert(all(islogical(openWells)));% avoid errors due to setting status to 1;
           W = W(openWells);
        end
    end
    dispif(vb, 'Time step length: %g day.\n', convertTo(dt(tstep), day))
    state0 = state;
    if useMrstSchedule && uniformSchedule
        state0.wellSol = initWellSolLocal(W, state, wellSol_init(openWells));
    else
        state0.wellSol = initWellSolLocal(W, state);
    end

    [state, its, conv] = solvefiADI(state0, dt(tstep), W, G, system);
    % check if any controls have been switched, and if so update W
    W = updateSwitchedControls(state.wellSol, W);
    wellSols{tstep} = state.wellSol;
    wellSols{tstep} = addWellInfo(wellSols{tstep}, W);

    iter(tstep) = its;
    if useMrstSchedule && uniformSchedule
        wellSol_init(openWells) = state.wellSol;
        ws = wellSol_zero;
        ws(openWells) = state.wellSol;
        state.wellSol = ws;
    end
    wellSols{tstep} = state.wellSol;
    wellSols{tstep} = addWellInfo(wellSols{tstep}, W);
    if outputStates
        states{tstep + 1} = state;
    end
    if ~isempty(opt.plotCallback)
        opt.plotCallback(G, state)
    end

    prevControl = control;
    if opt.writeOutput && schedule.step.repStep(tstep)
        repStep = repStep + 1;
        save(outNm(repStep), 'state');
    end
    convergence = [convergence; conv]; %#ok
    dispif(~opt.Verbose, 'Step %4g of %4g (Used %3g iterations)\n', ...
           tstep, numel(dt), its);
end

timend = toc(timero);
dispif(vb, ['************Simulation done: %7.2f seconds ', ...
            '********************\n'], timend)
varargout{1} = wellSols;

if opt.writeOutput
   save(fullfile(opt.outputDir, 'wellSols'), 'wellSols');
end

if outputStates
    varargout{2} = states;
end

if outputIter
    varargout{3} = iter;
end

if outputConvergence
    varargout{4} = convergence;
end

end

%--------------------------------------------------------------------------

function wellSol = addWellInfo(wellSol, W)
   %nm = fieldnames(W);
   nm = {'name', 'sign'};
   for k = 1:numel(nm)
      for wnum = 1:numel(W)
         wellSol(wnum).(nm{k}) = W(wnum).(nm{k});
      end
   end
end
