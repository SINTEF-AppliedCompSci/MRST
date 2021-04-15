function varargout = runAdjointADI(G, rock, fluid, schedule, objective, system, varargin)
% Compute adjoint gradients for a schedule using the fully implicit ad solvers
%
% SYNOPSIS:
%   grad = runAdjointADI(G, rock,s, f, schedule,
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
%  fluid      - Fluid as defined by initDeckADIFluid.
%
%  schedule   - Schedule (usually found in the deck.SCHEDULE field from
%               the return value of readEclipseDeck from the deckformat
%               module). This fully defines the well configurations for all
%               timesteps.
%
%  objective  - Function handle to the objective function with compute
%               partials enabled. The interface should be of the format
%               obj = @(tstep)NPVOW(G, wellSols, schedule, ...
%                                 'ComputePartials', true, 'tStep', tstep);
%               Where wellSols correspond to a foward simulation of the
%               same schedule.
%
%  system     - System configuration as defined by initADISystem.
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
%                 writeOutput is enabled. Defaults to 'adjoint'.
%
% simOutputName - The name used for outputName for the corresponding
%                 forward simulation. Is required if ForwardStates is not
%                 supplied.
%
% ForwardStates - A cell array of states corresponding to the forward run
%                 of the same schedule. Should be given if the states are
%                 not to be read from file using simOutputName.
%
% ControlVariables - Indices of the control variables. Default to the last
%                 variable in the implicit system, i.e. typically well
%                 closure variables.
%
% scaling       - Struct containing scaling factors .rate and .pressure.
%                 Defaults to 1 and 1.
% RETURNS:
%
%  grad         - gradient vector with respect to control step.
%
%  gradFull {OPTIONAL} - gradient vector with respect to time step.
%
% COMMENTS:
%
%
% SEE ALSO:
%

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


   directory = fullfile(fileparts(mfilename('fullpath')), 'cache');

   opt = struct('Verbose',             mrstVerbose, ...
                'writeOutput',         false, ...
                'directory',           directory, ...
                'outputName',          'adjoint', ...
                'outputNameFunc',      [], ...
                'simOutputName',       'state', ...
                'simOutputNameFunc',   [], ...
                'ForwardStates',       [], ...
                'ControlVariables',    [], ...
                'scaling',             [], ...
                'timeStepGradient',    false);

   opt = merge_options(opt, varargin{:});



   vb = opt.Verbose;
   states = opt.ForwardStates;

   if ~isempty(opt.scaling)
      scalFacs = opt.scaling;
   else
      scalFacs.rate = 1; scalFacs.pressure = 1;
   end

   %--------------------------------------------------------------------------
   dts = schedule.step.val;
   tm = cumsum(dts);
   dispif(vb, '*****************************************************************\n')
   dispif(vb, '**** Starting adjoint simulation: %5.0f steps, %5.0f days *******\n', numel(dts), tm(end)/day)
   dispif(vb, '*****************************************************************\n')
   %--------------------------------------------------------------------------

   if opt.writeOutput
      % delete existing output
      % delete([opt.outputName, '*.mat']);
      % output file-names
      if isempty(opt.outputNameFunc)
         outNm  = @(tstep)fullfile(opt.directory, [opt.outputName, sprintf('%05.0f', tstep)]);
      else
         outNm  = @(tstep)fullfile(opt.directory, opt.outputNameFunc(tstep));
      end
   end

   if isempty(opt.simOutputNameFunc)
      inNm  = @(tstep)fullfile(opt.directory, [opt.simOutputName, sprintf('%05.0f', tstep)]);
   else
      inNm  = @(tstep)fullfile(opt.directory, opt.simOutputNameFunc(tstep));
   end

   nsteps = numel(dts);

   prevControl = inf;
   adjVec    = [];
   eqs_p       = [];
   gradFull = cell(1, numel(dts));

   % Load state - either from passed states in memory or on disk
   % representation.
   state = loadState(states, inNm, numel(dts));
   % We strip wellSols for closed wells,
   state.wellSol = state.wellSol(vertcat(state.wellSol.status) == 1);

   timero = tic;
   useMrstSchedule = isfield(schedule.control(1), 'W');
   for tstep = numel(dts):-1:1
      dispif(vb, 'Time step: %5.0f\n', tstep); timeri = tic;
      control = schedule.step.control(tstep);
      if control~=prevControl
        if (control==0)
           W = processWellsLocal(G, rock, [], 'createDefaultWell', true);
        else
            if ~useMrstSchedule
                W = processWells(G, rock, schedule.control(control), 'Verbose', opt.Verbose, ...
                                      'DepthReorder', true);
            else
                W = schedule.control(control).W;
            end
        end
        openWells = vertcat(W.status);
        W = W(openWells);
      end

      state_m = loadState(states, inNm, tstep-1);
      % We strip wellSols for closed wells,
      state_m.wellSol = state_m.wellSol(vertcat(state_m.wellSol.status) == 1);

      % For adjoint computation, we do not allow control switching.
      system.well.allowControlSwitching = false;

      % Sanity check: wellSols and W have same type. We also set equal val. Note that
      % val should actually not be used in adjoint simulation.
      for wn = 1:numel(W)
         assert(strcmp(W(wn).type, state.wellSol(wn).type), ['In runAdjointADI, well types do not ' ...
                             'match']);
         W(wn).val = state.wellSol(wn).val;
      end


      eqs   = system.getEquations(state_m, state  , dts(tstep), G, W, system, fluid, 'scaling', ...
                                           scalFacs,  'stepOptions', ...
                                           system.stepOptions, 'iteration', inf);
      if tstep < nsteps
         eqs_p = system.getEquations(state, state_p, dts(tstep+1), G, W_p, system, fluid, ...
                                              'reverseMode', true, 'scaling', scalFacs, 'stepOptions', ...
                                              system.stepOptions, 'iteration', inf);
      end

      [adjVec, ii] = solveAdjointEqsADI(eqs, eqs_p, adjVec, objective(tstep), system);
      if isempty(opt.ControlVariables)
         gradFull{tstep} = zeros(numel(openWells), 1);
         gradFull{tstep}(openWells) = -adjVec(ii(end,1):ii(end,2));
      else
         gradFull{tstep} = -adjVec(mcolon(ii(opt.ControlVariables,1),ii(opt.ControlVariables,2)));
      end

      state_p = state;
      state   = state_m;
      W_p = W;
      prevControl = control;
      dispif(vb, 'Done %6.2f\n', toc(timeri))
      if opt.writeOutput
         save(outNm(tstep), 'adjVec');
      end
   end

   grad = cell(1, numel(schedule.control));
   for k = 1:numel(schedule.control)
      ck = (schedule.step.control == k);
      grad{k} = sum(cell2mat(gradFull(ck)), 2);
   end

   varargout{1} = grad;

   if opt.timeStepGradient
      varargout{2} = gradFull;
   end

   dispif(vb, '************Simulation done: %7.2f seconds ********************\n', toc(timero))
end
%--------------------------------------------------------------------------

function state = loadState(states, inNm, tstep)
   if ~isempty(states)
      % State has been provided, no need to look for it.
      state = states{tstep+1};
      return
   end
   out = load(inNm(tstep));
   fn  = fieldnames(out);
   if numel(fn) ~= 1
      error('Unexpected format for simulation output')
   else
      state = out.(fn{:});
   end
end
