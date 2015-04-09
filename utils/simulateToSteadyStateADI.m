function [state, varargout] = simulateToSteadyStateADI(state0, G, ...
   rock, fluid, varargin)
% Run simulation to a steady state solution is found using fully implicit 
% adi solver.
% 
% SYNOPSIS:
%   state         = simulateToSteadyStateADI(state0, G, rock, ...
%                       fluid, 'pn', pv, ...)
%   [state, meta] = simulateToSteadyStateADI(...)
% 
% PARAMETERS:
%   state0  - initial state
%   G       - grid structure
%   rock    - rock structure
%   fluid   - fluid structure
% 
% OPTIONAL PARAMETERS
%   'pn'/pv     - List of 'key'/value pairs defining optional parameters.
%                 The supported options are as follows:
%   bc          - fixed boundary conditions
%   Gp          - periodic grid structure
%   bcp         - periodic boundary conditions
%   dtInit      - initial timestep
%   dtIncFac    - timestep multiplyer for increasing timestep
%   dtMax       - maximum timestep allowed
%   dtIncTolSat - tolerence factor on the change in saturation for when 
%                 to increase timestep
%   absTolPres  - absolute pressure tolerence
%   absTolSat   - absolute saturation tolerence
%   absTolFlux  - absolute flux tolerence
%   absTolPoly  - absolute polymer tolerence
%   nIterMax    - maximum number of iterations
% 
% RETURNS:
%   state   - steady state solution
%   meta    - (OPTIONAL) meta data structure
% 

opt = struct( ...
   'bc',           [],       ... % fixed boundary conditions TODO not tried
   'Gp',           [],       ... % periodic grid
   'bcp',          [],       ... % periodic boundary conditions
   'dtInit',       10*day,   ... % initial dt
   'dtIncFac',     2.5,      ... % dt increasing factor
   'dtRedFac',     0.5,      ... % dt reduction factor
   'dtMax',        300*year, ... % maximum dt allowed
   'dtIncTolSat',  300,      ... % tol for absDiffSat to reduce dt
   'absTolPres',   10,       ... % absolute pressure tolerence
   'absTolSat',    1e-4,     ... % absolute saturation tolerence, 1e-4
   'absTolFlux',   1e-8,     ... % absolute flux tolerence, 1e-8
   'absTolPoly',   1e-5,     ... % absolute polymer tolerence, 1e-5
   'nIterMax',     20,       ... % max number of time iterations
   'adiOpt',       [],       ... % passed on to initADISystemBC
   'plotCallback', [] ...
   );

opt = merge_options(opt, varargin{:});

% % Log
% if exist('logf','file') == 2
%    wlog = @(s, varargin) logf(s, varargin{:});
% else
%    wlog = @(s, varargin) 0; % do nothing
% end

% Verbose option
vb = mrstVerbose;

% Output meta data
doReturnMeta = nargout > 1;

% Create dummy injection well
%schedule.control = getDummyWellSchedule();

if isempty(opt.adiOpt)
   opt.adiOpt = {};
end

% Initialize pressure
if ~isempty(opt.bc)
   state0 = initPressure(state0, G, fluid, rock, opt.bc);
end

% Set up ADI system
if isfield(state0, 'c')
   polymer = true;
   adiInput = {'oil','water','polymer'};
   system = initADISystemBC(adiInput, G, rock, fluid, 'bc', opt.bc, ....
      'Gp', opt.Gp, opt.adiOpt{:});
   system.stepFunction = @(state0, state, meta, dt, W, G, system) ...
      stepOWPolymerBC(state0, state, meta, dt, G, W, system, fluid, ...
      'bc', opt.bc, 'bcp', opt.bcp);
%    system.getEquations = @eqsfiOWPolymerExplicitWellsBC;
else
   polymer = false;
   adiInput = {'oil','water'};
   system = initADISystemBC(adiInput, G, rock, fluid, 'bc', opt.bc, ...
      'Gp', opt.Gp, opt.adiOpt{:});
   system.stepFunction = @(state0, state, meta, dt, W, G, system) ...
      stepOWBC(state0, state, meta, dt, G, W, system, fluid, ...
      'bc', opt.bc, 'bcp', opt.bcp);
%    system.getEquations = @eqsfiOWExplicitWellsBC;
end

% % Non-linear convergence tolerance
% systemPolymer.nonlinear.tolMB = 1e-5; % default: 1e-3
% systemPolymer.nonlinear.tolCNV = 1e-9; % default: 1e-7

% Set up solution loop
dt            = opt.dtInit;
steadyState   = false;
nIter         = 1;
nIterMax      = opt.nIterMax;
nAdiSteps     = 1; % number of timesteps in adi code
prevState     = state0;

% Allocate vectors
timesteps = zeros(nIterMax, 1);
if polymer
   diffVec = zeros(nIterMax, 4);
else
   diffVec = zeros(nIterMax, 3);
end

if isempty(opt.Gp)
   Gact = G;
else
   Gact = opt.Gp;
end

% Loop to steady state or max iterations
while ~steadyState && nIter <= nIterMax
   
   % Store timestep:
   timesteps(nIter) = dt;
   
   % Create schedule
   schedule.step.val = dt.*ones(nAdiSteps, 1);
   
   % Run simulation
   [~, outputStates, ~, conv] = ...
      runScheduleADIBC(prevState, Gact, rock, system, schedule);
   curState = outputStates{end};
   
   % Call plot callback
   if ~isempty(opt.plotCallback)
      opt.plotCallback(G, curState, prevState);
   end
   
   % Check if Newton iterations failed
   if ~conv.converged %(newtonIts == system.nonlinear.maxIterations)
      
      % Cut timestep in two and try again
      dt = dt * opt.dtRedFac;
      
   else
      
      foundSteadyState = false;
      
      
      absDiffSat = norm(curState.s - prevState.s, inf) * (year/dt);
      
      if nIter > 1
         
         % Compute changes and check if steady state is reached

         absDiffPres = norm(curState.pressure - prevState.pressure, inf) * ...
            (year/dt);
         absDiffFlux = norm(curState.flux - prevState.flux, inf) * (year/dt); 
         if polymer
            absDiffPoly = norm(curState.c - prevState.c, inf) * (year/dt);
         end

         diffVec(nIter,1:3) = [absDiffSat, absDiffPres, absDiffFlux];
         if polymer
            diffVec(nIter,4) = absDiffPoly;
         end

         dispif(vb, ['absDiffSat=%1.2e, absDiffPres=%1.2e, ' ...
            'absDiffFlux=%1.2e'], diffVec(nIter,1:3));
         if polymer
            dispif(vb, 'absDiffPoly=%1.2e\n', diffVec(nIter,4));
         else
            dispif(vb, '\n');
         end


         if absDiffSat < opt.absTolSat && absDiffPres < opt.absTolPres ...
               && absDiffFlux < opt.absTolFlux
            if polymer
               if absDiffPoly < opt.absTolPoly
                  foundSteadyState = true;
               end
            else
               foundSteadyState = true;
            end
         end
         
      end
      
      if foundSteadyState
         
         steadyState = true;
         dispif(vb, 'Steady state found after %d iterations.\n', nIter);
         
      else
         
         % Steady state is not reached yet. Move to next iteration.
         
         % Check timestep
         if absDiffSat < opt.dtIncTolSat
            % TODO Note: Halvors code 'simulateToSteadyStatePeriodic' also
            % checkes for change in flux.
            if dt < opt.dtMax
               dt = min(opt.dtMax, opt.dtIncFac*dt);
               dispif(vb, 'Timestep changed to %0.3e days\n', dt/day);
            end
         end
         
         % Initialize next iteration
         prevState = curState;
         nIter = nIter + 1;
         
      end
      
   end
   
end

nIter = min(nIter, nIterMax);

if ~steadyState
   warning('Max number of iterations reached. Steady state not found.');
end

% Output state is last state
state = curState;

% Return metadata if requested
if doReturnMeta
   meta.steadyState = steadyState;
   meta.iterations  = nIter;
   meta.timesteps   = timesteps(1:nIter);
   meta.absDiffSat  = diffVec(1:nIter, 1);
   meta.absDiffPres = diffVec(1:nIter, 2);
   meta.absDiffFlux = diffVec(1:nIter, 3);
   if polymer
      meta.absDiffPoly = diffVec(1:nIter, 4);
   end
   
   input.dtInit      = opt.dtInit;
   input.dtIncFac    = opt.dtIncFac;
   input.dtRedFac    = opt.dtRedFac;
   input.dtMax       = opt.dtMax;
   input.dtIncTolSat = opt.dtIncTolSat;
   input.absTolPres  = opt.absTolPres;
   input.absTolSat   = opt.absTolSat;
   input.absTolFlux  = opt.absTolFlux;
   input.absTolPoly  = opt.absTolPoly;
   input.nIterMax    = opt.nIterMax;
   input.adiOpt      = opt.adiOpt;
   meta.input        = input;
   
   varargout{1}      = meta;
end


end


% function c = getDummyWellSchedule()
% % Creates and returns a dummy well schedule. The schedule contains a water
% % injection well with zero injection rate.
% 
% f = {'WELSPECS', 'COMPDAT', 'WCONHIST', 'WCONINJ', 'WCONINJE', ...
%     'WCONINJH', 'WCONPROD', 'GCONPROD', 'GCONINJE', 'GRUPTREE', ...
%     'WGRUPCON', 'WPOLYMER', 'DRSDT'};
% c = cell2struct(cell(1,length(f)), f, 2);
% 
% c.WELSPECS = ...
%    {'DUMMY' 'G' 1 1 0 'WATER' 0 'STD' 'SHUT' 'NO' 0 'SEG' 0};
% c.COMPDAT  = ...
%    {'DUMMY' 1 1 1 1 'OPEN' -1 0 1e-10 0 0 'Default' 'Z' -1};
% c.WCONINJE = ...
%    {'DUMMY' 'WATER' 'OPEN' 'RATE' 0 Inf NaN Inf 0 0};
% c.WPOLYMER = ...
%    {'DUMMY' 0 0 'Default' 'Default'};
% c.DRSDT    = {Inf  'ALL'};
% 
% end






function state0 = initPressure(state0, G, fluid, rock, bc)
% Get initial pressure guess from incompTPFA

% Fix fluid
fluid.properties = @(state) fluidPropertiesField(fluid, state);
fluid.saturation = @(state) state.s;

if ~isfield(fluid, 'krO')
   fluid.krO = fluid.krOW;
end
fluid.relperm = @(s,state) [fluid.krW(state.s(:,1)) ...
                            fluid.krO(state.s(:,2)) ];

T  = computeTrans(G, rock, 'Verbose', true);
   
% Compute pressures
resSol = initResSol(G, state0.pressure, state0.s);
resSol = incompTPFA(resSol, G, T, fluid, 'bc', bc);

% Create initial state
state0.pressure = resSol.pressure;


end




function [mu, rho] = fluidPropertiesField(fluid, state)

pW = state.pressure - fluid.pcOW(state.s(:,1));
bO = fluid.bO(state.pressure);
bWmean = mean(fluid.bW(pW));
bOmean = mean(bO);

if isfield(fluid, 'BOxmuO')
   mu  = [fluid.muW(pW), fluid.BOxmuO(state.pressure).*bO];
else
   mu  = [fluid.muW(pW), fluid.muO(state.pressure)];
end
rho = [fluid.rhoWS*bWmean,  fluid.rhoOS*bOmean];

end


