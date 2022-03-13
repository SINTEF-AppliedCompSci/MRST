function [state, varargout] = simulateToSteadyStateADI(G, ...
    rock, fluid, sW0, varargin)
% Run simulation to a steady state solution is found using fully implicit
% adi solver.
%
% SYNOPSIS:
%   state           = simulateToSteadyStateADI(G, rock, fluid, state0, ...
%                         'pn', pv, ...)
%   [state, report] = simulateToSteadyStateADI(...)
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

opt = struct( ...
    'bc',             [],       ... % fixed boundary conditions TODO not tried
    'bcp',            [],       ... % periodic boundary conditions
    'dtInit',         10*day,   ... % initial dt
    'dtIncFac',       2.5,      ... % dt increasing factor
    'dtRedFac',       0.5,      ... % dt reduction factor
    'dtMax',          300*year, ... % maximum dt allowed
    'dtIncTolSat',    300,      ... % tol for absDiffSat to reduce dt
    'absTolPres',     10,       ... % absolute pressure tolerence
    'absTolSat',      1e-4,     ... % absolute saturation tolerence, 1e-4
    'absTolFlux',     1e-8,     ... % absolute flux tolerence, 1e-8
    'absTolPoly',     1e-5,     ... % absolute polymer tolerence, 1e-5
    'nIterMax',       20,       ... % max number of time iterations
    'adiOpt',         [],       ... % passed on to initADISystemBC
    'plotCallback',   [],       ...
    'verbose',        false,    ...
    'detailedOutput', false     ...
    );

opt = merge_options(opt, varargin{:});

% The following MRST modules are required
require ad-blackoil

% Verbose option
vb = opt.verbose;

returnReport = nargout > 1;

%polymer = isfield(fluid, 'ads'); % TEMP SOLUTION
polymer = false; % TEMP SOLUTION

% In the steady-state upscaling, we assume the fluid is incompressible.
% This needs to be specified to the equations.
fluid.isIncomp = true;

% Set up model and initial state.
if polymer
    error('Steady-state polymer model not implemented!')
    model  = TwoPhaseOilWaterModel(G, rock, fluid);
else
    model  = TwoPhaseOilWaterModel_BCP(G, rock, fluid);
end
model = model.validateModel();
state0 = initResSol(G, 50*barsa, [sW0 1-sW0]);
state0.wellSol = initWellSolAD([], model, state0);


%% Simulate some timesteps


solver = NonLinearSolver();

% Set up solution loop
dt            = opt.dtInit;
steadyState   = false;
nIter         = 0;
nIterMax      = opt.nIterMax;
prevState     = state0;

% Allocate vectors
timesteps = zeros(nIterMax, 1);
if polymer
    diffVec = zeros(nIterMax, 4);
else
    diffVec = zeros(nIterMax, 3);
end

% Loop to steady state or max iterations
while ~steadyState && nIter < nIterMax
    
    nIter = nIter + 1;
    
    % Store timestep:
    timesteps(nIter) = dt;
    
    % Solve
    [curState, report] = solver.solveTimestep(prevState, dt, model, ...
        'bc', opt.bc, 'bcp', opt.bcp);
    
    % Check if Newton iterations failed
    if ~report.Converged
        
        % Cut timestep in two and try again
        dt = dt * opt.dtRedFac;
        continue;
        
    end
    
    foundSteadyState = false;
    
    diffnorm = @(curr, prev) norm(curr-prev, inf)*(year/dt);
	absDiffSat  = diffnorm(curState.s, prevState.s);
    
    
    if vb && nIter==1 % Print heading in verbose mode
        fprintf('  Saturation  Pressure     Flux');
        if polymer, fprintf('      Polymer'), end
        fprintf('\n');
        n = 3; if polymer, n = 4; end
        fprintf(repmat('      -    ',1,n));
    end
    
    % Run minimum two iterations
    if nIter > 1
        
        % Compute convergence norms and check if steady state is reached
        absDiffPres = diffnorm(curState.pressure, prevState.pressure);
        absDiffFlux = diffnorm(curState.flux, prevState.flux);
        if polymer
            absDiffPoly = diffnorm(curState.c, prevState.c);
        end
        
        % Make vector of convergence norms
        diffVec(nIter,1:3) = [absDiffSat, absDiffPres, absDiffFlux];
        if polymer
            diffVec(nIter,4) = absDiffPoly;
        end
        
        % Display status if in verbose mode
        if vb
            fprintf('   %1.2e   %1.2e   %1.2e', diffVec(nIter,1:3));
            if polymer
                fprintf('   %1.2e', diffVec(nIter,4));
            end
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
        dispif(vb, '\n');
        dispif(vb, '   Steady state found after %d iterations.\n', nIter);
        
    else
        
        % Steady state is not reached yet. Move to next iteration.
        
        % Check timestep
        if absDiffSat < opt.dtIncTolSat
            % TODO Note: Halvors code 'simulateToSteadyStatePeriodic' also
            % checkes for change in flux.
            if dt < opt.dtMax
                dt = min(opt.dtMax, opt.dtIncFac*dt);
                dispif(vb, '  -> Timestep: %0.3e days', dt/day);
            end
        end
        
        dispif(vb,'\n');
        
        % Initialize next iteration
        prevState = curState;
        
    end
    
end

nIter = min(nIter, nIterMax);

if ~steadyState
    warning('Max number of iterations reached. Steady state not found.');
end

% Output state is last state
state = curState;

% Return metadata if requested
if returnReport
    report.Converged   = steadyState;
    report.Iterations  = nIter;
    report.Timesteps   = timesteps(1:nIter);
    report.Residuals.Saturation  = diffVec(1:nIter, 1);
    report.Residuals.Pressure    = diffVec(1:nIter, 2);
    report.Residuals.Flux        = diffVec(1:nIter, 3);
    if polymer
        report.Residuals.Polymer = diffVec(1:nIter, 4);
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
    report.input      = input;
    
    varargout{1}      = report;
end


end
