function [state, its, convergence] = ...
      solvefiADI(state0, dt, W, G, system, varargin)
%Solve single timestep of automatic differentiation system
%
% SYNOPSIS:
%   state = solvefiADI(state0, dt, W, G, system)
%   state = solvefiADI(state0, dt, W, G, system, ...)
%
% PARAMETERS:
%   state0 - Initial state as defined by initResSol
%
%   dt     - Size of timestep.
%
%   W      - Well configuration. Must have valid sign for all producers and
%            injectors.
%
%   G      - Grid. See grid_structure.
%
%   system - Valid system as defined by initADISystem. This function
%            defines all options for convergence, how to set up equations
%            and parameters for the non-linear solver.
%
%   ...    - Additional parameters passed on to 'system.stepFunction'.
%            OPTIONAL.
%
% RETURNS:
%   state  - Updated state after convergence or max iterations
%            being reached.
%
% SEE ALSO:
%   stepOW, stepBlackOil, stepOWPolymer, runScheduleADI

% COMMENTS:
%
%   The fully implicit AD solvers are written to minimize duplication of
%   coding effort. Primarly, this consists of a splitting of general logic
%   for solving schedules and timesteps and the setup of equations and
%   their solution. In general, a system is created which defines the
%   phases present along with handles to the stepFunction for that specific
%   type of system.
%
%   When solving for a single timestep, the stepFunction is called until
%   the stepFunction reports convergence or stopping in the meta struct it
%   returns.
%
%   The step functions can be completely arbitrary, but for the current
%   solvers the general structure are:
%
%    - solvefiADI calls stepFunction
%    - The stepFunction calls eqsfi<systemtype>, solves a single Newton
%      step to accuracies defined by the system input using some kind of
%      linear solver (Direct solver / some CPR preconditioner
%      configuration) and then updates the state.
%    - solvefiADI checks wells and returns the final state.
%
%   Depending on the simulation, solvefiADI may be exposed directly in a
%   for loop or hidden within a wrapper such as runScheduleADI.

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


    % Solve equations using a general iterative process defined by
    % stepFunction. This is typically Newton iterations.

    step = @(state, meta) ...
       system.stepFunction(state0, state, meta, dt, W, G, ...
                           system, varargin{:});

    meta = struct('converged'   , false, ...
                  'stopped'     , false, ...
                  'res_history', [], ...
                  'linsolver_diverged', false, ...
                  'wellschanged', false, ...
                  'relax'       , system.nonlinear.relaxation, ...
                  'stagnate'    , false, ...
                  'iteration'   , 0);

    state         = state0;
    %state.wellSol = initWellSolLocal(W, state);

    timer = tic;

    while ~ (meta.converged || meta.stopped),
        % Save iteration number in meta info
        meta.iteration = meta.iteration + 1;

        [state, meta] = step(state, meta);

        if meta.stagnate,
            warning('newt:stagnate', 'Non-linear solver stagnated...')
        end
    end
    if(isfield(system,'updateFinal'))
            state=system.updateFinal(state, state0);
    end
    
    if ~isempty(meta.res_history)
       convergence.residuals = meta.res_history(1:meta.iteration, :);
    else
       convergence.residuals = [];
    end
    convergence.converged = meta.converged;

    if meta.stopped
       if meta.linsolver_diverged
          warning('newt:linsolvdiv', ['Linear solver diverged']);
       else
          warning('newt:maxit', ...
                  ['Non-linear solver did not converge, stopped ', ...
                   'by max iterations...']);
       end
    end

    dispif(mrstVerbose, 'Completed %d iterations in %1.2f s\n', ...
           meta.iteration, toc(timer));

    its = meta.iteration;
end
