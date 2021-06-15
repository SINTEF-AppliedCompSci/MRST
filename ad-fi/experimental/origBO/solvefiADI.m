function [state, its] = solvefiADI(state0, dt, W, G, system)
% Solve a automatic differentiation system for a single timestep.
%
% SYNOPSIS:
%   state = solvefiADI(state, dt, W, G);
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
% RETURNS:
%   state  - Updated state after convergence or max iterations being reached.
%
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
%    step to accuracies defined by the system input using some kind of
%    linear solver (Direct solver / some CPR preconditioner configuration)
%    and then updates the state.
%    - solvefiADI checks wells and returns the final state.
%
%   Depending on the simulation, solvefiADI may be exposed directly in a
%   for loop or hidden within a wrapper such as runScheduleADI.
%
% SEE ALSO:
%   stepOW, stepBlackOil, stepOWPolymer, runScheduleADI

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

    meta.converged = false;
    meta.stopped = false;
    meta.wellschanged = false;
    meta.relax = system.nonlinear.relaxation;

    state0.wellSol = initWellSolLocal(W, state0);

    timer = tic;


    meta.iteration = 0;
    while ~meta.converged && ~meta.stopped
        if meta.iteration == 0
            state = state0;
        end
        % Save iteration number in meta info
        meta.iteration = meta.iteration + 1;
        [state, meta] = system.stepFunction(state0, state, meta, dt, W, G, system);
        if system.nonlinear.changeWells || system.nonlinear.bhpcontrols
            [W, meta, state, changed] = checkWellLimits(W, meta, state, system.nonlinear.bhpcontrols);
            if changed
                meta.iteration = meta.iteration - 1;
                meta.converged = false;
            end
        end

    end
    if(isfield(state,'smax'))
        assert(isfield(state,'smin'));
        if(isfield(system,'updateFinal'))
            state=system.updateFinal(state, state0);
        end
    end
    if meta.stopped
        warning('newt:maxit', 'Non-linear solver did not converge, stopped by max iterations...');
    end
    dispif(mrstVerbose, 'Completed %d iterations in %1.2f s\n', meta.iteration, toc(timer));
    its = meta.iteration;
end

function [W, meta, state, changed] = checkWellLimits(W, meta, state, useBHP)
    changed = false;
%     if meta.iteration == 1
%         return
%     end
    for i = 1 : numel(W)
        w = W(i);
        ws = state.wellSol(i);
        if strcmpi(w.type, 'bhp')
            continue
        end
        if ((w.sign ==  1) && ws.pressure > w.bhpLimit) || ...
           ((w.sign == -1) && ws.pressure < w.bhpLimit) || useBHP
            W(i).type = 'bhp';
            W(i).val = w.bhpLimit;
            state.wellSol(i).bhp = w.bhpLimit;
%             fprintf('\t%d Limit %2.2g, was %2.2g\n', w.sign, w.bhpLimit, ws.pressure);
            fprintf('Changing well %s!\n', W(i).name);
            changed = true;
        end
    end
end

