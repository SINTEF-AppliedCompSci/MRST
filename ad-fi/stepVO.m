function [state, meta] = stepVO(state0, state, meta, dt, G, W, system, fluid, varargin)
% Do a single step of a nonlinear solve for a Black-Oil system (with gas dissolution)
% This function should in general not be called directly and is as such not
% documented with regards to input/output: See solvefiADI for an
% explanation of how the ad-fi solvers are implemented.

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


   opt = struct('Verbose', mrstVerbose);
   opt = merge_options(opt, varargin{:});

   s = system.s;

   % If we are at the first iteration, some metafields must be created.
   if ~isfield(meta, 'history')
      meta.history = [];
   end

   if ~isfield(meta, 'gmresflag')
      meta.gmresflag = 0;
      meta.gmresits = [0 0];
   end

   % Get the equations given current states
   [eqs, state] = eqsfiVO(state0, state, dt, G, W, s, fluid, system, ...
                'stepOptions', system.stepOptions, ...
                'iteration', meta.iteration);

   [meta, residuals] = getResiduals(meta, eqs, system, meta.gmresflag);

   % Check convergence
   [converged, CNV, MB] = getConvergence(state, eqs, fluid, system, dt);
   %  add well convergence: qWs, qOs, qGs order of flow
   wellConverged = and(all(residuals(4:6)<1/day), residuals(7)<1*barsa );
   if ~system.nonlinear.use_ecltol
       % Override ecl convergence, but keep CNV/MB numbers for printout
       converged = all(residuals <= system.nonlinear.tol);
   end
   if meta.iteration == 1, converged = false;end
   meta.converged = converged&&wellConverged;

   if opt.Verbose
      eqnnames = {'Oil', 'Water', 'Gas', 'qWs', 'qOs', 'qGs', 'control'};
      printResidual(residuals, meta.gmresits, eqnnames, meta.iteration, CNV, MB);
   end

   if meta.converged
       return
   end


   if system.nonlinear.cpr && isempty(system.podbasis)
      % Solve the equations using CPR preconditioner
      %p  = state.pressure; rs = state.rs; rss = fluid.rsSat(p);
      %bW = fluid.bW(p); bO = fluid.bO(p, rs, rs >= rss); bG = fluid.bG(p);
      p  = mean(state.pressure); rs = fluid.rsSat(p);
      bW = fluid.bW(p); bO = fluid.bO(p, rs, true); bG = fluid.bG(p);
      sc = [1./bO, 1./bW, 1./bG];
      %sc = [1./bW, 1./bO-rs./bG, 1./bG];
      [dx, gmresits, gmresflag] = ...
         cprGeneric(eqs, system,                                      ...
                    'ellipSolve', system.nonlinear.cprEllipticSolver, ...
                    'cprType',    system.nonlinear.cprType,           ...
                    'relTol',     system.nonlinear.cprRelTol,         ...
                    'eqScale',       sc);
   else
      dx = SolveEqsADI(eqs, system.podbasis);
      gmresits = [0 0];
      gmresflag = 0;
   end
   meta.gmresits = gmresits;

   %
   searchfail = true;
   if system.nonlinear.linesearch
      stepOptions = system.stepOptions;
      stepOptions.solveWellEqs = false;
      getEqs = @(state) eqsfiBlackOil(state0, state, dt, G, W, s, fluid, 'stepOptions', stepOptions, 'history', meta.history, 'resOnly', true);
      upState = @(dx, explTrms) updateState(W, state, dx, fluid, system);
      [state, dx, searchfail] = linesearchADI(state, dx, system, getEqs, upState, true);
   end

   % Update reservoir conditions once a delta has been found.
   if searchfail
      dispif(mrstVerbose && system.nonlinear.linesearch, 'Linesearch failed!\n')
      [dx, meta] = stabilizeNewton(dx, meta, system);
      % If the line search failed, uncritically accept the first step and
      % pray the other measures (relaxation / dampening) handle the error.
      state = updateStateVO(W, state, dx, fluid, system);
   end




   meta.stopped = meta.iteration == system.nonlinear.maxIterations && ~converged;
%   meta.history = history;
   meta.gmresflag = gmresflag;
end


