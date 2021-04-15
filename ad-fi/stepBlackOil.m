function [state, meta] = stepBlackOil(state0, state, meta, dt, G, W, system, fluid, varargin)
% Do a single step of a nonlinear solve for a Black-Oil system (with gas dissolution)
% This function should in general not be called directly and is as such not
% documented with regards to input/output: See solvefiADI for an
% explanation of how the ad-fi solvers are implemented.

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
   [eqs, state, history] = eqsfiBlackOil(state0, state, dt, G, W, system, fluid, ...
       'stepOptions', system.stepOptions, 'history', meta.history, ...
       'iteration', meta.iteration);


   [meta, residuals] = getResiduals(meta, eqs, system, meta.gmresflag);

   % Check convergence
   [converged, CNV, MB] = getConvergence(state, eqs, fluid, system, dt);
   %  add well convergence: qWs, qOs, qGs order of flow
   wellConverged = and(all(residuals(5:7)<1/day), residuals(8)<1*barsa );
   if ~system.nonlinear.use_ecltol
       % Override ecl convergence, but keep CNV/MB numbers for printout
       converged = all(residuals <= system.nonlinear.tol);
   end
   if meta.iteration == 1, converged = false;end
   meta.converged = converged&&wellConverged;

   if opt.Verbose
      eqnnames = {'Oil', 'Water', 'Gas', 'Dis.Gas', 'qWs', 'qOs', 'qGs', 'control'};
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
      getEqs = @(state) eqsfiBlackOil(state0, state, dt, G, W, system, fluid, 'stepOptions', stepOptions, 'history', meta.history, 'resOnly', true);
      upState = @(dx, explTrms) updateState(W, state, dx, fluid, system);
      [state, dx, searchfail] = linesearchADI(state, dx, system, getEqs, upState, true);
   end

   % Update reservoir conditions once a delta has been found.
   if searchfail
      dispif(mrstVerbose && system.nonlinear.linesearch, 'Linesearch failed!\n')
      [dx, meta] = stabilizeNewton(dx, meta, system);
      % If the line search failed, uncritically accept the first step and
      % pray the other measures (relaxation / dampening) handle the error.
      [state, nInc] = updateStateBO(W, state, dx, fluid, system);
   end




   meta.stopped = meta.iteration == system.nonlinear.maxIterations && ~converged;
   meta.history = history;
   meta.gmresflag = gmresflag;



end

%{
% %--------------------------------------------------------------------------
% function [state, nInc] = updateState(W, state, dx, f, explTrms)
% %maxSatStep = .25;
% explTrms = [];
% appW = false;
% appG = true;
%
% maxSatCh = .2;
% %maxRsCh  = inf;
%
% % d = findDampening(state, W, dx{8});
% % dx = cellfun(@(x) x*d, dx, 'UniformOutput', false);
%
% dp  = dx{1};
% dp = sign(dp).*min(abs(dp), 20*barsa);
%
% dsw = dx{2};
% inx = abs(dsw) > maxSatCh;
% dsw(inx) = sign(dsw(inx))*maxSatCh;
%
% dsg = dx{3};
% inx = abs(dsg) > maxSatCh;
% dsg(inx) = sign(dsg(inx))*maxSatCh;
%
% drs = dx{4};
%
% cap = @(x) max(x, 0);
%
% nInc = max( [norm(dp,'inf')/norm(state.pressure, 'inf'), ...
%              norm(dsw,'inf'), ...
%              norm(dsg,'inf'), ...
%              norm(drs,'inf')/norm(state.rs, 'inf')] );
%
% state0 = state;
% state.pressure = cap(state.pressure + dp);
%
% rs = cap(state.rs + drs);
%
% epsilon = sqrt(eps);
% if appW
%     warning('Not implemented')
% else
%     sw = state.s(:,1) + dsw;
% end
%
% if appG
%     sg0 = state0.s(:,3);
%     isSat0  = sg0 > 0;
%     keepSat = sg0 > epsilon*(1+eps);
%
%     rs0    = state0.rs;
%     rsSat0 = f.rsSat(state0.pressure);
%     keepUSat = rs0 < rsSat0*(1-epsilon)*(1-eps);
%
%     sg  = sg0 + dsg;
%     rsSat = f.rsSat(state.pressure);
%     isSat = or(sg > 0, rs > rsSat);
%
%     %setToEps  = and(isSat<isSat0, keepSat);
%     setToEps  = and(sg<0, keepSat);
%     setToZero = and(and(sg<0,sg0>0), ~keepSat);
%
%     setToMaxEps = and(rs>rsSat, keepUSat);
%     setToMax    = and(isSat>isSat0, ~keepUSat);
%
%     sg(setToEps)  = epsilon;
%     rs(setToEps)  = rsSat(setToEps);
%
%     sg(setToZero) = 0;
%     rs(setToZero) = rsSat(setToZero)*(1-epsilon);
%
%     sg(setToMaxEps) = 0;
%     rs(setToMaxEps) = rsSat(setToMaxEps)*(1-epsilon);
%
%     sg(setToMax) = epsilon;
%     rs(setToMax) = rsSat(setToMax);
%     rs(isSat)    = rsSat(isSat);%nei!!!
%
% else
%     sg = state.s(:,3) + dsg;
% end
%
%
% state.rs = rs;
% state.s  = [sw, 1-sw-sg, sg];
%
% %--------------
% %
% % dpBHP = dx{8}; %dpBHP = sign(dpBHP).*min(abs(dpBHP), 25*barsa);%maybe put some limit on this?
% % dqWs  = dx{5};
% % dqOs  = dx{6};
% % dqGs  = dx{7};
% %
% %
% %
% % for w = 1:numel(state.wellSol)
% %     if 0
% %     state.wellSol(w) = updateWellSol(W(w), state.wellSol(w),...
% %                                     dpBHP(w), ...
% %                                     dqWs(w), ...
% %                                     dqOs(w), ...
% %                                     dqGs(w));
% %     else
% %     state.wellSol(w).bhp = state.wellSol(w).bhp + dpBHP(w);
% %     state.wellSol(w).qWs      = state.wellSol(w).qWs + dqWs(w);
% %     state.wellSol(w).qOs      = state.wellSol(w).qOs + dqOs(w);
% %     state.wellSol(w).qGs      = state.wellSol(w).qGs + dqGs(w);
% %     end
% %end
%
% % Add explicit terms
% if isfield(explTrms, 'wellFlux')
%     cc = 0;
%     winx = [0; cumsum(arrayfun(@(x)numel(x.cells), W))];
%     for wnr = 1:numel(W)
%         state.wellSol(wnr).flux = explTrms.wellFlux((winx(wnr)+1) : winx(wnr+1),:);
%         if any(state.wellSol(wnr).flux.*W(wnr).sign < 0)
%             cc = cc +1;
%         end
%     end
% %    if cc > 0
% %        warning(['Crossflow appeared in ', num2str(cc), ' wells']);
% %    end
% end
% if isfield(explTrms, 'conPr')
%     winx = [0; cumsum(arrayfun(@(x)numel(x.cells), W))];
%     for wnr = 1:numel(W)
%         state.wellSol(wnr).cpr = explTrms.conPr((winx(wnr)+1) : winx(wnr+1),:);
%     end
% end
% end
%}

%{
function wsol = updateWellSol(w, wsol, dp, dqW, dqO, dqG)

   dmod = 1;
   pressure = wsol.pressure + dp;
   if ((w.sign ==  1) && pressure > w.bhpLimit) ||...
                ((w.sign == -1) && pressure < w.bhpLimit)
      overstep = w.sign*(pressure - w.bhpLimit);
      dmod = 1 - overstep/dp;
      pressure = w.bhpLimit + w.sign*1*barsa;
   end
   wsol.pressure = pressure;
   % Ensure that we take a step equal to
   wsol.qWs = wsol.qWs + dmod*dqW;
   wsol.qOs = wsol.qOs + dmod*dqO;
   wsol.qGs = wsol.qGs + dmod*dqG;
end
%}

%{
function d = findDampening(state, W, dp)

   d = 1;
   for i = 1:numel(W)
      pressure = state.wellSol(i).bhp + dp(i);
      w = W(i);
      if ~strcmp(w.type, 'bhp')
         assert(((w.sign ==  1) && state.wellSol(i).bhp - sqrt(eps) <= w.bhpLimit) || ((w.sign == -1) && state.wellSol(i).bhp + sqrt(eps) >= w.bhpLimit))
      end
      if ((w.sign ==  1) && pressure > w.bhpLimit) ||...
                   ((w.sign == -1) && pressure < w.bhpLimit)
         overstep = w.sign*(pressure - w.bhpLimit);

         dmod = overstep/dp(i);
         %             if( abs(dmod) == 1)
         %                 dmod = .1;
         %             end
         d = min(1 - dmod, d);
         if d < sqrt(eps)
            d = 1e-3;
         end
      end
   end
end
%}
