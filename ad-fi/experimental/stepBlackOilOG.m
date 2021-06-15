function [state, meta] = stepBlackOilOG(state0, state, meta, dt, G, W, system, fluid, varargin)
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

% Get the equations given current states
%[eqs, history, explTrms] = eqsfiBlackOilExplicitWellsOG(state0, state, dt, G, W, s, fluid, 'history', meta.history);
[eqs, history, explTrms] =system.getEquations(state0, state, dt, G, W, s, fluid,'history', meta.history);
if system.nonlinear.cpr && isempty(system.podbasis)
    % Solve the equations using CPR preconditioner
    [dx, gmresits, gmresflag] = cprGeneric(eqs, system,...
                                'ellipSolve', system.nonlinear.cprEllipticSolver,...
                                'cprType',    system.nonlinear.cprType,...
                                'relTol',     system.nonlinear.cprRelTol);
else
    dx = SolveEqsADI(eqs, system.podbasis);
    gmresits = [0 0];
    gmresflag = 0;
end

% dx_old = meta.dx;
% meta.dx = dx;

[meta, residuals] = getResiduals(meta, eqs, system, gmresflag);

%

if system.nonlinear.linesearch
    getEqs = @(state) eqsfiBlackOilExplicitWellsOG(state0, state, dt, G, W, s, fluid, 'history', meta.history, 'resOnly', true);
    upState = @(dx, explTrms) updateState(W, state, dx, fluid, explTrms);
    [state, dx, searchfail] = linesearchADI(state, dx, system, getEqs, upState, true);
else
   searchfail = false;
end

% Update reservoir conditions once a delta has been found.
if ~searchfail
    [dx, meta] = stabilizeNewton(dx, meta, system);
    % If the line search failed, uncritically accept the first step and
    % pray the other measures (relaxation / dampening) handle the error.
    [state, nInc] = updateState(W, state, dx, fluid, explTrms);
end

% Check convergence
[converged CNV MB] = getConvergence(state, eqs, fluid, system, dt);

meta.converged = converged;
meta.stopped = meta.iteration == system.nonlinear.maxIterations && ~converged;
meta.history = history;

if opt.Verbose
    eqnnames = {'Oil', 'Gas', 'Dis.Gas', 'qOs', 'qGs', 'pBHP'};
    printResidual(residuals, gmresits, eqnnames, meta.iteration, CNV, MB);
end

end

%--------------------------------------------------------------------------
% This is a 'cleaned' version of this function.  The original version, with
% varous comment-outs and original indenting, is preserved below
% (commented-out) for future reference.
function [state, nInc] = updateState(W, state, dx, f, explTrms)

    maxSatCh = .2;
    dp       = dx{1};
    dsg      = dx{2};
    drs      = dx{3};
    inx      = abs(dsg) > maxSatCh;
    dsg(inx) = sign(dsg(inx))*maxSatCh;
    cap      = @(x) max(x, 0);

    nInc = max( [norm(dp,'inf')/norm(state.pressure, 'inf'), ...
                 norm(dsg,'inf'), ...
                 norm(drs,'inf')/norm(state.rs, 'inf')] );

    state.pressure = cap(state.pressure + dp);
    % state.pressure(state.pressure<200*barsa)=201*barsa; % @@ This line really
                                                          % doesn't seem right,
                                                          % so I commented it out
                                                          % for now.
    state.rs = cap(state.rs + drs);
    sg       = state.s(:,2) + dsg;

    sg=max(0,sg);
    sg=min(1,sg);
    state.s  = [1-sg, sg];

    %--------------
    % NB: eqnnames = {'Oil', 'Gas', 'Dis.Gas', 'qOs', 'qGs', 'pBHP'};
    dqOs  = dx{4};
    dqGs  = dx{5};
    dpBHP = dx{6}; %maybe put some limit on this?

    for w = 1:numel(state.wellSol)
        state.wellSol(w).bhp = state.wellSol(w).bhp + dpBHP(w);
        state.wellSol(w).qOs      = state.wellSol(w).qOs + dqOs(w);
        state.wellSol(w).qGs      = state.wellSol(w).qGs + dqGs(w);
    end

    % Add explicit terms
    if isfield(explTrms, 'wellFlux')
        winx = [0; cumsum(arrayfun(@(x)numel(x.cells), W))];
        for wnr = 1:numel(W)
            state.wellSol(wnr).flux = explTrms.wellFlux((winx(wnr)+1) : winx(wnr+1),:);
        end
    end
    if isfield(explTrms, 'conPr')
        winx = [0; cumsum(arrayfun(@(x)numel(x.cells), W))];
        for wnr = 1:numel(W)
            state.wellSol(wnr).cpr = explTrms.conPr((winx(wnr)+1) : winx(wnr+1),:);
        end
    end
end


% %--------------------------------------------------------------------------
% function [state, nInc] = updateState(W, state, dx, f, explTrms)
% %maxSatStep = .25;


% maxSatCh = .2;
% %maxRsCh  = inf;

% % d = findDampening(state, W, dx{8});
% % dx = cellfun(@(x) x*d, dx, 'UniformOutput', false);

% dp  = dx{1};

% %dsw = dx{2};
% %inx = abs(dsw) > maxSatCh;
% %dsw(inx) = sign(dsw(inx))*maxSatCh;

% dsg = dx{2};
% inx = abs(dsg) > maxSatCh;
% dsg(inx) = sign(dsg(inx))*maxSatCh;

% drs = dx{3};

% cap = @(x) max(x, 0);

% nInc = max( [norm(dp,'inf')/norm(state.pressure, 'inf'), ...
%              norm(dsg,'inf'), ...
%              norm(drs,'inf')/norm(state.rs, 'inf')] );


% %dp(dp>30*barsa)=30*barsa;
% %dp(dp<30*barsa)=-30*barsa;

% state.pressure = cap(state.pressure + dp);
% %state.pressure = (state.pressure + dp);
% state.pressure(state.pressure<200*barsa)=201*barsa;
% rs = cap(state.rs + drs);


% %{
% %This has serios troubles with time dependent resolution
% state0 = state;
% epsilon = sqrt(eps);
% appG = true;
% if appG
%     sg0 = state0.s(:,2);
%     isSat0  = sg0 > 0;
%     keepSat = sg0 > epsilon*(1+eps);

%     rs0    = state0.rs;
%     rsSat0 = f.rsSat(state0.pressure);
%     keepUSat = rs0 < rsSat0*(1-epsilon)*(1-eps);

%     sg  = sg0 + dsg;
%     rsSat = f.rsSat(state.pressure);
%     isSat = or(sg > 0, rs > rsSat);

%     setToEps  = and(isSat<isSat0, keepSat);
%     setToZero = and(isSat<isSat0, ~keepSat);

%     setToMaxEps = and(isSat>isSat0, keepUSat);
%     setToMax    = and(isSat>isSat0, ~keepUSat);

%     sg(setToEps)  = epsilon;
%     rs(setToEps)  = rsSat(setToEps);

%     sg(setToZero) = 0;
%     rs(setToZero) = rsSat(setToZero)*(1-epsilon);

%     sg(setToMaxEps) = 0;
%     rs(setToMaxEps) = rsSat(setToMaxEps)*(1-epsilon);

%     sg(setToMax) = epsilon;
%     rs(setToMax) = rsSat(setToMax);

% else
% %}
%     sg = state.s(:,2) + dsg;
% %end


% state.rs = rs;
% if(isfield(f,'cutValues') & false)
%     state.s  = [1-sg, sg];
%     state=f.cutValues(state);
% else
%     sg=max(0,sg);
%     sg=min(1,sg);
%     state.s  = [1-sg, sg];
% end


% %--------------
% %eqnnames = {'Oil', 'Gas', 'Dis.Gas', 'qOs', 'qGs', 'pBHP'};
% dpBHP = dx{6}; %maybe put some limit on this?
% dqOs  = dx{4};
% dqGs  = dx{5};



% for w = 1:numel(state.wellSol)
%     state.wellSol(w).bhp = state.wellSol(w).bhp + dpBHP(w);
%     state.wellSol(w).qOs      = state.wellSol(w).qOs + dqOs(w);
%     state.wellSol(w).qGs      = state.wellSol(w).qGs + dqGs(w);
% end

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

% function wsol = updateWellSol(w, wsol, dp, dqO, dqG)

%     dmod = 1;
%     pressure = wsol.pressure + dp;
%     if ((w.sign ==  1) && pressure > w.bhpLimit) ||...
%        ((w.sign == -1) && pressure < w.bhpLimit)
%         overstep = w.sign*(pressure - w.bhpLimit);
%         dmod = 1 - overstep/dp;
%         pressure = w.bhpLimit + w.sign*1*barsa;
%     end
%     wsol.pressure = pressure;
%     % Ensure that we take a step equal to
%     wsol.qOs = wsol.qOs + dmod*dqO;
%     wsol.qGs = wsol.qGs + dmod*dqG;
% end

% function d = findDampening(state, W, dp)

%     d = 1;
%     for i = 1:numel(W)
%         pressure = state.wellSol(i).bhp + dp(i);
%         w = W(i);
%         if ~strcmp(w.type, 'bhp')
%             assert(((w.sign ==  1) && state.wellSol(i).bhp - sqrt(eps) <= w.bhpLimit) || ((w.sign == -1) && state.wellSol(i).bhp + sqrt(eps) >= w.bhpLimit))
%         end
%         if ((w.sign ==  1) && pressure > w.bhpLimit) ||...
%            ((w.sign == -1) && pressure < w.bhpLimit)
%             overstep = w.sign*(pressure - w.bhpLimit);

%             dmod = overstep/dp(i);
% %             if( abs(dmod) == 1)
% %                 dmod = .1;
% %             end
%             d = min(1 - dmod, d);
%             if d < sqrt(eps)
%                 d = 1e-3;
%             end
%         end
%     end
% end
