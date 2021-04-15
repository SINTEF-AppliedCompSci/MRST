function [state, meta] = stepOG(state0, state, meta, dt, G, W, system, fluid, varargin)
% Do a single step of a nonlinear solve for a Oil-Water system
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


opt = struct('Verbose', mrstVerbose,'bc',[],'ecl_conv',false);
opt = merge_options(opt, varargin{:});
s = system.s;

% if ~isempty(system.podbasis)
%     solve = @(eqs) SolveEqsADIPOD(eqs, opt.podbasis);
% else
%eqs = eqsfiOGExplicitWells(state0, state, dt, G, W, s, fluid);
eqs =system.getEquations(state0, state, dt, G, W, s, fluid,'bc',opt.bc);
if system.nonlinear.cpr && isempty(system.podbasis)
    if 0
    [dx, gmresits, gmresflag] =cprTwoPhaseExplicitWells(eqs,...
                    'ellipSolve', system.nonlinear.cprEllipticSolver,...
                    'cprType',    system.nonlinear.cprType,...
                    'relTol',     system.nonlinear.cprRelTol);
    else
    [dx, gmresits, gmresflag] = cprGeneric(eqs, system,...
                    'ellipSolve', system.nonlinear.cprEllipticSolver,...
                    'cprType',    system.nonlinear.cprType,...
                    'relTol',     system.nonlinear.cprRelTol);
    end
else
    dx = SolveEqsADI(eqs, system.podbasis);
    gmresits = [0 0];
    gmresflag = 0;
end



% dx = solve(eqs);

% [state, nInc] = updateState(state, dx);
[meta, residuals] = getResiduals(meta, eqs, system, gmresflag);

searchfail = true;
if system.nonlinear.linesearch
    getEqs = @(state) system.getEquations(state0, state, dt, G, W, s, fluid, 'resOnly', true,'bc',opt.bc);
    upState = @(dx) updateState(state, dx);
    [state, dx, searchfail] = linesearchADI(state, dx, system, getEqs, upState, false);
end

% Update reservoir conditions once a delta has been found.
if searchfail
    [dx, meta] = stabilizeNewton(dx, meta, system);
    % If the line search failed, uncritically accept the first step and
    % pray the other measures (relaxation / dampening) handle the error.
    [state, nInc] = updateState(state, dx);
end

[converged CNV MB] = getConvergence(state, eqs, fluid, system, dt);
meta.converged = all(residuals< system.nonlinear.tol);
if(~opt.ecl_conv)
    converged=meta.converged;
end


%meta.converged = converged;
meta.stopped = meta.iteration == system.nonlinear.maxIterations && ~converged;

if opt.Verbose
    residuals = cellfun(@(x) norm(x.val, 'inf'), eqs);
    eqnnames = {'Oil', 'Gas',  'qOs', 'qGs', 'pBHP'};
    printResidual(residuals, [], eqnnames, meta.iteration, CNV, MB);
end
end

%--------------------------------------------------------------------------

function [state, nInc] = updateState(state, dx)
maxSatStep = .2;

% if ~isempty(phi)
%     for i = 1:numel(dx)
%         dx{i} = phi.basis{i}*dx{i};
%     end
% end

dp = dx{1};
ds = dx{2};
nInc = max( norm(dp,'inf')/norm(state.pressure, 'inf'), ...
            norm(ds,'inf')/norm(state.s(:,1), 'inf') );

maxch = norm(ds, 'inf');
step = min(1, maxSatStep./maxch);


state.pressure = state.pressure + step*dp;
sg = state.s(:,2) + step*ds;
% Cap values
sg = min(sg, 1); sg = max(sg, 0);

state.s = [1-sg, sg];

dqGs  = step*dx{3};
dqOs  = step*dx{4};
dpBHP = step*dx{5};

for w = 1:numel(state.wellSol)
    state.wellSol(w).bhp = state.wellSol(w).bhp + dpBHP(w);
    state.wellSol(w).qGs = state.wellSol(w).qGs + dqGs(w);
    state.wellSol(w).qOs = state.wellSol(w).qOs + dqOs(w);
end
end
