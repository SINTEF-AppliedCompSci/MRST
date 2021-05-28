function [state, meta] = stepOWT(state0, state, meta, dt, G, W, system, fluid, varargin)
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


opt = struct('Verbose', mrstVerbose, 'temperature', false, 'minerals', false);
opt = merge_options(opt, varargin{:});
s = system.s;

% if ~isempty(system.podbasis)
%     solve = @(eqs) SolveEqsADIPOD(eqs, opt.podbasis);
% else
[eqs, state] = system.getEquations(state0, state, dt, G, W, s, fluid, ...
    'temperature', opt.temperature, ...
    'minerals', opt.minerals, ...
    'iteration', meta.iteration);

if system.nonlinear.cpr && isempty(system.podbasis)
    p  = mean(state0.pressure);
    bW = fluid.bW(p); bO = fluid.bO(p);
    sc = [1./bO, 1./bW];

    vargs = { 'ellipSolve', system.nonlinear.cprEllipticSolver, ...
              'cprType'   , system.nonlinear.cprType          , ...
              'relTol'    , system.nonlinear.cprRelTol        , ...
              'eqScale'   , sc};

    [dx, gmresits, solver_diverged] = cprGeneric(eqs, system, vargs{:});

else

    dx = SolveEqsADI(eqs, system.podbasis);
    solver_diverged = false;
    gmresits = [0 0];
end



% dx = solve(eqs);

% [state, nInc] = updateState(state, dx);


searchfail = true;
if system.nonlinear.linesearch
    getEqs = @(state) system.getEquations(state0, state, dt, G, W, s, fluid, 'resOnly', true,...
        'temperature', opt.temperature,...
        'minerals',opt.minerals);
    upState = @(dx) updateState(state, dx, opt);
    [state, dx, searchfail] = linesearchADI(state, dx, system, getEqs, upState, false);
end

% Update reservoir conditions once a delta has been found.
if searchfail
    [dx, meta] = stabilizeNewton(dx, meta, system);
    % If the line search failed, uncritically accept the first step and
    % pray the other measures (relaxation / dampening) handle the error.
    [state, nInc] = updateState(state, dx, opt);
end


[meta, residuals] = getResiduals(meta, eqs, system, solver_diverged);


%%{
    %residuals = cellfun(@(x) norm(x.val, 'inf'), eqs);
    meta.converged = all(residuals < system.nonlinear.tol);    
%}


if opt.Verbose
    
    if meta.iteration == 1
        eqnnames = {'Oil', 'Water', 'T', 'qOs', 'qWs', 'pBHP'};
        if(opt.temperature)
            eqnnames{end+1} = 'T';
        end
        
        if(opt.minerals)
            eqnnames{end+1} = 'MI';
        end
        fprintf('%-9s', eqnnames{:})
        fprintf('\n');
    end
    fprintf('%8.2e ', residuals);
    fprintf('\n');
end


meta.stopped = meta.iteration == system.nonlinear.maxIterations && ~meta.converged;

end

%--------------------------------------------------------------------------

function [state, nInc] = updateState(state, dx, opt)
dsMax = .2;
dpMax = 0.3;
% if ~isempty(phi)
%     for i = 1:numel(dx)
%         dx{i} = phi.basis{i}*dx{i};
%     end
% end

dp = dx{1};
ds = dx{2};
nInc = max( norm(dp,'inf')/norm(state.pressure, 'inf'), ...
            norm(ds,'inf')/norm(state.s(:,1), 'inf') );

%maxch = norm(ds, 'inf');
%step = min(1, maxSatStep./maxch);

ds = sign(ds).*min(abs(ds), dsMax);
%dp = sign(dp).*min(abs(dp), abs(dpMax.*state.pressure));
dp = sign(dp).*min(abs(dp), abs(10*barsa));

state.pressure = state.pressure + dp;
sw = state.s(:,1) + ds;
% Cap values
sw = min(sw, 1); sw = max(sw, 0);

state.s = [sw, 1-sw];
dT=dx{3};
dT=sign(dT).*min(abs(dT),10);
%if(opt.temperature)
    state.T=state.T+dT;
%end
dqWs    = dx{4};
dqOs    = dx{5};
dpBHP   = dx{6};
if(prod(size(dpBHP))>0)
    dpBHP = sign(dpBHP).*min(abs(dpBHP), abs(dpMax.*vertcat(state.wellSol.bhp)));
end
for w = 1:numel(state.wellSol)
    state.wellSol(w).bhp      = state.wellSol(w).bhp + dpBHP(w);
    state.wellSol(w).qWs      = state.wellSol(w).qWs + dqWs(w);
    state.wellSol(w).qOs      = state.wellSol(w).qOs + dqOs(w);
end
end
