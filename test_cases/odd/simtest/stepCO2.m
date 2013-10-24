function [state, meta] = stepCO2(state0, state, meta, dt, W, G, system, fluid)

[eqs, info] = system.getEquations(state0, state, dt, G, W, system.s, fluid);
dx  = SolveEqsADI(eqs, []);

%% store residual history (required by the 'solvefiADI' function)
residuals = cellfun(@(x) norm(x.val, 'inf'), eqs);
if ~isfield(meta, 'res_history')
    meta.res_history = zeros(system.nonlinear.maxIterations, numel(residuals));
end
meta.res_history(meta.iteration,:) = residuals;
    
%% determine if converged
[converged CNV MB] = getCO2Convergence(state, eqs, fluid, system, dt); 
                                                                       
meta.converged = converged;
meta.stopped = meta.iteration == system.nonlinear.maxIterations && ~ ...
    converged;

%% updating state

maxHStep = .2 * min(G.cells.H); %@@ this should perhaps be refined
dp = dx{1};
dh = dx{2};
maxch = norm(dh, 'inf');
step = min(1, maxHStep./maxch);

state.pressure = state.pressure + step*dp;
state.h = state.h + step * dh;

% cap height values
state.h(state.h < 0) = 0;

hix = find(state.h > G.cells.H);
state.h(hix) = G.cells.H(hix);

% additional info from 'getEquations' that should be passed on to the user
state.info = info;


%@@ should we put well solution stuff in here (as is done in stepOG).
