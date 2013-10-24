function [converged CNV MB] = getCO2Convergence(state, eqs, fluid, system, dt)

tol_mb = system.nonlinear.tolMB;
tol_cnv = system.nonlinear.tolCNV;
pvsum = sum(system.s.pv);

% We assume eqs{1} represents CO2 and eqs{2} represents water.

W_res = eqs{2}.val;   % residual of water
C_res =eqs{1}.val;    % residual of CO2

assert(isa(fluid.CO2.rho, 'function_handle'));
assert(isa(fluid.water.rho, 'double'));

temperature = system.s.top_temp + (system.s.temp_grad * state.h / 1000);
rhoC = fluid.CO2.rho(state.pressure, temperature);
rhoW = fluid.water.rho * ones(size(state.pressure));

CNVW = dt * max(abs(W_res.*rhoW)./system.s.pv);
CNVC = dt * max(abs(C_res.*rhoC)./system.s.pv);
%CNVC = fluid.rhoC * dt * max(abs(C_res)./system.s.pv);

% Check if material balance for each phase fulfils residual convergence
% criterion
MB = abs([rhoW' * sum(W_res), rhoC' * C_res]); 
%MB = abs([fluid.rhoW * sum(W_res), fluid.rhoC * sum(C_res)]);
converged_MB = all(MB < tol_mb*pvsum/dt);

%Check maximum normalized residuals (maximum mass error)
CNV = [CNVW, CNVC];
converged_CNV = all(CNV < tol_cnv);

converged = converged_MB & converged_CNV;

