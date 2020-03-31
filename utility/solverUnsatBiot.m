function [p, p_m, u, iter] = solverUnsatBiot(G, p_n, u_n, modelEqs, ...
    time_param, solver_param, sourceFlow, sourceMech)
% Newton solver for the equations of unsaturated poroelasticity
%
% SYNOPSIS:
%   function [p, p_m, u, iter] = solverUnsatBiot(G, p_n, u_n, modelEqs, ...
%       time_param, solver_param, sourceFlow, sourceMech)
%
% PARAMETERS:
%   G             - Structure, grid structure from MRST
%   p_n           - Vector, pressure evaluated at the last time level
%   u_n           - Vector, displacement evaluated at the last time level
%   modelEqs      - Structure, containing the model equations
%   time_param    - Structure, containing the time parameters
%   solver_param  - Structure, containing the solver parameters
%   sourceFlow    - Vector, source term for the flow equation  
%   sourceMech    - Vector, source term for the mechanics equation  
%
% RETURNS:
%   p_ad          - Vector, updated pressure
%   p_m           - Vector, pressure evaluated at the last iteration
%   u             - Vector, updated displacement
%   iter          - Scalar, number of iterations needed for convergence
%
% See also solverRE.

%{
Copyright 2018-2019, University of Bergen.

This file is part of the fv-unsat module.

fv-unsat is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

fv-unsat is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file.  If not, see <http://www.gnu.org/licenses/>.
%} 

res = 100; % initializing residual
iter = 1; % maximum number of iterations
Nd = G.griddim; % grid dimension
Nc = G.cells.num; % number of cells

% Initializing AD-variables
p_ad = initVariablesADI(p_n);
u_ad = initVariablesADI(u_n);

% Newton loop
while (res > solver_param.tol) && (iter <= solver_param.maxIter)
    
    p_m = p_ad.val; % current iteration level (m-index)
    
    % Calling equations
    eq1 = modelEqs.uEq1(u_ad);
    eq2 = modelEqs.uEq2(p_ad, p_n, sourceMech);
    eq3 = modelEqs.pEq1(p_n, u_ad, u_n);
    eq4 = modelEqs.pEq2(p_ad, p_n, p_m, time_param.tau, sourceFlow);
    
    J = [eq1.jac{1} eq2.jac{1}; eq3.jac{1} eq4.jac{1}]; % Jacobian
    R = [eq1.val + eq2.val; eq3.val + eq4.val]; % Residual
    Y = J \ -R; % solve linear system
    u_ad.val = u_ad.val + Y(1:Nd*Nc); % update u
    p_ad.val = p_ad.val + Y(Nd*Nc+1:end); % update p
    res = norm(R); % compute tolerance
    
    % Checking convergence...
    if res <= solver_param.tol
        fprintf('Time: %.2f \t Iter: %d \t Error: %.2e \n', ...
            time_param.time, iter, res);
    elseif iter >= solver_param.maxIter
        error('Solver failed to converge. Try decreasing tol or increasing maxIter.');
    else
        iter = iter + 1;
    end
    
end

p = p_ad.val; % updating pressure value
u = u_ad.val; % updating displacement value