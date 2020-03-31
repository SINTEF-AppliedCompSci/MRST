function [psi, psi_m, iter] = solverRE(psi_n, modelEqs, time_param, ...
    solver_param, source)
% Newton solver for Richards' equation
%
% SYNOPSIS:
%   function [psi, psi_m, iter] = solverRE(psi_n, modelEqs, time_param, ... 
%       solver_param, source)
%
% PARAMETERS:
%   psi_n        - Vector, pressure head evaluated at the last time level
%   modelEqs     - Structure, containing the model equations
%   time_param   - Structure, containing the time parameters
%   solver_param - Structure, containing the solver parameters
%   source       - Vector, source term for the flow equation  
%
% RETURNS:
%   psi          - Vector, updated pressure
%   psi_m        - Vector, pressure head evaluated at the last iteration
%   iter         - Scalar, number of iterations needed for convergence
%
% See also solverUnsatBiot.

%{
Copyright 2018-2020, University of Bergen.

This file is part of the fv-re-biot module.

fv-re-biot is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

fv-re-biot is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file.  If not, see <http://www.gnu.org/licenses/>.
%} 

res = 100; % residual value
iter = 1; % iterations

% Initialiazing AD-variable
psi_ad = initVariablesADI(psi_n); 

% Newton loop
while (res > solver_param.tol) && (iter <= solver_param.maxIter)
    
    psi_m = psi_ad.val; % current iteration level (m-index)
    eq = modelEqs.psiEq(psi_ad, psi_n, psi_m, time_param.tau, ...
        source); % call equation from model
    R = eq.val; % residual
    J = eq.jac{1}; % Jacobian
    Y = J\-R; % solve linear system
    psi_ad.val  = psi_ad.val + Y; % update
    res = norm(R); % compute tolerance
    
    % Checking convergence
    if res <= solver_param.tol
        fprintf('Time: %.2f \t Iter: %d \t Error: %.2e \n',...
            time_param.time, iter, res);
    elseif iter >= solver_param.maxIter
        error('Solver failed to converge. Try decreasing tol or increasing maxIter.');
    else
        iter = iter + 1;
    end
end

psi = psi_ad.val; % return updated pressure head