function [psi, psi_m, iter] = solverRE(psi_n, psiEq, tau, source, ...
    currentTime, tol, maxIter)
% Newton solver for the Richards' equation
%
% SYNOPSIS:
%   function [psi, psi_m, iter] = solverRE(psi_n, psiEq, tau, source, ..
%   currentTime, tol, maxIter)
%
% PARAMETERS:
%   psi_n        - Vector, pressure head evaluated at the last time level.
%   psiEq        - Function, pressure head equation, i.e. mass balance eq.
%   tau          - Scalar, time step
%   source       - Vector, source term for the flow equation  
%   currentTime  - Scalar, current simulation time
%   tol          - Scalar, convergence tolerance
%   maxIter      - Scalar, maximum number of iterations
%
%  RETURNS:
%   psi          - AD-object, updated pressure head AD-variable.
%   psi_m        - Vector, pressure head evaluated at the last iteration.
%   iter         - Scalar, number of iterations needed for convergence
%

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

res = 100;      % residual value
iter = 1;       % iterations

% Initaliazing AD-variable
psi_ad = initVariablesADI(psi_n); 

% Newton loop
while (res > tol) && (iter <= maxIter)
    
    psi_m = psi_ad.val; % current iteration level (m-index)
    eq = psiEq(psi_ad, psi_n, psi_m, tau, source); % call eq
    R = eq.val;       % residual
    J = eq.jac{1};    % Jacobian
    Y = J\-R;         % solve linear system
    psi_ad.val  = psi_ad.val + Y; % update
    res = norm(R); % compute tolerance
    
    % Checking convergence
    if res <= tol
        fprintf('Time: %.2f \t Iter: %d \t Error: %.2e \n',...
            currentTime, iter, res);
    elseif iter >= maxIter
        error('Solver failed to converge. Try decreasing tol or increasing maxIter.');
    else
        iter = iter + 1;
    end
end

psi = psi_ad.val; % updated pressure head