function [h_ad,h_m0,iter,timeCum,tf] = newtonSolver(h_ad,h_eq,timeCum,dt,maxTolPresHead,maxIter)
% Newton-AD (Automatic Differentiation) solver
%
% SYNOPSIS:
%   [h_ad,h_m0,iter,timeCum,tf] = newtonSolver(h_ad,h_eq,timeCum,dt,maxTolPresHead,maxIter)
%
% PARAMETERS:
%   h_ad            - AD variable, hydraulic head (primary variable)
%   h_eq            - Function, hydraulic head discrete equation
%   timeCum         - Scalar, cumulative time
%   dt              - Scalar, time step
%   maxTolPresHead  - Scalar, maximum tolerance of pressure head ~ 1 [cm]
%   maxIter         - Scalar, maximum iteration number ~ 10 
%
%  RETURNS:
%   h_ad            - AD variable, updated hydraulic head AD varible
%   h_m0            - Vector, containing pressure head of the last
%                     iteration level
%   iter            - Scalar, number of iterations
%   timeCum         - Scalar, updated cumulative time
%   tf              - Scalar, computational time

ti = tic();                     % saves current time (for CPU time)

% Newton parameters
absTolPresHead = 100;           % [cm] initializing maximum relative tolerance
iter = 1;                       % [iter] initializing iterations
h_n0 = h_ad.val;                % [cm] current time step h (n-index)
timeCum = timeCum + dt;         % [s] cumulative time

% Newton loop
while (absTolPresHead > maxTolPresHead) && (iter < maxIter)
    
    h_m0 = h_ad.val;                % current iteration step h (m-index)
    eq = h_eq(h_ad,h_n0,h_m0,dt);   % calling equation
    R = eq.val;                     % determing residual
    J = eq.jac{1};                  % determing Jacobian
    y = J\-R;                       % Newton update
    h_ad.val  = h_ad.val + y;       % root for the k-th step
    absTolPresHead = max(abs(h_ad.val - h_m0)); % checking convergence
    
    % Printing convergence info...
    if (absTolPresHead <= maxTolPresHead) && (iter <= maxIter)
        fprintf('Time: %.3f [s] \t Iter: %d \t Error: %d [cm] \n',...
            timeCum,iter,absTolPresHead);
    elseif (iter > maxIter)
        error('Newton method did not converge!');
    end
    
    iter = iter + 1;                % iteration ++
end

tf = toc(ti);

end

