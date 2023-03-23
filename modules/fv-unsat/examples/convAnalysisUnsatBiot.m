%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               fv-unsat                                  %              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example 2: Convergence analysis for the unsaturated poroelastic equations
% Author: Jhabriel Varela. E-mail: Jhabriel.Varela@uib.no. 

%{
Copyright 2018-2020, University of Bergen.

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

%% Clearing workspace and cleaning console
clear; clc();

%% Performing the converge test

% We are interested in the convergence rate for the pressure, displacement,
% fluxes and tractions. To obtained them, we perform the analysis using 5
% levels of spatial and time refinements. 

numCells   = [5, 10, 20, 40, 80, 160]'; % number of cells
timeLevels = [10, 10, 10, 10, 10, 10]'; % time levels

h = 1 ./ numCells; % spatial mesh size
tau = 1 ./ timeLevels; % time step

error_p = zeros(length(numCells), 1); % L2 error pressure
error_u = zeros(length(numCells), 1); % L2 error displacement
error_Q = zeros(length(numCells), 1); % L2 error flux
error_T = zeros(length(numCells), 1); % L2 error traction

% Looping to obtain the errors
for ii=1:length(error_p)
    fprintf('\nPerforming simulation for h = %f and dt = %f\n', ...
        h(ii), tau(ii));
    [error_p(ii), error_u(ii), error_Q(ii), error_T(ii)] = ...
        convergenceUnsatBiot(numCells(ii), timeLevels(ii));
end

% Computing the error reduction
Red_p = error_p(1:end-1) ./ error_p(2:end);
Red_u = error_u(1:end-1) ./ error_u(2:end);
Red_Q = error_Q(1:end-1) ./ error_Q(2:end);
Red_T = error_T(1:end-1) ./ error_T(2:end);

% Numerical convergence rates
rate_p = log2(Red_p);
rate_u = log2(Red_u);
rate_Q = log2(Red_Q);
rate_T = log2(Red_T);

%% Printing results
table(h, tau, ...
    error_p, [nan; Red_p], [nan; rate_p], ...
    error_u, [nan; Red_u], [nan; rate_u],...
'VariableNames', {...
    'Grid_size', 'Time_step', ...
    'Error_p', 'Red_p', 'Rate_p', ...
    'Error_u', 'Red_u', 'Rate_u',...
    })

table(h, tau, ...
    error_Q, [nan; Red_Q], [nan; rate_Q], ...
    error_T, [nan; Red_T], [nan; rate_T],...
'VariableNames', {...
    'Grid_size', 'Time_step', ...
    'Error_Q', 'Red_Q', 'Rate_Q', ...
    'Error_T', 'Red_T', 'Rate_T',...
    })
