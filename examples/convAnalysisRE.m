%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               fv-unsat                                  %                                                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example 1: Convergence analysis for the Richards' equation
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

% We are interested in the convergence rate for the pressure head and 
% the fluxes. In order to obtain them, we perform the analysis using five
% levels of spatial refinements for a fixed time step.

numCells = [10, 20, 40, 80, 160]'; % number of cells
timeLevels = [10, 10, 10, 10, 10]'; % time levels

h = 1 ./ numCells; % spatial mesh size
tau = 1 ./ timeLevels; % time step

errorPsi = zeros(length(numCells), 1); % L2 error psi
errorFlux = zeros(length(numCells), 1); % L2 error flux

% Looping to obtain the errors
for ii=1:length(errorPsi)
    fprintf('\nPerforming simulation for h = %f and dt = %f\n', ...
        h(ii), tau(ii));
    [errorPsi(ii), errorFlux(ii)] = ...
        convergenceRE(numCells(ii), timeLevels(ii));
end

% Computing the error reduction
RedPsi = errorPsi(1:end-1) ./ errorPsi(2:end);
RedFlux = errorFlux(1:end-1) ./ errorFlux(2:end);

% Numerical convergence rates
RatePsi = log2(RedPsi);
RateFlux = log2(RedFlux);

%% Printing results
table(h, tau, ...
    errorPsi, [nan; RedPsi], [nan; RatePsi], ...
    errorFlux,  [nan; RedFlux], [nan; RateFlux], ...
'VariableNames', {...
    'Grid_size', 'Time_step', ...
    'Error_psi', 'Red_psi', 'Rate_psi', ...
    'Error_flux', 'Red_flux', 'Rate_flux', ...
    })
