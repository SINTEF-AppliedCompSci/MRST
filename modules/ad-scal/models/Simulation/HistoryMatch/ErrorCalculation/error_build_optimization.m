function error = error_build_optimization(simulation, observation, mcmc_state, error)
%
% DESCRIPTION: calculates the error between the experimental measurements
%              and the simulation prediction - basically the objective
%              function formula is defined here
%
% SYNOPSIS:
%   error = error_build_optimization(simulation, observation, mcmc_state, error)
%
% PARAMETERS:
%   - simulation: simulation predictions
%   - observation: experimental measurements (data to match)
%   - mcmc_state: boolean, defining if we are in a MCMC simulation (true)
%   or another type of history matching (false)
%   - error: relative error for the data (input for the MCMC simulations
%   otherwise use 0)
%
% RETURNS:
%   error - the error between the experimental measurements and the
%   simulation predictions
%
% ----------------------------------
% (c) 2020-2022
% Siroos Azizmohammadi
% Omidreza Amrollahinasab
% Montanuniversität Leoben, Austria
% Chair of Reservoir Engineering
% https://dpe.ac.at/
% ----------------------------------
%
%%
if mcmc_state
    error = -0.5 * sum( rmmissing( ( (simulation - observation)...
        ./ (error / 100 * observation ) ).^2 ), 'all' ) / numel(observation);
    % to test the mcmc objective with HM
%     error = 0.5 * sum( rmmissing( ( (simulation - observation)...
%         ./ (error / 100 * observation ) ).^2 ), 'all' ) / numel(observation);
else
    error_temp = (simulation - observation) ./ observation;
    error = rms(rmmissing( error_temp(not(isinf(error_temp))) ),'all');
end