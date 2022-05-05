function prod_error = calculate_prod_error(model)
%
% DESCRIPTION: calculates the error between the simulated and experimental
%              measured production
%
% SYNOPSIS:
%   prod_error = calculate_prod_error(model)
%
% PARAMETERS:
%   - model: struct on which the history match is running
%
% RETURNS:
%   prod_error - production error with the experimental measurements
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
% get the data from model
cumTime = model.dynamic.params.cumScheduleSteps; 
calc_production = model.history_match.Qp_net;
observed_t_production = model.experiment.observation.prod.table{:,1};
observed_production = model.experiment.observation.prod.table{:,2};
observed_production(observed_production == 0) = NaN;

% make the data vector to the same size
interped_production = interp1(cumTime,calc_production,observed_t_production,'linear','extrap');

% an idea to normalize the data
% prod_combined = [observed_production; interped_production]; 
% norm_prod_combined = normalize(prod_combined,'range');
% observed_production = norm_prod_combined(1:height(observed_production));
% interped_production = norm_prod_combined(height(observed_production)+1:end);

if isfield(model.history_match, 'mcmc')
    mcmc_state = model.history_match.mcmc;
    error = model.history_match.prod_error;
    assert(model.history_match.prod_error ~= 0, 'The production error for mcmc must not be zero')
else
    mcmc_state = false;
    error = 0;
    assert(model.history_match.prod_error == 0, 'The production error for history match must be zero')
end

prod_error = error_build_optimization(interped_production, observed_production, ...
    mcmc_state, error) * model.history_match.prod_weight;
    
% helper plot to check the data
% figure; plot(observed_t_production, [observed_production interped_production]*1e6)

