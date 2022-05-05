function swavg_error_SS = calculate_swavg_error(model)
%
% DESCRIPTION: calculates the error between the simulated and experimental
%              measured average water saturation
%
% SYNOPSIS:
%   swavg_error_SS = calculate_swavg_error(model)
%
% PARAMETERS:
%   - model: struct on which the history match is running
%
% RETURNS:
%   swavg_error_SS - average water saturation error with the 
%   experimental measurements
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
calcSwavg = model.history_match.SwAvg;
observed_t_swavg = model.experiment.observation.swavg.table{:,1};
observedSwAvg = model.experiment.observation.swavg.table{:,2};
observedSwAvg(observedSwAvg == 0) = NaN;

% make the data vector to the same size
interped_Swavg = interp1(cumTime,calcSwavg,observed_t_swavg,'linear','extrap');


if isfield(model.history_match, 'mcmc')
    mcmc_state = model.history_match.mcmc;
    error = model.history_match.swavg_error;
    assert(model.history_match.swavg_error ~= 0, 'The sw average error for mcmc must not be zero')
else
    mcmc_state = false;
    error = 0;
    assert(model.history_match.swavg_error == 0, 'The sw average error for history match must be zero')
end

swavg_error_SS = error_build_optimization(interped_Swavg, observedSwAvg, ....
    mcmc_state, error) * model.history_match.swavg_weight;

% helper plot to check the data
% figure; plot(observed_t_swavg, [observedSwAvg interped_Swavg])
