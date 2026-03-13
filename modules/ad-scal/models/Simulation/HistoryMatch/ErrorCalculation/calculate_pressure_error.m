function pressure_error = calculate_pressure_error(model)
%
% DESCRIPTION: calculates the error between the simulated and experimental
%              measured pressure
%
% SYNOPSIS:
%   pressure_error = calculate_pressure_error(model)
%
% PARAMETERS:
%   - model: struct on which the history match is running
%
% RETURNS:
%   pressure_error - pressure error with the experimental measurements
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
calcp = model.history_match.calculatedp;
observed_t_p = model.experiment.observation.pressure.table{:,1};
observedp = model.experiment.observation.pressure.table{:,2};
observedp(observedp == 0) = NaN;

% make the data vector to the same size
interped_calculatedp = interp1(cumTime,calcp,observed_t_p,'linear','extrap');

% an idea to normalize the data
% pressure_combined = [observedp; interped_calculatedp]; 
% norm_pressure_combined = normalize(pressure_combined,'range');
% observedp = norm_pressure_combined(1:height(observedp));
% interped_calculatedp = norm_pressure_combined(height(observedp)+1:end);


if isfield(model.history_match, 'mcmc')
    mcmc_state = model.history_match.mcmc;
    error = model.history_match.pdiff_error;
    assert(model.history_match.pdiff_error ~= 0, 'The pressure error for mcmc must not be zero')
else
    assert(model.history_match.pdiff_error == 0, 'The pressure error for history match must be zero')
    mcmc_state = false;
    error = 0;
end

pressure_error = error_build_optimization(interped_calculatedp, observedp, ....
    mcmc_state, error) * model.history_match.pdiff_weight;

if isfield(model.experiment.observation, "pressure_mid")
    calcp_mid = model.dynamic.params.pDiff_mid;
    observed_t_p_mid = model.experiment.observation.pressure_mid.table{:,1};
    observedp_mid = model.experiment.observation.pressure_mid.table{:,2};
    observedp_mid(observedp_mid == 0) = NaN;
    interped_calculatedp_mid = interp1(cumTime, calcp_mid, observed_t_p_mid,'linear','extrap');
    pressure_error = pressure_error + error_build_optimization(interped_calculatedp_mid, observedp_mid, ....
    mcmc_state, error) * model.history_match.pdiff_weight;
end

% helper plot to check the data
% figure; plot(observed_t_p, [observedp interped_calculatedp]/1e5)
