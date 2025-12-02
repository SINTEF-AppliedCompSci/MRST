function satprof_error = calculate_sat_profile_error(model)
%
% DESCRIPTION: calculates the error between the simulated and experimental
%              measured water saturation profile
%
% SYNOPSIS:
%   satprof_error = calculate_sat_profile_error(model)
%
% PARAMETERS:
%   - model: struct on which the history match is running
%
% RETURNS:
%   satprof_error - water saturaiton profile error with the
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
calcSatProfile = model.history_match.Sw_profile(2:end,2:end);
calcsatProfile_location = model.history_match.Sw_profile(1,2:end);
observed_sat_profile_t = model.experiment.observation.satProfile.table(2:end,1);
observed_sat_profile_location = model.experiment.observation.satProfile.table(1,2:end);
observed_sat_profile = model.experiment.observation.satProfile.table(2:end,2:end);
observed_sat_profile(observed_sat_profile == 0) = NaN;

% make the data vector to the same size
interped_Satprof = interp2(calcsatProfile_location,cumTime,calcSatProfile,...
    observed_sat_profile_location,observed_sat_profile_t,'linear');


if isfield(model.history_match, 'mcmc')
    mcmc_state = model.history_match.mcmc;
    error = model.history_match.sat_profile_error;
    assert(model.history_match.sat_profile_error ~= 0, 'The saturation profile error for mcmc must not be zero')
else
    mcmc_state = false;
    error = 0;
    assert(model.history_match.sat_profile_error == 0, 'The saturation profile error for history match must be zero')
end

satprof_error = error_build_optimization(interped_Satprof, observed_sat_profile,...
    mcmc_state, error) * model.history_match.sat_profile_weight;