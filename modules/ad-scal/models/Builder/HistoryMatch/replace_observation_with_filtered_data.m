function model = replace_observation_with_filtered_data(model, filtered_data)
%
% DESCRIPTION: replaces the results of the median filter from 
%              PlotObservation_pre_hm module in the model struct to be used during the
%              history matching process
%
% SYNOPSIS:
%   model = replace_observation_with_filtered_data(model, filtered_data)
%
% PARAMETERS:
%   model - struct containing following fields:
%   - experiment
%   - filtered_data: median filter applied on the experimental data
%
% RETURNS:
%   model - struct containing following fields:
%   - experiment: replaces the experimental measurements with the 
%   median filter results
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
if isfield(filtered_data, 'pressure')
    model.experiment.observation.pressure.table{:,2} = ...
        filtered_data.pressure;
    nan_mask = isnan(filtered_data.pressure);
    model.experiment.observation.pressure.table(nan_mask,:) = [];
end

if isfield(filtered_data, 'sw_avg')
    model.experiment.observation.sw_avg.table{:,2} = ...
        filtered_data.sw_avg;
    nan_mask = isnan(filtered_data.sw_avg);
    model.experiment.observation.sw_avg.table(nan_mask,:) = [];
end

if isfield(filtered_data, 'prod')
    model.experiment.observation.prod.table{:,2} = ...
        filtered_data.prod;
    nan_mask = isnan(filtered_data.prod);
    model.experiment.observation.prod.table(nan_mask,:) = [];
end