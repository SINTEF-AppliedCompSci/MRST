function model = replace_observation_with_filtered_data(model, filtered_data)

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