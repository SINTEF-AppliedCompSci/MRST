function model = run_history_match_homorogeneous(x, model, type)

% build the saturation functions
params_no_kr = get_params_no_kr(model);
model = Create_pc_history_match(x, model, params_no_kr);
model = Create_kr_history_match(x, model);

% run homogeneous simulation
model = CreateFluid(model);
if strcmpi(type, 'ss')
    model = change_config_to_SS(model);
end
model = Run(model, false);
