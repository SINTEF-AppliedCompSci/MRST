function model = run_history_match_heterogeneous(x, model, type)
%
% DESCRIPTION: run the simulation case for history matching purposes of the
%              heterogeneous case
%
% SYNOPSIS:
%   model = run_history_match_homogeneous(x, model, type)
%
% PARAMETERS:
%   - x: iteration parameters
%   - model: struct on which the history match is running
%   - type: string used to turn of the gravity for the flooding 
%   experiments - example: 'ss' 
%
% RETURNS:
%   model - struct with the simulation results inside
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
% build the saturation functions
params_no_kr = get_params_no_kr(model);
model = Create_pc_history_match(x, model, params_no_kr);
model = Create_kr_history_match(x, model);
alpha = get_alpha(x, model);

% run homogeneous simulation
model = CreateFluid(model);
if strcmpi(type, 'ss')
    model = change_config_to_SS(model);
end
model = Run(model, false);

% run heterogeneous simulation
f_factor = f_factor_calculator(model, alpha, false, model.experiment.rock.het_index_mask);
model.experiment.rock.heterogeneous = true;
model = CreateGrid_heterogeneous(model);
model = CreateRock_heterogeneous(model, f_factor);
model = CreateFluid_heterogenous(model, f_factor);
model = Run(model, false);
