function total_error = objectivefun_sync_multi_obj(x, model_flooding, ...
    model_cent)
%
% DESCRIPTION: defines the objective function including the simultaneous
%              case, and homogeneous case for multi objective optimization 
%              using genetic algorithm
%
% SYNOPSIS:
%   total_error = objectivefun_sync_multi_obj(x, model_flooding, ...
%     model_cent)
%
% PARAMETERS:
%   - x: the parameters on which we are iterating
%   model_flooding - primary model
%   model_cent - optional, secondary model in case of a simultaneous
%   history matching
%
% RETURNS:
%   total_error - error with the experimental measurements
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
obj_fun_type = model_flooding.history_match.obj_fun;
try
    model_flooding = run_history_match_homogeneous(x, model_flooding, 'ss');

    if isfield(model_flooding.experiment.observation,'pressure')
        if isfield(model_flooding.experiment.observation.pressure,'table')
            pressure_error_SS = calculate_pressure_error(model_flooding);
        end
    end
    if isfield(model_flooding.experiment.observation,'swavg')
        if isfield(model_flooding.experiment.observation.swavg,'table')
            swavg_error_SS = calculate_swavg_error(model_flooding);
        end
    end
    if isfield(model_flooding.experiment.observation,'prod')
        if isfield(model_flooding.experiment.observation.prod,'table')
            prod_error_SS = calculate_prod_error(model_flooding);
        end
    end
    if isfield(model_flooding.experiment.observation,'satProfile')
        if isfield(model_flooding.experiment.observation.satProfile,'table')
            satprof_error_SS = calculate_sat_profile_error(model_flooding);
        end
    end
    if strcmp(obj_fun_type,'Simultaneous')
        model_cent = run_history_match_homogeneous(x, model_cent, 'cent');
        
        if isfield(model_cent.experiment.observation,'swavg')
            if isfield(model_cent.experiment.observation.swavg,'table')
                swavg_error_cent = calculate_swavg_error(model_cent);
            end
        end
        if isfield(model_cent.experiment.observation,'prod')
            if isfield(model_cent.experiment.observation.prod,'table')
                prod_error_cent = calculate_prod_error(model_cent);
            end
        end
    end
    total_error = [];
    if exist('pressure_error_SS','var'), total_error = [total_error; pressure_error_SS]; end
    if exist('swavg_error_SS','var'), total_error = [total_error; swavg_error_SS]; end
    if exist('prod_error_SS','var'), total_error = [total_error; prod_error_SS]; end
    if exist('satprof_error_SS','var'), total_error = [total_error; satprof_error_SS]; end
    if exist('swavg_error_cent','var'), total_error = [total_error; swavg_error_cent]; end
    if exist('prod_error_cent','var'), total_error = [total_error; prod_error_cent]; end
catch  
    obs_no = length(fieldnames(model_flooding.experiment.observation));
    if isfield(model_cent, 'experiment')
       obs_no = obs_no + length(fieldnames(model_cent.experiment.observation));
    end
    total_error = ones(obs_no,1)*1e20;
end

