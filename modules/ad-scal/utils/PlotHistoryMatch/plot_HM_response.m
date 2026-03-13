function plot_HM_response(model, experiment, varargin)
%
% DESCRIPTION: after the history matching is done, plots the simulation 
%              predictions using the history match resutls and compares 
%              them with the experimental measurements in the model 
%
% SYNOPSIS:
%   plot_HM_response(model, experiment, varargin)
%
% PARAMETERS:
%   - model - struct after the history matching is done
%   - experiment - in case of a simultaneous history matching use 'ss' for
%   flooding experiments and 'cent' for centrifuge experiments - can be
%   left empty for non-simultaneous history matched
%   - optional: in case of heterogeneous history matching, the indexes of
%   the saturation profiles used
%
% RETURNS:
%   comparasion between simulation prediction, using history match results,
%   and the experimental measurements
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
    x = model.history_match.x;
    obj_fun = model.history_match.obj_fun;
    if strcmpi(obj_fun,'simultaneous') && strcmpi(experiment,'ss')
        type = 'ss';
    else
        type = '';
    end
    if model.experiment.rock.heterogeneous
        model = run_history_match_heterogeneous(x, model, type);
    else
        model = run_history_match_homogeneous(x, model, type);
    end

    calcp = model.history_match.calculatedp;
    calcSwavg = model.history_match.SwAvg;
    calc_production = model.history_match.Qp_net;
    calcSatProfile = model.history_match.Sw_profile(2:end,2:end);
    calcsatProfile_location = model.history_match.Sw_profile(1,2:end);
    cumTime = model.dynamic.params.cumScheduleSteps;
    
    tile_width = 1;
    tile_height = 1;
    obs_no = length(fieldnames(model.experiment.observation));
    if obs_no > 1; tile_width = 2; end
    if obs_no > 2; tile_height = 2; end
    figure('Name','History match predictions','NumberTitle','off')
    tiledlayout(tile_width, tile_height);
    if isfield(model.experiment.observation,'pressure')
        if isfield(model.experiment.observation.pressure,'table')
            observed_t_p = model.experiment.observation.pressure.table{:,1};
            observedp = model.experiment.observation.pressure.table{:,2};
            interped_calculatedp = interp1(cumTime,calcp,observed_t_p,'linear','extrap');
            nexttile
            hold on
            observed_t_p = observed_t_p / 3600;
            observedp = observedp * 1e-5; interped_calculatedp = interped_calculatedp * 1e-5;
            plot(observed_t_p,observedp, 'ro', 'MarkerSize',2)
            plot(observed_t_p,interped_calculatedp, 'b-', 'LineWidth',2)
            legend('Measured Pressure', 'History match pressure prediction',...
                'location','best')
            title('Pressure Difference')
            ylabel('Pressure (bar)')
            xlabel('Simulation Time (hour)')
            grid('on')
        end
    end
    if isfield(model.experiment.observation,'swavg')
        if isfield(model.experiment.observation.swavg,'table')
            observed_t_swavg = model.experiment.observation.swavg.table{:,1};
            observedSwAvg = model.experiment.observation.swavg.table{:,2};
            interped_Swavg = interp1(cumTime,calcSwavg,observed_t_swavg,'linear','extrap');
            nexttile
            hold on 
            observed_t_swavg = observed_t_swavg / 3600;
            plot(observed_t_swavg,observedSwAvg, 'ro', 'MarkerSize',2)
            plot(observed_t_swavg,interped_Swavg, 'b-', 'LineWidth',2)
            legend('Measured Saturation', 'History match saturation prediction',...
                'location','best')
            title('Average Water Saturation')
            ylabel('Average Water Saturation (fraction)')
            xlabel('Simulation Time (hour)')
            grid('on')
        end
    end
    if isfield(model.experiment.observation,'satProfile')
        if isfield(model.experiment.observation.satProfile,'table')
            if nargin > 2
                index_mask = varargin{1};
            end

            observed_sat_profile_t = model.experiment.observation.satProfile.table(2:end,1);
            observed_fractional_flows = time_to_fractional_flow(observed_sat_profile_t,model);
            observed_sat_profile_location = model.experiment.observation.satProfile.table(1,2:end);
            observed_sat_profile = model.experiment.observation.satProfile.table(2:end,2:end);
            interped_Satprof = interp2(calcsatProfile_location,cumTime,calcSatProfile,...
                observed_sat_profile_location,observed_sat_profile_t,'linear');
            nexttile
            if nargin > 2
                index_mask = varargin{1};
            else
                index_mask = 1:height(observed_sat_profile);
            end
            p_obs = plot(observed_sat_profile_location * 100, observed_sat_profile(index_mask,:),...
                'o','DisplayName','Measured Saturation Profile', 'MarkerSize',4);
            hold on
            p_sim = plot(observed_sat_profile_location * 100, interped_Satprof(index_mask,:),...
                'DisplayName','History match saturation profile prediction','LineWidth',2);
            if nargin > 2
                for i = 1:length(index_mask)
                    p_sim(i).Color = p_obs(i).Color;
                end
                leg = legend(num2cell(string(round(observed_sat_profile_t(index_mask)...
                /3600',2))),'Location','bestoutside');
                title(leg,'Time (hour)')
            else
                for i = 1:length(observed_fractional_flows)
                    p_sim(i).Color = p_obs(i).Color;
                end
                legend(repmat(string(observed_fractional_flows),2,1),...
                    'Location','bestoutside')
            end
            xlabel('Position (cm)')
            ylabel('Water Saturation')
            grid('on')
        end
    end
    if isfield(model.experiment.observation,'prod')
        if isfield(model.experiment.observation.prod,'table')
            observed_t_production = model.experiment.observation.prod.table{:,1};
            observed_production = model.experiment.observation.prod.table{:,2};
            interped_production = interp1(cumTime,calc_production,observed_t_production,'linear','extrap');
            nexttile
            hold on 
            observed_t_production = observed_t_production / 3600;
            observed_production = observed_production * 1e6; interped_production = interped_production * 1e6;
            plot(observed_t_production,observed_production, 'ro', 'MarkerSize',2)
            plot(observed_t_production,interped_production, 'b-', 'LineWidth',2)
            legend('Measured Production', 'History match production prediction',...
                'location','best')
            title('Production')
            ylabel('Production (cm^3)')
            xlabel('Simulation Time (hour)')
            grid('on')
        end
    end
end

function factional_flow = time_to_fractional_flow(time,model)
    table_size = size(model.experiment.schedule.table);
    factional_flow = [];
    for j = 1:length(time)
        for i = 1:table_size(1)-1
            if time(j) > model.experiment.schedule.table{i,1} && time(j) ...
                    <= model.experiment.schedule.table{i+1,1}
                factional_flow = [factional_flow;...
                    model.experiment.schedule.table{i,4}];
            end
        end
    end
end
