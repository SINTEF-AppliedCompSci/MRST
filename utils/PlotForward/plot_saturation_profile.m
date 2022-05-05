function plot_saturation_profile (model)
%
% DESCRIPTION: plot the saturation profile after the simulation is done and
% compares it with the measured data, at the end of each flooding/rotation
% schedule
%
% SYNOPSIS:
%   plot_saturation_profile (model)
%
% PARAMETERS:
%   model - struct containing following fields:
%   - experiment: saturation functions used for forward modeling
%   - history_match: parameters used in gradient based and MCMC
%   - dynamic: information saved during the simulation
%
% RETURNS:
%   plot of the saturation profiles after simulation is completed
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
    core_length = model.experiment.geometry.length.value / Convert('cm');
    calcsatProfile = model.history_match.Sw_profile(2:end,2:end);
    calcsatProfile_location = model.history_match.Sw_profile(1,2:end);
    cumTime = model.dynamic.params.cumScheduleSteps;
    if(isfield(model.experiment.observation,'satProfile')) && model.experiment.observation.satProfile.include
        observed_sat_profile_t = model.experiment.observation.satProfile.table(2:end,1);
        observed_fractional_flows = time_to_fractional_flow(observed_sat_profile_t,model);
        observed_sat_profile_location = model.experiment.observation.satProfile.table(1,2:end);
        observed_sat_profile = model.experiment.observation.satProfile.table(2:end,2:end);
        interped_Satprof = interp2(calcsatProfile_location,cumTime,calcsatProfile,...
            observed_sat_profile_location,observed_sat_profile_t,'linear');
        assert(not(any(isnan(interped_Satprof),'all')), 'There are observations in the saturation profile which are not simulated!')
        figure('Name','Saturation Profile','NumberTitle','off','WindowStyle','docked')
        p_obs = plot(observed_sat_profile_location / Convert('cm'), observed_sat_profile,...
            's','DisplayName','Observed Saturation Profile','LineWidth',2, 'MarkerSize',3);
        hold on
        p_sim = plot(observed_sat_profile_location / Convert('cm'), interped_Satprof, ...
            'DisplayName','Simulated Saturation Profile','LineWidth',2, 'MarkerSize',3);
        for i = 1:length(observed_fractional_flows)
            p_sim(i).Color = p_obs(i).Color;
        end
        legend(repmat(string(observed_fractional_flows),2,1),'Location','bestoutside')
        xlabel('Position (cm)')
        ylabel('Water Saturation')
        ylim([0 1]); grid on; xlim([0 core_length])
    else
        endofsched_sat_profile_t = model.history_match.endofsched_satprofile(2:end,1);
        endofsched_fractional_flows = time_to_fractional_flow(endofsched_sat_profile_t,model);
        endofsched_sat_profile_location = model.history_match.endofsched_satprofile(1,2:end);
        endofsched_sat_profile = model.history_match.endofsched_satprofile(2:end,2:end);
        color = rand(length(endofsched_fractional_flows),3);
        figure('Name','Saturation Profile','NumberTitle','off','WindowStyle','docked')
        if model.experiment.rock.heterogeneous
            size_satprof = size(endofsched_sat_profile);
            for i = 1:size_satprof(1)
                stairs(endofsched_sat_profile_location / Convert('cm'),...
                    endofsched_sat_profile(i,:),'-', 'DisplayName','Saturation Profile',...
                    'LineWidth',1, 'MarkerSize',2);
                hold on
            end
        else
            p_obs = plot(endofsched_sat_profile_location / Convert('cm'),...
                endofsched_sat_profile,'-', 'DisplayName','Saturation Profile',...
                'LineWidth',1, 'MarkerSize',2);
            for i = 1:length(endofsched_fractional_flows)
                p_obs(i).Color = color(i,:);
            end
        end
        legend(repmat(string(endofsched_fractional_flows),1,1),'Location','bestoutside')
        xlabel('Position (cm)')
        ylabel('Water Saturation')
        ylim([0 1]); grid on; xlim([0 core_length])
    end
end

function factional_flow = time_to_fractional_flow(time,model)
% Calculates the fractional flow for each schedule row
    table_size = size(model.experiment.schedule.table);
    factional_flow = [];
    for j = 1:length(time)
        for i = 1:table_size(1)-1
            if time(j) > model.experiment.schedule.table{i,1} && time(j) ...
                    <= model.experiment.schedule.table{i+1,1}
                factional_flow = [factional_flow; ...
                    model.experiment.schedule.table{i,4}];
            end
        end
    end
end