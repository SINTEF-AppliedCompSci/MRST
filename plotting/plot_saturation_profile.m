function plot_saturation_profile (model)
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
        assert(not(any(isnan(interped_Satprof),'all')), 'Saturation profile observation time or location is out of simulation range')
        color = rand(length(observed_fractional_flows),3);
        figure('Name','Saturation Profile','NumberTitle','off')
        p_obs = plot(observed_sat_profile_location / Convert('cm'), observed_sat_profile,...
            's','DisplayName','Observed Saturation Profile','LineWidth',2, 'MarkerSize',3);
        for i = 1:length(observed_fractional_flows)
            p_obs(i).Color = color(i,:);
        end
        hold on
        p_sim = plot(observed_sat_profile_location / Convert('cm'), interped_Satprof, ...
            'DisplayName','Simulated Saturation Profile','LineWidth',2, 'MarkerSize',3);
        for i = 1:length(observed_fractional_flows)
            p_sim(i).Color = color(i,:);
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
        figure('Name','Saturation Profile','NumberTitle','off')
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