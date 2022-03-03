function plot_HM_response_multi_obj(model,experiment,index,observations)
% function plot_HM_response(model,experiment)
%     experiment = 'CENT';
    solutions = model.history_match.x(index,:);
    fval = model.history_match.fval(index,:);
    obj_fun = model.history_match.obj_fun;
    if strcmpi(obj_fun,'Simultaneous')
        if strcmpi(experiment,'CENT')
            model = change_config_to_cent(model);
        else
            model = change_config_to_SS(model);
        end
    end
    % f = waitbar(0, 'Starting');
    n = height(solutions);
    rng default
    color = rand(n,3);
    for j=1:n
        x = solutions(j,:);
        model = CreateSatFun_HM(x,model);
        model = CreateFluid(model);
        model = Run(model);
        calcp = model.history_match.calculatedp;
        calcSwavg = model.history_match.SwAvg;
        calc_production = model.history_match.Qp_net;
        calcSatProfile = model.history_match.Sw_profile(2:end,2:end);
        calcsatProfile_location = model.history_match.Sw_profile(1,2:end);
        cumTime = model.dynamic.params.cumScheduleSteps;
        if any(strcmpi(observations,'pressure'))
            if isfield(model.experiment.observation,'pressure')
                if isfield(model.experiment.observation.pressure,'table')
                    observed_t_p = model.experiment.observation.pressure.table{:,1};
                    observedp = model.experiment.observation.pressure.table{:,2};
                    interped_calculatedp = interp1(cumTime,calcp,observed_t_p,'linear','extrap');
                    if j==1
                        figure('Name','Pressure Difference','NumberTitle','off')
                        ax_pressure = gca;
                        hold(ax_pressure,'on')
                        plot(ax_pressure,observed_t_p / 3600,observedp * 1e-5, 'ko', 'MarkerSize',2, 'DisplayName', "Measured Pressure")
                        title('Pressure Difference')
                        ylabel('Pressure (bar)')
                        xlabel('Simulation Time (hour)')
                        grid('on')
                        legend show
                    end
                    plot(ax_pressure,observed_t_p / 3600 ,interped_calculatedp * 1e-5, '-', 'LineWidth',2,'DisplayName',...
                        "error = " + string(fval(j,1)) + " index = " + string(j) + " total error = " ...
                        + string(sum(fval(j,:),2)), 'color', color(j,:))
                end
            end
        end
        if any(strcmpi(observations,'prod'))
            if isfield(model.experiment.observation,'prod')
                if isfield(model.experiment.observation.prod,'table')
                    observed_t_production = model.experiment.observation.prod.table{:,1};
                    observed_production = model.experiment.observation.prod.table{:,2};
                    interped_production = interp1(cumTime,calc_production,observed_t_production,'linear','extrap');
                    if j==1
                        figure('Name','Production','NumberTitle','off')
                        ax_pod = gca;
                        hold(ax_pod,'on')
                        plot(ax_pod, observed_t_production/3600,observed_production*1e6, 'ko', 'MarkerSize',2, 'DisplayName', "Measured Production")
                        title('Production')
                        ylabel('Production (cm^3)')
                        xlabel('Simulation Time (hour)')
                        grid('on')
                        legend show
                    end
                    plot(ax_pod, observed_t_production/3600,interped_production*1e6, '-', 'LineWidth',2,'DisplayName',...
                        "error = " + string(fval(j,2)) + " index = " + string(j), 'color', color(j,:))
                end
            end
        end
        if any(strcmpi(observations,'swavg'))
            if isfield(model.experiment.observation,'swavg')
                if isfield(model.experiment.observation.swavg,'table')
                    observed_t_swavg = model.experiment.observation.swavg.table{:,1};
                    observedSwAvg = model.experiment.observation.swavg.table{:,2};
                    interped_Swavg = interp1(cumTime,calcSwavg,observed_t_swavg,'linear','extrap');
                    if j==1
                        figure('Name','Average Water Saturation','NumberTitle','off')
                        ax_swavg = gca;
                        hold(ax_swavg,'on')
                        plot(ax_swavg, observed_t_swavg/3600,observedSwAvg, 'ko', 'MarkerSize',2, 'DisplayName', "Measured Saturation")
                        title('Average Water Saturation')
                        ylabel('Average Water Saturation (fraction)')
                        xlabel('Simulation Time (hour)')
                        grid('on')
                        legend show
                    end
                    plot(ax_swavg, observed_t_swavg/3600,interped_Swavg, '-', 'LineWidth',2,'DisplayName',...
                        "error = " + string(fval(j,3)) + " index = " + string(j), 'color', color(j,:))
                end
            end
        end
        if any(strcmpi(observations,'satProfile')) 
            if isfield(model.experiment.observation,'satProfile')
                if isfield(model.experiment.observation.satProfile,'table')
                    observed_sat_profile_t = model.experiment.observation.satProfile.table(2:end,1);
                    observed_fractional_flows = time_to_fractional_flow(observed_sat_profile_t,model);
                    observed_sat_profile_location = model.experiment.observation.satProfile.table(1,2:end);
                    observed_sat_profile = model.experiment.observation.satProfile.table(2:end,2:end);
                    interped_Satprof = interp2(calcsatProfile_location,cumTime,calcSatProfile,...
                        observed_sat_profile_location,observed_sat_profile_t,'linear');
                    figure('Name','Saturation Profile','NumberTitle','off')
                    ax_sat_prof = gca;
                    hold on
                    color_satprof = rand(length(observed_fractional_flows),3);
                    p_obs = plot(ax_sat_prof, observed_sat_profile_location * 100, observed_sat_profile, 'o','DisplayName','Measured Saturation Profile', 'MarkerSize',4);
                    for i = 1:length(observed_fractional_flows)
                        p_obs(i).Color = color_satprof(i,:);
                    end
                    p_sim = plot(ax_sat_prof, observed_sat_profile_location * 100, interped_Satprof, 'DisplayName','Simulated Saturation Profile','LineWidth',2);
                    for i = 1:length(observed_fractional_flows)
                        p_sim(i).Color = color_satprof(i,:);
                    end
                    title("error = " + string(fval(j,2)) + " index = " + string(j))
                    xlabel('Position (cm)')
                    ylabel('Water Saturation')
                    legend(repmat(string(observed_fractional_flows),2,1),'Location','bestoutside')
                    grid('on')
                    ylim([0 1])
                end
            end
        end
        % waitbar(j/n, f, sprintf('Progress: %d %%', floor(j/n*100)));
        pause(1e-6)
    end
    % close(f)
end

function factional_flow = time_to_fractional_flow(time,model)
    table_size = size(model.experiment.schedule.table);
    factional_flow = [];
    for j = 1:length(time)
        for i = 1:table_size(1)-1
            if time(j) > model.experiment.schedule.table{i,1} && time(j) <= model.experiment.schedule.table{i+1,1}
                factional_flow = [factional_flow; model.experiment.schedule.table{i,4}];
            end
        end
    end
end
