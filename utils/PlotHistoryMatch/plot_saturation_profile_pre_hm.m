function plot_saturation_profile_pre_hm (model)
%
% DESCRIPTION: this function comes after the PlotObservation_pre_hm module
%              to also plot and check the saturation profiles before 
%              starting the history match
%
% SYNOPSIS:
%   plot_saturation_profile_pre_hm (model)
%
% PARAMETERS:
%   model - struct containing following fields:
%   - experiment
%
% RETURNS:
%    - plot of the experimentally measured saturation profiles
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
    if(isfield(model.experiment.observation,'satProfile')) && model.experiment.observation.satProfile.include
        observed_sat_profile_t = model.experiment.observation.satProfile.table(2:end,1);
        observed_fractional_flows = time_to_fractional_flow(observed_sat_profile_t,model);
        observed_sat_profile_location = model.experiment.observation.satProfile.table(1,2:end);
        observed_sat_profile = model.experiment.observation.satProfile.table(2:end,2:end);
        nexttile;
        p = plot(observed_sat_profile_location / Convert('cm'), observed_sat_profile,...
            '-','DisplayName','Observed Saturation Profile','LineWidth',2, 'MarkerSize',3);
%         legend(repmat(string( round(observed_sat_profile_t / 3600,2)),1,1),'Location','bestoutside')
        xlabel('Position (cm)'); title('Saturation Profile');
        ylabel('Water Saturation'); 
        ylim([0 1]); grid off; xlim([0 core_length])
        colors = linspace(0.8, 0 ,length(observed_sat_profile_t))';
        colors = repmat(colors,1,3);
        for i = 1:length(observed_sat_profile_t)
            p(i).Color = colors(i,:);
            p(i).DataTipTemplate.DataTipRows(1).Label = 'Position (cm)';
            p(i).DataTipTemplate.DataTipRows(2).Label = 'Water Saturation';
            row = dataTipTextRow('Time (hour)', repmat(observed_sat_profile_t (i)...
                / Convert('hr'),1,length(observed_sat_profile)));
            p(i).DataTipTemplate.DataTipRows(end+1) = row;
        end

        nexttile;
%         plot3(repmat(observed_sat_profile_location / Convert('cm'),height(observed_sat_profile),1),...
%             repmat(observed_sat_profile_t / Convert('hr'),1,length(observed_sat_profile)),...
%             observed_sat_profile)
        s = surf(observed_sat_profile_location / Convert('cm'), ...
            observed_sat_profile_t / Convert('hr'),...
            observed_sat_profile); s.EdgeColor = 'none';
        xlabel('Position (cm)'); grid on 
        ylabel('Time (hour)'); title('Saturation profile')
        zlabel('Water Saturation'); zlim([0 1]);
        colorbar;
    end
end

function factional_flow = time_to_fractional_flow(time,model)
% calculates the fractional flow at each schedule row
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