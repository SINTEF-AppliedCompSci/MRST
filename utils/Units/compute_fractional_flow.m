function [sw, fw] = compute_fractional_flow (model)

    sw_avg = model.history_match.endofsched_swavg;
    factional_flow = model.experiment.schedule.procedure{:,5};
    
% older code    
%     time = model.history_match.endofsched_t;
%     table_size = size(model.experiment.schedule.procedure);
%     factional_flow = [];
%     for j = 1:length(time)
%         for i = 1:table_size(1)-1
%             if time(j) > model.experiment.schedule.procedure{i,2} && time(j) ...
%                     <= model.experiment.schedule.procedure{i+1,2}
%                 factional_flow = [factional_flow;...
%                     model.experiment.schedule.procedure{i,5}];
%             end
%         end
%     end

    sw = sw_avg;
    fw = factional_flow;
    
%     figure;
%     plot(sw_inlet, fw_inlet);
end