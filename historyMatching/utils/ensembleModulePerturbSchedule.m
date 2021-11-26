function schedule = ensembleModulePerturbSchedule(schedule)
    for n=2:max(schedule.step.control)
        schedule.control(n) = schedule.control(1);
        for i=1:numel(schedule.control(n).W)
            W = schedule.control(n).W(i);
            switch W.type
                case 'rate'
                    W.val = (.75 + .5*rand)*W.val;
                case 'bhp'
                    %if rand < 0.2
                    %    W.status = false;
                    % else
                        W.val = (.95 + 0.1*rand)*W.val;
                    %end
            end
            schedule.control(n).W(i) = W;
        end
    end


end

