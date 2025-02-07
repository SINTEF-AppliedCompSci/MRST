function lastSteps = getCyclesLastSteps(schedule)
    %% Register the last step at the end of each injection, production, or idle period
    lastSteps.charge = [];
    lastSteps.discharge = [];
    lastSteps.shut = [];
    lastSteps.cushion = [];
    stageNames = fieldnames(lastSteps); 

    for i = 2:(length(schedule.step.control) - 1)
        currentwell = schedule.control(schedule.step.control(i)).W.name; 
        nextwell = schedule.control(schedule.step.control(i + 1)).W.name; 
        
        if ~strcmp(currentwell, nextwell)
            lastSteps.(currentwell)(end + 1) = i; 
        end
    end
    if ~isempty(currentwell)
        lastSteps.(currentwell)(end) = length(schedule.step.control); % Register the final step
    end

    % Display the last steps of each cycle for each stage
    disp('Last steps of each cycle for each stage:');
    for j = 1:length(stageNames)
        stage = stageNames{j};
        disp([stage, ': ', num2str(lastSteps.(stage))]);
    end
end
