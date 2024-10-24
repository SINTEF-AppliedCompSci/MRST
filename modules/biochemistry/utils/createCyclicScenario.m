function schedule = createCyclicScenario(rampupTime, nCycles, buildUpLength, restLength, injectionLength, productionLength, idleLength, W)
    % Create time steps
    TotalTime = buildUpLength + restLength + nCycles*(injectionLength+productionLength +idleLength);
    deltaT = rampupTimesteps(TotalTime, rampupTime, 0);
    schedule = simpleSchedule(deltaT);
    
    % Initialize control steps
    schedule.step.control = zeros(size(deltaT));
    tmp = cell(5,1);
    schedule.control = struct('W',tmp);
    % conversion
    buildUpSteps = round(buildUpLength/rampupTime);
    restSteps = round(restLength/rampupTime);
    injectionSteps = round(injectionLength/rampupTime);
    productionSteps = round(productionLength/rampupTime);
    idleSteps = round(idleLength/rampupTime);
    sumSTeps = sum(buildUpSteps+restSteps+(injectionSteps+productionSteps+idleSteps)*nCycles);
    if sumSTeps<length(schedule.step.val)
        buildUpSteps = buildUpSteps + (sumSTeps-length(schedule.step.val));
    end
    % Assign control 1 for the build-up phase
    schedule.step.control(1:buildUpSteps) = 1;
    schedule.control(1).W = W(1);
    % Assign cyclic controls for n cycles
    currentStep = buildUpSteps + 1;
    schedule.step.control(currentStep:currentStep+restSteps-1) = 2;
    schedule.control(2).W = W(2);
    currentStep =  currentStep + restSteps;
    for cycle = 1:nCycles
        schedule.step.control(currentStep:currentStep+injectionSteps-1) = 3;
        schedule.control(3).W = W(3);
        currentStep = currentStep + injectionSteps;
        schedule.step.control(currentStep:currentStep+productionSteps-1) = 4;
        schedule.control(4).W = W(4);
        currentStep = currentStep + productionSteps;
        schedule.step.control(currentStep:currentStep+idleSteps-1) = 5;
        schedule.control(5).W = W(2);
        currentStep = currentStep + idleSteps;
    end
    
    % Ensure the schedule does not exceed the total number of steps
    schedule.step.control = schedule.step.control(1:length(deltaT));
end