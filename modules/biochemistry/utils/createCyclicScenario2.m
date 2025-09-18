function [schedule, TotalTime, buildUpSteps, restSteps, injectionSteps, idleSteps, productionSteps, idle1Steps] = ...
    createCyclicScenario2(rampupTime, nCycles, buildUpLength, restLength, injectionLength, ...
                         idleLength, productionLength, idleFinalLength, W)
% Create a cyclic injection schedule for reservoir simulation
%
% SYNOPSIS:
%   [schedule, TotalTime, buildUpSteps, ...] = createCyclicScenario2(...)
%
% PARAMETERS:
%   rampupTime    - Initial timestep size (seconds)
%   nCycles       - Number of injection cycles
%   buildUpLength - Duration of initial buildup phase (seconds)
%   restLength    - Duration of rest period after buildup (seconds)
%   injectionLength - Duration of injection phase (seconds)
%   idleLength    - Duration of idle period (seconds)
%   productionLength - Duration of production phase (seconds)
%   idleFinalLength   - Duration of final idle period (seconds)
%   W             - Well structures for different phases
%
% RETURNS:
%   schedule      - MRST schedule structure
%   TotalTime     - Total simulation time (seconds)
%   buildUpSteps  - Number of timesteps in buildup phase
%   ...           - Number of timesteps in other phases
%
% DESCRIPTION:
%   Creates a schedule with cyclic injection/production periods for reservoir
%   simulation. The schedule consists of:
%   1. Initial buildup phase
%   2. Rest period
%   3. Multiple cycles of:
%      a. Injection
%      b. Idle
%      c. Production
%      d. Idle

    % Calculate total simulation time
    TotalTime = buildUpLength + restLength + ...
                nCycles*(injectionLength + idleLength + productionLength + idleFinalLength);
    
    % Create timesteps
    deltaT = rampupTimesteps(TotalTime, rampupTime, 0);
    schedule = simpleSchedule(deltaT);
    
    % Initialize control structure
    schedule.step.control = zeros(size(deltaT));
    schedule.control = struct('W', cell(6,1));
    
    % Calculate number of steps for each phase
    buildUpSteps    = round(buildUpLength / rampupTime);
    restSteps       = round(restLength / rampupTime);
    injectionSteps  = round(injectionLength / rampupTime);
    idleSteps       = round(idleLength / rampupTime);
    productionSteps = round(productionLength / rampupTime);
    idle1Steps      = round(idleFinalLength / rampupTime);
    
    % Adjust steps if total doesn't match
    sumSteps = buildUpSteps + restSteps + ...
               nCycles*(injectionSteps + idleSteps + productionSteps + idle1Steps);
    
    if sumSteps < length(schedule.step.val)
        buildUpSteps = buildUpSteps + (length(schedule.step.val) - sumSteps);
    end
    
    % Assign controls to each phase
    
    % 1. Build-up phase
    schedule.step.control(1:buildUpSteps) = 1;
    schedule.control(1).W = W(1);
    
    % 2. Rest period
    currentStep = buildUpSteps + 1;
    schedule.step.control(currentStep:currentStep + restSteps - 1) = 2;
    schedule.control(2).W = W(2);
    currentStep = currentStep + restSteps;
    
    % 3. Cyclic phases
    for cycle = 1:nCycles
        % a. Injection
        schedule.step.control(currentStep:currentStep + injectionSteps - 1) = 3;
        schedule.control(3).W = W(3);
        currentStep = currentStep + injectionSteps;
        
        % b. Idle
        schedule.step.control(currentStep:currentStep + idleSteps - 1) = 4;
        schedule.control(4).W = W(2);
        currentStep = currentStep + idleSteps;
       
        % c. Production
        schedule.step.control(currentStep:currentStep + productionSteps - 1) = 5;
        schedule.control(5).W = W(4);
        currentStep = currentStep + productionSteps;
        
        % d. Idle
        schedule.step.control(currentStep:currentStep + idle1Steps - 1) = 6;
        schedule.control(6).W = W(2);
        currentStep = currentStep + idle1Steps;
    end
    
    % Ensure we don't exceed allocated steps
    schedule.step.control = schedule.step.control(1:length(deltaT));
end

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}