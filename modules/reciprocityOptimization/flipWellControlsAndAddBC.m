function newSchedule = flipWellControlsAndAddBC(schedule, state, bc)
% Modify an MRST schedule by flipping well control types and adding BC controls.
%
% Inputs:
%   schedule - Original MRST schedule structure containing well controls.
%   state    - MRST state structure with pressure and rate values for each well.
%   bc       - Boundary condition structure, provided as a cell array with 
%              boundary conditions for each timestep.
%
% Output:
%   newSchedule - Modified schedule with flipped well controls and BCs.
%
% The function creates a new schedule with:
%   - Flipped well controls for each well (rate <-> pressure)
%   - Boundary conditions applied to each step if `bc` is not empty.

% Initialize wells with the same structure as input
nSteps = numel(schedule.step.val);
nw = numel(schedule.control(1).W);  % Number of wells per control

% Loop through each step to adjust controls and add BCs
for step = 1:nSteps    
    % Initialize flipped well control
    W = schedule.control(schedule.step.control(step)).W;
    for w = 1:nw
        % current state values
        currentBhp = state{step}.wellSol(w).bhp;
        % Ensure any empty flow rates are set to zero to avoid sum errors
        currentqWs = getFieldOrDefault(state{step}.wellSol(w), 'qWs', 0);
        currentqOs = getFieldOrDefault(state{step}.wellSol(w), 'qOs', 0);
        currentqGs = getFieldOrDefault(state{step}.wellSol(w), 'qGs', 0);
        sumq = currentqWs +currentqOs + currentqGs;
        % Flip control type and set new value
        if (strcmp(W(w).type, 'rate') && (step>0))
            W(w).type = 'bhp';               % Switch to BHP control
            W(w).val = currentBhp;           % Set BHP based on state pressure
            W(w).lims  = [];           % Set BHP based on state pressure
            W(w).lims.wrat= currentqWs;
%             W(w).lims.grat= currentqGs;
%             W(w).lims.lrat= -Inf;
%             W(w).lims.bhp= currentBhp+1*barsa();
%             W(w).lims.thp= 0;
            W(w).compi = [abs(currentqWs./sumq) abs(currentqOs./sumq)];
            W(w).sign = 1;
        elseif (strcmp(W(w).type, 'bhp') && (step>0))
            W(w).type = 'rate';                         
            % Switch to rate control
            W(w).val = sumq; % Sum of flow rates
            W(w).lims  = [];           % Set BHP based on state pressure
%             W(w).lims.orat= currentqOs;
%             W(w).lims.wrat= currentqWs;
%             W(w).lims.lrat= -Inf;
%               W(w).lims.bhp= currentBhp*0.95;
%             W(w).lims.thp= 0;
            W(w).compi = [abs(currentqWs./sumq) abs(currentqOs./sumq)];
            W(w).sign = -1; 
        end
    end
    
    % Update control in new schedule
    if ~isempty(bc)
        newSchedule{step} = simpleSchedule(schedule.step.val(step),'W', W, 'bc', bc{step});
    else
        newSchedule{step} = simpleSchedule(schedule.step.val(step),'W', W);
    end
end

newSchedule = combineSchedules(newSchedule{:},'makeConsistent',false);
end

function value = getFieldOrDefault(structure, fieldName, defaultValue)
% Helper function to retrieve field value or default if field is empty or missing
    if isfield(structure, fieldName) && ~isempty(structure.(fieldName))
        value = structure.(fieldName);
    else
        value = defaultValue;
    end
end
