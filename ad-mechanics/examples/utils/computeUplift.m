function obj = computeUplift(model, states, schedule, topnode, varargin)
% Compute the average of the vertical displacement at the top of the domain 
% This function is used in runAdjointExample

    opt = struct('tStep', [], 'ComputePartials', false);
    opt = merge_options(opt, varargin{:});
    
    
    numSteps = numel(schedule.step.val);
    lastStep = numSteps;
    tSteps = opt.tStep;
    if isempty(tSteps) 
        % do all
        tSteps = (1 : numSteps)';
    else
        numSteps = 1;
    end
    obj = repmat({[]}, numSteps, 1);

    for step = 1 : numSteps
        obj{step} = computeUpliftForState(model, states{tSteps(step)}, topnode, ...
                                          opt.ComputePartials);
        if tSteps(step) ~= lastStep
            obj{step} = double2ADI(0, obj{step});
        end
    end
    
end