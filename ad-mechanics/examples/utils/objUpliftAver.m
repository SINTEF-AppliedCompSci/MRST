function obj = objUpliftAver(model, states, schedule, topnode, varargin)
% Compute the average of the vertical displacement at the top of the domain
% This function is used in runAdjointExample

    opt = struct('tStep', [], ...
                 'ComputePartials', false, ...
                 'exponent', 100, ...
                 'normalizationConstant', 1);

    opt = merge_options(opt, varargin{:});

    p = opt.exponent;

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

    C = opt.normalizationConstant;
    for step = 1 : numSteps
        dt = schedule.step.val(tSteps(step));
        obj{step} = computeUpliftForState(model, states{tSteps(step)}, topnode, ...
                                          'ComputePartials', ...
                                          opt.ComputePartials);
        obj{step} = exp(p*log(eps + C*max(0, obj{step})))*dt;
    end

end