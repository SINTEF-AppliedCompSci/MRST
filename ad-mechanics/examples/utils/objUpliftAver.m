function obj = objUpliftAver(model, states, schedule, topnode, varargin)
%
%
% SYNOPSIS:
%   function obj = objUpliftAver(model, states, schedule, topnode, varargin)
%
% DESCRIPTION: Compute a time average of the vertical displacement at the
% node index by topnode. This function is used in runAdjointExample
%
% PARAMETERS:
%   model    - poroelastic model (fully coupled type e.g. MechWaterModel or MechOilWaterModel)
%   states   - cell array of state, as returned by simulateScheduleAD
%   schedule - schedule which was used as argument in simulateScheduleAD
%   topnode  - index of the node where the uplift will be computed
%   varargin - 
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   computePartials  - if true, the derivatives are computed
%   exponent         - We use the l^p norm as an average where p=exponent
%
% RETURNS:
%   obj - 
%
% EXAMPLE: runAdjointExample
%
% SEE ALSO: computeUpliftForState
%


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