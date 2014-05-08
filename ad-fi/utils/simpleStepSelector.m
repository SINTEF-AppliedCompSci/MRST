function [dt, dt_history] = simpleStepSelector(dt_history, dt_current, its, varargin)
    opt = struct('targetIts', 5, ...
        'dt_min', 0.001*day, ...
        'dt_max', inf, ...
        'stepModifier', 1.5);

    opt = merge_options(opt, varargin{:});


    if its == opt.targetIts || its == 0
        modifier = 1;
    elseif its > opt.targetIts
        modifier = 1/opt.stepModifier;
    else
        modifier = opt.stepModifier;
    end

    dt = modifier*dt_current;

    dt_history = vertcat(dt_history, [dt, its]);
    dt = max(dt, opt.dt_min);
    dt = min(dt, opt.dt_max);
end
