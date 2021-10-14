function schedule = perturbedSimpleSchedule(dt, varargin)
opt = struct('pressureFac', 0.05, ...
             'rateFac',     0.05, ...
             'perturbStep', []);
[opt, extra] = merge_options(opt, varargin{:});
schedule = simpleSchedule(dt, extra{:});

pfac = opt.pressureFac;
rfac = opt.rateFac;
p   = opt.perturbStep;
if isempty(p)
    p = (1:numel(dt))';
end
np = max(p);
assert(numel(p) == numel(dt))

schedule.step.control = p;
schedule.control = repmat(schedule.control, [np, 1]);
W  = schedule.control(1).W; 
nw = numel(W);
isPressure = arrayfun(@(w)strcmp(w.type, 'bhp'), W);

for kw = 1:nw
    if isPressure(kw)
        fac = pfac;
    else
        fac = rfac;
    end
    for kp = 1:np
        schedule.control(kp).W(kw).val = W(kw).val*(1 + fac*(rand-.5));
    end
end
end
