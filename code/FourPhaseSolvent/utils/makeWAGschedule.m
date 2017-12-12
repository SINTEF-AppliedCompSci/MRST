function schedule = makeWAGschedule(W, nCycles, varargin)

opt = struct('gas_end'  , 0.5   , ...
             'time'     , 1*year, ...
             'dt'       , 30*day, ...
             'useRampup', false );
         
opt = merge_options(opt, varargin{:});

%%

time = opt.time;
dt      = opt.dt;

tCycle  = time/nCycles;
gas_end = opt.gas_end;

if opt.useRampup
    dtG1 = rampupTimesteps(gas_end*tCycle    , dt*gas_end, 8);
    dtW1 = rampupTimesteps((1-gas_end)*tCycle, dt, 4);
    dtG  = rampupTimesteps(gas_end*tCycle    , dt, 2);
    dtW  = rampupTimesteps((1-gas_end)*tCycle, dt, 2);
else
    [dtG1, dtG] = deal(rampupTimesteps(gas_end*tCycle    , dt, 0));
    [dtW1, dtW] = deal(rampupTimesteps((1-gas_end)*tCycle, dt, 0));
end

step.val     = [dtG1; dtW1; repmat([dtG; dtW], nCycles-1, 1)];
step.control = [1*ones(numel(dtG1),1); 2*ones(numel(dtW1),1); ...
                repmat([1*ones(numel(dtG),1); 2*ones(numel(dtW),1)], nCycles-1, 1)];

[WG, WW] = deal(W);
for wNo = 1:numel(W)
    WG(wNo).compi = [0, 0, 0, 1];
    WW(wNo).compi = [1, 0, 0, 0];
end

control(1).W = WG;
control(2).W = WW;

schedule.control = control;
schedule.step    = step;

end

