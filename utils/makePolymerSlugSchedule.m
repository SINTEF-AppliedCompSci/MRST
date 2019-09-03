function schedule = makePolymerSlugSchedule(W, fluid, varargin)

    opt = struct('dt'     , 30*day        , ...
                 'nRampup', 0             , ...
                 'wTime'  , [1,1]*0.5*year, ...
                 'pTime'  , 0.5*year      );
    opt = merge_options(opt, varargin{:});

    dtw1        = rampupTimesteps2(opt.wTime(1), opt.dt, opt.nRampup);
    schedule_w1 = simpleSchedule(dtw1, 'W', W);
    dtw2        = rampupTimesteps2(opt.wTime(2), opt.dt, opt.nRampup);
    schedule_w2 = simpleSchedule(dtw2, 'W', W);
    dtp         = rampupTimesteps2(opt.pTime, opt.dt, opt.nRampup);
    schedule_p  = simpleSchedule(dtp, 'W', W);
    schedule_p.step.control = schedule_p.step.control + 1;
    
    schedule              = schedule_w1;
    schedule.step.val     = [schedule_w1.step.val; ...
                             schedule_p.step.val; ...
                             schedule_w2.step.val];
    schedule.step.control = [schedule_w1.step.control; ...
                             schedule_p.step.control; ...
                             schedule_w2.step.control];
    schedule.control(2)       = schedule.control(1);
    [schedule.control(2).W.c] = deal(fluid.cmax);

end