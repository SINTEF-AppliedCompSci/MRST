function schedule = makeStorageSchedule(W, chargeIx, varargin)

    month = year/12;
    opt = struct( ...
        'time'            , [4, 4, 4]*month            , ...
        'dt'              , [7, 7, 7]*day              , ...
        'temperature'     , convertFromCelcius([90,10]), ...
        'rate'            , [50, 50]*litre/second      , ...
        'reverseDischarge', true                       , ...
        'numCycles'       , 1                          , ...
        'bhpLim'          , [1, 15]*barsa                ...
    );
    opt = merge_options(opt, varargin{:});
    
    if ~islogical(chargeIx)
        tmp           = chargeIx;
        chargeIx      = false(numel(W), 1);
        chargeIx(tmp) = true;
    end
    
    [WCharge   , dtCharge   ] = makeChargeStage(W, chargeIx, opt);
    [WDischarge, dtDischarge] = makeDischargeStage(W, chargeIx, opt);
    
    schedule = simpleSchedule(1);
    schedule.control(1).W = WCharge;
    schedule.control(2).W = WDischarge;
    
    nCharge    = numel(dtCharge);
    nDischarge = numel(dtDischarge);
    schedule.step.val = repmat([dtCharge; dtDischarge], opt.numCycles, 1);
    
    cno = rldecode([1,2], [nCharge, nDischarge], 2)';
    schedule.step.control = repmat(cno, opt.numCycles, 1);
                 
end

%-------------------------------------------------------------------------%
function [W, dt] = makeChargeStage(W, chargeIx, opt)
    
    injectors = chargeIx;
    [W(injectors).type] = deal('rate');
    [W(injectors).val ] = deal(opt.rate(1));
    [W(injectors).sign] = deal(1);

    [W(~injectors).type] = deal('bhp');
    [W(~injectors).val ] = deal(opt.bhpLim(1));
    [W(~injectors).sign] = deal(-1);
    
    [W.T] = deal(opt.temperature(1));
    
    dt = rampupTimesteps(opt.time(1), opt.dt(1));
    
end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [W, dt] = makeDischargeStage(W, chargeIx, opt)

    if opt.reverseDischarge
        injectors = ~chargeIx;
    else
        injectors = chargeIx;
    end
    
    [W(injectors).type] = deal('rate');
    [W(injectors).val ] = deal(opt.rate(2));
    [W(injectors).sign] = deal(1);

    [W(~injectors).type] = deal('bhp');
    [W(~injectors).val ] = deal(opt.bhpLim(1));
    [W(~injectors).sign] = deal(-1);
    
    [W.T] = deal(opt.temperature(2));
    
    dt = rampupTimesteps(opt.time(1), opt.dt(1));

end
%-------------------------------------------------------------------------%