function schedule = makeStorageSchedule(W, chargeIx, varargin)
% Undocumented Helper Function

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
    
    limits = defineLimits(opt);
    
    [WCharge   , dtCharge   ] = makeChargeStage(W, chargeIx, limits, opt);
    [WDischarge, dtDischarge] = makeDischargeStage(W, chargeIx, limits, opt);
    
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
function limits = defineLimits(opt)

    limits = struct('bhp', opt.bhpLim(2));

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [W, dt] = makeChargeStage(W, chargeIx, limits, opt)
    
    injectors = chargeIx;
    [W(injectors).type] = deal('rate');
    [W(injectors).val ] = deal(opt.rate(1)/nnz(injectors));
    [W(injectors).sign] = deal(1);
    [W(injectors).lims] = deal(limits);

    [W(~injectors).type] = deal('bhp');
    [W(~injectors).val ] = deal(opt.bhpLim(1));
    [W(~injectors).sign] = deal(-1);
    
    [W.T] = deal(opt.temperature(1));
    
    dt = rampupTimesteps(opt.time(1), opt.dt(1));
    
end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [W, dt] = makeDischargeStage(W, chargeIx, limits, opt)

    if opt.reverseDischarge
        injectors = ~chargeIx;
    else
        injectors = chargeIx;
    end
    
    [W(injectors).type] = deal('rate');
    [W(injectors).val ] = deal(opt.rate(2)/nnz(injectors));
    [W(injectors).sign] = deal(1);
    [W(injectors).lims] = deal(limits);

    [W(~injectors).type] = deal('bhp');
    [W(~injectors).val ] = deal(opt.bhpLim(1));
    [W(~injectors).sign] = deal(-1);
    
    [W.T] = deal(opt.temperature(2));
    
    dt = rampupTimesteps(opt.time(1), opt.dt(1));
end
