function schedule = makeWAGschedule(W, nCycles, varargin)
% Make schedule for water-alternating gas (WAG) injection.
%
% SYNOPSIS:
%   schedule = makeWAGschedule(W, nCycles)
%
% DESCRIPTION:
%   This function makes a WAG schedule based on a well structure, to be
%   used with the BlackOilSolventModel.
%
% REQUIRED PARAMETERS:
%   W         - Well structure, properly initialized with e.g. addWell().
%               The function assumes this is initialized with desired
%               rates/bhp.
%   nCycles   - Number of cycles.
%
% OPTIONAL PARAMETERS:
%   'time'      - Total duration of the WAG injection period. Defaults to 1
%                 year.
%   'dt'          Target time step size. Defaluts to 30 days.
%   'gas_end'   - Duration of gas injection, 0 < gas_end < 1. The gas
%                 injection period of each cycle will last for
%                 time/ncycle*gas_end. Defaults to 0.5.
%   'useRampUp' - Use rampup each time the well control changes to ease the
%                 nonlinear solution process. Defaults to false.
%
% RETURNS:
%   schedule - WAG injection schedule.
%
% SEE ALSO:
%   simpleSchedule, BlackOilSolventModel

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

opt = struct('time'     , 1*year, ...
             'dt'       , 30*day, ...
             'gas_end'  , 0.5   , ...
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

