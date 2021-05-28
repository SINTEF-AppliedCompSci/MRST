function schedule = updateSchedule(schedule, vals, varargin)
%Undocumented Utility Function

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

opt = struct('mode', 'replace', 'updateInx', []);
opt = merge_options(opt, varargin{:});

for k = 1:numel(schedule.control)
    v = vals{k};
    if ~isempty(opt.updateInx)
        inx = opt.updateInx{k};
    else
        inx = true(numel(v),1);
    end

    inj    = schedule.control(k).WCONINJE;
    numInj = size(inj, 1);
    updInj = find(inx(1:numInj));
    numUpdInj = numel(updInj);
    vi     = v(1:numUpdInj);
    for ki = 1:numUpdInj
        inj(updInj(ki),:) = updateInjector(inj(updInj(ki),:), vi(ki), opt.mode);
    end
    schedule.control(k).WCONINJE = inj;

    prod    = schedule.control(k).WCONPROD;
    updProd = find(inx(numInj+1:end));
    numUpdProd = numel(updProd);
    vp      = v(numUpdInj+1:end);
    for kp = 1:numUpdProd
        prod(updProd(kp),:) = updateProducer(prod(updProd(kp),:), vp(kp), opt.mode);
    end
    schedule.control(k).WCONPROD = prod;
end
end

%--------------------------------------------------------------------------

function inj = updateInjector(inj, val, mode)
cntrMode = inj{4};
switch cntrMode
    case 'RATE'
        inx = 5;
    case 'RESV'
        inx = 6;
    case 'BHP'
        inx = 7;
    otherwise
        error(['Cannot handle control mode: ', cntrMode]);
end
if strcmp(mode, 'replace')
    inj{inx} = val;
elseif strcmp(mode, 'add')
    inj{inx} = inj{inx}+val;
end
end

function prod = updateProducer(prod, val, mode)
cntrMode = prod{3};
switch cntrMode
    case 'ORAT'
        inx = 4;
    case 'WRAT'
        inx = 5;
    case 'LRAT'
        inx = 6;
    case 'RESV'
        inx = 8;
    case 'BHP'
        inx = 9;
    otherwise
        error(['Cannot handle control mode: ', cntrMode]);
end
if strcmp(mode, 'replace')
    prod{inx} = val;
elseif strcmp(mode, 'add')
    prod{inx} = prod{inx}+val;
end
end


