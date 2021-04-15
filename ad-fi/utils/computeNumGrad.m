function grad = computeNumGrad(state, G, rock, system, schedule, obj, varargin)
% Compute numerical gradient w.r.t wells

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

opt = struct('Verbose', mrstVerbose, 'scaling', []);
opt = merge_options(opt, varargin{:});

if ~isempty(opt.scaling)
    scalFacs = opt.scaling;
else
    scalFacs.rate = 1; scalFacs.pressure = 1;
end

%initial run:
wellSols = runScheduleADI(state, G, rock, system, schedule, 'Verbose', opt.Verbose);
val0 = sum( cell2mat(obj(wellSols)));



schedule0 = schedule;
grad = cell(1, numel(schedule.control));
%perturbed runs
for cn = 1:numel(schedule.control)
    inj   = schedule0.control(cn).WCONINJE;
    nInj  = size(inj,1);

    prod  = schedule0.control(cn).WCONPROD;
    nProd = size(prod, 1);

    nWell = nInj + nProd;
    grad{cn} = zeros(nWell, 1);


    dispif(opt.Verbose, 'Solving for control %d of %d, with %d injectors and %d producers\n', cn, numel(schedule.control), nInj, nProd);
    for k = 1:nInj
        dispif(opt.Verbose, 'Processing injector %d of %d\n', k, nInj);
        isBHP = strcmp(inj{k,4}, 'BHP');
        if isBHP
            e = scalFacs.pressure*1e-7;
        else
            e = scalFacs.rate*1e-7;
        end
%         -ws(k)*
        schedule = schedule0;
        upd = zeroVals(schedule);
        upd{cn}(k) = e;
        schedule = updateSchedule(schedule, upd, 'mode', 'add');%, 'updateInx',  cellfun(@(x) x~= 0, upd, 'UniformOutput', false));

        wellSols = runScheduleADI(state, G, rock, system, schedule, 'writeOutput', false, 'Verbose', opt.Verbose);
        valk = sum( cell2mat(obj(wellSols)));
        grad{cn}(k) = (valk-val0)/e;
    end

    for k = 1:nProd
        dispif(opt.Verbose, 'Processing producer %d of %d\n', k, nProd);
        isBHP = strcmp(prod{k,3}, 'BHP');
        if isBHP
            ws = 1;
            e = scalFacs.pressure*1e-7;
        else
            ws = -1;
            e = scalFacs.rate*1e-7;
        end
%         -ws(nInj+k)
        schedule = schedule0;
        upd = zeroVals(schedule);
        upd{cn}(nInj+k) = e;
        schedule = updateSchedule(schedule, upd, 'mode', 'add');%, 'updateInx',  cellfun(@(x) x~= 0, upd, 'UniformOutput', false));
        wellSols = runScheduleADI(state, G, rock, system, schedule, 'writeOutput', false, 'Verbose', opt.Verbose);
        valk = sum( cell2mat(obj(wellSols)));
        grad{cn}(nInj+k) = ws*(valk-val0)/e;
    end
end
end

function zz = zeroVals(schedule)
zz = cell(1, numel(schedule.control));
for k = 1:numel(zz)
    zz{k} = zeros( size(schedule.control(k).WCONINJE, 1) + ...
                   size(schedule.control(k).WCONPROD, 1) , 1);
end
end


