function sols = makeWellSolsConsistent(sols)
% check if sols are states or wellSols

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

if isfield(sols{1}, 'wellSol')
    wellSols = cellfun( @(x)x.wellSol, sols,'UniformOutput', false);
    stateInput = true;
else
    wellSols = sols;
    stateInput = false;
end
% make template wellSol-struct
tmp = wellSols{end};
zfields = {'qWs', 'qOs', 'qGs', 'bhp', 'resv'};
efields = {'type', 'val', 'sign', 'compi'};
for k = 1:numel(tmp)
    nc = numel(tmp(k).cells);
    for fn = 1:numel(zfields)
        tmp(k).(zfields{fn}) = 0;
    end
    for fn = 1:numel(efields)
        tmp(k).(efields{fn}) = 0;
    end
    tmp(k).status  = false;
    tmp(k).cstatus = false(nc, 1);
    tmp(k).flux    = zeros(nc, 1);
    tmp(k).cqs     = zeros(nc, 3);
end

nr = numel(wellSols);
ws = repmat({tmp}, [nr, 1]);
for k = 1:numel(ws)
    w  = wellSols{k};
    nm = {ws{k}.name};
    for k2 = 1:numel(w)
        ix = strcmp(w(k2).name, nm);
        if numel(ws{k}(ix).cells) == numel(w(k2).cells)
            ws{k}(ix) = w(k2);
        else % set each field
            for f = [zfields, efields]
                ws{k}(ix).(f{:}) = w(k2).(f{:});
            end
            ws{k}(ix).status = true;
            na = numel(w(k2).cells);
            ws{k}(ix).cstatus(1:na) = w(k2).cstatus;
            ws{k}(ix).flux(1:na)    = w(k2).flux;
            ws{k}(ix).cqs(1:na, :)  = w(k2).cqs;
        end
    end
end
% traverse forth and back to get sign-field for all steps, give warning if 
% sign can't be found
wsign = {ws{1}.sign};
for k = 1:nr
    ix  = [ws{k}.sign] ~= 0;
    [wsign{ix}] = deal(ws{k}(ix).sign); 
end

for k = nr:-1:1
    ix  = [ws{k}.sign] == 0;
    [ws{k}(ix).sign] = deal(wsign{ix});
    [wsign{~ix}]     = deal(ws{k}(~ix).sign);
end
    
if stateInput
    for k = 1:numel(sols)
        sols{k}.wellSol = ws{k};
    end
else
    sols = ws;
end

end



    
