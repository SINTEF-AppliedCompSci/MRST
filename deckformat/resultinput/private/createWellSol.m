function wellSol = createWellSol(rstrt, repStep, ijk2act)
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

[w, ih]   = getRestartWellInfo(rstrt, repStep);
unit = {'metric', 'field', 'lab'};
unit = unit{ih.unit};
u    = getUnits(unit);

% only include open wells
%horzcat(w.stat)>0
%wix  = find(horzcat(w.stat)>0);
%wsix = 0;
if numel(w)==0
    wellSol = [];
end
for k = 1 : numel(w)
    %wsix = wsix +1;
    % only include open cons
    %cix = w(k).cstat > 0;
    ws = struct('name',     w(k).name, ...
        'cells',    ijk2act(w(k).cijk), ...
        'type',     '', ...
        'val',      [], ...
        'WI',       convertFrom(w(k).cwi, u.trans), ...
        'compi',    [], ...
        'refDepth', convertFrom(w(k).depth, u.len), ...
        'dZ',       convertFrom(w(k).cdepth-w(k).depth, u.len), ...
        'dir',      w(k).cdir, ...
        'sign',     0, ...
        'status',   w(k).stat > 0, ...
        'cstatus',  w(k).cstat > 0, ...
        'qWs',      convertFrom(w(k).qWs, u.ql), ...
        'qOs',      convertFrom(w(k).qOs, u.ql), ...
        'qGs',      convertFrom(w(k).qGs, u.qg), ...
        'bhp',      convertFrom(w(k).bhp, u.p),  ...
        'resv',     convertFrom(w(k).qr,  u.qr), ...
        'flux',     convertFrom(w(k).cqr, u.qr), ...
        'press',    convertFrom(w(k).press, u.p), ...
        'cqs',      [convertFrom(w(k).cqs(:,2), u.ql), ...
                     convertFrom(w(k).cqs(:,1), u.ql), ...
        	         convertFrom(w(k).cqs(:,3), u.qg)] );
    % if open, deduce type, val, compi and sign
    if w(k).stat > 0
        types  = {'orat', 'wrat', 'grat', 'lrat',         'resv',  'unknown', 'bhp'};
        vals   = [ws.qOs,  ws.qWs, ws.qGs, ws.qOs+ws.qWs,  ws.resv, nan,       ws.bhp];
        compis = {[0 1 0], [1 0 0], [0 0 1]};
        if w(k).type == 1 % producer
            sign  = -1;
            compi = [0 1 0];
            if w(k).cntr <= numel(types) && w(k).cntr > 0
                type = types{w(k).cntr};
                val  = vals(w(k).cntr);
            else
                type = '';
                val  = nan;
            end
        else              % injector
            sign = 1;
            if w(k).type <= 4
                compi = compis{w(k).type-1};
            else
                compi = vals(1:3)/sum(vals(1:3));
            end
            if w(k).cntr <= numel(types) && w(k).cntr > 0
                val = vals(w(k).cntr);
                if w(k).cntr <= 3
                    type = 'rate';
                else
                    type = types{w(k).cntr};
                end
            else
                type = '';
                val  = nan;
            end
        end
        ws.type  = type;
        ws.val   = val;
        ws.compi = compi;
        ws.sign  = sign;
    end
    wellSol(k,1) = ws; %#ok
end
end

%% ------------------------------------------------------------------------
function u = getUnits(unit)
switch unit
    case 'metric'
        u.len = meter;
        u.p   = barsa;
        u.ql  = meter^3/day;
        u.qg  = meter^3/day;
        u.qr  = meter^3/day;
        u.trans = centi*poise * meter^3 / (day * barsa);
    case 'field'
        u.len = ft;
        u.p   = psia;
        u.ql  = stb/day;
        u.qg  = 1000*ft^3/day;
        u.qr  = stb/day;
        u.trans = centi*poise * stb / (day * psia);
    otherwise
        error(['Unit ', unit, ' not supported']);
end
end
