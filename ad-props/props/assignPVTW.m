function f = assignPVTW(f, pvtw, reg)
    [f.cW, f.muWr, f.bW, f.muW] = getFunctions(pvtw, reg);
end

function [cW, muWr, bW_fn, muW_fn] = getFunctions(pvtw, reg)
    assert(size(pvtw, 1) == reg.pvt);
    cW   = pvtw(:, 3)';
    muWr = pvtw(:, 4)';
    [bW_fn, muW_fn] = deal(cell(1, reg.pvt));
    
    for pvtnum = 1:size(pvtw, 1)
        pwr  = pvtw(pvtnum, 1); % ref pres
        bwr  = pvtw(pvtnum, 2); % ref fvf
        vbw  = pvtw(pvtnum,5); % viscosibility
        muwr = muWr(pvtnum);
        bW_fn{pvtnum}= getFunction(@(pw) bW(pw, cW(pvtnum), bwr, pwr), reg);
        if vbw > 0
            muW_fn{pvtnum} = getFunction(@(pw) muW(pw, pwr, muwr, vbw), reg);
        else
            muW_fn{pvtnum} = @(pw) repmat(muwr, numelValue(pw), 1);
        end
    end
end

function f = getFunction(fn, reg)
    ps = reg.prange;
    if isempty(ps)
        f = @(p) fn(p);
    else
        fs = fn(ps);
        f = @(p) reg.interp1d_uniform(ps, fs, p);
    end
end

function v = bW(pw, cw, bwr, pwr)
    X = cw.*(pw-pwr);
    v = (1+X+X.^2/2)./bwr;
end

function v = muW(pw, pwr, muwr, vbw)
    Y = -vbw.*(pw-pwr);
    v = muwr./(1+Y+Y.^2/2);
end


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
