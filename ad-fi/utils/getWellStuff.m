function [Tw, dzw, Rw, wc, perf2well, pInx, iInxW, iInxO] = getWellStuff(W)
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

if isempty(W)
    [Tw, dzw, Rw, wc, perf2well, pInx, iInxW] = deal([]);
    return
end

nPerf = cellfun(@numel, {W.cells})';
nPerfTot = sum(nPerf);
nw    = numel(W);
perf2well = rldecode((1:nw)', nPerf);

Rw = sparse((1:nPerfTot)', rldecode((1:nw)', nPerf), 1, nPerfTot, nw);

%------------------------------------------

Tw    = vertcat(W.WI);
dzw   = vertcat(W.dZ);
wc    = vertcat(W.cells);
inj   = vertcat(W.sign)==1;
compi = vertcat(W.compi);
iInx  = rldecode(inj, nPerf);
pInx  = ~iInx;
iInx   = find(iInx);
pInx   = find(pInx);
iInxW  = iInx(compi(perf2well(iInx),1)==1);
iInxO  = iInx(compi(perf2well(iInx),2)==1);
%------------------------------

end
