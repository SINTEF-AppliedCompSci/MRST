function unode = computeNodeDisp(ucell, tbls)
%Undocumented Utility Function

%{
Copyright 2020 University of Bergen and SINTEF Digital, Mathematics & Cybernetics.

This file is part of the MPSA-W module for the MATLAB Reservoir Simulation Toolbox (MRST).

The MPSA-W module is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The MPSA-W module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with the MPSA-W module.  If not, see <http://www.gnu.org/licenses/>.
%}


    celltbl        = tbls.celltbl;
    cellcoltbl     = tbls.cellcoltbl;
    nodetbl        = tbls.nodetbl;
    nodecoltbl     = tbls.nodecoltbl;
    cellnodetbl    = tbls.cellnodetbl;
    cellnodecoltbl = tbls.cellnodecoltbl;
    
    map = TensorMap();
    map.fromTbl = cellnodetbl;
    map.toTbl = nodetbl;
    map.mergefds = {'nodes'};
    map = map.setup();
    
    ncpernode = map.eval(ones(cellnodetbl.num, 1));
    coef = 1./ncpernode;
    
    map = TensorMap();
    map.fromTbl = nodetbl;
    map.toTbl = cellnodetbl;
    map.mergefds = {'nodes'};
    map = map.setup();
    
    coef = map.eval(coef);
    
    prod = TensorProd();
    prod.tbl1 = cellnodetbl;
    prod.tbl2 = cellcoltbl;
    prod.tbl3 = nodecoltbl;
    prod.reducefds = {'cells'};
    prod = prod.setup();
    
    unode = prod.eval(coef, ucell);
end
