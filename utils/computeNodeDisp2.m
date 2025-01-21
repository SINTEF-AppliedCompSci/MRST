function unode = computeNodeDisp2(ucell, tbls, mappings, varargin)
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

    opt = struct('useVirtual', false);
    opt = merge_options(opt, varargin{:});

    useVirtual = opt.useVirtual;
    
    celltbl        = tbls.celltbl;
    cellvectbl     = tbls.cellvectbl;
    nodetbl        = tbls.nodetbl;
    nodevectbl     = tbls.nodevectbl;
    cellnodetbl    = tbls.cellnodetbl;
    cellnodevectbl = tbls.cellnodevectbl;
    vectbl         = tbls.vectbl;
    
    map = TensorMap();
    map.fromTbl = cellnodetbl;
    map.toTbl = nodetbl;
    map.mergefds = {'nodes'};
    map = map.setup();
    
    ncpernode = map.eval(ones(cellnodetbl.num, 1));
    coef = 1./ncpernode;
    
    map = TensorMap();
    map.fromTbl  = nodetbl;
    map.toTbl    = cellnodetbl;
    map.mergefds = {'nodes'};

    if useVirtual
        
        map.pivottbl = cellnodetbl;
        map.dispind1 = mappings.node_from_cellnode;
        map.dispind2 = (1 : cellnodetbl.num)';

        map.issetup = true;
        
    else
        map = map.setup();
    end
    
    coef = map.eval(coef);
    
    prod = TensorProd();
    prod.tbl1 = cellnodetbl;
    prod.tbl2 = cellvectbl;
    prod.tbl3 = nodevectbl;
    prod.reducefds = {'cells'};

    
    if useVirtual

        cellnodevectbl = tbls.cellnodevectbl;
        
        prod.pivottbl = cellnodevectbl;

        [vec, i] = ind2sub([vectbl.num, cellnodetbl.num], (1 : cellnodevectbl.num)');
        
        prod.dispind1 = i;
        prod.dispind2 = sub2ind([vectbl.num, celltbl.num], vec, mappings.cell_from_cellnode(i));
        prod.dispind3 = sub2ind([vectbl.num, nodetbl.num], vec, mappings.node_from_cellnode(i));
        
        prod.issetup = true;
        
    else
        
        prod = prod.setup();
        
    end
    
    
    unode = prod.eval(coef, ucell);
end
