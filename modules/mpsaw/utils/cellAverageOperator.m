function aver = cellAverageOperator(tbls, mappings)
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


    cellnodecolrowtbl = tbls.cellnodecolrowtbl;
    cellcolrowtbl = tbls.cellcolrowtbl;
    cellnodetbl = tbls.cellnodetbl;
    celltbl = tbls.celltbl;
    coltbl = tbls.coltbl;
    
    cell_from_cellnode = mappings.cell_from_cellnode;
    
    % shortcuts
    d_num = coltbl.num;
    c_num = celltbl.num;
    cn_num = cellnodetbl.num;
    cncr_num = cellnodecolrowtbl.num;
    
    % Compute number of nodes per cells
    map = TensorMap();
    map.fromTbl = cellnodetbl;
    map.toTbl = celltbl;
    map.mergefds = {'cells'};
    map = map.setup();
    
    nnodepercell = map.eval(ones(cellnodetbl.num, 1));
    
    % Compute cell average stress
    prod = TensorProd();
    prod.tbl1 = celltbl;
    prod.tbl2 = cellnodecolrowtbl;
    prod.tbl3 = cellcolrowtbl;
    prod.mergefds = {'cells'};
    
    prod.pivottbl = cellnodecolrowtbl;
    [r, c, i] = ind2sub([d_num, d_num, cn_num], (1 : cncr_num)');
    prod.dispind1 = cell_from_cellnode(i);
    prod.dispind2 = (1 : cncr_num)';
    prod.dispind3 = sub2ind([d_num, d_num, c_num], r, c, cell_from_cellnode(i));
    prod.issetup = true;

    aver_T = SparseTensor();
    aver_T = aver_T.setFromTensorProd(1./nnodepercell, prod);

    aver = aver_T.getMatrix();
end
