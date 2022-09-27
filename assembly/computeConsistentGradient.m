function g = computeConsistentGradient(G, eta, tbls, mappings, varargin)
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


    opt = struct('bcetazero', false);
    opt = merge_options(opt, varargin{:});
    bcetazero = opt.bcetazero;
    
    cellnodefacecents = computeNodeFaceCentroids(G, eta, tbls, 'bcetazero', opt.bcetazero);

    coltbl             = tbls.coltbl;
    cellnodetbl        = tbls.cellnodetbl;
    cellnodecoltbl     = tbls.cellnodecoltbl;
    cellnodefacetbl    = tbls.cellnodefacetbl;
    cellnodefacecoltbl = tbls.cellnodefacecoltbl;
    
    cellnode_from_cellnodeface = mappings.cellnode_from_cellnodeface;
    
    d_num    = coltbl.num;
    cn_num   = cellnodetbl.num;
    cnc_num  = cellnodecoltbl.num; 
    cnfc_num = cellnodefacecoltbl.num;
    cnf_num  = cellnodefacetbl.num;
    
    [c, i] = ind2sub([d_num, cnf_num], (1 : cnfc_num)');
    ind1 = i;
    ind2 = sub2ind([d_num, cn_num], c, cellnode_from_cellnodeface(i));

    assert(cnc_num == cnf_num, ['This implementation of mpsaw cannot handle ' ...
                        'this grid']);

    A = sparse(ind1, ind2, cellnodefacecents, cnc_num, cnc_num);

    opt.invertBlocks = 'mex';
    bi = blockInverter(opt);

    sz = repmat(coltbl.num, cellnodetbl.num, 1);
    invA = bi(A, sz);

    ind = sub2ind([cnf_num, cnc_num], ind2, ind1);
    
    g = invA(ind);
end
