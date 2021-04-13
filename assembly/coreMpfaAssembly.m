function [matrices, bcvals, extra] = coreMpfaAssembly(G, K, bcdirichlet, tbls, mappings, opts)
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

    
    bcetazero = opts.bcetazero;
    eta = opts.eta;
    dooptimize = opts.dooptimize;
    
    cellcolrowtbl         = tbls.cellcolrowtbl;
    cellnodecolrowtbl     = tbls.cellnodecolrowtbl;
    cellnodeface2coltbl   = tbls.cellnodeface2coltbl;
    cellnodeface2tbl      = tbls.cellnodeface2tbl;
    cellnodefacecolrowtbl = tbls.cellnodefacecolrowtbl;
    cellnodefacecoltbl    = tbls.cellnodefacecoltbl;
    cellnodefacetbl       = tbls.cellnodefacetbl;
    celltbl               = tbls.celltbl;
    coltbl                = tbls.coltbl;
    nodeface2tbl          = tbls.nodeface2tbl;
    nodefacecoltbl        = tbls.nodefacecoltbl;
    nodefacetbl           = tbls.nodefacetbl;
    
    if dooptimize
        % fetch the index mappings to set explictly the tensor products or tensor mappings
        cell_from_cellnodeface     = mappings.cell_from_cellnodeface;
        nodeface_from_cellnodeface = mappings.nodeface_from_cellnodeface;
        cellnodeface_1_from_cellnodeface2 = mappings.cellnodeface_1_from_cellnodeface2;
        cellnodeface_2_from_cellnodeface2 = mappings.cellnodeface_2_from_cellnodeface2;
        nodeface_1_from_nodeface2 = mappings.nodeface_1_from_nodeface2;
        nodeface_2_from_nodeface2 = mappings.nodeface_2_from_nodeface2;
    end
    
    % Some shortcuts
    c_num     = celltbl.num;
    cnf_num   = cellnodefacetbl.num;
    nf_num    = nodefacetbl.num;
    cnfcr_num = cellnodefacecolrowtbl.num;
    d_num     = coltbl.num;

    % g belongs to cellnodefacecoltbl;
    g = computeConsistentGradient(G, eta, tbls, mappings, 'bcetazero', bcetazero);

    % facetNormals belongs to cellnodefacecoltbl;
    normals = computeFacetNormals(G, cellnodefacetbl);

    % K belongs to cellcolrowtbl

    prod = TensorProd();
    prod.tbl1 = cellcolrowtbl;
    prod.tbl2 = cellnodefacecoltbl;
    prod.tbl3 = cellnodefacecoltbl;
    prod.replacefds2 = {{'coldim', 'rowdim'}};
    prod.mergefds = {'cells'};
    prod.reducefds = {'rowdim'};
    
    if dooptimize
        prod.pivottbl = cellnodefacecolrowtbl;
        [r, c, i] = ind2sub([d_num, d_num, cnf_num], (1 : cnfcr_num)');
        prod.dispind1 = sub2ind([d_num, d_num, c_num], r, c, cell_from_cellnodeface(i));
        prod.dispind2 = sub2ind([d_num, cnf_num], r, i);
        prod.dispind3 = sub2ind([d_num, cnf_num], c, i);
        prod.issetup = true;
    else
        prod = prod.setup();
    end
    
    Kg = prod.eval(K, g);
    
    prod = TensorProd();
    prod.tbl1 = cellnodefacecoltbl;
    prod.tbl2 = cellnodefacecoltbl;
    prod.tbl3 = cellnodeface2tbl;
    prod.replacefds1 = {{'faces', 'faces1'}};
    prod.replacefds2 = {{'faces', 'faces2'}};
    prod.mergefds = {'cells', 'nodes'};
    prod.reducefds = {'coldim'};
    
    if dooptimize
        cnf2_num = cellnodeface2tbl.num;
        cnf2c_num = cellnodeface2coltbl.num;
        prod.pivottbl = cellnodeface2coltbl;
        [c, i] = ind2sub([d_num, cnf2_num], (1 : cnf2c_num)');
        prod.dispind1 = sub2ind([d_num, cnf_num], c, cellnodeface_1_from_cellnodeface2(i));
        prod.dispind2 = sub2ind([d_num, cnf_num], c, cellnodeface_2_from_cellnodeface2(i));
        prod.dispind3 = i;
        prod.issetup = true;
    else
        prod = prod.setup();
    end
    
    nKg = prod.eval(normals, Kg);
    
    
    % Setup A11 matrix (facenode dof -> facenode dof)
    
    map = TensorMap();
    map.fromTbl = cellnodeface2tbl;
    map.toTbl = nodeface2tbl;
    map.mergefds = {'nodes', 'faces1', 'faces2'};
    map = map.setup(); % not optimized (use generic setup function)
    
    A11 = map.eval(nKg);
    
    prod = TensorProd();
    prod.tbl1 = nodeface2tbl;
    prod.tbl2 = nodefacetbl;
    prod.tbl3 = nodefacetbl;
    prod.replacefds1 = {{'faces1', 'faces'}};
    prod.replacefds2 = {{'faces', 'faces2'}};
    prod.mergefds = {'nodes'};
    prod.reducefds = {'faces2'};
    
    if dooptimize
        prod.pivottbl = nodeface2tbl;
        i = (1 : nodeface2tbl.num)';
        prod.dispind1 = i;
        prod.dispind2 = nodeface_2_from_nodeface2(i);
        prod.dispind3 = nodeface_1_from_nodeface2(i);
        prod.issetup = true;
    else
        prod = prod.setup();
    end
    
    A11_T = SparseTensor();
    A11_T = A11_T.setFromTensorProd(A11, prod);
    A11 = A11_T.getMatrix();
    
    [~, sz] = rlencode(nodefacetbl.get('nodes'), 1);
    opt.invertBlocks = 'mex';
    bi = blockInverter(opt);
    invA11 = bi(A11, sz);

    
    % Setup A12 matrix (cell dof -> facenode dof)
   
    map = TensorMap();
    map.fromTbl = cellnodeface2tbl;
    map.toTbl = cellnodefacetbl;
    map.replaceFromTblfds = {{'faces1', 'faces'}};
    map.mergefds = {'cells', 'nodes', 'faces'};
    
    if dooptimize
        map.pivottbl = cellnodeface2tbl;
        cnf2_num = cellnodeface2tbl.num;
        i = (1 : cnf2_num)';
        map.dispind1 = i;
        map.dispind2 = cellnodeface_1_from_cellnodeface2(i);
        map.issetup = true;
    else
        map = map.setup();
    end
    
    % beware minus sign here
    A12 = - map.eval(nKg);
    
    prod = TensorProd();
    prod.tbl1 = cellnodefacetbl;
    prod.tbl2 = celltbl;
    prod.tbl3 = nodefacetbl;
    prod.reducefds = {'cells'};
    
    if dooptimize
        prod.pivottbl = cellnodefacetbl;
        i = (1 : cellnodefacetbl.num)';
        prod.dispind1 = i;
        prod.dispind2 = cell_from_cellnodeface(i);
        prod.dispind3 = nodeface_from_cellnodeface(i);
        prod.issetup = true;
    else
        prod = prod.setup();
    end
    
    A12_T = SparseTensor();
    A12_T = A12_T.setFromTensorProd(A12, prod);
    A12 = A12_T.getMatrix();

    
    % Setup A21 matrix (facenode dof -> cell dof)

    map = TensorMap();
    map.fromTbl = cellnodeface2tbl;
    map.toTbl = cellnodefacetbl;
    map.replaceFromTblfds = {{'faces2', 'faces'}};
    map.mergefds = {'cells', 'nodes', 'faces'};
    
    if dooptimize
        map.pivottbl = cellnodeface2tbl;
        i = (1 : cellnodeface2tbl.num)';
        map.dispind1 = i;
        map.dispind2 = cellnodeface_2_from_cellnodeface2(i);
        map.issetup = true;
    else
        map = map.setup();
    end
    
    % beware minus sign here
    A21 = -map.eval(nKg);
    
    prod = TensorProd();
    prod.tbl1 = cellnodefacetbl;
    prod.tbl2 = nodefacetbl;
    prod.tbl3 = celltbl;
    prod.reducefds = {'faces', 'nodes'};
    prod = prod.setup();
    
    if dooptimize
        prod.pivottbl = cellnodefacetbl;
        i = (1 : cellnodefacetbl.num)';
        prod.dispind1 = i;
        prod.dispind2 = nodeface_from_cellnodeface(i);
        prod.dispind3 = cell_from_cellnodeface(i);
    else
        prod = prod.setup();
    end
    
    A21_T = SparseTensor();
    A21_T = A21_T.setFromTensorProd(A21, prod);
    A21 = A21_T.getMatrix();
    
    
    % Setup A22 matrix (cell dof -> cell dof)

    map = TensorMap();
    map.fromTbl = cellnodeface2tbl;
    map.toTbl = celltbl;
    map.mergefds = {'cells'};

    if dooptimize
        map.pivottbl = cellnodeface2tbl;
        i = (1 : cellnodeface2tbl.num);
        map.dispind1 = i;
        map.dispind2 = cell_from_cellnodeface(cellnodeface_1_from_cellnodeface2(i));
        map.issetup = true;
    else
        map = map.setup();
    end
    
    A22 = map.eval(nKg);
    
    prod = TensorProd();
    prod.tbl1 = celltbl;
    prod.tbl2 = celltbl;
    prod.tbl3 = celltbl;
    prod.mergefds = {'cells'};
    
    if dooptimize
        prod.pivottbl = celltbl;
        i = (1 : celltbl.num);
        prod.dispind1 = i;
        prod.dispind2 = i;
        prod.dispind3 = i;
        prod.issetup = true;
    else
        prod = prod.setup();
    end
    
    A22_T = SparseTensor();
    A22_T = A22_T.setFromTensorProd(A22, prod);
    A22 = A22_T.getMatrix();
    
    if ~isempty(bcdirichlet)
        [D, bcvals] = setupMpfaNodeFaceBc(bcdirichlet, tbls);
        coef = max(abs(A11), [], 'all');
        D = coef*D;
        bcvals = coef*bcvals;
    else
        D = ones(size(A11, 1), 0);
        bcvals = [];
    end

    matrices = struct('invA11', invA11, ...
                      'A11', A11, ...
                      'A12', A12, ...
                      'A21', A21, ...
                      'A22', A22, ...
                      'D', D);
    
    extra.nKg = nKg;
    
end
