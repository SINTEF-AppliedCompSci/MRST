function [matrices, bcvals, extra] = coreMpsaAssembly(G, C, bc, nnodesperface, tbls, mappings, opts)
% Assembly of MPSA-weak
%
% Reference paper:
% Finite volume methods for elasticity with weak symmetry
% Keilegavlen, Eirik and Nordbotten, Jan Martin
% International Journal for Numerical Methods in Engineering
% 2017

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


    % the solution is given by the system
    %
    % A = [[A11, A12, -D];
    %      [A21, A22,  0];
    %      [D' , 0  ,  0]];
    %
    % u = [u (displacement at nodefacecoltbl);
    %      u (displacement at cellcoltbl);
    %      lagmult (forces in the linear directions at the boundary)];
    %
    % f = [extforce (force at nodefacecoltbl);
    %      force    (volumetric force at cellcoltbl);
    %      bcvals   (for the linear form at the boundary)];
    %
    % A*u = f
    %
    % Note: extforce is sparse and should only give contribution at facets
    % that are at the boundary
    %
    % By construction of the method, the matrix A11 is block-diagonal. Hence,
    % we invert it directly and reduce to a cell-centered scheme.
        
    bcetazero = opts.bcetazero;
    eta = opts.eta;
    
    coltbl      = tbls.coltbl;
    colrowtbl   = tbls.colrowtbl;
    col2row2tbl = tbls.col2row2tbl;
    celltbl     = tbls.celltbl;
    nodetbl     = tbls.nodetbl;
    cellnodetbl = tbls.cellnodetbl;
    nodefacetbl = tbls.nodefacetbl;
    cellcoltbl  = tbls.cellcoltbl;
    nodecoltbl  = tbls.nodecoltbl;
    
    nodefacecoltbl        = tbls.nodefacecoltbl;
    cellnodefacetbl       = tbls.cellnodefacetbl;
    cellnodecoltbl        = tbls.cellnodecoltbl;
    cellnodecolrowtbl     = tbls.cellnodecolrowtbl;
    cellnodefacecoltbl    = tbls.cellnodefacecoltbl;
    cellnodefacecolrowtbl = tbls.cellnodefacecolrowtbl;
    nodecolrowtbl         = tbls.nodecolrowtbl;
    cellcol2row2tbl       = tbls.cellcol2row2tbl;
    cellnodecol2row2tbl   = tbls.cellnodecol2row2tbl;

    cell_from_cellnode         = mappings.cell_from_cellnode;
    node_from_cellnode         = mappings.node_from_cellnode;
    cellnode_from_cellnodeface = mappings.cellnode_from_cellnodeface;
    nodeface_from_cellnodeface = mappings.nodeface_from_cellnodeface;
    
    % Some shortcuts
    c_num     = celltbl.num;
    n_num     = nodetbl.num;
    cnf_num   = cellnodefacetbl.num;
    cnfc_num  = cellnodefacecoltbl.num;
    cn_num    = cellnodetbl.num;
    cncr_num  = cellnodecolrowtbl.num;
    nf_num    = nodefacetbl.num;
    nfc_num   = nodefacecoltbl.num;
    cnfcr_num = cellnodefacecolrowtbl.num;
    d_num     = coltbl.num;
    cc2r2_num = cellcol2row2tbl.num; %shortcut
    c2r2_num  = col2row2tbl.num; %shortcut
        
    dim = coltbl.num;
    
    % Construction of tensor g (as defined in paper eq 4.1.2)
    
    g = computeConsistentGradient(G, eta, tbls, mappings, 'bcetazero', bcetazero);
    
    % Construction of the gradient operator
    %

    % Construction of gradnodeface_op : nodefacecoltbl -> cellnodecolrowtbl
    %
    % The nodefacecol part of the grad operator from nodefacecoltbl to cellnodecolrowtbl is obtained for any u in
    % nodefacecoltbl by using v = prod.eval(g, u) where prod is defined below and this is how we set the correspinding
    % tensor.
    %
    prod = TensorProd();
    prod.tbl1 = cellnodefacecoltbl;
    prod.tbl2 = nodefacecoltbl;
    prod.tbl3 = cellnodecolrowtbl;
    prod.replacefds2 = {'coldim', 'rowdim'};
    prod.reducefds   = {'faces'};
    prod.mergefds    = {'nodes'};

    prod.pivottbl = cellnodefacecolrowtbl;
    [r, c, i] = ind2sub([d_num, d_num, cnf_num], (1 : cnfcr_num)');
    
    prod.dispind1 = sub2ind([d_num, cnf_num], c, i);
    prod.dispind2 = sub2ind([d_num, cnf_num], r, nodeface_from_cellnodeface(i));
    prod.dispind3 = sub2ind([d_num, d_num, cn_num], r, c, cellnode_from_cellnodeface(i));
    prod.issetup = true;
    
    gradnodeface_T = SparseTensor('matlabsparse', true);
    gradnodeface_T = gradnodeface_T.setFromTensorProd(g, prod);

    % Construction of gradcell_T : cellcoltbl -> cellnodecolrowtbl
    %
    % The cellcol part of the grad operator from cellcoltbl to cellnodecolrowtbl is
    % obtained for any u in cellcoltbl by using v = prod.eval(greduced, u)
    % where greduced and prod are defined below
    %
    map = TensorMap();
    map.fromTbl = cellnodefacecoltbl;
    map.toTbl = cellnodecoltbl;
    map.mergefds = {'cells', 'nodes', 'coldim'};

    map.pivottbl = cellnodefacecoltbl;
    map.dispind1 = (1 : cnfc_num)';
    [c, i] = ind2sub([d_num, cnf_num], (1 : cnfc_num)');
    map.dispind2 = sub2ind([d_num, cn_num], c, cellnode_from_cellnodeface(i));
    map.issetup = true;

    % note minus sign
    greduced = - map.eval(g);
    
    prod = TensorProd();
    prod.tbl1 = cellnodecoltbl;
    prod.tbl2 = cellcoltbl;
    prod.tbl3 = cellnodecolrowtbl;
    prod.replacefds2 = {'coldim', 'rowdim'};
    prod.mergefds = {'cells'};

    prod.pivottbl = cellnodecolrowtbl;
    [r, c, i] = ind2sub([d_num, d_num, cn_num], (1 : cncr_num)');
    prod.dispind1 = sub2ind([d_num, cn_num], c, i);
    prod.dispind2 = sub2ind([d_num, c_num], r, cell_from_cellnode(i));
    prod.dispind3 = (1 : cncr_num);
    prod.issetup = true;

    % prod = prod.setup();
    
    gradcell_T = SparseTensor('matlabsparse', true);
    gradcell_T = gradcell_T.setFromTensorProd(greduced, prod);

    % Construction of the divergence operator
    %
    % setup the facet normals

    facetNormals = computeFacetNormals(G, cellnodefacetbl);
    
    % divnodeface_T : cellnodecolrowtbl -> nodefacecoltbl
    %
    % The nodefacecol part of the divergence operator from cellnodecolrowtbl to
    % nodefacecoltbl is obtained for any u in cellnodecolrowtbl by evaluating the
    % expression divnodeface_T.eval(d, u) where d and divnodeface_T are defined
    % below
    %
    d = facetNormals;
    prod = TensorProd();
    prod.tbl1 = cellnodefacecoltbl;
    prod.tbl2 = cellnodecolrowtbl;
    prod.replacefds1 = {'coldim', 'rowdim'};
    prod.replacefds2 = {'coldim', 'rowdim', 'interchange'};
    prod.reducefds = {'rowdim', 'cells'};
    prod.mergefds = {'nodes'};
    prod.tbl3 = nodefacecoltbl;

    prod.pivottbl = cellnodefacecolrowtbl;
    [r, c, i] = ind2sub([d_num, d_num, cnf_num], (1 : cnfcr_num)');
    prod.dispind1 = sub2ind([d_num, cnf_num], r, i);
    prod.dispind2 = sub2ind([d_num, d_num, cn_num], c, r, cellnode_from_cellnodeface(i));
    prod.dispind3 = sub2ind([d_num, nf_num], c, nodeface_from_cellnodeface(i));
    prod.issetup = true;

    divnodeface_T = SparseTensor('matlabsparse', true);
    divnodeface_T = divnodeface_T.setFromTensorProd(d, prod);

    % divcell_T : cellnodecolrowtbl -> cellcoltbl
    %
    % the cellcol part of the divergence operator from cellnodecolrowtbl to
    % cellcoltbl is obtained for any u in cellnodecolrowtbl by evaluating the
    % expression divcell_T.eval(dreduced, u) where dreduced and divcell_T
    % are defined below
    %

    fds = {'cells', 'nodes', 'coldim'};
    % note the minus sign below (see formula in paper)
    map = TensorMap();
    map.fromTbl = cellnodefacecoltbl;
    map.toTbl = cellnodecoltbl;
    map.mergefds = {'cells', 'nodes', 'coldim'};
    
    map.pivottbl = cellnodefacecoltbl;
    map.dispind1 = (1 : cnfc_num)';
    [c, i] = ind2sub([d_num, cnf_num], (1 : cnfc_num)');
    map.dispind2 = sub2ind([d_num, cn_num], c, cellnode_from_cellnodeface(i));
    map.issetup = true;

    dreduced = - map.eval(facetNormals);

    prod = TensorProd();
    prod.tbl1 = cellnodecoltbl;
    prod.tbl2 = cellnodecolrowtbl;
    prod.tbl3 = cellcoltbl;
    prod.replacefds1 = {'coldim', 'rowdim'};
    prod.replacefds2 = {'coldim', 'rowdim', 'interchange'};
    prod.reducefds   = {'rowdim', 'nodes'};
    prod.mergefds    = {'cells'};

    prod.pivottbl = cellnodecolrowtbl;
    [r, c, i] = ind2sub([d_num, d_num, cn_num], (1 : cncr_num)');
    prod.dispind1 = sub2ind([d_num, cn_num], r, i);
    prod.dispind2 = sub2ind([d_num, d_num, cn_num], c, r, i);
    prod.dispind3 = sub2ind([d_num, c_num], c, cell_from_cellnode(i));
    prod.issetup = true;
    
    divcell_T = SparseTensor('matlabsparse', true);
    divcell_T = divcell_T.setFromTensorProd(dreduced, prod);

    % Construction of transpose operator for matrices at nodes (that are
    % elements of nodecolrowtbl)
    %
    %  trans_T: nodecolrowtbl -> nodecolrowtbl

    clear symcol2row2tbl;
    symcol2row2tbl.coldim2 = colrowtbl.get('coldim');
    symcol2row2tbl.rowdim2 = colrowtbl.get('rowdim');
    symcol2row2tbl.coldim1 = colrowtbl.get('rowdim');
    symcol2row2tbl.rowdim1 = colrowtbl.get('coldim');
    symcol2row2tbl = IndexArray(symcol2row2tbl);

    prod = TensorProd();
    prod.tbl1 = symcol2row2tbl;
    prod.tbl2 = nodecolrowtbl;
    prod.tbl3 = nodecolrowtbl;
    prod.replacefds1 = {{'coldim1', 'coldim'}, ...
                        {'rowdim1', 'rowdim'}};
    prod.replacefds2 = {{'coldim', 'coldim2'}, ...
                        {'rowdim', 'rowdim2'}};
    prod.reducefds = {'coldim2', 'rowdim2'};

    symnodecol2row2tbl = crossIndexArray(nodetbl, symcol2row2tbl, {});
    nc2r2_num = symnodecol2row2tbl.num; % shortcut

    % (note the definition of symcol2row2tbl above)
    prod.pivottbl = symnodecol2row2tbl;
    [r, c, i] = ind2sub([d_num, d_num, n_num], (1 : nc2r2_num)');
    c2 = c;
    r2 = r;
    c1 = r;
    r1 = c;
    prod.dispind1 = sub2ind([d_num, d_num], r, c);
    prod.dispind2 = sub2ind([d_num, d_num, n_num], r2, c2, i);
    prod.dispind3 = sub2ind([d_num, d_num, n_num], r1, c1, i);
    prod.issetup = true;

    trans_T = SparseTensor('matlabsparse', true);
    trans_T = trans_T.setFromTensorProd(ones(symcol2row2tbl.num, 1), prod);

    % Construction of nodal average for cellnode tensor
    %
    % transnodeaverage_T : cellnodecolrowtbl -> nodecolrowtbl
    %
    % (later this operator is dispatched to cells)
    %

    % Compute number of cell per node
    map = TensorMap();
    map.fromTbl = cellnodetbl;
    map.toTbl = nodetbl;
    map.mergefds = {'nodes'};
    map = map.setup();
    
    nnodepercell = map.eval(ones(cellnodetbl.num, 1));
    
    map = TensorMap();
    map.fromTbl = nodetbl;
    map.toTbl = cellnodetbl;
    map.mergefds = {'nodes'};
    map = map.setup();
    
    coef = map.eval(1./nnodepercell);

    % we eliminitate the places (at the boundaries) where the local reconstruction
    % is ill-posed: nodes with one cell in 2d (corners of a Cartesian grid) and
    % nodes with less the two nodes in 3d (edges of a Cartesian grid);

    switch dim
      case 2
        maxnnodepercell = 1;
      case 3
        maxnnodepercell = 2;
    end

    clear fixnodetbl
    fixnodetbl.nodes = find(nnodepercell <= maxnnodepercell);
    fixnodetbl = IndexArray(fixnodetbl);

    coef(coef >= 1/maxnnodepercell) = 0;

    prod = TensorProd();
    prod.tbl1 = cellnodetbl;
    prod.tbl2 = cellnodecolrowtbl;
    prod.tbl3 = nodecolrowtbl;
    prod.reducefds = {'cells'};
    prod.mergefds = {'nodes'};

    prod.pivottbl = cellnodecolrowtbl;
    [r, c, i] = ind2sub([d_num, d_num, cn_num], (1 : cncr_num)');
    prod.dispind1 = i;
    prod.dispind2 = (1 : cncr_num)';
    prod.dispind3 = sub2ind([d_num, d_num, n_num], r, c, node_from_cellnode(i));
    prod.issetup = true;

    % prod = prod.setup();
    
    nodeaverage_T = SparseTensor('matlabsparse', true);
    nodeaverage_T = nodeaverage_T.setFromTensorProd(coef, prod);

    transnodeaverage_T = trans_T*nodeaverage_T;

    % We need to dispatch this tensor to cellnodecolrowtbl.
    % Now we have
    % transnodeaverage_T : cellnodecolrowtbl -> cellnodecolrowtbl

    map = TensorMap();
    map.fromTbl = nodecolrowtbl;
    map.toTbl = cellnodecolrowtbl;
    map.mergefds = {'nodes', 'coldim', 'rowdim'};
    
    map.pivottbl = cellnodecolrowtbl;
    [r, c, i] = ind2sub([d_num, d_num, cn_num], (1 : cncr_num)');
    map.dispind1 = sub2ind([d_num, d_num, n_num], r, c, node_from_cellnode(i));
    map.dispind2 = (1 : cncr_num)';
    map.issetup = true;
    % map = map.setup();
    
    celldispatch_T = SparseTensor('matlabsparse', true);
    celldispatch_T = celldispatch_T.setFromTensorMap(map);

    transnodeaverage_T = celldispatch_T*transnodeaverage_T;

    % We need to multiply by 2 at the place where we discarded the symmetry requirement

    coef = ones(nodetbl.num, 1);
    coef(fixnodetbl.get('nodes')) = 2;

    map = TensorMap();
    map.fromTbl = nodetbl;
    map.toTbl = cellnodecolrowtbl;
    map.mergefds = {'nodes'};

    map.pivottbl = cellnodecolrowtbl;
    [r, c, i] = ind2sub([d_num, d_num, cn_num], (1 : cncr_num)');
    map.dispind1 = node_from_cellnode(i);
    map.dispind2 = (1 : cncr_num)';
    map.issetup = true;
    % map = map.setup();
    
    coef = map.eval(coef);

    prod = TensorProd();
    prod.tbl1 = cellnodecolrowtbl;
    prod.tbl2 = cellnodecolrowtbl;
    prod.mergefds = {'cells', 'nodes', 'coldim', 'rowdim'};
    prod.tbl3 = cellnodecolrowtbl;

    prod.pivottbl = cellnodecolrowtbl;
    cncr_num = cellnodecolrowtbl.num; %shortcut
    prod.dispind1 = (1 : cncr_num)';
    prod.dispind2 = (1 : cncr_num)';
    prod.dispind3 = (1 : cncr_num)';
    prod.issetup = true;
    % prod = prod.setup();

    bcfix_T = SparseTensor('matlabsparse', true);
    bcfix_T = bcfix_T.setFromTensorProd(coef, prod);

    % Construction of the stiffness operator
    %
    % C_T : cellnodecolrowtbl -> cellnodecolrowtbl
    %
    
    map = TensorMap();
    map.fromTbl = cellcol2row2tbl;
    map.toTbl = cellnodecol2row2tbl;
    map.mergefds = {'cells', 'coldim1', 'coldim2', 'rowdim1', 'rowdim2'};

    map.pivottbl = cellnodecol2row2tbl;
    cnc2r2_num = cellnodecol2row2tbl.num; %shortcut
    [c2r2, i] = ind2sub([c2r2_num, cn_num], (1 : cnc2r2_num)');
    map.dispind1 = sub2ind([c2r2_num, c_num], c2r2, cell_from_cellnode(i));
    map.dispind2 = (1 : cnc2r2_num)';
    map.issetup = true;
    % map = map.setup();

    C = map.eval(C);

    prod = TensorProd();
    prod.tbl1 = cellnodecol2row2tbl;
    prod.tbl2 = cellnodecolrowtbl;
    prod.replacefds1 = {{'coldim1', 'coldim'}, {'rowdim1', 'rowdim'}};
    prod.replacefds2 = {{'coldim', 'coldim2'}, {'rowdim', 'rowdim2'}};
    prod.mergefds = {'cells', 'nodes'};
    prod.reducefds = {'coldim2', 'rowdim2'};
    prod.tbl3 = cellnodecolrowtbl;

    prod.pivottbl = cellnodecol2row2tbl;
    d = d_num; %shortcut
    [r2, c2, r1, c1, i] = ind2sub([d, d, d, d, cn_num], (1 : cnc2r2_num)');
    prod.dispind1 = (1 : cnc2r2_num)';
    prod.dispind2 = sub2ind([d, d, cn_num], r1, c1, i);
    prod.dispind3 = sub2ind([d, d, cn_num], r2, c2, i);
    prod.issetup = true;
    % prod = prod.setup();
    
    C_T = SparseTensor('matlabsparse', true);
    C_T = C_T.setFromTensorProd(C, prod);

    % Assembly

    Cgradnodeface_T = bcfix_T*C_T*gradnodeface_T;
    transaverCgradnodeface_T = transnodeaverage_T*Cgradnodeface_T;

    combCgradnodeface_T = 0.5*(Cgradnodeface_T + transaverCgradnodeface_T);

    Cgradcell_T = bcfix_T*C_T*gradcell_T;
    transaverCgradcell_T = transnodeaverage_T*Cgradcell_T;

    combCgradcell_T = 0.5*(Cgradcell_T + transaverCgradcell_T);

    C1 = combCgradnodeface_T.getMatrix();
    C2 = combCgradcell_T.getMatrix();
    
    A11 = divnodeface_T*combCgradnodeface_T;
    A12 = divnodeface_T*combCgradcell_T;
    A21 = divcell_T*combCgradnodeface_T;
    A22 = divcell_T*combCgradcell_T;

    A11 = A11.getMatrix();
    A12 = A12.getMatrix();
    A21 = A21.getMatrix();
    A22 = A22.getMatrix();

    [nodes, sz] = rlencode(nodefacecoltbl.get('nodes'), 1);
    opt.invertBlocks = 'mex';
    bi = blockInverter(opt);
    invA11 = bi(A11, sz);

    % Matrix for boundary conditions
    [D, bcvals] = setupMpsaNodeFaceBc(bc, G, nnodesperface, tbls);
    
    %
    % The divergence operator (integrated over the volume)
    % is given by 
    %
    %  div[c] = sum (m[f,s] u_[f,n,i] n[c,f,i])
    %
    % where u:solution, n:normal, m:area
    % indices : c:cell, f:face, n:node.
    
    % The facetNormals are already weighted with respect to area
    
    prod = TensorProd();
    prod.tbl1 = cellnodefacecoltbl;
    prod.tbl2 = nodefacecoltbl;
    prod.tbl3 = celltbl;
    prod.reducefds = {'faces', 'nodes', 'coldim'};
    % prod = prod.setup();
    
    prod.pivottbl = cellnodefacecoltbl;
    prod.dispind1 = (1 : cnfc_num)';
    [c, i] = ind2sub([d_num, cnf_num], (1 : cnfc_num)');
    prod.dispind2 = sub2ind([d_num, nf_num], c, nodeface_from_cellnodeface(i));
    prod.dispind3 = cell_from_cellnode(cellnode_from_cellnodeface(i));
    prod.issetup = true;
    
    div_T = SparseTensor;
    div_T = div_T.setFromTensorProd(facetNormals, prod);
    div = div_T.getMatrix();
    
    matrices = struct('A11'   , A11   , ...
                      'A12'   , A12   , ...
                      'A21'   , A21   , ...
                      'A22'   , A22   , ...
                      'D'     , D     , ...
                      'invA11', invA11, ...
                      'C1'    , C1    , ...
                      'C2'    , C2    , ...
                      'div'   , div);
    
    extra = struct('g', g);
end
