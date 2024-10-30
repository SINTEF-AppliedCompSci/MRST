function [matrices, bcvals, extra] = coreMpsaAssembly2(G, C, bc, nnodesperface, tbls, mappings, opts)
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
    
    vectbl       = tbls.vectbl;
    vec12tbl     = tbls.vec12tbl;
    vec1212tbl   = tbls.vec1212tbl;
    celltbl      = tbls.celltbl;
    nodetbl      = tbls.nodetbl;
    cellnodetbl  = tbls.cellnodetbl;
    nodefacetbl  = tbls.nodefacetbl;
    cellvectbl   = tbls.cellvectbl;
    cellvec12tbl = tbls.cellvec12tbl;
    nodevectbl   = tbls.nodevectbl;

    nodefacevectbl         = tbls.nodefacevectbl;
    cellnodefacetbl        = tbls.cellnodefacetbl;
    cellnodevectbl         = tbls.cellnodevectbl;
    cellnodevec12tbl       = tbls.cellnodevec12tbl;
    cellnodefacevectbl     = tbls.cellnodefacevectbl;
    cellnodefacevec12tbl   = tbls.cellnodefacevec12tbl;
    cellnodeface12tbl      = tbls.cellnodeface12tbl;
    cellnodeface12vec12tbl = tbls.cellnodeface12vec12tbl;
    nodeface12tbl          = tbls.nodeface12tbl;
    nodeface12vec12tbl     = tbls.nodeface12vec12tbl;    
    nodevec12tbl           = tbls.nodevec12tbl;
    nodefacevec12tbl       = tbls.nodefacevec12tbl;
    cellvec1212tbl         = tbls.cellvec1212tbl;
    cellnodevec1212tbl     = tbls.cellnodevec1212tbl;
    
    % Construction of tensor g (as defined in paper eq 4.1.2)
    % g belongs to cellnodefacevectbl
    g = computeConsistentGradient2(G, eta, tbls, mappings, 'bcetazero', bcetazero);    
    
    %% Computation of Cg
    
    prod = TensorProd();
    prod.tbl1        = cellvec1212tbl;
    prod.tbl2        = cellnodefacevectbl;
    prod.replacefds1 = {{'vec22', 'vec2'}};
    prod.replacefds2 = {{'vec', 'vec21'}};
    prod.reducefds   = {'vec21'};
    prod.mergefds    = {'cells'};
    prod = prod.setup();
    
    cellnodefacevec122tbl = prod.tbl3;
    Cg = prod.eval(C, g); % Cg is in cellnodefacevec122tbl

    % Compute number of cell per node, nnodepercell in nodetbl
    map = TensorMap();
    map.fromTbl  = cellnodetbl;
    map.toTbl    = nodetbl;
    map.mergefds = {'nodes'};
    map = map.setup();
    
    nnodepercell = map.eval(ones(cellnodetbl.num, 1));

    switch vectbl.num
      case 2
        maxnnodepercell = 1;
      case 3
        maxnnodepercell = 2;
    end

    fixnodes = find(nnodepercell <= maxnnodepercell);
    
    fixbc1 = ones(nodetbl.num, 1);
    fixbc1(fixnodes) = 2;
    
    fixbc2 = ones(nodetbl.num, 1);
    fixbc2(fixnodes) = 0;

    %% Setup Cg1
    
    prod = TensorProd();
    prod.tbl1 = nodetbl;
    prod.tbl2 = cellnodefacevec122tbl;
    prod.tbl3 = cellnodefacevec122tbl;
    prod.mergefds = {'nodes'};
    prod = prod.setup();
    
    Cg1 = prod.eval(fixbc1, Cg);

    %% Setup Cg2

    coef = (1./nnodepercell).*fixbc2;
    Cg2 = prod.eval(coef, Cg);

    % Take cell cell average around node
    
    prod = TensorProd();
    prod.tbl1        = cellnodetbl;
    prod.tbl2        = cellnodefacevec122tbl;
    prod.mergefds    = {'nodes'};
    prod.replacefds2 = {{'cells', 'cells2'}};
    prod.reducefds2  = {'cells2'};
    prod.mergefds    = {'nodes'};
    prod = prod.setup();

    extcellnodefacevec122tbl = prod.tbl3;
    
    delta = ones(cellnodetbl.num, 1);
    Cg2 = prod.eval(delta, Cg2);

    % Take transpose

    map = TensorMap();
    map.fromTbl           = extcellnodefacevec122tbl;
    map.toTbl             = extcellnodefacevec122tbl;
    map.replaceFromTblfds = {{'vec11', 'vec12', 'interchange'}};
    map.mergefds          = {'cells', 'nodes', 'faces', 'vec11', 'vec12', 'vec2'};
    map = map.setup();

    Cg2 = map.eval(Cg2);

    %% We reset Cg as the average of Cg1 and Cg2

    % map = TensorMap();
    % map.fromTbl  = cellnodefacevec122tbl;
    % map.toTbl    = extcellnodefacevec122tbl;
    % map.mergefds = {'cells', 'nodes', 'faces', 'vec11', 'vec12', 'vec2'};
    % map = map.setup();

    % Cg1 = map.eval(Cg1);
    
    % Cg = 0.5*(Cg1 + Cg2);
    
    % cellnodefacevec122tbl = extcellnodefacevec122tbl;
    
    %% We multiply Cg by the normals to obtain nCg

    % Setup normals at facets (get vector in cellnodefacevectbl)
    
    facetNormals = computeFacetNormals2(G, cellnodefacetbl);
    
    prod = TensorProd();
    prod.tbl1        = cellnodefacevectbl;
    prod.tbl2        = cellnodefacevec122tbl;
    prod.replacefds1 = {{'faces', 'faces1'}, {'vec', 'redvec'}};
    prod.replacefds2 = {{'faces', 'faces2'}, {'vec11', 'vec1'}, {'vec12', 'redvec'}};
    prod.reducefds   = {'redvec'};
    prod.mergefds    = {'cells', 'nodes'};
    prod = prod.setup();

    cellnodeface12vec12tbl = prod.tbl3;

    nCg = prod.eval(facetNormals, Cg);

    %% Setup of the matrices A11, A12, A21, A22

    % Setup of A11
    
    prod = TensorProd();
    prod.tbl1        = cellnodeface12vec12tbl;
    prod.tbl2        = nodefacevectbl;
    prod.tbl3        = nodefacevectbl;
    prod.replacefds1 = {{'faces1', 'faces'}, {'vec1', 'vec'}};
    prod.replacefds2 = {{'faces', 'faces2'}, {'vec', 'vec2'}};
    prod.reducefds   = {'faces2', 'vec2'};
    prod.reducefds1  = {'cells'};
    prod.mergefds    = {'nodes'};
    prod = prod.setup();
    
    A11 = prod.setupMatrix(nCg);

    keyboard
    
    % Setup of A12
    
    prod = TensorProd();
    prod.tbl1        = cellnodeface12vec12tbl;
    prod.tbl2        = cellvectbl;
    prod.tbl3        = nodefacevectbl;
    prod.replacefds1 = {{'faces1', 'faces'}, {'vec1', 'vec'}};
    prod.replacefds2 = {{'vec', 'vec2'}};
    prod.reducefds   = {'vec2', 'cells'};
    prod.reducefds1  = {'faces2'};
    prod = prod.setup();

    % note the sign
    A12 = -prod.setupMatrix(nCg);
    
    % Setup of A21
    
    prod = TensorProd();
    prod.tbl1        = cellnodeface12vec12tbl;
    prod.tbl2        = nodefacevectbl;
    prod.tbl3        = cellvectbl;
    prod.replacefds1 = {{'vec1', 'vec'}};
    prod.replacefds2 = {{'vec', 'vec2'}, {'faces', 'faces2'}};
    prod.reducefds   = {'nodes', 'faces2', 'vec2'};
    prod.reducefds1  = {'faces1'};
    prod = prod.setup();

    % note the sign
    A21 = -prod.setupMatrix(nCg);
    
    % Setup of A22
    
    prod = TensorProd();
    prod.tbl1        = cellnodeface12vec12tbl;
    prod.tbl2        = cellvectbl;
    prod.tbl3        = cellvectbl;
    prod.replacefds1 = {{'vec1', 'vec'}};
    prod.replacefds2 = {{'vec', 'vec2'}};
    prod.reducefds   = {'vec2'};
    prod.reducefds1  = {'faces1', 'faces2', 'nodes'};
    prod.mergefds    = {'cells'};
    prod = prod.setup();

    A22 = prod.setupMatrix(nCg);

    keyboard
    
    [nodes, sz] = rlencode(nodefacevectbl.get('nodes'), 1);
    opt.invertBlocks = 'm';
    bi = blockInverter(opt);
    invA11 = bi(A11, sz);

    % Matrix for boundary conditions
    [D, bcvals] = setupMpsaNodeFaceBc2(bc, G, nnodesperface, tbls);

    
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
