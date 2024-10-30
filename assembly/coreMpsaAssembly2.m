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
    eta       = opts.eta;
    
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
    prod.replacefds1 = {{'vec21', 'vec2'}};
    prod.replacefds2 = {{'vec', 'vec22'}};
    prod.reducefds   = {'vec22'};
    prod.mergefds    = {'cells'};
    prod = prod.setup();
    
    cellnodefacevec122tbl = prod.tbl3;
    Cg = prod.eval(C, g); % Cg is in cellnodefacevec122tbl

    %% Setup transpose CgT

    map = TensorMap();
    map.fromTbl           = cellnodefacevec122tbl;
    map.toTbl             = cellnodefacevec122tbl;
    map.replaceFromTblfds = {{'vec11', 'vec12', 'interchange'}};
    map.mergefds          = {'cells', 'nodes', 'faces', 'vec11', 'vec12', 'vec2'};
    map = map.setup();

    CgT = map.eval(Cg); 
    
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

    gen = CrossIndexArrayGenerator();
    gen.tbl1        = cellnodetbl;
    gen.tbl2        = cellnodetbl;
    gen.replacefds1 = {{'cells', 'cells1'}};
    gen.replacefds2 = {{'cells', 'cells2'}};
    gen.mergefds    = {'nodes'};

    cell12nodetbl = gen.eval();

    icell12nodetbl = cellnodetbl;
    icell12nodetbl = replacefield(icell12nodetbl, {{'cells', 'cells1'}});
    icell12nodetbl = icell12nodetbl.addInd('cells2', cellnodetbl.get('cells'));

    map = TensorMap();
    map.fromTbl  = icell12nodetbl;
    map.toTbl    = cell12nodetbl;
    map.mergefds = {'cells1', 'cells2', 'nodes'};
    map = map.setup();

    a = map.eval(ones(icell12nodetbl.num, 1));

    map = TensorMap();
    map.fromTbl           = cellnodetbl;
    map.toTbl             = cell12nodetbl;
    map.replaceFromTblfds = {{'cells', 'cells1'}};
    map.mergefds          = {'cells1', 'nodes'};
    map = map.setup();

    aT = map.eval(ones(cellnodetbl.num, 1));

    prod = TensorProd();
    prod.tbl1     = nodetbl;
    prod.tbl2     = cell12nodetbl;
    prod.tbl3     = cell12nodetbl;
    prod.mergefds = {'nodes'};
    prod = prod.setup();
    
    a  = 0.5*prod.eval(fixbc1, a);
    aT = 0.5*prod.eval(fixbc2./nnodepercell, aT);
    
    %% Setup of S
    
    prod = TensorProd();
    prod.tbl1 = cell12nodetbl;
    prod.tbl2 = cellnodefacevec122tbl;
    prod.replacefds2 = {{'cells', 'cells2'}};
    prod.mergefds = {'cells2', 'nodes'};
    prod = prod.setup();

    cell12nodefacevec122tbl = prod.tbl3;
    S = prod.eval(a, Cg) + prod.eval(aT, CgT);

    
    %% We multiply S by the normals to obtain nS

    % Setup normals at facets (get vector in cellnodefacevectbl)
    
    facetNormals = computeFacetNormals2(G, cellnodefacetbl);
    
    prod = TensorProd();
    prod.tbl1        = cellnodefacevectbl;
    prod.tbl2        = cell12nodefacevec122tbl;
    prod.replacefds1 = {{'faces', 'faces1'}, {'vec', 'redvec'}, {'cells', 'cells1'}};
    prod.replacefds2 = {{'faces', 'faces2'}, {'vec11', 'vec1'}, {'vec12', 'redvec'}};
    prod.reducefds   = {'redvec'};
    prod.mergefds    = {'cells1', 'nodes'};
    prod = prod.setup();

    cell12nodeface12vec12tbl = prod.tbl3;

    nS = prod.eval(facetNormals, S);

    %% Setup of the matrices A11, A12, A21, A22

    % Setup of A11
    
    prod = TensorProd();
    prod.tbl1        = cell12nodeface12vec12tbl;
    prod.tbl2        = nodefacevectbl;
    prod.tbl3        = nodefacevectbl;
    prod.replacefds1 = {{'faces1', 'faces'}, {'vec1', 'vec'}};
    prod.replacefds2 = {{'faces', 'faces2'}, {'vec', 'vec2'}};
    prod.reducefds   = {'faces2', 'vec2'};
    prod.reducefds1  = {'cells1', 'cells2'};
    prod.mergefds    = {'nodes'};
    prod = prod.setup();
    
    A11 = prod.setupMatrix(nS);

    % Setup of A12
    
    prod = TensorProd();
    prod.tbl1        = cell12nodeface12vec12tbl;
    prod.tbl2        = cellvectbl;
    prod.tbl3        = nodefacevectbl;
    prod.replacefds1 = {{'faces1', 'faces'}, {'vec1', 'vec'}};
    prod.replacefds2 = {{'vec', 'vec2'}, {'cells', 'cells2'}};
    prod.reducefds   = {'cells2', 'vec2'};
    prod.reducefds1  = {'cells1', 'faces2'};
    prod = prod.setup();

    % note the sign
    A12 = - prod.setupMatrix(nS);
    
    % Setup of A21
    
    prod = TensorProd();
    prod.tbl1        = cell12nodeface12vec12tbl;
    prod.tbl2        = nodefacevectbl;
    prod.tbl3        = cellvectbl;
    prod.replacefds1 = {{'vec1', 'vec'}, {'cells1', 'cells'}};
    prod.replacefds2 = {{'vec', 'vec2'}, {'faces', 'faces2'}};
    prod.reducefds   = {'nodes', 'faces2', 'vec2'};
    prod.reducefds1  = {'faces1', 'cells2'};
    prod = prod.setup();

    % note the sign
    A21 = - prod.setupMatrix(nS);
    
    % Setup of A22
    
    prod = TensorProd();
    prod.tbl1        = cell12nodeface12vec12tbl;
    prod.tbl2        = cellvectbl;
    prod.tbl3        = cellvectbl;
    prod.replacefds1 = {{'vec1', 'vec'}, {'cells1', 'cells'}};
    prod.replacefds2 = {{'vec', 'vec2'}, {'cells', 'cells2'}};
    prod.reducefds   = {'vec2', 'cells2'};
    prod.reducefds1  = {'faces1', 'faces2', 'nodes'};
    prod = prod.setup();

    A22 = prod.setupMatrix(nS);

    [nodes, sz] = rlencode(nodefacevectbl.get('nodes'), 1);
    opt.invertBlocks = 'm';
    bi = blockInverter(opt);
    invA11 = bi(A11, sz);

    % Matrix for boundary conditions
    [D, bcvals] = setupMpsaNodeFaceBc2(bc, G, nnodesperface, tbls);
    
    matrices = struct('A11'   , A11   , ...
                      'A12'   , A12   , ...
                      'A21'   , A21   , ...
                      'A22'   , A22   , ...
                      'D'     , D     , ...
                      'invA11', invA11);
    
    extra = struct('g', g);
    
end
