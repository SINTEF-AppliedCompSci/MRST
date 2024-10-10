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
    
    g = computeConsistentGradient2(G, eta, tbls, mappings, 'bcetazero', bcetazero);

    prod = TensorProd();
    prod.tbl1 = cellnodefacevectbl;
    prod.tbl2 = nodefacevectbl;
    prod.tbl3 = cellnodevec12tbl;
    prod.replacefds1 = {{'vec', 'vec1'}};
    prod.replacefds2 = {{'vec', 'vec2'}};
    prod.reducefds = {'faces'};
    prod.mergefds = {'nodes'};
    prod = prod.setup();
    
    gmat = prod.setupMatrix(g);

    prod = TensorProd();
    prod.tbl1 = cellvec1212tbl;
    prod.tbl2 = cellnodevec12tbl;
    prod.tbl3 = cellnodevec12tbl;
    prod.replacefds1 = {{'vec11', 'vec1'}, {'vec12', 'vec2'}, {'vec21', 'redvec1'}, {'vec22', 'redvec2'}};
    prod.replacefds2 = {{'vec1', 'redvec1'}, {'vec2', 'redvec2'}};
    prod.reducefds = {'redvec1', 'redvec2'};
    prod.mergefds = {'cells'};
    prod = prod.setup();

    Cmat = prod.setupMatrix(C);

    facetNormals = computeFacetNormals2(G, cellnodefacetbl);
    
    prod = TensorProd();
    prod.tbl1 = cellnodefacevectbl;
    prod.tbl2 = cellnodevec12tbl;
    prod.tbl3 = nodefacevectbl;
    prod.replacefds1 = {{'vec', 'vec1'}};
    prod.replacefds2 = {{'vec2', 'vec'}};
    prod.reducefds = {'cells', 'vec1'};
    prod.mergefds = {'nodes'};
    prod = prod.setup();

    divmat = prod.setupMatrix(facetNormals);

    Cgmat = Cmat*gmat;
    
    gen = CrossIndexArrayGenerator();
    gen.tbl1 = cellnodefacevec12tbl;
    gen.tbl2 = vectbl;
    gen.replacefds1 = {{'vec1', 'vec11'}, {'vec2', 'vec12'}};
    gen.replacefds2 = {{'vec', 'vec2'}};

    cellnodefacevec122tbl = gen.eval();

    tconv = TensorConvert;
    tconv.fromTbl           = nodefacevectbl;
    tconv.toTbl             = cellnodevec12tbl;
    tconv.pivottbl          = cellnodefacevec122tbl;
    tconv.replaceFromTblfds = {{'vec', 'vec2'}};
    tconv.replaceToTblfds   = {{'vec1', 'vec11'}, {'vec2', 'vec12'}};
    tconv = tconv.setup();

    Cg2 = tconv.convert(Cgmat);

    % gen = CrossIndexArrayGenerator();
    % gen.tbl1 = nodefacevec12tbl;
    % gen.tbl2 = vectbl;
    % gen.replacefds1 = {{'vec1', 'vec11'}, {'vec2', 'vec12'}};
    % gen.replacefds2 = {{'vec', 'vec2'}};

    % nodefacevec122tbl = gen.eval();
    
    prod = TensorProd();
    prod.tbl1 = cellvec1212tbl;
    prod.tbl2 = cellnodefacevectbl;
    prod.tbl3 = cellnodefacevec122tbl;
    prod.replacefds1 = {{'vec22', 'vec2'}, {'vec21', 'redvec'}};
    % prod.replacefds1 = {{'vec22', 'redvec'}, {'vec21', 'vec2'}};
    prod.replacefds2 = {{'vec', 'redvec'}};
    prod.reducefds   = {'redvec'};
    prod.mergefds    = {'cells'};
    prod = prod.setup();

    Cg = prod.eval(C, g);

    prod = TensorProd();
    prod.tbl1 = cellnodefacevec122tbl;
    prod.tbl2 = nodefacevectbl;
    prod.tbl3 = cellnodevec12tbl;
    prod.replacefds1 = {{'vec2', 'redvec'}, {'vec11', 'vec1'}, {'vec12', 'vec2'}};
    prod.replacefds2 = {{'vec', 'redvec'}};
    prod.reducefds = {'faces', 'redvec'};
    prod.mergefds = {'nodes'};
    prod = prod.setup();

    Cgmat2 = prod.setupMatrix(Cg);

    prod = TensorProd();
    prod.tbl1        = cellnodefacevec122tbl;
    prod.tbl2        = cellnodefacevectbl;
    prod.tbl3        = nodeface12vec12tbl;
    prod.replacefds1 = {{'vec11', 'redvec'}, {'vec12', 'vec1'}, {'faces', 'faces2'}};
    prod.replacefds2 = {{'vec', 'redvec'}, {'faces', 'faces1'}};
    prod.reducefds   = {'redvec', 'cells'};
    prod.mergefds    = {'nodes'};
    prod = prod.setup();
    
    nCg1 = prod.eval(Cg, facetNormals);

    prod = TensorProd();
    prod.tbl1        = nodeface12vec12tbl;
    prod.tbl2        = nodefacevectbl;
    prod.tbl3        = nodefacevectbl;
    prod.replacefds1 = {{'faces1', 'faces'}, {'vec1', 'vec'}};
    prod.replacefds2 = {{'faces', 'faces2'}, {'vec', 'vec2'}};
    prod.mergefds    = {'nodes'};
    prod.reducefds   = {'faces2', 'vec2'};
    prod = prod.setup();

    A11 = prod.setupMatrix(nCg1);
    
    
    % Compute number of cell per node
    map = TensorMap();
    map.fromTbl = cellnodetbl;
    map.toTbl = nodetbl;
    map.mergefds = {'nodes'};
    map = map.setup();
    
    ncellpernode = map.eval(ones(cellnodetbl.num, 1));

    switch vectbl.num
      case 2
        maxncellpernode = 1;
      case 3
        maxncellpernode = 2;
    end

    fixnodes = find(ncellpernode <= maxncellpernode);
    
    fixbc1 = 0.5*ones(nodetbl.num, 1);
    fixbc1(fixnodes) = 1;
    
    fixbc2 = 0.5*ones(nodetbl.num, 1);
    fixbc2(fixnodes) = 0;

    prod = TensorProd();
    prod.tbl1 = nodetbl;
    prod.tbl2 = cellnodefacevec122tbl;
    prod.tbl3 = cellnodefacevec122tbl;
    prod.mergefds = {'nodes'};
    prod = prod.setup();
    
    Cg1 = prod.eval(fixbc1, Cg);

    coef = (1./ncellpernode).*fixbc2;

    map = TensorMap();
    map.fromTbl = nodetbl;
    map.toTbl = cellnodetbl;
    map.mergefds = {'nodes'};
    map = map.setup();

    coef = map.eval(coef);

    prod = TensorProd();
    prod.tbl1 = cellnodetbl;
    prod.tbl2 = cellnodefacevec122tbl;
    prod.reducefds = {'cells'};
    prod.mergefds = {'nodes'};
    prod = prod.setup();

    nodefacevec122tbl = prod.tbl3;
    
    Cg2 = prod.eval(coef, Cg);

    map = TensorMap();
    map.fromTbl = nodefacevec122tbl;
    map.toTbl = cellnodefacevec122tbl;
    map.replaceFromTblfds = {{'vec11', 'vec12', 'interchange'}};
    map.mergefds = {'nodes', 'faces', 'vec11', 'vec12', 'vec2'};
    map = map.setup();
    
    Cg2 = map.eval(Cg2);

    Cg = Cg1 + Cg2;
    
    % setup normals at facets (get vector in cellnodefacevectbl)
    facetNormals = computeFacetNormals2(G, cellnodefacetbl);

    prod = TensorProd();
    prod.tbl1        = cellnodefacevec122tbl;
    prod.tbl2        = cellnodefacevectbl;
    prod.replacefds1 = {{'vec11', 'redvec'}, {'vec12', 'vec1'}, {'faces', 'faces2'}};
    prod.replacefds2 = {{'vec', 'redvec'}, {'faces', 'faces1'}};
    prod.reducefds   = {'redvec'};
    prod.mergefds    = {'cells', 'nodes'};
    prod = prod.setup();

    cellnodeface12vec12tbl = prod.tbl3;
    
    nCg = prod.eval(Cg, facetNormals);
    
    % setup of A11

    % nodeface12vec12tbl = projIndexArray(cellnodeface12vec12tbl, {'nodes', 'faces1', 'faces2', 'vec1', 'vec2'});
    
    map = TensorMap();
    map.fromTbl  = cellnodeface12vec12tbl;
    map.toTbl    = nodeface12vec12tbl;
    map.mergefds = {'nodes', 'faces1', 'faces2', 'vec1', 'vec2'};
    map = map.setup();

    A11 = map.eval(nCg);

    % setup of A12

    map = TensorMap();
    map.fromTbl           = cellnodeface12vec12tbl;
    map.toTbl             = cellnodefacevec12tbl;
    map.replaceFromTblfds = {{'faces1', 'faces'}};
    map.mergefds          = {'cells', 'nodes', 'faces', 'vec1', 'vec2'};
    map = map.setup();

    A12 = map.eval(-nCg);

    % setup of A21

    map = TensorMap();
    map.fromTbl           = cellnodeface12vec12tbl;
    map.toTbl             = cellnodefacevec12tbl;
    map.replaceFromTblfds = {{'faces2', 'faces'}};
    map.mergefds          = {'cells', 'nodes', 'faces', 'vec1', 'vec2'};
    map = map.setup();

    A21 = map.eval(-nCg);

    % setup of A22

    map = TensorMap();
    map.fromTbl  = cellnodeface12vec12tbl;
    map.toTbl    = cellvec12tbl;
    map.mergefds = {'cells', 'vec1', 'vec2'};
    map = map.setup();

    A22 = map.eval(nCg);

    prod = TensorProd();
    prod.tbl1        = nodeface12vec12tbl;
    prod.tbl2        = nodefacevectbl;
    prod.tbl3        = nodefacevectbl;
    prod.replacefds1 = {{'faces1', 'faces'}, {'vec1', 'vec'}};
    prod.replacefds2 = {{'faces', 'faces2'}, {'vec', 'vec2'}};
    prod.mergefds    = {'nodes'};
    prod.reducefds   = {'faces2', 'vec2'};
    prod = prod.setup();

    A11 = prod.setupMatrix(A11);

    [nodes, sz] = rlencode(nodefacevectbl.get('nodes'), 1);
    opt.invertBlocks = 'm';
    bi = blockInverter(opt);
    invA11 = bi(A11, sz);

    prod = TensorProd();
    prod.tbl1        = cellnodefacevec12tbl;
    prod.tbl2        = cellvectbl;
    prod.tbl3        = nodefacevectbl;
    prod.replacefds1 = {{'vec1', 'vec'}};
    prod.replacefds2 = {{'vec', 'vec2'}};
    prod.reducefds   = {'cells', 'vec2'};
    prod = prod.setup();

    A12 = prod.setupMatrix(A12);

    prod = TensorProd();
    prod.tbl1        = cellnodefacevec12tbl;
    prod.tbl2        = nodefacevectbl;
    prod.tbl3        = cellvectbl;
    prod.replacefds1 = {{'vec1', 'vec'}};
    prod.replacefds2 = {{'vec', 'vec2'}};
    prod.reducefds   = {'nodes', 'faces', 'vec2'};
    prod = prod.setup();

    A21 = prod.setupMatrix(A21);
    
    prod = TensorProd();
    prod.tbl1        = cellvec12tbl;
    prod.tbl2        = cellvectbl;
    prod.tbl3        = cellvectbl;
    prod.replacefds1 = {{'vec1', 'vec'}};
    prod.replacefds2 = {{'vec', 'vec2'}};
    prod.mergefds    = {'cells'};
    prod.reducefds   = {'vec2'};
    prod = prod.setup();
    
    A22 = prod.setupMatrix(A22);

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
