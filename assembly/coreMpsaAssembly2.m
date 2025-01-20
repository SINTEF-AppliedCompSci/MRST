function output =  coreMpsaAssembly2(G, C, bc, nnodesperface, tbls, mappings, opts, varargin)
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


    opt = struct('useVirtual', false);
    opt = merge_options(opt, varargin{:});

    useVirtual = opt.useVirtual;
    
    bcetazero = opts.bcetazero;
    eta       = opts.eta;
    
    vectbl      = tbls.vectbl;
    vec12tbl    = tbls.vec12tbl;
    vec1212tbl  = tbls.vec1212tbl;
    celltbl     = tbls.celltbl;
    nodetbl     = tbls.nodetbl;
    cellnodetbl = tbls.cellnodetbl;
    nodefacetbl = tbls.nodefacetbl;
    cellvectbl  = tbls.cellvectbl;
    
    dim = vectbl.num;
    
    nodefacevectbl           = tbls.nodefacevectbl;
    cellnodevec12tbl         = tbls.cellnodevec12tbl;
    cellnodefacetbl          = tbls.cellnodefacetbl;
    cellnodefacevectbl       = tbls.cellnodefacevectbl;
    cellvec1212tbl           = tbls.cellvec1212tbl;
    cellnodefacevec122tbl    = tbls.cellnodefacevec122tbl;
    cell12nodetbl            = tbls.cell12nodetbl;
    cell12nodefacevec122tbl  = tbls.cell12nodefacevec122tbl;
    cell12nodeface12vec12tbl = tbls.cell12nodeface12vec12tbl;

    if useVirtual
        
        d_num   = vectbl.num;
        c_num   = celltbl.num;
        cnf_num = cellnodefacetbl.num;

        cell_from_cellnode                      = mappings.cell_from_cellnode;
        face_from_cellface                      = mappings.face_from_cellface;
        cell_from_cellface                      = mappings.cell_from_cellface;
        node_from_cellnode                      = mappings.node_from_cellnode;
        cell_from_cellnodeface                  = mappings.cell_from_cellnodeface;
        cellnode_from_cellnodeface              = mappings.cellnode_from_cellnodeface;
        nodeface_from_cellnodeface              = mappings.nodeface_from_cellnodeface;
        
    end
    
    % Construction of tensor g (as defined in paper eq 4.1.2)
    % g belongs to cellnodefacevectbl
    g = computeConsistentGradient2(G, eta, tbls, mappings, 'bcetazero', bcetazero, 'useVirtual', useVirtual);    

    %% Computation of Cg
    
    prod = TensorProd();
    prod.tbl1        = cellvec1212tbl;
    prod.tbl2        = cellnodefacevectbl;
    prod.tbl3        = cellnodefacevec122tbl;
    prod.replacefds1 = {{'vec21', 'vec2'}};
    prod.replacefds2 = {{'vec', 'vec22'}};
    prod.reducefds   = {'vec22'};
    prod.mergefds    = {'cells'};

    if useVirtual
        cellnodefacevec1212tbl = crossIndexArray(cellnodefacetbl, vec1212tbl, {}, 'optpureproduct', true, 'virtual', true);
        prod.pivottbl = cellnodefacevec1212tbl;
        [vec22, vec21, vec12, vec11, i] = ind2sub([d_num, d_num, d_num, d_num, cnf_num], (1 : cellnodefacevec1212tbl.num)');

        prod.dispind1 = sub2ind([d_num, d_num, d_num, d_num, c_num], vec22, vec21, vec12, vec11 , cell_from_cellnodeface(i));
        prod.dispind2 = sub2ind([d_num, cnf_num], vec22, i);
        % order in cellnodefacevec122tbl is cellnodefac, vec11, vec12, vec2
        prod.dispind3 = sub2ind([d_num, d_num, d_num, cnf_num], vec21, vec12, vec11 , i);
        
        prod.issetup = true;
        
    else
        prod = prod.setup();
    end
    
    Cg = prod.eval(C, g); % Cg is in cellnodefacevec122tbl

    %% Setup transpose CgT

    map = TensorMap();
    map.fromTbl           = cellnodefacevec122tbl;
    map.toTbl             = cellnodefacevec122tbl;
    map.replaceFromTblfds = {{'vec11', 'vec12', 'interchange'}};
    map.mergefds          = {'cells', 'nodes', 'faces', 'vec11', 'vec12', 'vec2'};

    if useVirtual
        map.pivottbl = cellnodefacevec122tbl;
        % order in cellnodefacevec122tbl is cellnodefac, vec11, vec12, vec2
        [vec2, vec12, vec11, i] = ind2sub([d_num, d_num, d_num, cnf_num], (1 : cellnodefacevec122tbl.num)');
        map.dispind1 = (1 : cellnodefacevec122tbl.num)';
        map.dispind2 = sub2ind([d_num, d_num, d_num, cnf_num], vec2, vec11, vec12 , i);
        map.issetup = true;
    else
        map = map.setup();
    end

    CgT = map.eval(Cg); 
    
    % Compute number of cell per node, nnodepercell in nodetbl
    map = TensorMap();
    map.fromTbl  = cellnodetbl;
    map.toTbl    = nodetbl;
    map.mergefds = {'nodes'};
    map = map.setup();


    if useVirtual
        map.pivottbl = cellnodetbl;
        map.dispind1 = (1 : cellnodetbl.num)';
        map.dispind2 = node_from_cellnode;
        map.issetup = true;
    else
        map = map.setup();
    end
        
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
    prod.tbl1        = cell12nodetbl;
    prod.tbl2        = cellnodefacevec122tbl;
    prod.tbl3        = cell12nodefacevec122tbl;
    prod.replacefds2 = {{'cells', 'cells2'}};
    prod.mergefds    = {'cells2', 'nodes'};
    
    if useVirtual

        prod.pivottbl = cell12nodefacevec122tbl;

        cell12nodefacetbl = tbls.cell12nodefacetbl;
        cell12nodetbl     = tbls.cell12nodetbl;
        
        N = cell12nodefacetbl.num;
        [vec2, vec12, vec11, i] = ind2sub([d_num, d_num, d_num, N], (1 : cell12nodefacevec122tbl.num)');

        prod.dispind1 = mappings.cell12node_from_cell12nodeface(i);
        
        N = cellnodefacetbl.num;
        prod.dispind2 = sub2ind([d_num, d_num, d_num,  N], vec2, vec12, vec11 , mappings.cell2nodeface_from_cell12nodeface(i));

        prod.dispind3 = (1 : cell12nodefacevec122tbl.num)';
        
        prod.issetup = true;
        
    else
        prod = prod.setup();
    end

    S = prod.eval(a, Cg) + prod.eval(aT, CgT); % S is in cell12nodefacevec122tbl
    
    %% We multiply S by the normals to obtain nS

    % Setup normals at facets (get vector in cellnodefacevectbl)
    
    facetNormals = computeFacetNormals2(G, cellnodefacetbl);
    
    prod = TensorProd();
    prod.tbl1        = cellnodefacevectbl;
    prod.tbl2        = cell12nodefacevec122tbl;
    prod.tbl3        = cell12nodeface12vec12tbl;
    prod.replacefds1 = {{'faces', 'faces1'}, {'vec', 'redvec'}, {'cells', 'cells1'}};
    prod.replacefds2 = {{'faces', 'faces2'}, {'vec11', 'vec1'}, {'vec12', 'redvec'}};
    prod.reducefds   = {'redvec'};
    prod.mergefds    = {'cells1', 'nodes'};


    vec122tbl = tbls.vec122tbl;

    cell12nodeface12vec122tbl = crossIndexArray(tbls.cell12nodeface12tbl, vec122tbl, {}, ...
                                                'optpureproduct', true, ...
                                                'virtual', true);        
    prod.pivottbl = cell12nodeface12vec122tbl;
    
    map1 = mappings.cell1nodeface1_from_cell12nodeface12;
    map2 = mappings.cell12nodeface2_from_cell12nodeface12;

    [v2, v12, v11, i] = ind2sub([dim, dim, dim, tbls.cell12nodeface12tbl.num], (1 : cell12nodeface12vec122tbl.num)');

    prod.dispind1 = sub2ind([dim, tbls.cellnodefacetbl.num], v12, map1(i));
    prod.dispind2 = sub2ind([dim, dim, dim, tbls.cell12nodefacetbl.num], v2, v12, v11, map2(i));
    prod.dispind3 = sub2ind([dim, dim, tbls.cell12nodeface12tbl.num], v2, v11, i);

    prod.issetup = true;

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

    if useVirtual

        prod.pivottbl = cell12nodeface12vec12tbl;

        N = tbls.cell12nodeface12tbl.num;
        [vec2, vec1, i] = ind2sub([d_num, d_num, N], (1 : cell12nodeface12vec12tbl.num)');

        prod.dispind1 = (1 : cell12nodeface12vec12tbl.num)';
        
        N = nodefacetbl.num;
        j = mappings.cell12nodeface2_from_cell12nodeface12(i);
        j = mappings.cell2nodeface_from_cell12nodeface(j);
        j = nodeface_from_cellnodeface(j);
        prod.dispind2 = sub2ind([d_num, N], vec2, j);

        N = nodefacetbl.num;
        j = mappings.cell1nodeface1_from_cell12nodeface12(i);
        j = nodeface_from_cellnodeface(j);
        prod.dispind3 = sub2ind([d_num, N], vec1, j);
        
        prod.issetup = true;
        
    else

        prod = prod.setup();
        
    end
    
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

    if useVirtual

        prod.pivottbl = cell12nodeface12vec12tbl;

        N = tbls.cell12nodeface12tbl.num;
        [vec2, vec1, i] = ind2sub([d_num, d_num, N], (1 : cell12nodeface12vec12tbl.num)');

        prod.dispind1 = (1 : cell12nodeface12vec12tbl.num)';
        
        N = celltbl.num;
        j = mappings.cell12nodeface2_from_cell12nodeface12(i);
        j = mappings.cell12node_from_cell12nodeface(j);
        j = mappings.cell2node_from_cell12node(j);
        j = cell_from_cellnode(j);
        prod.dispind2 = sub2ind([d_num, N], vec2, j);

        N = nodefacetbl.num;
        j = mappings.cell1nodeface1_from_cell12nodeface12(i);
        j = mappings.nodeface_from_cellnodeface(j);
        prod.dispind3 = sub2ind([d_num, N], vec1, j);
        
        prod.issetup = true;
        
    else

        prod = prod.setup();
        
    end
    
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

    if useVirtual

        prod.pivottbl = cell12nodeface12vec12tbl;

        N = tbls.cell12nodeface12tbl.num;
        [vec2, vec1, i] = ind2sub([d_num, d_num, N], (1 : cell12nodeface12vec12tbl.num)');

        prod.dispind1 = (1 : cell12nodeface12vec12tbl.num)';

        N = nodefacetbl.num;
        j = mappings.cell12nodeface2_from_cell12nodeface12(i);
        j = mappings.cell2nodeface_from_cell12nodeface(j);
        j = mappings.nodeface_from_cellnodeface(j);
        prod.dispind2 = sub2ind([d_num, N], vec2, j);

        N = celltbl.num;
        j = mappings.cell1nodeface1_from_cell12nodeface12(i);
        j = mappings.cell_from_cellnodeface(j);
        prod.dispind3 = sub2ind([d_num, N], vec1, j);
        
        prod.issetup = true;
        
    else

        prod = prod.setup();
        
    end
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
    if useVirtual

        prod.pivottbl = cell12nodeface12vec12tbl;

        N = tbls.cell12nodeface12tbl.num;
        [vec2, vec1, i] = ind2sub([d_num, d_num, N], (1 : cell12nodeface12vec12tbl.num)');

        prod.dispind1 = (1 : cell12nodeface12vec12tbl.num)';

        N = celltbl.num;
        j = mappings.cell12nodeface2_from_cell12nodeface12(i);
        j = mappings.cell2nodeface_from_cell12nodeface(j);
        j = mappings.cell_from_cellnodeface(j);
        prod.dispind2 = sub2ind([d_num, N], vec2, j);

        N = celltbl.num;
        j = mappings.cell1nodeface1_from_cell12nodeface12(i);
        j = mappings.cell_from_cellnodeface(j);
        prod.dispind3 = sub2ind([d_num, N], vec1, j);
        
        prod.issetup = true;
        
    else

        prod = prod.setup();
        
    end
    
    A22 = prod.setupMatrix(nS);

    [nodes, sz] = rlencode(nodefacevectbl.get('nodes'), 1);
    opt.invertBlocks = 'm';
    bi = blockInverter(opt);
    invA11 = bi(A11, sz);

    % Matrix for boundary conditions
    
    [D, bcvals] = setupMpsaNodeFaceBc2(bc, G, nnodesperface, tbls);

    % Matrix for the stress computation

    prod = TensorProd();
    prod.tbl1        = cell12nodefacevec122tbl;
    prod.tbl2        = nodefacevectbl;
    prod.tbl3        = cellnodevec12tbl;
    prod.replacefds1 = {{'cells1', 'cells'}, {'vec2', 'redvec'}, {'vec11', 'vec1'}, {'vec12', 'vec2'}};
    prod.replacefds2 = {{'vec', 'redvec'}};
    prod.reducefds   = {'redvec', 'faces'};
    prod.reducefds1  = {'cells2'};
    prod.mergefds    = {'nodes'};

    if useVirtual

        prod.pivottbl = cell12nodefacevec122tbl;

        N = tbls.cell12nodefacetbl.num;
        [vec2, vec12, vec11, i] = ind2sub([d_num, d_num, d_num, N], (1 : cell12nodefacevec122tbl.num)');

        prod.dispind1 = (1 : cell12nodefacevec122tbl.num)';

        N = nodefacetbl.num;
        j = mappings.cell2nodeface_from_cell12nodeface(i);
        j = mappings.nodeface_from_cellnodeface(j);
        prod.dispind2 = sub2ind([d_num, N], vec2, j);

        N = cellnodetbl.num;
        j = mappings.cell1node_from_cell12nodeface(i);
        prod.dispind3 = sub2ind([d_num, d_num, N], vec12, vec11, j);
        
        prod.issetup = true;
        
    else

        prod = prod.setup();
        
    end

    C1 = prod.setupMatrix(S);

    prod = TensorProd();
    prod.tbl1        = cell12nodefacevec122tbl;
    prod.tbl2        = cellvectbl;
    prod.tbl3        = cellnodevec12tbl;
    prod.replacefds1 = {{'cells1', 'cells'}, {'vec2', 'redvec'}, {'vec11', 'vec1'}, {'vec12', 'vec2'}};
    prod.replacefds2 = {{'vec', 'vec2'}, {'cells', 'cells2'}};
    prod.reducefds   = {'vec2', 'cells2'};
    prod.reducefds1  = {'faces'};

    if useVirtual

        prod.pivottbl = cell12nodefacevec122tbl;

        N = tbls.cell12nodefacetbl.num;
        [vec2, vec12, vec11, i] = ind2sub([d_num, d_num, d_num, N], (1 : cell12nodefacevec122tbl.num)');

        prod.dispind1 = (1 : cell12nodefacevec122tbl.num)';

        N = celltbl.num;
        j = mappings.cell2nodeface_from_cell12nodeface(i);
        j = mappings.cell_from_cellnodeface(j);
        prod.dispind2 = sub2ind([d_num, N], vec2, j);

        N = cellnodetbl.num;
        j = mappings.cell1node_from_cell12nodeface(i);
        prod.dispind3 = sub2ind([d_num, d_num, N], vec12, vec11, j);
        
        prod.issetup = true;
        
    else

        prod = prod.setup();
        
    end

    % note the sign
    C2 = -prod.setupMatrix(S);
    
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
    prod.tbl1 = cellnodefacevectbl;
    prod.tbl2 = nodefacevectbl;
    prod.tbl3 = celltbl;
    prod.reducefds = {'faces', 'nodes', 'vec'};

    if useVirtual
        
        prod.pivottbl = cellnodefacevectbl;

        N = tbls.cellnodefacetbl.num;
        [vec, i] = ind2sub([d_num, N], (1 : cellnodefacevectbl.num)');

        prod.dispind1 = (1 : cellnodefacevectbl.num)';

        N = nodefacevectbl.num;
        j = mappings.nodeface_from_cellnodeface(i);
        prod.dispind2 = sub2ind([d_num, N], vec, j);

        prod.dispind3 = mappings.cell_from_cellnodeface(i);
        
        prod.issetup = true;
        
    else
        prod = prod.setup();
    end
    
    div = prod.setupMatrix(facetNormals);

    matrices = struct('A11'   , A11   , ...
                      'A12'   , A12   , ...
                      'A21'   , A21   , ...
                      'A22'   , A22   , ...
                      'D'     , D     , ...
                      'invA11', invA11, ...
                      'C1'    , C1    , ...
                      'C2'    , C2    , ...
                      'div'   , div );
    
    extra = struct('g', g);

    tbls.cell12nodeface12vec12tbl = cell12nodeface12vec12tbl;
    
    indexarrays = struct('nS', nS);

    output = struct('tbls'       , tbls       , ...
                    'indexarrays', indexarrays, ...
                    'matrices'   , matrices   , ...
                    'bcvals'     , bcvals     , ...
                    'extra'      , extra);

end
