function vagstruct = computeVagTrans(G, rock)
%
%
% SYNOPSIS:
%   function vagstruct = computeVagTrans(G, rock)
%
% DESCRIPTION: Computes the VAG transmissibilities. The implementation makes
% extensive use of the indexing table mechanism, introduced utils/tbl_utils in
% module mrst-core. It is not clear yet whether indexing tables scale well to
% large problem. Therefore, the function computeVagTrans may be slow.
%
% PARAMETERS:
%   G    - Grid
%   rock - Rock structure 
%
% RETURNS:
%   vagstruct - structure containing
%
%                  vagstruct.A - vector containing the transmissibility.
%                                The index is given by the table cellnode2tbl
%
%                  vagstruct.cellnode2tbl - Indexing table of the cell-node-node
%                                           triplets (for each cell, all pairs of nodes that belong to
%                                           that cell)
%
%                  vagstruct.cellnodetbl  - Indexing table of the cell-node
%                                           connectivities (for each cell,
%                                           the nodes that belong to the cell)
%
% SEE ALSO:
%   `computeMimeticIP`.

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    %% Setup cell, face, node and vertices indexing tables
        
    % NOTE : By vertices, we denote a point that can either be a cell centroid, a
    % face centroid or a node position. The vertices have their own numbering
    % (from 1 to np, see below). 
    assert(G.griddim == 3, 'only 3D code for now!');
    
    nc = G.cells.num;
    nf = G.faces.num;
    nn = G.nodes.num;

    celltbl.cells    = (1 : nc)';
    celltbl.vertices = (1 : nc)';
    celltbl.num      = nc;

    facetbl.faces    = (1 : nf)';
    facetbl.vertices = nc + (1 : nf)';
    facetbl.num      = nf;

    nodetbl.nodes    = (1 : nn)';
    nodetbl.vertices = nc + nf + (1 : nn)';
    nodetbl.num      = nn;

    np = nc + nf + nn;
    verttbl.vertices = (1 : np)';
    verttbl.num      = np;
    vertfds = {'vertices'};
    
    %% Setup vertcoltbl and vertcent
    
    % vertcent contains the coordinates of the vertices and is indexed using
    % the indexing table vertcoltbl
     
    coltbl.coldim   = (1 : 3)';
    coltbl.num      = 3;


    % We collect the cell, face and node centroids given by the grid structure in
    % vectors that follows the indexing of the indexing tables cellcoltbl,
    % facecoltbl and nodecoltbl.
    cellcent = G.cells.centroids;
    cellcent = reshape(cellcent', [], 1);
    facecent = G.faces.centroids;
    facecent = reshape(facecent', [], 1);
    nodecent = G.nodes.coords;
    nodecent = reshape(nodecent', [], 1);

    vertcoltbl = generateSubspace(verttbl, coltbl, []);
    cellcoltbl = generateSubspace(celltbl, coltbl, []);
    facecoltbl = generateSubspace(facetbl, coltbl, []);
    nodecoltbl = generateSubspace(nodetbl, coltbl, []);

    [~, indstruct] = generateSubspace(cellcoltbl, vertcoltbl, {'vertices', 'coldim'});
    vertcent = tblmap1to2(cellcent, indstruct);
    [~, indstruct] = generateSubspace(facecoltbl, vertcoltbl, {'vertices', 'coldim'});
    vertcent = vertcent + tblmap1to2(facecent, indstruct);
    [~, indstruct] = generateSubspace(nodecoltbl, vertcoltbl, {'vertices', 'coldim'});
    vertcent = vertcent + tblmap1to2(nodecent, indstruct);

    %% Setup tetratbl
    
    % There corresponds a tetrahedra for each triplet (cell, face, edge), where the
    % face belongs to the face and the edge to the face. The indexing table
    % tetratbl describe all these tetrahedra. The vertices of the tetrahedra are
    % given by the cell centroid, the face centroid, and the two nodes that
    % constitute the edge.
    

    nc = G.cells.num;
    clear cellfacetbl;
    cellfacetbl.cells = rldecode((1 : nc)', diff(G.cells.facePos));
    cellfacetbl.faces = G.cells.faces(:, 1);
    cellfacetbl.num   = numel(cellfacetbl.cells);

    nf = G.faces.num;
    clear faceedgetbl;
    if ~isfield(G.faces, 'edgePos')
        G = computeVEMGeometry(G);
    end
    faceedgetbl.faces = rldecode((1 : nf)', diff(G.faces.edgePos));
    faceedgetbl.edges = G.faces.edges;
    faceedgetbl.num   = numel(faceedgetbl.faces);

    [~, cellfaceedgetbl] = setupTableMapping(cellfacetbl, faceedgetbl, ...
                                                          {'faces'});
    

    %% Setup tetraverttbl
    
    % We set the table tetraverttbl which, for each tetrahedra, gives all the
    % vertices that constitute the tetrahedra.
    %
    % In addition, the vertices are given a index indverts (1, 2 or 3) depending on
    % which entity they correspond to (cell, face or node).

    tetratbl = cellfaceedgetbl;
    tetrafds = {'cells', 'faces', 'edges'};

    [~, tetracoltbl] = setupTableMapping(coltbl, tetratbl, []);
    tetracolfds = {'cells', 'faces', 'edges', 'coldim'};

    ne = G.edges.num;
    clear edgetbl;
    edgetbl.edges = (1 : ne)';
    edgetbl.nodes = G.edges.nodes(G.edges.nodePos(1 : (end - 1)));
    edgetbl.num   = ne;

    [~, tetratbl2] = setupTableMapping(tetratbl, edgetbl, {'edges'});
    tetratbl2 = replacefield(tetratbl2, {'nodes', 'nodes1'});

    clear edgetbl;
    edgetbl.edges = (1 : ne)';
    edgetbl.nodes = G.edges.nodes(G.edges.nodePos(1 : (end - 1)) + 1);
    edgetbl.num   = ne;

    [~, tetratbl2] = setupTableMapping(tetratbl2, edgetbl, {'edges'});
    tetratbl2 = replacefield(tetratbl2, {'nodes', 'nodes2'});

    exttetrafds = {'cells', 'faces', 'edges', 'nodes1', 'nodes2'};

    verts = cell(4, 1);
    indverts = cell(4, 1);
    nt = tetratbl.num;
    verts{1}    = tblmap(celltbl.vertices, celltbl, tetratbl2, {'cells'});
    indverts{1} = ones(nt, 1);
    verts{2}    = tblmap(facetbl.vertices, facetbl, tetratbl2, {'faces'});
    indverts{2} = 2*ones(nt, 1);
    verts{3}    = tblmap(nodetbl.vertices, nodetbl, tetratbl2, {{'nodes', 'nodes1'}});
    indverts{3} = 3*ones(nt, 1);
    verts{4}    = tblmap(nodetbl.vertices, nodetbl, tetratbl2, {{'nodes', 'nodes2'}});
    indverts{4} = 3*ones(nt, 1);

    clear tetraverttbl;
    tetraverttbl.vertices = vertcat(verts{:});
    tetraverttbl.indverts = vertcat(indverts{:});
    tetraverttbl.cells    = repmat(tetratbl2.cells, 4, 1);
    tetraverttbl.faces    = repmat(tetratbl2.faces, 4, 1);
    tetraverttbl.edges    = repmat(tetratbl2.edges, 4, 1);
    tetraverttbl.num = numel(tetraverttbl.vertices);

    tetravertindtbl = tetraverttbl;
    tetraverttbl = rmfield(tetraverttbl, 'indverts');

    tetravertfds = {'cells', 'faces', 'edges', 'vertices'};

    [~, tetravertcoltbl] = setupTableMapping(tetraverttbl, coltbl, []);

    tetravertcolfds = {'cells', 'faces', 'edges', 'vertices', 'coldim'};

    %% Compute the normals for tetraverttbl

    % For each vertex in a tetrahedra, we are going to compute the normal of the
    % face that is opposite to the tetrahedra. The normals will be stored in a
    % vector that belongs that follows the indexing table tetravertcoltbl
    
    
    % We first set up the table tetravertopptbl, which contains for each
    % vertex of each tetrahedra, the three opposite vertices.

    
    % start the construction of tetraveropptbl
    [~, tetravertopptbl] = setupTableMapping(tetraverttbl, tetraverttbl, {'cells', ...
                        'faces', 'edges'}, 'crossextend', {{'vertices', {'vertices', ...
                        'oppvertices'}}});

    tetravertoppfds = {'vertices', 'oppvertices', 'cells', 'faces', 'edges'};
    A = convertTableToArray(tetravertopptbl, tetravertoppfds);
    A = A(A(:, 1) ~= A(:, 2), :);

    tetravertopptbl = convertArrayToTable(A, tetravertoppfds);
    A = convertTableToArray(tetravertopptbl, {'cells', 'faces', 'edges', 'vertices', ...
                        'oppvertices'});
    A = sortrows(A);

    loctbl.oppvertinds = [1; 2; 3];
    loctbl.num = 3;
    a = repmat(loctbl.oppvertinds, tetraverttbl.num, 1);

    A = [A, a];

    tetravertopptbl = convertArrayToTable(A, {'cells', 'faces', 'edges', 'vertices', ...
                        'oppvertices', 'oppvertinds'});
    % end of the construction of tetraveropptbl
    
    % We construct the table tetravertoppcoltbl, which contains a vector for
    % each of the opposite vertices
    tetravertoppcoltbl = generateSubspace(tetravertopptbl, coltbl, []);

    % tetraoppvertcent gives the coordinate of each of the opposite vertices
    tetraoppvertcent = tblmap(vertcent, vertcoltbl, tetravertoppcoltbl, {{'vertices', 'oppvertices'}, 'coldim'});
    
    % We set up the vectors u1 and u2. Given an opposite face, let us denote the
    % vertices by P1, P2 and P3. Then, we have u1 = P1 - P2 and u2 = P1 -
    % P3. The normal of the opposite face will then be obtained by computing
    % the cross-product of u1 and u2.

    prod = TensorProd();
    prod.tbl1 = loctbl;
    prod.tbl2 = tetravertoppcoltbl;
    prod.reducefds = {'oppvertinds'};
    prod.prodtbl   = tetravertcoltbl;
    prod = prod.setup();
    
    u1 = prod.evalProd([1; -1; 0], tetraoppvertcent);
    u2 = prod.evalProd([1; 0; -1], tetraoppvertcent);
    
    % We start computing the cross product of u1 and u2
    colrowcrosstbl.coldim   = [2; 3; 1; 3; 1; 2];
    colrowcrosstbl.rowdim   = [3; 2; 3; 1; 2; 1];
    colrowcrosstbl.crossdim = [1; 1; 2; 2; 3; 3];
    colrowcrosstbl.num = numel(colrowcrosstbl.coldim);

    crossvec = [1; -1; -1; 1; 1; -1];
    [tetravertcolrowcrosstbl, indstruct] = generateSubspace(colrowcrosstbl, tetraverttbl, []);
    crossvec = tbldispatch1(crossvec, indstruct);
    
    prod = TensorProd();
    prod.tbl1 = tetravertcoltbl;
    prod.tbl2 = tetravertcolrowcrosstbl;
    prod.reducefds = {'coldim'};
    prod.mergefds  = tetravertfds;
    prod = prod.setup();
    
    oppfacenormals = prod.evalProd(u1, crossvec);
    
    
    tetravertrowcrosstbl = prod.prodtbl;
    tetravertcolcrosstbl = replacefield(tetravertrowcrosstbl, {'rowdim', ...
                        'coldim'});
    
    prod = TensorProd();
    prod.tbl1 = tetravertcolcrosstbl;
    prod.tbl2 = tetravertcoltbl;
    prod.reducefds = {'coldim'};
    prod.mergefds = tetravertfds;
    prod = prod.setup();
    
    oppfacenormals = prod.evalProd(oppfacenormals, u2);
    
    tetravertcrosstbl = prod.prodtbl;
    tetravertcrosstbl = replacefield(tetravertcrosstbl, {'crossdim', 'coldim'});

    tetravertcolfds = gettblfds(tetravertcoltbl);
    [~, indstruct] = generateSubspace(tetravertcrosstbl, tetravertcoltbl, tetravertcolfds);
    oppfacenormals = tblmap1to2(oppfacenormals, indstruct);


    %% Computation of the approximate of the gradient in the tetrahedra
    
    % For each vertex, it is given by the normal of the opposite face (let us denote
    % it N) divided by a coefficient equal to the scalar product of N with (V -
    %  FV), where V is the vector pointing to the vertex and FV is the vector
    %  pointing to the centroid of the opposite face.
    
    % We compute the centroid of the opposite face
    reducemap = setupTableMapping(tetravertoppcoltbl, tetravertcoltbl, tetravertcolfds);
    tetraoppfacecent = 1/3*(reducemap*tetraoppvertcent);

    map = setupTableMapping(vertcoltbl, tetravertcoltbl, {'vertices', 'coldim'});
    tetravertcent = map*vertcent;

    coef = (oppfacenormals.*(tetravertcent - tetraoppfacecent));
    reducemap = setupTableMapping(tetravertcoltbl, tetraverttbl, tetravertfds);
    coef = reducemap*coef;
    coef = 1./coef;

    dispatchmap = setupTableMapping(tetraverttbl, tetravertcoltbl, tetravertfds);
    coef = dispatchmap*coef;

    grad = coef.*oppfacenormals;


    %% We set up the permeability tensor

    % The permeability tensor is stored in perm which is indexed according to the
    % indexing table cellcolrowtbl
    rowtbl.rowdim = (1 : 3)';
    rowtbl.num = 3;

    [~, colrowtbl] = setupTableMapping(coltbl, rowtbl, {});
    [~, cellcolrowtbl] = setupTableMapping(celltbl, colrowtbl, []);
    cellcolrowtbl = sortTable(cellcolrowtbl, {'cells', 'coldim', 'rowdim'});

    permmat = permTensor(rock, G.griddim);
    permmat = permmat(celltbl.cells, :);
    perm = reshape(permmat', [], 1);
    
    
    % We dispatch perm to each vertex of each tetrahedra and multiply it with
    % grad. We obtain Kgrad, which is indexed according to indexing table
    % tetravertcoltbl
    
    [~, tetravertcolrowtbl] = setupTableMapping(tetravertcoltbl, rowtbl, []);
    dispatchmap = setupTableMapping(cellcolrowtbl, tetravertcolrowtbl, {'cells', ...
                        'coldim', 'rowdim'});
    K = dispatchmap*perm;

    %% Setup of Amat
    % The tensor Amat is stored according indexing table tetravert2tbl. For each
    % tetrahedra, it is a two-dimensional tensor for pairs of values at the
    % vertices of the tetrahedra.
    
    dispatchmap = setupTableMapping(tetravertcoltbl, tetravertcolrowtbl, tetravertcolfds);
    Kgrad = K.*(dispatchmap*grad);
    tetravertrowtbl = replacefield(tetravertcoltbl, {'coldim', 'rowdim'});
    reducemap = setupTableMapping(tetravertcolrowtbl, tetravertrowtbl, {tetravertfds{:}, ...
                        'rowdim'});
    Kgrad = reducemap*Kgrad;

    
    [~, tetravert2tbl] = setupTableMapping(tetraverttbl, tetraverttbl, {'cells', ...
                        'faces', 'edges'}, 'crossextend', {{'vertices', {'vertices1', ...
                        'vertices2'}}});
    tetravert2fds = {'cells', 'faces', 'edges', 'vertices1', 'vertices2'};
    [~, tetravert2coltbl] = setupTableMapping(coltbl, tetravert2tbl, []);
    tetravert2colfds = {'cells', 'faces', 'edges', 'vertices1', 'vertices2', 'coldim'};

    dispatchmap = setupTableMapping(tetravertcoltbl, tetravert2coltbl, {'cells', ...
                        'faces', 'edges', 'coldim', {'vertices', 'vertices1'}});
    Kgrad = dispatchmap*Kgrad;
    dispatchmap = setupTableMapping(tetravertcoltbl, tetravert2coltbl, {'cells', ...
                        'faces', 'edges', 'coldim', {'vertices', 'vertices2'}});
    grad = dispatchmap*grad;
    Amat = Kgrad.*grad;

    reducemapping = setupTableMapping(tetravert2coltbl, tetravert2tbl, tetravert2fds);
    Amat = reducemapping*Amat;

    %% Computation of the volume of each tetrahedra
    
    % Let us denote by P1, P2, P3 and P4 the vertices of a tetrahedra, we set u1 =
    % P2 - P1, u2 = P3 - P1 and u3 = P4 - P1. The volume of the tetrahedra is
    % given by 1/6 multiplied with the determinant of [u1, u2, u3].
    
    A = convertTableToArray(tetraverttbl, {'cells', 'faces', 'edges', 'vertices'});
    A = sortrows(A);

    clear loctbl
    loctbl.locind = [1; 2; 3; 4];
    loctbl.num = 4;

    locinds = repmat(loctbl.locind, tetratbl.num, 1);

    A = [A, locinds];

    tetraverttbl2 = convertArrayToTable(A, {'cells', 'faces', 'edges', 'vertices', ...
                        'locind'});
    [~, tetravertcoltbl2] = setupTableMapping(tetraverttbl2, coltbl, []);

    map = setupTableMapping(vertcoltbl, tetravertcoltbl2, {'vertices', 'coldim'});
    tetravertcent = map*vertcent;

    mult1 = [-1; 1; 0; 0];
    dispatchmap = setupTableMapping(loctbl, tetravertcoltbl2, {'locind'});
    mult1 = dispatchmap*mult1;
    u1 = mult1.*tetravertcent;
    reducemap = setupTableMapping(tetravertcoltbl2, tetracoltbl, {'cells', ...
                        'faces', 'edges', 'coldim'});
    u1 = reducemap*u1;

    mult2 = [-1; 0; 1; 0];
    mult2 = dispatchmap*mult2;
    u2 = mult2.*tetravertcent;
    u2 = reducemap*u2;

    mult3 = [-1; 0; 0; 1];
    mult3 = dispatchmap*mult3;
    u3 = mult3.*tetravertcent;
    u3 = reducemap*u3;

    determinanttbl.dim1 = [1; 1; 2; 2; 3; 3];
    determinanttbl.dim2 = [2; 3; 1; 3; 1; 2];
    determinanttbl.dim3 = [3; 2; 3; 1; 2; 1];
    determinanttbl.num = numel(determinanttbl.dim1);
    [~, tetradettbl] = setupTableMapping(tetratbl, determinanttbl, []);
    detmultiplier = [1; -1; -1; 1; 1; -1];
    map1 = setupTableMapping(tetracoltbl, tetradettbl, {tetrafds{:}, {'coldim', 'dim1'}});
    map2 = setupTableMapping(tetracoltbl, tetradettbl, {tetrafds{:}, {'coldim', 'dim2'}});
    map3 = setupTableMapping(tetracoltbl, tetradettbl, {tetrafds{:}, {'coldim', 'dim3'}});
    map4 = setupTableMapping(determinanttbl, tetradettbl, {'dim1', 'dim2', ...
                        'dim3'});
    vol = (map4*detmultiplier).*(map1*u1).*(map2*u2).*(map3*u3);
    reducemap = setupTableMapping(tetradettbl, tetratbl, tetrafds);
    vol = reducemap*vol;
    vol = 1/6*abs(vol);

    % We multiply the tensor Amat by the volume
    
    dispatchmap = setupTableMapping(tetratbl, tetravert2tbl, tetrafds);
    vol = dispatchmap*vol;
    
    Amat = vol.*Amat;


    %% Summation of contribution over the tetrahedra in a cell
    
    % We create indexing table cellverttbl which contains all the vertices for a
    % cell.
    cellverttbl = projTable(tetraverttbl, {'cells', 'vertices'});

    % We create table cellvert2tbl which contains the two-dimensional tensor
    % acting on pair of vertices 
    [~, cellvert2tbl] = setupTableMapping(cellverttbl, cellverttbl, {'cells'}, 'crossextend', {{'vertices', ...
                        {'vertices1', 'vertices2' }}});

    reducemapping = setupTableMapping(tetravert2tbl, cellvert2tbl, {'cells', 'vertices1', ...
                        'vertices2'});
    % The tensor Amat is now indexed according to cellvert2tbl:
    Amat = reducemapping*Amat;

    %% Barycentric condensation

    % setup cellcnverttbl: For each cell, the vertices of the cell and of the nodes
    % that belong to the cell.
    cellnodetbl.cells = rldecode((1 : nc)', diff(G.cells.nodePos));
    cellnodetbl.nodes = G.cells.nodes;
    cellnodetbl.num   = numel(cellnodetbl.cells);

    [~, cellnodetbl] = setupTableMapping(cellnodetbl, nodetbl, {'nodes'});
    a1 = convertTableToArray(cellnodetbl, {'cells', 'vertices'});
    a1 = unique(a1(:, [1, 2]), 'rows'); % all the vertices from the nodes for
                                        % given cell
    a2 = convertTableToArray(celltbl, {'cells', 'vertices'}); % vertices of the
                                                              % cell for given
                                                              % cells
    a = [a1; a2];
    cellcnverttbl = convertArrayToTable(a, {'cells', 'cnvertices'});
    cellcnvertfds = {'cells', 'cnvertices'};

    % We setup cellvertcnvert: This table is used to describe the barycentric
    % condensation from vertices to cnvertices. It is used to store a 2D tensor
    % of pair of values at cnvertices and vertices, which we call the
    % barycentric mapping.
    [~, cellvertcnverttbl] = setupTableMapping(cellverttbl, cellcnverttbl, ...
                                                            {'cells'});
    cellvertcnvertfds = {'cells', 'cnvertices', 'vertices'};

    % We store in barcoef the barycentric mapping
    vertind   = cellvertcnverttbl.vertices;
    cnvertind = cellvertcnverttbl.cnvertices;
    barcoef = zeros(cellvertcnverttbl.num, 1);
    % For the cnvertices the mapping is just the identity
    barcoef(vertind == cnvertind) = 1; % todo: barcoef is sparse. can we use a different approach?

    % We setup nodecoef: for a face, it is equal to 1/(number of nodes that belongs to the face)

    facenodetbl.faces = rldecode((1 : nf)', diff(G.faces.nodePos));
    facenodetbl.nodes = G.faces.nodes;
    facenodetbl.num = numel(facenodetbl.faces);
    nnodeperface = diff(G.faces.nodePos);
    % Here we assume that facetbl.faces = (1 : nf)'
    dispatchmap = setupTableMapping(facetbl, facenodetbl, {'faces'});
    nodecoef = 1./(dispatchmap*nnodeperface);

    [~, facenodetbl] = setupTableMapping(facenodetbl, nodetbl, {'nodes'});
    facenodetbl = replacefield(facenodetbl, {'vertices', 'nodevertices'});
    [~, facenodetbl] = setupTableMapping(facenodetbl, facetbl, {'faces'});
    facenodetbl = replacefield(facenodetbl, {'vertices', 'facevertices'});

    facenodefds = {'faces', 'nodes', 'nodevertices', 'facevertices'};

    [cellfacenodetbl, indstruct] = generateSubspace(cellfacetbl, facenodetbl, {'faces'});
    nodecoef = tbldispatch2(nodecoef, indstruct);
    % dispatchmap = setupTableMapping(facenodetbl, cellfacenodetbl, {'faces', 'nodes'});
    % nodecoef = dispatchmap*nodecoef;
    
    cellfacenodefds = {'cells', 'faces', 'nodes', 'nodevertices', 'facevertices'};

    map = setupTableMapping(cellfacenodetbl, cellvertcnverttbl, {'cells', ...
                        {'facevertices', 'vertices'}, {'nodevertices', ...
                        'cnvertices'}});
    barcoef = barcoef + map*nodecoef;

    % We apply barcoef to vertices1 and vertices2 in the Amat matrix.
    
    prod = TensorProd();
    prod.tbl1 = cellvertcnverttbl;
    prod.tbl2 = cellvert2tbl;
    prod.replacefds1 = {{'vertices', 'vertices1'}, {'cnvertices', ...
                        'cnvertices1'}};
    prod.reducefds   = {'vertices1'};
    prod.mergefds    = {'cells'};
    prod = prod.setup();
    
    Amat = prod.evalProd(barcoef, Amat);
    
    cellvert2cnvert1tbl = prod.prodtbl;

    prod = TensorProd();
    prod.tbl1 = cellvert2cnvert1tbl;
    prod.tbl2 = cellvertcnverttbl;
    prod.replacefds2 = {{'vertices', 'vertices2'}, {'cnvertices', ...
                        'cnvertices2'}};
    prod.reducefds   = {'vertices2'};
    prod.mergefds    = {'cells'};
    prod = prod.setup();
    
    Amat = prod.evalProd(Amat, barcoef);
    
    cellcnvert2tbl = prod.prodtbl;
    

    %% Finally, we collect the elements for only the node vertices

    [~, cellnverttbl] = setupTableMapping(cellcnverttbl, nodetbl, {{'cnvertices', ...
                        'vertices'}});
    cellnverttbl = replacefield(cellnverttbl, {'cnvertices', 'nvertices'});

    clear crossextends
    crossextends{1} = {'nvertices', {'nvertices1', 'nvertices2'}};
    crossextends{2} = {'nodes', {'nodes1', 'nodes2'}};
    [~, cellnvert2tbl] = setupTableMapping(cellnverttbl, cellnverttbl, {'cells'}, ...
                                                         'crossextend', crossextends);

    map = setupTableMapping(cellcnvert2tbl, cellnvert2tbl, {'cells', {'cnvertices1', ...
                        'nvertices1'}, {'cnvertices2', 'nvertices2'}});
    A = map*Amat;

    cellnode2tbl = cellnvert2tbl;
    cellnode2tbl = rmfield(cellnode2tbl, 'nvertices1');
    cellnode2tbl = rmfield(cellnode2tbl, 'nvertices2');
    
    cellnodetbl = cellnverttbl;
    cellnodetbl = rmfield(cellnodetbl, 'nvertices');
    
    vagstruct = struct('A'           , A           , ...
                       'cellnode2tbl', cellnode2tbl, ...
                       'cellnodetbl' , cellnodetbl);
    
end
