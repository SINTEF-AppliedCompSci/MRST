function vagstruct = computeVagTrans2(G, rock)
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
%
% Reference :
%  title={Benchmark 3D: the VAG scheme},
%  author={Eymard, Robert and Guichard, Cindy and Herbin, Raphaele},
%  booktitle={Finite Volumes for Complex Applications VI Problems \& Perspectives},
%  pages={1013--1022},
%  year={2011},
%  publisher={Springer}


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
    
    vagstruct = [];
    
    nc = G.cells.num;
    nf = G.faces.num;
    nn = G.nodes.num;

    celltbl.cells    = (1 : nc)';
    celltbl.vertices = (1 : nc)';
    celltbl = IndexArray(celltbl);

    facetbl.faces    = (1 : nf)';
    facetbl.vertices = nc + (1 : nf)';
    facetbl = IndexArray(facetbl);

    nodetbl.nodes    = (1 : nn)';
    nodetbl.vertices = nc + nf + (1 : nn)';
    nodetbl = IndexArray(nodetbl);

    np = nc + nf + nn;
    verttbl.vertices = (1 : np)';
    verttbl = IndexArray(verttbl);
    vertfds = {'vertices'};
    
    %% Setup vertcoltbl and vertcent
    
    % vertcent contains the coordinates of the vertices and is indexed using
    % the indexing table vertcoltbl
     
    coltbl.coldim   = (1 : 3)';
    coltbl = IndexArray(coltbl);


    % We collect the cell, face and node centroids given by the grid structure in
    % vectors that follows the indexing of the indexing tables cellcoltbl,
    % facecoltbl and nodecoltbl.
    cellcent = G.cells.centroids;
    cellcent = reshape(cellcent', [], 1);
    facecent = G.faces.centroids;
    facecent = reshape(facecent', [], 1);
    nodecent = G.nodes.coords;
    nodecent = reshape(nodecent', [], 1);

    vertcoltbl = crossIndexArray(verttbl, coltbl, [], 'optpureproduct', true);
    cellcoltbl = crossIndexArray(celltbl, coltbl, [], 'optpureproduct', true);
    facecoltbl = crossIndexArray(facetbl, coltbl, [], 'optpureproduct', true);
    nodecoltbl = crossIndexArray(nodetbl, coltbl, [], 'optpureproduct', true);

    
    map = TensorMap();
    map.fromTbl = cellcoltbl;
    map.toTbl = vertcoltbl;
    map.mergefds = {'vertices', 'coldim'};
    map = map.setup();
    vertcent = map.eval(cellcent);
    
    map = TensorMap();
    map.fromTbl = facecoltbl;
    map.toTbl = vertcoltbl;
    map.mergefds = {'vertices', 'coldim'};
    map = map.setup();
    vertcent = vertcent + map.eval(facecent);
    
    map = TensorMap();
    map.fromTbl = nodecoltbl;
    map.toTbl = vertcoltbl;
    map.mergefds = {'vertices', 'coldim'};
    map = map.setup();
    vertcent = vertcent + map.eval(nodecent);
        

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
    cellfacetbl = IndexArray(cellfacetbl);

    nf = G.faces.num;
    clear faceedgetbl;
    if ~isfield(G.faces, 'edgePos')
        G = computeVEMGeometry(G);
    end
    faceedgetbl.faces = rldecode((1 : nf)', diff(G.faces.edgePos));
    faceedgetbl.edges = G.faces.edges;
    faceedgetbl = IndexArray(faceedgetbl);

    cellfaceedgetbl = crossIndexArray(cellfacetbl, faceedgetbl, {'faces'});
    
    
    %% Setup tetraverttbl
    
    % We set the table tetraverttbl which, for each tetrahedra, gives all the three vertices
    % given by the vertex of the face and the two vertices of the edge.
    % 
    % A tetrahedra is given by triplet (cell, face, edge).
    %

    tetratbl = cellfaceedgetbl;
    tetratbl = sortIndexArray(tetratbl, {'cells', 'faces', 'edges'});
    tetrafds = {'cells', 'faces', 'edges'};

    tetracoltbl = crossIndexArray(tetratbl, coltbl, [], 'optpureproduct', true);
    tetracolfds = {'cells', 'faces', 'edges', 'coldim'};

    ne = G.edges.num;
    clear edgetbl;
    edgetbl.edges = (1 : ne)';
    edgetbl.nodes = G.edges.nodes(G.edges.nodePos(1 : (end - 1)));
    edgetbl = IndexArray(edgetbl);

    tetratbl2 = crossIndexArray(tetratbl, edgetbl, {'edges'});
    tetratbl2 = replacefield(tetratbl2, {'nodes', 'nodes1'});

    clear edgetbl;
    edgetbl.edges = (1 : ne)';
    edgetbl.nodes = G.edges.nodes(G.edges.nodePos(1 : (end - 1)) + 1);
    edgetbl = IndexArray(edgetbl);

    tetratbl2 = crossIndexArray(tetratbl2, edgetbl, {'edges'});
    tetratbl2 = replacefield(tetratbl2, {'nodes', 'nodes2'});
    
    exttetrafds = {'cells', 'faces', 'edges', 'nodes1', 'nodes2'};
    
    tetraverttbls = cell(3, 1);
    
    tetraverttbls{1} = crossIndexArray(tetratbl2, facetbl, {'faces'});
    tetraverttbls{2} = crossIndexArray(tetratbl2, nodetbl, {{'nodes1', ...
                        'nodes'}});
    tetraverttbls{3} = crossIndexArray(tetratbl2, nodetbl, {{'nodes2', ...
                        'nodes'}});
    
    tetravertinds = [];
    for i = 1 : 3
        % note that we assume given row ordering of the indices
        % {'cells', 'faces', 'edges', 'nodes1', 'nodes2', 'vertices'};
        tetravertinds = vertcat(tetravertinds, tetraverttbls{i}.inds);
    end
    tetravertfdnames = {'cells', 'faces', 'edges', 'nodes1', 'nodes2', 'vertices'};

    tetraverttbl = IndexArray([]);
    tetraverttbl = tetraverttbl.setup(tetravertfdnames, tetravertinds);
    tetraverttbl = tetraverttbl.removefield({'nodes1', 'nodes2'});
    
    tetraverttbl = sortIndexArray(tetraverttbl, {'cells', 'faces', 'edges', 'vertices'});
    
    tetravertcoltbl = crossIndexArray(tetraverttbl, coltbl, [], 'optpureproduct', ...
                                      true);

    tetravertcolfds = {'cells', 'faces', 'edges', 'vertices', 'coldim'};

    map = TensorMap();
    map.fromTbl = vertcoltbl;
    map.toTbl = tetravertcoltbl;
    map.mergefds = {'vertices', 'coldim'};
    map = map.setup();
    
    tetravertcent = map.eval(vertcent);
    
    map = TensorMap();
    cellcoltbl2 = replacefield(cellcoltbl, {{'vertices', ''}});
    map.fromTbl = cellcoltbl2;
    map.toTbl = tetravertcoltbl;
    map.mergefds = {'cells', 'coldim'};
    map = map.setup();
    
    tetracellcent = map.eval(cellcent);
    
    tetravect = tetravertcent - tetracellcent;

    dotest = true;
    if dotest
        onetetratbl.cells = 1;
        onetetratbl.faces = 1;
        onetetratbl.edges = 2;
        onetetratbl = IndexArray(onetetratbl);

        onetetraverttbl = crossIndexArray(onetetratbl, tetraverttbl, {'cells', 'faces', 'edges'});
        onetetraverttbl = sortIndexArray(onetetraverttbl, {'cells', 'faces', 'edges', 'vertices'});
        onetetravertcoltbl = crossIndexArray(onetetraverttbl, coltbl, {}, 'optpureproduct', true);
        
        map = TensorMap();
        map.fromTbl = tetravertcoltbl;
        map.toTbl = onetetravertcoltbl;
        map.mergefds = {'cells', 'faces', 'edges', 'vertices', 'coldim'};
        map = map.setup();
        
        onetetravect = map.eval(tetravect);
        onetetravertcent = map.eval(tetravertcent);
        onetetracellcent = map.eval(tetracellcent);
    
    end
    
    % prepare for the blockwise inversion
    map = TensorMap();
    map.fromTbl = tetraverttbl;
    map.toTbl = tetravertcoltbl;
    map.mergefds = {'cells', 'faces', 'edges', 'vertices'};
    ind1 = getDispatchInd(map);
    
    map = TensorMap();
    map.fromTbl = tetracoltbl;
    map.toTbl = tetravertcoltbl;
    map.mergefds = {'cells', 'faces', 'edges', 'coldim'};
    ind2 = getDispatchInd(map);    
    
    A = sparse(ind1, ind2, tetravect, tetraverttbl.num, tetracoltbl.num);
    
    bi = @invertDiagonalBlocksMex;
    sz = repmat(coltbl.num, tetratbl.num, 1);
    invA = bi(A, sz);
    
    return
    

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


    % We setup a 4D tensor barcoef2, which is the tensor product of the 2D
    % tensor barcoef.
    crossextend1 = {'vertices'  , {'vertices1'  , 'vertices2'}};
    crossextend2 = {'cnvertices', {'cnvertices1', 'cnvertices2'}};
    crossextend = {crossextend1, crossextend2};
    [~, cellvert2cnvert2tbl] = setupTableMapping(cellvertcnverttbl, ...
                                                 cellvertcnverttbl, {'cells'}, ...
                                                 'crossextend', crossextend);

    dispatchmap = setupTableMapping(cellvert2tbl, cellvert2cnvert2tbl, {'cells', 'vertices1', ...
                        'vertices2'});
    Amat = dispatchmap*Amat;

    dispatchmap1 = setupTableMapping(cellvertcnverttbl, cellvert2cnvert2tbl, {'cells', {'vertices', 'vertices1'}, ...
                        {'cnvertices', 'cnvertices1'}});
    barcoef1 = dispatchmap1*barcoef;
    dispatchmap2 = setupTableMapping(cellvertcnverttbl, cellvert2cnvert2tbl, {'cells', {'vertices', 'vertices2'}, ...
                        {'cnvertices', 'cnvertices2'}});
    barcoef2 = dispatchmap2*barcoef;

    Amat = Amat.*barcoef1.*barcoef2;

    [~, cellcnvert2tbl] = setupTableMapping(cellcnverttbl, cellcnverttbl, {'cells'}, 'crossextend', ...
                                                          {{'cnvertices', ...
                        {'cnvertices1', 'cnvertices2'}}});
    reducemap = setupTableMapping(cellvert2cnvert2tbl, cellcnvert2tbl, {'cells', ...
                        'cnvertices1', 'cnvertices2'});

    Amat = reducemap*Amat;


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

function iA = invertDiagonalBlocksMex(A, sz)
   sz     = int32(sz);
   blocks = matrixBlocksFromSparse(A, sz);
   iA     = blockDiagMatrix(invv(blocks, sz), sz);
end

function A = blockDiagMatrix(V, sz)
   [I, J] = blockDiagIndex(sz);
   A      = sparse(I, J, V);
end
