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
%
% Reference :
%  title={Benchmark 3D: the VAG scheme},
%  author={Eymard, Robert and Guichard, Cindy and Herbin, Raphaele},
%  booktitle={Finite Volumes for Complex Applications VI Problems \& Perspectives},
%  pages={1013--1022},
%  year={2011},
%  publisher={Springer}


%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

    
    %% Setup cell, face, node and vertices IndexArray
        
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
    rowtbl = replacefield(coltbl, {{'coldim', 'rowdim'}});
    colrowtbl = crossIndexArray(coltbl, rowtbl, {}, 'optpureproduct', true);
    
    celltbl2 = replacefield(celltbl, {{'vertices', ''}});
    cellcoltbl2 = crossIndexArray(celltbl2, coltbl, {}, 'optpureproduct', true);
    cellcolrowtbl = crossIndexArray(celltbl2, colrowtbl, {}, 'optpureproduct', true);
    
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

    nc = G.cells.num;
    cellfacetbl.cells = rldecode((1 : nc)', diff(G.cells.facePos));
    cellfacetbl.faces = G.cells.faces(:, 1);
    cellfacetbl = IndexArray(cellfacetbl);

    nf = G.faces.num;
    facenodetbl.faces = rldecode((1 : nf)', diff(G.faces.nodePos));
    facenodetbl.nodes = G.faces.nodes;
    facenodetbl = IndexArray(facenodetbl); 

    dorecomputefacecentroids = true;
    if dorecomputefacecentroids
        % We recompute the face centroids
        nnodeperface = diff(G.faces.nodePos);
        
        map = TensorMap();
        map.fromTbl = facetbl;
        map.toTbl = facenodetbl;
        map.mergefds = {'faces'};
        map = map.setup();
        nodecoef = 1./map.eval(nnodeperface);
        
        prod = TensorProd();
        prod.tbl1 = facenodetbl;
        nodecoltbl2 = replacefield(nodecoltbl, {{'vertices', ''}});
        facecoltbl2 = replacefield(facecoltbl, {{'vertices', ''}});
        prod.tbl2 = nodecoltbl2;
        prod.tbl3 = facecoltbl2;
        prod.mergefds = {'nodes'};
        prod = prod.setup();
        
        facecent = prod.eval(nodecoef, nodecent);
    end
    
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
    


    cellfacenodetbl = crossIndexArray(cellfacetbl, facenodetbl, {'faces'});
    
    cellnodetbl = projIndexArray(cellfacenodetbl, {'cells', 'nodes'});
    cellnodetbl = sortIndexArray(cellnodetbl, {'cells', 'nodes'});
    
    cellnode2tbl = crossIndexArray(cellnodetbl, cellnodetbl, {'cells'}, 'crossextend', {{'nodes', {'nodes1', 'nodes2'}}});
    
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
    tetraverttbl = tetraverttbl.removeInd({'nodes1', 'nodes2'});

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
    map.fromTbl = cellcoltbl2;
    map.toTbl = tetravertcoltbl;
    map.mergefds = {'cells', 'coldim'};
    map = map.setup();
    
    tetracellcent = map.eval(cellcent);

    tetravect = tetravertcent - tetracellcent;
    
    dotest = false;
    if dotest
        % investigate one tetrahedra
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
    
    % Prepare for the blockwise inversion
    prod = TensorProd();
    prod.tbl1 = tetravertcoltbl;
    prod.tbl2 = tetracoltbl;
    prod.tbl3 = tetraverttbl;
    prod.mergefds = {'cells', 'faces', 'edges'};
    prod.reducefds = {'coldim'};
    
    [ind1, ind2] = prod.getDispatchInd();
    
    tv_num = tetraverttbl.num;
    tc_num = tetracoltbl.num;
    A = sparse(ind1, ind2, tetravect, tetraverttbl.num, tetracoltbl.num);
    
    bi = @invertDiagonalBlocksMex;
    sz = repmat(coltbl.num, tetratbl.num, 1);
    invA = bi(A, sz);
    
    ind = sub2ind([tc_num, tv_num], ind2, ind1);
    grad = invA(ind);

    dotest = false;
    if dotest
        % check if we have compute the inverse matrix of A
        rowtbl = replacefield(coltbl, {{'coldim', 'rowdim'}});
        tetracolrowtbl = crossIndexArray(tetracoltbl, rowtbl, {}, 'optpureproduct', true);
    prod = TensorProd();
    prod.tbl1 = tetravertcoltbl;
    prod.tbl2 = tetravertcoltbl;
        prod.tbl3 = tetracolrowtbl;
        prod.replacefds1 = {{'coldim', 'rowdim'}};
    prod.mergefds = {'cells', 'faces', 'edges'};
        prod.reducefds = {'vertices'};
    prod = prod.setup();
    
        id = prod.eval(grad, tetravect);
    
    prod = TensorProd();
        prod.tbl1 = tetracolrowtbl;
        prod.tbl2 = tetracoltbl;
        prod.tbl3 = tetracoltbl;
        prod.replacefds1 = {{'coldim', 'rowdim', 'interchange'}};
        prod.replacefds2 = {{'coldim', 'rowdim'}};
        prod.mergefds = {'cells', 'faces', 'edges'};
    prod.reducefds   = {'rowdim'};
    prod = prod.setup();
    
        id_T = SparseTensor();
        id_T = id_T.setFromTensorProd(id, prod);
        B = id_T.getMatrix();

        useit = id>1e-10;
        spy(B); % looks ok
    
    end

    permmat = permTensor(rock, G.griddim);
    K = reshape(permmat', [], 1);
    
    prod = TensorProd();
    prod.tbl1 = tetravertcoltbl;
    prod.tbl2 = cellcolrowtbl;
    prod.tbl3 = tetravertcoltbl;
    prod.replacefds1 = {{'coldim', 'rowdim'}};
    prod.replacefds2 = {{'coldim', 'rowdim', 'interchange'}};
    prod.mergefds = {'cells'};
    prod.reducefds = {'rowdim'};
    prod = prod.setup();
    
    Kgrad = prod.eval(grad, K);

    tetravert2tbl = crossIndexArray(tetraverttbl, tetraverttbl, {'cells', 'faces', 'edges'}, 'crossextend', {{'vertices', {'vertices1', 'vertices2'}}});
    
    prod = TensorProd();
    prod.tbl1 = tetravertcoltbl;
    prod.tbl2 = tetravertcoltbl;
    prod.tbl3 = tetravert2tbl;
    prod.replacefds1 = {{'vertices', 'vertices2'}};
    prod.replacefds2 = {{'vertices', 'vertices1'}};
    prod.mergefds = {'cells', 'faces', 'edges'};
    prod.reducefds = {'coldim'};
    prod = prod.setup();
    
    gradKgrad = prod.eval(Kgrad, grad);
    
    % We compute the volumes of the tetrahedras. tetraverttbl is sorted. We use this
    % now.
    nt = tetratbl.num;
    vertlocinds = repmat((1 : coltbl.num)', nt, 1);
    tetravertloctbl = tetraverttbl.addInd('vertloc', vertlocinds);
    tetravertloccoltbl = crossIndexArray(tetravertloctbl, coltbl, {}, 'optpureproduct', true);
    
    tetraloccoltbls = cell(coltbl.num, 1);
    for i = 1 : coltbl.num
        tetraloccoltbls{i}.vertloc = i;
        tetraloccoltbls{i} = IndexArray(tetraloccoltbls{i});
        tetraloccoltbls{i} = crossIndexArray(tetraloccoltbls{i}, tetravertloccoltbl, {'vertloc'});

        map = TensorMap();
        map.fromTbl = tetravertcoltbl;
        map.toTbl = tetraloccoltbls{i};
        map.mergefds = {'cells', 'faces', 'edges', 'vertices', 'coldim'};
        map = map.setup();
        tetravects{i} = map.eval(tetravect);
    end

    u = tetravects{1};
    v = tetravects{2};
    w = tetravects{3};
    
    % We compute cross-product of u and v
    
    crossinds = [1; 3; 2; 3; 2; 1; 2; 1; 3];
    colrowcrosstbl = colrowtbl.addInd('crossdim', crossinds);
    crossvect = [0; 1; -1; -1; 0; 1; 1; -1; 0];
    
    tetracolrowtbl = crossIndexArray(tetratbl, colrowtbl, {}, 'optpureproduct', true);
    
    prod = TensorProd();
    prod.tbl1 = colrowcrosstbl;
    prod.tbl2 = tetracoltbl;
    prod.tbl3 = tetracolrowtbl;
    prod.replacefds1 = {{'rowdim', 'ind'},  {'crossdim', 'rowdim'}};
    prod.replacefds2 = {{'coldim', 'ind'}};
    prod.reducefds = {'ind'};
    prod = prod.setup();
    
    uvcross = prod.eval(crossvect, u);
    
    prod = TensorProd();
    prod.tbl1 = tetracolrowtbl;
    prod.tbl2 = tetracoltbl;
    prod.tbl3 = tetracoltbl;
    prod.replacefds1 = {{'coldim', 'ind'},  {'rowdim', 'coldim'}};
    prod.replacefds2 = {{'coldim', 'ind'}};
    prod.mergefds = {'cells', 'faces', 'edges'};
    prod.reducefds = {'ind'};
    prod = prod.setup();
    
    uvcross = prod.eval(uvcross, v);
    
    % The volume is equal to the absolute value of the inner-product of w and
    % uvcross.
    
    prod = TensorProd();
    prod.tbl1 = tetracoltbl;
    prod.tbl2 = tetracoltbl;
    prod.tbl3 = tetratbl;
    prod.mergefds = {'cells', 'faces', 'edges'};
    prod.reducefds = {'coldim'};
    prod = prod.setup();
    
    vol = prod.eval(w, uvcross);
    vol = abs(vol);
    
    map = TensorMap();
    map.fromTbl = tetratbl;
    map.toTbl = tetravert2tbl;
    map.mergefds = {'cells', 'faces', 'edges'};
    map = map.setup();
    
    vol = map.eval(vol);
    
    gradKgrad = vol.*gradKgrad;
    
    cellverttbl = projIndexArray(tetraverttbl, {'cells', 'vertices'});
    cellvert2tbl = crossIndexArray(cellverttbl, cellverttbl, {'cells'}, 'crossextend', {{'vertices', {'vertices1', 'vertices2'}}});

    map = TensorMap();
    map.fromTbl = tetravert2tbl;
    map.toTbl = cellvert2tbl;
    map.mergefds = {'cells', 'vertices1', 'vertices2'};
    map = map.setup();
    
    gradKgrad = map.eval(gradKgrad);
    
    %% Barycentric reduction
    
    nnodeperface = diff(G.faces.nodePos);

    map = TensorMap();
    map.fromTbl = facetbl;
    map.toTbl = facenodetbl;
    map.mergefds = {'faces'};
    map = map.setup();
    nodecoef = 1./map.eval(nnodeperface);

    dotest = true;
    if dotest 
        % test consistency of face centroids
        prod = TensorProd();
        prod.tbl1 = facenodetbl;
        nodecoltbl2 = replacefield(nodecoltbl, {{'vertices', ''}});
        facecoltbl2 = replacefield(facecoltbl, {{'vertices', ''}});
        prod.tbl2 = nodecoltbl2;
        prod.tbl3 = facecoltbl2;
        prod.mergefds = {'nodes'};
    prod = prod.setup();
    
        newfacecent = prod.eval(nodecoef, nodecent);
        max(abs(newfacecent - facecent)); % they did not match exactly. We recompute them instead.
    
    end
    
    facenodeverttbl = crossIndexArray(facenodetbl, facetbl, {'faces'});
    nodeverttbl1 = projIndexArray(facenodeverttbl, {'nodes', 'vertices'});
    
    nodeverttbl2 = nodetbl;
    
    % Sanity check
    assert(all(strcmp(nodeverttbl1.fdnames, nodeverttbl2.fdnames)), 'inconsistent field names');
    inds = [nodeverttbl1.inds; nodeverttbl2.inds];
    fdnames =  nodeverttbl1.fdnames;
    nodeverttbl = IndexArray([]);
    nodeverttbl = nodeverttbl.setup(fdnames, inds);

    map = TensorMap();
    map.fromTbl = facenodeverttbl;
    map.toTbl = nodeverttbl;
    map.mergefds = {'nodes', 'vertices'};
    map = map.setup();

    barred1 = map.eval(nodecoef);

    map = TensorMap();
    map.fromTbl = nodetbl;
    map.toTbl = nodeverttbl;
    map.mergefds = {'nodes', 'vertices'};
    map = map.setup();
    
    barred2 = map.eval(ones(nodetbl.num, 1));

    barred = barred1 + barred2;

    prod = TensorProd();
    prod.tbl1 = nodeverttbl;
    prod.tbl2 = cellvert2tbl;
    prod.replacefds1 = {{'vertices', 'vertices2'}};
    prod.replacefds2 = {{'vertices1', 'vertices'}};
    prod.reducefds = {'vertices2'};
    prod = prod.setup();
    
    cellnodeverttbl = prod.tbl3;
    gradKgrad = prod.eval(barred, gradKgrad);
    
    prod = TensorProd();
    prod.tbl1 = nodeverttbl;
    prod.tbl2 = cellnodeverttbl;
    prod.tbl3 = cellnode2tbl;
    prod.replacefds1 = {{'nodes', 'nodes1'}};
    prod.replacefds2 = {{'nodes', 'nodes2'}};
    prod.reducefds = {'vertices'};
    prod = prod.setup();
    
    gradKgrad = prod.eval(barred, gradKgrad);
    
    dotest = false;
    
    if dotest
        % plot sparsity of matrix gradKgrad
        prod = TensorProd();
        prod.tbl1 = cellnode2tbl;
        prod.tbl2 = cellnodetbl;
        prod.tbl3 = cellnodetbl;
        prod.replacefds1 = {{'nodes1', 'nodes'}};
        prod.replacefds2 = {{'nodes', 'nodes2'}};
        prod.mergefds = {'cells'};
        prod.reducefds = {'nodes2'};
        
        [ind1, ind2] = prod.getDispatchInd();
        
        cn_num = cellnodetbl.num;
        A = sparse(ind1, ind2, gradKgrad, cn_num, cn_num);
        
        spy(A); 
    
    end
    
    tbls = struct('celltbl', celltbl, ...
                  'nodetbl', nodetbl, ...
                  'facetbl', facetbl, ...
                  'facenodetbl', facenodetbl, ...
                  'cellnodetbl', cellnodetbl, ...
                  'cellnode2tbl', cellnode2tbl);
    
    vagstruct = struct('A'   , gradKgrad, ...
                       'tbls', tbls);

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
