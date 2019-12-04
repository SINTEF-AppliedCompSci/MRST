function relation = handleMatchingFaces(G1, cells1, bdnodes1, G2)
% Compute intersection relation between layered boundaries of subgrids. 
% Basically, the faces on the layered boundaries are matching, and only the
% common areas are obtained from the cells and boundary nodes of G1. This
% function is tailored to grids of the near-wellbore model.
%
% SYNOPSIS:
%   relation = handleMatchingFaces(G1, cells1, bdnodes1, G2)
%
% PARAMETERS:
%  G1,G2     - Layered grid structures, G2 is loacted inside G1
%  cells1    - Cells in G1 which will be replaced by G2
%  bdnodes1  - Boundary nodes of cells1
%
% RETURNS:
% relation - Face intersection relation, n x 3 matrix
%            column 1     - Face of G1
%            column 2     - Face of G2
%            column 3     - Areas of intersection subfaces
%
% SEE ALSO:
%   `assembleGrids`, `handleNonMatchingFaces`

    g2 = G2.surfGrid;

    fun_faces = @(G, c)G.cells.faces...
        (G.cells.facePos(c) : G.cells.facePos(c+1)-1, 1);
    fun_nodes = @(G, f)G.faces.nodes...
        (G.faces.nodePos(f) : G.faces.nodePos(f+1)-1, 1);

    % Find boundary faces 2, 2D
    faces2 = find(~all(g2.faces.neighbors, 2));
    nodes2 = arrayfun(@(f)fun_nodes(g2,f), faces2, 'UniformOutput', false);

    bdnodes2 = g2.nodes.boundary;
    bdnodesTwoPts = [bdnodes2; bdnodes2(1)];

    bdfaces2 = zeros(length(bdnodes2), 1);
    for i = 1 : length(bdfaces2)
        idx = cellfun(@(n)all( ismember(n, bdnodesTwoPts(i:i+1, :)) ), nodes2);
        bdfaces2(i) = faces2(idx);
    end

    % Find boundary faces 1
    nRef = G2.layers.refinement;
    layerIdx = cumsum([1; nRef]);
    nxyfaces = nnz(G2.faces.surfaces==0) / G2.layers.num;

    relation = cell(length(cells1), 1);
    for k = 1 : length(cells1)
        faces1 = arrayfun(@(c)fun_faces(G1,c), cells1{k}, 'UniformOutput', false);
        tab = tabulate_NWM(cell2mat(faces1));
        faces1 = tab( tab(:,2) == 1, 1);
        nodes1 = arrayfun(@(f)fun_nodes(G1,f), faces1, 'UniformOutput', false);
        bdnFourPts = [bdnodes1{k}, bdnodes1{k+1}];
        bdnFourPts = [bdnFourPts; bdnFourPts(1,:)];
        ca0 = cell(length(bdnodes1{k}), 1);
        for i = 1 : length(ca0)
            idx = cellfun(@(n)all( ismember(n, bdnFourPts(i:i+1, :)) ), nodes1);
            bdface1 = faces1(idx);
            faces1_i = bdface1;
            faces2_i = bdfaces2(i) + ( (layerIdx(k) : layerIdx(k+1)-1) - 1)' * nxyfaces;
            areas_i  = G2.faces.areas( faces2_i );
            ca0{i} = [faces1_i * ones(size(faces2_i)), faces2_i, areas_i];
        end
        relation{k} = cell2mat(ca0);
    end
    relation = cell2mat(relation);
end