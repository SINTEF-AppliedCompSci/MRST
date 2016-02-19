function dofVec = VEM2D_loc_v2(G, K, k)

nodeNum = G.cells.nodePos(K):G.cells.nodePos(K+1)-1;
nodes = G.cells.nodes(nodeNum);
if size(nodes,1) == 1
    nodes = nodes';
end
nN = size(nodes,1);

if k ==1
    dofVec = nodes';

elseif k == 2
                            %   Cell edges and edge midpoint coordinates.
    edgeNum = G.cells.facePos(K):G.cells.facePos(K+1)-1;
    edges = G.cells.edges(edgeNum);
    if size(egdes,1) == 1
        edges = edges';
    end
    nE = size(edges,1);
    dofVec = [nodes', edges' + nN, K + nN + nE];
end

end