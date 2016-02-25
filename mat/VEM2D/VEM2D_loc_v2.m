function dofVec = VEM2D_loc_v2(G, K, k)

edgeNum = G.cells.facePos(K):G.cells.facePos(K+1)-1;
edges = G.cells.faces(edgeNum);
if size(edges,1) == 1
    edges = edges';
end
nE = size(edges,1);
edgeSign = (-ones(nE,1)).^(G.faces.neighbors(edges,1) ~= K); 

nodeNum = mcolon(G.faces.nodePos(edges),G.faces.nodePos(edges+1)-1);
nodes = G.faces.nodes(nodeNum);
if size(nodes,1) == 1
    nodes = nodes';
end
nodes = reshape(nodes,2,[])';
nodes(edgeSign == -1,:) ...
        = nodes(edgeSign == -1,2:-1:1);
nodes = nodes(:,1);

if k ==1
    dofVec = nodes';

elseif k == 2
                            %   Cell edges and edge midpoint coordinates.   
    dofVec = [nodes', edges' + G.nodes.num, K + G.nodes.num + G.faces.num];
end

end