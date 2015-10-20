function [nodes, points] = nodeData(G,K)
    
nodeNum = G.cells.nodePos(K):G.cells.nodePos(K+1)-1;
nodes = G.cells.nodes(nodeNum);
points = G.nodes.coords(nodes,:);

end
