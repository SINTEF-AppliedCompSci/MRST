function plotVEM_v2(G,U)

    nodeNum = mcolon(G.faces.nodePos(1:end-1),G.faces.nodePos(2:end)-1);
    nodes = G.faces.nodes(nodeNum);
    X = G.nodes.coords(nodes,:);

    U = U(nodes);
    
    plot3(X(:,1),X(:,2),U);
end
    