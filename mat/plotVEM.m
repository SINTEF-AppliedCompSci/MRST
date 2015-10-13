function plotVEM(G, u, type)

Nc = G.cells.num;

if strcmp(type, 'dof')
    baricenters = zeros(Nc,2);
    for c = 1:Nc
        nodeNum = G.cells.nodePos(c):G.cells.nodePos(c+1)-1;        
        nodes = G.cells.nodes(nodeNum);
        X = G.nodes.coords(nodes,:);
        [~, baricenters(c,:)] = baric(X);
    end
    X = [G.nodes.coords ; G.faces.centroids ; baricenters];
    plot3(X(:,1), X(:,2), u, 'ob')
end

hold on
for c = 1:Nc
    nodeNum = G.cells.nodePos(c):G.cells.nodePos(c+1)-1;
    nodes = G.cells.nodes(nodeNum);
    X = G.nodes.coords(nodes,:);
    if strcmp(type, 'dof')
        plot3(X([1:end, 1], 1), X([1:end, 1],2), u(nodes([1:end, 1])), 'k');
    else
        fill3(X([1:end, 1], 1), X([1:end, 1],2), u(nodes([1:end, 1])), u(nodes([1:end, 1])));
    end
end
hold off
view(3)
end