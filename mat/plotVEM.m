function plotVEM(G, u, type)

Nc = G.cells.num;
Nn = G.nodes.num;

if strcmp(type, 'dof')
    Xb = zeros(Nc, 2);
end

hold on
for c = 1:Nc
    nodeNum = G.cells.nodePos(c):G.cells.nodePos(c+1)-1;        
    nodes = G.cells.nodes(nodeNum);
    X = G.nodes.coords(nodes,:);
    if strcmp(type, 'dof')
        [~, Xb(c,:)] = baric(X);
        plot3(X([1:end, 1], 1), X([1:end, 1],2), u(nodes([1:end, 1])), 'k');
    else
        faceNum = G.cells.facePos(c):G.cells.facePos(c+1) -1;
        faces = G.cells.faces(faceNum);
        Xc = G.faces.centroids(faces,:);
        XdofC = zeros(2*size(X,1) + 1,2);
        XdofC(1:2:end,:) = X([1:end,1], :);
        XdofC(2:2:end-1,:) = Xc;
        dofs = zeros(size(XdofC,1),1);
        if size(nodes,1) == 1
            nodes = nodes';
        end
        if size(faces,1) == 1
            faces = faces';
        end
        dofs(1:2:end) = nodes([1:end, 1]);
        dofs(2:2:end-1,:) = faces + Nn;
        fill3(XdofC(:, 1), XdofC(:,2), u(dofs), 'w');        
        %plot3(XdofC(:, 1), XdofC(:,2), u(dofs), 'k', 'LineWidth', 0.01);
    end
end

if strcmp(type, 'dof')
    X = G.nodes.coords;
    Xc = G.faces.centroids;
    Xdof = [X; Xc; Xb];
    plot3(Xdof(:,1), Xdof(:,2), u, '.b')
end

hold off
view(3)
axis([0 1 0 1 min(u), max(u)]);
fontSize = 18;
xlabel('x', 'FontSize', fontSize); ylabel('y', 'FontSize', fontSize); zlabel('u_h(x,y)', 'FontSize', fontSize);
end