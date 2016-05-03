clc; clear all; close all;


fN = 2;

switch fN
    case 1

    addpath('../VEM2D/')
    n = 4;
    G = unitSquare([n,n],[1,1]);
    G = computeGeometry(G);
    G = mrstGridWithFullMappings(G);

    dx = 1/(n-1);
    dy = 1/(n-1);

    addFac = .15;

    intCells = find(G.cells.centroids(:,1) > dx   & ...
                    G.cells.centroids(:,1) < 1-dx & ...
                    G.cells.centroids(:,2) > dy   & ...
                    G.cells.centroids(:,2) < 1-dy       );                 

    nK = numel(intCells);
    K = intCells(round(rand(1,1)*(nK-1) + 1));
    Kc = G.cells.centroids(K,:);

    nodeNum = G.cells.nodePos(K):G.cells.nodePos(K+1)-1;
    nodes   = G.cells.nodes(nodeNum);
    X       = G.nodes.coords(nodes,:);
    nN = numel(nodes);
    
    faces   = G.cells.faces(G.cells.facePos(K):G.cells.facePos(K+1)-1);
    Ec = G.faces.centroids(faces,:);
    nE = numel(faces);

    Kc = G.cells.centroids(K,:);
    
    %%

    plotGrid(G, K, 'facealpha', .2)
    hold on;
    for i = 1:nN
        h1 = plot(X(:,1), X(:,2), 'ok', 'MarkerFaceColor', 'r', 'markersize', 6);
    end
    for i = 1:nE
        h2 = plot(Ec(:,1), Ec(:,2), 'sk', 'MarkerFaceColor', 'b', 'markersize', 6);
    end
    
    h3 = plot(Kc(:,1), Kc(:,2), 'dk', 'MarkerFaceColor', 'w', 'markersize', 6);
    
    
    h = legend([h1, h2, h3], '$\mathcal{V}^K$', '$\mathcal{E}^K$',...
                             '$\mathcal{P}^K$');
    set(h, 'interpreter', 'latex')
    axis equal off;
    

    %%

    cut = 4;
    ps = get(gcf, 'Position');
    ratio = (ps(4)-ps(2)) / (ps(3)-ps(1));
    paperWidth = 10;
    paperHeight = paperWidth*ratio - cut;
    set(gcf, 'paperunits', 'centimeters');
    set(gcf, 'papersize', [paperWidth paperHeight]);
    set(gcf, 'PaperPosition', [0    0   paperWidth paperHeight]);

    print(gcf, '-dpdf', '../../tex/thesis/fig/dofExample2D.pdf');
    
    
    case 2
    
     %%
        
    addpath('../../../pebiGridding/voronoi3D/')
    addpath('../VEM3D')
    addpath('../')

    %%
    
    close all;

    G = voronoiCube(50,[1,1,1]);

    G = computeVEM3DGeometry(G);

    intCells = find(G.cells.centroids(:,1) > .2 & ...
                    G.cells.centroids(:,1) < .8 & ...
                    G.cells.centroids(:,2) > .2 & ...
                    G.cells.centroids(:,2) < .8 & ...
                    G.cells.centroids(:,3) > .2 & ...
                    G.cells.centroids(:,3) < .8);

    nK = numel(intCells);
    K = intCells(round(rand(1,1)*(nK-1)+1));

    faces        = G.cells.faces(G.cells.facePos(K):G.cells.facePos(K+1)-1);
    faceNormals  = G.faces.normals(faces,:);
    nF = numel(faces);
    faceSigns    = (-ones(nF,1)).^(G.faces.neighbors(faces,1) ~= K);
    faceNormals  = bsxfun(@times, faceNormals,faceSigns);
    faces        = faces(sum(bsxfun(@times, faceNormals, [0,0,1]), 2) < 0);
    
    edges = G.cells.edges(G.cells.edgePos(K):G.cells.edgePos(K+1)-1);
    allEdges = edges;
    faceEdges    = G.faces.edges(mcolon(G.faces.edgePos(faces),G.faces.edgePos(faces+1)-1));
    edges = edges(ismember(edges,faceEdges));
    nodes = G.cells.nodes(G.cells.nodePos(K):G.cells.nodePos(K+1)-1);
    faceNodes    = G.faces.nodes(mcolon(G.faces.nodePos(faces),G.faces.nodePos(faces+1)-1));
    nodes = nodes(ismember(nodes,faceNodes));
    
    Kc = G.cells.centroids(K,:);

    X = G.nodes.coords(nodes,:);
    Ec = G.edges.centroids(edges,:);
    Fc = G.faces.centroids(faces,:);

    nN = numel(nodes);
    nE = numel(edges);
    nF = numel(faces);

    plotGrid(G, K, 'facealpha', .2, 'edgecolor', 'none')
    hold on;
 
    for i = 1:nE
        edgeNodes = G.edges.nodes(G.edges.nodePos(edges(i)):G.edges.nodePos(edges(i)+1)-1);
        Xe = G.nodes.coords(edgeNodes,:);
        plot3(Xe(:,1), Xe(:,2), Xe(:,3),'k');
    end
    
    remEdges = allEdges(~ismember(allEdges,edges));
    nRE = numel(remEdges);
    for i = 1:nRE
        edgeNodes = G.edges.nodes(G.edges.nodePos(remEdges(i)):G.edges.nodePos(remEdges(i)+1)-1);
        Xe = G.nodes.coords(edgeNodes,:);
        h = plot3(Xe(:,1), Xe(:,2), Xe(:,3));
        set(h, 'color', [0.5 0.5 0.5])
    end
    
    
    addFac = .1;
    for i = 1:nN
        h1 = plot3(X(:,1), X(:,2), X(:,3), 'ok', 'MarkerFaceColor', 'r', 'markersize', 4);
    end
    for i = 1:nE
        h2 = plot3(Ec(:,1), Ec(:,2), Ec(:,3), 'sk', 'MarkerFaceColor', 'b', 'markersize', 4);
    end
    for i = 1:nF
        h3 = plot3(Fc(:,1), Fc(:,2), Fc(:,3), 'dk', 'MarkerFaceColor', 'w', 'markersize', 4);
    end
    h4 = plot3(Kc(:,1), Kc(:,2), Kc(:,3), 'pk', 'MarkerFaceColor', 'm', 'markersize', 4);
    
    
    h = legend([h1, h2, h3, h4], '$\mathcal{V}^K$', '$\mathcal{E}^K$',...
                                 '$\mathcal{F}^K$', '$\mathcal{P}^K$');
    set(h, 'interpreter', 'latex')
    axis equal off;


    %%

    cut = 4;
    ps = get(gcf, 'Position');
    ratio = (ps(4)-ps(2)) / (ps(3)-ps(1));
    paperWidth = 10;
    paperHeight = paperWidth*ratio - cut;
    set(gcf, 'paperunits', 'centimeters');
    set(gcf, 'papersize', [paperWidth paperHeight]);
    set(gcf, 'PaperPosition', [0    0   paperWidth paperHeight]);
    print(gcf, '-dpdf', '../../tex/thesis/fig/dofExample3D.pdf');
end