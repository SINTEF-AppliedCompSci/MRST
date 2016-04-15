function G = computeVEM3DGeometry(G)

    fprintf('Solving Poisson equation on grid with %d cells \n\n', G.cells.num);
    fprintf('Computing VEM geometry ...\n');
    
    tic;

    G = computeGeometry(G);
    G = mrstGridWithFullMappings(G);
    
    fprintf('... computing edge data\n');

    nodeNum = mcolon(G.edges.nodePos(1:end-1),G.edges.nodePos(2:end)-1);
    nodes = G.edges.nodes(nodeNum);
    edgeVec   = G.nodes.coords(nodes(2:2:end),:) -  ...
                 G.nodes.coords(nodes(1:2:end-1),:);
    lengths   = sqrt(sum(edgeVec.^2,2));
    centroids = (G.nodes.coords(nodes(2:2:end),:) +  ...
                 G.nodes.coords(nodes(1:2:end-1),:))./2;
    clear nodeNum nodes edgeVec


    faceNormals = G.faces.normals;
    edgeNum = mcolon(G.faces.edgePos(1:end-1),G.faces.edgePos(2:end)-1);
    edges = G.faces.edges(edgeNum);
    signs = G.faces.edgeSign(edgeNum);
    nodeNum = mcolon(G.edges.nodePos(edges),G.edges.nodePos(edges+1)-1);
    nodes = G.edges.nodes(nodeNum);
    edgeVec = G.nodes.coords(nodes(2:2:end),:)-G.nodes.coords(nodes(1:2:end-1),:);
    edgeVec = edgeVec.*repmat(signs,1,3);
    normals = cross(edgeVec, rldecode(faceNormals, diff(G.faces.edgePos), 1));
    normals = normals./repmat(sqrt(sum(normals.^2,2)),1,3);

    clear faceNormals edgeNum edges signs nodeNum nodes edgeVec

    nodes = G.cells.nodes;
    edgeNodes = repmat(reshape(G.edges.nodes,2,[])',G.cells.num,1);
    edgeNodes = mat2cell(edgeNodes, G.edges.num*ones(1,G.cells.num),2);
    edges = repmat((1:G.edges.num)',G.cells.num,1);
    edges = mat2cell(edges, G.edges.num*ones(1,G.cells.num),1);
    nodes = mat2cell(nodes,diff(G.cells.nodePos),1);
    cellEdges = cellfun(@(X,Y,Z) Z(sum(ismember(X,Y),2) == 2), ...
           edgeNodes, nodes, edges,  'UniformOutput', false);
    edgePos = [1;cumsum(cellfun(@(X) size(X,1), cellEdges))+1];
    cellEdges = cell2mat(cellEdges);

    clear nodes edgeNodes edges

    G.edges.('lengths')     = lengths;
    G.edges.('centroids')   = centroids;
    G.faces.('edgeNormals') = normals;
    G.cells.('edges')       = cellEdges;
    G.cells.('edgePos')     = edgePos;

    fprintf('... computing diameters\n');
    
    cellDiameters = zeros(G.cells.num,1);
    for i = 1:G.cells.num
        nodeNum = G.cells.nodePos(i):G.cells.nodePos(i+1)-1;
        nodes = G.cells.nodes(nodeNum);
        X = G.nodes.coords(nodes,:);
        cellDiameters(i) = cellDiameter(X);
    end
    
    G.cells.('diameters')   = cellDiameters;
    
    faceDiameters = zeros(G.faces.num,1);
    for i = 1:G.faces.num
        nodeNum = G.faces.nodePos(i):G.faces.nodePos(i+1)-1;
        nodes = G.faces.nodes(nodeNum);
        X = G.nodes.coords(nodes,:);
        faceDiameters(i) = cellDiameter(X);
    end

    G.faces.('diameters')   = faceDiameters;
    
    stop = toc;
    
    fprintf('Done in %f seconds.\n\n', stop);
    
end

