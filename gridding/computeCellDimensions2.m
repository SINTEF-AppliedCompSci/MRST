function G = computeCellDimensions2(G)

    % Get extra grid topology
    G = createAugmentedGrid(G);

    % Get node coordinates
    xn  = G.nodes.coords(G.cells.nodes,:);
    ncn = diff(G.cells.nodePos);
    
    % Get minimum and maximum cell coordinates
    [G.cells.xMin, G.cells.xMax] = getMinMax(xn, ncn);

    % Compute cell bounding box dimensions
    G.cells.dx = G.cells.xMax - G.cells.xMin;
    
    % Get face coordinates
    xn  = G.nodes.coords(G.faces.nodes,:);
    nfn = diff(G.faces.nodePos);
    if G.griddim == 3
        % Create local coordinate systems on each face
        G.faces.coordSys = faceCoordSys(G);
%         
%         % Map to face coordinate system
%         node2face = rldecode((1:G.faces.num)', nfn, 1);
%         xnf       = zeros(sum(nfn), 2); 
%         for dNo = 1:2
%             vec = G.faces.coordSys{dNo}(node2face,:);
%             xnf(:, dNo) = sum(xn.*vec,2);
%         end
%         xn = xnf;
%         
    end
        
    % Get minimum and maximum cell coordinates
    [G.faces.xMin, G.faces.xMax] = getMinMax(xn, nfn);
    G.faces.dx = G.faces.xMax - G.faces.xMin;
        
     
end

function coordSys = faceCoordSys(G)
    
    % Get nodes of last edge of each face
    nfe   = diff(G.faces.edgePos);
    edges = G.faces.edges(cumsum(nfe));
    nodes = G.edges.nodes(mcolon(G.edges.nodePos(edges), ...
                                 G.edges.nodePos(edges+1)-1));

    % First coordinate axis points along last edge
    vec1 = G.nodes.coords(nodes(2:2:end),:) ...
         - G.nodes.coords(nodes(1:2:end-1),:);
    vec1 = vec1./sqrt(sum(vec1.^2,2));
    
    % Second vector constructed as cross product of vec1 and face normal
    n    = G.faces.normals./G.faces.areas;
    vec2 = cross(vec1, n, 2);
    
    coordSys = {vec1, vec2};

end

