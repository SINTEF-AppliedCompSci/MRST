function w  = createGridMappings(g)

%{ 
Copyright 2009-2014 SINTEF ICT, Applied Mathematics
%} 
%% Create mapping from sub-half-face to cell, node, face, half-face and
%% sub-face
    
    cellno   = rldecode(1:g.cells.num, diff(g.cells.facePos), 2) .';
    col      = 1 + (cellno == g.faces.neighbors(g.cells.faces(:,1), 2));
    nhfaces  = g.cells.facePos(end)-1;
    hfaces   = accumarray([g.cells.faces(:,1), col], 1:nhfaces);
    hfaces   = rldecode(hfaces, diff(g.faces.nodePos));
    cells    = rldecode(g.faces.neighbors, diff(g.faces.nodePos));
    nodes    = repmat(g.faces.nodes, [2,1]);
    faces    = repmat(rldecode(1:g.faces.num, diff(g.faces.nodePos),2)', [2,1]);
    i        = cells~=0;
    w        = [cells(i), nodes(i), hfaces(i), faces(i)];
    w        = double(sortrows(w));
    
end
