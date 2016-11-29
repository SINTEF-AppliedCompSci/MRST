function faces = adjacentFacesForFace(G, f, cts)
    if nargin < 3
        cts = 2;
    end
    active = false(G.nodes.num, 1);
    active(gridFaceNodes(G, f)) = true;

    % Count number of nodes shared per face
    nPos = rldecode(1:G.faces.num, diff(G.faces.nodePos), 2) .';
    counts = accumarray(nPos, active(G.faces.nodes(:, 1)));
    % Take as active anyone with more than cts overlapping nodes
    faces = find(counts >= cts);
    faces = faces(faces~=f);
end