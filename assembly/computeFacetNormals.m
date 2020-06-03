function facetNormals =  computeFacetNormals(G, cellnodefacetbl)
% compute facet normals (area weighted)
    fno = cellnodefacetbl.get('faces');
    cno = cellnodefacetbl.get('cells');
    
    numnodes = double(diff(G.faces.nodePos));
    numnodes = numnodes(fno);
    
    facetNormals = G.faces.normals(fno, :);
    facetNormals = bsxfun(@ldivide, numnodes, facetNormals);

    sgn = 2*(cno == G.faces.neighbors(fno, 1)) - 1;
    facetNormals = sgn.*facetNormals; % Outward normals with respect to cell
                                      % in cellnodefacetbl.
    
    facetNormals = reshape(facetNormals', [], 1);
    
end

