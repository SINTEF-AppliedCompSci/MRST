function G = addBoundingBoxFields(G, varargin)
% Add minimal bounding box for each cell/face
% Useful for quickly picking subset of candidate cells/faces before more 
% costly geometry is invoked
%
% SEE ALSO:
%  computeVerticalGridIntersection, computeTraversedCells,
%  addBoundingBoxFields
opt = struct('cells', false, ...
             'faces', true);
opt = merge_options(opt, varargin{:});

% cells (cellNodes fails for e.g., Olympus)
if opt.cells
    cno  = cellNodes(G);
    nn   = accumarray(cno(:,1), ones(size(cno,1), 1));
    npos = cumsum([1;nn]);
    % max and min ix after sorting
    minIx = npos(1:end-1);
    maxIx = npos(2:end)-1;
    bbox = nan(G.cells.num, G.griddim);
    for d = 1:G.griddim
        tmp = sortrows([cno(:,1), G.nodes.coords(cno(:,3),d)]);
        bbox(:,d) = tmp(maxIx, 2) - tmp(minIx, 2);
    end
    G.cells.bbox = bbox;
end

% faces
if opt.faces
    npos = G.faces.nodePos;
    fno  = rldecode((1:G.faces.num)', diff(npos));
    % max and min ix after sorting
    minIx = npos(1:end-1);
    maxIx = npos(2:end)-1;
    bbox = nan(G.faces.num, G.griddim);
    for d = 1:G.griddim
        tmp = sortrows([fno(:,1), G.nodes.coords(G.faces.nodes, d)]);
        bbox(:,d) = tmp(maxIx, 2) - tmp(minIx, 2);
    end
    G.faces.bbox = bbox;
end
end
