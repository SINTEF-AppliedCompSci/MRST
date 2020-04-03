function [G, bcfaces] = rotateGrid(G, angle)
    
    
    dim = G.griddim;
    assert(dim == 2, 'only 2D implementation for now');
    nodes = G.nodes.coords;
    
    bcfaces = getbcfaces(G);
    
    nodetbl.nodes = (1 : G.nodes.num)';
    nodetbl = IndexArray(nodetbl);
    
    vectbl.vecdim = (1 : dim)';
    vectbl = IndexArray(vectbl);
    
    % Index array for vector at nodes
    nodevectbl = crossIndexArray(nodetbl, vectbl, {});
    
    nodecents = G.nodes.coords;
    nodecents = reshape(nodecents', [], 1);
    
    
    vec2tbl = crossIndexArray(vectbl, vectbl, {}, 'crossextend', {{'vecdim', ...
                        {'vecdim1', 'vecdim2'}}});
    rotmat = [cos(angle); -sin(angle); sin(angle); cos(angle)];
    

    prod = TensorProd();
    prod.tbl1 = vec2tbl;
    prod.tbl2 = nodevectbl;
    prod.tbl3 = nodevectbl;
    prod.replacefds1 = {{'vecdim1', 'vecdim'}};
    prod.replacefds2 = {{'vecdim', 'vecdim2'}};
    prod.reducefds = {'vecdim2'};
    
    prod = prod.setup();
    
    nodecents = prod.eval(rotmat, nodecents);
    
    nodecents = reshape(nodecents, 2, [])';
    
    G.nodes.coords = nodecents;
    
end

function bcfaces = getbcfaces(G)
    
    G = computeGeometry(G);
    
    xvals = G.faces.centroids(:, 1);
    yvals = G.faces.centroids(:, 2);
    
    xmin = min(xvals);
    bcfaces.xmin = find(xvals - xmin < eps);
    xmax = max(xvals);
    bcfaces.xmax = find(xmax - xvals < eps);    
    
    ymin = min(yvals);
    bcfaces.ymin = find(yvals - ymin < eps);
    ymax = max(yvals);
    bcfaces.ymax = find(ymax - yvals < eps);    
    
end
