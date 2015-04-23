function [cellcenters, facecenters] = mapCenters(CG, blockCenters, faceCenters)
    if nargin < 3
        faceCenters = CG.faces.centroids;
        if nargin < 2
            blockCenters = CG.cells.centroids;
        end
    end
    
    
    G = CG.parent;
    cellcenters = zeros(CG.cells.num, 1);
    for i = 1:CG.cells.num
        local = find(CG.partition == i);
        cellcenters(i) = doMap(i, local, blockCenters, G.cells.centroids);
    end
    
    if nargout > 1
        facecenters = zeros(CG.faces.num, 1);
        for i = 1:CG.faces.num
            local = CG.faces.fconn(CG.faces.connPos(i):CG.faces.connPos(i+1)-1);
            facecenters(i) = doMap(i, local, faceCenters, G.faces.centroids);
        end
    end
end


function l = doMap(coarseind, localfine, coarsepts, finepts)
    % Add a tiny bit of bias to the selection so that we avoid
    % consistent mappings on cartesian grids where the center is exactly
    % between two cells
    dist = bsxfun(@(x, y) x-y, finepts(localfine,:), ...
                               coarsepts(coarseind,:) - 10*sqrt(eps));
    [v, ind] = min(sqrt(sum(dist.^2, 2))); %#ok
    l = localfine(ind);
end
