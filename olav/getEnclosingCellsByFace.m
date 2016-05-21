function [cells, isOutside] = getEnclosingCellsByFace(G, pts, varargin)
    opt = struct('candidates',[]);
    opt = merge_options(opt,varargin{:});
    bf = boundaryFaces(G);
    
    % Take all cells centers + boundary faces

    tripts = [G.cells.centroids; G.faces.centroids(bf, :)];
    if isempty(opt.candidates)
        cellNo = [(1:G.cells.num)'; sum(G.faces.neighbors(bf, :), 2)];
    end
    T = delaunayTriangulation(tripts);

    loc = T.nearestNeighbor(pts);
    
    % Map outside points to the inside
    isOutside = isnan(T.pointLocation(pts));
    pts(isOutside, :) = tripts(loc(isOutside), :);
    
    npts = size(pts, 1);
    
    cells = zeros(npts, 1);
    n = bsxfun(@rdivide, G.faces.normals, G.faces.areas);
    for i = 1:npts
        pt = pts(i, :);
        ok = false;
        tri = any(T.ConnectivityList == loc(i), 2);
        conn = T.ConnectivityList(tri, :);
        candidates = cellNo(unique(conn(:)));
        
        for j = 1:numel(candidates)
            c = candidates(j);
            faces = gridCellFaces(G, c);
            sgn = 1 - 2*(G.faces.neighbors(faces, 2) == c);
            fc = G.faces.centroids(faces, :);
            
            normals = bsxfun(@times, n(faces, :), sgn);
            
            outBound = bsxfun(@minus, fc, pt);
            
            d = sum(normals.*outBound, 2);
            if all(d >= 0)
                ok = true;
                break
            end
        end
        if ~ok
            warning('Unable to find position for point %d', i)
            isOutside(i) = true;
        end
        cells(i) = c;
    end
end
