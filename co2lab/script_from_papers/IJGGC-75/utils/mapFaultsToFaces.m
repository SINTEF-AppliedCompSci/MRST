function [faces, active] = mapFaultsToFaces(G, faults, varargin)
    opt = struct('refineFactor', 1000);
    opt = merge_options(opt, varargin{:});
    
    vb = mrstVerbose();
    
    % Define roughly cell size for refining the line segments
    dc = max(G.cells.centroids) - min(G.cells.centroids);
    cellsize = min(dc(1:2)./G.cartDims(1:2));
    % Triangulate the grid
    dispif(vb, 'Triangulating...\n');
    timer = ticif(vb);
    % Add virtual points at the boundary of the domain
    bf = boundaryFaces(G);    
    N = G.faces.neighbors;
    for i = 1:numel(bf)
        if N(bf(i), 1) == 0
            N(bf(i), 1) = i + G.cells.num;
        else
            N(bf(i), 2) = i + G.cells.num;
        end
    end
    % Pad the neighborship with additional zeros for the new boundary
    padding = [zeros(numel(bf), 1), (1:numel(bf))' + G.cells.num];
    N = [N; padding];
    % Points to triangulate are all cells + boundary faces
    allpts = [G.cells.centroids; G.faces.centroids(bf, :)];
    nf = size(N, 1);
    % Actual triangulation
    T = delaunayTriangulation(allpts);
    tocif(vb, timer);
    
    % Find the triangles that each face belongs to. A face is associated
    % with a triangle if two of the three points in the triangle correspond
    % to the cells on each side of the interface.
    dispif(vb, 'Finding triangle mapping...\n');
    timer = ticif(vb);
    C = T.ConnectivityList;
    ntri = size(C, 1);
    % 10 more than sufficient for regular grids
    triangles = nan(nf, 10);
    for f = 1:nf
        n1 = N(f, 1);
        n2 = N(f, 2);
        if n1 == 0 
            tt = any(C == n2, 2);
        elseif n2 == 0
            tt = any(C == n1, 2);
        else
            tt = any(C == n1, 2) & any(C == n2, 2);
        end
        tt = find(tt);
        triangles(f, 1:numel(tt)) = tt;
    end
    % Strip extra columns added as a safeguard.
    triangles = triangles(:, ~all(isnan(triangles), 1));
    triangles(isnan(triangles)) = 0;
    
    n_el = size(triangles, 2) - sum(triangles == 0, 2);
    tocif(vb, timer);
    
    % Categorize faces. A face is categorized if all triangles associated
    % with that face is intersected by the line segment. This gives us an
    % even division of the faults (stair-stepping and avoiding closed of
    % cells unless the faults are very close).
    timer = ticif(vb);
    dispif(vb, 'Categorizing faces...\n');
    faces = cell(numel(faults), 1);
    active = false(numel(faults), 1);
    for i = 1:numel(faults)
        fpts = refine(faults{i}, cellsize, opt);        
        pl = T.pointLocation(fpts);
        pl(isnan(pl)) = [];
        if isempty(pl)
            continue
        end
        active(i) = true;
        tmp = false(ntri+1, 1);
        tmp(pl+1) = true;
        
        faces{i} = find(sum(tmp(triangles+1), 2) >= n_el);
        faces{i} = faces{i}(faces{i} <= G.faces.num);
    end
    dispif(vb, 'All done!\n');
    tocif(vb, timer);
end

function p = refine(pts, cellsize, opt)
    % Refine faultline according to approx. grid size, length of segment
    % and optional fudge parameter
    l = 0;
    d = zeros(size(pts, 1), 1);
    for i = 2:size(pts, 1)
        lseg = norm(pts(i-1, :) - pts(i, :), 2);
        l = l + lseg;
        d(i) = l;
    end
    d = d./d(end);
    npt = opt.refineFactor*l/cellsize;
    d_new = linspace(0, 1, npt)';
    p = interp1(d, pts,  d_new, 'linear');
end