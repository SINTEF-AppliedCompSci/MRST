function plotCubature(G, cubature, elements, varargin)

    opt = struct('plotBoundingBox', true, 'plotTri', true, 'smin', 5, 'smax', 20, 'markerStyle', 'ok');
    [opt, extra] = merge_options(opt, varargin{:});
    
    if G.griddim - cubature.dim == 1
        faces = elements;
        cells = reshape(G.faces.neighbors(faces,:), [], 1);
        [~, x, w, cellNo, faceNo] = cubature.getCubature(faces, 'face');
    else
        cells = elements;
        [~, x, w, cellNo] = cubature.getCubature(cells, 'cell');
    end

    hold on
    plotGrid(G, cells, extra{:}); axis equal off
    if opt.plotBoundingBox
        for i = 1:numel(cells)
            rectangle('Position', [G.cells.centroids(cells(i),:) - G.cells.dx(cells(i),:)/2, G.cells.dx(cells(i),:);]);
        end
    end
    
    if any(strcmpi(class(cubature), {'TriangleCubature', 'TetrahedronCubature'})) && opt.plotTri
        tri = cubature.triangulation.triPos(cellNo):cubature.triangulation.triPos(cellNo+1)-1;
        if G.griddim == 2
            trimesh(cubature.triangulation.ConnectivityList(tri,:), ...
                    cubature.triangulation.Points(:,1), ...
                    cubature.triangulation.Points(:,2), ...
                    'linestyle', '-', 'color', 'k')
        else
            trimesh(cubature.triangulation.ConnectivityList(tri,:), ...
                    cubature.triangulation.Points(:,1), ...
                    cubature.triangulation.Points(:,2), ...
                    cubature.triangulation.Points(:,3), ...
                    'linestyle', '-', 'color', 'k')
        end
    end

    alpha = (w - min(w))/(max(w) - min(w));
    alpha(isnan(alpha)) = 1;
    if isempty(opt.smax)
        opt.smax = opt.smin*max(w)/min(w);
    end
    scaling = opt.smin*(1-alpha) + opt.smax*alpha;
    for i = 1:numel(w)
        plot(x(i,1), x(i,2), opt.markerStyle, 'markerSize', scaling(i), 'markerFaceColor', 'k')
    end

end