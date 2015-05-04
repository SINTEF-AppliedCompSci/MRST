function CG = addCoarseCenterPoints(CG, varargin)
    opt = struct('adjustCenters', true, ...
                 'centerOverride',    [], ...
                 'adjustDims',        1:CG.parent.griddim, ...
                 'edgeBoundaryCenters', true);
    opt = merge_options(opt, varargin{:});
    
    G = CG.parent;
    if opt.adjustCenters
        blockPts = zeros(CG.cells.num, G.griddim);
        for i = 1:CG.cells.num
            fa = gridCellFaces(CG, i);
            blockPts(i, :) = geometricMedian(CG.faces.centroids(fa, :));
        end
    else
        blockPts = CG.cells.centroids;
    end
    
    if opt.edgeBoundaryCenters
        bf = boundaryFaces(CG);
        for i = 1:numel(bf)
            f = bf(i);
            c = sum(CG.faces.neighbors(f, :), 2);
            
            dims = opt.adjustDims;
            cc = blockPts(c, dims);
            fc = CG.faces.centroids(f, dims);
            N = CG.faces.normals(f, dims);
            d = norm(cc - fc, 2);
            N = N./norm(N, 2);
            blockPts(c, dims) = cc - N*d;
        end
    end
    
    if ~isempty(opt.centerOverride)
        assert(size(opt.centerOverride, 1) == CG.cells.num);
        assert(size(opt.centerOverride, 2) == G.griddim);
        ok = ~any(isnan(opt.centerOverride), 2);
        
        blockPts(ok, :) = opt.centerOverride(ok, :);
    end
    
    [CG.cells.centers, CG.faces.centers] = mapCenters(CG, blockPts, CG.faces.centroids);
end




function pt = geometricMedian(pts)
    pt = mean(pts);
    len = @(v) sqrt(sum(v.^2, 2));
    for i = 1:10
        d = len(bsxfun(@minus, pts, pt));
        pt = sum(bsxfun(@rdivide, pts, d), 1)./sum(1./d, 1);
    end
end
