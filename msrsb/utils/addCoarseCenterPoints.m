function CG = addCoarseCenterPoints(CG, varargin)
% Adds center points to each block. A center is typically the fine cell
% closest to the geometrical centroid of the block.

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
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
            cc = blockPts(c, :);
            fc = CG.faces.centroids(f, :);
            N = CG.faces.normals(f, :);
            d = norm(cc - fc, 2);
            N = N./norm(N, 2);
            new = cc - N*d;
            blockPts(c, dims) = new(1, dims);
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
