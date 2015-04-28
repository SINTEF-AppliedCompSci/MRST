function [CG, interaction, triangulations] = storeInteractionRegion(CG, varargin)
%Store interaction region for coarse grid
%
% SYNOPSIS:
%   CG = storeInteractionRegion(CG);
%
% DESCRIPTION:
%   Compute and store the interaction regions for coarse blocks. Each
%   coares block will be assigned a set of fine cells for which a MsRSB
%   basis function will be supported. The underlying functions use
%   delaunay triangulation to get the support for unstructured models.
%
% REQUIRED PARAMETERS:
%   CG     - Desired coarse grid.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   TODO
%
% RETURNS:
%  CG
%
% EXAMPLE:
%    
%
% SEE ALSO:
%   storeIteractionRegionCart

%{
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

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
                 'skipSingleCellBlocks', false, ...
                 'simpleCellGrouping', false, ...
                 'ensureConnected',    false, ...
                 'centerOverride',    [], ...
                 'largeBasis',        false, ...
                 'useFaces',         true, ...
                 'edgeBoundaryCenters', true, ...
                 'localTriangulation',   true, ...
                 'useMultipoint', true);
    opt = merge_options(opt, varargin{:});
    % Store interaction regions in coarse grid    
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
            
            cc = blockPts(c, :);
            fc = CG.faces.centroids(f, :);
            N = CG.faces.normals(f, :);
            d = norm(cc - fc, 2);
            N = N./norm(N, 2);
            blockPts(c, :) = cc - N*d;
        end
    end
    
    if ~isempty(opt.centerOverride)
        assert(size(opt.centerOverride, 1) == CG.cells.num);
        assert(size(opt.centerOverride, 2) == G.griddim);
        ok = ~any(isnan(opt.centerOverride), 2);
        
        blockPts(ok, :) = opt.centerOverride(ok, :);
    end
    
    [CG.cells.centers, CG.faces.centers] = mapCenters(CG, blockPts, CG.faces.centroids);
    % Extend face points outside of domain so that we have edge that reach
    % outside the domain
    facePts = extrudeFaceCentroids(CG, G.faces.centroids(CG.faces.centers, :), ...
                                       G.cells.centroids(CG.cells.centers, :));
    
    centers = CG.cells.centers;

    interaction = cell(CG.cells.num, 1);
    assigned = false(G.cells.num, 1);
    
    faceCell = rldecode(1 : G.cells.num, diff(G.cells.facePos), 2) .';
    faceNodePos = rldecode(1 : G.faces.num, diff(G.faces.nodePos), 2) .';
    
    if ~opt.localTriangulation
        if opt.useFaces
        	points = vertcat(facePts, G.cells.centroids(centers, :));
        else
            points = G.cells.centroids(centers, :);
        end
        tri = delaunayTriangulation(points);
        triangulations = tri;
    else
        triangulations = cell(CG.cells.num, 1);
    end
    for i = 1:CG.cells.num
        fprintf('Handled coarse block %d / %d\n', i, CG.cells.num);
        
        if sum(CG.partition == i) == 1 && opt.skipSingleCellBlocks
            c = find(CG.partition == i);
            if 1
                n = c;
            else
                n = G.faces.neighbors(any(G.faces.neighbors == c, 2), :);
                n = unique(n(:));
                n = n(n~=0);
            end
            interaction{i} = n;
            assigned(n) = true;
            continue;
        end
        [coarseCells, coarseFaces] = coarseNeighbors(CG, i, opt.useMultipoint, faceNodePos);
        
        coarseCells = [coarseCells; i]; %#ok not really the case
        
        if opt.localTriangulation
            if opt.useFaces
                points = vertcat(facePts(coarseFaces, :), ...
                                 G.cells.centroids(centers(coarseCells), :));
            else
                points = G.cells.centroids(centers(coarseCells), :);
            end

            tri = delaunayTriangulation(points);
            if opt.useFaces
                % Use convex hull
                hull = tri.convexHull();
                hull = unique(hull(:));
                tri = delaunayTriangulation(tri.Points(hull, :));
            end
            
            if nargout > 2
                triangulations{i} = tri;
            end
            if ~isempty(tri.ConnectivityList)
                evaluateInternal = @(pts) tri.pointLocation(pts);
            else
                evaluateInternal = @(pts) nan(size(pts, 1), 1);
            end
            
        else
            indexInTriArray = [coarseFaces; ...
                              (coarseCells) + CG.faces.num];
            evaluateInternal = @(pts) globalEvaluateTriInside(pts, tri, indexInTriArray);
        end
        
        isCoarseNeigh = false(CG.cells.num, 1);
        isCoarseNeigh(coarseCells) = true;

        localCells = find(isCoarseNeigh(CG.partition));
        if opt.simpleCellGrouping
            % Simpler definition that is sometimes faster
            % inside = ~isnan(tri.pointLocation(G.cells.centroids));
            pts = G.cells.centroids(localCells, :);
            % Add a bit of noise to avoid triangulation having decision
            % problems for square grids
            pts = bsxfun(@plus, pts, sqrt(eps));
            int = ~isnan(evaluateInternal(pts));
            
            inside = false(G.cells.num, 1);
            inside(localCells(int)) = true;
        else
            [localFaces, fmap] = gridCellFaces(G, localCells);

            faceIsInside = evaluateInternal(G.faces.centroids(localFaces, :));
            localCellNo = rldecode((1:numel(localCells))', diff(fmap));

            % ok (i) = true indicates that localCells(i) has all faces
            % on the inside of triangulation
            ok = ~isnan(accumarray(localCellNo, faceIsInside));

            inside = false(G.cells.num, 1);
            inside(localCells(ok)) = true;
        end
        
        % Include all cells inside block as "inside"
        inside(CG.partition == i) = true;
        [bf, bc] = boundaryFaces(G, find(CG.partition == i));
        inside(bc) = true;
        if opt.ensureConnected
            p = processPartition(G, double(inside) + 1);
            inside = p == p(centers(i));
            cells = find(inside);
        else
            cells = find(inside);
        end
        % Ensure that any cells that belong to coarse blocks further away
        % than our neighborship definition are not included.
        ok = false(CG.cells.num, 1);
        ok(coarseCells) = true;
        ok(i) = true;
        
        isLocal = ok(CG.partition(cells));
        
        interaction{i} = cells(isLocal);
        if opt.largeBasis
            interaction{i} = find(ok(CG.partition));
        end
        assigned(interaction{i}) = true;
    end
    
    % Corner case: If a coarse block interaction covers all cells of
    % another coarse block, remove all this support immediately.
    if ~opt.largeBasis && 0
        disp('Glorb!')
        for i = 1:CG.cells.num
            c = interaction{i};
            pc = CG.partition(c);
            p = unique(pc);

            keep = true(numel(c), 1);
            for j = 1:numel(p)
                if p(j) == i
                    continue;
                end
                current = pc == p(j);
                if sum(current) == sum(CG.partition == p(j))
                    keep(current) = false;
                end
            end
            interaction{i} = c(keep);
        end
    end
    
    
    % If cells are not member of any triangulation (somewhere on the
    % boundary, typically) these are assigned to the interaction region of
    % the coarse block they belong to.
    bad = find(~assigned);
    for i = 1:numel(bad)
        bf = bad(i);
        bc = CG.partition(bad(i));
        interaction{bc} = [interaction{bc}; bf];
    end
    
    CG.cells.interaction = interaction;
end

function faceCentroids = extrudeFaceCentroids(CG, faceCentroids, cellCentroids)
    if nargin == 1
        faceCentroids = CG.faces.centroids;
    end
    
    if nargin < 3
        cellCentroids = CG.cells.centroids;
    end
    bf = any(CG.faces.neighbors == 0, 2);

    faceNo = sum(CG.faces.neighbors(bf, :), 2);

    fc = faceCentroids(bf, :);
    cc = cellCentroids(faceNo, :);
    faceCentroids(bf, :) = fc + 2*(fc - cc);

end


function pt = geometricMedian(pts)
    pt = mean(pts);
    len = @(v) sqrt(sum(v.^2, 2));
    for i = 1:10
        d = len(bsxfun(@minus, pts, pt));
        pt = sum(bsxfun(@rdivide, pts, d), 1)./sum(1./d, 1);
    end
end

function ok = globalEvaluateTriInside(pts, tri, indexInTriArray)
    tmp = [nan, 1];
    triangles = tri.pointLocation(pts);
    bad = isnan(triangles);
    triangles(bad) = 1;
    internal = all(ismember(tri.ConnectivityList(triangles, :), indexInTriArray), 2);
    internal(bad) = false;
    ok = tmp(internal + 1);
end
