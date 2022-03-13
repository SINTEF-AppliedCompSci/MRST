function [CG, interaction, triangulations] = storeInteractionRegion(CG, varargin)
%Store interaction region for coarse grid
%
% SYNOPSIS:
%   CG = storeInteractionRegion(CG);
%
% DESCRIPTION:
%   Compute and store the interaction/support regions for coarse blocks.
%   Each coares block will be assigned a set of fine cells for which a
%   MsRSB basis function will be supported. The underlying functions use
%   delaunay triangulation to get the support for unstructured models.
%
% REQUIRED PARAMETERS:
%   CG     - Desired coarse grid.
%
% OPTIONAL PARAMETERS:
%    localTriangulation - Use a local triangulation for each block. Usually
%                         gives the best results.
% RETURNS:
%  CG      - Coarsegrid with additional field CG.cells.interaction
%
% SEE ALSO:
%   `storeIteractionRegionCart`

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
                 'skipSingleCellBlocks', false, ...
                 'simpleCellGrouping', false, ...
                 'ensureConnected',    false, ...
                 'centerOverride',    [], ...
                 'largeBasis',        false, ...
                 'useFaces',         true, ...
                 'edgeBoundaryCenters', true, ...
                 'localTriangulation',   true, ...
                 'useMultipoint', isfield(CG.parent.faces, 'nodePos'));
    opt = merge_options(opt, varargin{:});
    % Store interaction regions in coarse grid    
    G = CG.parent;
    
    
    if ~isfield(CG.faces, 'centers') || ~isfield(CG.cells, 'centers')
        CG = addCoarseCenterPoints(CG, ...
                                   'adjustCenters', opt.adjustCenters, ...
                                   'centerOverride', opt.centerOverride, ...
                                   'edgeBoundaryCenters', opt.edgeBoundaryCenters);
    end
    % Extend face points outside of domain so that we have edge that reach
    % outside the domain
    facePts = G.faces.centroids(CG.faces.centers, :);
    facePts = extrudeFaceCentroids(CG, facePts);
    
    centers = CG.cells.centers;

    interaction = cell(CG.cells.num, 1);
    [assigned, isCenter] = deal(false(G.cells.num, 1));
    isCenter(centers) = true;
    if opt.useMultipoint
        faceNodePos = rldecode(1 : G.faces.num, diff(G.faces.nodePos), 2) .';
    else
        faceNodePos = [];
    end
    
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
            interaction{i} = c;
            assigned(c) = true;
            continue;
        end
        [coarseCells, coarseFaces] = coarseNeighbors(CG, i, opt.useMultipoint, faceNodePos);
        
        coarseCells = [coarseCells; i]; %#ok
        
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
                if ~isempty(tri.ConnectivityList)
                    hull = tri.convexHull();
                    hull = unique(hull(:));
                    tri = delaunayTriangulation(tri.Points(hull, :));
                end
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
        inside = false(G.cells.num, 1);
        if opt.simpleCellGrouping
            % Simpler definition that is sometimes faster
            % inside = ~isnan(tri.pointLocation(G.cells.centroids));
            pts = G.cells.centroids(localCells, :);
            % Add a bit of noise to avoid triangulation having decision
            % problems for square grids
            pts = bsxfun(@plus, pts, sqrt(eps));
            int = ~isnan(evaluateInternal(pts));
            inside(localCells(int)) = true;
        else
            [localFaces, fmap] = gridCellFaces(G, localCells);

            faceIsInside = evaluateInternal(G.faces.centroids(localFaces, :));
            localCellNo = rldecode((1:numel(localCells))', diff(fmap));

            % ok (i) = true indicates that localCells(i) has all faces
            % on the inside of triangulation
            ok = ~isnan(accumarray(localCellNo, faceIsInside));
            inside(localCells(ok)) = true;
        end
        
        % Include all cells inside block as "inside"
        inside(CG.partition == i) = true;
        [bf, bc] = boundaryFaces(G, find(CG.partition == i));
        inside(bc) = true;
        if opt.ensureConnected
            p = processPartition(G, double(inside) + 1);
            inside = p == p(centers(i));
        end
        inside(CG.cells.centers) = false;
        inside(CG.cells.centers(i)) = true;
        cells = find(inside);
        % Ensure that any cells that belong to coarse blocks further away
        % than our neighborship definition are not included.
        ok = false(CG.cells.num, 1);
        ok(coarseCells) = true;
        ok(i) = true;
        
        isLocal = ok(CG.partition(cells));
        % Avoid the centers of other blocks
        isCenter(centers(i)) = false;
        isLocal = isLocal & ~isCenter(cells);
        isCenter(centers(i)) = true;
        
        interaction{i} = cells(isLocal);
        if opt.largeBasis
            interaction{i} = find(ok(CG.partition));
        end
        assigned(interaction{i}) = true;
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

function faceCentroids = extrudeFaceCentroids(CG, faceCentroids)
    if nargin == 1
        faceCentroids = CG.faces.centroids;
    end

    bf = any(CG.faces.neighbors == 0, 2);

    faceNo = sum(CG.faces.neighbors(bf, :), 2);

    fc = faceCentroids(bf, :);
    dist = CG.faces.centroids(bf, :) - CG.cells.centroids(faceNo, :);
    faceCentroids(bf, :) = fc + 2*dist;

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
