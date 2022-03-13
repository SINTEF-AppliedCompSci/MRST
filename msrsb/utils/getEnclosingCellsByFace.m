function [cells, isOutside] = getEnclosingCellsByFace(G, pts, varargin)
% Get cells enclosing a set of points using dot product with center-face
% vectors.

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
    opt = struct('candidates', []);
    opt = merge_options(opt, varargin{:});

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
        if isempty(opt.candidates)
            candidates = cellNo(unique(conn(:)));
        else
            candidates = opt.candidates;
        end
        
        for j = 1:numel(candidates)
            c = candidates(j);
            faces = gridCellFaces(G, c);
            sgn = 1 - 2*(G.faces.neighbors(faces, 2) == c);
            fc = G.faces.centroids(faces, :);
            % Check sign of dot product
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
