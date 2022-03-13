function [dist, f, isInside] = getDistanceToBoundary(G, p, vec, faceIx)
% Compute distance from a given point in given direction to external boundary
%
% SYNOPSIS:
%   [d, fix, pnt] = getDistanceToBoundary(G, p, v, faceIx)
%
% REQUIRED PARAMETERS:
%    G      - Grid structure
%    p      - 1x3 coodinate of point
%    v      - 1x3 direction
%    faceIx - (optional) boundary faces of region to consider
%             default: faceIx = boundaryFaces(G);
%
% RETURNS:
%    d        - distance to external boundary from p in direction v
%               For multiple intersections, the maximal distance is returned.
%               If point lies outside external boundary and v points away from
%               the region, a search in the opposite direction is performed, and 
%               returs a negative d corresponing to smallest in magnitide distance.
%    fix      - index of corresponding intersected face
%    isInside - if requiested, tests if point p is inside grid. NOTE: will
%               return false if p is inside hole in domain 

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

if norm(vec) < sqrt(eps)*max(norm(p), 1)
    [dist, f, isInside] = deal(inf, [], []);
    return;
end
     
if nargin < 4
    faceIx = boundaryFaces(G);
end

if ~isfield(G.faces, 'bbox')
    warning('Precompute bounding boxes for grid faces for improved performance');
    G = addBoundingBoxFields(G);
end
if ~isfield(G, 'bbox')
    G.bbox = [min(G.nodes.coords); max(G.nodes.coords)];
end

% get sufficent length (make sure we get all intersections)
maxLength = norm(diff(G.bbox)) + norm(p-mean(G.bbox));

vec = vec/norm(vec);
[dist, f, isInside] = deal([]);
dist = inf;
s_max = [];
for direction = [1, -1]
    seg = [p; p+direction*maxLength*vec];
    
    %hold on
    %if direction == 1, col = '-ob'; else col = '-or';end
    %plot3(seg(:,1), seg(:,2), seg(:,3), col,'LineWidth', 2);
    
    fix = getFacesCloseToSegment(G, seg, 'faceIx', faceIx);
    if ~isempty(fix)
        % get face triangulation
        [v1, v2, v3, faceNo] = getTriangles(G, fix);
        
        % edges
        e1 = v2-v1;
        e2 = v3-v1;
        
        % some geometry (see moller-trumbore)
        t = seg(1,:)-v1;
        d = repmat(diff(seg), size(t, 1), 1);
        
        sg = cross(d, e2, 2);
        q = cross(t, e1, 2);
        den = dot(sg ,e1, 2);
        
        s = dot(q, e2, 2)./den;
        u = dot(sg, t,  2)./den;
        v = dot(q, d,  2)./den;
        
        % tolerance should be >0 to avoid missing intersections
        tol = sqrt(eps);
        % need 0 <= s <= 1 to be along segment
        chk = s >= -tol & s <= 1+tol;
        
        % intersection-test in barycentric coord
        chk = chk & u >= -tol & v >= -tol & (u+v) <= 1+tol;
        
        [s, faceNo, t_ix] = deal(s(chk), faceNo(chk), find(chk));
        if any(chk)
            [s_max, ix_max] = max(direction*s);
            if direction == 1
                s_max = max(s_max, 0);
            else
                s_max = min(s_max, 0);
            end
            if nargout == 3
                [s_min, ix_min] = min(abs(s));
                % check sign of this intersection
                isInside  = true;
                if abs(s_min) > tol
                    % see computeTraversedCells
                    nt  = cross(e1(t_ix(ix_min),:),e2(t_ix(ix_min),:), 2);
                    n   = G.faces.normals(faceNo(ix_min),:);
                    nt  = sign(dot(n,nt,2))*nt;
                    sgn = sign(nt*(diff(seg)'));
                    nix = 2;
                    if sgn< 1
                        nix = 1;
                    end
                    isInside = G.faces.neighbors(faceNo(ix_min), nix) == 0;
                end
            end
            break;
        end
    end
end
if ~isempty(s_max)
    dist   = maxLength*s_max;
    f = faceNo(ix_max);
else
    warning('No intersections found')
end
end


%--------------------------------------------------------------------------

function [v1, v2, v3, faceNo] = getTriangles(G, f)
isSubset = true;
if isempty(f)
    f = (1:G.faces.num);
    isSubset = false;
end

if ~isSubset % full grid
    faceNo  = rldecode(1:G.faces.num, diff(G.faces.nodePos), 2) .';
    p       = G.faces.nodePos;
    next    = (2:size(G.faces.nodes, 1)+1) .';
    next(p(2 : end) - 1) = p(1 : end-1);
    
    v1 = G.nodes.coords(G.faces.nodes,:);
    v2 = G.nodes.coords(G.faces.nodes(next),:);
    v3 = G.faces.centroids(faceNo,:);
else
    [np1, np2] = deal(G.faces.nodePos(f), G.faces.nodePos(f+1));
    faceNo = rldecode(f(:), np2-np1);
    
    nodes = G.faces.nodes(mcolon(np1, np2-1));
    nnode = np2-np1;
    locpos = [1; cumsum(nnode)+1];
    next   = (2:numel(nodes)+1) .';
    next(locpos(2 : end) - 1) = locpos(1 : end-1);
    
    % triangle vertices
    v1 = G.nodes.coords(nodes,:);
    v2 = G.nodes.coords(nodes(next),:);
    v3 = G.faces.centroids(faceNo,:);
end
end