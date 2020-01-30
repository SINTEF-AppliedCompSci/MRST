function intersections = surfaceIntersections3D(surfaces)
% Find the intersection of surfaces defined by polygons. 
%
% SYNOPSIS:
%   intersections = surfaceIntersections3D(surfaces)
%
% PARAMETERS
%   surfaces        - cell array of surfaces. Each element is a vector of
%                     the vertices of the surface
% RETURNS:
%   intersections   - cell-array of the intersections of the surfaces. Each
%                     element is a cell array containing the start and end
%                     point of the intersection line and the indices of the
%                     two surfaces that generated this intersection
%
% EXAMPLE:
%   f1 = [1,3,2; 4,3,2; 4,3,4; 1,3, 4];
%   f2 = [2,2,3.3; 5,2,3.3; 5,4,3.3; 2,4, 3.3];
%   surfaces = {f1, f2};
%   intersections = surfaceIntersections3D(surfaces);
%   % The intersection line is:
%   intersections{1}{1}
%   % The indices of the surfaces are
%   intersections{1}{2:3}
%
% SEE ALSO
%   lineGrid3D, surfaceGrid3D, volumeGrid3D, compositePebiGrid2D, pebi, surfaceSites2D, lineSites2D.

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}  

intersections = {};
for f1 = 1:numel(surfaces)
    % first find intersection between fault i and all other faults
    for f2 = f1+1:numel(surfaces)
       int_pts = polygonIntersection(surfaces{f1}, surfaces{f2});
       if size(int_pts, 1) > 0
           intersections= {intersections{:}, {int_pts, f1, f2}};
       end
    end
end

% then split intersection of intersections

for f = 1:numel(surfaces)
    f_pts = surfaces{f};
    center = mean(f_pts, 1);
    f_center = bsxfun(@minus, f_pts, center);
    R = rotationMatrixFromPlane(f_center);
    f_xy = f_center * R';
    assert(sum(abs(f_xy(:, 3)))<1e-6);
    f_xy = f_xy(:,1:2);
    internal_pts_xy = zeros(0, 2);
    int_in_plane = false(numel(intersections),1);
    for j = 1:numel(intersections)
        if intersections{j}{2}==f || intersections{j}{3}==f
%              if intersections{j}{2} < f || intersections{j}{3} < f
%                  % did already calculate this intersection
%                  continue
%              end
            int_in_plane(j) = true;
            pts_i = intersections{j}{1};
            pts_i_center = bsxfun(@minus, pts_i, center);
            pts_i_xy = pts_i_center * R';
            assert(all(abs(pts_i_xy(:,3))<1e-6))
            pts_i_xy = pts_i_xy(:, 1:2);
            assert(size(pts_i_xy,1)==2)
            internal_pts_xy = [internal_pts_xy; pts_i_xy];
        end
    end
    num_lines = size(internal_pts_xy,1)/2;
    internal_pts_xy = mat2cell(internal_pts_xy, 2 * ones(num_lines,1), 2)';
    [split_pts_xy, ~, ~, ic] = splitAtInt2D(internal_pts_xy, {});
    %split_pts_xy = cell2mat(split_pts_xy)';
    split_pts = cellfun(@(c) [c, zeros(2, 1)], split_pts_xy, 'un',0);
    split_pts = cellfun(@(c) c * R + center, split_pts, 'un', 0);
    
    % remove old intersection lines and add new
    int_vals_in_plane = intersections(int_in_plane);
    new_intersections = int_vals_in_plane(ic);
    if numel(new_intersections) == 0
        continue
    end
    new_intersections = cellfun(@(c1,c2) [c1, c2(2:3)], split_pts, new_intersections, 'un', 0);
    intersections = [intersections(~int_in_plane), new_intersections];


end
end
