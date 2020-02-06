function [grids_2] = surfaceGrid3D(surfaces, grids_1, intersections, ds, gamma)
% Find the sites and grids of surfaces such that the 2D PEBI-grids conform to
% the surface intersections.
%
% SYNOPSIS:
%   grids_2 = surfaceGrid3D(faults, grids_1, intersections, ds, gamma)
%
% PARAMETERS
%   surfaces        - cell array of surface. Each element is a vector of
%                     the vertices of the surface
%   grids_1         - cell array of 1D grids as returned from
%                     lineGrid3D
%   intersecitons   - The intersection of surfaces as returned from
%                     surfaceIntersections3D
%   ds              - Cell size of the 2D cells
%   gamma           - Distance surface sites are placed from the 1D grids
%
% RETURNS:
%   grids_2         - cell-array where each element is a 2D grid. In
%                     addition to a normal MRST-grid the field
%                     G.cells.sites is added, which is an array of the
%                     sites used to generate the grid. Likewise, we add a
%                     Boolean array G.cells.resSite that indicates if the
%                     the site is a background reservoir site or not
%
% EXAMPLE:
%   f1 = [1,3,2; 4,3,2; 4,3,4; 1,3, 4];
%   f2 = [2,2,3.3; 5,2,3.3; 5,4,3.3; 2,4, 3.3];
%   fracs = {f1, f2};
%   intersections = surfaceIntersections3D(fracs);
%   grids_1 = lineGrid3D(intersections, 0.2);
%   grids_2 = fautlSites(fracs, grids_1, intersections, 0.2, 0.2/6)
%   figure(); clf; hold on
%   for i =1:numel(grids_2)
%       plotGrid(grids_2{i})
%   end
% SEE ALSO
%   lineGrid3D, surfaceIntersections3D, volumeGrid3D, compositePebiGrid2D, pebi, surfaceSites2D, lineSites2D.

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}  
%sites_2 = {};
grids_2 = {};
for f = 1:numel(surfaces)
    f_pts = surfaces{f};
    center = mean(f_pts, 1);
    f_center = bsxfun(@minus, f_pts, center);
    R = rotationMatrixFromPlane(f_center);
    f_xy = f_center * R';
    assert(sum(abs(f_xy(:, 3)))<1e-6);
    f_xy = f_xy(:,1:2);
        
    internal_pts = [];
    tip_sites = [];
    grids_1_projected = {};
    for j = 1:numel(intersections)
        if intersections{j}{2}==f || intersections{j}{3}==f 
            pts_1 = grids_1{j}{1};
            vert_1 = grids_1{j}{2};
            if size(pts_1,1) == 0
                continue
            end
            pts_1_center = bsxfun(@minus, pts_1, center);
            pts_1_xy = pts_1_center * R';
            vert_1_center = bsxfun(@minus, vert_1, center);
            vert_1_xy = vert_1_center * R';
            assert(all(abs(pts_1_xy(:,3))<1e-6))
            assert(all(abs(vert_1_xy(:,3))<1e-6))
            
            pts_1_xy = pts_1_xy(:, 1:2);
            vert_1_xy = vert_1_xy(:, 1:2);
            grids_1_projected = {grids_1_projected{:}, {pts_1_xy, vert_1_xy}};
            
            tangent = pts_1_xy(end,:) - pts_1_xy(1, :);            
            tangent = tangent / sqrt(sum(tangent.^2));
            normal = [tangent(2), -tangent(1)];
            
            pts_2_xy = bsxfun(@plus, pts_1_xy, normal * gamma);
            pts_2_xy_neg = bsxfun(@minus, pts_1_xy, normal * gamma);

            internal_pts = [internal_pts; pts_2_xy; pts_2_xy_neg];
            
            % Now add an extra point at each end
            vert_s = vert_1_xy(1,:);
            vert_e= vert_1_xy(end,:);
            kappa_s = sum((vert_s - pts_1_xy(1,:)).^2, 2);
            kappa_e = sum((vert_e - pts_1_xy(end,:)).^2, 2);
            R_s = sqrt(kappa_s + gamma^2);
            R_e = sqrt(kappa_e + gamma^2);
            
            start_pts = vert_s - tangent*R_s;
            end_pts = vert_e + tangent*R_e;
            tip_sites = [tip_sites; start_pts; end_pts];

        end
    end

    % throw out tip sites otside fracture domain
    is_innside = true(1, size(tip_sites, 1));
    for i = 1:size(tip_sites, 1)
        is_innside(i) = inPolygon(f_xy, tip_sites(i,:));    
    end
    tip_sites = tip_sites(is_innside, :);
    tip_sites = faultSufCondFrom1D(tip_sites, grids_1_projected, gamma - 1e-6);
    rectangle = [min(f_xy); max(f_xy)];
	corners   = f_xy;
	fd        = @dpoly;
	vararg    = [f_xy; f_xy(1,:)];
    h = @(p, varargin) ds * ones(size(p, 1), 1);
    fixedPts = [internal_pts; tip_sites; corners];  
    [pts, ~, sorting] = distmesh2d(fd, h, ds, rectangle, fixedPts, false, vararg);
    % Distmesh change the order of all points. We undo this sorting.
    isInt = false(max(sorting),1); isInt(1:size(internal_pts,1)) = true;
    isInt = isInt(sorting);
    [~,If]   = sort(sorting(isInt));
    isRes   = ~isInt;

    fPts = pts(isInt,:);
    fPts = fPts(If,:);
        
    pts = faultSufCondFrom1D(pts(isRes, :), grids_1_projected, gamma-1e-6);
    %radius = 0.8 * ds * ones(size(fPts,1),1);
    pts = [fPts; pts];
    
    G = clippedPebi2D(pts, f_xy);
    G.nodes.coords = [G.nodes.coords, zeros(size(G.nodes.coords,1),1)];
    G.nodes.coords = G.nodes.coords*R + center;

    pts = [pts, zeros(size(pts,1),1)];
    %sites_2 = [sites_2; {pts * R + center, f}];
    G.cells.sites = pts * R + center;
    G.cells.resSite = true(size(pts,1),1);
    G.cells.resSite(1:numel(If))=false;
    grids_2 = [grids_2, G];
end
end


function [pts] = faultSufCondFrom1D(pts, grids_1, gamma)
for i = 1:numel(grids_1)
    c = grids_1{i}{1};
    c = [c; c(end,:)];
    v = grids_1{i}{2};
    kappaSqr = sum((c - v).^2, 2);
    R = sqrt(kappaSqr + gamma^2);
    pts = removeConflictPoints(pts, v, R);
end
end
