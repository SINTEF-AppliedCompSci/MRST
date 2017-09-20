function grids_2 = fault_sites(faults, grids_1, intersections, ds)
grids_2 = {};
for f = 1:numel(faults)
    f_pts = faults{f};
    center = mean(f_pts, 1);
    f_center = bsxfun(@minus, f_pts, center);
    R = rotation_matrix_from_plane(f_center);
    f_xy = f_center * R';
    assert(sum(abs(f_xy(:, 3)))<1e-6);
    f_xy = f_xy(:,1:2);
        
    internal_pts = [];
    for j = 1:numel(intersections)
        if intersections{j}{2}==f || intersections{j}{3}==f 
            pts_1 = grids_1{j};
            pts_1_center = bsxfun(@minus, pts_1, center);
            pts_1_xy = pts_1_center * R';
            assert(all(abs(pts_1_xy(:,3))<1e-6))
            pts_1_xy = pts_1_xy(:, 1:2);
            tangent = pts_1_xy(1,:) - pts_1_xy(end, :);
            normal = [tangent(2), -tangent(1)];
            normal = normal / sqrt(sum(normal.^2));
            pts_2_xy = bsxfun(@plus, pts_1_xy, normal * ds / 2);
            pts_2_xy_neg = bsxfun(@minus, pts_1_xy, normal * ds / 2);

            internal_pts = [internal_pts; pts_2_xy; pts_2_xy_neg];
       end
    end

    rectangle = [min(f_xy); max(f_xy)];
	corners   = f_xy;
	fd        = @dpoly;
	vararg    = [f_xy; f_xy(1,:)];
    h = @(p, varargin) ds * ones(size(p, 1), 1);
    fixedPts = [internal_pts; corners];  
    [pts, ~, sorting] = distmesh2d(fd, h, ds, rectangle, fixedPts, vararg);
    % Distmesh change the order of all points. We undo this sorting.
    isInt = false(max(sorting),1); isInt(1:size(internal_pts,1)) = true;
    isInt = isInt(sorting);
    [~,If]   = sort(sorting(isInt));
    isRes   = ~isInt;

    fPts = pts(isInt,:);
    fPts = fPts(If,:);
        
    radius = 0.8 * ds * ones(size(fPts,1),1);
    pts = removeConflictPoints2(pts(isRes, :), fPts, radius);
    pts = [fPts; pts];
    
%     %% Debugging
%     G = clippedPebi2D(pts, f_xy);
%     figure()
%     hold on
%     plotGrid(G)
%     plot(internal_pts(:,1), internal_pts(:,2), '.','markersize', 15)
%     axis equal
%     %% Debugging end
    pts = [pts, zeros(size(pts,1),1)];
    grids_2 = [grids_2; {pts * R + center, f}];
end

end