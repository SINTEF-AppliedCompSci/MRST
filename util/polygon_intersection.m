function int = polygon_intersection(poly_1, poly_2)
    tol = 1e-6;
    normal = normal_from_points(poly_1);
    center_1 = mean(poly_1, 1);
    poly_1 = bsxfun(@minus, poly_1, center_1);
    poly_2 = bsxfun(@minus, poly_2, center_1);
    R = rotation_matrix_from_plane(poly_1);
    poly_1 = poly_1 * R';
    assert(sum(abs(poly_1(:, 3)))<1e-6);
    poly_1 = poly_1(:,1:2);
    poly_2 = poly_2 * R';
    
    if all(poly_2(:, 3) > tol) || all(poly_2(:, 3) < -tol)
        int = zeros(0, 3);
        return
    elseif all(abs(poly_2(:, 3)) < tol)
        % We do not treat parallel polygons. Hope for the best
        int = zeros(0, 3);
        return
    end
    ind = [1:size(poly_2,1), 1];
    AB = [];
    for e = 1:numel(ind) - 1
        v1 = poly_2(ind(e),:);
        v2 = poly_2(ind(e + 1),:);
        I = edge_intersecting_plane([0,0,1], [0,0,0], v1, v2);
        if size(I, 1)==0
            continue
        end
        AB = [AB; I];
    end
    assert(size(AB,1)==2)
    assert(all(abs(AB(:,3))<tol))
    AB = AB(:, 1:2);
    A_innside = in_polygon(poly_1, AB(1,:));
    B_innside = in_polygon(poly_1, AB(2,:));
    AB_line = reshape(AB', 1,4);
    poly_1 = [poly_1; poly_1(1,:)];
    poly_1 = [poly_1(1:end-1,:),poly_1(2:end,  :)];
    [X, Y, segInt] = lineLineInt(AB_line, poly_1);
    
    if A_innside && B_innside
        int = AB;
    elseif A_innside && ~B_innside
        int = [AB(1,:); X(segInt), Y(segInt)];
    elseif ~A_innside && B_innside 
        int = [AB(2,:); X(segInt), Y(segInt)];
    elseif ~A_innside && ~B_innside
        int = [X(segInt), Y(segInt)];
    else
        assert False
    end
        
    assert(size(int, 1) == 2)
    
    int = [int, zeros(size(int,1),1)];
    
    int = int * R + center_1;
    
end
   
function inside = in_polygon(poly, ptn)
    tol = 1e-6;
    poly = poly - ptn;
    ind = [1:size(poly,1), 1];
    prev_side = 0;
    for e = 1:numel(ind) - 1
        v1 = poly(ind(e), :);
        v2 = poly(ind(e + 1), :);
        v1Xv2 = v2(1)*v1(2) - v1(1)*v2(2);
        current_side = sign(v1Xv2);
        if abs(v1Xv2) <= tol
            inside = false;
            return
        elseif prev_side == 0
            prev_side = current_side;
        elseif prev_side ~= current_side
            inside = false;
            return
        end
    end
    inside = true;
end

function val = x_product(p1, p2)
    val= p1(1)*p2(2) - p1(2)*p2(1);
end

function I = edge_intersecting_plane(normal, x0, v1, v2)
    sgn1 = sign(dot(v1, normal));
    sgn2 = sign(dot(v2, normal));
    if sgn1 == sgn2
        I = zeros(0,3);
        return
    end
    d1 = abs((v1 - x0)*normal');
    d2 = abs((v2 - x0)*normal');
    
    t = d1 / (d1 + d2);
    I = v1 + (v2 - v1) * t;    
end

function mat = rotation_matrix_from_plane(poly)
    normal = normal_from_points(poly);
    z_basis = [0,0,1];
    theta = acos(dot(normal, z_basis));
    v = cross(normal, z_basis);
    if sum(abs(v))<1e-8
        mat = eye(3);
        return
    end
    v = v/sqrt(sum(v.^2));
    mat_x = [     0, -v(3),  v(2);...
               v(3),     0, -v(1);...
              -v(2),  v(1),     0];
   mat_t = [v(1)^2, v(1)*v(2), v(1)*v(3); ...
            v(1)*v(2), v(2)^2, v(2)*v(3); ...
            v(1)*v(3), v(2)*v(3), v(3)^3];
   mat = cos(theta) * eye(3) + sin(theta) * mat_x + (1-cos(theta))*mat_t;
end

function normal = normal_from_points(p)
    assert(size(p,1)>=3)
    p1 = p(1, :);
    p2 = p(2,:);
    p_mean = mean(p, 1);
    v1 = p2 - p1;
    v2 = p_mean - p1;
    
    normal = cross(v1, v2);
    if sum(abs(normal))<1e-8
        normal = normal_from_points(p(:, 2:end));
        return
    end
    normal = normal / sqrt(sum(normal.^2));
end