function int = polygonIntersection(poly_1, poly_2)
    tol = 1e-6;
    center_1 = mean(poly_1, 1);
    poly_1 = bsxfun(@minus, poly_1, center_1);
    poly_2 = bsxfun(@minus, poly_2, center_1);
    R = rotationMatrixFromPlane(poly_1);
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
    % Calculate the intersection of each line segment and the plane
    for e = 1:numel(ind) - 1
        v1 = poly_2(ind(e),:);
        v2 = poly_2(ind(e + 1),:);
        I = edge_intersecting_plane([0,0,1], [0,0,0], v1, v2);
        if size(I, 1)==0
            continue
        end
        AB = [AB; I];
    end
    % If a polygon node lies exactly at the plane both edges connected to
    % the node will create an intersection. Remove duplicates
    AB = round(AB * 1e10)/1e10;
    AB = unique(AB, 'rows');
    % The if statement above ensures that the polygon crosses the plane.
    % Our convex polygon should therefore have exactly two intersections
    % with this plane.
    assert(size(AB,1)==2)
    assert(all(abs(AB(:,3))<tol))
    AB = AB(:, 1:2);
    
    % Now test if the two intersections are innside the first polygon
    A_innside = inPolygon(poly_1, AB(1,:));
    B_innside = inPolygon(poly_1, AB(2,:));
    AB_line = reshape(AB', 1,4);
    poly_1 = [poly_1; poly_1(1,:)];
    poly_1 = [poly_1(1:end-1,:),poly_1(2:end,  :)];
    [X, Y, segInt] = lineLineInt(AB_line, poly_1);
    
    if A_innside && B_innside
        int = AB;
    elseif A_innside && ~B_innside
        int = [AB(1,:); X(segInt), Y(segInt)];
    elseif ~A_innside && B_innside 
        int = [AB(2,:); X(segInt)', Y(segInt)'];
    elseif ~A_innside && ~B_innside
        int = [X(segInt)', Y(segInt)'];
    else
        assert False
    end
    % remove duplicates
    int = round(int * 1e10)/1e10;
    int = unique(int, 'rows');   
    assert(size(int, 1) == 2 | size(int, 1) == 0)
    
    int = [int, zeros(size(int,1),1)];
    
    int = int * R + center_1;
    
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



