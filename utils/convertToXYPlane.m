function [pts1, pts2, T, R, nor_z] = convertToXYPlane(pts1, n1, pts2, varargin)
% Convert the points p from fully 3D plane to horizontal xy plane. The 
% fully 3D plane is specified by pts1(n1(1), :), pts1(n1(2), :), and 
% pts1(n1(3), :).
% New z-axis: Along normals of the 3D plane
% New x-axis: pts1(n1(2), :) - pts1(n1(1), :)
% All points of pts1 and pts2 will be transformed.
% Optional: 'normalZ', provide z-normal of the plane

    opt = struct('normalZ', []);
    opt = merge_options(opt, varargin{:});
    
    % Get normal z
    if ~isempty(opt.normalZ)
        nor_z = opt.normalZ;
    else
        pFacez = pts1(n1(1:3), :);
        nor_z  = faceNormals(pFacez);
    end

    % Get normal x
    nor_x = pts1(n1(2), :) - pts1(n1(1), :);
    nor_x = nor_x / norm(nor_x, 2);

    % Get normal y
    pFacey = [nor_x; nor_z; [0, 0, 0]];
    nor_y  = faceNormals(pFacey);

    % Base point
    % x0 = min( pts1(n1, 1) );
    % y0 = min( pts1(n1, 2) );
    % z0 = min( pts1(n1, 3) );
    
    % Set base point yo the origin
    [x0, y0, z0] = deal(0);

    % Shift
    T = eye(4);
    T(4, [1, 2, 3]) = [-x0, -y0, -z0];

    % Rotate
    R = [nor_x', nor_y',nor_z'];
    R(4,4) = 1;

    % Transform
    pts1_extend = [pts1, ones(size(pts1,1), 1)];
    pts1 = pts1_extend * T * R;
    pts1 = pts1(:, 1:3);

    pts2_extend = [pts2, ones(size(pts2,1), 1)];
    pts2 = pts2_extend * T * R;
    pts2 = pts2(:, 1:3);
end

function nor = faceNormals(p)
    p1 = p(1, :);
    p2 = p(2, :);
    p3 = p(3, :);

    x1 = p1(1); x2 = p2(1); x3 = p3(1);
    y1 = p1(2); y2 = p2(2); y3 = p3(2);
    z1 = p1(3); z2 = p2(3); z3 = p3(3);

    a = (y2-y1)*(z3-z1) - (y3-y1)*(z2-z1);
    b = (z2-z1)*(x3-x1) - (z3-z1)*(x2-x1);
    c = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1);

    nor = [a, b, c];
    nor = nor / norm(nor, 2);
end