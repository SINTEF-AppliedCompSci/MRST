function mat = rotationMatrixFromPlane(poly)
    normal = normalFromPoints(poly);
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