function mat = rotationMatrixFromVector(theta, vec)
    if sum(abs(vec))<1e-8
        mat = eye(3);
        return
    end
    v = vec/sqrt(sum(vec.^2));
    mat_x = [     0, -v(3),  v(2);...
               v(3),     0, -v(1);...
              -v(2),  v(1),     0];
   mat_t = [v(1)^2, v(1)*v(2), v(1)*v(3); ...
            v(1)*v(2), v(2)^2, v(2)*v(3); ...
            v(1)*v(3), v(2)*v(3), v(3)^3];
   mat = cos(theta) * eye(3) + sin(theta) * mat_x + (1-cos(theta))*mat_t;
end