function int_u = faceIntegral(u, cubature, faces)
    % Integrate integrand u over cells using a given cubature
    % int_u(i) = (int_{face(i)} u ds)/|face(i)|
    if nargin < 3 || isinf(faces)
        % Empty cells means all internal connections of the grid
        faces = find(disc.internalConn)';
    end
    % Get cubature for all cells, transform coordinates to ref space
    W = cubature.getCubature(faces, 'face');
    int_u = W*u;
end