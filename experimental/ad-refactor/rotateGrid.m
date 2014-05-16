function G = rotateGrid(G, theta, center, axis_dir)

    X = G.nodes.coords';
    
    G.nodes.coords = translate(rotate(translate(X, -center), theta, axis_dir), center)';

end

% ----------------------------------------------------------------------------
function res = translate(X, v)
% Translate each column vector in X by v
    res = bsxfun(@plus, X, v);
end

% ----------------------------------------------------------------------------

function res = rotate(X, theta, axis_dir)
% Rotate each column vector in X theta radians around the axis passing
% through the origin with direction 'axis_dir'.
    
% constructing rotation quaternion
    axis_dir = sin(theta/2) * (axis_dir ./ norm(axis_dir));
    q    = [cos(theta/2);  reshape(axis_dir, [], 1)];
    qinv = [cos(theta/2); -reshape(axis_dir, [], 1)];
    
% making quaternions of X, performing the transformation, and extracting
% coordinates 
    res = bsxfun(@quatmul, q, [zeros(1, size(X,2)); X]);
    res = bsxfun(@quatmul, res, qinv);
    
    res = res(2:4,:);
end

function q = quatmul(l, r)
    
    q = zeros(4, 1);    % (w, x, y, z)

    q(1) = l(1) * r(1) - (l(2:4)' * r(2:4));
    q(2) = l(1) * r(2) + l(2) * r(1) + l(3) * r(4) - l(4) * r(3);
    q(3) = l(1) * r(3) - l(2) * r(4) + l(3) * r(1) + l(4) * r(2);
    q(4) = l(1) * r(4) + l(2) * r(3) - l(3) * r(2) + l(4) * r(1);
    
end


