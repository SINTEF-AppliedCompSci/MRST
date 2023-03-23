function [H] = rodriguesRotation(G, phi, u, x0)
    W = [0,   -u(3),  u(2);
         u(3),    0, -u(1);
        -u(2), u(1),     0];
    R = eye(3) + sin(phi) * W + (2 * sin(phi / 2)^2) * W^2;
    H = G;
    for i = 1:G.nodes.num
        H.nodes.coords(i, :) = (R * (G.nodes.coords(i, :)' - x0'))' + x0;
    end
end