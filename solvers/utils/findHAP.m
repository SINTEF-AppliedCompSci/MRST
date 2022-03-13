function interpFace = findHAP(G, rock)
% Find harmonic averaging points for 2D and 3D grids. Considering both
% Dirichlet and Neumann boundary conditions

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

% interpFace.coords: coordinates of interpolating points
% interpFace.weights: interpolating weights
% interpFace.fraction: the fraction of cells whose centroid is
% outside the convex hull
    dispif(mrstVerbose, 'findHAP\n');

    [K, ~, c] = permTensor(rock, G.griddim);
    interpFace.coords = zeros(G.faces.num, G.griddim);
    interpFace.weights = zeros(G.faces.num, 2);
    interpFace.fraction = 0;

    % Interior faces
    interior = all(G.faces.neighbors ~= 0, 2);
    c1 = G.faces.neighbors(interior, 1);
    c2 = G.faces.neighbors(interior, 2);

    n = G.faces.normals(interior, :);
    K1 = K(c1, :);
    K2 = K(c2, :);
    w1 = bsxfun(@times, K1, n(:, c));
    w2 = bsxfun(@times, K2, n(:, c));
    if G.griddim == 2
        sumfcn = @(w) [sum(w(:, 1:2), 2), sum(w(:, 3:4), 2)];
    else
        sumfcn = @(w) [sum(w(:, 1:3), 2), sum(w(:, 4:6), 2), sum(w(:, 7:9), 2)];
    end
    w1 = sumfcn(w1);
    w2 = sumfcn(w2);

    x1 = G.cells.centroids(c1, :);
    x2 = G.cells.centroids(c2, :);
    xf = G.faces.centroids(interior,:);
    xA = x1 + dot(xf-x1, n, 2) ./ dot(w1, n, 2) .* w1;
    xB = x2 + dot(xf-x2, n, 2) ./ dot(w2, n, 2) .* w2;
    w1 = vecnorm(w1, 2, 2) ./ vecnorm(xA-x1, 2, 2);
    w2 = vecnorm(w2, 2, 2) ./ vecnorm(xB-x2, 2, 2);

    interpFace.coords(interior, :) = (xA .* w1 + xB .* w2) ./ (w1 + w2);
    interpFace.weights(interior, 1) = w1;
    interpFace.weights(interior, 2) = w2;
    interpFace.weights(interior, :) = interpFace.weights(interior, :) ./ (w1 + w2);

    % Boundary faces as hom Neumann
    bdry = ~interior;
    c1 = G.faces.neighbors(bdry, 1);
    c2 = G.faces.neighbors(bdry, 2);
    cc = max(c1, c2);
    n = G.faces.normals(bdry, :);
    K1 = K(cc, :);
    w1 = bsxfun(@times, K1, n(:, c));
    w1 = sumfcn(w1);
    x1 = G.cells.centroids(cc, :);
    xf = G.faces.centroids(bdry, :);
    xA = x1 + dot(xf-x1, n, 2) ./ dot(w1, n, 2) .* w1;

    interpFace.coords(bdry, :) = xA;
    ex1 = G.faces.neighbors(:, 1) ~= 0;
    interpFace.weights(bdry & ex1, 1) = 1;
    interpFace.weights(bdry & ~ex1, 2) = 1;

    % Count the number of cells whose centroid is outside the convex hull-----
    counter = zeros(G.cells.num, 1);
    for i = 1:G.cells.num
        xc = G.cells.centroids(i, :);
        theFaces = G.cells.faces(G.cells.facePos(i):G.cells.facePos(i+1)-1);
        hap = interpFace.coords(theFaces, :);
        ind = convhull(hap);
        switch G.griddim
            case 2
                xv = hap(ind, 1);
                yv = hap(ind, 2);
                counter(i) = inpolygon(xc(1), xc(2), xv, yv);
            case 3
                counter(i) = mex_inhull(xc, hap, ind, -1e-5);
        end
    end

    interpFace.counter = logical(counter);
    interpFace.fraction = 1 - sum(counter) / G.cells.num;
end
