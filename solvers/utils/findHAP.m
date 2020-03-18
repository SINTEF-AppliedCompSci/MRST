function interpFace = findHAP(G, rock, bc)
% Find harmonic averaging points for 2D and 3D grids. Considering both
% Dirichlet and Neumann boundary conditions

% interpFace.coords: coordinates of interpolating points
% interpFace.weights: interpolating weights
% interpFace.fraction: the fraction of cells whose centroid is
% outside the convex hull
    
    dispif(mrstVerbose, 'findHAP\n');
    K = permTensor(rock, G.griddim);
    K = reshape(K', G.griddim, G.griddim, []);
    interpFace.coords = zeros(G.faces.num, G.griddim);
    interpFace.weights = zeros(G.faces.num, 2);
    interpFace.fraction = 0;

    % Find harmonic averaging point--------------------------------------------
    for i_face = 1 : G.faces.num
        
        c1 = G.faces.neighbors(i_face, 1);
        c2 = G.faces.neighbors(i_face, 2);
        xf = G.faces.centroids(i_face, :)';
        if (all([c1, c2] ~= 0))
            K1 = K(:, :, c1);
            K2 = K(:, :, c2);
            fn = G.faces.normals(i_face, :)';
            w1 = K1 * fn;
            w2 = K2 * fn;
            x1 = G.cells.centroids(c1, :)';
            x2 = G.cells.centroids(c2, :)';
            xA = x1 + dot(xf-x1, fn) / dot(w1, fn) * w1;
            xB = x2 + dot(xf-x2, fn) / dot(w2, fn) * w2;
            w1 = norm(w1) / norm(xA-x1);
            w2 = norm(w2) / norm(xB-x2);
            interpFace.coords(i_face, :) = (w1 * xA + w2 * xB)' / (w1 + w2);
            interpFace.weights(i_face, 1) = w1 / (w1 + w2);
            interpFace.weights(i_face, 2) = w2 / (w1 + w2);
        else
            ind = find(bc.face == i_face, 1);
            %ind,bc.type{ind}
            if (strcmpi(bc.type{ind}, 'pressure'))
                interpFace.coords(i_face, :) = xf';
                interpFace.weights(i_face, (c2 == 0)+1) = bc.value{ind}(xf);
            else
                c = max(c1, c2);
                K1 = K(:, :, c);
                fn = G.faces.normals(i_face, :)';
                w1 = K1 * fn;
                x1 = G.cells.centroids(c, :)';
                xA = x1 + dot(xf-x1, fn) / dot(w1, fn) * w1;
                interpFace.coords(i_face, :) = xA';
                a = norm(w1) / norm(x1-xA);
                gN = bc.value{ind}(xf) * G.faces.areas(i_face);
                interpFace.weights(i_face, (c1 == 0)+1) = 1;
                interpFace.weights(i_face, (c2 == 0)+1) = -gN / a;
                %interpFace.weights(i_face,(c2==0)+1)=bc.value{ind}(xf); %-gN/a;
            end
        end
    end

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
    interpFace.fraction = 1 - sum(counter) / G.cells.num;
end