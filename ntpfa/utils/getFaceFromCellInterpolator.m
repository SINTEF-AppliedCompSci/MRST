function [P, faceNo, cellNo, w] = getFaceFromCellInterpolator(G, rock)
    % Get nodes
    faceNo = G.cells.faces(:, 1);
    pos = G.cells.facePos;
    cellNo = rldecode((1:G.cells.num)', diff(pos));
    
    L = (G.cells.centroids(cellNo, :) - G.faces.centroids(faceNo, :));
%     for i = 1:size(L, 1)
%         l = G.faces.normals(faceNo(i), :)*getK(rock, cellNo(i), G.griddim);
%         L(i, :) = L(i, :)*norm(l);
%     end
    % One over distance weighting
    w = 1./sqrt(sum(L.^2, 2));
    
    % Divide by sum of weights for partition of unity
    sumw = accumarray(faceNo, w);
    w = w./sumw(faceNo);
    
    % n_n by n_c matrix constructing values on nodes by cells
    P = sparse(faceNo, cellNo, w, G.faces.num, G.cells.num);
end

function K = getK(rock, cell, dim)
    k = rock.perm(cell, :);
    switch numel(k)
        case 1
            % Scalar perm
            K = k;
        case dim
            % Diagonal tensor
            K = diag(k);
        case 3*(dim - 1)
            % Full symmetric tensor
            if dim == 2
                K = [k(1), k(2); ...
                     k(2), k(3)];
            else
                K = [k(1), k(2), k(3); ...
                     k(2), k(4), k(5); ...
                     k(3), k(5), k(6)];
            end
        otherwise
            error('What sorcery is this?!');
    end
end
