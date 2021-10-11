function [P, faceNo, cellNo, w] = getFaceFromCellInterpolator(G, rock)
%Undocumented Utility Function

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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
