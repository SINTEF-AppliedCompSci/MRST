function [interpFace] = correctHAP(G, interpFace, myRatio)
% Correct ill-placed harmonic averaging points. If the number of input arguments
% is 2, then the correction algorithm is applied only when some cell centroids
% lie outside their associated convex hull; if the number of input arguments is
% 3, then the last input argument myRatio is applied to all the harmonic
% averaging points.
%
%  G - Grid structure of MRST
%  interpFace - harmonic averaging point interplation without correction
%  myRatio - user specified ratio

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

    dispif(mrstVerbose, 'correctHAP\n');            
    
    % the locations of the original harmonic averaging points
    HAP = interpFace.coords; 
    interpFace.corrected = interpFace.fraction;

    if (nargin == 2 || (nargin == 3 && isempty(myRatio)))
        if (interpFace.fraction > 0)
            if (G.griddim == 2)
                R = 0.5 * G.faces.areas;
            else
                R = sqrt(G.faces.areas./pi);
            end

            % Find all cells with centroids not in the convex hull
            out = find_cells(G, 1:G.cells.num, interpFace);
            
            % Correct
            cells = 1:G.cells.num;
            cells = cells(out);
            while sum(out)
                for c = cells
                    while out(c)
                        interpFace = correct(G, c, interpFace, R, HAP);
                        out(c) = find_cells(G, c, interpFace);
                    end
                    faces = G.cells.faces(G.cells.facePos(c):G.cells.facePos(c+1)-1);
                    neighs = G.faces.neighbors(faces, :);
                    neighs = neighs(:);
                    neighs = neighs(neighs > 0);
                    neighs = neighs(neighs ~= c);
                    idx = find_cells(G, neighs, interpFace);
                    neighs = neighs(idx);
                    cells = [cells, neighs'];
                    out(neighs) = 1;
                end
            end
        end
    elseif (nargin == 3)
        if (G.griddim == 2)
            R = 0.5 * G.faces.areas;
        else
            R = sqrt(G.faces.areas./pi);
        end
        R = R * myRatio;
        xf = G.faces.centroids;
        hap = interpFace.coords;
        d = hap - xf;
        d = sqrt(dot(d, d, 2));
        ind = find(d > R);
        dispif(mrstVerbose, 'correcting %d HAPs using myRatio=%1.1f\n', numel(ind), myRatio);
        interpFace.corrected = numel(ind)/G.faces.num;
        for i = 1:numel(ind)
            interpFace = correctHAP_local(G, ind(i), interpFace, HAP(ind(i), :)', myRatio);
        end
    else
        error('Wrong number of inputs')
    end
end


function out = find_cells(G, cells, interpFace)

    in = zeros(numel(cells), 1);
    for i = 1:numel(cells)
        c = cells(i);
        faces = G.cells.faces(G.cells.facePos(c):G.cells.facePos(c+1)-1);
        hap = interpFace.coords(faces, :);
        ind = convhull(hap);
        in(i) = mex_inhull(G.cells.centroids(c, :), hap, ind, -1e-5);
    end
    out = ~in;
end


function interpFace = correct(G, cells, interpFace, R, HAP)

    for c = cells % Row vector
        faces = G.cells.faces(G.cells.facePos(c):G.cells.facePos(c+1)-1);
        d = interpFace.coords(faces, :) - G.faces.centroids(faces, :);
        d = vecnorm(d, 2, 2);
        [maxRatio, ind] = max(d./R(faces));
        y_sigma = HAP(faces(ind), :)';
        interpFace = correctHAP_local(G, faces(ind), interpFace, y_sigma, 0.9*maxRatio);
    end

end


function interpFace = correctHAP_local(G, i_face, interpFace, y_sigma, myRatio)
% Correct harmonic averaging point for i_face based on given myRatio
    if (myRatio > 0)
        if (G.griddim == 2)
            R = 0.5 * G.faces.areas(i_face) * myRatio;
        elseif (G.griddim == 3)
            R = myRatio * sqrt(G.faces.areas(i_face)/pi);
        end
        xm = G.faces.centroids(i_face, :)';
        interpFace.coords(i_face, :) = (xm + R * (y_sigma - xm) / norm(y_sigma-xm))';
    else
        interpFace.coords(i_face, :) = G.faces.centroids(i_face, :);
    end
end
