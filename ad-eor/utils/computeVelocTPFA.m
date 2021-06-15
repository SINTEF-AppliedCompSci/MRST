function veloc = computeVelocTPFA(G, intInx)
%
%
% SYNOPSIS:
%   function veloc = computeVelocTPFA(G, intInx)
%
% DESCRIPTION: 
%   Setup the operator which computes an approximation of the square
%   of the velocity for each cell given fluxes on the faces. 
%
%   We use a velocity reconstruction of the type  v_c  = 1/V * sum_{f} (x_f -x_c)u_f
%   where
%     v_c : Approximated alue of the velocity at the cell center
%     V   : Cell volume
%     f   : Face
%     x_f : Centroid of the face f
%     x_c : Centroid of the cell
%     u_f : flux at the face f
%    
%   Such reconstruction is exact for linear functions and first order accurate,
%   when a mimetic discretization is used or, in the case of TPFA, if the grid is
%   K-orthogonal. Note that it gives very large errors for cells that contain
%   well, as the pressure in such cells is only badly approximated by linear
%   functions.
%
%
% PARAMETERS:
%   G      - Grid structure
%   intInx - Logical vector giving the internal faces.
%
% RETURNS:
%   sqVeloc - Function of the form u=sqVeloc(v), which returns cell-valued
%   square of the velocity for given face-valued fluxes v.
%
% EXAMPLE:
%
% SEE ALSO: computeSqVelocTPFA
%

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

    cellNo = rldecode(1 : G.cells.num, diff(G.cells.facePos), 2).';
    dim    = size(G.nodes.coords, 2);

    % Define mapping from internal faces to half faces.
    nhf = size(G.cells.faces, 1); % number of half faces
    nf  = G.faces.num;            % number of faces
    nif = nnz(intInx);            % number of internal faces
    nc  = G.cells.num;            % number of cells
    fromIntfacesToFaces = sparse(find(intInx), 1 : nif, 1, nf, nif);
    sgn = 2*(cellNo == G.faces.neighbors(G.cells.faces(:,1), 1)) - 1;
    fromFacesToHalffaces = sparse(1 : nhf, G.cells.faces(:, 1), sgn, nhf, nf);
    fromIntfacesToHalffaces =  fromFacesToHalffaces*fromIntfacesToFaces;

    vol = G.cells.volumes;

    C = G.faces.centroids(G.cells.faces(:, 1), :) - G.cells.centroids(cellNo, :);

    sumHalffaces = sparse(cellNo, 1 : nhf, 1, nc, nhf);
    veloc = cell(dim, 1);
    for i = 1 : dim
        veloc{i} = @(v) 1./vol.*(sumHalffaces*(C(:, i).*(fromIntfacesToHalffaces*v)));
    end

end