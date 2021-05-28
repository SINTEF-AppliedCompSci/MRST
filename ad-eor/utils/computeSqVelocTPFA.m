function sqVeloc = computeSqVelocTPFA(G, intInx)
%
%
% SYNOPSIS:
%   function sqVeloc = computeSqVelocTPFA(G, intInx)
%
% DESCRIPTION: Setup the operator which computes an approximation of the square
% of the velocity for each cell given fluxes on the faces. The following
% method is used
%
% 1) For each face, we compute a square velocity by taking the square of the
% flux after having divided it by the area. We obtain a scalar-valued variable
% per face
%
% 2) From a scalar variable, a vector can be reconstructed on each face, using
% the coordinate of the vector joining the cell centroid to the face
% centroid. Let c denotes this vector, after renormalization, (example:
% c=[1;0;0]' for a face of a Cartesian 3D grid) then, from a scalar value f on
% the face, we construct the vector-valued variable given by f*c.
%
% 3) We use the reconstruction defined in 2) to construct vector-valued
% square velocity from the scalar-valued square velocity defined in 1).
%
% 4) We take the average of the vector-valued obtained in 3) . We use a
% weight in the average so that faces that are far from the cell centroid
% contribute less than those that are close. See implementation for how this
% is done.
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
% SEE ALSO: computeVelocTPFA
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

    sumHalffaces = sparse(cellNo, 1 : nhf, 1, nc, nhf);
    wSumHalffaces = sumHalffaces*sparse(1 : nhf, 1 : nhf, 1./(G.faces.areas(G.cells.faces(:, 1)).^2), ...
                                        nhf, nhf);

    C = G.faces.centroids(G.cells.faces(:, 1), :) - G.cells.centroids(cellNo, :);
    C = abs(C);
    Csum = sumHalffaces*C;
    C = C./Csum(cellNo, :);

    for i = 1: dim
        D{i} = wSumHalffaces*sparse(1 : nhf, 1 : nhf, C(:, i), nhf, nhf);
    end

    sqVeloc = @(v) computeSqVeloc(v, D, fromIntfacesToHalffaces);

end

function sqVeloc = computeSqVeloc(v, D, fromIntfacesToHalffaces)

    hf_v = fromIntfacesToHalffaces*v;
    hf_sq_v = hf_v.^2;

    dim = size(D, 2);
    sqVeloc = 0;
    for i = 1 : dim
        sqVeloc = sqVeloc + D{i}*hf_sq_v;
    end

end