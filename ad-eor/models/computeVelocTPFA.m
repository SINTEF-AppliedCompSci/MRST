function veloc = computeVelocTPFA(G, intInx)
%
%
% SYNOPSIS:
%   function veloc = computeVelocTPFA(G, intInx)
%
% DESCRIPTION: Setup the operator which computes an approximation of the square
% of the velocity for each cell given fluxes on the faces. 
%
% We use a velocity reconstruction of the type  v_c  = 1/V * sum_{f} (x_f -x_c)u_f
% where
%    v_c : Approximated alue of the velocity at the cell center
%    V   : Cell volume
%    f   : Face
%    x_f : Centroid of the face f
%    x_c : Centroid of the cell
%    u_f : flux at the face f
%    
% Such reconstruction is exact for linear functions and first order accurate,
% when a mimetic discretization is used or, in the case of TPFA, if the grid is
% K-orthogonal. Note that it gives very large errors for cells that contain
% well, as the pressure in such cells is only badly approximated by linear
% functions.
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