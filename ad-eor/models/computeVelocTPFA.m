function veloc = computeVelocTPFA(G, intInx)

% Compute approximation of velocity for TPFA
%
% Note:
%
% * The approximation of acceptable when flow is closed to linear. In particular, it is very bad
% when there is well cells.
%
% * We assume that there is no boundary flow.
%
% intInx is a logical vector giving the internal faces.

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