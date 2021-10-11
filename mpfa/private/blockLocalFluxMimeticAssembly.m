function [B, tbls] = blockLocalFluxMimeticAssembly(G, rock, globtbls, nodes, varargin)
% Computes the scalar product matrix for the cell-face-node degrees of freedom
% for the block given by nodes
%
% SYNOPSIS:
%   function [B, tbls] = blockLocalFluxMimeticAssembly(G, rock, globtbls, nodes, varargin)
%
% DESCRIPTION: See `computeLocalFluxMimetic`
%
% PARAMETERS:
%   G        - Grid
%   rock     - Rock data structure (see description in `PermTensor`)
%   globtbls - Global IndexArrays, for the whole domain
%   nodes    - A set of nodes that is used to set up the block. 
%              The matrix B will contain all the face-node degrees of freedom for which the node belongs to the block.
%   varargin - See below
%
% KEYWORD ARGUMENTS:
%   verbose       - true if verbose
%   ip_compmethod - Option sent to `computeLocalFluxMimetic`
%   eta           - Option sent to `computeLocalFluxMimetic`
%
% RETURNS:
%   B    - Scalar product matrix for the block
%   tbls - Structures with IndexArrays
%
% EXAMPLE:
%
% SEE ALSO:
% `computeLocalFluxMimetic`, `private/blockComputeMultiPointTrans`, `private/blockComputeNeumannMultiPointTrans`
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

    opt = struct('verbose'      , mrstVerbose, ...
                 'ip_compmethod', 'general'  , ...
                 'eta'          , 0);
    opt = merge_options(opt, varargin{:});

    nodetbl.nodes = nodes;
    nodetbl = IndexArray(nodetbl);

    nc  = G.cells.num;
    nf  = G.faces.num;
    dim = G.griddim;
    
    cellfacetbl     = globtbls.cellfacetbl;
    facenodetbl     = globtbls.facenodetbl;
    cellfacenodetbl = globtbls.cellfacenodetbl;
    cellnodetbl     = globtbls.cellnodetbl;

    % Restrict cellnodetbl to the nodes blocks (nodes that are sent as argument)
    cellnodetbl = crossIndexArray(cellnodetbl, nodetbl, {'nodes'});
    cellnodetbl = sortIndexArray(cellnodetbl, {'cells', 'nodes'});
    
    if cellnodetbl.num == 0
        % may happens when grid contains nodes that do not belong to any faces.
        B = [];
        tbls = [];
        return
    end
    % Restrict cellfacenodetbl to the node blocks. Each entry of cellfacenodetbl
    % determine a unique facet in a corner
    cellfacenodetbl = crossIndexArray(cellnodetbl, cellfacenodetbl, {'nodes', ...
                        'cells'});
    % We order cellnodeface in cell-node-face order. This is done to optimize the
    % for-end loop below.
    cellfacenodetbl = sortIndexArray(cellfacenodetbl, {'cells', 'nodes', 'faces'});
    cellfacenodetbl = cellfacenodetbl.addLocInd('cnfind');

    % We set up the facenode table. It is reordered so that we obtain a diagonal
    % matrix in the local indices.
    facenodetbl = projIndexArray(cellfacenodetbl, {'nodes', 'faces'});
    facenodetbl = sortIndexArray(facenodetbl, {'nodes', 'faces'});
    facenodetbl = facenodetbl.addLocInd('fnind');

    % We set up cellface table
    cellfacetbl = projIndexArray(cellfacenodetbl, {'cells', 'faces'});

    % col and row tables for matrices in the spatial dimension (dim).
    coltbl.coldim = (1 : dim)';
    coltbl = IndexArray(coltbl);
    rowtbl = coltbl;
    rowtbl = replacefield(rowtbl, {'coldim', 'rowdim'});

    % We set up cell table (only cells that belong to the block).
    celltbl = projIndexArray(cellnodetbl, {'cells'});
    celltbl = celltbl.addLocInd('cind');

    % Nodal scalar product is stored in vector B. The table mattbl
    % specifies how B is stored: a matrix for each "corner" (cell-node
    % pair).
    crossextend = {'faces', {'faces1', 'faces2'}};
    mattbl= crossIndexArray(cellfacenodetbl, cellfacenodetbl, {'cells', 'nodes'}, ...
                            'crossextend', {crossextend});
    % We order mattbl in cell-node-face1-face2 order. This is done to optimize
    % for-end loop below
    mattbl = sortIndexArray(mattbl, {'cells', 'nodes', 'faces1', 'faces2'});

    B = zeros(mattbl.num, 1);

    if opt.verbose
        fprintf('assemble facet normals ...\n');
    end

    % We compute the product of the permeability tensor K and the normal at each
    % node-face. The normal is weighted by area divided number of nodes the face
    % is dealt with.
    fno = cellfacenodetbl.get('faces');
    cno = cellfacenodetbl.get('cells');
    numnodes = double(diff(G.faces.nodePos));
    numnodes = numnodes(fno);
    facetNormals = G.faces.normals(fno, :);
    facetNormals = bsxfun(@ldivide, numnodes, facetNormals);

    sgn = 2*(cno == G.faces.neighbors(fno, 1)) - 1;
    facetNormals = sgn.*facetNormals; % Outward normals with respect to cell
                                      % in cellnodeface.

    cellnodefacecoltbl = crossIndexArray(cellfacenodetbl, coltbl, {}, 'optpureproduct', true);
    cellnodefacecoltbl = sortIndexArray(cellnodefacecoltbl, {'cells', 'nodes', ...
                        'faces', 'coldim', 'cnfind'});
    
    facetNormals = reshape(facetNormals', [], 1);

    if opt.verbose
        fprintf('assemble facet K*normals ...\n');
    end
    % Assemble facePermNormals which corresponds to $Kn$ where n are the *outward*
    % (weighted) normals at the facets
    [fullpermmat, r, c] = permTensor(rock, G.griddim);
    cells = celltbl.get('cells');
    permmat = fullpermmat(cells, :);
    perm = reshape(permmat', [], 1);
    % Setup cellcolrow table for the vector perm
    colrowtbl = crossIndexArray(coltbl, rowtbl, {}, 'optpureproduct', true);
    cellcolrowtbl = crossIndexArray(celltbl, colrowtbl, {}, 'optpureproduct', true);
    % cellcolrowtbl = sortIndexArray(cellcolrowtbl, {'cells', 'coldim', 'rowdim'});
    cellcolrowtbl = cellcolrowtbl.addLocInd('ccrind');

    prod = TensorProd();
    prod.tbl1 = cellcolrowtbl;
    prod.tbl2 = cellnodefacecoltbl;
    prod.tbl3 = cellnodefacecoltbl;
    prod.replacefds1 = {{'coldim', 'rowdim', 'interchange'}};
    prod.replacefds2 = {{'coldim', 'rowdim'}};
    prod.reducefds = {'rowdim'};
    prod.mergefds = {'cells'};
    prod = prod.setup();
    
    Kn = prod.eval(perm, facetNormals);
    
    % store Kn in matrix form in facePermNormals.
    % we know that 'coldim' is last in the ordering of the multi-index. Hence,
    facePermNormals = reshape(Kn, coltbl.num, [])';

    % Some shortcuts
    cno = cellfacenodetbl.get('cells');
    fno = cellfacenodetbl.get('faces');
    nno = cellfacenodetbl.get('nodes');
    % Default option (opt.eta = 0): Use original face centroids and cell centroids,
    % NOT actual subface centroids. This corresponds to an MPFA method
    % (O-method) R = G.faces.centroids(fno,:) - G.cells.centroids(cno,:);
    cellFacetVec = G.faces.centroids(fno,:) - G.cells.centroids(cno,:) + ...
        opt.eta*(G.nodes.coords(nno,:) - G.faces.centroids(fno,:));

    % Set up areas and volumes
    areas = G.faces.areas(fno);
    vols  = G.cells.volumes(cno);

    map = TensorMap();
    map.fromTbl = cellfacenodetbl;
    map.toTbl = cellnodetbl;
    map.mergefds = {'cells', 'nodes'};
    map = map.setup();
    
    nfaces = map.eval(ones(cellfacenodetbl.num, 1));
    
    map = TensorMap();
    map.fromTbl = cellnodetbl;
    map.toTbl = cellfacenodetbl;
    map.mergefds = {'cells', 'nodes'};
    map = map.setup();

    nfaces = map.eval(nfaces);

    cnf_i = 1; % start indice for the cellfacenodetbl index
    mat_i = 1; % start indice for the mattbl index

    for i = 1 : cellnodetbl.num

        % if opt.verbose
            % t0 = tic();
        % end

        nface = nfaces(cnf_i);
        cnfind = cnf_i : (cnf_i + (nface - 1));

        N     = facePermNormals(cnfind, :);
        R     = cellFacetVec(cnfind, :);
        a     = areas(cnfind); % areas of the faces the facets belong to
        v     = vols(cnf_i); % volume of the current cell
        faces = cellfacenodetbl.get('faces');
        faces = faces(cnfind);

        cellno = cellfacenodetbl.get('cells');
        cellno = cellno(cnf_i);
        K = reshape(fullpermmat(cellno, :), [dim, dim]);

        % Assemble local nodal scalar product
        switch opt.ip_compmethod
          case 'general'
            locM = node_ip(a, v, full(N), full(R), K);
          case 'nicecorner'
            locM = node_ip2(full(N), full(R));
          case 'directinverse'
            volE = G.cells.volumes(cellno);
            locM = node_ip3(full(N), K, G.griddim, volE);
          otherwise
            error('ip_compmethod not recognized');
        end

        locM = reshape(locM', [], 1);

        matind = mat_i : (mat_i + (nface*nface - 1));
        B(matind) = B(matind) + locM;

        % increment start indices
        cnf_i = cnf_i + nface;
        mat_i = mat_i + nface*nface;

        % if opt.verbose
            % t0 = toc(t0);
            % fprintf('assembly cellnode %d / %d took %g sec\n', i, cellnodetbl.num, ...
                    % t0);
        % end
    end

    if opt.verbose
        fprintf('include face sign ...\n');
    end
    
    cells = mattbl.get('cells');
    faces1 = mattbl.get('faces1');
    faces2 = mattbl.get('faces2');
    
    sgn1 = 2*(cells == G.faces.neighbors(faces1, 1)) - 1;
    sgn2 = 2*(cells == G.faces.neighbors(faces2, 1)) - 1;
    B = B.*sgn1.*sgn2;

    if opt.verbose
        fprintf('Condensate on nodes ...\n');
    end
    % Condensate on nodes (sum up cell contributions for given node).
    crossextfds = {{'faces', {'faces1', 'faces2'}}, {'fnind', {'fnind1', ...
                        'fnind2'}}};
    face2nodetbl = crossIndexArray(facenodetbl, facenodetbl, {'nodes'}, ...
                                                   'crossextend', crossextfds);
    
    map = TensorMap();
    map.fromTbl = mattbl;
    map.toTbl = face2nodetbl;
    map.mergefds = {'nodes', 'faces1', 'faces2'};
    map = map.setup();
    
    B = map.eval(B);

    if opt.verbose
        fprintf('Set up matrix ...\n');
    end

    tbls = struct('face2nodetbl'   , face2nodetbl   , ...
                  'cellfacenodetbl', cellfacenodetbl, ...
                  'cellfacetbl'    , cellfacetbl    , ...
                  'cellnodetbl'    , cellnodetbl    , ...
                  'facenodetbl'    , facenodetbl    , ...
                  'celltbl'        , celltbl);

end

function M = node_ip(a, v, N, R, K)
% a : areas of the facets
% v : volume of the corner (for the moment we use volume of cell)
% N : permeability*(facets' normals), corresponds to $\tilde N_c$ in paper of
%     Lipnikov et al (2009)
% R : vector of cell's to facets' centroids, corresponds to $R_c$ in paper of
%     Lipnikov et al (2009)

    k = 6 * sum(diag(K)) / size(K, 2);
    scalfact = k*sum(a)/numel(a);
    N = (1/scalfact)*N;
    [U, D, V] = svd(N);
    D = scalfact*D;
    fnum = size(N, 1);
    dim = size(N, 2);

    Dp = D(1 : dim, 1 : dim);
    d = diag(Dp);
    assert(prod(d)>0, 'cannot assemble mpfa, need extra fix');
    invd = 1./d;
    H = [diag(invd), zeros(dim, fnum - dim)];
    M = R*V*H*U';

    % Add stabilization term S
    S = blkdiag(zeros(dim), eye(fnum - dim));
    U = diag(a)*U;
    t = 6 * v^(5/3) * sum(diag(K)) / size(K, 2);
    regM = (1/t)*U*S*U';
    % fprintf('norm main: %g, norm reg: %g\n', norm(M), norm(regM));
    M = M + regM;
    % fprintf('condition number: %g\n', condest(M));
end

function M = node_ip2(N, R)
    M = R*inv(N);
end

function M = node_ip3(N, K, d, volE)
% volE : volume of cell
% d    : Spatial dimension (2 or 3)
% K    : permeability tensor
    invN = inv(N);
    switch d
      case 2
        mE = 3;
      case 3
        mE = 4;
      otherwise
        error('wrong spatial dimension (should be equal to 2 or 3)');
    end

    M = (volE/mE)*(invN'*K*invN);
end
