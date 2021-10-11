function [B, tbls] = computeLocalFluxMimetic(G, rock, varargin)
% Computes the scalar product matrix for the cell-face-node degrees of freedom
%
% SYNOPSIS:
%   function [B, tbls] = computeLocalFluxMimetic(G, rock, varargin)
%
% DESCRIPTION:
% 
%  Reference paper : 
%      title     = {Local flux mimetic finite difference methods},
%      author    = {Lipnikov, Konstantin and Shashkov, Mikhail and Yotov, Ivan},
%      journal   = {Numerische Mathematik},
%      volume    = {112},
%      number    = {1},
%      pages     = {115--152},
%      year      = {2009},
%      publisher = {Springer}
% PARAMETERS:
%   G        - Grid
%   rock     - Rock data structure (see description in `PermTensor`)
%   varargin - See below
%
% KEYWORD ARGUMENTS:
%   verbose       - true if verbose
%   ip_compmethod - Method that is used to compute the scalar product at a corner.
%                   The possible options for ip_compmethod are
%                     'general'       : general case, no special requirements on corner
%                     'nicecorner'    : case where at each corner, number of faces is equal to G.griddim
%                     'directinverse' : case where at each corner, number of faces is equal to
%                                       G.griddim AND eta is value such that N*R is diagonal (see Lipnikov reference paper)
%   eta           - scalar parameter which determines position of continuity point on face-node 
%                   (see definition of cellFacetVec in code, default value is zero)
%
% RETURNS:
%   B    - Scalar product matrix
%   tbls - Structures with IndexArrays
%
% SEE ALSO: 
%   `computeMultiPointTrans`.

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
    
    % Some  aliases 
    nc  = G.cells.num;
    nf  = G.faces.num;
    dim = G.griddim;

    coltbl.coldim = (1 : dim)';
    coltbl = IndexArray(coltbl);
    rowtbl = coltbl;
    rowtbl = replacefield(rowtbl, {'coldim', 'rowdim'});

    celltbl.cells = (1 : nc)';
    celltbl = IndexArray((celltbl));
    
    facetbl.faces = (1 : nf)';
    facetbl = IndexArray((facetbl));
    
    cellfacetbl.cells = rldecode((1 : nc)', diff(G.cells.facePos)); 
    cellfacetbl.faces = G.cells.faces(:, 1);
    cellfacetbl = IndexArray(cellfacetbl);

    facenodetbl.faces = rldecode((1 : nf)', diff(G.faces.nodePos)); 
    facenodetbl.nodes = G.faces.nodes;
    facenodetbl = IndexArray(facenodetbl);
    % We setup the face-node table and it is ordered along ascending node numbers so
    % that we will have a block structure for the nodal scalar product.
    facenodetbl = sortIndexArray(facenodetbl, {'nodes', 'faces'});
    
    cellnodefacetbl = crossIndexArray(cellfacetbl, facenodetbl, {'faces'});

    % We setup the cell-face-node table, cellnodefacetbl. Each entry determine a
    % unique facet in a corner
    % We order cellnodeface in cell-node-face order. This is done to optimize
    % for-end loop below.
    cellnodefacetbl = sortIndexArray(cellnodefacetbl, {'cells', 'nodes', 'faces'});
    cellnodefacetbl = cellnodefacetbl.addLocInd('cnfind');    
    
    % We setup the cell-node table, cellnodetbl. Each entry determine a unique
    % corner
    cellnodetbl = projIndexArray(cellnodefacetbl, {'cells', 'nodes'});
    % ordering to optimize for-end loop below
    cellnodetbl = sortIndexArray(cellnodetbl, {'cells', 'nodes'});
    
    % Nodal scalar product is stored in vector nodeM
    % mattbl is the table which specifies how nodeM is stored: a matrix for
    % each "corner" (cell-node pair).
    crossextend = {'faces', {'faces1', 'faces2'}};
    mattbl = crossIndexArray(cellnodefacetbl, cellnodefacetbl, {'cells', 'nodes'}, ...
                        'crossextend', {crossextend});
    % We order mattbl in cell-node-face1-face2 order
    % This is done to optimize for-end loop below
    mattbl = sortIndexArray(mattbl, {'cells', 'nodes', 'faces1', 'faces2'});
    
    nodeM = zeros(mattbl.num, 1);

    if opt.verbose
        fprintf('assemble facet normals ...\n');
    end
    
    fno = cellnodefacetbl.get('faces');
    cno = cellnodefacetbl.get('cells');
    numnodes = double(diff(G.faces.nodePos));
    numnodes = numnodes(fno);
    facetNormals = G.faces.normals(fno, :);
    facetNormals = bsxfun(@ldivide, numnodes, facetNormals);
    
    sgn = 2*(cno == G.faces.neighbors(fno, 1)) - 1;
    facetNormals = sgn.*facetNormals; % Outward normals with respect to cell
                                      % in cellnodeface.
    
    cellnodefacecoltbl = crossIndexArray(cellnodefacetbl, coltbl, {}, 'optpureproduct', true);
    cellnodefacecoltbl = sortIndexArray(cellnodefacecoltbl, {'cells', 'nodes', ...
                        'faces', 'coldim', 'cnfind'});
    facetNormals = reshape(facetNormals', [], 1);
    
    if opt.verbose
        fprintf('assemble facet K*normals ...\n');
    end
    % Assemble facePermNormals which corresponds to $Kn$ where n are the *outward*
    % normals at the facets.
    [perm, r, c] = permTensor(rock, G.griddim);
    permmat = perm;
    perm = reshape(permmat', [], 1);
    % setup cellcolrow table for the vector perm
    colrowtbl = crossIndexArray(coltbl, rowtbl, {});
    cellcolrowtbl = crossIndexArray(celltbl, colrowtbl, {}, 'optpureproduct', true);
    cellcolrowtbl = sortIndexArray(cellcolrowtbl, {'cells', 'coldim', ...
                        'rowdim'});
    cellcolrowtbl = cellcolrowtbl.addLocInd('ccrind');
    
    % Multiply perm with facetNormals
    prod = TensorProd();
    prod.tbl1 = cellcolrowtbl;
    prod.tbl2 = cellnodefacecoltbl;
    prod.tbl3 = cellnodefacecoltbl;
    prod.replacefds1 = {{'coldim', 'temp'}, {'rowdim', 'coldim'}, {'temp', 'rowdim'}};
    prod.replacefds2 = {'coldim', 'rowdim'};
    prod.mergefds = {'cells'};
    prod.reducefds = {'rowdim'};
    prod = prod.setup();
    
    Kn = prod.eval(perm, facetNormals);
   
    % store Kn in matrix form in facePermNormals.
    % Note that the indices are, by construction above, sorted.
    facePermNormals = reshape(Kn, coltbl.num, [])';
    
    % Some shortcuts
    cno = cellnodefacetbl.get('cells');
    fno = cellnodefacetbl.get('faces');
    nno = cellnodefacetbl.get('nodes');
    % Default option (opt.eta = 0): Use original face centroids and cell centroids,
    % NOT actual subface centroids. This corresponds to an MPFA method
    % (O-method) R = G.faces.centroids(fno,:) - G.cells.centroids(cno,:);
    cellFacetVec = G.faces.centroids(fno,:) - G.cells.centroids(cno,:) + ...
        opt.eta*(G.nodes.coords(nno,:) - G.faces.centroids(fno,:));
    
    % set up areas and volumes
    areas = G.faces.areas(fno);
    vols  = G.cells.volumes(cno);
    
    % number of faces per cell-nodes.
    map = TensorMap();
    map.fromTbl = cellnodefacetbl;
    map.toTbl = cellnodetbl;
    map.mergefds = {'cells', 'nodes'};
    map = map.setup();
    
    nfaces = map.eval(ones(cellnodefacetbl.num, 1));

    % we setup nfaces indexed along cellnodefacetbl
    map = TensorMap();
    map.fromTbl = cellnodetbl;
    map.toTbl = cellnodefacetbl;
    map.mergefds = {'cells', 'nodes'};
    map = map.setup();
    
    nfaces = map.eval(nfaces); 
    
    cnf_i = 1; % start indice for the cellnodefacetbl index
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
        faces = cellnodefacetbl.get('faces');
        faces = faces(cnfind);
        
        cellno = cellnodefacetbl.get('cells');
        cellno = cellno(cnf_i);
        
        K = reshape(permmat(cellno, :), [dim, dim]);
        
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
        nodeM(matind) = nodeM(matind) + locM;
        
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
    nodeM = nodeM.*sgn1.*sgn2;   
    
    if opt.verbose
        fprintf('Condensate on nodes ...\n');
    end    
    % Condensate on nodes (sum up cell contributions for given node).
    redmattbl = projIndexArray(mattbl, {'nodes', 'faces1', 'faces2'});
    map = TensorMap();
    map.fromTbl = mattbl;
    map.toTbl = redmattbl;
    map.mergefds = {'nodes', 'faces1', 'faces2'};
    map = map.setup();
    
    nodeM = map.eval(nodeM);
    
    if opt.verbose
        fprintf('Set up matrix ...\n');
    end    
    % Setup matrix
    % First set up facet indices in the redmattbl table
    fnind = (1 : facenodetbl.num)';
    
    map = TensorMap();
    map.fromTbl = facenodetbl;
    map.toTbl = redmattbl;
    map.replaceFromTblfds = {{'faces', 'faces1'}};
    map.mergefds = {'nodes', 'faces1'};
    
    facesind1 = map.getDispatchInd();
    
    map = TensorMap();
    map.fromTbl = facenodetbl;
    map.toTbl = redmattbl;
    map.replaceFromTblfds = {{'faces', 'faces2'}};
    map.mergefds = {'nodes', 'faces2'};
    
    facesind2 = map.getDispatchInd();
    
    % Assembly of B
    B = sparse(facesind1, ...
               facesind2, ...
               nodeM, ...
               facenodetbl.num, ...
               facenodetbl.num);
    
    tbls = struct('celltbl'        , celltbl        , ...
                  'facetbl'        , facetbl        , ...
                  'cellnodefacetbl', cellnodefacetbl, ...
                  'cellfacetbl'    , cellfacetbl    , ...
                  'cellnodetbl'    , cellnodetbl    , ...
                  'facenodetbl'    , facenodetbl);
    
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
