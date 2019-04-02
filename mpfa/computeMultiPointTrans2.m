function mpfastruct = computeMultiPointTrans2(G, rock, varargin)
%Compute multi-point transmissibilities.
%
% SYNOPSIS:
%   T = computeMultiPointTrans(G, rock)
%
% REQUIRED PARAMETERS:
%   G       - Grid data structure as described by grid_structure.
%
%   rock    - Rock data structure with valid field 'perm'.  The
%             permeability is assumed to be in measured in units of
%             metres squared (m^2).  Use function 'darcy' to convert from
%             (milli)darcies to m^2, e.g.,
%
%                 perm = convertFrom(perm, milli*darcy)
%
%             if the permeability is provided in units of millidarcies.
%
%             The field rock.perm may have ONE column for a scalar
%             permeability in each cell, TWO/THREE columns for a diagonal
%             permeability in each cell (in 2/3 D) and THREE/SIX columns
%             for a symmetric full tensor permeability.  In the latter
%             case, each cell gets the permeability tensor
%
%                 K_i = [ k1  k2 ]      in two space dimensions
%                       [ k2  k3 ]
%
%                 K_i = [ k1  k2  k3 ]  in three space dimensions
%                       [ k2  k4  k5 ]
%                       [ k3  k5  k6 ]
%
% OPTIONAL PARAMETERS
%   verbose   - Whether or not to emit informational messages throughout the
%               computational process.  Default value depending on the
%               settings of function 'mrstVerbose'.
%
%   facetrans -
%
%   invertBlocks -
%               Method by which to invert a sequence of small matrices that
%               arise in the discretisation.  String.  Must be one of
%                  - MATLAB -- Use an function implemented purely in MATLAB
%                              (the default).
%
%                  - MEX    -- Use two C-accelerated MEX functions to
%                              extract and invert, respectively, the blocks
%                              along the diagonal of a sparse matrix.  This
%                              method is often faster by a significant
%                              margin, but relies on being able to build
%                              the required MEX functions.
%
% RETURNS:
%   iB 
%   b          : defines the mpfa boundary faces
%
% COMMENTS:
%   PLEASE NOTE: Face normals have length equal to face areas.
%
% SEE ALSO:
%   `incompMPFA`, `mrstVerbose`.

%{
Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.

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


   opt = struct('verbose',      mrstVerbose,   ...
                'facetrans',    zeros([0, 2]), ...
                'invertBlocks', 'matlab',...
                'eta',0);

   opt = merge_options(opt, varargin{:});
   opt.invertBlocks = blockInverter(opt);

   if opt.verbose
      fprintf('Computing inner product on sub-half-faces ...\t');
      t0 = tic;
   else
      t0 = [];
   end

   [B, tbls] = robustComputeLocalFluxMimeticIP(G, rock, opt);
   
   tocif(opt.verbose, t0);
   
   if opt.verbose
      fprintf('Computing inverse mixed innerproduct ...\t');
      t0 = tic;
   end
   % The face-node degrees of freedom, as specified by the facenodetbl table, are
   % ordered by nodes first (see implementation below). It means in particular
   % that the matrix B is (or should be!) already block diagonal.
   facenodetbl = tbls.facenodetbl;
   [~, sz] = rlencode(facenodetbl.nodes); 
   iB   = opt.invertBlocks(B, sz);

   tocif(opt.verbose, t0);

   if opt.verbose
      fprintf('Computing multi-point transmissibilities ...\t');
      t0 = tic;
   end

   % Assembly of the divergence operator, from face-node faces to cell. We consider
   % only the internal faces for the moment (corresponds to Neumann boundary
   % condition).
   
   cellnodefacetbl = tbls.cellnodefacetbl;
   fno = cellnodefacetbl.faces; %alias
   cno = cellnodefacetbl.cells; %alias
   sgn = 2*(cno == G.faces.neighbors(fno, 1)) - 1;
   div = sparse(cellnodefacetbl.cells, ...
                (1 : cellnodefacetbl.num)', ...
                sgn, ...
                G.cells.num, ...
                cellnodefacetbl.num);
   % reduce from cell-face-node to face-node (equivalent to removing hybridization)
   op = setupTableMapping(facenodetbl, cellnodefacetbl, {'faces', 'nodes'});
   div = div*op;
   
   extfaces = (G.faces.neighbors(:, 1) == 0) |(G.faces.neighbors(:, 2) == 0);
   faceexttbl.faces = find(extfaces);
   [~, facenodeexttbl, bcmap, ~]= setupTableMapping(facenodetbl, faceexttbl, {'faces'});
   nbcdof = size(bcmap, 1);
   np = G.cells.num;
   nfacenode = facenodetbl.num;
   T = sparse(nfacenode, np + nbcdof);
   T(:, 1 : np) = iB*div';
   T(:, np + (1 : nbcdof)) = -iB*(bcmap');
   
   e_div = sparse(np + nbcdof, nfacenode);
   e_div(1 : np, :) = div;
   op = setupTableMapping(cellnodefacetbl, facenodeexttbl, {'faces', 'nodes'});
   e_sgn = op*sgn;
   tmpind = find(bcmap'*ones(nbcdof, 1));
   e_div(np + (1 : nbcdof), :) = sparse(1 : nbcdof, tmpind, e_sgn, nbcdof, nfacenode);
   
   tbls.faceexttbl = faceexttbl;
   tbls.facenodeexttbl = facenodeexttbl;
   
   mpfastruct = struct('T'    , T    , ...
                       'e_div', e_div, ...
                       'tbls' , tbls);
   
end

%--------------------------------------------------------------------------

function [B, tbls] = robustComputeLocalFluxMimeticIP(G, rock, opt)

    % Some short aliases 
    nc = G.cells.num;
    nf = G.faces.num;
    dim = G.griddim;
   
    cellfacetbl.cells = rldecode((1 : nc)', diff(G.cells.facePos)); 
    cellfacetbl.faces = G.cells.faces(:, 1);

    faces = rldecode((1 : nf)', diff(G.faces.nodePos)); 
    nodes = G.faces.nodes;
    % We setup the face-node table and it is ordered along ascending node numbers so
    % that we will have a block structure for the nodal scalar product.
    collatefacenode = [nodes, faces];
    collatefacenode = sortrows(collatefacenode, 1);
    facenodetbl.nodes = collatefacenode(:, 1);
    facenodetbl.faces = collatefacenode(:, 2);
    facenodetbl.num = numel(facenodetbl.faces);
    
    [op, cellnodefacetbl] = setupTableMapping(cellfacetbl, facenodetbl, {'faces'});

    % We setup the cell-face-node table, cellnodefacetbl. Each entry determine a
    % unique facet in a corner
    % We order cellfacenode in cell-node-face order. This is node to optimize
    % for-end loop below.
    orderingmat = [cellnodefacetbl.cells, cellnodefacetbl.nodes, ...
                   cellnodefacetbl.faces];
    orderingmat = sortrows(orderingmat);
    cellnodefacetbl.cells = orderingmat(:, 1);
    cellnodefacetbl.nodes = orderingmat(:, 2);
    cellnodefacetbl.faces = orderingmat(:, 3);
    
    % Somer shortcuts
    cno = cellnodefacetbl.cells;
    fno = cellnodefacetbl.faces;
    nno = cellnodefacetbl.nodes;
    
    % We setup the cell-node table, cellnodetbl. Each entry determine
    % a unique corner
    cfn = cellnodefacetbl;
    cellnode = [cfn.nodes, cfn.cells];
    cellnode = unique(cellnode, 'rows');
    cellnodetbl.nodes = cellnode(:, 1);
    cellnodetbl.cells = cellnode(:, 2);
    cellnodetbl.num   = numel(cellnodetbl.nodes);
    
    % Nodal scalar product is stored in vector nodeM
    % mattbl is the table which specifies how nodeM is stored: a matrix for
    % each "corner" (cell-node pair).
    duplicate = {'faces', {'faces1', 'faces2'}};
    [~, mattbl] = setupTableMapping(cellnodefacetbl, cellnodefacetbl, {'cells', 'nodes'}, ...
                                         'duplicate', {duplicate});
    % We order mattbl in cell-node-face1-face2 order
    % This is done to optimize for-end loop below
    orderingmat = [mattbl.cells, mattbl.nodes, ...
                   mattbl.faces1, mattbl.faces2];
    orderingmat = sortrows(orderingmat);
    mattbl.cells  = orderingmat(:, 1);
    mattbl.nodes  = orderingmat(:, 2);
    mattbl.faces1 = orderingmat(:, 3);    
    mattbl.faces2 = orderingmat(:, 4);    
    
    nodeM = zeros(mattbl.num, 1);
    
    numnodes = double(diff(G.faces.nodePos));
    numnodes = numnodes(fno);
    facetNormals = G.faces.normals(fno, :);
    facetNormals = bsxfun(@ldivide, numnodes, facetNormals);
    
    sgn = 2*(cno == G.faces.neighbors(fno, 1)) - 1;
    facetNormals = sgn.*facetNormals; % outward normals.
    
    % Assemble facePermNormals which corresponds to $Kn$ where n are the *outward*
    % normals at the facets.
    [perm, r, c] = permTensor(rock, G.griddim);
    Kn = bsxfun(@times, perm(cno, :), facetNormals(:, c));
    Kn = rldecode(Kn, dim*ones(numel(cno, 1)));
    coef = sparse(r, (1 : dim*dim)', ones(dim*dim, 1), dim, dim*dim);
    coef = repmat(coef, numel(cno), 1);
    Kn = coef.*Kn;
    Kn = sum(Kn, 2);
    Kn = reshape(Kn, dim, [])';
    
    facePermNormals = Kn;
    
    % Default option (opt.eta = 0): Use original face centroids and cell centroids,
    % NOT actual subface centroids. This corresponds to an MPFA method
    % (O-method) R = G.faces.centroids(fno,:) - G.cells.centroids(cno,:);
    cellFacetVec = G.faces.centroids(fno,:) - G.cells.centroids(cno,:) + ...
        opt.eta*(G.nodes.coords(nno,:) - G.faces.centroids(fno,:));
    
    % set up areas and volumes
    areas = G.faces.areas(fno);
    vols  = G.cells.volumes(cno);
    
    map = setupTableMapping(cellnodetbl, cellnodefacetbl, {'cells', 'nodes'}); 
    nfaces = diag(map'*map);
    
    cnf_i = 1; % start indice for the cellnodefacetbl index
    mat_i = 1; % start indice for the mattbl index
    
    for i = 1 : cellnodetbl.num
        
        nface = nfaces(i);
        cnfind = cnf_i : (cnf_i + (nface - 1));
        
        N     = facePermNormals(cnfind, :); 
        R     = cellFacetVec(cnfind, :);
        a     = areas(cnfind);
        v     = vols(cnfind);
        faces = cellnodefacetbl.faces(cnfind);
        
        cell = cellnodetbl.cells(i);
        node = cellnodetbl.nodes(i);
        
        K = reshape(perm(cell, :), [dim, dim]);
        
        % Assemble local nodal scalar product ( function node_ip2 below handle case when
        % N is invertible)
        locM = node_ip(a, v, ...
                       full(N), ...
                       full(R), ...
                       K);
        locM = reshape(locM, [], 1);
        
        matind = mat_i : (mat_i + (nface*nface - 1));
        nodeM(matind) = nodeM(matind) + locM;
        
        % increment start indices
        cnf_i = cnf_i + nface;
        mat_i = mat_i + nface*nface;        
        
    end

    sgn1 = 2*(mattbl.cells == G.faces.neighbors(mattbl.faces1, 1)) - 1;
    sgn2 = 2*(mattbl.cells == G.faces.neighbors(mattbl.faces2, 1)) - 1;
    nodeM = nodeM.*sgn1.*sgn2;   

    % Condensate on nodes (sum up cell contributions for give node).
    op = setupTableMapping(facenodetbl, facenodetbl, {'nodes'});
    [colind, rowind] = find(op);
    redmattbl.nodes  = facenodetbl.nodes(colind);
    redmattbl.faces1 = facenodetbl.faces(colind);
    redmattbl.faces2 = facenodetbl.faces(rowind);
    redmattbl.num    = numel(redmattbl.nodes);
    
    op = setupTableMapping(mattbl, redmattbl, {'nodes', 'faces1', 'faces2'});
    nodeM = op*nodeM;
    
    % Setup matrix
    % First set up facet indices in the redmattbl table
    mat1tbl.nodes  = facenodetbl.nodes;
    mat1tbl.faces1 = facenodetbl.faces;
    mat1tbl.ind    = (1 : facenodetbl.num)';
    op = setupTableMapping(redmattbl, mat1tbl, {'nodes', 'faces1'});
    [colind, rowind] = find(op);
    facetind1 = zeros(redmattbl.num, 1);
    facetind1(rowind) = mat1tbl.ind(colind);
    redmattbl.facetind1 = facetind1;
    
    mat2tbl.nodes  = facenodetbl.nodes;
    mat2tbl.faces2 = facenodetbl.faces;
    mat2tbl.ind    = (1 : facenodetbl.num)';
    op = setupTableMapping(redmattbl, mat2tbl, {'nodes', 'faces2'});
    [colind, rowind] = find(op);
    facetind2 = zeros(redmattbl.num, 1);
    facetind2(rowind) = mat2tbl.ind(colind);
    redmattbl.facetind2 = facetind2;
    
    % Assembly of B
    B = sparse(redmattbl.facetind1, ...
               redmattbl.facetind2, ...
               nodeM, ...
               facenodetbl.num, ...
               facenodetbl.num);
    
    tbls = struct('cellnodefacetbl', cellnodefacetbl, ...
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
    
    [U, D, V] = svd(N);
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
    t = 6 * sum(diag(K)) / size(K, 2);
    M = M + (t/v)*U*S*U';

end

function M = node_ip2(a, v, N, R, K)
    M = R*inv(N);
end

