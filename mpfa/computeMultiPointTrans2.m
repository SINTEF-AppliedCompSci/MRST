function [T, T_noflow] = computeMultiPointTrans2(g, rock, varargin)
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
%   T - structure with fields
%       'cf_trans'    : transmisibility calculate K\grad on mpfa faces from cell pressures and boundary pressures
%       'e_div'       : calculate div on cells and mpfa fluxes at boundary from mpfa fluxes
%       'cf_trans_g'  : calculate gravity contribution form gravity
%                       differences from mpfa half faces (not implemented yet)
%       'd1'          : map from mpfa faces to faces
%       'R'           : the continuity points for for calculating gravity
%                       contributions (not implemented yet).
%       'cno'         : cell numbers for mpfa faces
%       'counts'      :
%       'sb'          : defines the mpfa boundary faces
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


% Written by Jostein R. Natvig, SINTEF ICT, 2009.

   opt = struct('verbose',      mrstVerbose,   ...
                'facetrans',    zeros([0, 2]), ...
                'invertBlocks', 'matlab',...
                'eta',0);

   opt = merge_options(opt, varargin{:});
   opt.invertBlocks = blockInverter(opt);


   if opt.verbose
      fprintf('Computing inner product on sub-half-faces ...\t');
      t0 = tic;
   end

   [B, tables] = robustComputeLocalFluxMimeticIP(g, rock, opt);

   if opt.verbose
      fprintf('Computing inverse mixed innerproduct ...\t');
      t0 = tic;
   end

   % Create matrices needed to compute transmissibilities
   s  = 2*(cno == g.faces.neighbors(fno,1)) - 1;
   D  = sparse(subhfno, subfno, 1); % Hybrid D matrix
   Do = sparse(subhfno, subfno, s); % Mixed D matrix
   C  = sparse(subhfno, cno, 1);

   % c1 adds up sub-half-face contributions for each half-face
   % c1 = sparse(subhfno, hfno, 1);

   % d1 adds up sub-face contributions for each face.
   counts = accumarray(fno, 1, [g.faces.num, 1]);
   d1     = sparse(subfno, fno, 1);

   % Note that c1'*D*d1 is equal to the regular mimetic D matrix.

   % Construct the inverse of the mixed B
   % In the local-flux mimetic method, the mixed mass matrix DoBDo is block
   % diagonal with n_i x n_i blocks due to the special form of the hybrid
   % mass matrix.  Here n_i is the number of cells adjacent to a node in
   % the grid.
   DoBDo = Do'*B*Do;

   % Invert DoBDo
   tmp      = unique([nno, subfno], 'rows');
   p        = tmp(:,2);
   P        = sparse(1:numel(p), p, 1, numel(p), numel(p));
   [sz, sz] = rlencode(tmp(:,1));
   iDoBDo   = P' * opt.invertBlocks(DoBDo(p, p), sz) * P;
   clear tmp sz P p

   tocif(opt.verbose, t0);

   if opt.verbose
      fprintf('Computing multi-point transmissibilities ...\t');
      t0 = tic;
   end


   % Compute multi-point transmissibilities
   % for all faces in terms of cell pressures and outer boundary pressures.
  
   %T = c1'*Dm*inv(Dm'*B*Dm)*Dm'*[C, -D(:,b)*d1(b,:)];
    %if(nargout==2)
  
   tocif(opt.verbose, t0);
   % c1'*D*d1 har samme struktur som vanlig D.
   % T er feil størrelse å returnere dersom gravitasjon skal håndteres
   % skikkelig.  Gravitasjonsleddet kommer inn som c1'*Dm*iDmBDm*Dm'*f.
   % Siden gravitasjonsbidraget for all subfaces er likt kan de sikkert
   % skrives om til c1'*Dm*iDmBDm*Dm'*F*g der F*G=f.
   %T=struct('T',T,'Tg',Tg,'hfhf',Do*iDoBDo*Do','c1',c1,'D',D,'d1',d1,'C',C,'Do',Do,'R',R,'cno',cno);
   %{
    %old structure
    b = full(sum(D, 1)) == 1;
    Tg=c1'*Do*iDoBDo*Do';
    %end
    T = Tg*[C, -D(:,b)*d1(b,:)];
    Tg=Tg*c1;
    T=struct('T',T,'Tg',Tg,'hfhf',Do*iDoBDo*Do','c1',c1,'D',D,'d1',d1,'C',C,'Do',Do,'R',R,'cno',cno);
   %}
   sb = full(sum(D, 1)) == 1;
   %cf_mtrans=Do'*Do*iDoBDo*Do'*[C, -D(:,sb)];
   cf_mtrans=iDoBDo*Do'*[C, -D(:,sb)];
   % define div operaor
   e_div =  [C, -D(:,sb)]'*Do;
   % multiply fluxes with harmonic mean of mobility this to avoid for re asssembly
   % to be equivalent coupled reservoir simulation the value of sum of upwind
   % mobility should be used.
   % cf_trans_g = Do'*Do*iDoBDo*Do';
   cf_trans_g = iDoBDo*Do';
   T = struct('cf_trans'  , cf_mtrans ,... % transmisibility calculate K\grad on mpfa faces from cell pressures and boundary pressures
              'e_div'     , e_div     ,... % calulate div on cells and mpfa fluxes at boundary from mpfa fluxes
              'cf_trans_g', cf_trans_g,... % calulate gravity contribution form gravity diferences from mpfa half faces
              'd1'        , d1        ,... % map from mpfa faces to faces
              'R'         , R         ,... % the continuity points fo for calculating gravity contricutions
              'cno'       , cno       ,... % cno for mpfa faces
              'counts'    , counts    , ...
              'sb'        , sb ... % defines the mpfa boundary faces
            );
   %
   if nargout > 1
       % the usefull trans for  other methods are
       % Trans =d1'*iDoBDo*Do';
       % reduced Trans for neumann baundary
       nc=size(C,2);
       iface=~sb;
       Trans=cf_mtrans; % iDoBDo*Do'*[C, -D(:,sb)];
       A=Trans(iface,1:nc);
       B=Trans(iface,nc+1:end);
       C=Trans(~iface,1:nc);
       D=Trans(~iface,nc+1:end);
       rTrans = A-B*(D\C);
       % reduce to internal normal

       %
       intfaces=~any(g.faces.neighbors==0,2);
       rTrans = d1(iface,intfaces)'*rTrans; % mpfa trans from cell pressure to internal fluxes
       N=g.faces.neighbors(intfaces,:);
   
       % Create transmissibility as a hidden, undocumented output
       % gravity contributaion 
       gTrans = iDoBDo*d1; % note that this maps from gravity differences over faces including outer faces.
       % gravity trans ???
       rgTrans = gTrans(iface,:) + B*(D\gTrans(~iface,:));
       rgTrans = d1(iface,intfaces)'*rgTrans;   
       T_noflow=struct('rTrans',rTrans,...%calculate K\grad from cell pressures assuming no flow boundary 
       'rgTrans',rgTrans','N',N);%calculate mpfa gravity contribution from "gravity difference between cells and cells to bounary faces" to internal face flux 
   end
end

%--------------------------------------------------------------------------

function [B, tbls] = robustComputeLocalFluxMimeticIP(G, rock, opt)

    nc = G.cells.num;
    cellfacetbl.cells = rldecode((1 : nc)', diff(G.cells.facePos)); 
    cellfacetbl.faces = G.cells.faces(:, 1);

    nf = G.faces.num;
    faces = rldecode((1 : nf)', diff(G.faces.nodePos)); 
    nodes = G.faces.nodes;
    % we setup the face-node table and it is ordered along ascending node numbers so
    % that we will have a block structure for the nodal scalar product.
    collatefacenode = [nodes, faces];
    collatefacenode = sortrow(collatefacenode, 1);
    facenodetbl.nodes = collatefacenode(:, 1);
    facenodetbl.faces = collatefacenode(:, 2);
    facenodetbl.num = numel(facenodetbl.faces);
    
    [op, colind, rowind] = setupTableMapping(cellfacetbl, facenodetbl, ...
                                                          'faces');

    cellfacenodetbl.cells = cellfacetbl.cells(rowind);
    cellfacenodetbl.faces = cellfacetbl.faces(rowind);
    cellfacenodetbl.nodes = facenodetbl.nodes(colind);
    cellfacenodetbl.num = numel(cellfacenodetbl.cells);

    nn = G.nodes.num;
    cfn = cellfacenodetbl;
    op = sparse(cfn.cells, cfn.nodes, 1, nc, nn);
    [cind, nind] = find(op);
    cellnodetbl.cells = cind;
    cellnodetbl.nodes = nind;

    % nodal scalar product is stored in vector nodeM
    % mattbl is the table which specifies how nodeM is stored: a matrix for
    % each "corner" (cell-node pair).
    [~, colind, rowind] = setupTableMapping(cellfacenodetbl, cellfacenodetbl, ...
                                                          'cells', 'nodes'); 
    mattbl.cells  = cellfacenodetbl.cells(colind);
    mattbl.nodes  = cellfacenodetbl.nodes(colind);
    mattbl.faces1 = cellfacenodetbl.faces(colind);
    mattbl.faces2 = cellfacenodetbl.faces(rowind);
    mattbl.num = numel(mattbl.cells);
    nodeM = zeros(mattbl.num, 1);
    
    dim = G.griddim;

    faceNormals = bsxfun(@divide, G.faces.normals, G.faces.areas);
    fno = cellfacenodetbl.faces;
    facetNormals = faceNormals(fno, :);
    
    % Assemble facePermNormals which corresponds to $Kn$ where n are the
    % normals at the facets.
    [K, r, c] = permTensor(rock, G.griddim);
    error('fix the sign of the normals');
    cno = cellfacenodetbl.cells;
    Kn = bsxfun(@times, K(cno, :), facetNormal(:, c));
    Kn = rldecode(Kn, dim*ones(numel(cno, 1)));
    coef = sparse(r, (1 : dim*dim)', ones(dim*dim, 1), dim, dim*dim);
    coef = repmat(coef, numel(cno), 1);
    Kn = coef.*Kn;
    Kn = sum(Kn, 2);
    Kn = reshape(Kn, dim, [])';
    
    facePermNormals = Kn;
    
    % Use original face centroids and cell centroids, NOT actual subface
    % centroids.  This corresponds to an MPFA method (O-method)
    % R = G.faces.centroids(fno,:) - G.cells.centroids(cno,:);
    cellFacetVec = G.faces.centroids(fno,:) - G.cells.centroids(cno,:) + ...
        opt.eta*(G.nodes.coords(nno,:) - G.faces.centroids(fno,:));
    
    % set up areas and volumes
    areas = G.faces.areas(fno);
    vols  = G.cells.volumes(cfo);
    
    op = setupTableMapping(cellfacenodetbl, cellnodetbl, 'cells', 'nodes'); 

    for i = 1 : cellnodes.num
        
        logfacets = logical(op(i, :)');
        
        N = facePermNormals(logfacets, :); 
        R = cellFacetVec(logfacets, :);
        a = areas(logfacets);
        v = vols(logfacets);
        
        cell = cellnode.cells(i);
        node = cellnode.nodes(i);
        faces  = cellfacenodetbl.faces(logfacets);
        
        % Assemble local nodal scalar product
        locM = node_ip(a, v, N, R);
        locM = reshape(M, [], 1);
    
        nfaces = nnz(logfacets);
        loctbl.faces1 = repmat(faces, nfaces, 1);
        loctbl.faces2 = rldecode(faces, nfaces*ones(nfaces, 1));
        loctblNum = numel(loctbl.faces1);
        loctbl.cells = cell*ones(loctblNum, 1);
        loctbl.nodes = node*ones(loctblNum, 1);
        
        [~, colind, ~] = setupTableMapping(loctbl, mattbl, 'cells', ...
                                                   'nodes', 'faces1', ...
                                                   'faces2');
        nodeM(colind) = locM;
    
    end
    
    
    [~, colind, rowind] = setupTableMapping(facenodetbl, facenodetbl, 'nodes');
    redmattbl.nodes  = facenodetbl.nodes(colind);
    redmattbl.faces1 = facenodetbl.faces(colind);
    redmattbl.faces2 = facenodetbl.faces(rowind);

    op = setupTableMapping(redmattbl, mattbl, 'nodes', 'faces1', 'faces2');
    B = op'*nodeM*op;

    tbls = struct('cellfacenodetbl', cellfacenodetbl, ...
                  'cellnodetbl'    , cellnodetbl    , ...
                  'facenodetbl'    , facenodetbl)
    
end

function M = node_ip(a, v, N, R)
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
    assert(prod(diagDp)>0, 'cannot assemble mpfa, need extra fix'); 
    invd = 1./d;
    H = [diag(invd), zeros(fnum - dim)];
    M = R*V*H*U';
    
    % Add stabilization term S
    S = blkdiag(zeros(dim), eyes(fnum - dim));
    U = diag(a)*U;
    t = 6 * sum(diag(K)) / size(K, 2);
    M = M + (t/v)*U*S*U';
    
end

function bi = blockInverter(opt)
   if ~ischar(opt.invertBlocks)
      dispif(opt.verbose, ...
            ['Unsupported option value of type ''%s'' in ', ...
             'option ''invertBlocks''. Reset to default ' , ...
             '(''matlab'')\n'], class(opt.invertBlocks));
      opt.invertBlocks = 'matlab';
   end

   switch lower(opt.invertBlocks)
      case {'matlab', 'm', 'builtin'}
         bi = @invertDiagonalBlocks;
      case {'mex', 'c', 'accelerated'}
         bi = @invertDiagonalBlocksMex;
      otherwise
         dispif(opt.verbose, ...
               ['Unsupported value ''%s'' in option ', ...
                '''invertBlocks''.\nMust be one of ''matlab'' or ', ...
                '''mex''.\nReset to default (''matlab'').'], ...
                opt.invertBlocks);

         bi = @invertDiagonalBlocks;
   end
end

%--------------------------------------------------------------------------

% Matlab code calling mex function invertSmallMatrices.c
function iA = invertDiagonalBlocksMex(A, sz)
   sz     = int32(sz);
   blocks = matrixBlocksFromSparse(A, sz);
   iA     = blockDiagMatrix(invv(blocks, sz), sz);
end

%--------------------------------------------------------------------------

% Pure Matlab code using inv
function iA = invertDiagonalBlocks(A, sz)
   V = zeros([sum(sz .^ 2), 1]);
   [p1, p2] = deal(0);

   for b = 1 : numel(sz)
      n  = sz(b);
      n2 = n * n;
      i  = p1 + (1 : n);

      V(p2 + (1 : n2)) = inv(full(A(i, i)));

      p1 = p1 + n;
      p2 = p2 + n2;
   end

   iA = blockDiagMatrix(V, sz);
end

%--------------------------------------------------------------------------

function A = blockDiagMatrix(V, sz)
   [I, J] = blockDiagIndex(sz);
   A      = sparse(I, J, V);
end
