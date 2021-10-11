function [T, T_noflow] = computeMultiPointTransLegacy(g, rock, varargin)
% Compute multi-point transmissibilities. This is the legacy version.
%
% SYNOPSIS:
%   T = computeMultiPointTransLegacy(G, rock)
%
% DESCRIPTION:
%
% We follow the local flux mimetic approach as described in this reference paper:
%
%      title     = {Local flux mimetic finite difference methods},
%      author    = {Lipnikov, Konstantin and Shashkov, Mikhail and Yotov, Ivan},
%      journal   = {Numerische Mathematik},
%      volume    = {112},
%      number    = {1},
%      pages     = {115--152},
%      year      = {2009},
%      publisher = {Springer}
% 
% The legacy version can only handle grid cells where corners have the same
% number of faces as the spatial dimension (this is always the case in 2D but
% not in 3D). 
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
%   T - half-transmissibilities for each local face of each grid cell
%       in the grid.  The number of half-transmissibilities equal the
%       number of rows in G.cells.faces.
%
% COMMENTS:
%   PLEASE NOTE: Face normals have length equal to face areas.
%
% SEE ALSO: `computeMultiPointTrans`, `private/computeMultiPointTransLegacy`
%   `incompMPFA`, `computeMultiPointTrans`, `computeMultiPointTransTensorAssembly`

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


% Written by Jostein R. Natvig, SINTEF ICT, 2009.

   opt = struct('verbose'     , mrstVerbose,   ...
                'facetrans'   , zeros([0, 2]), ...
                'invertBlocks', 'matlab',...
                'eta'         , 0);
   
   opt = merge_options(opt, varargin{:});
   opt.invertBlocks = blockInverter(opt);

   if opt.verbose
      fprintf('Computing mappings between for subfaces ...\t');
      t0 = tic;
   else
      t0 = [];
   end

   % Enumerate sub-faces and sub-half-faces.
   % Create mappings from {cells, nodes, faces} to {subfaces, subhalffaces}.

   [cno, nno, hfno, fno, subfno, subhfno] = createMapping(g);

   tocif(opt.verbose, t0);

   if opt.verbose
      fprintf('Computing inner product on sub-half-faces ...\t');
      t0 = tic;
   end

   [B,R] = computeLocalFluxMimeticIP(g, rock, cno, fno, nno, subhfno, opt);
   B = processFaceTrans(B, g, opt.facetrans(:,1), opt.facetrans(:,2), fno);

   tocif(opt.verbose, t0);

   if opt.verbose
      fprintf('Computing inverse mixed innerproduct ...\t');
      t0 = tic;
   end

   % Create matrices needed to compute transmissibilities
   s     = 2*(cno==g.faces.neighbors(fno,1))-1;
   D     = sparse(subhfno, subfno, 1); % Hybrid D matrix
   Do    = sparse(subhfno, subfno, s); % Mixed  D matrix
   C     = sparse(subhfno, cno, 1);

   % c1 adds up sub-half-face contributions for each half-face
   % c1     = sparse(subhfno, hfno, 1);

   % d1 adds up sub-face contributions for each face.
   counts = accumarray(fno, 1, [g.faces.num, 1]);
   d1     = sparse(subfno, fno, 1);

   % Note that c1'*D*d1 is equal to the reglar mimetic D matrix.

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
   % multiply fluxes with harmonic mean of mobility
   % this to avoid for re-asssembly
   % to be equivalent coupled reservoir simulation the value of
   % sum of upwind mobility should be used.
   % cf_trans_g=Do'*Do*iDoBDo*Do';
   cf_trans_g=iDoBDo*Do';
   T=struct('cf_trans'  , cf_mtrans , ... % transmisibility calculate K\grad on mpfa faces from cell pressures and boundary pressures
            'e_div'     , e_div     , ... % calculate div on cells and mpfa fluxes at boundary from mpfa fluxes
            'cf_trans_g', cf_trans_g, ... % calculate gravity contribution from gravity differences from mpfa half-faces
            'd1'        , d1        , ... % map from mpfa faces to faces
            'R'         , R         , ... % the continuity points for calculating gravity contributions
            'cno'       , cno       , ... % cno for mpfa faces
            'counts'    , counts    , ...
            'sb'        , sb          ... % defines the mpfa boundary faces
            );
   %
   if nargout > 1
       % the usefull trans for  other methods are
       %Trans =d1'*iDoBDo*Do';
       % resdused Trans for neumann baundary
       nc=size(C,2);
       iface=~sb;
       Trans=cf_mtrans;%iDoBDo*Do'*[C, -D(:,sb)];
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

function [cno, nno, hfno, fno, subfno, subhfno] = createMapping(g)
% Create mapping from sub - half - face to cell, node, face, half - face and
% sub - face
    cellno  = rldecode(1 : g.cells.num, diff(g.cells.facePos), 2) .'; 
    col     = 1 + (cellno == g.faces.neighbors(g.cells.faces(:, 1), 2)); 
    nhfaces = g.cells.facePos(end) - 1; 
    hfaces  = accumarray([g.cells.faces(:, 1), col], 1 : nhfaces); 
    hfaces  = rldecode(hfaces, diff(g.faces.nodePos)); 


    cells    = rldecode(g.faces.neighbors, diff(g.faces.nodePos)); 
    nodes    = repmat(g.faces.nodes, [2, 1]); 
    faces    = repmat(rldecode(1:g.faces.num, diff(g.faces.nodePos), 2)', [2, 1]); 
    subfaces = repmat((1:size(g.faces.nodes, 1))', [2, 1]); 
    i        = cells ~= 0; 
    w        = [cells(i), nodes(i), hfaces(i), faces(i), subfaces(i)]; 
    w        = double(sortrows(w)); 

    cno     = w(:, 1); 
    nno     = w(:, 2); 
    hfno    = w(:, 3); 
    fno     = w(:, 4); 
    subfno  = w(:, 5); 
    subhfno = (1:numel(cno))'; 
end

function [B, Rvec] = computeLocalFluxMimeticIP(g, rock, cno, fno, nno, subhfno, opt)
   [a, blocksz] = rlencode([cno,nno]);
   dims         = size(g.nodes.coords, 2);
   assert(all(blocksz==dims));
   % Expand centroid differences, normals, signs and permeabilities to
   % block-diagonal matrices such that B may be constructed by matrix-matrix
   % multiplications:
   [i, j] = ndgrid(subhfno, 1:dims);
   j      = bsxfun(@plus, j, rldecode(cumsum(blocksz)-blocksz(1), blocksz));

   % Use original face centroids and cell centroids, NOT actual subface
   % centroids.  This corresponds to an MPFA method (O-method)
   %R     = g.faces.centroids(fno,:) - g.cells.centroids(cno,:);
   Rvec = g.faces.centroids(fno,:) - g.cells.centroids(cno,:)+...
          opt.eta*(g.nodes.coords(nno,:)-g.faces.centroids(fno,:));
   R    = sparse(i,j,Rvec);

   % Subface sign == face sign
   sgn   = 2*(cno == g.faces.neighbors(fno,1)) -1;

   % Use fraction of face normal as subface normal
   numnodes = double(diff(g.faces.nodePos));
   N     = g.faces.normals  (fno,:).*sgn(:,ones(1, dims))./ ...
           numnodes(fno, ones(1,dims));
   N     = sparse(i,j,N);

   k     = permTensor(rock, dims);

   assert (size(k,1) == g.cells.num, ...
          ['Permeability must be defined in active cells only.\n', ...
           'Got %d tensors, expected %d (== number of cells).'],   ...
           size(k,1), g.cells.num);

   K     = sparse(i,j,reshape(k(a(:,1),:)', dims, [])');
   clear k a blocksz sgn

   % Construct B : Invert diagonal blocks of size sz
   d     = size(g.nodes.coords,2);       % == dims
   sz    = repmat(d, [size(R, 1)/d, 1]);
   B     = R * opt.invertBlocks(N*K, sz);
end

function B = processFaceTrans(B, g, f, t, fno)
% Modify B to account for face transmissibility (t) for regular faces (f)
% (a) Distribute 2*t to each half-face such that the harmonic mean
%     equals t.  (Use 1*t for boundary faces)
%
% (b) As an approximation, divide facetrans by number of corners in face
%     to distribute face transmissibility on each sub-face.
%
% (c) Modify B <-- B + diag(1/(2*t)).  This correspond to harmomic
%     mean when the mimetic innerproduct B is diagonal.
%

   % (a)
   invFaceTrans     = zeros(g.faces.num, 1);
   invFaceTrans (f) = t.*sum(g.faces.neighbors(f,:) ~= 0, 2);
   invFaceTrans (f) = 1./invFaceTrans (f);

   % (b)
   numNodes        = diff(g.faces.nodePos);
   invSubFaceTrans = invFaceTrans(fno).*double(numNodes(fno));

   % (c)
   B               = B + spdiags(invSubFaceTrans, 0, size(B, 1), size(B, 2));
end

%--------------------------------------------------------------------------

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
