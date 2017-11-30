function T = computeMultiPointTrans_DFM(g, rock, varargin)
%Compute multi-point transmissibilities.
%
% SYNOPSIS:
%   T = computeMultiPointTrans_DFM(G, rock)
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
%   hybrid -    If true, the interaction region is splitted for all nodes
%               shared by more then one hybrid cell. And the continuity
%               point is moved half an aperture in the normal direction.
%
% RETURNS:
%   T - half-transmissibilities for each local face of each grid cell
%       in the grid.  The number of half-transmissibilities equal the
%       number of rows in G.cells.faces.
%
% COMMENTS:
%   PLEASE NOTE: Face normals have length equal to face areas.
%
% SEE ALSO:
%   `incompMPFA`, `mrstVerbose`.

%{
Copyright 2009, 2010, 2011 SINTEF ICT, Applied Mathematics.

Portions Copyright 2011-2012 University of Bergen.

This file is part of DFM module of The MATLAB Reservoir Simulation Toolbox
(MRST).

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

   opt = struct('verbose',   mrstVerbose,   ...
                'facetrans', zeros([0, 2]),...
                'eta',0,'hybrid',false);
   opt = merge_options(opt, varargin{:});

   if opt.verbose,
      fprintf('Computing mappings between for subfaces ...\t');
      t0 = tic;
   else
      t0 = [];
   end

   %%
   % Enumerate sub-faces and sub-half-faces.
   % Create mappings from {cells, nodes, faces} to {subfaces, subhalffaces}.

   [g,cno, nno, hfno, fno, subfno, subhfno] = createMapping(g,opt);

   tocif(opt.verbose, t0);

   if opt.verbose,
      fprintf('Computing inner product on sub-half-faces ...\t');
      t0 = tic;
   end

   [B] = computeLocalFluxMimeticIP(g, rock, cno, fno, nno, subfno,subhfno,opt);
   B = processFaceTrans(B, g, opt.facetrans(:,1), opt.facetrans(:,2), fno,opt);

   tocif(opt.verbose, t0);

   if opt.verbose,
      fprintf('Computing inverse mixed innerproduct ...\t');
      t0 = tic;
   end

   %% Create matrices needed to compute transmissibilities
   s     = 2*(cno==g.faces.neighbors(fno,1))-1;
   D     = sparse(subhfno, subfno, 1); % Hybrid D matrix
   Do    = sparse(subhfno, subfno, s); % Mixed  D matrix
   C     = sparse(subhfno, cno, 1);

   % c1 adds up sub-half-face contributions for each half-face
   c1     = sparse(subhfno, hfno, 1);

   % d1 adds up sub-face contributions for each face.
   d1     = sparse(subfno, fno, 1);

   % Note that c1'*D*d1 is equal to the reglar mimetic D matrix.

   %% Construct the inverse of the mixed B
   % In the local-flux mimetic method, the mixed mass matrix DoBDo is block
   % diagonal with n_i x n_i blocks due to the special form of the hybrid
   % mass matrix.  Here n_i is the number of cells adjacent to a node in
   % the grid.
   DoBDo = Do'*B*Do;

   %% Invert DoBDo
   tmp      = unique([nno, subfno], 'rows');
   p        = tmp(:,2);
   P        = sparse(1:numel(p), p, 1, numel(p), numel(p));
   [sz, sz] = rlencode(tmp(:,1));
   iDoBDo   = P'*invertDiagonalBlocks(DoBDo(p, p), sz)*P;
   %iDoBDo= P'*DoBDo(p,p)*P;
   clear tmp sz P p

   tocif(opt.verbose, t0);

   if opt.verbose,
      fprintf('Computing multi-point transmissibilities ...\t');
      t0 = tic;
   end


   %% Compute multi-point transmissibilities
   % for all faces in terms of cell pressures and outer boundary pressures.
   b = full(sum(D, 1)) == 1;
   %T = c1'*Dm*inv(Dm'*B*Dm)*Dm'*[C, -D(:,b)*d1(b,:)];
    %if(nargout==2)
   Tg=c1'*Do*iDoBDo*Do';
    %end
   T = Tg*[C, -D(:,b)*d1(b,:)];
   %T(abs(T)>0 & abs(T)<max(max(T))*eps)=0;
   Tg=Tg*c1;
   tocif(opt.verbose, t0);
   % c1'*D*d1 har samme struktur som vanlig D.
   % T er feil størrelse å returnere dersom gravitasjon skal håndteres
   % skikkelig.  Gravitasjonsleddet kommer inn som c1'*Dm*iDmBDm*Dm'*f.
   % Siden gravitasjonsbidraget for all subfaces er likt kan de sikkert
   % skrives om til c1'*Dm*iDmBDm*Dm'*F*g der F*G=f.
   T=struct('T',T,'Tg',Tg);
   %T(:,all(T==0, 1))=[];
end

function [g,cno, nno, hfno, fno, subfno, subhfno] = createMapping(g,opt)
   %% Create mapping from sub-half-face to cell, node, face, half-face and
   %% sub-face

   % Recover the cell numbers corresponding to g.cell.faces (hereafter
   % cellfaces)
   cellno   = rldecode(1:g.cells.num, diff(g.cells.facePos), 2) .';

   % Which column in g.faces neighbors is the cell stored in
   col      = 1 + (cellno == g.faces.neighbors(g.cells.faces(:,1), 2));

   % Number of cell faces
   nhfaces  = g.cells.facePos(end)-1;

   % Pair the cellfaces to form connections (resembling g.faces.neighbors,
   % but using cellface index instead of facecells)
   hfaces   = accumarray([g.cells.faces(:,1), col], 1:nhfaces);

   % Clone the cellface connections into sub faces, one for each node on the
   % face
   hfaces   = rldecode(hfaces, diff(g.faces.nodePos));

   % Clone facecell connections
   cells    =  rldecode(g.faces.neighbors, diff(g.faces.nodePos));

   % Clone face indices into.
   % Also copy the information to account for both columns in cells
   faces    =  repmat(rldecode(1:g.faces.num, diff(g.faces.nodePos),2)', [2,1]);

   % Copy facenodes
   nodes    =  repmat(g.faces.nodes, [2,1]);

   % Create a numbering for subfaces, and copy that
   subfaces =  repmat((1:size(g.faces.nodes,1))', [2,1]);

   % Find facecell connections that are not on the boundary
   i        =  cells ~= 0;

   w        =  [cells(i), nodes(i), hfaces(i), faces(i), subfaces(i)];
   w        =  double(sortrows(w));

   % Split the interaction region consisting of hybrid cells.
   if opt.hybrid
       [g, w] = splittInteractionRegion(g,w);
   end

   cno     = w(:,1);
   nno     = w(:,2);
   hfno    = w(:,3);
   fno     = w(:,4);
   subfno  = w(:,5);
   subhfno = (1:numel(cno))';

end

function [B] = computeLocalFluxMimeticIP(g, rock, cno, fno, nno, subfno,subhfno,opt)
[a, blocksz] = rlencode([cno,nno]);
dims         = size(g.nodes.coords, 2);
assert(all(blocksz==dims));

%% Expand centroid differences, normals, signs and permeabilities to
%  block-diagonal matrices such that B may be constructed by matrix-matrix
%  multiplications:
[i,j] = ndgrid(subhfno, 1:dims);
j     = bsxfun(@plus, j, rldecode(cumsum(blocksz)-blocksz(1), blocksz));

% Use original face centroids and cell centroids, NOT actual subface
% centroids.  This corresponds to an MPFA method (O-method)

% TODO: generalize for eta [0-1]
assert(opt.eta == 0,'only vallid for eta = 0')
eta=0*ones(length(fno),1);


iP = g.faces.centroids(fno,:)+bsxfun(@times,eta,(g.nodes.coords(nno,:)-g.faces.centroids(fno,:)));
R  = iP - g.cells.centroids(cno,:);

% Subface sign == face sign
sgn   = 2*(cno == g.faces.neighbors(fno,1)) -1;

% Use fraction of face normal as subface normal
numnodes = double(diff(g.faces.nodePos));

% Hybrid faces consist of one node in the grid, but two node in the interaction region
% i.e multiply by 2.
if(opt.hybrid)
    numnodes(g.faces.hybrid==1) = numnodes(g.faces.hybrid==1)*2;
    %numnodes(fno(isHybHybEnds)) = numnodes(fno(isHybHybEnds))-1;
end

N     = g.faces.normals  (fno,:).*sgn(:,ones(1, dims))./ ...
    numnodes(fno, ones(1,dims));

% For the matrix-fracture transmissibility the interaction region is
% expanded with half an aperture in the normal direction outward from the
% fracture
if(opt.hybrid)

    % Find hybrid faces
    isHybface=g.faces.hybrid(fno)==1;

    hybrids=g.cells.hybrid(cno) & ~isHybface;

    [~,hereMap]=unique(subfno,'first');

    next2Hybrids=hereMap(subfno(hybrids));

    % the unit normal vector of the fracture faces
    n = bsxfun(@times,g.faces.normals(fno(hybrids),:).*sgn(hybrids,ones(1, dims)),1./sqrt(sum(g.faces.normals(fno(hybrids),:).*g.faces.normals(fno(hybrids),:),2)));

    % the aperture
    apertures = g.cells.volumes(cno(hybrids))./g.faces.areas(fno(hybrids));

    % expand the interaction region
    expantion = bsxfun(@times,n,apertures./2);
    R(hybrids,:) = R(hybrids,:)+expantion;
    R(next2Hybrids,:) = R(next2Hybrids,:)+expantion;
end

R     = sparse(i,j,R);
N     = sparse(i,j,N);

k     = permTensor(rock, dims);
K     = sparse(i,j,reshape(k(a(:,1),:)', dims, [])');
clear k a blocksz sgn

%% Construct B : Invert diagonal blocks of size sz
d     = size(g.nodes.coords,2);       % == dims
sz    = repmat(d, [size(R, 1)/d, 1]);
B     = R*invertDiagonalBlocks(N*K, sz);

end

function B = processFaceTrans(B, g, f, t, fno,opt)
%% Modify B to account for face transmissibility (t) for regular faces (f)
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

   % Hybrid faces consist of one node in the grid, but two node in the interaction region
   % i.e multiply by 2.
   if(opt.hybrid)
       numNodes(g.faces.hybrid==1)=numNodes(g.faces.hybrid==1)*2;
       % numNodes(fno(isHybHybEnds))=numNodes(fno(isHybHybEnds))-1;
   end

   invSubFaceTrans = invFaceTrans(fno).*double(numNodes(fno));

   % (c)
   B               = B + spdiags(invSubFaceTrans, 0, size(B, 1), size(B, 2));
end


% %%{
% %% Matlab code calling mex function invertSmallMatrices.c
% function iA = invertDiagonalBlocks(A, sz)
%    sz = int32(sz);
%    blocks = matrixBlocksFromSparse(A, sz);
%
%    [i, j] = blockDiagIndex(sz, sz);
%    iA  = sparse(i,j, invv(blocks, sz));
% end
% %}

%
%% Pure Matlab code using inv
function iA = invertDiagonalBlocks(A, varargin)
   [a,b,c,d]=dmperm(A+speye(size(A)));
   if all(diff(c)==1),
      % For diagonal A,
      iA = inv(A);
   else
      % For block-diagonal A, we identify each diagonal block and invert
      % using full matrix version of inv.  This is more efficient than a
      % sparse inv call.
      sz = sum(diff(c).^2);
      I = zeros(sz, 1);
      J = I;
      V = I;
      pos = 0;
      for k=1:numel(c)-1,
         i     = a( c(k) : c(k+1)-1 );
         j     = b( d(k) : d(k+1)-1 );
         ip    = inv(full(A(i,j)));
         sz    = numel(i)^2;

         I(pos+1:pos+sz)  = repmat(i, [numel(j), 1])';
         J(pos+1:pos+sz)  = repmat(j, [numel(i), 1]);
         V(pos+1:pos+sz)  = ip;

         pos   = pos + sz;
      end
      clear A
      iA = sparse(I,J,V);
   end
end
%


function [g,new_w] = splittInteractionRegion(g,w)
% Splits the interaction region when more then one hybrid cell
% meets at a vertex.

% Faces with tag of at least this value are considered hybrid. Note that
% this puts a constraint on how the hybrid cells are identified
MIN_HYBRID_TAG_VAL = 1;

cno    = w(:,1);
nno    = w(:,2);
hfno   = w(:,3);
fno    = w(:,4);
subfno = w(:,5);

new_w  = w;
clear w;

% the hybrid cells
isHybridCell = g.cells.hybrid(cno) == 1;

% store face tags
faceTags = double(g.faces.tags(fno));

% subfno contains one or two indices, depending on whether the subface has
% a neighbor or not. Obtain a mapping between the two neighbors (if they
% exists)
[~, hereMap] = unique(subfno,'first');
[~,thereMap] = unique(subfno,'last');

% Find all faces that are both on a fracture and belong to a hybrid cell
fracFace = find(faceTags >= MIN_HYBRID_TAG_VAL & isHybridCell);

% Pick one fracture interface for each interaction region,
[~,I] = unique([nno,faceTags >= MIN_HYBRID_TAG_VAL],'rows');
here = I(faceTags(I) >= MIN_HYBRID_TAG_VAL);

% Nodes located on a fracture have multiple interaction regions, since the
% faces are split when computing transmissibilities. oneStep finds one such
% region by traversing faces surrounding this vertex, and adding new faces.
% Note that since this is the first region around this node, no new node is
% created.
[~,subhintf_h,fracIntf_red] = oneStep(subfno,nno,cno,fno,g,here,fracFace,hereMap,thereMap);

% Hybrid faces (between two hybrid cells) belongs to both (multiple)
% regions, even though they only have one node. Thus the face has not been
% split yet, but now is the time. Prepare for copying.
isHybrid = g.faces.hybrid(fno) == 1;
subhintf_h_copy = setxor(find(isHybrid),subhintf_h);

% number of nodes
new_N_num = g.nodes.num;

% number of subfaces
subfno_num = max(subfno);

% keep track of the number of steps, we dont support more then 100 regions.
% If there is still more fracture faces left something is wrong
count = 0;

% Loop over all fracture faces that have no interaction region
% TODO: 100 iterations have always been sufficient, but a more nuanced
% stopping criterion in case of bugs may be useful
while (~isempty(fracIntf_red) && count<100)

    % Pick next fracture face from the list of those that have not yet been
    % assigned an interaction region
    [~,I] = unique(nno(fracIntf_red));
    here = fracIntf_red(I);

    % Find the region for this fracture face
    [thisRegion,subhintf_h,fracIntf_red] = oneStep(subfno,nno,cno,fno,g,here,fracIntf_red,hereMap,thereMap);

    % The new region must be assigned a unique node number
    [old_N,I,J] = unique(nno(thisRegion));
    new_N = (new_N_num+1:new_N_num+length(I))';
    new_N_num = new_N_num+length(I);
    new_w(thisRegion,2) = new_N(J);
    Map = accumarray(old_N,new_N,[g.nodes.num,1]);
    g.nodes.coords = [g.nodes.coords;g.nodes.coords(old_N,:)];
    g.nodes.num = new_N_num;

    % copy the corresponing hybrid faces.
    intc = intersect(subhintf_h,subhintf_h_copy);
    new_w(intc,2) = Map(nno(intc));
    subhintf_h = setdiff(subhintf_h,intc);
    subhintf_h_copy = setdiff(subhintf_h_copy,intc);
    new_nno = Map(nno(subhintf_h));
    new_cno = cno(subhintf_h);
    new_hfno = hfno(subhintf_h);
    new_fno = fno(subhintf_h);
    new_subfno = (subfno_num+1:subfno_num+length(subhintf_h))';
    subfno_num = subfno_num+length(subhintf_h);
    new_w = [new_w;[new_cno,new_nno,new_hfno,new_fno,new_subfno]];
    count = count+1;
end

% Error if the loop reached the upper limit of iterations
if ~isempty(here) && count == 100
    error('Error when splitting hybrid interaction regions for MPFA')
end

% return a sorted list.
new_w = sortrows(new_w);

end


% Look for cells in the next region.
function [thisRegion,subhintf_h,fracIntf_red] = ...
                    oneStep(subfno,nno,cno,fno,g,here,fracIntf,hereMap,thereMap)

% keep track of this interaction region
thisRegion = [];

% Assign those hybrid faces that belongs to this region
subhintf_h = getHybridSubIntf(here,nno,cno,fno,g);

% Find the index of the sub face as seen from the other side of the
% connection. In this case, we move from the hybrid cell to 'normal' cell
% on the other side of the fracture face.
there = setxor(union(hereMap(subfno(here)),thereMap(subfno(here))),here);

% remove here from the fracture interface list; we know that this is the
% fracture
fracIntf_red = setdiff(fracIntf,here);

% update this interaction region
thisRegion = [thisRegion; here; there];

% walk from fracture interface to fracture interfaces as
% long as there are interfaces left, we don't want to walk more then
% 100 steps. If we are not back after 100 steps something is probably wrong.
count = 0;
while (~isempty(there) && count < 100)

    % Find the local sub index of all faces of the cell of there
    % (cno(there)), which also share the vertex with there.
    i = find(ismember([nno,cno],[nno(there),cno(there)],'rows'));

    % We have already processed the 'there' face, remove this.
    % Here now represents the next faces that we will include in the
    % interaction region.
    here = setdiff(i,there);

    % Find the face neighbors of here.
    there = setxor(union(hereMap(subfno(here)),thereMap(subfno(here))),here);

    % Find those there interfaces that are at a fracture.
    frac_intf_there = intersect(fracIntf_red,there);

    % Add new hybrid interfaces associated with frac_intf_there
    subhintf_h = [subhintf_h ; getHybridSubIntf(frac_intf_there,nno,cno,fno,g)];

    % update thisRegion with the here and theres of this step
    thisRegion = [thisRegion;here;there];

    % There faces that are on fractures are removed from the global list of
    % fractures without interaction regions
    fracIntf_red = setdiff(fracIntf_red,frac_intf_there);

    % Remove fracture faces from the list of seeds used in the next
    % iteration
    there = setdiff(there,frac_intf_there);

    count = count+1;
end

% Error if the loop reached the upper limit of iterations
if ~isempty(here) && count == 100
    error('Error when splitting hybrid interaction regions for MPFA')
end

end


function subhintf_h = getHybridSubIntf(here,nno,cno,fno,g)
% finds the hybrid subfaces associated with the region

% find the local sub index of the subface sharing node and cell with here
i = find(ismember([nno,cno],[nno(here),cno(here)],'rows'));

% pick out the new subhintf i.e not here,
subhintf_h = setxor(i,here);

% we are interested in the hybrid faces only
subhintf_h = subhintf_h(g.faces.hybrid(fno(setxor(i,here)))==1);
end


