function out = mpfa(G,rock,activeNodes,varargin)
% Discretize Darcy's law by a multipoint flux approximation
%
% Parameters:
% G - mrst grid structure, see http://www.sintef.no/projectweb/mrst/
%     To get an overview of the data format, type 'help grid_structure'
% rock - mrst rock structure. Should contain a field perm. To get an
%    overview of the data format, type 'help makeRock'
% Active nodes: Which nodes should the stencil be computed for. If empty,
%   all nodes in the grid is considered. 
%
%   NOTE: The function computes and stores in memory the inverse of a 
%   block-diagonal matrix
%   (which expresses gradients on the sub-cells as a function of pressures
%   in the surrounding cells). The size of each sub-matrix will be the
%   number of sub-cells sharing the vertex * number of dimensions, so each
%   system will be 8x8 for Cartesian 2D, 24x24 for Cartesian 3D. For large
%   3D grids, in particular with simplices, this may be heavy in terms of
%   memory needs. For this reason, the activeNodes option may be an atractive 
%   option. For a work-around that splits the discretization into several
%   calls to mpfa, and thus reduces memory consumption, see mpfa_subgrid.
%
% Optional parametrs, defined as keyword - value pairs:
% 'eta' - Location of continuity point on the half edges, defined
%   according to Aavatsmark 2002. Between 0 and 1. Default value is 0, for
%   simplices, the value should be 1/3 (will give symmetric method, see
%   Klausen, Radu Eigestad 2008).
% bc - boundary conditions, as defined by the MRST function addBC. If none
%    are provided, homogeneous Neumann conditions will be assigned.
% invertBlocks - method for inverting block diagonal systems. Should be
%    either 'matlab' or 'mex'. The former is pure matlab, which will be
%    slow for large problems. Mex is substantially faster, but requires
%    that matlab has access to a mex / c compiler.
% returnHalfEdgeValues - whether the return values for stresses etc should
%   be for sub edges/faces (or summed to the whole face). Default is the
%   full faces. Can be used for a partial update of the discretization
%   scheme.
%
%{
Copyright 2015-2016, University of Bergen.

This file is part of FVBiot.

FVBiot is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

FVBiot is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file.  If not, see <http://www.gnu.org/licenses/>.
%} 


opt=struct('eta', 0,...
    'bc',[],...
    'invertBlocks','mex', ...
    'returnHalfEdgeValues',0);

if nargin < 3 || numel(activeNodes) == 0
    activeNodes = 1:G.nodes.num;
end

opt=merge_options(opt,varargin{:});
opt.invertBlocks = blockInverter(opt);

% Bookkeeping
Nd = G.griddim;
Nc = G.cells.num;

% Number of nodes per face
nFaceNodes = diff(G.faces.nodePos);

% Boundary conditions
[isNeumann,isDirichlet] = classifyBoundaryFaces(G,opt.bc);

% Subcell topology; give numbers to sub-cells, sub-faces etc.
[cnoAll, nnoAll, ~, fnoAll, subfnoAll, subhfnoAll] = createSubcellMapping(G);

% If not all nodes are active, we need to differ between active and passive
% sub-faces, sub-cells etc.
localTopology = ismember(nnoAll, activeNodes);
cnoGlob = cnoAll(localTopology);
fnoGlob = fnoAll(localTopology);
subfnoGlob = subfnoAll(localTopology);
[~,~,nno] = unique(nnoAll(localTopology));
[tmpcno,~,cno] = unique(cnoAll(localTopology));
[~,~,fno] = unique(fnoAll(localTopology));
[~,~,subfno] = unique(subfnoAll(localTopology));
[~,~,subhfno] = unique(subhfnoAll(localTopology));

% Obtain permeability tensor
K = permTensor(rock,Nd); 
K = K(tmpcno,:);

%% Discretization

% Contribution of gradients to flux continuity (-nK)
[nK,cellNodeBlocks,subcind] = multiplyNormalVectorAndTensor(G,K,cno,fnoGlob,nno,subhfno);

% Distance from cell centers to face continuity points
pCont = computeDistFaceCell(G,cnoGlob,fnoGlob,nno,subhfno,opt.eta);

% Darcy's law, only computed from one side of the face
[~,uniqueSubfno] = unique(subfno);
darcy = -nK(uniqueSubfno,:);

hf2f = sparse(fno(uniqueSubfno),subfno(uniqueSubfno),1);
sgn   = 2*(cnoGlob == G.faces.neighbors(fnoGlob,1)) -1;
% nK and pCont was computed from both sides of a subface. Now combine them
pairOverSubfaces = sparse(subfno,subhfno,sgn);
nK = pairOverSubfaces * nK; pCont = pairOverSubfaces * pCont;

clear pairOverSubfaces

% Right hand side of the equations (corresponding to cell center
% contributions)
pContCC = sparse(subfno,cno,sgn);
nKCC = sparse(max(subfno),max(cno)); 

% Done with matching expression from the two sides of the faces, now reduce
% the fields to one per sub-face (e.g. from one side of the faces, but not
% the other).
nno = nno(uniqueSubfno);
subfno = subfno(uniqueSubfno);
% Also find the global index of unique faces
fnoGlob = fnoGlob(uniqueSubfno);
subfnoGlob = subfnoGlob(uniqueSubfno);

nsubfno = max(subfno);
sgn   = 2*(cnoGlob(uniqueSubfno) == G.faces.neighbors(fnoGlob,1)) -1;

%% Exclude equations from boundary faces
j = find(~isNeumann(fnoGlob));
i = 1:numel(j);
excludeNeumann = sparse(i,j,1,numel(i),nsubfno);

j = find(~isDirichlet(fnoGlob));
i = 1 : numel(j);
excludeDirichlet = sparse(i,j,1,numel(i),nsubfno);
nK = excludeDirichlet * nK;
pCont = excludeNeumann * pCont;

clear i j 

% Exclude boundary faces from right hand side
nKCC = excludeDirichlet * nKCC;
pContCC = excludeNeumann * pContCC;

%% Write the equations for the gradients on block diagonal form, and invert

% Map of rows to block-diag structure
% Do not consider excluded boundary equations
nnoFlux = excludeDirichlet * nno;
nnoPressure = excludeNeumann * nno;
[n,map] = sort([nnoFlux; nnoPressure]);
rows2blkDiag = sparse(1:numel(map),map,1);

clear map nnoFlux nnoPressure

% Map of columns to block-diag structure
[~,sortNode] = sort(cellNodeBlocks(:,2));
subcind = reshape(subcind(sortNode,:)',[],1);
cols2blkDiag = sparse(1:numel(subcind),subcind,1)';

% Set up block matrix for gradient contribution
grad = rows2blkDiag * [nK ; pCont ] * cols2blkDiag;

% Size of sub-systems.
sz = accumarray(n,1);

% Invert by chosen method
igrad = opt.invertBlocks(grad,sz);

clear grad n sz pCont subcind sortNode cellNodeBlocks
%% Compute basis functions for internal degrees of freedom

% Assign unit pressures in cell centers. When
% mupliplying with igrad, this will give gradients in terms of the
% pressures, that is, the basis function.

% Rearrange rows of rhs, and scale the elements
cc = rows2blkDiag * [nKCC ; pContCC];

subgrad = cols2blkDiag * igrad * (-cc); % Cell center values moved to right hand side

loc2globFace = sparse(subfnoGlob,1:numel(subfnoGlob),1,max(subfnoAll),numel(subfnoGlob));
cnoGlob = unique(cnoGlob);
loc2globCell = sparse(1:numel(cnoGlob),cnoGlob,1,numel(cnoGlob),max(cnoAll));

% The flux discretization
out.F = loc2globFace * darcy * subgrad * loc2globCell;

clear cc subgrad

%% Boundary conditions

% Similar to the cell centers, we compute sub-cell pressure gradients as a
% function of boundary conditions (Neumann or Dirichlet). This gives a form
% of basis function for the boundary conditions.

nNeu = sum(isNeumann(fnoGlob));
nDir = sum(isDirichlet(fnoGlob));
nbnd = nNeu + nDir;

% We need to operate with two sets of indices for Neumann and Dirichlet
% sub-faces: One represents their global indices, while the other gives
% their location in the local (block-diagonal) equation. The difference is
% that the latter need to take into account that the equation numbering
% is impacted by the removal of rows in the system on the boundaries (via
% the mappings excludeDirichlet and excludeNeumann)
neu_ind_eq = find(excludeDirichlet * isNeumann(fnoGlob));
neu_ind_glob = find(isNeumann(fnoGlob));
dir_ind_eq = find(excludeNeumann * isDirichlet(fnoGlob));
dir_ind_glob = find(isDirichlet(fnoGlob));

% The right hand sides should consider the equation based numbering
ccNeu = sparse(neu_ind_eq, 1:numel(neu_ind_eq), sgn(neu_ind_glob) ...
    ./nFaceNodes(fnoGlob(neu_ind_eq)), size(nK, 1), nNeu+nDir);
ccDir = sparse(dir_ind_eq,nNeu+(1:numel(dir_ind_eq)),(sgn(dir_ind_glob).^1), size(pContCC,1), nNeu+nDir);

% The matrix [ccNeu; ccDir] have one column for each subface on the
% boundary. Map this into one column for each face (including interior
% ones) in the grid.
% Here we need to use the global indices
subBndFace2AllFaces = sparse(1:nbnd,[neu_ind_glob; dir_ind_glob],1,nbnd,...
                             max(subfno)) * sparse(subfno,fnoGlob,1,max(subfno),G.faces.num);

% BoundFlux is here fluxes expressed through                         
out.boundFlux = loc2globFace * darcy * cols2blkDiag * igrad * ...
                  ( rows2blkDiag * [ccNeu ; ccDir]) * subBndFace2AllFaces; ...

% Map from sub-face to face value if desired
if opt.returnHalfEdgeValues
    out.hf2f = hf2f;
else
    out.F = hf2f * out.F;
    out.boundFlux = hf2f * out.boundFlux;
    out.div = scalarDivergence(G);
    out.A = out.div * out.F;
end

clear darcy igrad nK cols2blkDiag rows2blkDiag hf2f 
clear subfno fno cno nno subhfno subcind scalingNK scalingPcont
clear excludeDirichlet excludeNeumann pContCC

