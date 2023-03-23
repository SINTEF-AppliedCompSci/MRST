function out = mpsa(G,C,activeNodes,varargin)
% Discretize linear elasticity and Biot by a multipoint stress
% approximation. 
%
% The implementation is modelled on mpfa.m, but everything becomes somewhat
% more complex for the vector unknown. 
%
% NOTE: For poro-mechanical problems (not pure elasticity), non-homogeneous
% boundary conditions have not yet been implemented. This applies to the
% fields.
%
% Parameters:
% G - mrst grid structure, see http://www.sintef.no/projectweb/mrst/
%     To get an overview of the data format, type 'help grid_structure'
% C - Constitutive relation. See shear_normal_stress.m for details.
% Active nodes: Which nodes should the stencil be computed for. If empty,
%   all nodes in the grid is considered. 
%
%   NOTE: The function computes and stores in memory the inverse of a 
%   block-diagonal matrix
%   (which expresses gradients on the sub-cells as a function of pressures
%   in the surrounding cells). The size of each sub-matrix will be the
%   number of sub-cells sharing the vertex * number of dimensions^2, so each
%   system will be 16x16 for Cartesian 2D, 72x72 for Cartesian 3D. For large
%   3D grids, in particular with simplices, this may be heavy in terms of
%   memory needs. For this reason, the activeNodes option may be an atractive 
%   option. For a work-around that splits the discretization into several
%   calls to mpfa, and thus reduces memory consumption, see mpsa_subgrid.
%
% Optional parametrs, defined as keyword - value pairs:
% 'eta' - Location of continuity point on the half edges, defined
%   according to Aavatsmark 2002. Between 0 and 1. Default value is 0, for
%   simplices, the value should be 1/3 (will give symmetric method, see
%   Klausen, Radu, Eigestad 2008).
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
% Returns:
%   out - a structure containing the following fields:
%     Discretization of elasticity:
%       A - discretization matrix
%       div - divergence operator
%       stress - discretized Hook's law
%       boundStress - Discretization of boundary conditions
%
%    Fields used for discretization of poro-mechanics 
%       divD - Discretized volumetric stress operator
%       gradP - Discrete pressure gradient for fluid pressure in Biot's
%          equations
%       stabDelta - Stabilitazion term used for the flow equation.
%
% Example: See mpsa_ex and biot_ex
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


opt=struct('eta', 1/3,...
    'bc',[],...
    'invertBlocks','matlab',...
    'verbose',false, ...
    'returnHalfEdgeValues',0);

if nargin < 3 || numel(activeNodes) == 0
    activeNodes = 1:G.nodes.num;
end

opt=merge_options(opt,varargin{:});
opt.invertBlocks = blockInverter(opt);

% Bookkeeping
Nd = G.griddim;
Nc = G.cells.num;

nFaceNodes = diff(G.faces.nodePos);

%% Boundary conditions
[isNeumann,isDirichlet] = classifyBoundaryFaces(G,opt.bc);


%% Subcell topology
[cnoAll, nnoAll, ~, fnoAll, subfnoAll, subhfnoAll] = createSubcellMapping(G);

localTopology = ismember(nnoAll, activeNodes);

cnoGlob = cnoAll(localTopology);
fnoGlob = fnoAll(localTopology);
subfnoGlob = subfnoAll(localTopology);

[~,~,nno] = unique(nnoAll(localTopology));
[tmpcno,~,cno] = unique(cnoAll(localTopology));
[~,~,fno] = unique(fnoAll(localTopology));
[~,~,subfno] = unique(subfnoAll(localTopology));
[~,~,subhfno] = unique(subhfnoAll(localTopology));


% Convert permeability tensor to matrix form for easy product with normal
% vector
C = C(tmpcno);
% Split the stiffness matrix into symmetric and asymmetric part.
[Csym,CAsym] = splitStiffnessMatrix(C);

clear C

%% Discretization

if opt.verbose
    t = tic;
    disp('Calculating nC products')
end
% Contribution of gradients to stress continuity (nC), both symmetric and
% asymmetric part
[nCsym,cellNodeBlocks,subcind] = multiplyNormalVectorStiffnessMatrix(G,Csym,cno,fnoGlob,nno,subhfno);
nCAsym = multiplyNormalVectorStiffnessMatrix(G,CAsym,cno,fnoGlob,nno,subhfno,1);

if opt.verbose
    toc(t)
    disp('done')
end

nsubhfno = numel(subhfno);
% Hook's law for normal vectors. This field is used only for Biot,
% specifically discretization of grad(p) (inclusion of fluid pressure in
% elasticity). For each face, it must be computed as seen from both
% neighboring cell, so do this before faces are uniquified (below)
[~,indf] = sort(reshape(repmat(subhfno,1,Nd),[],1));
hookNormal = sparse(1:(Nd*nsubhfno),indf,1) * (cat(1,nCsym{:}) + cat(1,nCAsym{:}));

nsubcells = size(cellNodeBlocks,1);

% Distance from cell centers to face continuity points
pCont = computeDistFaceCell(G,cnoGlob,fnoGlob,nno,subhfno,opt.eta);

% Unique subface numbers
[~,uniqueSubfno] = unique(subfno);

hfi = bsxfun(@minus,repmat(subfno(uniqueSubfno) * Nd,1,Nd),fliplr(1:Nd)-1);
fi = bsxfun(@minus,repmat(fno(uniqueSubfno) * Nd,1,Nd),fliplr(1:Nd)-1);
hf2f = sparse(fi,hfi,1);
hookSymCell = cell(Nd,1);
hookASymCell = hookSymCell;
for iter1 = 1 : Nd
    hookSymCell{iter1} = nCsym{iter1}(uniqueSubfno,:);
    hookASymCell{iter1} = nCAsym{iter1}(uniqueSubfno,:);
end

clear nCAsym

% nK and pCont was computed from both sides of a subface. Now combine them
sgn   = 2*(cnoGlob == G.faces.neighbors(fnoGlob,1)) -1;
pairOverSubfaces = sparse(subfno,subhfno,sgn);

% Match stress equations from the two adjacent cells
for iter1 = 1 : Nd
    nCsym{iter1} = pairOverSubfaces * nCsym{iter1}; 
end

pCont = pairOverSubfaces * pCont;

clear pairOverSubfaces

% Right hand side of the equations (corresponding to cell center
% contributions)
sgn   = 2 * (cnoGlob == G.faces.neighbors(fnoGlob,1)) -1;

% Empty contribution from cell center displacements to face stresses (only
% strains count)
for iter1 = 1 : Nd
    nKCC{iter1} = sparse(max(subfno),max(cno)*Nd); 
end

pContCC = sparse(subfno,cno,sgn);

% Right hand side for fluid pressure gradient and stabilization term
rhsNormals = cell(Nd,1);
for iter1 = 1 : Nd
    rhsNormals{iter1} = sparse(subfno,cno,G.faces.normals(fnoGlob,iter1) ...
                        .* sgn ./nFaceNodes(fnoGlob));
end

i = reshape(bsxfun(@minus,repmat(Nd*cno,1,Nd),fliplr(1:Nd)-1)',[],1);
j = (1:numel(subhfno)*Nd)';
v = reshape(repmat(sgn,1,Nd)',[],1);
divGradP = sparse(i,j,v);

clear i j v

% Done with matching expression from the two sides of the faces
nno = nno(uniqueSubfno);
subfno = subfno(uniqueSubfno);
nsubfno = max(subfno);
fnoGlob = fnoGlob(uniqueSubfno);
subfnoGlob = subfnoGlob(uniqueSubfno);

sgn   = 2*(cnoGlob(uniqueSubfno) == G.faces.neighbors(fnoGlob,1)) -1;


% Now we can glue together Hook's law (symmetric and asymmetric parts) for
% the cell faces
hookCell = cell(Nd,1);
for iter1 = 1 : Nd
    % Could have done this with cellfun
    hookCell{iter1} = (hookSymCell{iter1} + hookASymCell{iter1});
end

% Final version of Hook's law. Multiply with gradients to get face stresses
[~,indf] = sort(repmat(subfno,Nd,1));
hook = sparse(1:Nd*numel(subfno),indf,1) * cat(1,hookCell{:});
clear hookSymCell hookASymCell hookCell indf

%% Exclude equations from boundary faces
j = find(~isNeumann(fnoGlob));
i = 1:numel(j);
excludeNeumann = sparse(i,j,1,numel(i),nsubfno);

j = find(~isDirichlet(fnoGlob));
i = 1 : numel(j);
excludeDirichlet = sparse(i,j,1,numel(i),nsubfno);

clear i j

for iter1 = 1 : Nd
    nCsym{iter1} = excludeDirichlet * nCsym{iter1};
    nKCC{iter1} = excludeDirichlet * nKCC{iter1};
end
pCont = excludeNeumann * pCont;

% Exclude boundary faces from right hand side
pContCC = excludeNeumann * pContCC;

%% Write the equations for the gradients on block diagonal form, and invert

% Map of rows to block-diag structure
% Do not consider excluded boundary equations
nnoStress = excludeDirichlet * nno;
nnoDispl = excludeNeumann * nno;
[n,map] = sort([repmat(nnoStress,Nd,1); repmat(nnoDispl,Nd,1)]);
rows2blkDiag = sparse(1:numel(map),map,1);

clear nnoStress nnoDispl map

% The columns in nC are number by covering all subcells in one cell at a
% time. To obtain a block structure, they must be centered around nodes
% instead
[~,sortNode] = sort(cellNodeBlocks(:,2));
subcind = reshape(subcind(sortNode,:)',[],1);
cols2blkDiag = sparse(subcind,1:numel(subcind),1);

nCsym = cat(1,nCsym{:}) * cols2blkDiag; 

% Displacement contiuity conditions are originally formed for a scalar
% equation. First duplicate the conditions to get one for each component
pCont = kron(eye(Nd),pCont);

[~,map] = sort(repmat(reshape(ones(Nd,1) * (1:nsubcells),[],1),Nd,1));
pCont = pCont(:,map) * cols2blkDiag;

clear map

grad = rows2blkDiag * [nCsym ; pCont ];

nNeuCond = size(nCsym,1);

% Set up block matrix for gradient contribution
sz = accumarray(n,1);

% Clear some variables before we compute the inverse gradient
clear nCsym pCont

if opt.verbose
    t = tic;
    disp('Inverting gradient')
end

% Invert by chosen method
igrad = cols2blkDiag * opt.invertBlocks(grad,sz) * rows2blkDiag;

if opt.verbose
    toc(t)
    disp('done')
end

clear pCont grad sz n cols2blkDiag rows2blkDiag

%% Compute basis functions for internal degrees of freedom

% Rearrange rows of rhs, and scale the elements
% cc = rows2blkDiag * [nKCC *scalingNK; pContCC* scalingPcont];
nKCC = cat(1,nKCC{:});
pContCC = kron(eye(Nd),pContCC);
[~,map] = sort(repmat((1:numel(unique(cnoGlob)))',Nd,1));
pContCC = pContCC(:,map);

nDirCond = size(pContCC,1);

rif = bsxfun(@minus, Nd*repmat(subfnoGlob,1,Nd), fliplr(1:Nd)-1);
cif = bsxfun(@minus, Nd*repmat((1:numel(subfnoGlob))',1,Nd), fliplr(1:Nd)-1);

loc2globFace = sparse(rif,cif,1,Nd*max(subfnoAll),Nd*numel(subfnoGlob));
cnoGlob = unique(cnoGlob);

ric = bsxfun(@minus, Nd*repmat((1:numel(cnoGlob))',1,Nd), fliplr(1:Nd)-1);
cic = bsxfun(@minus, Nd*repmat(cnoGlob,1,Nd), fliplr(1:Nd)-1);
loc2globCell = sparse(ric, cic,1, Nd*numel(cnoGlob),Nd*max(cnoAll));


cc =  [nKCC ; pContCC];
out.stress =  loc2globFace * hook *  igrad * (-cc) * loc2globCell; 

clear pContCC nKCC

%% Boundary conditions
nNeuBnd = sum(isNeumann(fnoGlob)) * Nd;
nDirBnd = sum(isDirichlet(fnoGlob)) * Nd;
nbnd = nNeuBnd + nDirBnd;

fnoLoc = repmat(fnoGlob,Nd,1);

% Neumann BC in terms of the block diagonal system (e.g. faces with
% Dirichlet conditions are removed)
neuCol_eq = find(excludeDirichlet * isNeumann(fnoGlob));
% Neumann BC in terms of global numebering of unknowns
neuCol_glob = find(isNeumann(fnoGlob));

neuCol_eq = reshape(bsxfun(@plus,neuCol_eq * ones(1,Nd),nNeuCond/Nd*((1:Nd)-1))',[],1);

% Values for implementation of Neumann boundary conditions. Refers to
% global field to get the signs correct
neuVal = reshape(repmat(sgn(neuCol_glob),1,Nd)',[],1) ...
    ./nFaceNodes(fnoLoc(neuCol_eq));
ccNeu = sparse(neuCol_eq, 1:numel(neuCol_eq), neuVal, nNeuCond, nbnd);
dirCol_eq = find(excludeNeumann * isDirichlet(fnoGlob));
dirCol_glob = find(isDirichlet(fnoGlob));

dirCol_eq = reshape(bsxfun(@plus,dirCol_eq * ones(1,Nd),nDirCond/Nd*((1:Nd)-1))',[],1);
dirVal = reshape(repmat(sgn(dirCol_glob),1,Nd)',[],1);
ccDir = sparse(dirCol_eq,nNeuBnd+(1:numel(dirCol_eq)),dirVal,nDirCond,nNeuBnd+nDirBnd);

clear neuCol neuCol2 neuVal dirCol dirCol2 dirVal

% Mapping from boundary counting to global counting. Uses global indices.
tmp = [neuCol_glob; dirCol_glob];
cols = reshape(bsxfun(@minus,repmat(tmp * Nd,1,Nd),fliplr(1:Nd)-1)',[],1);
allFace2boundFace = sparse(1:nbnd,cols,1,nbnd,max(subfnoAll)*Nd);
map = @(i) reshape(bsxfun(@minus,repmat(i*Nd,1,Nd),fliplr(1:Nd)-1)',[],1);
face2halfFace =  sparse(map(subfnoGlob),map(fnoGlob),1,Nd*max(subfnoAll),Nd*max(fnoAll));
out.boundStress = loc2globFace * hook * igrad * [ccNeu ; ccDir] * allFace2boundFace * face2halfFace;

clear hook face2halfFace allFace2boundFace tmp cols map

%% Coupling and stabilization terms for Biot

% Right hand side for the pressure gradient term 
for iter1 = 1 : Nd
    rhsNormals{iter1} = excludeDirichlet * rhsNormals{iter1};
end

rhsNormals = -[cat(1,rhsNormals{:}); sparse(nsubfno * Nd -nNeuBnd,numel(cnoGlob))];

loc2globCellsScalar = sparse(1:numel(cnoGlob),cnoGlob,1,numel(cnoGlob),Nc);

% Discretization of the term -grad(p) (impact of flow on elasticity)
out.gradP = loc2globCell' * divGradP * hookNormal * igrad * rhsNormals * ...
                    loc2globCellsScalar;

clear hookNormal

% Take divergence of mechanics and gradient to get divergence of
% displacement and stabilization term, respectively
if Nd == 2
    trace = [1 4];
else
    trace = [1 5 9];
end
[i,j] = ndgrid((1:nsubcells)',trace);
j = bsxfun(@plus,j,cumsum(Nd^2 * ones(nsubcells,1)) - Nd^2);

[~,cn] = gridCellNodes(G,cnoGlob);
cn = diff(cn);
cvol = G.cells.volumes(cnoGlob) ./ cn;

v = repmat(cvol(cellNodeBlocks(:,1)),1,Nd);
% Scale with subcell areas, and also add subcells together
divOp = sparse(cellNodeBlocks(:,1),(1:nsubcells)',1)* sparse(i,j,v);

% Impact of elasticity on flow equation (div(u))
out.divD = loc2globCellsScalar' * divOp * igrad * (-cc) * loc2globCell;
% Stabilization term, this is crucial when time step and compressibility
% goes to zero
out.stabDelta = -loc2globCellsScalar' * divOp * igrad * rhsNormals * loc2globCellsScalar;
clear igrad gradPsub i j cvol cn v rhsNormals

%%
% NOTE: Boundary conditions for the Biot terms have not been implemented.
% Should be straightforward extension of flow / elasticity.
%%

if opt.returnHalfEdgeValues
    out.hf2f = hf2f;
else
    out.stress = hf2f * out.stress;
    out.boundStress = hf2f * out.boundStress;
    out.div = vectorDivergence(G);
    out.A = out.div * out.stress;
end

