function G = repairNormals(G)
% temporary fix for non-XY 2D-grids, should be merged into computeGeometry!
% Also computes 2D cell-normals which is useful for dealing with non-planar
% grids

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

assert(G.griddim == 2);
compGeom = ~any(strcmp(G.type, 'computeGeometry'));
if size(G.nodes.coords,2)==2
    G.nodes.coords(:,3) = 0;
    compGeom = true;
end
if compGeom
    G = computeGeometry(G);
end

dispif(mrstVerbose, 'Repairing 2D grid geometry\n')

nc = G.cells.num;
% first do cell normals (all equal if planar grid)
% face-nodes come in pairs, so get cellnodes directly
nfn = diff(G.faces.nodePos);
assert(all(nfn==2));
[fno, fno_next] = deal(G.faces.nodes(1:2:end), G.faces.nodes(2:2:end));
cNo = rldecode((1:nc)', diff(G.cells.facePos));
acc = sparse(cNo, (1:numel(cNo))', 1);
cellFaces = G.cells.faces(:,1);
isPos = cNo == G.faces.neighbors(cellFaces,1);  
cellNodes = fno(G.cells.faces(:,1));
next      = fno_next(G.cells.faces(:,1));
[cellNodes(~isPos), next(~isPos)] = deal( next(~isPos), cellNodes(~isPos));
[~,nNodes] = rlencode(cNo);

% do wrt average centroid
cent  = bsxfun(@rdivide, acc*G.nodes.coords(cellNodes,:), nNodes);

[a,b] = deal(G.nodes.coords(cellNodes,:), G.nodes.coords(next,:));
nt    = cross(b-a, cent(cNo,:)-a, 2)./2;
G.cells.normals = acc*nt;
G.cells.volumes = sqrt(sum(G.cells.normals.^ 2, 2));

% get cell-normals for each face (constant for planar grids). Use average
% normal for neighboring cells
neig  = G.faces.neighbors;
nNeig = sum(neig>0, 2);
n     = [[0 0 0]; G.cells.normals./G.cells.volumes];
neig  = neig +1;
ncf   = bsxfun(@rdivide, n(neig(:,1),:)+n(neig(:,2),:), nNeig);
% re-normalization of normals may or may be not considered depending on application
[a,b] = deal(G.nodes.coords(fno,:), G.nodes.coords(fno_next,:));
G.faces.normals  = cross(b-a, ncf, 2);
G.faces.areas    = sqrt(sum(G.faces.normals.^ 2, 2));

G.type = [G.type, 'repairNormals'];
end
