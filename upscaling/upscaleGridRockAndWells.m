function [rock_cg,CG,W_cg] = upscaleGridRockAndWells(G,rock,coarseDim,W)
% [rock_cg,CG,W_cg] = upscaleGridRockAndWells(G,rock,coarseDim)
% function to upscale and make coarse grid compatible with many of our
% solvers

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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


% make partion
p = partitionUI(G, coarseDim);
% generate caurse grid
CG = generateCoarseGrid(G, p);
[nsub, sub] = subFaces(G, CG);
% add geometry to CG to be able to use solvers
% volume and volume weighted cells.centroids
CG.cells.centroids=sparse(CG.partition,1:G.cells.num,G.cells.volumes)*[G.cells.centroids,ones(G.cells.num,1)];
CG.cells.volumes=CG.cells.centroids(:,end);
CG.cells.centroids=bsxfun(@rdivide,CG.cells.centroids(:,1:end-1),CG.cells.centroids(:,end));
CG.faces.centroids=sparse(rldecode(1:numel(nsub),nsub,2),1:numel(sub),G.faces.areas(sub))*[G.faces.centroids(sub,:),ones(numel(sub),1)];
% area and area weighted faces.centroids
CG.faces.areas=CG.faces.centroids(:,end);
CG.faces.centroids=bsxfun(@rdivide,CG.faces.centroids(:,1:end-1),CG.faces.centroids(:,end));
CG.faces.normals=sparse(rldecode(1:numel(nsub),nsub,2),1:numel(sub),1)*G.faces.normals(sub,:);
% griddim
CG.griddim=G.griddim;
% fake some fields needed for solvers
CG.cells.indexMap=1:CG.cells.num;
CG.cartDims=coarseDim;
% upscale perm
rock_cg.perm=upscalePerm(G,CG,rock);
% upscale poro
rock_cg.poro=accumarray(CG.partition,poreVolume(G,rock))./accumarray(CG.partition,G.cells.volumes);
if isfield(rock, 'cr')
   rock_cg.cr=rock.cr;
end
% make
% only work for default arguments, may make more than one perforation
% per cell
W_cg =[];
for i=1:numel(W)
   W_cg = addWell(W_cg, CG, rock_cg, CG.partition(W(i).cells), ...
                  'Comp_i', W(i).compi, 'type', W(i).type, ...
                  'val', W(i).val, 'InnerProduct', 'ip_simple');
end
end
