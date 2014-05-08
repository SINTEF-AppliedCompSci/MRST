function testGrid(g)
% Perform simple checks for the sanity of a top-surface grid
%
% SYNOPSOS:
%   G = testGrid(G)
%
% PARAMETERS:
%   G       - Grid structure holding top-surface grid
%
% COMMENTS:
%    Performs the following simple sanity checks:
%     * Positive and correctly normalized face normals
%     * Positive and correctly normalized cell normals
%     * Positive cell volumes
%     * Positive face areas
%

%{
Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.

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

%check normals
internal= (sum(g.faces.neighbors>0,2)==2);
facecells = g.faces.neighbors(internal,:);
test_normal= g.cells.centroids(facecells(:,2),:)-g.cells.centroids(facecells(:,1),:);
wrong_face_normals=find(sum(test_normal.*g.faces.normals(internal,:),2)<0);
%test_normal.*g.faces.normals(internal,:)
if(numel(wrong_face_normals)>0)
   disp('Something suspicious about normals of faces')
end
wrong_cell_normals=find(g.cells.normals(:,3)<0)
if(numel(wrong_cell_normals)>0)
   disp('Something suspicious about normals of cells')
end
if (sum(g.cells.volumes>0)~= g.cells.num)
   disp('Something suspicious about cell volumes')
end
if (sum(g.faces.areas>0) ~= g.faces.num)
   disp('Something suspicious about face areas')
end
return
sum((sqrt(sum(g.faces.normals.^2,2))-g.faces.areas)./g.faces.areas<1e-10)
sum((sqrt(sum(g.cells.normals.^2,2))-g.cells.volumes)./g.cells.volumes>1e-6)
%% make also test for that the cell centroid is betwen maximum of face
%% centoids
