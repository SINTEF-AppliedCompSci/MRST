function ok=checkGrid(G)
%Apply Basic Consistency Checks to MRST Grid Geometry
%
% SYNOPSIS:
%   status = checkGrid(G)
%
% PARAMETERS:
%   G - MRST Grid Structure with associate geometry.
%
% RETURNS:
%   status - Whether or not the geometry of the input grid satisfies basic
%            consistency checks such as normals not pointing in opposite
%            directions from the centroid vectors and all cells having
%            positive bulk volumes and all interfaces having positive
%            areas.
%
% SEE ALSO:
%   `grid_structure`, `computeGeometry`.

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

ok=true;
internal=sum(G.faces.neighbors~=0,2)>1;
c=G.cells.centroids(G.faces.neighbors(internal,2),:)-G.cells.centroids(G.faces.neighbors(internal,1),:);
n=G.faces.normals(internal,:);

if ~all(sum(n.*c,2)>0)
    ok=false;
    disp('Something wrong between neigbours and sgn of face normal')
end

cellno=rldecode([1:G.cells.num]',diff(G.cells.facePos));
cc=G.faces.centroids(G.cells.faces(:,1),:)-G.cells.centroids(cellno,:);
ncc=bsxfun(@times,G.faces.normals(G.cells.faces(:,1),:),(1-2*(G.faces.neighbors(G.cells.faces(:,1),2)==cellno)));
if ~all(sum(ncc.*cc,2)>0)
    ok=false;
    disp('Something wrong between cellface normals')
end

%{
vv=volumeByGaussGreens(G);
if ~all(vv>0)
    ok=false;
   disp('Greans gauss volumes wrong')
end
%}

if (~all(G.cells.volumes>0))
    ok = false;
    disp('There are som negative cell volumes');
end

if (~all(G.faces.areas>0))
    ok =false;
    disp('There are som negative face areas');
end

end

