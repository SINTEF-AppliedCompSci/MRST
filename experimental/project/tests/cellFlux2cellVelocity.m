function cellvel = cellFlux2cellVelocity(G, cellflux)
%Undocumented internal helper in Vertical Equilibrium module.

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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

   %ff = mcolon(G.cells.facePos(1:end-1),G.cells.facePos(2:end));
   %hf = G.cells.faces(G.cells.facePos(1:end-1),G.cells.facePos(2:end),1);
   hf=G.cells.faces(:,1);
   cells = rldecode([1:G.cells.num]',diff(G.cells.facePos));
   hfc = G.faces.centroids(hf,:)-G.cells.centroids(cells,:);
   tmp = bsxfun(@times,hfc,cellflux)';
   %cellvelpos = [repmat(cells,size(tmp,1),1),repmat([1:size(tmp,1)]',size(tmp,2),1)];
   %cellvel=accumarray(cellvelpos,tmp(:));
   cellvel(:,1)=accumarray(cells,tmp(1,:));
   cellvel(:,2)=accumarray(cells,tmp(2,:));
   cellvel = bsxfun(@rdivide,cellvel,G.cells.volumes);

end
