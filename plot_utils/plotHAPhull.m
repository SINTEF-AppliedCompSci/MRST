function [] = plotHAPhull(G,interpFace,mycell)
%plot harmonic averaging points and its convex hull associated with a given
%cell 'mycell'
%   Detailed explanation goes here

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

if(G.griddim==2)
    x=G.cells.centroids(mycell,1);
    y=G.cells.centroids(mycell,2);
    myFaces=G.cells.faces(G.cells.facePos(mycell):G.cells.facePos(mycell+1)-1);
    hap=interpFace.coords(myFaces,:);
    ind=convhull(hap);
    plotGrid(G,mycell,'facecolor','none');hold on
    plot(x,y,'r.','markersize',30);
    plot(hap(:,1),hap(:,2),'b.','markersize',30);
    plot(hap(ind,1),hap(ind,2),'m-');axis equal tight
else
    x=G.cells.centroids(mycell,1);
    y=G.cells.centroids(mycell,2);
    z=G.cells.centroids(mycell,3);
    myFaces=G.cells.faces(G.cells.facePos(mycell):G.cells.facePos(mycell+1)-1);
    hap=interpFace.coords(myFaces,:);
    ind=convhull(hap);

    plotGrid(G,mycell,'facealpha',0.3);hold on;view(3);
    plot3(x,y,z,'r.','markersize',30);
    plot3(hap(:,1),hap(:,2),hap(:,3),'b.','markersize',30);
    trisurf(ind,hap(:,1),hap(:,2),hap(:,3),'facecolor','g','facealpha',0.5);
end
end
