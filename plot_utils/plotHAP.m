function plotHAP(G,interpFace,varargin)
%plot interpolating points on faces for the given interpFace in 2D

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

xf=G.faces.centroids;
ys=interpFace.coords;
d=sqrt(dot(xf-ys,xf-ys,2));
ind=d>0.5*G.faces.areas;
plotGrid(G,'facecolor','none');axis image; hold on
x=interpFace.coords(:,1);
y=interpFace.coords(:,2);
plot(x(ind),y(ind),'r.',varargin{:});
plot(x(~ind),y(~ind),'b.',varargin{:});
end

% x=interpFace.coords(theFaces,1);
% y=interpFace.coords(theFaces,2);
% z=interpFace.coords(theFaces,3);
% xc=G.cells.centroids(mycell,1);
% yc=G.cells.centroids(mycell,2);
% zc=G.cells.centroids(mycell,3);
% ind=convhull(x,y,z);
% figure, plotGrid(G,mycell);alpha(0.2);hold on;
% plot3(xc,yc,zc,'r.','markersize',20);
% plot3(x,y,z,'b.','markersize',20);
% view(3);grid on
% hand=trisurf(ind,x,y,z,'facecolor','y');
% alpha(hand,0.5)
%
%
% x=interpFace.coords(theFaces,1);
% y=interpFace.coords(theFaces,2);
% xc=G.cells.centroids(mycell,1);
% yc=G.cells.centroids(mycell,2);
% ind=convhull(x,y);
% figure, plotGrid(G,mycell);hold on;
% plot(xc,yc,'r.','markersize',20);
% plot(x,y,'b.','markersize',20);
% grid on
% plot(x(ind),y(ind))
