doprint = true;

%% Example of circumcirles
clf;
p = [0 0; 0.7 0.7; -0.2 1; -0.7 .8];
t = DelaunayTri(p);
c = circumcenters(t, [1;2]);
plot(p(:,1),p(:,2),'o','MarkerSize',8,'LineWidth',1);
hold on;
triplot(t,'LineWidth',1);
plot(c(:,1),c(:,2),'*r',c(:,1),c(:,2),'or'); 
theta = linspace(0,2*pi,100); xc = cos(theta); yc = sin(theta);
r1 = norm(p(1,:)-c(1,:));
r2 = norm(p(1,:)-c(2,:));
plot(c(1,1)+xc*r1, c(1,2)+yc*r1,'--');
plot(c(2,1)+xc*r2, c(2,2)+yc*r2,'--');
hold off
axis equal off;
if doprint
   print -deps2 tri-circumcircle;
end

%{
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

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
