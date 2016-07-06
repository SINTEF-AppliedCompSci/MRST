%
[X,Y] = meshgrid(-1:2); Z = 0*X;

figure('Position',[10 400 1400 300]);set(gcf,'PaperPositionMode','auto');
%
subplot(1,4,1);
[x,y]=meshgrid(-1:.1:1,0:.1:1);
z = 1+x; z(x>0)=1-x(x>0);
plot3(X,Y,Z,'k-',X',Y',Z,'k-');
hold on
plot3([0 1 1 0 0], [0 0 1 1 0], [0 0 0 0 0],'k-', 'LineWidth',2);
surf(x,y,z); hold off; axis equal tight off
view(-45,60); zoom(1.3); set(gca,'Projection','perspective')
%
subplot(1,4,2);
[x,y]=meshgrid(0:.1:2,0:.1:1);
z = x; z(x>1)=2-x(x>1);
plot3(X,Y,Z,'k-',X',Y',Z,'k-');
hold on
plot3([0 1 1 0 0], [0 0 1 1 0], [0 0 0 0 0],'k-', 'LineWidth',2);
surf(x,y,z); hold off; axis equal tight off
view(-45,60); zoom(1.3); set(gca,'Projection','perspective')
%
subplot(1,4,3);
[x,y]=meshgrid(0:0.1:1,-1:.1:1);
z = 1+y; z(y>0)=1-y(y>0);
plot3(X,Y,Z,'k-',X',Y',Z,'k-');
hold on
plot3([0 1 1 0 0], [0 0 1 1 0], [0 0 0 0 0],'k-', 'LineWidth',2);
surf(x,y,z); hold off; axis equal tight off
view(-45,60); zoom(1.3); set(gca,'Projection','perspective')
%
subplot(1,4,4);
[x,y]=meshgrid(0:0.1:1,0:.1:2);
z = y; z(y>1)=2-y(y>1);
plot3(X,Y,Z,'k-',X',Y',Z,'k-');
hold on
plot3([0 1 1 0 0], [0 0 1 1 0], [0 0 0 0 0],'k-', 'LineWidth',2);
surf(x,y,z); hold off; axis equal tight off
view(-45,60); zoom(1.3); set(gca,'Projection','perspective')

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
