%% Examples of Voronoi / Perpendicular Bisector (PEBI) grids
% In this script, we go through several examples of polyhedral grids in 3D.

%% Plot of Voronoi grid + pillars
[x,y] = meshgrid((0:2)*2*cos(pi/6),0:2);
x = [x(:); x(:)+cos(pi/6)];
y = [y(:); y(:)+sin(pi/6)];
a = [0 4.3 0 2.5 0 1.5];

clf
plot3(x,y,zeros(size(x)),'o'); view(40,30); axis tight off;
T = triangleGrid([x(:),y(:)]);
T.nodes.coords(:,3) = zeros(size(T.nodes.coords(:,1)));
plotGrid(T,'FaceColor','none'); view(40,30);
axis(a); axis equal tight off;
set(gca,'zdir','normal');

%%
G = pebi(T);
G.nodes.coords(:,3) = zeros(size(G.nodes.coords(:,1)));
plotGrid(G,'FaceColor','none','EdgeColor','r'); set(gca,'zdir','normal');

%%
newplot
plotGrid(G, 'FaceColor',[.9 .4 .4]); axis(a); axis off
hold on

G1 = G;
G1.nodes.coords(:,3) = 1.5;
G1.nodes.coords(:,1:2) = 0.925*G.nodes.coords(:,1:2);
plot3([G.nodes.coords(:,1) G1.nodes.coords(:,1)]', ...
   [G.nodes.coords(:,2) G1.nodes.coords(:,2)]', ...
   [G.nodes.coords(:,3) G1.nodes.coords(:,3)]', '-k','LineWidth',1);
set(gca,'zdir','normal'); view(40,30), 
axis(a), axis equal tight off

%%
G1 = G;
G1.nodes.coords(:,3) = 1;
G1.nodes.coords(:,1:2) = 0.95*G1.nodes.coords(:,1:2);
plotGrid(G1, 'FaceColor',[.4 .4 .9]); set(gca,'zdir','normal');
hold off

%% Make a good radial grid
newplot
P = [];
for r = exp(-3.5:0.25:0),
   [x,y,z] = cylinder(r,10); P = [P [x(1,:); y(1,:)]];          %#ok<AGROW>
end
P = unique(P','rows');
G = makeLayeredGrid(pebi(triangleGrid(P)), 5);
plotGrid(G,'FaceColor',[.8 .8 .8]); view(30,50), axis tight off

%% Make a grid draped over peaks
[x,y] = meshgrid((0:6)*2*cos(pi/6),0:7);
x = [x(:); x(:)+cos(pi/6)]; x=(x - mean(x(:)))/2;
y = [y(:); y(:)+sin(pi/6)]; y=(y - mean(y(:)))/2;
G = pebi(triangleGrid([x(:),y(:)]));
G.nodes.coords(:,3) = -peaks(G.nodes.coords(:,1),G.nodes.coords(:,2));
clf,
plotGrid(G); view(25,50), axis tight off

%{
Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.

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
