mrstModule add vemmech coarsegrid msrsb

%%

n = 50;
G = pebiGrid(1/n, [1,1]);

G.nodes.coords(:,1) = G.nodes.coords(:,1)*3;

plotGrid(G)
axis equal tight

%%

G = computeGeometry(G);
G = computeCellDimensions2(G);

%%

% Cartesian coarse grid
G_cart = cartGrid([100, 100]);
t = pi/8;
R = [cos(t), -sin(t); sin(t), cos(t)];
G_cart.nodes.coords = G_cart.nodes.coords*R;
G_cart = computeGeometry(G_cart);
p_cart = partitionUI(G_cart, [10, 10]);
p_cart = sampleFromBox(G, reshape(p_cart, G_cart.cartDims));
GC = generateCoarseGrid(G, p_cart);
GC = coarsenGeometry(GC);
GC = addCoarseCenterPoints(GC);
GC = coarsenCellDimensions(GC);

%%

rock = makeRock(G, 1,1);
T = computeTrans(G, rock);
p = partitionMETIS(G, T, 50);
GC = generateCoarseGrid(G, p);
GC = coarsenGeometry(GC);
GC = addCoarseCenterPoints(GC);
GC = coarsenCellDimensions(GC);

%%

close all

cNo = 30;
plotGrid(GC)
plotGrid(GC, cNo, 'facec', 'r');
hold on
x = [GC.cells.xMin(cNo,:);GC.cells.xMax(cNo,:)];
plot([x(1,1), x(2,1), x(2,1), x(1,1), x(1,1)], ...
     [x(1,2), x(1,2), x(2,2), x(2,2), x(1,2)], 'b', 'linew', 2);
 
f = GC.cells.faces(GC.cells.facePos(cNo):GC.cells.facePos(cNo+1)-1);
for fNo = 1:numel(f)
    x = GC.nodes.coords(GC.faces.nodePos(f(fNo)):GC.faces.nodePos(f(fNo)+1)-1,:);
    plot(x(:,1), x(:,2), '--g', 'linew', 2);
end
axis equal tight

%%

close all

cNo = 45;
plotGrid(GC)
plotGrid(GC, cNo, 'facec', 'r');
hold on
x = [GC.cells.xMin(cNo,:);GC.cells.xMax(cNo,:)];
plot([x(1,1), x(2,1), x(2,1), x(1,1), x(1,1)], ...
     [x(1,2), x(1,2), x(2,2), x(2,2), x(1,2)], 'b', 'linew', 2);
 
f = GC.cells.faces(GC.cells.facePos(cNo):GC.cells.facePos(cNo+1)-1);
for fNo = 1:numel(f)
    ff = f(fNo);
    sgn = 1 - 2*(GC.faces.neighbors(ff,1) ~= cNo);
    n = GC.faces.normals(ff,:)./GC.faces.areas(ff).*sgn;
    v = [n(2), -n(1)];
    x = [GC.faces.centroids(ff,:) + GC.faces.areas(ff,:)/2*v;
         GC.faces.centroids(ff,:) - GC.faces.areas(ff,:)/2*v];
%      x = GC.faces.centroids(ff,:);
    plot(x(:,1), x(:,2), '--g', 'linew', 2);
    
    quiver(GC.faces.centroids(ff,1), GC.faces.centroids(ff,2), n(1),n(2), 0.05, 'linew', 2)
    
%     x = [GC.faces.centroids(ff,:);
%          GC.faces.centroids(ff,:) + GC.faces.areas(ff,:).*n];
% %      x = GC.faces.centroids(ff,:);
%     plot(x(:,1), x(:,2), '--g', 'linew', 2);
    
    
end
axis equal tight

%%

n = 5;
G = computeGeometry(cartGrid([n,n,n], [1,1,1]));
t = pi/8;
R = [cos(t), -sin(t), 0; sin(t), cos(t), 0; 0, 0, 1];
G.nodes.coords = G.nodes.coords*R';

close all
plotGrid(G)
axis equal tight

%%

G = computeGeometry(G);
G = computeCellDimensions2(G);

%%

close all

cNo = 11;
plotGrid(G, 'facec', 'none')
plotGrid(G, cNo, 'facec', 'r');
hold on
x = [G.cells.xMin(cNo,:);G.cells.xMax(cNo,:)];

G.cells.dx(cNo,:)

faces = G.cells.faces(G.cells.facePos(cNo):G.cells.facePos(cNo+1)-1);
G.faces.dx

view(3)
axis equal tight

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
