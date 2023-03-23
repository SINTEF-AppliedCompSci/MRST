%% Tessellations of 2D space
% We consider a set of different examples of tessellations of increasing
% complexity. The first is just a standard Cartesian grid based on a
% regular quadrilaterals. The next three are n-polygonal extensions of a
% regular triangular tessellation.
colormap(.6*parula(128)+.4*ones(128,3));

%% Example 1: Cartesian grid
% This grid is formed by laying out regular quadrilaterals
[nx,ny] = deal(15,10);
[x,y] = meshgrid(linspace(0,1,nx+1),linspace(0,1,ny+1));
p = [x(:) y(:)];
n = (nx+1)*(ny+1);
I = reshape(1:n,ny+1,nx+1);
T = [
    reshape(I(1:end-1,1:end-1),[],1)';
    reshape(I(1:end-1,2:end  ),[],1)';
    reshape(I(2:end,  2:end  ),[],1)';
    reshape(I(2:end,  1:end-1),[],1)'
    ]';

G = tessellationGrid(p, T);
clf, plotGrid(G);



%% Example 2: convex/concave hexagonal tiles
% We first generate the two hexagonal tiles. Then, we identify three
% symmetry lines (that together make up a triangle for each tile) and use
% these to fit the tiles together in space.

% Construct two symmetric convex/concave hexagonal tiles
[dx, dy, dPhi] = deal(cos(pi/3),sin(pi/3), pi*15/180);
v = pi/180*[0 120 240]';
dv = [cos(v-dPhi) sin(v-dPhi) cos(v+dPhi) sin(v+dPhi)]/2;
P = [ 0           0           0           0;
      0+dv(1,1)   0+dv(1,2)   0+dv(1,3)   0+dv(1,4);
      1           0           1           0;
      1+dv(2,1)   0+dv(2,2)   1+dv(2,3)   0+dv(2,4);
      dx          dy          dx          dy;
      dx+dv(3,1)  dy+dv(3,2)  dx+dv(3,3)  dy+dv(3,4)];
P1 = P(:,1:2);
P2 = [P([1 6:-1:2],3) -P([1 6:-1:2],4)];
clf
plotCellData(tessellationGrid([P1; P2], reshape(1:12,6,2)'),(1:2)',...
    'EdgeColor','k');

% Add help triangles that give the symmetry-directions we will use to fit
% the tiles together 
P = [0 0; 1 0; dx dy];
Gt = triangleGrid([P; P([1 3 2],1), -P([1 3 2],2)],reshape(1:6,3,2)');
plotGrid(Gt,'FaceColor','none'); axis equal off;
plotGrid(Gt,'FaceColor','none','LineWidth',2);
axis equal off
% print -depsc2 convexConcave-1.eps;

%%
% To make a full tiling, we make a basic pattern consisting of four
% triangles and then use this pattern to tile as much of space as we want
[p,t,n] = deal([],[],0);
T  = reshape(1:24,6,4)';
for j=0:1
    for i=0:3
        p = [p; bsxfun(@plus,P1,[i 2*j*dy])];                          %#ok<AGROW>
        p = [p; bsxfun(@plus,P2,[i-dx (2*j+1)*dy])];                   %#ok<AGROW>
        p = [p; bsxfun(@plus,P1,[i-dx (2*j+1)*dy])];                   %#ok<AGROW>
        p = [p; bsxfun(@plus,P2,[i 2*(j+1)*dy])];                      %#ok<AGROW>
        t = [t; T+n]; n=n+24;                                          %#ok<AGROW>
    end
end

[p,~,ic] = unique(round(p*1e5)/1e5,'rows');
G = tessellationGrid(p, ic(t));
i=repmat((1:2)',G.cells.num/2,1);
clf; plotCellData(G,i(:));
%plotFaces(G,find(any(G.faces.neighbors==0,2)),'EdgeColor','r','LineWidth',2);
axis equal tight off; drawnow
% print -depsc2 convexConcave-2.eps;


%% Example 3: nonagonal tiles
% Using the same approach as in the previous example, we can make a set of
% nonagonal tiles
[dx, dy, dPhi] = deal(cos(pi/3),sin(pi/3), pi*40/180);
v = pi/180*[0 180 120 -60 240 60]';
dv = [cos(v-dPhi) sin(v-dPhi) cos(v+dPhi) sin(v+dPhi)]/3.5;
P = [...
    0           0           0           0;
    0+dv(1,1)   0+dv(1,2)   0+dv(1,3)   0+dv(1,4);
    1+dv(2,1)   0+dv(2,2)   1+dv(2,3)   0+dv(2,4);
    1           0           1           0;
    1+dv(3,1)   0+dv(3,2)   1+dv(3,3)   0+dv(3,4);
    dx+dv(4,1)  dy+dv(4,2)  dx+dv(4,3)  dy+dv(4,4);
    dx          dy          dx          dy;
    dx+dv(5,1)  dy+dv(5,2)  dx+dv(5,3)  dy+dv(5,4);
     0+dv(6,1)  0+dv(6,2)   0+dv(6,3)   0+dv(6,4);
    ];
m  = size(P,1);
P1 = P(:,1:2);
P2 = [P([1 m:-1:2],3) -P([1 m:-1:2],4)];
T  = reshape(1:4*m,m,4)';

[p,t,n] = deal([],[],0);
for j=0:2
    for i=0:4
        p = [p; bsxfun(@plus,P1,[i 2*j*dy])];                              %#ok<AGROW>
        p = [p; bsxfun(@plus,P2,[i-dx (2*j+1)*dy])];                       %#ok<AGROW>
        p = [p; bsxfun(@plus,P1,[i-dx (2*j+1)*dy])];                       %#ok<AGROW>
        p = [p; bsxfun(@plus,P2,[i 2*(j+1)*dy])];                          %#ok<AGROW>
        t = [t; T+n]; n=n+4*m;                                             %#ok<AGROW>
    end
end

[p,ia,ic] = unique(round(p*1e5)/1e5,'rows');
G = tessellationGrid(p, ic(t));
i=repmat((1:2)',G.cells.num/2,1);
clf
plotCellData(G,i(:));
%plotFaces(G,find(any(G.faces.neighbors==0,2)),'LineWidth',2);
axis tight off;



%% Example 4: pentadecagonal tiles
% In this example, we use a slightly different approach. We first create a
% regular triangular tiling, and then we go through the line segments of
% the triangles and add one by one point

% Basic patterns consisting of uniform triangles
[dx, dy] = deal(cos(pi/3),sin(pi/3));
P1 = [0 0; 1 0; dx dy];
P2 = [P1([1 3 2],1) -P1([1 3 2],2)];
[p,t,n] = deal([],[],0);
for j=0:1
    for i=0:3
        p = [p; bsxfun(@plus,P1,[i 2*j*dy])];                              %#ok<AGROW>
        p = [p; bsxfun(@plus,P2,[i-dx (2*j+1)*dy])];                       %#ok<AGROW>
        p = [p; bsxfun(@plus,P1,[i-dx (2*j+1)*dy])];                       %#ok<AGROW>
        p = [p; bsxfun(@plus,P2,[i 2*(j+1)*dy])];                          %#ok<AGROW>
        t = [t; [1:3; 4:6; 7:9; 10:12]+n]; n=n+12;                         %#ok<AGROW>
    end
end

%{ 
 % In case you just want to show two polygons as in the book
 p = [P1; P2];
 t = reshape(1:6,3,2)';
%}

%%
% We then extract the end-points on each edge and compute the angle the
% corresponding line makes with the x-axis. We will use this information to
% perturb the points. By using this orientation of the lines, we can easily
% make sure that the perturbations are introduced in the correct direction.
dPhi   = pi*35/180;
e      = reshape(t(:,[1 2 2 3 3 1])',2,[])';
v      = p(e(:,2),:) - p(e(:,1),:);
phi    = atan2(v(:,2),v(:,1));
P      = p;

T         = reshape(e',6,[])';
T         = T(:,[1 3 5]);
[Pp,~,ic] = unique(round(P*1e5)/1e5, 'rows','stable');
G         = tessellationGrid(Pp,ic(T)); 
i=repmat((1:2)',G.cells.num/2,1);
clf
plotCellData(G,i(:));
axis equal off
% print -depsc2 trigon.eps;

%%
% Add the first set of points, making an irregular hexagonal tile
pn     = p(e(:,1),:) + 1/3*[cos(phi-dPhi)  sin(phi-dPhi)];
e      = e(:,[1 1 2]);
e(:,2) = size(P,1)+(1:size(pn,1)).';
P      = [P; pn];

T         = reshape(e',9,[])';
T         = T(:,[1:2 4:5 7:8]);
[Pp,~,ic] = unique(round(P*1e5)/1e5, 'rows','stable');
G         = tessellationGrid(Pp,ic(T)); 
i=repmat((1:2)',G.cells.num/2,1);
clf
plotCellData(G,i(:));
axis equal off
% print -depsc2 hexagon.eps;

%%
% Reshape into irregular nonagons
pn     = p(e(:,1),:) + 1/3*[cos(phi-.2*dPhi)  sin(phi-.2*dPhi)];
e      = e(:,[1:3 3]);
e(:,3) = size(P,1)+(1:size(pn,1)).';
P      = [P; pn];

T         = reshape(e',12,[])';
T         = T(:,[1:3 5:7 9:11]);
[Pp,~,ic] = unique(round(P*1e5)/1e5, 'rows','stable');
G         = tessellationGrid(Pp,ic(T)); 
i=repmat((1:2)',G.cells.num/2,1);
clf
plotCellData(G,i(:));
axis equal off
% print -depsc2 nonagon.eps;

%%
% Reshape into irregular dodecagons
phi    = atan2(-v(:,2),-v(:,1));
pn     = p(e(:,4),:) + 1/3*[cos(phi-.2*dPhi)  sin(phi-.2*dPhi)];
e      = e(:,[1:4 4]);
e(:,4) = size(P,1)+(1:size(pn,1)).';
P      = [P; pn];

T         = reshape(e',15,[])';
T         = T(:,[1:4 6:9 10:14]);
[Pp,~,ic] = unique(round(P*1e5)/1e5, 'rows','stable');
G         = tessellationGrid(Pp,ic(T)); 
i=repmat((1:2)',G.cells.num/2,1);
clf
plotCellData(G,i(:));
axis equal off
% print -depsc2 dodecagon.eps;

%% Reshape into irregular pentadecagon
e      = e(:,[1:5 5]);
pn     = p(e(:,5),:) + 1/3*[cos(phi-dPhi) sin(phi-dPhi)];
e(:,5) = size(P,1)+(1:size(pn,1)).';
P      = [P; pn];

T         = reshape(e',18,[])';
T         = T(:,[1:5 7:11 13:17]);
[P,~,ic] = unique(round(P*1e5)/1e5, 'rows','stable');
G         = tessellationGrid(P,ic(T)); 
i=repmat((1:2)',G.cells.num/2,1);
clf
plotCellData(G,i(:));
axis equal off;
% print -depsc2 pentadecagon.eps;

%%
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
