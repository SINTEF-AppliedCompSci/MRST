%% Triangulation of a set of random points
clf
p = rand(10,2);
t = delaunayn(p);
G = triangleGrid(p,t);
plotGrid(G,'FaceColor','none');
hold on; plot(p(:,1),p(:,2),'o'); hold off; axis off;

%% Triangulation of a rectangular mesh
[x,y] = meshgrid(1:10,1:8);
t = delaunay(x(:),y(:));
G = triangleGrid([x(:) y(:)],t);
plot(x(:),y(:),'o','MarkerSize',8);
plotGrid(G,'FaceColor','none');
axis([.9 10 0.9 8]); axis off;

%%
t = delaunayn([x(:) y(:)]);
G = triangleGrid([x(:) y(:)], t);
plot(x(:),y(:),'o','MarkerSize',8);
plotGrid(G,'FaceColor','none');
axis([.9 10 0.9 8]); axis off;

%% Triangulation of a set of perturbed mesh points
clf;
N=7; M=5; K=3;
[x,y,z] = ndgrid(0:N,0:M,0:K);
x(2:N,2:M,:) = x(2:N,2:M,:) + 0.3*randn(N-1,M-1,K+1);
y(2:N,2:M,:) = y(2:N,2:M,:) + 0.3*randn(N-1,M-1,K+1);
p = [x(:) y(:) z(:)];
t = delaunayn(p);
G = tetrahedralGrid(p,t);
plotGrid(G, 'FaceColor',[.8 .8 .8]); view(-40,60); axis tight off;

%% Seamount: A standard example from Matlab
clf
load seamount
t = delaunay(x,y);
G = triangleGrid([x(:) y(:)], t);
plot(x(:),y(:),'o');
plotGrid(G,'FaceColor','none'); axis off

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
