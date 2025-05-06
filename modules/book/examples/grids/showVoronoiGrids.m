%% Illustrate the correspondence between Voronoi and Delaunay
N = 7; M=5;
[x,y]=ndgrid(0:N,0:M);
x(2:N,2:M) = x(2:N,2:M) + 0.3*randn(N-1,M-1);
y(2:N,2:M) = y(2:N,2:M) + 0.2*randn(N-1,M-1);
%
subplot(2,2,1);
plot(x,y,'ok');
axis equal tight off;
%
subplot(2,2,2);
t = delaunay(x,y);
triplot(t,x,y);
hold on; plot(x,y,'ok'); hold off;
axis equal tight off;
%
subplot(2,2,3);
voronoi(x,y);
axis equal tight off;
p=get(gca,'position'); p(2)=p(2)+0.1; set(gca,'position',p);
%
subplot(2,2,4);
h = voronoi(x,y,'b'); set(h,'LineWidth',1);
hold on; triplot(t,x,y,'r','LineWidth',0.25); hold off;
axis equal tight off;
p=get(gca,'position'); p(2)=p(2)+0.1; set(gca,'position',p);

%% Examples of Voronoi grids
clf;
%
% Regular grid
[x,y] = meshgrid(0:10,0:8);
voronoi(x,y);
axis([1 8 1 5]); axis equal off;
%%
% Regular grid
[x,y] = meshgrid(0:10,0:8);
x = [x(:); x(:)+0.5];
y = [y(:); y(:)+0.5];
voronoi(x,y);
axis([1 8 1 5]); axis equal off;

%%
% Honeycomb grids
[x,y] = meshgrid((0:10)*2*cos(pi/6),0:8);
x = [x(:); x(:)+cos(pi/6)];
y = [y(:); y(:)+sin(pi/6)];
voronoi(x,y);
axis([1 8 1 5]); axis equal off;

%% Triangulation of a set of perturbed mesh points
clf;
N=7; M=5;
[x,y] = ndgrid(0:N,0:M);
x(2:N,2:M) = x(2:N,2:M) + 0.3*randn(N-1,M-1);
y(2:N,2:M) = y(2:N,2:M) + 0.3*randn(N-1,M-1);
p = [x(:) y(:)];
T = triangleGrid(p,delaunayn(p));
V = makeLayeredGrid(pebi(T), 3);
plotGrid(V, 'FaceColor',[.8 .8 .8]); view(-40,60); axis tight off;

%% Seamount: A standard example from Matlab
clf
load seamount
V = pebi( triangleGrid([x(:) y(:)], delaunay(x,y)));
plotGrid(V,'FaceColor','none'); axis off;

%% A honeycombed grid
clf
[x,y] = meshgrid((0:4)*2*cos(pi/6),0:3);
x = [x(:); x(:)+cos(pi/6)];
y = [y(:); y(:)+sin(pi/6)];
G = triangleGrid([x(:),y(:)],delaunay(x(:),y(:)));
plotGrid(pebi(G), 'FaceColor','none'); axis equal off

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
