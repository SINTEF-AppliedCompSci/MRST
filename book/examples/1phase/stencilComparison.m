%% Triangular grid
load seamount
T = triangleGrid([x(:) y(:)], delaunay(x,y));
[Tmin,Tmax] = deal(min(T.nodes.coords), max(T.nodes.coords));
T.nodes.coords = bsxfun(@times, ...
   bsxfun(@minus, T.nodes.coords, Tmin), 1000./(Tmax - Tmin));
T = computeGeometry(T);
clear x y z Tmin Tmax;

%% Cartesian grids
G = computeGeometry(cartGrid([25 25], [1000 1000]));
inside = isPointInsideGrid(T, G.cells.centroids);
G = removeCells(G, ~inside);

Gr = computeGeometry(cartGrid([250 250], [1000 1000]));
inside = isPointInsideGrid(T, Gr.cells.centroids);
Gr = removeCells(Gr, ~inside);

%% Radial grid
P = [];
for r = exp([-3.5:.2:0, 0, .1]),
   [x,y] = cylinder(r,25); P = [P [x(1,:); y(1,:)]]; %#ok<AGROW>
end
P = unique([P'; 0 0],'rows');
[Pmin,Pmax] = deal(min(P), max(P));
P = bsxfun(@minus, bsxfun(@times, ...
   bsxfun(@minus, P, Pmin), 1200./(Pmax-Pmin)), [150 100]);
inside = isPointInsideGrid(T, P);
V = pebi( triangleGrid(P(inside,:)) );
V = computeGeometry(V);
clear P* x y;


%% Simulation loop
mrstModule add incomp diagnostics
state = cell(4,1);
src   = cell(4,1);
bc    = cell(4,1);
A     = cell(4,1);
tof   = cell(4,1);
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);
g     = {G, T, V, Gr};
for i=1:4
   rock.poro = repmat(0.2, g{i}.cells.num, 1);
   rock.perm = repmat(100*milli*darcy, g{i}.cells.num, 1);
   hT = simpleComputeTrans(g{i}, rock);
   pv = sum(poreVolume(g{i}, rock));

   tmp = (g{i}.cells.centroids - repmat([450, 500],g{i}.cells.num,1)).^2;
   [~,ind] = min(sum(tmp,2));
   src{i} = addSource(src{i}, ind, -.02*pv/year);

   f = boundaryFaces(g{i});
   bc{i} = addBC([], f, 'pressure', 50*barsa);

   state{i} = incompTPFA(initResSol(g{i},0,1), ...
      g{i}, hT, fluid, 'src', src{i}, 'bc', bc{i}, 'MatrixOutput', true);

   [tof{i},A{i}] = computeTimeOfFlight(state{i}, g{i}, rock,...
      'src', src{i},'bc',bc{i}, 'reverse', true);
end

%% Plot solutions
figure(1), clf, set(gcf,'Position', [400 420 925 400]);
ttext = {'Cartesian','Triangular','Radial','Reference'};
for i=1:4
   subplot(1,4,i),
   plotCellData(g{i},state{i}.pressure/barsa,'EdgeColor','k', 'EdgeAlpha', .1);
   plotGrid(g{i},src{i}.cell, 'FaceColor', 'w');
   title([ttext{i} ': ' num2str(g{i}.cells.num) ' cells']);
   caxis([40 50]); axis tight off
end
set(get(gca,'Children'),'EdgeColor','none');
h=colorbar('Location','South');
set(h,'position',[0.13 0.01 0.8 0.03],'YTick',1,'YTickLabel','[bar]');
set(gcf,'PaperPositionMode', 'auto');
% print -dpng stencil-p.png;

%% Plot radial solutions
col = 'rbgk';
ms  = [12 8 8 2];
figure(2), clf, hold on
for i=[4 1:3]
  d = g{i}.cells.centroids - repmat(g{i}.cells.centroids(src{i}.cell,:),g{i}.cells.num,1);
  r = (sum(d.^2,2)).^.5;
  plot(r, state{i}.pressure/barsa, [col(i) '.'],'MarkerSize',ms(i));
end
axis([0 300 40 50]);
h=legend(ttext{[4 1:3]},'Location','SouthEast'); set(h,'FontSize',14);
chld = get(h,'Children');
set(chld(1:3:end),'MarkerSize',20);
% print -depsc2 stencil-rad.eps;

%% Plot matrix structures: TPFA matrix
figure(3); clf, set(gcf,'Position', [400 420 925 400]);
for i=1:3
   subplot(1,3,i),
   spy(state{i}.A);
   title(ttext{i});
end
set(gcf,'PaperPositionMode', 'auto');
% print -depsc2 stencil-A.eps;

%% Plot solutions
figure(4), clf, set(gcf,'Position', [400 420 925 400]);
ttext = {'Cartesian','Triangular','Radial','Reference'};
for i=1:4
   subplot(1,4,i),
   plotCellData(g{i},tof{i}/year,'EdgeColor','k', 'EdgeAlpha', .1);
   plotGrid(g{i},src{i}.cell, 'FaceColor', 'w');
   title([ttext{i} ': ' num2str(g{i}.cells.num) ' cells']);
   axis tight off
end
set(get(gca,'Children'),'EdgeColor','none');
h=colorbar('Location','South');
set(h,'position',[0.13 0.01 0.8 0.03],'YTick',1,'YTickLabel','[year]');
set(gcf,'PaperPositionMode', 'auto');
% print -dpng stencil-tof.png;

%% Plot matrix structures: TPFA matrix
figure(5); clf, set(gcf,'Position', [400 420 925 500]);
for i=1:3
   subplot(2,3,i),
   spy(A{i});
   title(ttext{i});
   subplot(2,3,i+3)
   [~,q] = sort(state{i}.pressure);
   spy(A{i}(q,q));
   l = triu(A{i}(q,q),1);
   if sum(l(:))
      disp(['Discretization matrix: ' ttext{i} ', *not* lower triangular']);
   else
      disp(['Discretization matrix: ' ttext{i} ', lower triangular']);
   end
end
set(gcf,'PaperPositionMode', 'auto');
% print -depsc2 stencil-A-tof.eps;

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
