%% Exercise 6.3.1
mrstModule add mimetic
edit mimeticExample1


%% Exercise 6.3.2
% Revisit the quarter five spot
mrstModule add mimetic incomp
[nx,ny] = deal(64);
G = cartGrid([nx,ny],[500,500]);
G.nodes.coords = twister(G.nodes.coords);
G = computeGeometry(G);
rock.perm = ones(G.cells.num, 1)*100*milli*darcy;
rock.poro = ones(G.cells.num, 1)*.2;

S = computeMimeticIP(G, rock);
% hT = computeTrans(G, rock);

gravity reset off
fluid = initSingleFluid('mu', 1*centi*poise,'rho', 1014*kilogram/meter^3);

pv  = sum(poreVolume(G,rock));
src = addSource([], 1, pv);
src = addSource(src, G.cells.num, -pv);

state = initResSol(G, 0.0, 1.0);
state = incompMimetic(state, G, S, fluid, 'src', src);
% state = incompTPFA(state, G, hT, fluid, 'src',src);

clf,
plotCellData(G, state.pressure);
plotGrid(G, src.cell, 'FaceColor', 'w');
axis equal tight; colormap(jet(128));

mrstModule add diagnostics
tof = computeTimeOfFlight(state, G, rock, 'src', src);
clf,
plotCellData(G, tof,'EdgeColor','none');
axis equal tight;
colormap(jet(16)); caxis([0,1]);

mrstModule add streamlines;
seed = (nx:nx-1:nx*ny).';
Sf = pollock(G, state, seed, 'substeps', 1);
Sb = pollock(G, state, seed, 'substeps', 1, 'reverse', true);
hf=streamline(Sf);
hb=streamline(Sb);
set([hf; hb],'Color',.7*[1 1 1]);
plotGrid(G,src.cell,'FaceColor','w');


%% Exercise 6.3.3
% Stencil comparison on the unstructured seamount grid
load seamount
T = triangleGrid([x(:) y(:)], delaunay(x,y));
[Tmin,Tmax] = deal(min(T.nodes.coords), max(T.nodes.coords));
T.nodes.coords = bsxfun(@times, ...
   bsxfun(@minus, T.nodes.coords, Tmin), 1000./(Tmax - Tmin));
T = computeGeometry(T);
clear x y z Tmin Tmax;

%
G = computeGeometry(cartGrid([25 25], [1000 1000]));
inside = isPointInsideGrid(T, G.cells.centroids);
G = removeCells(G, ~inside);

Gr = computeGeometry(cartGrid([250 250], [1000 1000]));
inside = isPointInsideGrid(T, Gr.cells.centroids);
Gr = removeCells(Gr, ~inside);

%
P = [];
for r = exp([-3.5:.2:0, 0, .1])
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

%
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
   S = computeMimeticIP(g{i}, rock);
   pv = sum(poreVolume(g{i}, rock));

   tmp = (g{i}.cells.centroids - repmat([450, 500],g{i}.cells.num,1)).^2;
   [~,ind] = min(sum(tmp,2));
   src{i} = addSource(src{i}, ind, -.02*pv/year);

   f = boundaryFaces(g{i});
   bc{i} = addBC([], f, 'pressure', 50*barsa);

   state{i} = incompMimetic(initResSol(g{i},0,1), ...
      g{i}, S, fluid, 'src', src{i}, 'bc', bc{i}, 'MatrixOutput', true);

   [tof{i},A{i}] = computeTimeOfFlight(state{i}, g{i}, rock,...
      'src', src{i},'bc',bc{i}, 'reverse', true);
end

% Plot pressure solutions
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

% Plot radial solutions
col = 'rbgk';
ms  = [12 8 8 2];
figure(2), clf, hold on
for i=[4 1:3]
  d = g{i}.cells.centroids - repmat(g{i}.cells.centroids(src{i}.cell,:),g{i}.cells.num,1);
  r = (sum(d.^2,2)).^.5;
  plot(r, state{i}.pressure/barsa, [col(i) '.'],'MarkerSize',ms(i));
end
axis([0 300 40 50]);
h=legend(ttext{[4 1:3]},'location','best'); set(h,'FontSize',14);
chld = get(h,'Children');
set(chld(1:3:end),'MarkerSize',20);

% Plot tof solutions
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


%% Exercise 6.3.4
% Compare TPFA and mimetic for the SAIGUP model
% Grid
grdecl = fullfile(getDatasetPath('SAIGUP'), 'SAIGUP.GRDECL');
if ~exist(grdecl, 'file')
   error('SAIGUP model data is not available.')
end
grdecl = readGRDECL(grdecl);
usys   = getUnitSystem('METRIC');
grdecl = convertInputUnits(grdecl, usys);
% G = processGRDECL(grdecl);
% G = computeGeometry(G);
mrstModule add libgeometry deckformat
G = mcomputeGeometry(processgrid(grdecl));

% Petrophysics
rock = grdecl2Rock(grdecl, G.cells.indexMap);
is_pos = rock.perm(:, 3) > 0;
rock.perm(~is_pos, 3) = 1e-6*min(rock.perm(is_pos, 3));

% Transmissibilities/innerproduct
hT = computeTrans(G, rock);
S  = computeMimeticIP(G, rock);

% Fluid
gravity reset on
fluid = initSingleFluid('mu',1*centi*poise,'rho', 1000*kilogram/meter^3);

% Wells: remember to change innerproduct for mimetic solver
nz = G.cartDims(3);
I = [ 3, 20,  3, 25,  3, 30,  5, 29];
J = [ 4,  3, 35, 35, 70, 70,113,113];
R = [ 1,  3,  3,  3,  2,  4,  2,  3]*500*meter^3/day;
[Wt,Wm] = deal([]);
for i = 1 : numel(I)
   Wt = verticalWell(Wt, G, rock, I(i), J(i), 1:nz, 'Type', 'rate', ...
      'Val', R(i), 'Radius', .1*meter, 'Comp_i', 1, ...
      'name', ['I$_{', int2str(i), '}$']);
   Wm = verticalWell(Wm, G, rock, I(i), J(i), 1:nz, 'Type', 'rate', ...
      'Val', R(i), 'Radius', .1*meter, 'Comp_i', 1, ...
      'name', ['I$_{', int2str(i), '}$'],'InnerProduct','ip_simple');
end
plotWell(G, Wt, 'height', 30, 'color', 'k');
in = numel(Wt);

I = [15, 12, 25, 21, 29, 12];
J = [25, 51, 51, 60, 95, 90];
for i = 1 : numel(I)
   Wt = verticalWell(Wt, G, rock, I(i), J(i), 1:nz, 'Type', 'bhp', ...
      'Val', 200*barsa(), 'Radius', .1*meter, ...
      'name', ['P$_{', int2str(i), '}$'], 'Comp_i',1);
   Wm = verticalWell(Wm, G, rock, I(i), J(i), 1:nz, 'Type', 'bhp', ...
      'Val', 200*barsa(), 'Radius', .1*meter, 'InnerProduct', 'ip_simple', ...
      'name', ['P$_{', int2str(i), '}$'], 'Comp_i',1);
end
plotWell(G,Wt(in+1:end),'height',30,'color','b');

% Assemble and solve system
state = initState(G, Wt, 350*barsa, 1);
tstate = incompTPFA(state, G, hT, fluid, 'wells', Wt);

mrstModule add linearsolvers
mstate = incompMimetic(state, G, S, fluid, 'wells', Wm, 'LinSolve', @callAMGCL);

% Plot pressure discrepancy
figure('position',[440 317 866 480]);
plotCellData(G, convertTo(tstate.pressure-mstate.pressure, barsa), ...
   'EdgeColor','k','EdgeAlpha',0.1);
plotWell(G, Wt(1:in),     'height', 100, 'color', 'b');
plotWell(G, Wt(in+1:end), 'height', 100, 'color', 'k');
axis off; view(-80,36)
h=colorbar; set(h,'Position',[0.88 0.15 0.03 0.67]);

% Plot discrepancy in time-of-flight in wells
mrstModule add diagnostics
ttf = computeTimeOfFlight(tstate, G, rock, 'wells', Wt)/year;
ttb = computeTimeOfFlight(tstate, G, rock, 'wells', Wt, 'reverse', true)/year;
mtf = computeTimeOfFlight(mstate, G, rock, 'wells', Wm)/year;
mtb = computeTimeOfFlight(mstate, G, rock, 'wells', Wm, 'reverse', true)/year;

ip=in+1:numel(Wt); 
for i=1:numel(ip) 
   subplot(2,3,i), plot(Wt(ip(i)).dZ,ttf(Wt(ip(i)).cells), '-o',...
      Wm(ip(i)).dZ, mtf(Wm(ip(i)).cells), '--s'); 
end

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
