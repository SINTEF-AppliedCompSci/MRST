%% Example: specification of boundary conditions
% In this example, we will show how to set boundary conditions to create
% pressure drop across a reservoir. For a rectangular model with
% homogeneous properties, this will create a linear pressure drop.
% The user can choose between different setups:
%
% # Rectangular model with Dirichlet boundary conditions on the left
% boundary and Neumann conditions on the right boundary, which will create
% a linear pressure drop
% # Rectangular model with hydrostatic boundary conditions and a fluid sink
% at the midpoint of the upper geological layer
% # Corner-point model with the same type of Dirichlet/Neumann boundary
% conditions on the left and right global boundaries
mrstModule add incomp

addpath('src');
setup = 3;

%% Define geometry
[nx,ny,nz] = deal(20, 20, 5);
[Lx,Ly,Lz] = deal(1000, 1000, 50);
switch setup
   case 1,
      gravity reset off
      G = cartGrid([nx ny nz], [Lx Ly Lz]);
   case 2,
      gravity reset on
      G = cartGrid([nx ny nz], [Lx Ly Lz]);
   case 3,
      gravity reset on
      G = processGRDECL(makeModel3([nx ny nz], [Lx Ly Lz/5]));
      G.nodes.coords(:,3) = 5*(G.nodes.coords(:,3) ...
                               - min(G.nodes.coords(:,3)));
end
G.nodes.coords(:,3) = G.nodes.coords(:,3) + 500;
G = computeGeometry(G);

%% Set rock and fluid data
rock.poro = repmat(.2, [G.cells.num, 1]);
rock.perm = repmat([1000, 300, 10].* milli*darcy(), [G.cells.num, 1]);
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);

%% Compute transmissibility and initialize reservoir state
hT = simpleComputeTrans(G, rock);
[mu,rho] = fluid.properties();
state = initResSol(G, G.cells.centroids(:,3)*rho*norm(gravity), 1.0);

%% Impose boundary conditions
% Our flow solvers automatically assume no-flow conditions on all outer
% (and inner) boundaries; other type of boundary conditions need to be
% specified explicitly.
[src,bc] = deal([]);
switch setup
   case 1,
      bc = fluxside(bc, G, 'EAST', 5e3*meter^3/day);
      bc = pside   (bc, G, 'WEST', 50*barsa);

      clf, plotGrid(G,'FaceColor', 'none'); view(3);
      plotFaces(G, bc.face(strcmp(bc.type,'flux')), 'b');
      plotFaces(G, bc.face(strcmp(bc.type,'pressure')), 'r');

   case 2,
      % Alternative 1: use psideh
      % bc = psideh(bc, G, 'EAST', fluid);
      % bc = psideh(bc, G, 'WEST', fluid);
      % bc = psideh(bc, G, 'SOUTH', fluid);
      % bc = psideh(bc, G, 'NORTH', fluid);
      %
      % Alternative 2: compute using face centroids
      % f = boundaryFaces(G);
      % f = f(abs(G.faces.normals(f,3))<eps);
      % bc = addBC(bc,f,'pressure',G.faces.centroids(f,3)*rho*norm(gravity));
      %
      % Alternative 3: sample from initialized state object
      f = boundaryFaces(G);
      f = f(abs(G.faces.normals(f,3))<eps);  % vertical faces only
      cif = sum(G.faces.neighbors(f,:),2);   % cells adjacent to face
      bc = addBC(bc, f, 'pressure', state.pressure(cif));

      clf, plotGrid(G,'FaceColor', 'none'); view(-40,40)
      plotFaces(G, f, state.pressure(cif)/barsa);

      ci = round(.5*(nx*ny-nx));
      ci = [ci; ci+nx*ny];
      src = addSource(src, ci, repmat(-1e3*meter^3/day,numel(ci),1));
      plotGrid(G, ci, 'FaceColor', 'y'); axis tight
      h=colorbar; set(h,'XTick', 1, 'XTickLabel','[bar]');
      
   case 3,
      clf,
      show = false(nz,1); show([1 nz])=true;
      k = repmat((1:nz),nx*ny,1); k = k(G.cells.indexMap);
      plotGrid(G,'FaceColor','none');
      plotGrid(G, show(k), 'FaceColor', [1 1 .7]);
      view(40,20);
      
      % First attempt: use pside/fluxside
      bcf = fluxside([], G, 'EAST', 5e3*meter^3/day);
      bcp = pside   ([], G, 'WEST', 50*barsa);
      hf  = plotFaces(G, bcf.face, 'b');
      hp  = plotFaces(G, bcp.face, 'r');

      % Second attempt: count faces
      bcf = fluxside([], G, 'EAST', 5e3*meter^3/day, 4:15, 1:5);
      bcp = pside   ([], G, 'WEST', 50*barsa,        7:17, []);
      delete([hf hp]);
      hf  = plotFaces(G, bcf.face, 'b');
      hp  = plotFaces(G, bcp.face, 'r');
      
      % Third attempt: do it the hard way
      f = boundaryFaces(G);
      f = f(abs(G.faces.normals(f,3))<eps);  % vertical faces only
      
      x = G.faces.centroids(f,1);
      [xm,xM] = deal(min(x), max(x));
      ff = f(x>xM-1e-5);
      bc = addBC(bc, ff, 'flux', (5e3*meter^3/day) ...
                 * G.faces.areas(ff)/ sum(G.faces.areas(ff)));
      fp = f(x<xm+1e-5);
      bc = addBC(bc, fp, 'pressure', repmat(50*barsa, numel(fp), 1));

      delete([hf hp]);
      plotFaces(G, ff, 'b'); plotFaces(G, fp, 'r');

end

%% Solve the linear system and plot result
state = simpleIncompTPFA(state, G, hT, fluid, 'bc', bc, 'src', src);
clf
plotCellData(G, convertTo(state.pressure, barsa()), 'EdgeColor', 'k');
xlabel('x'), ylabel('y'), zlabel('Depth');
view(3); axis tight;
h=colorbar; set(h,'XTick', 1, 'XTickLabel','[bar]');

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
