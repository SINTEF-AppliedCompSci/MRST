%% Grid-orientation effect
% This example demonstrates the basic use of the MPFA-O pressure solver. To
% this end, we consider a square [0,30]x[0,30] m^2 with uniform
% permeability and a prescribed pressure drop from the left to the right
% boundary, which gives an analytical solution p(x,y)=2-x/30. The domain is
% represented on a skewed grid, and we compare the MPFA-O solution with
% solutions computed by the TPFA method and a mimetic method.

mrstModule add incomp mimetic mpfa

%% Set up simulation model
gravity off
G = cartGrid([1, 1],[1 1]);
G.nodes.coords = twister(G.nodes.coords);
G = computeGeometry(G);
X=G.nodes.coords(:,1);Y=G.nodes.coords(:,2);
X=[X;G.cells.centroids(:,1)+0.2];Y=[Y;G.cells.centroids(:,2)]
p     = [X(:), Y(:)];
t     = delaunayn(p);
G     = triangleGrid(p, t);
G = computeGeometry(G);
rock.perm = ones(G.cells.num, 1);

xfaces=find(abs(G.faces.centroids(:,1))<1e-4)
yfaces=find(abs(G.faces.centroids(:,1)-1)<1e-4)
bc=addBC([],xfaces,'pressure',2)
bc=addBC(bc,yfaces,'pressure',1)
%bc  = pside([], G, 'left',  2);
%bc  = pside(bc, G, 'right', 1);

fluid = initSingleFluid('mu' ,    1     , ...
                        'rho', 0);

%% Mimetic method
fprintf('Mimetic method\t... ')
tic
S = computeMimeticIP(G, rock);
xr1 = incompMimetic(initResSol(G, 0, 0), G, S, fluid, 'bc', bc);
toc

%% MPFA-O method
fprintf('MPFA-O method\t... ')
tic
T1  = computeMultiPointTrans(G, rock,'eta',1/3);
xr2 = incompMPFA(initResSol(G, 0, 0), G, T1, fluid, ...
                 'bc', bc,'MatrixOutput',true);
toc

%% TPFA method
fprintf('TPFA Method\t... ')
tic
T2  = computeTrans(G, rock);
xr3 = incompTPFA(initResSol(G, 0, 0), G, T2, fluid, ...
                 'bc', bc,'MatrixOutput',true);
toc

%% Plot solutions
hf2cn      = gridCellNo(G);
flux_int   = @(x) accumarray(hf2cn, abs(x.flux(G.cells.faces(:,1))));
plot_var   = @(v) plotCellData(G, v, 'EdgeColor','none');
plot_press = @(x) plot_var(x.pressure(1:G.cells.num));
plot_flux  = @(x) plot_var(convertTo(flux_int(x), meter^3/day));

clf, set(gcf,'Position',[300 250 1000 500]);

subplot(2,3,1),
plot_flux(xr1); cax = caxis; axis equal tight, title('Mimetic')

subplot(2,3,2),
plot_flux(xr2); caxis(cax);  axis equal tight, title('MPFA-O')

subplot(2,3,3),
plot_flux(xr3); caxis(cax);  axis equal tight, title('TPFA')
colorbar('Position',[.92 .58 .02 .34])

subplot(2,3,4),
plot_press(xr1); cax = caxis; axis equal tight, title('Mimetic')

subplot(2,3,5),
plot_press(xr2); caxis(cax);  axis equal tight, title('MPFA-O')

subplot(2,3,6),
plot_press(xr3); caxis(cax);  axis equal tight, title('TPFA')
colorbar('Position',[.92 .11 .02 .34])

 
%% Compute discrepancies in flux and errors in pressure
p.pressure = 2 - G.cells.centroids(:,1);
err        = @(q1, q2) norm(q1 - q2, inf);
err_press  = @(x1, x2) err(x1.pressure(1:G.cells.num), ...
                           x2.pressure(1:G.cells.num));
err_flux   = @(x1, x2) err(flux_int(x1), flux_int(x2));

fprintf(['\nFlux Difference:\n', ...
         '\to Mimetic/MPFA-O : %.15e\n',    ...
         '\to Mimetic/TPFA   : %.15e\n',    ...
         '\to MPFA-O /TPFA   : %.15e\n\n'], ...
        err_flux(xr1, xr2), err_flux(xr1, xr3), err_flux(xr2, xr3));

fprintf(['Cell Pressure Error:\n', ...
         '\to Mimetic        : %.15e\n',  ...
         '\to MPFA-O         : %.15e\n',  ...
         '\to TPFA           : %.15e\n'], ...
        err_press(xr1, p), err_press(xr2, p), err_press(xr3, p));

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
        %%
max(max(xr2.A-xr2.A'))/max(abs(xr2.A(:)))
     