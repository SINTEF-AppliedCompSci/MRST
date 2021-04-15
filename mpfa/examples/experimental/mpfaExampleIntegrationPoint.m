%% Example 2: Grid-orientation effect
% The multipoint flux-approximation (MPFA-O) scheme is developed to be
% consistent on grids that are not necessarily K-orthogonal. This example
% demonstrates the basic use of the MPFA-O pressure solver by applying it
% to a single-phase flow problem posed on a square [0,30]x[0,30] m^2 with
% homogeneous and isotropic permeability. To discretize the problem, we
% introduce a skewed, curvilinear grid in which the majority of the grid
% cells are not K-orthogonal. The classical TPFA scheme can therefore be
% expected to give significant grid-orientation effects. To investigate
% this, we consider a single-phase flow problem with a prescribed pressure
% drop from the left to the right boundary, which gives an analytical
% solution p(x,y)=2-x/30. We compare the MPFA-O solution with solutions
% computed by the TPFA method and a mimetic method.

mrstModule add incomp mimetic mpfa

%% Set up simulation model
eta=0;
gravity off;
gravity('reset',[0.3 0.3]);
gravity on;
g_vec=gravity();g_vec=g_vec(1:2)';
L = [1 1];
dim= [10 10];
G = cartGrid(dim,L);
G.nodes.coords = twister(G.nodes.coords);
G = computeGeometry(G);
X=G.nodes.coords(:,1);Y=G.nodes.coords(:,2);
if(G.cells.num==1)
    X=[X;G.cells.centroids(:,1)+0.2];Y=[Y;G.cells.centroids(:,2)];
end
%%{
ref     = [X(:), Y(:)];
t     = delaunayn(ref);
G     = triangleGrid(ref, t);
G = computeGeometry(G);
%}
perm=1;
rock.perm = perm*ones(G.cells.num, 1);
fluid = initSingleFluid('mu' ,    10     , ...
                        'rho', 1);
[~,rho]=fluid.properties();                    
xfaces = find(abs(G.faces.centroids(:,1))<1e-4);
yfaces = find(abs(G.faces.centroids(:,1)-1)<1e-4);
bc_left=2;bc_right=1;
bc = addBC([],xfaces,'pressure',bc_left+rho*G.faces.centroids(xfaces,:)*g_vec);
bc = addBC(bc,yfaces,'pressure',bc_right+rho*G.faces.centroids(yfaces,:)*g_vec);

%% MPFA-O method
fprintf('MPFA-O method\t... ')
tic
T1  = computeMultiPointTrans(G, rock,'eta',eta);
xr1 = incompMPFA(initResSol(G, 0, 0), G, T1, fluid, ...
                 'bc', bc,'MatrixOutput',true);
toc

%% Mimetic method
fprintf('Mimetic method\t... ')
tic
S = computeMimeticIP(G, rock);
xr2 = incompMimetic(initResSol(G, 0, 0), G, S, fluid, 'bc', bc);
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
plot_flux(xr1); cax=caxis();  axis equal tight, title('Flux: MPFA-O')

subplot(2,3,2),
plot_flux(xr2); caxis(cax); axis equal tight, title('Flux: Mimetic')

subplot(2,3,3),
plot_flux(xr3); caxis(cax);  axis equal tight, title('Flux: TPFA')
colorbar('Position',[.92 .58 .02 .34])

subplot(2,3,4),
plot_press(xr1); cax=caxis();  axis equal tight, title('Pressure: MPFA-O')

subplot(2,3,5),
plot_press(xr2); caxis(cax); axis equal tight, title('Pressure: Mimetic')

subplot(2,3,6),
plot_press(xr3); caxis(cax);  axis equal tight, title('Pressure: TPFA')
colorbar('Position',[.92 .11 .02 .34])
figure(1)

%% Compute discrepancies in flux and errors in pressure
[mu,rho]=fluid.properties();
ref = initResSol(G, 0);
ref.pressure = (bc_left -(bc_left - bc_right)*G.cells.centroids(:,1)/L(1)) + rho*G.cells.centroids*g_vec;
v = (perm/mu)*(bc_right-bc_left)/L(1);
ref.flux = G.faces.normals(:,1)*v;
err        = @(q1, q2) norm(q1 - q2, inf);
err_press  = @(x1, x2) err(x1.pressure(1:G.cells.num), ...
                           x2.pressure(1:G.cells.num));
err_flux   = @(x1, x2) err(flux_int(x1), flux_int(x2));

fprintf(['\nInternal flux error:\n', ...
         '\to Mimetic        : %.15e\n',    ...
         '\to MPFA-O         : %.15e\n',    ...
         '\to TPFA           : %.15e\n\n'], ...
        err_flux(xr1, ref), err_flux(xr2, ref), err_flux(xr3, ref));
    
    

fprintf(['Cell Pressure Error:\n', ...
         '\to Mimetic        : %.15e\n',  ...
         '\to MPFA-O         : %.15e\n',  ...
         '\to TPFA           : %.15e\n',  ...
         ], ...
         err_press(xr1, ref), err_press(xr2, ref), err_press(xr3, ref));

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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
     