%{
Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.

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

% testTopSurfaceGrid
close all

gravity([0 0 9.8])
gravity on;
test_grid=false;
int_vert=false;%true is not properly implemented for the pressuresolver
int_vert_poro=false;
semi_implicit=false;
injected_volume=4e6;h0=10;
rate=1e6;
% atestgrid for testing
if test_grid
   %nx=120;ny=200;nz=1;
   nx=5;ny=5;nz=1;
   dims= [nx,ny,nz];
   physdims=[6000,10000,25];
   g = cartGrid(dims,physdims); %, [1 1 1]);
   g.nodes.coords(:,2)=g.nodes.coords(:,2)-10000;
   %r2=(g.nodes.coords(:,1).^2+g.nodes.coords(:,2).^2);
   g.nodes.coords(:,3)=800+g.nodes.coords(:,3)-0.5*g.nodes.coords(:,1);
   g = computeGeometry(g);
   rock.perm = ones(g.cells.num,1)*100*milli*darcy;
   rock.poro = 0.3*ones(g.cells.num,1);
   grdecl = [];
   g_top = topSurfaceGrid(g, 'grdecl', grdecl);
else
   %% check makeSome2DGrid for manipulations 
   % now the two upper layers are removed
   [g_top,g,rock,grdecl]=makeSome2DGrids('sleipner');
   %% if one will experiment with flat surface
   %g_top.cells.z=700*ones(size(g_top.cells.z));
   %g_top.faces.z=700*ones(size(g_top.faces.z));
   %g_top.nodes.z=700*ones(size(g_top.faces.z));
   %rock2D = averageRock(rock, g_top);
   %g_top.nodes.z=700*ones(size(g_top.nodes.z));
   %% if one will use constant permeability
   %rock.perm = ones(g.cells.num,1)*mean(rock2D.perm);
   %rock.poro = ones(g.cells.num,1)*mean(rock2D.poro);
   %rock2D = averageRock(rock, g_top);
   clear g;
   %load norne
end
%%
% define the cartesian coordinates for nicer plotting
X=reshape(g_top.cells.centroids(:,1),g_top.cartDims(1),g_top.cartDims(2));
Y=reshape(g_top.cells.centroids(:,2),g_top.cartDims(1),g_top.cartDims(2));
Z=reshape(g_top.cells.z,g_top.cartDims(1),g_top.cartDims(2));
rock2D = averageRock(rock, g_top);
PERM=reshape(rock2D.perm,g_top.cartDims(1),g_top.cartDims(2));
%return
%%
% time data for injection rate at the moment not used
a=xlsread('data/sleipner/Injection rates Layer 9.xls');
time_values=a(12:end,3);
injectiondata=a(12:end,4:6);
density_surface=a(6,7);
density_reserviour=a(6,4);
%return
% 
figure;
plotCellData(g_top, rock2D.perm)
%%


%% define fluid - use same flux function for entire grid,
% make table for interpolation  

minH = min(g_top.cells.H);

% compute one datapoint in each column (different height in each column)
z = linspace(0, minH, g_top.cells.num)';
f = integrateVertically(rock.perm, z, g_top)./rock2D.perm;

% pick data points for making table:
ix = mcolon(1, numel(z), 10);
f = f(ix); z = z(ix);

fluid = initVEFluid(g_top, 'mu', [0.1 0.4]*centi*poise, ...
                            'val_z', z, 'val_f', f, 'rho', [600 1000],'sr',0.0);
      
% plot interpolated flux against integrated flux
state.h = rand(g_top.cells.num,1)*minH;

kr_H = integrateVertically(rock.perm(:,1), inf, g_top);
kr =integrateVertically(rock.perm(:,1), state.h, g_top);
mob2 = [kr./fluid.mu(1), (kr_H-kr)./fluid.mu(2)];
mob = fluid.mob(state);
ix = 1:100;
  
figure;    
plot(mob(ix,1).*rock2D.perm(ix),'g')
hold on;
plot(mob2(ix,1), 'k--'); 
figure
plot(mob(ix,2).*rock2D.perm(ix),'b')
hold on
 
 plot(mob2(ix,2), 'r--')
 

 
return;

