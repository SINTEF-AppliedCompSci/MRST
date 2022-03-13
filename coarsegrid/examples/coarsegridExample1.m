%% Introduction to Coarse Grids in MRST
% In MRST, a coarse grid always refers to a grid that is defined as
% a partition of another grid, which is referred to as the 'fine' grid. The
% coarsegrid module defines a basic grid structure for such coarse grids
% and supplies simple tools for partitioning fine grids. In this example,
% we will show you the basics of the coarse-grid structure and how to
% define partitions of simple 2D Cartesian grids.

mrstModule add coarsegrid
plotPartition = @(G, p) plotCellData(G, p, 'EdgeColor','w','EdgeAlpha',.2);

%% 2x2 partition of a 4x4 Cartesian grid
% To define a coarse grid, we must first compute a partition. To this end,
% we use partitionCartGrid, which exploits the logical structure and creates
% a uniform partitioning in logical space.
G = computeGeometry(cartGrid([4,4]));
p = partitionCartGrid(G.cartDims, [2,2]);
plotPartition(G, p);
colormap(jet), axis tight off

%% Generate coarse-grid structure
% We can then call the routine that builds the structure representing the
% coarse grid. As a naming convention, we usually call this grid CG. The
% coarse grid can be passed on to any plotting routine in MRST. By calling
% coarsenGeometry on a grid which has been returned from processGeometry,
% we can get coarse centroids, volumes and so on.
CG = generateCoarseGrid(G, p);
CG = coarsenGeometry(CG);

cla
plotCellData(CG, (1:CG.cells.num)', 'EdgeColor','w','EdgeAlpha',.5);
plotFaces(CG, (1:CG.faces.num)', 'FaceColor','none','LineWidth', 2);

%% Show cell/block indices
% In its basic form, the structure only represents topological information
% that specifies the relationship between blocks and block interfaces, etc.
% The structure also contains information of the underlying fine grid. Let
% us start by plotting cell/block indices
tg = text(G.cells.centroids(:,1), G.cells.centroids(:,2), ...
   num2str((1:G.cells.num)'),'FontSize',16, 'HorizontalAlignment','center');
tcg = text(CG.cells.centroids(:,1), CG.cells.centroids(:,2), ...
   num2str((1:CG.cells.num)'),'FontSize',24, 'HorizontalAlignment','center');
axis off;
set(tcg,'BackgroundColor','w','EdgeColor','none');
colormap(.5*jet+.5*ones(size(jet)));

%% Show face indices of fine/coarse grids
delete([tg; tcg]);
tg = text(G.faces.centroids(:,1), G.faces.centroids(:,2), ...
   num2str((1:G.faces.num)'),'FontSize',14, 'HorizontalAlignment','center');
tcg = text(CG.faces.centroids(:,1), CG.faces.centroids(:,2), ...
   num2str((1:CG.faces.num)'),'FontSize',24, 'HorizontalAlignment','center');
set(tcg,'BackgroundColor','w','EdgeColor','none');


%% Solvers using coarse grids
% The coarse grid will in most aspects behave as the standard grids and can
% be used in many of the solvers in MRST. Utilities such as 'coarsenBC'
% makes it easier to solve problems at several scales.

mrstModule add incomp;

% Make a somewhat larger grid
G  = cartGrid([20, 20], [1 1]);
G  = computeGeometry(G);
p  = partitionCartGrid(G.cartDims,[5 2]);
CG = generateCoarseGrid(G, p);
CG = coarsenGeometry(CG);

% Trivial fluid
fluid = initSingleFluid('mu', 1, 'rho', 1);

% Uniform permeability and a pressure drop
rock.perm = ones(G.cells.num, 1);
T = computeTrans(G, rock);
state = initState(G, [], 0);
bc = pside([], G, 'left', 100*barsa);
bc = pside(bc, G, 'right', 0*barsa);

% Solve for fine scale
state = incompTPFA(state, G, T, fluid, 'bc', bc);

% Create the same trivial permeability for the coarse grid. This could be
% replaced by for example an upscaling routine from the upscaling module if
% the permeability was non-uniform.
rock.perm = ones(CG.cells.num, 1);
T_coarse = computeTrans(CG, rock);
state_coarse = initState(CG, [], 0);

% Take the existing boundary condition and sample it in the coarse grid.
bc_coarse = coarsenBC(CG, bc);
state_coarse = incompTPFA(state_coarse, CG, T_coarse, fluid, 'bc', bc_coarse);

% Plot the solutions
subplot(2,1,1)
plotCellData(G, state.pressure,'EdgeColor','k','EdgeAlpha',.2)
title('Fine-scale solution')
subplot(2,1,2)
plotCellData(CG, state_coarse.pressure,'EdgeColor','none')
plotFaces(CG,(1:CG.faces.num)','FaceColor','none');
title('Coarse-scale solution')

% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
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
