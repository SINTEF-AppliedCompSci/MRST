% testTopSurfaceGrid

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


if false
g = cartGrid([3 1 3]); %, [1 1 1]);
%g = cartGrid([2 2 1]); %, [1 1 1]);
rock.perm = [ones(ceil(g.cells.num/2),1); 2*ones(floor(g.cells.num/2),1)]; 
rock.poro = 0.3*ones(g.cells.num,1);

remove = [5 6 3 ];

g = removeCells(g, remove);
%g.nodes.coords = twister(g.nodes.coords);

g = computeGeometry(g);

%rock.perm = 100*rand(g.cells.num,1)*milli*darcy;
%rock.perm = (1:g.cells.num)'; %*milli*darcy;


rock.perm(remove,:) = [];
rock.poro(remove,:) = [];




%figure
%plotGrid(g); view(3)
grdecl = [];
else
 load nograv
 grdecl = [];
 
%  cells = find(g.cells.centroids(:,1)>1000 | ...
%    (g.cells.centroids(:,2)>5200 | g.cells.centroids(:,2)<5000));
%                
% g = removeCells(g, cells);
% g = computeGeometry(g);
% 
% rock.perm(cells, :) = []; %100*milli*darcy;
% rock.poro(~cells, :) = [];

 
end

%rock.perm(:) = 100*milli*darcy;

g_top = topSurfaceGrid(g, 'grdecl', grdecl);

%g_top2 = topSurfaceGrid(g);


fluid = initVEFluid(g_top, 'mu', [0.1 0.4]*centi*poise, 'rho', [600 1000]);
resSol.h = 2*ones(g_top.cells.num,1);

 
kr_H = integrateVertically(rock.perm(:,1), resSol.h, ...
                           g_top);
   
rock2D = averageRock(rock, g_top);

mob = fluid.mob(resSol);

kr = rock2D.perm.*(mob(:,1)*fluid.mu(1));



figure;
subplot(1, 3, 1)
plotGrid(g_top)
title('g_top')
subplot(1, 3, 2)
cn_top = cellNodes(g_top);
%plotCellData(g_top, accumarray(cn_top(:,1), g_top.nodes.z(cn_top(:,3)))./accumarray(cn_top(:,1),1));
plotCellData(g_top, g_top.cells.z)
colorbar
title('g_top.cells.z')
ca = caxis;

subplot(1, 3, 3)
plotCellData(g, g.cells.centroids(:,3))
colorbar
title('g.cells.centroids')
%caxis(ca);


figure;

subplot(1, 2, 1)
plotCellData(g_top, convertTo(rock2D.perm, darcy))
colorbar
ca2 = caxis;
title('perm g_top')
subplot(1, 2, 2)

plotCellData(g, convertTo(rock.perm(:,1), darcy))
colorbar
view(3)
title('perm g')
%caxis(ca2);





%g_top.nodes.z

% g_top.columns.cells
% g_top.columns.kPos
% g_top.columns.dz

%g_top.cells.normals
