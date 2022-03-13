function test_newplotting()
%Demo for new plotting routines
%
% SYNOPSIS:
%   testPlotting()
%
% NEW FEATURES
%   - Coarse grid plotting
%   - Plotting of faces in both 2D and 3D
%   - Plotting of corner nodes for coarse and fine grids in 2D and 3D.
%
% The new plotting routines feature seamless support for plotting of coarse
% grids and fine grids in the same package. Plotting of sub-grid gemoetry
% of coarse faces happens behind the scenes by extracting necessary data
% from the (new) coarse grid field 'parent'.
% If the 'parent' field is present in the input grid, the plotting options
% are interpreted as intended modify coarse grid edges, faces and vertices.
% To plot fine grid edges, faces or vertices, a grid without the 'parent'
% field must be supplied.
%

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

% Written by Jostein R. Natvig, SINTEF Applied Mathematics.

%% test coarse grid plotting
   test2D()
   test3D()
end

function test2D()
   %% 2D grid
   G = twister(cartGrid([32,32]));
   p = partitionCartGrid(G.cartDims, G.cartDims/8);
   %cg = generateCoarseGrid(G, p);

   pf = [0;partitionCartGrid(G.cartDims, [8, 8])];
   pf = max(pf(G.faces.neighbors+1), 2);
   cg = generateCoarseGrid(G, p, pf);

   plot1(cg);
   plot2(cg);
   plot3(cg);
   plot4(cg);

   plot2d1(cg);
   plot2d2(cg);
end
function test3D()
   %% 3D grid
   G = twister(cartGrid([32,32,32]));
   p = partitionCartGrid(G.cartDims, [2,2,2]);

   pf = [0;partitionCartGrid(G.cartDims, [4,2,2])];
   pf = max(pf(G.faces.neighbors+1), [], 2);
   cg = generateCoarseGrid(G, p, pf);
   %%{
   plot1(cg);
   plot2(cg);
   plot3(cg);
   plot4(cg);
   %}
   plot2d1(cg);
   plot2d2(cg);
end

%% 2D TEST for plotPatches
%% 2D TEST for plotFaceOutline

%% 2D TESTS for plotFaces
function plot2d1(cg)
   f = figure;
   plotFaces(cg, 1:cg.faces.num,'linestyle','--', 'marker','o', ...
                'markerfacecolor','b');
   axisTight off
   view(cg.griddim)
   print('-dpng','-r500','plot-2d-1.png');
   close(f)
end

%% 2D TESTS for plotGrid
function plot2d2(cg)
   f = figure;
   plotGrid(cg.parent, find(cg.partition==1), 'facec','none', ...
               'edgea', 0.1, 'edgec', 'b');
   plotGrid(cg, 1, 'facea', 0.5, 'edgea', 0.2, 'edgec', 'r');
   axisTight off
   view(cg.griddim)
   print('-dpng','-r500','plot-2d-2.png');
   close(f);
end

%% 2D TESTS for plotCellData

%% 2D TESTS for plotFaceData (?)

%% 3D TEST for plotPatches
%% 3D TEST for plotFaceOutline

%% 3D TESTS for plotFaces
function plot1(cg)
   % Plot all faces if 'faces' argument IS NOT given.
   f = figure;
   plotFaces(cg,        'edgea',0.3,'edgec','r');
   plotFaces(cg.parent, 'edgea',0.05,'edgec','g');
   view(cg.griddim); axisTight; axis off
   axisTight off
   print('-dpng','-r500','plot1.png');
   close(f);
end
function plot2(cg)
   % Plot a selection of faces if 'faces' argument IS given
   f = figure;
   plotFaces(cg,        1:4,  'edgea',0.3,'edgec','r');
   plotFaces(cg.parent, 1:40, 'edgea',0.2,'edgec','b');
   if cg.griddim==2, view(2), else view(130, 40); end
   axisTight off
   print('-dpng','-r500','plot2.png');
   close(f);
end


%% 3D TESTS for plotGrid
function plot3(cg)
   % Plot all (outer) edges if coarse grid, no æfacesæargument given.
   f = figure;
   plotGrid(cg,        'facec', 'g',   'facea', 0.3, 'edgec', 'k', ...
                          'linewidth', 2);
   plotGrid(cg.parent, find(cg.partition==1), 'facec','none', ...
                          'edgea', 0.2, 'edgec', 'b');
   if cg.griddim==2, view(2), else view(40, 45); end
   axisTight off
   print('-dpng','-r500','plot3.png');
   close(f);
end


%% 3D TESTS for plotCellData
function plot4(cg)
   f = figure;
   plotCellData(cg, rand(cg.cells.num, 1), 1:cg.cells.num-1);
   if cg.griddim==2, view(2), else view(110, 45); end
   camlight
   axisTight off
   print('-dpng','-r500','plot4.png');
   close(f);
end

function axisTight(varargin)
   axis auto tight;
   a = axis;
   d = a(2:2:end)-a(1:2:end);
   d = reshape([-d(:),d(:)]', size(a));
   axis(a + 0.01*d);
   axis(varargin{:});
end

