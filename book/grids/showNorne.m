%% The Norne Field
% Norne is an oil and gas field lies located in the Norwegian Sea. The
% reservoir is found in Jurrasic sandstone at a depth of 2500 meter below
% sea level. Operator Statoil and partners (ENI and Petoro) have agreed
% with NTNU to release large amounts of subsurface data from the Norne
% field for research and education purposes.  The Norne Benchmark datasets
% are hosted and supported by the Center for Integrated Operations in the
% Petroleum Industry (IO Center) at NTNU:
% http://www.ipt.ntnu.no/~norne/wiki/doku.php 
%
% Here, we will use the simulation model released as part of "Package 2:
% Full Field model" (2013) as an example of a real reservoir.

mrstModule add mrst-experimental

%% Read model
grdecl = fullfile(ROOTDIR, 'examples', 'data', 'norne', 'GSmodel.grdecl');
if ~exist(grdecl, 'file'),
   error('Model data is not available.')
end
grdecl = readGRDECL(grdecl);


%% Visualize the whole model including inactive cells
% Save ACTNUM until later and then override the ACTNUM values to obtain a
% grid that contains all cells
actnum        = grdecl.ACTNUM;
grdecl.ACTNUM = ones(prod(grdecl.cartDims),1);
G             = processGRDECL(grdecl,'checkgrid', false);

%%
% Having obtained the grid, we first plot the outline of the whole model
% and highlight all faults
newplot;
plotGrid(G,'FaceColor','none','EdgeAlpha',0.1);
plotFaces(G,find(G.faces.tag>0), ...
   'FaceColor','red','FaceAlpha',0.2, 'EdgeColor','r','EdgeAlpha',0.1);
axis off; view(-155,80); zoom(1.7);

%%
% We then distinguish the active and inactive cells
cla;
plotGrid(G,find(~actnum(G.cells.indexMap)), ...
   'FaceColor','none','EdgeAlpha',0.1);
plotGrid(G,find( actnum(G.cells.indexMap)), ...
   'FaceColor','y', 'EdgeAlpha',0.1);

%% Consider a subgrid of the model
% Prepare to extract a subgrid - show the part that will be extracted
clear ijk;
[ijk{1:3}] = ind2sub(G.cartDims, G.cells.indexMap(:));
ijk        = [ijk{:}];

[I,J,K] = meshgrid(6:15, 80:100, 1:22);
cellNo = find (ismember(ijk, [I(:), J(:), K(:)], 'rows'));
a = actnum(G.cells.indexMap); in = find(a(cellNo));
plotGrid(G,cellNo(in),'FaceColor','m');
view(-186,68), zoom(1.7)

[I,J,K] = meshgrid(9:11, 65:67, 1:22);
cellNo = find (ismember(ijk, [I(:), J(:), K(:)], 'rows'));
plotGrid(G,cellNo,'FaceColor','r');

[I,J,K] = meshgrid(29:31, 55:75, 1:22);
cellNo = find (ismember(ijk, [I(:), J(:), K(:)], 'rows'));
plotGrid(G,cellNo,'FaceColor','b');

%%
% Extract the subgrid and compute pillars
grdecl.ACTNUM = actnum;
cut_grdecl = cutGrdecl(grdecl, [6 15; 80 100; 1 22]);
g = computeGeometry(processGRDECL(cut_grdecl));

figure, clf
plotCellData(g,g.cells.volumes,'EdgeColor','k','EdgeAlpha',.1);
colormap(jet), brighten(.5)
axis tight off, view(210,15); zoom(1.2);

%%
cut_grdecl = cutGrdecl(grdecl, [9 11; 65 67; 1 22]);
g = processGRDECL(cut_grdecl);

figure,
[X,Y,Z] = buildCornerPtPillars(cut_grdecl,'Scale',true);
[x,y,z] = buildCornerPtNodes(cut_grdecl);
plot3(X',Y',Z','-k',x(:),y(:),z(:),'or');
plotGrid(g);
axis tight off, view(180,20);

%%
cut_grdecl = cutGrdecl(grdecl, [29 31; 55 75; 1 22]);
g = processGRDECL(cut_grdecl);

figure,
plotGrid(g);
plotFaces(g,find(g.faces.tag>0),'FaceColor','b','EdgeColor','b','FaceAlpha',.5);
axis tight off, view(180,27); zoom(1.7)

%{
Copyright 2009-2014 SINTEF ICT, Applied Mathematics.

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
