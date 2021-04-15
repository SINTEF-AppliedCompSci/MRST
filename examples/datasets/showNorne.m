%% Introduction to the Norne Model
% Norne is an oil and gas field lies located in the Norwegian Sea. The
% reservoir is found in Jurrasic sandstone at a depth of 2500 meter below
% sea level. Operator Statoil and partners (ENI and Petoro) have agreed
% with NTNU to release large amounts of subsurface data from the Norne
% field for research and education purposes.  The
% <http://www.ipt.ntnu.no/~norne/wiki/doku.php Norne Benchmark> datasets
% are hosted and supported by the Center for Integrated Operations in the
% Petroleum Industry (IO Center) at NTNU. Recently, the
% <http://www.opm-project.org OPM Initiative> released the full simulation
% model as an open data set on <https://github.com/OPM/opm-data GitHub>.
%
% Here, we will go through the grid and petrophysical properties from the
% simulation model.

mrstModule add deckformat

%% Read and process the model
% As the Norne dataset is available from the OPM project's public GitHub
% repository, we can download a suitable subset of the simulation model and
% process that subset.  Function |makeNorneSubsetAvailable| checks if the
% corner-point geometry and associate petrophysical properties (porosity,
% permeability and net-to-gross factors) is already available and downloads
% this subset if it is not.  Similarly, function |makeNorneGRDECL| creates
% a simple .GRDECL file with a known name that includes the datafiles in
% the correct order.  We refer to example |showSAIGUP| for a more in-depth
% discussion of how to read and process such input data.

if ~ (makeNorneSubsetAvailable() && makeNorneGRDECL())
   error('Unable to obtain simulation model subset');
end

grdecl = fullfile(getDatasetPath('norne'), 'NORNE.GRDECL');
grdecl = readGRDECL(grdecl);
usys   = getUnitSystem('METRIC');
grdecl = convertInputUnits(grdecl, usys);

%% Visualize the whole model including inactive cells
% Save ACTNUM until later and then override the ACTNUM values to obtain a
% grid that contains all cells
actnum        = grdecl.ACTNUM;
grdecl.ACTNUM = ones(prod(grdecl.cartDims),1);
G             = processGRDECL(grdecl,'checkgrid', false);

%%
% Having obtained the grid, we first plot the outline of the whole model
% and highlight all faults
clf
pargs = {'EdgeAlpha'; 0.1; 'EdgeColor'; 'k'};
plotGrid(G,'FaceColor','none', pargs{:});
plotFaces(G,find(G.faces.tag>0), ...
   'FaceColor','red','EdgeColor','none');
axis off tight; view(-170,80); zoom(1.3); camdolly(0,-.05,0); 
set(gca,'Clipping','off')

%%
% We then distinguish the active and inactive cells
cla;
plotGrid(G,find(~actnum(G.cells.indexMap)), 'FaceColor','none',pargs{:});
plotGrid(G,find( actnum(G.cells.indexMap)), 'FaceColor','y', pargs{:});

%% Show various subsets of the model
% Prepare to extract a subgrid - show the part that will be extracted
clear ijk;
[ijk{1:3}] = ind2sub(G.cartDims, G.cells.indexMap(:));
ijk        = [ijk{:}];
[I,J,K] = meshgrid(6:15, 80:100, 1:22);
cellNo = find (ismember(ijk, [I(:), J(:), K(:)], 'rows'));
a = actnum(G.cells.indexMap); in = find(a(cellNo));
plotGrid(G,cellNo(in),'FaceColor','m');
view(-186,68)

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

figure
plotCellData(g,g.cells.volumes, pargs{:});
colormap(jet), brighten(.5)
axis tight off, view(210,15); zoom(1.2);
set(gca,'Clipping','off')

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
axis tight off, view(180,27); zoom(1.2); camdolly(0,-0.2,0)
set(gca,'Clipping','off')

%% Extract the active part of the model
% To get the whole grid, we need to reprocess the input data. This gives a
% grid with two components: the main reservoir and a detached stack of
% twelve cells, which we ignore henceforth.
G  = processGRDECL(grdecl);
clf
plotGrid(G(1),pargs{:});
plotGrid(G(2),pargs{:},'FaceColor','r');
view(-110,50), axis tight off, set(gca,'DataAspect',[20 11 1]); zoom(1.7);
set(gca,'Clipping','off')

%% Outline petrophysical properties
% The petrophysical properties are included in the simulation model subset
% represented by function |makeNorneSubsetAvailable|.  Consequently, the
% necessary data has been read into the |grdecl| structure.  We can then
% extract the petrophysical properties. Notice that here, the rock
close all, G = G(1);
rock = grdecl2Rock(grdecl, G.cells.indexMap);

%% 
% Porosity
% First, we show the porosities as they are generated on a Cartesian grid
% that corresponds to geological 'time zero', i.e., an imaginary time at
% which all the depositional layers were stacked as horizontally on a
% Cartesian grid.
p = reshape(grdecl.PORO, G.cartDims);
clf
slice(p, 1, 1, 1); view(-135,30), shading flat,
axis equal off, set(gca,'ydir','reverse','zdir','reverse')
colorbarHist(p(p(:)>0), [.05 .35],'South',100);

%% 
% Then we show the porosities mapped onto the structural grid
clf
plotCellData(G,rock.poro, pargs{:});
axis tight off; set(gca,'DataAspect',[2 1 0.1]); 
view(-110,55); zoom(1.8);
colorbarHist(rock.poro, [.05 .35],'South',100);
set(gca,'Clipping','off')

%%
% show also the net-to-gross
clf
plotCellData(G,rock.ntg, pargs{:});
axis tight off; set(gca,'DataAspect',[1 1 0.1]); 
view(-110,55); zoom(1.8);
colorbarHist(rock.ntg,[.05 1],'South',100);
set(gca,'Clipping','off')

%%
% Permeability
% The permeability is generally a tridiagonal tensor K = diag(Kx, Ky, Kz).
% For the Norne model, the data file only specifies Kx, and then various
% manipulations are done to get the correct vertical permeability.
clf
p = log10(rock.perm(:,3));
plotCellData(G,p, pargs{:});
axis tight off; set(gca,'DataAspect',[1 1 0.1]);
view(-110,55); zoom(1.8);
set(gca,'Clipping','off')

% Manipulate the colorbar to get the ticks we want
cs = [1 10 100 1000 10000];
caxis(log10([min(cs) max(cs)]*milli*darcy));
hc = colorbarHist(p(~isinf(p)),caxis,'South',100);
set(hc, 'YTick', 0.5, 'YTickLabel','mD', ...
    'XTick', log10(cs*milli*darcy), 'XTickLabel', num2str(cs'));

%%
clf
p = log10(rock.perm(:,3));
plotCellData(G,p, pargs{:});
axis tight off; set(gca,'DataAspect',[1 1 0.1]);
view(-110,55); zoom(1.8);
set(gca,'Clipping','off')

% Manipulate the colorbar to get the ticks we want
cs = [0.1 1 10 100 1000];
caxis(log10([.1 2500]*milli*darcy));
hc = colorbarHist(p(~isinf(p)),caxis,'South',100);
set(hc, 'YTick', 0.5, 'YTickLabel','mD', ...
    'XTick', log10(cs*milli*darcy), 'XTickLabel', num2str(cs'));

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
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
