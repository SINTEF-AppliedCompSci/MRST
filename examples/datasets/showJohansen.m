%% Visualizing the Johansen Data Set
% The Johansen formation is a candidate site for large-scale CO2 storage
% offshore the south-west coast of Norway. The
% <http://www.sintef.no/MatMoRa MatMoRA project> has developed a set of
% <http://www.sintef.no/Projectweb/MatMorA/Downloads/Johansen/ geological
% models> based on available seismic and well data. Herein, we
% will inspect one instance of the model in more detail.

%% Faults and active/inactive cells
% We start by reading the model from a file in the Eclipse format (GRDECL),
% picking the sector model with five vertical layers in the Johansen
% formation and with five shale layers above and one below.
%
% The file contains both active and inactive cells. We take the chance that
% the inactive cells do not contain garbage and visualize the whole model,
% marking with red color all faults found during processing.
dpath = getDatasetPath('johansen');
sector = fullfile(dpath, 'NPD5');
filename = [sector, '.grdecl'];

grdecl        = readGRDECL(filename);  clear filename
actnum        = grdecl.ACTNUM;
grdecl.ACTNUM = ones(prod(grdecl.cartDims),1);
G             = processGRDECL(grdecl, 'checkgrid', false);

%%
% Plot the results
clf, subplot('position',[0.025 0.025 0.95 0.95]);
h = plotGrid(G,'FaceColor',.85*[1 1 1], 'FaceAlpha',.4,'EdgeColor','none');
plotFaces(G,find(G.faces.tag>0),'FaceColor','r','EdgeColor','none');
axis tight off; view(-145,60);
set(gca,'ZDir','normal'), camlight headlight, set(gca,'Zdir','reverse');

%%
% Next we mark the active part of the model
delete(h);
plotGrid(G,find(~actnum(G.cells.indexMap)), ...
         'FaceColor',.85*[1 1 1], 'FaceAlpha',.4,'EdgeColor','none');
plotGrid(G,find(actnum(G.cells.indexMap)), ...
         'FaceColor', 'b', 'FaceAlpha', 0.4,'EdgeColor','none');
view(15,40);
%set(gca,'ZDir','normal'), camlight headlight, set(gca,'Zdir','reverse');

%% Height map
% It is only meaningful to show a height map of the active cells.
% Therefore, to inspect only the active model, we reset the ACTNUM field to
% its original values and recreate the grid. Now, inactive cells will be
% ignored and we therefore get a different unstructured grid.
grdecl.ACTNUM = actnum; clear actnum;
G = processGRDECL(grdecl); clear grdecl;
G = computeGeometry(G);

% Plotting a height map of the field using the z-component of the centroids
% of the cells
clf,
plotCellData(G,G.cells.centroids(:,3),'EdgeColor','k','EdgeAlpha',0.1);
colorbar, view(3), axis tight off, view(-20,40), zoom(1.2)
set(gca,'Clipping','off')

%% Porosity
% The porosity data are given with one value for each cell in the model. We
% read all values and then pick only the values corresponding to active
% cells in the model.
clf
p = reshape(load([sector, '_Porosity.txt'])', prod(G.cartDims), []);
poro = p(G.cells.indexMap); clear p
plotCellData(G, poro,'EdgeColor','k','EdgeAlpha',0.1);
colorbarHist(poro,[0.09 0.31],'West',50); view(-45,15), axis tight off, zoom(1.2)
set(gca,'Clipping','off')

%%
% From the plot, it seems like the formation has been pinched out and only
% contains the shale layers in the front part of the model. We verify this
% by plotting a filtered porosity field in which all values smaller than or
% equal 0.1 have been taken out.
clf
view(-15,40)
plotGrid(G,'FaceColor','none','EdgeAlpha',0.1);
plotCellData(G, poro, poro>0.1, 'EdgeColor','k','EdgeAlpha',0.1);
colorbarHist(poro(poro>.1),[0.09 0.31],'West',50);
axis tight off, zoom(1.45), camdolly(0,.2,0);
set(gca,'Clipping','off')

%% Permeability
% The permeability is given as a scalar field (Kx) similarly as the
% porosity. The tensor is given as K = diag(Kx, Kx, 0.1Kx) and we therefore
% only plot the x-component, Kx, using a logarithmic color scale.
clf
K = reshape(load([sector, '_Permeability.txt']')', prod(G.cartDims), []);
perm = bsxfun(@times, [1 1 0.1], K(G.cells.indexMap)).*milli*darcy; clear K;
rock = makeRock(G, perm, poro);
p  = log10(rock.perm(:,1));
plotCellData(G,p,'EdgeColor','k','EdgeAlpha',0.1);
view(-45,15), axis tight off, zoom(1.2)
set(gca,'Clipping','off')

% Manipulate the colorbar to get the ticks we want
cs = [0.01 0.1 1 10 100 1000];
caxis(log10([min(cs) max(cs)]*milli*darcy));
h = colorbarHist(p,caxis,'West',100);
set(h, 'XTick', 0.5, 'XTickLabel','mD', ...
   'YTick', log10(cs*milli*darcy), 'YTickLabel', num2str(cs'));

%%
% To show more of the permeability structure, we strip away the shale
% layers, starting with the layers with lowest permeability on top.
clf
idx = p>log10(0.01*milli*darcy);
plotGrid(G,'FaceColor','none','EdgeAlpha',0.1);
plotCellData(G, p, idx, 'EdgeColor','k', 'EdgeAlpha', 0.1);
view(-20,35), axis tight off, zoom(1.45), camdolly(0,.12,0);
set(gca,'Clipping','off')

cs = [0.1 1 10 100 1000];
cx = log10([.05 1010]*milli*darcy); caxis(cx);
h = colorbarHist(p(idx),cx,'West',50);
set(h, 'XTick', 0.5, 'XTickLabel','mD', ...
   'YTick', log10(cs*milli*darcy), 'YTickLabel', num2str(cs'));
         
%%
% Then we also take away the lower shale layer and plot the permeability
% using a linear color scale.
clf
idx = p>log10(0.1*milli*darcy);
plotGrid(G,'FaceColor','none','EdgeAlpha',0.1);
plotCellData(G, p, idx, 'EdgeColor','k', 'EdgeAlpha', 0.1);
view(-20,35), axis tight off, zoom(1.45), camdolly(0,.12,0);
set(gca,'Clipping','off')

cs = [10 100 1000];
cx = log10([9 1010]*milli*darcy); caxis(cx);
h = colorbarHist(p(idx),cx,'West',50); caxis(cx);
set(h, 'XTick', 0.5, 'XTickLabel','mD', ...
   'YTick', log10(cs*milli*darcy), 'YTickLabel', num2str(cs'));

%% Well
% Finally, we read the well data and plot the injection well at the correct
% position.
w = load([sector, '_Well.txt']);
W = verticalWell([], G, rock,  w(1,1), w(1,2), w(1,3):w(1,4),  ...
                 'InnerProduct', 'ip_tpf', ...
                 'Radius', 0.1, 'name', 'I');
plotWell(G,W,'height',500,'color','b');


%% Copyright notice

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
