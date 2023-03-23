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
sector = fullfile(getDatasetPath('johansen'), 'NPD5');
filename = [sector, '.grdecl'];
G = processGRDECL(readGRDECL(filename));


%% Porosity
% The porosity data are given with one value for each cell in the model. We
% read all values and then pick only the values corresponding to active
% cells in the model.
args = {'EdgeAlpha'; 0.1; 'EdgeColor'; 'k'};
p    = load([sector, '_Porosity.txt'])'; p = p(G.cells.indexMap);
h    = plotCellData(G, p, args{:}); view(-45,15),
axis tight off, zoom(1.15), caxis([0.1 0.3]), colorbar; 
set(gca,'Clipping','off')

%%
% From the plot, it seems like the formation has been pinched out and only
% contains the shale layers in the front part of the model. We verify this
% by plotting a filtered porosity field in which all values smaller than or
% equal 0.1 have been taken out.
delete(h), view(-15,40)
plotGrid(G,'FaceColor','none', args{:}); 
plotCellData(G, p, find(p>0.1), args{:})

%% Permeability
% The permeability is given as a scalar field (Kx) similarly as the
% porosity. The tensor is given as K = diag(Kx, Kx, 0.1Kx) and we therefore
% only plot the x-component, Kx, using a logarithmic color scale.
clf
K = load([sector '_Permeability.txt'])'; K=K(G.cells.indexMap);
plotCellData(G, log10(K), args{:});
view(-45,15), axis tight off, zoom(1.15)
set(gca,'Clipping','off')

[hc,hh] = colorbarHist(K, [0.001 1000], 'East', 100, true);
p = get(hh,'Position'); set(hh,'Position',p.*[1 1 1.5 1]);
set(get(hh,'Children'),'FaceColor',[.6 .6 .6]);
set(hc,'Fontsize',16,'FontWeight','bold', 'XTick', 0.5, 'XTickLabel','mD', ...
    'YTickLabel', {10.^(-3:3)});

%%
% To show more of the permeability structure, we strip away the shale
% layers, starting with the layers with lowest permeability on top.
clf,
plotGrid(G,'FaceColor','none',args{:});
plotCellData(G, log10(K), find(K>0.01), args{:});
view(-60,40), axis tight off, zoom(1.15)
set(gca,'Clipping','off')

[hc,hh] = colorbarHist(K(K>0.01), [0.01 1000], 'East', 100, true);
set(get(hh,'Children'),'FaceColor',[.6 .6 .6]);
set(hc,'Fontsize',16,'FontWeight','bold', 'XTick', 0.5, 'XTickLabel','mD', ...
    'YTickLabel', {10.^(-2:3)});

%%
% Then we also take away the lower shale layer and plot the permeability
% using a linear color scale.
clf; plotGrid(G,'FaceColor','none',args{:});
h = plotCellData(G, K, find(K>0.1), args{:});
view(-60,40), axis tight off, zoom(1.15)
set(gca,'Clipping','off')

[hc,hh] = colorbarHist(K(K>0.1), [0 900], 'East', 100, false);
set(get(hh,'Children'),'FaceColor',[.6 .6 .6]);
set(hc,'Fontsize',16,'FontWeight','bold', 'XTick', 0.5, 'XTickLabel','mD');

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
