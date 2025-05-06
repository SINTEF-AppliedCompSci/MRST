%% Visualizing a realization of the SAIGUP project
% The <http://www.fault-analysis-group.ucd.ie/Projects/SAIGUP.html>SAIGUP
% project is a systematic assessment of uncertainty in reserves
% and production estimates within an objectively defined geological
% parameterisation encompassing the majority of European clastic oil
% reservoirs. A broad suite of shallow marine sedimentological reservoir
% types are indexed to continuously varying 3D anisotropy and heterogeneity
% levels. Structural complexity ranges from unfaulted to compartmentalised,
% and fault types from transmissive to sealing. Several geostatistical
% realisations each for the geologically diverse reservoir types covering
% the pre-defined parameter-space are up-scaled, faulted and simulated with
% an appropriate production strategy for an approximately 20 year period.
% Herein, we will inspect in detail one instance of the model, which can be
% downloaded from the <http://www.sintef.no/Projectweb/MRST>MRST webpage
%


%% Show the structural model
% We start by reading the model from a file in the Eclipse format (GRDECL),
% The file contains both active and inactive cells
grdecl = fullfile(getDatasetPath('SAIGUP'), 'SAIGUP.GRDECL');
grdecl = readGRDECL(grdecl);

%%
% MRST uses the strict SI conventions in all of its internal caluclations.
% The SAIGUP model, however, is provided using the ECLIPSE 'METRIC'
% conventions (permeabilities in mD and so on).  We use the functions
% <matlab:doc('getUnitSystem') getUnitSystem> and
% <matlab:doc('convertInputUnits') convertInputUnits> to assist in
% converting the input data to MRST's internal unit conventions.
usys   = getUnitSystem('METRIC');
grdecl = convertInputUnits(grdecl, usys);

%% Define geometry and rock properties
% We generate a space-filling geometry using the
% <matlab:doc('processGRDECL') processGRDECL> function and then compute a
% few geometric primitives (cell volumes, centroids, etc.) by means of the
% <matlab:doc('computeGeometry') computeGeometry> function.
G = processGRDECL(grdecl);
G = computeGeometry(G);

figure
args = {'EdgeAlpha'; 0.1; 'EdgeColor'; 'k'};
plotGrid(G,'FaceColor','none', args{:});
plotFaces(G,find(G.faces.tag>0),'FaceColor','r', args{:});
axis tight off; set(gca,'DataAspect',[1 1 0.1]); 
view(-65,60); zoom(1.1); camdolly(0,-0.2,0)
set(gca,'Clipping','off')

% The media (rock) properties can be extracted by means of the
% <matlab:doc('grdecl2Rock') grdecl2Rock> function.  
rock = grdecl2Rock(grdecl, G.cells.indexMap);


%% Porosity
% The porosity data are given with one value for each cell in the model.
% First, we show the porosities as they are generated on a Cartesian grid
% that corresponds to geological 'time zero', i.e., an imaginary time at
% which all the depositional layers were stacked as horizontally on a
% Cartesian grid.
p = reshape(grdecl.PORO, G.cartDims);
clf
slice(p, 1, 1, 1); view(-135,30), shading flat,
axis equal off, set(gca,'ydir','reverse','zdir','reverse')
colorbar('horiz'); caxis([0.01 0.3]);

%% 
% Likewise, we show a histogram of the porosity
hist(p(:),100);
h=get(gca,'Children'); 
set(h(1),'EdgeColor',[0 0 0.4],'FaceColor','none')
h = legend('Porosity'); set(h,'FontSize',16);

%% 
% Then we show the porosities mapped onto the structural grid
clf
plotCellData(G,rock.poro, args{:});
axis tight off; set(gca,'DataAspect',[1 1 0.1]); 
view(-65,55); zoom(1.4); camdolly(0,-0.2,0)
set(gca,'Clipping','off')
colorbar('horiz'); caxis([0.1 0.3])

%%
% show also the net-to-gross
clf
plotCellData(G,rock.ntg, args{:});
axis tight off; set(gca,'DataAspect',[1 1 0.1]); 
view(-65,60); zoom(1.4); camdolly(0,-0.2,0)
set(gca,'Clipping','off')
colorbar('horiz'); %caxis([0.1 0.3])

%% 
% and the satnum
clf
SN = grdecl.SATNUM(G.cells.indexMap); j = jet(60);
plotCellData(G,SN, args{:});
axis tight off; set(gca,'DataAspect',[1 1 0.1]); 
view(-65,60); zoom(1.4); camdolly(0,-0.2,0)
set(gca,'Clipping','off')
colorbar('horiz'); caxis([0.5 6.5]), colormap(j(1:10:end,:))

%% 
% and a plot where we split them up
clf
plotGrid(G,'FaceColor','none', args{:});
axis tight off; set(gca,'DataAspect',[1 1 0.1]);
view(-65,60); zoom(1.4); camdolly(0,-0.2,0)
set(gca,'Clipping','off')
caxis([0.5 6.5]), colormap(j(1:10:end,:))

h1 = plotCellData(G,SN, find(SN==1), args{:});
h2 = plotCellData(G,SN, find(SN==5), args{:});

%%
delete([h1,h2])
h1 = plotCellData(G,SN, find(SN==2), args{:});
h2 = plotCellData(G,SN, find(SN==4), args{:});

%%
delete([h1,h2])
h1 = plotCellData(G,SN, find(SN==3), args{:});
h2 = plotCellData(G,SN, find(SN==6), args{:});

%% Permeability
% The permeability is given as a tridiagonal tensor K = diag(Kx, Ky, Kz)
% and we therefore plot each of the horizontal and vertical components
% using a logarithmic color scale.
figure;
h = plotCellData(G,log10(rock.perm(:,1)), args{:});
axis tight off; set(gca,'DataAspect',[1 1 0.1]);
view(-65,55); zoom(1.4); camdolly(0,-0.2,0)
set(gca,'Clipping','off')

% Manipulate the colorbar to get the ticks we want
hc = colorbar('horiz');
cs = [0.001 0.01 0.1 1 10 100 1000 10000];
caxis(log10([min(cs) max(cs)]*milli*darcy));
set(hc, 'YTick', 0.5, 'YTickLabel','mD', ...
   'XTick', log10(cs*milli*darcy), 'XTickLabel', num2str(cs'));

%%
% Vertical permeability
delete(h)
plotCellData(G,log10(rock.perm(:,3)), args{:});

%% 
% Likewise, we show a histogram of the permeability
clf
Kx = rock.perm(:,1)./(milli*darcy);
Kz = rock.perm(:,3)./(milli*darcy);
hist(log10(Kx),100);
hold on
hist(log10(Kz(Kz>0)),100);
hold off
h=get(gca,'Children'); 
set(h(1),'EdgeColor',[0 0 0.4],'FaceColor','none')
set(h(2),'EdgeColor',[0.6 0 0],'FaceColor','none')
h = legend('Horizontal','Vertical'); set(h,'FontSize',16);

%%
% and one histogram per rock type
figure('Position',[0 60 900 350]);
col = jet(60); col=col(1:10:end,:);
for i=1:6 
   subplot(2,3,i); 
   hist(log10(Kx(SN==i)), 100);
   h = findobj(gca,'Type','patch');
   set(h,'FaceColor',col(i,:),'EdgeColor',col(i,:))
   set(gca,'XLim',[-3 4],'XTick',-2:2:4,'XTickLabel',num2str(10.^(-2:2:4)'));
end

%% Finally, show the multiplicator values
figure
Mx = grdecl.MULTX(G.cells.indexMap);
My = grdecl.MULTY(G.cells.indexMap);
Mz = grdecl.MULTZ(G.cells.indexMap);
plotGrid(G,'FaceColor','none','EdgeAlpha',0.1);
plotCellData(G,Mz,find(Mz<1),args{:});
axis tight off; set(gca,'DataAspect',[1 1 0.1]);
view(-65,55); zoom(1.4); camdolly(0,-0.2,0)
set(gca,'Clipping','off')
colorbar('horiz');

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
