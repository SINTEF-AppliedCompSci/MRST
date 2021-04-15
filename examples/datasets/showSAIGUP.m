%% Introduction to the SAIGUP Model
% In this example, we will examine a model from the project Sensitivity
% Analysis of the Impact of Geological Uncertainties on Production
% forecasting in clastic hydrocarbon reservoirs <http://www.nr.no/saigup
% (SAIGUP)>. The model has faults, inactive cells, and disconnected
% components, but no pinchout.
%
% We will show how to read, process, and visualize the reservoir model and
% discuss some of the petrophysical properties.

%% Read and process the model
% We start by reading the model from a file in the Eclipse format (GRDECL)
% that can be read using the <matlab:doc('readGRDECL') readGRDECL>
% function. (If the model is not available the
% <matlab:doc('getDatasetPath') getDatasetPath> function will download and
% install it for you).
grdecl = fullfile(getDatasetPath('SAIGUP'), 'SAIGUP.GRDECL');
grdecl = readGRDECL(grdecl);

%%
% MRST uses the strict SI conventions in all of its internal calculations.
% The SAIGUP model, however, is provided using the ECLIPSE 'METRIC'
% conventions (permeabilities in mD and so on).  We use the functions
% <matlab:doc('getUnitSystem') getUnitSystem> and
% <matlab:doc('convertInputUnits') convertInputUnits> to assist in
% converting the input data to MRST's internal unit conventions.
usys   = getUnitSystem('METRIC');
grdecl = convertInputUnits(grdecl, usys);

grdecl    %#ok  (intentional display)
%%
% From the output of |readGRDECL|, we see that the file contains four fields:
%
% * The dimension of the underlying logical Cartesian grid (keyword
% SPECGRID, equal 40x120x20)
% * The coordinates of the pillars (keyword COORD, 6x41x121 values)
% * The coordinates along the pillars (keyword ZCORN, 8x40x120x20 values)
% * The flag for active/inactive cells (keyword ACTNUM, 40x120x20 values)
%
% Since the keyword ACTNUM is present, the model is likely to contain both
% active and inactive cells. To be able to plot both the active and the
% inactive cells, we need to override the ACTNUM field when processing the
% input, because if not, the inactive cells will be ignored when the
% unstructured grid is built.
%
% We generate a space-filling geometry using the
% <matlab:doc('processGRDECL') processGRDECL> function
%
actnum        = grdecl.ACTNUM;
grdecl.ACTNUM = ones(prod(grdecl.cartDims),1);
G             = processGRDECL(grdecl, 'Verbose', true, 'checkgrid', false);

%%
% WARNING: inactive cells often contain garbage data and should generally
% not be inspected in this manner. Here, most inactive cells are defined in
% a reasonable way. By not performing basic sanity checks on the resulting
% grid (option |'checkgrid'=false|), we manage to process the grid and
% produce reasonable graphical output. In general, however, we strongly
% advice that |'checkgrid'| remain set in its default state of |true|.
%
% To simplify the processing, a single layer of artificial cells is added
% above the top and below the bottom of the model, but not touching the
% model.  (See the tutorial <gridTutorialCornerPoint.html#1 "Read, Display,
% and Manipulate"> for more details). In the following, we therefore work
% with a 40x120x22 model.
%
% In the first phase, we process all faces with normals in the logical
% i-direction. There should be 40x120x22=105600, out of which 96778 are
% not degenerate or at a fault. In the next phase, we process the faults
% and split faces to obtain a matching grid. Here there are faults at 521
% pairs of pillars and the splitting of these results in 27752 new faces.
% If each face were split in two, we would have obtained
% 521x(20x2+2)=21882, which means that some of the faces have been split
% into at least three subfaces. The process is then repeated in the logical
% j-direction.
%
% The processing assumes that there are no faults in the logical
% k-direction and therefore processes only regular connections. In absence
% of inactive or pinched cells, there should be (20+1+4)x120x40=120000
% faces (where +4 is due to the artificial cells) in the k-direction. The
% result of the grid processing is a new structure G, outlined below

G   %#ok  (intentional display)

%% Inspect the whole model
% Having obtained the grid in the correct unstructured format, we first
% plot the outline of the whole model and highlight all faults. This model
% consist of two separated grids so that numel(G)=2
newplot; subplot('position',[0.025 0.025 0.95 0.95]);
for i=1:numel(G)
   plotGrid(G(i),'FaceColor','none','EdgeColor',[0.65 0.65 0.65], ...
           'EdgeAlpha',0.2);
   plotFaces(G(i),find(G(i).faces.tag>0),'FaceColor','red','EdgeAlpha',0.1);
end
axis off; view(-10,40); camdolly(0,-0.2,0)

%%
% Then we distinguish the active and inactive cells using the |'FaceColor'|
% property set to |'none'| for the inactive cells and to |'y'| for the
% active cells. We notice that that only |G(1)| has active cells, this is
% indicated with the warning.
cla;
for i=1:numel(G)
   hi = plotGrid(G(i),find(~actnum(G(i).cells.indexMap)), ...
        'FaceColor','none','EdgeColor',[0.65 0.65 0.65],'EdgeAlpha',0.2);
   ha = plotGrid(G(i),find( actnum(G(i).cells.indexMap)), ...
                'FaceColor','y','EdgeAlpha',0.1);
end
axis off; view(-15,20); camdolly(0,-0.2,0); axis normal

%% Inspect the active model
% To inspect only the active model, we reset the ACTNUM field to its
% original values and recreate the grid. Now, inactive cells will be
% ignored and we therefore get a different unstructured grid.
% If we include the actnum, G from proccessGRDECL has only one component.
grdecl.ACTNUM = actnum; clear actnum;
G = processGRDECL(grdecl);

clf, 
pargs = {'EdgeAlpha'; 0.1; 'EdgeColor'; 'k'};
plotGrid(G, pargs{:}); view(-100,40); axis tight off

%%
% Once the active part of the model is read, we can compute geometric
% primitives like cell and face centroids, cell areas, and face centroids
% and normals by means of the <matlab:doc('computeGeometry')
% computeGeometry> function.

G = computeGeometry(G);

clf
plotCellData(G,G.cells.centroids(:,3), pargs{:});
view(-100,40); axis tight off

%%
% The media (rock) properties can be extracted by means of the
% <matlab:doc('grdecl2Rock') grdecl2Rock> function.  This function inspects
% the keywords present in the input data and constructs a |rock| data
% structure that holds the permeability tensor $K$ (in the structure
% field |perm|), possibly porosity values $\phi$ (in the |poro| field).
%
% The first input argument is the raw input data as constructed by function
% |readGRDECL|.  The second, optional, input argument is a list of active
% cells represented as a mapping from the grid's cells to the global grid
% cells from the original Nx-by-Ny-by-Nz box description. WARNING: the
% |grdecl2Rock| routine returns permeability values in units [md] and to
% use the rock structure for simulation, we therefore need to convert the
% permeabilities to SI units. This was done by the function
% |convertInputUnits| when we first read the data.
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
colorbarHist(grdecl.PORO, [.01 .3],'South',100);

%% 
% Then we show the porosities mapped onto the structural grid
clf
plotCellData(G,rock.poro, pargs{:});
axis tight off; set(gca,'DataAspect',[1 1 0.1]); 
view(-65,55);
colorbarHist(rock.poro, [.01 .3],'South',100);

%%
% show also the net-to-gross
clf
plotCellData(G,rock.ntg, pargs{:});
axis tight off; set(gca,'DataAspect',[1 1 0.1]); 
view(-65,60);
colorbarHist(rock.ntg,[.2 1],'South',100);

%% 
% and the satnum
clf
SN = grdecl.SATNUM(G.cells.indexMap); jj = jet(60);
plotCellData(G,SN, pargs{:});
axis tight off; set(gca,'DataAspect',[1 1 0.1]); 
view(-65,60);
colorbarHist(SN, [.5 6.5], 'South', 0:.5:6.5);
colormap(jj(1:10:end,:))

%% 
% and a plot where we split them up
clf
plotGrid(G,'FaceColor','none', pargs{:});
axis tight off; set(gca,'DataAspect',[1 1 0.1]);
view(-65,60);
caxis([0.5 6.5]), colormap(jj(1:10:end,:))

h1 = plotCellData(G,SN, find(SN==1), pargs{:});
h2 = plotCellData(G,SN, find(SN==5), pargs{:});

%%
delete([h1,h2])
h1 = plotCellData(G,SN, find(SN==2), pargs{:});
h2 = plotCellData(G,SN, find(SN==4), pargs{:});

%%
delete([h1,h2])
h1 = plotCellData(G,SN, find(SN==3), pargs{:});
h2 = plotCellData(G,SN, find(SN==6), pargs{:});

%% Permeability
% The permeability is given as a tridiagonal tensor K = diag(Kx, Ky, Kz)
% and we therefore plot each of the horizontal and vertical components
% using a logarithmic color scale. The vertical permeability has a large
% number of zero values. These are shown as dark blue color, but are not
% included in the histogram
close
for i=[1,3]
    figure;
    p = log10(rock.perm(:,i));
    plotCellData(G,p, pargs{:});
    axis tight off; set(gca,'DataAspect',[1 1 0.1]);
    view(-65,55); 
    % Manipulate the colorbar to get the ticks we want
    cs = [0.001 0.01 0.1 1 10 100 1000 10000];
    caxis(log10([min(cs) max(cs)]*milli*darcy));
    hc = colorbarHist(p(~isinf(p)),caxis,'South',100);
    set(hc, 'YTick', 0.5, 'YTickLabel','mD', ...
        'XTick', log10(cs*milli*darcy), 'XTickLabel', num2str(cs'));
    colormap(jet)
end

%% 
% Likewise, we show a combined histogram of the horizontal/vertical
% permeability (with all zero values excluded).
figure
Kx = convertTo(rock.perm(:,1),milli*darcy);
Kz = convertTo(rock.perm(:,3),milli*darcy);
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
for i=1:6, 
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
plotCellData(G,Mz,Mz<1,pargs{:});
axis tight off; set(gca,'DataAspect',[1 1 0.1]);
view(-65,55);
colorbarHist(Mz(Mz<1),[0 1],'South',100);

%% Copyright notice

displayEndOfDemoMessage(mfilename)

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
