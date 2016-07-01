%% Introduction to the Norne Model
% Norne is an oil and gas field lies located in the Norwegian Sea. The
% reservoir is found in Jurrasic sandstone at a depth of 2500 meter below
% sea level. Operator Statoil and partners (ENI and Petoro) have agreed
% with NTNU to release large amounts of subsurface data from the Norne
% field for research and education purposes.  TheH
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
% For simplicity, we assume that all the pertinent sections of the model
% have been downloaded and concatenated into a single input model. To get
% the grid, it is sufficient to concatenate the two files
% |IRAP_1005.GRDECL| and |ACTNUM_0704.prop| into a single file, which we
% here have called |OPM.GRDECL|. Refer to |showSAIGUP| for a more
% in-depth discussion of how to read and process such input data.
grdecl = fullfile(getDatasetPath('norne'), 'OPM.GRDECL');
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
   'FaceColor','red','FaceAlpha',0.2, 'EdgeColor','r','EdgeAlpha',0.1);
axis off; view(-155,80); zoom(1.7);

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

figure
plotCellData(g,g.cells.volumes, pargs{:});
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


%% Extract the active part of the model
% To get the whole grid, we need to reprocess the input data. This gives a
% grid with two components: the main reservoir and a detached stack of
% twelve cells, which we ignore henceforth.
G  = processGRDECL(grdecl);
clf
plotGrid(G(1),pargs{:});
plotGrid(G(2),pargs{:},'FaceColor','r');
view(-110,50), axis tight off, set(gca,'DataAspect',[20 11 1]); zoom(1.7);

%% Outline petrophysical properties
% To also include the petrophysical data, we need to either also
% concatenate the files |NTG_0704.prop|, |PERM_0704.prop|, and
% |PORO_0704.prop| into the same input file as the grid specification, or
% alternatively work with the full data file and remove the parts that we
% are not interested in. Here, we assume that this has been done so that
% the necessary data has been read into the |grdecl| structure.  We can
% then extract the petrophysical properties. Notice that here, the rock 
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

%%
% show also the net-to-gross
clf
plotCellData(G,rock.ntg, pargs{:});
axis tight off; set(gca,'DataAspect',[1 1 0.1]); 
view(-110,55); zoom(1.8);
colorbarHist(rock.ntg,[.05 1],'South',100);


%%
% Permeability
% The permeability is generally a tridiagonal tensor K = diag(Kx, Ky, Kz).
% For the Norne model, the data file only specifies Kx, and then various
% manipulations are done to get the correct vertical permeability. These
% are not incorporated herein, and hence we only show the horizontal
% component of the permeability.
clf
p = log10(rock.perm(:,1));
plotCellData(G,p, pargs{:});
axis tight off; set(gca,'DataAspect',[1 1 0.1]);
view(-110,55); zoom(1.8);
    
% Manipulate the colorbar to get the ticks we want
cs = [1 10 100 1000 10000];
caxis(log10([min(cs) max(cs)]*milli*darcy));
hc = colorbarHist(p(~isinf(p)),caxis,'South',100);
set(hc, 'YTick', 0.5, 'YTickLabel','mD', ...
    'XTick', log10(cs*milli*darcy), 'XTickLabel', num2str(cs'));

%%
displayEndOfDemoMessage(mfilename)

% #COPYRIGHT_EXAMPLE#
