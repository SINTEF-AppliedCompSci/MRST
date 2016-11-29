%% CO2 Storage Atlas for the Norwegian North Sea
% The CO2 lab has functionality to construct various sand bodies based
% on public dataset supplied by the Norwegian Petroleum Directorate as part
% of the recent CO2 storage atlas,
% http://www.npd.no/en/Publications/Reports/CO2-Storage-Atlas-/.
% In this example, we will show the corresponding formations on a map,
% visualize the sand bodies, and compute the potential for structural
% trapping.
mrstModule add co2lab;

%% Visualize all formations
% In this example, we visualizing all formations that can be constructed
% along with a map of Norway and point plots of all production wells in the
% Norwegian continental shelf.
%
% The well data comes from the Norwegian Petroleum Directorate and can be
% found in more detail at http://factpages.npd.no/factpages/.
%
% The map of Norway comes from The Norwegian Mapping and Cadastre Authority
% and can be found at http://www.kartverket.no/. Note that the map is only
% provided for scale and rough positioning - no claims are made regarding
% the accuracy in relation the subsea reservoirs.
 
% Again we read all the fine grids and process them to get proper grids
% with coordinates. The coarsening is useful here as we will only need the
% basic outline for the upcoming figure.
grdecls = getAtlasGrid(getNorthSeaNames(),'coarsening', 5);
ng = numel(grdecls);

grids = cell(ng,1);
for i = 1:ng
    gd = processGRDECL(grdecls{i});
    grids{i} = computeGeometry(gd(1));
end

clf;
hold on

for i=1:ng;
    G = grids{i};
    % We want to colorize each grid differently
    data = repmat(i, G.cells.num, 1);
    plotCellData(grids{i}, data, 'facea', .5, 'edgea', .025, 'edgec', 'k');
end

legend(cellfun(@(x) x.name, grdecls, 'UniformOutput', false), 'Location', 'EastOutside')

box on
view(2)
set(gcf,'Color',[0.8 0.8 0.8]);set(gca,'Color',[0.8 0.8 0.8]);
set(gca,'XColor',[0,0,0])
set(gca,'YColor',[0,0,0]);set(gca,'LineWidth',1)

ax = axis();
colormap(colorcube(ng));

%
% Load and plot a map 
load(fullfile(getDatasetPath('co2atlas'), 'norway.mat'));
for k=1:length(norway), 
    line(norway{k}(:,1) + 6.8e5, norway{k}(:,2)); 
end;
axis(ax)

hold on
load(fullfile(getDatasetPath('co2atlas'), 'welldata.mat'));
plot(welldata(:,2), welldata(:,1), '.k', 'MarkerSize', 5)
hold off


%% Redo the visualization in 3D with higher resolution
% To show the depth of the various formations, we redo the plot in 3D.
grdecls = getAtlasGrid(getNorthSeaNames(),'coarsening', 3);
ng = numel(grdecls);

grids = cell(ng,1);
for i = 1:ng
    gd = processGRDECL(grdecls{i});
    grids{i} = computeGeometry(gd(1));
end

clf;
hold on
for i=1:ng;
    G = grids{i};
    % We want to colorize each grid differently
    data = repmat(i, G.cells.num, 1);
    plotCellData(grids{i}, data, 'facea', .8, 'edgea', .01, 'edgec', 'k');
end
axis tight; box on; view(-300,50)
hold off
%legend(cellfun(@(x) x.name, grdecls, 'UniformOutput', false), 'Location', 'EastOutside')


%% Full 3D visualization of all grids
% Not all formations in the data set supply both a height map of the top
% surface and a map of the formation thickness, which are both needed to
% construct a volumetric sandbody. Next, we visualize the fourteen
% different sandbodies that can be constructed based on pairs of depth and
% thickness maps. 
moduleCheck('libgeometry','opm_gridprocessing');
viewMat = [-58 50; -90 65; -120 30; -120 30; -110 60; -110 60; ...
   -105 55; -105 -10; -95 40; -40 15; -90 80; -55 80; ...
   105 55; 60 20; -100 55; 90 45];
for i=1:ng;
   grdecl = getAtlasGrid(grdecls{i}.name, 'coarsening', 1);
   try
      G = processgrid(grdecl{1});
      G = mcomputeGeometry(G(1));
   catch
      G = processGRDECL(grdecl{1});
      G = computeGeometry(G(1));
   end
   clf;
   plotGrid(G,'FaceColor', 'y', 'EdgeAlpha', .05);
   view(viewMat(i,:)); axis tight off
   light('Position',[-1 -1 -1],'Style','infinite');lighting phong
   
   drawnow
   %print('-dpng', ['grids/' grdecls{i}.name '.png']);
end

%% Basic capacity estimates for the CO2 atlas data
% Finally, we compute the potential volumes available for structural
% trapping. To better report progress, we first load a low resolution
% version to get names of all aquifers. Then we load and process the
% full-resolution versions using both the cell-based and node-based
% methods.
moduleCheck('coarsegrid','matlab_bgl');
grdecls = getAtlasGrid(getNorthSeaNames(),'coarsening',10);
ng = numel(grdecls);
res = cell(ng,1);
for i=1:ng
   fprintf('------------------------------------------------\n');
   fprintf('Processing %s ....\n', grdecls{i}.name);
   grdecl  = getAtlasGrid(grdecls{i}.name, 'coarsening', 1);
   try
      G = mprocessGRDECL(grdecl{1});
      G = mcomputeGeometry(G(1));
   catch
      G = processGRDECL(grdecl{1});
      G = computeGeometry(G(1));
   end
   Gt  = topSurfaceGrid(G);
   tan = trapAnalysis(Gt, false);
   tac = trapAnalysis(Gt, true);
   
   res{i}.name      = grdecls{i}.name;
   res{i}.cells     = Gt.cells.num;
   res{i}.zmin      = min(Gt.cells.z);
   res{i}.zmax      = max(Gt.cells.z);
   res{i}.volume    = sum(G.cells.volumes);
   res{i}.ctrapvols = volumesOfTraps(Gt,tac, []);
   res{i}.ccapacity = sum(res{i}.ctrapvols);
   res{i}.ntrapvols = volumesOfTraps(Gt,tan, []);
   res{i}.ncapacity = sum(res{i}.ntrapvols);
   fprintf('done\n');
end

%% Show table of volumes
fprintf('\n\n%-13s& Cells  & Min  & Max  & Volume  & Capacity & Percent &  Capacity& Percent\\\\\n', 'Name');
for i=1:ng
   fprintf('%-13s& %6d & %4.0f & %4.0f & %4.2e &  %4.2e & %5.2f &   %4.2e & %5.2f \\\\\n',...
      res{i}.name, res{i}.cells, res{i}.zmin, res{i}.zmax, res{i}.volume, ...
      res{i}.ncapacity, res{i}.ncapacity/res{i}.volume*100, ...
      res{i}.ccapacity, res{i}.ccapacity/res{i}.volume*100);
end


