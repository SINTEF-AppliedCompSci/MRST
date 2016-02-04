%% CO2 Storage Atlas: the Barents Sea
% The Norwegian Petroleum Directorate (NPD) has developed an atlas that
% gives an overview over areas in the Norwegian Contential Shelf (NCS)
% where CO2 can be stored safely in the subsurface for a long time, see
% http://www.npd.no/en/Publications/Reports/CO2-Storage-Atlas-/.
%
% While the NPD has released geographical data for the formations in the
% North Sea, the Barents Sea formation datasets were obtained via email,
% and are not currently downloaded on their website.
%
% The following analyzes the Barents Sea datasets, in particular the
% Hammerfest Aquifer, which is made up of three stacked formations, and
% where the Snohvit injection site is located.

%% Read data
% This example assumes that you have already obtained the Barents Sea
% datasets, and have been placed in the same directory as the North Sea
% datasets (which can be downloaded using downloadDataSets('atlas')).

try
   require co2lab
catch %#ok<CTCH>
   mrstModule add co2lab
end

moduleCheck deckformat

fprintf('Loading atlas data (this may take a few minutes)..\n');
[grdecls, rawdata] = getAtlasGrid( [getBarentsSeaNames() getNorwegianSeaNames() getNorthSeaNames()] ); %#ok
fprintf('done\n');

%% Description of raw data
% Show the raw data. Each dataset contains four fields:
% - Name, which is the name of the formation
% - Variant: Either thickness or height, indicating wether the dataset
% represents height data or thickness data.
% - Data: The actual datasets as a matrix.
% - Meta: Metadata. The most interesting field here is the
%   xllcorner/yllcorner variable which indicates the position in ED50
%   datum space.

fprintf('\nRaw data:\n')
fprintf('----------------------------------------------------------------\n');
for i=1:numel(rawdata);
    rd = rawdata{i};
    fprintf('Dataset %-2i is %-12s (%-9s). Resolution: %4i meters\n', ...
            i, rd.name, rd.variant,  rd.meta.cellsize)
end
fprintf('----------------------------------------------------------------\n');

% Store names for convenience
names = cellfun(@(x) x.name, rawdata, 'UniformOutput', false)';

%% Show the data directly: Bjarmeland formation
% The datasets are perfectly usable for visualization on their own. To see
% this, we find the datasets corresponding to the Bjarmeland formation and
% plot both the thickness and the heightmap.
%
% Note that the resolution of the thickness map is different from the
% resolution of the height map, in the case of Bjarmeland. These maps are
% processed using interpolation in order to create the appropriate model
% grids.
%
% % Note that the datasets are not entirely equal: Some sections are not
% % included in the thickness map and vice versa. In addition to this, the
% % coordinates are not always overlapping, making interpolation neccessary.
% %
% % Some formations are only provided as thickness maps; These are processed
% % by sampling the relevant part of the Jurassic formation for top surface
% % structure.

bjarmeland_rd = rawdata(strcmpi(names, 'Bjarmelandfm'));
clf;
for i = 1:numel(bjarmeland_rd)
    rd = bjarmeland_rd{i};
    
    subplot(2,1,i)
    surf(rd.data)
    
    title([rd.name ' ' rd.variant])
    shading interp
    view(0,90)
    axis tight off
end


%% Visualize all the formations
% We then visualize the formations present in the Barnets Sea along with a
% map of Norway and point plots of all production wells in the Norwegian
% Continental Shelf.
%
% The well data comes from the Norwegian Petroleum Directorate and can be
% found in more detail at http://factpages.npd.no/factpages/.
%
% The map of Norway comes from The Norwegian Mapping and Cadastre Authority
% and can be found at http://www.kartverket.no/. Note that the map is only
% provided for scale and rough positioning - no claims are made regarding
% the accuracy in relation the subsea reservoirs.
%
% To visualize the formations, we load a 5x5 coarsened version of each data
% set and use this to create a simple volumetric model that represents
% approximately the outline of each formation. More details about how to
% create 3D grid models are given in the script 'modelsFromAtlas.m'

names = [getBarentsSeaNames() getNorwegianSeaNames() getNorthSeaNames()];

% Remove certain formation names:
names = names(~strcmpi(names,'Nordmelafm'));
names = names(~strcmpi(names,'Rorfm'));
names = names(~strcmpi(names,'Notfm'));
names = names(~strcmpi(names,'Knurrfm'));       % @@ can be constructed??
names = names(~strcmpi(names,'Fruholmenfm'));   % @@
names = names(~strcmpi(names,'Cookfm'));
names = names(~strcmpi(names,'Dunlingp'));
names = names(~strcmpi(names,'Paleocene'));

grdecls = getAtlasGrid( names ,'coarsening', 5);
ng = numel(grdecls);

grids = cell(ng,1);
for i = 1:ng
    gd = processGRDECL(grdecls{i});
    grids{i} = computeGeometry(gd(1));
end
clf;

% Create color map with distinguishable colors for each grid
assert(exist('colorspace')~=0, 'Ensure colorspace exists and is on path.')
func = @(x) colorspace('RGB->Lab',x);
mymap = distinguishable_colors(ng,'w',func);
% mymap = [
%          0         0    1.0000
%          0    1.0000         0
%     1.0000         0         0
%          0         0    0.5000
%          0    0.5000         0
%     0.7500         0         0         
%          0    1.0000    1.0000
%     1.0000         0    1.0000
%     1.0000    1.0000         0
%     0.8500    0.8500    0.8500
%     0.4000    0.4000    0.4000    
%          0         0         0
%     0.2500    0.5000    1.0000
%     0.5000    0.2500    0.7500
%     0.7500    0.2500    0.5000
%     0.7500    0.5000    0.2500
%     0.2500    1.0000    0.5000
%     0.5000    1.0000    0.2500
%          ];
hold on

for i=1:ng;
    G = grids{i};
    bf = boundaryFaces(G);
    ind = G.faces.normals(bf,3)==0;
    plotFaces(G, bf(ind), 'FaceColor', 'none', 'EdgeColor', mymap(i,:), 'LineWidth',2);
%    % We want to colorize each grid differently
%    data = repmat(i, G.cells.num, 1);
%    plotCellData(grids{i}, data, 'facea', .3, 'edgea', .05, 'edgec', 'k');
end

legend(cellfun(@(x) x.name, grdecls, 'UniformOutput', false), ...
   'Location', 'EastOutside')
box on
view(2)
set(gcf,'Color',[0.8 0.8 0.8]);set(gca,'Color',[0.8 0.8 0.8]);
set(gca,'XColor',[0,0,0])
set(gca,'YColor',[0,0,0]);set(gca,'LineWidth',1)

%ax = axis();
colormap hsv

% Load and plot a map 
%load(fullfile(mrstPath('co2lab'), 'data', 'atlas', 'norway.mat'));
load(fullfile(getDatasetPath('co2atlas'), 'norway.mat'));
for k=1:length(norway), 
    line(0.85*(norway{k}(:,1)) + 0.37e6, norway{k}(:,2)); 
end;
%axis(ax)
axis equal
axis tight

% Note: well data coordinates do not line up with map of Norway.
% hold on
% %load(fullfile(mrstPath('co2lab'), 'data', 'atlas', 'welldata.mat'));
% load(fullfile(getDatasetPath('co2atlas'), 'welldata.mat'));
% plot(welldata(:,2), welldata(:,1), '.k', 'MarkerSize', 5)


%% Redo the visualization in 3D with higher resolution
% To show the depth of the various Barents Sea formations, we redo the plot
% in 3D.
figure;
p = get(gcf,'Position'); set(gcf,'Position', p + [-300 0 300 0]);
grdecls = getAtlasGrid( [getBarentsSeaNames()] ,'coarsening', 3);
%grdecls = getAtlasGrid( {'Stofm','Nordmelafm','Tubaenfm'} ,'coarsening', 3);
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
legend(cellfun(@(x) x.name, grdecls, 'UniformOutput', false), ...
   'Location', 'EastOutside')
