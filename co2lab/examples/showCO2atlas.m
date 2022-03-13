%% CO2 Storage Atlas: the Norwegian North Sea
% The Norwegian Petroleum Directorate (NPD) has developed an atlas that
% gives an overview over areas in the Norwegian part of the North Sea where
% CO2 can be stored safely in the subsurface for a long time, see
% http://www.npd.no/en/Publications/Reports/CO2-Storage-Atlas-/. As part of
% the atlas, NPD has released geographical data for many of the formations
% that are described in the atlas. The data are given in a GIS formate
% (shape- and rasterfiles) and can be downloaded from their webpage.
%
% One of the purposes of the 'co2lab' in MRST is to provide simple access
% to public datasets. In this example, we describe how you can
% use the functionality in 'co2lab' to download and display the
% data sets that accompany the atlas. For the sake of convenience, we have
% converted the files to an ASCII format more suitable for our
% applications. These files can be inspected using any text editor. We also
% process the files and get both the raw datasets and data structures
% suitable for constructing volumetric corner-point grids of several sand
% bodies. At the end of the example, we analyse the different formations
% and compute the capacity for structural trapping in each formation.

%% Read data
% This example assumes that you have already downloaded the North Sea
% datasets from the NPD webpages. If you have not done so, you can use the
% following command: downloadDataSets('atlas')
mrstModule add co2lab deckformat

fprintf('Loading atlas data (this may take a few minutes)..');
[grdecls, rawdata] = getAtlasGrid(getNorthSeaNames()); %#ok
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

%% Show the data directly: Utsira formation
% The datasets are perfectly usable for visualization on their own. To see
% this, we find the datasets corresponding to the Utsira formation and plot
% both the thickness and the heightmap.
%
% Note that the datasets are not entirely equal: Some sections are not
% included in the thickness map and vice versa. In addition to this, the
% coordinates are not always overlapping, making interpolation neccessary.
%
% Some formations are only provided as thickness maps; These are processed
% by sampling the relevant part of the Jurassic formation for top surface
% structure.

utsira_rd = rawdata(strcmpi(names, 'Utsirafm'));
clf;
for i = 1:numel(utsira_rd)
    urd = utsira_rd{i};
    
    subplot(2,1,i)
    surf(urd.data)
    
    title([urd.name ' ' urd.variant])
    shading interp
    view(0,90)
    axis tight off
end


%% Visualize all the formations
% We then visualize the formations present in the North Sea along with a
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

grdecls = getAtlasGrid(getNorthSeaNames(),'coarsening', 5);
ng = numel(grdecls);

grids = cell(ng,1);
for i = 1:ng
    gd = processGRDECL(grdecls{i});
    grids{i} = computeGeometry(gd(1));
end
clf;
mymap = [
         0         0    1.0000
         0    1.0000         0
    1.0000         0         0
         0         0    0.5000
         0    0.5000         0
    0.7500         0         0         
         0    1.0000    1.0000
    1.0000         0    1.0000
    1.0000    1.0000         0
    0.8500    0.8500    0.8500
    0.4000    0.4000    0.4000    
         0         0         0
    0.2500    0.5000    1.0000
    0.5000    0.2500    0.7500
    0.7500    0.2500    0.5000
    0.7500    0.5000    0.2500
    0.2500    1.0000    0.5000
    0.5000    1.0000    0.2500
         ];
hold on

for i=1:ng;
    G = grids{i};
    bf = boundaryFaces(G);
    ind = G.faces.normals(bf,3)==0;
    plotFaces(G,bf(ind), 'FaceColor', 'none', 'EdgeColor', mymap(i,:), 'LineWidth',2);
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

ax = axis();
colormap hsv

% Load and plot a map 
%load(fullfile(mrstPath('co2lab'), 'data', 'atlas', 'norway.mat'));
load(fullfile(getDatasetPath('co2atlas'), 'norway.mat'));
for k=1:length(norway), 
    line(norway{k}(:,1) + 6.8e5, norway{k}(:,2)); 
end;
axis(ax)

hold on
%load(fullfile(mrstPath('co2lab'), 'data', 'atlas', 'welldata.mat'));
load(fullfile(getDatasetPath('co2atlas'), 'welldata.mat'));
plot(welldata(:,2), welldata(:,1), '.k', 'MarkerSize', 5)


%% Redo the visualization in 3D with higher resolution
% To show the depth of the various North Sea formations, we redo the plot
% in 3D.
figure;
p = get(gcf,'Position'); set(gcf,'Position', p + [-300 0 300 0]);
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
legend(cellfun(@(x) x.name, grdecls, 'UniformOutput', false), ...
   'Location', 'EastOutside')

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
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
