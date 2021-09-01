function [plumes,depths,plist] = getSleipnerPlumeOutlines2019()
% Utility function to read in plume outlines from textfile in Sleipner
% 2019 benchmark model.
%
% SYNOPSIS:
%   [plumes,depths,plist] = getSleipnerPlumeOutlines2019()
% 
% RETURNS:
%   plumes - cell array of size 1 x 9 where each entry corresponds to a
%             different layer and contains a structure with entries:
%    
%             CO2plumeOutline: array with a matrix of x,y,z coords for each
%                                    polyline in that layer
%             polynos: original polyline numbers for points in dataset
%             layer: Layer identifier
%             cp: x,y,z coords of centre of each polyline.
%        
%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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


[plumes,depths] = getPlumeOutlinesAllLayers();


for i = 1:numel(plumes)
    for j = 1:numel(plumes{i}.CO2plumeOutline)        
        coords = [plumes{i}.CO2plumeOutline{j}(:,1) ...
            plumes{i}.CO2plumeOutline{j}(:,2) plumes{i}.CO2plumeOutline{j}(:,3)+depths(i,1)];       
        plumes{i}.cp{j} = mean(coords);
    end
end

% Polylines for each layer, in order of size (largest
% to smallest).

plist{1} = [1];
plist{2} = [1];
plist{3} = [1];
plist{4} = [1,5,6];
plist{5} = [1];
plist{6} = [4,1,5];
plist{7} = [1];
plist{8} = [2,1];
plist{9} = [2,1];




function [plumes,depths] = getPlumeOutlinesAllLayers()
% Read in plume outlines for 2019 Sleipner layered model
datafolder = fullfile(getDatasetPath('sleipner2019'),...
    'sleipner_plumes_boundaries',...
    'Sleipner_Plumes_Boundaries','data')

filelist = dir(fullfile(datafolder,'L*'))
for k = 1 : length(filelist)
  fullFileName = fullfile(datafolder, filelist(k).name);
  fprintf(1, 'Now reading %s\n', fullFileName);
  dat = dlmread(fullFileName,'', 16,0);
  polynos = unique(dat(:,4));
  for j = 1:numel(polynos)
    inds = dat(:,4) == polynos(j);
    plumes{k}.CO2plumeOutline{j}  = dat(inds,1:3);
  end
  plumes{k}.polynos = dat(:,4);
  plumes{k}.layer = filelist(k).name;
end

% Read in layer depths

depths = dlmread('SleipnerLayerDepths.txt','',1,1);

