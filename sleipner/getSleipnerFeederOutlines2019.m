function [feeders] = getSleipnerFeederOutlines2019()
% Utility function to read in feeder outlines from textfile in Sleipner
% 2019 benchmark model.
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
datafolder = fullfile(getDatasetPath('sleipner2019'),...
    'feeders','data');

filelist = dir(fullfile(datafolder,'*_feeder_*'));

% Layers that each feeder is present in according to 
% Sleipner 2019 Reference Model - Polygons_feeders_chimney.pptx
% Add layer 9 as there are 10 layers in our model including the thick shale
% layer.
layernos{1} = 1:1:9;
layernos{2} = 5;
layernos{3} = 7;

feeders = cell(length(filelist),1);

for k = 1 : length(filelist)
  fullFileName = fullfile(datafolder, filelist(k).name);
  dat = dlmread(fullFileName,'', 16,0);
  polynos = unique(dat(:,4));
  for j = 1:numel(polynos)
    inds = dat(:,4) == polynos(j);
    feeders{k}.outline{j}  = dat(inds,1:3);
  end
  feeders{k}.layernos = layernos{k};
end


