function downloadDataSets(varargin)
% Download data sets for CO2 laboratory
%
% SYNOPSIS:
%   downloadDataSets()
%   downloadDataSets('all')
%   downloadDataSets({'atlas', 'johansen', ..})
%
% DESCRIPTION:
%   This script contains functionality to download and install all data
%   sets used by the CO2 laboratory that are publicly available and do not
%   require registration by the user:
%
%   * Name: 'atlas'
%     The Norwegian Petroleum Directorate (NPD) has developed an atlas that
%     gives an overview over areas in the Norwegian part of the North Sea
%     where CO2 can be stored safely in the subsurface for a long time, see
%     http://www.npd.no/en/Publications/Reports/CO2-Storage-Atlas-/. As
%     part of the atlas, NPD has released geographical data for many of the
%     formations that are described in the atlas. The data are given in a
%     GIS formate (shape- and rasterfiles)
%     Size of data set: ~14 MB
%
%   * Name: 'igems'
%     The IGEMS project (http://www.nr.no/nb/IGEMS) studied how top surface
%     morphology influences the CO2 storage capacity. Alternative
%     top-surface morphologies are created stochastically by combining
%     different stratigraphic scenarios with different structural
%     scenarios. For more information about the data set, see
%     http://files.nr.no/igems/data.pdf.
%     Size of data set: ~1.1 GB
%
%   * Name: 'johansen'
%     The Johansen formation is a candidate site for large-scale CO2
%     storage offshore the south-west coast of Norway. The Johansen data
%     set, consists of four different sector models and a full-field model.
%     Here, we only download the NPD5 sector and the full-field models. For
%     more details, see:
%     http://www.sintef.no/Projectweb/MatMorA/Downloads/Johansen/
%     Size of data set: ~13 MB
%
%   * Name: 'slopingaquifer'
%     A conceptual model from the IGEMS study. This synthetic model was
%     created early in the IGEMS project (http://www.nr.no/igems) as an
%     example of a large 3D grid model with variations in the topsurface
%     topography.
%     Size of data set: ~3 MB
%
% If the function is called with no argument, it goes through the list and
% asks the user for permission to download each dataset. When given a
% specific name or 'all', the function downloads the datasets without
% asking for permission.
%
% In addition to the above, a few examples use the Sleipner data set from
% ieaghg.org, which requires user registration and hence must be downloaded
% manually; see "sleipner/README" for more details.

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

if nargin > 0
    ask = false;
    given = @(x) any(strcmpi(varargin{1}, x));
    if given('all')
        getJohansen = true;
        getSloping = true;
        getAtlas = true;
        getIGEMS = true;
    else
        getJohansen = given('johansen');
        getSloping = given('slopingaquifer');
        getAtlas = given('atlas');
        getIGEMS = given('igems');
    end
else
    ask = true;
    getJohansen = false;
    getSloping = false;
    getAtlas = false;
    getIGEMS = false;
end

% Johansen model
if getJohansen || (ask && userConsent('Download the Johansen data set (~13 MB)'))
   fprintf(1, '  Downloading Johansen sector model (NPD5) ..');
   unzip('https://www.sintef.no/project/MatMoRA/Johansen/NPD5.zip', ...
      fullfile(mrstPath('co2lab'), 'data', 'johansen'));
   fprintf(1,'done\n');

   fprintf(1,'  Downloading Johansen full-field model..');
   unzip('https://www.sintef.no/project/MatMoRA/Johansen/FULLFIELD_Eclipse.zip', ...
      fullfile(mrstPath('co2lab'), 'data', 'johansen'));
    fprintf(1,'done\n');  
end

% IGEMS conceptual model
if getSloping || (ask && userConsent('Download synthetic sloping aquifer model (~3 MB)'))
   fprintf(1,'  Downloading synthetic sloping aquifer model ..');
   untar('https://www.sintef.no/project/MRST/IGEMS.tar.gz',...
      fullfile(mrstPath('co2lab'), 'data'));
   fprintf(1,'done\n');
end

% CO2 Atlas
if getAtlas || (ask && userConsent('Download CO2 Atlas (~14 MB)'))
   fprintf(1,'  Downloading CO2 Atlas (this may take a few minutes) ..');
   untar('https://www.sintef.no/project/MRST/CO2_atlas.tar.gz',...
      fullfile(mrstPath('co2lab'), 'data'));
   untar('https://www.sintef.no/project/MRST/mapAndWells.tar.gz', ...
      fullfile(mrstPath('co2lab'), 'data', 'atlas'));
   fprintf(1,'done\n');
end

% IGEMS data set
if getIGEMS || (ask && userConsent('Download IGEMS data set (~1.1 GB)'))
   % Download and check surfaces
   fprintf(1,'  Downloading IGEMS surfaces (this may take several minutes) ..');
   datadir = fullfile(mrstPath('co2lab'), 'data', 'igems');
   unzip('https://files.nr.no/igems/surfaces.zip', datadir);
   dirs = dir(fullfile(datadir,'surfaces'));
   for i=1:numel(dirs)
      if sum(isletter(dirs(i).name))==0, continue, end
      if numel(dir(fullfile(datadir,'surfaces',dirs(i).name,'*.irap')))<100
         disp(['Installation of surfaces/' dirs(i).name ' seems wrong']);
      end
   end
   fprintf(1,'done\n');
   
   % Download and check eclipse data sets
   fprintf(1, '  Downloading IGEMS grids (this may take several minutes) ..');
   unzip('https://files.nr.no/igems/eclipsegrids.zip', datadir);
   if numel(dir(fullfile(datadir,'eclipsegrids','*.GRDECL')))~=15
      disp('Installation of eclipsegrids/ seems wrong');
   end
   fprintf(1,'done\n');
end
